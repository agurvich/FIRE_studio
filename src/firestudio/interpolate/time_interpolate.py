import os
import itertools
import multiprocessing

import numpy as np

from abg_python.interpolate.time_interpolate_utils import find_bordering_snapnums
from abg_python.galaxy.gal_utils import ManyGalaxy

from .time_helper import split_into_n_approx_equal_chunks,single_threaded_control_flow
from ..studios.gas_studio import GasStudio
from ..studios.star_studio import StarStudio
from ..studios.FIRE_studio import FIREStudio
from ..studios.composition import Composition

class TimeInterpolationHandler(object):

    def __repr__(self,verbose=True):
        return ('TimeInterpolationHandler(%d unique snapshots (%d-%d) - %d frames)'%(
            len(self.unique_snaps),self.unique_snaps.min(),self.unique_snaps.max(),len(self.snap_pairs))+
        verbose*('\n'+repr(self.snap_pairs[-10:]).replace('array([[','[ [...],\n[').replace(' ','')[:-1]))

    def __init__(
        self,
        dGyr_or_times_gyr,
        snap_times_gyr=None,
        dGyr_tmin=None,
        dGyr_tmax=None,
        coord_interp_mode='spherical'):
        """ """

        ## need to check/decide how one would convert to Myr if necessary. I think only
        ##  the snap interpolation would break.

        if snap_times_gyr is None:
            print("Using default snapshot_times.txt for 600 snapshots")
            snapnums,scale_factors,redshifts,snap_times_gyr,dgyrs = np.genfromtxt(
                os.path.join(
                    os.path.dirname(__file__),'snapshot_times.txt')).T

        dGyr_or_times_gyr = np.array(dGyr_or_times_gyr)

        if len(dGyr_or_times_gyr.shape) == 0:
            ## we were passed a time spacing
            times_gyr,snap_pairs,snap_pair_times = find_bordering_snapnums(
                    snap_times_gyr,
                    dGyr=dGyr_or_times_gyr,
                    tmin=dGyr_tmin,
                    tmax=dGyr_tmax)

        elif len(dGyr_or_times_gyr.shape) == 1:
            times_gyr,snap_pairs,snap_pair_times = find_bordering_snapnums(
                snap_times_gyr,
                times_gyr=dGyr_or_times_gyr)
        else:
            raise ValueError("Could not interpret dGyr_or_times_gyr",dGyr_or_times_gyr)

        self.times_gyr = times_gyr
        self.snap_pairs = snap_pairs
        self.snap_pair_times = snap_pair_times
        self.unique_snaps = np.unique(self.snap_pairs)
        self.nframes = len(self.times_gyr)

        ## how should we interpolate the coordinates? in spherical coords? cylindrical?
        self.coord_interp_mode=coord_interp_mode

        ## initialize dummy dictionaries
        self.scene_kwargss = [{} for i in range(self.nframes)]


    def interpolateAndRender(
        self,
        galaxy_kwargs, ## only 1 dict, shared by all frames
        scene_kwargss=None, ## 1 dict per frame
        studio_kwargss=None,
        render_kwargss=None, ## only 1 dict, shared by all frames
        savefigs=None,
        which_studios=None,
        multi_threads=1,
        timestamp=0, ## offset by 0 Myr, pass None for no timestamp
        check_exists=True
        ):
        """ """

        if which_studios is None: which_studios = [GasStudio]

        ## handle default arguments
        if studio_kwargss is None: studio_kwargss = [{} for i in range(len(which_studios))]
        if render_kwargss is None: render_kwargss = [{} for i in range(len(which_studios))]
        if scene_kwargss is None: scene_kwargss = self.scene_kwargss

        ## prepended to frame_%0{log10(N)//1+1}d.png
        if savefigs is None: savefigs = [None for i in range(len(which_studios))]  

        for savefig,studio_kwargs in zip(savefigs,studio_kwargss):
            studio_kwargs['savefig'] = savefig

        for i,which_studio in enumerate(which_studios):
            ## determine which studio we should initialize inside the worker_function
            if which_studio is None: which_studios[i] = GasStudio
            elif (which_studio is not GasStudio and 
                which_studio is not StarStudio and
                which_studio is not FIREStudio and
                which_studio is not Composition): 
                raise TypeError("%s is not GasStudio, StarStudio, or FIREStudio"%repr(which_studio))


        format_str = '_frame_%0'+'%dd.png'%(np.ceil(np.log10(self.nframes)))
        ## initialize array of savefig values
        for i in range(self.nframes):
            ## determine minimum number of leading zeros
            scene_kwargss[i]['savefig_suffix'] = format_str%i
        
        times_gyr = self.times_gyr
        snap_pairs = self.snap_pairs
        snap_pair_times = self.snap_pair_times

        if hasattr(self,'keyframes'):
            times_gyr = times_gyr[self.keyframes]
            snap_pairs = snap_pairs[self.keyframes]
            snap_pair_times = snap_pair_times[self.keyframes]
            scene_kwargss = np.array(scene_kwargss)[self.keyframes]

        ## anything to do with the galaxy
        if 'final_orientation' in galaxy_kwargs:
            orientation_flag = galaxy_kwargs.pop('final_orientation')

        ## create a many galaxy instance
        many_galaxy = ManyGalaxy(
            galaxy_kwargs['name'],
            suite_name=galaxy_kwargs['suite_name'] if 'suite_name' in galaxy_kwargs else 'metal_diffusion')

        if orientation_flag:
            theta,phi = many_galaxy.get_final_orientation()
            galaxy_kwargs['force_theta_TB'] = theta
            galaxy_kwargs['force_phi_TB'] = phi

        if check_exists:
            ## address png caching here that way we can load balance appropriately
            ##  for multiprocessing
            frames_to_do = []

            for i,scene_kwargs in enumerate(scene_kwargss):
                for studio_kwargs in studio_kwargss:
                    if studio_kwargs['savefig'] is None:
                        this_fname = None
                    else:
                        this_fname = os.path.join(
                            many_galaxy.datadir,
                            'firestudio',
                            studio_kwargs['savefig']+scene_kwargs['savefig_suffix'])
                    if this_fname is None or not os.path.isfile(this_fname): 
                        frames_to_do.append(i)
                        break

            if len(frames_to_do) == 0: return
            times_gyr = times_gyr[frames_to_do]
            snap_pairs = snap_pairs[frames_to_do]
            snap_pair_times = snap_pair_times[frames_to_do]
            scene_kwargss = np.array(scene_kwargss)[frames_to_do]


        if multi_threads == 1:
            ## collect positional arguments for worker_function
            return single_threaded_control_flow(
                which_studios,
                times_gyr,
                snap_pairs,
                snap_pair_times,
                galaxy_kwargs, ## 1 dictionary
                scene_kwargss, ## Ntimesteps many dictionaries, shared for all studios
                studio_kwargss, ## Nstudios
                render_kwargss, ## Nstudios
                many_galaxy.datadir,
                timestamp,
                self.coord_interp_mode)
            
        elif multi_threads > 1:
            ## split the pairs of snapshots into approximately equal chunks
            ##  prioritizing  matching pairs of snapshots
            mps_indices = split_into_n_approx_equal_chunks(snap_pairs,multi_threads)

            split_times_gyr = np.array_split(times_gyr,mps_indices)
            split_snap_pairs = np.array_split(snap_pairs,mps_indices)
            split_snap_pair_times = np.array_split(snap_pair_times,mps_indices)
            scene_kwargss = np.array_split(scene_kwargss,mps_indices)
            
            ## collect positional arguments for worker_function
            argss = zip(
                itertools.repeat(which_studios), ## Nstudios
                split_times_gyr,
                split_snap_pairs,
                split_snap_pair_times,
                itertools.repeat(galaxy_kwargs),
                scene_kwargss, ## Ntimesteps
                itertools.repeat(studio_kwargss), ## Nstudios
                itertools.repeat(render_kwargss), ## Nstudios
                itertools.repeat(many_galaxy.datadir),
                itertools.repeat(timestamp),
                itertools.repeat(self.coord_interp_mode))

            with multiprocessing.Pool(multi_threads) as my_pool:
                return_value = my_pool.starmap(single_threaded_control_flow,argss)
            return np.hstack(return_value)
        else:
            raise ValueError("Specify a number of threads >=1, not",multi_threads)