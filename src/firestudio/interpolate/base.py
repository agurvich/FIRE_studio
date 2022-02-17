import os 
import copy
import itertools
import multiprocessing
import time

import numpy as np

from abg_python.plot_utils import plt,ffmpeg_frames
from abg_python.galaxy.gal_utils import ManyGalaxy
from abg_python.galaxy.gal_utils import Galaxy

from ..studios.gas_studio import GasStudio
from ..studios.star_studio import StarStudio 
from ..studios.FIRE_studio import FIREStudio
from ..studios.simple_studio import SimpleStudio
from ..studios.composition import Composition

from .time_helper import single_threaded_control_flow

class BaseInterpolate(object):
    def interpolateAndRender(
        self,
        galaxy_kwargs, ## only 1 dict, shared by all frames
        scene_kwargss=None, ## 1 dict per frame
        studio_kwargss=None,
        render_kwargss=None, ## only 1 dict, shared by all frames
        which_studios=None,
        keyframes=False,
        multi_threads=1,
        timestamp=0, ## offset by 0 Myr, pass None for no timestamp
        check_exists=True,
        add_composition=False
        ):
        """ """

        ## should we rotate to match the final orientation?
        if 'final_orientation' in galaxy_kwargs:
            orientation_flag = galaxy_kwargs.pop('final_orientation')

        ## create a many galaxy instance to get the datadir
        ##  and potentially the final orientation of the simulation
        many_galaxy = ManyGalaxy(
            galaxy_kwargs['name'],
            suite_name=galaxy_kwargs['suite_name'] if 'suite_name' in galaxy_kwargs else 'metal_diffusion')

        ## get the face-on rotation angles for the z=0 snapshot
        ##  and render the movie in that static frame
        ##  otherwise we'll use the face-on frame at each redshift
        ##  (which may introduce weird artifacts at early times)
        if orientation_flag:
            theta,phi = many_galaxy.get_final_orientation()
            galaxy_kwargs['force_theta_TB'] = theta
            galaxy_kwargs['force_phi_TB'] = phi

        ## take input dictionaries and make sure input is valid
        ##  checks:
        ## *1: whether FIREstudio studios have been passed in which_studios
        ## *2: whether arrays are the correct length i.e.
        ##    len(which_studios) == len(studio_kwargss) == len(render_kwargss)
        ## *3: masks which frames need to be rendered if we've asked only 
        ##    keyframes to be rendered or if we're checking if frames already exist
        ## *4: adds an extra "Composition" studio which is a tiled mosaic of all the
        ##    studios in which_studios.
        ## *5: names frames if savefig is not present in studio_kwargss[i]. defaults to 
        ##    outputting a frame for each studio. this can be overridden by setting
        ##    studio_kwargs['savefig'] = None
        ## *6: adds '_frame%d.png' as savefig suffixes, this is what numbers each frame
        (which_studios,
        galaxy_kwargs,
        scene_kwargss,
        studio_kwargss,
        render_kwargss,) = self.sanitize_input(
            which_studios,
            galaxy_kwargs,
            scene_kwargss,
            studio_kwargss,
            render_kwargss,
            keyframes,
            many_galaxy.datadir if check_exists else None,
            add_composition=add_composition)
        
        ## check if there's actually any work to be done
        if len(scene_kwargss) == 0:
            print('all frames already rendered, exiting...')
            render_ffmpeg_frames(studio_kwargss,galaxy_kwargs,self.nframes,self.fps)
            return [None]

        if 'snap_pair' in scene_kwargss[0]:
            snap_pairs = [scene_kwargs['snap_pair'] for scene_kwargs in scene_kwargss]
        else: snap_pairs = None

        ## single threaded, the vanilla experience
        if multi_threads == 1:
            ## collect positional arguments for worker_function
            return_value = single_threaded_control_flow(
                which_studios,
                galaxy_kwargs, ## 1 dictionary
                scene_kwargss, ## Ntimesteps many dictionaries, shared for all studios
                studio_kwargss, ## Nstudios
                render_kwargss, ## Nstudios
                many_galaxy.datadir,
                timestamp,
                self.coord_interp_mode)
            
        ## multi-threaded, sometimes can have some h5py conflicts
        ##  so definitely not *completely* thread-safe but usually 
        ##  one can get around this by just resubmitting the job
        elif multi_threads > 1:
            ## split the pairs of snapshots into approximately equal chunks
            ##  prioritizing  matching pairs of snapshots
            mps_indices = split_into_n_approx_equal_chunks(snap_pairs,multi_threads)

            scene_kwargss = np.array_split(scene_kwargss,mps_indices)
            
            ## collect positional arguments for worker_function
            argss = zip(
                itertools.repeat(which_studios), ## Nstudios
                itertools.repeat(galaxy_kwargs), ## 1 dictionary
                scene_kwargss, ## Ntimesteps
                itertools.repeat(studio_kwargss), ## Nstudios
                itertools.repeat(render_kwargss), ## Nstudios
                itertools.repeat(many_galaxy.datadir),
                itertools.repeat(timestamp),
                itertools.repeat(self.coord_interp_mode))

            with multiprocessing.Pool(multi_threads) as my_pool:
                return_value = my_pool.starmap(single_threaded_control_flow,argss)
            return_value = np.hstack(return_value)
        else:
            raise ValueError("Specify a number of threads >=1, not",multi_threads)

        ## use ffmpeg to produce an mp4 of the frames
        render_ffmpeg_frames(studio_kwargss,galaxy_kwargs,self.nframes,self.fps)

        return return_value

    def sanitize_input(
        self,
        which_studios,
        galaxy_kwargs,
        scene_kwargss,
        studio_kwargss,
        render_kwargss,
        keyframes,
        datadir,
        add_composition=False):
        """ 
        ## take input dictionaries and make sure input is valid
        ##  checks:
        ## *1: whether FIREstudio studios have been passed in which_studios
        ## *2: whether arrays are the correct length i.e.
        ##    len(which_studios) == len(studio_kwargss) == len(render_kwargss)
        ## *3: masks which frames need to be rendered if we've asked only 
        ##    keyframes to be rendered or if we're checking if frames already exist
        ## *4: adds an extra "Composition" studio which is a tiled mosaic of all the
        ##    studios in which_studios.
        ## *5: names frames if savefig is not present in studio_kwargss[i]. defaults to 
        ##    outputting a frame for each studio. this can be overridden by setting
        ##    studio_kwargs['savefig'] = None
        ## *6: adds '_frame%d.png' as savefig suffixes, this is what numbers each frame
        ## *7: adds minimum required snapshot keys to load from disk based on which_studios
        """

        ## in the rare case when 
        if which_studios is None: which_studios = [SimpleStudio]

        if 'keys_to_extract' not in galaxy_kwargs:
            keys_to_extract = []
            for which_studio in which_studios:
                keys_to_extract+=which_studio.required_snapdict_keys
            galaxy_kwargs['keys_to_extract'] = list(np.unique(keys_to_extract))

        ## prepended to frame_%0{log10(N)//1+1}d.png
        for which_studio,studio_kwargs in zip(which_studios,studio_kwargss):
            if 'savefig' not in studio_kwargs: studio_kwargs['savefig'] = which_studio.__name__

        ## handle default arguments
        if studio_kwargss is None: studio_kwargss = [{} for i in range(len(which_studios))]
        if render_kwargss is None: render_kwargss = [{} for i in range(len(which_studios))]
        if scene_kwargss is None: scene_kwargss = self.scene_kwargss

        ## check contents of which_studios
        for which_studio in which_studios:
            if (which_studio is not GasStudio and 
                which_studio is not StarStudio and
                which_studio is not FIREStudio and
                which_studio is not SimpleStudio and
                which_studio is not Composition): 
                raise TypeError("%s is not GasStudio, StarStudio, or FIREStudio"%repr(which_studio))

        ## check lengths of arrays
        nstudios = len(which_studios)
        if len(studio_kwargss) != nstudios: 
            raise ValueError(
                'Mismatched lengths: %d studios - %d studio_kwargs'%(nstudios,len(studio_kwargss)))

        if len(render_kwargss) != nstudios: 
            raise ValueError(
                'Mismatched lengths: %d studios - %d render_kwargs'%(nstudios,len(render_kwargss)))


        ## if we don't have a dictionary for each requested frame we'll just repeat the last frame
        ##  until we hit the requested duration
        ndiff =  self.nframes - len(scene_kwargss)
        ## repeat the last frame until we reach the total duration \_(ツ)_/
        scene_kwargss = np.array(scene_kwargss + [copy.copy(scene_kwargss[-1]) for i in range(ndiff)],dtype=object)
        if ndiff >0:
            print('Repeating the last frame %d times because %d frames were requested.'%(ndiff,self.nframes))

        ## add a Composition studio (which is a tiled mosaic of all the studios in which_studios)
        ##  if requested
        if add_composition:
            if type(add_composition) != str:
                c_savefig = '_'.join([which_studio.__name__ for which_studio in which_studios])
            else:
                c_savefig=add_composition
            c_studio_kwargs = {
                'studios_tuple':tuple(which_studios),
                'subplots_kwargs':{
                    'wspace':0,'hspace':0,
                    'left':0,'right':1,
                    'bottom':0,'top':1},
                'savefig':c_savefig,
                'studio_kwargss':copy.deepcopy(studio_kwargss),
                'ncols':min(len(which_studios),4)
                #'size_inches':(12,6),
                }

            ## tell the interpolator to initialize a composition
            which_studios = which_studios+[Composition]
            ## add the kwargs that should be passed to the Composition at initialization
            studio_kwargss = studio_kwargss+[c_studio_kwargs]
            ## join the render kwargs, they'll be ignored by the studios that don't need them
            c_render_kwargs = {}
            for this_render_kwargs in render_kwargss: c_render_kwargs.update(this_render_kwargs)

            ## compositions should only be run immediately after their
            ##  components
            c_render_kwargs['use_metadata'] = True
            c_render_kwargs['assert_cached'] = True
            render_kwargss = render_kwargss + [c_render_kwargs]

        ## mask for only keyframes if requested
        if keyframes: scene_kwargss = scene_kwargss[self.keyframes]

        ## _frame%d to savefig_suffixes to each scene_kwargs dictionary
        ##  this is what gives each frame its number, important this goes
        ##  before the masking below (but after the keyframe masking above)
        scene_kwargss = add_savefig_suffix_to_dicts(scene_kwargss,self.nframes)

        ## if we were passed a directory let's assume we want to check 
        ##  if the frames exist in it
        if datadir is not None:
            frames_to_do = png_frame_cache(
                scene_kwargss,
                studio_kwargss,
                datadir)
        else: frames_to_do = np.ones(scene_kwargss.shape[0],dtype=bool)

        return which_studios,galaxy_kwargs,scene_kwargss[frames_to_do],studio_kwargss,render_kwargss

def add_savefig_suffix_to_dicts(scene_kwargss,nframes):

    format_str = '_frame_%0'+'%dd.png'%(np.ceil(np.log10(nframes)))
    ## initialize array of savefig values
    for i in range(nframes):
        ## determine minimum number of leading zeros
        scene_kwargss[i]['savefig_suffix'] = format_str%i
    return scene_kwargss
        
def png_frame_cache(scene_kwargss,studio_kwargss,datadir):
    ## address png caching here that way we can load balance appropriately
    ##  for multiprocessing
    frames_to_do = []

    for i,scene_kwargs in enumerate(scene_kwargss):
        for studio_kwargs in studio_kwargss:
            if studio_kwargs['savefig'] is None: continue
            else: 
                this_fname = os.path.join(
                    datadir,
                    'firestudio',
                    studio_kwargs['savefig']+scene_kwargs['savefig_suffix'])
            if not os.path.isfile(this_fname): 
                #print(this_fname,'missing')
                frames_to_do.append(i)
                break

    return frames_to_do

def render_ffmpeg_frames(studio_kwargss,galaxy_kwargs,nframes,fps):
    for studio_kwargs in studio_kwargss:
        if 'keys_to_extract' in galaxy_kwargs: galaxy_kwargs.pop('keys_to_extract')

        if 'snapnum' not in galaxy_kwargs: galaxy_kwargs['snapnum'] = None
        galaxy = Galaxy(**galaxy_kwargs)

        this_savefig = studio_kwargs['savefig']
        ## in case we rendered multiple types of frames
        ##  we'll loop through a list of savefigs (even if there's only one)
        if this_savefig is not None:
            format_str = '%s'%this_savefig + '_frame_%0'+'%dd.png'%(np.ceil(np.log10(nframes)))
            ## ffmpeg the frames
            ffmpeg_frames(
                os.path.join(galaxy.datadir,'firestudio'),
                [format_str],
                savename=galaxy_kwargs['name'],
                framerate=fps,
                extension='.mp4')

#### MPS LOAD BALANCING
def split_into_n_approx_equal_chunks(snap_pairs,nchunks):

    
    if snap_pairs is None: mps_chunks = nchunks
    else:
        nrenders = len(snap_pairs)
        ## split matching pairs into groups
        indices = split_pairs(snap_pairs[:])
        splitted = np.array_split(snap_pairs,indices)
        nchunks = min(len(indices),nchunks)
        
        ## determine how many of each matching pair there are
        n_renders_per_split = [len(this_split) for this_split in splitted]
        
        mps_chunks = []
        this_chunk = 0
        # 1 would bias early, do this if the remainder would dump a lot of pairs on the last
        ##  process, i.e. when the remainder > nchunks/2
        per_chunk = nrenders//nchunks + (nrenders%nchunks > nchunks//2)
        
        for i in range(len(indices)+1):
            #print(this_chunk,per_chunk)
            if this_chunk >= per_chunk:
                mps_chunks+=[indices[i-1]]
                this_chunk=0
            this_chunk+=n_renders_per_split[i]
        
        mps_chunks = np.array(mps_chunks)#-indices[0]
        mps_chunks=list(mps_chunks)
        print('split into:','%d threads'%(len(mps_chunks)+1),np.diff([0]+mps_chunks+[len(snap_pairs)]))
        time.sleep(5)
    return mps_chunks

def split_pairs(snap_pairs):
    indices = []
    changed = split_head(snap_pairs)
    cur_index = 0
    while changed > 0:
        changed = split_head(snap_pairs[cur_index:])
        cur_index+=changed
        indices+=[cur_index]

    return indices[:-1]#np.array_split(snap_pairs,indices[:-1])
    
def split_head(snap_pairs):
    split_index = np.argmax(np.logical_not(np.all(snap_pairs==snap_pairs[0],axis=1)))
    return split_index