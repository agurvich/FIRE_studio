import itertools
import multiprocessing
import os

import numpy as np
import matplotlib.pyplot as plt

from abg_python.galaxy.gal_utils import Galaxy

from abg_python.interpolate.time_interpolate_utils import find_bordering_snapnums,single_threaded_control_flow
from .time_helper import split_into_n_approx_equal_chunks

class TimeInterpolationHandler(object):

    def __repr__(self):
        return ('TimeInterpolationHandler with %d snapshot pairs:\n'%len(self.snap_pairs)+
        repr(self.snap_pairs[-10:]).replace('array([[','[ [...],\n[').replace(' ','')[:-1])

    def __init__(
        self,
        dGyr_or_times_gyr,
        snap_times_gyr=None,
        dGyr_tmin=None,
        dGyr_tmax=None):
        """ """

        ## need to check/decide how one would convert to Myr if necessary. I think only
        ##  the snap interpolation would break.

        if snap_times_gyr is None:
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

    def interpolate_on_snap_pairs(
        self,
        function_to_call_on_interpolated_dataframe,
        galaxy_kwargs=None,
        multi_threads=1,
        extra_keys_to_extract=None):
        """ """

        if galaxy_kwargs is None:
            galaxy_kwargs = {}

        if multi_threads == 1:

            return single_threaded_control_flow(
                function_to_call_on_interpolated_dataframe,
                self.times_gyr,
                self.snap_pairs,
                self.snap_pair_times,
                galaxy_kwargs,
                extra_keys_to_extract)
            
        elif multi_threads > 1:
            ## split the pairs of snapshots into approximately equal chunks
            ##  prioritizing  matching pairs of snapshots
            mps_indices = split_into_n_approx_equal_chunks(self.snap_pairs,multi_threads)
            split_times_gyr = np.array_split(self.times_gyr,mps_indices)
            split_snap_pairs = np.array_split(self.snap_pairs,mps_indices)
            split_snap_pair_times = np.array_split(self.snap_pair_times,mps_indices)
            
            argss = zip(
                itertools.repeat(function_to_call_on_interpolated_dataframe),
                split_times_gyr,
                split_snap_pairs,
                split_snap_pair_times,
                itertools.repeat(galaxy_kwargs),
                itertools.repeat(extra_keys_to_extract))

            with multiprocessing.Pool(multi_threads) as my_pool:
                return_value = my_pool.starmap(single_threaded_control_flow,argss)
            return np.hstack(return_value)
        else:
            raise ValueError("Specify a number of threads >=1, not",multi_threads)
        
def main():

    savename = 'm10q_res250'
    snapnum = 600 
    suite_name = 'metal_diffusion'

    ## load up a galaxy object to get its snapnums and snap_gyrs attribtues
    prev_galaxy = Galaxy(
        savename,
        snapnum,
        suite_name=suite_name)

    time_handler = TimeInterpolationHandler(
        prev_galaxy.snapnums,
        dGyr=2,
        snap_times=prev_galaxy.snap_gyrs)

    galaxy_kwargs = {
        'name':savename,
        'suite_name':suite_name}

    time_handler.interpolate_on_snap_pairs(
        my_func,
        galaxy_kwargs=galaxy_kwargs,
        multi_threads=1)

from firestudio.studios.gas_studio import GasStudio
def my_func(interp_snap):
    ## decide what we want to pass to the GasStudio
    render_kwargs = {
        'weight_name':'Masses',
        'quantity_name':'Temperature',
        'min_quantity':2,
        'max_quantity':7,
        'quantity_adjustment_function':np.log10,
        'quick':False,
        'loud':False,
        #'save_meta':False,
        #'use_metadata':False,
        #'min_weight':-0.5,
        #'max_weight':3,
        #'weight_adjustment_function':lambda x: np.log10(x/(30**2/1200**2)) + 10 - 6, ## msun/pc^2,
        }

    my_gasStudio = GasStudio(
        interp_snap['studio_datadir'],
        interp_snap['next_snapnum'], ## attribute this data to the next_snapnum's projection file
        interp_snap['name'],
        gas_snapdict=interp_snap,
        master_loud=False)
    
    ## differentiate this time to << Myr precision
    my_gasStudio.this_setup_id += "_time%.5f"%interp_snap['this_time'] 

    ## create a new figure for this guy
    fig = plt.figure()
    my_gasStudio.render(plt.gca(),**render_kwargs)
    return fig


if __name__ == '__main__':
    main()