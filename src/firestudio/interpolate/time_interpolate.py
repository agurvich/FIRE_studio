import os

import numpy as np

from abg_python.interpolate.time_interpolate_utils import find_bordering_snapnums

from .base import BaseInterpolate

class TimeInterpolationHandler(BaseInterpolate):

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
        self.scene_kwargss = [dict(
            time=self.times_gyr[i],
            snap_pair=self.snap_pairs[i],
            snap_pair_time=self.snap_pair_times[i])
            for i in range(self.nframes)] 