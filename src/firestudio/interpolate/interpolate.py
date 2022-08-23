import numpy as np
import os
import copy

from .time_interpolate import TimeInterpolationHandler
from .scene_interpolate import SceneInterpolationHandler
from .base import BaseInterpolate

class InterpolationHandler(BaseInterpolate):
    def __repr__(self):
        return "InterpolationHandler(%s - %s)"%(
            self.time_handler.__repr__(verbose=False),
            self.scene_handler.__repr__())

    def __init__(
        self,
        total_duration_sec,
        sim_time_begin=None,
        sim_time_end=None,
        fps=24,
        snapshot_times=None,
        coord_interp_mode='spherical',
        **scene_kwargs):
        
        ## need to interpolate camera orientation or other scene properties
        ##  the scene handler will have to be called interactively, I think. 
        ##  it gets so complicated trying to add stuff all at the beginning
        self.scene_handler = SceneInterpolationHandler(total_duration_sec,fps,**scene_kwargs)

        ## need to interpolate in time
        if sim_time_begin and sim_time_end is not None:
            self.time_handler = TimeInterpolationHandler(
                np.linspace(sim_time_begin,sim_time_end,int(total_duration_sec*fps)),
                snapshot_times,
                coord_interp_mode=coord_interp_mode)
            self.nframes = self.time_handler.nframes
        else: 
            self.time_handler = None
            self.nframes = int(fps*total_duration_sec)
        
        self.scene_handler.nframes = self.nframes
        self.fps = fps

    def interpolateAndRender(
        self,
        galaxy_kwargs, ## only 1 dict, shared by all frames
        studio_kwargss=None, ## only 1 dicts, shared by all frames
        render_kwargss=None, ## only 1 dict, shared by all frames
        which_studios=None,
        multi_threads=1,
        keyframes=False,
        check_exists=True,
        timestamp=0, ## offset in Gyr for timestamp, None = no timestamp
        add_composition=False,
        shared_memory=False,
        time_slice=None
        ):

        ## handle simple case of moving camera at fixed time
        if self.time_handler is not None: 
            if keyframes: self.time_handler.keyframes = self.scene_handler.keyframes
            elif hasattr(self.time_handler,'keyframes'): del self.time_handler.keyframes

            ndiff =  self.nframes - len(self.scene_handler.frame_kwargss)
            scene_kwargs = self.scene_handler.frame_kwargss + [
                copy.copy(self.scene_handler.frame_kwargss[-1]) 
                for i in range(ndiff)]

            ## merge dictionaries with priority such that
            ## this_scene_kwargs > this_time_kwargs > studio_kwargs
            scene_kwargss = [{**this_time_kwargs,**this_scene_kwargs} for 
                this_time_kwargs,this_scene_kwargs in 
                zip(self.time_handler.scene_kwargss,scene_kwargs)]
            self.snap_pairs = self.time_handler.snap_pairs
            self.coord_interp_mode = self.time_handler.coord_interp_mode
        else:
            if 'snapnum' not in galaxy_kwargs: raise KeyError("galaxy_kwargs must contain snapnum.")
            scene_kwargss = self.scene_handler.scene_kwargss
            self.snap_pairs = None
            self.coord_interp_mode = None


        return_value = super().interpolateAndRender(
            galaxy_kwargs, ## only 1 dict, shared by all frames
            scene_kwargss, ## nframe dicts, 1 for each frame
            studio_kwargss, ## only 1 dict, shared by all frames
            render_kwargss, ## only 1 dict, shared by all frames
            which_studios,
            keyframes=keyframes,
            multi_threads=multi_threads if not shared_memory else 'shared',
            timestamp=timestamp,
            check_exists=check_exists,
            add_composition=add_composition,
            time_slice=time_slice)

        ## point many threads at a single shared memory buffer to render multiple orientations 
        ##  of the same snapshot simultaneously, super powerful!
        if shared_memory and return_value is not None: 
            return_value = self.scene_handler.interpolateAndRenderMultiprocessing(multi_threads,*return_value)
    
        return return_value