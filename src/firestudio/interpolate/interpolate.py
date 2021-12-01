import numpy as np
import os

from abg_python.plot_utils import plt

from ..studios.gas_studio import GasStudio
from ..studios.star_studio import StarStudio 

from .time_interpolate import TimeInterpolationHandler
from .scene_interpolate import SceneInterpolationHandler

class InterpolationHandler(object):
    def __repr__(self):
        return "InterpolationHandler(%s - %s)"%(
            self.time_handler.__repr__(verbose=False),
            self.scene_handler.__repr__())

    def __init__(
        self,
        total_duration_sec,
        sim_time,
        sim_time_end=None,
        fps=15,
        snapshot_times=None,
        **scene_kwargs):
        
        self.nframes = int(total_duration_sec*fps)

        ## need to interpolate in time
        if sim_time_end is not None:
            self.time_handler = TimeInterpolationHandler(
                np.linspace(sim_time,sim_time_end,self.nframes),
                snapshot_times)
        else: self.time_handler = None

        ## need to interpolate camera orientation or other scene properties
        ##  the scene handler will have to be called interactively, I think. 
        ##  it gets so complicated trying to add stuff all at the beginning
        self.scene_handler = SceneInterpolationHandler(total_duration_sec,fps,**scene_kwargs)

    def interpolateAndRender(
        self,
        galaxy_kwargs, ## only 1 dict, shared by all frames
        studio_kwargs=None, ## only 1 dicts, shared by all frames
        render_kwargs=None, ## only 1 dict, shared by all frames
        savefig=True,
        which_studio=None,
        multi_threads=1,
        keyframes=False):

        if keyframes: self.time_handler.keyframes = self.scene_handler.keyframes
        elif hasattr(self.time_handler,'keyframes'): del self.time_handler.keyframes

        ndiff =  self.nframes - len(self.scene_handler.frame_kwargss)
        scene_kwargs = self.scene_handler.frame_kwargss + [self.scene_handler.frame_kwargss[-1]]*ndiff

        ## merge dictionaries with priority such that
        ## studio_kwargs < this_time_kwargs < this_scene_kwargs
        frame_kwargss = [{**this_time_kwargs,**this_scene_kwargs} for 
            this_time_kwargs,this_scene_kwargs in 
            zip(self.time_handler.frame_kwargss,scene_kwargs)]

        ## handle simple case of moving camera at fixed time
        if self.time_handler is None: 
            raise NotImplementedError("Need to handle case of just moving camera at fixed time")
            snapnum = np.argmin(snapshot_times-sim_time)**2
            ## now we handle a single snapnum 
            ##  ....
            return self.scene_handler.interpolateAndRender()
            if multi_threads > 1: 
                return self.scene_handler.interpolateAndRenderMultiprocessing(multi_threads=multi_threads)

        ## handle complex case of moving camera and incrementing time
        else:
            #for this_kwargs,time_kwargs in zip(frame_kwargss,self.time_handler.frame_kwargss):
                #this_kwargs.update(time_kwargs)

            return self.time_handler.interpolateAndRender(
                galaxy_kwargs, ## only 1 dict, shared by all frames
                frame_kwargss=frame_kwargss, ## nframe dicts, 1 for each frame
                studio_kwargs=studio_kwargs, ## only 1 dict, shared by all frames
                render_kwargs=render_kwargs, ## only 1 dict, shared by all frames
                savefig=savefig,
                which_studio=which_studio,
                multi_threads=multi_threads)

def worker_function(
    which_studio,
    this_snapdict,
    this_star_snapdict=None,
    studio_kwargs=None,
    add_render_kwargs=None):

    if studio_kwargs is None: studio_kwargs = {}
    if add_render_kwargs is None: add_render_kwargs = {}

    ## decide what we want to pass to the GasStudio
    if which_studio is GasStudio: render_kwargs = {
        'weight_name':'Masses',
        'quantity_name':'Temperature',
        'min_quantity':2,
        'max_quantity':7,
        'quantity_adjustment_function':np.log10,
        #'save_meta':False,
        #'use_metadata':False,
        #'min_weight':-0.5,
        #'max_weight':3,
        #'weight_adjustment_function':lambda x: np.log10(x/(30**2/1200**2)) + 10 - 6, ## msun/pc^2,
        }
    elif which_studio is StarStudio: render_kwargs = {}
    else: raise TypeError("%s is not GasStudio or StarStudio"%repr(which_studio))

    render_kwargs.update(add_render_kwargs)

    my_studio = which_studio(
        os.path.join(this_snapdict['datadir'],'firestudio'),
        this_snapdict['snapnum'], ## attribute this data to the next_snapnum's projection file
        this_snapdict['name'],
        gas_snapdict=this_snapdict,
        star_snapdict=this_star_snapdict,
        master_loud=False,
        **studio_kwargs)
    
    ## differentiate this time to << Myr precision
    my_studio.this_setup_id += "_time%.5f"%this_snapdict['this_time'] 

    ## create a new figure for this guy
    fig,ax = plt.subplots(nrows=1,ncols=1)
    my_studio.render(ax,**render_kwargs)
    if studio_kwargs['savefig'] is not None: plt.close(fig)
    else: return fig