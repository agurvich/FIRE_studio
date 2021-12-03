import numpy as np
import os
import copy

from abg_python.plot_utils import plt,ffmpeg_frames
from abg_python.galaxy.gal_utils import Galaxy

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
        sim_time_begin=None,
        sim_time_end=None,
        fps=24,
        snapshot_times=None,
        **scene_kwargs):
        
        self.nframes = int(total_duration_sec*fps)

        ## need to interpolate in time
        if sim_time_begin and sim_time_end is not None:
            self.time_handler = TimeInterpolationHandler(
                np.linspace(sim_time_begin,sim_time_end,self.nframes),
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
        savefig='frame',
        which_studio=None,
        multi_threads=1,
        keyframes=False):

        ## handle simple case of moving camera at fixed time
        if self.time_handler is None: 
            if 'snapnum' not in galaxy_kwargs: raise KeyError("galaxy_kwargs must contain snapnum.")

            if multi_threads > 1: 
                return_value = self.scene_handler.interpolateAndRenderMultiprocessing(
                    galaxy_kwargs,
                    studio_kwargs,
                    render_kwargs,
                    savefig,
                    which_studio,
                    multi_threads,
                    keyframes)
            else:
                return_value = self.scene_handler.interpolateAndRender(
                    galaxy_kwargs,
                    studio_kwargs,
                    render_kwargs,
                    savefig,
                    which_studio,
                    keyframes)

        ## handle complex case of moving camera and incrementing time
        else:
            if keyframes: self.time_handler.keyframes = self.scene_handler.keyframes
            elif hasattr(self.time_handler,'keyframes'): del self.time_handler.keyframes

            ndiff =  self.nframes - len(self.scene_handler.frame_kwargss)
            scene_kwargs = self.scene_handler.frame_kwargss + [copy.copy(self.scene_handler.frame_kwargss[-1]) for i in range(ndiff)]

            ## merge dictionaries with priority such that
            ## studio_kwargs < this_time_kwargs < this_scene_kwargs
            frame_kwargss = [{**this_time_kwargs,**this_scene_kwargs} for 
                this_time_kwargs,this_scene_kwargs in 
                zip(self.time_handler.frame_kwargss,scene_kwargs)]

            return_value = self.time_handler.interpolateAndRender(
                galaxy_kwargs, ## only 1 dict, shared by all frames
                frame_kwargss=frame_kwargss, ## nframe dicts, 1 for each frame
                studio_kwargs=studio_kwargs, ## only 1 dict, shared by all frames
                render_kwargs=render_kwargs, ## only 1 dict, shared by all frames
                savefig=savefig,
                which_studio=which_studio,
                multi_threads=multi_threads)

        if savefig is not None:
            if 'keys_to_extract' in galaxy_kwargs: galaxy_kwargs.pop('keys_to_extract')

            galaxy = Galaxy(**galaxy_kwargs)

            format_str = '%s'%savefig + '_%0'+'%dd.png'%(np.ceil(np.log10(self.nframes)))
            ## ffmpeg the frames
            ffmpeg_frames(
                os.path.join(galaxy.datadir,'firestudio'),
                [format_str],
                savename=galaxy_kwargs['name'],
                framerate=self.scene_handler.fps,
                extension='.mp4')

        return return_value

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
    if 'this_time' in this_snapdict: my_studio.this_setup_id += "_time%.5f"%this_snapdict['this_time'] 

    ## create a new figure for this guy
    fig,ax = plt.subplots(nrows=1,ncols=1)
    my_studio.render(ax,**render_kwargs)
    if studio_kwargs['savefig'] is not None: plt.close(fig)
    else: return fig