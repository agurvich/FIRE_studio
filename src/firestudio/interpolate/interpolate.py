import numpy as np
import os
import copy
import warnings

from abg_python.plot_utils import plt,ffmpeg_frames
from abg_python.galaxy.gal_utils import Galaxy

from ..studios.gas_studio import GasStudio
from ..studios.star_studio import StarStudio 
from ..studios.FIRE_studio import FIREStudio
from ..studios.composition import Composition

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
        time_slice=None,
        coord_interp_mode='spherical',
        **scene_kwargs):
        
        if time_slice is None: time_slice = slice(0,None)

        ## need to interpolate camera orientation or other scene properties
        ##  the scene handler will have to be called interactively, I think. 
        ##  it gets so complicated trying to add stuff all at the beginning
        self.scene_handler = SceneInterpolationHandler(total_duration_sec,fps,**scene_kwargs)

        ## need to interpolate in time
        if sim_time_begin and sim_time_end is not None:
            self.time_handler = TimeInterpolationHandler(
                np.linspace(sim_time_begin,sim_time_end,total_duration_sec*fps)[time_slice],
                snapshot_times,
                coord_interp_mode=coord_interp_mode)
            self.nframes = self.time_handler.nframes
        else: 
            self.time_handler = None
            self.nframes = fps*total_duration_sec
        
        self.scene_handler.nframes = self.nframes

    def interpolateAndRender(
        self,
        galaxy_kwargs, ## only 1 dict, shared by all frames
        studio_kwargss=None, ## only 1 dicts, shared by all frames
        render_kwargss=None, ## only 1 dict, shared by all frames
        savefigs=None,
        which_studios=None,
        multi_threads=1,
        keyframes=False,
        check_exists=True,
        timestamp=0, ## offset in Gyr for timestamp, None = no timestamp
        add_composition=False,
        ):

        ## will default to gas studio
        if which_studios is None: which_studios = [GasStudio]
        if savefigs is None: savefigs = [which_studio.__name__ for which_studio in which_studios]

        if add_composition:
            c_studio_kwargs = {
                'studios_tuple':tuple(which_studios),
                'subplots_kwargs':{
                    'wspace':0,'hspace':0,
                    'left':0,'right':1,
                    'bottom':0,'top':1},
                'savefig':'_'.join([which_studio.__name__ for which_studio in which_studios]),
                'studio_kwargss':copy.deepcopy(studio_kwargss),
                'size_inches':(12,6),
                }

            ## tell the interpolator to initialize a composition
            which_studios = which_studios+[Composition]
            ## add the kwargs that should be passed to the Composition at initialization
            studio_kwargss = studio_kwargss+[c_studio_kwargs]
            ## join the render kwargs, they'll be ignored by the studios that don't need them
            c_render_kwargs = {}
            for this_render_kwargs in render_kwargss: c_render_kwargs.update(this_render_kwargs)
            render_kwargss = render_kwargss + [c_render_kwargs]
            ## and let's go ahead and add the savefig string to the list of outputs
            savefigs = savefigs + [c_studio_kwargs['savefig']]

        ## handle simple case of moving camera at fixed time
        if self.time_handler is None: 
            raise NotImplementedError("Haven't updated scene interpolation for multiple studios")
            if 'snapnum' not in galaxy_kwargs: raise KeyError("galaxy_kwargs must contain snapnum.")

            if multi_threads > 1: 
                return_value = self.scene_handler.interpolateAndRenderMultiprocessing(
                    galaxy_kwargs,
                    studio_kwargs,
                    render_kwargs,
                    savefigs,
                    which_studio,
                    multi_threads,
                    keyframes,
                    check_exists)
            else:
                return_value = self.scene_handler.interpolateAndRender(
                    galaxy_kwargs,
                    studio_kwargs,
                    render_kwargs,
                    savefigs,
                    which_studio,
                    keyframes,
                    check_exists)

        ## handle complex case of moving camera and incrementing time
        else:
            if keyframes: self.time_handler.keyframes = self.scene_handler.keyframes
            elif hasattr(self.time_handler,'keyframes'): del self.time_handler.keyframes

            ndiff =  self.nframes - len(self.scene_handler.frame_kwargss)
            scene_kwargs = self.scene_handler.frame_kwargss + [copy.copy(self.scene_handler.frame_kwargss[-1]) for i in range(ndiff)]

            ## merge dictionaries with priority such that
            ## studio_kwargs < this_time_kwargs < this_scene_kwargs
            scene_kwargss = [{**this_time_kwargs,**this_scene_kwargs} for 
                this_time_kwargs,this_scene_kwargs in 
                zip(self.time_handler.scene_kwargss,scene_kwargs)]

            return_value = self.time_handler.interpolateAndRender(
                galaxy_kwargs, ## only 1 dict, shared by all frames
                scene_kwargss=scene_kwargss, ## nframe dicts, 1 for each frame
                studio_kwargss=studio_kwargss, ## only 1 dict, shared by all frames
                render_kwargss=render_kwargss, ## only 1 dict, shared by all frames
                savefigs=savefigs,
                which_studios=which_studios,
                multi_threads=multi_threads,
                check_exists=check_exists,
                timestamp=timestamp)

        if savefigs is not None:
            if 'keys_to_extract' in galaxy_kwargs: galaxy_kwargs.pop('keys_to_extract')

            if 'snapnum' not in galaxy_kwargs: galaxy_kwargs['snapnum'] = None
            galaxy = Galaxy(**galaxy_kwargs)

            ## in case we rendered multiple types of frames
            ##  we'll loop through a list of savefigs (even if there's only one)
            for this_savefig in savefigs:
                if this_savefig is not None:
                    format_str = '%s'%this_savefig + '_frame_%0'+'%dd.png'%(np.ceil(np.log10(self.nframes)))
                    ## ffmpeg the frames
                    ffmpeg_frames(
                        os.path.join(galaxy.datadir,'firestudio'),
                        [format_str],
                        savename=galaxy_kwargs['name'],
                        framerate=self.scene_handler.fps,
                        extension='.mp4')

        return return_value


def multi_worker_function(
    which_studios,
    this_snapdict,
    this_star_snapdict=None,
    scene_kwargs=None,
    studio_kwargss=None,
    add_render_kwargss=None):

    if len(which_studios) != len(studio_kwargss) != len(add_render_kwargss):
        raise ValueError(
            "Must have matching length lists of: "+
            "studios, studio kwargs, and render kwargs")

    for which_studio,studio_kwargs,add_render_kwargs in zip(
        which_studios,studio_kwargss,add_render_kwargss):
        try:
            worker_function(
                which_studio,
                this_snapdict,
                this_star_snapdict,
                {**scene_kwargs,**studio_kwargs},
                add_render_kwargs)
        except AssertionError: raise
        except:
            print(this_snapdict['snapnum'],'failed')
            raise

def worker_function(
    which_studio,
    this_snapdict,
    this_star_snapdict=None,
    studio_kwargs=None,
    add_render_kwargs=None):

    if studio_kwargs is None: studio_kwargs = {'savefig':None}
    if add_render_kwargs is None: add_render_kwargs = {}

    if studio_kwargs['savefig'] is not None and 'savefig_suffix' in studio_kwargs:
        this_savefig = studio_kwargs['savefig'] + studio_kwargs.pop('savefig_suffix')
    else: this_savefig = None

    ## decide what we want to pass to the GasStudio
    if which_studio is GasStudio: render_kwargs = {
        'weight_name':'Masses',
        'quantity_name':'Temperature',
        'min_quantity':2,
        'max_quantity':7,
        'quantity_adjustment_function':np.log10,
        #'min_weight':-0.5,
        #'max_weight':3,
        }
    elif which_studio is StarStudio: render_kwargs = {}
    elif which_studio is FIREStudio: render_kwargs = {}
    elif which_studio is Composition: 
        studios_tuple = studio_kwargs.pop('studios_tuple')
        which_studio_actual = which_studio
        which_studio = lambda *args,**kwargs: which_studio_actual(
            studios_tuple,*args,**kwargs)
        render_kwargs = {}
    else: raise TypeError("%s is not GasStudio or StarStudio"%repr(which_studio))

    render_kwargs.update(add_render_kwargs)

    #print(which_studio)
    #print(studio_kwargs)
    #print(render_kwargs)
    my_studio = which_studio(
        os.path.join(this_snapdict['datadir'],'firestudio'),
        this_snapdict['snapnum'], ## attribute this data to the next_snapnum's projection file
        this_snapdict['name'],
        gas_snapdict=this_snapdict,
        star_snapdict=this_star_snapdict,
        master_loud=False,
        setup_id_append="_time%.5f"%this_snapdict['this_time'],
        **{**studio_kwargs,'savefig':this_savefig})

    if which_studio is GasStudio and render_kwargs['weight_name'] == 'Masses':
        render_kwargs['weight_adjustment_function'] = lambda x: np.log10(x/my_studio.Acell) + 10 - 6 ## msun/pc^2,
    

    ## ignore warnings to reduce console spam
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ## overwrite loud rendering to reduce console spam
        ax,im = my_studio.render(None,**{**render_kwargs,'loud':False})
    axs = np.array([ax]).reshape(-1)
    fig = axs[0].get_figure()

    if my_studio.savefig is not None: 
        plt.close(fig)
        print(this_snapdict['snapnum'],my_studio.savefig)
    else: return fig