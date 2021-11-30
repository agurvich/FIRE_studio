import gc
import multiprocessing
import itertools

import numpy as np

from abg_python.parallel.multiproc_utils import copySnapshotNamesToMPSharedMemory
from abg_python.galaxy.gal_utils import Galaxy

from ..utils.camera_utils import Camera
from ..studios.gas_studio import GasStudio
from ..studios.star_studio import StarStudio

studio_kwargs = {
    'quaternion':(1,0,0,0),
    'frame_half_thickness':15, ## half-thickness of image in z direction
    'aspect_ratio':1, ## shape of image, y/x TODO figure out if this is necessary to pass?
    'pixels':1200, ## pixels in x direction, resolution of image
    'scale_line_length':5, ## length of the scale line in kpc
    'fontsize':12,  ## font size of scale bar and figure label
    'font_color':(1,1,1,1), ## font color of scale bar and figure label

    #'figure_label':'', ## string to be put in upper right corner

    #'scale_bar':True,  ## flag to plot length scale bar in lower left corner
    #'figure_label_side':'right', ## corner to put label in
    #'noaxis':True, ## turns off axis ticks

    #'savefig':None, ## save the image as a png if passed a string
    #'snapdir':None,
    #'snapnum':None,
    #'sim_name':None
    }

star_kwargs = {
    'maxden' : None, ## 
    'dynrange' : None}, ## controls the saturation of the image in a non-obvious way

## set_ImageParams only controls basic aspects of the colorbar which
##  don't make sense to try and interpolate
gas_kwargs = {}

default_kwargs ={}
for kwargs in [studio_kwargs,star_kwargs,gas_kwargs]: default_kwargs.update(kwargs)

class SceneInterpolationHandler(object):
    def __repr__(self):
        return (
            "SceneInterpolationHandler(%d/%d frames (%d keyframes) - %s)"%(
                len(self.frame_kwargss),
                self.nframes,
                self.nkeyframes,
                repr(list(self.frame_kwargss[0].keys()))))
    
    def __getitem__(self,key):
        return self.kwargs[key]

    def __init__(
        self,
        total_duration_sec,
        fps=15,
        **kwargs):

        self.fps = fps
        self.total_duration_sec = total_duration_sec
        self.nframes = int(self.total_duration_sec*self.fps)

        kwargs = self.parse_kwargs(**kwargs)

        ## initialize the first keyframe. we'll interpolate from it on the first call to 
        ##  self.add_keyframe and we'll copy it if, for whatever reason, we make it to 
        ##  interpolater.render
        self.frame_kwargss = [kwargs]
        self.nkeyframes = 1

    def parse_kwargs(self,**kwargs):
        camera_kwargs = {'camera_pos':[0,0,15],}
        for kwarg in list(kwargs.keys()): 
            ## overwrite the camera_kwargs with anything passed through
            if 'camera_' in kwarg: camera_kwargs[kwarg] = kwargs.pop(kwarg)

            elif (kwarg in default_kwargs): pass
            else: raise KeyError(
                'Invalid key: %s - try one of:\n%s'%(
                    kwarg,
                    repr(list(default_kwargs.keys()))))

        if 'quaternion' not in kwargs:
            ## will want to initialize the camera in the studio instance w/ 
            ##  quaternion that is interpolated from the quaternions taken
            ##  from the cameras initialized by the interpolater here
            if 'camera' not in kwargs or kwargs['camera'] is None:
                kwargs['camera'] = Camera(**camera_kwargs)
        
            kwargs['quaternion'] = kwargs.pop('camera').quaternion


        return kwargs

    def add_keyframe(
        self,
        time_since_last_keyframe_sec,
        time_clip=True,
        loud=True,
        **kwargs):

        new_kwargs = self.parse_kwargs(**kwargs)
        prev_kwargs = self.frame_kwargss[-1]

        nsteps = int(np.ceil(time_since_last_keyframe_sec*self.fps))

        ## handle case when we are asked to interpolate past
        ##  the total duration of the movie
        if len(self.frame_kwargss) + nsteps > self.nframes:
            message = ("time since last keyframe too large,"+
            " this segment (%d + %d frames) would exceed total duration: %d frames (%d sec)"%(
                len(self.frame_kwargss),
                nsteps,
                self.nframes,
                self.total_duration_sec))

            if not time_clip: raise ValueError(message+". Use time_clip=True to avoid this error message.")
            else: 
                nsteps = self.nframes - len(self.frame_kwargss)
                if loud: print(message+'... clipping to %d frames instead'%nsteps)
        
        ## make sure each dictionary has one-to-one
        ##  corresponding keys: 
        for new_kwarg in new_kwargs:
            ## handle invalid kwarg
            if new_kwarg not in default_kwargs: raise KeyError(
                'kwarg %s must be one of:\n%s'%(new_kwarg,repr(default_kwargs.keys())))

            ## handle case where a new kwarg is not in the previous
            if new_kwarg not in prev_kwargs:
                ## *explicitly* set the previous values for each frame to be the default
                for sub_prev_kwargs in self.frame_kwargss:
                    sub_prev_kwargs[new_kwarg] = default_kwargs[new_kwarg]
        
        ## handle case where we are not changing a previously specified kwarg
        for prev_kwarg in prev_kwargs:
            if prev_kwarg not in new_kwargs:
                ##  just copy the old value over
                new_kwargs[prev_kwarg] = prev_kwargs[prev_kwarg]

        if nsteps == 1: self.frame_kwargss.append(new_kwargs)
        ## start at i = 1 to avoid repeating frames
        for i in range(1,nsteps+1):
            this_kwargs = {}
            for kwarg in new_kwargs:
                pval = prev_kwargs[kwarg]
                nval = new_kwargs[kwarg]
                ## TODO should have some kind of interpolation function
                ##  so we don't have to do just linear
                ##  then again we can always string together keyframes
                ##  to get complex interpolations
                this_kwargs[kwarg] = pval + i*(nval-pval)/(nsteps-1)
            self.frame_kwargss.append(this_kwargs)

        self.nkeyframes+=1
        if loud: print(self)
    
    def interpolateAndRender(
        self,
        galaxy_kwargs,
        studio_kwargs=None,
        render_kwargs=None,
        savefig=True,
        which_studio=None,
        multi_threads=None):

        ## put here to avoid circular import
        from .interpolate import worker_function

        if galaxy_kwargs is None: raise ValueError(
            'galaxy_kwargs must be a dictionary with,'+
            ' at minimum, name (i.e. m12i_res7100) and snapnum (i.e. 600).')

        galaxy = Galaxy(**galaxy_kwargs)
        galaxy.extractMainHalo()

        if multi_threads is not None: raise ValueError("Use interpolateAndRenderMultiprocessing instead.")

        ## handle default arguments
        if studio_kwargs is None: studio_kwargs = {}
        if render_kwargs is None: render_kwargs = {}

        frame_kwargss = [{**this_frame_kwargs,**studio_kwargs} for this_frame_kwargs in self.frame_kwargss]

        ## determine which studio we should initialize inside the worker_function
        if which_studio is None: which_studio = GasStudio
        elif which_studio is not GasStudio and which_studio is not StarStudio: 
            raise TypeError("%s is not GasStudio or StarStudio"%repr(which_studio))

        ## initialize array of savefig values
        savefigs = ['frame_%04d.png'%i if savefig else False 
            for i in range(len(self.nframes))]

        ## collect positional arguments for worker_function
        argss = zip(
            itertools.repeat(which_studio),
            itertools.repeat(galaxy.sub_snap),
            itertools.repeat(galaxy.sub_star_snap),
            frame_kwargss,
            itertools.repeat(render_kwargs),
            savefigs)

        these_figs = [worker_function(*args) for args in argss]

        return these_figs

    def interpolateAndRenderMultiprocessing(
        self,
        galaxy_kwargs,
        studio_kwargs=None,
        render_kwargs=None,
        savefig=True,
        which_studio=None,
        multi_threads=1):

        if galaxy_kwargs is None: raise ValueError(
            'galaxy_kwargs must be a dictionary with,'+
            ' at minimum, name (i.e. m12i_res7100) and snapnum (i.e. 600).')

        galaxy = Galaxy(**galaxy_kwargs)
        galaxy.extractMainHalo()

        global_snapdict_name = 'gas_snapshot_%03d'%galaxy.snapnum
        global_star_snapdict_name = 'star_snapshot_%03d'%galaxy.snapnum

        if multi_threads is None: multi_threads = multiprocessing.cpu_count()-1

        ## handle default arguments
        if studio_kwargs is None: studio_kwargs = {}
        if render_kwargs is None: render_kwargs = {}

        frame_kwargss = [{**this_frame_kwargs,**studio_kwargs} for this_frame_kwargs in self.frame_kwargss]

        ## determine which studio we should initialize inside the worker_function
        if which_studio is None: which_studio = GasStudio
        elif which_studio is not GasStudio and which_studio is not StarStudio: 
            raise TypeError("%s is not GasStudio or StarStudio"%repr(which_studio))

        ## initialize array of savefig values
        savefigs = ['frame_%04d.png'%i if savefig else False 
            for i in range(len(self.nframes))]

        ## collect positional arguments for worker_function
        argss = zip(
            itertools.repeat(which_studio),
            itertools.repeat(global_snapdict_name),
            itertools.repeat(global_star_snapdict_name),
            frame_kwargss,
            itertools.repeat(render_kwargs),
            savefigs)

        ## initialize dictionary that will point to shared memory buffers
        wrapper_dict = {}
        try:
            ## use as few references so i have to clean up fewer below lol
            wrapper_dict,shm_buffers = copySnapshotNamesToMPSharedMemory(
                ['Coordinates',
                'Masses',
                'Temperature',
                'SmoothingLength'],
                galaxy.sub_snap,
                finally_flag=True,
                loud=True)

            del galaxy
            globals()[global_snapdict_name] = wrapper_dict
            globals()[global_star_snapdict_name] = None

            ## don't remove these lines, they perform some form of dark arts
            ##  that helps the garbage collector its due
            ## attempt to wrangle shared memory buffer and avoid memory leak 
            locals().keys()
            globals().keys()
            gc.collect()

            ## attempt to wrangle shared memory buffer and avoid memory leak 
            for obj in gc.get_objects():
                if isinstance(obj,Galaxy):
                    print(obj,'will be copied to child processes and is probably large.')
            
            with multiprocessing.Pool(multi_threads) as my_pool:
                these_figs = my_pool.starmap(worker_function_wrapper,argss)

            ## attempt to wrangle shared memory buffer and avoid memory leak 
            del my_pool
            locals().keys()
            globals().keys()
            gc.collect()
        except: raise
        finally:
            ## TODO clean up anything that contains a reference to a shared
            ##  memory object. globals() must be purged before the shm_buffers
            ##  are unlinked or python will crash.
            globals().pop(global_snapdict_name)
            globals().pop(global_star_snapdict_name)
            del wrapper_dict
            for shm_buffer in shm_buffers:
                ## handle case where multiprocessing isn't used
                if shm_buffer is not None:
                    shm_buffer.close()
                    shm_buffer.unlink()

            del shm_buffers

        return these_figs

def worker_function_wrapper(
    which_studio,
    global_snapdict_name, ## important to access shared memory :\
    global_star_snapdict_name,
    studio_kwargs=None,
    add_render_kwargs=None,
    savefig=False):

    ## put here to avoid circular import
    from .interpolate import worker_function

    ## read the unique global name for the relevant snapshot dictionary
    ##  TODO: could I handle time interpolation right here by checking if 
    ##  if I was passed multiple snapdict names... then I could compute
    ##  current_time_gyr and make a new combination snapshotdictionary 
    ##  that was interpolated.
    ##  TODO: think more about if this is how I want to do this if multiprocessing 
    ##  is turned off which should be the default mode tbh.
    this_snapdict = globals()[global_snapdict_name]
    this_star_snapdict = globals()[global_star_snapdict_name]
    worker_function(which_studio,this_snapdict,this_star_snapdict,studio_kwargs,add_render_kwargs,savefig)