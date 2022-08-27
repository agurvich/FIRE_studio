import gc
import multiprocessing
import itertools

import numpy as np

from abg_python.parallel.multiproc_utils import copySnapshotNamesToMPSharedMemory
from abg_python.galaxy.gal_utils import Galaxy

from .time_interpolate import TimeInterpolationHandler
from .time_helper import get_single_snap,get_load_flags,multi_worker_function,render_ffmpeg_frames

studio_kwargs = {
    'quaternion':(1,0,0,0),
    'camera_pos':(0,0,15),
    'camera_focus':(0,0,0),
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

class SceneInterpolationHandler(TimeInterpolationHandler):

    snap_pairs = None

    def __repr__(self):
        return (
            "SceneInterpolationHandler(%d/%d frames (%d keyframes) - %s)"%(
                len(self.scene_kwargss),
                self.nframes,
                len(self.keyframes),
                repr(list(self.scene_kwargss[0].keys()))))
    
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
        self.scene_kwargss = [kwargs]
        self.keyframes = [0]

    def parse_kwargs(self,**kwargs):
        if 'camera' in kwargs:
            camera = kwargs.pop('camera')
            for key in ['quaternion','camera_pos','camera_focus']:
                kwargs[key] = getattr(camera,key)

        for kwarg in list(kwargs.keys()): 
            if (kwarg in default_kwargs): pass
            else: raise KeyError(
                'Invalid key: %s - try one of:\n%s'%(
                    kwarg,
                    repr(list(default_kwargs.keys()))))

        return kwargs

    def add_keyframe(
        self,
        time_since_last_keyframe_sec,
        time_clip=True,
        loud=True,
        **kwargs):

        new_kwargs = self.parse_kwargs(**kwargs)
        prev_kwargs = self.scene_kwargss[-1]

        nsteps = int(np.ceil(time_since_last_keyframe_sec*self.fps))

        ## handle case when we are asked to interpolate past
        ##  the total duration of the movie
        if len(self.scene_kwargss) + nsteps > self.nframes:
            message = ("time since last keyframe too large,"+
            " this segment (%d + %d frames) would exceed total duration: %d frames (%.1f sec)"%(
                len(self.scene_kwargss),
                nsteps,
                self.nframes,
                self.total_duration_sec))

            if not time_clip: raise ValueError(message+". Use time_clip=True to avoid this error message.")
            else: 
                nsteps = self.nframes - len(self.scene_kwargss)
                if loud: print(message+'... clipping to %d frames instead'%nsteps)
        
        ## handle case where we are not changing a previously specified kwarg
        for prev_kwarg in prev_kwargs:
            if prev_kwarg not in new_kwargs:
                ##  just copy the old value over
                new_kwargs[prev_kwarg] = prev_kwargs[prev_kwarg]

        ## make sure each dictionary has one-to-one
        ##  corresponding keys: 
        for new_kwarg in new_kwargs:
            ## handle invalid kwarg
            if new_kwarg not in default_kwargs: raise KeyError(
                'kwarg %s must be one of:\n%s'%(new_kwarg,repr(default_kwargs.keys())))

            ## handle case where a new kwarg is not in the previous
            if new_kwarg not in prev_kwargs:
                ## *explicitly* set the previous values for each frame to be the default
                for sub_prev_kwargs in self.scene_kwargss:
                    sub_prev_kwargs[new_kwarg] = default_kwargs[new_kwarg]
        
        if nsteps == 1: self.scene_kwargss.append(new_kwargs)
        ## start at i = 1 to avoid repeating frames
        for i in range(1,nsteps+1):
            this_kwargs = {}
            for kwarg in new_kwargs:
                pval = prev_kwargs[kwarg]
                nval = new_kwargs[kwarg]
                ## convert args that are lists/tuples to arrays
                if kwarg in ['quaternion','camera_pos','camera_focus','font_color']:
                    pval = np.array(pval)
                    nval = np.array(nval)
                ## TODO should have some kind of interpolation function
                ##  so we don't have to do just linear
                ##  then again we can always string together keyframes
                ##  to get complex interpolations
                this_kwargs[kwarg] = pval + i*(nval-pval)/(nsteps)
            self.scene_kwargss.append(this_kwargs)

        ## note the index of this keyframe
        self.keyframes.append(len(self.scene_kwargss)-1)

        if loud: print(self)
    
    def interpolateAndRenderMultiprocessing(
        self,
        multi_threads,
        galaxy_kwargs,
        scene_kwargss=None,
        studio_kwargss=None,
        render_kwargss=None,
        which_studios=None,
        fixed_star_hsml=0.028
        ):


        if 'keys_to_extract' in galaxy_kwargs.keys(): keys_to_extract = galaxy_kwargs.pop('keys_to_extract')
        else: keys_to_extract = []

        load_gas,load_star = get_load_flags(which_studios,render_kwargss)

        gas_snapdict,star_snapdict = get_single_snap(
            load_gas,
            load_star,
            keys_to_extract=keys_to_extract,
            **galaxy_kwargs)

        global_snapdict_name = 'gas_snapshot_%03d'%galaxy_kwargs['snapnum']
        global_star_snapdict_name = 'star_snapshot_%03d'%galaxy_kwargs['snapnum']

        ## if we were bold enough to extract everything, copy nothing to the child processes.
        ##  that'll teach us!
        #if keys_to_extract is None: 
            ### todo, why not just use all the keys if they're going to go to a shared memory buffer?
            #raise KeyError("Use keys_to_extract to specify field keys you need for rendering,"+
            #" they're going to be put into a shared memory buffer so we will *not* pass all keys by default.")

        if multi_threads is None: multi_threads = multiprocessing.cpu_count()-1

        ## collect positional arguments for worker_function
        argss = zip(
            itertools.repeat(which_studios),
            itertools.repeat(global_snapdict_name),
            itertools.repeat(global_star_snapdict_name),
            scene_kwargss,
            itertools.repeat(studio_kwargss),
            itertools.repeat(render_kwargss))

        ## initialize dictionary that will point to shared memory buffers
        gas_wrapper_dict = {}
        star_wrapper_dict = {}

        try:
            if load_gas:
                ## use as few references so i have to clean up fewer below lol
                gas_wrapper_dict,gas_shm_buffers = copySnapshotNamesToMPSharedMemory(
                    ['Coordinates',
                    'Masses',
                    'SmoothingLength']+keys_to_extract,
                    gas_snapdict,
                    finally_flag=True,
                    loud=True)
            else: gas_shm_buffers = [None]

            if load_star:

                if 'SmoothingLength' not in star_snapdict: 
                    star_snapdict['SmoothingLength'] = np.repeat(fixed_star_hsml,star_snapdict['Coordinates'].shape[0])

                ## NOTE the lack of smoothing lengths might mess this up if a bunch of processes all
                ##  try and compute smoothing lengths and write to the same file :\
                star_wrapper_dict,star_shm_buffers = copySnapshotNamesToMPSharedMemory(
                    ['Coordinates',
                    'Masses',
                    'SmoothingLength',
                    'AgeGyr']+keys_to_extract,
                    star_snapdict,
                    finally_flag=True,
                    loud=True)
                
            else: star_shm_buffers = [None]
                

            for key in ['name','datadir','snapnum']: gas_wrapper_dict[key] = gas_snapdict[key]
            for key in ['name','datadir','snapnum']: star_wrapper_dict[key] = star_snapdict[key]

            del gas_snapdict,star_snapdict
            globals()[global_snapdict_name] = gas_wrapper_dict
            globals()[global_star_snapdict_name] = star_wrapper_dict

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
                these_figs = my_pool.starmap(multi_worker_function_wrapper,argss)

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
            del gas_wrapper_dict
            del star_wrapper_dict
            for shm_buffer in gas_shm_buffers:
                ## handle case where multiprocessing isn't used
                if shm_buffer is not None:
                    shm_buffer.close()
                    shm_buffer.unlink()

            del gas_shm_buffers
            for shm_buffer in star_shm_buffers:
                ## handle case where multiprocessing isn't used
                if shm_buffer is not None:
                    shm_buffer.close()
                    try: shm_buffer.unlink()
                    except FileNotFoundError: pass
            del star_shm_buffers

        ## use ffmpeg to produce an mp4 of the frames
        render_ffmpeg_frames(studio_kwargss,galaxy_kwargs,self.nframes,self.fps)

        return these_figs

def multi_worker_function_wrapper(
    which_studios,
    global_snapdict_name, ## important to access shared memory :\
    global_star_snapdict_name,
    scene_kwargss,
    studio_kwargss,
    render_kwargss):

    ## read the unique global name for the relevant snapshot dictionary
    ##  TODO: could I handle time interpolation right here by checking if 
    ##  if I was passed multiple snapdict names... then I could compute
    ##  current_time_gyr and make a new combination snapshotdictionary 
    ##  that was interpolated.
    multi_worker_function(
        which_studios,
        globals()[global_snapdict_name],
        globals()[global_star_snapdict_name],
        scene_kwargss,
        studio_kwargss,
        render_kwargss)