import itertools 

import matplotlib.pyplot as plt
import numpy as np 
import os

import time
import gc

import multiprocessing
import itertools

from abg_python.multiproc_utils import copySnapshotNamesToMPSharedMemory
from abg_python.galaxy.gal_utils import Galaxy
import abg_python.all_utils as all_utils

from firestudio.studios.gas_studio import GasStudio
from firestudio.studios.studio import Studio
from firestudio.studios.star_studio import StarStudio

def interpolationHelper(duration,previous_chain=None,framerate=15,this_class=None,**kwargs):
    if this_class is None:
        this_class = Studio
        
    if previous_chain is None:
        previous_chain = {'nsteps_tot':0,'nlinks':0,'nstepss':[]}
    
    nsteps = duration*framerate
    for kwarg in kwargs:
        if kwarg+'s' not in previous_chain:
            if previous_chain['nlinks'] > 0:
                try:
                    previous_chain[kwarg+'s'] = ([
                        this_class.set_ImageParams.default_kwargs[kwarg]] *
                        (previous_chain['nlinks']+1))
                except KeyError:
                    print("Looked in",
                        this_class,
                        "pass the class you'd like to look in"+
                        " using the this_class kwarg")
            else:
                previous_chain[kwarg+'s'] = [kwargs[kwarg][0]]
                kwargs[kwarg]=kwargs[kwarg][1]
                
        previous_chain[kwarg+'s']+=[kwargs[kwarg]]
    
    for pkwarg in previous_chain:
        if pkwarg[:-1] in kwargs or pkwarg in ['nlinks','nsteps_tot','nstepss']:
            continue
        previous_chain[pkwarg]+=[previous_chain[pkwarg][-1]]
    
    previous_chain['nstepss']+=[nsteps]
    previous_chain['nlinks']+=1
    previous_chain['nsteps_tot']+=nsteps
    return previous_chain
    
class Interpolater(object):
    
    def __init__(
        self,
        nstepss,
        keyframes_only=False,
        **kwargs):


        self.interp_kwargs = [
            'frame_half_widths',
            'frame_half_thicknesss',
            'frame_centers',
            'thetas',
            'phis',
            'psis',
            'aspect_ratios',
            'snapnums']


        ## perform each interpolation in a single step
        ##  producing only the interpolation key frames
        if keyframes_only:
            nstepss = [[1] for i in range(len(nstepss))]

        self.nstepss = nstepss

        for kwarg in kwargs.keys():
            if kwarg == 'aspect_ratios':
                raise NotImplementedError(
                    "not sure aspect_ratios should be allowed (yet)")
            elif kwarg == 'snapnums':
                raise NotImplementedError(
                    "need to sort lists by snapnum and connect to"+
                    " a time interpolater class")
            elif kwarg not in self.interp_kwargs:
                raise KeyError(kwarg+" not allowed",self.interp_kwargs)
            else:
                value = kwargs[kwarg]

                ## passes quality checks?
                if len(value) != (len(nstepss)+1):
                    raise ValueError(
                        kwarg+
                        " has wrong length this: %d req: %d"%(len(value),len(nstepss)+1))

                interpds = []
                if kwarg != 'frame_centers':
                    for li in range(len(nstepss)):
                        ri = li+1
                        interpds.append(np.linspace(
                            kwargs[kwarg][li],
                            kwargs[kwarg][ri],
                            nstepss[li]+1,
                            endpoint=True))  
                else:
                    frame_centers = kwargs[kwarg]
                    for li in range(len(nstepss)):
                        ri = li+1
                        xs = np.linspace(
                            frame_centers[li][0],
                            frame_centers[ri][0],
                            nstepss[li]+1,
                            endpoint=True)  
                        ys = np.linspace(
                            frame_centers[li][1],
                            frame_centers[ri][1],
                            nstepss[li]+1,
                            endpoint=True)  
                        zs = np.linspace(
                            frame_centers[li][2],
                            frame_centers[ri][2],
                            nstepss[li]+1,
                            endpoint=True)
                        interpds.append(np.array([xs,ys,zs]).T)

                for li in range(len(nstepss)):
                    ## li > 0 business is to avoid repeating
                    ##  edges
                    ## i.e. frame_half_widths = [30,5,30]
                    ##  should go as [30,20,10,5,10,20,30]
                    ##  not [30,20,10,5,5,10,20,30]
                    if li > 0 :
                        interpds[li]=interpds[li][1:]

            ## bind to the object
            setattr(
                self,
                kwarg,
                interpds)

        ## see above comment for li > 0 business
        #for li in range(len(nstepss)):
        nstepss[0]+=1

    def interpolateAndRender(self,frame_offset=0,galaxy=None):

        if galaxy is None:
            snapdir = "/scratch/projects/xsede/GalaxiesOnFIRE/metal_diffusion/m12i_res7100/output"
            snapnum = 600 
            galaxy = Galaxy(
                'm12i_res7100',
                snapdir,
                600,
                datadir='/scratch/04210/tg835099/data/metal_diffusion')
            galaxy.extractMainHalo()
        else:
            snapnum = galaxy.snapnum

        ## let's put the FIREstudio projections into a sub-directory of our Galaxy class instance
        studio_datadir = os.path.join(os.path.dirname(galaxy.datadir),'firestudio')

        global render_kwargs
        render_kwargs = {
            'weight_name':'Masses',
            'quantity_name':'Temperature',
            #min_quantity=2,
            #max_quantity=7,
            'quantity_adjustment_function':np.log10,
            'quick':True,
            'min_weight':-0.5,
            'max_weight':3,
            'weight_adjustment_function':lambda x: np.log10(x/(30**2/1200**2)) + 10 - 6, ## msun/pc^2,
            'cmap':'afmhot'}

         ## pass in snapshot dictionary
        wrapper_dict = {}
        global_this_snapdict_name = 'gas_snapshot_%03d'%snapnum
        ## use as few references so i have to clean up fewer below lol
        wrapper_dict = galaxy.sub_snap
        globals()[global_this_snapdict_name] = wrapper_dict
        ## don't remove these lines, they perform some form of dark arts
        ##  that helps the garbage collector its due

        frame_num = frame_offset
        for interp_i,nsteps in enumerate(self.nstepss):
            im_param_kwargs = []
            for step in range(nsteps):
                this_im_param_kwargs = {'frame_num':frame_num}
                frame_num+=1
                ## load the kwargs for this interpolation frame
                for interp_kwarg in self.interp_kwargs:
                    if hasattr(self,interp_kwarg):
                        this_im_param_kwargs[interp_kwarg[:-1]] = getattr(self,interp_kwarg)[interp_i][step]
                im_param_kwargs.append(this_im_param_kwargs)

            args = zip(
                itertools.repeat(GasStudio), ## this class,gas studio or star studio
                itertools.repeat(studio_datadir), ## datadir
                itertools.repeat(snapnum), ## snapnum for projection file to be saved to
                itertools.repeat(global_this_snapdict_name), ## what to look up in globals() for gas
                im_param_kwargs) ##

            these_axs = [worker_function(*arg) for arg in args]


        return these_axs

    def interpolateAndRenderMultiprocessing(self,frame_offset=0,galaxy=None,nproc=None):

        if nproc is None:
            nproc = multiprocessing.cpu_count()-1

        if galaxy is None:
            snapdir = "/scratch/projects/xsede/GalaxiesOnFIRE/metal_diffusion/m12i_res7100/output"
            snapnum = 600 
            galaxy = Galaxy(
                'm12i_res7100',
                snapdir,
                600,
                datadir='/scratch/04210/tg835099/data/metal_diffusion')


            galaxy.extractMainHalo()
        else:
            snapnum = galaxy.snapnum

        ## let's put the FIREstudio projections into a sub-directory of our Galaxy class instance
        studio_datadir = os.path.join(os.path.dirname(galaxy.datadir),'firestudio')

        global render_kwargs
        render_kwargs = {
            'weight_name':'Masses',
            'quantity_name':'Temperature',
            #min_quantity=2,
            #max_quantity=7,
            'quantity_adjustment_function':np.log10,
            'quick':True,
            'min_weight':-0.5,
            'max_weight':3,
            'weight_adjustment_function':lambda x: np.log10(x/(30**2/1200**2)) + 10 - 6, ## msun/pc^2,
            'cmap':'afmhot'}

        ## pass in snapshot dictionary
        wrapper_dict = {}
        try:
            global_this_snapdict_name = 'gas_snapshot_%03d'%snapnum
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
            globals()[global_this_snapdict_name] = wrapper_dict
            ## don't remove these lines, they perform some form of dark arts
            ##  that helps the garbage collector its due
            locals().keys()
            globals().keys()
            gc.collect()

            for obj in gc.get_objects():
                if isinstance(obj,Galaxy):
                    print(obj,'will be copied to child processes and is probably large.')
            
            frame_num = frame_offset
            for interp_i,nsteps in enumerate(self.nstepss):
                im_param_kwargs = []
                for step in range(nsteps):
                    this_im_param_kwargs = {'frame_num':frame_num}
                    frame_num+=1
                    ## load the kwargs for this interpolation frame
                    for interp_kwarg in self.interp_kwargs:
                        if hasattr(self,interp_kwarg):
                            this_im_param_kwargs[interp_kwarg[:-1]] = getattr(self,interp_kwarg)[interp_i][step]
                    im_param_kwargs.append(this_im_param_kwargs)

                args = zip(
                    itertools.repeat(GasStudio), ## this class,gas studio or star studio
                    itertools.repeat(studio_datadir), ## datadir
                    itertools.repeat(snapnum), ## snapnum for projection file to be saved to
                    itertools.repeat(global_this_snapdict_name), ## what to look up in globals() for gas
                    im_param_kwargs) ##

                this_nprocs = min(nproc,len(im_param_kwargs))
                print("Starting a pool of %d processes..."%this_nprocs)
                with multiprocessing.Pool(this_nprocs) as my_pool:
                    my_pool.starmap(worker_function,args)
                print("pool finished.")
                del my_pool
                del args
                del im_param_kwargs
                locals().keys()
                globals().keys()
                gc.collect()
        except:
            raise
        finally:
            ## TODO clean up anything that contains a reference to a shared
            ##  memory object. globals() must be purged before the shm_buffers
            ##  are unlinked or python will crash.
            globals().pop(global_this_snapdict_name)
            del wrapper_dict
            for shm_buffer in shm_buffers:
                ## handle case where multiprocessing isn't used
                if shm_buffer is not None:
                    shm_buffer.close()
                    shm_buffer.unlink()

            del shm_buffers

        return None 

def worker_function(
    this_class,
    datadir,
    snapnum,
    global_this_snapdict_name,
    im_param_kwargs):

    ## read the unique global name for the relevant snapshot dictionary
    ##  TODO: could I handle time interpolation right here by checking if 
    ##  if I was passed multiple snapdict names... then I could compute
    ##  current_time_gyr and make a new combination snapshotdictionary 
    ##  that was interpolated.
    ##  TODO: think more about if this is how I want to do this if multiprocessing 
    ##  is turned off which should be the default mode tbh.
    this_snapdict = globals()[global_this_snapdict_name]

    frame_num = im_param_kwargs.pop('frame_num')
    ## initialize the GasStudio instance
    my_studio = this_class(
        datadir,
        snapnum,
        datadir,
        gas_snapdict=this_snapdict, ## pass in snapshot dictionary
        star_snapdict=None, ## TODO: have a flag to pass in snapshot dictionary
        **im_param_kwargs,
        loud=False,
        master_loud=False)

    ax, pixel_map =  my_studio.render(
        **render_kwargs) 

    del my_studio

    ax.axis('on')
    fig = ax.get_figure()
    fig.savefig('frame_%03d'%frame_num)
    plt.close(fig)
    del fig
    del ax
    return None


def main():
    
    init = time.time()
    #my_interp = Interpolater(
        #[40,20,40],
        #frame_half_widths=[30,5,5,30],
        #frame_centers = [[0,0,0],[5,5,0],[5,0,0],[0,0,0]],
        #thetas=[0,90,90,0],
        #frame_half_thicknesss=[15,5,5,15])
    my_interp = Interpolater(
        [300,150,300],
        thetas=[0,-90,-60,-60],
        psis=[0,720,720+360,720+360+720])
    my_interp.interpolateAndRenderMultiprocessing()
    print(time.time()-init,'s elapsed')

if __name__ == '__main__':
    main()


#chain = interpolationHelper(2,theta=(0,90))
#chain = interpolationHelper(2,chain,theta=45,psi=720)
#interpolationHelper(2,chain,frame_half_width=10)
#{'nsteps_tot': 90,
 #'nlinks': 3,
 #'nstepss': [30, 30, 30],
 #'thetas': [0, 90, 45, 45],
 #'psis': [0, 0, 720, 720],
 #'frame_half_widths': [15, 15, 15, 10]}
