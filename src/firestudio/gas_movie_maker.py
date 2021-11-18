import numpy as np 
import os,sys,getopt
import copy
import multiprocessing
import itertools
import gc

from .studios.gas_studio import GasStudio

def renderGalaxy(
    ax,
    snapdir,snapnum,
    datadir,
    frame_half_width,frame_half_thickness,
    edgeon = 0,
    **kwargs):
    ## copy the docstring
    GasStudio.__doc__
    
    return render(
        snapdir,snapnum,
        datadir,
        frame_half_width,frame_half_thickness,
        edgeon,kwargs,ax)

def render(
    snapdir,
    snapnum,
    datadir,
    frame_half_width,
    frame_half_thickness,
    edgeon,
    min_weight,
    max_weight,
    min_quantity,
    max_quantity,
    ax,
    kwargs):

    firestudio_datadir = os.path.join(datadir,'firestudio')

    ## assumes that the simulation name immediately precedes "output"
    sim_name = snapdir.split(os.path.sep)[-2]

    fname = 'snapshot_%03d.png'%snapnum
    if edgeon:
        fname = 'edgeon_'+fname

    fname = os.path.join(
            firestudio_datadir,
            'plots',
            fname)

    if os.path.isfile(fname):
        return

    if ax is None:
        ax = plt.gca()
    fig = ax.get_figure()

    gasStudio = GasStudio(
        sim_name=sim_name,
        snapdir=snapdir,
        snapnum=snapnum,
        datadir=firestudio_datadir,
        frame_half_width=frame_half_width,
        frame_half_thickness=frame_half_thickness,
        **kwargs)

    if edgeon:
        gasStudio.set_ImageParams(
            theta = 90,
            aspect_ratio = frame_half_thickness/frame_half_width)

    try:
        gasStudio.render(
            ax,
            min_weight=min_weight,
            max_weight=max_weight,
            min_quantity=min_quantity,
            max_quantity=max_quantity,
            weight_adjustment_function= lambda x: np.log10(x/gasStudio.Acell)+10-6, ## log10(msun/pc^2)
            quantity_adjustment_function=np.log10)
    except:
        gasStudio.load_SnapshotData(use_saved_subsnapshots=True)
        gasStudio.render(
            ax,
            min_weight=min_weight,
            max_weight=max_weight,
            min_quantity=min_quantity,
            max_quantity=max_quantity,
            weight_adjustment_function= lambda x: np.log10(x/gasStudio.Acell)+10-6, ## log10(msun/pc^2)
            quantity_adjustment_function=np.log10)

    fig.savefig(
        fname,
        dpi=150,
        pad_inches=0,
        bbox_inches='tight')

    ## don't remove these lines, they perform
    ##  some kind of dark arts that helps the
    ##  garbage collecter collect
    del gasStudio
    locals().keys()
    globals().keys()
    gc.collect()

    return ax
    
def main(
    snapdir,
    snapstart,
    snapmax,
    datadir,
    frame_half_width,
    frame_half_thickness,
    edgeon=0,
    min_weight=-0.5,
    max_weight=1.6,
    min_quantity=2,
    max_quantity=7, 
    multiproc=1,
    **kwargs):

    argss = zip(
        itertools.repeat(snapdir),
        range(snapstart,snapmax+1),
        itertools.repeat(datadir),
        itertools.repeat(frame_half_width),
        itertools.repeat(frame_half_thickness),
        itertools.repeat(edgeon),
        itertools.repeat(min_weight),
        itertools.repeat(max_weight),
        itertools.repeat(min_quantity),
        itertools.repeat(max_quantity),
        itertools.repeat(None),
        itertools.repeat(kwargs))

    if multiproc > 1:
        ## map a wrapper to a pool of processes
        with multiprocessing.Pool(multiproc) as my_pool:
            my_pool.starmap(render,argss)
    else:
        ## just do a for loop
        for args in argss:
            render(*args)

if __name__=='__main__':
    import matplotlib.pyplot as plt
    argv = sys.argv[1:]
    opts,args = getopt.getopt(argv,'',[
        'frame_half_width=',
        'frame_half_thickness=',
        'pixels=',
        'edgeon=',
        'noaxis=',

        'min_den=','max_den=',
        'min_temp=','max_temp=', 
        'single_image=',
        'use_colorbar=', 
        'cmap=',
        'cbar_label=',
        'take_log_of_quantity=',

        'snapdir=',
        'datadir=',
        'snapstart=',
        'snapmax=',
        'multiproc=',])

    #options:
    #--min/max_den/temp: bottom/top of color scales for density/temperature
    #--cmap : string of colormap name to use
    #--single_image : string of quantity name that you want to make a one-color mass weighted map of
    #--use_colorbar : flag to create a colorbar
    #--cbar_label :  flag to label the colorbar
    #--take_log_of_quantity : flag to take the log of the quantity you are making a map of

    for i,opt in enumerate(opts):
        if opt[1]=='':
            opts[i]=('mode',opt[0].replace('-'))
        else:
            key = opt[0].replace('-','')
            if key in ['snapdir','datadir']:
                value= opt[1]
            else:
                # turn arguments from strings to whatever
                try:
                    value = eval(opt[1])
                except:
                    value = opt[1]
            opts[i]=(key,value)
    main(**dict(opts))
