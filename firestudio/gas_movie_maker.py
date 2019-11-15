import numpy as np 
import os,sys,getopt
import copy
import multiprocessing
import itertools

from firestudio.studios.gas_studio import GasStudio

def renderGalaxy(
    ax,
    snapdir,snapnum,
    datadir,
    frame_half_width,frame_depth,
    edgeon = 0,
    **kwargs):
    ## copy the docstring
    GasStudio.__doc__
    
    return render(
        snapdir,snapnum,
        datadir,
        frame_half_width,frame_depth,
        edgeon,kwargs,ax)

def renderWrapper(args):
    return render(*args)

def render(
    snapdir,snapnum,
    datadir,
    frame_half_width,frame_depth,
    edgeon,
    kwargs,
    ax=None,
    image_names = None):

    ## this is stupid and should change
    image_names = [
        'columnDensityMap',
        'massWeightedTemperatureMap',
        'two_color']

    gasStudio = GasStudio(
        snapdir=snapdir,
        snapnum=snapnum,
        datadir=datadir,
        frame_half_width=frame_half_width,
        frame_depth=frame_depth,
        **kwargs)

    if edgeon:
        gasStudio.renderFaceAppendEdgeViews(image_names)
    else:
        if ax is None:
            ax = plt.gca()
            ## set figure size to square
            ax.get_figure().set_size_inches(6,6)
        gasStudio.render(ax,image_names)
        if ax is None:
            plt.clf()
    return ax
    
def main(
    snapdir,
    snapstart,snapmax,
    datadir,
    frame_half_width,
    frame_depth,
    edgeon=0,
    **kwargs):

    if 'multiproc' in kwargs and kwargs['multiproc']:
        ## map a wrapper to a pool of processes
        argss = itertools.izip(
            itertools.repeat(snapdir),
            range(snapstart,snapmax+1),
            itertools.repeat(datadir),
            itertools.repeat(frame_half_width),
            itertools.repeat(frame_depth),
            itertools.repeat(edgeon),
            itertools.repeat(kwargs),
            itertools.repeat(None))
        my_pool = multiprocessing.Pool(int(kwargs['multiproc']))
        my_pool.map(renderWrapper,argss)
    else:
        ## just do a for loop
        for snapnum in range(snapstart,snapmax+1):
            render(
                snapdir,snapnum,
                datadir,
                frame_half_width,frame_depth,
                edgeon,
                kwargs,
                None)

if __name__=='__main__':
    import matplotlib.pyplot as plt
    from firestudio.studios.studio import shared_kwargs
    argv = sys.argv[1:]
    opts,args = getopt.getopt(argv,'',[
        'min_den=','max_den=',
        'min_temp=','max_temp=', 
        'single_image=',
        'use_colorbar=', 
        'cmap=',
        'cbar_label=',
        'take_log_of_quantity=',
    ]+shared_kwargs)

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
