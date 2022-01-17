import numpy as np 
import os,sys,getopt
import itertools
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt

from .studios.star_studio import StarStudio

def renderStarGalaxy(
    ax,
    snapdir,snapnum,
    datadir,
    frame_half_width,frame_depth,
    edgeon=False,
    **kwargs):
    ## copy the docstring
    StarStudio.__doc__

    ## make the call to renderWrapper
    return render(snapdir,snapnum,datadir,frame_half_width,frame_depth,edgeon,kwargs,ax)

######## Meat and Potatoes render looping functions
def renderWrapper(args):
    return render(*args)

def render(
    snapdir,snapnum,
    datadir,
    frame_half_width,frame_depth,
    edgeon,
    kwargs,
    ax=None):

    image_names = ['out_u','oul_g','out_r','hubble']
    starStudio = StarStudio(
        snapdir=snapdir,
        snapnum=snapnum,
        datadir=datadir,
        frame_half_width=frame_half_width,
        frame_depth=frame_depth,
        **kwargs)

    if edgeon:
        starStudio.renderFaceAppendEdgeViews(image_names)
    else:
        if ax is None:
            ax = plt.gca()
            ## set figure size to square
            ax.get_figure().set_size_inches(6,6)
        starStudio.render(ax,image_names)
        if ax is None:
            plt.clf()

    return starStudio
    
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
            render(snapdir,snapnum,datadir,frame_half_width,frame_depth,edgeon,kwargs,None)

if __name__=='__main__':
    import matplotlib.pyplot as plt
    from firestudio.studios.studio import shared_kwargs
    argv = sys.argv[1:]
    opts,args = getopt.getopt(argv,'',[
        'dynrange=','maxden=','color_scheme_nasa='] + shared_kwargs)    



    # options:
    # --dynrange : TODO unknown
    # --maxden : TODO unknown
    # --color_scheme_nasa: True - flag to use nasa colors (vs. SDSS if false) 

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
