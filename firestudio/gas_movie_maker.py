import numpy as np 
import os,sys,getopt
import copy



import matplotlib 
matplotlib.use('Agg') 

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from firestudio.utils.gas_utils.gas_utils import addPrettyGalaxyToAx

from abg_python.snapshot_utils import openSnapshot
from abg_python.cosmoExtractor import rotateVectorsZY,diskFilterDictionary
from abg_python.cosmo_utils import load_AHF

import multiprocessing

def renderGalaxy(
    ax,
    snapdir,snapnum,
    savefig=1,noaxis=0,
    datadir=None,
    extract_galaxy=0,
    overwrite=0,
    edgeon = 0,
    **kwargs):
    """
    Input:
        ax - matplotlib axis object to draw to
        snapdir - location that the snapshots live in
        snapnum - snapshot number

    Optional:

        extract_galaxy=False - flag to use abg_python.cosmoExtractor to extract main halo

        overwrite=0 - flag to overwrite the intermediate grid files
        datadir=None - place to save intermediate files (None -> snapdir)

    Mandatory kwargs to be passed along:

        snapdict=None - snapshot dictionary of gas particles to use,
            will ignore extract_galaxy if present

            Required elements in snapdict:
                Coordinates - coordinates of particles to be projected, in kpc
                Masses - masses of particles to be projected, in 1e10 msun
                Quantity - quantity of particles to be mass weighted/projected
        
                BoxSize - c routine needs it, probably fine to pass in a large number

        frame_center - origin of image in data space 
        frame_half_width - half-width of image in data space
        frame_depth - half-depth of image in data space 

    Optional kwargs to be passed along: 
        quantity_name='Temperature' - the name of the quantity that you're mass weighting
            should match whatever array you're passing in as quantity

        theta=0- euler rotation angle
        phi=0- euler rotation angle
        psi=0 - euler rotation angle
        aspect_ratio - the 'shape' of the image (y/x)
        pixels=1200 - the resolution of image (pixels x pixels)
        min_den=-0.4 - the minimum of the density color scale
        max_den=1.6 - the maximum of the density color scale
        min_quantity=2 - the minimum of the temperature color scale
        max_quantity=7 - the maximum of the temperature color scale

        h5prefix='' - a string that you can prepend to the projection filename if desired
        this_setup_id=None - string that defines the projection setup, None by default means
            it defaults to a gross combination of frame params + angles
        cmap='viridis' - string name for cmap to use 
        scale_bar=1 - should you plot a scale bar in the bottom left corner
        figure_label=None - what string should you put in the top right corner? 
        fontsize=None - fontsize for all text in frame
        single_image=None - string, if it's "Density" it will plot a column 
            density projection, if it's anything else it will be a mass weighted
            `quantity_name` projection. None will be a "two-colour" projection
            with hue determined by `quantity_name` and saturation by density
        use_colorbar=False - flag for whether to plot a colorbar at all
        cbar_label=None - string that shows up next to colorbar, good for specifying 
            units, otherwise will default to just quantity_name.title()
        take_log_of_quantity=True - should we plot the log of the resulting quantity map?
    """


def renderWrapper(snapnum,image_names):
    gasStudio = GasStudio(
        snapdir=glob_snapdir,
        snapnum=snapnum,
        **glob_kwargs)

    if glob_edgeon:
        gasStudio.renderFaceAppendEdgeViews(image_names)
    else:
        ax = plt.gca()
        ## set figure size to square
        ax.get_figure().set_size_inches(6,6)
        gasStudio.render(ax,image_names)
        plt.clf()
    
def main(snapdir,snapstart,snapmax,edgeon=0,**kwargs):
    global glob_kwargs,glob_snapdir,glob_edgeon
    glob_kwargs = kwargs
    glob_snapdir = snapdir
    glob_edgeon = edgeon
        
    if 'multiproc' in kwargs and kwargs['multiproc']:
        ## map a wrapper to a pool of processes
        my_pool = multiprocessing.Pool(int(kwargs['multiproc']))
        my_pool.map(renderWrapper,range(snapstart,snapmax))
    else:
        ## just do a for loop
        for snapnum in range(snapstart,snapmax+1):
            renderWrapper(snapnum)

if __name__=='__main__':
    argv = sys.argv[1:]
    opts,args = getopt.getopt(argv,'',[
        'snapdir=',
        'snapstart=','snapmax=',
        'pixels=','frame_half_width=','frame_depth=',
        'theta=','phi=','psi=',
        'edgeon=',
        'min_den=','max_den=',
        'min_temp=','max_temp=',
        'datadir=',
        'noaxis=',
        'multiproc=',
        'extract_galaxy=',
        'ahf_path=',
        'figure_label=',
        'cmap=',
        'single_image=',
        'use_colorbar=', 
        'cbar_label=',
        'take_log_of_quantity=',
        'overwrite=',
    ])

    #options:
    #--snapdir: place where snapshots live
    #--snapstart : which snapshot to start the loop at
    #--snapmax : which snapshot to end the loop at
    #--frame_half_width : half width of frame in kpc
    #--frame_depth : half depth of frame in kpc
    #--datadir: place to output frames to

    #--theta,phi,psi : euler angles for rotation
    #--edgeon : flag for sticking a 90 degree edge on rotation underneath 
    #--pixels : how many pixels in each direction to use, defaults to 1200
    #--min/max_den/temp: bottom/top of color scales for density/temperature
    #--noaxis : flag for removing axis and whitespace for just the pretty part
    #--multiproc : how many processes should be run simultaneously, keep in mind memory constraints
    #--extract_galaxy=False : flag to use abg_python.cosmoExtractor to extract main halo
    #--ahf_path : path relative to snapdir where the halo files are stored
    #--figure_label: text to put in the upper right hand corner

    #--cmap : string of colormap name to use
    #--single_image : string of quantity name that you want to make a one-color mass weighted map of
    #--use_colorbar : flag to create a colorbar
    #--cbar_label :  flag to label the colorbar
    #--take_log_of_quantity : flag to take the log of the quantity you are making a map of

    #--overwrite: flag to  overwrite the cached projection if it exists


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

