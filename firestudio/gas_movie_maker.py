import numpy as np 
import os,sys,getopt
import copy


import matplotlib 
matplotlib.use('Agg') 

import matplotlib.pyplot as plt

from firestudio.utils.gas_utils.gas_utils import addPrettyGalaxyToAx

from abg_python.snapshot_utils import openSnapshot
from abg_python.cosmoExtractor import rotateVectorsZY,diskFilterDictionary
from abg_python.cosmo_utils import load_AHF

import multiprocessing

def renderGalaxy(
    ax,
    snapdir,snapnum,
    savefig=1,noaxis=0,mode='r',
    datadir=None,
    extract_galaxy=0,
    overwrite=0,
    **kwargs):
    """
    Input:
        ax - matplotlib axis object to draw to
        snapdir - location that the snapshots live in
        snapnum - snapshot number

    Optional:
        savefig=1 - flag to save figure to datadir (default snapdir, but can be a kwarg)
        noaxis=0 - flag to turn off axis (1=off 0=on)
        mode='r' - 'r' for reading from intermediate files, anything else to ignore intermediate files
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
        take_log_of_quantity=True - should we save the log of the quantity being plotted
            to the intermediate hdf5 file (to be subsequently plotted?)
    """

    ## default value for pixels, handled this way so we don't have to pass it around
    if 'pixels' not in kwargs:
        kwargs['pixels']=1200

    ## copy the dictionary so we don't mess anything up 
    copydict = copy.copy(kwargs)
    print(list(copydict.keys()),'keys passed to FIRE_studio')

    ## pass along input to next routines through copydict
    copydict['overwrite']=overwrite

    ## default to snapdir, pass datadir to next routines through copydict
    if 'datadir' is None:
        datadir = snapdir
    copydict['datadir']=datadir

    try:
        ## if we're being told to overwrite we shouldn't use the previous intermediate files
        assert not overwrite
        print("Trying to use a previous projection...")

        ## pass a dummy frame_center if not given one
        if 'frame_center' not in copydict:
            copydict['frame_center']=np.zeros(3)

        ax = addPrettyGalaxyToAx(
            ax,snapdir,snapnum, 
            **copydict)

    except (IOError,AssertionError,TypeError):
        print("Failed to use a previous projection")
        ## unpack the snapshot keys from snapdict or open the snapshot 
        ##  itself to load data into the copydict
        copydict.update(addSnapKeys(snapdir,snapnum,extract_galaxy,**copydict))
    
        ax = addPrettyGalaxyToAx(
            ax,snapdir,snapnum,
            **copydict)

    ## set the figure size appropriately, extended for the edgeon view
    if 'edgeon' in kwargs and kwargs['edgeon']:
        ax.get_figure().set_size_inches(6,8)
    else:
        ax.get_figure().set_size_inches(6,6)

    ## turn off the axis if asked
    if noaxis:
        ax.axis('off')

    ## save the figure if asked
    if savefig:
        savefig_args={} 
        if noaxis:
            ## remove whitespace around the axis, apparently the x/y origin is offset in pixel 
            ## space and so the bounding box doesn't actually reflect the left/bottom edge of the 
            ## axis
            ax.xaxis.set_major_locator(plt.NullLocator())
            ax.yaxis.set_major_locator(plt.NullLocator())
            savefig_args['bbox_inches']='tight'
            savefig_args['pad_inches']=0

        pixels = 1200 if 'pixels' not in kwargs else kwargs['pixels'] 
        image_name = "frame_%03d_%dkpc.png" % (snapnum, 2*kwargs['frame_half_width'])

        ax.get_figure().savefig(
            os.path.join(datadir,'Plots','GasTwoColour',image_name),dpi=300,
            **savefig_args)

    return ax 

def addSnapKeys(
    snapdir,snapnum,
    extract_galaxy=False,
    snapdict=None,
    frame_center=None,
    **kwargs):
    quantity_name = kwargs['quantity_name'] if ('quantity_name' in kwargs) else 'Temperature'

    ## have to go get the snapdict and add it to the copydict
    if snapdict is not None:
        pass
    elif not extract_galaxy:
        ## isolated galaxy huh? good choice. 
        ##  don't worry, cosmological will be overwritten in openSnapshot if 
        ##  HubbleParam != 1, none of the coordinates will be offset
        ##  and the galaxy won't be rotated or extracted
        snapdict=openSnapshot(snapdir,snapnum,ptype=0,cosmological=0)
    else:
        ## cosmological snapshot it is then... 
        snapdict=openSnapshot(snapdir,snapnum,ptype=0,cosmological=1)
        scom,rvir,vesc,rstar_half = load_AHF(
            snapdir,snapnum,
            snapdict['Redshift'],
            ahf_path=kwargs['ahf_path'] if 'ahf_path' in kwargs else None)

        ## filter all the keys in the snapdict as necessary to extract a spherical volume
        ##  centered on scom (the halo center), using 5*rstar_half  (thanks Zach for galaxy definition)
        diskFilterDictionary(
            None,snapdict,
            radius=rstar_half*5,scom=scom)

    Coordinates = snapdict['Coordinates']
    Masses = snapdict['Masses']

    ## temperature is computed in openSnapshot
    Quantity = snapdict[quantity_name]

    mydict = {
        'Coordinates':Coordinates,'Masses':Masses,'Quantity':Quantity,
        'BoxSize':snapdict['BoxSize'],'frame_center' : frame_center
    }
    return mydict

def multiProcRender(snapnum):
    ax = plt.gca()
    renderGalaxy(ax,glob_snapdir,snapnum,**glob_kwargs)
    plt.clf()
    
def main(snapdir,snapstart,snapmax,**kwargs):
    if 'multiproc' in kwargs and kwargs['multiproc']:
        ## map a wrapper to a pool of processes
        global glob_kwargs,glob_snapdir
        glob_kwargs = kwargs
        glob_snapdir=snapdir
        my_pool = multiprocessing.Pool(int(kwargs['multiproc']))
        my_pool.map(multiProcRender,range(snapstart,snapmax))
    else:
        ## just do a for loop
        for snapnum in range(snapstart,snapmax):
            ax = plt.gca()
            renderGalaxy(ax,snapdir,snapnum,**kwargs)
            plt.clf()       

if __name__=='__main__':
    argv = sys.argv[1:]
    opts,args = getopt.getopt(argv,'rs',[
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
    ])

    #options:
    # -r/s = use readsnap or use single snapshot loader
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

