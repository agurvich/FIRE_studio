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

def addSnapKeys(
    snapdir,snapnum,
    extract_galaxy=False,
    snapdict=None,
    frame_center=None,
    **kwargs):

    ## have to go get the snapdict and add it to the copydict
    if snapdict is not None:
        pass
    elif not extract_galaxy:
        ## isolated galaxy huh? good choice. 
        ##  don't worry, cosmological will be overwritten in openSnapshot if 
        ##  HubbleParam != 1, none of the coordinates will be updated
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

    pos_all = snapdict['Coordinates']
    mass_all = snapdict['Masses']

    ## temperature is computed in openSnapshot
    temperature_all = snapdict['Temperature']

    mydict = {
        'pos_all':pos_all,'mass_all':mass_all,'temperature_all':temperature_all,
        'HubbleParam':snapdict['HubbleParam'],'time_Myr':snapdict['Time'],
        'BoxSize':snapdict['BoxSize'],'frame_center' : frame_center
    }
    return mydict

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
            datadir=None - place to save intermediate files (None -> snapdir)
            extract_galaxy=False - flag to use abg_python.cosmoExtractor to extract main halo
            overwrite=0 -  flag to overwrite intermediate files
        Available kwargs:
            snapdict=None - snapshot dictionary of gas particles to use, will ignore extract_galaxy if present

            -- make image --
            theta=0 - euler rotation angle
            phi=0 - euler rotation angle
            psi=0 - euler rotation angle
            pixels=1200 - the resolution of image (pixels x pixels)
            min_den=-0.4 - the minimum of the log(density) color scale
            max_den=1.6 - the maximum of the log(density) color scale
            min_temp=2 - the minimum of the log(temperature) color scale
            max_temp=7 - the maximum of the log(temperature) color scale
            edgeon=0 - flag to create and plot an edgeon view

            frame_center=None - origin of image in data space, if None will use [0,0,0]
            frame_half_width=None - half-width of image in data space, if None will use ? 
            frame_depth=None - half-depth of image in data space, if None will use ? 

            fontsize=None - size of font on the scale bar and the text label
            scale_bar=1 - flag to include a scale bar
            figure_label="" - text to display in the upper right corner, defaults to current time?
    """
    ## copy the dictionary so we don't mess anything up 
    copydict = copy.copy(kwargs)
    print copydict.keys(),'keys passed'

    ## pass along input to next routines through copydict
    copydict['overwrite']=overwrite

    ## default to snapdir, pass datadir to next routines through copydict
    if 'datadir' is None:
        datadir = snapdir
    copydict['datadir']=datadir

    try:
        ## if we're being told to overwrite we shouldn't use the previous intermediate files
        assert not overwrite
        print "Trying to use a previous projection..."

        ## pass a dummy frame_center, it's not used but will raise an error
        if 'frame_center' not in copydict:
            copydict['frame_center']=np.zeros(3)

        ax = addPrettyGalaxyToAx(
            ax,snapdir,snapnum, 
            **copydict)

    except (IOError,AssertionError,TypeError):
        print "Failed to use a previous projection"
        ## add the snapshot keys to the copydict
        copydict.update(addSnapKeys(snapdir,snapnum,extract_galaxy,**copydict))
    
        ax = addPrettyGalaxyToAx(
            ax,snapdir,snapnum,
            **copydict)

    ## set the figure size appropriately, extended for the edgeon view
    if 'edgeon' in kwargs and kwargs['edgeon']:
        ax.get_figure().set_size_inches(6,8)
    else:
	# TODO 8x8, 6x6 is normal
	print "SETTING THE SIZE"
        ax.get_figure().set_size_inches(8,8)

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
            os.path.join(datadir,'Plots','GasTwoColour',image_name),dpi=387,
            **savefig_args)

    return ax 

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
        for snapnum in xrange(snapstart,snapmax):
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
        'min_den=','max_den=','min_temp=','max_temp=','datadir=',
        'noaxis=',
        'multiproc=',
        'extract_galaxy=',
        'ahf_path=',
	'figure_label='])

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
                value = eval(opt[1])
            opts[i]=(key,value)
    main(**dict(opts))

