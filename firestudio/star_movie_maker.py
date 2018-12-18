import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np 
import h5py
import os,sys,getopt

from firestudio.utils.stellar_utils import raytrace_projection,load_stellar_hsml
import firestudio.utils.stellar_utils.make_threeband_image as makethreepic

from abg_python.all_utils import filterDictionary
from abg_python.snapshot_utils import openSnapshot
from abg_python.cosmo_utils import load_AHF
from abg_python.cosmoExtractor import diskFilterDictionary as orientFaceon



## Stellar light attenuation projection
def calc_stellar_hsml(x,y,z):
    starhsml = load_stellar_hsml.get_particle_hsml
    h_all = starhsml(x,y,z)
    return h_all

def raytrace_ugr_attenuation(
    x,y,z,
    mstar,ages,metals,
    h_star, 
    gx,gy,gz,
    mgas, gas_metals,
    h_gas,
    pixels = 1200,
    xlim = None, ylim = None, zlim = None
    ):

    ## setup boundaries to cut-out gas particles that lay outside
    ## range
    if xlim is None:
        xlim = [np.min(x),np.max(x)]
    if ylim is None:
        ylim = [np.min(y),np.max(y)]
    if zlim is None:
        zlim = [np.min(z),np.max(z)]

    stellar_raytrace = raytrace_projection.stellar_raytrace


#   band=BAND_ID; # default=bolometric
#   j = [  0,  6,  7,  8,  9, 10, 11, 12, 13,  1,   2,   3,   4,   5] # ordering I'm used to
#   i = [  0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10,  11,  12,  13] # ordering of this
#   band_standardordering = band
#   band = j[band]
#   if (band > 13): 
#       print 'BAND_ID must be < 13'; 
#       return 0;
#   
#   b=['Bolometric', \
#   'Sloan u','Sloan g','Sloan r','Sloan i','Sloan z', \
#   'Johnsons U','Johnsons B', 'Johnsons V','Johnsons R','Johnsons I', \
#   'Cousins J','Cousins H','Cousins K']
    ## pick color bands by their IDs, see above
    BAND_IDS=[9,10,11]
    #gas_out,out_u,out_g,out_r = stellar_raytrace(
    return stellar_raytrace(
        BAND_IDS,
        x,y,z,
        mstar,ages,metals,
        h_star,
        gx,gy,gz,
        mgas, gas_metals,
        h_gas,
        pixels=pixels,
        xlim=xlim,ylim=ylim,zlim=zlim,
        ADD_BASE_METALLICITY=0.001*0.02,ADD_BASE_AGE=0.0003,
        IMF_SALPETER=0,IMF_CHABRIER=1
    )

def make_threeband_image(ax,
    out_r,out_g,out_u,
    maxden = 1e-3, dynrange = 1.001,
    pixels = 1200):

    image24, massmap = makethreepic.make_threeband_image_process_bandmaps(
        out_r,out_g,out_u,
        maxden=maxden,dynrange=dynrange,pixels=pixels,
        color_scheme_nasa=1,color_scheme_sdss=0)

    ## for some reason it's rotated 90 degrees...? kind of like transposed but different
    ##  need to take the rows of the output image and make them the columns, iteratively,
    ##  for now... 
    #image24=np.rot90(image24,k=1,axes=(0,1))
    return np.transpose(image24,axes=(1,0,2))
    
def get_h_star(xs,ys,zs,savefile=None,overwrite=0):
    try:
        assert not overwrite
        if savefile is not None:
            with h5py.File(savefile,'r') as handle:
                return np.array(handle['StarHs'])
        else:
            raise IOError
    except (IOError,KeyError,AssertionError):
        h_star =  calc_stellar_hsml(xs,ys,zs)
        if savefile is not None:
            ## save it for next time!
            mode = 'a' if not overwrite else 'w'
            with h5py.File(savefile,mode) as handle:
                handle['StarHs']=h_star
        return h_star

def get_bands_out(
    xs,ys,zs,
    mstar,ages,metals,
    h_star,
    gxs,gys,gzs,
    mgas,gas_metals,h_gas,
    pixels = 1200,
    savefile=None
    ):
    try:
        if savefile is not None:
            with h5py.File(savefile,'r') as handle:
                return np.array(handle['out_u']),np.array(handle['out_g']),np.array(handle['out_r'])
        else:
            raise IOError
    except (IOError,KeyError):
        gas_out,out_u,out_g,out_r = raytrace_ugr_attenuation(
            xs,ys,zs,
            mstar,ages,metals,
            h_star,
            gxs,gys,gzs,
            mgas,gas_metals,h_gas,
            pixels=pixels
        ) 
        if savefile is not None:
            ## save it for next time!
            with h5py.File(savefile,'a') as handle:
                handle['out_u']=out_u
                handle['out_g']=out_g
                handle['out_r']=out_r

        return out_u,out_g,out_r

def add_scale_bar(final_image,frame_half_width,npix,scale_line_length):
    # Convert to pixels
    length_per_pixel = 2.0*frame_half_width/ npix
    scale_line_length_px = int(scale_line_length / length_per_pixel)

    # Position in terms of image array indices
    scale_line_x_start = int(0.05 * npix)
    scale_line_x_end = min(scale_line_x_start + scale_line_length_px,npix)
    scale_line_y = int(0.02 * npix)

    # Go through pixels for scale bar, setting them to white
    for x_index in xrange(scale_line_x_start, scale_line_x_end):
        final_image[scale_line_y, x_index, 0] = 1
        final_image[scale_line_y, x_index, 1] = 1
        final_image[scale_line_y, x_index, 2] = 1
        final_image[scale_line_y + 1, x_index, 0] = 1
        final_image[scale_line_y + 1, x_index, 1] = 1
        final_image[scale_line_y + 1, x_index, 2] = 1
        final_image[scale_line_y + 2, x_index, 0] = 1
        final_image[scale_line_y + 2, x_index, 1] = 1
        final_image[scale_line_y + 2, x_index, 2] = 1
        final_image[scale_line_y + 3, x_index, 0] = 1
        final_image[scale_line_y + 3, x_index, 1] = 1
        final_image[scale_line_y + 3, x_index, 2] = 1
        final_image[scale_line_y + 4, x_index, 0] = 1
        final_image[scale_line_y + 4, x_index, 1] = 1
        final_image[scale_line_y + 4, x_index, 2] = 1
        final_image[scale_line_y + 5, x_index, 0] = 1
        final_image[scale_line_y + 5, x_index, 1] = 1
        final_image[scale_line_y + 5, x_index, 2] = 1

    return final_image

def get_indices(
    snap,
    frame_center,
    frame_half_width,frame_depth=None):

    ## offset, unpack, and square coordinates
    xs2,ys2,zs2 = (snap['Coordinates']-frame_center).T**2

    ## compare to half-width squared
    xindices = xs2 <= (frame_half_width)**2
    yindices = ys2 <= (frame_half_width)**2
    zindices = zs2 <= (frame_depth)**2

    return np.logical_and(np.logical_and(xindices,yindices),zindices)

def renderStarGalaxy(
    ax,
    snapdir,snapnum,
    datadir,
    frame_half_width,frame_depth,
    frame_center=np.zeros(3),
    savefig=1,noaxis=0,
    savefile=None,mode='r',overwrite=0,
    fontsize=None,scale_bar=1,
    extract_galaxy = 1,
    dynrange = 10,
    maxden = 1e-2,
    pixels=1200,
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
            extract_galaxy=True - flag to use abg_python.cosmoExtractor to extract main halo
            overwrite=0 -  flag to overwrite intermediate files
            frame_half_width=None - half-width of image in data space, if None, will use width of data provided
            frame_depth=None - half-depth of image in data space, if None will use frame_half_width
            frame_center=np.zeros(3) - origin of image in data space 
            fontsize=None - size of font on the scale bar and the text label
            scale_bar=1 - flag to include a scale bar
            dynrange=10 - something to do with the dynamic range in the image
            maxden=1e-2 - something to do with the colors in the image
        Available kwargs:
            -- make image --
            pixels=1200 - the resolution of image (pixels x pixels)
            #TODO
            theta=0 - euler rotation angle
            phi=0 - euler rotation angle
            psi=0 - euler rotation angle
            edgeon=0 - flag to create and plot an edgeon view
    """
    ## handle default arguments
    datadir = snapdir if datadir is None else datadir
    if savefile is None:
        savefile = os.path.join(datadir,'Plots','Projections','hubble_image_%03d.hdf5'%snapnum)
    ## force no savefile to be created, need to manually pass this
    elif savefile is '':
        savefile = None

    ## we've already been passed open snapshot data
    if ('star_snap' in kwargs) and ('snap' in kwargs):
        star_snap = kwargs['star_snap']
        snap = kwargs['snap']
    else:
        snap,star_snap = loadSnapshotData(extract_galaxy,snapdir,snapnum)
    
    ## now find only the particles within the viewbox
    ##  apply indices to all the arrays in the snapdict
    indices = get_indices(star_snap,frame_center,frame_half_width,frame_depth)
    star_snap = filterDictionary(star_snap,indices)

    indices = get_indices(star_snap,frame_center,frame_half_width,frame_depth)
    snap = filterDictionary(snap,indices)
    image24 = processSnapshots(
        ax,snap,star_snap,dynrange,maxden,savefile,pixels=pixels,mode=mode,overwrite=overwrite)
    addScaleBarAndFigureLabel(ax,image24,frame_half_width,fontsize,pixels,**kwargs)

def processSnapshots(ax,snap,star_snap,dynrange,maxden,savefile,pixels=1200, mode='r',overwrite=0):
    ## unpack relevant information
    xs,ys,zs = star_snap['Coordinates'].T
    mstar,ages, metals = star_snap['Masses'],star_snap['AgeGyr'],star_snap['Metallicity'][:,0]
    gxs, gys, gzs, = snap['Coordinates'].T
    mgas,gas_metals, h_gas = snap['Masses'],snap['Metallicity'][:,0],snap['SmoothingLength']

    h_star = get_h_star(xs,ys,zs,savefile if mode =='r' else None,overwrite=overwrite)
    
    out_u,out_g,out_r=get_bands_out(
        xs,ys,zs,
        mstar,ages,metals,
        h_star,
        gxs,gys,gzs,
        mgas,gas_metals,h_gas,
        pixels = pixels,
        savefile=savefile if mode =='r' else None
        )

    image24=make_threeband_image(
        ax,
        out_r,out_g,out_u,
        dynrange=dynrange,
        maxden=maxden)

    return image24

def addScaleBarAndFigureLabel(
    ax,image24,frame_half_width,fontsize,
    pixels,scale_bar=1,noaxis=1,savefig=0,**kwargs):
    if 2*frame_half_width > 15 : 
        scale_line_length = 5
        scale_label_text = r"$\mathbf{5 \, \rm{kpc}}$"

    elif 2*frame_half_width > 1.5 : 
        scale_line_length = 1.0 
        scale_label_text = r"$\mathbf{1 \, \rm{kpc}}$"

    else:
        scale_line_length = .1
        scale_label_text = r"$\mathbf{100 \, \rm{pc}}$"

    ## handle default
    fontsize=8 if fontsize is None else fontsize
    scale_label_position = 0.06 

    if scale_bar:
        image24=add_scale_bar(image24,frame_half_width,pixels,scale_line_length)
        label2 = ax.text(scale_label_position,
            0.03, scale_label_text, fontweight = 'bold', transform = ax.transAxes)
        label2.set_color('white')
        label2.set_fontsize(fontsize*0.75)

    ax.imshow(image24,interpolation='bicubic',origin='lower')#,aspect='normal')
    ax.get_figure().set_size_inches(6,6)
    if noaxis:
        ax.axis('off')
    if savefig:
        savename = os.path.join(datadir,'hubble_image_%03d'%snapnum)
        savefig_args={} 
        if noaxis:
            ## remove whitespace around the axis, apparently the x/y origin is offset in pixel 
            ##  space and so the bounding box doesn't actually reflect the left/bottom edge of the 
            ##  axis
            ax.xaxis.set_major_locator(plt.NullLocator())
            ax.yaxis.set_major_locator(plt.NullLocator())
            savefig_args['bbox_inches']='tight'
            savefig_args['pad_inches']=0

        ax.get_figure().savefig(savename,
            **savefig_args)
        
    def loadSnapshotData(extract_galaxy,snapdir,snapnum):
        try:
            assert extract_galaxy
        except:
            raise Exception("Fullbox render from the command line isn't supported yet!")

        print "Assuming all stars are part type 4, hope this isn't an isolated galaxy!"
        ## load star particles
        star_snap = openSnapshot(
            snapdir,snapnum,4,
            cosmological=1)
            #,keys_to_extract=['Coordinates','Masses','AgeGyr','Metallicity','Velocities'])

        ## load gas particles
        snap = openSnapshot(
            snapdir,snapnum,0,
            cosmological=1)
            #,keys_to_extract=[
            #   'Coordinates','Masses','Metallicity','Velocities',
            #   'SmoothingLength','InternalEnergy','ElectronAbundance',
            #   'Density'])

        ## load the center of the halo
        scom,rvir,vesc = load_AHF(snapdir,snapnum,snap['Redshift'],
            ahf_path=kwargs['ahf_path'] if 'ahf_path' in kwargs else None,
            extra_names_to_read = [])

        ## overwrite star_snap/snap with a subset that contains only within 0.2 rvir
        ##  of the halo center-- or within a frame_half_width radius, whichever is larger
        ##  orient_stars = 0 -> find angular momentum vector of gas and make that new z axis
        radius = 0.2*rvir if frame_half_width is None else np.max([0.2*rvir,frame_half_width]) 
        snap,star_snap = orientFaceon(star_snap,snap,
            radius = radius, ## even if frame_half_width is None this works, huh..
            scom = scom, orient_stars = 0)
        return snap,star_snap

def multiProcRender(snapnum):
    ax = plt.gca()
    renderStarGalaxy(ax,glob_snapdir,snapnum,**glob_kwargs)
    plt.clf()

def main(snapdir,snapstart,snapmax,datadir=None,**kwargs):

    ## handle datadir creation before possible multiproc
    datadir = snapdir if datadir is None else datadir
    kwargs['datadir']=datadir
    ## creates plots/twocolour and plots/projections/den and plots/projections/temp
    makeOutputDirectories(datadir)

    if 'multiproc' in kwargs and kwargs['multiproc']:
        ## map a wrapper to a pool of processes
        global glob_kwargs,glob_snapdir
        glob_kwargs = kwargs
        glob_snapdir=snapdir
        my_pool = multiprocessing.Pool(int(kwargs['multiproc']))
        my_pool.map(multiProcRender,range(snapstart,snapmax+1))
    else:
        ## just do a for loop
        for snapnum in xrange(snapstart,snapmax+1):
            ax = plt.gca()
            renderStarGalaxy(ax,snapdir,snapnum,**kwargs)
            plt.clf()       

if __name__=='__main__':
    argv = sys.argv[1:]
    opts,args = getopt.getopt(argv,'rs',[
        'snapdir=',
        'snapstart=','snapmax=',
        'pixels=','frame_half_width=','frame_depth=',
        'theta=','phi=','psi=',
        'dynrange=',
        'maxden=',
        'edgeon=',
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
    # -r/s = use readsnap or use single snapshot loader
    #--snapdir: place where snapshots live
    #--snapstart : which snapshot to start the loop at
    #--snapmax : which snapshot to end the loop at
    #--frame_half_width : half width of frame in kpc
    #--frame_depth : half depth of frame in kpc
    #--datadir: place to output frames to

    #--theta,phi,psi : euler angles for rotation

    #--dynrange : TODO unknown
    #--maxden : TODO unknown

    #--edgeon : flag for sticking a 90 degree edge on rotation underneath 
    #--pixels : how many pixels in each direction to use, defaults to 1200
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


