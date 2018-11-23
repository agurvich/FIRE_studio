from firestudio.utils.gas_utils.plotTwoColorGrid import plot_image_grid
from firestudio.utils.gas_utils.projectDensityAndQuantity import compute_image_grid 

from abg_python.snapshot_utils import openSnapshot

import matplotlib 
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
import numpy as np 
import os
import h5py

## system functions
def makeOutputDirectories(datadir):
    ## make the path to store all plots in the datadir
    path = os.path.join(datadir,'Plots')
    if not os.path.isdir(path):
        os.mkdir(path)

    ## make the path to store final images
    if not os.path.isdir(os.path.join(path,'GasTwoColour')):
        os.mkdir(os.path.join(path,'GasTwoColour'))

    ## make the paths to store projection hdf5 files
    path = os.path.join(path,'Projections')
    if not os.path.isdir(path):
        os.mkdir(path)

def checkProjectionFile(
    projection_file,
    pixels, frame_half_width,frame_depth,
    frame_center,
    theta=0,phi=0,psi=0,
    aspect_ratio=1,
    **kwargs):
    try:
        with h5py.File(projection_file,'r') as handle:
            for group in handle.keys():
                this_group = handle[group]
                flag = True
                print(this_group['aspect_ratio'].value,aspect_ratio)
                key,variable = 'aspect_ratio',aspect_ratio
                print(np.round(this_group[key].value,decimals=2) , np.round(variable,decimals=2))
                for key,variable in zip(
                    ['npix_x','frame_half_width','frame_depth',
                    'frame_center','theta','phi','psi','aspect_ratio'],
                    [ pixels , frame_half_width , frame_depth ,
                     frame_center , theta , phi , psi , aspect_ratio ]):

                    ## read the value in the hdf5 file and compare to variable
                    if key not in ['npix_x']:
                        ## key is not an integer so we have to round it somehow
                        flag = flag and np.all(
                            np.round(this_group[key].value,decimals=2) == np.round(variable,decimals=2))
                    else:
                        ## key is an integer
                        flag = flag and this_group[key].value == variable
                ## found the one we wanted
                if flag:
                    return 1 
        return 0 
    except IOError:
        return 0
    
## physics functions
def projectColumnDensityAndQuantity(
    snapdir,snapnum,
    projection_dir,**kwargs):
    """ 
    Input: 
        snapdir - directory the snapshots live in
        snapnum - snapshot number currently rendering
        projection_dir - place to save projections
        available kwargs:
            theta=0- euler rotation angle
            phi=0- euler rotation angle
            psi=0 - euler rotation angle
            pixels=1200 - the resolution of image (pixels x pixels)

            frame_center - origin of image in data space
            frame_half_width - half-width of image in data space
            frame_depth - half-depth of image in data space
    """

    print('projecting the image grids')
    den_grid=compute_image_grid(
        projection_dir = projection_dir,snapnum = snapnum,
        **kwargs)

def addPrettyGalaxyToAx(
    ax,
    snapdir,snapnum,
    overwrite=0,datadir=None,
    **kwargs):
    """
    Input:
        ax - matplotlib axis object to draw to
        snapdir - location that the snapshots live in
        snapnum - snapshot number

    Optional:
        overwrite=0 - flag to overwrite the intermediate grid files
        datadir=None - directory to output the the intermediate grid files and output png

    Mandatory kwargs to be passed along:
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

    print('Drawing %s:%d'%(snapdir,snapnum)+' to:%s'%datadir)
    makeOutputDirectories(datadir)

    ## where to find/save column density/quantity maps-- this could get crowded!
    projection_dir=os.path.join(datadir,'Plots','Projections')

    ## what are we going to call the intermediate filename? 
    ##	will it have a prefix? 
    h5prefix='' if 'h5prefix' not in kwargs else kwargs['h5prefix']
    h5name=h5prefix+"proj_maps_%03d.hdf5" % snapnum

    ## check if we've already projected this setup and saved it to intermediate file
    this_setup_in_projection_file = checkProjectionFile(
	os.path.join(projection_dir,h5name),**kwargs)
    
    if overwrite or not this_setup_in_projection_file:
        ## compute the projections
        projectColumnDensityAndQuantity(
	    snapdir,snapnum,
	    projection_dir,**kwargs)

    print('plotting image grid')
    print(list(kwargs.keys()),'passed to plot_image_grid')
    plot_image_grid(
        ax,
        projection_dir = projection_dir,
        snapnum = snapnum,
        **kwargs)

    return ax
