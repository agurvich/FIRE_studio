from firestudio.utils.gas_utils.movie_maker_gasTwoColour import plot_image_grid
from firestudio.utils.gas_utils.movie_maker_gasDensity_v2 import compute_image_grid as compute_density_grid
from firestudio.utils.gas_utils.movie_maker_gasTemperature_v2 import compute_image_grid as compute_temp_grid

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
    **kwargs):
    try:
	with h5py.File(projection_file,'r') as handle:
	    for group in handle.keys():
		this_group = handle[group]
		flag = 1
		for key,variable in zip(
		    ['npix_x','frame_half_width','frame_depth','frame_center','theta','phi','psi'],
		    [ pixels , frame_half_width , frame_depth , frame_center , theta , phi , psi ]):
    
		    ## read the value in the hdf5 file and compare to variable
		    if key not in ['npix_x']:
			## key is not an integer so we have to round it somehow
			flag = flag and np.all(this_group[key].value == np.round(variable,decimals=2))
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

    raise Exception("NO!")
    print('projecting the density grid')
    den_grid=compute_density_grid(
        projection_dir = projection_dir,snapnum = snapnum,
        **kwargs)

def addPrettyGalaxyToAx(ax,snapdir,snapnum,
    overwrite=0,datadir=None,**kwargs):
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
	
	BoxSize -

	frame_center - origin of image in data space 
	frame_half_width - half-width of image in data space
	frame_depth - half-depth of image in data space 

    Optional kwargs to be passed along: 
	quantity_name='temperature' - the name of the quantity that you're mass weighting
	    should match whatever array you're passing in as quantity

	theta=0- euler rotation angle
	phi=0- euler rotation angle
	psi=0 - euler rotation angle
	pixels=1200 - the resolution of image (pixels x pixels)
	min_den=-0.4 - the minimum of the density color scale
	max_den=1.6 - the maximum of the density color scale
	min_temp=2 - the minimum of the temperature color scale
	max_temp=7 - the maximum of the temperature color scale

	h5prefix='' - a string that you can prepend to the projection filename if desired
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
    print kwargs.keys(),'passed to plot_image_grid'
    plot_image_grid(ax,snapnum,projection_dir,
        **kwargs)

    return ax
