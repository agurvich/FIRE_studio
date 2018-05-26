from firestudio.utils.gas_utils.movie_maker_gasTwoColour import plot_image_grid
from firestudio.utils.gas_utils.movie_maker_gasDensity_v2 import compute_image_grid as compute_density_grid
from firestudio.utils.gas_utils.movie_maker_gasTemperature_v2 import compute_image_grid as compute_temp_grid

from abg_python.snapshot_utils import openSnapshot

import matplotlib 
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
import numpy as np 
import os

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
    if not os.path.isdir(os.path.join(path,'Den')):
        os.mkdir(os.path.join(path,'Den'))
    if not os.path.isdir(os.path.join(path,'Temp')):
        os.mkdir(os.path.join(path,'Temp'))

## physics functions
def projectDenTemp(snapdir,snapnum,
    dprojects,tprojects,
    flag_cosmological,
    **kwargs):
    """ 
    Input: 
        snapdir - directory the snapshots live in
        snapnum - snapshot number currently rendering
        dprojects - place to save density projections
        tprojects - place to save temperature projections
        available kwargs:
            theta=0- euler rotation angle
            phi=0- euler rotation angle
            psi=0 - euler rotation angle
            pixels=1200 - the resolution of image (pixels x pixels)
            min_den=-0.4 - the minimum of the density color scale
            max_den=1.6 - the maximum of the density color scale
            min_temp=2 - the minimum of the temperature color scale
            max_temp=7 - the maximum of the temperature color scale
            edgeon - create and plot an edgeon view

            frame_center - origin of image in data space
            frame_width - half-width of image in data space
            frame_depth - half-depth of image in data space
    """

    print 'projecting the density grid'
    den_grid=compute_density_grid(
        output_dir = dprojects,isnap = snapnum,
        **kwargs)

    print 'projecting the temperature grid'
    temp_grid=compute_temp_grid(
        output_dir = tprojects,isnap = snapnum,
        **kwargs)

def addPrettyGalaxyToAx(ax,snapdir,snapnum,
    overwrite=0,datadir=None,**kwargs):
    """
    Input:
        ax - matplotlib axis object to draw to
        snapdir - location that the snapshots live in
        snapnum - snapshot number
    Optional:
        overwrite - flag to overwrite the intermediate grid files
        datadir - directory to output the the intermediate grid files and output png

        kwargs able to be passed along: 
            theta=0- euler rotation angle
            phi=0- euler rotation angle
            psi=0 - euler rotation angle
            pixels=1200 - the resolution of image (pixels x pixels)
            min_den=-0.4 - the minimum of the density color scale
            max_den=1.6 - the maximum of the density color scale
            min_temp=2 - the minimum of the temperature color scale
            max_temp=7 - the maximum of the temperature color scale
            edgeon - create and plot an edgeon view

            frame_center - origin of image in data space 
            frame_width - half-width of image in data space
            frame_depth - half-depth of image in data space 

    """

    print 'Drawing',snapdir,'to:',datadir
    makeOutputDirectories(datadir)

    #where to find/save gas/temperature density grids-- this could get crowded!
    dprojects=os.path.join(datadir,'Plots','Projections','Den/')
    tprojects=os.path.join(datadir,'Plots','Projections','Temp/')

    if overwrite:
        ## compute the projections
        projectDenTemp(snapdir,snapnum,dprojects,tprojects,
            flag_cosmological = 'm12i' in snapdir,**kwargs)

    print 'plotting image grid'
    plot_image_grid(ax,snapnum,dprojects,tprojects,
        **kwargs)

    return ax
