from movie_maker_gasTwoColour import plot_image_grid
from movie_maker_gasDensity_v2 import compute_image_grid as compute_density_grid
from movie_maker_gasTemperature_v2 import compute_image_grid as compute_temp_grid


import matplotlib.pyplot as plt
import numpy as np 
import os

from readsnap import readsnap


## system functions
def makeOutputDirectories(datadir):
    os.mkdir(os.path.join(datadir,'Plots'))
    os.mkdir(os.path.join(datadir,'Plots','GasTwoColour'))
    os.mkdir(os.path.join(datadir,'Plots','Projections'))
    os.mkdir(os.path.join(datadir,'Plots','Projections','Den'))
    os.mkdir(os.path.join(datadir,'Plots','Projections','Temp'))

## math functions 

## physics functions
def getTemperature(U_code,y_helium,ElectronAbundance):
    """U_codes = res['u']
        y_heliums = res['z'][:,1]
        ElectronAbundance=res['ne']"""
    U_cgs = U_code*1e10
    gamma=5/3.
    kB=1.38e-16 #erg /K
    m_proton=1.67e-24 # g
    mu = (1.0 + 4*y_helium) / (1+y_helium+ElectronAbundance) 
    mean_molecular_weight=mu*m_proton
    return mean_molecular_weight * (gamma-1) * U_cgs / kB

def projectDenTemp(snapdir,snapnum,
    dprojects,tprojects,
    frame_center,frame_width,frame_depth,
    flag_cosmological,
    **kwargs):
    """ 
    Input: 
        snapdir - directory the snapshots live in
        snapnum - snapshot number currently rendering
        dprojects - place to save density projections
        tprojects - place to save temperature projections
        frame_center - origin of image in data space
        frame_width - half-width of image in data space
        frame_depth - half-depth of image in data space
        available kwargs:
            theta=0- euler rotation angle
            phi=0- euler rotation angle
            psi=0 - euler rotation angle
            pixels=1200 - the resolution of image (pixels x pixels)
            min_den=-0.4 - the minimum of the density color scale
            max_den=1.6 - the maximum of the density color scale
            min_temp=2 - the minimum of the temperature color scale
            max_temp=7 - the maximum of the temperature color scale
    """

    keys = ['pos_all','mass_all','temperature_all',
        'HubbleParam','time_Myr','BoxSize']

    needToRead = False
    for key in keys:
        if key not in kwargs:
            needToRead = True

    if needToRead:
        print "Reading the snapshot values from",snapdir,snapnum
        res = readsnap(snapdir,snapnum,0,cosmological=flag_cosmological)
        kwargs['pos_all'] = res['p']
        kwargs['mass_all'] = res['m']

        kwargs['temperature_all']=getTemperature(res['u'],res['z'][:,1],res['ne'])

        kwargs['HubbleParam'] = res['hubble']
        kwargs['time_Myr'] = res['time'] # could be scale factor if flag_cosmological
        kwargs['BoxSize'] = res['boxsize']
        print 'Done reading the snapshot values'
    else:
        print 'You provided values as kwargs! nice.'

    print 'projecting the density grid'
    den_grid=compute_density_grid(
        frame_center=frame_center,frame_width=frame_width,frame_depth=frame_depth,
        output_dir = dprojects,isnap = snapnum,
        **kwargs)

    print 'projecting the temperature grid'
    temp_grid=compute_temp_grid(
        frame_center=frame_center,frame_width=frame_width,frame_depth=frame_depth,
        output_dir = tprojects,isnap = snapnum,
        **kwargs)

def addPrettyGalaxyToAx(ax,snapdir,snapnum,
    frame_center,frame_width,frame_depth,
    overwrite=0,savefig=0,datadir=None,**kwargs):
    """
    Input:
        ax - matplotlib axis object to draw to
        snapdir - directory where the snapshots live
        snapnum - snapshot number to draw
        frame_center - origin of image in data space 
        frame_width - half-width of image in data space
        frame_depth - half-depth of image in data space 
        overwrite - flag to overwrite the intermediate grid files
        savefig - flag to save a png of the image
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
    """

    if datadir ==None:
        datadir = snapdir

    print 'Drawing',snapdir,'to:',datadir
    if 'Plots' not in os.listdir(datadir):
        makeOutputDirectories(datadir)

    #where to find/save gas/temperature density grids-- this could get crowded!
    dprojects=os.path.join(datadir,'Plots','Projections','Den/')
    tprojects=os.path.join(datadir,'Plots','Projections','Temp/')

    try:
        if overwrite:
            ## pretend we couldn't find the precomputed image grids
            raise IOError
        print 'plotting image grid'
        plot_image_grid(ax,snapnum,dprojects,tprojects,
            frame_center,frame_width,**kwargs)

    except IOError:
        print "Haven't computed image grids"
        projectDenTemp(snapdir,snapnum,dprojects,tprojects,
            frame_center,frame_width,frame_depth,flag_cosmological = 'm12i' in snapdir,**kwargs)
        print 'plotting image grid'
        plot_image_grid(ax,snapnum,dprojects,tprojects,
            frame_center,frame_width,**kwargs)

    if savefig:
        ax.get_figure().set_size_inches(6,6)
        pixels = 1200 if 'pixels' not in kwargs else kwargs['pixels'] 
        image_name = "frame_%3d_%.2fkpc.png" % (snapnum, 2*frame_width)
        ax.get_figure().savefig(os.path.join(datadir,'Plots','GasTwoColour',image_name))
    
    return ax
