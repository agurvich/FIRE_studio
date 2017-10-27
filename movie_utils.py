from movie_maker_gasTwoColour import plot_image_grid
from movie_maker_gasDensity_v2 import compute_image_grid as compute_density_grid
from movie_maker_gasTemperature_v2 import compute_image_grid as compute_temp_grid


import matplotlib.pyplot as plt
import numpy as np 
import os

from distinct_colours import get_distinct

from sne_utils import findClusterExtent,drawSN
from all_utils import getTemperature

from readsnap import readsnap


def makeOutputDirectories(datadir):
    os.mkdir(os.path.join(datadir,'Plots'))
    os.mkdir(os.path.join(datadir,'Plots','GasTwoColour'))
    os.mkdir(os.path.join(datadir,'Plots','Projections'))
    os.mkdir(os.path.join(datadir,'Plots','Projections','Den'))
    os.mkdir(os.path.join(datadir,'Plots','Projections','Temp'))


def plotClusterAtSnap(ax,gal,snap,cluster,datadir,savename,snappath):
    short = 1
    to_plot=[]

    timepath=os.path.join(datadir,"snaptimes"+"_shortsnap.txt")
    snaps,times=np.genfromtxt(timepath,unpack=1)
    snaptime = times[snaps==snap]

    for sn in cluster:
        if sn.launch_time < snaptime:
            to_plot+=[sn]
    if len(to_plot):
        center,extent,(xs,ys,rcools)=findClusterExtent(cluster,real=snappath)

        #add background
        new_ax=addPrettyGalaxyToAx(ax,snap,gal,center,frame_width=extent,frame_depth=extent,
            snappath=snappath,datadir=datadir,savefig=0,savename=savename,overwrite=0
            ,min_den=-0.35,short=short)

        for j,sn in enumerate(to_plot): 
            drawSN(plt,ax,sn,snaptime,loc=(xs[j],ys[j]))
        new_ax.set_xlim([center[0]-extent,center[0]+extent])
        new_ax.set_ylim([center[1]-extent,center[1]+extent])
        new_ax.get_figure().savefig(os.path.join(datadir,'movies','frames',savename))
        plt.close()

    return ax

def projectDenTemp(snapdir,snapnum,savename,
    dprojects,tprojects,
    frame_center,frame_width,frame_depth,
    flag_cosmological,
    **kwargs):
    """ 
    Input: 
        snapdir - directory the snapshots live in
        snapnum - snapshot number currently rendering
        savename - name of .hdf5 file
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
        savename="%s_%d"%(savename,snapnum),output_dir = dprojects,isnap = snapnum,
        **kwargs)

    print 'projecting the temperature grid'
    temp_grid=compute_temp_grid(
        frame_center=frame_center,frame_width=frame_width,frame_depth=frame_depth,
        savename="%s_%d"%(savename,snapnum),output_dir = tprojects,isnap = snapnum,
        **kwargs)

def addPrettyGalaxyToAx(ax,snapdir,snapnum,savename,
    frame_center,frame_width,frame_depth,
    overwrite=0,savefig=0,datadir=None,**kwargs):
    """
    Input:
        ax - matplotlib axis object to draw to
        snapdir - directory where the snapshots live
        snapnum - snapshot number to draw
        savename - what to call the .hdf5 file
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

    print 'Drawing',savename,'to:',datadir
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
            frame_center,frame_width,savename=savename+'.hdf5',**kwargs)

    except IOError:
        print "Haven't computed image grids"
        projectDenTemp(snapdir,snapnum,savename+'.hdf5',dprojects,tprojects,
            frame_center,frame_width,frame_depth,flag_cosmological = 'm12i' in savename,**kwargs)
        print 'plotting image grid'
        plot_image_grid(ax,snapnum,dprojects,tprojects,
            frame_center,frame_width,savename=savename+'.hdf5',**kwargs)

    if savefig:
        ax.get_figure().set_size_inches(6,6)
        pixels = 1200 if 'pixels' not in kwargs else kwargs['pixels'] 
        image_name = "frame_%3d_%.2fkpc.png" % (snapnum, 2*frame_width)
        ax.get_figure().savefig(os.path.join(datadir,'Plots','GasTwoColour',image_name))
    
    return ax

def main(savename,snapnum=0,pixels=1200,frame_width=10,cosmo=0,**kwargs):
    cosmo == 'm12i' in savename
    if cosmo:
        raise Exception("Not Implemented!")
    else:
        frame_center=(0,0)

    plt.plot(frame_center[0],frame_center[1]) 
    ax=plt.gca()

    return addPrettyGalaxyToAx(ax,snapnum,savename,frame_center,frame_width,frame_depth,pixels)


if __name__=='__main__':
    argv = sys.argv[1:]
    opts,args = getopt.getopt(argv,'',['pixels=','gal=','snap=',
        'short=','frame_width='])

    #options:
    #--pixels : how many pixels in each direction to use, defaults to 1200
    #--gal : which simulation to use, defaults to m1e11f20
    #--snap : which snapshot
    #--short : flag for short snapshots starting at 250
    #--frame_width : size of frame in kpc

    for i,opt in enumerate(opts):
        if opt[1]=='':
            opts[i]=('mode',opt[0].replace('-'))
        else:
            opts[i]=(opt[0].replace('-',''),opt[1])
    main(**dict(opts))
