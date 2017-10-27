import numpy as np 
import os,sys,getopt
import copy
import matplotlib.pyplot as plt

from utils.cosmoExtractor import findGalaxyAndOrient,rotateVectorsZY
from utils.movie_utils import addPrettyGalaxyToAx,getTemperature
from utils.readsnap import readsnap

def loadDataFromSnapshot(snapdir,snapnum,mode,frame_width,frame_depth):
    if 'r' in mode:
        print "using readsnap to load in data"
        ## assumes the only sort of multi-part snapshot you would have is in cosmological units
        ## if you have a multipart snapshot that isn't in cosmological units pray that h=1, or change
        ## this flag
        res = readsnap(snapdir,snapnum,0,cosmological = 1)

        pos_all = res['p']
        mass_all = res['m']
        temperature_all = getTemperature(res['u'],res['z'][:,1],res['ne'])
        
        thetay,thetaz,frame_center,gindices = findGalaxyAndOrient(snapdir,snapnum,pos_all,res['rho'],
            frame_width,frame_depth)
        
        ## filter and free up memory
        pos = rotateVectorsZY(thetay,thetaz,pos_all[gindices]-frame_center)
        frame_center = np.zeros(3) # plot at center of mass
        del pos_all
        mass = mass_all[gindices]
        del mass_all
        temperature = temperature_all[gindices]
        del temperature_all

        mydict = {
            'pos_all':pos,'mass_all':mass,'temperature_all':temperature,
            'HubbleParam':res['hubble'],'time_Myr':res['time'],
            'BoxSize':res['boxsize'],'frame_center' : frame_center
        }
        
    elif 's' in mode:
        print 'using h5py to load in data'
        mydict={}
        raise Exception("Unimplemented!")
    return mydict


def renderGalaxy(ax,snapdir,snapnum,savefig=1,mode='r',**kwargs):
    # copy the dictionary so we don't mess anything up 
    copydict = copy.copy(kwargs)
    # add in filtered galaxy data
    copydict.update(loadDataFromSnapshot(snapdir,snapnum,mode,kwargs['frame_width'],kwargs['frame_depth']))

    # make projection map and add to canvas
    print copydict.keys(),'passed'
    addPrettyGalaxyToAx(
        ax,snapdir,snapnum,
        savefig=savefig,**copydict)

def main(snapdir,snapstart,snapmax,**kwargs):
    for snapnum in xrange(snapstart,snapmax):
        renderGalaxy(plt.gca(),snapdir,snapnum,**kwargs)
        plt.clf()
        

if __name__=='__main__':
    argv = sys.argv[1:]
    opts,args = getopt.getopt(argv,'rs',[
        'snapdir=',
        'snapstart=','snapmax=',
        'pixels=','frame_width=','frame_depth=',
        'theta=','phi=','psi=',
        'min_den=','max_den=','min_temp=','max_temp=','datadir='])

    #options:
    # -r/s = use readsnap or use single snapshot loader
    #--pixels : how many pixels in each direction to use, defaults to 1200
    #--snapstart : which snapshot to start the loop at
    #--snapmax : which snapshot to end the loop at
    #--frame_width : half width of frame in kpc
    #--frame_depth : half depth of frame in kpc
    #--theta,phi,psi : euler angles for rotation
    #--min/max_den/temp: bottom/top of color scales for density/temperature
    #--datadir: place to output frames to

    for i,opt in enumerate(opts):
        if opt[1]=='':
            opts[i]=('mode',opt[0].replace('-'))
        else:
            key = opt[0].replace('-','')
            if key in ['snapdir']:
                value= opt[1]
            else:
                # turn arguments from strings to whatever
                value = eval(opt[1])
            opts[i]=(key,value)
    main(**dict(opts))

