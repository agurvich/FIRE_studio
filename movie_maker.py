import numpy as np 
import os,sys,getopt
import copy
import matplotlib.pyplot as plt

from cosmoExtractor import extractDiskFromArrays,rotateVectorsZY
from movie_utils import addPrettyGalaxyToAx
from readsnap import readsnap


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

def findGalaxyAndOrient(snapdir,snapnum,gaspos,gasdens,frame_width,frame_depth):
    ## load in stars to find their center of mass
    star_res = readsnap(snapdir,snapnum,4,cosmological=1)

    args = {
        'srs':star_res['p']
        ,'svs':star_res['v']
        ,'smasses':star_res['m']
        ,'rs':gaspos
        ,'rhos':gasdens
        ,'radius':2**0.5*frame_width#kpc
        ,'cylinder':frame_depth#kpc
        }

    thetay,thetaz,scom,vscom,gindices,sindices = extractDiskFromArrays(**args)
    del star_res

    return thetay,thetaz,scom,gindices

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
    savename = snapdir[-6:]+".hdf5" #ignored for now...
    # copy the dictionary so we don't mess anything up 
    copydict = copy.copy(kwargs)
    # add in filtered galaxy data
    copydict.update(loadDataFromSnapshot(snapdir,snapnum,mode,kwargs['frame_width'],kwargs['frame_depth']))

    # make projection map and add to canvas
    print copydict.keys(),'passed'
    addPrettyGalaxyToAx(
        ax,snapdir,snapnum,savename=savename,
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
        'min_den=','max_den=','min_temp=','max_temp='])

    #options:
    # -r/s = use readsnap or use single snapshot loader
    #--pixels : how many pixels in each direction to use, defaults to 1200
    #--snapstart : which snapshot to start the loop at
    #--snapmax : which snapshot to end the loop at
    #--frame_width : half width of frame in kpc
    #--frame_depth : half depth of frame in kpc
    #--theta,phi,psi : euler angles for rotation
    #--min/max_den/temp: bottom/top of color scales for density/temperature
    

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

