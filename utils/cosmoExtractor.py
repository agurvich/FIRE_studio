##### This script will output *only* the gas particles near a galaxy from a cosmological
##### simulation into a "subsnapshot" that has been centered/rotated to be edge on

import h5py,sys,getopt,os
import numpy as np
#from readsnap import readsnap

def makeOutputDir(savename):
    projectDir="/home/abg6257/projects/subsnaps"
    datadir=os.path.join(projectDir,savename)
    if savename not in os.listdir(os.path.join(projectDir)):
        print 'making directory subsnaps/%s in'%savename
        os.mkdir(datadir)
    return datadir

###### Math Functions
def rotationMatrixY(theta):
    return np.array([
            [np.cos(theta),0,-np.sin(theta)],
            [0,1,0],
            [np.sin(theta),0,np.cos(theta)]
        ])

def rotationMatrixX(theta):
    return np.array([
            [1,0,0],
            [0,np.cos(theta),np.sin(theta)],
            [0,-np.sin(theta),np.cos(theta)]
        ])

def rotationMatrixZ(theta):
    return np.array([
            [np.cos(theta),np.sin(theta),0],
            [-np.sin(theta),np.cos(theta),0],
            [0,0,1]
        ])

def rotateVectorsZY(thetay,thetaz,vectors):
    rotatedCoords=rotateVectors(rotationMatrixZ(thetaz),vectors)
    rotatedCoords=rotateVectors(rotationMatrixY(thetay),rotatedCoords)
    return rotatedCoords

def rotateVectors(rotationMatrix,vectors):
    return np.dot(rotationMatrix,vectors.T).T

def getVcom(masses,velocities):
    return np.sum(masses[:,None]*velocities,axis=0)/np.sum(masses)

def getAngularMomentum(vectors,masses,velocities):
    return np.sum(np.cross(vectors,masses[:,None]*velocities),axis=0)

def extractSphericalVolumeIndices(rs,rcom,radius2):
    return np.sum((rs - rcom)**2.,axis=1) < radius2

def extractCylindricalVolumeIndices(rs,rcom,radius,height):
    xyindices = np.sum((rs-rcom)[:,:2]**2.,axis=1)<radius**2.
    #absolute value is expensive?
    zindices = (rs[:,2]-rcom[2])**2. < height**2.
    indices = np.logical_and(xyindices,zindices)
    return indices

def iterativeCoM(coords,masses,n=4,r0=np.array([0,0,0])):
    rcom = r0
    for i in xrange(n):
        sindices= extractSphericalVolumeIndices(coords,rcom,1000**2/10**i)
        #sindices= extractSphericalVolumeIndices(coords,rcom,10000**2/100**i)
        rcom = np.sum(coords[sindices]*masses[sindices][:,None],axis=0)/np.sum(masses[sindices])
    return rcom

def getThetas(angMom):
    thetay = np.arctan2(np.sqrt(angMom[0]**2+angMom[1]**2),(angMom[2]))
    thetaz = np.arctan2(angMom[1],angMom[0])
    print "want to rotate by",thetay*180/np.pi,thetaz*180/np.pi
    return thetay,thetaz

def extractDisk(snap,radius=30,cylinder=15):
    """Takes an open snapshot and returns the orientation of the disk,
        its location, and the index arrays for the stars/gas particles 
        that make it up
        Input:
            snap - an h5py.File instance
            radius - in kpc, the radius to which you want to include particles
        Output:
            thetay/thetaz - orientation angles
            scom/vscom - location of the center of mass of the star particles and its velocity
            gindinces/sindices - index arrays for gas/star particles within 30kpc"""
    HubbleParam = snap['Header'].attrs['HubbleParam']

    ## load star particle info
    srs = np.array(snap['PartType4/Coordinates'],dtype=np.float64)/HubbleParam #kpc
    svs = np.array(snap['PartType4/Velocities'],dtype=np.float64) #km/s
    smasses = np.array(snap['PartType4/Masses'],dtype=np.float64)/HubbleParam*1e10 #solar masses

    ## load gas particle info
    rs = np.array(snap['PartType0/Coordinates'],dtype=np.float64)/HubbleParam #kpc
    rhos = np.array(snap['PartType0/Density'],dtype=np.float64)*HubbleParam**2.# mass / len^3 (code)
    return extractDiskFromArrays(srs,svs,smasses,rs,rhos,radius,cylinder=cylinder)

def extractDiskFromArrays(srs,svs,smasses,rs,rhos,radius,cylinder=0):
    """Takes arrays from a snapshot and returns the information required
        from extractDisk. Useful to separate so that external protocols can 
        call it (multiple snapshots, for instance)
        Input: 
            srs/svs/smasses - positions,velocities, and masses of star particles
            rs/rhos - positions and densities of gas particles 
            radius - radius to extract particles from
            cylinder - cylindrical height, if 0 extracts a sphere instead"""
    ## find com using iterative shells
    scom = np.sum(smasses[:,None]*srs,axis=0)/np.sum(smasses)
    scom = iterativeCoM(srs,smasses,r0=scom)#rs[np.argmax(rhos)])

    ## extract particles within radius
    sindices= extractSphericalVolumeIndices(srs,scom,radius**2)
    gindices= extractSphericalVolumeIndices(rs,scom,radius**2)

    ## find com velocity of disk
    vscom = getVcom(smasses[sindices],svs[sindices])

    ## find angular momentum vector
    angMom = getAngularMomentum((srs-scom)[sindices],
        smasses[sindices],(svs-vscom)[sindices])
    ## find angles necessary to rotate coordinates to align with angMom
    thetay,thetaz = getThetas(angMom)

    if cylinder:
        new_com = rotateVectorsZY(thetay,thetaz,scom)
        new_rs=rotateVectorsZY(thetay,thetaz,rs-scom)
        new_srs=rotateVectorsZY(thetay,thetaz,srs-scom)
        sindices= extractCylindricalVolumeIndices(new_srs,np.array([0,0,0]),radius,cylinder)
        gindices= extractCylindricalVolumeIndices(new_rs,np.array([0,0,0]),radius,cylinder)

    return thetay,thetaz,scom,vscom,gindices,sindices

def getHaloIndices(snap):
    thetay,thetaz,scom,vscom,gindices,sindices=extractDisk(snap)
    rs = np.array(snap['PartType0/Coordinates'])
    rotatedCoords=rotateVectorsZY(thetay,thetaz,rs[gindices]-scom)
    return gindices,rotatedCoords

def outputHDF5(datadir,snapdir,snapnum):
    with h5py.File(snapdir+"/snapshot_%03d.hdf5"%snapnum,'r') as snap:
        gindices,rotatedCoords=getHaloIndices(snap)
        writeSnap(snap,gindices,rotatedCoords,datadir,snapnum)

def writeSnap(snap,gindices,rotatedCoords,datadir,snapnum):
    HubbleParam=snap['Header'].attrs['HubbleParam']
    with h5py.File(datadir+"/snapshot_%03d.hdf5"%snapnum,'w') as buildsnap:
        ## write header info
        buildsnap.create_group("Header")
        buildsnap['Header'].attrs['Time']=snap['Header'].attrs['Time']
        buildsnap['Header'].attrs['HubbleParam']=1
        buildsnap['Header'].attrs['BoxSize']=snap['Header'].attrs['BoxSize']

        ## write gas particle info in non-cosmo units
        buildsnap['PartType0/Coordinates']=rotatedCoords
        buildsnap['PartType0/Masses']=np.array(snap['PartType0/Masses'],dtype=np.float64)[gindices]/HubbleParam

        ## write "cosmo-unit-less" gas particle info
        buildsnap['PartType0/ElectronAbundance']=np.array(snap['PartType0/ElectronAbundance'],dtype=np.float64)[gindices]
        buildsnap['PartType0/Metallicity']=np.array(snap['PartType0/Metallicity'],dtype=np.float64)[gindices]
        buildsnap['PartType0/InternalEnergy']=np.array(snap['PartType0/InternalEnergy'],dtype=np.float64)[gindices]

def main(low,high,savename,**kwargs):
    low,high=int(low),int(high)
    snapdir = "/home/abg6257/projects/isoDisk/%s/output"%savename
    datadir = makeOutputDir(savename)
    for snapnum in xrange(low,high):
        outputHDF5(datadir,snapdir,snapnum)
    print "all done!"
    
if __name__=='__main__':
    argv = sys.argv[1:]
    opts,args = getopt.getopt(argv,'',['low=','high=',"savename="])
    #options:
    #-r : read in previous starobjects.npz, default behavior
    #-n : ignore any preexisting starobjects.npz 
    #--snap : snapshot number 
    #--gal  : galaxy prefix 
    for i,opt in enumerate(opts):
        if opt[1]=='':
            opts[i]=('mode',opt[0].replace('-',''))
        else:
            opts[i]=(opt[0].replace('-',''),opt[1])
    main(**dict(opts))    
