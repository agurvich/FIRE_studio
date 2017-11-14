import h5py,sys,getopt,os
import numpy as np
from readsnap import readsnap

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

## geometry functions
def extractSphericalVolumeIndices(rs,rcom,radius2):
    return np.sum((rs - rcom)**2.,axis=1) < radius2

def extractCylindricalVolumeIndices(rs,rcom,radius,height):
    xyindices = np.sum((rs-rcom)[:,:2]**2.,axis=1)<radius**2.
    #absolute value is expensive?
    zindices = (rs[:,2]-rcom[2])**2. < height**2.
    indices = np.logical_and(xyindices,zindices)
    return indices

## physics functions
def getVcom(masses,velocities):
    return np.sum(masses[:,None]*velocities,axis=0)/np.sum(masses)

def getAngularMomentum(vectors,masses,velocities):
    return np.sum(np.cross(vectors,masses[:,None]*velocities),axis=0)

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

## main function
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

