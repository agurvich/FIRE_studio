#########################################################
# Ported IDL script for use on the cluster              #
# For making 2d images of GADGET sims                   #
# Created 02/04/2010 by RAC                             #
# Last modified: N/A                                    #
#########################################################

#!/usr/bin/env python

# Import system and vtk routines
import os
import sys
import string
import numpy as np 
#import tables
import h5py
import ctypes

def compute_image_grid(
    pos_all,mass_all,temperature_all,
    BoxSize,
    frame_center,frame_half_width,frame_depth,
    output_dir,snapnum,
    theta=0,phi=0,psi=0,
    pixels=1200,
    min_den=-1.0,max_den=1.2,
    edgeon=0,
    h5filename='',**kwargs):
    if edgeon==1:
	raise Exception("Unimplemented edgeon!")

    print("extra kwargs in compute_den:",kwargs.keys())
    print(' rotation = (',theta,',',phi,',',psi,')')

    ## Set image size 
    npix_x   = pixels #1200 by default
    npix_y   = pixels #1200 by default

    ## calculate edges of the frame
    Xmin,Ymin = -frame_half_width + frame_center[:2]
    Xmax,Ymax = frame_half_width + frame_center[:2]
    Zmin,Zmax = -frame_depth+frame_center[2],frame_depth+frame_center[2]

    ## extract a cube of particles that are in relevant area
    print('extracting cube')
    ind_box = ((pos_all[:,0] > Xmin) & (pos_all[:,0] < Xmax) &
               (pos_all[:,1] > Ymin) & (pos_all[:,1] < Ymax) &
               (pos_all[:,2] > Zmin) & (pos_all[:,2] < Zmax))

    pos = pos_all[ind_box].astype(np.float32)
    mass = mass_all[ind_box].astype(np.float32)
    quantity = temperature_all[ind_box].astype(np.float32)

    print('-done')

    ## rotate by euler angles if necessary
    pos = rotateEuler(theta,phi,psi,pos,frame_center)

    ## make the actual C call
    columnDensityMap, massWeightedQuantityMap = getImageGrid(
	BoxSize,
	Xmin,Xmax,
	Ymin,Ymax,
	Zmin,Zmax,
	npix_x,npix_y,
	pos,mass,quantity)

    ## write the output to an .hdf5 file
    writeImageGrid(
        snapnum,
        output_dir,
	columnDensityMap, massWeightedQuantityMap,
        npix_x,frame_half_width,frame_depth,
	frame_center,
	theta,phi,psi,
	h5filename=h5filename)

def rotateEuler(
    theta,phi,psi,
    pos,frame_center):

    ## if need to rotate at all really -__-
    if theta==0 and phi==0 and psi==0:
        return pos
    # rotate particles by angle derived from frame number
    theta_rad = np.pi*theta/ 1.8e2
    phi_rad   = np.pi*phi  / 1.8e2
    psi_rad   = np.pi*psi  / 1.8e2

    # construct rotation matrix
    rot_matrix      = new_matrix(3,3)
    print('theta = ',theta_rad)
    print('phi   = ',phi_rad)
    print('psi   = ',psi_rad)

    ## explicitly define the euler rotation matrix 
    rot_matrix = np.array([
	[np.cos(phi_rad)*np.cos(psi_rad), #xx
	    -np.cos(phi_rad)*np.sin(psi_rad), #xy
	    np.sin(phi_rad)], #xz
	[np.cos(theta_rad)*np.sin(psi_rad) + np.sin(theta_rad)*np.sin(phi_rad)*np.cos(psi_rad),#yx
	    np.cos(theta_rad)*np.cos(psi_rad) - np.sin(theta_rad)*np.sin(phi_rad)*np.sin(psi_rad),#yy
	    -np.sin(theta_rad)*np.cos(phi_rad)],#yz
	[np.sin(theta_rad)*np.sin(psi_rad) - np.cos(theta_rad)*np.sin(phi_rad)*np.cos(psi_rad),#zx
	    np.sin(theta_rad)*np.cos(psi_rad) - np.cos(theta_rad)*np.sin(phi_rad)*np.sin(psi_rad),#zy
	    np.cos(theta_rad)*np.cos(phi_rad)]#zz
	])
	    
    ## translate the particles so we are rotating about the center of our frame
    pos-=frame_center

    ## rotate each of the vectors in the pos array
    pos = np.dot(rot_matrix,pos.T).T

    # translate particles back
    pos+=frame_center
    
    return pos_rot
    
def getImageGrid(
    BoxSize,
    Xmin,Xmax,
    Ymin,Ymax,
    Zmin,Zmax,
    npix_x,npix_y,
    pos,mass,quantity):

    ## set c-routine variables
    desngb   = 32
    Axis1    = 0
    Axis2    = 1
    Axis3    = 2

    Hmax     = 0.5*(Xmax-Xmin)

    n_smooth = pos.shape[0]

    ## output array for sum along the line of sight
    totalMassMap = np.zeros(shape = (npix_x,npix_y),dtype=np.float32)

    ## output array for average along the line of sight
    massWeightedQuantityMap = np.zeros(shape = (npix_x,npix_y),dtype=np.float32)
    
    ## create hsml output array
    hsml = np.zeros(mass.shape[0])
    
    c_f_p      = ctypes.POINTER(ctypes.c_float)
    pos_p      = pos.ctypes.data_as(c_f_p)
    hsml_p     = hsml.ctypes.data_as(c_f_p)
    mass_p     = mass.ctypes.data_as(c_f_p)
    quantity_p = quantity.ctypes.data_as(c_f_p)
    w_f_p    = totalMassMap.ctypes.data_as(c_f_p)
    q_f_p    = massWeightedQuantityMap.ctypes.data_as(c_f_p)

    print('------------------------------------------')
    curpath = os.path.realpath(__file__)
    curpath = curpath[:len("utils")+curpath.index("utils")] #split off this filename
    c_obj = ctypes.CDLL(os.path.join(curpath,'gas_utils','HsmlAndProject_cubicSpline/HsmlAndProject.so'))
    
   #print(n_smooth)
   #print(pos_p)
   #print(hsml_p)
   #print(mass_p)
   #print(quantity_p)
   #print(Xmin,Xmax)
   #print(Ymin,Ymax)
   #print(Zmin,Zmax)
   #print(npix_x,npix_y)
   #print(desngb)
   #print(Axis1,Axis2,Axis3)
   #print(Hmax,BoxSize)

    c_obj.findHsmlAndProject(
	ctypes.c_int(n_smooth), ## number of particles
	pos_p,hsml_p,mass_p,quantity_p, ## position, mass, and "quantity" of particles
        ctypes.c_float(Xmin),ctypes.c_float(Xmax), ## xmin/xmax
	ctypes.c_float(Ymin),ctypes.c_float(Ymax), ## ymin/ymax
	ctypes.c_float(Zmin),ctypes.c_float(Zmax), ## zmin/zmax
        ctypes.c_int(npix_x),ctypes.c_int(npix_y), ## npixels
	ctypes.c_int(desngb), ## neighbor depth
        ctypes.c_int(Axis1),ctypes.c_int(Axis2),ctypes.c_int(Axis3), ## axes...?
	ctypes.c_float(Hmax),ctypes.c_double(BoxSize), ## maximum smoothing length and size of box
	w_f_p,q_f_p) ## pointers to output cell-mass and cell-mass-weighted-quantity
    print('------------------------------------------')

    # normalise by area to get SFC density (column density)
    ## NOTE does the image HAVE to be square? should look at the C routine
    Acell = (Xmax-Xmin)/npix_x * (Ymax-Ymin)/npix_y
    columnDensityMap = totalMassMap/(Acell) # 10^10 Msun / kpc^-2 
    
    # convert into Msun/pc^2
    unitmass_in_g = 1.9890000e+43 
    solar_mass    = 1.9890000e+33
    conv_fac = (unitmass_in_g/solar_mass) / (1.0e3)**2 ## Msun/pc^2
    columnDensityMap *= conv_fac
    columnDensityMap = np.log10(columnDensityMap)
    print(
	'log10 minmax(columnDensityMap)',
	np.min(columnDensityMap),
	np.max(columnDensityMap))

    # massWeightedQuantityMap contains the mass-weighted temperature, in K, most often
    massWeightedQuantityMap = np.log10(massWeightedQuantityMap)
    print(
	'log10 minmax(massWeightedQuantityMap)',
	np.min(massWeightedQuantityMap),
	np.min(massWeightedQuantityMap))
   
    return columnDensityMap,massWeightedQuantityMap

def writeImageGrid(
    snapnum,
    output_dir,
    columnDensityMap, massWeightedQuantityMap,
    npix_x,frame_half_width,frame_depth,
    frame_center,
    theta,phi,psi,
    h5filename='',
    this_setup_id=None):

    ## Write the image grids to HDF5 file 
    h5filename += "proj_maps_%d.hdf5" % snapnum
    h5name=os.path.join(output_dir,h5filename)

    ## what should we call this setup? need a unique identifier
    ## let the user give it a
    ##	custom name through kwargs later on TODO
    this_setup_id = (
	"npix%d_width%.2fkpc_depth%.2fkpc_x%.2f_y%.2f_z%.2f_theta%.2f_phi%.2f_psi%.2f"%(
	    npix_x, 2*frame_half_width,frame_depth,
	    frame_center[0],frame_center[1],frame_center[2],
	    theta,phi,psi)
	if this_setup_id is None else this_setup_id)

    with h5py.File(h5name, "w") as h5file:
	this_group = h5file.create_group(this_setup_id)
	
	## save the maps themselves
        this_group['columnDensityMap']=columnDensityMap
        this_group['massWeightedQuantityMap']=massWeightedQuantityMap
	
	## save the meta data
	this_group['npix_x']=npix_x
	this_group['frame_center']=frame_center
	this_group['frame_half_width']=frame_half_width
	this_group['frame_depth']=frame_depth
	this_group['theta']=theta
	this_group['phi']=phi
	this_group['psi']=psi

	## TODO should I put in metadata that allows you to recreate
	##  frames without access to the relevant snapshots?
	##  e.g. current time/redshift 
