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
import copy

def compute_image_grid(
    Coordinates,Masses,Quantity,
    BoxSize,
    frame_half_width,frame_depth,
    frame_center,
    projection_dir,snapnum,
    quantity_name='Temperature',
    this_setup_id = None,
    theta=0,phi=0,psi=0,
    pixels=1200,
    edgeon=0,
    h5prefix='',
    take_log_of_quantity=True,
    Hsml = None,
    **kwargs):
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
    ind_box = ((Coordinates[:,0] > Xmin) & (Coordinates[:,0] < Xmax) &
               (Coordinates[:,1] > Ymin) & (Coordinates[:,1] < Ymax) &
               (Coordinates[:,2] > Zmin) & (Coordinates[:,2] < Zmax))

    pos = Coordinates[ind_box].astype(np.float32)
    mass = Masses[ind_box].astype(np.float32)
    quantity = Quantity[ind_box].astype(np.float32)
    frame_center = frame_center.astype(np.float32)
    if Hsml is None:
        hsml = Hsml
    else:
        hsml = Hsml[ind_box].astype(np.float32)

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
	pos,mass,quantity,
        take_log_of_quantity,
        hsml = hsml)

    ## write the output to an .hdf5 file
    writeImageGrid(
        snapnum,
        projection_dir,
	columnDensityMap,
	massWeightedQuantityMap, quantity_name,
        npix_x,frame_half_width,frame_depth,
	frame_center,
	theta,phi,psi,
	h5prefix=h5prefix)

def rotateEuler(
    theta,phi,psi,
    pos,frame_center):

    ## if need to rotate at all really -__-
    if theta==0 and phi==0 and psi==0:
        return pos
    theta=phi=psi=1e-8
    # rotate particles by angle derived from frame number
    theta_rad = np.pi*theta/ 1.8e2
    phi_rad   = np.pi*phi  / 1.8e2
    psi_rad   = np.pi*psi  / 1.8e2

    # construct rotation matrix
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
	]).astype(np.float32)
	    

    ## rotate each of the vectors in the pos array about the frame center
    pos_rot = np.dot(rot_matrix,(pos-frame_center).T).T

    # translate particles back
    pos_rot+=frame_center
    

    assert np.allclose(pos,pos_rot)
    print "the arrays are identical"
    return pos
    return pos_rot.astype(np.float32)

def rotateEuler(theta,phi,psi,pos,frame_center):
    def new_matrix(m,n):
        # Create zero matrix
        matrix = [[0 for row in range(n)] for col in range(m)]
        return matrix

    pos-=frame_center
    ## if need to rotate at all really -__-
    if theta==0 and phi==0 and psi==0:
        return pos
    # rotate particles by angle derived from frame number
    pi        = 3.14159265e0
    theta_rad = pi*theta/ 1.8e2
    phi_rad   = pi*phi  / 1.8e2
    psi_rad   = pi*psi  / 1.8e2

    # construct rotation matrix
    rot_matrix      = new_matrix(3,3)
    #print 'theta = ',theta_rad
    #print 'phi   = ',phi_rad
    #print 'psi   = ',psi_rad
    rot_matrix[0][0] =  np.cos(phi_rad)*np.cos(psi_rad)
    rot_matrix[0][1] = -np.cos(phi_rad)*np.sin(psi_rad)
    rot_matrix[0][2] =  np.sin(phi_rad)
    rot_matrix[1][0] =  np.cos(theta_rad)*np.sin(psi_rad) + np.sin(theta_rad)*np.sin(phi_rad)*np.cos(psi_rad)
    rot_matrix[1][1] =  np.cos(theta_rad)*np.cos(psi_rad) - np.sin(theta_rad)*np.sin(phi_rad)*np.sin(psi_rad)
    rot_matrix[1][2] = -np.sin(theta_rad)*np.cos(phi_rad)
    rot_matrix[2][0] =  np.sin(theta_rad)*np.sin(psi_rad) - np.cos(theta_rad)*np.sin(phi_rad)*np.cos(psi_rad)
    rot_matrix[2][1] =  np.sin(theta_rad)*np.cos(psi_rad) - np.cos(theta_rad)*np.sin(phi_rad)*np.sin(psi_rad)
    rot_matrix[2][2] =  np.cos(theta_rad)*np.cos(phi_rad)
    #print rot_matrix

    ## this needs to be done before passing to rotateeuler
    # translate particles so centre = origin
    #pos[:,0] -= (Xmin + (L_x/2.))
    #pos[:,1] -= (Ymin + (L_y/2.))
    #pos[:,2] -= (Zmin + (L_z/2.))

    ## global variable inside real render 
    n_box = pos.shape[0]

    # rotate about each axis with a matrix operation
    pos_rot = np.ndarray((n_box, 3), dtype=np.float32)
    for ipart in range(n_box):
        pos_matrix = new_matrix(3,1)
        pos_matrix[0][0] = pos[ipart,0]
        pos_matrix[1][0] = pos[ipart,1]
        pos_matrix[2][0] = pos[ipart,2]
        rotated = np.matmul(rot_matrix,pos_matrix)
        pos_rot[ipart,0] = rotated[0][0]
        pos_rot[ipart,1] = rotated[1][0]
        pos_rot[ipart,2] = rotated[2][0]
    
    ## occurs outside rotateeuler
    # translate particles back
    #pos[:,0] += (Xmin + (L_x/2.))
    #pos[:,1] += (Ymin + (L_y/2.))
    #pos[:,2] += (Zmin + (L_z/2.))
    
    pos_rot+=frame_center
    return pos_rot
    
def getImageGrid(
    BoxSize,
    Xmin,Xmax,
    Ymin,Ymax,
    Zmin,Zmax,
    npix_x,npix_y,
    pos,mass,quantity,
    take_log_of_quantity,
    hsml=None):

    ## set c-routine variables
    desngb   = 32
    Axis1    = 0
    Axis2    = 1
    Axis3    = 2

    Hmax     = 0.5*(Xmax-Xmin)

    n_smooth = pos.shape[0]
    print(type(n_smooth),'nsmooth')

    ## output array for sum along the line of sight
    totalMassMap = np.zeros(shape = (npix_x,npix_y),dtype=np.float32)

    ## output array for average along the line of sight
    massWeightedQuantityMap = np.zeros(shape = (npix_x,npix_y),dtype=np.float32)
    
    ## create hsml output array
    if hsml is None:
        hsml = np.zeros(mass.shape[0],dtype=np.float32)
    else:
        print("Using provided smoothing lengths")
    
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
    print(pos,type(pos),type(pos[0]),type(pos[0][0]),np.shape(pos))

    c_obj.findHsmlAndProject(
	ctypes.c_int(n_smooth), ## number of particles
	pos_p,hsml_p,mass_p,quantity_p, ## position, mass, and "quantity" of particles
        ctypes.c_float(Xmin.astype(np.float32)),ctypes.c_float(Xmax.astype(np.float32)), ## xmin/xmax
	ctypes.c_float(Ymin.astype(np.float32)),ctypes.c_float(Ymax.astype(np.float32)), ## ymin/ymax
	ctypes.c_float(Zmin.astype(np.float32)),ctypes.c_float(Zmax.astype(np.float32)), ## zmin/zmax
        ctypes.c_int(npix_x),ctypes.c_int(npix_y), ## npixels
	ctypes.c_int(desngb), ## neighbor depth
        ctypes.c_int(Axis1),ctypes.c_int(Axis2),ctypes.c_int(Axis3), ## axes...?
	ctypes.c_float(Hmax),ctypes.c_double(BoxSize), ## maximum smoothing length and size of box
	w_f_p,q_f_p) ## pointers to output cell-mass and cell-mass-weighted-quantity
    print('------------------------------------------')
    print(pos,type(pos),type(pos[0]),type(pos[0][0]),np.shape(pos))

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

    # massWeightedQuantityMap contains the mass-weighted quantity
    if take_log_of_quantity:
        massWeightedQuantityMap = np.log10(massWeightedQuantityMap)
    print(
	'log10 minmax(massWeightedQuantityMap)',
	np.min(massWeightedQuantityMap),
	np.min(massWeightedQuantityMap))
   
    return columnDensityMap,massWeightedQuantityMap

def writeImageGrid(
    snapnum,
    projection_dir,
    columnDensityMap,
    massWeightedQuantityMap,quantity_name,
    npix_x,frame_half_width,frame_depth,
    frame_center,
    theta,phi,psi,
    h5prefix='',
    this_setup_id=None):

    ## Write the image grids to HDF5 file 
    h5prefix += "proj_maps_%d.hdf5" % snapnum
    h5name=os.path.join(projection_dir,h5prefix)

    ## what should we call this setup? need a unique identifier
    ## let the user give it a
    ##	custom name through kwargs later on TODO
    this_setup_id = (
	"npix%d_width%.2fkpc_depth%.2fkpc_x%.2f_y%.2f_z%.2f_theta%.2f_phi%.2f_psi%.2f"%(
	    npix_x, 2*frame_half_width,frame_depth,
	    frame_center[0],frame_center[1],frame_center[2],
	    theta,phi,psi)
	if this_setup_id is None else this_setup_id)

    with h5py.File(h5name, "a") as h5file:
        if this_setup_id not in list(h5file.keys()):
            this_group = h5file.create_group(this_setup_id)
            ## save the maps themselves
            this_group['columnDensityMap']=columnDensityMap
            this_group['massWeighted%sMap'%quantity_name.title()]=massWeightedQuantityMap

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

        else:
            ## appending another quantity, cool!
            this_group = h5file[this_setup_id]
            assert "massWeighted%sMap"%quantity_name.title() not in this_group.keys() 
            ## save this new quantity
            this_group['massWeighted%sMap'%quantity_name.title()]=massWeightedQuantityMap
