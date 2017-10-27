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
import getopt
from io import *
import numpy as np 
#import tables
import h5py
from numpy import zeros
from numpy import where
from numpy import ravel
from numpy import shape
from numpy import array
from numpy import ndarray
from numpy import float32
from numpy import arange
from numpy import uint8
from numpy import uint16
from numpy import cos
from numpy import sin
from numpy import reshape
from numpy import log10 as np_log10
from numpy.ctypeslib import ndpointer
from ctypes import *
import matplotlib
#matplotlib.use('Agg') 
import math
import pylab

class colour_table_class:
    
    index  = zeros(256,dtype=uint8)
    red    = zeros(256,dtype=uint8)
    green  = zeros(256,dtype=uint8)
    blue   = zeros(256,dtype=uint8)

    def __init__(self):
        return

def round_to_nearest_integer(x):
    delta_x = x - int(x)
    if delta_x < 0.5:
        return int(x)
    else:
        return int(x) + 1

def read_colour_table(file):

    ff = open(file)
    lines = ff.readlines()
    ff.close()
    
    colour_table = colour_table_class()

    for i in lines:
        if i[0]!='#':
            idx,r,g,b = string.split(i)
            idx = uint8(idx)
            r   = uint8(r)
            g   = uint8(g)
            b   = uint8(b)
            colour_table.index[idx] = idx
            colour_table.red[idx]   = r
            colour_table.green[idx] = g
            colour_table.blue[idx]  = b

    return colour_table

def new_matrix(m,n):
    # Create zero matrix
    matrix = [[0 for row in range(n)] for col in range(m)]
    return matrix
           

def matmul(matrix1,matrix2):
    # Matrix multiplication
    if len(matrix1[0]) != len(matrix2):
        # Check matrix dimensions
        print 'Matrices must be m*n and n*p to multiply!'
        print len(matrix1[0]),len(matrix2)
    else:
        # Multiply if correct dimensions
        my_matrix = new_matrix(len(matrix1),len(matrix2[0]))
        for i in range(len(matrix1)):
            for j in range(len(matrix2[0])):
                for k in range(len(matrix2)):
                    my_matrix[i][j] += matrix1[i][k]*matrix2[k][j]
        return my_matrix

def rotateEuler(theta,phi,psi,pos):
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
    print 'theta = ',theta_rad
    print 'phi   = ',phi_rad
    print 'psi   = ',psi_rad
    rot_matrix[0][0] =  cos(phi_rad)*cos(psi_rad)
    rot_matrix[0][1] = -cos(phi_rad)*sin(psi_rad)
    rot_matrix[0][2] =  sin(phi_rad)
    rot_matrix[1][0] =  cos(theta_rad)*sin(psi_rad) + sin(theta_rad)*sin(phi_rad)*cos(psi_rad)
    rot_matrix[1][1] =  cos(theta_rad)*cos(psi_rad) - sin(theta_rad)*sin(phi_rad)*sin(psi_rad)
    rot_matrix[1][2] = -sin(theta_rad)*cos(phi_rad)
    rot_matrix[2][0] =  sin(theta_rad)*sin(psi_rad) - cos(theta_rad)*sin(phi_rad)*cos(psi_rad)
    rot_matrix[2][1] =  sin(theta_rad)*cos(psi_rad) - cos(theta_rad)*sin(phi_rad)*sin(psi_rad)
    rot_matrix[2][2] =  cos(theta_rad)*cos(phi_rad)
    print rot_matrix

    # translate particles so centre = origin
    pos[:,0] -= (Xmin + (L_x/2.))
    pos[:,1] -= (Ymin + (L_y/2.))
    pos[:,2] -= (Zmin + (L_z/2.))

    # rotate about each axis with a matrix operation
    pos_rot = ndarray((n_box, 3), dtype=float32)
    for ipart in range(n_box):
        pos_matrix = new_matrix(3,1)
        pos_matrix[0][0] = pos[ipart,0]
        pos_matrix[1][0] = pos[ipart,1]
        pos_matrix[2][0] = pos[ipart,2]
        rotated = matmul(rot_matrix,pos_matrix)
        pos_rot[ipart,0] = rotated[0][0]
        pos_rot[ipart,1] = rotated[1][0]
        pos_rot[ipart,2] = rotated[2][0]
    
    # translate particles back
    pos[:,0] += (Xmin + (L_x/2.))
    pos[:,1] += (Ymin + (L_y/2.))
    pos[:,2] += (Zmin + (L_z/2.))
    
    return pos_rot

def writeImageGrid(isnap,output_dir,savename,
    ResultQ,ResultQ_edge,
    image_length,npix_x,time_Myr,HubbleParam,edgeon=0):
    array_name = "ResultQ" 
    # Write the image grid to HDF5 file 
    h5filename = "gasTemp_proj_%.3d_%.2fkpc.hdf5" % (isnap, image_length)

    output_name = "%s_faceOn" % (array_name, )

    h5name=output_dir + h5filename#savename
    with h5py.File(h5name, "w") as h5file:
        #h5file.createArray(h5file.root, output_name, ResultQ)
        h5file[output_name]=ResultQ

        if edgeon:
            # Write edge on image grid to HDF5 file. 
            output_name = "%s_edgeOn" % (array_name, ) 
            #h5file.createArray(h5file.root, output_name, ResultQ_edge)
            h5file[output_name]=ResultQ_edge
        try: 
            #h5file.createArray(h5file.root, "image_length", np.array([image_length]))
            #h5file.createArray(h5file.root, "npix_x", np.array([npix_x]))
            #h5file.createArray(h5file.root, "Time_Myr", np.array([time_Myr]))
            #h5file.createArray(h5file.root, "HubbleParam", np.array([HubbleParam]))

            h5file['image_length']=np.array([image_length])
            h5file['npix_x']=np.array([npix_x])
            h5file['Time_Myr']=np.array([time_Myr])
            h5file['HubbleParam']=np.array([HubbleParam])

        except:
            print "name already exists"
            raise Exception("overwriting!!")
        h5file.close()


def compute_image_grid(pos_all,mass_all,temperature_all,HubbleParam,time_Myr,BoxSize,
    frame_center,frame_width,frame_depth,
    savename,output_dir,isnap,
    edgeon=1,
    theta=0,phi=0,psi=0,
    pixels=1200,min_temp=2,max_temp=7,**kwargs):

    print "extra kwargs in compute_temp:",kwargs.keys()
    print ' rotation = (',theta,',',phi,',',psi,')'
  
    # Set paths

    # Set parameters
    npix_x   = pixels#1200
    npix_y   = pixels#1200
    desngb   = 32
    Axis1    = 0
    Axis2    = 1
    Axis3    = 2

    """
    # cgs units
    u_all *= 1.0e10   # cgs                                                                                                                                                                   

    XH = 0.76
    yhelium = (1.0 - XH) / (4.0 * XH)
    mu_all = (1.0 + (4.0 * yhelium)) / (1.0 + yhelium + elec_all)

    temperature_all = (2.0 / 3.0) * mu_all * 1.67e-24 * u_all / 1.38e-16
    """

    npart    = len(mass_all)
    hsml_all = ndarray(npart,dtype=float32)
    hsml_all[:] = 0.0

    image_length = 2*frame_width

    global Xmin,Xmax,Ymin,Ymax,Zmin,Zmax
    global L_x,L_y,L_z
    L_x     = image_length * HubbleParam	   #Size of image in x (in kpc / h)
    L_y     = image_length * HubbleParam 	   #Size of image in y
    L_z     = image_length * HubbleParam / 2.0
    Xmin    = -image_length * HubbleParam / 2.0	+frame_center[0]
    Xmax    = image_length * HubbleParam / 2.0	+frame_center[0] 
    Ymin    = -image_length * HubbleParam / 2.0 +frame_center[1]
    Ymax    = image_length * HubbleParam / 2.0 	+frame_center[1]
    Zmin    = -frame_depth + frame_center[2]#-image_length * HubbleParam / 4.0 
    Zmax    = frame_depth + frame_center[2]#image_length * HubbleParam / 4.0 

    print 'extracting cube'
    
    ind_box = ((pos_all[:,0] > Xmin) & (pos_all[:,0] < Xmax) &
               (pos_all[:,1] > Ymin) & (pos_all[:,1] < Ymax) &
               (pos_all[:,2] > Zmin) & (pos_all[:,2] < Zmax))
    global n_box
    n_box = sum(ind_box)
    
    print 'n_box = ',n_box	#This is the number of particles in the box.
    
    pos = ndarray((n_box, 3), dtype=float32)
    pos[:,0] = pos_all[ind_box,0]
    pos[:,1] = pos_all[ind_box,1]
    pos[:,2] = pos_all[ind_box,2]
    mass     = mass_all[ind_box].astype(float32) 
    hsml     = hsml_all[ind_box].astype(float32)
    temperature = temperature_all[ind_box].astype(float32)
    quantity = temperature 
    Hmax     = 0.5*L_y	
    n_smooth = n_box

    print '-done'

    # rotate particles by angle derived from frame number
    ## rotate by euler angles if necessary
    pos = rotateEuler(theta,phi,psi,pos)

    ResultW = ndarray((npix_x,npix_y),dtype=float32)
    ResultQ = ndarray((npix_x,npix_y),dtype=float32)
    ResultW[:,:] = 0.0
    ResultQ[:,:] = 0.0
    
    c_f_p      = POINTER(c_float)


    c_f_p      = POINTER(c_float)
    pos_p      = pos.ctypes.data_as(c_f_p)
    hsml_p     = hsml.ctypes.data_as(c_f_p)
    mass_p     = mass.ctypes.data_as(c_f_p)
    quantity_p = quantity.ctypes.data_as(c_f_p)
    w_f_p    = ResultW.ctypes.data_as(c_f_p)
    q_f_p    = ResultQ.ctypes.data_as(c_f_p)


    print '------------------------------------------'
    c_obj = CDLL('/home/abg6257/starformation/movie_maker/HsmlAndProject_cubicSpline/HsmlAndProject.so')
    c_obj.findHsmlAndProject(c_int(n_smooth),pos_p,hsml_p,mass_p,quantity_p,\
                             c_float(Xmin),c_float(Xmax),c_float(Ymin),c_float(Ymax),c_float(Zmin),c_float(Zmax),\
                             c_int(npix_x),c_int(npix_y),c_int(desngb),\
                             c_int(Axis1),c_int(Axis2),c_int(Axis3),c_float(Hmax),c_double(BoxSize),w_f_p,q_f_p)
    print '------------------------------------------'

    # ResultQ contains the mass-weighted temperature, in K. 

    ResultQ = np_log10(ResultQ)
   
    print 'log10 minmax(ResultQ)',min(ravel(ResultQ)),max(ravel(ResultQ))


    # EDGE ON orientation.    
    # Halve number of y pixels
    npix_y/= 2

    ResultW_edgeOn= ndarray((npix_x,npix_y),dtype=float32)
    ResultQ_edgeOn= ndarray((npix_x,npix_y),dtype=float32)
    ResultW_edgeOn[:,:] = 0.0
    ResultQ_edgeOn[:,:] = 0.0

    if not edgeon:
        # free up some memory
        print '-releasing original snapshot memory'
        del pos_all
        del mass_all
        del hsml_all
        del temperature_all 
        print '-done'
    else:
        ind_box = ((pos_all[:,0] > Xmin) & (pos_all[:,0] < Xmax) &
                   (pos_all[:,1] > Ymin) & (pos_all[:,1] < Ymax) &
                   (pos_all[:,2] > Zmin) & (pos_all[:,2] < Zmax))
        n_box = sum(ind_box)
        
        print 'n_box = ',n_box	#This is the number of particles in the box.
        
        pos = ndarray((n_box, 3), dtype=float32)
        pos[:,0] = pos_all[ind_box,0]
        pos[:,1] = pos_all[ind_box,1]
        pos[:,2] = pos_all[ind_box,2]
        mass     = mass_all[ind_box].astype(float32) 
        hsml     = hsml_all[ind_box].astype(float32)
        temperature = temperature_all[ind_box].astype(float32)
        quantity = temperature 
        Hmax     = 0.5*L_y	
        n_smooth = n_box

        print '-done'

        # free up some memory
        print '-releasing original snapshot memory'
        del pos_all
        del mass_all
        del hsml_all
        del temperature_all 
        print '-done'

        # rotate particles by 90 degrees 
        pos = rotateEuler(-90,0,0,pos)

        c_f_p      = POINTER(c_float)
        pos_p      = pos.ctypes.data_as(c_f_p)
        hsml_p     = hsml.ctypes.data_as(c_f_p)
        mass_p     = mass.ctypes.data_as(c_f_p)
        quantity_p = quantity.ctypes.data_as(c_f_p)
        w_f_p    = ResultW_edgeOn.ctypes.data_as(c_f_p)
        q_f_p    = ResultQ_edgeOn.ctypes.data_as(c_f_p)

        print '------------------------------------------'
        c_obj = CDLL('/home/abg6257/starformation/movie_maker/HsmlAndProject_cubicSpline/HsmlAndProject.so')
        c_obj.findHsmlAndProject(c_int(n_smooth),pos_p,hsml_p,mass_p,quantity_p,\
                                 c_float(Xmin),c_float(Xmax),c_float(Zmin),c_float(Zmax),c_float(Ymin),c_float(Ymax),\
                                 c_int(npix_x),c_int(npix_y),c_int(desngb),\
                                 c_int(Axis1),c_int(Axis2),c_int(Axis3),c_float(Hmax),c_double(BoxSize),w_f_p,q_f_p)
        print '------------------------------------------'

        ResultQ_edgeOn = np_log10(ResultQ_edgeOn)
       
        print 'log10 minmax(ResultQ_edgeOn)',min(ravel(ResultQ_edgeOn)),max(ravel(ResultQ_edgeOn))

    writeImageGrid(isnap,output_dir,savename,ResultQ,ResultQ_edgeOn,image_length,npix_x,time_Myr,HubbleParam,edgeon=edgeon)
    
    return 

def plot_image_grid(isnap, sim_name):

    # Set paths
    sim_dir     = './%s/' % (sim_name, )	#Make this the data directory.
    
    my_cbar_label = r"$\log \, T \, (\rm{K})$"
    data_dir = sim_dir + 'Plots/GasTemperature/' 
    output_dir  = sim_dir + 'Plots/GasTemperature/standard/' 

    # Linear colormap (rainbow) 
    cols = [(0,0,0)]
    for y in np.linspace(0,1, 254):
        # 2 < log T < 6 
#        if y < 0.5:
#            x = (4.0 / 3.0) * (y ** 2.0) 
#        else:
#            x = (4.0 * y / 3.0) - (1.0 / 3.0)

        # 2 < log T < 7
        if y < 0.4:
            x = (125.0 / 24.0) * (y ** 3.0)
        else:
            x = ((10.0 / 9.0) * y) - (1.0 / 9.0) 
            
        rcol = (0.472-0.567*x+4.05*x**2)/(1.+8.72*x-19.17*x**2+14.1*x**3)
        gcol = 0.108932-1.22635*x+27.284*x**2-98.577*x**3+163.3*x**4-131.395*x**5+40.634*x**6
        bcol = 1./(1.97+3.54*x-68.5*x**2+243*x**3-297*x**4+125*x**5)
        cols.append((rcol, gcol, bcol))
    cols.append((rcol, gcol, bcol))
    colour_map = matplotlib.colors.LinearSegmentedColormap.from_list("PaulT_rainbow", cols)
    colour_map.set_under('k')

    # Set parameters
    npix_x   = 1200 
    npix_y   = 1200 
    image_length = 8.0
    scale_label_position = 0.06 
    scale_line_length = 1.0 
    scale_label_text = r"$\mathbf{1 \, \rm{kpc}}$"
    cbar_ticks = [2.0, 3.0, 4.0, 5.0, 6.0, 7.0] 
    tf_min = min_temp
    tf_max = max_temp

    print 'tf_min = ',tf_min
    print 'tf_max = ',tf_max

    h5filename = "gasTemperature_projection_%.3d_%.2fkpc_%dpx.hdf5" % (isnap, image_length, npix_x)
    h5file = tables.openFile(data_dir + h5filename, "r")
    
    # Face on image 
    ResultQ = float32(h5file.root.ResultQ_faceOn.read())
    try:
        HubbleParam = h5file.root.HubbleParam.read()[0] 
    except:
        HubbleParam=1
    time_Myr = h5file.root.Time_Myr.read()[0] 

    Xmin    = -image_length * HubbleParam / 2.0	   
    Xmax    = image_length * HubbleParam / 2.0	   
    Ymin    = -image_length * HubbleParam / 2.0    
    Ymax    = image_length * HubbleParam / 2.0 	   
    
#    ResultQ = ResultQ - tf_min 
#    ResultQ = ResultQ / (tf_max - tf_min)
#    ResultQ[ResultQ < 0.0] = 0.0
#    ResultQ[ResultQ > 1.0] = 1.0
#    ResultQ = ResultQ*255.0

#    print 'Image range (8bit): ',min(ravel(ResultQ)),max(ravel(ResultQ))

#    ResultQ = ResultQ.astype(uint16)    

    # translate ResultQ into a 3-channel array using the colour table
    # Also, the array 'image' has twice as many x-pixels, to accommodate
    # the other orientation. 
    image = ndarray((1.5*npix_y,npix_x),dtype=float32)
    for i in range(0,npix_y):
        for j in range(0,npix_x):
            value = ResultQ[j,i]
            image[i + (0.5 * npix_y),j] = value
            image[i + (0.5 * npix_y),j] = value 
            image[i + (0.5 * npix_y),j] = value 
            
    # Edge on image 
    npix_y /= 2 

    ResultQ = float32(h5file.root.ResultQ_edgeOn.read())
    
#    ResultQ = ResultQ - tf_min			
#    ResultQ = ResultQ / (tf_max - tf_min)
#    ResultQ[ResultQ < 0.0] = 0.0
#    ResultQ[ResultQ > 1.0] = 1.0
#    ResultQ = ResultQ*255.0

#    print 'Image range (8bit): ',min(ravel(ResultQ)),max(ravel(ResultQ))

#    ResultQ = ResultQ.astype(uint16) 
            
    for i in range(0,npix_y):
        for j in range(0,npix_x):
            value = ResultQ[j,i]
            image[i,j] = value 
            image[i,j] = value
            image[i,j] = value
            
    fig = pylab.figure(figsize = (npix_x / 600.0, 3.0 * npix_y / 600.0), dpi=600, frameon=False)
    ax = pylab.Axes(fig, [0,0,1,1])
    ax.set_axis_off()
    fig.add_axes(ax)

    # We need to create the scale bar

    # Length of bar in code units
    scale_line_length *= HubbleParam
    # Convert to pixels
    length_per_pixel = (Xmax - Xmin) / npix_x
    scale_line_length_px = int(scale_line_length / length_per_pixel)

    # Position in terms of image array indices
    scale_line_x_start = int(0.05 * npix_x)
    scale_line_x_end = scale_line_x_start + scale_line_length_px
    scale_line_y = int(0.2 * npix_y)

    # Go through pixels for scale bar, setting them to white
    for x_index in xrange(scale_line_x_start, scale_line_x_end):
        image[scale_line_y, x_index] = 0
        image[scale_line_y + 1, x_index] = 0
        image[scale_line_y + 2, x_index] = 0
        image[scale_line_y + 3, x_index] = 0
        image[scale_line_y + 4, x_index] = 0
        image[scale_line_y + 5, x_index] = 0

    figure_label2 = r"$\rm{UVBthin}$"
    
    imgplot = ax.imshow(image, extent = (Xmin,Xmax,1.5*Ymin,1.5*Ymax),origin = 'lower', aspect = 'auto', cmap = colour_map)
    imgplot.set_clim(tf_min, tf_max)
    figure_label = r"$%03d \, \rm{Myr}$" % (round_to_nearest_integer(time_Myr), )
    label = pylab.text(0.70, 0.92, figure_label, fontsize = 8, transform = ax.transAxes)
    label.set_color('k')
    label2 = pylab.text(scale_label_position, 0.03, scale_label_text, fontweight = 'bold', transform = ax.transAxes)
    label2.set_color('k')
    label2.set_fontsize(6)
    #label3 = pylab.text(0.10, 0.92, figure_label2, fontsize = 8, transform = ax.transAxes)
    #label3.set_color('k')

    # colour bar
    ax_c = pylab.axes([0.15, 0.415, 0.7, 0.015])
    norm = matplotlib.colors.Normalize(vmin = tf_min, vmax = tf_max)
    cbar = matplotlib.colorbar.ColorbarBase(ax_c, cmap = colour_map, norm = norm, orientation = 'horizontal') 
    cbar.set_ticks(cbar_ticks)
    cbar_tick_labels = []
    for i in cbar_ticks:
        next_label = r"$\mathbf{%.1f}$" % (i, )
        cbar_tick_labels.append(next_label)
    cbar.set_ticklabels(cbar_tick_labels)
    ax_c.spines["bottom"].set_linewidth(1.0)
    ax_c.spines["top"].set_linewidth(1.0)
    ax_c.spines["left"].set_linewidth(1.0)
    ax_c.spines["right"].set_linewidth(1.0)
    ax_c.xaxis.set_tick_params(width=0.5, length=1.2)
    for t in ax_c.get_xticklabels():
        t.set_fontsize(5)
        t.set_fontweight('bold')
        t.set_color('k')
        t.set_y(t.get_position()[1] + 0.5)
    cbar.outline.set_linewidth(0.4)
    cbar.set_label(my_cbar_label, color = 'k', fontsize=6, fontweight='bold', labelpad = 0.5)

    print cbar.get_clim()

    '''ax_c = pylab.axes([0.05, 0.2, 0.9, 0.02])
    norm = matplotlib.colors.Normalize(vmin = tf_min, vmax = tf_max)
    cbar = matplotlib.colorbar.ColorbarBase(ax_c, cmap = colour_map, norm = norm, orientation = 'horizontal') 
    cbar.set_ticks(cbar_ticks)
    cbar_tick_labels = []
    for i in cbar_ticks:
        next_label = r"$\mathbf{%.1f}$" % (i, )
        cbar_tick_labels.append(next_label)
    cbar.set_ticklabels(cbar_tick_labels)
    ax_c.spines["bottom"].set_linewidth(1.0)
    ax_c.spines["top"].set_linewidth(1.0)
    ax_c.spines["left"].set_linewidth(1.0)
    ax_c.spines["right"].set_linewidth(1.0)
    ax_c.xaxis.set_tick_params(width=0.5, length=1.5)
    for t in ax_c.get_xticklabels():
        t.set_fontsize(5)
        t.set_fontweight('bold')
        t.set_color('k')
        t.set_y(t.get_position()[1] + 0.5)
    cbar.outline.set_linewidth(0.4)
    cbar.set_label(my_cbar_label, color = 'k', fontsize=6, fontweight='bold', labelpad = 0.5)'''

    extent = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    exec "pylab.savefig('%s/frame%04d.png', dpi = 600, bbox_inches=extent)" % (output_dir, isnap)
    print '-done'
    print ''
    pylab.close()

    h5file.close() 

    return 

'''def main():
    isnap_low = int(sys.argv[1]) 
    isnap_hi = int(sys.argv[2]) 
    sim_name = sys.argv[3]   # of the form 'data_dir' 
    chimes_mode = int(sys.argv[4]) # 0 - standard GIZMO cooling 
                                   # 1 - CHIMES 
    plot_mode = int(sys.argv[5])   # 0 - compute image grid 
                                   # 1 - plot image grid 

    for isnap in xrange(isnap_low, isnap_hi):
        if plot_mode == 0: 
            # compute temperature maps. 
            compute_image_grid(isnap, sim_name, chimes_mode) 
        elif plot_mode == 1: 
            plot_image_grid(isnap, sim_name) 
        else:
            print "ERROR: plot mode not recognised. " 

    return 

main() '''
