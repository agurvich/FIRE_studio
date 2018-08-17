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

class colour_table_class:
    
    index  = np.zeros(256,dtype=np.uint8)
    red    = np.zeros(256,dtype=np.uint8)
    green  = np.zeros(256,dtype=np.uint8)
    blue   = np.zeros(256,dtype=np.uint8)

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
        print('Matrices must be m*n and n*p to multiply!')
        print(len(matrix1[0]),len(matrix2))
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
    print('theta = ',theta_rad)
    print('phi   = ',phi_rad)
    print('psi   = ',psi_rad)
    rot_matrix[0][0] =  cos(phi_rad)*cos(psi_rad)
    rot_matrix[0][1] = -cos(phi_rad)*sin(psi_rad)
    rot_matrix[0][2] =  sin(phi_rad)
    rot_matrix[1][0] =  cos(theta_rad)*sin(psi_rad) + sin(theta_rad)*sin(phi_rad)*cos(psi_rad)
    rot_matrix[1][1] =  cos(theta_rad)*cos(psi_rad) - sin(theta_rad)*sin(phi_rad)*sin(psi_rad)
    rot_matrix[1][2] = -sin(theta_rad)*cos(phi_rad)
    rot_matrix[2][0] =  sin(theta_rad)*sin(psi_rad) - cos(theta_rad)*sin(phi_rad)*cos(psi_rad)
    rot_matrix[2][1] =  sin(theta_rad)*cos(psi_rad) - cos(theta_rad)*sin(phi_rad)*sin(psi_rad)
    rot_matrix[2][2] =  cos(theta_rad)*cos(phi_rad)
    print(rot_matrix)

    # translate particles so centre = origin
    pos[:,0] -= (Xmin + (L_x/2.))
    pos[:,1] -= (Ymin + (L_y/2.))
    pos[:,2] -= (Zmin + (L_z/2.))

    # rotate about each axis with a matrix operation
    pos_rot = np.ndarray((n_box, 3), dtype=np.float32)
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

def writeImageGrid(isnap,output_dir,
    ResultW,ResultW_edge,
    image_length,npix_x,time_Myr,edgeon=0,h5filename=''):
    array_name = "ResultW_tot" 
    # Write the image grid to HDF5 file 
    h5filename += "gas_proj_%.3d_%.2fkpc.hdf5" % (isnap, image_length)

    output_name = "%s_faceOn" % (array_name, )

    h5name=output_dir + h5filename
    with h5py.File(h5name, "w") as h5file:
        h5file[output_name]=ResultW

        if edgeon:
            output_name = "%s_edgeOn" % (array_name, ) 
            h5file[output_name]=ResultW_edge

        try: 
            h5file['image_length']=np.array([image_length])
            h5file['npix_x']=np.array([npix_x])
            h5file['Time_Myr']=np.array([time_Myr])
        except:
            print("name already exists")
            raise Exception("overwriting!!")
        h5file.close()
        
    
def getImageGrid(BoxSize,npix_x,npix_y,pos,mass,hsml,quantity):
    ResultW = np.ndarray((npix_x,npix_y),dtype=np.float32)
    ResultQ = np.ndarray((npix_x,npix_y),dtype=np.float32)
    ResultW[:,:] = 0.0
    ResultQ[:,:] = 0.0
    
    c_f_p      = ctypes.POINTER(ctypes.c_float)
    pos_p      = pos.ctypes.data_as(c_f_p)
    hsml_p     = hsml.ctypes.data_as(c_f_p)
    mass_p     = mass.ctypes.data_as(c_f_p)
    quantity_p = quantity.ctypes.data_as(c_f_p)
    w_f_p    = ResultW.ctypes.data_as(c_f_p)
    q_f_p    = ResultQ.ctypes.data_as(c_f_p)

    print('------------------------------------------')
    curpath = os.path.realpath(__file__)
    curpath = curpath[:len("utils")+curpath.index("utils")] #split off this filename
    c_obj = ctypes.CDLL(os.path.join(curpath,'gas_utils','HsmlAndProject_cubicSpline/HsmlAndProject.so'))
    c_obj.findHsmlAndProject(
	ctypes.c_int(n_smooth),
	pos_p,hsml_p,mass_p,quantity_p,
        ctypes.c_float(Xmin),ctypes.c_float(Xmax),
	ctypes.c_float(Ymin),ctypes.c_float(Ymax),
	ctypes.c_float(Zmin),ctypes.c_float(Zmax),
        ctypes.c_int(npix_x),ctypes.c_int(npix_y),
	ctypes.c_int(desngb),
        ctypes.c_int(Axis1),ctypes.c_int(Axis2),ctypes.c_int(Axis3),
	ctypes.c_float(Hmax),ctypes.c_double(BoxSize),w_f_p,q_f_p)
    print('------------------------------------------')

    # normalise by area
    Lcell    = (L_x)/float(npix_x)
    ResultW /= (Lcell*Lcell)             # 10^10 Msun / kpc^-2 
    
    
    # convert into Msun/pc^2
    unitmass_in_g = 1.9890000e+43 
    solar_mass    = 1.9890000e+33
    conv_fac = (unitmass_in_g/solar_mass) / (1.0e6)
    ResultW *= conv_fac

    ResultW = np.log10(ResultW)
   
    print('log10 minmax(ResultW)',min(np.ravel(ResultW)),max(np.ravel(ResultW)))
    return ResultW

def compute_image_grid(pos_all,mass_all,time_Myr,BoxSize,
    frame_center,frame_half_width,frame_depth,
    output_dir,isnap,
    edgeon=1,
    theta=0,phi=0,psi=0,
    pixels=1200,min_den=-1.0,max_den=1.2,h5filename='',**kwargs):

    print("extra kwargs in compute_den:",kwargs.keys())
    print(' rotation = (',theta,',',phi,',',psi,')')

    # Set parameters
    npix_x   = pixels#1200
    npix_y   = pixels#1200
    global desngb,Axis1,Axis2,Axis3
    desngb   = 32
    Axis1    = 0
    Axis2    = 1
    Axis3    = 2

    npart    = len(mass_all)
    hsml_all = np.ndarray(npart,dtype=np.float32)
    hsml_all[:] = 0.0

    image_length = 2*frame_half_width
    image_depth = 2*frame_depth

    global L_x,L_y,L_z
    L_x = image_length  #Size of image in x (in kpc)
    L_y = image_length  #Size of image in y (in kpc)
    L_z = image_depth   #Size of image in z (in kpc)
    global Xmin,Xmax,Ymin,Ymax,Zmin,Zmax
    Xmin = -image_length/2.0 + frame_center[0]
    Xmax = image_length/2.0 + frame_center[0]
    Ymin = -image_length/2.0 + frame_center[1]
    Ymax = image_length/2.0 + frame_center[1]
    Zmin = -frame_depth + frame_center[2]
    Zmax = frame_depth + frame_center[2] 

    print('extracting cube')
    
    ind_box = ((pos_all[:,0] > Xmin) & (pos_all[:,0] < Xmax) &
               (pos_all[:,1] > Ymin) & (pos_all[:,1] < Ymax) &
               (pos_all[:,2] > Zmin) & (pos_all[:,2] < Zmax))
    global n_box
    n_box = np.sum(ind_box)
    
    print('n_box = ',n_box)	#This is the number of particles in the box.
    
    pos = np.ndarray((n_box, 3), dtype=np.float32)
    pos[:,0] = pos_all[ind_box,0]
    pos[:,1] = pos_all[ind_box,1]
    pos[:,2] = pos_all[ind_box,2]

    mass     = mass_all[ind_box].astype(np.float32)
    hsml     = hsml_all[ind_box].astype(np.float32)
    quantity = mass
    global Hmax
    Hmax     = 0.5*L_y	
    global n_smooth
    n_smooth = n_box

    print('-done')

    ## rotate by euler angles if necessary
    pos = rotateEuler(theta,phi,psi,pos)

    ResultW = getImageGrid(BoxSize,npix_x,npix_y,pos,mass,hsml,quantity)

    # EDGE ON orientation.    
    # Halve number of y pixels
    npix_y   /= 2

    if not edgeon:
        # free up some memory
        print('-releasing original snapshot memory')
        del pos_all
        del mass_all
        del hsml_all
        print('-done')

        ## set them just to 0
        ResultW_edge = np.ndarray((npix_x,npix_y),dtype=np.float32)
        ResultW_edge[:,:] = 0.0
        
    else:
        ind_box = ((pos_all[:,0] > Xmin) & (pos_all[:,0] < Xmax) &
                   (pos_all[:,1] > Ymin) & (pos_all[:,1] < Ymax) &
                   (pos_all[:,2] > Zmin) & (pos_all[:,2] < Zmax))
        n_box = sum(ind_box)
        
        print('n_box = ',n_box)#This is the number of particles in the box.
        
        pos = np.ndarray((n_box, 3), dtype=np.float32)
        pos[:,0] = pos_all[ind_box,0]
        pos[:,1] = pos_all[ind_box,1]
        pos[:,2] = pos_all[ind_box,2]
        mass     = mass_all[ind_box].astype(np.float32)
        hsml     = hsml_all[ind_box].astype(np.float32)
        quantity = mass
        Hmax     = 0.5*L_y	
        n_smooth = n_box

        print('-done')

        # free up some memory
        print('-releasing original snapshot memory')
        del pos_all
        del mass_all
        del hsml_all
        print('-done')

        # rotate particles by 90 degrees to do edge on view
        pos = rotateEuler(-90,0,0,pos)
        Ymin,Ymax,Zmin,Zmax=Zmin,Zmax,Ymin,Ymax

        ResultW_edge = getImageGrid(BoxSize,npix_x,npix_y,pos,mass,hsml,quantity)
        print('log10 minmax(ResultW_edge)',min(np.ravel(ResultW_edge)),max(np.ravel(ResultW_edge)))
            
    writeImageGrid(
        isnap,
        output_dir,ResultW,ResultW_edge,
        image_length,npix_x,time_Myr,edgeon=edgeon,h5filename=h5filename)


####### ignore beyond here
"""
from matplotlib import colors, colorbar 
def plot_image_grid(isnap, sim_name, species_index):

    # Set paths
    sim_dir     = './%s/' % (sim_name, )	# Make this the data directory.
    
    if species_index == 999:
        array_name = "ResultW_tot" 
        my_cbar_label = r"$\mathbf{\log \, \Sigma_{\rm{gas}} \, (\rm{M}_{\odot} \, \rm{pc}^{-2})}$"
        output_dir  = sim_dir + 'Plots/GasDensity/totGas/'    
    elif species_index == 1:
        array_name = "ResultW_HI" 
        my_cbar_label = r"$\mathbf{\log \, \Sigma_{\rm{HI}} \, (\rm{M}_{\odot} \, \rm{pc}^{-2})}$"
        output_dir  = sim_dir + 'Plots/GasDensity/HI/'    
    elif species_index == 2:
        array_name = "ResultW_HII"
        my_cbar_label = r"$\mathbf{\log \, \Sigma_{\rm{HII}} \, (\rm{M}_{\odot} \, \rm{pc}^{-2})}$"
        output_dir  = sim_dir + 'Plots/GasDensity/HII/' 
    elif species_index == 8:
        array_name = "ResultW_CII"
        my_cbar_label = r"$\mathbf{\log \, \Sigma_{\rm{CII}} \, (\rm{M}_{\odot} \, \rm{pc}^{-2})}$"
        output_dir  = sim_dir + 'Plots/GasDensity/CII/' 
    elif species_index == 137: 
        array_name = "ResultW_H2" 
        my_cbar_label = r"$\mathbf{\log \, \Sigma_{\rm{H_{2}}} \, (\rm{M}_{\odot} \, \rm{pc}^{-2})}$"
        output_dir  = sim_dir + 'Plots/GasDensity/H2/' 
    elif species_index == 148: 
        array_name = "ResultW_CO" 
        my_cbar_label = r"$\mathbf{\log \, \Sigma_{\rm{CO}} \, (\rm{M}_{\odot} \, \rm{pc}^{-2})}$"
        output_dir  = sim_dir + 'Plots/GasDensity/CO/' 
    elif species_index == 300:
        array_name = "ResultW_Htot" 
        my_cbar_label = r"$\mathbf{\log \, \Sigma_{\rm{H,} \, \rm{tot}} \, (\rm{M}_{\odot} \, \rm{pc}^{-2})}$"
        output_dir  = sim_dir + 'Plots/GasDensity/Htot/' 
    else: 
        print("ERROR: species_index not recognised. ")

    data_dir = sim_dir + 'Plots/GasDensity/' 

    ct_file = '/home/ajr882/analysis/colour_tables/REDTEMPERATURE.ct'
    colour_table = read_colour_table(ct_file)
             
    cols = []
    for x in range(256):
        rcol = float(colour_table.red[x]) / 255.
        gcol = float(colour_table.green[x]) / 255.
        bcol = float(colour_table.blue[x]) / 255.
        cols.append((rcol, gcol, bcol))
    colour_map = colors.LinearSegmentedColormap.from_list("PaulT_rainbow", cols)
    colour_map.set_under('k')

    # Set parameters
    npix_x   = 1200 
    npix_y   = 1200 
    image_length = 8.0
    scale_label_position = 0.06 
    scale_line_length = 1.0 
    scale_label_text = r"$\mathbf{1 \, \rm{kpc}}$"

    tf_min = -0.4 
    tf_max = 1.6 

    cbar_ticks = np.arange(tf_min,tf_max+0.4,0.4)#[-0.4, 0.0, 0.4, 0.8, 1.2, 1.6] 


    print('tf_min = ',tf_min)
    print('tf_max = ',tf_max)

    h5filename = "gas_projection_%03d_%.2fkpc_%dpx.hdf5" % (isnap, image_length, npix_x)
    h5file = tables.openFile(data_dir + h5filename, "r")
    
    # Face on image 
    exec "ResultW = h5file.root.%s_faceOn.read()" % (array_name, )
    
    Xmin    = -image_length/2.0	   
    Xmax    = image_length/2.0	   
    Ymin    = -image_length/2.0    
    Ymax    = image_length/2.0 	   
    
    ResultW = ResultW - tf_min 
    ResultW = ResultW / (tf_max - tf_min)
    ResultW[ResultW < 0.0] = 0.0
    ResultW[ResultW > 1.0] = 1.0
    ResultW = ResultW*255.0

    print('Image range (8bit): ',min(np.ravel(ResultW)),max(np.ravel(ResultW)))

    ResultW = ResultW.astype(uint16)    

    # translate ResultW into a 3-channel array using the colour table
    # Also, the array 'image' has twice as many x-pixels, to accommodate
    # the other orientation. 
    image = np.ndarray((1.5*npix_y,npix_x,3),dtype=uint8)
    for i in range(0,npix_y):
        for j in range(0,npix_x):
            value = ResultW[j,i]
            image[i + (0.5 * npix_y),j,0] = colour_table.red[value]
            image[i + (0.5 * npix_y),j,1] = colour_table.green[value]
            image[i + (0.5 * npix_y),j,2] = colour_table.blue[value]
            
    # Edge on image 
    npix_y /= 2 

    ResultW = h5file['%s_edgeOn'% (array_name, )] 
    
    ResultW = ResultW - tf_min			
    ResultW = ResultW / (tf_max - tf_min)
    ResultW[ResultW < 0.0] = 0.0
    ResultW[ResultW > 1.0] = 1.0
    ResultW = ResultW*255.0

    print('Image range (8bit): ',min(np.ravel(ResultW)),max(np.ravel(ResultW)))

    ResultW = ResultW.astype(uint16) 
            
    for i in range(0,npix_y):
        for j in range(0,npix_x):
            value = ResultW[j,i]
            image[i,j,0] = colour_table.red[value]
            image[i,j,1] = colour_table.green[value]
            image[i,j,2] = colour_table.blue[value]
            
    import pylab
    fig = pylab.figure(figsize = (npix_x / 600.0, 3.0 * npix_y / 600.0), dpi=600, frameon=False)
    ax = pylab.Axes(fig, [0,0,1,1])
    ax.set_axis_off()
    fig.add_axes(ax)

    # We need to create the scale bar

    # Convert to pixels
    length_per_pixel = (Xmax - Xmin) / npix_x
    scale_line_length_px = int(scale_line_length / length_per_pixel)

    # Position in terms of image array indices
    scale_line_x_start = int(0.05 * npix_x)
    scale_line_x_end = scale_line_x_start + scale_line_length_px
    scale_line_y = int(0.2 * npix_y)

    # Go through pixels for scale bar, setting them to white
    for x_index in range(scale_line_x_start, scale_line_x_end):
        image[scale_line_y, x_index, 0] = 255
        image[scale_line_y, x_index, 1] = 255
        image[scale_line_y, x_index, 2] = 255
        image[scale_line_y + 1, x_index, 0] = 255
        image[scale_line_y + 1, x_index, 1] = 255
        image[scale_line_y + 1, x_index, 2] = 255
        image[scale_line_y + 2, x_index, 0] = 255
        image[scale_line_y + 2, x_index, 1] = 255
        image[scale_line_y + 2, x_index, 2] = 255
        image[scale_line_y + 3, x_index, 0] = 255
        image[scale_line_y + 3, x_index, 1] = 255
        image[scale_line_y + 3, x_index, 2] = 255
        image[scale_line_y + 4, x_index, 0] = 255
        image[scale_line_y + 4, x_index, 1] = 255
        image[scale_line_y + 4, x_index, 2] = 255
        image[scale_line_y + 5, x_index, 0] = 255
        image[scale_line_y + 5, x_index, 1] = 255
        image[scale_line_y + 5, x_index, 2] = 255

    #figure_label2 = r"$\rm{UVBthin}$"
        
    imgplot = ax.imshow(image, extent = (Xmin,Xmax,1.5*Ymin,1.5*Ymax),origin = 'lower', aspect = 'auto')
    figure_label = r"$%03d \, \rm{Myr}$" % (round_to_nearest_integer(time_Myr), )
    label = pylab.text(0.70, 0.92, figure_label, fontsize = 8, transform = ax.transAxes)
    label.set_color('white')
    label2 = pylab.text(scale_label_position, 0.03, scale_label_text, fontweight = 'bold', transform = ax.transAxes)
    label2.set_color('white')
    label2.set_fontsize(6)
    #label3 = pylab.text(0.10, 0.92, figure_label2, fontsize = 8, transform = ax.transAxes)
    #label3.set_color('white')

    # colour bar
    ax_c = pylab.axes([0.15, 0.415, 0.7, 0.015])
    norm = colors.Normalize(vmin = tf_min, vmax = tf_max)
    cbar = colorbar.ColorbarBase(ax_c, cmap = colour_map, norm = norm, orientation = 'horizontal') 
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
        t.set_color('w')
        t.set_y(t.get_position()[1] + 0.5)
    cbar.outline.set_linewidth(0.4)
    cbar.set_label(my_cbar_label, color = 'w', fontsize=6, fontweight='bold', labelpad = 0.5)

    print(cbar.get_clim())

    extent = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    exec "pylab.savefig('%s/frame%04d.png', dpi = 600, bbox_inches=extent)" % (output_dir, isnap)
    
    print('-done')
    pylab.close() 

    h5file.close() 

    return 
"""

'''def main():
    isnap_low = int(sys.argv[1]) 
    isnap_hi = int(sys.argv[2]) 
    sim_name = sys.argv[3]   # of the form 'data_dir' 
    chimes_flag = int(sys.argv[4])    # 0 - Standard GIZMO cooling 
                                      # 1 - CHIMES 
    plot_mode = int(sys.argv[5])   # 0 - compute image grid; given species only. 
                                   # 1 - plot image grid. 
    species_index = int(sys.argv[6]) 

    for isnap in range(isnap_low, isnap_hi):
        if plot_mode == 0:
            compute_image_grid(isnap, sim_name, eqm_mode, species_index)
        elif plot_mode == 1: 
            plot_image_grid(isnap, sim_name, eqm_mode, species_index)
        else:
            print "ERROR: plot mode not recognised. "

    return 


main()'''
