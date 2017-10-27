import matplotlib 
matplotlib.use('Agg') 
import pylab 
import numpy as np 
import my_colour_maps as mcm 
import sys 
import h5py

def round_to_nearest_integer(x):
    delta_x = x - int(x)
    if delta_x < 0.5:
        return int(x)
    else:
        return int(x) + 1

def plot_image_grid(ax,isnap,dprojects,tprojects,
    frame_center,frame_width,pixels=1200,
    min_den=-1.0,max_den=1.2,min_temp=2,max_temp=7,edgeon=0,**kwargs): 
    print "extra kwargs in plot_2color_image:",kwargs.keys()
    # Set paths
    
    array_name = "ResultW_tot" 

    data_dir_rho = dprojects #sim_dir + 'Plots/GasDensity/' 
    data_dir_T = tprojects #sim_dir + 'Plots/GasTemperature/'

    # Set parameters
    npix_x   = pixels#1200 
    npix_y   = pixels#1200 
    image_length =2*frame_width #8.0
    scale_label_position = 0.06 

    if image_length > 1.5 : 
        scale_line_length = 1.0 
        scale_label_text = r"$\mathbf{1 \, \rm{kpc}}$"
    else:
        scale_line_length = .1
        scale_label_text = r"$\mathbf{100 \, \rm{pc}}$"

    tf_min_rho = min_den
    tf_max_rho = max_den


    print 'tf_min_rho = ',tf_min_rho
    print 'tf_max_rho = ',tf_max_rho

    # First, read in density image arrays 
    h5filename = "gas_proj_%3d_%.2fkpc.hdf5" % (isnap, image_length)
    h5name=h5filename
    with h5py.File(data_dir_rho + h5name, "r") as h5file:
        exec "ResultW_rho = np.array(h5file['%s_faceOn'])" % (array_name, )
        try:
            HubbleParam = h5file['HubbleParam'][0]
        except:
            HubbleParam = 1

        try:
            time_Myr = h5file['Time_Myr'][0]
        except:
            print dir(h5file.root)
            print h5file.keys()
            time_Myr=h5file.root
            raise Exception("STOP!")

        Xmin    = -image_length * HubbleParam / 2.0	+frame_center[0]
        Xmax    = image_length * HubbleParam / 2.0	+frame_center[0]
        Ymin    = -image_length * HubbleParam / 2.0 +frame_center[1]
        Ymax    = image_length * HubbleParam / 2.0 	+frame_center[1]
        
        ResultW_rho = ResultW_rho - tf_min_rho 
        ResultW_rho = ResultW_rho / (tf_max_rho - tf_min_rho)
        ResultW_rho[ResultW_rho < 0.0] = 0.0
        ResultW_rho[ResultW_rho > 1.0] = 1.0
        ResultW_rho = ResultW_rho*255.0

        print 'Image range (8bit): ',min(np.ravel(ResultW_rho)),max(np.ravel(ResultW_rho))

        ResultW_rho = ResultW_rho.astype(np.uint16) 

        # Create the image array for both panels, and read 
        # ResultW_rho into it. 
        if not edgeon:
            image_rho = np.ndarray((npix_y,npix_x),dtype=np.uint16)
            for i in range(0,npix_y):
                for j in range(0,npix_x):
                    value = ResultW_rho[j,i]
                    image_rho[i,j] = value
        else:
            image_rho = np.ndarray((1.5*npix_y,npix_x),dtype=np.uint16)
            for i in range(0,npix_y):
                for j in range(0,npix_x):
                    value = ResultW_rho[j,i]
                    image_rho[i+ (0.5 * npix_y),j] = value
                


        h5file.close() 
        if edgeon:
            # Edge on image 
            npix_y /= 2 

            exec "ResultW_rho = h5file['%s_edgeOn']" % (array_name, )
            
            ResultW_rho = ResultW_rho - tf_min_rho			
            ResultW_rho = ResultW_rho / (tf_max_rho - tf_min_rho)
            ResultW_rho[ResultW_rho < 0.0] = 0.0
            ResultW_rho[ResultW_rho > 1.0] = 1.0
            ResultW_rho = ResultW_rho*255.0

            print 'Image range (8bit): ',min(np.ravel(ResultW_rho)),max(np.ravel(ResultW_rho))

            ResultW_rho = ResultW_rho.astype(np.uint16) 
                    
            for i in range(0,npix_y):
                for j in range(0,npix_x):
                    value = ResultW_rho[j,i]
                    image_rho[i,j] = value
                    


    # Now read in temperature image array 
    tf_min_T = min_temp
    tf_max_T = max_temp

    print 'tf_min_T = ',tf_min_T
    print 'tf_max_T = ',tf_max_T

    h5filename = "gasTemp_proj_%3d_%.2fkpc.hdf5" % (isnap, image_length)
    h5name=h5filename
    with h5py.File(data_dir_T + h5name, "r") as h5file:
        #npix_y *= 2 
        ResultQ_T = np.array(h5file['ResultQ_faceOn'])
        
        print 'Image range (temp): ',min(np.ravel(ResultQ_T)),max(np.ravel(ResultQ_T))
        ResultQ_T = ResultQ_T - tf_min_T 
        ResultQ_T = ResultQ_T / (tf_max_T - tf_min_T)
        ResultQ_T[ResultQ_T < 0.0] = 0.0
        ResultQ_T[ResultQ_T > 1.0] = 1.0
        ResultQ_T = ResultQ_T*255.0

        print 'Image range (8bit): ',min(np.ravel(ResultQ_T)),max(np.ravel(ResultQ_T))

        ResultQ_T = ResultQ_T.astype(np.uint16)    

        if not edgeon:
            image_T = np.ndarray((npix_y,npix_x),dtype=np.uint16)
            for i in range(0,npix_y):
                for j in range(0,npix_x):
                    value = ResultQ_T[j,i]
                    image_T[i,j] = value 
        else:
            image_T = np.ndarray((1.5*npix_y,npix_x),dtype=np.uint16)
            for i in range(0,npix_y):
                for j in range(0,npix_x):
                    value = ResultQ_T[j,i]
                    image_rho[i+ (0.5 * npix_y),j] = value
        # Edge on image 
        if edgeon:
            npix_y /= 2 
            ResultQ_T = h5file['ResultQ_edgeOn']
            
            ResultQ_T = ResultQ_T - tf_min_T 
            ResultQ_T = ResultQ_T / (tf_max_T - tf_min_T) 
            ResultQ_T[ResultQ_T < 0.0] = 0.0
            ResultQ_T[ResultQ_T > 1.0] = 1.0
            ResultQ_T = ResultQ_T*255.0

            print 'Image range (8bit): ',min(np.ravel(ResultQ_T)),max(np.ravel(ResultQ_T))

            ResultQ_T = ResultQ_T.astype(np.uint16) 
                    
            for i in range(0,npix_y):
                for j in range(0,npix_x):
                    value = ResultQ_T[j,i]
                    image_T[i,j] = value 


    # Now take the rho and T images, and combine them 
    # to produce the final image array. 
    final_image = mcm.produce_viridis_hsv_image(image_T, image_rho) 
            
    """
    fig = pylab.figure(figsize = (npix_x / 600.0, npix_y / 600.0), dpi=600, frameon=False)
    ax = pylab.Axes(fig, [0,0,1,1])
    ax.set_axis_off()
    fig = pylab.gcf()
    fig.add_axes(ax)
    """

    # Length of bar in code units
    scale_line_length *= HubbleParam
    # Convert to pixels
    length_per_pixel = (Xmax - Xmin) / npix_x
    scale_line_length_px = int(scale_line_length / length_per_pixel)

    # Position in terms of image array indices
    scale_line_x_start = int(0.05 * npix_x)
    scale_line_x_end = min(scale_line_x_start + scale_line_length_px,npix_x)
    scale_line_y = int(0.02 * npix_y)

    # Go through pixels for scale bar, setting them to white
    for x_index in xrange(scale_line_x_start, scale_line_x_end):
        final_image[scale_line_y, x_index, 0] = 1
        final_image[scale_line_y, x_index, 1] = 1
        final_image[scale_line_y, x_index, 2] = 1
        final_image[scale_line_y + 1, x_index, 0] = 1
        final_image[scale_line_y + 1, x_index, 1] = 1
        final_image[scale_line_y + 1, x_index, 2] = 1
        final_image[scale_line_y + 2, x_index, 0] = 1
        final_image[scale_line_y + 2, x_index, 1] = 1
        final_image[scale_line_y + 2, x_index, 2] = 1
        final_image[scale_line_y + 3, x_index, 0] = 1
        final_image[scale_line_y + 3, x_index, 1] = 1
        final_image[scale_line_y + 3, x_index, 2] = 1
        final_image[scale_line_y + 4, x_index, 0] = 1
        final_image[scale_line_y + 4, x_index, 1] = 1
        final_image[scale_line_y + 4, x_index, 2] = 1
        final_image[scale_line_y + 5, x_index, 0] = 1
        final_image[scale_line_y + 5, x_index, 1] = 1
        final_image[scale_line_y + 5, x_index, 2] = 1

    #figure_label2 = r"$\rm{UVBthin}$"
    imgplot = ax.imshow(final_image, 
        extent = (Xmin,Xmax,Ymin,Ymax),origin = 'lower', aspect = 'auto')
    figure_label = r"$%03d \, \rm{Myr}$" % (round_to_nearest_integer(time_Myr), )
    figure_label = r"$%.2f \, \rm{Myr}$" % (time_Myr)
    label = pylab.text(0.70, 0.92, figure_label, fontsize = 8, transform = ax.transAxes)
    label.set_color('white')
    label2 = pylab.text(scale_label_position, 0.03, scale_label_text, fontweight = 'bold', transform = ax.transAxes)
    label2.set_color('white')
    label2.set_fontsize(6)
    #label3 = pylab.text(0.10, 0.92, figure_label2, fontsize = 8, transform = ax.transAxes)
    #label3.set_color('white')

    # colour bar
    use_colorbar=0
    if use_colorbar:

        colour_map = matplotlib.colors.LinearSegmentedColormap.from_list("PaulT_rainbow", cols)
        colour_map.set_under('k')

        ax_c = pylab.axes([0.15, 0.415, 0.7, 0.015])
        norm = matplotlib.colors.Normalize(vmin = tf_min_T, vmax = tf_max_T)
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
            t.set_color('w')
            t.set_y(t.get_position()[1] + 0.5)
        cbar.outline.set_linewidth(0.4)
        cbar.set_label(my_cbar_label, color = 'w', fontsize=6, fontweight='bold', labelpad = 0.5)

        print cbar.get_clim()

    """
    extent = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    exec "pylab.savefig('%s/frame%04d.png', dpi = 600, bbox_inches=extent)" % (output_dir, isnap)
    pylab.close() 
    """

    return ax

def main(): 
    isnap_low = int(sys.argv[1]) 
    isnap_hi = int(sys.argv[2]) 
    sim_name = sys.argv[3]   # of the form 'data_dir' 

    for isnap in xrange(isnap_low, isnap_hi):
        plot_image_grid(isnap)

    return 

if __name__=='__main__':
    main() 
