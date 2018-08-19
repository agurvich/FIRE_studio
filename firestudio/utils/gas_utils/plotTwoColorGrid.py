import matplotlib 
matplotlib.use('Agg') 
import pylab 
import numpy as np 
import firestudio.utils.gas_utils.my_colour_maps as mcm 
import sys 
import h5py
import os

def plot_image_grid(
    ax,
    frame_half_width,frame_depth,
    frame_center,
    projection_dir,snapnum,
    quantity_name='Temperature',
    this_setup_id = None,
    pixels=1200,
    theta=0,phi=0,psi=0,
    cmap='viridis',
    min_den=-1.0,max_den=1.2,
    min_quantity=2,max_quantity=7,
    h5prefix='',
    scale_bar = 1,
    figure_label=None,
    fontsize=None,
    **kwargs): 

    print("extra kwargs in plot_2color_image:",list(kwargs.keys()))

    ## Set parameters
    npix_x   = pixels  
    npix_y   = pixels  
    image_length =2*frame_half_width # kpc
    scale_label_position = 0.06 

    ## set scale bar length
    if image_length > 15 : 
        scale_line_length = 5
        scale_label_text = r"$\mathbf{5 \, \rm{kpc}}$"

    elif image_length > 1.5 : 
        scale_line_length = 1.0 
        scale_label_text = r"$\mathbf{1 \, \rm{kpc}}$"

    else:
        scale_line_length = .1
        scale_label_text = r"$\mathbf{100 \, \rm{pc}}$"

    print('min_den = ',min_den)
    print('max_den = ',max_den)

    ## uniquely identify this projection setup
    this_setup_id = (
	"npix%d_width%.2fkpc_depth%.2fkpc_x%.2f_y%.2f_z%.2f_theta%.2f_phi%.2f_psi%.2f"%(
	    npix_x, 2*frame_half_width,frame_depth,
	    frame_center[0],frame_center[1],frame_center[2],
	    theta,phi,psi)
	if this_setup_id is None else this_setup_id)

    ## read in the column density and quantity maps
    h5name=h5prefix+"proj_maps_%03d.hdf5" % snapnum
    with h5py.File(os.path.join(projection_dir,h5name), "r") as handle:
        this_group=handle[this_setup_id]
        columnDensityMap = np.array(this_group['columnDensityMap'])
        massWeightedQuantityMap = np.array(this_group['massWeighted%sMap'%quantity_name.title()])

    Xmin,Ymin = -frame_half_width+frame_center[:2]
    Xmax,Ymax = frame_half_width+frame_center[:2]
        
    print('Image range (rho): ',np.max(columnDensityMap),np.max(columnDensityMap))
    columnDensityMap = columnDensityMap - min_den 
    columnDensityMap = columnDensityMap / (max_den - min_den)
    
    ## clip anything outside the range
    columnDensityMap[columnDensityMap < 0.0] = 0.0
    columnDensityMap[columnDensityMap > 1.0] = 1.0
    columnDensityMap = columnDensityMap*255.0

    print('Image range (8bit): ',np.min(columnDensityMap),np.max(columnDensityMap))

    ## cast to integer to use as indices for cmap array
    columnDensityMap = columnDensityMap.astype(np.uint16) 
    image_rho = columnDensityMap.T

    print('min_%s = '%quantity_name,min_quantity)
    print('max_%s = '%quantity_name,max_quantity)
        
    print('Image range (%s): '%quantity_name,np.min(massWeightedQuantityMap),np.max(massWeightedQuantityMap))
    massWeightedQuantityMap = massWeightedQuantityMap - min_quantity 
    massWeightedQuantityMap = massWeightedQuantityMap / (max_quantity - min_quantity)
    
    ## clip anything outside the range
    massWeightedQuantityMap[massWeightedQuantityMap < 0.0] = 0.0
    massWeightedQuantityMap[massWeightedQuantityMap > 1.0] = 1.0
    massWeightedQuantityMap = massWeightedQuantityMap*255.0
    print('Image range (8bit): ',np.min(massWeightedQuantityMap),np.max(massWeightedQuantityMap))

    ## cast to integer to use as indices for cmap array 
    massWeightedQuantityMap = massWeightedQuantityMap.astype(np.uint16)    
    image_T = massWeightedQuantityMap.T
    
    ## Now take the rho and T images, and combine them 
    ##	to produce the final image array. 
    final_image = mcm.produce_cmap_hsv_image(image_T, image_rho,cmap=cmap) 

    ## fill the pixels of the the scale bar with white
    if scale_bar:
        # Convert to pixel space
        length_per_pixel = (Xmax - Xmin) / npix_x
        scale_line_length_px = int(scale_line_length / length_per_pixel)

        # Position in terms of image array indices
        scale_line_x_start = int(0.05 * npix_x)
        scale_line_x_end = min(scale_line_x_start + scale_line_length_px,npix_x)
        scale_line_y = int(0.02 * npix_y)

        # Go through pixels for scale bar, setting them to white
        for x_index in range(scale_line_x_start, scale_line_x_end):
            final_image[scale_line_y:scale_line_y+6, x_index,:3] = 1

    ## main imshow call
    imgplot = ax.imshow(
	final_image, 
        extent = (Xmin,Xmax,Ymin,Ymax),
	origin = 'lower', aspect = 'auto')

    ## handle any text additions
    fontsize=8 if fontsize is None else fontsize 
    if figure_label is not None:
    ## plot the  figure label in the top right corner
        label = pylab.text(
	    0.95, 0.92,
	    figure_label,
	    fontsize = fontsize,
	    transform = ax.transAxes,
	    ha='right')
        label.set_color('white')

    if scale_bar: 
	## plot the scale bar label
        label2 = pylab.text(
	    scale_label_position,0.03,
	    scale_label_text,
	    fontweight = 'bold',
	    fontsize=fontsize*0.75,
	    transform = ax.transAxes)
        label2.set_color('white')

    ## colour bar
    use_colorbar=0
    if use_colorbar:
        raise Exception('Unimplemented!')
        colour_map = matplotlib.colors.LinearSegmentedColormap.from_list("PaulT_rainbow", cols)
        colour_map.set_under('k')

        ax_c = pylab.axes([0.15, 0.415, 0.7, 0.015])
        norm = matplotlib.colors.Normalize(vmin = min_quantity, vmax = max_quantity)
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

        print(cbar.get_clim())

    """
    extent = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    exec "pylab.savefig('%s/frame%04d.png', dpi = 600, bbox_inches=extent)" % (output_dir, snapnum)
    pylab.close() 
    """

    return ax

"""
        if edgeon:
            # Edge on image 
            npix_y /= 2

            exec "columnDensityMap = np.array(h5file['%s_edgeOn'])" % (array_name, )
            
            columnDensityMap = columnDensityMap - min_den			
            columnDensityMap = columnDensityMap / (max_den - min_den)
            columnDensityMap[columnDensityMap < 0.0] = 0.0
            columnDensityMap[columnDensityMap > 1.0] = 1.0
            columnDensityMap = columnDensityMap*255.0

            print 'Image range (8bit): ',min(np.ravel(columnDensityMap)),max(np.ravel(columnDensityMap))

            columnDensityMap = columnDensityMap.astype(np.uint16) 
                    
            for i in range(0,npix_y):
                for j in range(0,npix_x):
                    value = columnDensityMap[j,i]
                    image_rho[i,j] = value

	"""



"""
        else:
            image_T = np.ndarray((3*npix_y/2,npix_x),dtype=np.uint16)
            for i in range(0,npix_y):
                new_i = i+ (npix_y/2)
                for j in range(0,npix_x):
                    value = massWeightedQuantityMap[j,i]
                    image_T[new_i,j] = value

        # Edge on image 
        if edgeon:
            npix_y /= 2 
            massWeightedQuantityMap = np.array(h5file['ResultQ_edgeOn'])
            
            massWeightedQuantityMap = massWeightedQuantityMap - min_quantity 
            massWeightedQuantityMap = massWeightedQuantityMap / (max_quantity - min_quantity) 
            massWeightedQuantityMap[massWeightedQuantityMap < 0.0] = 0.0
            massWeightedQuantityMap[massWeightedQuantityMap > 1.0] = 1.0
            massWeightedQuantityMap = massWeightedQuantityMap*255.0

            print 'Image range (8bit): ',min(np.ravel(massWeightedQuantityMap)),max(np.ravel(massWeightedQuantityMap))

            massWeightedQuantityMap = massWeightedQuantityMap.astype(np.uint16) 
                    
            for i in range(0,npix_y):
                for j in range(0,npix_x):
                    value = massWeightedQuantityMap[j,i]
                    image_T[i,j] = value 

    """
