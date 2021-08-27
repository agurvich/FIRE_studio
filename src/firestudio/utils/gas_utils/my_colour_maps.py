from matplotlib.colors import rgb_to_hsv, hsv_to_rgb 
import numpy as np 

import matplotlib.pyplot as plt
try:

    import palettable
except:
    print("palettable colormaps are not installed")

try:
    from pfh_colormaps import load_my_custom_color_tables
    load_my_custom_color_tables()
except:
    print("don't have phil's colormaps")

def get_cmap(cmap_name):
    try:
        cmap = plt.get_cmap(cmap_name)
    except:
        ## perhaps i was passed a palettable cmap path
        cmap_name = cmap_name.split(".")
        cmap = palettable
        for name in cmap_name:
            cmap = getattr(cmap,name)
        cmap = cmap.mpl_colormap
    return cmap

def produce_colmap(cmap_name):
    cmap = get_cmap(cmap_name)
    ## discretize the colormap into 256 parts...
    return [list(cmap(i/255.)[:3]) for i in range(0,256)]

def produce_cmap_hsv_image(image_1, image_2,cmap='viridis'): 
    # image_1 and image_2 are arrays of pixels 
    # with integer values in the range 0 to 255. 
    # These will be mapped onto hue and brightness, 
    # respectively, with saturation fixed at 1. 
    cols = produce_colmap(cmap)
    cols = np.array(cols) 

    cols_2d = np.ones((len(cols), 1, 3)) 
    cols_2d[:, 0, 0] = cols[:, 0] 
    cols_2d[:, 0, 1] = cols[:, 1] 
    cols_2d[:, 0, 2] = cols[:, 2] 
    cols_hsv = rgb_to_hsv(cols_2d) 
    hue_viridis = cols_hsv[:, 0, 0] 

    npix_x = len(image_1) 
    npix_y = len(image_1[0]) 
    if image_2 is not None:
        output_image_hsv = np.zeros((npix_x, npix_y, 3)) 
        for i in range(npix_x): 
            for j in range(npix_y): 
                output_image_hsv[i, j, 0] = hue_viridis[image_1[i, j]] 
                output_image_hsv[i, j, 1] = 1.0 
                output_image_hsv[i, j, 2] = float(image_2[i, j]) / 255.0 

        output_image_rgb = hsv_to_rgb(output_image_hsv) 
    else:
        output_image_rgb = np.zeros((npix_x, npix_y, 3)) 
        for i in range(npix_x): 
            for j in range(npix_y): 
                output_image_rgb[i,j]=cols[image_1[i,j]]
                
    return output_image_rgb 

