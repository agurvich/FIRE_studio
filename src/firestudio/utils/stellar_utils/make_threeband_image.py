import numpy as np
import matplotlib

def checklen(x): return np.array(x).shape[0]

def int_round(x): return np.int(np.round(x))

def ok_scan(arr,xmax=1.0e10,pos=0):
    if (pos==0): return np.logical_and(
        np.isfinite(arr),
        np.fabs(arr)<=xmax)
    if (pos==1): return np.logical_and(
        np.logical_and(
            np.isfinite(arr),
            np.fabs(arr)<=xmax),
        arr > 0.)

def clip_256(x,max=255,min=2):
    x = x*(max-min) + min
    x[x >= max]=max
    x[x <= min]=min
    x[np.isnan(x)]=min
    x /= 256.
    return x

def single_vec_sorted(x,reverse=False):
    return sorted(np.reshape(x,x.size),reverse=reverse)

def make_threeband_image_process_bandmaps(
    r,g,b,
    maxden,
    dynrange,
    color_scheme_nasa=1,
    color_scheme_sdss=0 , 
    filterset=None,
    QUIET=False):

    ## handle default arguments
    if filterset is None: filterset = ['r','g','b']

    ## now clip the maps and determine saturation levels
    cmap_m = np.zeros((checklen(r[:,0]),checklen(r[0,:]),3),dtype='f')
    cmap_m[:,:,0] = r
    cmap_m[:,:,1] = g
    cmap_m[:,:,2] = b

    ## now do the color processing on the maps
    maxnorm = maxden
    minnorm = (maxden/dynrange)
    if not QUIET: print(
        'maxden == ',maxnorm,
        'dynrange == ',dynrange,
        'minden == ',minnorm)

    rgb_avg = (r+g+b)/3.

    if (color_scheme_sdss==1):

        f_i = np.log10(rgb_avg/minnorm) / np.log10(maxnorm/minnorm)
        f_i[rgb_avg>=maxnorm] = 1.
        f_i[rgb_avg<=minnorm] = 0.

        q,alpha = 9.,0.3
        f_i = np.arcsinh( alpha * q * (rgb_avg/minnorm) ) / q
        wt = f_i/rgb_avg
        r *= wt
        g *= wt
        b *= wt

    if (color_scheme_nasa==1):
        r = np.log(r/minnorm) / np.log(maxnorm/minnorm)
        g = np.log(g/minnorm) / np.log(maxnorm/minnorm)
        b = np.log(b/minnorm) / np.log(maxnorm/minnorm)

    ## rescale to saturation limit
    bad = (rgb_avg<=0.)

    maxrgb = 0.
    if (checklen(rgb_avg[bad])>0):
        r[bad]=0.
        g[bad]=0.
        b[bad]=0.

    f_saturated=0.0004  ## fraction of pixels that should be saturated
    f_saturated=0.0001  ## fraction of pixels that should be saturated
    x0 = int_round( f_saturated * (np.float(checklen(r)) - 1.) )

    for rgb_v in [r,g,b]: 
        rgbm = single_vec_sorted(rgb_v,reverse=True)
        if(rgbm[x0]>maxrgb): maxrgb=rgbm[x0]

    if (maxrgb > 1.): 
        r/=maxrgb
        g/=maxrgb
        b/=maxrgb

    ## rescale to 256-colors to clip the extremes (rescales back to 0-1):
    r,g,b = clip_256(r),clip_256(g),clip_256(b)

    ## initialize the output image
    image24 = np.zeros(
        (checklen(r[:,0]),checklen(r[0,:]),3),
        dtype='f')

    ## fill the output image with rgb values
    image24[:,:,0] = r
    image24[:,:,1] = g
    image24[:,:,2] = b
    
    ## ok have r, g, b -- really just three re-scaled maps. no reason they 
    ##   have to map to r, g, b colors: use the filter set given to map them ::
    image24_new = 0.*image24
    for i in range(3):
        im = image24[:,:,i]
        if filterset[i] == 'r': image24_new[:,:,0] = im
        elif filterset[i] == 'g': image24_new[:,:,1] = im
        elif filterset[i] == 'b': image24_new[:,:,2] = im
        else:
            my_cmap = matplotlib.cm.get_cmap(filterset[i])
            rgb_im = my_cmap(im)
            image24_new += rgb_im[:,:,0:3] ## dropping the alpha channel here!

    image24 = image24_new

    return image24, cmap_m ## return both processed image and massmap

def layer_band_images(ims, maps):
    
    nx=ims[:,0,0].size; ny=ims[0,:,0].size;
    im_new = np.zeros((nx,ny,3))
    map_cum = np.zeros((nx,ny))

    cmaps = ['Blues',
        'heat_green',
        'heat_red']
    
    for i in range(ims.shape[2]):
        im = ims[:,:,i]
        map = maps[:,:,i]

        #im_0=im/np.max(im); # if want more saturated images
        ####alpha_im=maps[:,:,i]/map_sum; ## deprected
        my_cmap=matplotlib.cm.get_cmap(cmaps[i]); ## load cmap
        #rgb_im=my_cmap(im_0); ## get rgba values of image as mapped by this cmap
        
        rgb_im = my_cmap(im)[:,:,0:3]
        for j in [0,1,2]:
            im_new[:,:,j] = (map_cum*im_new[:,:,j] + map*rgb_im[:,:,j]) / (map_cum+map)
        map_cum += map
        
        #rgb_im[:,:,3]=alpha_im; ## replace alpha channel for this image
        #rgb_im[:,:,3]=0.*alpha_im+1.;
        #matplotlib.pyplot.imshow(rgb_im); ## plot it, with appropriate (new) alpha channel

    return im_new
