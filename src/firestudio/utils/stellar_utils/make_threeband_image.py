import numpy as np
import matplotlib

def checklen(x):
    return np.array(x).size;

def int_round(x):
    return np.int(np.round(x));

def ok_scan(input,xmax=1.0e10,pos=0):
    if (pos==1):
        return (np.isnan(input)==False) & (abs(input)<=xmax) & (input > 0.);
    else:
        return (np.isnan(input)==False) & (abs(input)<=xmax);

def clip_256(x,max=255,min=2):
    x = x*(max-min) + min;
    x[x >= max]=max;
    x[x <= min]=min;
    x[np.isnan(x)]=min;
    x /= 256.;
    return x;

def single_vec_sorted(x,reverse=False):
    return sorted(np.reshape(x,x.size),reverse=reverse); 

def make_threeband_image_process_bandmaps(
    r,g,b,
    dont_make_image=0, maxden=0, dynrange=0, pixels=720, 
    color_scheme_nasa=1, color_scheme_sdss=0 , 
    filterset = ['r','g','b'] ):

    ## now clip the maps and determine saturation levels
    cmap_m=np.zeros((checklen(r[:,0]),checklen(r[0,:]),3),dtype='f');
    cmap_m[:,:,0]=r; cmap_m[:,:,1]=g; cmap_m[:,:,2]=b; 
    if (dont_make_image==1): return cmap_m;

    if (maxden<=0):
        f_saturated=0.005 ## fraction of pixels that should be saturated 
        x0=int_round( f_saturated * (np.float(checklen(r)) - 1.) );
        for rgb_v in [r,g,b]: 
            rgbm=single_vec_sorted(rgb_v,reverse=True);
            if(rgbm[x0]>maxden): maxden=rgbm[x0]

    if (dynrange<=0):
        f_zeroed=0.1 ## fraction of pixels that should be black 		
        x0=int_round( f_zeroed * (np.float(checklen(r)) - 1.) ); minden=np.max(r); 
        rgbm=single_vec_sorted(r+g+b,reverse=False); 
        if(rgbm[x0]<minden): minden=rgbm[x0];
        #for rgb_v in [r,g,b]: 
        #    rgbm=single_vec_sorted(rgb_v,reverse=False); 
        #    if(rgbm[x0]<minden): minden=rgbm[x0];
        if (minden<=0):
            minden = np.min(np.concatenate((r[r>0.],g[g>0.],b[b>0.])));
        dynrange = maxden/minden;

    ## now do the color processing on the maps
    maxnorm=maxden; minnorm=(maxden/dynrange);
    print('maxden == ',maxnorm,' dynrange == ',dynrange,' minden == ',minnorm)

    i = (r+g+b)/3.
    f_i = np.log10(i/minnorm) / np.log10(maxnorm/minnorm);
    f_i[i>=maxnorm]=1.; f_i[i<=minnorm]=0.

    if (color_scheme_sdss==1):
        q=9.; alpha=0.3;
        f_i = np.arcsinh( alpha * q * (i/minnorm) ) / q; 
        wt=f_i/i; r*=wt; g*=wt; b*=wt;
    if (color_scheme_nasa==1):
        r = np.log(r/minnorm) / np.log(maxnorm/minnorm);
        g = np.log(g/minnorm) / np.log(maxnorm/minnorm);
        b = np.log(b/minnorm) / np.log(maxnorm/minnorm);

    ## rescale to saturation limit
    bad=(i<=0.); maxrgb=0.;
    if (checklen(i[bad])>0): r[bad]=0.; g[bad]=0.; b[bad]=0.;
    f_saturated=0.0004  ## fraction of pixels that should be saturated
    f_saturated=0.0001  ## fraction of pixels that should be saturated
    x0=int_round( f_saturated * (np.float(checklen(r)) - 1.) ); 
    for rgb_v in [r,g,b]: 
        rgbm=single_vec_sorted(rgb_v,reverse=True); 
        if(rgbm[x0]>maxrgb): maxrgb=rgbm[x0]
    if (maxrgb > 1.): r/=maxrgb; g/=maxrgb; b/=maxrgb;
    ## rescale to 256-colors to clip the extremes (rescales back to 0-1):
    max_c=255; min_c=2;
    r=clip_256(r,max=max_c,min=min_c);
    g=clip_256(g,max=max_c,min=min_c);
    b=clip_256(b,max=max_c,min=min_c);
    image24=np.zeros((checklen(r[:,0]),checklen(r[0,:]),3),dtype='f');
    image24[:,:,0]=r; image24[:,:,1]=g; image24[:,:,2]=b; 
    
    ## ok have r, g, b -- really just three re-scaled maps. no reason they 
    ##   have to map to r, g, b colors: use the filter set given to map them ::
    image24_new = 0.*image24
    for i in [0,1,2]:
        im=image24[:,:,i]
        if filterset[i]=='r': image24_new[:,:,0] = im
        if filterset[i]=='g': image24_new[:,:,1] = im
        if filterset[i]=='b': image24_new[:,:,2] = im
        if (filterset[i] != 'r') & (filterset[i] != 'g') & (filterset[i] != 'b'):
            my_cmap = matplotlib.cm.get_cmap(filterset[i])
            rgb_im = my_cmap(im)
            image24_new += rgb_im[:,:,0:3] ## dropping the alpha channel here!
    image24 = image24_new

    return image24, cmap_m; ## return both processed image and massmap