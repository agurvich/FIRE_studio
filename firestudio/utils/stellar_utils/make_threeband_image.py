import numpy as np
import math
import firestudio.utils.stellar_utils.contour_makepic as cmakepic
import firestudio.utils.stellar_utils.colors as viscolors
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

## some of my favorite regular cmaps below
def pick_custom_cmap(i):
    cm='hot'
    if(i==0): cm='heat_purple'
    if(i==1): cm='heat_green'
    if(i==2): cm='heat_blue'
    if(i==3): cm='heat_yellow'
    if(i==4): cm='heat_red'
    if(i==5): cm='heat_orange'
    if(i==6): cm='heat_redyellow'
    if(i==7): cm='pink' # whites out fairly quickly
    if(i==8): cm='bone' # pretty nice for white-ish (X-ray like colors)
    if(i==9): cm='copper' # orange-y
    if(i==10): cm='gray' # basic greyscale
    
    if(i==11): cm='spring'
    if(i==12): cm='summer'
    if(i==13): cm='winter'
    if(i==14): cm='autumn'
    if(i==15): cm='gist_earth'
    if(i==16): cm='Blues_r'
    if(i==17): cm='Greens_r'
    if(i==18): cm='Oranges_r'
    if(i==19): cm='Purples_r'
    if(i==20): cm='RdPu_r'
    if(i==21): cm='Reds_r'
    

    if(i==0): cm='heat_orange'
    if(i==0): cm='heat_redyellow'
    if(i==1): cm='heat_green'
    if(i==2): cm='heat_purple'
    return cm;


def layer_band_images(ims, maps):
    matplotlib.pyplot.imshow(0.*ims[:,:,0]); ## zero out background
    
    viscolors.load_my_custom_color_tables();
    nx=ims[:,0,0].size; ny=ims[0,:,0].size;
    im_new = np.zeros((nx,ny,3))
    map_cum = np.zeros((nx,ny))
    
    for i in range(ims.shape[2]):
        im = ims[:,:,i]
        map = maps[:,:,i]

        #im_0=im/np.max(im); # if want more saturated images
        ####alpha_im=maps[:,:,i]/map_sum; ## deprected
        cm=pick_custom_cmap(i);
        my_cmap=matplotlib.cm.get_cmap(cm); ## load cmap
        #rgb_im=my_cmap(im_0); ## get rgba values of image as mapped by this cmap
        
        rgb_im = my_cmap(im)[:,:,0:3]
        for j in [0,1,2]:
            im_new[:,:,j] = (map_cum*im_new[:,:,j] + map*rgb_im[:,:,j]) / (map_cum+map)
        map_cum += map
        
        #rgb_im[:,:,3]=alpha_im; ## replace alpha channel for this image
        #rgb_im[:,:,3]=0.*alpha_im+1.;
        #matplotlib.pyplot.imshow(rgb_im); ## plot it, with appropriate (new) alpha channel

    return im_new
    

def make_threeband_image_process_bandmaps(r,g,b, \
    dont_make_image=0, maxden=0, dynrange=0, pixels=720, \
    color_scheme_nasa=1, color_scheme_sdss=0 , \
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
    viscolors.load_my_custom_color_tables();
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



def make_threeband_image( x, y, lums, hsml=0, xrange=0, yrange=0, \
    dont_make_image=0, maxden=0, dynrange=0, pixels=720, \
    color_scheme_nasa=1, color_scheme_sdss=0 ):

    if (len(lums[0,:]) != len(x)) or (len(lums[:,0]) != 3):
        print(' expect error: lums must be an array of (3,n) ')
	
    ## set x/y ranges and clip the particle distribution 
    if checklen(xrange)<=1:
        x0=np.median(x); xm=np.median(np.fabs(x-x0)); xrange=[-xm+x0,xm+x0];
    if checklen(yrange)<=1:
        y0=np.median(y); ym=np.median(np.fabs(y-y0)); xrange=[-ym+y0,ym+y0];
    h=0.*x+1.; 
    if (checklen(hsml)==checklen(x)): h=hsml
    wt0=lums[0,:]; wt1=lums[1,:]; wt2=lums[2,:]

    tolfac=0.1 ; xr=xrange ; yr=yrange
    ok=ok_scan(x) & ok_scan(y) & ok_scan(h,pos=1) & ok_scan(wt0,pos=1,xmax=1.0e40) & \
        ok_scan(wt1,pos=1,xmax=1.0e40) & ok_scan(wt2,pos=1,xmax=1.0e40) & \
        (x >= xr[0]-(xr[1]-xr[0])*tolfac) & (x <= xr[1]+(xr[1]-xr[0])*tolfac) & \
        (y >= yr[0]-(yr[1]-yr[0])*tolfac) & (y <= yr[1]+(yr[1]-yr[0])*tolfac) 
    n_ok=checklen(x[ok]);
    x_c=0.5*(np.max(xr)+np.min(xr)); xlen=0.5*(np.max(xr)-np.min(xr)); xx=x-x_c;
    y_c=0.5*(np.max(yr)+np.min(yr)); ylen=0.5*(np.max(yr)-np.min(yr)); yy=y-y_c;
    zz=0.*xx; aspect_ratio=ylen/xlen; 
    
    ## make the image maps in each band
    u_band_map, r_band_map, k_band_map, dummy_pic = \
        cmakepic.contour_makepic( xx[ok], yy[ok], zz[ok], h[ok], \
        wt0[ok], weight2=wt1[ok], weight3=wt2[ok], \
        xlen=xlen, set_aspect_ratio=aspect_ratio, pixels=pixels );

    r=k_band_map #*4.9
    g=r_band_map #*5.7
    b=u_band_map #*7.8
    image24,cmap_m = make_threeband_image_process_bandmaps(r,g,b, \
        dont_make_image=dont_make_image, maxden=maxden, dynrange=dynrange, pixels=pixels, \
        color_scheme_nasa=color_scheme_nasa, color_scheme_sdss=color_scheme_sdss );
        
    ## insert the actual plotting here

    return image24,cmap_m; ## return both processed image and massmap
