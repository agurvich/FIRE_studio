## very simply projection plotting for SPH data :: 
##   note that clipping should be done *before* sending here, to ensure good behavior

import numpy as np
import ctypes
import colorsys
import firestudio.utils.stellar_utils.utilities as util
import math
import firestudio.utils.stellar_utils.colors as colors

def fcor(x):
    return np.array(x,dtype='f',ndmin=1)
def vfloat(x):
    return x.ctypes.data_as(ctypes.POINTER(ctypes.c_float));
def cfloat(x):
    return ctypes.c_float(x);
def checklen(x):
    return len(np.array(x,ndmin=1));

def contour_makepic( x, y, z, hsml, weight, 
        weight2=0, weight3=0, 
        xlen = 1, 
        pixels = 720, set_aspect_ratio = 1.0, 
        set_maxden = 1.0e-1, ## (gadget units, 10^10 msun/kpc^2 = 10^4 msun/pc^2)
        set_dynrng = 1.0e4, 
        invert_colorscale = 0, 
        set_percent_maxden = 0, set_percent_minden = 0 ):
        
    ## set some basic initial values
    xypixels=pixels; xpixels=xypixels; ypixels = np.around(float(xpixels)*set_aspect_ratio).astype(int);
    ylen=xlen*set_aspect_ratio;
    ma = set_maxden; mi = ma / set_dynrng;
    xmin=-xlen; xmax=xlen; ymin=-ylen; ymax=ylen;

    ## load the routine we need
    exec_call=util.return_python_routines_cdir()+'/SmoothedProjPFH/allnsmooth.so'
    smooth_routine=ctypes.cdll[exec_call];
    ## make sure values to be passed are in the right format
    N=checklen(x); x=fcor(x); y=fcor(y); M1=fcor(weight); M2=fcor(weight2); M3=fcor(weight3); H=fcor(hsml)
    xpixels=np.int(xpixels); ypixels=np.int(ypixels)
    ## check for whether the optional extra weights are set
    NM=1; 
    if(checklen(M2)==checklen(M1)): 
        NM=2; 
        if(checklen(M3)==checklen(M1)):
            NM=3;
        else:
            M3=np.copy(M1);
    else:
        M2=np.copy(M1);
        M3=np.copy(M1);
    ## initialize the output vector to recieve the results
    XYpix=xpixels*ypixels; MAP=ctypes.c_float*XYpix; MAP1=MAP(); MAP2=MAP(); MAP3=MAP();
    ## main call to the imaging routine
    smooth_routine.project_and_smooth( \
        ctypes.c_int(N), \
        vfloat(x), vfloat(y), vfloat(H), \
        ctypes.c_int(NM), \
        vfloat(M1), vfloat(M2), vfloat(M3), \
        cfloat(xmin), cfloat(xmax), cfloat(ymin), cfloat(ymax), \
        ctypes.c_int(xpixels), ctypes.c_int(ypixels), \
        ctypes.byref(MAP1), ctypes.byref(MAP2), ctypes.byref(MAP3) );
    ## now put the output arrays into a useful format 
    MassMap1=np.ctypeslib.as_array(MAP1).reshape([xpixels,ypixels]);
    MassMap2=np.ctypeslib.as_array(MAP2).reshape([xpixels,ypixels]);
    MassMap3=np.ctypeslib.as_array(MAP3).reshape([xpixels,ypixels]);
 
    # set boundaries and do some clipping
    MassMap = np.copy(MassMap1);
    print("MassMap : max: ", np.max(MassMap), "   min: ", np.min(MassMap))
    if (set_percent_maxden !=0) or (set_percent_minden !=0):
        print('percent max/min = ',set_percent_maxden,set_percent_minden)
        Msort=np.sort(MassMap);
        if (set_percent_maxden != 0): ma=Msort[set_percent_maxden*float(checklen(MassMap)-1)];
        mi=ma/set_dynrng;
        if (set_percent_minden != 0): mi=Msort[set_percent_minden*float(checklen(MassMap)-1)];
        if (set_percent_maxden == 0): ma=mi*set_dynrng;
        ok=(Msort > 0.) & (np.isnan(Msort)==False)
        if (mi <= 0) or (np.isnan(mi)): mi=np.min(Msort[ok]);
        if (ma <= 0) or (np.isnan(ma)): ma=np.max(Msort[ok]);
    print("Clipping at   ma= ", ma, " mi= ", mi)
    MassMap[MassMap < mi]=mi; MassMap[MassMap > ma]=ma;

    # now set colors
    cols=255. # number of colors
    Pic = (np.log(MassMap/mi)/np.log(ma/mi)) * (cols-3.) + 2.
    Pic[Pic > 255.]=255.; Pic[Pic < 1]=1.;
    backgrd = np.where((Pic<=2) | (np.isnan(Pic)))
    if (invert_colorscale==0): 
        Pic=256-Pic;
        Pic[backgrd] = 1; # white
    else:
        Pic[backgrd] = 0; # black
    
    if (NM>1):
        return MassMap1,MassMap2,MassMap3, Pic
    else:
        return MassMap, Pic




## best to send to a postscript file, but can do it as a part of that, or 
##   in any fashion desired. -- need to give it an xrange, otherwise fine!
def simple_makepic( x, y, 
    xrange=[-1.,1.], yrange=0, weights=0, hsml=0, 
    dont_make_plots=0, color_temperature=0, temp_weights=0, 
    set_temp_max=0, set_temp_min=0, 
    pixels = 720, invert_colorscale = 0, 
    set_maxden = 1.0e-1, set_dynrng = 1.0e4, 
    set_percent_maxden = 0, set_percent_minden = 0 ):
    
    
    xrange=np.array(xrange);
    if (checklen(yrange) <= 1): yrange=1.0*xrange;
    yrange=1.0*xrange;
    xmin=xrange[0]; xmax=xrange[1]; ymin=yrange[0]; ymax=yrange[1];
    xlen=0.5*(xmax-xmin); ylen=0.5*(ymax-ymin)
    x_c=0.5*(xmin+xmax); y_c=0.5*(ymin+ymax); xx=x-x_c; yy=y-y_c; zz=0.*xx;

    if (checklen(weights) <= 1): weights=0.*x+1.; weights /= math.fsum(weights);
    if (checklen(hsml) <= 1): hsml=0.*x+xlen/100.
    if ((color_temperature==0) or (checklen(temp_weights)<=1)): temp_weights = 0.*x+1.;
    temp_weights *= weights; ## include mass weighting in color averaging ## 
    aspect_ratio = ylen/xlen; 

    MassMap_1,MassMap_2,MassMap_3,Pic_1 = \
        contour_makepic( xx, yy, zz, hsml, weights, weight2=temp_weights, set_aspect_ratio=aspect_ratio, \
        xlen=xlen,set_maxden=set_maxden,set_dynrng=set_dynrng,pixels=pixels, \
        invert_colorscale=invert_colorscale, \
        set_percent_maxden=set_percent_maxden,set_percent_minden=set_percent_minden);
    bad = (MassMap_1 <= 0.); 
    M1Min=np.min(MassMap_1[MassMap_1 > 0.]); MassMap_1[bad]=M1Min;

    if (color_temperature==1):
        tempMap = MassMap_2/MassMap_1
        M2Min=np.min(MassMap_2[MassMap_2 > 0.]); 
        tempMap[bad] = (MassMap_2[bad] + M2Min) / M1Min;
        print('MinMax: ',np.min(MassMap_1),np.max(MassMap_1),np.min(MassMap_2),np.max(MassMap_2),np.min(tempMap),np.max(tempMap))
        Pic2_RGB = colors.temperature_map_color_index(Pic_1,tempMap, \
            set_temp_max=set_temp_max, set_temp_min=set_temp_min, huem100=0, invertreverse=1);

    if (dont_make_plots==1):
        if (color_temperature==1): return MassMap_1,Pic_1,MassMap_2,Pic2_RGB;
        return MassMap_1,Pic_1;


    if (color_temperature==1): return MassMap_1,Pic_1,MassMap_2,Pic2_RGB;
    return MassMap_1,Pic_1;


def test():
    N0=512
    x=np.random.rand(N0*N0)
    y=np.random.rand(N0*N0)
    m=0.*x+1./float(N0*N0)
    h=0.*x+0.01
    a,b,c,d=simple_makepic(x,y,xrange=[0.,1.],weights=m,hsml=h,set_maxden=1.0e-1,COLOR_TEMPERATURE=1)
    
    return a,b
