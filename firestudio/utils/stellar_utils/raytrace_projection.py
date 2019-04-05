import os
import numpy as np
import math
import ctypes
import firestudio.utils.stellar_utils.utilities as util
import scipy

def checklen(x):
    return len(np.array(x,ndmin=1));
def int_round(x):
    return np.int(np.round(x));
def ok_scan(input,xmax=1.0e30,pos=0):
    if (pos==0):
        return (np.isnan(input)==False) & (np.isfinite(input)) & (np.fabs(input)<=xmax);
    if (pos==1):
        return (np.isnan(input)==False) & (np.isfinite(input)) & (np.fabs(input)<=xmax) & (input > 0.);
    if (pos==2):
        return (np.isnan(input)==False) & (np.isfinite(input)) & (np.fabs(input)<=xmax) & (input >= 0.);
def fcor(x):
    return np.array(x,dtype='f',ndmin=1)
def vfloat(x):
    return x.ctypes.data_as(ctypes.POINTER(ctypes.c_float));


## 
## routine to use 'raytrace_projection_compute' to make mock gas images, 
##   with three color channels for different temperature ranges
## 
def gas_raytrace_temperature( TEMPERATURE_CUTS, \
        gas_x, gas_y, gas_z, gas_temperature, gas_mass, gas_hsml, \
        xlim=0, ylim=0, zlim=0, pixels=720, 
        KAPPA_UNITS=2.08854068444, kernel_width=0.05, use_log_t=1 , \
        isosurfaces=0 , add_temperature_weights=0 ):
    
    wtfn = gas_mass #* np.sqrt(gas_temperature/1.0e4) 
    if (add_temperature_weights==1): wtfn *= np.sqrt(1. + gas_temperature/1.0e4) 
    # weighting by sqrt(temp) makes the contributions more similar by temperature bins

    tcuts=TEMPERATURE_CUTS;
    tval=gas_temperature
    if(use_log_t==1):
        tcuts=np.log10(tcuts);
        tval=np.log10(tval);

    if np.array(kernel_width).size > 1:
        w = kernel_width
    else:
        w = kernel_width + np.zeros(3)
    ## continuous smoothing with gaussians for the temperature:
    if (isosurfaces==1):
        wt1 = np.exp(-(tval-tcuts[0])*(tval-tcuts[0])/(2.*w[0]*w[0]));
        wt2 = np.exp(-(tval-tcuts[1])*(tval-tcuts[1])/(2.*w[1]*w[1]));
        wt3 = np.exp(-(tval-tcuts[2])*(tval-tcuts[2])/(2.*w[2]*w[2]));
    else: ## isosurfaces==0, so do total in integral ranges set by temperature_cuts  
        wt1 = 0.5*(1.0-scipy.special.erf((tval-tcuts[0])/(np.sqrt(2.)*w[0])));
        wt3 = 0.5*(1.0-scipy.special.erf((tcuts[1]-tval)/(np.sqrt(2.)*w[1])));
        wt2 = 1.-wt1-wt3; wt2[wt2<0.]=0.;

    wt1*= wtfn; wt2*=wtfn; wt3*=wtfn;
    kappa = 200. * (1.+np.zeros((3)));
    kappa *= KAPPA_UNITS;
    print('KAPPA == ',kappa)

    return raytrace_projection_compute(gas_x,gas_y,gas_z,gas_hsml,gas_mass,\
        wt1,wt2,wt3,kappa[0],kappa[1],kappa[2],\
        xlim=xlim,ylim=ylim,zlim=zlim,pixels=pixels,TRIM_PARTICLES=1);




## 
## routine to use 'raytrace_projection_compute' to make mock stellar images, 
##   treating the starlight as sources and accounting for gas extinction via ray-tracing
## 
## important to set KAPPA_UNITS appropriately: code loads opacity (kappa) for 
##   the bands of interest in cgs (cm^2/g), must be converted to match units of input 
##   mass and size. the default it to assume gadget units (M=10^10 M_sun, l=kpc)
##
def stellar_raytrace( BAND_IDS, \
        stellar_x, stellar_y, stellar_z, \
        stellar_mass, stellar_age, stellar_metallicity, stellar_hsml, \
        gas_x, gas_y, gas_z, gas_mass, gas_metallicity, gas_hsml, \
        xlim=0, ylim=0, zlim=0, pixels=720, 
        KAPPA_UNITS=2.08854068444, \
        IMF_CHABRIER=1, IMF_SALPETER=0 , \
        ADD_BASE_METALLICITY=0.0, ADD_BASE_AGE=0.0 ):
        
    Nbands=len(np.array(BAND_IDS)); Nstars=len(np.array(stellar_mass)); Ngas=len(np.array(gas_mass));
    if (Nbands != 3): print("stellar_raytrace needs 3 bands, you gave"),Nbands; return -1,-1,-1,-1;
    ## check if stellar metallicity is a matrix
    if (len(stellar_metallicity.shape)>1): stellar_metallicity=stellar_metallicity[:,0];
    if (len(gas_metallicity.shape)>1): gas_metallicity=gas_metallicity[:,0];

    ## get opacities and luminosities at frequencies we need:
    stellar_metallicity[stellar_metallicity>0] += ADD_BASE_METALLICITY;
    gas_metallicity[gas_metallicity>0] += ADD_BASE_METALLICITY;
    stellar_age += ADD_BASE_AGE;
    kappa=np.zeros([Nbands]); lums=np.zeros([Nbands,Nstars]);
    for i_band in range(Nbands):
        nu_eff = util.colors_table(np.array([1.0]),np.array([1.0]), \
            BAND_ID=BAND_IDS[i_band],RETURN_NU_EFF=1);
        kappa[i_band] = util.opacity_per_solar_metallicity(nu_eff);
        l_m_ssp = util.colors_table( stellar_age, stellar_metallicity/0.02, \
            BAND_ID=BAND_IDS[i_band], CHABRIER_IMF=IMF_CHABRIER, SALPETER_IMF=IMF_SALPETER, CRUDE=1, \
            UNITS_SOLAR_IN_BAND=1); ## this is such that solar-type colors appear white
        l_m_ssp[l_m_ssp >= 300.] = 300. ## just to prevent crazy values here 
        l_m_ssp[l_m_ssp <= 0.] = 0. ## just to prevent crazy values here 
        lums[i_band,:] = stellar_mass * l_m_ssp
    gas_lum=np.zeros(Ngas); ## gas has no 'source term' for this calculation
    stellar_mass_attenuation = np.zeros(Nstars); ## stars have no 'attenuation term'
    gas_mass_metal = gas_mass * (gas_metallicity/0.02);
    kappa *= KAPPA_UNITS;
    
    ## combine the relevant arrays so it can all be fed into the ray-tracing
    x=np.concatenate([stellar_x,gas_x]); y=np.concatenate([stellar_y,gas_y]); z=np.concatenate([stellar_z,gas_z]);
    mass=np.concatenate([stellar_mass_attenuation,gas_mass_metal]);
    hsml=np.concatenate([stellar_hsml,gas_hsml]);
    wt1=np.concatenate([lums[0,:],gas_lum]); wt2=np.concatenate([lums[1,:],gas_lum]); wt3=np.concatenate([lums[2,:],gas_lum]);
    k1=kappa[0]; k2=kappa[1]; k3=kappa[2];
        
    return raytrace_projection_compute(x,y,z,hsml,mass,wt1,wt2,wt3,k1,k2,k3,\
        xlim=xlim,ylim=ylim,zlim=zlim,pixels=pixels,TRIM_PARTICLES=1);


##
##  Wrapper for raytrace_rgb, program which does a simply line-of-sight projection 
##    with multi-color source and self-extinction along the sightline: here called 
##    from python in its most general form: from the c-code itself:
##
##  int raytrace_rgb(
##    int N_xy, // number of input particles/positions
##    float *x, float *y, // positions (assumed already sorted in z)
##    float *hsml, // smoothing lengths for each
##    float *Mass, // total weight for 'extinction' part of calculation
##    float *wt1, float *wt2, float *wt3, // weights for 'luminosities'
##    float KAPPA1, float KAPPA2, float KAPPA3, // opacities for each channel
##    float Xmin, float Xmax, float Ymin, float Ymax, // boundaries of output grid
##    int Xpixels, int Ypixels, // dimensions of grid
##    float *OUT0, float *OUT1, float *OUT2, float*OUT3 ) // output vectors with final weights
##
def raytrace_projection_compute( x, y, z, hsml, mass, wt1, wt2, wt3, \
    kappa_1, kappa_2, kappa_3, xlim=0, ylim=0, zlim=0, pixels=720, \
    TRIM_PARTICLES=1 ):

    ## define bounaries
    if(checklen(xlim)<=1): xlim=[np.min(x),np.max(x)];
    if(checklen(ylim)<=1): ylim=[np.min(y),np.max(y)];
    if(checklen(zlim)<=1): zlim=[np.min(z),np.max(z)];
    xr=xlim; yr=ylim; zr=zlim;
    x00=0.5*(xr[1]+xr[0]); y00=0.5*(yr[1]+yr[0]); z00=0.5*(zr[1]+zr[0]); 
    tolfac = 1.0e10;
    if (TRIM_PARTICLES==1): tolfac = 0.05; 

    ## clip to particles inside those
    xlen=0.5*(xr[1]-xr[0]); ylen=0.5*(yr[1]-yr[0]); zlen=0.5*(zr[1]-zr[0]);
    x-=x00; y-=y00; z-=z00; dx=xlen*(1.+tolfac*2.); dy=ylen*(1.+tolfac*2.); dz=zlen*(1.+tolfac*2.);
    ok=ok_scan(x,xmax=dx) & ok_scan(y,xmax=dy) & ok_scan(z,xmax=dz) & \
        ok_scan(hsml,pos=1) & ok_scan(mass+wt1+wt2+wt3,pos=1)
        #& ok_scan(mass) & ok_scan(wt1) & ok_scan(wt2) & ok_scan(wt3);
    x=x[ok]; y=y[ok]; z=z[ok]; hsml=hsml[ok]; mass=mass[ok]; wt1=wt1[ok]; wt2=wt2[ok]; wt3=wt3[ok];
    N_p=checklen(x); xmin=-xlen; xmax=xlen; ymin=-ylen; ymax=ylen;
    if(N_p<=1): 
        print(' UH-OH: EXPECT ERROR NOW, there are no valid source/gas particles to send!'); return -1,-1,-1,-1;

    ## now sort these in z (this is critical!)
    s=np.argsort(z);
    x=x[s]; y=y[s]; z=z[s]; hsml=hsml[s]; mass=mass[s]; wt1=wt1[s]; wt2=wt2[s]; wt3=wt3[s];
    ## cast new copies to ensure the correct formatting when fed to the c-routine:
    x=fcor(x); y=fcor(y); z=fcor(z); hsml=fcor(hsml); mass=fcor(mass); wt1=fcor(wt1); wt2=fcor(wt2); wt3=fcor(wt3);

    ## load the routine we need
    curpath = os.path.realpath(__file__)
    curpath = curpath[:len("utils")+curpath.index("utils")] #split off this filename
    exec_call=os.path.join(curpath,'stellar_utils/c_libraries/RayTrace_RGB/raytrace_rgb.so')
    routine=ctypes.cdll[exec_call];
    
    ## cast the variables to store the results
    aspect_ratio=ylen/xlen; Xpixels=int_round(pixels); Ypixels=int_round(aspect_ratio*np.float(Xpixels));
    N_pixels=Xpixels*Ypixels; out_cast=ctypes.c_float*N_pixels; 
    out_0=out_cast(); out_1=out_cast(); out_2=out_cast(); out_3=out_cast();

    ## main call to the calculation routine
    routine.raytrace_rgb( ctypes.c_int(N_p), \
        vfloat(x), vfloat(y), vfloat(hsml), vfloat(mass), \
        vfloat(wt1), vfloat(wt2), vfloat(wt3), \
        ctypes.c_float(kappa_1), ctypes.c_float(kappa_2), ctypes.c_float(kappa_3), \
        ctypes.c_float(xmin), ctypes.c_float(xmax), ctypes.c_float(ymin), ctypes.c_float(ymax), \
        ctypes.c_int(Xpixels), ctypes.c_int(Ypixels), \
        ctypes.byref(out_0), ctypes.byref(out_1), ctypes.byref(out_2), ctypes.byref(out_3) );

    print(np.sum(out_1),np.sum(out_2),np.sum(out_3),'outputs')
    ## now put the output arrays into a useful format 
    out_0 = np.copy(np.ctypeslib.as_array(out_0));
    out_1 = np.copy(np.ctypeslib.as_array(out_1));
    out_2 = np.copy(np.ctypeslib.as_array(out_2));
    out_3 = np.copy(np.ctypeslib.as_array(out_3));
    out_0 = out_0.reshape([Xpixels,Ypixels]);
    out_1 = out_1.reshape([Xpixels,Ypixels]);
    out_2 = out_2.reshape([Xpixels,Ypixels]);
    out_3 = out_3.reshape([Xpixels,Ypixels]);

    return out_0, out_1, out_2, out_3;
