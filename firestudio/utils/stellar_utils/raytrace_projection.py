import os
import numpy as np
import math
import ctypes

from firestudio.utils.stellar_utils.colors_sps.colors_table import colors_table
from firestudio.utils.stellar_utils.attenuation.cross_section import opacity_per_solar_metallicity

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
def read_band_lums_from_tables(
    BAND_IDS, 
    stellar_mass,stellar_age,stellar_metallicity,
    ## flag to return luminosity in each band requested without projecting
    nu_effs=None,
    lums=None,
    QUIET=False,
    IMF_CHABRIER=1,
    IMF_SALPETER=0):
        
    if nu_effs is None:
        nu_effs = [None,None,None]

    ## count particles we're using
    Nstars=len(np.array(stellar_mass))

    ## count how many bands we're attenuating
    Nbands=len(np.array(BAND_IDS))

    ## require that we attenuate 3 bands to combine since attenuation
    ##  routine is hardcoded to accept 3 weights
    if (Nbands != 3): 
        print("stellar_raytrace needs 3 bands, you gave",Nbands)
        return -1,-1,-1;

    ## check if stellar metallicity is a matrix
    ##  i.e. mass fraction of many species. If so,
    ##  take the total metallicity 
    if (len(stellar_metallicity.shape)>1): 
        stellar_metallicity=stellar_metallicity[:,0];

    ## get opacities and luminosities at frequencies we need:
    kappa=np.zeros([Nbands])
    if lums is None:
        lums=np.zeros([Nbands,Nstars])
    else:
        ## verify shape of lums is correct
        if lums.shape != (Nbands,Nstars):
            raise ValueError(
                "Shape (%d,%d) of lums does not match (3,%d)"%(
                lums.shape[0],lums.shape[1],Nstars))

    for i_band in range(Nbands):
        if nu_effs[i_band] is None:
            if not np.all(lums[i_band]==0):
                raise ValueError(
                    "Non-zero lums passed in axis %d"%i_band+
                    " without corresponding nu_eff")
            ## find the frequency associated with this band
            nu_effs[i_band] = colors_table(
                np.array([1.0]), ## dummy value
                np.array([1.0]), ## dummy value
                BAND_ID=BAND_IDS[i_band], ## band index
                RETURN_NU_EFF=1,
                QUIET=True) ## flag to return effective NU in this band


        ## calculate the kappa in this band using:
        ##  Thompson scattering + 
        ##  Pei (1992) + -- 304 < lambda[Angstroms] < 2e7
        ##  Morrison & McCammon (1983) -- 1.2 < lambda[Angstroms] < 413
        kappa[i_band] = opacity_per_solar_metallicity(
            nu_effs[i_band])

        these_lums = lums[i_band]
        ## if lums were not passed in for this band
        if np.all( these_lums == 0):
            ## lookup the luminosity/mass in this band
            ##  given stellar ages and metallicities
            these_lums[:] = colors_table(
                stellar_age, ## ages in Gyr
                stellar_metallicity/0.02,  ## metallicity in solar
                BAND_ID=BAND_IDS[i_band], ## band index
                CHABRIER_IMF=IMF_CHABRIER, ## imf flags
                SALPETER_IMF=IMF_SALPETER, ## imf flags
                CRUDE=1, ## map particles to nearest table entry rather than interpolate
                UNITS_SOLAR_IN_BAND=1, ## return ((L_star)_band / L_sun) / M_sun
                QUIET=QUIET
                ) 

        #these_lums[these_lums >= 300.] = 300. ## just to prevent crazy values here 
        #these_lums[these_lums <= 0.] = 0. ## just to prevent crazy values here 
        lums[i_band] = stellar_mass * these_lums 

    return kappa,lums

def stellar_raytrace(
    stellar_x,stellar_y,stellar_z, 
    stellar_mass,stellar_age,stellar_metallicity,
    stellar_hsml, 
    gas_x,gas_y,gas_z,
    gas_mass,gas_metallicity,
    gas_hsml, 
    kappa,
    lums,
    xlim=0,
    ylim=0,
    zlim=0,
    pixels=720, 
    KAPPA_UNITS=2.08854068444, ## cm^2/g -> kpc^2/mcode
    QUIET=False):

    ## check if stellar metallicity is a matrix
    ##  i.e. mass fraction of many species. If so,
    ##  take the total metallicity 
    if (len(stellar_metallicity.shape)>1): 
        stellar_metallicity=stellar_metallicity[:,0];
    if (len(gas_metallicity.shape)>1): 
        gas_metallicity=gas_metallicity[:,0];

    ## apply minimum metallicity and ages
    #stellar_metallicity[stellar_metallicity>0] += ADD_BASE_METALLICITY; ## TODO why >0?
    #gas_metallicity[gas_metallicity>0] += ADD_BASE_METALLICITY;
    #stellar_age += ADD_BASE_AGE;

    ## count particles we're using
    Nstars=len(np.array(stellar_mass))
    Ngas=len(np.array(gas_mass))

    ## dummy values to use for source and attenuation terms 
    gas_lum=np.zeros(Ngas) ## gas has no 'source term' for this calculation
    stellar_mass_attenuation = np.zeros(Nstars) ## stars have no 'attenuation term'

    ## convert units
    gas_mass_metal = gas_mass * (gas_metallicity/0.02)
    kappa *= KAPPA_UNITS
    
    ## combine the relevant arrays so it can all be fed into the ray-tracing
    ##  positions
    x=np.concatenate([stellar_x,gas_x])
    y=np.concatenate([stellar_y,gas_y])
    z=np.concatenate([stellar_z,gas_z])
    
    ##  masses for attenuation purposes
    mass=np.concatenate([stellar_mass_attenuation,gas_mass_metal])
    ##  smoothing lengths 
    hsml=np.concatenate([stellar_hsml,gas_hsml])

    ##  source terms in each band
    wt1=np.concatenate([lums[0,:],gas_lum])
    wt2=np.concatenate([lums[1,:],gas_lum])
    wt3=np.concatenate([lums[2,:],gas_lum])

    ## opacity in each band
    k1,k2,k3=kappa
        
    if not QUIET:
        print("Projecting with attenuation...")
        print('total lum before attenuation in each band (Lsun/1e10):',np.sum(lums,axis=1))
        print('opacity in each band:',kappa)
        print('total gas mass:',np.sum(gas_mass_metal))
    return raytrace_projection_compute(
        x,y,z,
        hsml,mass,
        wt1,wt2,wt3,
        k1,k2,k3,
        xlim=xlim,ylim=ylim,zlim=zlim,
        pixels=pixels,
        TRIM_PARTICLES=1)
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
def raytrace_projection_compute(
    x,y,z,
    hsml,mass,
    wt1,wt2,wt3,
    kappa_1,kappa_2,kappa_3,
    xlim=0,ylim=0,zlim=0,
    pixels=720,
    TRIM_PARTICLES=1):

    ## define bounaries
    if(checklen(xlim)<=1): 
        xlim=[np.min(x),np.max(x)]
    if(checklen(ylim)<=1): 
        ylim=[np.min(y),np.max(y)]
    if(checklen(zlim)<=1): 
        zlim=[np.min(z),np.max(z)]

    ## midpoint of box
    x00=0.5*(xlim[1]+xlim[0])
    y00=0.5*(ylim[1]+ylim[0])
    z00=0.5*(zlim[1]+zlim[0])

    ## re-center particles
    x-=x00
    y-=y00
    z-=z00

    ## half-width of box
    xlen=0.5*(xlim[1]-xlim[0])
    ylen=0.5*(ylim[1]-ylim[0])
    zlen=0.5*(zlim[1]-zlim[0])

    tolfac = 1.0e10 ## dummy, all particles will be in box
    if (TRIM_PARTICLES==1): 
        tolfac = 0.05

    ## determine which particles are in box
    dx=xlen*(1.+tolfac*2.)
    dy=ylen*(1.+tolfac*2.)
    dz=zlen*(1.+tolfac*2.)

    ## produce an 'in-box' "ok" mask
    ok=(ok_scan(x,xmax=dx) & 
        ok_scan(y,xmax=dy) & 
        ok_scan(z,xmax=dz) & 
        ok_scan(hsml,pos=1) & 
        ok_scan(mass+wt1+wt2+wt3,pos=1))

    ## apply "ok" mask
    x=x[ok]
    y=y[ok]
    z=z[ok]
    hsml=hsml[ok]
    mass=mass[ok]
    wt1=wt1[ok]
    wt2=wt2[ok]
    wt3=wt3[ok]

    ## limits of box
    xmin=-xlen
    xmax=xlen
    ymin=-ylen
    ymax=ylen

    N_p=checklen(x)
    if(N_p<=1): 
        print(
            'UH-OH: EXPECT ERROR NOW',
            'there are no valid source/gas particles to send!')
        return -1,-1,-1,-1;

    ## now sort these in z (this is critical!)
    ##  get sort indices
    s=np.argsort(z);

    ##  apply sort indices
    x,y,z=x[s],y[s],z[s]
    mass=mass[s]
    hsml=hsml[s]
    wt1,wt2,wt3=wt1[s],wt2[s],wt3[s]

    ## cast new copies to ensure the correct formatting when fed to the c-routine:
    ##  cast to single precision
    x,y,z=fcor(x),fcor(y),fcor(z)
    mass=fcor(mass)
    hsml=fcor(hsml)
    wt1,wt2,wt3=fcor(wt1),fcor(wt2),fcor(wt3)

    ## load the routine we need
    curpath = os.path.realpath(__file__)
    curpath = curpath[:len("utils")+curpath.index("utils")] #split off this filename
    exec_call=os.path.join(curpath,'stellar_utils/c_libraries/RayTrace_RGB/raytrace_rgb.so')
    routine=ctypes.cdll[exec_call];
    
    ## cast the variables to store the results
    aspect_ratio=ylen/xlen
    Xpixels=int_round(pixels)
    Ypixels=int_round(aspect_ratio*np.float(Xpixels))
    N_pixels=Xpixels*Ypixels

    ## create output array pointers
    out_cast=ctypes.c_float*N_pixels
    out_0=out_cast() ## mass map
    out_1=out_cast() ## band 1 
    out_2=out_cast() ## band 2
    out_3=out_cast() ## band 3

    ## main call to the attenuation routine in C
    routine.raytrace_rgb( 
        ctypes.c_int(N_p), ## number of star + gas particles
        vfloat(x),vfloat(y), ## x-y positions of star + gas particles
        vfloat(hsml),  ## smoothing lengths of star + gas particles
        vfloat(mass), ## attenuation masses of star + gas particles, stars are 0 
        ## emission in each band of star+gas particles, gas is 0 
        vfloat(wt1),vfloat(wt2),vfloat(wt3),  
        ## opacity in each band
        ctypes.c_float(kappa_1),ctypes.c_float(kappa_2),ctypes.c_float(kappa_3), 
        ## x-y limits of the image
        ctypes.c_float(xmin),ctypes.c_float(xmax),ctypes.c_float(ymin),ctypes.c_float(ymax), 
        ctypes.c_int(Xpixels),ctypes.c_int(Ypixels), ## output shape
        ctypes.byref(out_0), ## mass map
        ctypes.byref(out_1),ctypes.byref(out_2),ctypes.byref(out_3) ) ## band maps

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

