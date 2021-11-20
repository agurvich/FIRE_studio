import os
import numpy as np
import ctypes

def checklen(x): return np.array(x,ndmin=1).shape[0]

def int_round(x): return np.int(np.round(x))

def ok_scan(arr,xmax=1.0e30,pos=0):
    if (pos==0): return np.logical_and(
        np.isfinite(arr),
        np.fabs(arr)<=xmax)
    if (pos==1): return np.logical_and(
        np.logical_and(
            np.isfinite(arr),
            np.fabs(arr)<=xmax),
        arr > 0.)

def farr(x): return np.array(x,dtype='f',ndmin=1)

def f_cpointer(x): return x.ctypes.data_as(ctypes.POINTER(ctypes.c_float))

def stellar_raytrace(
    stellar_x,stellar_y,stellar_z, 
    stellar_hsml, 
    gas_x,gas_y,gas_z,
    gas_mass,gas_metallicity,
    gas_hsml, 
    kappa,lums,
    xlim=0,ylim=0,zlim=0,
    pixels=720, 
    QUIET=False):
    """
    routine to use 'raytrace_projection_compute' to make mock stellar images, 
    treating the starlight as sources and accounting for gas extinction via ray-tracing
    """

    ## check if metallicity is a matrix
    ##  i.e. mass fraction of many species. If so,
    ##  take the total metallicity which should be the first column
    if (len(gas_metallicity.shape)>1): gas_metallicity=gas_metallicity[:,0]

    ## count particles we're using
    Nstars = np.array(stellar_x).shape[0]
    Ngas = np.array(gas_mass).shape[0]

    ## dummy values to use for source and attenuation terms 
    gas_lum = np.zeros(Ngas) ## gas has no 'source term' for this calculation
    stellar_mass_attenuation = np.zeros(Nstars) ## stars have no 'attenuation term'

    ## convert from metal mass fraction in solar units to actual metal mass
    gas_mass_metal = gas_mass * (gas_metallicity/0.02)
    
    ## combine the relevant arrays so it can all be fed into the ray-tracing
    ##  positions
    x = np.concatenate([stellar_x,gas_x])
    y = np.concatenate([stellar_y,gas_y])
    z = np.concatenate([stellar_z,gas_z])
    
    ##  masses for attenuation purposes
    mass = np.concatenate([stellar_mass_attenuation,gas_mass_metal])
    ##  smoothing lengths 
    hsml = np.concatenate([stellar_hsml,gas_hsml])

    ##  source terms in each band
    wt1 = np.concatenate([lums[0,:],gas_lum])
    wt2 = np.concatenate([lums[1,:],gas_lum])
    wt3 = np.concatenate([lums[2,:],gas_lum])

    ## unpack the opacity in each band
    k1,k2,k3 = kappa
        
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
        xlim,ylim,zlim,
        pixels)

def raytrace_projection_compute(
    x,y,z,
    hsml,mass,
    wt1,wt2,wt3,
    kappa_1,kappa_2,kappa_3,
    xlim=0,ylim=0,zlim=0,
    pixels=720,
    TRIM_PARTICLES=1):
    """
    Wrapper for raytrace_rgb, program which does a simply line-of-sight projection 
    with multi-color source and self-extinction along the sightline
    """


    ## define bounaries
    if(checklen(xlim)<=1): xlim=[np.min(x),np.max(x)]
    if(checklen(ylim)<=1): ylim=[np.min(y),np.max(y)]
    if(checklen(zlim)<=1): zlim=[np.min(z),np.max(z)]

    ## midpoint of box
    x00 = 0.5*(xlim[1]+xlim[0])
    y00 = 0.5*(ylim[1]+ylim[0])
    z00 = 0.5*(zlim[1]+zlim[0])

    ## re-center particles
    x -= x00
    y -= y00
    z -= z00

    ## half-width of box
    xlen = 0.5*(xlim[1]-xlim[0])
    ylen = 0.5*(ylim[1]-ylim[0])
    zlen = 0.5*(zlim[1]-zlim[0])

    if (TRIM_PARTICLES==1): tolfac = 0.05
    else: tolfac = 1.0e10 ## dummy, all particles will be in box

    ## determine which particles are in box
    dx = xlen*(1.+tolfac*2.)
    dy = ylen*(1.+tolfac*2.)
    dz = zlen*(1.+tolfac*2.)

    ## produce an 'in-box' "ok" mask
    ok_mask = (
        ok_scan(x,xmax=dx) & 
        ok_scan(y,xmax=dy) & 
        ok_scan(z,xmax=dz) & 
        ok_scan(hsml,pos=1) & 
        ok_scan(mass+wt1+wt2+wt3,pos=1))

    ## apply "ok" mask
    x = x[ok_mask]
    y = y[ok_mask]
    z = z[ok_mask]
    hsml = hsml[ok_mask]
    mass = mass[ok_mask]
    wt1 = wt1[ok_mask]
    wt2 = wt2[ok_mask]
    wt3 = wt3[ok_mask]

    ## limits of box
    xmin=-xlen
    xmax=xlen
    ymin=-ylen
    ymax=ylen

    num_points = checklen(x)
    if num_points <= 1: 
        print(
            'UH-OH: EXPECT ERROR NOW',
            'there are no valid source/gas particles to send!')
        return -1,-1,-1,-1

    ## now sort these in z (this is critical!)
    ##  get sort indices
    sort_indices = np.argsort(z)

    ##  apply sort indices
    x,y,z = x[sort_indices],y[sort_indices],z[sort_indices]
    mass = mass[sort_indices]
    hsml = hsml[sort_indices]
    wt1,wt2,wt3 = wt1[sort_indices],wt2[sort_indices],wt3[sort_indices]

    ## cast new copies to ensure the correct formatting when fed to the c-routine:
    ##  cast to single precision
    x,y,z = farr(x),farr(y),farr(z)
    mass = farr(mass)
    hsml = farr(hsml)
    wt1,wt2,wt3 = farr(wt1),farr(wt2),farr(wt3)

    ## load the routine we need
    curpath = os.path.realpath(__file__)
    curpath = curpath[:len("utils")+curpath.index("utils")] #split off this filename
    exec_call=os.path.join(curpath,'C_routines/RayTrace_RGB/raytrace_rgb.so')
    routine=ctypes.cdll[exec_call]
    
    ## cast the variables to store the results
    aspect_ratio = ylen/xlen
    Xpixels = int_round(pixels)
    Ypixels = int_round(aspect_ratio*np.float(Xpixels))
    N_pixels = Xpixels*Ypixels

    ## create output array pointers
    out_cast=ctypes.c_float*N_pixels
    out_0=out_cast() ## mass map
    out_1=out_cast() ## band 1 
    out_2=out_cast() ## band 2
    out_3=out_cast() ## band 3

    ## main call to the attenuation routine in C
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
    routine.raytrace_rgb( 
        ctypes.c_int(num_points), ## number of star + gas particles
        f_cpointer(x),f_cpointer(y), ## x-y positions of star + gas particles
        f_cpointer(hsml),  ## smoothing lengths of star + gas particles
        f_cpointer(mass), ## attenuation masses of star + gas particles, stars are 0 
        ## emission in each band of star+gas particles, gas is 0 
        f_cpointer(wt1),f_cpointer(wt2),f_cpointer(wt3),  
        ## opacity in each band
        ctypes.c_float(kappa_1),ctypes.c_float(kappa_2),ctypes.c_float(kappa_3), 
        ## x-y limits of the image
        ctypes.c_float(xmin),ctypes.c_float(xmax),ctypes.c_float(ymin),ctypes.c_float(ymax), 
        ctypes.c_int(Xpixels),ctypes.c_int(Ypixels), ## output shape
        ctypes.byref(out_0), ## mass map
        ctypes.byref(out_1),ctypes.byref(out_2),ctypes.byref(out_3) ) ## band maps

    ## now put the output arrays into a useful format 
    out_0 = np.copy(np.ctypeslib.as_array(out_0))
    out_1 = np.copy(np.ctypeslib.as_array(out_1))
    out_2 = np.copy(np.ctypeslib.as_array(out_2))
    out_3 = np.copy(np.ctypeslib.as_array(out_3))
    out_0 = out_0.reshape([Xpixels,Ypixels])
    out_1 = out_1.reshape([Xpixels,Ypixels])
    out_2 = out_2.reshape([Xpixels,Ypixels])
    out_3 = out_3.reshape([Xpixels,Ypixels])

    return out_0, out_1, out_2, out_3