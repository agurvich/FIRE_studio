## returns the attenuation at a given frequency nu (in Hz) for a 
##    given column density (log(NH/cm^-2)), and metallicity (in solar units)
##
##		i.e. for some luminosity L, observed = L * attenuate(nu,alog10(NH),Z/0.02)
##
##	flags: can change the reddening curve adopted 
##		SMC (default) = smc-like reddening curve (strong UV & mid-IR extinction)
##			observations find this most appropriate for quasars: Hopkins et al. 2004
##		MW  = MW-like (weaker UV & mid-IR; strong 2200 angstrom bump) 
##		LMC = LMC-like (basically halfway between the two)
##			(gas to dust ratio in each case is MW-like average scaled 
##				linearly with metallicity)
##
##	can call for specific bands (i.e. net attenuation integrated over the band):
##		BB = B-band
##      IR = mid-IR (15 microns), 
##		SX = soft X-ray (0.5-2 keV), 
##		HX = hard X-ray (2-10 keV)
##
##  the dust extinction curves are from Pei et al. 1992.; gas-to-dust ratios 
##	  from Bouchet et al. 1985. 
##
##  above the Lyman edge, photoionization cross sections computed following 
##   Morrison & McCammon 1983. the code follows Brant's decomposition of this into 
##   a metals-free component and a metals component that scales linearly with 
##   the input metallicity
##
##  compton scattering is treated as non-relativistic achromatic Thompson scattering
##	 at low energies (<~ 4 keV), and above this the transmission curves from 
##   Matt, Pompilio, & La Franca 1995 have been extracted to interpolate over a wide 
##   range of frequencies and NH. These are decomposed into a H-dominated continuum 
##   component that basically scales in a metal-free manner, and an iron flourescence 
##   component (between ~5.8-7.4 keV) that emits in a manner scaling with metallicity. 
##   the iron flourescence is fine for most applications, but can give strange 
##   (although fixed s.t. non-divergent and still monotonic) behavior when you're  
##   simultaneously at very low but non-zero metallicities (<~ 0.1 solar) 
##   and very high (logNH > 25-25.5) column densities -- in this regime the 
##   flourescence should really be calculated self-consistently. 
##
##

import numpy as np
import math
import ctypes

def vdouble(x):
    return x.ctypes.data_as(ctypes.POINTER(ctypes.c_double));


def attenuate( nu_in_Hz, log_NH, metallicity_in_solar, \
	SMC=0, LMC=0, MW=0, BB=0, IR=0, SX=0, HX=0):

    import utilities as util

    ## location of shared library
    exec_call=util.return_python_routines_homedir()+'/attenuation/attenuate_py.so'
    lib=ctypes.cdll[exec_call];

    ## default to SMC-like reddening
    dust_key = 2
    if (MW==1):  dust_key = 0
    if (LMC==1): dust_key = 1
    if (SMC==1): dust_key = 2

    ## frequency
    if (BB==1): nu_in_Hz = -1.0
    if (IR==1): nu_in_Hz = -2.0
    if (SX==1): nu_in_Hz = -3.0
    if (HX==1): nu_in_Hz = -4.0	

    NH = 10.**np.array(log_NH,ndmin=1,dtype='d')
    nu_in_Hz=np.array(nu_in_Hz,ndmin=1,dtype='d');
    metallicity_in_solar=np.array(metallicity_in_solar,ndmin=1,dtype='d');
    N_nu = len(nu_in_Hz)
    N_NH = len(NH)
    N_metal = len(metallicity_in_solar)
    if (N_metal <= 1): metallicity_in_solar = 0.*NH + metallicity_in_solar[0]
    atten = np.zeros(N_nu*N_NH,dtype='d')
    NH = np.array(NH,ndmin=1,dtype='d')
    metallicity_in_solar = np.array(metallicity_in_solar,ndmin=1,dtype='d')
    out_cast=ctypes.c_double*(N_nu*N_NH); atten=out_cast();
    
    ## actual function call
    lib.main( ctypes.c_int(N_nu), vdouble(nu_in_Hz), \
        ctypes.c_int(N_NH), vdouble(NH), \
        vdouble(metallicity_in_solar), \
        ctypes.c_int(dust_key), \
        ctypes.byref(atten) );

    ## now put the output arrays into a useful format 
    atten = np.copy(np.ctypeslib.as_array(atten));

    atten_f = np.zeros((N_nu,N_NH),dtype='d')
    for i in range(N_nu):
        atten_f[i,:] = atten[i*N_NH:(i+1)*N_NH]
    atten=atten_f
    atten[atten==0]=1.0e-40; 
    atten[np.isnan(atten)]=1.0e-40;

    return atten

