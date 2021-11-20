import os
import numpy as np
import math
import scipy.ndimage.interpolation as interpolate
import struct

def read_band_lums_from_tables(
    BAND_IDS, 
    stellar_mass,stellar_age,stellar_metallicity,
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
        return -1,-1,-1

    ## check if stellar metallicity is a matrix
    ##  i.e. mass fraction of many species. If so,
    ##  take the total metallicity 
    if (len(stellar_metallicity.shape)>1): 
        stellar_metallicity=stellar_metallicity[:,0]

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

    return lums,nu_effs

def colors_table(
    age_in_Gyr,
    metallicity_in_solar_units, 
    BAND_ID=0,
    SALPETER_IMF=0,CHABRIER_IMF=1,
    QUIET=0,CRUDE=0, 
    RETURN_NU_EFF=0,RETURN_LAMBDA_EFF=0,
    UNITS_SOLAR_IN_BAND=0):

    age_in_Gyr=np.array(age_in_Gyr,ndmin=1)
    metallicity_in_solar_units=np.array(metallicity_in_solar_units,ndmin=1)

    ## remap BAND_ID to a different ordering @_@
    ## ordering I'm used to
    #j = [  0,  6,  7,  8,  9, 10, 11, 12, 13,  1,   2,   3,   4,   5] 
    #BAND_ID = j[BAND_ID]

    band_names=['Bolometric',
        'Sloan u','Sloan g','Sloan r','Sloan i','Sloan z', 
        'Johnsons U','Johnsons B', 'Johnsons V','Johnsons R','Johnsons I', 
        'Cousins J','Cousins H','Cousins K']


    lam_eff=np.array([
            ## Bolometric (?)
                    1.e-5,  
            ## SDSS u       g       r      i      z
                    3551. , 4686. , 6165., 7481., 8931.,
            ##      U       B       V      R
                    3600. , 4400. , 5556., 6940., 
            ##      I      J       H       K
                    8700., 12150., 16540., 21790.])

    if not QUIET: 
        print('Calculating L/M in ' +
            band_names[BAND_ID] + 
            ' (BAND_ID=%d,l=%d A)'%(BAND_ID,lam_eff[BAND_ID]))
    
    if RETURN_NU_EFF or RETURN_LAMBDA_EFF:
    
        nu_eff = 2.998e18 / lam_eff ## c/lambda
        if RETURN_NU_EFF: 
            return nu_eff[BAND_ID]
        if RETURN_LAMBDA_EFF: 
            return lam_eff[BAND_ID]
    
    ## find directory of colors.dat file for a given IMF
    curpath = os.path.realpath(__file__)
    ##  split off this filename
    curpath = curpath[:len("utils")+curpath.index("utils")] 
    ##  directory in which the data binaries are stored
    froot = os.path.join(curpath,'stellar_utils','colors_sps/') 

    if (CHABRIER_IMF==1): fname=froot+'colors.chabrier.dat'
    if (SALPETER_IMF==1): fname=froot+'colors.salpeter.dat'

    ## read color data from binary 
    ##  TODO consider rewriting these tables to disk as text...
    with open(fname,'rb') as lut:
        lut_dat = lut.read()
        Nl,Na,Nz = struct.unpack('3i',lut_dat[0:12])
        z_grid = np.array(struct.unpack(str(Nz)+'d',lut_dat[12:12+8*Nz]))
        age_grid = np.array(struct.unpack(str(Na)+'d',lut_dat[12+8*Nz:12+8*Nz+8*Na]))
        l_all_l = np.array(struct.unpack(str(Nl*Na*Nz)+'d',lut_dat[12+8*Nz+8*Na:12+8*Nz+8*Na+8*Nl*Na*Nz]))
        l_all = np.transpose(l_all_l.reshape(Nz,Na,Nl))
    
    ## pick the luminosity in this band
    l_band = np.zeros((Na,Nz),dtype=np.float64)
    for iz in range(Nz): 
        l_band[:,iz]=l_all[BAND_ID,:,iz]
    ## TODO i think this is the same as l_band = l_all[band]...
    
    # allow for extreme metallicities (extrapolate linearly past table)
    push_metals = 1
    if (push_metals==1):
        Nz = Nz + 1
        z_ext = [1000.0]
        z_grid = np.concatenate([z_grid,z_ext])
        lb1 = l_band[:,Nz-3]
        lb2 = l_band[:,Nz-2]
        lbx = np.zeros((Na,Nz),dtype=np.float64)
        lbx[:,0:Nz-1] = l_band
        lbx[:,Nz-1] = ( (lb2 - lb1) / 
            np.log10(z_grid[Nz-2]/z_grid[Nz-3]) * 
            np.log10(z_grid[Nz-1]/z_grid[Nz-2]) ) 
        l_band = lbx

    # get the x-axis (age) locations of input points
    ia_pts=np.interp(
        np.log10(age_in_Gyr)+9.0,
        age_grid,
        np.arange(0,Na,1))

    # this returns the boundary values for points outside of them (no extrapolation)
    #f=interp.interp1d(age_grid,np.arange(0,Na,1),kind='linear') 
    #ia_pts=f(np.log10(age_in_Gyr)+9.0)

    # get the y-axis (metallicity) locations of input points
    zsun = 0.02
    iz_pts=np.interp(
        np.log10(metallicity_in_solar_units*zsun),
        np.log10(z_grid),
        np.arange(0,Nz,1))

    #f=interp.interp1d(np.log10(z_grid),np.arange(0,Nz,1),kind='linear') 
    #iz_pts=f(np.log10(metallicity_in_solar_units*zsun))

    ## map requested points to their location on the grid.
    ##  TODO i think these are index maps, in which case 
    ##  np.digitize would do this
    if (CRUDE==1):
        ia_pts=np.around(ia_pts).astype(int)
        iz_pts=np.around(iz_pts).astype(int)
        l_b=l_band[ia_pts,iz_pts]
    else:
        l_b = interpolate.map_coordinates(l_band, (ia_pts,iz_pts), order=1)
    l_b = 10.**l_b
	
    # at this point, output is currently L/M in L_sun_IN_THE_BAND_OF_INTEREST/M_sun, 
    # but we want our default to be L/M in units of L_bolometric/M_sun = 3.9e33/2.0e33, so 
    #   need to get rid fo the L_sun_IN_THE_BAND_OF_INTEREST/L_bolometric
    if not UNITS_SOLAR_IN_BAND:
        l_b = renormalize_band_luminosity(l_b,BAND_ID)

    return l_b

def renormalize_band_luminosity(l_b,BAND_ID):
    # AB system solar luminosities used for determining L_sun in absolute units for each of these
    N_BANDS=14
    mag_sun_ab = np.zeros(N_BANDS,dtype=float)
    mag_sun_ab[0] = 4.74  
    mag_sun_ab[1] = 6.75  #SDSS u (unprimed AB)
    mag_sun_ab[2] = 5.33  #SDSS g (unprimed AB)
    mag_sun_ab[3] = 4.67  #SDSS r (unprimed AB)
    mag_sun_ab[4] = 4.48  #SDSS i (unprimed AB)
    mag_sun_ab[5] = 4.42  #SDSS z (unprimed AB)
    mag_sun_ab[6] = 6.34  #U (BESSEL)
    mag_sun_ab[7] = 5.33  #B (BESSEL)
    mag_sun_ab[8] = 4.81  #V (BESSEL)
    mag_sun_ab[9] = 4.65  #R (KPNO)
    mag_sun_ab[10] = 4.55 #I (KPNO)
    mag_sun_ab[11] = 4.57 #J (BESSEL)
    mag_sun_ab[12] = 4.71 #H (BESSEL)
    mag_sun_ab[13] = 5.19 #K (BESSEL)

    # Effective wavelengths of the bands (in Angstroms), to compute nuLnu<->Lnu
    # UBVRIJHK from https://cass.ucsd.edu/archive/physics/ph162/mags.html
    # SDSS ugriz from http://www.sdss.org/dr4/instruments/imager/index.html#filters
    lambda_eff=np.array([
        ## Bolometric (?)
                1.e-5,  
        ## SDSS u       g       r      i      z
                3551. , 4686. , 6165., 7481., 8931.,
        ##      U       B       V      R
                3600. , 4400. , 5556., 6940., 
        ##      I      J       H       K
                8700., 12150., 16540., 21790.])

    c_light = 2.998e10 # speed of light in cm/s
    nu_eff  = c_light / (lambda_eff * 1.0e-8) # converts to nu_eff in Hz

    ten_pc   = 10.e0 * 3.086e18 # 10 pc in cm
    log_S_nu = -(mag_sun_ab + 48.6)/2.5 # zero point definition for ab magnitudes
    S_nu     = 10.**log_S_nu # get the S_nu at 10 pc which defines M_AB
    lnu_sun_band = S_nu * (4.*math.pi*ten_pc*ten_pc) # multiply by distance modulus 
    nulnu_sun_band = lnu_sun_band * nu_eff # multiply by nu_eff to get nu*L_nu
    l_bol_sun = nulnu_sun_band[0]

    l_b *= nulnu_sun_band[BAND_ID] / l_bol_sun

    return l_b
