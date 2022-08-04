import os
import numpy as np
import math
import scipy.ndimage.interpolation as interpolate
import struct

def read_band_lums_from_tables(
    BAND_NAMES,
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
    Nbands=len(np.array(BAND_NAMES))

    ## require that we attenuate 3 bands to combine since attenuation
    ##  routine is hardcoded to accept 3 weights
    if (Nbands != 3): 
        print("stellar_raytrace needs 3 bands, you gave ",Nbands)
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
                band_name=BAND_NAMES[i_band], ## band index
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
                band_name=BAND_NAMES[i_band], ## band index
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
    band_name='Bolometric',
    SALPETER_IMF=0,CHABRIER_IMF=1,
    QUIET=0,CRUDE=0, 
    RETURN_NU_EFF=0,RETURN_LAMBDA_EFF=0,
    UNITS_SOLAR_IN_BAND=0):

    age_in_Gyr=np.array(age_in_Gyr,ndmin=1)
    metallicity_in_solar_units=np.array(metallicity_in_solar_units,ndmin=1)

    # band_names pulled from fsps.get_filters, to get your favorite band name:
    #   either: a) import fsps, run fsps.find_filter(''), b) look at https://dfm.io/python-fsps/current/filters/

    # MEO commentary: this code only accounts for attenuation, so don't be an idiot making maps in like the Spitzer
    #                 filters and expect that it's getting the reprocessed dust radiation right.  I have included all
    #                 these filters since python FSPS includes them and it could somehow be useful, not because they
    #                 all make sense to use, or are vetted in any way.
    if not QUIET: print('Calculating L/M in ' + band_name)

    ## find directory of colors.dat file for a given IMF
    curpath = os.path.realpath(__file__)
    ##  split off this filename
    curpath = curpath[:len("utils")+curpath.index("utils")] 
    ##  directory in which the data binaries are stored
    froot = os.path.join(curpath,'stellar_utils','colors_sps/') 

    # NOTE: right now just using Chabrier FSPS pops (since I'm lazy, and that's what FIRE uses)
    #if (CHABRIER_IMF==1): fname=froot+'colors.chabrier.dat'
    if (SALPETER_IMF==1):
        print("NUR EIN KLEINES PROBLEM: SALPETER_IMF NOT IMPLEMENTED! SCHAUEN SIE AN DIE FSPS DOKUMENTATION.")

    # Find/Build FSPS Z-age grid of magnitudes for selected band
    min_age_Gyr = np.min(age_in_Gyr)
    max_age_Gyr = np.max(age_in_Gyr)

    min_metallicity = np.max([np.min(metallicity_in_solar_units), 0.02 * 1e-4]) # don't let it go < 1e-4 Solar
    max_metallicity = np.max(metallicity_in_solar_units)

    NZbins = 36
    forceFSPSbuild = 1
    logZ_pts_array = np.linspace(-2.5, 1., num=NZbins)  # not sure how many Z pts we need?, in log solar Z units

    # check if FSPS tables already exist, open or create
    fsps_colors = h5py.File(froot+"fsps_colors.hdf5", 'a')
    if band_name in fsps_colors.keys():
        if (QUIET == 0): print("FSPS Table Found (BAND: %s)!!" % band_name)
        # file exists, and the filter/band has a table!
        lband = np.array(fsps_colors[band_name]['lband'])
        (Na, Nz) = lband.shape
        age_pts_array = fsps_colors["SSP_ages"]  # age grid of SSPs, in Gyr

        if (NZbins == Nz):
            # correct metallacity array shape! assuming that the table is valid
            forceFSPSbuild = 0
            if (QUIET == 0): print("FSPS Table valid shape!!")

    if forceFSPSbuild:
        # actually need FSPS here
        try: import fsps
        except: print("No FSPS module! Download and install pythonFSPS (https://dfm.io/python-fsps)...")
        # file either didn't exist or the filter/band didn't have a table!
        if not band_name in fsps.find_filter(''):
            if not band_name == 'Bolometric':
                print("BAND: %s not supported by (this? build of) FSPS. I'm going to break now..." % band_name)
                return -1
        if (QUIET == 0): print("No valid FSPS Tables Found (BAND: %s)!! Doing it ourselves" % band_name)
        if (QUIET == 0): print("Initializing FSPS SSPs...")
        if (QUIET == 0): print("Computing SSP for logZ=", logZ_pts_array[0])
        sp = fsps.StellarPopulation(zcontinuous=1, logzsol=logZ_pts_array[0], tage=0.0)
        sp.params['imf_type'] = 1 # Chabrier 2003 IMF (god help us if you need another IMF)

        age_pts_array = 10. ** (sp.ssp_ages - 9.)  # age grid of SSPs, in Gyr
        Na = len(age_pts_array)
        Nz = len(logZ_pts_array)

        bol_bool = 0
        if band_name == 'Bolometric':
            bol_bool = 1
            band_name = 'sdss_u' # dummy filter, since 'Bolometric' isn't a filter in FSPS
        lband = np.zeros((Na, Nz), dtype=np.float64);
        lband[:, 0] = 10.**sp.get_mags(bands=[band_name])[:, 0]
        if bol_bool: lband[:, 0] = sp.log_lbol
        # build array of mags for band in age/Z space
        for zz in range(1, len(logZ_pts_array)):
            if (QUIET == 0): print("Computing SSP for logZ=", logZ_pts_array[zz])
            sp.params['logzsol'] = logZ_pts_array[zz]
            lband[:, zz] = sp.get_mags(bands=[band_name])[:, 0] # in M_ab/Msun units
            if bol_bool: lband[:, zz] = 10.**sp.log_lbol[:, 0]

        # convert to solar Luminosity units
        if not bol_bool: lband = 10. ** ((lband - fsps.get_filter(band_name).msun_ab) / -2.5) # now in Lsun/Msun units!
        # Now saving the dataset under BAND_NAME, closing the file.
        fsps_colors.create_group(band_name)
        fsps_colors.create_dataset(band_name+'/lband', data=lband)
        if bol_bool:
            fsps_colors.create_dataset(band_name+'/msun_ab', data=4.74)
            fsps_colors.create_dataset(band_name+'/lambda_eff', data=4243.93)
        else:
            fsps_colors.create_dataset(band_name + '/msun_ab', data=fsps.get_filter(band_name).msun_ab)
            fsps_colors.create_dataset(band_name + '/lambda_eff', data=fsps.get_filter(band_name).lambda_eff)
        if not ("SSP_ages" in fsps_colors.keys()):
            fsps_colors.create_dataset("SSP_ages", data=age_pts_array)
        if not ("logZ_pts" in fsps_colors.keys()):
            fsps_colors.create_dataset("logZ_pts", data=logZ_pts_array)


        if (QUIET == 0):
            if bol_bool: print("FSPS Tables Built (BAND: Bolometric)!!")
            else: print("FSPS Tables Built (BAND: %s)!!" % band_name)
        # done with FSPS
    # load out filter properties
    lam_eff = float(np.array(fsps_colors[band_name]['lambda_eff']))
    mag_sun_ab = float(np.array(fsps_colors[band_name]['msun_ab']))

    fsps_colors.close()

    # I know that having this after the 'load lband' block means it's less efficient.. but it's better conceptually
    if RETURN_NU_EFF or RETURN_LAMBDA_EFF:
        nu_eff = 2.998e18 / lam_eff ## c/lambda
        if RETURN_NU_EFF: return nu_eff
        if RETURN_LAMBDA_EFF: return lam_eff
# yes
    # allow for extreme metallicities (extrapolate linearly past table)
    push_metals = 1;
    # we're working with Lsun/Msun NOT log(Lsun/Msun) here.. might fuss things up
    if (push_metals == 1):
        Nz = Nz + 1;
        z_ext = [3.];
        logZ_pts_array = np.concatenate([logZ_pts_array, z_ext])
        lb1 = lband[:, Nz - 3]
        lb2 = lband[:, Nz - 2]
        lbx = np.zeros((Na, Nz), dtype=np.float64)
        lbx[:, 0:Nz - 1] = lband
        lbx[:, Nz - 1] = (lb2 - lb1) / (logZ_pts_array[Nz - 2] - logZ_pts_array[Nz - 3]) * \
                         (logZ_pts_array[Nz - 1] - logZ_pts_array[Nz - 2])  # might be really wrong.. i was copy-pasting
        lband = lbx;

    # get the x-axis (age) locations of input points
    ia_pts=np.interp(age_in_Gyr,age_pts_array,np.arange(0,Na,1));

    # this returns the boundary values for points outside of them (no extrapolation)
    #f=interp.interp1d(age_grid,np.arange(0,Na,1),kind='linear') 
    #ia_pts=f(np.log10(age_in_Gyr)+9.0)

    # get the y-axis (metallicity) locations of input points
    iz_pts=np.interp(np.log10(metallicity_in_solar_units), logZ_pts_array,np.arange(0,Nz,1));

    #f=interp.interp1d(np.log10(z_grid),np.arange(0,Nz,1),kind='linear') 
    #iz_pts=f(np.log10(metallicity_in_solar_units*zsun))

    ## map requested points to their location on the grid.
    ##  TODO i think these are index maps, in which case 
    ##  np.digitize would do this
    if (CRUDE==1):
        ia_pts=np.around(ia_pts).astype(int)
        iz_pts=np.around(iz_pts).astype(int)
        l_b=lband[ia_pts,iz_pts]
    else:
        l_b = interpolate.map_coordinates(lband, (ia_pts,iz_pts), order=1)
	
    # at this point, output is currently L/M in L_sun_IN_THE_BAND_OF_INTEREST/M_sun, 
    # but we want our default to be L/M in units of L_bolometric/M_sun = 3.9e33/2.0e33, so 
    #   need to get rid fo the L_sun_IN_THE_BAND_OF_INTEREST/L_bolometric
    if not UNITS_SOLAR_IN_BAND:
        nu_eff = 2.998e18 / lam_eff  # converts to nu_eff in Hz

        ten_pc = 10.e0 * 3.086e18  # 10 pc in cm
        log_S_nu = -(mag_sun_ab + 48.6) / 2.5  # zero point definition for ab magnitudes
        S_nu = 10. ** log_S_nu  # get the S_nu at 10 pc which defines M_AB
        lnu_sun_band = S_nu * (4. * math.pi * ten_pc * ten_pc)  # multiply by distance modulus
        nulnu_sun_band = lnu_sun_band * nu_eff  # multiply by nu_eff to get nu*L_nu

        l_bol_sun = 10. ** (-(4.74 + 48.6) / 2.5) * (4. * math.pi * ten_pc * ten_pc) * 2.998e18 / 1.e-5 # done by hand

        l_b *= nulnu_sun_band / l_bol_sun

    return l_b