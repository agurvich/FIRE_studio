import numpy as np

#-------------------------------------------------------------------------------------
# Pei (1992) MW,LMC,SMC dust laws used to calculate differential extinction
#   for dust at most frequencies
# (INPUT LAMBDA IN MICRONS)
#     NOTE: Pei's model is based on observations from lambda > 1000 Angstroms
#       (f < 3.0e15 Hz), but the grain models extrapolate these fits at least to
#		lambda = 0.001 micrometer (= 10 Angstroms, or f = 3.0e17 Hz), where 
#		the extinction drops off rapidly. So it's probably safe to use for all freq.
#-------------------------------------------------------------------------------------
def pei_dustparam( lambda_microns, MW=0, LMC=0, SMC=1 ):

    ## see Table 4 of Pei92
    if (MW==1):
        ##  BKG    FUV   2175 A 9.7 um 18 um  FIR
        a = [165. , 14. , 0.045, 0.002, 0.002, 0.012]
        l = [0.047, 0.08, 0.22 , 9.7  , 18.  , 25.  ]
        b = [90.  , 4.00, -1.95, -1.95, -1.80, 0.00 ]
        n = [2.0  , 6.5 , 2.0  , 2.0  , 2.0  , 2.0  ]
        R_V = 3.08 ## Table 2
    if (LMC==1):
        ##  BKG    FUV   2175 A 9.7 um 18 um  FIR
        a = [175. , 19. , 0.023, 0.005, 0.006, 0.020]
        l = [0.046, 0.08, 0.22 , 9.7  , 18.  , 25.  ]
        b = [90.  , 5.50, -1.95, -1.95, -1.80, 0.00 ]
        n = [2.0  , 4.5 , 2.0  , 2.0  , 2.0  , 2.0  ]
        R_V = 3.16 ## Table 2
    if (SMC==1):
        ##  BKG    FUV   2175 A 9.7 um 18 um  FIR
        a = [185. , 27. , 0.005, 0.010, 0.012, 0.030]
        l = [0.042, 0.08, 0.22 , 9.7  , 18.  , 25.  ]
        b = [90.  , 5.50, -1.95, -1.95, -1.80, 0.00 ]
        n = [2.0  , 4.0 , 2.0  , 2.0  , 2.0  , 2.0  ]
        R_V = 2.93 ## Table 2

    xsi = 0.*lambda_microns;
    for i in range(len(a)):
        xsi += a[i] / ( (lambda_microns/l[i])**(n[i]) + (l[i]/lambda_microns)**(n[i]) + b[i] )
    R_lam = (1.0 + R_V) * xsi;

    return xsi;


#-------------------------------------------------------------------------------------
# For 0.03 keV < E < 10 keV  
#   (7.2e15 < nu[Hz] < 2.4e18  or   1.2 < lambda[Angstroms] < 413)
#   we use the photoelectric absorption cross sections of 
#   Morrison & McCammon (1983)
#     NOTE: these assume solar abundances and no ionization, 
#             the appropriate number probably scales linearly with both
#				- keyword ABUNDANCE gives ratio to solar/MW metallicity
#   (this is all for the COMPTON THIN regime)
#-------------------------------------------------------------------------------------
def morrison_photoelec( frequency ):
    # (can restructure this to allow vector i/o if the time for the run is too long)
    f_003keV = 7.253e15
    f_10keV  = 2.418e18
    keV_per_Hz = 4.13608e-18    
    fkeV = frequency * keV_per_Hz	# convert frequency from Hz to keV

    # Now set the appropriate polynomial terms from Table 2 
    #   of Morrison & McCammon (for a given frequency range)
    if  (fkeV < 0.100):
        c = [	17.3,	608.1,	-2150.0	]
    elif(fkeV < 0.284):
        c = [	34.6,	267.9,	 -476.1	]
    elif(fkeV < 0.400):
        c = [	78.1,	 18.8,	    4.3	]
    elif(fkeV < 0.532):
        c = [	71.4,	 66.8,	  -51.4	]
    elif(fkeV < 0.707): 
        c = [	95.5,	145.8,	  -61.1	]
    elif(fkeV < 0.867): 
        c = [   308.9,  -380.6,	  294.0	]
    elif(fkeV < 1.303): 
        c = [   120.6,	169.3,	  -47.7	]
    elif(fkeV < 1.840): 
        c = [   141.3,	146.8,	  -31.5	]
    elif(fkeV < 2.471): 
        c = [   202.7,	104.7,	  -17.0	]
    elif(fkeV < 3.210): 
        c = [   342.7,	 18.7,	    0.0	]
    elif(fkeV < 4.038): 
        c = [   352.2,	 18.7,	    0.0	]
    elif(fkeV < 7.111): 
        c = [   433.9,	 -2.4,	  -0.75	]
    elif(fkeV < 8.331): 
        c = [   629.0,	 30.9,      0.0	]
    elif(fkeV < 10.00): 
        c = [   701.2,	 25.2,	    0.0	]
    else:
        c = [   701.2,	 25.2,	    0.0	]

    # Use these coefficients to calculate the cross section (in 10^-24 cm^2) per hydrogen atom
    return (c[0] + c[1]*fkeV + c[2]*(fkeV*fkeV)) / (fkeV*fkeV*fkeV)  * 1.0e-24;



#--------------------------------------------------------------------------
# Function to calculate the cross section per H atom at a given 
#    frequency f (in Hz)
#--------------------------------------------------------------------------
def cross_section( f0, METALLICITY_OVER_SOLAR=1.0 ):
    f0=np.array(f0);
    SIGMA = 0.0*np.array(f0);

    # For 0.03 keV < E < 10 keV  
    #   (7.2e15 < nu[Hz] < 2.4e18  or   1.2 < lambda[Angstroms] < 413)
    #   we use the photoelectric absorption cross sections of 
    #   Morrison & McCammon (1983)
    #     NOTE: these assume solar abundances and no ionization, 
    #             the appropriate number probably scales linearly with both
    #   (this is all for the COMPTON THIN regime)
    #
    f_003keV = 7.253e15
    f_H_edge = 1.362/3.0 * f_003keV
    f_10keV  = 2.418e18
    if (SIGMA.size > 1):
        for i in range(SIGMA.size): 
            if (f_H_edge <= f0[i]): SIGMA[i] += morrison_photoelec(f0[i]);
    else:
        if (f_H_edge <= f0): SIGMA += morrison_photoelec(f0);
        SIGMA *= METALLICITY_OVER_SOLAR;
        # (make sure this uses the total column, not just the 
        #   neutral column, as ions still contribute), otherwise divide by mean (~1/3,say)

    # For optical-IR regions, we use the Pei numerical approximations below.
    #
    # xsi = tau(lambda)/tau(B) is the ratio of extinction at lambda to the 
    #    extinction in the B-band. 
    # k = 10^21 (tau_B / NH)   (NH in cm^2) gives the dimensionless gas-to-dust
    #    ratio, with k=0.78 for MW, k=0.16 for LMC, k=0.08 for SMC.
    #    k is INDEPENDENT of the grain properties, and seems to scale rougly
    #    linearly with metallicity
    # so, for now, assume solar metallicity and k = k_MW = 0.78. we can rescale later.
    #
    # tau_B = ( NH / (10^21 cm^-2) ) * k --> SIGMA_B = k*10^-21  cm^2
    # tau_lambda = xsi * tau_B --> SIGMA = xsi * SIGMA_B
    #
    # k = 0.78 for the MW
    # k = 0.08 for the SMC, approximately in line with the MW/LMC/SMC metallicity 
    #  sequence, so we take a k_MW then scaled by the metallicity
    #
    k = 0.78 * METALLICITY_OVER_SOLAR
    if (SIGMA.size > 1):
        ok=(f0 < f_003keV)
        lambda_microns = 2.998e14 / f0[ok];   # convert frequency [Hz] to wavelength [microns]
        xsi = pei_dustparam( lambda_microns, SMC=1 );
        SIGMA[ok] += xsi * k * 1.0e-21;
    else:
        if (f0 < f_003keV):
            lambda_microns = 2.998e14 / f0;   # convert frequency [Hz] to wavelength [microns]
            xsi = pei_dustparam( lambda_microns, SMC=1 );
            SIGMA += xsi * k * 1.0e-21;

    # No double-counting, but I have checked it in detail, and there is really almost
    #    absolutely no difference -- even up to NH~10^24-25, it's like a factor of 1.1
    #    or so between including both these factors and not, and only in a very, very 
    #    narrow frequency range

    # + Thompson scattering cross sections (NR Compton scattering)
    #    full compton scattering dsigma/dOmega = (1/2)*r0^2*(1+cos^2(theta)) 
    #    ( sigma_thompson = (8pi/3)*r0^2 )

    sigma_thompson = 6.65e-25	## cm^-2
    SIGMA += sigma_thompson * 0.5  
    # b/c half scattered into viewing angle, half out (*very* roughly) (reflection)
                                           
    return SIGMA


## simple routine for opacity scaled to solar metallicity at frequencies f0
def opacity_per_solar_metallicity( f0 ):
    mu = 0.75 ## assume dense gas
    return cross_section( f0, METALLICITY_OVER_SOLAR=1.0 ) / ( mu * 1.67e-24 );
    
