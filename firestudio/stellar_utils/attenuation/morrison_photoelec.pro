;-------------------------------------------------------------------------------------
; For 0.03 keV < E < 10 keV  
;   (7.2e15 < nu[Hz] < 2.4e18  or   1.2 < lambda[Angstroms] < 413)
;   we use the photoelectric absorption cross sections of 
;   Morrison & McCammon (1983)
;     NOTE: these assume solar abundances and no ionization, 
;             the appropriate number probably scales linearly with both
;				- keyword ABUNDANCE gives ratio to solar/MW metallicity
;   (this is all for the COMPTON THIN regime)
;-------------------------------------------------------------------------------------
pro morrison_photoelec, frequency, sigma, ABUNDANCE=ABUNDANCE
  ; (can restructure this to allow vector i/o if the time for the run is too long)
  
  f_003keV = 7.253d15
  f_10keV  = 2.418d18

  keV_per_Hz = 4.13608d-18    
  fkeV = frequency * keV_per_Hz		; convert frequency from Hz to keV
    
  ; Now set the appropriate polynomial terms from Table 2 
  ;   of Morrison & McCammon (for a given frequency range)
    CASE 1 OF
      (fkeV LT 0.100): c = [	17.3,	608.1,	-2150.0	]
      (fkeV LT 0.284): c = [	34.6,	267.9,	 -476.1	]
      (fkeV LT 0.400): c = [	78.1,	 18.8,	    4.3	]
      (fkeV LT 0.532): c = [	71.4,	 66.8,	  -51.4	]
      (fkeV LT 0.707): c = [	95.5,	145.8,	  -61.1	]
      (fkeV LT 0.867): c = [   308.9,  -380.6,	  294.0	]
      (fkeV LT 1.303): c = [   120.6,	169.3,	  -47.7	]
      (fkeV LT 1.840): c = [   141.3,	146.8,	  -31.5	]
      (fkeV LT 2.471): c = [   202.7,	104.7,	  -17.0	]
      (fkeV LT 3.210): c = [   342.7,	 18.7,	    0.0	]
      (fkeV LT 4.038): c = [   352.2,	 18.7,	    0.0	]
      (fkeV LT 7.111): c = [   433.9,	 -2.4,	  -0.75	]
      (fkeV LT 8.331): c = [   629.0,	 30.9,      0.0	]
      (fkeV LT 10.00): c = [   701.2,	 25.2,	    0.0	]
      ELSE: c = [   701.2,	 25.2,	    0.0	]
    ENDCASE
    
  ; Use these coefficients to calculate the cross section
  ;   (in 10^-24 cm^2) per hydrogen atom
  sigma = (c[0] + c[1]*fkeV + c[2]*(fkeV*fkeV)) / (fkeV*fkeV*fkeV)  * 1.0d-24
  
  if (KEYWORD_SET(ABUNDANCE)) then sigma = sigma * ABUNDANCE
  
  return
end
