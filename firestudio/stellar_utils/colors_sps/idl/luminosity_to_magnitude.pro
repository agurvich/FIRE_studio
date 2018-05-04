;;
;;
;; function to take the input luminosity of a given band and return the 
;;    corresponding magnitude, for the appropriate bands and units as 
;;    listed below
;;
;; UNITS_... keywords tell it the units of L, whether solar L in the 
;;    band, bolometric L_sun (nuLnu/3.9d33), cgs (nuLnu/[erg/s])
;;		default : bolometric L_sun
;;
;; BAND_... keywords tell it the relevant band, of UBVRIJHKugriz (SDSS 
;;	  un-primed filters), or /BOLOMETRIC, or BAND_NUMBER uses the number of the 
;;    band below instead of making you explicitly name it
;;		default : bolometric
;;
;; VEGA or AB keywords have it return either vega or ab magnitudes
;;		default : VEGA for UBVRIJHK (Johnsons UBVRI, Cousins JHK), AB for ugriz  
;;
;; L_NU or NU_L_NU tell it whether the input is in specific luminosity L_NU (e.g. 
;;   erg/s/Hz), or NU_L_NU (e.g. erg/s)
;;		default : NU_L_NU
;;
;;
;;
function luminosity_to_magnitude, L, $
	UNITS_SOLAR_BOL=UNITS_SOLAR_BOL, UNITS_SOLAR_BAND=UNITS_SOLAR_BAND, $
	UNITS_CGS=UNITS_CGS, $
	NU_L_NU=NU_L_NU, L_NU=L_NU, $
	BAND_U=BAND_U,BAND_B=BAND_B,BAND_V=BAND_V,BAND_R=BAND_R,BAND_I=BAND_I,$
	BAND_J=BAND_J,BAND_H=BAND_H,BAND_K=BAND_K,BAND_SDSS_u=BAND_SDSS_u, $
	BAND_SDSS_g=BAND_SDSS_g,BAND_SDSS_r=BAND_SDSS_r,BAND_SDSS_i=BAND_SDSS_i, $
	BAND_SDSS_z=BAND_SDSS_z,BOLOMETRIC=BOLOMETRIC, BAND_NUMBER=BAND_NUMBER, $
	VEGA=VEGA, AB=AB

;Luminosity index legend
; 0 = bolometric luminosity
; 1 = Johnsons U
; 2 = Johnsons B
; 3 = Johnsons V
; 4 = Johnsons R
; 5 = Johnsons I
; 6 = Cousins J
; 7 = Cousins H
; 8 = Cousins K
; 9 = Sloan u
;10 = Sloan g
;11 = Sloan r
;12 = Sloan i
;13 = Sloan z

N_BANDS = 14

;VEGA system
;from www.ucolick.org/~cnaw/sun.html
;following Fukugita et al. 1995, PASP, 105, 945
mag_sun_vega = fltarr(N_BANDS)
mag_sun_vega(0) = 4.74;  bolometric from Allen's Astrophysical Quantities
mag_sun_vega(1) = 5.56;  //U (BESSEL)
mag_sun_vega(2) = 5.45;  //B (BESSEL)
mag_sun_vega(3) = 4.80;  //V (BESSEL)
mag_sun_vega(4) = 4.46;  //R (KPNO)
mag_sun_vega(5) = 4.10;  //I (KPNO)
mag_sun_vega(6) = 3.66;  //J (BESSEL)
mag_sun_vega(7) = 3.32;  //H (BESSEL)
mag_sun_vega(8) = 3.28;  //K (BESSEL)
mag_sun_vega(9) = 5.82;  //SDSS u (unprimed Vega)
mag_sun_vega(10) = 5.44; //SDSS g (unprimed Vega)
mag_sun_vega(11) = 4.52; //SDSS r (unprimed Vega)
mag_sun_vega(12) = 4.11; //SDSS i (unprimed Vega)
mag_sun_vega(13) = 3.89; //SDSS z (unprimed Vega)


;AB system
mag_sun_ab = fltarr(N_BANDS)
mag_sun_ab(0) = 4.74;  
	l_bol_sun = 3.9d33; // bolometric solar in erg/s
mag_sun_ab(1) = 6.34;  //U (BESSEL)
mag_sun_ab(2) = 5.33;  //B (BESSEL)
mag_sun_ab(3) = 4.81;  //V (BESSEL)
mag_sun_ab(4) = 4.65;  //R (KPNO)
mag_sun_ab(5) = 4.55;  //I (KPNO)
mag_sun_ab(6) = 4.57;  //J (BESSEL)
mag_sun_ab(7) = 4.71;  //H (BESSEL)
mag_sun_ab(8) = 5.19;  //K (BESSEL)
mag_sun_ab(9) = 6.75;  //SDSS u (unprimed AB)
mag_sun_ab(10) = 5.33; //SDSS g (unprimed AB)
mag_sun_ab(11) = 4.67; //SDSS r (unprimed AB)
mag_sun_ab(12) = 4.48; //SDSS i (unprimed AB)
mag_sun_ab(13) = 4.42; //SDSS z (unprimed AB)


;Effective wavelengths of the bands (in Angstroms), to compute nuLnu<->Lnu
; UBVRIJHK from http://cassfos02.ucsd.edu/physics/ph162/mags.html
; SDSS ugriz from http://www.sdss.org/dr4/instruments/imager/index.html#filters
lambda_eff = fltarr(N_BANDS)
lambda_eff(0) = 1.0;  bolometric, no nu
lambda_eff(1) = 3600.0;  //U
lambda_eff(2) = 4400.0;  //B
lambda_eff(3) = 5556.0;  //V
lambda_eff(4) = 6940.0;  //R
lambda_eff(5) = 8700.0;  //I
lambda_eff(6) = 12150.;  //J
lambda_eff(7) = 16540.;  //H
lambda_eff(8) = 21790.;  //K
lambda_eff(9)  = 3551.;  //SDSS u
lambda_eff(10) = 4686.;  //SDSS g
lambda_eff(11) = 6165.;  //SDSS r
lambda_eff(12) = 7481.;  //SDSS i
lambda_eff(13) = 8931.;  //SDSS z
c_light = 2.998d10; // speed of light in cm/s
nu_eff  = c_light/(lambda_eff * 1.0d-8); // converts to nu_eff in Hz


i_BAND = 0
	if (keyword_set(BAND_NUMBER)) then i_BAND=BAND_NUMBER
	if (keyword_set(BAND_U))  then i_BAND=1
	if (keyword_set(BAND_B))  then i_BAND=2
	if (keyword_set(BAND_V))  then i_BAND=3
	if (keyword_set(BAND_R))  then i_BAND=4
	if (keyword_set(BAND_I))  then i_BAND=5
	if (keyword_set(BAND_J))  then i_BAND=6
	if (keyword_set(BAND_H))  then i_BAND=7
	if (keyword_set(BAND_K))  then i_BAND=8
	if (keyword_set(BAND_SDSS_u)) then i_BAND=9
	if (keyword_set(BAND_SDSS_g)) then i_BAND=10
	if (keyword_set(BAND_SDSS_r)) then i_BAND=11
	if (keyword_set(BAND_SDSS_i)) then i_BAND=12
	if (keyword_set(BAND_SDSS_z)) then i_BAND=13
	if (keyword_set(BOLOMETRIC))  then i_BAND=0
	

;; default to Vega for bolometric & UBVRIJHK, and AB for ugriz
vega_key = 1
if ((i_BAND GT 8) AND (i_BAND LE 13)) then vega_key = 0
if (keyword_set(VEGA)) then vega_key=1
if (keyword_set(AB)) then vega_key=0
magnitude_zero_point = mag_sun_vega(i_BAND)
if (vega_key EQ 0) then magnitude_zero_point = mag_sun_ab(i_BAND)


;; use the AB magnitudes to convert to an actual L_nu of the sun in each band
lnu_sun_band = fltarr(N_BANDS)
ten_pc   = 10.0 * 3.086d18; //10 pc in cm
log_S_nu = -(mag_sun_ab + 48.6)/2.5;	//zero point definition for ab magnitudes
S_nu     = 10^(DOUBLE(log_S_nu)); // get the S_nu at 10 pc which defines M_AB
lnu_sun_band = S_nu * (4.0*!PI*ten_pc*ten_pc); // multiply by distance modulus 
nulnu_sun_band = lnu_sun_band * nu_eff; // multiply by nu_eff to get nu*L_nu
;; correct the bolometric
lnu_sun_band(0) = l_bol_sun;
nulnu_sun_band(0) = l_bol_sun;


;; alright, now have lnu of the sun in each band (the appropriate normalization
;;   for either magnitude system), can compare with the input luminosity
;;
nulnu_given = L;
	if (keyword_set(NU_L_NU)) then nulnu_given = L;
	if (keyword_set(L_NU))    then nulnu_given = nu_eff(i_BAND) * L;


;;default to assume in units of solar bolometric (if nu*L_nu):  
	l_in_solar_in_band = nulnu_given*(l_bol_sun/nulnu_sun_band(i_BAND));
;; or L_nu(sun) in the band (if given L_nu):
	if (keyword_set(L_NU)) then l_in_solar_in_band = L;

;; convert to the appropriate units
if (keyword_set(UNITS_SOLAR_BAND)) then l_in_solar_in_band = nulnu_given; //given in solar in band
	if (keyword_set(UNITS_SOLAR_BAND) AND keyword_set(L_NU)) then l_in_solar_in_band = L;
if (keyword_set(UNITS_SOLAR_BOL))  then l_in_solar_in_band = nulnu_given * (l_bol_sun/nulnu_sun_band(i_BAND));
if (keyword_set(UNITS_CGS))        then l_in_solar_in_band = nulnu_given / nulnu_sun_band(i_BAND);


return, magnitude_zero_point - 2.5*alog10(l_in_solar_in_band);

end
