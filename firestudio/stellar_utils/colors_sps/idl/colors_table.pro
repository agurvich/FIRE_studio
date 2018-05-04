function colors_table, age_in_Gyr, metallicity_in_solar_units, BAND_ID=BAND_ID, $
  SALPETER_IMF=SALPETER_IMF,CHABRIER_IMF=CHABRIER_IMF, QUIET=QUIET, CRUDE=CRUDE, $
  RETURN_NU_EFF=RETURN_NU_EFF,RETURN_LAMBDA_EFF=RETURN_LAMBDA_EFF, $
  UNITS_SOLAR_IN_BAND=UNITS_SOLAR_IN_BAND


  ;; default = bolometric 
  band = 0 ; bolometric
   if (keyword_set(BAND_ID)) then band=BAND_ID
   j = [  0,  6,  7,  8,  9, 10, 11, 12, 13,  1,   2,   3,   4,   5] ;; ordering I'm used to
   i = [  0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10,  11,  12,  13] ;; ordering of this
     band_standardordering = band
     band = j[band]
   if (band GT 13) then begin
	print,'BAND_ID must be set < 13'
	return, 0
   endif
   b = ['Bolometric','Sloan u','Sloan g','Sloan r','Sloan i','Sloan z','Johnsons U','Johnsons B', $
	'Johnsons V','Johnsons R','Johnsons I','Cousins J','Cousins H','Cousins K']
   if keyword_set(RETURN_NU_EFF) or keyword_set(RETURN_LAM_EFF) then begin
   nu_eff=[0.,10^[		14.92,	14.79,		14.678,	14.590,		14.527,		14.93,		14.84,	$
   		14.732,		14.635,		14.537,		14.39,		14.26,		14.14	]]
   lam_eff=[0., 		3541.,	4653.,		6147.,	7461.,		8904.,	3600.,		4400.,	$
   		5556.,		6940.,		8700.,		12150.,		16540.,		21790.]
   nu_eff=2.998d18/lam_eff
   endif
   if (keyword_set(QUIET) eq 0) then print, 'Calculating M/L in ',b[band],' (',band,',',band_standardordering,')'

   froot = return_idl_routines_homedir(1)+'/colors_bc03_code/' ;; directory in which the data binaries are stored
    fnm_c= 'colors.chabrier.dat'
    fnm_s= 'colors.salpeter.dat'
   ;; default to Chabrier IMF
   fname = froot+fnm_c
     if (keyword_set(CHABRIER_IMF)) then fname = froot+fnm_c
     if (keyword_set(SALPETER_IMF)) then fname = froot+fnm_s
   
   OPENR,1,fname,/SWAP_IF_BIG_ENDIAN
   Nl = READ_BINARY(1,DATA_DIMS=0,DATA_TYPE=3)
   Na = READ_BINARY(1,DATA_DIMS=0,DATA_TYPE=3)
   Nz = READ_BINARY(1,DATA_DIMS=0,DATA_TYPE=3)
   z_grid   = READ_BINARY(1,DATA_DIMS=Nz,DATA_TYPE=5)
   age_grid = READ_BINARY(1,DATA_DIMS=Na,DATA_TYPE=5)
   l_all    = READ_BINARY(1,DATA_DIMS=[Nl,Na,Nz],DATA_TYPE=5)
   CLOSE,1

  l_band = dblarr(Na,Nz)
  for iz=0,Nz-1 do l_band(*,iz)=l_all(band,*,iz)

  ;; allow for extreme metallicities (i.e. extrapolate 
  ;;   linearly past limit above)
  PUSHMETALS = 1
  if PUSHMETALS then begin
   Nz = Nz + 1
   z_ext  = 1.0
   z_grid = [z_grid, z_ext]
    lb1 = l_band(*,Nz-3)
    lb2 = l_band(*,Nz-2)
    lbx = dblarr(Na,Nz)
    lbx[*,0:Nz-2] = l_band
    lbx[*,Nz-1] = (lb2 - lb1)/(alog10(z_grid[Nz-2]/z_grid[Nz-3]))*alog10(z_grid[Nz-1]/z_grid[Nz-2])
    l_band = lbx
  endif
  
  ;; get the x-axis (age) locations of input points
    iage = findgen(Na)
    ia_pts = INTERPOL(iage,age_grid,alog10(age_in_Gyr)+9.0)
     bad = where(ia_pts GT Na, bad_i)
     if (bad_i GT 0) then ia_pts[bad] = Na
  ;; get the y-axis (metallicity) locations of input points
    imet = findgen(Nz)
    zsun   = 0.02
    iz_pts = INTERPOL(imet,alog10(z_grid),alog10(metallicity_in_solar_units*zsun))
     bad = where(iz_pts GT Nz, bad_i)
     if (bad_i GT 0) then iz_pts[bad] = Nz

	if (keyword_set(CRUDE)) then begin
	  ia_pts=round(ia_pts) & iz_pts=round(iz_pts)	
	  l_b=l_band[ia_pts,iz_pts]
	endif else begin
	  l_b = INTERPOLATE(l_band,ia_pts,iz_pts)
	endelse
	
	if keyword_set(RETURN_NU_EFF) then RETURN_NU_EFF=nu_eff[band]
	if keyword_set(RETURN_LAMBDA_EFF) then RETURN_NU_EFF=lam_eff[band]
        l_b=10.d0^(DOUBLE(l_b))
 
  ;; output is currently L/M in L_sun_IN_THE_BAND_OF_INTEREST/M_sun, 
   ;; but we want our default to be L/M in units of L_bolometric/M_sun = 3.9e33/2.0e33, so 
   ;;   need to get rid fo the L_sun_IN_THE_BAND_OF_INTEREST/L_bolometric


;AB system solar luminosities used for determining L_sun in absolute units for each of these
N_BANDS=14
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
lambda_eff(0) = 4243.93;  bolometric, no nu
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

  ten_pc   = 10.d0 * 3.086d18; //10 pc in cm
  log_S_nu = -(mag_sun_ab + 48.6d0)/2.5d0; //zero point definition for ab magnitudes
  S_nu     = 10.d0^(DOUBLE(log_S_nu)); // get the S_nu at 10 pc which defines M_AB
  lnu_sun_band = S_nu * (4.d0*!PI*ten_pc*ten_pc); // multiply by distance modulus 
  nulnu_sun_band = lnu_sun_band * nu_eff; // multiply by nu_eff to get nu*L_nu
  l_bol_sun=nulnu_sun_band[0]

  if not keyword_set(UNITS_SOLAR_IN_BAND) then $
  l_b = l_b * nulnu_sun_band[band_standardordering]/l_bol_sun  

  return, l_b
end

