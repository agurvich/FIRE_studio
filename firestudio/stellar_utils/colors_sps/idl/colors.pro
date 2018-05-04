function colors,N,mass,age,metallicity,model


luminosities = fltarr(14,N)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;	colors.pro
;
;	TAKES:
;	
;	mass        in Solar Masses (no h)
;	age         in Gyr          (no h)
;	metallicity in Solar units  (Z_sun = 0.02)
;	model:
;		0 == Salpeter IMF, 0.1, 100.0
;		1 == Chabrier IMF, 0.1, 100.0
;		default = Salpeter
;
;
;	USES:
;		Padova 1994 stellar tracks
;
;
;	RETURNS:
;
;		Luminosities in Solar Luminosities (no h)
;		NOTE: returns solar luminosities in the chosen
;		band!
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

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

;VEGA system
;from www.ucolick.org/~cnaw/sun.html
;following Fukugita et al. 1995, PASP, 105, 945
s_UBVRIJHK = fltarr(8)
s_UBVRIJHK(0) = 5.56;  //U (BESSEL)
s_UBVRIJHK(1) = 5.45;  //B (BESSEL)
s_UBVRIJHK(2) = 4.80;  //V (BESSEL)
s_UBVRIJHK(3) = 4.46;  //R (KPNO)
s_UBVRIJHK(4) = 4.10;  //I (KPNO)
s_UBVRIJHK(5) = 3.66;  //J (BESSEL)
s_UBVRIJHK(6) = 3.32;  //H (BESSEL)
s_UBVRIJHK(7) = 3.28;  //K (BESSEL)

;AB  magnitudes
;from www.ucolick.org/~cnaw/sun.html
;following Fukugita et al. 1995, PASP, 105, 945
s_ugrizJHK = fltarr(5)
s_ugrizJHK(0) = 6.75; //u SDSS
s_ugrizJHK(1) = 5.33; //g SDSS
s_ugrizJHK(2) = 4.67; //r SDSS
s_ugrizJHK(3) = 4.48; //i SDSS
s_ugrizJHK(4) = 4.42; //z SDSS

solar_mags = fltarr(14)
solar_mags(0) = 4.74 ;bolometric from Allen's Astrophysical Quantities, p.
solar_mags(1) = s_UBVRIJHK(0)
solar_mags(2) = s_UBVRIJHK(1)
solar_mags(3) = s_UBVRIJHK(2)
solar_mags(4) = s_UBVRIJHK(3)
solar_mags(5) = s_UBVRIJHK(4)
solar_mags(6) = s_UBVRIJHK(5)
solar_mags(7) = s_UBVRIJHK(6)
solar_mags(8) = s_UBVRIJHK(7)
solar_mags(9) = s_ugrizJHK(0)
solar_mags(10) = s_ugrizJHK(1)
solar_mags(11) = s_ugrizJHK(2)
solar_mags(12) = s_ugrizJHK(3)
solar_mags(13) = s_ugrizJHK(4)

fname=return_idl_c_dir()+'/colors/colors'
S = CALL_EXTERNAL(fname, $
	'main', $
	long(N), $
	age, $           ;age lin gyr
	metallicity, $   ;metallicity in solar units
	luminosities, $  ;luminosities in abs magnitudes
	long(model), $	 ;model==0 -> Salpeter, model==1 Chabrier
	/F_VALUE)

for i=0,13 do luminosities(i,*) = (10.0^(-0.4*(luminosities(i,*)-solar_mags(i))))*mass(*)

return,luminosities ;in Solar Luminosities in each band (no factors of h)

end
