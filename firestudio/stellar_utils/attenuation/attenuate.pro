;; returns the attenuation at a given frequency nu (in Hz) for a 
;;    given column density (log(NH/cm^-2)), and metallicity (in solar units)
;;
;;		i.e. for some luminosity L, observed = L * attenuate(nu,alog10(NH),Z/0.02)
;;
;;
;;	flags: can change the reddening curve adopted 
;;		SMC (default) = smc-like reddening curve (strong UV & mid-IR extinction)
;;			observations find this most appropriate for quasars: Hopkins et al. 2004
;;		MW  = MW-like (weaker UV & mid-IR; strong 2200 angstrom bump :: 
;;		LMC = LMC-like (basically halfway between the two)
;;			(gas to dust ratio in each case is MW-like average scaled 
;;				linearly with metallicity)
;;
;;	can call for specific bands (i.e. net attenuation integrated over the band):
;;		BB = B-band, IR = mid-IR (15 microns), 
;;		SX = soft X-ray (0.5-2 keV), 
;;		HX = hard X-ray (2-10 keV)
;;
;;  the dust extinction curves are from Pei et al. 1992. gas-to-dust ratios 
;;	  from Bouchet et al. 1985. 
;;
;;  above the Lyman edge, photoionization cross sections computed following 
;;   Morrison & McCammon 1983. the code follows Brant's decomposition of this into 
;;   a metals-free component and a metals component that scales linearly with 
;;   the input metallicity
;;
;;  compton scattering is treated as non-relativistic achromatic Thompson scattering
;;	 at low energies (<~ 4 keV), and above this the transmission curves from 
;;   Matt, Pompilio, & La Franca 1995 have been extracted to interpolate over a wide 
;;   range of frequencies and NH. These are decomposed into a H-dominated continuum 
;;   component that basically scales in a metal-free manner, and an iron flourescence 
;;   component (between ~5.8-7.4 keV) that emits in a manner scaling with metallicity. 
;;   the iron flourescence is fine for most applications, but can give strange 
;;   (although fixed s.t. non-divergent and still monotonic) behavior when you're  
;;   simultaneously at very low but non-zero metallicities (<~ 0.1 solar) 
;;   and very high (logNH > 25-25.5) column densities -- in this regime the 
;;   flourescence should really be calculated self-consistently. 
;;
;;
;;
function attenuate, nu_in_Hz, log_NH, metallicity_in_solar, $
	SMC=SMC, LMC=LMC, MW=MW, $
	BB=BB, IR=IR, SX=SX, HX=HX


	exec_call=return_idl_routines_homedir(0)+'/attenuation/attenuate.so'

	
	;; default to SMC-like reddening
	dust_key = 2L
		if (keyword_set(MW))  then dust_key = 0L
		if (keyword_set(LMC)) then dust_key = 1L
		if (keyword_set(SMC)) then dust_key = 2L


	nu_in_Hz     = DOUBLE(nu_in_Hz)
		if (keyword_set(BB)) then nu_in_Hz = -1.0d0
		if (keyword_set(IR)) then nu_in_Hz = -2.0d0
		if (keyword_set(SX)) then nu_in_Hz = -3.0d0
		if (keyword_set(HX)) then nu_in_Hz = -4.0d0	
		N_nu  = LONG(n_elements(nu_in_Hz))


	NH = 10^(DOUBLE(log_NH))
		N_NH = LONG(n_elements(NH))
	metallicity_in_solar=DOUBLE(metallicity_in_solar)
		N_metal=LONG(n_elements(metallicity_in_solar))
		if (N_metal LE 1) then metallicity_in_solar=0.0d0*NH + DOUBLE(metallicity_in_solar(0))

	atten = DOUBLE(fltarr(N_nu*N_NH))
	S = CALL_EXTERNAL(exec_call, $
           'main',N_nu,nu_in_Hz,N_NH,NH,metallicity_in_solar,dust_key,atten)

	atten_f = DOUBLE(fltarr(N_nu,N_NH))
	for i=0,N_nu-1 do begin
		atten_f[i,*]=atten[i*N_NH:(i+1)*N_NH-1]
	endfor
	atten=atten_f

	bad = where((FINITE(atten) NE 1) OR (FINITE(atten,/NAN) EQ 1) OR (atten EQ 0.),n_bad)
	if (n_bad GT 0) then atten[bad] = 1.0d-40

;;	for i=0,n_lum-1 do begin
;;		plot,alog10(nu_in_Hz),l_band_all[*,i]-log_l_bol[i],ystyle=1,yrange=[-4.,0.]
;;	endfor


	return, atten
end

