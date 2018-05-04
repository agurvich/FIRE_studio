;;; routine to return luminosities (in solar) in each of several bands 
;;    based on the Bruzual & Charlot 2003 stellar models
;;

function get_bc03_lum, metallicity, age_in_Gyr, VECTOR=VECTOR

	exec_call = return_idl_routines_homedir(0)+'/colors_bc03_code/colors'

	s_UBVRIJHK = fltarr(8)
	s_UBVRIJHK(0) = 5.66; //U
	s_UBVRIJHK(1) = 5.47; //B
	s_UBVRIJHK(2) = 4.82; //V
	s_UBVRIJHK(3) = 4.28; //R
	s_UBVRIJHK(4) = 3.94; //I
	s_UBVRIJHK(5) = 3.64; //J ?
	s_UBVRIJHK(6) = 3.44; //H ?
	s_UBVRIJHK(7) = 3.33; //K

	s_ugrizJHK = fltarr(8)
	s_ugrizJHK(0) = 6.2789; //u
	s_ugrizJHK(1) = 4.9489; //g
	s_ugrizJHK(2) = 4.44964; //r
	s_ugrizJHK(3) = 4.34644; //i
	s_ugrizJHK(4) = 4.3592; //z
	s_ugrizJHK(5) = 3.64; //J ?
	s_ugrizJHK(6) = 3.44; //H ?
	s_ugrizJHK(7) = 3.33; //K

	solar_mags = fltarr(14)
	solar_mags(0) = 4.74 ;bolometric
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

	
	l_stars = fltarr(14)	;luminosities
	l_stars_temp = fltarr(14,2)
	age_input = age_in_Gyr
	if (n_elements(age_in_Gyr) GT 1) then age_input = age_in_Gyr[0]
	met_input = metallicity
	if (n_elements(metallicity) GT 1) then met_input = metallicity[0]

    	n_temp = 1L
    	z_temp = 0.0E
    	m_temp = 1.0E
    	t_temp = float(age_input)
    	zm_temp= float(met_input)
    	
    	if (keyword_set(VECTOR)) then begin
    		n_temp = long(n_elements(metallicity))
    		m_temp = fltarr(n_elements(metallicity))+1.0E
    		t_temp = float(age_in_Gyr)
    		zm_temp= float(metallicity)
    		l_stars = fltarr(14,n_elements(metallicity))
    		l_stars_temp = fltarr(14,n_elements(metallicity))
    	endif
    	
		S = CALL_EXTERNAL(exec_call, $
			'colors', $
			n_temp, $
			z_temp, $
			m_temp, $
			t_temp, $
			zm_temp, $
			l_stars_temp)
	
		for i=0,13 do begin
			l_stars(i) = (10.0^(-0.4*(l_stars_temp(i,0)-solar_mags(i))))
		endfor
		
    	if (keyword_set(VECTOR)) then begin
		for i=0,13 do begin
			l_stars(i,*) = (10.0^(-0.4*(l_stars_temp(i,*)-solar_mags(i))))
		endfor
		endif	

	return, l_stars
end
