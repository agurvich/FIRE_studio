;; routine to return the solar absolute magnitude in each band 
;;    that the colors code gives 

function get_solar_mags

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

	return, solar_mags
end
