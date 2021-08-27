#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

double morrison_photoeletric_absorption(double x);
double morrison_photoeletric_absorption_nometals(double x);
double return_tau(double log_NH, double nu);
double attenuation(double nu, double NH, double metallicity_over_solar, int dust_key);
double pei_dust_extinction(double lambda_in_microns, int extinction_curve_key);
double g_bf(double x, double Te, int n);
double photoelectric_absorption(double nu, double metallicity_over_solar);
double compton_transmission(double nu_in_keV, double log_NH, double metallicity_in_solar);
double get_sx_atten_lookup(double NH, double z_in_solar);


int main(int argc, void *argv[])
{
  	if (argc < 5) {
  	  fprintf(stderr, "Expected >=5 arguments (found %d)\n",argc); 
  	  exit(0); 
  	}
	
	double nu, *nu_vec, NH, metallicity_over_solar, *NH_vec, *metal_vec, *atten_vec;
	int dust_key,N_nu,N_NH,i_nu,i_NH;
	
	N_nu = *(int *)argv[0];
	nu_vec = (double *)argv[1];
	N_NH = *(int *)argv[2];
	NH_vec = (double *)argv[3];
	metal_vec = (double *)argv[4];
	dust_key = *(int *)argv[5];
	atten_vec = (double *)argv[6];

	for (i_nu=0;i_nu<N_nu;i_nu++) {
		nu = nu_vec[i_nu];
		for (i_NH=0;i_NH<N_NH;i_NH++) {
			NH = NH_vec[i_NH];
			metallicity_over_solar=metal_vec[i_NH];
			atten_vec[i_nu*N_NH+i_NH]=attenuation(nu,NH,metallicity_over_solar,dust_key);
		}
	}
	return 1;
}

// shorthand LUT to get the SX-band attenuation as a function of NH and metallicity
//   -- note, this will depend on the shape of the soft X-ray spectrum (here HRH06), 
//   and for the marconi et al. 2004 spectrum weakly on L_bol (due to their interpolating
//   with alpha_ox all the way to 1keV, which is probably a bad approroximation
//   based on the Elvis et al. X-ray compilations
//
double get_sx_atten_lookup(double NH, double z_in_solar)
{
double log_nh[41]={
20.00, 20.10, 20.20, 20.30, 20.40, 20.50, 20.60, 20.70, 20.80, 20.90, 21.00, 21.10, 21.20, 21.30, 
21.40, 21.50, 21.60, 21.70, 21.80, 21.90, 22.00, 22.10, 22.20, 22.30, 22.40, 22.50, 22.60, 22.70, 
22.80, 22.90, 23.00, 23.10, 23.20, 23.30, 23.40, 23.50, 23.60, 23.70, 23.80, 23.90, 24.00};
double z_lut[7]={0.00,  0.10,  0.16,  0.25,  0.40,  0.63,  1.00};
double atten_0[41]={
0.00,  0.00,  0.01,  0.01, 0.01,  0.01,  0.01,  0.02,  0.02,  0.03,  0.04,  0.04,  0.05,  0.07, 
0.08,  0.10,  0.12, 0.14,  0.17,  0.20,  0.23,  0.27,  0.31,  0.35,  0.40,  0.45,  0.50,  0.57, 
0.64,  0.72, 0.80,  0.91,  1.02,  1.16,  1.32,  1.50,  1.72,  1.97,  2.28,  2.65,  3.10};
double atten_1[41]={
0.00, 0.01, 0.01,  0.01,  0.01,  0.01,  0.02,  0.02,  0.03,  0.04,  0.04,  0.05,  0.07,  0.08,
0.10, 0.12,  0.15,  0.18,  0.21,  0.25,  0.29,  0.33,  0.38,  0.44,  0.50,  0.57,  0.64, 
0.73, 0.83, 0.94,  1.07,  1.23,  1.41,  1.62,  1.87,  2.17,  2.54,  2.98,  3.52,  4.19,  5.00};
double atten_2[41]={
0.01,  0.01, 0.01,  0.01,  0.01,  0.02,  0.02,  0.03,  0.03,  0.04,  0.05,  0.06,  0.08, 0.09,
0.11,  0.14, 0.17,  0.20,  0.23,  0.27,  0.32,  0.37,  0.42,  0.49,  0.56,  0.63, 0.72, 
0.82,  0.93,  1.07, 1.22,  1.40,  1.62,  1.88,  2.18,  2.56,  3.01,  3.56,  4.24, 5.07,  6.09};
double atten_3[41]={
0.01,  0.01,  0.01, 0.01,  0.02,  0.02,  0.02,  0.03,  0.04,  0.05,  0.06, 0.07,  0.09,  0.11,
0.13,  0.16,  0.19, 0.23,  0.27,  0.32,  0.37,  0.42,  0.49,  0.56, 0.64,  0.73,  0.83, 
0.95,  1.09,  1.26,  1.45, 1.67,  1.94,  2.27,  2.66,  3.14,  3.73, 4.46,  5.35,  6.44,  7.80};
double atten_4[41]={
0.01,  0.01,  0.01,  0.01, 0.02,  0.02,  0.03,  0.04,  0.05, 0.06,  0.07,  0.09,  0.11,  0.13,
0.16,  0.19,  0.23,  0.27, 0.32,  0.38,  0.44,  0.50, 0.58,  0.67,  0.76,  0.87,  1.00, 
1.15,  1.33,  1.54,  1.78,  2.08, 2.44,  2.87,  3.40, 4.06,  4.86,  5.85,  7.07,  8.59, 10.47};
double atten_5[41]={
0.01,  0.01,  0.02,  0.02,  0.02, 0.03, 0.04, 0.05,  0.06,  0.07,  0.09,  0.11,  0.14,  0.17,
0.20,  0.24,  0.29,  0.34,  0.40, 0.46, 0.54, 0.62,  0.71,  0.82,  0.94,  1.08,  1.25,  1.45,  
1.68,  1.96,  2.29,  2.70,  3.19, 3.80, 4.55, 5.47,  6.61,  8.02,  9.76, 11.93, 14.63};  
double atten_6[41]={
0.01,  0.02,  0.02,  0.03,  0.03, 0.04, 0.05,  0.06,  0.08,  0.10,  0.12,  0.15,  0.18, 0.22, 
0.26,  0.31,  0.37,  0.43, 0.51,  0.59, 0.68,  0.78,  0.90,  1.04,  1.20,  1.39,  1.61, 1.87, 
2.19,  2.58,  3.05, 3.63,  4.34,  5.22, 6.31,  7.65,  9.31, 11.38, 13.95, 17.15, 21.15};

int n0 = (int )((log10(NH)-20.00)/0.10);
	if (n0 > 39) n0=39;
	if (n0 < 0)  n0=0;
	int np = n0+1;
double x_interp,f0,f1,f;
	x_interp = (log10(NH)-log_nh[n0])/0.10;

if (z_in_solar < 0.0) z_in_solar=0.;
// at low metallicity, scale linearly with metallicity
if (z_in_solar <= 0.1) {
	f0 = atten_0[n0] + (atten_0[np]-atten_0[n0])*x_interp;
	f1 = atten_1[n0] + (atten_1[np]-atten_1[n0])*x_interp;
	f = pow(10., f0 + (z_in_solar/0.1)*(f1-f0));
}
if ((z_in_solar > 0.10)&&(z_in_solar<=0.16)) {
	f0 = atten_1[n0] + (atten_1[np]-atten_1[n0])*x_interp;
	f1 = atten_2[n0] + (atten_2[np]-atten_2[n0])*x_interp;
	f = pow(10., f0 + ((z_in_solar-0.10)/(0.16-0.10))*(f1-f0));
}
if ((z_in_solar > 0.16)&&(z_in_solar<=0.25)) {
	f0 = atten_2[n0] + (atten_2[np]-atten_2[n0])*x_interp;
	f1 = atten_3[n0] + (atten_3[np]-atten_3[n0])*x_interp;
	f = pow(10., f0 + ((z_in_solar-0.16)/(0.25-0.16))*(f1-f0));
}
if ((z_in_solar > 0.25)&&(z_in_solar<=0.40)) {
	f0 = atten_3[n0] + (atten_3[np]-atten_3[n0])*x_interp;
	f1 = atten_4[n0] + (atten_4[np]-atten_4[n0])*x_interp;
	f = pow(10., f0 + ((z_in_solar-0.25)/(0.40-0.25))*(f1-f0));
}
if ((z_in_solar > 0.40)&&(z_in_solar<=0.63)) {
	f0 = atten_4[n0] + (atten_4[np]-atten_4[n0])*x_interp;
	f1 = atten_5[n0] + (atten_5[np]-atten_5[n0])*x_interp;
	f = pow(10., f0 + ((z_in_solar-0.40)/(0.63-0.40))*(f1-f0));
}
if ((z_in_solar > 0.63)&&(z_in_solar<=1.00)) {
	f0 = atten_5[n0] + (atten_5[np]-atten_5[n0])*x_interp;
	f1 = atten_6[n0] + (atten_6[np]-atten_6[n0])*x_interp;
	f = pow(10., f0 + ((z_in_solar-0.63)/(1.00-0.63))*(f1-f0));
}
if (z_in_solar > 1.0) {
	// in this regime, metals-dominated, can simply treat as an effective 
	//   (NH-dependent) cross-section
	// shift by log(z)
	n0 = (int )((log10(NH) + log10(z_in_solar) - 20.00)/0.10);
		if (n0 > 39) n0=39;
		if (n0 < 0)  n0=0;
		np = n0+1;
		x_interp = (log10(NH) + log10(z_in_solar) - log_nh[n0])/0.10;
		f = pow(10.,atten_6[n0] + (atten_6[np]-atten_6[n0])*x_interp);
}
f = 1./f;
if ((f>1.0)||(f<0.0)) f = 1.0;
return f;
}
// same, for hard X-rays; 
double get_hx_atten_lookup(double NH, double z_in_solar)
{
double log_nh[41]={
21.50, 21.60, 21.70, 21.80, 21.90, 22.00, 22.10, 22.20, 22.30, 22.40, 22.50, 22.60, 22.70,
22.80, 22.90, 23.00, 23.10, 23.20, 23.30, 23.40, 23.50, 23.60, 23.70, 23.80, 23.90, 24.00,
24.10, 24.20, 24.30, 24.40, 24.50, 24.60, 24.70, 24.80, 24.90, 25.00, 25.10, 25.20, 25.30,
25.40, 25.50};
double z_lut[7]={0.00,  0.10,  0.16,  0.25,  0.40,  0.63,  1.00};
double atten_0[41]={
0.001, 0.001, 0.002, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.009, 0.011, 0.014, 0.018, 0.023,
0.028, 0.035, 0.047, 0.059, 0.073, 0.089, 0.106, 0.142, 0.179, 0.219, 0.266, 0.315, 0.400, 0.490,
0.596, 0.715, 0.820, 1.086, 1.346, 1.650, 1.995, 2.346, 2.891, 3.566, 4.400, 5.431, 6.703};
double atten_1[41]={
0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.009, 0.012, 0.014, 0.018, 0.023, 0.028, 0.035, 0.044,
0.055, 0.068, 0.086, 0.107, 0.131, 0.158, 0.187, 0.238, 0.290, 0.348, 0.414, 0.483, 0.598, 0.720,
0.861, 1.020, 1.169, 1.512, 1.856, 2.250, 2.681, 3.110, 3.788, 4.620, 5.641, 6.897, 8.448};
double atten_2[41]={
0.003, 0.004, 0.005, 0.006, 0.008, 0.009, 0.012, 0.015, 0.019, 0.023, 0.029, 0.036, 0.045, 0.056,
0.069, 0.085, 0.108, 0.132, 0.161, 0.192, 0.226, 0.283, 0.342, 0.407, 0.481, 0.558, 0.687, 0.823,
0.981, 1.159, 1.330, 1.707, 2.091, 2.527, 3.003, 3.483, 4.233, 5.150, 6.277, 7.665, 9.384};
double atten_3[41]={
0.004, 0.005, 0.006, 0.008, 0.010, 0.013, 0.016, 0.020, 0.025, 0.031, 0.039, 0.048, 0.060, 0.074,
0.091, 0.111, 0.138, 0.168, 0.202, 0.240, 0.280, 0.345, 0.412, 0.487, 0.572, 0.662, 0.809, 0.967,
1.149, 1.355, 1.561, 1.986, 2.428, 2.927, 3.471, 4.034, 4.892, 5.943, 7.237, 8.836, 10.822};
double atten_4[41]={
0.006, 0.007, 0.009, 0.011, 0.014, 0.018, 0.022, 0.028, 0.035, 0.044, 0.054, 0.067, 0.082, 0.101,
0.123, 0.148, 0.182, 0.218, 0.259, 0.304, 0.352, 0.427, 0.506, 0.595, 0.696, 0.805, 0.978, 1.167,
1.386, 1.634, 1.894, 2.388, 2.914, 3.505, 4.153, 4.847, 5.873, 7.135, 8.694, 10.631, 13.044};
double atten_5[41]={
0.008, 0.011, 0.013, 0.017, 0.021, 0.026, 0.033, 0.040, 0.050, 0.062, 0.077, 0.094, 0.114, 0.139,
0.167, 0.199, 0.240, 0.284, 0.333, 0.387, 0.446, 0.536, 0.632, 0.741, 0.866, 1.004, 1.216, 1.451,
1.724, 2.036, 2.376, 2.970, 3.616, 4.344, 5.154, 6.057, 7.348, 8.946, 10.931, 13.406, 16.498};
double atten_6[41]={
0.013, 0.016, 0.020, 0.025, 0.031, 0.039, 0.048, 0.059, 0.073, 0.090, 0.109, 0.132, 0.160, 0.191,
0.226, 0.266, 0.316, 0.370, 0.430, 0.497, 0.572, 0.683, 0.804, 0.943, 1.104, 1.287, 1.556, 1.859,
2.213, 2.620, 3.077, 3.814, 4.637, 5.573, 6.644, 7.881, 9.597, 11.732, 14.396, 17.726, 21.893};


int n0 = (int )((log10(NH)-21.50)/0.10);
	if (n0 > 39) n0=39;
	if (n0 < 0)  n0=0;
	int np = n0+1;
double x_interp,f0,f1,f;
	x_interp = (log10(NH)-log_nh[n0])/0.10;

if (z_in_solar < 0.0) z_in_solar=0.;
// at low metallicity, scale linearly with metallicity
if (z_in_solar <= 0.1) {
	f0 = atten_0[n0] + (atten_0[np]-atten_0[n0])*x_interp;
	f1 = atten_1[n0] + (atten_1[np]-atten_1[n0])*x_interp;
	f = pow(10., f0 + (z_in_solar/0.1)*(f1-f0));
}
if ((z_in_solar > 0.10)&&(z_in_solar<=0.16)) {
	f0 = atten_1[n0] + (atten_1[np]-atten_1[n0])*x_interp;
	f1 = atten_2[n0] + (atten_2[np]-atten_2[n0])*x_interp;
	f = pow(10., f0 + ((z_in_solar-0.10)/(0.16-0.10))*(f1-f0));
}
if ((z_in_solar > 0.16)&&(z_in_solar<=0.25)) {
	f0 = atten_2[n0] + (atten_2[np]-atten_2[n0])*x_interp;
	f1 = atten_3[n0] + (atten_3[np]-atten_3[n0])*x_interp;
	f = pow(10., f0 + ((z_in_solar-0.16)/(0.25-0.16))*(f1-f0));
}
if ((z_in_solar > 0.25)&&(z_in_solar<=0.40)) {
	f0 = atten_3[n0] + (atten_3[np]-atten_3[n0])*x_interp;
	f1 = atten_4[n0] + (atten_4[np]-atten_4[n0])*x_interp;
	f = pow(10., f0 + ((z_in_solar-0.25)/(0.40-0.25))*(f1-f0));
}
if ((z_in_solar > 0.40)&&(z_in_solar<=0.63)) {
	f0 = atten_4[n0] + (atten_4[np]-atten_4[n0])*x_interp;
	f1 = atten_5[n0] + (atten_5[np]-atten_5[n0])*x_interp;
	f = pow(10., f0 + ((z_in_solar-0.40)/(0.63-0.40))*(f1-f0));
}
if ((z_in_solar > 0.63)&&(z_in_solar<=1.00)) {
	f0 = atten_5[n0] + (atten_5[np]-atten_5[n0])*x_interp;
	f1 = atten_6[n0] + (atten_6[np]-atten_6[n0])*x_interp;
	f = pow(10., f0 + ((z_in_solar-0.63)/(1.00-0.63))*(f1-f0));
}
if (z_in_solar > 1.0) {
	// in this regime, metals-dominated, can simply treat as an effective 
	//   (NH-dependent) cross-section
	// shift by log(z)
	n0 = (int )((log10(NH) + log10(z_in_solar) - 21.50)/0.10);
		if (n0 > 39) n0=39;
		if (n0 < 0)  n0=0;
		np = n0+1;
		x_interp = (log10(NH) + log10(z_in_solar) - log_nh[n0])/0.10;
		f = pow(10.,atten_6[n0] + (atten_6[np]-atten_6[n0])*x_interp);
}
f = 1./f;
if ((f>1.0)||(f<0.0)) f = 1.0;
return f;
}




//
// returns the attenuation for a given nu in Hz, column density in cm^-2
//   and metallicity in solar units. dust_key = 0 (MW), 1 (LMC), 2 (SMC)
//
double attenuation(double nu, double NH, double metallicity_over_solar, int dust_key)
{
	double sigma = 0.;
	double keV_in_Hz = 2.418e17;
	double c_light = 2.998e8;
	double micron  = 1.0e-6;
	double atten   = 1.;

	/* first, check for the usual band 'shortcuts' */
	if (nu==-1.0) return attenuation(c_light/(4.4e-7),NH,metallicity_over_solar,dust_key);
	if (nu==-2.0) return attenuation(c_light/(1.5e-5),NH,metallicity_over_solar,dust_key);
	if (nu==-3.0) return get_sx_atten_lookup(NH,metallicity_over_solar);
	if (nu==-4.0) return get_hx_atten_lookup(NH,metallicity_over_solar);


/*
  ; For optical-IR regions, we use the Pei numerical approximations below.
  ;
  ; xsi = tau(lambda)/tau(B) is the ratio of extinction at lambda to the 
  ;    extinction in the B-band. 
  ; k = 10^21 (tau_B / NH)   (NH in cm^2) gives the dimensionless gas-to-dust
  ;    ratio, with k=0.78 for MW, k=0.16 for LMC, k=0.08 for SMC.
  ;    k is INDEPENDENT of the grain properties, and seems to scale rougly
  ;    linearly with metallicity
  ; so, for now, assume solar metallicity and k = k_MW = 0.78. we can rescale later.
  ;
  ; tau_B = ( NH / (10^21 cm^-2) ) * k --> SIGMA_B = k*10^-21  cm^2
  ; tau_lambda = xsi * tau_B --> SIGMA = xsi * SIGMA_B
  ;
  ; k = 0.78 for the MW
  ; k = 0.08 for the SMC, approximately in line with the MW/LMC/SMC metallicity 
  ;  sequence, so we take a k_MW then scaled by the metallicity
*/
	double k_dust_to_gas = 0.78 * metallicity_over_solar;
	double lambda_microns = c_light / nu / micron;
	if (nu < 0.03*keV_in_Hz) 
		sigma += pei_dust_extinction(lambda_microns,dust_key) * k_dust_to_gas * 1.0e-21;


/*
  ; For 0.03 keV < E < 10 keV  
  ;   (7.2e15 < nu[Hz] < 2.4e18  or   1.2 < lambda[Angstroms] < 413)
  ;   we use the photoelectric absorption cross sections of 
  ;   Morrison & McCammon (1983)
  ;     NOTE: these assume solar abundances and no ionization, 
  ;             the appropriate number probably scales linearly with both
  ;   (this is all for the COMPTON THIN regime)
*/
	sigma += photoelectric_absorption(nu,metallicity_over_solar);

	atten *= exp(-NH*sigma);


/*
  ; Floor in cross-section set by non-relativistic (achromatic) Thompson scattering, 
  ;   but often higher (& get iron emission) -- 
  ;   calculated self-consistently for the induced ionization state of the 
  ;   gas from Matt, Pompilio, & La Franca; their extracted transmission curves 
  ;   applied here. The z scaling should be weak. Technically defines an attenuation 
  ;   factor, need to 
*/
	//sigma += 6.65e-25;
	atten *= compton_transmission(nu/keV_in_Hz,log10(NH),metallicity_over_solar);


return atten;
}







// pei et al. fitted extinction curves for MW, LMC, and SMC
//    extinction_curve_key = 0 (MW), 1 (LMC), 2 (SMC)
double pei_dust_extinction(double lambda_in_microns, int extinction_curve_key)
{
int MW_key=0;
int LMC_key=0;
int SMC_key=1;
	if (extinction_curve_key==0) {MW_key=1; LMC_key=0; SMC_key=0;}
	if (extinction_curve_key==1) {MW_key=0; LMC_key=1; SMC_key=0;}
	if (extinction_curve_key==2) {MW_key=0; LMC_key=0; SMC_key=1;}
int i;
double xsi = 0.0*lambda_in_microns;
if (MW_key==1) {
  double a[6] = {165., 14., 0.045, 0.002, 0.002, 0.012};
  double l[6] = {0.047, 0.08, 0.22, 9.7, 18., 25.};
  double b[6] = {90., 4.00, -1.95, -1.95, -1.80, 0.00};
  double n[6] = {2.0, 6.5, 2.0, 2.0, 2.0, 2.0};
  double R_V = 3.08;
  for(i=0;i<6;i++) xsi += a[i] / ( pow(lambda_in_microns/l[i],n[i]) + pow(l[i]/lambda_in_microns,n[i]) + b[i] );
}
if (LMC_key==1) {
  double a[6] = {175., 19., 0.023, 0.005, 0.006, 0.020};
  double l[6] = {0.046, 0.08, 0.22, 9.7, 18., 25.};
  double b[6] = {90., 5.50, -1.95, -1.95, -1.80, 0.00};
  double n[6] = {2.0, 4.5, 2.0, 2.0, 2.0, 2.0};
  double R_V = 3.16;
  for(i=0;i<6;i++) xsi += a[i] / ( pow(lambda_in_microns/l[i],n[i]) + pow(l[i]/lambda_in_microns,n[i]) + b[i] );
}
if (SMC_key==1) {
  double a[6] = {185., 27., 0.005, 0.010, 0.012, 0.030};
  double l[6] = {0.042, 0.08, 0.22, 9.7, 18., 25.};
  double b[6] = {90., 5.50, -1.95, -1.95, -1.80, 0.00};
  double n[6] = {2.0, 4.0, 2.0, 2.0, 2.0, 2.0};
  double R_V = 2.93;
  for(i=0;i<6;i++) xsi += a[i] / ( pow(lambda_in_microns/l[i],n[i]) + pow(l[i]/lambda_in_microns,n[i]) + b[i] );
}
//double R_lam = (1.0 + R_V) * xsi;
return xsi;
}

// uses the fully relativistic compton transmission curves calculated in 
//   Matt, Pompilio, & La Franca, extracted & MM83 component factored out. 
//   optionally, can include iron flourescence, although less well defined for 
//   low (but still non-zero) metallicities, in which case it should really be 
//   calculated self-consistently
//
double compton_transmission(double nu_in_keV, double log_NH, double metallicity_in_solar)
{
double sigma_thompson = 6.65e-25;
double floor = exp(-sigma_thompson * pow(10.0,log_NH));

if (nu_in_keV < 4.00) return floor;

// first do the no-metals attenuation 
double x_nu[20]={
4.00,   5.04,   6.34,   7.98,  10.05,  12.65,  15.92,  20.05,  25.24,  31.77, 
40.00,  50.36,  63.40,  79.81, 100.48, 126.49, 159.24, 200.47, 252.38, 317.73};
double y_1_23[20]={
0.9399, 0.9457, 0.9513, 0.9550, 0.9590, 0.9686, 0.9788, 0.9766, 0.9836, 0.9808,
0.9782, 0.9760, 0.9752, 0.9745, 0.9771, 0.9752, 0.9754, 0.9732, 0.9707, 0.9700};
double y_3_23[20]={
0.8191, 0.8274, 0.8469, 0.8574, 0.8724, 0.8838, 0.8966, 0.9042, 0.9136, 0.9150,
0.9181, 0.9176, 0.9145, 0.9193, 0.9216, 0.9228, 0.9191, 0.9140, 0.9117, 0.9124};
double y_1_24[20]={
0.5143, 0.5143, 0.5349, 0.5576, 0.6134, 0.6608, 0.6903, 0.7096, 0.7268, 0.7430,
0.7504, 0.7439, 0.7443, 0.7558, 0.7645, 0.7530, 0.7466, 0.7347, 0.7304, 0.7333};
double y_3_24[20]={
0.1360, 0.1360, 0.1395, 0.1584, 0.2231, 0.2991, 0.3902, 0.4454, 0.4809, 0.5012,
0.4977, 0.4735, 0.4530, 0.4546, 0.4530, 0.4300, 0.4146, 0.3980, 0.3945, 0.3946};
double y_1_25[20]={
0.0013, 0.0013, 0.0013, 0.0030, 0.0078, 0.0699, 0.1322, 0.1766, 0.1795, 0.1548,
0.1222, 0.0955, 0.0793, 0.0729, 0.0663, 0.0574, 0.0516, 0.0493, 0.0489, 0.0482};

int n0 = (int )(log10(nu_in_keV/4.00)/0.1);
	if (n0 > 18) n0=18;
	int n0_p = n0 + 1;
double xinterp = (nu_in_keV-x_nu[n0])/(x_nu[n0_p]-x_nu[n0]);
double z_1_23 = pow(10.,log10(y_1_23[n0]) + log10(y_1_23[n0_p]/y_1_23[n0])*xinterp);
double z_3_23 = pow(10.,log10(y_3_23[n0]) + log10(y_3_23[n0_p]/y_3_23[n0])*xinterp);
double z_1_24 = pow(10.,log10(y_1_24[n0]) + log10(y_1_24[n0_p]/y_1_24[n0])*xinterp);
double z_3_24 = pow(10.,log10(y_3_24[n0]) + log10(y_3_24[n0_p]/y_3_24[n0])*xinterp);
double z_1_25 = pow(10.,log10(y_1_25[n0]) + log10(y_1_25[n0_p]/y_1_25[n0])*xinterp);
double atten=0;
if (log_NH <= 23.0) atten = pow(z_1_23,pow(10.0,log_NH-23.0));
if ((log_NH>23.0)&&(log_NH<=23.5)) 
	atten = pow(10.,log10(z_1_23) + log10(z_3_23/z_1_23)*(log_NH-23.0)/(0.5));
if ((log_NH>23.5)&&(log_NH<=24.0)) 
	atten = pow(10.,log10(z_3_23) + log10(z_1_24/z_3_23)*(log_NH-23.5)/(0.5));
if ((log_NH>24.0)&&(log_NH<=24.5)) 
	atten = pow(10.,log10(z_1_24) + log10(z_3_24/z_1_24)*(log_NH-24.0)/(0.5));
if ((log_NH>24.5)&&(log_NH<=25.0))
	atten = pow(10.,log10(z_3_24) + log10(z_1_25/z_3_24)*(log_NH-24.5)/(0.5));
if ((log_NH>25.0))
	atten = pow(10.,log10(z_1_25) * pow(10.,log_NH-25.0));
if (atten <= floor) atten=floor;
if ((nu_in_keV >= x_nu[15])&&(atten>1.0)) atten = 1.0;


// now add the iron line emission -- scales linearly w. metallicity
if ((nu_in_keV > 5.80)&&(nu_in_keV < 7.40)) {
double atten_line = 0.0;
double x_nu[41]={
5.80, 5.84, 5.88, 5.92, 5.96, 6.00, 6.04, 6.08, 6.12, 6.16, 6.20, 6.24, 6.28, 6.32,
6.36, 6.40, 6.44, 6.48, 6.52, 6.56, 6.60, 6.64, 6.68, 6.72, 6.76, 6.80, 6.84, 6.88,
6.92, 6.96, 7.00, 7.04, 7.08, 7.12, 7.16, 7.20, 7.24, 7.28, 7.32, 7.36, 7.40};
double y_nu[41]={
0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00017,
0.00287, 0.00558, 0.00911, 0.01293, 0.04275, 0.20076, 0.38741, 0.23514, 0.08959, 0.00177,
0.00221, 0.00265, 0.00294, 0.00314, 0.00308, 0.00242, 0.00176, 0.00308, 0.00455, 0.00564,
0.00378, 0.00194, 0.02293, 0.05314, 0.06960, 0.03874, 0.00899, 0.00000, 0.00000, 0.00000,
0.00000};
// nearly constant fraction of continuum comes out the in the iron flourescence -- 
//   that's the better approximation to be using, as opposed to comparing with the 
//   scaling relative to the attenuation, which gives runaway results since the two 
//   are only weakly correlated
int n0 = (int )((nu_in_keV - 5.80)/0.04);
	if (n0 > 39) n0=39;
	int n0_p = n0 + 1;
double xinterp = (nu_in_keV-x_nu[n0])/(x_nu[n0_p]-x_nu[n0]);
atten_line = y_nu[n0] + (y_nu[n0_p]-y_nu[n0])*xinterp;
atten_line *= metallicity_in_solar;
if (atten_line > 1.0) atten_line = 1.0;
// continuum needs to be absorbed out for the line to appear (e.g. MM83)
if (atten >= atten_line) {
	return atten;
} else {
	return atten_line;	// switch to atten_line here if want to return the line:: 
					//   not clear that's really self-consistent (since implicit in MM83), 
					//   so for now turn this off in the code
}
}

return atten;
}



double photoelectric_absorption(double nu, double metallicity_over_solar)
{
	double nu_keV = 2.418e17;
	double nu_Rydberg = (2.998e8)/(912.0e-10);
	if (nu <= nu_Rydberg) return 0.0;
	double sigma_nometals = morrison_photoeletric_absorption_nometals(nu/nu_keV);
	return sigma_nometals + metallicity_over_solar * 
			(morrison_photoeletric_absorption(nu/nu_keV) - sigma_nometals);
}
double morrison_photoeletric_absorption_nometals(double x)
{
	//Morrison and McCammon 1983
	double x_keV = (2.418e17)/(3.0e14); // 1/keV_in_microns
	double c0, c1, c2;
		c0 =    17.3; c1 =   608.1; c2 = -2150.0;
	double s03 = (1.0e-24)*(c0+c1*0.03+c2*0.03*0.03)/(0.03*0.03*0.03);
	double s10 = (1.0e-24)*(c0+c1*0.10+c2*0.10*0.10)/(0.10*0.10*0.10);

	if (x<0.03) 			
		return s03*pow(x/0.03,-2.43);
	if ((x>=0.03)&&(x<=0.1))	
		return (1.0e-24)*(c0+c1*x+c2*x*x)/(x*x*x);
	if (x>0.1)
		return s10*g_bf(x*x_keV,0.0,1)*pow(x/0.1,-3.0)/g_bf(0.1*x_keV,0.0,1);
}
double g_bf(double x, double Te, int n)
{
	double x_g;
	double nu_g = 3.28989e15;//Hz
	double c = 3.0e14;//microns / sec
	double u_n = 0.0;
	double n_d = 0.0;

	n_d = ((double) n);
	x_g = nu_g/c;
	u_n = n_d*n_d*(x/x_g);

	//from allen's astrophysical quantities
	return 1.0 + (0.1728*(u_n-1.0)/(n_d*pow(u_n+1.0,2.0/3.0))) - (0.0496*(u_n*u_n + (4.0/3.0)*u_n + 1.0)/(n_d*pow(u_n+1.0,4.0/3.0)));
}
double morrison_photoeletric_absorption(double x)	// x is nu in keV
{    
//  ; set the appropriate polynomial terms from Table 2 
//  ;   of Morrison & McCammon 1983 (for a given frequency range)
	double c0, c1, c2;
	if(x<0.03)
	{
		c0 =    17.3;
		c1 =   608.1;
		c2 = -2150.0;
		return (1.0e-24)*(c0+c1*0.03+c2*0.03*0.03)/(0.03*0.03*0.03)*pow(x/0.03,-2.43);
	}
	else if((x>=0.03)&&(x<0.1))
	{
		c0 =    17.3;
		c1 =   608.1;
		c2 = -2150.0;

	}else if((x>=0.1)&&(x<0.284))
	{
		c0 =    34.6;
		c1 =   267.9;
		c2 =  -476.1;

	}else if((x>=0.284)&&(x<0.4))
	{
		c0 =    78.1;
		c1 =    18.8;
		c2 =     4.3;

	}else if((x>=0.4)&&(x<0.532))
	{
		c0 =    71.4;
		c1 =    66.8;
		c2 =   -51.4;

	}else if((x>=0.532)&&(x<0.707))
	{
		c0 =    95.5;
		c1 =   145.8;
		c2 =   -61.1;

	}else if((x>=0.707)&&(x<0.867))
	{
		c0 =   308.9;
		c1 =  -380.6;
		c2 =   294.0;

	}else if((x>=0.867)&&(x<1.303))
	{
		c0 =   120.6;
		c1 =   169.3;
		c2 =   -47.7;

	}else if((x>=1.303)&&(x<1.840))
	{
		c0 =   141.3;
		c1 =   146.8;
		c2 =   -31.5;

	}else if((x>=1.840)&&(x<2.471))
	{
		c0 =   202.7;
		c1 =   104.7;
		c2 =   -17.0;

	}else if((x>=2.471)&&(x<3.210))
	{
		c0 =   342.7;
		c1 =    18.7;
		c2 =     0.0;

	}else if((x>=3.210)&&(x<4.038))
	{
		c0 =   352.2;
		c1 =    18.7;
		c2 =     0.0;

	}else if((x>=4.038)&&(x<7.111))
	{
		c0 =   433.9;
		c1 =    -2.4;
		c2 =     0.75;

	}else if((x>=7.111)&&(x<8.331))
	{
		c0 =   629.0;
		c1 =    30.9;
		c2 =     0.0;

	}else if((x>=8.331)&&(x<10.00))
	{
		c0 =   701.2;
		c1 =    25.2;
		c2 =     0.0;

	}else{	
		// extrapolate the > 10 keV results to higher frequencies
		c0 =   701.2;
		c1 =    25.2;
		c2 =     0.0;
	}
// Use these coefficients to calculate the cross section per hydrogen atom
	return (1.0e-24)*(c0+c1*x+c2*x*x)/(x*x*x); //cm^2
}
