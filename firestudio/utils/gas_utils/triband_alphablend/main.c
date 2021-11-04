#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* 
    this routine takes a set of points and does an approximate projection through them, 
        weighting the quantity by the weight. Returns sum(weight) in pixel to OUT0 and 
        sum(weight*quantity) in pixel to OUT1.
      
*/

/* extremely fast approximation function for the exponential, 
    useful here since fractional accuracy errors are smaller than the kernel sources anyways */
inline double fast_exp(double y) {
    double d;
    *((int*)(&d) + 0) = 0;
    *((int*)(&d) + 1) = (int)(1512775 * y + 1072632447);
    return d;
}

void neighborloop(
  long imin, long imax, // row indices
  long jmin, long jmax, // column indices
  float * x, float * y, // particle positions
  double * x_i, double * y_j, // cell positions
  // pass-through scalars
  long n, // particle index
  int use_kernel,
  double h, double h2, double h2_i,
  double dx_n, double dy_n,
  double * Kernel,
  double kernel_spacing_inv,
  long Ypixels,
  // output
  double * wt_sum,
  float * R, float * G, float * B, float * alphas,
  float * OUT0, float * OUT1, float * OUT2){

  double x2_n, r2_n, wk;
  long i,j,k;
  int accumulate;
  accumulate = (*wt_sum > 0);

  // avoid re-computing background alpha inside the loop
  float bkg_alpha=(1-alphas[n]);

  // loop through all pixels this particle covers
  for(i=imin;i<imax;i++){
    dx_n = x[n]-x_i[i]; 
    if (fabs(dx_n) < h){
      x2_n = dx_n*dx_n;
      // loop through all pixels this particle covers
      for(j=jmin;j<jmax;j++){
        dy_n = y[n]-y_j[j]; 
        r2_n = x2_n + dy_n*dy_n; 
        if (r2_n < h2){
          // use the kernel lookup table compiled above for the weighting //
          r2_n *= h2_i*kernel_spacing_inv; // ok now have the separation in units of hsml, then kernel_table // 
          k = (long)r2_n;

          wk = 1; // initialize wk assuming we will not use the kernel
          if (use_kernel) wk = h2_i * (Kernel[k] + (Kernel[k+1]-Kernel[k])*(r2_n-k)); // ok that's the weighted result

          k = j + Ypixels*i; // j runs 0-Ypixels-1, so this provides the necessary indexing //

          if (!accumulate) *wt_sum += wk;
          else{
            // renormalize by sum of wk
            wk = wk / *wt_sum;

            // pre-multiply alpha blend in 3 bands
            OUT0[k] = R[n]*wk + bkg_alpha*OUT0[k];
            OUT1[k] = G[n]*wk + bkg_alpha*OUT1[k];
            OUT2[k] = B[n]*wk + bkg_alpha*OUT2[k];

          }// if accumulate... else
        }// if r2_n < h2
      }// for j=jmin;j<jmax,j++
    }// if fbas(dx_n) < h
  }// for i=imin;i<imax;i++
}// void neighborloop

// changed to better fit python wrapper, not IDL //
int triband_alphablend(
    int N_xy, // number of input particles/positions
    int use_kernel,
    float* x, float* y, // positions (assumed already sorted in z)
    float* hsml, // smoothing lengths for each
    float* R, float* G, float* B, float* alphas,
    float Xmin, float Xmax, float Ymin, float Ymax, // boundaries of output grid
    int Xpixels, int Ypixels, // dimensions of grid
    float* OUT0, float* OUT1, float* OUT2) // output vectors for weightMap and weightWeightedQuantityMap
{
  // print out the input parameters // 
  printf("N_xy=%d...",N_xy); 
  printf("Xmin=%f...Xmax=%f...Ymin=%f...Ymax=%f...Xpixels=%d...Ypixels=%d\n",
    Xmin,Xmax,Ymin,Ymax,Xpixels,Ypixels);

  double dx, dy, dx_i, dy_i, dx_n, dy_n, i_x_flt, i_y_flt, d_ij, h, hmin;
  double h2, x2_n, y2_n, r2_n, h2_i, wk, hkernel_over_hsml_to_use, *Kernel; 
  double * wt_sum;

  double dpi=3.1415926535897932384626433832795;
  long n,i,j,k,imin,imax,jmin,jmax,N_KERNEL_TABLE;
  
  dx = (Xmax - Xmin)/((double)Xpixels);
  dy = (Ymax - Ymin)/((double)Ypixels);
  dx_i = 1./dx; dy_i = 1./dy;
  hmin = 0.5*sqrt(dx*dx+dy*dy); // ensures at least one cell 'sees' the particle // 
  
  // pre-define cell positions so we save a step in the loop //
  double x_i[Xpixels]; for(i=0;i<Xpixels;i++) x_i[i]=Xmin+dx*((double)i+0.5);
  double y_j[Ypixels]; for(i=0;i<Ypixels;i++) y_j[i]=Ymin+dy*((double)i+0.5);
  
  // build a kernel lookup table to save on the calculation below // 
  hkernel_over_hsml_to_use = 1.0; 
  //hkernel_over_hsml_to_use = 2.0; 
  // default = 1 (integrate out to kernel), but low-density regions are represented 
  //   more accurately if this is larger (~2); code expense increases though!
  N_KERNEL_TABLE = 1000;
  Kernel = calloc(N_KERNEL_TABLE+1, sizeof(double)); 
  dx_n=(hkernel_over_hsml_to_use*hkernel_over_hsml_to_use)/((double)N_KERNEL_TABLE); r2_n=0.;
  double kernel_spacing_inv = 1./dx_n;
  for(n=0;n<N_KERNEL_TABLE;n++)
  {
   h = sqrt(r2_n); // radius (to save sqrt operation we're interpolating in r^2 //
   // approximate gaussian for projected, integrate kernel: //
        //wk = (135./(14.*dpi)) * exp(-(135./14.)*h*h); // quintic spline we're using now   
        //wk = (135./(14.*dpi)) * exp(-(135./14.)*h*h/(1.+h)); 
        // this has more extended tails, designed to reduce artifacts in low-density regions
        //  (where we would really want to use a proper volume-render)
        // wk = 1.91 * exp(-5.56*h*h); // cubic spline 
   // cubic spline kernel
        h2=(1.-h); if(h<=0.5) {wk=(1.-6.*h*h*h2);} else {wk=2.*h2*h2*h2;} wk*=8/dpi; //wk*=40./(7.0*dpi);
   Kernel[n] = wk;
   r2_n += dx_n; // radius at this point in the table
  }
  Kernel[N_KERNEL_TABLE]=0;
  
  // zero out the output vectors before the main sum
  for(n=0;n<Xpixels*Ypixels;n++)
  {
    OUT0[n]=0.0; OUT1[n]=0.0;
  }
  
  // loop over particles // 
  for(n=0;n<N_xy;n++)
  {
    h = hsml[n]; if(h<hmin) h=hmin; // assume 'intrinsic' h is smeared by some fraction of pixel
    h2=h*h;
    h2_i = 1./(h*h); // here we need the 'real' h (not the expanded search) // 
    h *= hkernel_over_hsml_to_use; // make search area larger for kernel //

    i_x_flt = (x[n] - Xmin) * dx_i;
    i_y_flt = (y[n] - Ymin) * dy_i;
    d_ij=h*dx_i; imin=(long)(i_x_flt-d_ij); imax=(long)(i_x_flt+d_ij)+1; if(imin<0) imin=0; if(imax>Xpixels-1) imax=Xpixels-1;
    d_ij=h*dy_i; jmin=(long)(i_y_flt-d_ij); jmax=(long)(i_y_flt+d_ij)+1; if(jmin<0) jmin=0; if(jmax>Ypixels-1) jmax=Ypixels-1;
    
    // initialize wt_sum assuming we won't use the kernel and want
    //  flat opacities over the surface
    *wt_sum = 1;

    if (use_kernel){
      // need to accumulate the total weight that will be deposited in order to renormalize
      *wt_sum = 0.;
      neighborloop(
        imin, imax,
        jmin, jmax,
        x, y,
        x_i, y_j,
    // pass-through scalars
        n,
        use_kernel,
        h, h2, h2_i,
        dx_n, dy_n,
        Kernel,
        kernel_spacing_inv,
        Ypixels,
    // output
        wt_sum, // = 0 initially -> accumulate = True
        NULL, NULL, NULL,NULL,
        NULL, NULL, NULL);
    }

    // now that wt_sum is set, actually deposit the mass
    neighborloop(
      imin, imax,
      jmin, jmax,
      x, y,
      x_i, y_j,
  // pass-through scalars
      n,
      use_kernel,
      h, h2, h2_i,
      dx_n, dy_n,
      Kernel,
      kernel_spacing_inv,
      Ypixels,
  // output
      wt_sum, // != 0 initially -> accumulate = False
      R, G, B, alphas,
      OUT0, OUT1, OUT2);

  } // for(n=0;n<N_xy;n++)
  return 1;
} // closes main program 


