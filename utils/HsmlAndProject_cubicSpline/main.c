#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "proto.h"
#include "allvars.h"

int findHsmlAndProject(int _NumPart, struct particle_data* _P, float* _Hsml, float* _Mass, float* _Quantity, float _Xmin, float _Xmax, float _Ymin, float _Ymax, float _Zmin, float _Zmax, int _Xpixels, int _Ypixels, int _DesDensNgb, int _Axis1, int _Axis2, int _Axis3, float _Hmax, double _BoxSize, float* _Value, float* _ValueQuantity);

int findHsmlAndProject(int _NumPart, struct particle_data* _P,
		       float* _Hsml, float* _Mass, float* _Quantity,
		       float _Xmin, float _Xmax, float _Ymin, float _Ymax, float _Zmin, float _Zmax,
		       int _Xpixels, int _Ypixels, int _DesDensNgb, int _Axis1, int _Axis2, int _Axis3,
		       float _Hmax, double _BoxSize, float* _Value, float* _ValueQuantity)

{
 
  NumPart = _NumPart;
  P = _P;
  Hsml = _Hsml;
  Mass = _Mass;
  Quantity = _Quantity;
  Xmin = _Xmin;
  Xmax = _Xmax;
  Ymin = _Ymin;
  Ymax = _Ymax;
  Zmin = _Zmin;
  Zmax = _Zmax;
  Xpixels = _Xpixels;
  Ypixels = _Ypixels;
  DesDensNgb =  _DesDensNgb;
  Axis1 = _Axis1;
  Axis2 = _Axis2;
  Axis3 = _Axis3;
  Hmax = _Hmax;
  BoxSize = _BoxSize;
  Value =         _Value;
  ValueQuantity = _ValueQuantity;

  /*fprintf(stdout, "NumPart = %d\n", NumPart);
  fprintf(stdout, "P = %p \n", P);
  fprintf(stdout, "Hsml = %p \n", Hsml);
  fprintf(stdout, "Mass = %p \n", Mass);
  fprintf(stdout, "Quantity = %p \n", Quantity);
  fprintf(stdout, "Xmin = %f, Xmax = %f \n", Xmin,Xmax);
  fprintf(stdout, "Ymin = %f, Ymax = %f \n", Ymin,Ymax);
  fprintf(stdout, "Zmin = %f, Zmax = %f \n", Zmin,Zmax);
  fprintf(stdout, "Xpixels = %d, Ypixels = %d \n", Xpixels,Ypixels);
  fprintf(stdout, "DesDensNgb = %d \n", DesDensNgb);
  fprintf(stdout, "Axis1 = %d, Axis2 = %d, Axis3 = %d\n", Axis1,Axis2,Axis3);
  fprintf(stdout, "Hmax = %f \n", Hmax);
  fprintf(stdout, "BoxSize = %g \n", BoxSize);
  fprintf(stdout, "Value = %p\n", Value);
  fprintf(stdout, "ValueQuantity = %p\n", ValueQuantity);

  printf("peano-hilbert order...\n");*/
  peano_hilbert_order();
  /*printf("done\n");*/

  if (NumPart < 100000) {
    //printf(" allocating memory for %d tree nodes\n",100*NumPart);
      tree_treeallocate(100*NumPart, NumPart);
  } else {
    //printf(" allocating memory for %d tree nodes\n",100*NumPart);
      tree_treeallocate(100*NumPart, NumPart);
    }
  //Softening = (Xmax-Xmin)/Xpixels/100.;
  Softening = (Xmax-Xmin)/pow(NumPart,0.3333)/10.0;

  //printf("build tree...\n");
  tree_treebuild();
  //printf("done.\n");

  //printf("finding neighbours...\n");
  determine_hsml();
  //printf("done.\n");

  tree_treefree();
    
  //printf("projecting\n");
  make_map();
  //printf("done\n");
    
  // find min/max of arrays
  /*float minval = 1.0e20;
  float maxval = -minval;
  float minvalQ = 1.0e20;
  float maxvalQ = -minvalQ;
  int nelements = Xpixels*Ypixels;
  int i;
  for (i = 0; i < nelements; i++)
    {
      if (Value[i] > maxval)
	{
	  maxval = Value[i];
	}
      else if (Value[i] < minval)
	{
	  minval = Value[i];
	}
      
      if (ValueQuantity[i] > maxvalQ)
	{
	  maxvalQ = ValueQuantity[i];
	}
      else if (ValueQuantity[i] < minvalQ)
	{
	  minvalQ = ValueQuantity[i];
	}    
      }    
  
    fprintf(stdout, "min(Value) = %f, max(Value) = %f\n", minval,maxval);
    fprintf(stdout, "min(ValueQuantity) = %f, max(ValueQuantity) = %f\n", minvalQ,maxvalQ);*/


  return 0;
}

void determine_hsml(void)
{
  int i, signal;
  double h;

  for(i = 0, signal = 0, h = 0; i < NumPart; i++)
    {
      if(i > (signal / 100.0) * NumPart)
        {
          /*printf("x");
	    fflush(stdout);*/
          signal++;
        }

      if(Hsml[i] == 0)
        Hsml[i] = h  = ngb_treefind(P[i].Pos, DesDensNgb, h * 1.1);
    }

  //printf("\n");
}


#ifdef PERIODIC
#define NGB_PERIODIC(x) (((x)>BoxHalf)?((x)-BoxSize):(((x)<-BoxHalf)?((x)+BoxSize):(x)))
#else
#define NGB_PERIODIC(x) (x)
#endif



void make_map(void)
{
  int i, j, n;
  int dx, dy, nx, ny;
  double h, r, u, wk;
  double pos[2];
  double LengthX;
  double LengthY;
  double r2, h2;
  double sum, hmin, hmax, x, y, xx, yy, xxx, yyy;
  double pixelsizeX, pixelsizeY;

  BoxHalf = 0.5 * BoxSize;

  for(i = 0; i < Xpixels; i++)
    for(j = 0; j < Ypixels; j++)
      {
        Value[i * Ypixels + j] = 0;
        ValueQuantity[i * Ypixels + j] = 0;
      }
 

  LengthX = Xmax-Xmin;
  LengthY = Ymax-Ymin;

  Xc = 0.5 * (Xmax + Xmin);
  Yc = 0.5 * (Ymax + Ymin);

  pixelsizeX = LengthX / Xpixels;
  pixelsizeY = LengthY / Ypixels;

  if(pixelsizeX < pixelsizeY)
    hmin = 1.001 * pixelsizeX / 2;
  else
    hmin = 1.001 * pixelsizeY / 2;


  hmax = Hmax;


  for(n = 0; n < NumPart; n++)
    {
      /*if((n % (NumPart / 100)) == 0)
	{
	  printf(".");
	  fflush(stdout);
	  }*/

      if(P[n].Pos[Axis3]< Zmin || P[n].Pos[Axis3] > Zmax)
        continue;

      pos[0]= P[n].Pos[Axis1]-Xmin;
      pos[1]= P[n].Pos[Axis2]-Ymin;


      h = Hsml[n];
      
      if(h < hmin)
        h = hmin;

      if(h > hmax)
        h = hmax;

#ifdef PERIODIC
      if((NGB_PERIODIC(Xc - P[n].Pos[Axis1]) + 0.5 * (Xmax - Xmin)) < -Hsml[n])
        continue;
      if((NGB_PERIODIC(Xc - P[n].Pos[Axis1]) - 0.5 * (Xmax - Xmin)) > Hsml[n])
        continue;
      
      if((NGB_PERIODIC(Yc - P[n].Pos[Axis2]) + 0.5 * (Ymax - Ymin)) < -Hsml[n])
        continue;
      if((NGB_PERIODIC(Yc - P[n].Pos[Axis2]) - 0.5 * (Ymax - Ymin)) > Hsml[n])
        continue;
#else
      if(pos[0] + h < 0 || pos[0] - h >  LengthX
         || pos[1] + h < 0 || pos[1] - h > LengthY)
        continue;
#endif
 
      h2 = h * h;

      nx = h / pixelsizeX + 1;
      ny = h / pixelsizeY + 1;

      /* x,y central pixel of region covered by the particle on the mesh */
      
      x = (floor(pos[0] / pixelsizeX) + 0.5) * pixelsizeX;
      y = (floor(pos[1] / pixelsizeY) + 0.5) * pixelsizeY;

      /* determine kernel normalizaton */      
      sum = 0;

      for(dx = -nx; dx <= nx; dx++)
        for(dy = -ny; dy <= ny; dy++)
          {
            xx = x + dx * pixelsizeX - pos[0];
            yy = y + dy * pixelsizeY - pos[1];
            r2 = xx * xx + yy * yy;
            
            if(r2 < h2)
              {
                r = sqrt(r2);
                u = r / h;

		// The following is for a Cubic Spline kernel. 
                if(u < 0.5)
                  wk = (2.546479089470 + 15.278874536822 * (u - 1) * u * u);
                else
                  wk = 5.092958178941 * (1.0 - u) * (1.0 - u) * (1.0 - u);
				
		// The following is for a Wendland C2 kernel. 
		//wk = 3.342253804929802 * pow(1.0 - u, 4.0) * (1.0 + (4.0 * u));
                
                sum += wk;
              }
          }
      
      if(sum < 1.0e-10)
        continue;

      for(dx = -nx; dx <= nx; dx++)
        for(dy = -ny; dy <= ny; dy++)
          {
            xxx = x + dx * pixelsizeX;
            yyy = y + dy * pixelsizeY;

#ifdef PERIODIC
            xxx = NGB_PERIODIC(xxx + Xmin - Xc) + Xc - Xmin;
            yyy = NGB_PERIODIC(yyy + Ymin - Yc) + Yc - Ymin;
#endif
           
            if(xxx >= 0 && yyy >= 0)
              {
                i = xxx / pixelsizeX;
                j = yyy / pixelsizeY;
                
                if(i >= 0 && i < Xpixels)
                  if(j >= 0 && j < Ypixels)
                    {
                      xx = x + dx * pixelsizeX - pos[0];
                      yy = y + dy * pixelsizeY - pos[1];
                      r2 = xx * xx + yy * yy;
                      
                      if(r2 < h2)
                        {
                          r = sqrt(r2);
                          u = r / h;

						  // Cubic Spline kernel. 
                          if(u < 0.5)
                            wk = (2.546479089470 + 15.278874536822 * (u - 1) * u * u);
                          else
                            wk = 5.092958178941 * (1.0 - u) * (1.0 - u) * (1.0 - u);
						  
						  // Wendland C2 kernel. 
						  //wk = 3.342253804929802 * pow(1.0 - u, 4.0) * (1.0 + (4.0 * u)); 
                          
                          Value[i * Ypixels + j] += Mass[n] * wk / sum;
                          ValueQuantity[i * Ypixels + j] += Mass[n]*Quantity[n]*wk / sum;
                        }
                    }
              }
          }
    }


  for(i = 0; i < Xpixels; i++)
    for(j = 0; j < Ypixels; j++)
      if(Value[i * Ypixels + j]>0)
        ValueQuantity[i * Ypixels + j] /= Value[i * Ypixels + j];


  //printf("\n");
}

