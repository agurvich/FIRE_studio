#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "proto.h"
#include "ngbtree3d.h"

void allocate_3d(int Nsize);
void set_particle_pointer(int Nsize, float* x, float* y, float* z);
void free_memory_3d(int Nsize);

struct particle_3d 
{
  float Pos[3];
} **P3d;


// revised call for python calling //
//int stellarhsml(int argc,void *argv[])
int stellarhsml(int N_in, float* x, float* y, float* z, int DesNgb, float Hmax, float* H_OUT)
{
  float h_guess, h2, xyz[3], dummy[3], h_guess_0;
  int i, ngbfound;
  allocate_3d(N_in);
  float *r2list;
  int *ngblist;

  //printf("N=%d\n",N_in); printf("Hmax=%g\n",Hmax); printf("DesNgb=%d\n",DesNgb);

  ngb3d_treeallocate(N_in, 2*N_in);
  set_particle_pointer(N_in, x,y,z);
  ngb3d_treebuild((float **)&P3d[1], N_in, 0, dummy, dummy);
  h_guess = Hmax/150.0e0; h_guess_0=h_guess;
  for(i=0;i<N_in;i++)
  {
	  xyz[0]=P3d[i+1]->Pos[0]+1.0e-10;
	  xyz[1]=P3d[i+1]->Pos[1]+1.0e-10;
	  xyz[2]=P3d[i+1]->Pos[2]+1.0e-10;
	  h2=ngb3d_treefind( xyz, DesNgb ,1.04*h_guess, &ngblist, &r2list, Hmax, &ngbfound); 

    /*
    if(!(i%10000))
    {
    printf("i=%d hmax=%g h_guess=%g h=%g xyz=%g|%g|%g ngb=%d \n",
        i,Hmax,h_guess,sqrt(h2),xyz[0],xyz[1],xyz[2],ngbfound); fflush(stdout);
    }
    */
      H_OUT[i] = sqrt(h2);
      h_guess = H_OUT[i]; // use this value for next guess, should speed things up // 
      //if (h_guess>10.*h_guess_0) h_guess=2.*h_guess_0;
    } 

  ngb3d_treefree();
  free_memory_3d(N_in);
  //printf("done\n");
  return 0;
}



void set_particle_pointer(int Nsize, float* x, float* y, float* z)
{
  int i;
  float *pos;
  for(i=1;i<=Nsize;i++)
    {
      P3d[i]->Pos[0] = x[i-1];
      P3d[i]->Pos[1] = y[i-1];
      P3d[i]->Pos[2] = z[i-1];
    }
}


void allocate_3d(int Nsize)
{
  //printf("allocating memory...\n");
  int i;
  if(Nsize>0)
    {
      if(!(P3d=malloc(Nsize*sizeof(struct particle_3d *))))
	{
	  //printf("failed to allocate memory. (A)\n");
	  exit(0);
	}
      P3d--;   /* start with offset 1 */
      if(!(P3d[1]=malloc(Nsize*sizeof(struct particle_3d))))
	{
	  //printf("failed to allocate memory. (B)\n");
	  exit(0);
	}
      for(i=2;i<=Nsize;i++)   /* initiliaze pointer table */
	P3d[i]=P3d[i-1]+1;
    }
  //printf("allocating memory...done\n");
}


void free_memory_3d(int Nsize)
{
  if(Nsize>0)
    {
      free(P3d[1]);
      P3d++;
      free(P3d);
    }
}
