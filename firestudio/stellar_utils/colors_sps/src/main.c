#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<gsl/gsl_interp.h>
#include"pl.h"

int colors(int argc, void *argv[])
{
	int   N;
	float z_obs;
	float zz;
	float *a_form;
	float *Age;
	float *Zmet;
	float *Mass;
	float *Value;
	float *vv;
	int i,j;

	int nm = 14;
	double *mags;

	N    =*(int *)argv[0];
	z_obs=*(float *)argv[1];
	Mass =(float*)argv[2];
	Age  =(float*)argv[3];
	Zmet =(float*)argv[4];
	Value=(float*)argv[5];

	//printf("....finding colors......");
	mags = (double *) malloc(nm*sizeof(double));

	AllocateMagnitudes();
	LoadMagnitudeData();

	for(i=0;i<N;i++)
	{
		//if(!(i%10000)) printf("%d..",i);
		GetMags(log10(Age[i]*1.0e9),Zmet[i],mags);
		for(j=0;j<nm;j++)
		{
			vv = Value + nm*i +j;
			*vv = mags[j];
		}
	}
	for(i=0;i<N;i++)
		vv = Value + nm*i;

	//printf("\n");
	FreeMagnitudes();
	free(mags);

	return 0;

}
