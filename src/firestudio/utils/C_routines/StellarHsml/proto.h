
#include "nrutil.h"

#include "ngbtree.h"
#include "ngbtree3d.h"

#include "forcetree.h"


int loadpositions(char *fname,int allocflag);


double selectb(unsigned long k, unsigned long n, double arr[],int ind[]);


int set_particles_into_zones(void);

double cooling_rate(double dens,double press);

void treebuild(void);
int  find_neighbors(double xyz[3],double r,int *nl,double *r2l,int maxngb);



void begrun(void);
void endrun(int ierr);

void read_initial_conditions(char *fname);

/* in eos.c */
double cool(double enp1);
void   enttab(int irho, int itemp);

void   eospre(void);
void   eoscor(void);

void   eos2(int iact,int *iter2,int *iterflag);
void   eos3(int iact,int *iter2,int *iterflag);

void   settab(void);

void   getion(double Density,double EnergySpecific,
	      double *temp,double *decol,double *photo,double *MassFracIons);

void   xflux(void);

double phh(double xnu);
double phhe(double xnu);
double phhe2(double xnu);
double pchh(double xnu);
double pchhe(double xnu);
double pchhe2(double xnu);
double xinu(double xnu);
double xinunu(double xnu);
/****/


void increment_file_suffix(void);

void frduni(void);
double z2t(double z,double OmegaMatter,double OmegaLambda,double HubbleToday);
double t2a(double Time,double OmegaMatter,double OmegaLambda,double HubbleToday,double TimeToday);
double ht0(double eta);


void glbout(void);

void init(void);
void read_parameter_file(char *fname);

void forwrd(void);

void bckwrd(void);

void oldnew(void);

void pltout(void);

void prnout(void);

void propag(void);



void restart(int mod);



/*** system.c ***/
double second(void);
double tremain(void);
double dmin(double,double);
double dmax(double,double);
int imin(int,int);
int imax(int,int);
int    nint(double);
/****/

void densty(int iact, int nsph, int *isph);
void trweos(void);


void grvfrc(void);
void trwfrc(void);



void trwngh(void);


void tstep(int startflag);







double zbrentMS(double (*func)(double), double x1, double x2, double tol,int *ierr);
double midinf(double (*funk)(double), double aa, double bb, int n);
void polint(double xa[], double ya[], int n, double x, double *y, double *dy);
double zriddr(double (*func)(double), double x1, double x2, double xacc);




void indexintx(unsigned long n, int arr[], unsigned long indx[]);
void indexshortx(unsigned long n, short arr[], int indx[]);
void indexx(unsigned long n, double arr[], int indx[]);
void rank(unsigned long n, int indx[], int irank[]);



double selip(unsigned long k, unsigned long n, double arr[]);

