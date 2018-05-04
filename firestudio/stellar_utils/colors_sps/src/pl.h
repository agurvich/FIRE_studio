#ifndef PL_BRANT
#define PL_BRANT
extern int nm;
extern int nz;
extern int na;
extern double *zz;
extern double *lage;
extern double *m;
void AllocateMagnitudes(void);
void FreeMagnitudes(void);
void GetMags(double lpa, double pz, double mags[]);
void LoadMagnitudeData(void);
#endif PL_BRANT
