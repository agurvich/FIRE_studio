

void ngb3d_treebuild(float **pospointer,int Npart,int MinMaxFlag,float XminOpt[3],float XmaxOpt[3]);

void ngb3d_treefree(void);
void ngb3d_treeallocate(int npart,int maxnodes);  /* usually maxnodes=2*npart is suffiecient */

float ngb3d_treefind(float xyz[3], int desngb, float hguess,int **ngblistback, float **r2listback,float hmax,int *ngbfound);


float ngb3d_treetest(float xyz[3], int desngb, float hguess,int **ngblistback, float **r2listback);


void allocate_ngblists(void);
void free_ngblists(void);


