

void ngb_treebuild(double **pospointer,int Npart,int MinMaxFlag,double XminOpt[3],double XmaxOpt[3]);

void ngb_treefree(void);
void ngb_treeallocate(int npart,int maxnodes);  /* usually maxnodes=2*npart is suffiecient */

double ngb_treefind(double xyz[3], int desngb, double hguess,int **ngblistback, double **r2listback);


double ngb_treetest(double xyz[3], int desngb, double hguess,int **ngblistback, double **r2listback);
