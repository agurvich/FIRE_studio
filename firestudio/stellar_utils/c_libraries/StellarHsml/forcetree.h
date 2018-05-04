


void force_treeallocate(int maxnodes);  /* usually maxnodes=2*npart is suffiecient */
void force_treefree(void);


void force_treebuild(double **pospointer,int Npart,
		     double thetamax,
		     double maxSoftening,
		     double *softeningtable);



void force_treeevaluate(double *targetpart);


void force_testforce(double *targetpart);
