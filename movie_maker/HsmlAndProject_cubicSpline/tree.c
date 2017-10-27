#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define  MAX_REAL_NUMBER  1.0e35

#include "allvars.h"
#include "proto.h"


static int last;		/* auxialiary variable used to set-up non-recursive walk */

#ifdef PERIODIC
#define NEAREST(x) (((x)>boxhalf)?((x)-boxsize):(((x)<-boxhalf)?((x)+boxsize):(x)))
static double boxhalf, boxsize;
#else

#define NEAREST(x) (x)

#endif





int tree_treebuild(void)
{
  int i, j, subnode = 0, parent = -1, numnodes;
  int nfree, th, nn;
  double lenhalf;
  struct NODE *nfreep;
#ifndef PERIODIC
  float xmin[3], xmax[3], len;
#endif


  Nodes = Nodes_base - NumPart;

  /* select first node */
  nfree = NumPart;
  nfreep = &Nodes[nfree];

  /* create an empty  root node  */

#ifdef PERIODIC
  boxsize = BoxSize;
  boxhalf = 0.5 * BoxSize;
  for(j = 0; j < 3; j++)
    nfreep->center[j] = 0.5 * BoxSize;
  nfreep->len = BoxSize;
#else

  for(j = 0; j < 3; j++)
    {
      xmin[j] = MAX_REAL_NUMBER;
      xmax[j] = -MAX_REAL_NUMBER;
    }

  for(i = 0; i < NumPart; i++)
    {
      for(j = 0; j < 3; j++)
	{
	  if(P[i].Pos[j] > xmax[j])
	    xmax[j] = P[i].Pos[j];
	  if(P[i].Pos[j] < xmin[j])
	    xmin[j] = P[i].Pos[j];
	}
    }

  /* determine maxmimum extension */

  len = xmax[0] - xmin[0];
  if((xmax[1] - xmin[1]) > len)
    len = xmax[1] - xmin[1];
  if((xmax[2] - xmin[2]) > len)
    len = xmax[2] - xmin[2];

  for(j = 0; j < 3; j++)
    nfreep->center[j] = 0.5* (xmin[j]+xmax[j]);
  nfreep->len = len;
#endif

  for(i = 0; i < 8; i++)
    nfreep->u.suns[i] = -1;

  numnodes = 1;
  nfreep++;
  nfree++;

  /* insert all particles */

  for(i = 0; i < NumPart; i++)
    {
      th = NumPart;
      
      while(1)
        {
          if(th >= NumPart)	/* we are dealing with an internal node */
            {
              subnode = 0;
              if(P[i].Pos[0] > Nodes[th].center[0])
                subnode += 1;
              if(P[i].Pos[1] > Nodes[th].center[1])
                subnode += 2;
              if(P[i].Pos[2] > Nodes[th].center[2])
                subnode += 4;
              
              nn = Nodes[th].u.suns[subnode];
              
              if(nn >= 0)	/* ok, something is in the daughter slot already, need to continue */
                {
                  parent = th;	/* note: subnode can still be used in the next step of the walk */
                  th = nn;
                }
              else
                {
                  /* here we have found an empty slot where we can 
                   * attach the new particle as a leaf 
                   */
                  Nodes[th].u.suns[subnode] = i;
                  break;	/* done for this particle */
                }
            }
          else
            {
              /* we try to insert into a leaf with a single particle
               * need to generate a new internal node at this point 
               */
              Nodes[parent].u.suns[subnode] = nfree;
	      
	      nfreep->len = 0.5 * Nodes[parent].len;
	      lenhalf = 0.25 * Nodes[parent].len;

		  if(subnode & 1)
		    nfreep->center[0] = Nodes[parent].center[0] + lenhalf;
		  else
		    nfreep->center[0] = Nodes[parent].center[0] - lenhalf;

		  if(subnode & 2)
		    nfreep->center[1] = Nodes[parent].center[1] + lenhalf;
		  else
		    nfreep->center[1] = Nodes[parent].center[1] - lenhalf;

		  if(subnode & 4)
		    nfreep->center[2] = Nodes[parent].center[2] + lenhalf;
		  else
		    nfreep->center[2] = Nodes[parent].center[2] - lenhalf;

		  nfreep->u.suns[0] = -1;
		  nfreep->u.suns[1] = -1;
		  nfreep->u.suns[2] = -1;
		  nfreep->u.suns[3] = -1;
		  nfreep->u.suns[4] = -1;
		  nfreep->u.suns[5] = -1;
		  nfreep->u.suns[6] = -1;
		  nfreep->u.suns[7] = -1;

		  subnode = 0;
		  if(P[th].Pos[0] > nfreep->center[0])
		    subnode += 1;
		  if(P[th].Pos[1] > nfreep->center[1])
		    subnode += 2;
		  if(P[th].Pos[2] > nfreep->center[2])
		    subnode += 4;

		  if(nfreep->len < 1.0e-3 * Softening)
		    {
		      /* seems like we're dealing with particles   
		       * at identical locations. randomize 
		       * subnode index (well below gravitational softening scale). 
		       */
		      subnode = (int) (8.0 * drand48());
		      if(subnode >= 8)
			subnode = 7;
		    }

		  nfreep->u.suns[subnode] = th;

		  th = nfree;	/* resume trying to insert the new particle at 
				   the newly created internal node */

		  numnodes++;
		  nfree++;
		  nfreep++;

		  if((numnodes) >= MaxNodes)
		    {
		      printf("maximum number %d of tree-nodes reached.\n", MaxNodes);
		      printf("for particle %d  %g %g %g\n", i, P[i].Pos[0], P[i].Pos[1], P[i].Pos[2]);
		      exit(1);
		    }
		}
	    }
    }

  /* now compute the multipole moments recursively */
  last = -1;
  tree_update_node_recursive(NumPart, -1, -1);

  if(last >= NumPart)
    Nodes[last].u.d.nextnode = -1;
  else
    Nextnode[last] = -1;

  /*
  printf("Have put %d particles into tree.\n", num);
  */

  return numnodes;
}

/* this routine computes the multipole moments for a given internal node and
 * all its subnodes using a recursive computation.  Note that the moments of
 * the daughter nodes are already stored in single precision. For very large
 * particle numbers, loss of precision may results for certain particle
 * distributions
 */
void tree_update_node_recursive(int no, int sib, int father)
{
  int j, jj, p, pp = 0, nextsib, mass, suns[8];
  double s[3];

  if(no >= NumPart)
    {
      for(j = 0; j < 8; j++)
	suns[j] = Nodes[no].u.suns[j];	/* this "backup" is necessary because the nextnode entry will
					   overwrite one element (union!) */
      if(last >= 0)
	{
	  if(last >= NumPart)
	    Nodes[last].u.d.nextnode = no;
	  else
	    Nextnode[last] = no;
	}

      last = no;

      mass = 0;
      s[0] = 0;
      s[1] = 0;
      s[2] = 0;

      for(j = 0; j < 8; j++)
	{
	  if((p = suns[j]) >= 0)
	    {
	      /* check if we have a sibling on the same level */
	      for(jj = j + 1; jj < 8; jj++)
		if((pp = suns[jj]) >= 0)
		  break;

	      if(jj < 8)	/* yes, we do */
		nextsib = pp;
	      else
		nextsib = sib;

	      tree_update_node_recursive(p, nextsib, no);

	      if(p >= NumPart)	/* an internal node or pseudo particle */
		{
		  mass += Nodes[p].u.d.mass;	/* we assume a fixed particle mass */
		  s[0] += Nodes[p].u.d.mass * Nodes[p].u.d.s[0];
		  s[1] += Nodes[p].u.d.mass * Nodes[p].u.d.s[1];
		  s[2] += Nodes[p].u.d.mass * Nodes[p].u.d.s[2];
		}
	      else		/* a particle */
		{
		  mass += 1;
		  s[0] += P[p].Pos[0];
		  s[1] += P[p].Pos[1];
		  s[2] += P[p].Pos[2];
		}
	    }
	}

      if(mass)
	{
	  s[0] /= mass;
	  s[1] /= mass;
	  s[2] /= mass;
	}
      else
	{
	  s[0] = Nodes[no].center[0];
	  s[1] = Nodes[no].center[1];
	  s[2] = Nodes[no].center[2];
	}

      Nodes[no].u.d.s[0] = s[0];
      Nodes[no].u.d.s[1] = s[1];
      Nodes[no].u.d.s[2] = s[2];
      Nodes[no].u.d.mass = mass;

      Nodes[no].u.d.sibling = sib;
      Nodes[no].u.d.father = father;
    }
  else				/* single particle or pseudo particle */
    {
      if(last >= 0)
	{
	  if(last >= NumPart)
	    Nodes[last].u.d.nextnode = no;
	  else
	    Nextnode[last] = no;
	}

      last = no;

      if(no < NumPart)		/* only set it for single particles */
	Father[no] = father;
    }
}







int ngb_compare_key(const void *a, const void *b)
{
  if(((struct r2data *) a)->r2 < (((struct r2data *) b)->r2))
    return -1;

  if(((struct r2data *) a)->r2 > (((struct r2data *) b)->r2))
    return +1;

  return 0;
}


float ngb_treefind(float xyz[3], int desngb, float hguess)
{
  int numngb, iter;
  float h2max;

  if(hguess == 0)
    {
      hguess = Hmax;
    }

  iter = 0;

  do
    {

      iter++;
      /*
         printf("hguess= %g\n", hguess);
       */
      numngb = ngb_treefind_variable(xyz, hguess);

      if(numngb < desngb)
	{
	  hguess *= 1.26;
	  continue;
	}

      if(numngb >= desngb)
	{
	  qsort(R2list, numngb, sizeof(struct r2data), ngb_compare_key);
	  h2max = R2list[desngb - 1].r2;
	  break;
	}

      hguess *= 1.26;

    }
  while(1);

  return sqrt(h2max);
}


int ngb_treefind_variable(float searchcenter[3], float hguess)
{
  int numngb, no, p;
  double dx, dy, dz, r2, h2;
  struct NODE *this;

  h2 = hguess * hguess;

  numngb = 0;
  no = NumPart;

  while(no >= 0)
    {
      if(no < NumPart)		/* single particle */
	{
	  p = no;
	  no = Nextnode[no];

	  dx = NEAREST(P[p].Pos[0] - searchcenter[0]);
	  if(dx < -hguess)
	    continue;
	  if(dx > hguess)
	    continue;

	  dy = NEAREST(P[p].Pos[1] - searchcenter[1]);
	  if(dy < -hguess)
	    continue;
	  if(dy > hguess)
	    continue;

	  dz = NEAREST(P[p].Pos[2] - searchcenter[2]);
	  if(dz < -hguess)
	    continue;
	  if(dz > hguess)
	    continue;

	  r2 = dx * dx + dy * dy + dz * dz;

	  if(r2 < h2)
	    {
	      R2list[numngb].r2 = r2;
	      R2list[numngb].index = p;
	      numngb++;
	    }
	}
      else
	{
	  this = &Nodes[no];

	  no = Nodes[no].u.d.sibling;	/* in case the node can be discarded */

	  if((NEAREST(this->center[0] - searchcenter[0]) + 0.5 * this->len) < -hguess)
	    continue;
	  if((NEAREST(this->center[0] - searchcenter[0]) - 0.5 * this->len) > hguess)
	    continue;
	  if((NEAREST(this->center[1] - searchcenter[1]) + 0.5 * this->len) < -hguess)
	    continue;
	  if((NEAREST(this->center[1] - searchcenter[1]) - 0.5 * this->len) > hguess)
	    continue;
	  if((NEAREST(this->center[2] - searchcenter[2]) + 0.5 * this->len) < -hguess)
	    continue;
	  if((NEAREST(this->center[2] - searchcenter[2]) - 0.5 * this->len) > hguess)
	    continue;

	  no = this->u.d.nextnode;	/* ok, we need to open the node */
	}
    }

  /*
     printf("numngb=%d\n", numngb);
   */
  return numngb;
}








/* this function allocates memory used for storage of the tree
 * and auxiliary arrays for tree-walk and link-lists.
 */
size_t tree_treeallocate(int maxnodes, int maxpart)	/* usually maxnodes=0.7*maxpart is sufficient */
{
  size_t bytes, allbytes = 0;

  MaxNodes = maxnodes;
  //printf("MaxNodes = %d\n", MaxNodes);

  if(!(Nodes_base = malloc(bytes = (MaxNodes + 1) * sizeof(struct NODE))))
    {
      printf("failed to allocate memory for %d tree-nodes (%g MB).\n", MaxNodes, bytes / (1024.0 * 1024.0));
      exit(3);
    }
  allbytes += bytes;

  if(!(Nextnode = malloc(bytes = maxpart * sizeof(int))))
    {
      printf("Failed to allocate %d spaces for 'Nextnode' array (%g MB)\n", maxpart,
	     bytes / (1024.0 * 1024.0));
      exit(4);
    }
  allbytes += bytes;

  if(!(Father = malloc(bytes = maxpart * sizeof(int))))
    {
      printf("Failed to allocate %d spaces for 'Father' array (%g MB)\n", maxpart, bytes / (1024.0 * 1024.0));
      exit(5);
    }
  allbytes += bytes;

  if(!(R2list = malloc(bytes = maxpart * sizeof(struct r2data))))
    {
      printf("failed to allocate memory for R2list\n");
      exit(3);
    }

  allbytes += bytes;

  return allbytes;
}


/* free the allocated memory
 */
void tree_treefree(void)
{
  free(R2list);
  free(Father);
  free(Nextnode);
  free(Nodes_base);
}
