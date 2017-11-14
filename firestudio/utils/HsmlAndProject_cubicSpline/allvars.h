
#ifndef ALLVARS_H
#define ALLVARS_H

#include <stdio.h>

extern int     NumPart;

extern double TreeAllocFactor;

extern int *Head, *Next, *Tail, *Len;

extern int *GroupLen, *GroupOffset;

extern FILE *Logfile;


extern double BoxSize, BoxHalf;

extern int    DesDensNgb;

extern int Axis1, Axis2, Axis3;

extern float *Hsml;
extern float *Quantity;
extern float *Mass;
extern float *Rho;

extern float Xmin, Ymin, Xmax, Ymax, Zmin, Zmax, Hmax;

extern float Xc, Yc;

extern int Xpixels, Ypixels;

extern float *Value, *ValueQuantity;

extern float Softening;


extern struct particle_data
{
  float Pos[3];			/*!< particle position at its current time */
}
*P;                            /*!< points to particles on this processor */





extern struct r2data
{
  float r2;
  int   index;
}
*R2list;






extern int    AnzNodes;
extern int    MaxNodes;



extern int *Nextnode;
extern int *Father;

extern struct NODE
{
  float len;			/*!< sidelength of treenode */
  float center[3];		/*!< geometrical center of node */
  union
  {
    int suns[8];		/*!< temporary pointers to daughter nodes */
    struct
    {
      float s[3];               /*!< center of mass of node */
      int mass;            /*!< mass of node */
      int cost;            /*!< counts the number of interactions in which this node is used */
      int sibling;         /*!< this gives the next node in the walk in case the current node can be used */
      int nextnode;        /*!< this gives the next node in case the current node needs to be opened */
      int father;          /*!< this gives the parent node of each node (or -1 if we have the root node) */
    }
    d;
  }
  u;
}
*Nodes_base,                    /*!< points to the actual memory allocted for the nodes */
*Nodes;                         /*!< this is a pointer used to access the nodes which is shifted such that Nodes[All.MaxPart] 
				  gives the first allocated node */


#endif








