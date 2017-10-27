#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"

#define  BITS_PER_DIMENSION 10	/* maximum is 10 to fit in 32-bit integer ! */

#define  MAX_REAL_NUMBER  1.0e35

int peano_hilbert_key(int x, int y, int z, int bits);
int compare_key(const void *a, const void *b);
void reorder_particles(void);


static struct peano_hilbert_data
{
  int key;
  int index;
}
 *mp;

static int *Id;

void peano_hilbert_order(void)
{
  int i, j;
  double xmin[3], xmax[3], len;
  double scalefac;


  //printf("establishing Peano-Hilbert order...\n");

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

  scalefac = 1.0 / len * ((1 << BITS_PER_DIMENSION) - 1);


  if(NumPart > 0)
    {
      mp = malloc(sizeof(struct peano_hilbert_data) * (NumPart));
      Id = malloc(sizeof(int) * (NumPart));

      for(i = 0; i < NumPart; i++)
	{
	  mp[i].index = i;
	  mp[i].key = peano_hilbert_key((P[i].Pos[0] - xmin[0]) * scalefac,
					(P[i].Pos[1] - xmin[1]) * scalefac,
					(P[i].Pos[2] - xmin[2]) * scalefac, BITS_PER_DIMENSION);
	}

      qsort(mp, NumPart, sizeof(struct peano_hilbert_data), compare_key);

      for(i = 0; i < NumPart; i++)
	Id[mp[i].index] = i;

      reorder_particles();

      free(Id);
      free(mp);
    }

  //printf("done.\n");
}


int compare_key(const void *a, const void *b)
{
  if(((struct peano_hilbert_data *) a)->key < (((struct peano_hilbert_data *) b)->key))
    return -1;

  if(((struct peano_hilbert_data *) a)->key > (((struct peano_hilbert_data *) b)->key))
    return +1;

  return 0;
}



void reorder_particles(void)
{
  int i;
  struct particle_data Psave, Psource;
  int idsource, idsave, dest;
  float mass_source, quantity_source, hsml_source, mass_save, quantity_save, hsml_save;   

  for(i = 0; i < NumPart; i++)
    {
      if(Id[i] != i)
	{
	  Psource = P[i];
	  idsource = Id[i];
          mass_source = Mass[i];
          quantity_source = Quantity[i];
          hsml_source = Hsml[i];


	  dest = Id[i];

	  do
	    {
	      Psave = P[dest];
	      idsave = Id[dest];
              mass_save = Mass[dest];
              quantity_save = Quantity[dest];
              hsml_save = Hsml[dest];


	      P[dest] = Psource;
	      Id[dest] = idsource;
              Mass[dest] = mass_source; 
              Quantity[dest] = quantity_source; 
              Hsml[dest] = hsml_source;

	      if(dest == i)
		break;

	      Psource = Psave;
	      idsource = idsave;
              mass_source = mass_save;
              quantity_source = quantity_save;
              hsml_source = hsml_save;

	      dest = idsource;
	    }
	  while(1);
	}
    }
}






static int quadrants[24][2][2][2] = {
  /* rotx=0, roty=0-3 */
  {{{0, 7}, {1, 6}}, {{3, 4}, {2, 5}}},
  {{{7, 4}, {6, 5}}, {{0, 3}, {1, 2}}},
  {{{4, 3}, {5, 2}}, {{7, 0}, {6, 1}}},
  {{{3, 0}, {2, 1}}, {{4, 7}, {5, 6}}},
  /* rotx=1, roty=0-3 */
  {{{1, 0}, {6, 7}}, {{2, 3}, {5, 4}}},
  {{{0, 3}, {7, 4}}, {{1, 2}, {6, 5}}},
  {{{3, 2}, {4, 5}}, {{0, 1}, {7, 6}}},
  {{{2, 1}, {5, 6}}, {{3, 0}, {4, 7}}},
  /* rotx=2, roty=0-3 */
  {{{6, 1}, {7, 0}}, {{5, 2}, {4, 3}}},
  {{{1, 2}, {0, 3}}, {{6, 5}, {7, 4}}},
  {{{2, 5}, {3, 4}}, {{1, 6}, {0, 7}}},
  {{{5, 6}, {4, 7}}, {{2, 1}, {3, 0}}},
  /* rotx=3, roty=0-3 */
  {{{7, 6}, {0, 1}}, {{4, 5}, {3, 2}}},
  {{{6, 5}, {1, 2}}, {{7, 4}, {0, 3}}},
  {{{5, 4}, {2, 3}}, {{6, 7}, {1, 0}}},
  {{{4, 7}, {3, 0}}, {{5, 6}, {2, 1}}},
  /* rotx=4, roty=0-3 */
  {{{6, 7}, {5, 4}}, {{1, 0}, {2, 3}}},
  {{{7, 0}, {4, 3}}, {{6, 1}, {5, 2}}},
  {{{0, 1}, {3, 2}}, {{7, 6}, {4, 5}}},
  {{{1, 6}, {2, 5}}, {{0, 7}, {3, 4}}},
  /* rotx=5, roty=0-3 */
  {{{2, 3}, {1, 0}}, {{5, 4}, {6, 7}}},
  {{{3, 4}, {0, 7}}, {{2, 5}, {1, 6}}},
  {{{4, 5}, {7, 6}}, {{3, 2}, {0, 1}}},
  {{{5, 2}, {6, 1}}, {{4, 3}, {7, 0}}}
};


static int rotxmap_table[24] = { 4, 5, 6, 7, 8, 9, 10, 11,
  12, 13, 14, 15, 0, 1, 2, 3, 17, 18, 19, 16, 23, 20, 21, 22
};

static int rotymap_table[24] = { 1, 2, 3, 0, 16, 17, 18, 19,
  11, 8, 9, 10, 22, 23, 20, 21, 14, 15, 12, 13, 4, 5, 6, 7
};

static int rotx_table[8] = { 3, 0, 0, 2, 2, 0, 0, 1 };
static int roty_table[8] = { 0, 1, 1, 2, 2, 3, 3, 0 };

static int sense_table[8] = { -1, -1, -1, +1, +1, -1, -1, -1 };


int peano_hilbert_key(int x, int y, int z, int bits)
{
  int i, key, quad, bitx, bity, bitz;
  int mask, rotation, rotx, roty, sense;

  mask = 1 << (bits - 1);
  key = 0;
  rotation = 0;
  sense = 1;


  for(i = 0; i < bits; i++, mask >>= 1)
    {
      bitx = (x & mask) ? 1 : 0;
      bity = (y & mask) ? 1 : 0;
      bitz = (z & mask) ? 1 : 0;

      quad = quadrants[rotation][bitx][bity][bitz];

      key <<= 3;
      key += (sense == 1) ? (quad) : (7 - quad);

      rotx = rotx_table[quad];
      roty = roty_table[quad];
      sense *= sense_table[quad];

      while(rotx > 0)
	{
	  rotation = rotxmap_table[rotation];
	  rotx--;
	}

      while(roty > 0)
	{
	  rotation = rotymap_table[rotation];
	  roty--;
	}
    }

  return key;
}
