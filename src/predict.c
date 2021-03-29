#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_math.h>

#include "allvars.h"
#include "proto.h"


/*! This function makes sure that all particle coordinates (Pos) are
 *  periodically mapped onto the interval [0, BoxSize].  After this function
 *  has been called, a new domain decomposition should be done, which will
 *  also force a new tree construction.
 */
#ifdef PERIODIC
void do_box_wrapping(void)
{
  int i, j;
  double boxsize[3];

  for(j = 0; j < 3; j++)
    boxsize[j] = All.BoxSize;

	#ifdef LONG_X
  	boxsize[0] *= LONG_X;
	#endif
	#ifdef LONG_Y
  	boxsize[1] *= LONG_Y;
	#endif
	#ifdef LONG_Z
  	boxsize[2] *= LONG_Z;
	#endif

  for(i = 0; i < NumPart; i++)
    for(j = 0; j < 3; j++)
    {
	    while(P[i].Pos[j] < 0)
	      P[i].Pos[j] += boxsize[j];

	    while(P[i].Pos[j] >= boxsize[j])
	      P[i].Pos[j] -= boxsize[j];
    }
}
#endif
