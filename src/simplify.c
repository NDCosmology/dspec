/************************************************
Title: simplify.c
Purpose: Keeps only the bare essential data to make
         writing it easier.
Notes:
************************************************/ 
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include "allvars.h"
#include "proto.h"



/***********************
 simplify_particle_data
***********************/
void simplify_particle_data(void)
{
   // Most of the elements of P are unneeded for this code,
   // so, to enable smaller communications between processors
   // as well as simplify my life, I'm just going to copy over
   // the essential members of P to a simplified version SP

   int i;
   int j;

   // Allocate memory
   if(!(SP = calloc(N_gas, sizeof(SIMPLE_PARTICLE_DATA))))
   {
      printf("Error, could not allocate memory for SP on task %d\n", ThisTask);
      MPI_Finalize();
      exit(EXIT_FAILURE);
   }

   // Copy over
   for(i = 0; i < N_gas; i++)
   {
      for(j = 0; j < 3; j++)
      {
         SP[i].Pos[j] = P[i].Pos[j];
         SP[i].Vel[j] = P[i].Vel[j];
      }

      SP[i].Mass = P[i].Mass;
      SP[i].Density = SphP[i].Density;
      SP[i].Hsml = SphP[i].Hsml;
      SP[i].Type = P[i].Type;
      SP[i].ID = P[i].ID;
   }

   // Clean up
   free(P);
   free(SphP);
}
