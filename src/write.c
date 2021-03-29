/************************************************
Title: write.c
Purpose: Writes the dm data (including rho and T)
         to a snapshot for use in exspec.
Notes:   * This function based on rsb.c
************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include "allvars.h"
#include "proto.h"



/***********************
       write_data
***********************/
void write_data(char *snapname)
{
   // Does as the name says

   FILE *fd;
   char dspec_file[100];
   char proc_num[4];
   char *split_str;
   int blksize;
   int pc = 0;
   int pc_new;
   int k;
   int n;
   int n_with_mass = 0; // This is very bad, but I don't feel like figuring it out (is the header on each
                        // processor different?) Also, since the particles are really dm particles, their 
                        // masses will always all be the same

   // ORIGINAL
   // Open file for writing
   sprintf(dspec_file, "%s-dspec.%d", snapname, ThisTask);

   // Split snapfilename. dspec is passed snapshot-fakegas_xxx
   /*split_str = strtok(snapname, "_.");

   // Copy first part of name over, the snapshot-fakegas part
   strcat(dspec_file, split_str);

   // Now append -dspec_ to the name
   strcat(dspec_file, "-dspec_");

   // Now append the snapshot number
   split_str = strtok(NULL, "_.");
   strcat(dspec_file, split_str);

   // Now append the processor number
   strcat(dspec_file, ".");
   sprintf(proc_num, "%d", ThisTask);
   strcat(dspec_file, proc_num);
*/
   if(!(fd = fopen(dspec_file, "wb")))
   {
      printf("Error, could not open file for writing dm data on task %d!\n", ThisTask);
      exit(EXIT_FAILURE);
   }

   // Write the header
   blksize = 256;

   // Change types back to DM (this shouldn't affect the write loops below)
   // Check to make sure there are gas particles first! (i.e. I didn't forget to
   // run the snapshot through fake_gas!!)
   if(N_gas == 0)
   {
      printf("Error, N_gas = 0. Forget to use fake_gas?\n");
      exit(EXIT_FAILURE);
   }

   header.npart[0] = 0;
   header.npart[1] = N_gas;
   header.npartTotal[1] = header.npartTotal[0];
   header.npartTotal[0] = 0;
   header.mass[1] = header.mass[0];
   header.mass[0] = 0.0; 

   // Change number of files per snapshot to be number of processors
   header.num_files = NTask;
 
   my_fwrite(&blksize, sizeof(int), 1, fd);
   my_fwrite(&header, sizeof(struct io_header), 1, fd);
   my_fwrite(&blksize, sizeof(int), 1, fd);

   // Positions
   blksize = N_gas * 3 * sizeof(float);
   my_fwrite(&blksize, sizeof(int), 1, fd);
   for(k = 0, pc_new = pc; k < 6; k++)
   {
      for(n = 0; n < header.npart[k]; n++)
      {
         my_fwrite(&SP[pc_new].Pos[0], sizeof(float), 3, fd);
         pc_new++;
      }
   }
   my_fwrite(&blksize, sizeof(int), 1, fd);
   
   // Velocities
   blksize = N_gas * 3 * sizeof(float);
   my_fwrite(&blksize, sizeof(int), 1, fd);

   double a3inv = 1.0 / (header.time * header.time * header.time);
   for(k = 0, pc_new = pc; k < 6; k++)
   {
      for(n = 0; n < header.npart[k]; n++)
      {
         if(All.ComovingIntegrationOn == 1)
         {
            // Convert velocity from p back to u
            SP[pc_new].Vel[0] *= sqrt(a3inv);
            SP[pc_new].Vel[1] *= sqrt(a3inv);
            SP[pc_new].Vel[2] *= sqrt(a3inv);
         }
         my_fwrite(&SP[pc_new].Vel[0], sizeof(float), 3, fd);
         pc_new++;
      }
   }
   my_fwrite(&blksize, sizeof(int), 1, fd);

   // IDs
   blksize = N_gas * sizeof(int);
   my_fwrite(&blksize, sizeof(int), 1, fd);
   for(k = 0, pc_new = pc; k < 6; k++)
   {
      for(n = 0; n < header.npart[k]; n++)
      {
         my_fwrite(&SP[pc_new].ID, sizeof(int), 1, fd);
         pc_new++;
      }
   }
   my_fwrite(&blksize, sizeof(int), 1, fd);
   
   // Masses
   blksize = n_with_mass * sizeof(float);

   if(n_with_mass > 0)
   {
      my_fwrite(&blksize, sizeof(int), 1, fd);
   }

   for(k = 0, pc_new = pc; k < 6; k++)
   {
      for(n = 0; n < header.npart[k]; n++)
      {
         if(header.mass[k] == 0)
         {
            my_fwrite(&SP[pc_new].Mass, sizeof(float), 1, fd);
         }
         pc_new++;  
      }
   }

   if(n_with_mass > 0)
   {
      my_fwrite(&blksize, sizeof(int), 1, fd);
   }

   // Sph properties
   // Density
   blksize = N_gas * sizeof(float);
   my_fwrite(&blksize, sizeof(int), 1, fd);
   for(k = 0, pc_new = pc; k < 6; k++)
   {
      for(n = 0; n < header.npart[k]; n++)
      {
         my_fwrite(&SP[pc_new].Density, sizeof(float), 1, fd);
         pc_new++;
      }
   }
   my_fwrite(&blksize, sizeof(int), 1, fd);

   // Hsml
   blksize = N_gas * sizeof(float);
   my_fwrite(&blksize, sizeof(int), 1, fd);
   for(k = 0, pc_new = pc; k < 6; k++)
   {
      for(n = 0; n < header.npart[k]; n++)
      {
         my_fwrite(&SP[pc_new].Hsml, sizeof(float), 1, fd);
         pc_new++;
      }
   }
   my_fwrite(&blksize, sizeof(int), 1, fd);

   // Close the file
   fclose(fd);
}
