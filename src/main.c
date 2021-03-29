#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"


/*! \file main.c
 *  \brief start of the program
 */

/*!
 *  This function initializes the MPI communication packages, and sets
 *  cpu-time counters to 0.  Then begrun() is called, which sets up
 *  the simulation either from IC's or from restart files.  Finally,
 *  run() is started, the main simulation loop, which iterates over
 *  the timesteps.
 */
int main(int argc, char **argv)
{
	/**************************** HELP SECTION ******************************/
	if( argv[1][0] == '-' && argv[1][1] == 'h')
	{
      #ifndef SCALARFIELD
		   printf("USAGE: mpirun -np 2 ./dspec snapshot_xxx\n");
		   exit(0);
      #else
         printf("USAGE: mpirun -np 2 ./dspec snapshot_xxx\n");
         exit(0);
      #endif
	}

   double t0, t1, t2;
   FILE *com;
   char com_file[100];

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
  MPI_Comm_size(MPI_COMM_WORLD, &NTask);


  // Set args
  sprintf(snap_filename, "%s", argv[1]);

   #ifdef SCALARFIELD
      sprintf(dark_file, "%s", argv[4]);
   #endif

  for(PTask = 0; NTask > (1 << PTask); PTask++);

  //RestartFlag = 2;    //For reading from a snapshot.
   // Jared added because it was failing to read in the snapshot because it expected there
   // to be densities and such when I don't have those yet
   RestartFlag = 0;

  All.CPU_TreeConstruction = All.CPU_TreeWalk = All.CPU_Gravity = All.CPU_Potential = 
		All.CPU_Domain = All.CPU_Snapshot = All.CPU_Total = All.CPU_CommSum = 
		All.CPU_Imbalance = All.CPU_Hydro = All.CPU_HydCompWalk = All.CPU_HydCommSumm = 
		All.CPU_HydImbalance = All.CPU_EnsureNgb = All.CPU_Predict = All.CPU_TimeLine = 
		All.CPU_PM = All.CPU_Peano = 0;

  CPUThisRun = 0;

  t0 = second();

  begrun();			/* set-up run  */

  t1 = second();
  CPUThisRun += timediff(t0, t1);
  All.CPU_Total += timediff(t0, t1);

   if(ThisTask == 0)
   {
      printf("Tree constructed!\n");  
      fflush(stdout);
   }

   // Print out the COM of the root node in the tree
   // This is probably wrong, but as far as I can tell, the
   // root node is Nodes[All.MaxPart] (see forcetree.c comments)
   /*if(ThisTask == 0)
   {
      sprintf(com_file, "%s-COM", argv[1]);
      if(!(com = fopen(com_file, "w")))
      {
         printf("Error, could not open file for writing COM!\n");
         MPI_Finalize();
         exit(EXIT_FAILURE);
      }

      fprintf(com, "%f \t %f \t %f\n", Nodes[All.MaxPart].u.d.s[0], 
             Nodes[All.MaxPart].u.d.s[1], Nodes[All.MaxPart].u.d.s[2]);

      fclose(com);
   }*/

   // Get which particles are in halos and the m_vir of those halos
   if(ThisTask == 0)
   {
      printf("Simplifying particle data...\n");
   }
   simplify_particle_data();

   // Write the data to a new snapshot for use in exspec (MAKE SURE TO CHANGE THE
   // TYPE OF THE DM PARTICLES BACK TO 1!!!!) 
   if(ThisTask == 0)
   {
      printf("Writing data...\n");
   }
   write_data(snap_filename);

  MPI_Finalize();		/* clean up & finalize MPI */

	t2 = second();

	printf("Total Time to Run : %fs\n", t2 - t0);

  return 0;
}
