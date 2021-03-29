#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"


/*! \file init.c
 *  \brief Code for initialisation of a simulation from initial conditions
 */


/*! This function reads the initial conditions, and allocates storage for the
 *  tree. Various variables of the particle data are initialised and An intial
 *  domain decomposition is performed. If SPH particles are present, the inial
 *  SPH smoothing lengths are determined.
 */
void init(void)
{
  int i, j;
  double a3;
	#ifdef COOLING
	#ifdef SFR
	#ifdef METALS
		int index = 0.0;
	#endif
	#endif
	#endif


  All.Time = All.TimeBegin;

  switch (All.ICFormat)
  {
    case 1:
			#if (MAKEGLASS > 1)
      	seed_glass();
			#else
      	read_ic(All.InitCondFile);
			#endif
      break;
    case 2:
    case 3:
      read_ic(All.InitCondFile);
      //read_in_snapshot();
      break;
    default:
      if(ThisTask == 0)
      	printf("ICFormat=%d not supported.\n", All.ICFormat);
      endrun(0);
  }
  error_check_header(); 

  All.Time = All.TimeBegin;
  All.Ti_Current = 0;

  if(All.ComovingIntegrationOn)
  {
    All.Timebase_interval = (log(All.TimeMax) - log(All.TimeBegin)) / TIMEBASE;
    a3 = All.Time * All.Time * All.Time;
  }
  else
  {
    All.Timebase_interval = (All.TimeMax - All.TimeBegin) / TIMEBASE;
    a3 = 1;
  }

  //set_softenings();

  All.NumCurrentTiStep = 0;	/* setup some counters */
  All.SnapshotFileCount = 0;
  if(RestartFlag == 2)
    All.SnapshotFileCount = atoi(All.InitCondFile + strlen(All.InitCondFile) - 3) + 1;

  All.TotNumOfForces = 0;
  All.NumForcesSinceLastDomainDecomp = 0;

  if(All.ComovingIntegrationOn)
    if(All.PeriodicBoundariesOn == 1)
      check_omega();

  All.TimeLastStatistics = All.TimeBegin - All.TimeBetStatistics;

  if(All.ComovingIntegrationOn)	/*  change to new velocity variable */
  {
    for(i = 0; i < NumPart; i++)
	    for(j = 0; j < 3; j++)
	      P[i].Vel[j] *= sqrt(All.Time) * All.Time;
  }


	#ifdef PMGRID
  	All.PM_Ti_endstep = All.PM_Ti_begstep = 0;
	#endif

	#ifdef FLEXSTEPS
  	All.PresentMinStep = TIMEBASE;
  	for(i = 0; i < NumPart; i++)	/*  start-up initialization */
  	{
  	  P[i].FlexStepGrp = (int) (TIMEBASE * get_random_number(P[i].ID));
  	}
	#endif

  ngb_treeallocate(MAX_NGB);

  force_treeallocate(All.TreeAllocFactor * All.MaxPart, All.MaxPart);

  All.NumForcesSinceLastDomainDecomp = 1 + All.TotNumPart * All.TreeDomainUpdateFrequency;

  Flag_FullStep = 1;		/* to ensure that Peano-Hilber order is done */

  domain_Decomposition();	/* do initial domain decomposition (gives equal numbers of particles) */

  ngb_treebuild();		/* will build tree */

   setup_smoothinglengths();


  TreeReconstructFlag = 1;

  /* at this point, the entropy variable normally contains the 
   * internal energy, read in from the initial conditions file, unless the file
   * explicitly signals that the initial conditions contain the entropy directly. 
   * Once the density has been computed, we can convert thermal energy to entropy.
   *
	#ifndef ISOTHERM_EQS
	  if(header.flag_entropy_instead_u == 0)
	    for(i = 0; i < N_gas; i++)
	      //SphP[i].Entropy = GAMMA_MINUS1 * SphP[i].Entropy / pow(SphP[i].Density / a3, GAMMA_MINUS1);
	#endif
  */
}


/*! This routine computes the mass content of the box and compares it to the
 *  specified value of Omega-matter.  If discrepant, the run is terminated.
 */
void check_omega(void)
{
  double mass = 0, masstot, omega;
  int i;

  for(i = 0; i < NumPart; i++)
    mass += P[i].Mass;

  MPI_Allreduce(&mass, &masstot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  omega =
    masstot / (All.BoxSize * All.BoxSize * All.BoxSize) / (3 * All.Hubble * All.Hubble / (8 * M_PI * All.G));

  if(fabs(omega - All.Omega0) > 1.0e-3)
  {
    if(ThisTask == 0)
	  {
	    printf("\n\nI've found something odd!\n");
	    printf
	    ("The mass content accounts only for Omega=%g,\nbut you specified Omega=%g in the parameterfile.\n",
	     omega, All.Omega0);
	    printf("\nI better stop.\n");

	    fflush(stdout);
	  }
    endrun(1);
  }
}



/***********************
 setup_smoothinglengths
***********************/
void setup_smoothinglengths(void)
{
  int i, no, p;

  //if(RestartFlag == 0)
  //{

    for(i = 0; i < N_gas; i++)
  	{
  	  no = Father[i];

  	  while(10 * All.DesNumNgb * P[i].Mass > Nodes[no].u.d.mass)
	    {
	      p = Nodes[no].u.d.father;

	      if(p < 0)
      		break;

	      no = p;
	    }
			#ifndef TWODIMS
	    	SphP[i].Hsml =
	      	     pow(3.0 / (4 * M_PI) * All.DesNumNgb * P[i].Mass / Nodes[no].u.d.mass, 1.0 / 3) * Nodes[no].len;
			#else
	    	SphP[i].Hsml =
	      	     pow(1.0 / (M_PI) * All.DesNumNgb * P[i].Mass / Nodes[no].u.d.mass, 1.0 / 2) * Nodes[no].len;
			#endif
	  }
  //}

  density();
}
