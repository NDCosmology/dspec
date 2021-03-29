#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <sys/types.h>
#include <unistd.h>
#include <gsl/gsl_rng.h>

#include "allvars.h"
#include "proto.h"


/*! \file begrun.c
 *  \brief initial set-up of a simulation run
 *
 *  This file contains various functions to initialize a simulation run. In
 *  particular, the parameterfile is read in and parsed, the initial
 *  conditions or restart files are read, and global variables are
 *  initialized to their proper values.
 */


/*! This function performs the initial set-up of the simulation. First, the
 *  parameterfile is set, then routines for setting units, reading
 *  ICs/restart-files are called, auxialiary memory is allocated, etc.
 */
void begrun(void)
{
  struct global_data_all_processes all;

  if(ThisTask == 0)
  {
    printf("\nThis is Gadget, version `%s'.\n", GADGETVERSION);
    printf("\nRunning on %d processors.\n", NTask);
  }

  read_parameter_file(ParameterFile);	/* ... read in parameters for this run */

  allocate_commbuffers();	/* ... allocate buffer-memory for particle 
	                      			   exchange during force computation */
  set_units();

	#if defined(PERIODIC) && (!defined(PMGRID) || defined(FORCETEST))
  	ewald_init();
	#endif

  //open_outputfiles();

  random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);
  gsl_rng_set(random_generator, 42);	/* start-up seed */

	#ifdef PMGRID
	  //long_range_init();
	#endif

  All.TimeLastRestartFile = CPUThisRun;

  if(RestartFlag == 0 || RestartFlag == 2)
  {
    set_random_numbers();

    init();			/* ... read in initial model */
  }
  else
  {
    all = All;		/* save global variables. (will be read from restart file) */

    /*restart(RestartFlag);	/ * ... read restart file. Note: This also resets 
				                         all variables in the struct `All'. However, during the
                                 run, some variables in the parameter file are allowed 
                                 to be changed, if desired. These need to copied in the
                                 way below.
				                         Note:  All.PartAllocFactor is treated in restart() 
                                        separately.  
				                  */

		//NOTE: May need to adjust for FeedbackRadius.
    All.MinSizeTimestep = all.MinSizeTimestep;
    All.MaxSizeTimestep = all.MaxSizeTimestep;
    All.BufferSize = all.BufferSize;
    All.BunchSizeForce = all.BunchSizeForce;
    All.BunchSizeDensity = all.BunchSizeDensity;
    All.BunchSizeHydro = all.BunchSizeHydro;
    All.BunchSizeDomain = all.BunchSizeDomain;

    All.TimeLimitCPU = all.TimeLimitCPU;
    All.ResubmitOn = all.ResubmitOn;
    All.TimeBetSnapshot = all.TimeBetSnapshot;
    All.TimeBetStatistics = all.TimeBetStatistics;
    All.CpuTimeBetRestartFile = all.CpuTimeBetRestartFile;
    All.ErrTolIntAccuracy = all.ErrTolIntAccuracy;
    All.MaxRMSDisplacementFac = all.MaxRMSDisplacementFac;

    All.ErrTolForceAcc = all.ErrTolForceAcc;

    All.TypeOfTimestepCriterion = all.TypeOfTimestepCriterion;
    All.TypeOfOpeningCriterion = all.TypeOfOpeningCriterion;
    All.NumFilesWrittenInParallel = all.NumFilesWrittenInParallel;
    All.TreeDomainUpdateFrequency = all.TreeDomainUpdateFrequency;

    All.SnapFormat = all.SnapFormat;
    All.NumFilesPerSnapshot = all.NumFilesPerSnapshot;
    All.MaxNumNgbDeviation = all.MaxNumNgbDeviation;
    All.ArtBulkViscConst = all.ArtBulkViscConst;


    All.OutputListOn = all.OutputListOn;
    All.CourantFac = all.CourantFac;

    All.OutputListLength = all.OutputListLength;
    memcpy(All.OutputListTimes, all.OutputListTimes, 
           sizeof(double) * All.OutputListLength);


    strcpy(All.ResubmitCommand, all.ResubmitCommand);
    strcpy(All.OutputListFilename, all.OutputListFilename);
    strcpy(All.OutputDir, all.OutputDir);
    strcpy(All.RestartFile, all.RestartFile);
    strcpy(All.EnergyFile, all.EnergyFile);
    strcpy(All.InfoFile, all.InfoFile);
    strcpy(All.CpuFile, all.CpuFile);
    strcpy(All.TimingsFile, all.TimingsFile);
    strcpy(All.SnapshotFileBase, all.SnapshotFileBase);

  }

  All.TimeLastRestartFile = CPUThisRun;
}




/*! Computes conversion factors between internal code units and the
 *  cgs-system.
 */
void set_units(void)
{
  double meanweight;

  All.UnitTime_in_s = All.UnitLength_in_cm / All.UnitVelocity_in_cm_per_s;
  All.UnitTime_in_Megayears = All.UnitTime_in_s / SEC_PER_MEGAYEAR;

  if(All.GravityConstantInternal == 0)
    All.G = GRAVITY / pow(All.UnitLength_in_cm, 3) * All.UnitMass_in_g * pow(All.UnitTime_in_s, 2);
  else
    All.G = All.GravityConstantInternal;

  All.UnitDensity_in_cgs = All.UnitMass_in_g / pow(All.UnitLength_in_cm, 3);
  All.UnitPressure_in_cgs = All.UnitMass_in_g / All.UnitLength_in_cm / pow(All.UnitTime_in_s, 2);
  All.UnitCoolingRate_in_cgs = All.UnitPressure_in_cgs / All.UnitTime_in_s;
  All.UnitEnergy_in_cgs = All.UnitMass_in_g * pow(All.UnitLength_in_cm, 2) / pow(All.UnitTime_in_s, 2);

  /* convert some physical input parameters to internal units */

  All.Hubble = HUBBLE * All.UnitTime_in_s;

  if(ThisTask == 0)
  {
    printf("\nHubble (internal units) = %g\n", All.Hubble);
    printf("G (internal units) = %g\n", All.G);
    printf("UnitMass_in_g = %g \n", All.UnitMass_in_g);
    printf("UnitTime_in_s = %g \n", All.UnitTime_in_s);
    printf("UnitVelocity_in_cm_per_s = %g \n", All.UnitVelocity_in_cm_per_s);
    printf("UnitDensity_in_cgs = %g \n", All.UnitDensity_in_cgs);
    printf("UnitEnergy_in_cgs = %g \n", All.UnitEnergy_in_cgs);
    printf("\n");
  }

	/*Photo-ionized even at low temps. Self shielding ignored. See
		p1214 in Schaye & Dalla Vecchia 2008 */
	#ifndef UVBACKGROUND
 		if(All.InitGasTemp > 1.0e4)   /* assuming FULL ionization */
 		  meanweight = 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));
		else                          /* assuming NEUTRAL GAS */
		  meanweight = 4 / (1 + 3 * HYDROGEN_MASSFRAC);
	#endif
	#ifdef UVBACKGROUND
 		meanweight = 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));
	#endif		
	

 // meanweight = 4.0 / (1 + 3 * HYDROGEN_MASSFRAC);	/* note: we assume neutral gas here */

	#ifdef ISOTHERM_EQS
	  All.MinEgySpec = 0;
	#else
	  All.MinEgySpec = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * All.MinGasTemp;
	  All.MinEgySpec *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;
	#endif

}



/*!  This function closes the global log-files.
 */
void close_outputfiles(void)
{
  if(ThisTask != 0)		/* only the root processor writes to the log files */
    return;

  fclose(FdCPU);
  fclose(FdInfo);
  fclose(FdEnergy);
  fclose(FdTimings);
	#ifdef FORCETEST
	  fclose(FdForceTest);
	#endif
}




/*! This function parses the parameterfile in a simple way.  Each paramater
 *  is defined by a keyword (`tag'), and can be either of type double, int,
 *  or character string.  The routine makes sure that each parameter
 *  appears exactly once in the parameterfile, otherwise error messages are
 *  produced that complain about the missing parameters.
 */
void read_parameter_file(char *fname)
{
	#define DOUBLE 1
	#define STRING 2
	#define INT 3
	#define MAXTAGS 300               //Maximum number of parameters

  int i, nt;
  char tag[MAXTAGS][50];						//Parameter field, ex: MinGasTemp
  int  errorFlag = 0;


  if(sizeof(long long) != 8)
  {
    if(ThisTask == 0)
	    printf("\nType `long long' is not 64 bit on this platform. Stopping.\n\n");
    endrun(0);
  }

  if(sizeof(int) != 4)
  {
    if(ThisTask == 0)
    	printf("\nType `int' is not 32 bit on this platform. Stopping.\n\n");
    endrun(0);
  }

  if(sizeof(float) != 4)
  {
    if(ThisTask == 0)
    	printf("\nType `float' is not 32 bit on this platform. Stopping.\n\n");
    endrun(0);
  }

  if(sizeof(double) != 8)
  {
    if(ThisTask == 0)
    	printf("\nType `double' is not 64 bit on this platform. Stopping.\n\n");
    endrun(0);
  }


  if(ThisTask == 0)		/* read parameter file on process 0 */
  {
    nt = 0;

    strcpy(All.InitCondFile, snap_filename);
    strcpy(All.OutputDir, ".");
    All.DesNumNgb       = 48;
    All.MaxNumNgbDeviation = 2;
    All.ICFormat        = 1;
    All.SnapFormat      = 1;
    All.NumFilesPerSnapshot = 1;
    All.BufferSize      = 30;
    All.PartAllocFactor = 3.0;
    All.TreeAllocFactor = 0.8;
    All.NumFilesWrittenInParallel = 1;
    All.PeriodicBoundariesOn = 1;
    All.MinHsml = 0;
    All.InitGasTemp = 1000;   //Used for read_ic.c conditional
    All.GravityConstantInternal = 0;  // I didn't see this explicitly set anywhere else,
                                       // so even though it's probably fine, b/c it probably
                                       // just gets set to 0 anyway, I decided to be explicit
                                       // about it.
    //All.ComovingIntegrationOn = 1;
    // Set units explicitly here because they aren't
    All.UnitMass_in_g = 1.989e43;
    All.UnitLength_in_cm = 3.085678e21;
    All.UnitVelocity_in_cm_per_s = 1e5;

    //Other fields of All.* are filled in read_ic.c


    if(errorFlag != 2)
	    for(i = 0; i < nt; i++)
	    {
	      if(*tag[i])
	      {
	      	printf("Error. I miss a value for tag '%s' in parameter file '%s'.\n", 
                 tag[i], fname);
		      errorFlag = 1;
	      }
	    }

    /*
    if(All.OutputListOn && errorFlag == 0)
	    errorFlag += read_outputlist(All.OutputListFilename);
    else
	    All.OutputListLength = 0;
    */
  }


  MPI_Bcast(&errorFlag, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if(errorFlag)
  {
    MPI_Finalize();
    exit(0);
  }

  /* now communicate the relevant parameters to the other processes */
  MPI_Bcast(&All, sizeof(struct global_data_all_processes), MPI_BYTE, 0, 
            MPI_COMM_WORLD);


  if(All.NumFilesWrittenInParallel < 1)
  {
    if(ThisTask == 0)
    	printf("NumFilesWrittenInParallel MUST be at least 1\n");
    endrun(0);
  }

  if(All.NumFilesWrittenInParallel > NTask)
  {
    if(ThisTask == 0)
    	printf("NumFilesWrittenInParallel MUST be smaller than number of processors\n");
    endrun(0);
  }

	#ifdef PERIODIC
	  if(All.PeriodicBoundariesOn == 0)
	  {
	    if(ThisTask == 0)
	  	{
	  	  printf("Code was compiled with periodic boundary conditions switched on.\n");
	  	  printf("You must set `PeriodicBoundariesOn=1', or recompile the code.\n");
	  	}
	    endrun(0);
	  }
	#else
	  if(All.PeriodicBoundariesOn == 1)
	  {
	    if(ThisTask == 0)
		  {
		    printf("Code was compiled with periodic boundary conditions switched off.\n");
		    printf("You must set `PeriodicBoundariesOn=0', or recompile the code.\n");
		  }
	    endrun(0);
	  }
	#endif


  if(All.TypeOfTimestepCriterion >= 1)
  {
    if(ThisTask == 0)
	  {
	    printf("The specified timestep criterion\n");
	    printf("is not valid\n");
	  }
    endrun(0);
  }

	#if defined(LONG_X) ||  defined(LONG_Y) || defined(LONG_Z)
		#ifndef NOGRAVITY
  		if(ThisTask == 0)
  		{
    		printf("Code was compiled with LONG_X/Y/Z, but not with NOGRAVITY.\n");
    		printf("Stretched periodic boxes are not implemented for gravity yet.\n");
  		}
  		endrun(0);
		#endif
	#endif

	#undef DOUBLE
	#undef STRING
	#undef INT
	#undef MAXTAGS
}





