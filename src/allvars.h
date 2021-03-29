/*! \file allvars.h
 *  \brief declares global variables.
 *
 *  This file declares all global variables. Further variables should be added here,
 *  and declared as 'extern'. The actual existence of these variables is provided by 
 *  the file 'allvars.c'. To produce 'allvars.c' from 'allvars.h', do the following:
 *
 *     - Erase all #define's, typedef's, and enum's
 *     - add #include "allvars.h", delete the #ifndef ALLVARS_H conditional
 *     - delete all keywords 'extern'
 *     - delete all struct definitions enclosed in {...}, e.g.
 *        "extern struct global_data_all_processes {....} All;"
 *        becomes "struct global_data_all_processes All;"

	K04 : Kobayashi 2004

 */

#ifndef ALLVARS_H
#define ALLVARS_H

#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <mpi.h>
#include "tags.h"

#define  GADGETVERSION   "2.0"   /*!< code version string */

#define  TIMEBASE        (1<<28) /*!< The simulated timespan is mapped onto the 
                                      integer interval [0,TIMESPAN], where TIMESPAN 
                                      needs to be a power of 2. Note that (1<<28) 
                                      corresponds to 2^28 */

#define  MAXTOPNODES     200000  /*!< Maximum number of nodes in the top-level tree 
                                      used for domain decomposition */


typedef  long long  peanokey;    /*!< defines the variable type used for Peano-Hilbert 
                                      keys */

#define  BITS_PER_DIMENSION 18	 /*!< Bits per dimension available for Peano-Hilbert 
                                      order. Note: If peanokey is defined as type int, 
                                      the allowed maximum is 10. If 64-bit integers are
                                      used, the maximum is 21 */

#define  PEANOCELLS (((peanokey)1)<<(3*BITS_PER_DIMENSION))  /*!< The number of 
                                                                 different Peano-Hilbert
                                                                 cells */


#define  RNDTABLE         3000   /*!< gives the length of a table with random numbers, 
                                      refreshed at every timestep. This is used to 
                                      allow application of random numbers to a specific
                                      particle in a way that is independent of the 
                                      number of processors used. */
#define  MAX_REAL_NUMBER  1e37
#define  MIN_REAL_NUMBER  1e-37

#define  MAXLEN_FILENAME  300    /*!< Maximum number of characters for filenames 
                                      (including the full path) */

#ifdef   ISOTHERM_EQS
	#define  GAMMA         (1.0)     /*!< index for isothermal gas */
#else
	#define  GAMMA         (5.0/3)   /*!< adiabatic index of simulated gas */
#endif

#define  GAMMA_MINUS1  (GAMMA-1)






/* Below Mass fractions are NOT independent. If you change one, you need to check
	 and change the others accordingly. See Solar_Abundances.txt and 
	 Solar_Abundances.xlsx for method
 */
#define  HYDROGEN_MASSFRAC         0.76        /* mass fraction of hydrogen, 
                                                 relevant only for radiative cooling */
#define  SOLAR_METALMASSFRAC       0.013765	/* = (M_metal / M_total)_sun Unitless. 
                                                 Derive by using Table 7.1 in  
                                                 Hazy : A brief introduction to CLOUDY 
                                                 C10. See Solar_Abundances.txt */
#define  PRIMORDIAL_METALMASSFRAC  1.3957e-5  /* = (M_metal /M_total)_primordial. 
																								 Corresponds to 
																								 log10(metal_frac/solar_metal_frac) = 
																								 -2.994. See Solar_Abundances.txt and
																									Solar_Abundances.xlsx*/
#define  INIT_METAL_TO_H_RATIO    1.0001e-3	/* Equivalent to "metals = -3" command
																							 in CLOUDY, the initial metal mass ratio
																							 to solar metal mass. Ensure w/in 
																							 interpolation bounds. Not quite the same
																							 as real metallicity.*/

#define IMF_MAX 100.0						//upper limit in stellar initial mass function
#define IMF_MIN 0.1							//lower limit in stellar initial mass function
#define NELEMENTS 25							/* Number of elements in SN yield tables.
																	 	 This is checked in begrun.c */
#define MAX_COOLING_ITER 100				/* Maximum number of iterations used in implicit 
																		 cooling integration */




/* Some physical constants in cgs units */

#define  GRAVITY                 6.672e-8   /* Gravitational constant (in cgs units) */
#define  SOLAR_MASS              1.989e33
#define  SOLAR_LUM               3.826e33
#define  RAD_CONST               7.565e-15
#define  AVOGADRO                6.0222e23
#define  BOLTZMANN               1.3806e-16			// erg/K
#define  GAS_CONST               8.31425e7
#define  C                       2.9979e10
#define  PLANCK                  6.6262e-27
#define  CM_PER_MPC              3.085678e24		
#define  CM_PER_PC               3.085678e18		
#define  PROTONMASS              1.6726e-24			// g
#define  ELECTRONMASS            9.10953e-28		// g
#define  THOMPSON                6.65245e-25
#define  ELECTRONCHARGE          4.8032e-10
#define  HUBBLE                  3.2407789e-18	/* in h/sec */

/* Some conversion factors */

#define  SEC_PER_MEGAYEAR  3.1557e13
#define  SEC_PER_YEAR      3.1557e7

#ifndef ASMTH
#define ASMTH 1.25                /*!< ASMTH gives the scale of the short-range/long-
                                       range force split in units of FFT-mesh cells */
#endif

#ifndef RCUT
#define RCUT  4.5                 /*!< RCUT gives the maximum distance (in units of the
                                       scale used for the force split) out to which 
                                       short-range forces are evaluated in the short-
                                       range tree walk. THis is 'r_s'in eqns. 20 & 21 */
#endif

#define MAX_NGB             20000  /*!< defines maximum length of neighbour list */

#define MAXLEN_OUTPUTLIST   500	   /*!< maxmimum number of entries in list of snapshot 
                                        output times */

#define DRIFT_TABLE_LENGTH  1000   /*!< length of the lookup table used to hold the 
                                        drift and kick factors */ 

#define MAXITER             150    /*!< maxmimum number of steps for SPH neighbour 
                                        iteration */


#ifdef DOUBLEPRECISION             /*!< If defined, the variable type FLOAT is set to 
                                        "double", otherwise to FLOAT */
	#define FLOAT double
#else
	#define FLOAT float
#endif

/* Below are coefficients for SPH spline kernel and its derivatives.
	 See eqn. 4 in GADGET2 paper.
*/
#ifndef  TWODIMS
#define  NUMDIMS 3                              /*!< For 3D-normalized kernel */
#define  KERNEL_COEFF_1  2.546479089470         // = 8/PI
#define  KERNEL_COEFF_2  15.278874536822				// = 8*6/PI
#define  KERNEL_COEFF_3  45.836623610466				// = 8*6*3/PI
#define  KERNEL_COEFF_4  30.557749073644				// = 8*6*2/PI
#define  KERNEL_COEFF_5  5.092958178941					// = 8*2/PI
#define  KERNEL_COEFF_6  (-15.278874536822)			// = -8*6/PI
#define  NORM_COEFF      4.188790204786         /* Coeff. for kernel normalization. 
																									 Note: 4.0/3 * PI = 4.188790204786 */ 
#else
#define  NUMDIMS 2                                 /*!< For 2D-normalized kernel */
#define  KERNEL_COEFF_1  (5.0/7*2.546479089470)    /*!< Coefficients for SPH spline 
                                                        kernel and its derivative */ 
#define  KERNEL_COEFF_2  (5.0/7*15.278874536822)
#define  KERNEL_COEFF_3  (5.0/7*45.836623610466)
#define  KERNEL_COEFF_4  (5.0/7*30.557749073644)
#define  KERNEL_COEFF_5  (5.0/7*5.092958178941)
#define  KERNEL_COEFF_6  (5.0/7*(-15.278874536822))
#define  NORM_COEFF      M_PI                      /*!< Coefficient for kernel 
                                                        normalization. */
#endif



extern int ThisTask;		        /*!< the rank of the local processor */
extern int NTask;               /*!< number of processors */
extern int PTask;	              /*!< smallest integer such that NTask <= 2^PTask */

extern int NumPart;		          /*!< number of particles on the LOCAL processor */
extern int N_gas;		            /*!< number of gas particles on the LOCAL processor  */
#ifdef SFR
extern int N_star;							/*!< number of star particles on the LOCAL processor */
extern int N_star_tot;					/*!< ALI: Global Number of star Particles */
#endif
extern long long Ntype[6];      /*!< total number of particles of each type */
extern int NtypeLocal[6];       /*!< local number of particles of each type 
																		 ALI: Note that this is not updated at the same
																					time as N_gas and NumPart. I would be careful
																					about using this */

extern int NumForceUpdate;      /*!< number of active particles on local processor in 
                                     current timestep  */
extern int NumSphUpdate;        /*!< number of active SPH particles on local processor 
                                     in current timestep  */

extern int NumGridUpdate;        /*!< number of active Star particles on local 
																			processor in current timestep  */

extern double CPUThisRun;	      /*!< Sums the CPU time for the process (current 
                                     submission only) */


extern int RestartFlag;         /*!< taken from command line used to start code. 0 is 
                                     normal start-up from initial conditions, 1 is 
                                     resuming a run from a set of restart files, while 2
                                     marks a restart from a snapshot file. */

extern char *Exportflag;        /*!< Buffer used for flagging whether a particle needs 
                                     to be exported to another process. I believe this
																		 is used to flag a processor which has a particle
																		 to export. */

extern int  *Ngblist;           /*!< Buffer to hold indices of neighbours retrieved by 
                                     the neighbour search routines */

extern int TreeReconstructFlag; /*!< Signals that a new tree needs to be constructed */

extern int Flag_FullStep;       /*!< This flag signals that the current step involves 
                                     all particles */


extern gsl_rng *random_generator; /*!< the employed random number generator of the GSL 
                                       library */

extern double RndTable[RNDTABLE]; /*!< Hold a table with random numbers, refreshed 
                                       every timestep */


extern double DomainCorner[3];    /*!< gives the lower left corner of simulation 
                                       volume */
extern double DomainCenter[3];    /*!< gives the center of simulation volume */
extern double DomainLen;          /*!< gives the (maximum) side-length of simulation 
                                       volume */
extern double DomainFac;          /*!< factor used for converting particle coordinates
                                       to a Peano-Hilbert mesh covering the simulation 
                                       volume */
extern int    DomainMyStart;      /*!< first domain mesh cell that resides on the local
                                       processor */
extern int    DomainMyLast;       /*!< last domain mesh cell that resides on the local 
                                       processor */
extern int    *DomainStartList;   /*!< a table that lists the first domain mesh cell 
                                       for all processors */
extern int    *DomainEndList;     /*!< a table that lists the last domain mesh cell for
                                       all processors */
extern double *DomainWork;        /*!< a table that gives the total "work" due to the 
                                       particles stored by each processor */
extern int    *DomainCount;       /*!< a table that gives the total number of particles
                                       held by each processor */
extern int    *DomainCountSph;    /*!< a table that gives the total number of SPH 
                                       particles held by each processor */

extern int    *DomainTask;        /*!< this table gives for each leaf of the top-level
                                       tree the processor it was assigned to */
extern int    *DomainNodeIndex;   /*!< this table gives for each leaf of the top-level
                                       tree the corresponding node of the gravitational
                                       tree */
extern FLOAT  *DomainTreeNodeLen; /*!< this table gives for each leaf of the top-level
                                       tree the side-length of the corresponding node of
                                       the gravitational tree */
extern FLOAT  *DomainHmax;        /*!< this table gives for each leaf of the top-level
                                       tree the maximum SPH smoothing length among the
                                       particles of the corresponding node of the 
                                       gravitational tree */

extern struct DomainNODE
{
  FLOAT s[3];                     /*!< center-of-mass coordinates */
  FLOAT vs[3];                    /*!< center-of-mass velocities */
  FLOAT mass;                     /*!< mass of node */
#ifdef UNEQUALSOFTENINGS
#ifndef ADAPTIVE_GRAVSOFT_FORGAS
  int   bitflags;                 /*!< this bit-field encodes the particle type with 
                                       the largest softening among the particles of the
                                       nodes, and whether there are particles with
                                       different softening in the node */
#else
  FLOAT maxsoft;                  /*!< hold the maximum gravitational softening of 
                                       particles in the node if the 
                                       ADAPTIVE_GRAVSOFT_FORGAS option is selected */
#endif
#endif
}
 *DomainMoment;                   /*!< this table stores for each node of the top-level
                                       tree corresponding node data from the 
                                       gravitational tree */

extern peanokey *DomainKeyBuf;    /*!< this points to a buffer used during the 
                                       exchange of particle data */

extern peanokey *Key;             /*!< a table used for storing Peano-Hilbert keys for
                                       particles */
extern peanokey *KeySorted;       /*!< holds a sorted table of Peano-Hilbert keys for
                                       all particles, used to construct top-level tree*/


extern int NTopnodes;             /*!< total number of nodes in top-level tree */
extern int NTopleaves;            /*!< number of leaves in top-level tree. Each leaf 
                                       can be assigned to a different processor */

extern struct topnode_data
{
  int Daughter;                   /*!< index of first daughter cell (out of 8) of 
                                       top-level node */
  int Pstart;                     /*!< for the present top-level node, this gives the 
                                       index of the first node in the concatenated list
                                       of topnodes collected from all processors */
  int Blocks;                     /*!< for the present top-level node, this gives the 
                                       number of corresponding nodes in the concatenated
                                       list of topnodes collected from all processors */
  int Leaf;                       /*!< if the node is a leaf, this gives its number when
                                       all leaves are traversed in Peano-Hilbert order
                                       ALI: I think this is the index used to map the 
                                            topnode to DomainNodeIndex[] which itself
                                            maps to Nodes[] */
  peanokey Size;                  /*!< number of Peano-Hilbert mesh-cells represented
                                       by top-level node */
  peanokey StartKey;              /*!< first Peano-Hilbert key in top-level node */
  long long Count;                /*!< counts the number of particles in this top-level
                                       node */
}
 *TopNodes;                       /*!< points to the root node of the top-level tree */


extern double TimeOfLastTreeConstruction; /*!< holds what it says, only used in 
                                               connection with FORCETEST */



/* variables for input/output, usually only used on process 0 */

extern char ParameterFile[MAXLEN_FILENAME];  /*!< file name of parameterfile used for 
                                                  starting the simulation */

extern FILE *FdInfo;       /*!< file handle for info.txt log-file. */
extern FILE *FdEnergy;     /*!< file handle for energy.txt log-file. */
extern FILE *FdTimings;    /*!< file handle for timings.txt log-file. */
extern FILE *FdCPU;        /*!< file handle for cpu.txt log-file. */

#ifdef FORCETEST
extern FILE *FdForceTest;  /*!< file handle for forcetest.txt log-file. */
#endif


extern double DriftTable[DRIFT_TABLE_LENGTH];     /*!< table for the cosmological 
                                                        drift factors */
extern double GravKickTable[DRIFT_TABLE_LENGTH];  /*!< table for the cosmological kick
                                                       factor for gravitational forces*/
extern double HydroKickTable[DRIFT_TABLE_LENGTH]; /*!< table for the cosmological kick 
                                                       factor for hydrodynmical forces*/

extern void *CommBuffer;   /*!< points to communication buffer, which is used in the 
                                domain decomposition, the parallel tree-force 
                                computation, the SPH routines, etc. */



/*! This structure contains data which is the SAME for all tasks (mostly code 
    parameters read from the parameter file).  Holding this data in a structure is 
    convenient for writing/reading the restart file, and it allows the introduction of 
    new global variables in a simple way. The only thing to do is to introduce
    them into this structure.
 */
extern struct global_data_all_processes
{
  long long TotNumPart;		/*!< total particle numbers (global value) */
  long long TotN_gas;  	  /*!< total gas particle number (global value) */
	#ifdef COOLING
	#ifdef SFR
		long long TotN_star;		/*!< ALI: total star particle number (global value) */	
	#endif
	#endif

  int MaxPart;			      /*!< This gives the maxmimum number of particles that can 
                               be stored on one processor. */
  int MaxPartSph;		      /*!< This gives the maxmimum number of SPH particles that can
                               be stored on one processor. */

  double BoxSize;         /*!< Boxsize in case periodic boundary conditions are used */

  int ICFormat;			      /*!< selects different versions of IC file-format */

  int SnapFormat;		      /*!< selects different versions of snapshot file-formats */

  int NumFilesPerSnapshot;      /*!< number of files in multi-file snapshot dumps */
  int NumFilesWrittenInParallel;/*!< maximum number of files that may be written 
                                     simultaneously when writing/reading restart-files,
                                     or when writing snapshot files */ 

  int BufferSize;		            /*!< size of communication buffer in MB */
  int BunchSizeForce;		        /*!< number of particles fitting into the buffer in the
                                     parallel tree-force algorithm  */
  int BunchSizeDensity;         /*!< number of particles fitting into the communication
                                     buffer in the density computation */
  int BunchSizeHydro;           /*!< number of particles fitting into the communication
                                     buffer in the SPH hydrodynamical force computation
                                     
																		 NOTE: That this is NOT the full length of the 
																		 CommBuffer.
																 */
		int BunchSizeGasSearch;			/* number of particles fitting into the communication
																	 buffer in the gas search compuation.
																	 NOTE: This is NOT the full length of CommBuffer
																*/
		int BunchSizeFeedback;			/* number of particles fitting into the communication
																	 buffer in the feedback compuation.
																	 NOTE: This is NOT the full length of CommBuffer
																*/
			int BunchSizeNormalization; /* Store Normalization constant, eqn 6 Stinson et al
																		 (2006) 
																	*/

  int BunchSizeDomain;          /*!< number of particles fitting into the communication
                                     buffer in the domain decomposition */

  double PartAllocFactor;	      /*!< in order to maintain work-load balance, the 
                                     particle load will usually NOT be balanced.  Each 
                                     processor allocates memory for PartAllocFactor 
                                     times the average number of particles to allow for
                                     that */

  double TreeAllocFactor;	      /*!< Each processor allocates a number of nodes which is
                                     TreeAllocFactor times the maximum(!) number of 
                                     particles.  Note: A typical local tree for N
				                             particles needs usually about ~0.65*N nodes. */

  /* some SPH parameters */

  double DesNumNgb;             /*!< Desired number of SPH neighbours */
  double MaxNumNgbDeviation;    /*!< Maximum allowed deviation neighbour number */

  double ArtBulkViscConst;      /*!< Sets the parameter $\alpha\$ of the artificial 
                                     viscosity. See Eqn 14 */
  double InitGasTemp;		        /*!< may be used to set the temperature in the IC's */
  double MinGasTemp;		        /*!< may be used to set a floor for the gas 
                                     temperature */
  double MinEgySpec;            /*!< the minimum allowed temperature expressed as energy
                                     per unit mass */

	double StarDesNumNgb;						//Desired number of gas neighbors for each star
	double StarMaxNumNgbDeviation;	//Max deviation from the desired Ngb


  /* some force counters  */

  long long TotNumOfForces;	    /*!< counts total number of force computations  */
  long long NumForcesSinceLastDomainDecomp;  /*!< count particle updates since last 
                                                  domain decomposition */


  /* system of units  */

  double G;                    /*!< Gravity-constant in internal units in 
                                    (kpc cm^2 / M_solar s^2)  */
  double UnitTime_in_s;   	   /*!< factor to convert internal time unit to seconds/h*/
  double UnitMass_in_g;        /*!< factor to convert internal mass unit to grams/h */
  double UnitVelocity_in_cm_per_s; /*!< factor to convert intqernal velocity unit to 
                                        cm/sec */
  double UnitLength_in_cm;     /*!< factor to convert internal length unit to cm/h */
  double UnitPressure_in_cgs;  /*!< factor to convert internal pressure unit to 
																		g*h^2/(cm * s^2)  (little 'h' still around!) */
  double UnitDensity_in_cgs;   /*!< factor to convert internal length unit to 
                                    g*h^2/cm^3  */
  double UnitCoolingRate_in_cgs;  /*!< factor to convert internal cooling rate to cgs 
                                        ergs*h^3/(cm^3 * s) units */
  double UnitEnergy_in_cgs;       /*!< factor to convert internal energy ergs/h */
  double UnitTime_in_Megayears;   /*!< factor to convert internal time to megayears/h*/
  double GravityConstantInternal; /*!< If set to zero in the parameterfile, the 
                                       internal value of the gravitational constant is
                                       set to the Newtonian value based on the system 
                                       of units specified.Otherwise the value provided
                                       is taken as internal gravity constant G. */


  /* Cosmological parameters */

  double Hubble;       /*!< Hubble-const @ z=0 (1/internal time) */
  double Omega0;       /*!< matter density in units of the critical density (at z=0)*/
  double OmegaLambda;  /*!< vaccum energy density relative to crictical density (at z=0) */
  double OmegaBaryon;  /*!< baryon density in units of the critical density (at z=0)*/
  double HubbleParam;  /*!< little `h', i.e. Hubble constant in units of 100 km/s/Mpc.
                            Only needed to get absolute physical values for cooling 
                            physics. I ignore it when ComovingIntegrationOn = 0 b/c it 
														is only useful in cosmological simulations */
  

  /* Code options */

  int ComovingIntegrationOn;  	/*!< flags that comoving integration is enabled. 
																		 When =0 I ignore 'h'. It doesn't make a lot
																		 of sense to use 'h' if the Universe isn't 
																		 expanding. */
  int PeriodicBoundariesOn;     /*!< flags that periodic boundaries are enabled */
  int ResubmitOn;               /*!< flags that automatic resubmission of job to queue 
                                     system is enabled */
  int TypeOfOpeningCriterion;   /*!< determines tree cell-opening criterion: 0 for 
                                     Barnes-Hut, 1 for relative criterion */
  int TypeOfTimestepCriterion;  /*!< gives type of timestep criterion (only 0 supported
                                     right now - unlike gadget-1.1) */
  int OutputListOn;             /*!< flags that output times are listed in a specified 
                                     file */


  /* Parameters determining output frequency */

  int SnapshotFileCount;        /*!< number of snapshot that is written next */
  double TimeBetSnapshot;       /*!< simulation time interval between snapshot files */
  double TimeOfFirstSnapshot;   /*!< simulation time of first snapshot files */
  double CpuTimeBetRestartFile; /*!< cpu-time between regularly generated restart 
                                     files */
  double TimeLastRestartFile;   /*!< cpu-time when last restart-file was written */
  double TimeBetStatistics;     /*!< simulation time interval between computations of 
                                     energy statistics */
  double TimeLastStatistics;    /*!< simulation time when the energy statistics was 
                                     computed the last time */
  int NumCurrentTiStep;         /*!< counts the number of system steps taken up to this
                                     point */


  /* Current time of the simulation, global step, and end of simulation */

  double Time;                  /*!< = 1/(1+z) current time of the simulation and 
																		 scale factor */
  double TimeBegin;             /*!< time of initial conditions of the simulation */
  double TimeStep;              /*!< difference between current times of previous and 
                                     current timestep */
  double TimeMax;	              /*!< marks the point of time until the simulation is to
                                     be evolved */


  /* variables for organizing discrete timeline */

  double Timebase_interval;     /*!< factor to convert from floating point time 
                                     interval to integer timeline. 
                                    = (log(All.TimeMax) - log(All.TimeMin))/TIMEBASE */
  int Ti_Current;               /*!< current time on integer timeline [0,TIMEBASE]*/ 
  int Ti_nextoutput;            /*!< next output time on integer timeline */
#ifdef FLEXSTEPS
  int PresentMinStep;           /*!< If FLEXSTEPS is used, particle timesteps are chosen
                                     as multiples of the present minimum timestep. */
  int PresentMaxStep;		        /*!< If FLEXSTEPS is used, this is the maximum timestep
                                     in timeline units, rounded down to the next power 
                                     2 division */
#endif
#ifdef PMGRID
  int PM_Ti_endstep;            /*!< begin of present long-range timestep */
  int PM_Ti_begstep;            /*!< end of present long-range timestep */
#endif


  /* Placement of PM grids */

#ifdef PMGRID
  double Asmth[2];              /*!< Gives the scale of the long-range/short-range 
                                     split (in mesh-cells), both for the coarse and the
                                     high-res mesh */
  double Rcut[2];               /*!< Gives the maximum radius for which the short-range
                                     force is evaluated with the tree (in mesh-cells), 
                                     both for the coarse and the high-res mesh */
  double Corner[2][3];          /*!< lower left corner of coarse and high-res PM-mesh */
  double UpperCorner[2][3];     /*!< upper right corner of coarse and high-res PM-mesh*/
  double Xmintot[2][3];         /*!< minimum particle coordinates both for coarse and 
                                     high-res PM-mesh */
  double Xmaxtot[2][3];         /*!< maximum particle coordinates both for coarse and 
                                     high-res PM-mesh */
  double TotalMeshSize[2];      /*!< total extension of coarse and high-res PM-mesh */
#endif


  /* Variables that keep track of cumulative CPU consumption */

  double TimeLimitCPU;          /*!< CPU time limit as defined in parameterfile */
  double CPU_TreeConstruction;  /*!< time spent for constructing the gravitational 
                                     tree */
  double CPU_TreeWalk;          /*!< actual time spent for pure tree-walks */
  double CPU_Gravity;           /*!< cumulative time used for gravity computation 
                                     (tree-algorithm only) */
  double CPU_Potential;         /*!< time used for computing gravitational potentials */
  double CPU_Domain;            /*!< cumulative time spent for domain decomposition */
  double CPU_Snapshot;          /*!< time used for writing snapshot files */
  double CPU_Total;             /*!< cumulative time spent for domain decomposition */
  double CPU_CommSum;           /*!< accumulated time used for communication, and for 
                                     collecting partial results, in tree-gravity */
  double CPU_Imbalance;         /*!< cumulative time lost accross all processors as work
                                     -load imbalance in gravitational tree */
  double CPU_HydCompWalk;       /*!< time used for actual SPH computations, including 
                                     neighbour search */
  double CPU_HydCommSumm;       /*!< cumulative time used for communication in SPH, and
                                     for collecting partial results */
  double CPU_HydImbalance;      /*!< cumulative time lost due to work-load imbalance in
                                     SPH */
	#ifdef FEEDBACK
		double CPU_FeedCompWalk;		/* These aren't used. Added Just in case */
		double CPU_FeedCommSumm;
		double CPU_FeedImbalance;
	
		#ifdef FEEDBACK_W_SPH_KERNAL
			double CPU_NormCompWalk;		/* These aren't used. Added Just in case */
			double CPU_NormCommSumm;
			double CPU_NormImbalance;
		#endif
	#endif 

  double CPU_Hydro;             /*!< cumulative time spent for SPH related 
                                     computations */
  double CPU_EnsureNgb;         /*!< time needed to iterate on correct neighbour 
                                     numbers */
  double CPU_Predict;           /*!< cumulative time to drift the system forward in 
                                     time, including dynamic tree updates */
  double CPU_TimeLine;          /*!< time used for determining new timesteps, and for 
                                     organizing the timestepping, including kicks of 
                                     active particles */
  double CPU_PM;                /*!< time used for long-range gravitational force */
  double CPU_Peano;             /*!< time required to establish Peano-Hilbert order */

  /* tree code opening criterion */

  double ErrTolTheta;		        /*!< BH tree opening angle */
  double ErrTolForceAcc;	      /*!< parameter for relative opening criterion in tree 
                                     walk */


  /* adjusts accuracy of time-integration */

  double ErrTolIntAccuracy;	    /*!< accuracy tolerance parameter \f$ \eta \f$ for 
                                     timestep criterion. The timestep is \f$ \Delta t =
                                     \sqrt{\frac{2 \eta eps}{a}} \f$ */

  double MinSizeTimestep;       /*!< minimum allowed timestep. Normally, the simulation
                                     terminates if the timestep determined by the 
                                     timestep criteria falls below this limit. */ 
  double MaxSizeTimestep;       /*!< maximum allowed timestep */

  double MaxRMSDisplacementFac; /*!< this determines a global timestep criterion for 
                                     cosmological simulations in comoving coordinates.
                                     To this end, the code computes the rms velocity
                                     of all particles, and limits the timestep such 
                                     that the rms displacement is a fraction of the 
                                     mean particle separation (determined from the 
                                     particle mass and the cosmological parameters). 
                                     This parameter specifies this fraction. */

  double CourantFac;		/*!< SPH-Courant factor */


  /* frequency of tree reconstruction/domain decomposition */

  double TreeDomainUpdateFrequency; /*!< controls frequency of domain decompositions  */


  /* Gravitational and hydrodynamical softening lengths (given in terms of an 
     `equivalent' Plummer softening length). Five groups of particles are supported 
      0="gas", 1="halo", 2="disk", 3="bulge", 4="stars", 5="bndry"
   */

  double MinGasHsmlFractional;  /*!< minimum allowed SPH smoothing length in units of 
                                     SPH gravitational softening length */
  double MinGasHsml;            /*!< minimum allowed SPH smoothing length */

	double MinStarHsmlFractional; /* min. permitted Gas particles in units of star 
																	 softening length. */

	double MinHsml;   				/* Minimum star smoothing length at a given redshift */

  double SofteningGas;      /*!< comoving gravitational softening lengths for type 0 */ 
  double SofteningHalo;     /*!< comoving gravitational softening lengths for type 1 */ 
  double SofteningDisk;     /*!< comoving gravitational softening lengths for type 2 */ 
  double SofteningBulge;    /*!< comoving gravitational softening lengths for type 3 */ 
  double SofteningStars;    /*!< comoving gravitational softening lengths for type 4 */ 
  double SofteningBndry;    /*!< comoving gravitational softening lengths for type 5 */ 

  double SofteningGasMaxPhys;   /*!< maximum physical softening length for type 0 */ 
  double SofteningHaloMaxPhys;  /*!< maximum physical softening length for type 1 */ 
  double SofteningDiskMaxPhys;  /*!< maximum physical softening length for type 2 */ 
  double SofteningBulgeMaxPhys; /*!< maximum physical softening length for type 3 */ 
  double SofteningStarsMaxPhys; /*!< maximum physical softening length for type 4 */ 
  double SofteningBndryMaxPhys; /*!< maximum physical softening length for type 5 */ 

  double SofteningTable[6];     /*!< current (comoving) gravitational softening lengths
                                     for each particle type */
  double ForceSoftening[6];     /*!< the same, but multiplied by a factor 2.8 - at that
                                     scale the force is Newtonian */


  double MassTable[6];          /*!< Table with particle masses for particle types with
                                     equal mass. If particle masses are all equal for 
                                     one type, the corresponding entry in MassTable 
                                     is set to this value, allowing the size of the 
                                     snapshot files to be reduced. */
	double FeedbackRadius;				//in Physical kpc.
	double ScaleFactorEnd;				/* Only useful when ComovingIntegrationOn=0. When 
																	 COOLING and UVBACKGROUND are enabled we need a way 
																	 to match the linear time to the redshift */
	double CosmTimeEnd;						/* Maps the ScaleFactorEnd and All.TimeMax to the 
																	 cosmological time so that it can be used when the 
																	 UVBACKGROUND is enabled AND 
																	 ComovingIntegrationOn = 0 */
	double CosmTimeBeg; 					/* Maps the All.TimeBeg to the cosmological time
																	 so that it can be used when the UVBACKGROUND is 
																	 enabled AND ComovingIntegrationOn = 0*/


  /* some filenames */

  char InitCondFile[MAXLEN_FILENAME];          /*!< filename of initial conditions */
  char OutputDir[MAXLEN_FILENAME];             /*!< output directory of the code */
  char SnapshotFileBase[MAXLEN_FILENAME];      /*!< basename to construct the names of 
                                                    snapshotf files */
  char EnergyFile[MAXLEN_FILENAME];            /*!< name of file with energy 
                                                    statistics*/
  char CpuFile[MAXLEN_FILENAME];               /*!< name of file with cpu-time 
                                                    statistics */
  char InfoFile[MAXLEN_FILENAME];              /*!< name of log-file with a list of the
                                                    timesteps taken */
  char TimingsFile[MAXLEN_FILENAME];           /*!< name of file with performance 
                                                    metrics of gravitational tree 
                                                    algorithm */
  char RestartFile[MAXLEN_FILENAME];           /*!< basename of restart-files */
  char ResubmitCommand[MAXLEN_FILENAME];       /*!< name of script-file that will be 
                                                    executed for automatic restart */
  char OutputListFilename[MAXLEN_FILENAME];    /*!< name of file with list of desired 
                                                    output times */

  double OutputListTimes[MAXLEN_OUTPUTLIST];   /*!< table with desired output times.
                                                    Specified in text file, specifies
                                                    the snapshot file times. For ex
                                                    see outputs_lcdm_gas.txt*/
  int OutputListLength;                        /*!< number of output times stored in 
																										the table of desired output times*/

}
 All;                                          /*!< a container variable for global 
                                                    variables that are equal on all 
                                                    processors */



/*! This structure holds all the information that is
 *  stored for each particle of the simulation.
 *
 *  NOTE: I believe that P[].Vel[] are really the conjugate momenta (from Hamilton's
 *        eqns) divided by the particle mass. Be VERY VERY careful when thinking about
 *        it, b/c it is evolved via eqn 1, 27, and 28. He defines the conjugate 
 *        momenta on p1107 under eqn 1... 
 *
 *                      p_i = a^2 * m_i * xdot_i   (updated in eqn 28)
 *                      P[].Vel = a^2 * xdot_i
 */
extern struct particle_data
{
  FLOAT Pos[3];			       /*!< particle position at its current time */
  FLOAT Mass;			         /*!< particle mass, in 10^10 M_sun/h */
  FLOAT Vel[3];		       	 /*!< particle velocity at its current time */
  FLOAT GravAccel[3];		   /*!< particle acceleration due to gravity */
	#ifdef PMGRID
  	FLOAT GravPM[3];		   /*!< particle acceleration due to long-range PM gravity 
                                force*/
	#endif
	#ifdef FORCETEST
  	FLOAT GravAccelDirect[3];/*!< particle acceleration when computed with direct 
    	                            summation */
	#endif
  FLOAT Potential;		     /*!< gravitational potential */
  FLOAT OldAcc;			       /*!< magnitude of old gravitational force. Used in relative 
                                opening criterion */
	#ifndef LONGIDS
  	unsigned int ID;		     /*!< particle identifier */
	#else
  	unsigned long long ID;   /*!< particle identifier */
	#endif

  int Type;		             /*!< flags particle type.  0=gas, 1=halo, 2=disk, 3=bulge, 
                                4=stars, 5=bndry */
  int Ti_endstep;          /*!< marks end of current timestep of particle on integer 
                                timeline */ 
  int Ti_begstep;          /*!< marks start of current timestep of particle on integer 
                                timeline */
	#ifdef FLEXSTEPS
  	int FlexStepGrp;		   /*!< a random 'offset' on the timeline to create a smooth 
                                groouping of particles */
	#endif
  float GravCost;		       /*!< weight factor used for balancing the work-load */
	#ifdef PSEUDOSYMMETRIC
  	float AphysOld;        /*!< magnitude of acceleration in last timestep.Used to make
                                a first order prediction of the change of acceleration 
                                expected in the future, thereby allowing to guess 
                                whether a decrease/increase of the timestep should 
                                occur in the timestep that is started. */
	#endif
	
  	FLOAT StarAge;       		//  Physical age of universe (mega-years) when star formed.

			long NumGasNearby;		/*	Num of gas particles within Feedback Radius of Star*/
			long NumGasFed;				//  Use to check that feedback is working.
			FLOAT TurnOffMass;		/*	Mass of oldest stars. Defined in eq. 26 
																Kobayashi 2004 (K04). in M_sun NOT M_sun/h */
	
			FLOAT MetalMass;	/* 	This is the (M_metal/M_gas) where  
																M_metal has solar abundance (See CLOUDY documentation, 
																Hazy Prt. 1, Table 7.1). 
																								::::::NOTE:::::
																We assume that ratio of (n_He / n_H) stays the same 
																even though this is not true in stellar evolution.
																We can do this b/c He contributes relatively little to 
																the overall cooling rate. See Solar_Abundances.txt
																*/

				FLOAT DtMass;				//Change in Mass due to feedback
				FLOAT StarHsml;			//Star smoothing length
				FLOAT StarHsmlLeft;	//Lower Bound on StarHsml in search
				FLOAT StarHsmlRight;//Upper Bound on StarHsml in search
				FLOAT FeedbackNorm; //FeedbackNorm or GasDensity are same? Maybe eliminate one?
				FLOAT GasDensity;		//Estimated Gas Density at Star Location using kernal
				FLOAT GasPressure;	//Estimated Pressure at Star Location using kernal
				FLOAT BlastRadius;	//in comoving inten. Units,see eqn 9 in Stinson et al 2006
        FLOAT origStarMass; //Original star mass;
				double R_II;        //Used to check that if Cooling should be turned off.
				double R_Ia;        //Type Ia SN rate.
				double R_SW;        //Stellar Wind Rate.
        double R_AGB;       //AGB and sAGB Rate

			#ifdef CARBON
				FLOAT CarbonMass;
			#endif

			#ifdef NITROGEN
				FLOAT NitrogenMass;
			#endif

			#ifdef OXYGEN
				FLOAT OxygenMass;
			#endif

			#ifdef FLORINE
				FLOAT FlorineMass;
			#endif

			#ifdef NEON
				FLOAT NeonMass;
			#endif

			#ifdef SODIUM
				FLOAT SodiumMass;
			#endif

			#ifdef MAGNESIUM
				FLOAT MagnesiumMass;
			#endif

			#ifdef ALUMINUM
				FLOAT AluminumMass;
			#endif

			#ifdef SILICON
				FLOAT SiliconMass;
			#endif

			#ifdef PHOSPHORUS
				FLOAT PhosphorusMass;
			#endif

			#ifdef SULFUR
				FLOAT SulfurMass;
			#endif

			#ifdef CHLORINE
				FLOAT ChlorineMass;
			#endif

			#ifdef ARGON
				FLOAT ArgonMass;
			#endif

			#ifdef POTASSIUM
				FLOAT PotassiumMass;
			#endif

			#ifdef CALCIUM
				FLOAT CalciumMass;
			#endif

			#ifdef SCANDIUM
				FLOAT ScandiumMass;
			#endif

			#ifdef TITANIUM
				FLOAT TitaniumMass;
			#endif

			#ifdef VANADIUM
				FLOAT VanadiumMass;
			#endif

			#ifdef CHROMIUM
				FLOAT ChromiumMass;
			#endif

			#ifdef MANGANESE
				FLOAT ManganeseMass;
			#endif

			#ifdef IRON
				FLOAT IronMass;
			#endif
	
			#ifdef COBALT
				FLOAT CobaltMass;
			#endif

			#ifdef NICKEL
				FLOAT NickelMass;
			#endif

			#ifdef COPPER
				FLOAT CopperMass;
			#endif

			#ifdef ZINC
				FLOAT ZincMass;
			#endif

      float temp;
      double m_vir;
      int in_halo;

}
 *StarP,					//Used in star_form.c to help move particles in P[]
 *P,              /*!< holds particle data on local processor */
 *DomainPartBuf;  /*!< buffer for particle data used in domain decomposition */


/* the following struture holds data that is stored for each SPH particle in addition 
   to the collisionless variables.
 */
extern struct sph_particle_data
{
  FLOAT Entropy;       /*!< current value of entropy (actually entropic function) 
                                of particle. Per p1108 in GADGET2 this is PHYSICAL
																entropy. */
  FLOAT Density;		       /*!< current COMOVING baryonic mass density of particle */
  FLOAT Hsml;			         /*!< current smoothing length */
  FLOAT Left;              /*!< lower bound in iterative smoothing length search */  
  FLOAT Right;             /*!< upper bound in iterative smoothing length search */ 
  FLOAT NumNgb;            /*!< weighted number of neighbours found */
  FLOAT Pressure;		       /*!<  = Entropy_phys * Density_comov^Gamma, is comoving*/
  FLOAT DtEntropy;         /*!< rate of change of PHYSICAL entropy (in internal code
																code units*/
  FLOAT HydroAccel[3]; 		 /*!< acceleration due to hydrodynamical force */
  FLOAT VelPred[3];		     /*!< predicted SPH particle velocity at the current time */
  FLOAT DivVel;			       /*!< local velocity divergence */
  FLOAT CurlVel;		       /*!< local velocity curl. ALWAYS postive, see density.c*/
  FLOAT Rot[3];		         /*!< local velocity curl */
  FLOAT DhsmlDensityFactor;/*!< correction factor needed in the equation of motion of 
                                the conservative entropy formulation of SPH. Related
																to f_i and grad_i(W_ij) from eqns 7 and 8. I can't 
																quite get the factors to work out, but its damn close*/
  FLOAT MaxSignalVel;      /*!< maximum "signal velocity" occuring for this particle 
																Eqn. 13 in GADGET2. Also used in dt_courant*/
  	FLOAT cooling_rate;		//In internal Units.

			int cooling_on;			//Flag to turn cooling off and on.
			FLOAT time_turned_off;	//Age of universe when flag was turned off (My).


}
 *SphP,                  	 /*!< holds SPH particle data on local processor */
 *DomainSphBuf;            /*!< buffer for SPH particle data in domain decomposition. 
																Points to previously allocated CommBuffer */





/*  Variables for Tree
 */

extern int MaxNodes;		   /*!< maximum allowed number of internal nodes */
extern int Numnodestree;	 /*!< number of (internal) nodes in each tree */

extern struct NODE
{
  FLOAT len;			         /*!< sidelength of treenode */
  FLOAT center[3];		     /*!< geometrical center of node */
	#ifdef ADAPTIVE_GRAVSOFT_FORGAS
  	FLOAT maxsoft;           /*!< hold the maximum gravitational softening of particles
																	in the node if the ADAPTIVE_GRAVSOFT_FORGAS option is
      	                          selected */
	#endif
  union
  {
    int suns[8];		       /*!< temporary pointers to daughter nodes */
    struct
    {
      FLOAT s[3];          /*!< center of mass of node */
      FLOAT mass;          /*!< mass of node */
      int bitflags;        /*!< a bit-field with various information on the node */
      int sibling;         /*!< this gives the next node in the walk in case the current
                                node can be used */
      int nextnode;        /*!< this gives the next node in case the current node needs
                                to be opened */
      int father;          /*!< this gives the parent node of each node (or -1 if we 
                                have the root node) */
    }
    d;
  }
  u;
}
 *Nodes_base,              /*!< points to the actual memory allocted for the nodes */
 *Nodes;                   /*!< this is a pointer used to access the nodes which is 
                                shifted such that Nodes[All.MaxPart] gives the first 
                                allocated node */


extern int *Nextnode;	     /*!< gives next node in tree walk */
extern int *Father;	       /*!< gives parent node in tree    */


extern struct extNODE      /*!< this structure holds additional tree-node information 
                                which is not needed in the actual gravity computation*/
{
  FLOAT hmax;			         /*!< maximum SPH smoothing length in node. Only used for gas
                                particles */
  FLOAT vs[3];			       /*!< center-of-mass velocity */
}
 *Extnodes_base,           /*!< points to the actual memory allocted for the extended 
                                node information */
 *Extnodes;                /*!< provides shifted access to extended node information, 
                                parallel to Nodes/Nodes_base */





/*! Header for the standard file format.

  !!!!!! IMPORTANT !!!!!!

  The ORDER of the header MATTERS! Because the header has a double which is 
  an 8 byte word, memory gets allocated in 8 byte chunks. So... Volker has
  ordered the 'ints' and the 'doubles' such that the memory for each variable is
  allocated contiguously. This means that you need 'ints' allocated in pairs
  of 2 so that they use up the entire 8 byte word.
 */
extern struct io_header
{
  int npart[6];            /*!< number of particles of each type in this file */
  double mass[6];          /*!< mass of particles of each type. If 0, then the masses 
                                are explicitly stored in the mass-block of the snapshot
                                file, otherwise they are omitted */
  double time;             /*!< time of snapshot file */
  double redshift;         /*!< redshift of snapshot file */
  int flag_sfr;            /*!< flags whether the simulation was including star 
                                formation */
  int flag_feedback;       /*!< flags whether feedback was included (obsolete) */
  unsigned int npartTotal[6];  /*!< total number of particles of each type in this 
                                    snapshot. This can be different from npart if one is
                                     dealing with a multi-file snapshot. */
  int flag_cooling;            /*!< flags whether cooling was included  */
  int num_files;               /*!< number of files in multi-file snapshot */
  double BoxSize;              /*!< box-size of simulation in case periodic boundaries 
                                    were used */
  double Omega0;               /*!< matter density in units of critical density */
  double OmegaLambda;          /*!< cosmological constant parameter */
  double HubbleParam;          /*!< Hubble parameter in units of 100 km/sec/Mpc */
  int flag_stellarage;         /*!< flags whether the file contains formation times of 
                                    star particles */
  int flag_metals;             /*!< flags whether the file contains metallicity values 
                                    for gas and star particles */
  unsigned int npartTotalHighWord[6];  /*!< High word of the total number of particles 
                                            of each type */
  int  flag_entropy_instead_u;         /*!< flags that IC-file contains entropy instead
                                            of u */
  //ADD NEW FLAGS BELOW HERE

  int flag_RII;               //TypeII SN rate per year
  int flag_RIa;               //TypeIa SN rate per year

	char flag_Carbon;
	char flag_Nitrogen;
	char flag_Oxygen;
	char flag_Florine;
	char flag_Neon;
	char flag_Sodium;
	char flag_Magnesium;
	char flag_Aluminum;
	char flag_Silicon;
	char flag_Phosphorus;
	char flag_Sulfur;
	char flag_Chlorine;
	char flag_Argon;
	char flag_Potassium;
	char flag_Calcium;
	char flag_Scandium;
	char flag_Titanium;
	char flag_Vanadium;
	char flag_Chromium;
	char flag_Manganese;
	char flag_Iron;
	char flag_Cobalt;
	char flag_Nickel;
	char flag_Copper;
	char flag_Zinc;

  char fill[27];	                     /*!< fills to 256 Bytes */
}
 header;                               /*!< holds header for snapshot files */


#define IO_NBLOCKS 39      /*!< total number of defined information blocks for snapshot
                                files. Must be equal to the number of entries in "enum 
                                iofields" */

enum iofields              /*!< this enumeration lists the defined output blocks in 
                                snapshot files. Not all of them need to be present. */
{ 
  IO_POS,				// =  0				= Position
  IO_VEL,				// =  1       = Velocity
  IO_ID,				// =  2       = Identifier
  IO_MASS,    	// =  3 			= Mass
  IO_U,					// =  4				= Internal Energy
  IO_RHO, 			// =  5				= Internal Density
  IO_HSML,			// =  6				= Smoothing Length
  IO_POT,				// =  7				= Potential Energy
  IO_ACCEL,			// =  8				= Acceleration
  IO_DTENTR,		// =  9				= Change in Entropy
  IO_TSTP,			// = 10				= Time Step
	//IO_SFR,				// = 11				= Star formation rate
  IO_RII,
  IO_RIA,
	IO_STELLARAGE,// = 12				= Stellar Age
	IO_CARBON,		// = 13				= Carbon Mass Frac
	IO_NITROGEN,	// = 14				= Nitrogen Mass Frac
	IO_OXYGEN,		// = 15				= Oxygen Mass Frac
	IO_FLORINE, 	// = 16       = Florine Mass Frac
	IO_NEON,			// = 17				= Neon Mass Frac
	IO_SODIUM,		// = 18				= Sodium Mass Frac
	IO_MAGNESIUM, // = 19			  = Magnesium Mass Frac
	IO_ALUMINUM,  // = 20				= Aluminum Mass Frac
	IO_SILICON,   // = 21				= Silicon Mass Frac
	IO_PHOSPHORUS,// = 22				= Phosphorus Mass Frac
	IO_SULFUR,		// = 23				= Sulfur Mass Frac
	IO_CHLORINE,  // = 24				= Chlorine Mass Frac
	IO_ARGON,		  // = 25				= Argon Mass Frac
	IO_POTASSIUM, // = 26				= Potassium Mass Frac
	IO_CALCIUM,	  // = 27				= Calcium Mass Frac
	IO_SCANDIUM,  // = 28				= Scandium Mass Frac
	IO_TITANIUM,  // = 29				= Titanium Mass Frac
	IO_VANADIUM,  // = 30				= Vanadium Mass Frac
	IO_CHROMIUM,  // = 31				= Chromium Mass Frac
	IO_MANGANESE, // = 32 			= Manganese Mass Frac
	IO_IRON,			// = 33				= Iron Mass Frac
	IO_COBALT,		// = 34				= Cobalt Mass Frac
	IO_NICKEL,		// = 35				= Nickel Mass Frac
	IO_COPPER,		// = 36				= Copper Mass Frac
	IO_ZINC,			// = 37				= Zinc Mass Frac
};


extern char Tab_IO_Labels[IO_NBLOCKS][4];   /*<! This table holds four-byte character 
                                                 tags used for fileformat 2 */


/* global state of system, used for global statistics
 */
extern struct state_of_system
{
  double Mass;
  double EnergyKin;
  double EnergyPot;
  double EnergyInt;
  double EnergyTot;
  double Momentum[4];
  double AngMomentum[4];
  double CenterOfMass[4];
  double MassComp[6];
  double EnergyKinComp[6];
  double EnergyPotComp[6];
  double EnergyIntComp[6];
  double EnergyTotComp[6];
  double MomentumComp[6][4]; 
  double AngMomentumComp[6][4]; 
  double CenterOfMassComp[6][4];
}
 SysState;                       /*<! Structure for storing some global statistics 
                                      about the simulation. */
 


/* Various structures for communication
 */
extern struct gravdata_in
{
  union
  {
    FLOAT Pos[3];
    FLOAT Acc[3];
    FLOAT Potential;
  }
  u;
	#ifdef UNEQUALSOFTENINGS
  	int Type;
		#ifdef ADAPTIVE_GRAVSOFT_FORGAS
		  FLOAT Soft;
		#endif
	#endif
  union
  {
    FLOAT OldAcc;
    int Ninteractions;
  }
  w;
}
 *GravDataIn,         /*!< holds particle data to be exported to other processors */
 *GravDataGet,        /*!< holds particle data imported from other processors */
 *GravDataResult,     /*!< holds the partial results computed for imported particles. 
                           Note: We use GravDataResult = GravDataGet, such that the 
                           result replaces the imported data */
 *GravDataOut;        /*!< holds partial results received from other processors. This 
                           will overwrite the GravDataIn array */

extern struct gravdata_index
{
  int Task;
  int Index;
  int SortIndex;
}
 *GravDataIndexTable;  /*!< the particles to be exported are grouped by task-number. 
                            This table allows the results to be disentangled again and 
                            to be assigned to the correct particle */



extern struct densdata_in
{
  FLOAT Pos[3];
  FLOAT Vel[3];
  FLOAT Hsml;
  int Index;
  int Task;
}
 *DensDataIn,          /*!< holds particle data for SPH density computation to be 
                            exported to other processors */
 *DensDataGet;         /*!< holds imported particle data for SPH density computation */

extern struct densdata_out
{
  FLOAT Rho;
  FLOAT Div, Rot[3];
  FLOAT DhsmlDensity;
  FLOAT Ngb;
}
 *DensDataResult,          /*!< stores the locally computed SPH density results for 
                                imported particles */
 *DensDataPartialResult;   /*!< imported partial SPH density results from other 
                                processors */



extern struct hydrodata_in
{
  FLOAT Pos[3];
  FLOAT Vel[3];
  FLOAT Hsml;
  FLOAT Mass;
  FLOAT Density;
  FLOAT Pressure;
  FLOAT F1;											//See GADGET2 eqn 17.
  FLOAT DhsmlDensityFactor;
  int   Timestep;
  int   Task;										//Task to export LOCAL i'th to
  int   Index;									//LOCAL particle index  (i)
}
 *HydroDataIn,             /*!< holds particle data for SPH hydro-force computation to 
                                be exported to other processors. Points to CommBuffer*/
 *HydroDataGet;            /*!< holds imported particle data for SPH hydro-force 
                                computation */

extern struct hydrodata_out
{
  FLOAT Acc[3];
  FLOAT DtEntropy;
  FLOAT MaxSignalVel;
	#ifdef COOLING
  	FLOAT cooling_rate;
	#endif
}
 *HydroDataResult,         /*!< stores the locally computed SPH hydro results for 
                                imported particles */
 *HydroDataPartialResult;  /*!< imported partial SPH hydro-force results from other 
                                processors */


/*
	This is used to export the location and hsml of star particles to get the 
	normalization constant from eqn 6 in Stinson et al. (2006)
*/
extern struct feedbacknormdata_in
{
	FLOAT Pos[3];
	FLOAT Hsml;
	int Task;
	int Index;

	#ifndef LONGIDS
		unsigned int ID;					
	#else
		unsigned long long ID;
	#endif
}
*FeedbackNormDataIn,
*FeedbackNormDataGet;


extern struct feedbacknormdata_out
{
	FLOAT Pos[3];
	FLOAT FeedbackNorm;				//Normalization Constant
	FLOAT GasDensity;					//Use SPH kernal to find
	FLOAT GasTemperature;				//Use SPH kernal to find
	int Task;

	#ifndef LONGIDS
		unsigned int ID;					
	#else
		unsigned long long ID;
	#endif
}
*FeedbackNormDataResult,
*FeedbackNormDataPartialResult;


/*
	This is used to export the location of star particles to search for the number of 
	nearby SphP[] on other processors.
*/
extern struct gaspartsearchdata_in
{
	FLOAT Pos[3];
	int Task;											/* Destination Task, not LOCAL Task */
	int Index;										/* Index in LOCAL Task */
	FLOAT Hsml;			      				//Star Smoothing Length

	#ifndef LONGIDS
		unsigned int ID;					
	#else
		unsigned long long ID;
	#endif
}
*GasPartSearchDataIn,
*GasPartSearchDataGet;


/*
	This is used to recieve the number of nearby SphP[] on other processors.
*/
extern struct gaspartsearchdata_out
{
	FLOAT Pos[3];
	int Task;
	#ifndef LONGIDS
		unsigned int ID;
	#else
		unsigned long long ID;
	#endif
	int NumGasNearby;
}
*GasPartSearchDataResult,
*GasPartSearchDataPartialResult;



extern struct feedbackdata_in
{
	FLOAT Pos[3];
	double DtEnergy;
	FLOAT DtMass;
	FLOAT R_II;
	FLOAT BlastRadius;
	FLOAT Dt_Z_Mass[NELEMENTS];	/* Should be a pointer, but it would be a pain to 
																 know the length point we don't know a priori how 
																 long NElements is 
													 		*/
	int Task;
	int Index;
	long NumGasNearby;					// Gas w/in FeedbackRadius
	#ifndef LONGIDS
		unsigned int ID;
	#else
		unsigned long long ID;
	#endif

	FLOAT Hsml;
	FLOAT FeedbackNorm;
}
*FeedbackDataIn,
*FeedbackDataGet;




/* This is unneccessary b/c there are no interactions from the gas to star particle.
	 So, the stars only need to give and not receive. HOWEVER, we'll use this to double
	 check that the correct number of gas particles were fed (i.e. NumGasNearby = */
extern struct feedbackdata_out
{
	int Task;
	long NumGasFed;
	#ifndef LONGIDS
		unsigned int ID;
	#else
		unsigned long long ID;
	#endif
}
*FeedbackDataResult,
*FeedbackDataPartialResult;


extern double phys_time;


/* Node Structures */
extern int ngrid;
extern char snap_filename[MAXLEN_FILENAME];


/* Grid Structures */
typedef struct CELL
{
  float Rho;
  float min[3];       //min bounds of volume element
  float max[3];       //max bounds of volume element
} CELL;

typedef struct GRID
{
  int ID;
  int isActive;
  float Rho;
  float Temperature;
  float Pos[3];        //Location of grid.
  int   Task;         //Task that computes this density.
  float Hsml;
  float HsmlLeft;
  float HsmlRight;
  int NumGasNearby;

} GRID;

extern CELL *** cell;
extern GRID *** grid;
extern GRID **  grid_1D;    //Points to elements in grid, makes easier to loop over
extern int      ngrid;      //number of grid elements on one side.
extern int      NumGrid;    //Number of grid elements... = (ngrid^3)
extern int      NumGridUpdate;

/* Parallel Structures */
extern int    ndim;         // NTask = (ndim)^3, ndim = 2, 3, 4



   // Simplified particle data struct
   typedef struct SIMPLE_PARTICLE_DATA
   {
      float Pos[3];
      float Vel[3];
      float Mass;
      float Density;
      float Hsml;
      int Type;
      int ID;
   } SIMPLE_PARTICLE_DATA;

   extern SIMPLE_PARTICLE_DATA *SP;
#endif
