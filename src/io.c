/****************************************************************************
Author: Volker Springel
Modified: Jared Coughlin and Ali Snedden
Date: 5/24/13

Purpose: This writes the snapshot files.

	NOTE: SPLASH expects the cooling parameters (Ne and Nh) to be written 
				BEFORE the smoothing length and the SFR to be written AFTER the 
				smoothing length. As for Stellar Age, that will require telling 
				SPLASH to expect an extra column. 

	NOTE: THis will only be tested for file type 1 and 1 file written in parallel.


****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <errno.h>

#ifdef HAVE_HDF5
#include <hdf5.h>
#endif

#include "allvars.h"
#include "proto.h"



/*! \file io.c
 *  \brief Routines for producing a snapshot file on disk.
 */

//static int n_type[6];												//Number of LOCAL particles per type
//static long long ntot_type_all[6];					//Number of TOTAL particles per type





/*! This function tells the size of one data entry in each of the blocks
 *  defined for the output file. If one wants to add a new output-block, this
 *  function should be augmented accordingly.
 */
int get_bytes_per_blockelement(enum iofields blocknr)
{
  int bytes_per_blockelement = 0;

  switch (blocknr)
  {
    case IO_POS:
    case IO_VEL:
    case IO_ACCEL:
      bytes_per_blockelement = 3 * sizeof(float);
      break;

    case IO_ID:
			#ifdef LONGIDS
  	    bytes_per_blockelement = sizeof(long long);
			#else
      	bytes_per_blockelement = sizeof(int);
			#endif
      break;

    case IO_MASS:
    case IO_U:
    case IO_RHO:
    case IO_HSML:
    case IO_POT:
    case IO_DTENTR:
    case IO_TSTP:
		//case IO_SFR:
		case IO_RII:
		case IO_RIA:
		case IO_STELLARAGE:
		case IO_CARBON:
		case IO_NITROGEN:
		case IO_OXYGEN:
		case IO_FLORINE:
		case IO_NEON:
		case IO_SODIUM:
		case IO_MAGNESIUM:
		case IO_ALUMINUM:
		case IO_SILICON:
		case IO_PHOSPHORUS:
		case IO_SULFUR:
		case IO_CHLORINE:
		case IO_ARGON:
		case IO_POTASSIUM:
		case IO_CALCIUM:
		case IO_SCANDIUM:
		case IO_TITANIUM:
		case IO_VANADIUM:
		case IO_CHROMIUM:
		case IO_MANGANESE:
		case IO_IRON:
		case IO_COBALT:
		case IO_NICKEL:
		case IO_COPPER:
		case IO_ZINC:
      bytes_per_blockelement = sizeof(float);
      break;
  }

  return bytes_per_blockelement;
}


/*! This function returns the type of the data contained in a given block of
 *  the output file. If one wants to add a new output-block, this function
 *  should be augmented accordingly.
 */
int get_datatype_in_block(enum iofields blocknr)
{
  int typekey;

  switch (blocknr)
  {
    case IO_ID:
			#ifdef LONGIDS
      	typekey = 2;		/* native long long */
			#else
      	typekey = 0;		/* native int */
			#endif
      break;

    default:
      typekey = 1;		/* native float */
      break;
  }

  return typekey;
}


/*! This function informs about the number of elements stored per particle for
 *  the given block of the output file. If one wants to add a new
 *  output-block, this function should be augmented accordingly.
 */
int get_values_per_blockelement(enum iofields blocknr)
{
  int values = 0;

  switch (blocknr)
  {
    case IO_POS:
    case IO_VEL:
    case IO_ACCEL:
      values = 3;
      break;

    case IO_ID:
    case IO_MASS:
    case IO_U:
    case IO_RHO:
    case IO_HSML:
    case IO_POT:
    case IO_DTENTR:
    case IO_TSTP:
		//case IO_SFR:
		case IO_RII:
		case IO_RIA:
		case IO_STELLARAGE:
		case IO_CARBON:
		case IO_NITROGEN:
		case IO_OXYGEN:
		case IO_FLORINE:
		case IO_NEON:
		case IO_SODIUM:
		case IO_MAGNESIUM:
		case IO_ALUMINUM:
		case IO_SILICON:
		case IO_PHOSPHORUS:
		case IO_SULFUR:
		case IO_CHLORINE:
		case IO_ARGON:
		case IO_POTASSIUM:
		case IO_CALCIUM:
		case IO_SCANDIUM:
		case IO_TITANIUM:
		case IO_VANADIUM:
		case IO_CHROMIUM:
		case IO_MANGANESE:
		case IO_IRON:
		case IO_COBALT:
		case IO_NICKEL:
		case IO_COPPER:
		case IO_ZINC:
      values = 1;
      break;
  }

  return values;
}


/*! This function determines how many particles there are in a given block,
 *  based on the information in the header-structure.  It also flags particle
 *  types that are present in the block in the typelist array. If one wants to
 *  add a new output-block, this function should be augmented accordingly.
 *
 *	Note that this is for the GLOBAL (i.e. NOT LOCAL) number of particles.	
 */
int get_particles_in_block(enum iofields blocknr, int *typelist)
{
  int i;							//Indice
	int nall;    				//Total number of particles
	int ntot_withmasses, ngas, nstars;

  nall = 0;
  ntot_withmasses = 0;

  for(i = 0; i < 6; i++)
  {
    typelist[i] = 0;

		//sum up total number of particles and flag particle types
    if(header.npart[i] > 0)
	  {
	    nall += header.npart[i];
	    typelist[i] = 1;
	  }

		//Get number of particles with P[].Mass
    if(All.MassTable[i] == 0)
	    ntot_withmasses += header.npart[i];
  }

  ngas = header.npart[0];
  nstars = header.npart[4];


  switch (blocknr)
  {
    case IO_POS:
    case IO_VEL:
    case IO_ACCEL:
    case IO_TSTP:
    case IO_ID:
    case IO_POT:
      return nall;
      break;

		/*IO_MASS not used for types that have All.MassTable[i] > 0*/
    case IO_MASS:
      for(i = 0; i < 6; i++)
	    {
	      typelist[i] = 0;
	      if(All.MassTable[i] == 0 && header.npart[i] > 0)
	        typelist[i] = 1;
	    }
      return ntot_withmasses;
      break;

    case IO_U:
    case IO_RHO:
    case IO_HSML:
    case IO_DTENTR:
		//case IO_SFR:
      for(i = 1; i < 6; i++)
	      typelist[i] = 0;
      return ngas;
      break;

		//ALI : Fix?
		case IO_STELLARAGE:
    case IO_RII:
    case IO_RIA:
			for(i=0; i<6; i++)
			{
				typelist[i] = 0;
				if(i == 4){
					typelist[i] = 1;
				}
			}
			return nstars;
			break;

		/*Elements.  These hold for gas and stars, so need to return ngas+nstars.*/
		case IO_CARBON:
		case IO_NITROGEN:
		case IO_OXYGEN:
		case IO_FLORINE:
		case IO_NEON:
		case IO_SODIUM:
		case IO_MAGNESIUM:
		case IO_ALUMINUM:
		case IO_SILICON:
		case IO_PHOSPHORUS:
		case IO_SULFUR:
		case IO_CHLORINE:
		case IO_ARGON:
		case IO_POTASSIUM:
		case IO_CALCIUM:
		case IO_SCANDIUM:
		case IO_TITANIUM:
		case IO_VANADIUM:
		case IO_CHROMIUM:
		case IO_MANGANESE:
		case IO_IRON:
		case IO_COBALT:
		case IO_NICKEL:
		case IO_COPPER:
		case IO_ZINC:
			for(i = 0; i < 6; i++)
			{
				typelist[i] = 0;
				if(i == 0 || i == 4)
				{
					typelist[i] = 1;
				}
			}
			return (ngas + nstars);
			break;

  }

  endrun(212);
  return 0;
}



/*! This function tells whether or not a given block in the output file is
 *  present, depending on the type of simulation run and the compile-time
 *  options. If one wants to add a new output-block, this function should be
 *  augmented accordingly.
 */
int blockpresent(enum iofields blocknr)
{
	#ifndef OUTPUTPOTENTIAL
	  if(blocknr == IO_POT)
  	  return 0;
	#endif

	#ifndef OUTPUTACCELERATION
  	if(blocknr == IO_ACCEL)
    	return 0;
	#endif

	#ifndef OUTPUTCHANGEOFENTROPY
  	if(blocknr == IO_DTENTR)
    	return 0;
	#endif

	#ifndef OUTPUTTIMESTEP
  	if(blocknr == IO_TSTP)
    	return 0;
	#endif

	/*#ifndef SFR
		if(blocknr == IO_SFR)
			return 0;
	#endif*/

	#ifndef FEEDBACK
		if(blocknr == IO_RII)
			return 0;
	#endif

	#ifndef FEEDBACK
		if(blocknr == IO_RIA)
			return 0;
	#endif

	#ifndef STELLARAGE
		if(blocknr == IO_STELLARAGE)
			return 0;
	#endif

	#ifndef CARBON
		if(blocknr == IO_CARBON)
			return 0;
	#endif

	
	#ifndef NITROGEN
		if(blocknr == IO_NITROGEN)
			return 0;
	#endif

	
	#ifndef OXYGEN
		if(blocknr == IO_OXYGEN)
			return 0;
	#endif

	
	#ifndef FLORINE
		if(blocknr == IO_FLORINE)
			return 0;
	#endif

	
	#ifndef NEON
		if(blocknr == IO_NEON)
			return 0;
	#endif

	
	#ifndef SODIUM
		if(blocknr == IO_SODIUM)
			return 0;
	#endif

	
	#ifndef MAGNESIUM
		if(blocknr == IO_MAGNESIUM)
			return 0;
	#endif

	
	#ifndef ALUMINUM
		if(blocknr == IO_ALUMINUM)
			return 0;
	#endif

	
	#ifndef SILICON
		if(blocknr == IO_SILICON)
			return 0;
	#endif

	
	#ifndef PHOSPHORUS
		if(blocknr == IO_PHOSPHORUS)
			return 0;
	#endif

	
	#ifndef SULFUR
		if(blocknr == IO_SULFUR)
			return 0;
	#endif

	
	#ifndef CHLORINE
		if(blocknr == IO_CHLORINE)
			return 0;
	#endif

	
	#ifndef ARGON
		if(blocknr == IO_ARGON)
			return 0;
	#endif

	
	#ifndef POTASSIUM
		if(blocknr == IO_POTASSIUM)
			return 0;
	#endif

	
	#ifndef CALCIUM
		if(blocknr == IO_CALCIUM)
			return 0;
	#endif

	
	#ifndef SCANDIUM
		if(blocknr == IO_SCANDIUM)
			return 0;
	#endif

	
	#ifndef TITANIUM
		if(blocknr == IO_TITANIUM)
			return 0;
	#endif

	
	#ifndef VANADIUM
		if(blocknr == IO_VANADIUM)
			return 0;
	#endif

	
	#ifndef CHROMIUM
		if(blocknr == IO_CHROMIUM)
			return 0;
	#endif

	
	#ifndef MANGANESE
		if(blocknr == IO_MANGANESE)
			return 0;
	#endif

	
	#ifndef IRON
		if(blocknr == IO_IRON)
			return 0;
	#endif

	
	#ifndef COBALT
		if(blocknr == IO_COBALT)
			return 0;
	#endif

	
	#ifndef NICKEL
		if(blocknr == IO_NICKEL)
			return 0;
	#endif

	
	#ifndef COPPER
		if(blocknr == IO_COPPER)
			return 0;
	#endif

	
	#ifndef ZINC
		if(blocknr == IO_ZINC)
			return 0;
	#endif

  return 1;			/* default: present */
}




/*! This function associates a short 4-character block name with each block
 *  number.  This is stored in front of each block for snapshot
 *  FileFormat=2. If one wants to add a new output-block, this function should
 *  be augmented accordingly.
 */
void fill_Tab_IO_Labels(void)
{
  enum iofields i;

  for(i = 0; i < IO_NBLOCKS; i++)
    switch (i)
    {
      case IO_POS:
	      strncpy(Tab_IO_Labels[IO_POS], "POS ", 4);
	      break;
      case IO_VEL:
	      strncpy(Tab_IO_Labels[IO_VEL], "VEL ", 4);
	      break;
      case IO_ID:
	      strncpy(Tab_IO_Labels[IO_ID], "ID  ", 4);
	      break;
      case IO_MASS:
	      strncpy(Tab_IO_Labels[IO_MASS], "MASS", 4);
	      break;
      case IO_U:
	      strncpy(Tab_IO_Labels[IO_U], "U   ", 4);
	      break;
      case IO_RHO:
	      strncpy(Tab_IO_Labels[IO_RHO], "RHO ", 4);
	      break;
      case IO_HSML:
	      strncpy(Tab_IO_Labels[IO_HSML], "HSML", 4);
	      break;
      case IO_POT:
	      strncpy(Tab_IO_Labels[IO_POT], "POT ", 4);
	      break;
      case IO_ACCEL:
	      strncpy(Tab_IO_Labels[IO_ACCEL], "ACCE", 4);
	      break;
      case IO_DTENTR:
	      strncpy(Tab_IO_Labels[IO_DTENTR], "ENDT", 4);
	      break;
      case IO_TSTP:
	      strncpy(Tab_IO_Labels[IO_TSTP], "TSTP", 4);
	      break;
			/*case IO_SFR:
				strncpy(Tab_IO_Labels[IO_SFR], "SFR ", 4);
				break;*/	
			case IO_RII:
				strncpy(Tab_IO_Labels[IO_RII], "RII ", 4);
				break;
			case IO_RIA:
				strncpy(Tab_IO_Labels[IO_RIA], "RIA ", 4);
				break;
			case IO_STELLARAGE:
				strncpy(Tab_IO_Labels[IO_STELLARAGE], "SAGE", 4);
				break;	
			case IO_CARBON:
				strncpy(Tab_IO_Labels[IO_CARBON], "CMF ", 4);
				break;
			case IO_NITROGEN:
				strncpy(Tab_IO_Labels[IO_NITROGEN], "NMF ", 4);
				break;
			case IO_OXYGEN:
				strncpy(Tab_IO_Labels[IO_OXYGEN], "OMF ", 4);
				break;
			case IO_FLORINE:
				strncpy(Tab_IO_Labels[IO_FLORINE], "FMF ", 4);
				break;
			case IO_NEON:
				strncpy(Tab_IO_Labels[IO_NEON], "NeMF", 4);
				break;
			case IO_SODIUM:
				strncpy(Tab_IO_Labels[IO_SODIUM], "NaMF", 4);
				break;
			case IO_MAGNESIUM:
				strncpy(Tab_IO_Labels[IO_MAGNESIUM], "MgMF", 4);
				break;
			case IO_ALUMINUM:
				strncpy(Tab_IO_Labels[IO_ALUMINUM], "AlMF", 4);
				break;
			case IO_SILICON:
				strncpy(Tab_IO_Labels[IO_SILICON], "SiMF", 4);
				break;
			case IO_PHOSPHORUS:
				strncpy(Tab_IO_Labels[IO_PHOSPHORUS], "PMF ", 4);
				break;
			case IO_SULFUR:
				strncpy(Tab_IO_Labels[IO_SULFUR], "SMF ", 4);
				break;
			case IO_CHLORINE:
				strncpy(Tab_IO_Labels[IO_CHLORINE], "ClMF", 4);
				break;
			case IO_ARGON:
				strncpy(Tab_IO_Labels[IO_ARGON], "ArMF", 4);
				break;
			case IO_POTASSIUM:
				strncpy(Tab_IO_Labels[IO_POTASSIUM], "KMF ", 4);
				break;
			case IO_CALCIUM:
				strncpy(Tab_IO_Labels[IO_CALCIUM], "CaMF", 4);
				break;
			case IO_SCANDIUM:
				strncpy(Tab_IO_Labels[IO_SCANDIUM], "ScMF", 4);
				break;
			case IO_TITANIUM:
				strncpy(Tab_IO_Labels[IO_TITANIUM], "TiMF", 4);
				break;
			case IO_VANADIUM:
				strncpy(Tab_IO_Labels[IO_VANADIUM], "VMF ", 4);
				break;
			case IO_CHROMIUM:
				strncpy(Tab_IO_Labels[IO_CHROMIUM], "CrMF", 4);
				break;
			case IO_MANGANESE:
				strncpy(Tab_IO_Labels[IO_MANGANESE], "MnMF", 4);
				break;
			case IO_IRON:
				strncpy(Tab_IO_Labels[IO_IRON], "FeMF", 4);
				break;
			case IO_COBALT:
				strncpy(Tab_IO_Labels[IO_COBALT], "CoMF", 4);
				break;
			case IO_NICKEL:
				strncpy(Tab_IO_Labels[IO_NICKEL], "NiMF", 4);
				break;
			case IO_COPPER:
				strncpy(Tab_IO_Labels[IO_COPPER], "CuMF", 4);
				break;
			case IO_ZINC:
				strncpy(Tab_IO_Labels[IO_ZINC], "ZnMF", 4);
				break;
    }
}

/*! This function returns a descriptive character string that describes the
 *  name of the block when the HDF5 file format is used.  If one wants to add
 *  a new output-block, this function should be augmented accordingly.
 */
void get_dataset_name(enum iofields blocknr, char *buf)
{

  strcpy(buf, "default");

  switch (blocknr)
  {
    case IO_POS:
      strcpy(buf, "Coordinates");
      break;
    case IO_VEL:
      strcpy(buf, "Velocities");
      break;
    case IO_ID:
      strcpy(buf, "ParticleIDs");
      break;
    case IO_MASS:
      strcpy(buf, "Masses");
      break;
    case IO_U:
      strcpy(buf, "InternalEnergy");
      break;
    case IO_RHO:
      strcpy(buf, "Density");
      break;
    case IO_HSML:
      strcpy(buf, "SmoothingLength");
      break;
    case IO_POT:
      strcpy(buf, "Potential");
      break;
    case IO_ACCEL:
      strcpy(buf, "Acceleration");
      break;
    case IO_DTENTR:
      strcpy(buf, "RateOfChangeOfEntropy");
      break;
    case IO_TSTP:
      strcpy(buf, "TimeStep");
      break;
		/*case IO_SFR:
			strcpy(buf, "StarFormRate");
			break;*/
		case IO_RII:
			strcpy(buf, "RII");
			break;
		case IO_RIA:
			strcpy(buf, "RIa");
			break;
		case IO_STELLARAGE:
			strcpy(buf, "StarAge");
			break;
		case IO_CARBON:
			strcpy(buf, "CarbonMass");
			break;
		case IO_NITROGEN:
			strcpy(buf, "NitrogenMass");
			break;
		case IO_OXYGEN:
			strcpy(buf, "OxygenMass");
			break;
		case IO_FLORINE:
			strcpy(buf, "FlorineMass");
			break;
		case IO_NEON:
			strcpy(buf, "NeonMass");
			break;
		case IO_SODIUM:
			strcpy(buf, "SodiumMass");
			break;
		case IO_MAGNESIUM:
			strcpy(buf, "MagnesiumMass");
			break;
		case IO_ALUMINUM:
			strcpy(buf, "AluminumMass");
			break;
		case IO_SILICON:
			strcpy(buf, "SiliconMass");
			break;
		case IO_PHOSPHORUS:
			strcpy(buf, "PhosphorusMass");
			break;
		case IO_SULFUR:
			strcpy(buf, "SulfurMass");
			break;
		case IO_CHLORINE:
			strcpy(buf, "ChlorineMass");
			break;
		case IO_ARGON:
			strcpy(buf, "ArgonMass");
			break;
		case IO_POTASSIUM:
			strcpy(buf, "PotassiumMass");
			break;
		case IO_CALCIUM:
			strcpy(buf, "CalciumMass");
		case IO_SCANDIUM:
			strcpy(buf, "ScandiumMass");
			break;
		case IO_TITANIUM:
			strcpy(buf, "TitaniumMass");
			break;
		case IO_VANADIUM:
			strcpy(buf, "VanadiumMass");
			break;
		case IO_CHROMIUM:
			strcpy(buf, "ChromiumMass");
			break;
		case IO_MANGANESE:
			strcpy(buf, "ManganeseMass");
			break;
		case IO_IRON:
			strcpy(buf, "IronMass");
			break;
		case IO_COBALT:
			strcpy(buf, "CobaltMass");
			break;
		case IO_NICKEL:
			strcpy(buf, "NickelMass");
			break;
		case IO_COPPER:
			strcpy(buf, "CopperMass");
			break;
		case IO_ZINC:
			strcpy(buf, "ZincMass");
			break;
  }
}





/*! This catches I/O errors occuring for my_fwrite(). In this case we
 *  better stop.
 */
size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE * stream)
{
  size_t nwritten;

  if((nwritten = fwrite(ptr, size, nmemb, stream)) != nmemb)
  {
    printf("I/O error (fwrite) on task=%d has occured: %s\n", ThisTask, strerror(errno));
    fflush(stdout);
    endrun(777);
  }
  return nwritten;
}


/*! This catches I/O errors occuring for fread(). In this case we
 *  better stop.
 */
size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE * stream)
{
  size_t nread;

  if((nread = fread(ptr, size, nmemb, stream)) != nmemb)
  {
    printf("I/O error (fread) on task=%d has occured: %s\n", ThisTask, strerror(errno));
    fflush(stdout);
    endrun(778);
  }
  return nread;
}
