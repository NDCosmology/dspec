/*************************************************************************************
Author: Ali Snedden 
Date: 	10/8/13

Dep. of Physics
University of Notre Dame

Purpose: 
	This is adapted from density.c. Below we want to implement mass, metal
	and energy feedback from supernovae.  This requires 3 searches tree walks.

	1. Iteratively search the tree and assign smoothing length to 
		 the star particles.

	2. Search the tree for SPH particles within the smoothing length.
		 Communicate neighbor particle info back to local particle / processor.
		 -> Calculate Normalization (eqn. 6 in Stinson et al. 2006)

	3. Search tree again. Distribute metals and energy.


NOTE: When you want to implement AGB stars, look at Wiersma et al. (2009).
  Its appendix provides a great place to start.



IMPORTANT DISCUSSION OF GSL LIBRARY:
  If the gsl_integration_qag() is used, at some point the following error is thrown
  "ERROR: roundoff error prevents tolerance from being achieved". I believe this
  is due to my flooring of the zMassFrac in TypeII_Zi_yield_mass_frac() when
  8 <= mass < 13. This is b/c I don't have data SN yield data for that range.

  So instead I use gsl_integration_glfixed(). It appears to give similarly decent
  results. I left in the lines related to gsl_integration_qag() to permit further
  testing in the future.
**************************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>


#include "allvars.h"
#include "proto.h"

//NOTE : If these values change, you'll need to change case IO_RII and IO_RIA in io.c


#ifdef PERIODIC
static double boxSize, boxHalf;

#ifdef LONG_X
static double boxSize_X, boxHalf_X;
#else
#define boxSize_X boxSize
#define boxHalf_X boxHalf
#endif
#ifdef LONG_Y
static double boxSize_Y, boxHalf_Y;
#else
#define boxSize_Y boxSize
#define boxHalf_Y boxHalf
#endif
#ifdef LONG_Z
static double boxSize_Z, boxHalf_Z;
#else
#define boxSize_Z boxSize
#define boxHalf_Z boxHalf
#endif
#endif



/***********************
        density
***********************/
void density(void)
{
  long long ntot;						//Total (global) number of SphP needing updated
	long long ntotleft;				//Number of SphP left to update
  int *noffset, *nbuffer, *nsend, *nsend_local;
	int *numlist; 						//List of SphP needing updated per processor. 
	int *ndonelist;
  int i, j, n;
	int ndone;								//Number of SphP done updating
	int npleft, maxfill, source, iter = 0;
  int level, ngrp, sendTask, recvTask, place, nexport;
  double dt_entr, tstart, tend, tstart_ngb = 0, tend_ngb = 0;
  double sumt, sumcomm, timengb, sumtimengb;
  double timecomp = 0, timeimbalance = 0, timecommsumm = 0, sumimbalance;
  MPI_Status status;

	#ifdef PERIODIC
  	boxSize = All.BoxSize;
  	boxHalf = 0.5 * All.BoxSize;
		#ifdef LONG_X
  		boxHalf_X = boxHalf * LONG_X;
  		boxSize_X = boxSize * LONG_X;
		#endif
		#ifdef LONG_Y
		  boxHalf_Y = boxHalf * LONG_Y;
		  boxSize_Y = boxSize * LONG_Y;
		#endif
		#ifdef LONG_Z
		  boxHalf_Z = boxHalf * LONG_Z;
		  boxSize_Z = boxSize * LONG_Z;
		#endif
	#endif

	//Allocate memory
  noffset = malloc(sizeof(int) * NTask);	/* offsets of bunches in common list */
  nbuffer = malloc(sizeof(int) * NTask);
  nsend_local = malloc(sizeof(int) * NTask);
  nsend = malloc(sizeof(int) * NTask * NTask);
  ndonelist = malloc(sizeof(int) * NTask);

	//Get number of SphP on local processor ready for update
  for(n = 0, NumSphUpdate = 0; n < N_gas; n++)
  {
    SphP[n].Left = SphP[n].Right = 0;

    if(P[n].Ti_endstep == All.Ti_Current)
	    NumSphUpdate++;
  }

	//Find total (global) number of SphP needing updated
  numlist = malloc(NTask * sizeof(int) * NTask);
  MPI_Allgather(&NumSphUpdate, 1, MPI_INT, numlist, 1, MPI_INT, MPI_COMM_WORLD);
  for(i = 0, ntot = 0; i < NTask; i++)
    ntot += numlist[i];
  free(numlist);



  /* we will repeat the whole thing for those particles where we didn't
     find enough neighbours */
  do
  {
    i = 0;							/* beginn with this index */
    ntotleft = ntot;		/* particles left for all tasks together */

    while(ntotleft > 0)
	  {
	    for(j = 0; j < NTask; j++)
	      nsend_local[j] = 0;

	    /* do local particles and prepare export list */
	    tstart = second();
	    for(nexport=0, ndone=0; i<N_gas && nexport< All.BunchSizeDensity - NTask; i++)
	      if(P[i].Ti_endstep == All.Ti_Current)
	      {
	      	ndone++;

					//Reset Exportflag. If true, particle needs exported 
	      	for(j = 0; j < NTask; j++)
	      	  Exportflag[j] = 0;

	      	density_evaluate(i, 0);			//Evaluate density

	      	for(j = 0; j < NTask; j++)
	    	  {
	    	    if(Exportflag[j])
		        {
	        		DensDataIn[nexport].Pos[0] = P[i].Pos[0];
	        		DensDataIn[nexport].Pos[1] = P[i].Pos[1];
	        		DensDataIn[nexport].Pos[2] = P[i].Pos[2];
	        		DensDataIn[nexport].Vel[0] = SphP[i].VelPred[0];
	        		DensDataIn[nexport].Vel[1] = SphP[i].VelPred[1];
	        		DensDataIn[nexport].Vel[2] = SphP[i].VelPred[2];
	        		DensDataIn[nexport].Hsml = SphP[i].Hsml;
	        		DensDataIn[nexport].Index = i;
	        		DensDataIn[nexport].Task = j;
	        		nexport++;
	        		nsend_local[j]++;
		        }
		      }
	      }
	    tend = second();
	    timecomp += timediff(tstart, tend);

	    qsort(DensDataIn, nexport, sizeof(struct densdata_in), dens_compare_key);

	    for(j = 1, noffset[0] = 0; j < NTask; j++)
	      noffset[j] = noffset[j - 1] + nsend_local[j - 1];

	    tstart = second();

	    MPI_Allgather(nsend_local, NTask, MPI_INT, nsend, NTask, MPI_INT, MPI_COMM_WORLD);

	    tend = second();
	    timeimbalance += timediff(tstart, tend);


	    /* now do the particles that need to be exported */

	    for(level = 1; level < (1 << PTask); level++)
	    {
	      tstart = second();
	      for(j = 0; j < NTask; j++)
		      nbuffer[j] = 0;
	      for(ngrp = level; ngrp < (1 << PTask); ngrp++)
		    {
		      maxfill = 0;
		      for(j = 0; j < NTask; j++)
		      {
		        if((j ^ ngrp) < NTask)
		        	if(maxfill < nbuffer[j] + nsend[(j ^ ngrp) * NTask + j])
		        	  maxfill = nbuffer[j] + nsend[(j ^ ngrp) * NTask + j];
		      }
		      if(maxfill >= All.BunchSizeDensity)
		        break;

		      sendTask = ThisTask;
		      recvTask = ThisTask ^ ngrp;

		      if(recvTask < NTask)
		      {
		        if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + 
               ThisTask] > 0)
		      	{
		      	  /* get the particles */
		      	  MPI_Sendrecv(&DensDataIn[noffset[recvTask]],
		         		       nsend_local[recvTask] * sizeof(struct densdata_in), MPI_BYTE,
		         		       recvTask, TAG_DENS_A,
		         		       &DensDataGet[nbuffer[ThisTask]],
		        		       nsend[recvTask * NTask + ThisTask] * sizeof(struct densdata_in),
		        		       MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, &status);
			      }
		      }

		      for(j = 0; j < NTask; j++)
		        if((j ^ ngrp) < NTask)
		          nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
		    }
	      tend = second();
	      timecommsumm += timediff(tstart, tend);


	      tstart = second();
	      for(j = 0; j < nbuffer[ThisTask]; j++)
	      	density_evaluate(j, 1);
	      tend = second();
	      timecomp += timediff(tstart, tend);

	      /* do a block to explicitly measure imbalance */
	      tstart = second();
	      MPI_Barrier(MPI_COMM_WORLD);
	      tend = second();
	      timeimbalance += timediff(tstart, tend);

	      /* get the result */
	      tstart = second();
	      for(j = 0; j < NTask; j++)
	      	nbuffer[j] = 0;
	      for(ngrp = level; ngrp < (1 << PTask); ngrp++)
	  	  {
	  	    maxfill = 0;
	  	    for(j = 0; j < NTask; j++)
		      {
		        if((j ^ ngrp) < NTask)
	        		if(maxfill < nbuffer[j] + nsend[(j ^ ngrp) * NTask + j])
	        		  maxfill = nbuffer[j] + nsend[(j ^ ngrp) * NTask + j];
		      }
	  	    if(maxfill >= All.BunchSizeDensity)
	  	      break;

	  	    sendTask = ThisTask;
		      recvTask = ThisTask ^ ngrp;

		      if(recvTask < NTask)
		      {
		        if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + 
               ThisTask] > 0)
			      {
			        /* send the results */
			        MPI_Sendrecv(&DensDataResult[nbuffer[ThisTask]],
				             nsend[recvTask * NTask + ThisTask] * sizeof(struct densdata_out),
				             MPI_BYTE, recvTask, TAG_DENS_B,
		      		       &DensDataPartialResult[noffset[recvTask]],
				             nsend_local[recvTask] * sizeof(struct densdata_out),
				             MPI_BYTE, recvTask, TAG_DENS_B, MPI_COMM_WORLD, &status);

			        /* add the result to the particles */
			        for(j = 0; j < nsend_local[recvTask]; j++)
			        {
			          source = j + noffset[recvTask];
			          place = DensDataIn[source].Index;

			          SphP[place].NumNgb += DensDataPartialResult[source].Ngb;
			          SphP[place].Density += DensDataPartialResult[source].Rho;
			          SphP[place].DivVel += DensDataPartialResult[source].Div;

			          SphP[place].DhsmlDensityFactor += 
                                    DensDataPartialResult[source].DhsmlDensity;

			          SphP[place].Rot[0] += DensDataPartialResult[source].Rot[0];
			          SphP[place].Rot[1] += DensDataPartialResult[source].Rot[1];
			          SphP[place].Rot[2] += DensDataPartialResult[source].Rot[2];
			        }
			      }
		      }

		      for(j = 0; j < NTask; j++)
		        if((j ^ ngrp) < NTask)
		          nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
		    }
	      tend = second();
	      timecommsumm += timediff(tstart, tend);

	      level = ngrp - 1;
	    }

	    MPI_Allgather(&ndone, 1, MPI_INT, ndonelist, 1, MPI_INT, MPI_COMM_WORLD);
	    for(j = 0; j < NTask; j++)
	      ntotleft -= ndonelist[j];
	  }



    /* do final operations on results */
    tstart = second();
    for(i = 0, npleft = 0; i < N_gas; i++)
  	{
  	  if(P[i].Ti_endstep == All.Ti_Current)
	    {
	      {
      		SphP[i].DhsmlDensityFactor =
         		  1 / (1 + SphP[i].Hsml * SphP[i].DhsmlDensityFactor / (NUMDIMS * SphP[i].Density));

      		SphP[i].CurlVel = sqrt(SphP[i].Rot[0] * SphP[i].Rot[0] +
    				       SphP[i].Rot[1] * SphP[i].Rot[1] +
    				       SphP[i].Rot[2] * SphP[i].Rot[2]) / SphP[i].Density;

      		SphP[i].DivVel /= SphP[i].Density;

      		dt_entr = (All.Ti_Current - (P[i].Ti_begstep + P[i].Ti_endstep) / 2) * 
                     All.Timebase_interval;

	      	SphP[i].Pressure =
	      	      (SphP[i].Entropy + SphP[i].DtEntropy * dt_entr) * pow(SphP[i].Density, 
                 GAMMA);
	      }


	      /* now check whether we had enough neighbours, adjust SphP[].Hsml accordingly */

	      if(SphP[i].NumNgb < (All.DesNumNgb - All.MaxNumNgbDeviation) ||
		        (SphP[i].NumNgb > (All.DesNumNgb + All.MaxNumNgbDeviation)
		        && SphP[i].Hsml > (1.01 * All.MinGasHsml)))
		    {
		      /* need to redo this particle */
		      npleft++;

		      if(SphP[i].Left > 0 && SphP[i].Right > 0)
		        if((SphP[i].Right - SphP[i].Left) < 1.0e-3 * SphP[i].Left)
		        {
			        /* this one should be ok */
			        npleft--;
			        P[i].Ti_endstep = -P[i].Ti_endstep - 1;	/* Mark as inactive */
			        continue;
		        }

		      if(SphP[i].NumNgb < (All.DesNumNgb - All.MaxNumNgbDeviation))
		        SphP[i].Left = dmax(SphP[i].Hsml, SphP[i].Left);
		      else
		      {
		        if(SphP[i].Right != 0)
			      {
			        if(SphP[i].Hsml < SphP[i].Right)
			          SphP[i].Right = SphP[i].Hsml;
			      }
		        else
			        SphP[i].Right = SphP[i].Hsml;
		      }

		      if(iter >= MAXITER - 10)
		      {
		        printf
			           ("i=%d task=%d ID=%d Hsml=%g Left=%g Right=%g Ngbs=%g Right-Left=%g\n   pos=(%g|%g|%g)\n",
			            i, ThisTask, (int) P[i].ID, SphP[i].Hsml, SphP[i].Left, SphP[i].Right,
			            (float) SphP[i].NumNgb, SphP[i].Right - SphP[i].Left, P[i].Pos[0], 
                  P[i].Pos[1], P[i].Pos[2]);
		        fflush(stdout);
		      }

		      if(SphP[i].Right > 0 && SphP[i].Left > 0)
		        SphP[i].Hsml = pow(0.5 * (pow(SphP[i].Left, 3) + pow(SphP[i].Right, 3)), 
                               1.0 / 3);		//Average volume between the two limits
		      else
		      {
		        if(SphP[i].Right == 0 && SphP[i].Left == 0)
		        	endrun(8188);	/* can't occur */

		        /***** this is to limit h not to over shoot *****/

		        double fac;

		        if(SphP[i].Right == 0 && SphP[i].Left > 0)
		      	{
		      	  if(P[i].Type == 0 && fabs(SphP[i].NumNgb - All.DesNumNgb) < 
                 0.5 * All.DesNumNgb)
			        {

			          fac = 1 - (SphP[i].NumNgb - All.DesNumNgb) / (NUMDIMS * SphP[i].NumNgb)
                          * SphP[i].DhsmlDensityFactor; 

			          if(fac > 1.26) 
				          fac = 1.26; 
			          if(fac < 1/1.26) 
				          fac = 1/1.26; 
			          SphP[i].Hsml *= fac; 

			          /*****   
                SphP[i].Hsml *= 1 - (SphP[i].NumNgb -
				                      All.DesNumNgb) / (NUMDIMS * SphP[i].NumNgb) * 
                              SphP[i].DhsmlDensityFactor; 
                *****/
			        }
			        else
			          SphP[i].Hsml *= 1.26;
			      }

		        if(SphP[i].Right > 0 && SphP[i].Left == 0)
			      {
			        if(P[i].Type == 0 && fabs(SphP[i].NumNgb - All.DesNumNgb) < 0.5 * All.DesNumNgb)
			        {

			          fac = 1 - (SphP[i].NumNgb - All.DesNumNgb) / (NUMDIMS * SphP[i].NumNgb) * SphP[i].DhsmlDensityFactor; 

			          if(fac > 1.26) 
			          	fac = 1.26; 
			          if(fac < 1/1.26) 
			          	fac = 1/1.26; 
			          SphP[i].Hsml *= fac; 

			          /*****			      SphP[i].Hsml *=
			               	1 - (SphP[i].NumNgb -
			               	All.DesNumNgb) / (NUMDIMS * SphP[i].NumNgb) * 
                      SphP[i].DhsmlDensityFactor; 
                *****/
			        }
			        else
			          SphP[i].Hsml /= 1.26;
			      }
		      }

		      if(SphP[i].Hsml < All.MinGasHsml)
		        SphP[i].Hsml = All.MinGasHsml;
		    }
	      else
		      P[i].Ti_endstep = -P[i].Ti_endstep - 1;	/* Mark as inactive */
	    }
	  }
    tend = second();
    timecomp += timediff(tstart, tend);


    numlist = malloc(NTask * sizeof(int) * NTask);
    MPI_Allgather(&npleft, 1, MPI_INT, numlist, 1, MPI_INT, MPI_COMM_WORLD);
    for(i = 0, ntot = 0; i < NTask; i++)
	    ntot += numlist[i];
    free(numlist);

    if(ntot > 0)
	  {
	    if(iter == 0)
	      tstart_ngb = second();

	    iter++;

	    if(iter > 0 && ThisTask == 0)
	    {
	      printf("ngb iteration %d: need to repeat for %d%09d particles.\n", 
               iter, (int) (ntot / 1000000000), (int) (ntot % 1000000000));
	      fflush(stdout);
	    }

	    if(iter > MAXITER)
	    {
	      printf("failed to converge in neighbour iteration in density()\n");
	      fflush(stdout);
	      endrun(1155);
	    }
	  }
    else
	    tend_ngb = second();
  }
  while(ntot > 0);


  /* mark as active again */
  for(i = 0; i < NumPart; i++)
    if(P[i].Ti_endstep < 0)
      P[i].Ti_endstep = -P[i].Ti_endstep - 1;

  free(ndonelist);
  free(nsend);
  free(nsend_local);
  free(nbuffer);
  free(noffset);


  /* collect some timing information */
  if(iter > 0)
    timengb = timediff(tstart_ngb, tend_ngb);
  else
    timengb = 0;

  MPI_Reduce(&timengb, &sumtimengb, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&timecomp, &sumt, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&timecommsumm, &sumcomm, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&timeimbalance, &sumimbalance, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
  {
    All.CPU_HydCompWalk += sumt / NTask;
    All.CPU_HydCommSumm += sumcomm / NTask;
    All.CPU_HydImbalance += sumimbalance / NTask;
    All.CPU_EnsureNgb += sumtimengb / NTask;
  }
}



/***********************
    density_evaluate
***********************/
void density_evaluate(int target, int mode)
{
  int j, n;
	int startnode;					// get node in tree walk 
	int numngb;							// Index used by Ngblist[]
	int numngb_inbox;				// Number of 'nearby' SphP in square box.
  double h;								// Smoothing length
	double h2;							// h^2
	double fac;
	double hinv;						// 1/h
	double hinv3;						// (1/h)^3
	double hinv4;						// (1/h)^4
  double rho, divv;
	double wk;							// Smoothing kernel, W(r,h) from eqn. 4
	double dwk;							// = d (W(r,h)) / dr
  double dx, dy, dz; 			// Displacement of SphP in x,y,z from target
	double r;								// distance SphP from target
	double r2;							// r^2
	double u;								// =r/h....used in eqn. 4, GADGET2 paper
	double mass_j;					// mass of SphP
  double dvx, dvy, dvz;
	double rotv[3];
  double weighted_numngb;	// Weighted number of neighbors (N_sph from eqn 6 rearranged)
	double dhsmlrho;				// eqn 8, I'm not sure,can't get factors to work
  FLOAT *pos;							// Position
	FLOAT *vel;							// Velocity

  if(mode == 0)
  {
    pos = P[target].Pos;
    vel = SphP[target].VelPred;
    h = SphP[target].Hsml;
  }
  else
  {
    pos = DensDataGet[target].Pos;
    vel = DensDataGet[target].Vel;
    h = DensDataGet[target].Hsml;
  }

	//Get smoothing length factors
  h2 = h * h;
  hinv = 1.0 / h;
	#ifndef  TWODIMS
  	hinv3 = hinv * hinv * hinv;
	#else
  	hinv3 = hinv * hinv / boxSize_Z;
	#endif
  hinv4 = hinv3 * hinv;

	//Initialize to 0
  rho = divv = rotv[0] = rotv[1] = rotv[2] = 0;
  weighted_numngb = 0;
  dhsmlrho = 0;

  startnode = All.MaxPart;
  numngb = 0;
  do																//while startnode > 0
  {
		//This writes the index of all nearby SphP to Ngblist[]
    numngb_inbox = ngb_treefind_variable(&pos[0], h, &startnode);

		//Itereate through List of nearby SphP
    for(n = 0; n < numngb_inbox; n++)
  	{
  	  j = Ngblist[n];								//Get nearby SphP's index

			//Get distance from target
  	  dx = pos[0] - P[j].Pos[0];
  	  dy = pos[1] - P[j].Pos[1];
  	  dz = pos[2] - P[j].Pos[2];

			/* now find the closest image in the given box size  */
			#ifdef PERIODIC			
  	  	if(dx > boxHalf_X)
  	  	  dx -= boxSize_X;
  	  	if(dx < -boxHalf_X)
  	  	  dx += boxSize_X;
  	  	if(dy > boxHalf_Y)
  	  	  dy -= boxSize_Y;
  	  	if(dy < -boxHalf_Y)
  	  	  dy += boxSize_Y;
  	  	if(dz > boxHalf_Z)
  	  	  dz -= boxSize_Z;
  	  	if(dz < -boxHalf_Z)
  	  	  dz += boxSize_Z;
			#endif
  	  r2 = dx * dx + dy * dy + dz * dz;


			/*SphP within smoothing length? Eliminates the corners of the
				box ngb_treefind_variable()'s returns to Ngblist[]. Reduces particles
				in a box to particles in a sphere. Loop will iterate over Particles themselves
				so ... r_ij = (0.0, 0.0, 0.0). */
  	  if(r2 < h2)
	    {
	      numngb++;

	      r = sqrt(r2);

	      u = r * hinv;

				//eqn. 4 in GADGET2 paper
	      if(u < 0.5)
    		{
    		  wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
    		  dwk = hinv4 * u * (KERNEL_COEFF_3 * u - KERNEL_COEFF_4);
    		}
	      else
		    {
		      wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);
		      dwk = hinv4 * KERNEL_COEFF_6 * (1.0 - u) * (1.0 - u);
		    }

	      mass_j = P[j].Mass;

	      rho += mass_j * wk;													//eqn 5 in GADGET2 paper

	      weighted_numngb += NORM_COEFF * wk / hinv3; //eqn 6, looks like N_sph

	      dhsmlrho += -mass_j * (NUMDIMS * hinv * wk + u * dwk); //eqn 8 (I can't derive)

				//Calculate local divergence and curl
	      if(r > 0)
		    {
		      fac = mass_j * dwk / r;

		      dvx = vel[0] - SphP[j].VelPred[0];
		      dvy = vel[1] - SphP[j].VelPred[1];
		      dvz = vel[2] - SphP[j].VelPred[2];

		      divv -= fac * (dx * dvx + dy * dvy + dz * dvz);

		      rotv[0] += fac * (dz * dvy - dy * dvz);
		      rotv[1] += fac * (dx * dvz - dz * dvx);
		      rotv[2] += fac * (dy * dvx - dx * dvy);
		    }
	    }
	  }
  }
  while(startnode >= 0);

	//write to SphP
  if(mode == 0)
  {
    SphP[target].NumNgb = weighted_numngb;
    SphP[target].Density = rho;
    SphP[target].DivVel = divv;
    SphP[target].DhsmlDensityFactor = dhsmlrho;
    SphP[target].Rot[0] = rotv[0];
    SphP[target].Rot[1] = rotv[1];
    SphP[target].Rot[2] = rotv[2];
  }
  else
  {
    DensDataResult[target].Rho = rho;
    DensDataResult[target].Div = divv;
    DensDataResult[target].Ngb = weighted_numngb;
    DensDataResult[target].DhsmlDensity = dhsmlrho;
    DensDataResult[target].Rot[0] = rotv[0];
    DensDataResult[target].Rot[1] = rotv[1];
    DensDataResult[target].Rot[2] = rotv[2];
  }
}



/*********************************************************************************
	Find the total number of particles within radius of P[target] 

	mode = 0   :  Local    particles, P[]
	mode = 1   :  NonLocal particles, GasPartSearchDataGet[]
**********************************************************************************/
long get_number_gas_part_nearby(FLOAT radius, int target, int mode)
{
	int startnode				= All.MaxPart;
	int n 							= 0;							//Ngblist[] index
	int j;
	long numngb					= 0; 							//Number of "nearby" gas part.
	long num_gas_nearby	= 0;							//Number of gas part w/in 'radius' of P[target]
	double dx 					= 0.0;
	double dy						= 0.0;
	double dz 					= 0.0;
	double r 						= 0.0;						// | P[target].Pos - P[j].Pos |
	FLOAT *Pos;

	/* In this case b/c boxHalf_* won't be set b/c they are set in 
		 get_number_of_gas_nearby which isn't called */
	#ifdef FEEDBACK_W_SPH_KERNAL
		#ifdef PERIODIC
  		boxSize = All.BoxSize;
  		boxHalf = 0.5 * All.BoxSize;
			#ifdef LONG_X
			  boxHalf_X = boxHalf * LONG_X;
			  boxSize_X = boxSize * LONG_X;
			#endif
			#ifdef LONG_Y
			  boxHalf_Y = boxHalf * LONG_Y;
			  boxSize_Y = boxSize * LONG_Y;
			#endif
			#ifdef LONG_Z
  			boxHalf_Z = boxHalf * LONG_Z;
		  	boxSize_Z = boxSize * LONG_Z;
			#endif
		#endif
	#endif


	if(mode == 0){
		Pos = grid_1D[target]->Pos;
	}
	else{
		Pos = GasPartSearchDataGet[target].Pos;
	}


	/* Get number of LOCAL particles within 'radius' of the P[target]'s location */
	do
	{
		numngb =	ngb_treefind_variable(&Pos[0], radius, &startnode);

		/* from Ngblist find particles w/in 'radius' */
		for(n = 0; n < numngb; n++)
		{
			j = Ngblist[n];										//Get index of 'nearby' gas part.

			dx = Pos[0] - P[j].Pos[0];	
			dy = Pos[1] - P[j].Pos[1];
			dz = Pos[2] - P[j].Pos[2];

			#ifdef PERIODIC			/*  find the closest image in the given box size  */
	    	if(dx > boxHalf_X) { 		dx -= boxSize_X; 	}
	   		if(dx < -boxHalf_X){		dx += boxSize_X;	}
	    	if(dy > boxHalf_Y) {		dy -= boxSize_Y;	}
	    	if(dy < -boxHalf_Y){		dy += boxSize_Y;	}
	    	if(dz > boxHalf_Z) {		dz -= boxSize_Z;	}
	    	if(dz < -boxHalf_Z){		dz += boxSize_Z;	}
			#endif
	    r = sqrt(dx * dx + dy * dy + dz * dz);
	
			if( r <= radius)
			{
				num_gas_nearby++;
			}
		}
	}
	while(startnode >= 0);

	return num_gas_nearby;
}



/*****************************************************************************
 This is a comparison kernel for a sort routine, which is used to group
 particles that are going to be exported to the same CPU.
******************************************************************************/
int gaspartsearch_compare_key(const void *a, const void *b)
{
  if(((struct gaspartsearchdata_in *) a)->Task < (((struct gaspartsearchdata_in *) b)->Task))
    return -1;
  if(((struct gaspartsearchdata_in *) a)->Task > (((struct gaspartsearchdata_in *) b)->Task))
    return +1;
  return 0;
}



/*! This routine is a comparison kernel used in a sort routine to group
 *  particles that are exported to the same processor.
 */
int dens_compare_key(const void *a, const void *b)
{
  if(((struct densdata_in *) a)->Task < (((struct densdata_in *) b)->Task))
    return -1;

  if(((struct densdata_in *) a)->Task > (((struct densdata_in *) b)->Task))
    return +1;

  return 0;
}
