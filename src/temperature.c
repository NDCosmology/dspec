/************************************************
Title: temperature
Purpose: Contains functions related to calculating
         the temperature of the gas particles.
Notes:  * Follows Bertone's thesis ch. 6.
        * She does high and low density particles 
          seperately. High density are those in halos
          and low density are those not in halos. 
        * She says the cutoff is roughly delta = 10,
          so I could either run amiga or use delta.
************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include "allvars.h"
#include "proto.h"



/***********************
        get_T
***********************/
void get_T(void)
{
   // This function calculates the temperature for
   // every particle, following Bertone's thesis ch. 6.
   // Basically, if the particle is in a halo, we use eq. 2.12.
   // If not in a halo, we use eq. 1.4

   int i;
   int NumPartTot;
   //double rho_bar_partial = 0.0;
   double rho_bar;
   double T0;
   double delta;
   double t_H = 2.06e17;
   double t_heat = 5.41e-11;
   double L_eps = 1.7e-20;
   double alpha = 1.5;
   double gamma = 1.6;
   double L_cc = -7.31e-30;
   double mu = 0.588;     // Reduced mass (value from Bertone thesis ch. 2, just under 2.12) 
   double m_p = 1.67e-27; // Proton mass (kg)
   double k_B = 1.38e-23; // Boltzmass constant (J/K) 
   double G = 6.67e-11;   // Newton's constant (m^3 kg^-1 s^-2)
   double H;              // Hubble parameter (s^-1)
   double H0;             // Present value of Hubble parameter (s^-1)
   double omega_b;        // Value of omega_baryon at current redshift. See 
                          // Notes from 12/8/15
   double T;              // Bertone thesis eq. 2.12
   double V_c;            // Bertone thesis eq. 2.11
   double R_vir;          // Bertone thesis eq. 2.10

   /*// Get the total number of particles in the sim
   MPI_Allreduce(&NumPart, &NumPartTot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

   // Set up for eos temp (particles not in a halo)
   // Get the mean density
   for(i = 0; i < N_gas; i++)
   {
      rho_bar_partial += SP[i].Density;
   }

   // Need to gather the rho_bar from all other processors, then send the result
   // to all processors
   MPI_Allreduce(&rho_bar_partial, &rho_bar, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);

   rho_bar = rho_bar / (float)NumPartTot;*/

   rho_bar = 9.37e-31;  // Average (critical) density of the universe according to Ali. This is
                        // physical, g/cm^3.  Need to convert to comoving GUDs for use here

   // Convert rho_bar to comoving 
   // For whatever reason, I'm pretty sure this needs to be at z = 0. That's
   // what Bertone used, and this value Ali gave me is also at z = 0, so scaling back
   // to whatever z the snapshot is at makes no sense.  Still comoving, but a = 1.
   //rho_bar = header.time * header.time * header.time * rho_bar;

   // Convert to GUD (See notes 2/3/16)
   rho_bar = rho_bar * 1.478e21 / (LITTLE_H * LITTLE_H);

   // rho_baryon = omega_baryon * rho_crit. The above is rho_crit
   rho_bar = rho_bar * OMEGA_BARYON;

   // Get Hubble param
   if(ThisTask == 0)
   {
      H = get_hubble();
   }

   // Broadcast H
   MPI_Bcast(&H, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

   // Get value of omega_baryon at current redshift
   H0 = 100.0 * LITTLE_H;

   // Convert to s^-1. 1 Mpc = 3.086e19 km
   H0 = H0 / 3.086e19;

   // I'm not sure about this. Is it supposed to be Omega_b0 or scaled to the snapshot's
   // redshift? Basically, does the equation for T0 below take into account the redshift scaling
   // to make this Omega_b0? I don't know. Looked at the deriv for T0 eq in Theuns98, and I can't
   // figure it out. Also, this equation really only holds fully into the matter dominated era. 
   omega_b = (OMEGA_BARYON * pow(1.0 + header.redshift, 3.0)) / pow(H / H0, 2.0);

   // Check this for nans and infs
   if(isnan(omega_b) == 1)
   {
      printf("Error, omega_b is a nan!\n");
      MPI_Finalize();
      exit(EXIT_FAILURE);
   }

   if(isfinite(omega_b) == 0)
   {
      printf("Error, omega_b is infinite!\n");
      MPI_Finalize();
      exit(EXIT_FAILURE);
   }
   // Get T0
   T0 = omega_b * LITTLE_H * t_H * L_eps * pow(1.0 + header.redshift, 1.5);
   T0 = T0 / (t_heat * (2.0 + alpha));
   T0 = T0 / (1.0 - (L_cc * t_H * pow(1.0 + header.redshift, 2.5) / (LITTLE_H * t_heat)));
   T0 = pow(T0, gamma - 1.0);

   // Check for infinities and nans
   if(isnan(T0) == 1)
   {
      printf("Error, T0 is a nan!\n");
      exit(EXIT_FAILURE);
   }

   if(isfinite(T0) == 0)
   {
      printf("Error, T0 is infinite!\n");
      exit(EXIT_FAILURE);
   }


   // Get the temps

   // NOTE: There's a potential problem doing it this way (that is, using delta
   // as the criteria). That is, if the particle is not in a halo, but has delta > 10, b/c
   // then SP[i].m_vir will be 0. I could just assign it to the closest halo based on the halo
   // centers? That's a pretty rough hack, though.    
   for(i = 0; i < N_gas; i++)
   {
      delta = SP[i].Density / rho_bar;
      if(SP[i].in_halo == 0)
      //if(delta <= 10.0)
      {
         // Not in halo
         //delta = 1.0 - (rho_bar / SP[i].Density); // THIS IS WRONG, I THINK. 
         // This (the above) is Bertone's definition on pg. 14 of her thesis
         // She does not refer to the absolute value, and this term will always 
         // be less than 10.  If we're in an underdense region, where rho < rho_bar, 
         // that term will be larger than 1, which makes 1 - that negative. 
         // Binney and Tremaine use delta = rho / rho_bar - 1. That makes more sense to me.
         // Croft98 uses delta = rho / rho_bar. Gnedin and Hu, which is the paper she says she gets her
         // formula from uses: rho / rho_bar, as well. At the end of the intro, he says, "The symbol rho is
         // used to denote mass density, as usual, as well as the mass density in units of the cosmic mean (
         // i.e. rho and 1 + delta used interchangably)."  If 1 + delta is the mass density in units of the cosmic
         // mean, that means that 1 + delta = rho / rho_bar.   
         SP[i].temp = T0 * pow(1.0 + delta, gamma - 1.0);
      }

      //else if((delta > 10.0) && (SP[i].in_halo == 1))
      else if(SP[i].in_halo == 1)
      {
         // In halo

         // Convert m_vir from M_sun / h to kg
         SP[i].m_vir = (SP[i].m_vir * 1.989e30) / LITTLE_H;

         // Check converted mass for infinity and nan
         if(isnan(SP[i].m_vir) == 1)
         {
            printf("Error, m_vir for particle %d is a nan!\n", i);
            exit(EXIT_FAILURE);
         }

         if(isfinite(SP[i].m_vir) == 0)
         {
            printf("Error, m_vir for particle %d is infinite!\n", i);
            exit(EXIT_FAILURE);
         }

         // Get R_vir (Bertone thesis eq. 2.10)
         R_vir = pow((G * SP[i].m_vir) / (100.0 * H * H), 1.0 / 3.0);

         // Get V_c (Bertone thesis eq. 2.11)
         V_c = sqrt((G * SP[i].m_vir) / R_vir);

         // Get T (Bertone thesis eq. 2.12)
         T = (mu * m_p * V_c * V_c) / (2.0 * k_B);

         SP[i].temp = T;
      }

      /*else if((delta > 10.0) && (SP[i].in_halo == 0))
      {
         // Instead of quitting, I should assign it to the closest halo. This would make the most sense
         // to do when reading halos files, instead of here. That way there, I could do away with the in_halo
         // flag all together.  Just assign every particle the mass of either the halo it's in if it's in one, or 
         // the closest halo if it's not in one.
         printf("Error, high density particle not in a halo, index: %d, task %d\n", i, ThisTask);
         MPI_Finalize();
         exit(EXIT_FAILURE);
      }*/

      // Check for inf and nan
      if(isnan(SP[i].temp) == 1)
      {
         printf("Error, temp for particle %d is a nan!\n", i);
         exit(EXIT_FAILURE);
      } 

      if(isfinite(SP[i].temp) == 0)
      {
         printf("Error, temp for particle %d is infinite!\n", i);
         exit(EXIT_FAILURE);
      }
   }
}
