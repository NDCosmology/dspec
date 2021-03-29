/************************************************
Title: read_snapshot.c
Purpose:
  Check for errors. 
************************************************/
#include <stdlib.h>
#include <stdio.h>
#include "allvars.h"
#include "proto.h"
#include <mpi.h>


void error_check_header(void)
{
  //Check that the compile time options match the file 
  if(header.flag_RII == 1){
 	  #ifndef FEEDBACK
		  fprintf(stderr,"ERROR : Recompile sph2grid including FEEDBACK\n");
      exit(1); 
	  #endif
  }else{
 	  #ifdef FEEDBACK
		  fprintf(stderr,"ERROR : Recompile sph2grid without FEEDBACK\n");
      exit(1); 
	  #endif
  }


  if(header.flag_RIa == 1){
 	  #ifndef FEEDBACK
		  fprintf(stderr,"ERROR : Recompile sph2grid including FEEDBACK\n");
      exit(1); 
	  #endif
  }else{
 	  #ifdef FEEDBACK
		  fprintf(stderr,"ERROR : Recompile sph2grid without FEEDBACK\n");
      exit(1); 
	  #endif
  }



  if(header.flag_stellarage == 1){
 	  #ifndef STELLARAGE
		  fprintf(stderr,"ERROR : Recompile sph2grid including STELLARAGE\n");
      exit(1); 
	  #endif
  }else{
 	  #ifdef STELLARAGE
		  fprintf(stderr,"ERROR : Recompile sph2grid without STELLARAGE\n");
      exit(1); 
	  #endif
  }


  if(header.flag_Carbon == '1'){
 	  #ifndef CARBON
		  fprintf(stderr,"ERROR : Recompile sph2grid including CARBON\n");
      exit(1); 
	  #endif
  }else{
 	  #ifdef CARBON
		  fprintf(stderr,"ERROR : Recompile sph2grid without CARBON\n");
      exit(1); 
	  #endif
  }


  if(header.flag_Nitrogen == '1'){
	  #ifndef NITROGEN
		  fprintf(stderr,"ERROR : Recompile sph2grid including NITROGEN\n");
      exit(1); 
	  #endif
  }else{
	  #ifdef NITROGEN
		  fprintf(stderr,"ERROR : Recompile sph2grid without NITROGEN\n");
      exit(1); 
	  #endif
  }


  if(header.flag_Oxygen == '1'){
	  #ifndef OXYGEN
		  fprintf(stderr,"ERROR : Recompile sph2grid including OXYGEN\n");
      exit(1); 
    #endif	
  }else{
	  #ifdef OXYGEN
		  fprintf(stderr,"ERROR : Recompile sph2grid without OXYGEN\n");
      exit(1); 
    #endif
  }


  if(header.flag_Florine == '1'){
	  #ifndef FLORINE
		  fprintf(stderr,"ERROR : Recompile sph2grid including FLORINE\n");
      exit(1); 
	  #endif
  }else{
	  #ifdef FLORINE
		  fprintf(stderr,"ERROR : Recompile sph2grid without FLORINE\n");
      exit(1); 
	  #endif
  }


  if(header.flag_Neon == '1'){
	  #ifndef NEON
		  fprintf(stderr,"ERROR : Recompile sph2grid including NEON\n");
      exit(1); 
	  #endif	
  }else{
	  #ifdef NEON
		  fprintf(stderr,"ERROR : Recompile sph2grid without NEON\n");
      exit(1); 
	  #endif	
  }


  if(header.flag_Sodium == '1'){
	  #ifndef SODIUM
		  fprintf(stderr,"ERROR : Recompile sph2grid including SODIUM\n");
      exit(1); 
	  #endif
  }else{
	  #ifdef SODIUM
		  fprintf(stderr,"ERROR : Recompile sph2grid without SODIUM\n");
      exit(1); 
	  #endif
  }

  if(header.flag_Magnesium == '1'){
	  #ifndef MAGNESIUM
		  fprintf(stderr,"ERROR : Recompile sph2grid including MAGNESIUM\n");
      exit(1); 
	  #endif
  }else{
	  #ifdef MAGNESIUM
		  fprintf(stderr,"ERROR : Recompile sph2grid without MAGNESIUM\n");
      exit(1); 
	  #endif
  }


  if(header.flag_Aluminum == '1'){
	  #ifndef ALUMINUM
		  fprintf(stderr,"ERROR : Recompile sph2grid including ALUMINUM\n");
      exit(1); 
	  #endif
  }else{
	  #ifdef ALUMINUM
		  fprintf(stderr,"ERROR : Recompile sph2grid without ALUMINUM\n");
      exit(1); 
	  #endif
  }

  if(header.flag_Silicon == '1'){
	  #ifndef SILICON
		  fprintf(stderr,"ERROR : Recompile sph2grid including SILICON\n");
      exit(1); 
	  #endif
  }else{
	  #ifdef SILICON
		  fprintf(stderr,"ERROR : Recompile sph2grid without SILICON\n");
      exit(1); 
	  #endif
  }


  if(header.flag_Phosphorus == '1'){
	  #ifndef PHOSPHORUS
		  fprintf(stderr,"ERROR : Recompile sph2grid including PHOSPHORUS\n");
      exit(1); 
	  #endif
  }else{
	  #ifdef PHOSPHORUS
		  fprintf(stderr,"ERROR : Recompile sph2grid including PHOSPHORUS\n");
      exit(1); 
	  #endif
  }


  if(header.flag_Sulfur == '1'){
	  #ifndef SULFUR
		  fprintf(stderr,"ERROR : Recompile sph2grid including SULFUR\n");
      exit(1); 
	  #endif
  }else{
	  #ifdef SULFUR
		  fprintf(stderr,"ERROR : Recompile sph2grid without SULFUR\n");
      exit(1); 
	  #endif
  }


  if(header.flag_Chlorine == '1'){
	  #ifndef CHLORINE
		  fprintf(stderr,"ERROR : Recompile sph2grid including CHLORINE\n");
      exit(1); 
	  #endif
  }else{
	  #ifdef CHLORINE
		  fprintf(stderr,"ERROR : Recompile sph2grid without CHLORINE\n");
      exit(1); 
	  #endif
  }


  if(header.flag_Argon == '1'){
	  #ifndef ARGON
		  fprintf(stderr,"ERROR : Recompile sph2grid including ARGON\n");
      exit(1); 
	  #endif
  }else{
	  #ifdef ARGON
		  fprintf(stderr,"ERROR : Recompile sph2grid without ARGON\n");
      exit(1); 
	  #endif
  }


  if(header.flag_Potassium == '1'){
	  #ifndef POTASSIUM	
		  fprintf(stderr,"ERROR : Recompile sph2grid including POTASSIUM\n");
      exit(1); 
	  #endif
  }else{
	  #ifdef POTASSIUM	
		  fprintf(stderr,"ERROR : Recompile sph2grid without POTASSIUM\n");
      exit(1); 
	  #endif
  }


  if(header.flag_Calcium == '1'){
	  #ifndef CALCIUM	
		  fprintf(stderr,"ERROR : Recompile sph2grid including CALCIUM\n");
      exit(1); 
	  #endif
  }else{
	  #ifdef CALCIUM	
		  fprintf(stderr,"ERROR : Recompile sph2grid without CALCIUM\n");
      exit(1); 
	  #endif
  }


  if(header.flag_Scandium == '1'){
	  #ifndef SCANDIUM	
		  fprintf(stderr,"ERROR : Recompile sph2grid including SCANDIUM\n");
      exit(1); 
	  #endif
  }else{
	  #ifdef SCANDIUM	
		  fprintf(stderr,"ERROR : Recompile sph2grid without SCANDIUM\n");
      exit(1); 
	  #endif
  }


  if(header.flag_Titanium == '1'){
		#ifndef TITANIUM
		  fprintf(stderr,"ERROR : Recompile sph2grid including TITANIUM\n");
      exit(1); 
	  #endif
  }else{
		#ifdef TITANIUM
		  fprintf(stderr,"ERROR : Recompile sph2grid without TITANIUM\n");
      exit(1); 
	  #endif
  }


  if(header.flag_Vanadium == '1'){
	  #ifndef VANADIUM
		  fprintf(stderr,"ERROR : Recompile sph2grid including VANADIUM\n");
      exit(1); 
	  #endif
  }else{
	  #ifdef VANADIUM
		  fprintf(stderr,"ERROR : Recompile sph2grid without VANADIUM\n");
      exit(1); 
	  #endif
  }


  if(header.flag_Chromium == '1'){
	  #ifndef CHROMIUM
		  fprintf(stderr,"ERROR : Recompile sph2grid including CHROMIUM\n");
      exit(1); 
    #endif 
  }else{
	  #ifdef CHROMIUM
		  fprintf(stderr,"ERROR : Recompile sph2grid without CHROMIUM\n");
      exit(1); 
    #endif 
  }


  if(header.flag_Manganese == '1'){
  	#ifndef MANGANESE
		  fprintf(stderr,"ERROR : Recompile sph2grid including MANGANESE\n");
      exit(1); 
	  #endif 
  }else{
  	#ifdef MANGANESE
		  fprintf(stderr,"ERROR : Recompile sph2grid without MANGANESE\n");
      exit(1); 
	  #endif 
  }


  if(header.flag_Iron == '1'){
	  #ifndef IRON
		  fprintf(stderr,"ERROR : Recompile sph2grid including IRON\n");
      exit(1); 
	  #endif
  }else{
	  #ifdef IRON
		  fprintf(stderr,"ERROR : Recompile sph2grid without IRON\n");
      exit(1); 
	  #endif
  }


  if(header.flag_Cobalt == '1'){
	  #ifndef COBALT
		  fprintf(stderr,"ERROR : Recompile sph2grid including COBALT\n");
      exit(1); 
	  #endif
  }else{
	  #ifdef COBALT
		  fprintf(stderr,"ERROR : Recompile sph2grid  without COBALT\n");
      exit(1); 
	  #endif
  }


  if(header.flag_Nickel == '1'){
	  #ifndef NICKEL
		  fprintf(stderr,"ERROR : Recompile sph2grid including NICKEL\n");
      exit(1); 
	  #endif
  }else{
	  #ifdef NICKEL
		  fprintf(stderr,"ERROR : Recompile sph2grid without NICKEL\n");
      exit(1); 
	  #endif
  }


  if(header.flag_Copper == '1'){
	  #ifndef COPPER
		  fprintf(stderr,"ERROR : Recompile sph2grid including COPPER\n");
      exit(1); 
	  #endif
  }else{
	  #ifdef COPPER
		  fprintf(stderr,"ERROR : Recompile sph2grid without COPPER\n");
      exit(1); 
	  #endif
  }


  if(header.flag_Zinc == '1'){
	  #ifndef ZINC
		  fprintf(stderr,"ERROR : Recompile sph2grid including ZINC\n");
      exit(1); 
	  #endif 
  }else{
	  #ifdef ZINC
		  fprintf(stderr,"ERROR : Recompile sph2grid without ZINC\n");
      exit(1); 
	  #endif 
  }
}

