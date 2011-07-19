/******************************************************************************/
/*******************                                        *******************/
/*******************       Custom Calculator for m23        *******************/
/*******************                                        *******************/
/******************************************************************************/
/*
Purpose
   Calculate new variables from existing MFIX variables.
Original
   Author: Kent Eschenberg eschenbe@psc.edu
   Organization: Pittsburgh Supercomputing Center, CMU
   Date: 2004/11/3 17:38
Current
   $Source$
   $Author$
   $Revision$
   $Date$
Routines
   calcVar2 ... create variable 2
   calcVar1 ... create variable 1
   calc ....... create one new variable
*/

#include "kLib.h"
#include "MfixReader.h"
#include "m2e.h"

static char *M2E_CALC_Id
   = "@(#) $Id$";

/***************************  Local variables  ********************************/

static char *var1Name = "Vel_max";
static char *var2Name = "EP2_g";

/******************************************************************************/
/*******************                                        *******************/
/*******************            Create variable 2           *******************/
/*******************                                        *******************/
/******************************************************************************/
/*
Purpose:
   Create the variable "EP2_g" as tne square of "EP_g". Support all 3
   data formats (cells, nodes by average, and nodes by cell centers).
*/

static int calcVar2(
   enum MfixReaderFormat format,
   int *created,
   int *veclen,
   int *useTime
)
{
   int *I, iX, iY, iZ, nX, nY, nZ, iT, rawVarAvail;
   FILE *file;
   float a, *F;
   static char *rawVarName = "EP_g";
   kIn( calcVar2 );
   file = 0;

/****
*****  Get information about the raw variable that will be used
****/

   getVariableInfo( rawVarName,format, &rawVarAvail, &nX, &nY, &nZ );

   if( !rawVarAvail ) {
      printf( "Not creating calculated variable %s\n", var2Name );
      printf( "because %s is not available in this data set\n", rawVarName );
      goto done;
   }

/****
*****  Open output file (the "1" is the vector length)
****/

   if( !openCalcFile( var2Name, format, 1, &file ) ) goto done;

/****
*****  Loop over time
****/

   for( iT=0; iT<nTimes; ++iT ) {

   /****
   *****  Get raw variable used to calculate new variable
   *****     For this example only EP_g is needed
   *****     The 0,0 are the species and phase, meaningless for EP_g
   *****     EP_g is always float data hence always returned in F
   ****/

      getVariable( rawVarName, 0, 0, format, iT, &I, &F );
      if( !kStatus ) goto wrapup;

   /****
   *****  Calculate new variable
   *****     Values change fastest along the X axis
   *****     If EP_g was a vector the component changes fastest.
   ****/

      for( iX=0; iX<nX; ++iX ) {
	 for( iY=0; iY<nY; ++iY ) {
	    for( iZ=0; iZ<nZ; ++iZ ) {
	       a = F[ iX + iY*nX + iZ*nX*nY ];
	       F[ iX + iY*nX + iZ*nX*nY ] = a * a;
	    }
	 }
      }


   /****  Write new variable
   *****     Comment " ... " should be no more than 72 characters
   ****/

      writeVariable( file, "EP_g squared", *veclen, format, F );
      if( !kStatus ) goto wrapup;
   }

/****
*****   Success! Tell main program
****/

   *created = 1;
   *veclen = 1;
   *useTime = 1;

wrapup:
   if( file != 0 ) fclose( file );

done:
   kOut;
   return( kStatus );
}

/******************************************************************************/
/*******************                                        *******************/
/*******************            Create variable 1           *******************/
/*******************                                        *******************/
/******************************************************************************/
/*
Purpose:
   Create the variable "VEL_max" as tne maximum of "VEL_g" and "VEL_s".
   Support all 3 data formats (cells, nodes by average, and nodes by cell
   centers).
*/

static int calcVar1(
   enum MfixReaderFormat format,
   int *created,
   int *veclen,
   int *useTime
)
{
   int i, *I, iX, iY, iZ, nX, nY, nZ, iT, rawVar1Avail, rawVar2Avail;
   int arraySizeBytes;
   FILE *file;
   float a, a0, a1, a2, b, b0, b1, b2, *A, *B;
   static char *rawVar1Name = "Vel_g";
   static char *rawVar2Name = "Vel_s";
   kIn( calcVar1 );
   file = 0;
   A = 0;

   getVariableInfo( rawVar1Name, format, &rawVar1Avail, &nX, &nY, &nZ );
   getVariableInfo( rawVar2Name, format, &rawVar2Avail, &nX, &nY, &nZ );

   if( (!rawVar1Avail) || (!rawVar2Avail) ) {
      printf( "Not creating calculated variable %s\n", var1Name );
      printf( "because %s and/or %s are not available in this data set\n",
	 rawVar1Name, rawVar2Name );
      goto done;
   }

   if( !openCalcFile( var1Name, format, 3, &file ) ) goto done;
   arraySizeBytes = 3 * nX * nY * nZ * sizeof(float);
   A = (float *)malloc( arraySizeBytes );

   for( iT=0; iT<nTimes; ++iT ) {
      getVariable( rawVar1Name, 0, 0, format, iT, &I, &B );
      if( !kStatus ) goto wrapup;
      memcpy( A, B, arraySizeBytes );  /* May overwrite array on next get */

      getVariable( rawVar2Name, 0, 0, format, iT, &I, &B );
      if( !kStatus ) goto wrapup;

      for( iX=0; iX<nX; ++iX ) {
	 for( iY=0; iY<nY; ++iY ) {
	    for( iZ=0; iZ<nZ; ++iZ ) {
	       i = 3 * ( iX + iY*nX + iZ*nX*nY );
	       a0 = A[i];
	       a1 = A[i+1];
	       a2 = A[i+2];
	       b0 = B[i];
	       b1 = B[i+1];
	       b2 = B[i+2];
	       a = sqrtf( a0*a0 + a1*a1 + a2*a2 );
	       b = sqrtf( b0*b0 + b1*b1 + b2*b2 );

	       if( b > a ) {
		  A[i] = b0;
		  A[i+1] = b1;
		  A[i+2] = b2;
	       }
	    }
	 }
      }

      writeVariable( file, "maximum vector", 3, format, A );
      if( !kStatus ) goto wrapup;
   }

   *created = 1;
   *veclen = 3;
   *useTime = 1;

wrapup:
   if( file != 0 ) fclose( file );

done:
   if(A) free(A);
   kOut;
   return( kStatus );
}

/******************************************************************************/
/*******************                                        *******************/
/*******************         Create one new variable        *******************/
/*******************                                        *******************/
/******************************************************************************/
/*
Purpose:
   Create one new variable.
Return:
   created
      0: variable was not created
      1: variable was created
   veclen
      1: variable is a scalar
      3: variable is a vector
   useTime
      0: variable is static
      1: variable defined for "nTimes" listed in "times" (see m2e.h)
*/

int calc(
   char *varName,
   enum MfixReaderFormat format,
   int *created,
   int *veclen,
   int *useTime
)
{
   kIn( calc );
   if( M2E_CALC_Id ) *created = 0;  /* Default is variable not created */

   if( !strcmp( varName, var1Name ) )
      calcVar1( format, created, veclen, useTime );
   else if( !strcmp( varName, var2Name ) )
      calcVar2( format, created, veclen, useTime );
   else
      printf( "Warning: %s is unknown and will not be created\n", varName );

   kOut;
   return( kStatus );
}

/******************************  End of calc.c  *******************************/
