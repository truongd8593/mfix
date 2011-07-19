/******************************************************************************/
/*******************                                        *******************/
/*******************      Kent's Library of Application     *******************/
/*******************           Support Components           *******************/
/*******************                                        *******************/
/******************************************************************************/
/*
Purpose
   Routine used by the kLib system. See kLib.h for more.
History
   Original entered by KEE 6-NOV-2001
   $Source$
   $Author$
   $Revision$
   $Date$
*/

#define kLibDeclare 1
#include "kLib.h"

static char *kLibId
   = "@(#) $Id$";

/********************************  Variables  *********************************/

static int kDepth = 0;             /* Subroutine nesting depth */
static int maxDepth = 10;          /* Maximum nesting depth shown */
static int inOutLevel = 5;         /* Make in/out log for kBugLevel >= this */

static int maxName = 32;   /* Maximim size of application and routine names */

static char *ind = "   ";   /* Indentation per depth */
static int indLen = 3;      /*  " length */

static char kLine[256];         /* Buffer used to build messages */
static int kLineBufLen = 180;   /* Max size of space in kLine for kBuf */

static int nTrace = 0;       /* Number of entries in traceback */
static int nTraceAlloc = 0;  /* Allocated space in "trace" */
static char **trace = 0;     /* List of routines for traceback */

static MsgHandler msgHandler = 0;    /* Optional message handler */

/*************************  Set the message handler  **************************/

void setMsgHandler( MsgHandler newHandler )
{
   msgHandler = newHandler;
   return;
}

/*********************************  Routine  **********************************/
/*
Code
   -3 means this is an error message
   -2 means this is a routine exit message
   -1 means this is a routine entry message
    otherwise this is an application debug message
Name
   Name of the routine from which the request has been sent.
kBuf
   Error and debug messages must have put something into this global buffer
   before calling this routine.
Safety
   Assume the worst. Check that the length of all strings that come from
   outside are within limits; if not, only use the part of the string that
   will fit. However, the program can still crash if a string pointer goes
   outside program memory; in that case, the use of "strlen" may fail as
   it walks down the the string trying to find the terminator.
Combine Strings
   Use "memcpy" and keep track of the current string length using "len". Do
   not bother with the terminating null until the entire output string has
   been built.
*/

void kReport( int code, char *name )
{
   int i, len, lenA, lenN, lenB, depth, printing;
   static char **traceOld, *errText = "error, ";

   if( kLibId ) i = 1;   /* Just to force kLibId into object code */

   if( code == -3 ) kStatus = 0;  /* Indicate an error */
   if( code == -2 ) --kDepth;     /* Update depth before print if leaving */

/****  Determine whether anything will be printed */

   if( code ==  -3 )
      printing = 1;
   else
      printing = ( ( code < 0 ) ? inOutLevel : code ) <= kBugLevel;

/****  Determine printing identation depth */

   depth = ( kDepth < 0 ) ? 0
      : ( ( kDepth > maxDepth ) ? maxDepth : kDepth);

/****  Prepare application name */

   if( ( code == -1 ) || ( printing && ( ( code == -3 ) || ( depth == 0 ) ) ) ){
      if( kAppName == 0 ) {
	 lenA = 0;
      }
      else {
	 lenA = strlen( kAppName );
	 if( lenA > maxName ) lenA = maxName;
      }
   }

/****  Prepare routine name */

   if( printing || ( code == -1 ) ) {
      if( name == 0 ) {
	 lenN = 0;
      }
      else {
	 lenN = strlen( name );
	 if( lenN > maxName ) lenN = maxName;
      }
   }

/****  Prepare message */

   if( printing && ( ( code == -3 ) || ( code >= 0 ) ) ) {
      if( kBuf == 0 ) {
	 lenB = 0;
      }
      else {
	 lenB = strlen( kBuf );
	 if( lenB > kBufLen ) lenB = kBufLen;
	 if( lenB > kLineBufLen ) lenB = kLineBufLen;
      }
   }

/****  Build error message */

   if( code == -3 ) {
      memcpy( kLine, kAppName, lenA );
      len = lenA;
      memcpy( kLine + len, "(", 1 );
      len += 1;
      memcpy( kLine + len, name, lenN );
      len += lenN;
      memcpy( kLine + len, "): ", 3 );
      len += 3;
      memcpy( kLine + len, errText, 7  );
      len += 7;
      memcpy( kLine + len, kBuf, lenB );
      len += lenB;
      kLine[ len ] = 0;
   }

/****  Build other messages */

   else if( printing ) {

   /****  Build start of string (application and/or routine names) */

      if( depth > 0 ) {
	 for( i=len=0; i<depth; ++i, len+=indLen )
	    memcpy( kLine + len, ind, indLen );

	 memcpy( kLine + len, name, lenN );
	 len += lenN;
	 memcpy( kLine + len, ": ", 2 );
	 len += 2;
      }
      else {
	 memcpy( kLine, kAppName, lenA );
	 len = lenA;
	 memcpy( kLine + len, "(", 1 );
	 len += 1;
	 memcpy( kLine + len, name, lenN );
	 len += lenN;
	 memcpy( kLine + len, "): ", 3 );
	 len += 3;
      }

   /****  Wrap up depending upon the type of message */

      switch( code ) {
	 case -2:
	    memcpy( kLine + len, "out", 3 );
	    len += 3;
	    break;
	 case -1:
	    memcpy( kLine + len, "in", 2 );
	    len += 2;
	    break;
	 default:
	    memcpy( kLine + len, kBuf, lenB );
	    len += lenB;
      }

      kLine[ len ] =0;
   }

   if( code == -1 ) ++kDepth;   /* Increase depth after message when enter */

   if( printing ) {
      if( msgHandler ) {
	 (*msgHandler)( kLine, (code==-3)?0:1 );
      }
      else {
         strcat( kLine, "\n" );
	 fputs( kLine, stderr );
      }
   }

/****  Update traceback stack */

   if( code == -1 ) {

   /****  Make room if needed */

      if( nTraceAlloc <= nTrace ) {
	 nTraceAlloc = nTrace + 20;
	 traceOld = trace;
	 trace = (char **)malloc( nTraceAlloc * sizeof(char *) );
	 memset( trace, 0, nTraceAlloc * sizeof(char *) );

	 if( traceOld != 0 ) {
	    memcpy( trace, traceOld, nTrace * sizeof(char *) );
	    free( traceOld );
	 }
      }

   /****  Compose location */

      memcpy( kLine, kAppName, lenA );
      len = lenA;
      memcpy( kLine + len, "(", 1 );
      len += 1;
      memcpy( kLine + len, name, lenN );
      len += lenN;
      memcpy( kLine + len, ")", 1 );
      len += 1;
      kLine[ len ] = 0;

   /****  Save location */

      trace[ nTrace ] = strdup( kLine );
      ++nTrace;
   }
   else if( code == -2 ) {
      --nTrace;
      if( nTrace < 0 ) nTrace = 0;
      if( trace[ nTrace ] != 0 ) free( trace[ nTrace ] );
   }

/****  Show traceback if available and just displayed an error */

   if( code == -3 ) {
      for( i = nTrace-2; i >= 0; --i ) {
	 sprintf( kLine, "   Called from %s", trace[i] );

	 if( msgHandler )
	    (*msgHandler)( kLine, 1 );
	 else
	    printf( "%s\n", kLine );
      }
   }

   return;
}

/******************************  End of kLib.c  *******************************/

