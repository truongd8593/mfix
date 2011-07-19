/******************************************************************************/
/*******************                                        *******************/
/*******************      Kent's Library of Application     *******************/
/*******************           Support Components           *******************/
/*******************                                        *******************/
/******************************************************************************/
/*
Purpose
   Include commonly used files and provide a library of macros and routines
   to support debugging and handle errors.
Original
   Author: Kent Eschenberg eschenbe@psc.edu
   Organization: Pittsburgh Supercomputing Center, CMU
   Date: 2001/11/6 08:00:00
Current
   $Source$
   $Author$
   $Revision$
   $Date$
Including
      Some variables in this file must be declared as "extern" everywhere
   except in one place. That one place is "kLib.c". It defines "kLibDeclare"
   before including this file. NO OTHER FILE SHOULD DO THIS!
      Other application header files included after this one may have a similar
   issue so the macros kExtern and kInit(v) are defined here. They default to
   be suitable for external variables. However, in one piece of code, the macro
   kAllocate should be defined; in that code, the macros kExtern and kInit are
   defined to actually create and initialize the variables.
Reserved Words
   This library uses several variables that you might accidentally use in the
   application. So avoid the routines, macros and
      kAppName
      kRouName
      kStatus
      kBuf
      kBugLevel
Debugging
   All calls to "kBug" are eliminated if "kBugDisable" is defined. If you
   prefer *slightly* faster code with no debug messages, define "kBugDisable"
   before including this file. If debug messages are enabled, the level of
   detail is controlled by the variable "kBugLevel". If this is zero, you get
   no debug messages. It can be changed at any time.
Using
   Put "kin( routinename );" at the top of every routine with the
   routine name as the argument. Put "kOut;" at the end of every routine.
   This will save the routine name for error and debug messages, and will
   print a notice that you have entered or left a routine if the debug level
   "kBugLevel" is 5 or above.

   At the top of the main routine, use "kMain( appname )" instead of kIn.
   This sets the application name that is then used in the other error and
   debug messages.

   Declare every routine to return an "int" and end every routine with
   "return( kStatus )". Thus, in the caller, you can easily check for errors.
   The value of "kStatus" is set to zero by "kIn" and to 1 by "kErr".

   Report an error with
         kErr1(s)   or
         kErr2(f,a)   or
         kErr3(f,a,b)
         etc.
   where
      s is a simple string
      f is a format suitable for "sprintf"
      a, b, etc. are of any type as long as they match the format you provided

   and the numeral in the name is the number of arguments. Log a debug message
   with
      kBug2(v,s)   or
      kBug3(v,f,a)   or
      kBug4(v,f,a,b)
      etc.
   where the arguments are the same as for "kErr" except for the
   addition of
      v the "verbosity" level
         0 means no messages
         1 means warnings
	 2-4 means report on major events
	 >=5 means entering and exiting of routines is reported
	 <=8 means verbose but human readable
	 >8 means very verbose and should be read with an analysis tool
	 10 is the maximum
   The message is printed only if "v" is equal to or above the current
   value of "kBugLevel". Since the latter defaults to zero, messages
   with v>0 are not printed by default. You could provide a way to set
   "kBugLevel" in your application or set it in your favorite debugger
   when you get to an area of the program that is having problems.

   The length of the message passed to kErr or kBug must not exceed the length
   of the message buffer "kBuf" (currently of length 180).
*/

#ifdef kAllocate
   #define kExtern
   #define kInit(v) = v
#else
   #define kExtern extern
   #define kInit(v)
#endif

#ifndef kLibIncluded
#define kLibIncluded 1

/*****************************  Common header files  **************************/

#include "stddef.h"
#include "stdlib.h"
#include "stdio.h"
#include "ctype.h"
#include "string.h"
#include "math.h"
#include "time.h"
#include "sys/time.h"

/********************************  Misc  **************************************/

#ifdef kLibDeclare
   int kStatus = 1;                 /* 0=error, 1=OK; reset to 1 by kIn */
   int kBugLevel = 0;               /* Level of detail for debug messages */
   char *kAppName = "application";  /* Application name; redefine in main */

   char kBuf[181];      /* Buffer used to build messages */
   int kBufLen =180;    /* Maximum number of char that will fit in buffer */
#else
   extern int kStatus;
   extern int kBugLevel;
   extern char *kAppName;

   extern char kBuf[181];
   extern int kBufLen;
#endif
/****  #ifdef kLibDeclare  ****/

void kReport( int code, char *name );  /* kLib internal routine */

/************************  Displaying messages  *******************************/
/*
   Messages are displayed by printing to stdout when no message handler
   has been registered. Otherwise, they will be sent to the handler. The
   type of message is indicated as follows:
      0 = error
      1 = information
*/

typedef void (*MsgHandler)( char *message, int type );
void setMsgHandler( MsgHandler newHandler );

/*******************************  Tracking  ***********************************/
/*
   "Tracking" is a process that sets the variable kRouName to be the name of
   the current routine. It is always set, even if debugging is off, since
   it is also used for error messages. These macros also log the entry or
   exit from a routine if debugging is enabled and the bug level is >= 5.
*/

#ifdef kBugDisable
   #define kMain( name ) \
      kAppName = #name; \
      kRouName = "main"; \
      kStatus = 1;
   #define kIn( name ) \
      static char *kRouName = #name; \
      kStatus = 1
   #define kOut
#else
   #define kMain( name ) \
      static char *kRouName = "main"; \
      kAppName = #name; \
      kStatus = 1;\
      kReport( -1, kRouName )
   #define kIn( name ) \
      static char *kRouName = #name; \
      kStatus = 1; \
      kReport( -1, kRouName )
   #define kOut \
       kReport( -2, kRouName )
#endif
/****  #ifdef kBugDisable  ****/

/********************************  Errors  ************************************/

#define kErr1( msg ) { \
   strcpy( kBuf, #msg ); \
   kReport( -3, kRouName ); \
}

#define kErr2( format, a ) { \
   sprintf( kBuf, #format, a ); \
   kReport( -3, kRouName ); \
}

#define kErr3( format, a1, a2 ) { \
   sprintf( kBuf, #format, a1, a2 ); \
   kReport( -3, kRouName ); \
}

#define kErr4( format, a1, a2, a3 ) { \
   sprintf( kBuf, #format, a1, a2, a3 ); \
   kReport( -3, kRouName ); \
}

#define kErr5( format, a1, a2, a3, a4 ) { \
   sprintf( kBuf, #format, a1, a2, a3, a4 ); \
   kReport( -3, kRouName ); \
}

#define kErr6( format, a1, a2, a3, a4, a5 ) { \
   sprintf( kBuf, #format, a1, a2, a3, a4, a5 ); \
   kReport( -3, kRouName ); \
}

#define kErr7( format, a1, a2, a3, a4, a5, a6 ) { \
   sprintf( kBuf, #format, a1, a2, a3, a4, a5, a6 ); \
   kReport( -3, kRouName ); \
}

#define kErr8( format, a1, a2, a3, a4, a5, a6, a7 ) { \
   sprintf( kBuf, #format, a1, a2, a3, a4, a5, a6, a7 ); \
   kReport( -3, kRouName ); \
}

#define kErr9( format, a1, a2, a3, a4, a5, a6, a7, a8 ) { \
   sprintf( kBuf, #format, a1, a2, a3, a4, a5, a6, a7, a8 ); \
   kReport( -3, kRouName ); \
}

#define kErr10( format, a1, a2, a3, a4, a5, a6, a7, a8, a9 ) { \
   sprintf( kBuf, #format, a1, a2, a3, a4, a5, a6, a7, a8, a9 ); \
   kReport( -3, kRouName ); \
}

/********************************  Debug  *************************************/

#ifdef kBugDisable
   #define kBug2( v, msg )
   #define kBug3( v, format, a )
   #define kBug4( v, format, a1, a2 )
   #define kBug5( v, format, a1, a2, a3 )
   #define kBug6( v, format, a1, a2, a3, a4 )
   #define kBug7( v, format, a1, a2, a3, a4, a5 )
   #define kBug8( v, format, a1, a2, a3, a4, a5, a6 )
   #define kBug9( v, format, a1, a2, a3, a4, a5, a6, a7 )
   #define kBug10( v, format, a1, a2, a3, a4, a5, a6, a7, a8 )
#else
   #define kBug2( v, msg ) { \
      strcpy( kBuf, #msg ); \
      kReport( v, kRouName ); \
   }

   #define kBug3( v, format, a ) { \
      sprintf( kBuf, #format, a ); \
      kReport( v, kRouName ); \
   }

   #define kBug4( v, format, a1, a2 ) { \
      sprintf( kBuf, #format, a1, a2 ); \
      kReport( v, kRouName ); \
   }

   #define kBug5( v, format, a1, a2, a3 ) { \
      sprintf( kBuf, #format, a1, a2, a3 ); \
      kReport( v, kRouName ); \
   }

   #define kBug6( v, format, a1, a2, a3, a4 ) { \
      sprintf( kBuf, #format, a1, a2, a3, a4 ); \
      kReport( v, kRouName ); \
   }

   #define kBug7( v, format, a1, a2, a3, a4, a5 ) { \
      sprintf( kBuf, #format, a1, a2, a3, a4, a5 ); \
      kReport( v, kRouName ); \
   }

   #define kBug8( v, format, a1, a2, a3, a4, a5, a6 ) { \
      sprintf( kBuf, #format, a1, a2, a3, a4, a5, a6 ); \
      kReport( v, kRouName ); \
   }

   #define kBug9( v, format, a1, a2, a3, a4, a5, a6, a7 ) { \
      sprintf( kBuf, #format, a1, a2, a3, a4, a5, a6, a7 ); \
      kReport( v, kRouName ); \
   }

   #define kBug10( v, format, a1, a2, a3, a4, a5, a6, a7, a8 ) { \
      sprintf( kBuf, #format, a1, a2, a3, a4, a5, a6, a7, a8 ); \
      kReport( v, kRouName ); \
   }
#endif
/****  #ifdef kBugDisable  ****/

#endif
/****  #ifndef kLibIncluded  ****/

/******************************  End of kLib.h  *******************************/
