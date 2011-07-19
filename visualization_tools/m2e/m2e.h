/******************************************************************************/
/*******************                                        *******************/
/*******************         Header file for m2e.c          *******************/
/*******************                                        *******************/
/******************************************************************************/
/*
Purpose
   Please see m2e.c.
Original
   Author: Kent Eschenberg eschenbe@psc.edu
   Organization: Pittsburgh Supercomputing Center, CMU
   Date: 2004/11/16
Current
   $Source$
   $Author$
   $Revision$
   $Date$
*/

#ifndef m2ehIncluded
#define m2ehIncluded 1

/***************************  Data from MFIX files  ***************************/

kExtern int handle kInit(0);

kExtern struct MfixReaderInfo *info kInit(0);    /* List of information items */
kExtern int nInfo kInit(0);                      /*  " number */

kExtern struct MfixReaderDataChoice *choices kInit(0);   /* List of variables */
kExtern int nChoices kInit(0);                           /*  " number */

kExtern float *xyz kInit(0);     /* Node coords (all X then all Y then all Z) */
kExtern float *xyzp kInit(0);   /*  " pseudo grid (cell centers) */

kExtern float *dataBuf kInit(0);  /* Data buffer */

kExtern int nXN, nYN, nZN, nXC, nYC, nZC, nXPN, nYPN, nZPN, nXPC, nYPC, nZPC;
kExtern int nCells[3], nNodes[3], nPCells[3], nPNodes[3];
kExtern int nXYN, nXYZN, nXYC, nXYZC, nXYPN, nXYZPN, nXYPC, nXYZPC;
kExtern int month, day, year, hours, minutes, seconds;
kExtern int versionMajor, versionMinor;
kExtern char *name, *description;

/**************************  Time management  *********************************/

kExtern float *times kInit(0);       /* Times */
kExtern int nTimes kInit(0);         /* Number of times */
kExtern int nTimesAlloc kInit(0);    /* Current number of spaces in times list */
kExtern time_t currentTime;     /* Time at which this program was run */

/*********************************  Routines  *********************************/

int writeVariable(
   FILE *file,
   char *description,
   int veclen,
   enum MfixReaderFormat format,
   float *dataF
);

int getVariable(
   char *name,
   int indexPrimary,
   int indexSecondary,
   enum MfixReaderFormat format,
   int iTime,
   int **dataInt,
   float **dataFloat
);

int getVariableInfo(
   char *varName,
   enum MfixReaderFormat format,
   int *available,
   int *nX,
   int *nY,
   int *nZ
);

int openCalcFile(
   char *varName,
   enum MfixReaderFormat format,
   int veclen,
   FILE **file
);

int calc(
   char *varName,
   enum MfixReaderFormat format,
   int *created,
   int *veclen,
   int *useTime
);

#endif

/*******************************  End of m2e.h  *******************************/
