/******************************************************************************/
/*******************                                        *******************/
/*******************       Routines to Read MFIX Files      *******************/
/*******************                                        *******************/
/******************************************************************************/
/*
Purpose
   Read a set of files created by the MFIX simulation. See "MfixReader.h" for
   more.
Original
   Author: Kent Eschenberg eschenbe@psc.edu
   Organization: Pittsburgh Supercomputing Center, CMU
   Date: 2001/11/6 08:00:00
Current
   $Source$
   $Author$
   $Revision$
   $Date$
Routines
   Manage buffers
      getReturnBuffer ................. get a buffer used to return data
      getInputBuf ..................... get the input buffer
      releaseBuf ...................... release one or all buffers
      cleanupBuf ...................... free unused buffers
      updateBuf ....................... update buffer size
   Utilities
      trimTrailingBlanks .............. replace trailing blanks with nulls
      findInstance .................... find an instance given the handle
      interpolate ..................... calculate data between saved times
      makeCenterCoord ................. calculate cell center coordinates
      makeNodetype .................... convert celltypes to nodetypes
      cell2NodeByAvg .................. convert cell data to node data by avg
   Low level input
      readRecord ...................... read one record
      readArray ....................... read an array
      skipRecords ..................... skip records
      getString ....................... get a string from a record
      getFloat ........................ get a float from a record
      getInt .......................... get an integer from a record
      readHeader ...................... read the file header
      openFile ........................ compose name then open file
   Read the data files
      skipToVariable .................. skip to a variable
      readVariable .................... read and filter one variable at one time
      selectTime ...................... select data time(s) to use
      readSpTimes ..................... gather time information for one sp file
      readSpHeader .................... read, check, store header from 1 sp file
      scanSpFiles ..................... partially read all sp files
      calcLocations ................... calculate locations in sp files
   Read and check the restart file
      readKEpsilon .................... read K Epsilon
      readNumberReactionRates ......... read number of reaction rate tables
      readNumberScalars ............... read number of scalars
      skipArrays2 ..................... skip more arrays
      getCellType ..................... get cell type
      readCellType .................... read cell types
      skipArrays1 ..................... skip a seris of arrays
      convertCoordinates .............. convert coordinate storage and system
      readRunInfo ..................... read run information
      readOneCoordinate ............... read one array of coordinates
      getCoordinates .................. read and filter coordinates
      readNumberSpecies ............... read number of species
      readMiscDim...................... read various dimensions
      readNumberSolidPhases ........... read number of solid phases
      updateInternals ................. apply filter and allocate input buffer
      readArrayDim .................... read array dimensions
      readRestartHeader ............... read, echo and check header
      readRestartFile ................. overall read logic
   Retreive information after open
      MFixReaderGetData ............... get data
      MfixReaderGetDataChoices ........ get list of available data
      MFixReaderGetCoordinates ........ get coordinates
      MfixReaderGetInfo ............... get information
   Options
      MfixReaderSetTimeRules .......... set the time selection rules
      MFixReaderSetLocationFilter ..... set the filter using locations
      MfixReaderSetDataTimeTolerance .. set the data time tolerance
      MfixReaderSetInputEndian ........ set the input file format
      MFixReaderSetRecordLength ....... set the input file record length
      MfixReaderSetFileName ........... set the input file name
   Main operations
      MfixReaderClose ................. close a set of files
      MFixReaderOpen .................. open the set of files
      MfixReaderDestroy ............... destroy a reader instance
      MfixReaderCreate ................ create a new reader instance
MFIX Version
   The MFIX file version, such 1.02 or 1.4, will be stored as two
   integer, such as (1,2) or (1,40), respectively.
Buffers
   The buffers returned by MfixReaderGetData are shared by all instances. 
   They should never be freed by the application. The application should
   assume that the buffers will change during any call to MfixReader. At
   the end of the routine MfixReaderGetData, all unused buffers will be
   freed except that at least one of each of the 5 sizes will be retained.
   The size used will be large enough for all instances. Additional buffers
   will lbe allocated as needed. There is only one input buffer.
Interpolation
   An interpolation rule and a difference threshold are used to determine
   whether interpolation over time is needed. If it is, linear interpolation
   is used. Each component of a vector is interpolated separately.
Cylindrical to Cartesian
   When the input file uses a cylindrical coordinate system, all values are
   converted to cartesian. Vectors are converted using the angle at the
   center of the cell which is the average of the angles of the nodes.
Cell to Node
   When node data is requested, it is formed by taking the average of the
   values at the cells around each node. Vector components are averaged
   separately. If the input file uses a cylindrical coordinate system, it
   is assumed that the first X coordinate is 0 and is on the central axis;
   and that the first and last angle are the same (i.e., 0 and 360 degrees).
   This assumption requires that the first cell in X and both the first and
   last cells in Z (boundary cells) be clipped.
Node Flags
   If requested, cell flags are converted to node flags following the rule
   that a node is "fluid" (a value of 1) if at least one cell around it is
   fluid; otherwise, it is "wall" (a value of 100). For example, all of the
   nodes around a single isolated non-fluid cell will be flagged fluid. This
   approach maps correctly into the valid/invalid flag system used by many
   visualization packages.
Record Length
   Currently MFIX only uses the first 512 bytes of a record, leaving the rest
   of a longer record unused. The "record length" parameter refers to the
   actual record length; if over 512, the gaps will be skipped.
*/

/*
Large File Support
   May work only with the GNU compilers and glibc.
   The following flag redefines some I/O features to use 64-bits.
   This code has been changed to also explicitly use fopen64 and fseeko64.
   As a result 32-bit systems to read files over 2GB.
   THE FOLLOWING FLAG MUST COME BEFORE INCLUDING OTHER FILES!
*/
#define _LARGEFILE64_SOURCE 1

#include "kLib.h"
#include "bio.h"
#include <fcntl.h>
#include <errno.h>
#include "MfixReader.h"

static char *MfixReaderId
   = "@(#) $Id$";

/************************  Externally visible variables  **********************/

char *MfixReaderCoordText[3] = { "unknown", "cartesian", "cylindrical" };
char *MfixReaderDataTypeText[3] = { "integer", "float", "string" };
char *MfixReaderAttachText[2] = { "cell", "node" };

char *MfixReaderFormatText[4] = { "no requirement", "cell data, no conversion",
   "node data, average method", "node data, cell-center method" };

char *MfixReaderInterpolationText[5] = { "never", "nearest", "prev", "next",
   "linear" };

char *MfixReaderExtrapolationText[2] = { "never", "nearest" };

/*********************************  Misc  *************************************/

#define nSpFiles 11   /* Number of sp files (sp1, sp2, ..., spb) */

static char *fileExt[ nSpFiles ] = { ".sp1", ".sp2", ".sp3", ".sp4", ".sp5",
   ".sp6", ".sp7", ".sp8", ".sp9", ".spa", ".spb" };

static float SMALL_TIME_DIFFERENCE = 0.000001;  /* Times are same if closer */

static int bioInited = 0;           /* Bio package has been initialized if 1 */
static char *coordUnitName = "cm";  /* Units for coordinates */

/* Variables needed only while reading the restart file */
static int tmpDimC, tmpDimIC, tmpDimBC, tmpDimIS;

/*******************************  Limits  *************************************/

static int majorMin = 1;      /* Acceptible major version numbers */
static int majorMax = 1;

static int minorMin = 0;      /* Acceptible minor version numbers */
static int minorMax = 60;

static int nXYZMax = 10000;   /* Maximum number of values along any one axes */

static int maxNTimes = 100000; /* Maximum number of times in one file */

/*********************  Information about the file set  ***********************/

#define nInfo 19   /* Number of information items */

static struct MfixReaderInfo infoTemplate[ nInfo ] = {
   { "dims",          MRD_integer, 3,0,0,0 }, /* Num of points along X, Y, Z */
   { "dimsoriginal",  MRD_integer, 3,0,0,0 }, /*  " before location filter */
   { "ntimes",        MRD_integer, 1,0,0,0 }, /* Number of time steps */
   { "timespan",      MRD_float,   2,0,0,0 }, /* Earliest & latest time (sec) */
   { "date",          MRD_integer, 3,0,0,0 }, /* Month, day, year */
   { "time",          MRD_integer, 3,0,0,0 }, /* Hour, minute, second */
   { "version",       MRD_integer, 2,0,0,0 }, /* Major, minor version */
   { "name",          MRD_string,  2,0,0,0 }, /* First, second run name */
   { "type",          MRD_string,  1,0,0,0 }, /* Type of run */
   { "description",   MRD_string,  1,0,0,0 },
   { "units",         MRD_string,  1,0,0,0 },
   { "coord",         MRD_integer, 1,0,0,0 }, /* Coordinate system */
   { "ngasspecies",   MRD_integer, 1,0,0,0 }, /* Called nmax(0) in MFIX */
   { "nsolidphases",  MRD_integer, 1,0,0,0 }, /* Called mmax in MFIX */
   { "nsolidspecies", MRD_integer, 0,0,0,0 }, /* Length nsolidphases; nmax() */
   { "rank",          MRD_integer, 1,0,0,0 }, /* 1/2/3 for 1D/2D/3D */
   { "nscalars",      MRD_integer, 1,0,0,0 }, /* Number of scalars */
   { "nreactionrates",MRD_integer, 1,0,0,0 }, /* Number of reaction rates */
   { "kepsilon",      MRD_integer, 1,0,0,0 }  /* Variable K_Epsilon */
};

/*******************************  Data choices  *******************************/

#define nChoices 16

/* Sp file number (0,1,...) for each data choice; -1 means the restart file */
static int data2file[ nChoices ]
   = { -1, 0, 1, 1, 2, 3, 4, 5, 5, 6, 6, 7, 8, 9, 10, 10 };

static struct MfixReaderDataChoice choicesTemplate[ nChoices ] = {
   /*[ 0 res ]*/ { "celltype", "cell types", 0, 0, 0, 0,
                   "enumeration", MRD_integer, MRA_cell, 1, 0, 0, 0 },
   /*[ 1 sp1 ]*/ { "EP_g", "void fraction", 0, 0, 0, 0,
                   "fraction", MRD_float, MRA_cell, 1, 0, 0, 0 },
   /*[ 2 sp2 ]*/ { "P_g", "gas pressure", 0, 0, 0, 0,
                   "dynes/cm^2", MRD_float, MRA_cell, 1, 0, 0, 0 },
   /*[ 3 sp2 ]*/ { "P_star", "fluid pressure", 0, 0, 0, 0,
                   "dynes/cm^2", MRD_float, MRA_cell, 1, 0, 0, 0 },
   /*[ 4 sp3 ]*/ { "Vel_g", "gas velocity vector", 0, 0, 0, 0,
                   "cm/sec", MRD_float, MRA_cell, 3, 0, 0, 0 },
   /*[ 5 sp4 ]*/ { "Vel_s", "solid velocity vector", "phase", 1, 0, 0,
                   "cm/sec", MRD_float, MRA_cell, 3, 0, 0, 0 },
   /*[ 6 sp5 ]*/ { "ROP_s", "solids density", "phase", 1, 0, 0,
                   "gm/cm^2", MRD_float, MRA_cell, 1, 0, 0, 0 },
   /*[ 7 sp6 ]*/ { "T_g", "fluid temperature", 0, 0, 0, 0,
                   "degrees Kelvin", MRD_float, MRA_cell, 1, 0, 0, 0 },
   /*[ 8 sp6 ]*/ { "T_s", "solid temperature", "phase", 1, 0, 0,
                   "degrees Kelvin", MRD_float, MRA_cell, 1, 0, 0, 0 },
   /*[ 9 sp7 ]*/ { "X_g", "gas species mass fraction", "species", 1, 0, 0,
                   "fraction", MRD_float, MRA_cell, 1, 0, 0, 0 },
   /*[10 sp7 ]*/ { "X_s", "solid species mass fraction", "phase", 1,"species",0,
                   "fraction", MRD_float, MRA_cell, 1, 0, 0, 0 },
   /*[11 sp8 ]*/ { "THETA_m", "granular temperature", "phase", 1, 0, 0,
                   "cm^2/sec^2", MRD_float, MRA_cell, 1, 0, 0, 0 },
   /*[12 sp9 ]*/ { "scalar", "scalar variable", "which", 1, 0, 0,
                   "unknown", MRD_float, MRA_cell, 1, 0, 0, 0 },
   /*[13 spa ]*/ { "r-rate", "reaction rates", "which", 1, 0, 0,
                   "unknown", MRD_float, MRA_cell, 1, 0, 0, 0 },
   /*[14 spb ]*/ { "K_Turb_G", "gas turbulance kinetic energy", 0, 0, 0, 0,
                   "cm^2/sec^2", MRD_float, MRA_cell, 1, 0, 0, 0 },
   /*[15 spb ]*/ { "E_Turb_G", "dissipation of K_Turb_G", 0, 0, 0, 0,
                   "cm^2/sec^2", MRD_float, MRA_cell, 1, 0, 0, 0 }
};

/***********************  Parameters for one instance  ************************/

struct FileSetStr {
   int handle;              /* NOT the index in "fileSets" */
   int opened;              /* 1=opened 0=closed */

   char *baseFileName;    /* Base name without extension */
   int big_endian_input;  /* Set to 1 when file stored in big-endian form */

   char *resFileName;     /* Restart file name */
   FILE *resFile;         /*  " file pointer (0 when closed) */

   char *spFileNames[nSpFiles];  /* Sp file names */
   FILE *spFiles[nSpFiles];      /*  " file pointer (0 when closed) */
   int spRecPerTime[nSpFiles];   /* Number of records per time step */

   int spNTimes[nSpFiles];       /* Number of time steps */
   float *spTimes[nSpFiles];     /* Times (seconds) */

   int spNTimesA[nSpFiles];      /* Number of acceptible time steps */
   float *spTimesA[nSpFiles];    /* Acceptible times (dups eliminated) */
   int *spITimesA[nSpFiles];     /* Index of corresponding acceptible time */

   int recl;                 /* Bytes per record */
   int reclStored;           /*  " stored in file (>=recl, extra not used) */
   int reclGap;              /* Gap size = reclStored - recl */

   int intPerRecord;         /* Number of 4-byte integers per record */
   int floatPerRecord;       /* Number of 4-byte floats per record */
   int doublePerRecord;      /* Number of 8-byte doubles per record */
   int recordPerFloatArray;  /* Number of records per float array of size nXYZ*/

   int skip[6];           /* Number of cells to skip at the ends */
   int rank;              /* 1/2/3 for 1D/2D/3D */
   int x0, x1;            /* 1st, last X cell index, inclusive, after filter */
   int y0, y1;            /*  " Y */
   int z0, z1;            /*  " Z */
   int nXC, nYC, nZC;     /* Number of cells along each axis */
   int nXN, nYN, nZN;     /*  " nodes */
   int nXCB, nYCB, nZCB;  /* Number of cells along each axis after filter */
   int nXNB, nYNB, nZNB;  /*  " nodes */
   int nXYZC;             /* Number of cells in 3D array, nX*nY*nZ */
   int nXYZCB;            /* Number of filtered cells in 3D array, nXB*nYB*nZB*/
   int nCells[3];         /* Number of original cells */
   int nCellsB[3];        /* Number of cells after filter */
   int nNodes[3];         /* Number of nodes */
   int nNodesB[3];        /* Number of nodes after filter */

   int nTimes;               /* Largest number of times in any one sp file */
   float timeSpan[2];        /* Earliest and latest time in any sp file */
   float tolerance;          /* Data time same as requested time if under */
   float dataTimeTolerance;  /* Ignores out-of-order times if under */

   enum MfixReaderInterpolation interpRule;
   enum MfixReaderExtrapolation extrapRule;

   char *runName[2], *runType, *desc, *units;
   int version[2], runDate[3], runTime[3];

   int nScalars;         /* Number of scalars stored in SP9 file */
   int nSolidPhases;     /* Number of solid phases (mmax)*/
   int *nSolidSpecies;   /* Number of species in solid phases (nmax[])*/
   int nGasSpecies;      /* Number of species in gas phases (nmax[0])*/
   int nReactionRates;   /* Number of reaction rates stored in SPA file */
   int KEpsilon;         /* Boolean K_Epsilon; SPB used only if True */

   struct MfixReaderInfo info[ nInfo ];             /* Used by GetInfo */
   struct MfixReaderDataChoice choices[ nChoices ]; /* Used by GetDataChoices */

   int varSkip[ nChoices ]; /* Number of variables in the file that preceed */
                            /* the variable */

   enum MfixReaderCoord coord;   /* Coordinate system */

   float *xyz;          /* Coord after filter; length (nXB+1)*(nYB+1)*(nZB+1) */
   float *xyzcc;        /*  " cell centers; length nXB*nYB*nZB */
   float *gridAngle;    /* Angle of cell centers */

   enum MfixReaderCellType *cellType;  /* Cell type (FLAG) */
};

typedef struct FileSetStr FileSet;

static FileSet **fileSets = 0;  /* List of filesets */
static int nFileSets = 0;              /*  " number */
static int nextHandle = 314159;        /* Next unique handle */

static FileSet *F;  /* Temporary pointer to current instance */

/*********************************  Buffers  **********************************/

static int currentRecord = 0;     /* Current record in file (1,2,...) */
static int reclDefault=512;       /* Default record length */
static int recl = 0;              /*  " current length */
static unsigned char *record = 0; /* Tmp storage for one record, size recl */

/*
The list of buffers contains buffers of 5 types: input, scalar cell, vector
cell, scalar node and vector node. Cell buffers are long enough to hold a float
(scalar) or 3 floats (vector) per cell. Node buffers are long enough to hold one
or 3 floats per node.
*/

struct BuffersStr {
   int used;     /* Buffer in use (boolean); always -1 for input buffer */
   int len;      /* Length of buffer (floats) */
   float *buf;   /* Buffer */
};

typedef struct BuffersStr Buffer;

static Buffer *buffers = 0;       /* List of buffer descriptions */
static int nBuffers = 0;        /*  " number */
static int nBuffersAlloc = 0;   /*  " number allocated */

static int nScaCellBuf = 0;   /* Number of scalar buffers for cells */
static int nVecCellBuf = 0;   /*  " vector cells */
static int nScaNodeBuf = 0;   /*  " scalar nodes */
static int nVecNodeBuf = 0;   /*  " vector nodes */

static int lenInputBuf = 0; /* Length of input buffer for cells before filter */

static int lenScaCellBuf = 0;   /* Length of scalar buffers for cells (bytes) */
static int lenVecCellBuf = 0;   /*  " vector cells */
static int lenScaNodeBuf = 0;   /*  " scalar nodes */
static int lenVecNodeBuf = 0;   /*  " vector nodes */

/******************************************************************************/
/*******************                                        *******************/
/*******************            Manage buffers              *******************/
/*******************                                        *******************/
/******************************************************************************/

/*******************
********************  Get a buffer used to return data
*******************/

static int getReturnBuffer( int len, float **buf )
{
   int i;
   Buffer *buf2;
   kIn( getReturnBuffer );

   for( i=0; i<nBuffers; ++i ) {
      if( buffers[i].used ) continue;
      if( buffers[i].len == len ) break;
   }

   if( i < nBuffers ) {
      *buf = buffers[i].buf;
      buffers[i].used = 1;
   }
   else {
      ++nBuffers;

      if( nBuffers > nBuffersAlloc ) {
	 buf2 = buffers;
	 nBuffersAlloc += 20;
	 buffers = (Buffer *)malloc( nBuffersAlloc * sizeof(Buffer) );

	 if( nBuffers > 1 ) {
	    memcpy( buffers, buf2, ( nBuffers - 1 ) * sizeof(Buffer) );
	    free( buf2 );
	 }
      }

      buffers[ nBuffers - 1 ].used = 1;
      buffers[ nBuffers - 1 ].len = len;
      *buf = (float *)malloc( len * sizeof(float) );
      buffers[ nBuffers - 1 ].buf = *buf;
   }

   kOut;
   return( kStatus );
}

/*******************
********************  Get the input buffer
*******************/

static int getInputBuf( float **buf )
{
   int i;
   kIn( getInputBuf );

   for( i=0; i<nBuffers; ++i ) {
      if( buffers[i].used < 0 ) break;
   }

   if( i >= nBuffers ) {
      kErr1( input buffer is missing );
      *buf = 0;
   }
   else {
      *buf = buffers[i].buf;
   }

   kOut;
   return( kStatus );
}

/*******************
********************  Release one or all buffers 
*******************/

static int releaseBuf( float *buf )
{
   int i;
   kIn( releaseBuf );

   if( buf ) {
      for( i=0; i<nBuffers; ++i ) {
	 if( buffers[i].buf != buf ) continue;
	 if( buffers[i].used > 0 ) buffers[i].used = 0;
	 break;
      }
   }
   else {
      for( i=0; i<nBuffers; ++i ) {
	 if( buffers[i].used > 0 ) buffers[i].used = 0;
      }
   }

   kOut;
   return( kStatus );
}

/*******************
********************  Free unused buffers
********************  Reduce number of each type to lesser of number used or 1
*******************/

static int cleanupBuf( void )
{
   int i, j, n, nU, len, lenT[4], nT[4], nTU[4];
   kIn( cleanupbuf );

/****  Determine number of buffers used of each type */

   lenT[0] = lenScaCellBuf;
   lenT[1] = lenVecCellBuf;
   lenT[2] = lenScaNodeBuf;
   lenT[3] = lenVecNodeBuf;

   memset( nT,  0, 4 * sizeof(int) );
   memset( nTU, 0, 4 * sizeof(int) );

   for( i=0; i<nBuffers; ++i ) {
      if( buffers[i].used < 1 ) continue;
      len = buffers[i].len;

      for( j=0; j<4; ++j ) {
	 if( len != lenT[i] ) continue;
	 ++nT[i];
	 if( buffers[i].used ) ++nTU[i];
	 break;
      }
   }

/****  Reduce number of buffers of each type to lower of used or 1 */

   for( j=0; j<4; ++j ) {
      len = lenT[j];
      n = nT[j];
      nU = nTU[j];

      while( ( n > 1 ) && ( n > nU ) ) {
	 for( i=0; i<nBuffers; ++i ) {
	    if( buffers[i].len != len ) continue;
	    if( !buffers[i].used ) break;
	 }

	 free( buffers[i].buf );

	 if( i < ( nBuffers - 1 ) )
	    memmove( buffers + i, buffers + i + 1,
	       ( nBuffers - i - 1 ) * sizeof(Buffer) );

	 --nBuffers;
	 --n;
      }
   }

   kOut;
   return( kStatus );
}

/*******************
********************  Update buffer size
*******************/

static int updateBuf( void )
{
   int i, n, nXYZC, nXYZCB, nXYZNB, change;
   FileSet *f;
   kIn( updateBuf );

/****  Recalculate buffer sizes using largest open grid */

   nXYZC = nXYZCB = nXYZNB = 0;

   for( i=0; i<nFileSets; ++i ) {
      f = fileSets[i];

      if( i == 0 ) {
	 nXYZC = f->nXYZC;
	 nXYZCB = f->nXYZCB;
	 nXYZNB = f->nXNB * f->nYNB * f->nZNB;
      }
      else {
	 n = f->nXNB * f->nYNB * f->nZNB;
	 if( nXYZC  < f->nXYZC  ) nXYZC  = f->nXYZC;
	 if( nXYZCB < f->nXYZCB ) nXYZCB = f->nXYZCB;
	 if( nXYZNB < n         ) nXYZNB = n;
      }
   }

/****  Change buffer size if needed */

   change = 0;

   if( nXYZC > lenInputBuf )
      change = 1;
   else if( nXYZC < ( lenInputBuf - 1000 ) )
      change = 1;
   else if( nXYZCB > lenScaCellBuf )
      change = 1;
   else if( nXYZCB < ( lenScaCellBuf - 1000 ) )
      change = 1;
   else if( nXYZNB > lenScaNodeBuf )
      change = 1;
   else if( nXYZNB < ( lenScaNodeBuf - 1000 ) )
      change = 1;

   if( !change ) goto done;

   /****  Change needed - first free old buffers */

   for( i=0; i<nBuffers; ++i ) free( buffers[i].buf );
   nBuffers = 0;

/****  Initialize main buffer management array */

   if( nBuffersAlloc < 1 ) {
      nBuffersAlloc = 20;
      buffers = (Buffer *)malloc( nBuffersAlloc * sizeof(Buffer) );
   }

/****  Initialize buffer system */

   if( nXYZC > 0 ) {
      lenInputBuf = nXYZC;
      lenScaCellBuf = nXYZCB;
      lenVecCellBuf = 3 * nXYZCB;
      lenScaNodeBuf = nXYZNB;
      lenVecNodeBuf = 3 * nXYZNB;
      nBuffers = 1;
      buffers[0].used = -1;
      buffers[0].len = lenInputBuf;
      buffers[0].buf = (float *)malloc( lenInputBuf * sizeof(float) );
   }
   else {
      nScaCellBuf = nVecCellBuf = nScaNodeBuf = nVecNodeBuf = 0;
   }

done:
   kOut;
   return( kStatus );
}

/******************************************************************************/
/*******************                                        *******************/
/*******************                 Utilities              *******************/
/*******************                                        *******************/
/******************************************************************************/

/*******************
********************  Replace trailing blanks with nulls
*******************/

static int trimTrailingBlanks( char *s, int n )
{
   int i;
   kIn( trimTrailingBlanks );

   s[n] = 0;

   for( i=n-1; i>=0; --i ) {
      if( s[i] != ' ' ) break;
      s[i] = 0;
   }

   kOut;
   return( kStatus );
}

/*******************
********************  Find an instance given the handle
*******************/

static int findInstance( int handle )
{
   int i;
   kIn( findInstance );

   for( i=0, F=0; i<nFileSets; ++i ) {
      kBug5( 7, comparing set %d handle %d to desired handle %d, i,
	 fileSets[i]->handle, handle );
      if( fileSets[i]->handle != handle ) continue;
      kBug2( 7, matched );
      F = fileSets[i];
      break;
   }

   if( i >= nFileSets ) {
      kErr2( handle %d is invalid, handle );
      kStatus = 0;
   }

   kOut;
   return( kStatus );
}

/*******************
********************  Calculate data between saved times
*******************/
/*
Purpose
   Perform linear interpolation between values at t0 (a0) and at t1
   (a1) to estimate values at t (a). Use the values from t0 or t1
   without interpolation if either is close enough to t.
Extrapolation
   The methods used here will not extrapolate, i.e., compute values for
   t where t is outside (t0,t1). Instead, the nearest value is used.
Formula
             t  - t0                            t  - t0           t  - t0
   at = a0 + ------- * ( a1 - a0 ) = a0 * ( 1 - -------  ) + a1 * -------
             t1 - t0                            t1 - t0           t1 - t0

                                          t  - t0             t  - t0
      = a0 * a + a1 * b   where   a = 1 - -------   and   b = -------
                                          t1 - t0             t1 - t0
*/

static int interpolate(
   float t,
   float t0,
   float t1,
   int veclen,
   float **A0,
   float **A1
)
{
   int i, n;
   float a, b, *a0, *a1, difference;
   kIn( interpolate );

/****  Handle extreme cases where interpolation is of little value */

   difference = t1 - t0;

   if( difference < SMALL_TIME_DIFFERENCE ) {
      releaseBuf( *A1 );
      goto done;
   }

   if( ( t - t0 ) < SMALL_TIME_DIFFERENCE ) {
      releaseBuf( *A1 );
      goto done;
   }

   if( ( t1 - t ) < SMALL_TIME_DIFFERENCE ) {
      releaseBuf( *A0 );
      *A0 = *A1;
      goto done;
   }

/****  Interpolate */

   a0 = *A0;
   a1 = *A1;
   b = ( t - t0 ) / difference;
   a = 1 - b;
   n = veclen * F->nXYZCB;

   for( i=0; i<n; ++i, ++a0, ++a1 ) *a0 = ( a * (*a0) ) + ( b * (*a1) );
   releaseBuf( *A1 );

done:
   kOut;
   return( kStatus );
}

/*******************
********************  Calculate cell center coordinates
*******************/

static int makeCenterCoord( void )
{
   int i001, i010, i011, i101, i110, i111, iX, iY, iZ, iC, nX, nY, nZ;
   float *xyz, *xyzcc;
   kIn( makeCenterCoord );

   nX = F->nXCB;
   nY = F->nYCB;
   nZ = F->nZCB;

   /* i000 = 0; */
   /* i100 = 3; */
   i010 = 3 * ( nX + 1 );
   i001 = i010 * ( nY + 1 );
   i011 = i001 + i010;

   i101 = i001 + 3;
   i110 = i010 + 3;
   i111 = i011 + 3;

   xyz = F->xyz;
   xyzcc = F->xyzcc = (float *)malloc( 3 * F->nXYZCB * sizeof(float) );

   for( iZ=0; iZ<nZ; ++iZ ) {
      for( iY=0; iY<nY; ++iY ) {
	 for( iX=0; iX<nX; ++iX ) {
	    for( iC=0; iC<3; ++iC, ++xyz, ++xyzcc )
	       *xyzcc = 0.125 * ( xyz[0] + xyz[3] + xyz[i001] + xyz[i101]
		  + xyz[i010] + xyz[i110] + xyz[i011] + xyz[i111] );
	 }

	 xyz += 3;
      }

      xyz += i010;
   }

   kOut;
   return( kStatus );
}

/*******************
********************  Convert celltypes to nodetypes
*******************/

static int makeNodetype( int veclen, int *A0, int **A1 )
{
   int iV, iX, iY, iXX, iYY, iX1, iX2, iY1, iY2, iW1, iW2, dX, dY, dZ;
   int iZ, iZZ, iZ1, iZ2, incZ, incW;
   int  nXN, nYN, nZN, nXC, nYC, nZC, nXYC, nXYN;
   int lastXN, lastYN, lastZN, lastXC, lastYC, lastZC, rank;
   int v, *a1;
   float *FA1;
   kIn( makeNodetype );

/****  Setup */

   rank = F->rank;

   nXC = F->nXCB;
   nYC = F->nYCB;
   nZC = F->nZCB;

   lastXC = nXC - 1;
   lastYC = nYC - 1;
   lastZC = nZC - 1;

   nXN = F->nXNB;
   nYN = F->nYNB;
   nZN = F->nZNB;

   lastXN = nXN - 1;
   lastYN = nYN - 1;
   lastZN = nZN - 1;

   nXYC = nXC * nYC;
   nXYN = nXN * nYN;

   dX = veclen;
   dY = veclen * nXC;
   dZ = veclen * nXYC;

/****  Get buffer for node data */

   if( veclen == 1 )
      getReturnBuffer( lenScaNodeBuf, &FA1 );
   else
      getReturnBuffer( lenVecNodeBuf, &FA1 );

   *A1 = (int *)FA1;

/****  Repeat for each vector component */

   for( iV=0; iV<veclen; ++iV ) {
      a1 = *A1 + iV;

   /****  Outler loops over iX,iY,iZ: calculate each node
   *****     Use three loops over XYZ
   *****     Calculate parameters for inner loops
   *****  Inner loops over iXX,iYY,iZZ: look at cells around current node
   *****     Use three loops over XYZ
   *****     In most cases there are 8 cells around a node
   *****     Number reduced at sides (4 cells), edges (2) and corners (1)
   *****     Number increased for cylindrical coordinates
   *****        Use cells on both sides of cut at the zero angle
   *****        Use all cells around center axis at a radius of zero
   ****/

   /****  Outer loop over Z */

      for( iZ=0; iZ<=lastZN; ++iZ ) {

      /****  Cartesian */

	 if( ( F->coord == MRC_cartesian ) || ( rank < 3 ) ) {
	    iZ1 = ( iZ == 0 ) ? 0 : iZ - 1;
	    iZ2 = ( iZ == lastZN ) ? lastZC : iZ;
	    incZ = 1;
	 }

      /****  Cylindrical at first or last angle */

	 else if( ( iZ == 0 ) || ( iZ == lastZN ) ) {
	    iZ1 = 0;
	    iZ2 = lastZC;
	    incZ = lastZC;
	    if( incZ < 0 ) incZ = 1;
	 }

      /****  Cylindrical inside */

	 else {
	    iZ1 = iZ - 1;
	    iZ2 = iZ;
	    incZ = 1;
	 }

      /****  Outer loop over Y */

	 for( iY=0; iY<=lastYN; ++iY ) {
	    iY1 = ( iY == 0 ) ? 0 : iY - 1;
	    iY2 = ( iY == lastYN ) ? lastYC : iY;

	 /****  Outer loop over X */

	    for( iX=0; iX<=lastXN; ++iX, a1+=veclen ) {
	       iX1 = ( iX == 0 ) ? 0 : iX - 1;
	       iX2 = ( iX == lastXN ) ? lastXC : iX;

	    /****  At center axis of cylinder - sum all the way around */

	       if( ( F->coord == MRC_cylindrical ) && ( iX == 0 ) ) {
		  iW1 = 0;
		  iW2 = lastZC;
		  incW = 1;
	       }

	    /****  Otherwise sum in Z as planned - linear or in angle */

	       else {
		  iW1 = iZ1;
		  iW2 = iZ2;
		  incW = incZ;
	       }

	    /****  Inner loops: examine cells surrounding current node */

	       v = MRT_wall;

	       for( iXX=iX1; iXX<=iX2; ++iXX ) {
		  for( iYY=iY1; iYY<=iY2; ++iYY ) {
		     for( iZZ=iW1; iZZ<=iW2; iZZ+=incW ) {
			if( A0[iV+iXX*dX+iYY*dY+iZZ*dZ] != MRT_fluid ) continue;
			v = MRT_fluid;
			goto innerLoopDone;
		     }
		  }
	       }

      innerLoopDone:
	       *a1 = v;
	    }
	 }
      }
   }

   kOut;
   return( kStatus );
}

/*******************
********************  Convert cell data to node data by the average method
*******************/

static int cell2NodeByAvg( int veclen, float *A0, float **A1 )
{
   int i, iV, iX, iY, iXX, iYY, iX1, iX2, iY1, iY2, iW1, iW2, dX, dY, dZ;
   int iZ, iZZ, iZ1, iZ2, incZ, incW;
   int ng, nXN, nYN, nZN, nXC, nYC, nZC, nXYC, nXYN;
   int lastXN, lastYN, lastZN, lastXC, lastYC, lastZC, rank;
   float g, *a1;
   enum MfixReaderCellType *cellType;
   kIn( cell2NodeByAvg );

/****  Setup */

   rank = F->rank;
   cellType = F->cellType;

   nXC = F->nXCB;
   nYC = F->nYCB;
   nZC = F->nZCB;

   lastXC = nXC - 1;
   lastYC = nYC - 1;
   lastZC = nZC - 1;

   nXN = F->nXNB;
   nYN = F->nYNB;
   nZN = F->nZNB;

   lastXN = nXN - 1;
   lastYN = nYN - 1;
   lastZN = nZN - 1;

   nXYC = nXC * nYC;
   nXYN = nXN * nYN;

   dX = veclen;
   dY = veclen * nXC;
   dZ = veclen * nXYC;

/****  Get buffer for node data */

   if( veclen == 1 )
      getReturnBuffer( lenScaNodeBuf, A1 );
   else
      getReturnBuffer( lenVecNodeBuf, A1 );

/****  Repeat for each vector component */

   for( iV=0; iV<veclen; ++iV ) {
      a1 = *A1 + iV;

   /****  Outler loops: calculate each node
   *****     Use three loops over XYZ
   *****     Calculate parameters for inner loops
   *****  Inner loops: add cells around current node
   *****     Use three loops over XYZ
   *****     In most cases there are 8 cells around a node
   *****     Number reduced at sides (4), edges (2) and corners (1)
   *****     Number increased for cylindrical coordinates
   *****        Use cells on both sides of cut at the zero angle
   *****        Use all cells around center axis at a radius of zero
   ****/

      for( iZ=0; iZ<=lastZN; ++iZ ) {

      /****  Cartesian */

	 if( ( F->coord == MRC_cartesian ) || ( rank < 3 ) ) {
	    iZ1 = ( iZ == 0 ) ? 0 : iZ - 1;
	    iZ2 = ( iZ == lastZN ) ? lastZC : iZ;
	    incZ = 1;
	 }

      /****  Cylindrical at first or last angle */

	 else if( ( iZ == 0 ) || ( iZ == lastZN ) ) {
	    iZ1 = 0;
	    iZ2 = lastZC;
	    incZ = lastZC;
	    if( incZ < 1 ) incZ = 1;
	 }

      /****  Cylindrical inside */

	 else {
	    iZ1 = iZ - 1;
	    iZ2 = iZ;
	    incZ = 1;
	 }

	 for( iY=0; iY<=lastYN; ++iY ) {
	    iY1 = ( iY == 0 ) ? 0 : iY - 1;
	    iY2 = ( iY == lastYN ) ? lastYC : iY;

	    for( iX=0; iX<=lastXN; ++iX, a1+=veclen ) {
	       iX1 = ( iX == 0 ) ? 0 : iX - 1;
	       iX2 = ( iX == lastXN ) ? lastXC : iX;

	    /****  At center axis of cylinder - sum all the way around */

	       if( ( F->coord == MRC_cylindrical ) && ( iX == 0 ) ) {
		  iW1 = 0;
		  iW2 = lastZC;
		  incW = 1;
	       }

	    /****  Otherwise sum in Z as planned - linear or in angle */

	       else {
		  iW1 = iZ1;
		  iW2 = iZ2;
		  incW = incZ;
	       }

	    /****  Sum the values from the surrounding cells */

	       g = 0.0;
	       ng = 0;

	       for( iXX=iX1; iXX<=iX2; ++iXX ) {
		  for( iYY=iY1; iYY<=iY2; ++iYY ) {
		     for( iZZ=iW1; iZZ<=iW2; iZZ+=incW ) {
			i = iV + iXX*dX + iYY*dY + iZZ*dZ;
			if( cellType[i] != MRT_fluid ) continue;
			g += A0[i];
			++ng;
		     }
		  }
	       }

	       if( ng < 2 )
		  *a1 = g;
	       else
		  *a1 = g / ng;
	    }
	 }
      }
   }

   releaseBuf( A0 );
   kOut;
   return( kStatus );
}

/******************************************************************************/
/*******************                                        *******************/
/*******************            Low level input             *******************/
/*******************                                        *******************/
/******************************************************************************/

/*******************
********************  Read one record
*******************/

static int readRecord( FILE *file, char *fileName, char *desc )
{
   kIn( readRecord );
   ++currentRecord;

   if( (int)fread( record, 1, F->recl, file ) != F->recl ) {
      kErr4( failed to read %s;record %d;file %s,
	 desc, currentRecord, fileName );
   }

   if( F->reclGap > 0 ) {
      if( fseeko64( file, (off64_t)(F->reclGap), SEEK_CUR ) ) {
	 kErr3( failed to skip to end of record;record %d;file %s,
	    currentRecord, fileName );
      }
   }

   kBug3( 6, current record now %d, currentRecord );

   kOut;
   return( kStatus );
}
/*******************
********************  Read an array
*******************/
/*
   Note: Always uses fseek to go to end of last record.
*/

static int readArray(
   FILE *file,                     /* File from which array is read */
   char *fileName,                 /* Name of file used in messages */
   char *desc,                     /* Description used in messages */
   enum MfixReaderDataType type,   /* Element type */
   int nTarget,                    /* Number of elements */
   void *target                    /* Data read from file (must be alloc)  */
 )
{
   int reclUsed;
   int iRec;
   int nRec;        /* Number of records that will be read */

   long long int fileLen;   /* Number of bytes that will be read */
   long long int targetLen; /* Number of bytes that will be put into target */
   off64_t nBytes;          /* Difference between fileLen and targetLen */

   unsigned char *b;
   kIn( readArray );

/****  Make sure arguments are valid */

   if( nTarget < 1 ) {
      kErr1( invalid value for nTarget );
      goto done;
   }

   if( target == 0 ) {
      kErr1( invalid value for target );
      goto done;
   }

/****  Calculate parameters */

   switch( type ) {
      case MRD_integer:
	 targetLen = nTarget * sizeof(int);
	 nRec = ( nTarget + F->intPerRecord - 1 ) / F->intPerRecord;
	 fileLen = nRec * F->intPerRecord * sizeof(int);
	 break;
      case MRD_float:
	 targetLen = nTarget * sizeof(float);
	 nRec = ( nTarget + F->floatPerRecord - 1 ) / F->floatPerRecord;
         fileLen = nRec * F->floatPerRecord * sizeof(float);
	 break;
      case MRD_double:
	 targetLen = nTarget * sizeof(double);
	 nRec = ( nTarget + F->doublePerRecord - 1 ) / F->doublePerRecord;
         fileLen = nRec * F->doublePerRecord * sizeof(double);
	 break;
      default:
	 kErr2( type %d is invalid, type );
	 goto done;
   }

   kBug6( 7, reading %d elements\n%d records\n%lld bytes\nsaving %lld bytes,
      nTarget, nRec, fileLen, targetLen );

/****  Read data and store in target */

   if( F->reclGap > 0 ) {

   /****  Read one record at a time, skipping unused end of record */

      b = (unsigned char *)target;
      reclUsed = F->recl;

      for( iRec=0; iRec<nRec; ++iRec ) {
	 if( !readRecord( file, fileName, desc ) ) goto done;
	 if( iRec == nRec-1 ) reclUsed = targetLen - iRec * F->recl;
	 memcpy( b, record, reclUsed );
	 b += reclUsed;
      }
   }
   else {

   /****  Read entire array then skip unused end of last record */

      if( (int)fread( target, 1, targetLen, file ) != (int)targetLen ) {
	 kErr4( failed to read %s;record %d;file %s,
	    desc, currentRecord, fileName );
	 goto done;
      }

      nBytes = (off64_t)fileLen - targetLen;

      if( nBytes > 0 ) {
	 if( fseeko64( file, nBytes, SEEK_CUR ) ) {
	    kErr4( could not skip %lld bytes;just read %s;file %s,
	       nBytes, desc, fileName );
	 }
      }

      currentRecord += nRec;
   }

/****  Convert binary to that used on this computer */

   switch( type ) {
      case MRD_integer:
	 BRdInt32Array( (unsigned char *)target, (int *)target,
	    (long long)nTarget );
	 break;
      case MRD_float:
	 BRdFloat32Array( (unsigned char *)target, (float *)target,
	    (long long)nTarget );
	 break;
      case MRD_double:
	 BRdFloat64Array( (unsigned char *)target, (double *)target,
	    (long long)nTarget );
	 break;
      case MRD_string:
	 break;
   }

done:
   kOut;
   return( kStatus );
}

/*******************
********************  Skip records
*******************/

static int skipRecords( FILE *file, char *fileName, int nRecords, int whence )
{
   int iWhence;
   off64_t offset;
   static char *whenceText[4] = { "unknown", "start", "current", "end" };
   kIn( skipRecords );

   if( nRecords < 1 ) {
      kErr2( cannot skip %d records, nRecords );
      goto done;
   }

   if( whence == SEEK_SET )
      iWhence = 1;
   else if( whence == SEEK_CUR )
      iWhence = 2;
   else if( whence == SEEK_END )
      iWhence = 3;
   else
      iWhence = 0;

   offset = ( (off64_t)nRecords ) * F->reclStored;

   kBug5( 6, nRec=%d nBytes=%lld iWhence=%d (1/2/3=set/cur/end),
      nRecords, (long long int)offset, iWhence );

   if( fseeko64( file, offset, whence ) != 0 ) {
      kErr5( failed skip %d re %s from %d;file %s,
         nRecords, whenceText[ iWhence ], currentRecord, fileName );
   }
   else {
      if( whence == SEEK_SET )
	 currentRecord = nRecords;
      else if( whence == SEEK_CUR )
	 currentRecord += nRecords;
      else
	 currentRecord = 0;
   }

   kBug3( 6, current record now %d, currentRecord );

done:
   kOut;
   return( kStatus );
}

/*******************
********************  Get a string from a record
*******************/

static int getString( char *fileName, char *desc, int off, int len, char **d )
{
   static char s[513];
   kIn( getString );

   if( ( len < 1 ) || ( len > 512 ) ) {
      kErr4( requested string of length %d for %s;file %s,
         len, desc, fileName );
      goto done;
   }

   memcpy( s, record+off, len );
   trimTrailingBlanks( s, len );
   *d = strdup(s);
   kBug5( 4, set %s to %s from file %s, desc, *d, fileName );

done:
   kOut;
   return( kStatus );
}

/*******************
********************  Get a float from a record
*******************/

static int getFloat( char *fileName, char *desc, int off, float *d )
{
   kIn( getFloat );

   *d = BRdFloat32( record + off );

   if( bio_error == 1 ) {
      kErr3( cannot decode %s;file %s, desc, fileName );
   }
   else {
      kBug5( 4, set %s to %f from file %s, desc, *d, fileName );
   }

   kOut;
   return( kStatus );
}

/*******************
********************  Get an integer from a record
*******************/

static int getInt( char *fileName, char *desc, int off, int *d )
{
   kIn( getInt );

   *d = BRdInt32( record + off );

   if( bio_error == 1 ) {
      kErr3( cannot decode %s;file %s, desc, fileName );
   }
   else {
      kBug5( 4, set %s to %d from file %s, desc, *d, fileName );
   }

   kOut;
   return( kStatus );
}

/*******************
********************  Read the file header
*******************/

static int readHeader(
   FILE *file,
   char *fileName,
   int *major,
   int *minor,
   char **runName,
   int *month,
   int *day,
   int *year,
   int *hour,
   int *minute,
   int *second,
   int *nextRec,
   int *numRec
)
{
   int lenMinorStr;
   char *s, minorStr[5];
   kIn( readHeader );

/****  Read version */

   if( readRecord( file, fileName, "version" ) == 0 ) goto done;

   s = strchr( (char *)record, '=' );

   if( s == 0 ) {
      kErr2( cannot find version;file %s, fileName );
      goto done;
   }

   ++s;

   if( sscanf( s, " %d.%3s", major, minorStr ) != 2 ) {
      kErr2( cannot decode version;file %s, fileName );
      goto done;
   }

   lenMinorStr = strlen( minorStr );

   if( ( lenMinorStr < 1 ) || ( lenMinorStr > 2 ) ) {
      kErr2( minor version number is < %s > but should be 1 or 2 digits,
	 minorStr );
      goto done;
   }

   sscanf( minorStr, "%d", minor );
   if( lenMinorStr == 1 ) *minor = 10 * (*minor);

/****  Read run information */

   if( readRecord( file, fileName, "date,time" ) == 0 ) goto done;

   if( runName != 0 ) {
      getString( fileName, "run name", 0, 60, runName );
      if( kStatus == 0 ) goto done;
   }

   if( getInt( fileName, "month",  60, month  ) == 0 ) goto done;
   if( getInt( fileName, "day",    64, day    ) == 0 ) goto done;
   if( getInt( fileName, "year",   68, year   ) == 0 ) goto done;
   if( getInt( fileName, "hour",   72, hour   ) == 0 ) goto done;
   if( getInt( fileName, "minute", 76, minute ) == 0 ) goto done;
   if( getInt( fileName, "second", 80, second ) == 0 ) goto done;

/****  Read number of records and records per time step */

   if( readRecord( file, fileName, "next rec" ) == 0 ) goto done;
   if( getInt( fileName, "next record", 0, nextRec ) == 0 ) goto done;
   if( getInt( fileName, "num rec", 4, numRec ) == 0 ) goto done;

done:
   kOut;
   return( kStatus );
}

/*******************
********************  Compose name then open file
*******************/
/*
Purpose
   Open a file using the given index. Indices from 0 to 10, inclusive,
   call for a sp data output file; other indices call for the restart file.
   Try all 6 combinations of upper- and lowercase name and extension. Do
   not try various cases on the directory - only the file name itself.
*/

static int openFile( int index )
{
   int i, extL, iCase, baseL, isSpFile;
   char *c, **name, *fileStart, *fileEnd, ext[6];
   FILE **file;
   kIn( openFile );

/****  Retreive pointer to place to save information about this file */

   isSpFile = ( index >= 0 ) && ( index <= 10 );

   if( isSpFile ) {
      file = &(F->spFiles[ index ]);
      name = &(F->spFileNames[ index ]);
   }
   else {
      file = &(F->resFile);
      name = &(F->resFileName);
   }

/****  Create storage for new name */

   if( *name != 0 ) free( *name );
   baseL = strlen( F->baseFileName );
   *name = malloc( 6 + baseL );

/****  Find the part of the base file name exclusive of the directory */

   fileEnd = F->baseFileName + baseL;
   c = strrchr( F->baseFileName, '/' );

   if( c == 0 )
      fileStart = F->baseFileName;
   else
      fileStart = c + 1;

/****  Prepare file name extension */

   if( isSpFile )
      strcpy( ext, fileExt[ index ] );
   else
      strcpy( ext, ".res" );

   extL = strlen( ext );

/****  Try all case combinations */

   for( iCase=0; iCase<6; ++iCase ) {
      strcpy( *name, F->baseFileName );

   /****  First 2 cases: use file name as given */

   /****  Second 2 cases: use all lowercase file name */

      if( ( iCase >= 2 ) && ( iCase <= 3 ) ) {
	 for( c=fileStart; c<fileEnd; ++c ) *c = tolower(*c);
      }

   /****  Third 2 cases: use all uppercase file name */

      else if( ( iCase >= 4 ) && ( iCase <= 5 ) ){
	 for( c=fileStart; c<fileEnd; ++c ) *c = toupper(*c);
      }

   /****  Odd numbered cases: use all lowercase file extension */

      if( iCase%2 == 1 ) {
	 for( i=0, c=ext; i<extL; ++i, ++c ) *c = toupper(*c);
      }

   /****  Even numbered cases: use all uppercase file extension */

      else {
	 for( i=0, c=ext; i<extL; ++i, ++c ) *c = tolower(*c);
      }

   /****  Combine file name and extension */

      strcat( *name, ext );
      kBug3( 7, trying %s, *name );

   /****  Exit loop if found a file */

      *file = fopen64( *name, "r" );

      if( *file == 0 ) {
         if( kBugLevel >= 4 ) perror( "fopen64 error:" );
         kBug3( 4, failed to open %s, *name );
         continue;
      }

      kBug3( 4, opened %s, *name );
      break;
   }

/****  Cleanup needed when all attempts failed */

   if( *file == 0 ) {
      free( *name );
      *name = 0;
      kStatus = 0;
   }

   kOut;
   return( kStatus );
}

/******************************************************************************/
/*******************                                        *******************/
/*******************          Read the data files           *******************/
/*******************                                        *******************/
/******************************************************************************/

/*******************
********************  Skip to a variable
********************

Purpose
   Before the array holding a variable can be read, the file position must
   be moved past
      (1) the file header (3 records)
      (2) previous times (if any)
      (3) the time header (1 record)
      (4) previous variables (if any)
      (5) unused vector components (if not using all of a vector)
*/

static int skipToVariable(
   FILE *file,
   char *fileName,
   int t,
   int recPerTime,
   int iC,
   int varSkip,
   int veclen,        /* Vector length (1 for scalars) */
   int component,     /* Desired component for vectors (0,1,2); else ignored */
   int pC,
   int *sC,
   int pN,
   int sN
)
{
   int iP, nSkip;
   kIn( skipToVariable );

/************  Calculate nSkip as number of arrays to be skipped  *************/

   /****  Skips previous variables */

   nSkip = varSkip;

   /****  Skip previous variants of this variable */

   switch( iC ) {
      case 5:            /* Vel_s     veclen=3 */
      case 6:            /* ROP_s     veclen=1 */
      case 8:            /* T_s       veclen=1 */
      case 9:            /* X_g       veclen=1 */
      case 11:           /* THETA_m   veclen=1 */
      case 13:           /* ReactionRates  veclen=1 ... add by pnicol 2009-12-16 */

	 nSkip += pN * veclen;
	 break;
      case 10:           /* X_s; veclen=1 */
	 if( pC>0 ) {

	 /****  Loop over all previous solid phases */

	    for( iP=0; iP<pN; ++iP ) {

	    /****  Add number of solid species for this phase */

	       if( sC != 0 )
		  nSkip += sC[ iP ];
	       else
		  nSkip += 1;
	    }

	 /****  Add previous solid species at current solid phase */

	    nSkip += sN;
	 }
	 break;
   }

   /****  Skip over previous vector components */

   if( veclen > 1 ) nSkip += component;

/****  Convert number of arrays to number of records */

   nSkip *= F->recordPerFloatArray;

/***  Add in number of header records to skip */

   nSkip += 4;

/****  Add in number of records to skip from previous time steps*/

   nSkip += t * recPerTime;

/****  Skip */

   if( skipRecords( file, fileName, nSkip, SEEK_SET ) == 0 ) {
      kBug4( 1, cannot skip to time index %d in file %s, t, fileName );
   }

   kOut;
   return( kStatus );
}

/*******************
********************  Read one variable at one time
********************

Purpose
   Read one variable at one time. Keep only the part of the array that
   falls within the bounds.
Position In File
   A previous call to skipToVariable moved the file pointer to the first byte
   that is needed. If only the 2nd or 3rd component of a vector is needed, the
   unneeded components will also have been skipped before coming here.
Vectors
   The argument veclen should be 1 or 3. If 1, the variable is a scalar or
   we wish to use only one component of a vector. If 3, we wish to return all
   3 components. The file stores all the values for one component at the
   current time followed in a like manner by the other components. In this case,
   the components will be reorganized so that all 3 components are together at
   each location in the array.
*/

static int readVariable(
   FILE *file,
   char *fileName,
   int veclen,
   float **B
)
{
   int iX, iXB, iY, iYB, iZ, iZB, iVeclen;
   int x0, y0, z0, yOffset, yOffsetB, zOffset, zOffsetB;
   int nXC, nYC, nZC, nXCB, nYCB, nZCB, nXYC, nXYCB, nXYZC;
   float *b, *BB, radius, angle, gridAngle;
   kIn( readVariable );
   getInputBuf( &BB );

/****  Get storage for array after bounds are applied */

   if( veclen == 1 )
      getReturnBuffer( lenScaCellBuf, B );
   else
      getReturnBuffer( lenVecCellBuf, B );

/****  Get each part of vector */

   x0 = F->x0;
   y0 = F->y0;
   z0 = F->z0;
   nXC = F->nXC;
   nYC = F->nYC;
   nZC = F->nZC;
   nXCB = F->nXCB;
   nYCB = F->nYCB;
   nZCB = F->nZCB;
   nXYC = nXC * nYC;
   nXYCB = nXCB * nYCB;
   nXYZC = F->nXYZC;

   b = *B;

   for( iVeclen=0; iVeclen<veclen; ++iVeclen ) {

   /****  Read entire 3D array of data into temporary storage */

      readArray( file, fileName, "data", MRD_float, nXYZC, (void *)BB );
      if( kStatus == 0 ) goto done;

   /****  Copy desired part from full array into bounded array */

      for( iZB=0, iZ=z0; iZB<nZCB; ++iZB, ++iZ ) {
	 zOffsetB = iZB * nXYCB * veclen + iVeclen;
	 zOffset = iZ * nXYC;

	 for( iYB=0, iY=y0; iYB<nYCB; ++iYB, ++iY ) {
	    yOffsetB = iYB * nXCB * veclen + zOffsetB;
	    yOffset = iY * nXC + zOffset;

	    for( iXB=0, iX=x0; iXB<nXCB; ++iXB, ++iX ) {
	       b[ iXB*veclen + yOffsetB ] = BB[ iX + yOffset ];
	    }
	 }
      }
   }

/****  When returning entire vector, convert cyl to cart if needed */

      if( ( veclen == 3 ) && ( F->coord == MRC_cylindrical ) && ( F->rank > 2)){
      for( iZ=0, b=*B; iZ<nZCB; ++iZ ) {
	 gridAngle = F->gridAngle[ iZ ];

	 for( iY=0; iY<nYCB; ++iY ) {
	    for( iX=0; iX<nXCB; ++iX, b+= 3 ) {
	       radius = b[0];
	       angle = gridAngle + b[2];

	       b[0] = radius * cosf( angle );
	       b[2] = - radius * sinf( angle );
	    }
	 }
      }
   }

done:
   kOut;
   return( kStatus );
}

/*******************
********************  Select data time(s) to use
********************

Return
   Returns kStatus=0 if there is a logic error. Returns -1 for *t0A if no
   suitable time can be found.
Variables t0 and t1
   These are indices into F->timesA for most of this routine. At the end, they
   are converted into indices into F->times so that they can be used to jump
   directly to the desired time in the data file.
*/

static int selectTime(
   int nTimesA,
   float *timesA,
   int *iTimesA,
   float dataTimeReq,
   int *t0A,
   int *t1A
)
{
   int t, t0, t1, last;
   float f, diff;
   kIn( selectTime );

   last = nTimesA - 1;

/****
 ****  Special case: only one time
 ****/

   if( nTimesA < 2 ) {
      t0 = t1 = 0;
      goto wrapup;
   }

/****
 ****  Apply extrapolaton rule
 ****/

   if( dataTimeReq < timesA[0] ) {
      if( ( dataTimeReq < ( timesA[0] -  F->tolerance ) )
	    && ( F->extrapRule == MRE_never ) ) {
	 t0 = t1 = -1;
	 goto wrapup;
      }

      t0 = t1 = 0;
      goto wrapup;
   }
   else {
      if( dataTimeReq > timesA[ last ] ) {
	 if( dataTimeReq > ( timesA[ last ] +  F->tolerance ) ) {
	    if( F->extrapRule == MRE_never ) {
	       t0 = t1 = -1;
	       goto wrapup;
	    }
	 }

	 t0 = t1 = last;
	 goto wrapup;
      }
   }

/****
 ****  Find time in file nearest to requested time
 ****/

   t0 = 0;
   diff = timesA[0] - dataTimeReq;
   if( diff < 0.0 ) diff = -diff;

   for( t=1; t<nTimesA; ++t ) {
      f = timesA[t] - dataTimeReq;
      if( f < 0.0 ) f = -f;

      if( diff > f ) {
	 diff = f;
	 t0 = t;
      }
   }

/****
 ****  Found time is within tolerance of requested time
 ****/

   if( diff <= F->tolerance ) {
      t1 = t0;
      goto wrapup;
   }

/****
 ****  Found time not within tolerance
 ****/

/****  Leave request unfilled */

   if( F->interpRule == MRI_never ) {
      t0 = t1 = -1;
      goto wrapup;
   }

/****  Use it anyway */

   if( F->interpRule == MRI_nearest_any ) {
      t1 = t0;
      goto wrapup;
   }

/****  Use nearest time at or before requested time */

   if( F->interpRule == MRI_nearest_prev ) {
      if( timesA[ t0 ] <= dataTimeReq ) {
	 t1 = t0;
      }
      else {
	 t1 = t0 = t0 - 1;
	 if( t0 < 0 ) {kErr1( major logic error 1 );}
      }

      goto wrapup;
   }

/****  Use nearest time at or after requested time */

   if( F->interpRule == MRI_nearest_next ) {
      if( timesA[ t0 ] >= dataTimeReq ) {
	 t1 = t0;
      }
      else {
	 t1 = t0 = t0 + 1;
	 if( t0 >= nTimesA ) {kErr1( major logic error 2 );}
      }

      goto wrapup;
   }

/****
 ****  Nothing left to do but interpolate
 ****/

   if( F->interpRule != MRI_linear ) {
      kErr1( unknown rule );
   }
   else if( timesA[ t0 ] < dataTimeReq ) {
      t1 = t0 + 1;
      if( t1 >= nTimesA ) { kErr1( major logic error 3 ); }
   }
   else {
      t1 = t0;
      t0 = t1 - 1;
      if( t0 < 0 ) { kErr1( major logic error 4 ); }
   }

/****
 ****  Pass selected time span back to caller
 ****/

 wrapup:
   if( kStatus == 0 ) t0 = t1 = -1;
   if( t0 >= 0 ) t0 = iTimesA[ t0 ];
   if( t1 >= 0 ) t1 = iTimesA[ t1 ];
   *t0A =t0;
   *t1A = t1;
   kOut;
   return( kStatus );
}

/*******************
********************  Gather time information for one sp file
*******************/

static int readSpTimes( int iFile, int nextRec )
{
   int i, iC, nSkip, nTimes, nTimesA, *iTimesA;
   float *times, *timesA, refTime, dataTimeTolerance;
   struct MfixReaderDataChoice *c;
   kIn( readSpTimes );

   times = 0;
   timesA = 0;
   iTimesA = 0;
   dataTimeTolerance = F->dataTimeTolerance;

/****  Calculate number of time steps in this file */

   nTimes = ( nextRec - 4 ) / F->spRecPerTime[ iFile ];

   if( ( nTimes < 1 ) || ( nTimes > maxNTimes ) ) {
      kErr4( expected 1 to %d time steps but found %d;file %s,
         maxNTimes, nTimes, F->spFileNames[ iFile ] );
      goto wrapup;
   }

   kBug4( 4, %d time steps in file %s, nTimes, F->spFileNames[ iFile ] );

/****  Allocate and initialize storage for times arrays */

   times = (float *)malloc( nTimes * sizeof(float) );
   memset( times, 0, nTimes * sizeof(float) );

   timesA = (float *)malloc( nTimes * sizeof(float) );
   memset( timesA, 0, nTimes * sizeof(float) );

   iTimesA = (int *)malloc( nTimes * sizeof(int) );
   for( i=0; i<nTimes; ++i ) iTimesA[i] = -1;

/****  Retreive times from file */

   for( i=0; i<nTimes; ++i ) {

   /****  Skip to time step */

      nSkip = 3 + i * F->spRecPerTime[ iFile ];
      skipRecords( F->spFiles[ iFile ], F->spFileNames[ iFile], nSkip,SEEK_SET);

      if( kStatus == 0 ) {
	 kBug4( 1, cannot skip to time index %d in file %s,
	    i, F->spFileNames[ iFile ] );
	 goto wrapup;
      }

   /****  Read, check, store time */

      readRecord( F->spFiles[ iFile ], F->spFileNames[ iFile ], "time header" );
      if( kStatus == 0 ) goto wrapup;

      getFloat( F->spFileNames[ iFile ], "time", 0, &(times[ i ]) );
      if( kStatus == 0 ) goto wrapup;
   }

/****
 ****  Determine which times to accept
 ****    After a restart, a few times from a previous run may be repeated.
 ****    Keep the time after the restart and throw away the previous one.
 ****    Put kept times in timesA and their index in iTimesA.
*/

   i = nTimes - 1;      /* Always keep last time */
   refTime = times[i];
   iTimesA[i] = i;

   for( i=nTimes-2; i>=0; --i ) {

   /****  Keep this time when order is correct */

      if( times[i] < refTime ) {
	 iTimesA[i] = i;
	 refTime = times[i];
      }

   /****  Order wrong - abort if too wrong */

      else {
	 if( times[i] > ( refTime + dataTimeTolerance ) ) {
	    kErr4( time %f occurs before time %f; file %s,
	       times[i], refTime, F->spFileNames[ iFile ] );
	    goto wrapup;
	 }
      }
   }

/****  Build second array of times with only kept times */

   for( i=nTimesA=0; i<nTimes; ++i ) {
      if( iTimesA[i] < 0 ) continue;
      iTimesA[ nTimesA ] = iTimesA[i];
      timesA[ nTimesA ] = times[i];
      ++nTimesA;
   }

/****  Save results in structure holding information per file */

   F->spNTimes[ iFile ] = nTimes;
   F->spTimes[ iFile ] = times;

   F->spNTimesA[ iFile ] = nTimesA;
   F->spTimesA[ iFile ] = timesA;
   F->spITimesA[ iFile ] = iTimesA;

/****  Also save results in structure holding information per variable */

   for( iC=1, c=(F->choices)+1; iC<nChoices; ++iC, ++c ) {
      if( data2file[ iC ] != iFile ) continue;
      if( (iC == 10) && (!c->primaryCount) ) continue;
      c->available = 1;
      c->nTimes = nTimes;
      c->times = times;
   }

   times = 0;
   timesA = 0;
   iTimesA = 0;

 wrapup:
   if( times != 0 ) free( times );
   if( timesA != 0 ) free( timesA );
   if( iTimesA != 0 ) free( iTimesA );
   kOut;
   return( kStatus );
}

/*******************
********************  Read, check, store header from 1 sp file
*******************/

static int readSpHeader( int iFile, int *numRec, int *nextRec )
{
   int major, minor, month, day, year, hour, minute, second;
   kIn( readSpHeader );

   readHeader( F->spFiles[ iFile ], F->spFileNames[ iFile ], &major, &minor,
      0, &month, &day, &year, &hour, &minute, &second, nextRec, numRec );

   if( kStatus == 0 ) goto done;

   if( (month!=F->runDate[0]) || (day!=F->runDate[1]) || (year!=F->runDate[2])){
      kErr8( date mismatch; %d/%d/%d in restart;%d/%d/%d in %s,
         F->runDate[0], F->runDate[1], F->runDate[2], month, day, year,
         F->spFileNames[ iFile ] );
      goto done;
   }

/****  Check records per time step */

   if( *numRec != F->spRecPerTime[ iFile ] ) {
      kErr4( uses %d but expected %d records/time;file %s,
	 *numRec, F->spRecPerTime[ iFile ], F->spFileNames[ iFile ] );
      goto done;
   }

   kBug4( 4, expect %d records per time step in file %s,
      F->spRecPerTime[ iFile ], F->spFileNames[ iFile ] );

done:
   kOut;
   return( kStatus );
}

/*******************
********************  Partially read sp files
*******************/

static int scanSpFiles( void )
{
   int n, m, iFile, numRec, nextRec, openedOne;
   float t, t0, t1;
   kIn( scanSpFiles );
   openedOne = 0;

/****  Try all sp files */

   for( iFile=0; iFile<nSpFiles; ++iFile ) {

      /* The scalars file (sp9) may exist even though it has no scalars */
      if( ( iFile == 8 ) && ( F->nScalars < 1 ) ) continue;

      /* The reaction rates file (spa) may exist even though it has no rates */
      if( ( iFile == 9 ) && ( F->nReactionRates < 1 ) ) continue;

      /* The SPB file may exist even though it has no data */
      if( ( iFile == 10 ) && ( F->KEpsilon < 1 ) ) continue;

      if( openFile( iFile ) == 0 ) {
	 kStatus = 1;   /* Individual failures are not an error */

	 kBug4( 1, cannot open sp file %d with base name %s,
	    iFile, F->baseFileName );

	 continue; 
      }

      if( readSpHeader( iFile, &numRec, &nextRec ) == 0 ) goto done;
      if( readSpTimes( iFile, nextRec ) == 0 ) goto done;
      openedOne = 1;
   }

/****  Must have at least one data file */

   if( openedOne == 0 ) {
      kErr1( no sp files could be found );
      goto done;
   }

/****  Determine overall time limits */

   for( iFile = 0; iFile < nSpFiles; ++iFile ) {
      if( F->spFiles[ iFile ] == 0 ) continue;
      m = F->spNTimesA[ iFile ];

      t = (F->spTimesA[ iFile ])[ m - 1 ];
      t1 = ( iFile == 0 ) ? t : ( ( t1<t ) ? t : t1 );

      t = (F->spTimesA[ iFile ])[0];
      t0 = ( iFile == 0 ) ? t : ( ( t0>t ) ? t : t0 );

      n = ( iFile == 0 ) ? m : ( ( n<m ) ? m : n );
   }

   F->timeSpan[0] = t0;
   F->timeSpan[1] = t1;
   F->nTimes = n;

done:
   kOut;
   return( kStatus );
}

/*******************
********************  Calculate locations in sp files
*******************/

static int calcLocations( void )
{
   int i, n, iC, iS, pC, *sC, iFile, veclen, *varSkip;
   struct MfixReaderDataChoice *c;
   kIn( calcLocations );

   varSkip = &(F->varSkip[0]);

/****  Determine number of variables that preceed a given variable in a file */

   varSkip[3] = varSkip[8] = varSkip[15] = 1;
   varSkip[10] = F->nGasSpecies;

/****  Determine number of records per time step for each file */

   for( iFile=0; iFile<nSpFiles; ++iFile ) {

   /****  Sum number of variable arrays per time step in this file */

      for( iC=1, n=0, c=(F->choices)+1; iC<nChoices; ++iC, ++c ) {
	 if( data2file[ iC ] != iFile ) continue;
	 veclen = c->veclen;
	 pC = c->primaryCount;
	 sC = c->secondaryCount;

	 if( pC < 1 ) {
	 /*
	    When the primary count for most variables is 0 it means
	    that there is actually 1. X_s is a special case: when its
	    primary count is 0 there is no data.
	 */

	    if( iC != 10 ) n += veclen;
	 }
	 else if( sC == 0 ) {
	    n += veclen * pC;
	 }
	 else {
	    for( iS=0; iS<pC; ++iS ) {
	       i = sC[ iS ];
	       if( i < 1 ) i = 1;
	       n += veclen * i;
	    }
	 }
      }

      kBug4( 4, file %d has %d arrays per time step, iFile, n );

   /****  Convert from number of variable arrays to number of records */

      n *= F->recordPerFloatArray;

   /****  One record is used for the time step header */

      ++n;

   /****  Store the result */

      F->spRecPerTime[ iFile ] = n;
   }

   kOut;
   return( kStatus );
}

/******************************************************************************/
/*******************                                        *******************/
/*******************     Read and check the restart file    *******************/
/*******************                                        *******************/
/******************************************************************************/

/*******************
********************  Read K_Epsilon
*******************/

static int readKEpsilon( void )
{
   kIn( readKEpsilon );

   if( F->version[1] >= 60 ) {

   /****  Read variable */

      if( !readRecord( F->resFile, F->resFileName, "K_Epsilon"  ) ) goto done;
      if( !getInt( F->resFileName, "K_Epsilon", 0, &(F->KEpsilon) )) goto done;

   /****  Check variable */
   
   /* 
      pan fix : 29-jul-2011 
      
      Some compilers are setting .true. to -1
      
   */
    if ( F->KEpsilon != 0 ) F->KEpsilon = 1; 

      if( ( F->KEpsilon < 0 ) || ( F->KEpsilon > 1 ) ) {
	 kErr2( K+Epsilon is %d; must be 0 to 1, F->KEpsilon );
	 goto done;
      }
   }

done:
   kOut;
   return( kStatus );
}

/*******************
********************  Read number of reaction rate tables
*******************/

static int readNumberReactionRates( void )
{
   kIn( readNumberReactionRates );

   if( F->version[1] >= 50 ) {

   /****  Read variable */

      if( !readRecord( F->resFile, F->resFileName, "nRR"  ) ) goto done;
      if( !getInt( F->resFileName, "nRR", 0, &(F->nReactionRates) ) ) goto done;

   /****  Check variable */

      if( ( F->nReactionRates < 0 ) || ( F->nReactionRates > 1000 ) ) {
	 kErr2( num reaction rates is %d; must be 0 to 1000, F->nReactionRates);
	 goto done;
      }

   /****  Set number of primary counts for reaction rates variable */

      F->choices[13].primaryCount = F->nReactionRates;
   }

done:
   kOut;
   return( kStatus );
}

/*******************
********************  Read number of scalars
*******************/

static int readNumberScalars( void )
{
   int nSkip, tmpDim;
   kIn( readNumberScalars );

   if( F->version[1] >= 30 ) {

   /****  Read variables */

      if( !readRecord( F->resFile, F->resFileName, "nScalars"  ) ) goto done;
      if( !getInt( F->resFileName, "nScalars", 0, &(F->nScalars) ) ) goto done;
      if( !getInt( F->resFileName, "tmpDim", 12, &tmpDim ) ) goto done;

   /****  Check variables */

      if( ( F->nScalars < 0 ) || ( F->nScalars > 1000 ) ) {
	 kErr2( NScalar is %d; must be 0 to 1000, F->nScalars );
	 goto done;
      }

      if( ( tmpDim < 1 ) || ( tmpDim > 1000 ) ) {
	 kErr2( DIM_tmp is %d; must be 1 to 1000, tmpDim );
	 goto done;
      }

   /****  Skip past array phase4scalar */

      nSkip = ( tmpDim + F->intPerRecord - 1 ) / F->intPerRecord;
      skipRecords( F->resFile, F->resFileName, nSkip, SEEK_CUR );

   /****  Set number of primary counts for scalars variable */

      F->choices[12].primaryCount = F->nScalars;  /* scalar */
   }

done:
   kOut;
   return( kStatus );
}

/*******************
********************  Skip more arrays
*******************/
/*
** Skip the part of the restart file after the flags array and before the
** record containing the number of scalars. Along the way, read the variable
** n_spx_res which is the dimensions of other arrays in this section. None of
** these variables or arrays are used anywhere else.
*/

static int skipArrays2( void )
{
   int ver;            /* Minor version of input file */
   int nSkip;          /* Number of records to skip */
   int iPR, dPR;       /* Integers, doubles per record */
   int nSolidPhases;   /* Number of solid phases */

   /* Variables from file; only used locally */
   int tmpSpxRes, tmpDimUsr;

   /* Num records used by arrays of corresponding length */
   int nDpDimIC, nDpDimBC, nDpDimIS, nIntDimIS, nDpDimUsr;

   kIn( skipArrays2 );

/****  Setup */

   tmpDimUsr = 5;

   ver          = F->version[1];
   iPR          = F->intPerRecord;
   dPR          = F->doublePerRecord;
   nSolidPhases = F->nSolidPhases;

   nDpDimIC  = ( tmpDimIC  + dPR - 1 ) / dPR;
   nDpDimBC  = ( tmpDimBC  + dPR - 1 ) / dPR;
   nDpDimIS  = ( tmpDimIS  + dPR - 1 ) / dPR;
   nIntDimIS = ( tmpDimIS  + iPR - 1 ) / iPR;
   nDpDimUsr = ( tmpDimUsr + dPR - 1 ) / dPR;

/* Skip down to the variable n_spx_res */

   nSkip = 0;

   if( ver >= 4 ) {
      nSkip += 8 * nDpDimIS;                 /* IS_X_W etc. */
      nSkip += 6 * nIntDimIS;                /* IS_I_W etc. */
      nSkip += tmpDimIS;                     /* IS_TYPE */

      if( ver >= 7 ) {
	 nSkip += nSolidPhases * nDpDimIS;   /* IS_VEL_S */

	 if( ver >= 8 ) {
	    nSkip += 1;                      /* CYCLIC_X etc. */
	    if( ver >= 9 ) nSkip += 1;       /* TIME etc. */
	 }
      }
   }

   kBug3( 4, skipping %d records, nSkip );
   skipRecords( F->resFile, F->resFileName, nSkip, SEEK_CUR );

/****  Set the variable n_spx_res */

   if( ver >= 50 ) {
      if( !readRecord( F->resFile, F->resFileName, "n_spx_res"  ) ) goto done;
      if( !getInt( F->resFileName, "n_spx_res", 0, &tmpSpxRes ) ) goto done;
   }
   else {
      tmpSpxRes = 9;
   }

/****  Skip down to the variable nScalar */

   nSkip = 0;

   if( ver >= 9 ) {
      nSkip += tmpSpxRes;                 /* SPX_DT */
      nSkip += 1 + nSolidPhases;          /* SPECIES_EQ */
      nSkip += nDpDimUsr * 7;             /* USR_DT etc. */
      nSkip += nDpDimIC * 2;              /* IC_P_STAR etc. */
      nSkip += nDpDimBC * 6;              /* BC_DT_0 etc. */
      nSkip += tmpDimUsr + tmpDimIC;      /* USR_FORMAT, IC_TYPE */

      if( ver >= 10 ) {
	 nSkip += 1;                      /* MU_GMAX */

	 if( ver >= 11 ) {
	    nSkip += 1;                   /* V_EX etc. */

	    if( ver >= 12 ) {
	       nSkip += 2;                                     /* P_REF etc. */
	       nSkip += nDpDimBC * 4 * ( 1 + nSolidPhases );   /* BC_HW_G etc.*/

	       if( ver >= 13 ) {
		  nSkip += 1;             /* MOMENTUM_X_EQ etc. */

		  if( ver >= 14 ) {
		     nSkip += 1;          /* DETECT_STALL */

		     if( ver >= 15 ) {
			nSkip += 1;       /* K_GO etc .*/

			nSkip += nDpDimIC * 2 * ( 1 + nSolidPhases );
		                          /* IC_GAMMA_RG etc. */

			if( ver >= 20 ) nSkip += 1;   /* NORM_G etc. */
		     }
		  }
	       }
	    }
	 }
      }
   }

   kBug3( 4, skipping %d records, nSkip );
   skipRecords( F->resFile, F->resFileName, nSkip, SEEK_CUR );

done:
   kOut;
   return( kStatus );
}

/*******************
********************  Get cell type
*******************/

static int getCellType(
   int nR,
   struct MfixReaderDataRequest *R
)
{
   int *A, iC, iR;
   char *variableName;
   struct MfixReaderDataRequest *r;
   struct MfixReaderDataChoice *c;
   kIn( getCellType );

/****  Find celltype */

   for( iC=0, c=F->choices; iC<nChoices; ++iC, ++c ) {
      if( c->available == 0 ) continue;      /* Skip if not available */
      if( data2file[ iC ] == -1 ) break;     /* Stop if celltype */
   }

   if( iC >= nChoices ) goto done;
   variableName = c->variableName;

/****  Look for requests for celltype */

   for( iR=0, r=R; iR<nR; ++iR, ++r ) {
      if( strcmp( variableName, r->variableName ) ) continue;

   /****  Copy to return buffer then convert to node data if needed */

      if( r->format == MRF_node_average ) {
	 makeNodetype( 1, (int *)(F->cellType), &A );
	 if( !kStatus ) goto done;
	 r->dims = &(F->nNodesB[0]);
	 r->dataInt = A;
      }
      else {
	 r->dims = &(F->nCellsB[0]);
	 r->dataInt = (int *)(F->cellType);
      }

   /****  Fill in rest of request */

      kBug4( 6, request %d for %s being read from restart, iR, variableName );
      r->units = c->units;
      r->dataType = c->dataType;
      r->attach = c->attach;
      r->dataTime[0] = r->dataTime[1] = 0.0;
      r->veclen = 1;
      r->interpolated = 0;
   }

done:
   kOut;
   return( kStatus );
}

/*******************
********************  Read cell types
*******************/

static int readCellType( void )
{
   int nUnknownCellType;
   int iX, iY, iZ, nPlane;
   int *b, *BB;
   float *BBF;
   enum MfixReaderCellType *t;
   kIn( readCellType );

   nPlane = F->nXC * F->nYC;
   getInputBuf( &BBF);
   BB = (int *)(BBF);

/****  Read entire 3D array of cell types into a temporary array */

   readArray( F->resFile, F->resFileName, "flags", MRD_integer, F->nXYZC,
      (void *)BB );
   if( kStatus == 0 ) goto done;

/****  Copy desired part of temporary array into permanent array */

   F->cellType = (enum MfixReaderCellType *)malloc(
      F->nXYZCB * sizeof(enum MfixReaderCellType) );

   for( iZ = F->z0, t = F->cellType; iZ <= F->z1; ++iZ ) {
      b = BB + iZ * nPlane + F->y0 * F->nXC + F->x0;

      for( iY = F->y0; iY <= F->y1; ++iY, t+=F->nXCB, b+=F->nXC  ) {
	 memcpy( t, b, F->nXCB * sizeof(int) );
      }
   }

/****  Ensure that flag is one of the know types */

   t = F->cellType;
   b = (int *)t;
   nUnknownCellType = 0;

   for( iZ=0; iZ < F->nZCB; ++iZ ) {
      for( iY=0; iY < F->nYCB; ++iY ) {
	 for( iX=0; iX < F->nXCB; ++iX, ++t, ++b ) {
	    if( ( (*b) < MRT_min ) || ( (*b) > MRT_max ) ) {
	       kBug6( 4, flag %d at (%d,%d,%d) unknown; change to fluid,
	          *b, iX, iY, iZ );
	       ++nUnknownCellType;
	       *t = MRT_fluid;
	    }
	 }
      }
   }

   if( nUnknownCellType > 0 ) {
      kBug3( 1, %d unknown flags changed to fluid, nUnknownCellType );
   }

   F->choices[0].available = 1;

done:
   kOut;
   return( kStatus );
}

/*******************
********************  Skip a series of arrays
*******************/

static int skipArrays1( void )
{
   int i;
   int ver;            /* Minor version number */
   int nSkip;          /* Number of records to skip */
   int iPR, dPR;       /* Integers, doubles per record */
   int nGasSpecies;    /* Number of gas species */
   int nSolidPhases;   /* Number of solid phases */

   /* Num records used by arrays of corresponding length */
   int nDpDimIC, nIntDimIC, nDpDimBC, nIntDimBC;

   kIn( skipArrays1 );

   ver          = F->version[1];
   iPR          = F->intPerRecord;
   dPR          = F->doublePerRecord;
   nGasSpecies  = F->nGasSpecies;
   nSolidPhases = F->nSolidPhases;

   nDpDimIC  = ( tmpDimIC + dPR - 1 ) / dPR;
   nIntDimIC = ( tmpDimIC + iPR - 1 ) / iPR;
   nDpDimBC  = ( tmpDimBC + dPR - 1 ) / dPR;
   nIntDimBC = ( tmpDimBC + iPR - 1 ) / iPR;

                   nSkip = 1;                              /* D_P etc. */
   if( ver >=  4 ) nSkip += (int)( ( nGasSpecies + dPR - 1 ) / dPR );  /* MW_G*/
   if( ver >=  4 ) nSkip += nSolidPhases;                  /* MW_S */
                   nSkip += 6 * nDpDimIC;                  /* IC_X_W etc. */
                   nSkip += 6 * nIntDimIC;                 /* IC_I_W etc. */
                   nSkip += 3 * nDpDimIC;                  /* IC_EP_G etc. */
   if( ver <  15 ) nSkip += 2 * nDpDimIC;                  /* IC_T_S etc. */
   if( ver >=  4 ) nSkip += nGasSpecies * nDpDimIC;        /* IC_X_G */
                   nSkip += 3 * nDpDimIC;                  /* IC_U_G etc. */
                   nSkip += 4 * nSolidPhases * nDpDimIC;   /* IC_ROP_S etc. */
   if( ver >= 15 ) nSkip += nSolidPhases * nDpDimIC;       /* IC_T_S */

   if( ver >= 4 ) {
      for( i=0; i<nSolidPhases; ++i ) {
	 nSkip += F->nSolidSpecies[i] * nDpDimIC;   /* IC_X_S */
      }
   }

                   nSkip += 6 * nDpDimBC;                  /* BC_X_W etc. */
                   nSkip += 6 * nIntDimBC;                 /* BC_I_W etc. */
                   nSkip += 3 * nDpDimBC;                  /* BC_EP_G etc. */
   if( ver <  15 ) nSkip += 2 * nDpDimBC;                  /* BC_T_S etc. */
   if( ver >=  4 ) nSkip += nGasSpecies * nDpDimBC;        /* BC_X_G */
                   nSkip += 7 * nDpDimBC;                  /* BC_U_G etc. */
		   nSkip += 5 * nSolidPhases * nDpDimBC;   /* BC_ROP_S etc. */
   if( ver >=  4 ) nSkip += nSolidPhases * nDpDimBC;       /* BC_W_S */
   if( ver >= 15 ) nSkip += nSolidPhases * nDpDimBC;       /* BC_T_S */

   if( ver >= 4 ) {   /* BC_X_S */
      for( i=0; i<nSolidPhases; ++i ) nSkip += F->nSolidSpecies[i] *nDpDimBC;
   }

   if( ver == 0 )   /* BC_TYPE */
      nSkip += 10;
   else
      nSkip += tmpDimBC;

   kBug3( 4, skipping %d records, nSkip );
   skipRecords( F->resFile, F->resFileName, nSkip, SEEK_CUR );
   kOut;
   return( kStatus );
}

/*******************
********************  Convert coordinate storage and system
*******************/

static int convertCoordinates( float *X, float *Y, float *Z )
{
   int i, n, iX, iY, iZ, jZ, nN, nXN, nYN, nZN, nXZN, nXYN, nYZN;
   float v, *U, *W, *xyz, *XYZ;
   float sina, cosa, radius, angle, *gridAngle;
   kIn( convertCoordinates );

   U = W = 0;

   nXN = F->nXNB;   /* The size of the array of vertices AFTER clipping */
   nYN = F->nYNB;
   nZN = F->nZNB;
   nXYN = nXN * nYN;
   nXZN = nXN * nZN;
   nYZN = nYN * nZN;
   nN =  nXN * nYZN;

/************************  Allocate coordinate array  *************************/

   XYZ = F->xyz = (float *)malloc( 3 * nN * sizeof(float) );

/*****************  Already in cartesian - convert to irregular  **************/

   if( ( F->coord == MRC_cartesian ) || ( F->rank < 3 ) ) {

   /****  Set X coordinates */

      for( iX=0; iX<nXN; ++iX ) {
	 v = X[ iX ]; 
	 xyz = XYZ + iX*3;
	 for( i=0; i<nYZN; ++i, xyz+=3*nXN ) *xyz = v;
      }

   /****  Set Y coordinates */

      for( iY=0; iY<nYN; ++iY ) {
	 v = Y[ iY ];

	 for( iZ=0; iZ<nZN; ++iZ ) {
	    xyz = 1 + XYZ + iY*3*nXN + iZ*3*nXYN;
	    for( iX=0; iX<nXN; ++iX, xyz+=3 ) *xyz = v;
	 }
      }

   /****  Set Z coodinates */

      for( iZ=0; iZ<nZN; ++iZ ) {
	 v = Z[ iZ ];

	 for( iY=0; iY<nYN; ++iY ) {
	    xyz = 2 + XYZ + iY*3*nXN + iZ*3*nXYN;
	    for( iX=0; iX<nXN; ++iX, xyz+=3 ) *xyz = v;
	 }
      }
   }

/*************  In cylindrical - convert to irregular cartesian  **************/

   else if( F->coord == MRC_cylindrical ) {

   /*
   The values in X,Z are currently radius,angle
   Convert every combination into cartesian and store (x,y) in U,W
   Also save the angles from Z in gridAngle for converting vector data later
   */

      n = nZN * sizeof(float);
      gridAngle = F->gridAngle = (float *)malloc(n);
      n *= nXN;
      U = (float *)malloc(n);
      W = (float *)malloc(n);

   /****  Outer loop is over angle */

      for( iZ=jZ=0; iZ<nZN; ++iZ, jZ+=nXN, ++gridAngle ) {
	 angle = *gridAngle = Z[ iZ ];
	 sina = - sinf( angle );
	 cosa = cosf( angle );

      /****  Inner loop is over radius */

	 for( iX=0; iX<nXN; ++iX ) {
	    radius = X[ iX ];
	    U[ iX + jZ ] = radius * cosa;
	    W[ iX + jZ ] = radius * sina;
	 }
      }

   /*
   Convert angle of nodes to angle of cells where
      cell_angle = average of node_angles at sides of cell
   The last float in F->gridAngle will no longer be used
   */

      for( iZ=0, gridAngle=F->gridAngle; iZ<F->nZCB; ++iZ, ++gridAngle )
	 gridAngle[0] = 0.5 * ( gridAngle[0] + gridAngle[1] );

   /****  Set X coordinates */

      for( iZ=0; iZ<nZN; ++iZ ) {
	 for( iX=0; iX<nXN; ++iX ) {
	    v = U[ iX + iZ * nXN ];
	    xyz = XYZ + iX*3 + iZ*3*nXYN;
	    for( iY=0; iY<nYN; ++iY, xyz+=3*nXN ) *xyz = v;
	 }
      }

   /****  Set Y coordinates */

      for( iY=0; iY<nYN; ++iY ) {
	 v = Y[ iY ];

	 for( iZ=0; iZ<nZN; ++iZ ) {
	    xyz = 1 + XYZ + iY*3*nXN + iZ*3*nXYN;
	    for( iX=0; iX<nXN; ++iX, xyz+=3 ) *xyz = v;
	 }
      }

   /****  Set Z coordinates */

      for( iZ=0; iZ<nZN; ++iZ ) {
	 for( iX=0; iX<nXN; ++iX ) {
	    v = W[ iX + iZ * nXN ];
	    xyz = 2 + XYZ + iX*3 + iZ*3*nXYN;
	    for( iY=0; iY<nYN; ++iY, xyz+=3*nXN ) *xyz = v;
	 }
      }
   }

/******************  Can not handle this coordinate system  *******************/

   else {
      kErr2( cannot handle coordinate type %d, F->coord );
   }

/*******************************  Cleanup  ************************************/

   if( X != 0 ) free(X);
   if( Y != 0 ) free(Y);
   if( Z != 0 ) free(Z);
   if( U != 0 ) free(U);
   if( W != 0 ) free(W);
   kOut;
   return( kStatus );
}

/*******************
********************  Read run information
*******************/

static int readRunInfo( void )
{
   char *coord;
   kIn( readRunInfo );
   coord = 0;

   if( readRecord( F->resFile, F->resFileName, "units"  ) == 0 ) goto wrapup;

   getString( F->resFileName, "2nd run name", 0, 60, &(F->runName[1]) );
   if( kStatus == 0 ) goto wrapup;

   getString( F->resFileName, "description", 60, 60, &(F->desc) );
   if( kStatus == 0 ) goto wrapup;

   if( getString( F->resFileName, "units", 120, 16,&(F->units))==0) goto wrapup;

   getString( F->resFileName, "run type", 136, 16, &(F->runType) );
   if( kStatus == 0 ) goto wrapup;

   if( getString( F->resFileName, "coord", 152, 16, &coord ) == 0 ) goto wrapup;
   F->coord = MRC_unknown;

   if( strstr( coord, "cartesian" ) != 0 ) {
      F->coord = MRC_cartesian;
   }
   if( strstr( coord, "Cartesian" ) != 0 ) {
      F->coord = MRC_cartesian;
   }
   if( strstr( coord, "CARTESIAN" ) != 0 ) {
      F->coord = MRC_cartesian;
   }
   else if( strstr( coord, "cylindrical" ) != 0 ) {
      F->coord = MRC_cylindrical;
   }
   else if( strstr( coord, "Cylindrical" ) != 0 ) {
      F->coord = MRC_cylindrical;
   }
   else if( strstr( coord, "CYLINDRICAL" ) != 0 ) {
      F->coord = MRC_cylindrical;
   }
   else {
      kBug4( 1, coord is %s (record %d); must be cartesian or cylindrical,
         coord, currentRecord );
      F->coord = MRC_unknown;
   }

   kBug4( 4, coordinate type is %d (%s),F->coord,MfixReaderCoordText[F->coord]);

wrapup:
   if( coord != 0 ) free( coord );
   kOut;
   return( kStatus );
}

/*******************
********************  Read one array of coordinates
*******************/

static int readOneCoordinate(
   double *D,    /* Temporary storage for full array of cells/coordinates */
   char *axis,   /* Name of axis being read (X, Y or Z) */
   int nA,       /* Number of cells before filter */
   int nAB,      /* Number of cells after filter */
   int off,      /* First cell to keep after filter (0,1,2,...,n-1) */
   float **s     /* Returned as coordinates after filter */
)
{
   int i;
   float *f;
   double *d;

   kIn( readOneCoordinate );
   kBug4( 4, reading %s starting on record %d, axis, currentRecord );

/*****  Read entire array of doubles into memory

   Put the cell widths into D starting at D[1], leaving the first location
   undefined for now. See next section for reason. The width of the first
   real cell is thus in D[1] and D[2]; of the second real cell in D[3]; etc.

****/

   readArray( F->resFile, F->resFileName, axis, MRD_double, nA, D+1 );
   if( kStatus == 0 ) goto done;

/****  Shortcut: if only one cell, there is only one coordinate */

   if( nAB < 2 ) {
      *s = (float *)malloc( sizeof(float) );
      **s = 0.0;
      goto done;
   }

/****  Convert from cell widths to node coordinates.

   The first coordinate, D[0], is the negative of width of the first real
      cell width, i.e., -D[1].
   The second coordinate, D[1], is always 0 (start of the first real cell).
   The third coordinate, D[2], is the the width of the first real cell, D[2].
   The fourth coordinate, D[3], is the sum of the width of the first 2 real
      cells, D[3] + D[2].
   The fifth coordinate, F[4], is the sum of the width of the first 3 real
      cells, D[4]+D[3] (at this point, D[3] is the location of the fourth
      coordinate which is in turn the sum of the width of the first 2 real
      cells).
*/

   D[0] = -D[1];
   for( i=1, d=D+1; i<=nA; ++i, ++d ) *d += *(d-1);

/****  Convert entire array to floats */

   for( i=0, d=D, f=(float *)D; i<=nA; ++i, ++f, ++d ) *f = (float)(*d);

/****  Copy coordinates after filter to coordinate storage */

   f = ((float *)D) + off;
   *s = (float *)malloc( ( nAB + 1 ) * sizeof(float) );
   memcpy( *s, f, ( nAB + 1 ) * sizeof(float) );

done:
   kOut;
   return( kStatus );
}

/*******************
********************  Read and filter coordinates
*******************/

static int getCoordinates( float **X, float **Y, float **Z )
{
   int n;
   double *D;
   kIn( getCoordinates );

   *X = *Y = *Z = 0;

   n = F->nXC;
   if( n < F->nYC ) n = F->nYC;
   if( n < F->nZC ) n = F->nZC;
   ++n;
   D = (double *)malloc( n * sizeof(double) );

   if( !readOneCoordinate( D, "X", F->nXC, F->nXCB, F->x0, X ) ) goto error;
   if( !readOneCoordinate( D, "Y", F->nYC, F->nYCB, F->y0, Y ) ) goto error;
   if( !readOneCoordinate( D, "Z", F->nZC, F->nZCB, F->z0, Z ) ) goto error;
   goto wrapup;

error:
   if( *X != 0 ) free( *X );
   if( *Y != 0 ) free( *Y );
   if( *Z != 0 ) free( *Z );

wrapup:
   if( D != 0 ) free(D);
   kOut;
   return( kStatus );
}

/*******************
********************  Read number of species
*******************/

static int readNumberSpecies( void )
{
   int i, n;
   kIn( readNumberSpecies );

   if( readRecord( F->resFile, F->resFileName, "species"  ) == 0 ) goto done;

/****  Obtain the number of gas species */

   getInt( F->resFileName, "n gas species", 0, &(F->nGasSpecies) );
   if( kStatus == 0 ) goto done;

   F->choices[ 9].primaryCount   = F->nGasSpecies;  /* X_g */

/****  Allocate an array for the number of solid species and then obtain */

   F->nSolidSpecies = (int *)malloc( F->nSolidPhases * sizeof(int) );

   for( i=n=0; i<F->nSolidPhases; ++i ) {
      getInt( F->resFileName, "n solid species", 4 + i*4,
	 &(F->nSolidSpecies[i]) );

      if( kStatus == 0 ) goto done;
      n += F->nSolidSpecies[i];
   }

/****
   When the number of solid species is zero for all solid phases the
   variable X_s is not written even when the number of solid phases is
   above zero. For this one case set the number of solid phases to
   zero.
****/

   if( !n ) F->choices[10].primaryCount = 0;
   F->info[ 14 ].iValue = F->nSolidSpecies;
   F->choices[10].secondaryCount = F->nSolidSpecies;  /* X_s */

 done:
   kOut;
   return( kStatus );
}

/*******************
********************  Read various dimensions
********************

Note: The record used here was read in readArrayDim.
*/

static int readMiscDim( void )
{
   kIn( readMiscDim );

   if( F->version[1] > 0 ) {
      if( getInt( F->resFileName, "dim IC", 4*15, &tmpDimIC ) == 0 ) goto done;
      if( getInt( F->resFileName, "dim BC", 4*16, &tmpDimBC ) == 0 ) goto done;

      if( ( tmpDimIC < 1 ) || ( tmpDimBC < 1 ) ) {
	 kErr3( IC/BC dim are %d/%d; must greater than 0, tmpDimIC, tmpDimBC );
	 goto done;
      }
   }
   else {
      tmpDimIC = tmpDimBC = 5;
   }

   if( F->version[1] >= 4 ) {
      if( getInt( F->resFileName, "dim C", 4*17, &tmpDimC ) == 0 ) goto done;

      if( tmpDimC < 1 ) {
	 kErr2( C dimension is %d; must greater than 0, tmpDimC );
	 goto done;
      }
   }
   else {
      tmpDimC = 5;
   }

   if( F->version[1] >= 5 ) {
      if( getInt( F->resFileName, "dim IS", 4*18, &tmpDimIS ) == 0 ) goto done;

      if( tmpDimIS < 1 ) {
	 kErr2( IS dimension is %d; must greater than 0, tmpDimIS );
	 goto done;
      }
   }
   else {
      tmpDimIS = 5;
   }

done:
   kOut;
   return( kStatus );
}

/*******************
********************  Read number of solid phases
********************

Note: The record used here was read in readArrayDim.
*/

static int readNumberSolidPhases( void )
{
   kIn( readNumberSolidPhases );

   getInt( F->resFileName, "number solid phases", 4*14, &(F->nSolidPhases) );
   if( kStatus == 0 ) goto done;

   if( ( F->nSolidPhases < 1 ) || ( F->nSolidPhases > 127 ) ) {
      kErr2( num solid phases is %d; must be 1 to 127, F->nSolidPhases );
      goto done;
   }

   F->info[ 14 ].nValues = F->nSolidPhases; /* Len array w/num solid species */

   F->choices[ 5].primaryCount = F->nSolidPhases;  /* Vel_s */
   F->choices[ 6].primaryCount = F->nSolidPhases;  /* ROP_s */
   F->choices[11].primaryCount = F->nSolidPhases;  /* THETA_m */

   if( F->version[1] > 15 )
      F->choices[8].primaryCount = F->nSolidPhases;  /*T_s*/

/*
   X_s
   Note: this will be changed to zero later if the number of solid
   species is zero for all phases
*/

   F->choices[10].primaryCount = F->nSolidPhases;

done:
   kOut;
   return( kStatus );
}

/*******************
********************  Apply filter and allocate input buffer
*******************/

static int updateInternals( void )
{
   kIn( updateInternals );

   F->nXN = F->nXC + 1;
   F->nYN = F->nYC + 1;
   F->nZN = F->nZC + 1;

   F->x0 = F->skip[0];
   F->y0 = F->skip[2];
   F->z0 = F->skip[4];
   F->x1 = F->nXC - ( 1 + F->skip[1] );
   F->y1 = F->nYC - ( 1 + F->skip[3] );
   F->z1 = F->nZC - ( 1 + F->skip[5] );

   F->x0 = ( F->x0 < 0 ) ? 0 : ( ( F->x0 >= F->nXC ) ? (F->nXC)-1 : F->x0 );
   F->y0 = ( F->y0 < 0 ) ? 0 : ( ( F->y0 >= F->nYC ) ? (F->nYC)-1 : F->y0 );
   F->z0 = ( F->z0 < 0 ) ? 0 : ( ( F->z0 >= F->nZC ) ? (F->nZC)-1 : F->z0 );
   F->x1 = ( F->x1 < F->x0 ) ? F->x0 : ( (F->x1 >= F->nXC) ? (F->nXC)-1 :F->x1);
   F->y1 = ( F->y1 < F->y0 ) ? F->y0 : ( (F->y1 >= F->nYC) ? (F->nYC)-1 :F->y1);
   F->z1 = ( F->z1 < F->z0 ) ? F->z0 : ( (F->z1 >= F->nZC) ? (F->nZC)-1 :F->z1);

   F->nXCB = F->x1 - F->x0 + 1;
   F->nYCB = F->y1 - F->y0 + 1;
   F->nZCB = F->z1 - F->z0 + 1;

   F->nXNB = ( F->nXCB < 2 ) ? 1 : F->nXCB + 1;
   F->nYNB = ( F->nYCB < 2 ) ? 1 : F->nYCB + 1;
   F->nZNB = ( F->nZCB < 2 ) ? 1 : F->nZCB + 1;

   F->nCellsB[0] = F->nXCB;
   F->nCellsB[1] = F->nYCB;
   F->nCellsB[2] = F->nZCB;

   F->nNodesB[0] = F->nXNB;
   F->nNodesB[1] = F->nYNB;
   F->nNodesB[2] = F->nZNB;

   F->nXYZCB = F->nXCB * F->nYCB * F->nZCB;

   kBug5( 4, filtered X coord: %d cells from %d to %d, F->nXCB, F->x0, F->x1 );
   kBug5( 4, filtered Y coord: %d cells from %d to %d, F->nYCB, F->y0, F->y1 );
   kBug5( 4, filtered Z coord: %d cells from %d to %d, F->nZCB, F->z0, F->z1 );

   updateBuf();

   kOut;
   return( kStatus );
}

/*******************
********************  Read array dimensions
*******************/

static int readArrayDim( void )
{
   int rank;
   kIn( readArrayDim );

/****  Get record (parts will be decoded in other routines) */

   if( readRecord( F->resFile, F->resFileName, "sizes"  ) == 0 ) goto done;

/****  Decode three axis dimensions */

   if( getInt( F->resFileName, "nx", 4*9,  &(F->nXC) ) == 0 ) goto done;
   if( getInt( F->resFileName, "ny", 4*10, &(F->nYC) ) == 0 ) goto done;
   if( getInt( F->resFileName, "nz", 4*11, &(F->nZC) ) == 0 ) goto done;

/****  Check individual dimensions */

   if( ( F->nXC < 1 ) || ( F->nXC == 2 ) || ( F->nXC > nXYZMax ) ) {
      kErr3( array size in X is %d; must be 1/3/4/.../%d, F->nXC, nXYZMax );
      goto done;
   }

   if( ( F->nYC < 1 ) || ( F->nYC == 2 ) || ( F->nYC > nXYZMax ) ) {
      kErr3( array size in Y is %d; must be 1/3/4/.../%d, F->nYC, nXYZMax );
      goto done;
   }

   if( ( F->nZC < 1 ) || ( F->nZC == 2 )  || ( F->nZC > nXYZMax ) ) {
      kErr3( array size in Z is %d; must be 1/3/4/.../%d, F->nZC, nXYZMax );
      goto done;
   }

/****  Check dimensions as a group while finding rank */

   if( F->nXC == 1 ) {
      if( ( F->nYC > 1 ) || ( F->nZC > 1 ) ) {
	 kErr1( nXC==1 but nYC>1 or nZC>1 );
	 goto done;
      }

      rank = 1;
   }
   else if( F->nYC == 1 ) {
      if( F->nZC > 1 ) {
	 kErr1( nYC==1 but nZC>1 );
	 goto done;
      }

      rank = 1;
   }
   else if( F->nZC == 1 ) {
      rank = 2;
   }
   else {
      rank = 3;
   }

/****  Store information in various places and determine node dimensions */

   F->rank = rank;
   F->nCells[0] = F->nXC;
   F->nCells[1] = F->nYC;
   F->nCells[2] = F->nZC;
   F->nNodes[0] = ( F->nXC < 2 ) ? 1 : F->nXC + 1;
   F->nNodes[1] = ( F->nYC < 2 ) ? 1 : F->nYC + 1;
   F->nNodes[2] = ( F->nZC < 2 ) ? 1 : F->nZC + 1;
   F->nXYZC = F->nXC * F->nYC * F->nZC;

/****  Calculate parameter dependent on grid size */

   F->recordPerFloatArray = ( F->nXYZC +F->floatPerRecord -1)/F->floatPerRecord;

done:
   kOut;
   return( kStatus );
}

/*******************
********************  Read, echo and check header
*******************/

static int readRestartHeader( void )
{
   int nextRec, numRec;
   kIn( readRestartHeader );

   readHeader( F->resFile, F->resFileName, F->version, (F->version)+1,
      &(F->runName[0]), &(F->runDate[0]), &(F->runDate[1]), &(F->runDate[2]),
      &(F->runTime[0]), &(F->runTime[1]), &(F->runTime[2]), &nextRec, &numRec );

   if( kStatus == 0 ) goto done;

   if( ( F->version[0] < majorMin ) || ( F->version[0] > majorMax ) ) {
      kErr4( major version is %d; must be %d to %d,
         F->version[0], majorMin, majorMax );
      goto done;
   }

   if( ( F->version[1] < minorMin ) || ( F->version[1] > minorMax ) ) {
      kErr4( minor version is %d; must be %d to %d,
         F->version[1], minorMin, minorMax );
      goto done;
   }

   if( F->version[1] <= 15 ) {
      F->choices[8].primaryName = 0;
      F->choices[8].primaryCount = 0;
   }

   kBug4( 4, version is %d.%d, F->version[0], F->version[1] );

done:
   kOut;
   return( kStatus );
}

/*******************
********************  Overall read logic
*******************/

static int readRestartFile( void )
{
   int n;
   float *X, *Y, *Z;
   kIn( readRestartFile );

   if( openFile(-1) == 0 ) {
       kErr2( cannot open restart file;base name is %s, F->baseFileName );
     goto wrapup;
   }

   if( readRestartHeader() == 0 ) goto wrapup;
   if( readArrayDim() == 0 ) goto wrapup;
   if( updateInternals() == 0 ) goto wrapup;
   if( readNumberSolidPhases() == 0 ) goto wrapup;
   if( readMiscDim() == 0 ) goto wrapup;

   if( F->version[1] >= 4 ) {
      n = ( tmpDimC + F->doublePerRecord - 1 ) / F->doublePerRecord;
      n += tmpDimC;
      skipRecords( F->resFile, F->resFileName, n, SEEK_CUR ); /*Skip C,C_NAME*/
   }

   if( !readNumberSpecies() ) goto wrapup;
   if( !getCoordinates( &X, &Y, &Z ) ) goto wrapup;
   if( !readRunInfo() ) goto wrapup;
   if( !convertCoordinates( X, Y, Z ) ) goto wrapup;
   if( !skipArrays1() ) goto wrapup;
   if( !readCellType() ) goto wrapup;
   if( !skipArrays2() ) goto wrapup;
   if( !readNumberScalars() ) goto wrapup;
   if( !readNumberReactionRates() ) goto wrapup;
   if( !readKEpsilon() ) goto wrapup;

wrapup:
   if( F->resFile != 0 ) fclose( F->resFile );
   F->resFile = 0;
   kOut;
   return( kStatus );
}

/******************************************************************************/
/*******************                                        *******************/
/*******************                  Get data              *******************/
/*******************                                        *******************/
/******************************************************************************/
/*
Purpose
   Get one or more arrays of data at a single time step. The location filter
   is applied.
Request Grouping
   The program will run more efficiently if you combine requests so that
   values stored in the same file are requested at the same time. For
   example, obtain U_g,V_g,W_g in the same request rather than in three
   seperate requests.
Vectors
   If the value of "component" in the request is 0, 1 or 2 and the variable is
   a vector, only one component is returned. Otherwise, all 3 components are
   returned with the component index changing fastest.
Return
   Check the returned data pointers to determine what data was available.
   Unavalable data is not considered an error. 
*/

int MfixReaderGetData(
   int handle,
   float dataTimeReq,   /* Desired time (seconds) */
   int nRequests,
   struct MfixReaderDataRequest *requests
)
{
   int t, t0, t1, iC, iR, pC, *sC, pN, sN, iFile, varSkip;
   int veclen;       /* Length of vector that will be returned */
   int component;    /* Component of a vector that will be returned */
   int recPerTime;   /* Number of records per time step */
   float *A0, *A1;
   char *vN, *fileName;
   FILE *file;
   struct MfixReaderDataChoice *c;
   struct MfixReaderDataRequest *r;
   kIn( MfixReaderGetData );

/****  Checks */

   if( findInstance( handle ) == 0 ) goto done;

   if( F->opened == 0 ) {
      kBug2( 1, ignoring; file set is not open );
      goto done;
   }

/****  Preparations */

   releaseBuf(0);
   bio_big_endian_input = F->big_endian_input;

   for( iR=0, r=requests; iR<nRequests; ++iR, ++r ) {
      r->units = 0;
      r->dataType = 0;
      r->attach = 0;
      r->veclen = 0;
      r->interpolated = 0;
      r->dataTime[0] = r->dataTime[1] = 0.0;
      r->dims = 0;
      r->dataInt = 0;
      r->dataFloat = 0;
   }

/****  Fill requests for cell types */

   if( getCellType( nRequests, requests ) == 0 ) goto done;

/****  Consider every variable */

   for( iC=1, c=1+F->choices; iC<nChoices; ++iC, ++c ) {
      if( c->available == 0 ) continue;      /* Skip if not available */
      iFile = data2file[ iC ];
      if( iFile == -1 ) continue;   /* Skip celltype */

      vN = c->variableName;
      pC = c->primaryCount;
      sC = c->secondaryCount;

   /****  Consider every request */

      for( iR=0, r=requests; iR<nRequests; ++iR, ++r ) {
	 if( strcmp( vN, r->variableName ) != 0 ) continue;

      /****  The current variable has been requested */
      /****  Determine whether requested variant is available */

         pN = r->primaryNum;
         sN = r->secondaryNum;

      /****  Variable has no variants - ignore variant part of request */

	 if( pC <= 0 ) {
	    pN = sN = 0;
	    break;
	 }

      /****  Requested primary variant not availalbe - ignore request */

         if( pN >= pC ) continue;

      /****  Variable has no secondary variant - ignore that part of request */

	 if( sC == 0 ) {
	    sN = 0;
	    break;
	 }

      /****  Requested secondary variant not available - ignore request */

         if( sN >= sC[ pN ] ) continue;
	 break;
      }

      if( iR >= nRequests ) continue;  /* Did any request match variable? */

   /****  Find suitable time(s) */

      if( selectTime( F->spNTimesA[ iFile ], F->spTimesA[ iFile ],
	     F->spITimesA[ iFile ], dataTimeReq, &t0, &t1 ) == 0 ) goto done;
      if( t0 < 0 ) continue;

   /****  Determine vector length */

      if( c->veclen != 3 ) {
	 veclen = 1;
	 component = 0;
      }
      else {
	 component = r->component;

	 if( ( component < 0 ) || ( component > 2 ) ) {
	    veclen = 3;
	    component = 0;
	 }
	 else {
	    veclen = 1;
	 }
      }

   /****  Get data at one time (or two if interpolating)  */

      file = F->spFiles[ iFile ];
      fileName = F->spFileNames[ iFile ];
      recPerTime = F->spRecPerTime[ iFile ];
      varSkip = F->varSkip[ iC ];

      for( t = t0; t <= t1; ++t ) {

         /****  Skip to first byte of data */
	 skipToVariable( file, fileName, t, recPerTime, iC, varSkip,
	    c->veclen, component, pC, sC, pN, sN );

	 if( kStatus ==0 ) goto done;

         /****  Read; crop as needed; reorganize vectors as xyz, xyz, ... */
	 if( !readVariable( file, fileName, veclen, &A1 ) ) goto done;

	 if( t == t0 ) A0 = A1;
      }

   /****  Interpolate if needed */

      if( t0 != t1 ) interpolate( dataTimeReq, c->times[t0], c->times[t1],
	 veclen, &A0, &A1 );

   /****  Convert to node data if needed */

      if( r->format == MRF_node_average ) {
	 if( !cell2NodeByAvg( veclen, A0, &A1 ) ) goto done;
	 r->dims = &(F->nNodesB[0]);
	 r->dataFloat = A1;
      }
      else {
	 r->dims = &(F->nCellsB[0]);
	 r->dataFloat = A0;
      }

   /****  Complete the request */

      r->units = c->units;
      r->dataType = c->dataType;
      r->attach = c->attach;
      r->dataTime[0] = c->times[ t0 ];
      r->dataTime[1] = c->times[ t1 ];
      r->veclen = veclen;
      r->interpolated = ( t0 != t1 );
   }

   cleanupBuf();

done:
   kOut;
   return( kStatus );
}

/******************************************************************************/
/*******************                                        *******************/
/*******************       Get list of available data       *******************/
/*******************                                        *******************/
/******************************************************************************/
/*
Purpose
   Obtain a list of the data items that are available along with a brief
   description of each item.
Varients
   Each data item may have several "varients" which are the same data item
   for different aspects of the simulation. In MFIX, the aspects that are
   used are the species (applied to some gas and some solid data) and the
   solid phase. In some cases, both varients are used. A particular
   varient is selected in a request passed to MfixReaderGetData using an
   integer that ranges from 1 to the maximum for that varient.
Primary and Secondary Varients
   When one varient is used, it is the primary varient. When two varients
   are used for the same data item, one is called the primary and one the
   secondary. The maximum secondary varient returned by
   MfixReaderGetDataChoices is that for all primary varients. A particular
   primary varient may not support all secondary varients.
Hint
   Variables using only the primary variant can be thought of as a
   single-dimensional array of variables; those that also use the secondary
   variant can be thought of as a two-dimensional array of variables where
   some of the locations in the 2D array are empty.
*/

int MfixReaderGetDataChoices(
   int handle,
   int *nChoicesArg,
   struct MfixReaderDataChoice **choices
)
{
   kIn( MfixReaderGetDataChoices );

   if( findInstance( handle ) == 0 ) goto done;

   if( F->opened == 0 ) {
      kBug2( 1, ignoring; file set is not open );
      goto done;
   }

   *nChoicesArg = nChoices;
   *choices = &(F->choices[0]);

done:
   kOut;
   return( kStatus );
}

/******************************************************************************/
/*******************                                        *******************/
/*******************              Get coordinates           *******************/
/*******************                                        *******************/
/******************************************************************************/

int MfixReaderGetCoordinates(
   int handle,
   enum MfixReaderFormat format,
   char **units,
   int *dims[3],
   float **xyz  /* Order: (xyz),(xyz),... from min to max, X changing fastest */
)
{
   kIn( MfixReaderGetCoordinates );

   if( findInstance( handle ) == 0 ) goto done;

   if( F->opened == 0 ) {
      kBug2( 1, ignoring; file set is not open );
      goto done;
   }

   *units = coordUnitName;

   if( format == MRF_node_center ) {
      if( F->xyzcc == 0 ) makeCenterCoord();
      *dims = &(F->nCellsB[0]);
      *xyz = F->xyzcc;
   }
   else {
      *dims = &(F->nNodesB[0]);
      *xyz = F->xyz;
   }

done:
   kOut;
   return( kStatus );
}

/******************************************************************************/
/*******************                                        *******************/
/*******************             Get information            *******************/
/*******************                                        *******************/
/******************************************************************************/
/*
Purpose
   Obtain information. Note that, for any given item, either integer(s)
   or string(s) are used. Here is a partial list of keywords:
      version - major, minor version
      date - month, day, year
      dims - number of nodes along X, Y, Z
      dimsoriginal - number of nodes before the locaton filter was applied
      ntime - number of time steps
      ngasspecies - called nmax(0) in MFIX
      nsolidphases - called mmax in MFIX
      nsolidspecies - an array of length nsolidphases; called nmax() in MFIX
*/

int MfixReaderGetInfo(
   int handle,
   int *nInfoArg,                /* Number of information items */
   struct MfixReaderInfo **info  /* Array of information items */
)
{
   kIn( MfixReaderGetInfo );

   if( findInstance( handle ) == 0 ) goto done;

   if( F->opened == 0 ) {
      kBug2( 1, ignoring; file set is not open );
      goto done;
   }

   *nInfoArg = nInfo;
   *info = &(F->info[0]);

done:
   kOut;
   return( kStatus );
}

/******************************************************************************/
/*******************                                        *******************/
/*******************      Set the time selection rules      *******************/
/*******************                                        *******************/
/******************************************************************************/

int MfixReaderSetTimeRules(
   int handle,
   float tolerance,
   enum MfixReaderExtrapolation extrapRule,
   enum MfixReaderInterpolation interpRule
)
{
   kIn( MfixReaderSetTimeRules );
   if( findInstance( handle ) == 0 ) goto done;
   if( tolerance < 0.0 ) tolerance = -tolerance;
   F->tolerance = tolerance;

   if( extrapRule != MRE_nearest ) extrapRule = MRE_never;
   F->extrapRule = extrapRule;

   if( ( interpRule < MRI_nearest_any ) || ( interpRule > MRI_linear ) )
      interpRule = MRI_never;
   F->interpRule = interpRule;

done:
   kOut;
   return( kStatus );
}

/******************************************************************************/
/*******************                                        *******************/
/*******************          Set the location filter       *******************/
/*******************                                        *******************/
/******************************************************************************/
/*
Purpose
   Set the location filter. The six integers provided indicate how many
   cells to remove at the ends of each axis. The first pair of integers are
   applied to the beginning and end of the X (or radius) axis; the second
   pair to the Y (or Z) axis; and the third to the Z (or angular) axis.
   The trailing end of the axis will be adjusted, if needed, to ensure that
   at least one value is retained along that axis. For example, a filter
   of {2,2,2,2,2,2} means to skip the two cells at the ends of each axis.
   The default is {0,0,0,0,0,0}.
*/

int MfixReaderSetLocationFilter( int handle, int skip[6] )
{
   kIn( MfixReaderSetLocationFilter );

   if( findInstance( handle ) == 0 ) goto done;

   if( F->opened != 0 ) {
      kBug2( 1, ignoring; file set is open );
      goto done;
   }

   memcpy( F->skip, skip, 6 * sizeof(int) );

done:
   kOut;
   return( kStatus );
}

/******************************************************************************/
/*******************                                        *******************/
/*******************      Set the data time tolerance       *******************/
/*******************                                        *******************/
/******************************************************************************/
/*
Purpose
   Set the tolerance (seconds) used to determine if an input file is corrupt.
   If a given time in a file is a smaller number (i.e., earlier) than the
   previous time in the file by an amount greater than this tolerance, the
   file is considered corrupt and the open of the file set will fail. Defaults
   to 0.000001.
*/

int MfixReaderSetDataTimeTolerance( int handle, float tol )
{
   kIn( MfixReaderSetDataTimeTolerance );

   if( findInstance( handle ) == 0 ) goto done;

   if( F->opened != 0 ) {
      kErr1( file set is open );
      goto done;
   }

   F->dataTimeTolerance = tol;

done:
   kOut;
   return( kStatus );
}

/******************************************************************************/
/*******************                                        *******************/
/*******************       Set the input file format        *******************/
/*******************                                        *******************/
/******************************************************************************/
/*
Purpose
   Set whether the input file stores binary numbers small-endian (endian
   equals 0) or big-endian (endian equals anything other than 0). The default
   (i.e., if this routine is not used) is big-endian.
*/

int MfixReaderSetInputEndian( int handle, int endian )
{
   kIn( MfixReaderSetInputEndian );

   if( findInstance( handle ) == 0 ) goto done;

   if( F->opened != 0 ) {
      kErr1( file set is open );
      goto done;
   }

   F->big_endian_input = ( endian == 0 ) ? 0 : 1;

done:
   kOut;
   return( kStatus );
}

/******************************************************************************/
/*******************                                        *******************/
/*******************    Set the input file record length    *******************/
/*******************                                        *******************/
/******************************************************************************/
/*
Purpose
   Set the number of bytes in a record. The default (i.e., if this routine
   is not used) is 512. Only the first 512 bytes are used.
*/

int MfixReaderSetRecordLength( int handle, int length )
{
   kIn( MfixReaderSetRecordLength );

   if( findInstance( handle ) == 0 ) goto done;

   if( F->opened != 0 ) {
      kErr1( file set is open );
      goto done;
   }

   if( ( length<512 ) || ( length>1048576 ) ) {
      kErr2( record length %d is invalid; must be 512 to 1048576, length );
      goto done;
   }

   F->reclStored = length;
   F->reclGap = length - F->recl;
      
done:
   kOut;
   return( kStatus );
}

/******************************************************************************/
/*******************                                        *******************/
/*******************         Set the input file name        *******************/
/*******************                                        *******************/
/******************************************************************************/
/*
Purpose
   Set the "base" or "root" name for the file set, Any extensions in the
   name provided (such as ".xxx") will be removed, Extensions will be added
   as needed: ".res", ".sp1", etc.
*/

int MfixReaderSetFileName( int handle, char *name )
{
   char *s;
   kIn( MfixReaderSetFileName );

   if( findInstance( handle ) == 0 ) goto done;

   if( F->opened != 0 ) {
      kBug2( 1, ignoring; file set is open );
      goto done;
   }

   if( name == 0 ) {
      kBug2( 1, ignoring null file name pointer );
      goto done;
   }

   if( *name == 0 ) {
      kBug2( 1, ignoring null file name );
      goto done;
   }

   if( F->baseFileName != 0 ) free( F->baseFileName );
   F->baseFileName = (char *)strdup( name );

   s = strstr( F->baseFileName, ".res" );

   s = strrchr( F->baseFileName, '.' );
   if( s != 0 ) *s = 0;

done:
   kOut;
   return( kStatus );
}

/******************************************************************************/
/*******************                                        *******************/
/*******************              Close a file set          *******************/
/*******************                                        *******************/
/******************************************************************************/
/*
Purpose
   Close a file set. The base file name and location filter are retained.
   Frees most arrays.
*/

int MfixReaderClose( int handle )
{
   int i;
   kIn( MfixReaderClose );

   if( findInstance( handle ) == 0 ) goto done;

   if( F->opened == 0 ) {
      kBug2( 1, ignoring; file set is not open );
      goto done;
   }

   F->opened = 0;

/*****************************  Close files  **********************************/

   if( F->resFile != 0 ) fclose( F->resFile );
   F->resFile = 0;

   for( i = 0; i < nSpFiles; ++i ) {
      if( F->spFiles[i] != 0 ) fclose( F->spFiles[i] );
      F->spFiles[i] = 0;
   }

/*****************************  Free file names  ******************************/

   if( F->resFileName != 0 ) free( F->resFileName );
   F->resFileName = 0;

   for( i = 0; i < nSpFiles; ++i ) {
      if( F->spFileNames[i] != 0 ) free( F->spFileNames[i] );
      F->spFileNames[i] = 0;
   }

/**************************  Free coordinates  ********************************/

   if( F->xyz   != 0 ) free( F->xyz   );
   if( F->xyzcc != 0 ) free( F->xyzcc );

   F->xyz = F->xyzcc = 0;
   F->nXC = F->nYC = F->nZC  = F->nXCB = F->nYCB = F->nZCB = 0;
   F->nXN = F->nYN = F->nZN  = F->nXNB = F->nYNB = F->nZNB = 0;
   F->nXYZC = F->nXYZCB = 0;
   F->x0 = F->x1 = F->y0 = F->y1 = F->z0 = F->z1 = 0;
   memset( F->nCells, 0, 3*sizeof(int) );
   memset( F->nCellsB, 0, 3*sizeof(int) );
   memset( F->nNodes, 0, 3*sizeof(int) );
   memset( F->nNodesB, 0, 3*sizeof(int) );
   F->coord = MRC_unknown;

/****************************  Free data arrays  ******************************/

   if( F->cellType != 0 ) free( F->cellType );
   F->cellType = 0;

/*************************  Free time information  ****************************/

   for( i = 0; i < nSpFiles; ++i ) {
      if( F->spTimes[i] != 0 ) free( F->spTimes[i] );
      F->spTimes[i] = 0;

      if( F->spTimesA[i] != 0 ) free( F->spTimesA[i] );
      F->spTimesA[i] = 0;

      if( F->spITimesA[i] != 0 ) free( F->spITimesA[i] );
      F->spITimesA[i] = 0;
   }

   F->nTimes = 0;

/***********************  Free and reset other items  *************************/

   for( i=0; i<nChoices; ++i ) {
      F->choices[i].available = 0;
      F->choices[i].nTimes = 0;
      F->choices[i].times = 0;
   }

   if( F->runName[0] != 0 ) free( F->runName[0] );
   if( F->runName[1] != 0 ) free( F->runName[1] );
   if( F->runType != 0 ) free( F->runType );
   if( F->desc != 0 ) free( F->desc );
   if( F->units != 0 ) free( F->units );
   F->runName[0] =  F->runName[1] = F->desc = F->units = 0;

   memset( F->version, 0, 2*sizeof(int) );
   memset( F->runDate, 0, 3*sizeof(int) );
   memset( F->runTime, 0, 3*sizeof(int) );

   if( F->nSolidSpecies != 0 ) free( F->nSolidSpecies );
   F->nSolidSpecies = 0;
   F->nGasSpecies = F->nSolidPhases = F->nScalars = F->nReactionRates = 0;
   F->KEpsilon = 0;

   updateBuf();

done:
   kOut;
   return( kStatus );
}

/******************************************************************************/
/*******************                                        *******************/
/*******************             Open a file set            *******************/
/*******************                                        *******************/
/******************************************************************************/
/*
Purpose
   Open a file set. The filename must be set. Optional settings such as
   the location filter, record length and file protocol must be set before
   this routine is called if they are to have an effect.
*/

int MfixReaderOpen( int handle )
{
   kIn( MfixReaderOpen );

/****  Ensure that this operation is valid */

   if( findInstance( handle ) == 0 ) goto done;

   if( F->opened != 0 ) {
      kBug2( 1, ignoring; file set is open );
      goto done;
   }

   if( F->baseFileName == 0 ) {
      kErr1( no file name );
      goto done;
   }

/****  Lock in settings */

   bio_big_endian_input = F->big_endian_input;
   F->intPerRecord      = F->recl / 4;
   F->floatPerRecord    = F->recl / 4;
   F->doublePerRecord   = F->recl / 8;

   kBug3( 4, %d=big_endian flag, F->big_endian_input );
   kBug6( 5, %d=recl; %d/%d/%d int/float/double per rec,
      F->recl, F->intPerRecord, F->floatPerRecord, F->doublePerRecord );
   kBug3( 5, %d=reclStored, F->reclStored );

/****  Initialize operations with records  */

   currentRecord = 0;

   if( recl < F->recl ) {
      if( record != 0 ) free( record );
      recl = F->recl;
      record = malloc( recl + 2 );
   }

/****  Read restart file and scan other files for times */

   if( readRestartFile() == 0 ) goto done;
   if( calcLocations() == 0 ) goto done;
   if( scanSpFiles() == 0 ) goto done;

   F->opened = 1;

done:
   kOut;
   return( kStatus );
}

/******************************************************************************/
/*******************                                        *******************/
/*******************        Destroy a reader instance       *******************/
/*******************                                        *******************/
/******************************************************************************/
/*
Purpose
   Destroy an instance of a file set reader. Calls MfixReaderClose if needed.
*/

int MfixReaderDestroy( int *handle )
{
   int i, j;
   kIn( MfixReaderDestroy );

   if( findInstance( *handle ) == 0 ) goto done;
   if( F->opened != 0 ) MfixReaderClose( *handle );
   if( F->baseFileName != 0 ) free( F->baseFileName );

   for( i=0; i<nFileSets; ++i ) {
      if( fileSets[i]->handle != *handle ) continue;
      if( i >= (nFileSets-1) ) break;
      for( j=i; j<(nFileSets-1); ++j ) fileSets[j] = fileSets[ j + 1 ];
      break;
   }

   --nFileSets;
   free( F );
   *handle = 0;

done:
   kOut;
   return( kStatus );
}

/******************************************************************************/
/*******************                                        *******************/
/*******************   Create a new file reader instance    *******************/
/*******************                                        *******************/
/******************************************************************************/
/*
Purpose
   Create an instance of a file set reader.
*/

int MfixReaderCreate( int *handle )
{
   FileSet **oldFileSets;
   kIn( MfixReaderCreate );

   if( MfixReaderId ) kStatus = kStatus;

/****  Initialize binary I/O package if needed */

   if( !bioInited ) {
      InitBIO();

      if( bio_error == 1 ) {
	 kErr1( the binary storage on this computer is not supported );
	 kStatus = 1;
	 goto done;
      }

      bioInited = 1;
   }

/****  Make room in list of instances */

   oldFileSets = fileSets;

   fileSets = (FileSet **)malloc( ( nFileSets + 1 )
      * sizeof(FileSet *) );

   if( nFileSets>0 ) {
      memcpy( fileSets, oldFileSets, nFileSets * sizeof(FileSet *) );
      free( oldFileSets );
   }

/****  Create new instance and add to list */

   F = (FileSet  *)malloc( sizeof(FileSet) );
   memset( F, 0, sizeof(FileSet) );
   fileSets[ nFileSets ] = F;

/****  Initialize new instance */

   F->handle = nextHandle;
   F->big_endian_input = 1;
   F->coord = MRC_unknown;
   F->recl = F->reclStored = reclDefault;
   F->reclGap = 0;
   F->skip[0] = F->skip[1] = F->skip[2] = F->skip[3] = F->skip[4] =F->skip[5]=1;
   F->tolerance = 0.000001;
   F->dataTimeTolerance = 0.000001;
   F->interpRule = MRI_linear;
   F->extrapRule = MRE_nearest;

   memcpy( &(F->info), &(infoTemplate), sizeof(infoTemplate) );
   F->info[0].iValue = &(F->nNodesB[0]);
   F->info[1].iValue = &(F->nNodes[0]);
   F->info[2].iValue = &(F->nTimes);
   F->info[3].fValue = &(F->timeSpan[0]);
   F->info[4].iValue = &(F->runDate[0]);
   F->info[5].iValue = &(F->runTime[0]);
   F->info[6].iValue = &(F->version[0]);
   F->info[7].cValue = &(F->runName[0]);
   F->info[8].cValue = &(F->runType);
   F->info[9].cValue = &(F->desc);
   F->info[10].cValue = &(F->units);
   F->info[11].iValue = (int *)( &(F->coord) );
   F->info[12].iValue = &(F->nGasSpecies);
   F->info[13].iValue = &(F->nSolidPhases);
   F->info[15].iValue = &(F->rank);
   F->info[16].iValue = &(F->nScalars);
   F->info[17].iValue = &(F->nReactionRates);
   F->info[18].iValue = &(F->KEpsilon);

   memcpy( &(F->choices), &(choicesTemplate), sizeof(choicesTemplate) );

/****  Wrapup */

   *handle = nextHandle;
   kBug3( 1, allocating instance with handle %d, *handle );
   ++nextHandle;
   ++nFileSets;

done:
   kOut;
   return( kStatus );
}

/***************************  End of MfixReader.c  ****************************/

