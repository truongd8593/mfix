/******************************************************************************/
/*******************                                        *******************/
/*******************          Translate MFIX files          *******************/
/*******************                   to                   *******************/
/*******************             Ensight 6 Files            *******************/
/*******************                                        *******************/
/******************************************************************************/
/*
Purpose
   Read a set of files in the MFIX format and then write the data
   into the Ensight G old format. See showHelp for options.
Original
   Author: Kent Eschenberg eschenbe@psc.edu
   Organization: Pittsburgh Supercomputing Center, CMU
   Date: 2002/1/29 10:41
Current
   $Source$
   $Author$
   $Revision$
   $Date$
Routines
   Miscellaneous
      showHelp .................. display help including command line options
      showInfo .................. show file information
      write80 ................... write a string of exactly 80 characters
      getVariableInfo ........... get information about a variable
   Manage Output Files
      makeFileName .............. make a variable's file name
      makeFileNameVar ........... make a var's file name given the var's name
      makeVarName ............... make a variable name
      openFile .................. open a file with checks
      openCalcFile .............. open a file for a calculated variable
      prepareNamesVariants ...... prepare to put variant indices in names
      prepareNames .............. prepare output names
   Write Information
      writeVariable ............. write one variable
      getVariable ............... get data for one variable
      writeOneVariant ........... write one variant at one time
      writeData ................. write all selected variables
      writeGeoDrawing ........... write one drawing to the geo file
      writeGeoPseudoGrid ........ write the geometry for the pseudo grid
      writeGeoGrid .............. write the geometry for the grid
      writeGeo .................. write the geometry file
      writeCase ................. write the case file
   Prepare information
      prepareCoordinates ........ prepare irregular coordinates
      prepareCelltypes .......... prepare cell types
      addTime ................... add one time to the list
      selectTimes ............... select times for data
      selectVariables ........... select variables
      addTube ................... add a void tube to the blanking
      prepareTest ............... prepare to create a test file
   Initialize
      getInfo ................... get file information
      checkCommandLine .......... check options from command line
      addVariable ............... add one variable to conversion list
      readCommandLine ........... parse the command line
      main ...................... main program
Output Files
   The names of the output files are constructed using part of the input
   file name and other information. The input file name, with the directory
   (if any) and extension (if any) removed is called <base>. The output
   files are always written to the directory that was current when the
   program is run. Output file names are:

      <base>.cas ............................. the Ensight case file
      <base>.geo ............................. the Ensight geometry file
      <base>_<var>_<p**><s**>.<tuple> ........ data files

   In all cases, the ** are replace with sufficient leading zeroes so that
   the number of characters are the same for all files for a given variable,
   or for all times for all variables.

      <var> ..... variable name such as EP_g
      <p**> ..... primary variant numbered from 0; not used if not needed
      <s**> ..... secondary variant numbered from 0; not used if not needed
      <tuple> ... sca/vec for scalars/vectors (all vectors have 3 elements)

   Letters used for "p" and "s" in the file names are the first letter of
   the variant name.
Vector Data in Cylindirical Coordinates
   The radius/Z/angle are in the cell's local coordinate system. This is
   converted to a global cartesian direction using the angle of the cell.
   That angle is estimated by taking the average angle of the nodes at the
   cell corners.
Test File Type 1
   A grid of 3x3x3 cells for 3 time steps is created. EP_g consists of all 0s,
   all 1s, and all 2s at the 3 times except that the center value is always 1.
   Vel_g consists of vectors all pointing towards +X, +Y or +Z at the 3 times.
   The input file name must be provided but is only used to generate the names
   of the output files.
Test File type 2
   A grid of 20x4x12 (xyz) cells for 3 time steps is created. EP_g is all
   0s for cells y=0 and y=1 then all 1s. Vel_g consists of vectors all
   pointing towards _X, _Y or _Z at the 3 times. The input file name must
   be provided but is only used to generate the names of the output files.
   If the -bound flag is also on, two "outside" regions or voids are
   created. Each is a tube of axa cells parallel to the Y axis. The two
   tubes are separated by b cells. a and b are taken from the command line
   after the -test2 flag and are forced into the range (1,4).
Pseudo Grid
   Cell data can be converted to node data by created a new grid, called the
   pseudo grid, where each vertex is at the center of a cell. This is called
   the "cell center" method for converted cell data to node data.
*/

#define kAllocate 1
#include "kLib.h"
#include "MfixReader.h"
#include "m2e.h"
#include "cell2geom.h"

static char *M2E_Id
   = "@(#) $Id$";

#define DEG2RAD ( 3.1415927 / 180.0 )

static int m2eVersionMajor = 3;
static int m2eVersionMinor = 7;

/************************  Command line options  ******************************/

/****
*****  Debug level; saved in kBugLevel
****/

static char *argBugLevel = "-bug";

/****
*****  Request to generate a test file
****/
/*
static char *argTest1 = "-test1";
static char *argTest2 = "-test2";
*/

static int testCellSize = 1;
static int testCellSpacing = 1;
static int testRequested = 0;

/****
*****  Request to only show information about the file
****/

static char *argInfo = "-i";
static int infoRequested = 0;

/****
*****  Input file endian form (0=little, 1=big)
****/

static char *argEndian = "-endian";
static int endian = 1;

/****
*****  Input file record length
****/

static char *argRecl = "-recl";
static int recl = 512;

/****
*****  Time reversal tolerance
****/

static char *argTimeTolerance = "-tol";
static float timeTolerance = 0.000001;

/****
*****  Time threshold
****/

static char *argTimeThreshold = "-tt";
static float timeThreshold = 0.001;

/****
*****  Time extrapolation method
****/

static char *argExtrapRule = "-te";
static enum MfixReaderExtrapolation extrapRule = MRE_nearest;

/****
*****  Time interpolation method
****/

static char *argInterpRule = "-tp";
static enum MfixReaderInterpolation interpRule = MRI_linear;

/****
*****  Selection of method for selecting times
****/

enum TimeMethod { TM_ALL, TM_ALL_BOUNDED, TM_STEP };
static char *timeMethodText[3] = { "all", "all bounded", "stepped" };
static char *argTimeMethod = "-tm";
static enum TimeMethod timeMethod = TM_ALL;

/****
*****  Selection of time span
****/

static char *argTimeBounds[2] = { "-t0", "-t1" };
static float timeBounds[2] = { 0.0, 1000000.0 };

/****
*****  Selection of time increment
****/

static char *argTimeIncrement = "-ti";
static float timeIncrement = 0.001;

/****
*****  Cells to skip
****/

static char *argSkip[6] = { "-x0", "-x1", "-y0", "-y1", "-z0", "-z1" };
static int skip[6] = { 0, 0, 0, 0, 0, 0 };

/****
*****  Variable format
****/

static char *argNodeAverage = "-na";
static char *argNodeCenter = "-nc";
static enum MfixReaderFormat varFormat = MRF_cell;

/****
*****  Overwrite of existing files
****/

static char *argOverwrite = "-o";
static int overwrite = 0;

/****
*****  Request for help
****/

static char *argHelp[4] = { "-h", "-H", "-help", "-HELP" };

/****
*****  Input file name (the extension is ignored so can be any in MFIX set)
****/

static char *inputFileName = 0;
static char *dirName = 0;
static char *baseName = 0;

/****
*****  Variables that will be written and their format
****/

struct VarRequestStr {
   char *name;
   enum MfixReaderFormat format;
};

typedef struct VarRequestStr VarRequest;

static char *argVariable = "-v";
static int nVarRequest = 0;       /* Number of variables */
static VarRequest *varRequest;    /* Variables */

static int usingPseudoGrid = 0;

/***************  Cell types and information derived from them  ***************/

static int needNodeFlags = 0;  /* Node flags needed */

static int *cellFlags = 0;             /* Cell flags, one per cell */
static int *pseudoGridCellFlags = 0;   /*  " for pseudo grid */
static int *nodeFlags = 0;             /* Node flags, one per vertex */

static int nCellsInside = 0;    /* Number of cells inside */
static int nPseudoInside = 0;   /* Number of pseudo-cells inside */

struct DrawTypesStr {
   char *name;   /* Short name for this type */
   int t0;       /* First integer representing this type, inclusive */
   int t1;       /* Last, inclusive */
};

static int nDrawTypes = 5;

static DrawType drawTypes[5] = {
   { "obstruction",     DR_ExclWallCorner, 100, 100 },
   { "velocity_inlet",  DR_None,           10,  10  },
   { "mass_inlet",      DR_None,           20,  20  },
   { "velocity_outlet", DR_None,           11,  11  },
   { "mass_outlet",     DR_None,           21,  21  }
};

static Drawing *drawings = 0;

/***************************  Part indices  ***********************************/

static int nextPartIndex = 1;
static int gridPartIndex = 0;
static int pseudoGridPartIndex = 0;

/*********************  Variable and file name generation  ********************/

static char *caseFileExtension = ".cas";
static char *geoFileExtension = ".geo";
static char *dataFileExtensions[2] = { ".sca", ".vec" };

static char *casePathName = 0;   /* Name of case output file incl directory */
static char *caseFileName = 0;   /*  " w/out directory */
static char *geoPathName = 0;    /* Name of geometry output file incl dir */
static char *geoFileName = 0;    /*  " w/out directory */

struct HelperStr {
   char *name;   /* Name of variable */
   int veclen;   /* 0=not used 1=scalar 2=vector */

   enum MfixReaderFormat format[3];  /* Desired formats for this variable */

   int useTime;        /* Variable changes over time if true */
   int nP;             /* Larger of 1, num 1st variant */
   int *nS;            /* Larger of 1, num 2nd variant per 1st */
   int nSmax;          /* Largest value in nS; always at least 1 */
   int nArgs;          /* Number of arguments used by fileFormat */
   char *fileFormat;   /* Format used to generate file name */
   char *nameFormat;   /* Format used to generate variable name */

/*
   Successfully wrote variable when True. Used as if [s][p] where p is
   the primary index (0,1,2,...,nP-1) and s is the secondard index
   (0,1,2,...,nS[iP]-1). The size is [nSmax][nP]. In some cases not all
   of the secondary indices are used.
*/
   int *success;
};

typedef struct HelperStr Helper;
static Helper *helper = 0;   /* 1 for each var in "C" and 1 for each calc var */
static int nHelper = 0;      /* Length of helper */

/******************************************************************************/
/**************                                                ****************/
/**************                  Miscellaneous                 ****************/
/**************                                                ****************/
/******************************************************************************/

/*******************
********************  Show help
*******************/

static int showHelp( void )
{
   kIn( showHelp );

   printf( "MFIX to Ensight Translator %d.%d\n",
      m2eVersionMajor, m2eVersionMinor );
   printf( "   This version supports input files of more than 2GB\n" );
   printf( "   Defaults shown after \"/\"\n" );
   printf( "   See also www.psc.edu/~eschenbe\n\n" );
   printf( "   Usage: m2e [options] file\n" );
   printf( "   file ......... the name of a file in the MFIX data set\n" );
   printf( "   -bug a ....... debug level/0 (0=none, 10=most)\n" );
   printf( "   -i ........... only show information/don't\n" );
   printf( "   -o ........... overwrite existing files/don't\n" );
   printf( "   -endian a .... input file byte order/1 (0=little,1=big)\n" );
   printf( "   -recl a ...... input file record length/512\n" );
   printf( "   -tol a ....... time reversal tolerance/0.000001\n" );
   printf( "   -tt a ........ threshold for equating times/0.001\n" );
   printf( "   -tm a ........ method for setting times/0\n" );
   printf( "                    0 - use all times\n" );
   printf( "                    1 - use all times within t0,t1\n" );
   printf( "                    2 - calculated using t0,t1,ti\n" );
   printf( "   -t0 a ........ first time/0\n" );
   printf( "   -t1 a ........ last time/1000000\n" );
   printf( "   -te a ........ time extrapolation mode/1\n" );
   printf( "                    0 - never\n" );
   printf( "                    1 - nearest\n" );
   printf( "   -ti a ........ time increment/0.001\n" );
   printf( "   -tp a ........ time interpolation mode/4\n" );
   printf( "                    0 - never\n" );
   printf( "                    1 - use the nearest time\n" );
   printf( "                    2 - use the nearest previous time\n" );
   printf( "                    3 - use the nearest following time\n" );
   printf( "                    4 - use linear interpolation\n" );
   printf( "   -x0 a ........ cells skipped at start of X axis/0\n" );
   printf( "   -x1 a ........ cells skipped at end of X axis/0\n" );
   printf( "   -y0 a ........ cells skipped at start of Y axis/0\n" );
   printf( "   -y1 a ........ cells skipped at end of Y axis/0\n" );
   printf( "   -z0 a ........ cells skipped at start of Z axis/0\n" );
   printf( "   -z1 a ........ cells skipped at end of Z axis/0\n");
   printf( "   -na .......... use average method to generate node data\n" );
   printf( "   -nc .......... use cell-center method to generate node data\n" );
   printf( "   -v a ......... variable to convert/all\n" );
   printf( "      celltype, EP_g, P_g, P_star, Vel_g, Vel_s\n" );
   printf( "      ROP_s, T_g, T_s, X_g, X_s, THETA_m, scalar\n" );
   printf( "      r-rate, K_Turb_g, E_Turb_G\n" );
   printf( "      plus any custom variables\n\n" );
/*
   printf( "   -test1 ....... generate test data (cube)\n" );
   printf( "   -test2 a b ... generate test data w/void size a " );
   printf(                    "and spacing b (voids)\n" );
*/

   kStatus = 0;  /* This causes the program to exit */
   kOut;
   return( kStatus );
}

/*******************
********************  Show file information
*******************/

static int showInfo( char *prefix, FILE *file )
{
   int i, j, n, nValues;
   char b[81];
   struct MfixReaderInfo *II;
   kIn( showInfo );

   for( i=0, II=info; i<nInfo; ++i, ++II ) {
      fprintf( file, "%skeyword %s ", prefix, II->keyword );
      nValues = II->nValues;

      switch( II->dataType ) {
	 case MRD_integer:
            fprintf( file, "(int): " );

	    for( j=0; j<nValues; ++j ) {
	       fprintf( file, "%d", (II->iValue)[j] );
	       if( j < (nValues-1) ) fprintf( file, "," );
	    }

	    fprintf( file, "\n" );
            break;

	 case MRD_float:
            fprintf( file, "(float): " );

	    for( j=0; j<nValues; ++j ) {
	       fprintf( file, "%f", (II->fValue)[j] );
	       if( j < (nValues-1) ) fprintf( file, "," );
	    }

	    fprintf( file, "\n" );
            break;

	 default:
            fprintf( file, "(string): " );

	    for( j=0; j<nValues; ++j ) {
	       n = strlen( (II->cValue)[j] );
	       if( n>32 ) n = 32;
	       strncpy( b, (II->cValue)[j], n );
	       b[n] = 0;

	       fprintf( file, "%s", b );
	       if( j < (nValues-1) ) fprintf( file, "," );
	    }

	    fprintf( file, "\n" );
            break;
      }
   }

   kOut;
   return( kStatus );
}

/*******************
********************  Write a string of exactly 80 characters
*******************/

static int write80( FILE *file, char *msg )
{
   int i, n;
   char buf[81];
   kIn( write80 );

   n = strlen( msg );
   if( n > 80 ) n = 80;
   strncpy( buf, msg, n );

   if( n < 80 ) {
      for( i=n; i<80; ++i ) buf[i] = ' ';
   }

   buf[80] = 0;
   fprintf( file, "%s", buf );

   kOut;
   return( kStatus );
}

/****
*****  Return information about a variable
****/

int getVariableInfo(
   char *varName,
   enum MfixReaderFormat format,
   int *available,
   int *nX,
   int *nY,
   int *nZ
)
{
   int iC;
   struct MfixReaderDataChoice *c;
   kIn( getVariableInfo );
   *available = 0;

   for( iC=0, c=choices; iC<nChoices; ++iC, ++c ) {
      if( strcmp( varName, c->variableName ) ) continue;
      if( !c->available ) break;

   /****  Found it - size depends upon format */

      *available = 1;

      switch( (int)format ) {
	 case MRF_node_average:
	    *nX = nXN;
	    *nY = nYN;
	    *nZ = nZN;
	    break;
	 case MRF_cell:
	    *nX = nXC;
	    *nY = nYC;
	    *nZ = nZC;
	    break;
	 case MRF_node_center:
	    *nX = nXPN;
	    *nY = nYPN;
	    *nZ = nZPN;
	    break;
      }
      
      break;
   }

   kOut;
   return( kStatus );
}

/******************************************************************************/
/*******************                                        *******************/
/*******************         Manage output files            *******************/
/*******************                                        *******************/
/******************************************************************************/

/*******************
********************  Make a variable's file name
*******************/

static int makeFileName( int iH, int iP, int iS, enum MfixReaderFormat f,
   int veclen, char *name )
{
   char *s, *nF, *ext;
   kIn( makeFileName );

   switch( (int)f ) {
      case MRF_cell: s = ""; break;
      case MRF_node_average: s = "_na"; break;
      case MRF_node_center: s = "_nc"; break;
   }

   nF = helper[ iH ].fileFormat;

   switch( helper[ iH ].nArgs ) {
      case 1: sprintf( name, nF, s );  break;
      case 2: sprintf( name, nF, iP, s ); break;
      case 3: sprintf( name, nF, iP, iS, s ); break;
   }

   ext = dataFileExtensions[ (veclen > 1) ];
   strcat( name, ext );

   kOut;
   return( kStatus );
}

/*******************
********************  Make a variable's file name given the variable's name
*******************/

static int makeFileNameVar( char *varName, int iP, int iS,
   enum MfixReaderFormat f, int veclen, char *fileName )
{
   int iH;
   Helper *h;
   kIn( makeFileNameVar );

   for( iH=0, h=helper; iH<nHelper; ++iH, ++h ) {
      if( !h->name ) continue;
      if( !strcmp( h->name, varName ) ) break;
   }

   if( iH >= nHelper ) {
      kErr2( variable %s is unknown, varName );
      goto done;
   }

   makeFileName( iH, iP, iS, f, veclen, fileName );

done:
   kOut;
   return( kStatus );
}

/*******************
********************  Make a variable name
*******************/

static int makeVarName( int iH, int iP, int iS, enum MfixReaderFormat f,
   char *name)
{
   char *s, *nF;
   kIn( makeVarName );

   switch( (int)f ) {
      case MRF_cell: s = ""; break;
      case MRF_node_average: s = "_na"; break;
      case MRF_node_center: s = "_nc"; break;
   }

   nF = helper[ iH ].nameFormat;

   switch( helper[ iH ].nArgs ) {
      case 1: sprintf( name, nF, s );  break;
      case 2: sprintf( name, nF, iP, s ); break;
      case 3: sprintf( name, nF, iP, iS, s ); break;
   }

   kOut;
   return( kStatus );
}

/*******************
********************  Open a file with checks
*******************/

static int openFile( char *fileName, FILE **file2 )
{
   FILE *file;
   kIn( openFile );

   if( overwrite == 0 ) {
      file = fopen( fileName, "r" );

      if( file != 0 ) {
	 kErr2( file %s already exists, fileName );
	 fclose( file );
	 file = 0;
	 goto wrapup;
      }
   }

   file = fopen( fileName, "w" );

   if( file == 0 ) {
      kErr2( cannot create %s, fileName );
      goto wrapup;
   }

   kBug3( 1, created output file %s, fileName );

wrapup:
   *file2 = file;
   kOut;
   return( kStatus );
}

/*******************
********************  Open a file for a calculated variable
*******************/

int openCalcFile(
   char *varName,
   enum MfixReaderFormat format,
   int veclen,
   FILE **file
)
{
   char fileName[256];
   kIn( openCalcFile );

   makeFileNameVar( varName, 0, 0, format, veclen, fileName );
   openFile( fileName, file );

   kOut;
   return( kStatus );
}

/*******************
********************  Prepare to put variant indices in names
*******************/

/*
   Calculate information used to build variable names using the
   variant indices. Note that indices start at 0 thus the largest
   that will be shown is (n-1), not n.
File Formats
   The format for the primary and secondary variant indices will be
   used in "sprintf" to generate the actual file or variable name. If a
   variant is not used, its format should be the null string. If it is
   used, it must generate the same number of digits for all cases, i.e.,
   ..., 08, 09, 10, 11, ... instead of ..., 8, 9, 10, 11, ... .
Loop Counts
   These are the same as that provided by the file reader except that
   the values for each will never be less than one. Thus, these values
   can be used to run loops over variants while being sure that the
   contents of the loop are executed at least once.
*/

static int prepareNamesVariants(
   struct MfixReaderDataChoice *c,
   char *pF,     /* Returns primary variant format */
   char *sF,     /*  " secondary */
   int *nP,      /* Returns primary variant loop count */
   int **nS,     /* Returns array of secondary variant loop counts */
   int *nSmax,   /* Returns maximum value in **nS (always at least 1) */
   int **success /* Returns array of success flags */
)
{
   int n;
   int iP;         /* Loop index for primary variant */
   int pC, *sC;    /* Primary, secondary variant count from file */

/*
   Maximum value for primary, secondary counts (0 = not used). Used to
   select the number of characters needed for this field in the file name.
*/
   int nPD, nSD;

   kIn( prepareNamesVariants );

/****  Get number of varients */

   pC = c->primaryCount;
   sC = c->secondaryCount;

/****
 ****  Perform analysis of primary and secondary variant counts
 ****     Determine the largest index that will be used (nPD, nSD)
 ****     Generate the upper loop bounds used in various places (nP,nS)
 ****/

   if( pC < 1 ) {
      nPD =  nSD = 0;

      *nP = 1;
      *nS = (int *)malloc( sizeof(int) );
      (*nS)[0] = 1;

      *nSmax = 1;
      *success = (int *)malloc( sizeof(int) );
      memset( (*success), 0, sizeof(int) );
   }
   else if( sC == 0 ) {
      nPD = pC;
      nSD = 0;

      *nP = pC;
      *nS = (int *)malloc( pC * sizeof(int) );
      for( iP=0; iP<pC; ++iP ) (*nS)[ iP ] = 1;

      *nSmax = 1;
      *success = (int *)malloc( pC * sizeof(int) );
      memset( (*success), 0, pC * sizeof(int) );
   }
   else {
      nPD = pC;
      nSD = 0;

      *nP = pC;
      *nS = (int *)malloc( pC * sizeof(int) );

      for( iP=0; iP<(*nP); ++iP ) {
	 n = sC[ iP ];
	 (*nS)[ iP ] = ( n<1 ) ? 1 : n;
         if( nSD < n ) nSD = n;
      }

      if( nSD < 1 ) nSD = 0;

      *nSmax = ( nSD < 1 ) ? 1 : nSD;
      n = pC * (*nSmax);
      *success = (int *)malloc( n * sizeof(int) );
      memset( (*success), 0, n * sizeof(int) );
   }

/****  Select format */

   if( nPD > 1000 ) {
      kErr3( number of primary variants for %s is %d; should be <= 1000,
         c->variableName, nPD );
      goto done;
   }
   else if( nPD > 100 ) {
      strcpy( pF, "%03d" );
   }
   else if( nPD > 10 ) {
      strcpy( pF, "%02d" );
   }
   else if( nPD > 0 ) {
      strcpy( pF, "%1d" );
   }
   else {
      *pF = 0;
   }

   if( nSD > 1000 ) {
      kErr3( number of secondary variants for %s is %d; should be <= 1000,
         c->variableName, nSD );
      goto done;
   }
   else if( nSD > 100 ) {
      strcpy( sF, "%03d" );
   }
   else if( nSD > 10 ) {
      strcpy( sF, "%02d" );
   }
   else if( nSD > 0 ) {
      strcpy( sF, "%1d" );
   }
   else {
      *sF = 0;
   }

done:
   kOut;
   return( kStatus );
}

/*******************
********************  Prepare output names
*******************/

static int prepareNames( void )
{
   int i, iC, iH, status;
   char *s, pTag, sTag, pF[10], sF[10], b[256];
   static char *format = "%s";
   Helper *h;
   struct MfixReaderDataChoice *c;
   kIn( prepareNames );

/****
*****  Split the input file name into its parts
****/

   s = strrchr( inputFileName, '/' );

   if( s == 0 ) {
      dirName = (char *)malloc(1);
      dirName[0] = 0;
      baseName = strdup( inputFileName );
   }
   else {
      *s = 0;
      dirName = (char *)malloc( 2 + strlen( inputFileName ) );
      strcpy( dirName, inputFileName );
      strcat( dirName, "/" );
      baseName = strdup( s+1 );

      if( strlen( baseName ) < 1 ) {
	 kErr2( input file name %s has no file name part, inputFileName );
	 goto done;
      }
   }

   s = strrchr( baseName, '.' );
   if( s != 0 ) *s = 0;

/****
*****  Form the case and geometry names
****/

   strcpy( b, baseName );           /* For now don't use dirName */
   strcat( b, caseFileExtension );
   casePathName = strdup(b);

   strcpy( b, baseName );
   strcat( b, caseFileExtension );
   caseFileName = strdup(b);

   strcpy( b, baseName );          /* For now don't use dirName */
   strcat( b, geoFileExtension );
   geoPathName = strdup(b);

   strcpy( b, baseName );
   strcat( b, geoFileExtension );
   geoFileName = strdup(b);

/****
*****  Create the formats used for each raw variable
****/

   for( iC=0, c=choices, h=helper; iC<nChoices; ++iC, ++c, ++h ) {
      for( i=0; i<3; ++i ) {
	 if( h->format[i] != MRF_none ) break;
      }

      if( i >= 3 ) continue;   /* Occurs when no formats selected */

   /****  Build the formats used for the primary and secondary variants */

      status = prepareNamesVariants( c, pF, sF, &(h->nP), &(h->nS),
	 &(h->nSmax), &(h->success) );

      if( !status ) goto done;

   /****  Select the variant tags as the first letter of the variant name */

      if( sF[0] != 0 ) {
	 pTag = (c->primaryName)[0];
	 sTag = (c->secondaryName)[0];
      }
      else if( pF[0] != 0 ) {
	 pTag = (c->primaryName)[0];
	 sTag = '?';
      }
      else {
	 pTag = sTag = '?';
      }

   /****  Build formats using both primary and secondary variants */

      if( sF[0] != 0 ) {
	 h->nArgs = 3;

	 sprintf( b, "%s_%c%s%c%s%s",
	    c->variableName, pTag, pF, sTag, sF, format );
	 h->nameFormat = strdup(b);

	 sprintf( b, "%s_%s_%c%s%c%s%s",
	    baseName, c->variableName, pTag, pF, sTag, sF, format );

	 h->fileFormat = strdup(b);
      }

   /****  Build formats using only primary variant */

      else if ( pF[0] != 0 ) {
	 h->nArgs = 2;

	 sprintf( b, "%s_%c%s%s", c->variableName, pTag, pF, format ); 
	 h->nameFormat = strdup(b);

	 sprintf( b, "%s_%s_%c%s%s",
	    baseName, c->variableName, pTag, pF, format );
	 h->fileFormat = strdup(b);
      }

   /****  Build formats using neither primary nor secondary variants */

      else {
	 h->nArgs = 1;

	 sprintf( b, "%s%s", c->variableName, format );
	 h->nameFormat = strdup(b);

	 sprintf( b, "%s_%s%s", baseName, c->variableName, format );
	 h->fileFormat = strdup(b);
      }
   }

/****
*****  Create the formats used for calculated variables, if any
****/

   if( nHelper > nChoices ) {
      for( iH=nChoices, h=helper+nChoices; iH<nHelper; ++iH, ++h ) {
	 for( i=0; i<3; ++i ) {
	    if( h->format[i] != MRF_none ) break;
	 }

	 if( i >= 3 ) continue;   /* Occurs when no formats selected */

      /****  Build formats using neither primary nor secondary variants */

	 h->nArgs = 1;

	 sprintf( b, "%s%s", h->name, format );
	 h->nameFormat = strdup(b);

	 sprintf( b, "%s_%s%s", baseName, h->name, format );
	 h->fileFormat = strdup(b);
      }
   }

 done:
   kOut;
   return( kStatus );
}

/******************************************************************************/
/*******************                                        *******************/
/*******************           Write information            *******************/
/*******************                                        *******************/
/******************************************************************************/

/*******************
********************  Write one variable
********************    Data array must include a value for all nodes or all
********************    cells even if some nodes or cells are unused. Data for
********************    unused nodes is written but not used by Ensight. Data
********************    for unused cells is not written.
*******************/

int writeVariable(
   FILE *file,
   char *description,
   int veclen,
   enum MfixReaderFormat format,
   float *dataF
)
{
   int i, *cf, iElement;
   float *t, *w, *W;
   kIn( writeVariable );

/****
*****  Header
****/

   write80( file, "BEGIN TIME STEP" );
   write80( file, description );
   write80( file, "part" );

   switch( (int)format ) {
      case MRF_node_average:
	 fwrite( &gridPartIndex, sizeof(int), 1, file );
	 write80( file, "coordinates" );
	 break;

      case MRF_cell:
	 fwrite( &gridPartIndex, sizeof(int), 1, file );

	 if( nZC > 1 )
	    write80( file, "hexa8" );
	 else
	    write80( file, "quad4" );
	 break;

      case MRF_node_center:
	 fwrite( &pseudoGridPartIndex, sizeof(int), 1, file );
	 write80( file, "coordinates" );
	 break;
   }

/****  Data is written from the pointer W. When writing scalars, W is the
*****  input array. When writing vectors, W is the temporary buffer "dataBuf"
*****  allocated elsewhere. To write a vector, the components are moved from
*****  the input array into dataBuf, one at a time, and written. In other
*****  words, Ensight wants all the values for the first component together
*****  followed by all the values for the second component and so on.
****/

   W = ( veclen == 1 ) ? dataF : dataBuf;

   switch( (int)format ) {

   /****  Write one scalar or vector per grid cell */
   /****     Do not write data to cells that are not used */

      case MRF_cell:
	 kBug2( 7, writing cell data );

	 for( iElement=0; iElement<veclen; ++iElement ) {
	    if( veclen == 3 ) {
	       for( i=0, w=W, t=dataF+iElement; i<nXYZC; ++i, ++w, t+=3 )
		  *w = *t;
	    }

	    for( i=0, w=W, cf=cellFlags; i<nXYZC; ++i, ++w, ++cf ) {
	       if( (*cf) != MRT_fluid ) continue;
	       fwrite( w, sizeof(float), 1, file );
	    }
	 }

	 break;

   /****  Write one scalar or vector per grid node */
   /****     Node data is needed even if the node is not used */

      case MRF_node_average:
	 kBug2( 7, writing node-average data );

	 for( iElement=0; iElement<veclen; ++iElement ) {
	    if( veclen == 3 ) {
	       for( i=0, w=W, t=dataF+iElement; i<nXYZN; ++i, ++w, t+=3 ) *w=*t;
	    }

	    fwrite( W, sizeof(float), nXYZN, file );
	 }

	 break;

   /****  Write one scalar or vector per pseudo grid node */
   /****     Node data is needed even if the node is not used */

      case MRF_node_center:
	 kBug2( 7, writing node-center data );

	 for( iElement=0; iElement<veclen; ++iElement ) {
	    if( veclen == 3 ) {
	       for( i=0, w=W, t=dataF+iElement; i<nXYZPN; ++i, ++w, t+=3 )*w=*t;
	    }

	    fwrite( W, sizeof(float), nXYZPN, file );
	 }

	 break;
   }

   write80( file, "END TIME STEP" );
   kOut;
   return( kStatus );
}

/*******************
********************  Get data for one variable
*******************/

int getVariable(
   char *name,
   int indexPrimary,
   int indexSecondary,
   enum MfixReaderFormat format,
   int iTime,
   int **dataInt,
   float **dataFloat
)
{
   struct MfixReaderDataRequest R;
   kIn( getVariable );

   R.dims = 0;
   R.variableName = name;
   R.primaryNum = indexPrimary;
   R.secondaryNum = indexSecondary;
   R.component = -1;
   R.format = format;

   if( !MfixReaderGetData( handle, times[ iTime ], 1, &R ) ) goto done;

   if( R.dataType == MRD_integer )
      *dataInt = R.dataInt;
   else if( R.dataType == MRD_float )
      *dataFloat = R.dataFloat;

done:
   kOut;
   return( kStatus );
}

/*******************
********************  Write one variant at one time
*******************/

static int writeOneVariant( FILE *file, int iC, int iT, int iP, int iS,
   enum MfixReaderFormat f, int *successful )
{
   int i, n, *S, nXYZvector, freeT;
   float *T, xv, yv, zv;
   char b[256];
   Helper *h;
   struct MfixReaderDataChoice *c;
   struct MfixReaderDataRequest R;
   kIn( writeOneVariant );

   kBug7( 7, iC=%d iT=%d iP=%d iS=%d format=%d, iC, iT, iP, iS, (int)f );

   freeT = 0;
   c = choices + iC;
   h = helper + iC;
   *successful = 0;

/****  Determine array size */

   switch( (int)f ) {
      case MRF_node_average:
	 nXYZvector = c->veclen * nXYZN;
	 break;
      case MRF_cell:
	 nXYZvector = c->veclen * nXYZC;
	 break;
      case MRF_node_center:
	 nXYZvector = c->veclen * nXYZPN;
	 break;
   }

/********************************  Get data  **********************************/

/****  Create request */

   R.dims = 0;
   R.variableName = c->variableName;
   R.primaryNum = iP;
   R.secondaryNum = iS;
   R.component = -1;
   R.format = f;

/****  Supply test data */

   if( testRequested ) {
      kBug3( 7, %d values from test, nXYZvector );

      R.veclen = c->veclen;
      R.dims = nCells;
      R.dataType = MRD_float;
      R.dataFloat = dataBuf;
      memset( dataBuf, 0, 3 * nXYZC * sizeof(float) );

      if( c->veclen == 1 ) {
	 if( testRequested == 1 ) {
	    dataBuf[13] = (float)iT;
	 }
	 else if( testRequested == 2 ) {
	    n = nXYC * ( nZC / 2 );
	    for( i=0; i<n; ++i ) dataBuf[i] = 1.0;
	 }
      }
      else {
	 switch( iT ) {
	    case 0: xv = 1.0; yv = zv = 0.0; break;
	    case 1: yv = 1.0; xv = zv = 0.0; break;
	    case 2: zv = 1.0; xv = yv = 0.0; break;
	 }

	 for( i=0; i<3*nXYZC; i+=3 ) {
	    dataBuf[i] = xv;
	    dataBuf[ i + 1 ] = yv;
	    dataBuf[ i + 2 ] = zv;
	 }
      }
   }

/****  Supply cell types already in memory */

   else if( !strcmp( R.variableName, "celltype" ) ) {
      R.veclen = 1;
      R.dataType = MRD_integer;

      switch( (int)R.format ) {
	 case MRF_cell:
	    R.dataInt = cellFlags;
	    kBug3( 7, %d values from cellflags, nXYZvector );
	    break;
	 case MRF_node_average:
	    R.dataInt = nodeFlags;
	    kBug3( 7, %d values from nodeflags, nXYZvector );
	    break;
	 case MRF_node_center:
	    R.dataInt = cellFlags;
	    kBug3( 7, %d values from cellflags, nXYZvector );
	    break;
      }
   }

/****  Read data */

   else {
      kBug3( 7, %d values from file, nXYZvector );
      if( MfixReaderGetData( handle, times[ iT ], 1, &R ) == 0 ) goto wrapup;
   }

/**********************  Convert to float if needed  **************************/

   if( R.dataType == MRD_integer ) {
      S = R.dataInt;

      if( !S ) {
	 printf( "Warning: m2e could not obtain data for %s\n",
	    R.variableName );

	 printf( "   primary=%d, secondary=%d, format=%d\n",
	    R.primaryNum, R.secondaryNum, R.format );

	 goto wrapup;
      }

      freeT = 1;
      T = (float *)malloc( nXYZvector * sizeof(float) );
      for( i=0; i<nXYZvector; ++i ) T[i] = (float)(S[i]);
   }
   else if( R.dataType == MRD_float ) {
      T = R.dataFloat;

      if( !T ) {
	 printf( "Warning: m2e could not obtain data for %s\n",
	    R.variableName );

	 printf( "   primary=%d, secondary=%d, format=%d\n",
	    R.primaryNum, R.secondaryNum, R.format );

	 goto wrapup;
      }
   }

/****************************  Create description  ****************************/

   n = strlen( c->variableDesc );
   if( n>72 ) n = 72;

   if( n<1 ) {
      strcpy( b, "unknown" );
   }
   else {
      strncpy( b, c->variableDesc, n );
      b[n] = 0;
   }

/*****************************  Write data  ***********************************/

   writeVariable( file, b, c->veclen, f, T );
   *successful = 1;

wrapup:
   if( freeT ) free(T);
   kOut;
   return( kStatus );
}

/*******************
********************  Write all selected variables
*******************/

static int writeData( void )
{
   int iC, iH, iP, nP, iS, nS, iT, nSmax, status, nTimes2, iFormat, created;
   int *success, successful;
   char fileName[256], msg[280];
   Helper *h;
   FILE *file;
   enum MfixReaderFormat f;
   struct MfixReaderDataChoice *c;
   kIn( writeData );

   file = 0;
   dataBuf = (float *)malloc( 3 * nXYZN * sizeof(float) );

/****  Write raw variables - loop over all available variables */

   for( iC=0, c=choices, h=helper; iC<nChoices; ++iC, ++c, ++h ) {
      success = h->success;
      nSmax = h->nSmax;

   /****  Loop over all possible formats */

      for( iFormat=0; iFormat<3; ++iFormat ) {

      /****  Skip unneeded formats */

	 f = h->format[ iFormat ];
	 if( f == MRF_none ) continue;

      /****  Prepare */

	 kBug4( 5, starting variable %s format %s, choices[ iC ].variableName,
	    MfixReaderFormatText[ (int)f ] );

	 nTimes2 = ( h->useTime ) ? nTimes : 1;
	 nP = h->nP;

      /****  Loop over primary variants */

	 for( iP=0; iP<nP; ++iP ) {
	    nS = (h->nS)[ iP ];

	 /****  Loop over secondary variants */

	    for( iS=0; iS<nS; ++iS ) {

	    /****
	    *****  Open an output file for each primary/secondary combo for
	    *****  each variable.
	    ****/

	       makeFileName( iC, iP, iS, f, h->veclen, fileName );
	       if( openFile( fileName, &file ) == 0 ) goto wrapup;

	    /****  Put all times into the same file */

	       for( iT=0; iT<nTimes2; ++iT ) {
		  status = writeOneVariant( file, iC, iT, iP, iS, f,
		     &successful );

		  if( !status ) goto wrapup;
		  if( !successful ) break;
	       }

	    /****  Close output file and record the success status */

	       fclose( file );
	       file = 0;
	       success[ iP*nSmax + iS ] = successful;

	    /****  On failure erase the output file */

	       if( !successful ) {
		  sprintf( msg, "rm %s", fileName );
		  system( msg );
	       }
	    }
	 }
      }
   }

/****  Write calculated variables if any */

   if( nHelper > nChoices ) {
      for( iH=nChoices, h=helper+nChoices; iH<nHelper; ++iH, ++h ) {
	 for( iFormat=0; iFormat<3; ++iFormat ) {

	 /****  Do nothing if variable not needed in this format */

	    f = h->format[ iFormat ];
	    if( f == MRF_none ) continue;

	 /****  Request the variable */
	 /****  This is the first place where we know the veclen and times */

	    calc( h->name, f, &created, &(h->veclen), &(h->useTime) );
	    if( !kStatus ) goto wrapup;

	 /****  Special case: veclen=0 means the var couldn't be calculated */

	    if( !created ) h->format[ iFormat ] = MRF_none;
	 }
      }
   }

wrapup:
   if( file != 0 ) fclose( file );
   kOut;
   return( kStatus );
}

/*******************
********************  Write one drawing to the geo file
*******************/

static int writeGeoDrawing( int index, FILE *file )
{
   int *elements, nElements;
   kIn( writeGeoDrawing );

   nElements = drawings[ index ].nE;
   elements = drawings[ index ].e;

   write80( file, "part" );
   fwrite( &nextPartIndex, sizeof(int), 1, file );
   ++nextPartIndex;
   write80( file, drawings[ index ].name );
   write80( file, "coordinates" );
   fwrite( &nXYZN, sizeof(int), 1, file );
   fwrite( xyz, sizeof(float), 3*nXYZN, file );
   write80( file, drawings[ index ].eName );
   fwrite( &nElements, sizeof(int), 1, file );
   fwrite( elements, sizeof(int), 4 * nElements, file );

   kOut;
   return( kStatus );
}

/*******************
********************  Write the geometry for the pseudo grid
*******************/

static int writeGeoPseudoGrid( FILE *file )
{
   int iX, iY, iZ, i[8], *pcf;
   kIn( writeGeoPseudoGrid );

   write80( file, "part" );
   fwrite( &pseudoGridPartIndex, sizeof(int), 1, file );
   write80( file, "cell-center data" );
   write80( file, "coordinates" );
   fwrite( &nXYZPN, sizeof(int), 1, file );
   fwrite( xyzp, sizeof(float), 3*nXYZPN, file );
   write80( file, "hexa8" );
   fwrite( &nPseudoInside, sizeof(int), 1, file );

   pcf = pseudoGridCellFlags;

   for( iZ=0; iZ<nZPC; ++iZ ) {
      for( iY=0; iY<nYPC; ++iY ) {
	 for( iX=0; iX<nXPC; ++iX, ++pcf ) {
	    if( (*pcf) != MRT_fluid ) continue;

            i[0] = iX + 2 + iY * nXPN + iZ * nXYPN;
	    i[1] = i[0] - 1;
	    i[2] = i[1] + nXYPN;
	    i[3] = i[2] + 1;
	    i[4] = i[0] + nXPN;
	    i[5] = i[1] + nXPN;
	    i[6] = i[2] + nXPN;
	    i[7] = i[3] + nXPN;

	    fwrite( i, sizeof(int), 8, file );
	 }
      }
   }

   kOut;
   return( kStatus );
}

/*******************
********************  Write the geometry for the grid
*******************/

static int writeGeoGrid( FILE *file )
{
   int n, iX, iY, iZ, i[8], *cf;
   kIn( writeGeoGrid );

   write80( file, "part" );
   fwrite( &gridPartIndex, sizeof(int), 1, file );
   write80( file, "grid for cell data" );
   write80( file, "coordinates" );
   fwrite( &nXYZN, sizeof(int), 1, file );
   fwrite( xyz, sizeof(float), 3*nXYZN, file );

   if( nZC > 1 ) {
      n = 8;
      write80( file, "hexa8" );
   }
   else {
      n = 4;
      write80( file, "quad4" );
   }

   fwrite( &nCellsInside, sizeof(int), 1, file );

   cf = cellFlags;

   for( iZ=0; iZ<nZC; ++iZ ) {
      for( iY=0; iY<nYC; ++iY ) {
	 for( iX=0; iX<nXC; ++iX, ++cf ) {
	    if( (*cf) != MRT_fluid ) continue;

	 /****  quad4 or -Z side of hexa8 */

	    i[0] = iX + 1 + iY * nXN + iZ * nXYN;
            i[1] = i[0] + 1;
	    i[2] = i[1] + nXN;
	    i[3] = i[2] - 1;

	 /****  +Z side of hexa8 */

	    if( n == 8 ) {
	       i[4] = i[0] + nXYN;
	       i[5] = i[1] + nXYN;
	       i[6] = i[2] + nXYN;
	       i[7] = i[3] + nXYN;
	    }

/*
            i[0] = iX + 2 + iY * nXN + iZ * nXYN;
	    i[1] = i[0] - 1;
	    i[2] = i[1] + nXYN;
	    i[3] = i[2] + 1;

	    if( n == 8 ) {
	       i[4] = i[0] + nXN;
	       i[5] = i[1] + nXN;
	       i[6] = i[2] + nXN;
	       i[7] = i[3] + nXN;
             }
*/

	    fwrite( i, sizeof(int), n, file );
	 }
      }
   }

   kOut;
   return( kStatus );
}

/*******************
********************  Write the geometry file
*******************/

static int writeGeo( void )
{
   int i, n;
   char buf[80];
   FILE *file;
   kIn( writeGeo );
   file = 0;

   gridPartIndex = nextPartIndex;
   ++nextPartIndex;

   if( usingPseudoGrid ) {
      pseudoGridPartIndex = nextPartIndex;
      ++nextPartIndex;
   }

   if( openFile( geoPathName, &file ) == 0 ) goto wrapup;
   write80( file, "C Binary" );
   write80( file, "BEGIN TIME STEP" );

   sprintf( buf, "See file %s for details", caseFileName );
   write80( file, buf );

   sprintf( buf, "Created %s", asctime( localtime( &currentTime ) ) );
   write80( file, buf );

   write80( file, "node id off" );
   write80( file, "element id off" );

   writeGeoGrid( file );
   if( usingPseudoGrid ) writeGeoPseudoGrid( file );

   if( drawings ) {
      n = drawings->nDrawings;
      for( i=0; i<n; ++i ) writeGeoDrawing( i, file );
   }

wrapup:
   if( file != 0 ) {
      write80( file, "END TIME STEP" );
      fclose( file );
   }

   kOut;
   return( kStatus );
}

/*******************
********************  Write the case file
*******************/

static int writeCase( void )
{
   int i, n, pn, iH, iC, iP, iS, nP, nS, iT, vec, lastT, iFormat, *success;
   char *varDecl, fileName[256], varName[256];
   FILE *file;
   Helper *h;
   struct MfixReaderDataChoice *c;
   Drawing *d;
   enum MfixReaderFormat f;
   kIn( writeCase );
   file = 0;

   if( openFile( casePathName, &file ) == 0 ) goto wrapup;

/***************  Write comments to document the data set  ********************/

   fprintf( file, "# Created by m2e (MFIX to Ensight) version %d.%d\n",
      m2eVersionMajor, m2eVersionMinor );
   fprintf( file, "#    %s", asctime( localtime( &currentTime ) ) );
   fprintf( file, "# Information about the input MFIX files\n" );
   fprintf( file, "#    For more on MFIX see www.mfix.org\n" );
   fprintf( file, "#    Input files base name was %s\n", baseName );

   if( testRequested ) {
      fprintf( file, "#    Generating test case type %d\n", testRequested );
   }
   else {
      fprintf( file, "#    Input files were in directory %s\n", dirName );
      showInfo( "#    ", file );
   }

   fprintf( file, "# Information about the translation\n" );
   fprintf( file, "#    For more on m2e contact eschenbe@psc.edu\n" );
   fprintf( file, "#    CVS id: %s\n", M2E_Id+10 );
   fprintf( file, "#    Input file endian: %d\n", endian );
   fprintf( file, "#    Input file recl: %d\n", recl );
   fprintf( file, "#    Time tolerance: %8.6f\n", timeTolerance );
   fprintf( file, "#    Time threshold: %8.6f\n", timeThreshold );

   fprintf( file, "#    Extrapolation rule: %s \n",
      MfixReaderExtrapolationText[ (int)extrapRule ] );

   fprintf( file, "#    Interpolation rule: %s \n",
      MfixReaderInterpolationText[ (int)interpRule ] );

   fprintf( file, "#    Time method: %s\n", timeMethodText[ timeMethod ] );
   fprintf( file, "#    Time span: %8.6f %8.6f\n", timeBounds[0],timeBounds[1]);
   fprintf( file, "#    Time increment: %8.6f\n", timeIncrement );

   fprintf( file, "#    Cells skipped: %d %d %d %d %d %d\n", skip[0], skip[1],
      skip[2], skip[3], skip[4], skip[5]);

   fprintf( file, "#    Overwrite: %d\n", overwrite );

   if( nVarRequest > 0 ) {
      fprintf( file, "#    Requested variables:\n" );

      for( i=0; i<nVarRequest; ++i )
	 fprintf( file, "#       %s, %s\n", varRequest[i].name,
	    MfixReaderFormatText[ (int)varRequest[i].format ] );
   }

   if( drawings ) {
      fprintf( file, "# Information about the generated parts\n" );
      d = drawings;
      n = d->nDrawings;
      pn = nextPartIndex + 1 + usingPseudoGrid;

      for( i=0; i<n; ++i, ++d, ++pn ) {
	 fprintf( file,
	    "#    part %d \"%s\" uses %d quads to surround %d cells\n",
	    pn, d->name, d->nE, d->nCells );
      }
   }

/***************  Write the format and geometry sections  *********************/

   fprintf( file, "FORMAT\ntype: ensight gold\n" );
   fprintf( file, "GEOMETRY\nmodel: 1 1 %s\n", geoFileName );

/*********************  Write the variable section  ***************************/

   fprintf( file, "VARIABLE\n" );

/****  Raw variables */

   for( iC=0, c=choices, h=helper; iC<nChoices; ++iC, ++c, ++h ) {
      for( iFormat=0; iFormat<3; ++iFormat ) {
	 f = h->format[ iFormat ];
	 if( f == MRF_none ) continue;
	 success = h->success;

      /****  Decide whether variable is a scalar or a vector */

	 if( h->veclen == 1 ) {
	    vec = 0;
	 }
	 else if( h->veclen == 3 ) {
	    vec = 1;
	 }

	 switch(f) {
	    case MRF_cell:
	       if( vec )
		  varDecl = "vector per element:";
	       else
		  varDecl = "scalar per element:";
               break;
	    default:
	       if( vec )
		  varDecl = "vector per node:";
	       else
		  varDecl = "scalar per node:";
	 }

      /****  Generate a line for each variant of this variable */

	 nP = h->nP;

	 for( iP=0; iP<nP; ++iP ) {
	    nS = (h->nS)[ iP ];

	    for( iS=0; iS<nS; ++iS ) {
	       if( !success[ iP*(h->nSmax) + iS ] ) continue;

	       makeVarName( iC, iP, iS, f, varName );
	       makeFileName( iC, iP, iS, f, h->veclen, fileName );

	       if( h->useTime )
		  fprintf( file, "%s 2 2 %12s %s\n", varDecl, varName,fileName);
	       else
		  fprintf( file, "%s 1 1 %12s %s\n", varDecl, varName,fileName);
	    }
	 }
      }
   }

/****  Calculated variables */

   if( nHelper > nChoices ) {
      for( iH=nChoices, h=helper+nChoices; iH<nHelper; ++iH, ++h ) {
	 if( h->veclen == 0 ) continue;

	 for( iFormat=0; iFormat<3; ++iFormat ) {
	    f = h->format[ iFormat ];
	    if( f == MRF_none ) continue;

	 /****  Decide whether variable is a scalar or a vector */

	    if( h->veclen == 1 ) {
	       vec = 0;
	    }
	    else if( h->veclen == 3 ) {
	       vec = 1;
	    }

	    switch(f) {
	       case MRF_cell:
		  if( vec )
		     varDecl = "vector per element:";
		  else
		     varDecl = "scalar per element:";
		  break;
	       default:
		  if( vec )
		     varDecl = "vector per node:";
		  else
		     varDecl = "scalar per node:";
	    }

	 /****  Generate a line for this variable */

	    makeVarName( iH, 0, 0, f, varName );
	    makeFileName( iH, 0, 0, f, h->veclen, fileName );

	    if( h->useTime )
	       fprintf( file, "%s 2 2 %12s %s\n", varDecl, varName, fileName );
	    else
	       fprintf( file, "%s 1 1 %12s %s\n", varDecl, varName, fileName );
	 }
      }
   }

/***********************  Write the time section  *****************************/

   fprintf( file, "TIME\n" );
   fprintf( file, "time set: 1\n" );
   fprintf( file, "number of steps: 1\n" );
   fprintf( file, "time values: 0.0\n" );

   fprintf( file, "time set: 2\n" );
   fprintf( file, "number of steps: %d\n", nTimes );

   for( iT=0; iT<nTimes; iT+=6 ) {
      lastT = iT + 5;
      if( lastT >= nTimes ) lastT = nTimes - 1;

      if( iT == 0 )
	 fprintf( file, "time values:" );
      else
	 fprintf( file, "            " );

      for( i=iT; i<=lastT; ++i ) fprintf( file, " %10.6f", times[i] );
      fprintf( file, "\n" );
   }

/************************  Write the fileset section  *************************/

   fprintf( file, "FILE\n" );
   fprintf( file, "file set: 1\n" );
   fprintf( file, "number of steps: 1\n" );
   fprintf( file, "file set: 2\n" );
   fprintf( file, "number of steps: %d\n", nTimes );

wrapup:
   if( file != 0 ) fclose( file );
   kOut;
   return( kStatus );
}

/******************************************************************************/
/*******************                                        *******************/
/*******************          Prepare information           *******************/
/*******************                                        *******************/
/******************************************************************************/

/*******************
********************  Prepare irregular coordinates
*******************/

static int prepareCoordinates( void )
{
   int i, j, iX, iY, iZ, jX, jY, jZ, *dims;
   float *g, *xyzTest;
   char *units;
   kIn( prepareCoordinates );

/**************************  Prepare basic grid  ******************************/

   xyz = (float *)malloc( 3 * nXYZN * sizeof(float) );

/****  Generate test grid or read real grid */

   if( testRequested ) {
      xyzTest = (float *)malloc( 3 * nXYZN * sizeof(float) );
      g = xyzTest;

      for( iZ=0; iZ<nZN; ++iZ ) {
	 for( iY=0; iY<nYN; ++iY ) {
	    for( iX=0; iX<nXN; ++iX, g+=3 ) {
	       g[0] = iX;
	       g[1] = iY;
	       g[2] = iZ;
	    }
	 }
      }

      g = xyzTest;
   }
   else {
      xyzTest = 0;
      MfixReaderGetCoordinates( handle, MRF_cell, &units, &dims, &g );
      if( !kStatus ) goto wrapup;
   }

/****  Rearrange for Ensight: all X then all Y then all Z */

   for( iZ=jZ=0; iZ<nZN; ++iZ, jZ+=nXYN ) {
      for( iY=0, jY=jZ; iY<nYN; ++iY, jY+=nXN ) {
	 for( iX=0, jX=jY; iX<nXN; ++iX, ++jX ) {
	    for( i=0, j=jX; i<3; ++i, ++g, j+=nXYZN ) xyz[j] = *g;
	 }
      }
   }

/************************  Prepare pseudo grid  *******************************/

   if( usingPseudoGrid ) {
      xyzp = (float *)malloc( 3 * nXYZPN * sizeof(float) );

      MfixReaderGetCoordinates( handle, MRF_node_center, &units, &dims, &g );
      if( !kStatus ) goto wrapup;

   /****  Rearrange for Ensight: all X then all Y then all Z */

      for( iZ=jZ=0; iZ<nZPN; ++iZ, jZ+=nXYPN ) {
	 for( iY=0, jY=jZ; iY<nYPN; ++iY, jY+=nXPN ) {
	    for( iX=0, jX=jY; iX<nXPN; ++iX, ++jX ) {
	       for( i=0, j=jX; i<3; ++i, ++g, j+=nXYZPN ) xyzp[j] = *g;
	    }
	 }
      }
   }

wrapup:
   if( xyzTest ) free( xyzTest );
   kOut;
   return( kStatus );
}

/*******************
********************  Prepare cell types
*******************/

static int prepareCelltypes( void )
{
   int i, iX, iY, iZ, *cf;
   struct MfixReaderDataRequest R;
   kIn( prepareCelltypes );

/*****************************  Read cell flags  ******************************/

   if( testRequested == 0 ) {
      R.dims = 0;
      R.variableName = "celltype";
      R.primaryNum = 0;
      R.secondaryNum = 0;
      R.component = -1;
      R.format = MRF_cell;

      if( MfixReaderGetData( handle, 0.0, 1, &R ) == 0 ) goto done;

      if( R.dims == 0 ) {
	 kErr1( could not obtain celltype );
	 goto done;
      }

      cellFlags = (int *)malloc( nXYZC * sizeof(int) );
      memcpy( cellFlags, R.dataInt, nXYZC * sizeof(int) );
   }

/***************************  Read node flags  ********************************/

   if( needNodeFlags ) {
      R.dims = 0;
      R.variableName = "celltype";
      R.primaryNum = 0;
      R.secondaryNum = 0;
      R.component = -1;
      R.format = MRF_node_average;

      if( MfixReaderGetData( handle, 0.0, 1, &R ) == 0 ) goto done;

      if( R.dims == 0 ) {
	 kErr1( could not obtain celltype at nodes );
	 goto done;
      }

      nodeFlags = (int *)malloc( nXYZN * sizeof(int) );
      memcpy( nodeFlags, R.dataInt, nXYZN * sizeof(int) );
   }

/*******************  Create pseudo grid cell flags  **************************/

   if( usingPseudoGrid ) {
      cf = pseudoGridCellFlags = (int *)malloc( nXYZPC * sizeof(int) );

      for( iZ=1; iZ<nZC; ++iZ ) {
	 for( iY=1; iY<nYC; ++iY ) {
	    for( iX=1; iX<nXC; ++iX, ++cf ) {
	       *cf = nodeFlags[ iX + iY*nXN + iZ*nXYN ];
	    }
	 }
      }
   }

/**********************  Calculate number of cells inside  ********************/

/****  Main grid */

   nCellsInside = nXYZC;
   cf = cellFlags;

   for( i=0; i<nXYZC; ++i, ++cf ) {
      if( *cf != MRT_fluid ) --nCellsInside;
   }

/****  Pseudo grid */

   if( usingPseudoGrid ) {
      nPseudoInside = nXYZPC;
      cf = pseudoGridCellFlags;

      for( i=0; i<nXYZPC; ++i, ++cf ) {
	 if( *cf != MRT_fluid ) --nPseudoInside;
      }
   }

done:
   kOut;
   return( kStatus );
}

/*******************
********************  Add one time to the list
*******************/

static int addTime( float newTime )
{
   int iT, iBest;
   float t, d, tBest, dBest, *timesOld;
   kIn( addTime );

/****  This is the first - merely stuff the new time into the list */

   if( nTimes < 1 ) {
      nTimesAlloc = 500;
      times = (float *)malloc( nTimesAlloc * sizeof(float) );
      times[0] = newTime;
      nTimes = 1;
      goto done;
   }

/****  Look for the time in the list that is closest to the new time */

   iBest = 0;
   tBest = times[0];
   dBest = tBest - newTime;
   if( dBest < 0.0 ) dBest = -dBest;
   if( dBest <= timeThreshold ) goto done;  /* No need to add to list */

   for( iT=1; iT<nTimes; ++iT ) {
      t = times[ iT ];
      d = t - newTime;
      if( d < 0.0 ) d = -d;
      if( d <= timeThreshold ) goto done;  /* No need to add to list */

      if( dBest > d ) {
	 dBest = d;
	 tBest = t;
	 iBest = iT;
      }
   }

/****  Extend the list to make room for the new time */

   if( nTimes >= nTimesAlloc ) {
      timesOld = times;
      nTimesAlloc += 500;
      times = (float *)malloc( nTimesAlloc * sizeof(float) );
      memcpy( times, timesOld, nTimes * sizeof(float) );
      free( timesOld );
   }

/****  Move old times out of the way, if needed */

   if( tBest <= newTime ) ++iBest;  /* iBest is now the location for the new */

   if( iBest < nTimes ) memmove( times+iBest+1, times+iBest,
      ( nTimes - iBest ) * sizeof(float) );

   times[ iBest ] = newTime;
   ++nTimes;

done:
   kOut;
   return( kStatus );
}

/*******************
********************  Select times for data
*******************/

static int selectTimes( void )
{
   int iC, iT, cNTimes, iFormat;
   float t, *cTimes;
   Helper *h;
   struct MfixReaderDataChoice *c;
   kIn( selectTimes );

/*********************  Create list of times by formula  **********************/

   if( timeMethod == TM_STEP ) {
      nTimes = 1.0 + timeThreshold
	 + ( ( timeBounds[1] - timeBounds[0] ) / timeIncrement );

      if( nTimes < 1 ) nTimes = 1;
      times = (float *)malloc( nTimes * sizeof(float) );

      for( iT=0, t=timeBounds[0]; iT<nTimes; ++iT, t+=timeIncrement )
	 times[ iT ] = t;
   }

/*************  Merge list of times from all used variable  *******************/

   else {

   /****  Merge times */

      for( iC=0, c=choices, h=helper; iC<nChoices; ++iC, ++c, ++h ) {
	 for( iFormat=0; iFormat<3; ++iFormat ) {
	    if( h->format[ iFormat ] != MRF_none ) break;
	 }

	 if( iFormat >= 3 ) continue;

	 cNTimes = c->nTimes;
	 cTimes = c->times;

	 for( iT=0; iT<cNTimes; ++iT ) {
	    addTime( cTimes[ iT ] );
	 }
      }

   /****  Apply bounds if needed */

      if( timeMethod == TM_ALL_BOUNDED ) {

      /****  Look for first acceptible time re lower time bound */

	 for( iT=0; iT<nTimes; ++iT ) {
	    if( timeBounds[0] <= times[ iT ] ) break;
	 }

	 if( iT >= nTimes ) {
	    kErr1( no times are within the requested time span );
	    goto done;
	 }

      /****  Need to eliminate times at the beginning of the list */

	 if( iT > 0 ) {
	    memmove( times, times + iT, ( nTimes - iT ) * sizeof(float) );
	    nTimes -= iT;
	 }

      /****  Look for last acceptible time re upper time bound */

	 for( iT=nTimes-1; iT>=0; --iT ) {
	    if( timeBounds[1] >= times[ iT ] ) break;
	 }

	 if( iT < 0 ) {
	    kErr1( no times are within the requested time span );
	    goto done;
	 }

      /****  This line eliminates times at the end of the list, if needed */

	    nTimes = iT + 1;
      }
   }

/****  Sanity check */

   if( nTimes < 1 ) {
      kErr1( no times selected );
   }

/****  Show times if debugging */

   if( kBugLevel > 1 ) {
      kBug2( 1, Times (sec) );
      for( iT=0; iT<nTimes; ++iT ) kBug4( 1, ... %d: %10.6f, iT, times[iT]);
   }

done:
   kOut;
   return( kStatus );
}

/*******************
********************  Extend helper list
*******************/

static int extendHelper( void )
{
   int i, iH;
   Helper *h, *hh;
   kIn( extendHelper );

   h = helper;

/****  Extend helper list */

   ++nHelper;
   helper = (Helper *)malloc( nHelper * sizeof(Helper) );

/****  Initialize list */

   memset( helper, 0, nHelper * sizeof(Helper) );

   for( iH=0, hh=helper; iH<nHelper; ++iH, ++hh ) {
      for( i=0; i<3; ++i ) hh->format[i] = MRF_none;
   }

/****  Copy data from old version of list then free it */

   memcpy( helper, h, ( nHelper - 1 ) * sizeof(Helper) );
   free(h);

   kOut;
   return( kStatus );
}

/*******************
********************  Select variables
*******************/

static int selectVariables( void )
{
   int i, iH, iC, iVar, nRawVar;
   char *variableName;
   Helper *h;
   struct MfixReaderDataChoice *c;
   VarRequest *v;
   kIn( selectVariables );

   nRawVar = 0;

/****  Initialize helper structure */

   nHelper = nChoices;
   helper = (Helper *)malloc( nHelper * sizeof(Helper) );
   memset( helper, 0, nHelper * sizeof(Helper) );

   for( iH=0, h=helper; iH<nHelper; ++iH, ++h ) {
      for( i=0; i<3; ++i ) h->format[i] = MRF_none;
   }

/*********************  Specific variables selected  **************************/

   if( nVarRequest > 0 ) {
      for( iVar=0, v=varRequest; iVar<nVarRequest; ++iVar, ++v ) {
	 variableName = v->name;

      /****  Look for requested variable in list of raw variables */

	 for( iC=0, c=choices, h=helper; iC<nChoices; ++iC, ++c, ++h ) {
	    if( strcmp( variableName, c->variableName ) == 0 ) break;
	 }

      /****  Not found - assume it is a calculated variable */

	 if( iC >= nChoices ) {
	    extendHelper();
	    h = helper + nHelper - 1;
	    h->veclen = 1;           /* May be changed by "calc" later */
	    h->useTime = 1;          /* May be changed by "calc" later */
	    h->name = variableName;
	 }

      /****  Was found - is raw variable available in current dataset? */

	 else if( !(c->available) ) {
	    kErr2( variable %s is unavailable, variableName );
	    goto done;
	 }

      /****  Was found - set raw variable info */

	 else {
	    ++nRawVar;
	    h = helper + iC;
	    h->name = variableName;
	    h->useTime = ( c->nTimes > 0 );

	    if( c->veclen == 1 ) {
	       h->veclen = 1;
	    }
	    else if( c->veclen == 3 ) {
	       h->veclen = 3;
	    }
	    else {
	       kErr3( variable %s has a vector length of %d but must be 1 or 3,
		  variableName, c->veclen );
	       goto done;
	    }
	 }

      /****  Use variable as cell data if not already doing so */

	 if( v->format == MRF_cell ) {
	    if( h->format[0] ) {
	       kErr2( variable %s (cell data) duplicated, variableName );
	       goto done;
	    }

	    h->format[0] = MRF_cell;
	 }

      /****  Use variable as averaged node data if not already doing so */

	 else if( v->format == MRF_node_average ) {
	    if( h->format[1] ) {
	       kErr2( variable %s (average node data) duplicated, variableName);
	       goto done;
	    }

	    h->format[1] = MRF_node_average;
	    if( !strcmp( variableName, "celltype" ) ) needNodeFlags = 1;
	 }

      /****  Use variable as centered node data if not already doing so */

	 else if( v->format == MRF_node_center ) {
	    if( h->format[2] ) {
	       kErr2( variable %s (center node data) duplicated, variableName );
	       goto done;
	    }

	    h->format[2] = MRF_node_center;
	    usingPseudoGrid = 1;
	 }
      }
   }

/****  No raw variables selected - select all using current format */

   if( nRawVar == 0 ) {
      switch( (int)varFormat ) {
	 case MRF_cell: i = 0; break;
	 case MRF_node_average: i = 1; break;
	 case MRF_node_center: i = 2; usingPseudoGrid = 1; break;
      }

      for( iC=0, c=choices, h=helper; iC<nChoices; ++iC, ++c, ++h ) {
	 if( !(c->available) ) continue;
	 h->format[i] = varFormat;
	 h->useTime = ( c->nTimes > 0 );

	 if( c->veclen == 1 ) {
	    h->veclen = 1;
	 }
	 else if( c->veclen == 3 ) {
	    h->veclen = 3;
	 }
	 else {
	    kErr3( variable %s has a vector length of %d but must be 1 or 3,
	       variableName, c->veclen );
	    goto done;
	 }

	 if( strcmp( c->variableName, "celltype" ) ) continue;
	 if( varFormat == MRF_node_average ) needNodeFlags = 1;
      }
   }

/****  Make a few checks needed for the pseudo-grid */

   if( usingPseudoGrid ) {
      needNodeFlags = 1;

      if( ( nXN < 3 ) || ( nYN < 3 ) || ( nZN < 3 ) ) {
	 kErr1( grid too small to create node data with the cell-center method);
	 goto done;
      }
   }

/****  Insure at least one variable */

   for( iH=0, h=helper; iH<nHelper; ++iH, ++h ) {
      for( i=0; i<3; ++i ) {
	 if( h->format[i] != MRF_none  ) goto done;
      }
   }

   kErr1( no variables were selected );

done:
   kOut;
   return( kStatus );
}

/*******************
********************  Add a void tube to the blanking
*******************/

static void addTube( int x0, int x1, int y0, int y1 )
{
   int iX, iY, iZ;
   kIn( addTube);

   for( iY=y0; iY<=y1; ++iY ) {
      for( iX=x0; iX<=x1; ++iX ) {
	 for( iZ=0; iZ<nZC; ++iZ ) {
	    cellFlags[ iX + iY * nXC + iZ * nXYC ] = MRT_wall;
	 }
      }
   }

   kOut;
   return;
}

/*******************
********************  Prepare to create a test file
*******************/

static int  prepareTest( void )
{
   static struct MfixReaderDataChoice choicesTest[2] = {
      { "EP_g", "void fraction", 0, 0, 0, 0,
	"fraction", MRD_float, MRA_cell, 1, 3, 0, 1 },
      { "Vel_g", "fluid velocity vector", 0, 0, 0, 0,
	"cm/sec", MRD_float, MRA_cell, 3, 3, 0, 1 },
   };

   int i, x0, x1, y0, y1;
   static float times[3] = { 0.0, 1.0, 2.0 };

   kIn( prepareTest );

   if( testRequested == 1 ) {
      description = "test data type 1";
      nNodes[0] = nXN = 4;
      nNodes[1] = nYN = 4;
      nNodes[2] = nZN = 4;
      nCells[0] = nXC = 3;
      nCells[1] = nYC = 3;
      nCells[2] = nZC = 3;
   }
   else if( testRequested == 2 ) {
      description = "test data type 2";
      nNodes[0] = nXN = 21;
      nNodes[1] = nYN = 13;
      nNodes[2] = nZN = 5;
      nCells[0] = nXC = 20;
      nCells[1] = nYC = 12;
      nCells[2] = nZC = 4;
   }

   nXYN = nXN * nYN;
   nXYZN = nXYN * nZN;
   nXYC = nXC * nYC;
   nXYZC = nXYC * nZC;

   month = day = year = 0;
   hours = minutes = seconds = 0;
   versionMajor = versionMinor = 0;
   name = "test";

   nChoices = 2;
   choices = choicesTest;
   choices[0].times = choices[1].times = times;

   cellFlags = (int *)malloc( nXYZC * sizeof(int) );
   for( i=0; i<nXYZC; ++i ) cellFlags[i] = MRT_fluid;

   if( testRequested == 2 ) {
      x1 = ( ( nXN - testCellSpacing ) / 2 ) - 1;
      x0 = x1 + 1 - testCellSize;
      y0 = ( nYN - testCellSize ) / 2;
      y1 = y0 + testCellSize - 1;
      addTube( x0, x1, y0, y1 );

      x0 = x1 + 1 + testCellSpacing;
      x1 = x0 + testCellSize - 1;
      addTube( x0, x1, y0, y1 );
   }

   kOut;
   return( kStatus );
}

/******************************************************************************/
/*******************                                        *******************/
/*******************              Initialize                *******************/
/*******************                                        *******************/
/******************************************************************************/

/*******************
********************  Get file information
*******************/

static int getInfo( void )
{
   int i;
   struct MfixReaderInfo *II;
   kIn( getInfo );

   if( MfixReaderGetInfo( handle, &nInfo, &info ) == 0 ) goto done;
   if( MfixReaderGetDataChoices( handle, &nChoices, &choices ) == 0 ) goto done;

   for( i=0, II=info; i<nInfo; ++i, ++II ) {
      if( strcmp( II->keyword, "dims" ) == 0 ) {
	 nNodes[0] = nXN = II->iValue[0];
	 nNodes[1] = nYN = II->iValue[1];
	 nNodes[2] = nZN = II->iValue[2];

	 nPNodes[0] = nXPN = nCells[0] = nXC = ( nXN < 2 ) ? 1 : nXN - 1;
	 nPNodes[1] = nYPN = nCells[1] = nYC = ( nYN < 2 ) ? 1 : nYN - 1;
	 nPNodes[2] = nZPN = nCells[2] = nZC = ( nZN < 2 ) ? 1 : nZN - 1;

	 nPCells[0] = nXPC = ( nXPN < 2 ) ? 1 : nXPN - 1;
	 nPCells[1] = nYPC = ( nYPN < 2 ) ? 1 : nYPN - 1;
	 nPCells[2] = nZPC = ( nZPN < 2 ) ? 1 : nZPN - 1;

	 nXYN = nXN * nYN;
	 nXYZN = nXYN * nZN;

	 nXYC = nXC * nYC;
	 nXYZC = nXYC * nZC;

	 nXYPN = nXPN * nYPN;
	 nXYZPN = nXYPN * nZPN;

	 nXYPC = nXPC * nYPC;
	 nXYZPC = nXYPC * nZPC;
      }
      else if( strcmp( II->keyword, "date" ) == 0 ) {
	 month = II->iValue[0];
	 day = II->iValue[1];
	 year = II->iValue[2];
      }
      else if( strcmp( II->keyword, "time" ) == 0 ) {
	 hours = II->iValue[0];
	 minutes = II->iValue[1];
	 seconds = II->iValue[2];
      }
      else if( strcmp( II->keyword, "version" ) == 0 ) {
	 versionMajor = II->iValue[0];
	 versionMinor = II->iValue[1];
      }
      else if( strcmp( II->keyword, "name" ) == 0 ) {
	 name = II->cValue[0];
      }
      else if( strcmp( II->keyword, "description" ) == 0 ) {
	 description = II->cValue[0];
      }
   }

done:
   kOut;
   return( kStatus );
}

/*******************
********************  Check options from command line
*******************/

static int checkCommandLine( void )
{
   int i;
   kIn( checkCommandLine );

   if( ( recl < 64 ) || ( recl > 65536 ) ) {
      kErr2( record length is %d but must be 64 to 65536, recl );
      goto done;
   }

   if( timeThreshold < 0.000001 ) {
      kBug4( 1, %s changed from %f to 0.001, argTimeThreshold, timeThreshold );
      timeThreshold = 0.001;
   }

   if( timeBounds[1] < timeBounds[0] ) {
      kBug5( 1, %s changed from %f to %f, argTimeBounds[1], timeBounds[1],
	 timeBounds[0] );
      timeBounds[1] = timeBounds[0];
   }

   if( timeIncrement < 0.000001 ) {
      kBug4( 1, %s changed from %f to 0.001, argTimeIncrement, timeIncrement );
      timeIncrement = 0.001;
   }

   for( i=0; i<6; ++i ) {
      if( skip[i] < 0 ) {
	 kBug4( 1, %s changed from %d to 0, argSkip[i], skip[i] );
	 skip[i] = 0;
      }
   }

   if( inputFileName == 0 ) {
      kErr1( failed to supply file name );
   }
   else if( strlen( inputFileName ) < 1 ) {
      kErr1( input file name is null );
   }

done:
   kOut;
   return( kStatus );
}

/*******************
********************  Add one variable to conversion list
*******************/

static int addVariable( char *name )
{
   VarRequest *v;
   kIn( addVariable );

   ++nVarRequest;
   v = varRequest;
   varRequest = (VarRequest *)malloc( nVarRequest * sizeof(VarRequest) );

   if( nVarRequest > 1 ) {
      memcpy( varRequest, v, ( nVarRequest - 1 ) * sizeof(VarRequest) );
      free(v);
   }

   v = varRequest + nVarRequest - 1;
   v->name = strdup( name );
   v->format = varFormat;

   kOut;
   return( kStatus );
}

/*******************
********************  Parse the command line
*******************/

static int readCommandLine( int nArg, char *args[] )
{
   int i, n, iArg, iLoop;
   float f;
   char *arg, *value;
   kIn( readCommandLine );

   iArg = 0;

  while( iArg < (nArg - 1) ) {
      ++iArg;
      arg = args[ iArg ];

   /****  Set flag for test data type 1 */
   /*
      if( strcmp( arg, argTest1 ) == 0 ) {
	 testRequested = 1;
	 kBug4( 2, %s just set to %d, arg, testRequested );
	 continue;
      }
   */
   /****  Set flag for test data type 2 */
   /*
      if( strcmp( arg, argTest2 ) == 0 ) {
	 if( iArg > (nArg-3) ) goto missing;
	 testRequested = 2;
	 kBug4( 2, %s just set to %d, arg, testRequested );
   */
         /****  Grab cell size */
   /*
	 ++iArg;
	 value = args[ iArg ];
	 n = sscanf( value, "%d", &i );

	 if( n != 1 ) {
	    kErr2( the value of the test2 cell size is %s; invalid, value );
	    goto done;
	 }

	 testCellSize = ( i < 1 ) ? 1 : ( ( i > 4 ) ? 4 : i ); 
   */
         /****  Grab cell spacing */
   /*
	 ++iArg;
	 value = args[ iArg ];
	 n = sscanf( value, "%d", &i );

	 if( n != 1 ) {
	    kErr2( the value of the test2 spacing is %s; invalid, value );
	    goto done;
	 }

	 testCellSpacing = ( i < 1 ) ? 1 : ( ( i > 4 ) ? 4 : i ); 
	 continue;
      }
      */

   /****  Set debug level */

      if( strcmp( arg, argBugLevel ) == 0 ) {
	 ++iArg;
	 if( iArg >= nArg ) goto missing;
	 value = args[ iArg ];
	 n = sscanf( value, "%d", &i );
	 
	 if( n != 1 ) {
	    kErr3( the value of %s is %s; invalid, arg, value );
	    goto done;
	 }
	 
	 kBugLevel = ( i < 0 ) ? 0 : ( ( i > 10 ) ? 10 : i );
	 kBug4( 2, %s just set to %d, arg, kBugLevel );
	 continue;
      }

   /****  Set input file endian */

      if( strcmp( arg, argEndian ) == 0 ) {
	 ++iArg;
	 if( iArg >= nArg ) goto missing;
	 value = args[ iArg ];
	 n = sscanf( value, "%d", &i );

	 if( n != 1 ) {
	    kErr3( the value of %s is %s; invalid, arg, value );
	    goto done;
	 }

	 endian = ( i == 0 ) ? 0 : 1;
	 kBug4( 2, %s just set to %d, arg, endian );
	 continue;
      }

   /****  Set input file record length */

      if( strcmp( arg, argRecl ) == 0 ) {
	 ++iArg;
	 if( iArg >= nArg ) goto missing;
	 value = args[ iArg ];
	 n = sscanf( value, "%d", &i );

	 if( n != 1 ) {
	    kErr3( the value of %s is %s; invalid, arg, value );
	    goto done;
	 }

	 recl = i;
	 kBug4( 2, %s just set to %d, arg, recl );
	 continue;
      }

   /****  Set time tolerance */

      if( strcmp( arg, argTimeTolerance ) == 0 ) {
	 ++iArg;
	 if( iArg >= nArg ) goto missing;
	 value = args[ iArg ];
	 n = sscanf( value, "%f", &f );

	 if( n != 1 ) {
	    kErr3( the value of %s is %s; invalid, arg, value );
	    goto done;
	 }

	 timeTolerance = f;
	 kBug4( 2, %s just set to %f, arg, timeThreshold );
	 continue;
      }

   /****  Set time threshold */

      if( strcmp( arg, argTimeThreshold ) == 0 ) {
	 ++iArg;
	 if( iArg >= nArg ) goto missing;
	 value = args[ iArg ];
	 n = sscanf( value, "%f", &f );

	 if( n != 1 ) {
	    kErr3( the value of %s is %s; invalid, arg, value );
	    goto done;
	 }

	 timeThreshold = f;
	 kBug4( 2, %s just set to %f, arg, timeThreshold );
	 continue;
      }

   /****  Set interpolation rule */

      if( strcmp( arg, argInterpRule ) == 0 ) {
	 ++iArg;
	 if( iArg >= nArg ) goto missing;
	 value = args[ iArg ];
	 n = sscanf( value, "%d", &i );
	 
	 if( n != 1 ) {
	    kErr3( the value of %s is %s; invalid, arg, value );
	    goto done;
	 }
	 
	 kBug4( 2, %s just set to %d, arg, i );

	 switch(i) {
	    case 4: interpRule = MRI_linear; break;
	    case 3: interpRule = MRI_nearest_next; break;
	    case 2: interpRule = MRI_nearest_prev; break;
	    case 1: interpRule = MRI_nearest_any; break;
	    default: interpRule = MRI_never;
	 }

	 continue;
      }

   /****  Set extrapolation rule */

      if( strcmp( arg, argExtrapRule ) == 0 ) {
	 ++iArg;
	 if( iArg >= nArg ) goto missing;
	 value = args[ iArg ];
	 n = sscanf( value, "%d", &i );
	 
	 if( n != 1 ) {
	    kErr3( the value of %s is %s; invalid, arg, value );
	    goto done;
	 }
	 
	 kBug4( 2, %s just set to %d, arg, i );

	 switch(i) {
	    case 1: extrapRule = MRE_nearest; break;
	    default: extrapRule = MRE_never;
	 }

	 continue;
      }

   /****  Set method for selecting time */

      if( strcmp( arg, argTimeMethod ) == 0 ) {
	 ++iArg;
	 if( iArg >= nArg ) goto missing;
	 value = args[ iArg ];
	 n = sscanf( value, "%d", &i );
	 
	 if( n != 1 ) {
	    kErr3( the value of %s is %s; invalid, arg, value );
	    goto done;
	 }
	 
	 kBug4( 2, %s just set to %d, arg, i );

	 switch(i) {
	    case 2: timeMethod = TM_STEP; break;
	    case 1: timeMethod = TM_ALL_BOUNDED; break;
	    default: timeMethod = TM_ALL;
	 }

	 continue;
      }

   /****  Set time bounds */

      for( iLoop=0; iLoop<2; ++iLoop ) {
	 if( strcmp( arg, argTimeBounds[ iLoop ] ) != 0 ) continue;
	 ++iArg;
	 if( iArg >= nArg ) goto missing;
	 value = args[ iArg ];
	 n = sscanf( value, "%f", &f );

	 if( n != 1 ) {
	    kErr3( the value of %s is %s; invalid, arg, value );
	    goto done;
	 }

	 timeBounds[ iLoop ] = f;
	 kBug4( 2, %s just set to %f, arg, timeBounds[ iLoop ] );
	 break;
      }

      if( iLoop<2 ) continue;

   /****  Set time increment */

      if( strcmp( arg, argTimeIncrement ) == 0 ) {
	 ++iArg;
	 if( iArg >= nArg ) goto missing;
	 value = args[ iArg ];
	 n = sscanf( value, "%f", &f );

	 if( n != 1 ) {
	    kErr3( the value of %s is %s; invalid, arg, value );
	    goto done;
	 }

	 timeIncrement = f;
	 kBug4( 2, %s just set to %f, arg, timeIncrement );
	 continue;
      }

   /****  Set subset size */

      for( iLoop=0; iLoop<6; ++iLoop ) {
	 if( strcmp( arg, argSkip[ iLoop ] ) != 0 ) continue;
	 ++iArg;
	 if( iArg >= nArg ) goto missing;
	 value = args[ iArg ];
	 n = sscanf( value, "%d", &i );

	 if( n != 1 ) {
	    kErr3( the value of %s is %s; invalid, arg, value );
	    goto done;
	 }

	 skip[ iLoop ] = i;
	 kBug4( 2, %s just set to %d, arg, skip[ iLoop ] );
	 break;
      }

      if( iLoop<6 ) continue;

   /****  Set data format */

      if( strcmp( arg, argNodeAverage ) == 0 ) {
	 varFormat = MRF_node_average;
	 kBug4( 2, %s just set to %s, arg, MfixReaderFormatText[1] );
	 continue;
      }

      if( strcmp( arg, argNodeCenter ) == 0 ) {
	 varFormat = MRF_node_center;
	 kBug4( 2, %s just set to %s, arg, MfixReaderFormatText[2] );
	 continue;
      }

   /****  Set information flag */

      if( strcmp( arg, argInfo ) == 0 ) {
	 infoRequested = 1;
	 kBug4( 2, %s just set to %d, arg, infoRequested );
	 continue;
      }

   /****  Set overwrite flag */

      if( strcmp( arg, argOverwrite ) == 0 ) {
	 overwrite = 1;
	 kBug4( 2, %s just set to %d, arg, overwrite );
	 continue;
      }

   /****  Select variables */

      if( strcmp( arg, argVariable ) == 0 ) {
	 ++iArg;
	 if( iArg >= nArg ) goto missing;
	 addVariable( args[ iArg ] );
	 kBug4( 2, %s just added request for %s, arg, args[ iArg ] );
	 continue;
      }

   /****  Provide help */

      for( iLoop=0; iLoop<4; ++iLoop ) {
	 if( strcmp( arg, argHelp[ iLoop ] ) != 0 ) continue;
	 showHelp();
	 goto done;
      }

    /****  Report bad command line flag */

      if( ( arg[0] == '-' ) || ( iArg < (nArg - 1) ) ) {
	 kErr2( unknown option %s; use -h for help, arg );
	 goto done;
      }

   /****  Set file name */

      if( iArg == (nArg - 1) ) {
	 inputFileName = arg;
	 kBug3( 2, filename set to %s, inputFileName );
	 continue;
      }
   }

   goto done;

missing:
   kErr2( the value for option %s is missing, arg );

done:
   kOut;
   return( kStatus );
}

/*******************
********************  Main
*******************/

int main( int nArg, char *args[] )
{
   kMain( M2E );
   currentTime = time(0);

   if( !readCommandLine( nArg, args ) ) goto done;
   if( !checkCommandLine() ) goto done;

   if( testRequested ) {
      prepareTest();
   }
   else {
      MfixReaderCreate( &handle );
      if( !kStatus ) goto destroyit;

      MfixReaderSetRecordLength( handle, recl );
      if( !kStatus ) goto destroyit;

      MfixReaderSetInputEndian( handle, endian );

      MfixReaderSetDataTimeTolerance( handle, timeTolerance );
      if( !kStatus ) goto destroyit;

      MfixReaderSetTimeRules( handle, timeThreshold, extrapRule, interpRule );
      if( !kStatus ) goto destroyit;

      MfixReaderSetFileName( handle, inputFileName );
      if( !kStatus ) goto destroyit;

      MfixReaderSetLocationFilter( handle, skip );
      if( !kStatus ) goto destroyit;

      MfixReaderOpen( handle );
      if( !kStatus ) goto destroyit;

      if( !getInfo() ) goto closeit;

      if( infoRequested ) {
	 printf( "General information from the data set:\n" );
	 showInfo( "   ", stdout );
	 goto closeit;
      }
   }

   if( !selectVariables() ) goto closeit;
   if( !selectTimes() ) goto closeit;
   if( !prepareNames() ) goto closeit;
   if( !prepareCoordinates() ) goto closeit;
   if( !prepareCelltypes() ) goto closeit;

   cell2geom( nCells, cellFlags, nDrawTypes, drawTypes, &drawings );
   if( !kStatus ) goto closeit;

   if( !writeGeo() ) goto closeit;
   if( !writeData() ) goto closeit;
   if( !writeCase() ) goto closeit;

closeit:
   MfixReaderClose( handle );

destroyit:
   MfixReaderDestroy( &handle );

done:
   kOut;
   return(0);
}

/******************************  End of m2e.c  ********************************/
