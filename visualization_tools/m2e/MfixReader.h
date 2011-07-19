/******************************************************************************/
/*******************                                        *******************/
/*******************       Routines to Read MFIX Files      *******************/
/*******************              Declarations              *******************/
/*******************                                        *******************/
/******************************************************************************/
/*
Purpose
   Read a set of files created by the MFIX simulation. If a parallel version
   of MFIX produced seperate files for each process, they must be combined
   before this routine can read them. The restart file (extension .res) is
   read along with one or more of the data files (extension .sp1, ..., sp9).
   The data files include the data for all time steps. However, the
   routine for reading the data files (MfixReaderGetData) reads only one time
   step each time it is called. 
Original
   Author: Kent Eschenberg eschenbe@psc.edu
   Organization: Pittsburgh Supercomputing Center, CMU
   Date: 2001/11/6 08:00:00
Current
   $Source$
   $Author$
   $Revision$
   $Date$
Usage
   You must first call "MfixReaderCreate". It allocates some memory to keep
   track of what your are doing and returns an integer "handle" that is
   used to refer to this file set. A "file set" includes the restart
   file plus all of the associated time-varying data files. You can
   have as many file sets open at once as you wish, limited only by
   your computer's memory and the number of open files allowed by your
   operating system.

   You should then call one or more of the "MfixReaderSet*" routines that
   control how the package operates. For example, you must call
   "MfixReaderSetFileName" to set the base name of the files. If you call
   "MfixReaderSetLocationFilter" to select a subset of the data, only the
   subset will be reported in later requests for sizes, coordinates and
   data.

   Not much happens until you call "MfixReaderOpen". This reads the restart
   file and opens all of the other files.

   You will probably wish to call one or more of the "MfixReaderGet*"
   routines to obtain the size of the data array, the number of times, and
   the coordinates.

   You are now ready to call "MfixReaderGetData" as many times as you wish
   where each call returns data at one time step.

   Finally, you should call "MfixReaderClose" to close the files, and
   "MfixReaderDestroy" to destroy the reader instance and release memory.
Memory
   Coordinate and data arrays will be allocated in this package. Calls to
   obtain this information will return pointers to these internal arrays.
   DO NOT FREE THE ARRAYS AT THESE POINTERS! The arrays will be freed when
   the instance is closed or destroyed. You may change the values in the
   data arrays since they will be reset at the next request for data.
   DO NOT CHANGE OTHER ARRAYS!
Compatibility
   This package should work on any system that stores integers, floats
   and doubles in 4, 4 and 8 bytes using the same format that was used
   to write the file. In other words, no attempt is made to convert
   from the file's binary format to the format of the machine running
   this code.
*/

#ifndef MfixReaderIncluded
#define MfixReaderIncluded 1

/*
Large File Support
   May work only with the GNU compilers and glibc.
   The following flag redefines some I/O features to use 64-bits.
   This code has been changed to also explicitly use fopen64 and fseeko64.
   As a result 32-bit systems to read files over 2GB.
   THE FOLLOWING FLAG MUST COME BEFORE INCLUDING OTHER FILES!
*/
#define _LARGEFILE64_SOURCE 1

/******************************************************************************/
/*******************                                        *******************/
/*******************      Enumerations and structures       *******************/
/*******************                                        *******************/
/******************************************************************************/

/****
 ****  Enumerations used to describe data and make requests
 ****
 ****/

enum MfixReaderCoord { MRC_unknown, MRC_cartesian, MRC_cylindrical };
extern char *MfixReaderCoordText[3];

enum MfixReaderCellType { 
   MRT_min=1,
   MRT_fluid=1,
   MRT_wall=100,
   MRT_inflow1=10,
   MRT_inflow2=20,
   MRT_outflow1=11,
   MRT_outflow2=21,
   MRT_max=200
};

enum MfixReaderDataType { MRD_integer, MRD_float, MRD_double, MRD_string };
extern char *MfixReaderDataTypeText[3];

enum MfixReaderAttach { MRA_cell, MRA_node };
extern char *MfixReaderAttachText[2];

enum MfixReaderFormat { MRF_none, MRF_cell, MRF_node_average, MRF_node_center };
extern char *MfixReaderFormatText[4];

/* Indices into array of data choices returned by *GetDataChoices */
enum MfixReaderVariables { MRV_celltype, MRV_EP_g, MRV_P_g, MRV_P_star, MRV_U_g,
   MRV_V_g,  MRV_W_g,  MRV_U_s,  MRV_V_s,  MRV_W_s,  MRV_ROP_s,  MRV_T_g,
   MRV_T_s,  MRV_X_g,  MRV_THETA_m,  MRV_scalar };

enum MfixReaderInterpolation { MRI_never, MRI_nearest_any, MRI_nearest_prev,
   MRI_nearest_next, MRI_linear };
extern char *MfixReaderInterpolationText[5];

enum MfixReaderExtrapolation { MRE_never, MRE_nearest };
extern char *MfixReaderExtrapolationText[2];

/****
 ****  A list of these structures is used to return information about the files
 ****
 ****/

struct MfixReaderInfo {
   char *keyword;
   enum MfixReaderDataType dataType;
   int nValues;                        /* Number of values (1,2,3,...) */
   int *iValue;                        /* Value(s) if integer else null */
   float *fValue;                      /* Value(s) if float else null */
   char **cValue;                      /* Value(s) if string else null */
};

/****
 ****  A list of these structures is used to return information about the data
 ****  variables in the files.
 ****
 ****  WARNING: The pointer "times" in the following structure is a copy of a
 ****  pointer stored elsewhere. Do NOT free the copy below.
 ****/

struct MfixReaderDataChoice {
   char *variableName;      /* Name of data variable */
   char *variableDesc;      /* Description of variable */
   char *primaryName;       /* Name of primary variant (null if not used) */
   int primaryCount;        /* Number of cases for primary variant */
   char *secondaryName;     /* Name of secondary variant (null if not used) */
   int *secondaryCount;     /* Num cases for secondary variant per primary */
   char *units;
   enum MfixReaderDataType dataType;
   enum MfixReaderAttach attach;
   int veclen;              /* Length of vector (1 if a scalar) */
   int nTimes;              /* Number of times found in file */
   float *times;            /* Times found in file (seconds) */
   int available;           /* Current status: 1=available, 0=not */
};

/****
 ****  A list of these structures are used to pass request(s) for data to
 ****  the library, and to return the data.  If data is not found, the value
 ****  returned for "dims" is null. When data is found, only one of the two
 ****  data pointers ("dataInt" or "dataFloat") is set; the other is null.
 ****  Component is is ignored for scalars; for vectors, a 0/1/2 results in
 ****  only the X,Y or Z component being returned and veclen being set to 1.
 ****  Any other value causes all the components to be returned with veclen
 ****  set to 3.
 ****/

struct MfixReaderDataRequest {
   char *variableName;
   int primaryNum;
   int secondaryNum;
   int component;                    /* If vector, comp 0,1,2; -1 for all */
   enum MfixReaderFormat format;
   char *units;                      /* Returned: units */
   enum MfixReaderDataType dataType; /* Returned: data type */
   enum MfixReaderAttach attach;     /* Returned: data attachment */
   int veclen;                       /* Returned: vector length */
   int interpolated;                 /* Returned: 1 if interpolated */
   float dataTime[2];                /* Returned: data times used for interp */
   int *dims;                        /* Returned: size of array (x,y,z) */
   int *dataInt;                     /* Returned: data if integer else null */
   float *dataFloat;                 /* Returned: data if float else null */
};

/******************************************************************************/
/*******************                                        *******************/
/*******************           Function prototypes          *******************/
/*******************                                        *******************/
/******************************************************************************/

/*********************  Create a new file reader instance  ********************/
/*
Purpose
   Create an instance of a file set reader.
*/

int MfixReaderCreate( int *handle );

/*********************  Destroy a file reader instance  ***********************/
/*
Purpose
   Destroy an instance of a file set reader. Calls MfixReaderClose if needed.
*/

int MfixReaderDestroy( int *handle );

/****************************  Open a file set  *******************************/
/*
Purpose
   Open a file set. The filename must be set. Optional settings such as
   the location filter, record length and file protocol must be set before
   this routine is called if they are to have an effect.
*/

int MfixReaderOpen( int handle );

/***************************  Close a file set  *******************************/
/*
Purpose
   Close a file set. The base file name and location filter are retained.
   Frees most arrays.
*/

int MfixReaderClose( int handle );

/**************************  Set the input file name  *************************/
/*
Purpose
   Set the "base" or "root" name for the file set, Any extensions in the
   name provided (such as ".xxx") will be removed, Extensions will be added
   as needed: ".res", ".sp1", etc.
*/

int MfixReaderSetFileName( int handle, char *name );

/********************  Set the input file record length  **********************/
/*
Purpose
   Set the number of bytes in a record. The default (i.e., if this routine
   is not used) is 512.
*/

int MfixReaderSetRecordLength( int handle, int length );

/***********************  Set the input file format  **************************/
/*
Purpose
   Set whether the input file stores binary numbers small-endian (endian
   equals 0) or big-endian (endian equals anything other than 0). The default
   (i.e., if this routine is not called) is big-endian.
*/

int MfixReaderSetInputEndian( int handle, int endian );

/*******************  Set the data time tolerance  ****************************/
/*
Purpose
   Set the tolerance (seconds) used to determine if an input file is corrupt.
   If a given time in a file is a smaller number (i.e., earlier) than the
   previous time in the file by an amount greater than this tolerance, the
   file is considered corrupt and the open of the file set will fail. Defaults
   to 0.000001.
*/

int MfixReaderSetDataTimeTolerance( int handle, float tol );

/***********************  Set the location filter  ****************************/
/*
Purpose
   Set the location filter. The six integers provided indicate how many
   cells to remove at the ends of each axis. The first pair of integers are
   applied to the beginning and end of the X (or radius) axis; the second
   pair to the Y (or Z) axis; and the third to the Z (or angular) axis.
   The trailing end of the axis will be adjusted, if needed, to ensure that
   at least one value is retained along that axis. For example, a filter
   of {2,2,2,2,2,2} means to skip the two cells at the ends of each axis.
   The default is 0,0,0,0,0,0}.
*/

int MfixReaderSetLocationFilter( int handle, int skip[6] );

/*************************  Set the time selection rules  *********************/
/*
Purpose
   Set the rules for selecting the data time(s) for a variable. Data times
   are the times for which the variable was saved in the MFIX files. The
   selection of one or two data times is performed when responding to a
   request for a particular variable at a particular time. These values
   can be changed at any time.
Extrapolation
   This rule is applied when the requested time is outside of the span
   of time covered by the data times. A value of MRE_nearest means to pick
   the nearest data time, i.e., the first or the last. A value of
   MRE_never means to signal "no data available".
Interpolation
   This rule is applied when the requested time is within the span of
   data times but is further than "tolerance" from any of the data times.
   A value of MRI_never means to signal "no data available". A value of
   MRI_nearest_any means to use the nearest data time, whether it
   is before or after the requested time. Values of MRI_nearest_prev or
   MRI_nearest_next use the nearest previous, or the nearest subsequent,
   data time. A value of MRI_linear means to use linear interpolation to
   estimate the variable using two data times on either side of the
   requested time.
*/

int MfixReaderSetTimeRules(
   int handle,
   float tolerance,
   enum MfixReaderExtrapolation extrapRule,
   enum MfixReaderInterpolation interpRule
);

/****************************  Get information  *******************************/
/*
Purpose
   Obtain information. Note that, for any given keyword, either integer(s)
   or string(s) are used. Here is a partial list of keywords:
      version - major, minor version
      date - month, day, year
      dims - number of points along X, Y, Z
      dimsoriginal - number of points before the locaton filter was applied
      ntimes - number of time steps
      timespan - earliest and latest time in seconds
      ngasspecies - called nmax(0) in MFIX
      nsolidphases - called mmax in MFIX
      nsolidspecies - an array of length nsolidphases; called nmax() in MFIX
*/

int MfixReaderGetInfo(
   int handle,
   int *nInfoArg,                /* Number of information items */
   struct MfixReaderInfo **info  /* Array of information items */
);

/**************************  Get coordinates  *********************************/
/*
Purpose
   Get the coordinates of the nodes in the cartesian system. When the file uses
   the cylindrical system, it is automatically converted to cartesion. The
   dimensions are after the location filter is applied. When the requested
   format is "MRF_node_center", the coordinates are for the cell centers and the
   dimensions are smaller by one. The coordinates are for the nodes for all
   other values for format.
*/

int MfixReaderGetCoordinates(
   int handle,
   enum MfixReaderFormat format,
   char **units,
   int *dims[3],
   float **xyz  /* Order: (xyz),(xyz),... from min to max, X changing fastest */
);

/************************  Get list of available data  ************************/
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
   secondary. The number of secondary varients is calculated by taking the
   maximum of all the secondary varients. For a particular primary varient,
   the number of secondary varients may be smaller than the maximum.
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
);

/*********************************  Get data  *********************************/
/*
Purpose
   Get one or more arrays of data at a single time step. The location filter,
   if set before the file set was opened, is applied before the arrays are
   returned.
Request Grouping
   The program will run more efficiently if you combine requests so that
   values stored in the same file are requested at the same time. For
   example, obtain U_g,V_g,W_g in the same request rather than in three
   seperate requests. The order does not matter - requesting W_g,V_g,U_g
   is just as efficient as requesting them in the opposite order.
Note
   MFIX data is always given at the cells.
Interpolation
   The current interpolation settings are used to decide whether
   interpolation is used and, if so, which style is used. If interpolation
   is not used, both values of "dataTime" in the request record are set
   to the time that was used.
Return
   The return status is set to 1 (an error) only if an error occurs; it is
   not set to 1 merely because the requested data does not exist. Check the
   returned data pointers in each request to discover what was found. 
*/

int MfixReaderGetData(
   int handle,
   float dataTime,   /* Desired time (seconds) */
   int nRequests,
   struct MfixReaderDataRequest *requests
);

#endif
/****  #ifndef MfixReaderIncluded  ****/

/***************************  End of MfixReader.h  ****************************/
