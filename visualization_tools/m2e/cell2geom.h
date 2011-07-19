/******************************************************************************/
/*******************                                        *******************/
/*******************      Header file for cell2geom.c       *******************/
/*******************                                        *******************/
/******************************************************************************/
/*
Purpose
   Please see cell2geom.c.
Original
   Author: Kent Eschenberg eschenbe@psc.edu
   Organization: Pittsburgh Supercomputing Center, CMU
   Date: 2004/03/0110:26pm
Current
   $Source$
   $Author$
   $Revision$
   $Date$
*/

/**************************  Cell types to be drawn  **************************/

enum DrawRuleEnum { DR_None, DR_ExclWallCorner };
typedef enum DrawRuleEnum DrawRule;

struct DrawTypeStr {
   char *name;      /* Short name for this type */
   DrawRule rule;   /* Rule used to select cells */
   int t0;          /* First integer representing this type, inclusive */
   int t1;          /* Last, inclusive */
};

typedef struct DrawTypeStr DrawType;

/*******************************  One drawing  **************************/

struct DrawingStr {
   int nDrawings;   /* Overall number of drawings */
   char *name;      /* Name */
   int nCells;      /* Number of cells used to form drawing */
   char *eName;     /* Element type (currently only "quad4")  */
   int nE;          /* Number of elements */
   int *e;          /* Arrays of origin-one indices into coordinates */
};

typedef struct DrawingStr Drawing;

/*************************  External entry to package  ************************/

int cell2geomFree( Drawing *d );

int cell2geom(
   int dim[3],      /* Number of cells in X, Y, Z */
   int *cells,      /* Cell types at (0,0,0), (1,0,0), ... */
   int nT,          /* Number of cell types to be drawn */
   DrawType *t,    /* Cell types to be drawn */
   Drawing **d     /* Return: drawings */
);

/****************************  End of cell2geom.h  ****************************/
