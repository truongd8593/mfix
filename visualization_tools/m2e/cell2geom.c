/******************************************************************************/
/*******************                                        *******************/
/*******************    Create drawings from cell types     *******************/
/*******************                                        *******************/
/******************************************************************************/
/*
Purpose
   Create drawings from cell types defined over a structured grid.
Original
   Author: Kent Eschenberg eschenbe@psc.edu
   Organization: Pittsburgh Supercomputing Center, CMU
   Date: 2004/03/0110:26pm
Current
   $Source$
   $Author$
   $Revision$
   $Date$
Routines
   Pick Cells
      plot .................... simple character plot used to debug cell picking
      keepZLines .............. keep Z lines
      keepXLines .............. keep X lines
      applyRange .............. pick all cells that match type range
      pickCells ............... pick cells
   Build List of Drawings
      addDrawing .............. add one drawing to list
      addQuad ................. add one quad to list
      makeOneDrawing .......... make one drawing
   External Routines
      cell2geomFree ........... free list of drawings
      cell2geom ............... create drawings
Input
   The grid size is defined using an array of 3 integers, none less than 1,
   each giving the number of nodes in X, Y or Z. The cell types are given
   as a 3D array of integers. The cell types for which drawings are desired
   are given as an array of integer pairs where all cells with values
   within the pair, inclusive, are part of that cell type.
Output
   For each cell type, a drawing will be returned if there are one or more
   cells of that type. Each drawing is an unordered array of rectangles
   representing the sides of the cells of the given type with interior sides
   excluded. Each rectangle is expressed as an array of four integers
   which are origin-one indices into the array of vertices. The integers
   are ordered so the the implicit normal points out of the cell.
*/

#include "kLib.h"
#include "cell2geom.h"
#include <unistd.h>

/*****************************  Local variables  ******************************/

static char *cell2geom_Id
   = "@(#) $Id$";

static int nXNodes, nYNodes, nZNodes, nXCells, nYCells, nZCells;
static int eXCells, eYCells, eZCell;
static int nXYNodes, nXYCells, nXYZCells;

static int *cells;       /* Application's cell flags */
static int *cellsBool;   /* Boolean flags (1=in current cell range) */
static int nCells;       /* Number of True flags in cellsBool */

static int *elements = 0;            /* Quad list for drawing being constructed */
static int nElements, nElementsAlloc;   /*  " number, number allocated */

static Drawing *drawings;   /* List of drawings being built */
static int nDrawings;       /*  " number */

/******************************************************************************/
/*******************                                        *******************/
/*******************              Pick cells                *******************/
/*******************                                        *******************/
/******************************************************************************/
/*
static int plot( char *msg, int iY )
{
   int a, iX, jX, jY, iZ, jZ;
   char buf[128];
   kIn( plot );

   if( ( iY != 92 ) && ( iY != 103 ) ) goto done;

   printf( "%s\n", msg );
   jY = iY * nXCells;

   for( iZ=0, jZ=jY; iZ<nZCells; ++iZ, jZ+=nXYCells ) {
      for( iX=0, jX=jZ; iX<nXCells; ++iX, ++jX ) {
	 a = cellsBool[ jX ];

	 if( !a )
	    buf[ iX ] = '-';
	 else
	    sprintf( buf + iX, "%X", a );
      }

      buf[ iX ] = 0;
      printf( "%s\n", buf );
   }

   sleep(1);

done:
   kOut;
   return( kStatus );
}
*/
/*******************
********************
*******************/

static int keepZLines( int jY )
{
   int iX, jX, kZ, lZ, iZ, jZ, iX0, jX0, iX1, jX1;
   kIn( keepZLines );

   /* At each Z find largest span in X bounded by non-null cells which */
   /* have at least one null cell between them and the outside of the grid. */

   for( iZ=0, jZ=jY; iZ<nZCells; ++iZ, jZ+=nXYCells ) {

   /****  From start of X axis search forward for first null cell */

      for( iX=0, jX=jZ; iX<nXCells; ++iX, ++jX ) {
	 if( !cellsBool[ jX ] ) break;
      }

      if( iX >= nXCells ) continue;
      iX0 = iX;
      jX0 = jX;

   /****  Continue until find first cell with count of 1 */

      for( iX=iX0+1, jX=jX0+1; iX<nXCells; ++iX, ++jX ) {
	 if( cellsBool[ jX ] == 1 ) break;
      }

      if( iX >= nXCells ) continue;

   /****  This is the start of the desired span in X */

      iX0 = iX;
      jX0 = jX;

   /****  From other end of X axis search backwards for first null cell */

      for( iX=eXCells, jX=jZ+eXCells; iX>iX0; --iX, --jX ) {
	 if( !cellsBool[ jX ] ) break;
      }

      if( iX <= iX0 ) continue;
      iX1 = iX;
      jX1 = jX;

   /****  Continue backwards until find first cell at 1 */
   /****  This search cannot fail but it may end when iX=iX0 */

      for( iX=iX1, jX=jX1; iX>iX0; --iX, --jX ) {
	 if( cellsBool[ jX ] == 1 ) break;
      }

   /****  This is the end of the desired X span */

      iX1 = iX;

   /****  For each X in (iX0,iX1) increase longevity of cells at 1 */
   /****  and strings of their similar Z neighbors */

      for( iX=iX0, jX=jX0; iX<=iX1; ++iX, ++jX ) {
	 if( cellsBool[ jX ] != 1 ) continue;

      /****  From [jX] work towards +Z end */

	 for( kZ=iZ, lZ=jX; kZ<nZCells; ++kZ, lZ+=nXYCells ) {
	    if( !cellsBool[ lZ ] ) break;
	    ++cellsBool[ lZ ];
	 }

      /****  From [jX] work towards -Z end */

	 for( kZ=iZ-1, lZ=jX-nXYCells; kZ>=0; --kZ, lZ-=nXYCells ) {
	    if( !cellsBool[ lZ ] ) break;
	    ++cellsBool[ lZ ];
	 }
      }
   }

   kOut;
   return( kStatus );
}

/*******************
********************
*******************/

static int keepXLines( int jY )
{
   int iX, jX, kX, lX, iZ, jZ, iZ0, jZ0, iZ1, jZ1;
   kIn( keepXLines );

   /* At each X find largest span in Z bounded by non-null cells which */
   /* have at least one null cell between them and the outside of the grid. */

   for( iX=0, jX=jY; iX<nXCells; ++iX, ++jX ) {

   /****  From start of Z axis search forward for first null cell */

      for( iZ=0, jZ=jX; iZ<nZCells; ++iZ, jZ+=nXYCells ) {
	 if( !cellsBool[ jZ ] ) break;
      }

      if( iZ >= nZCells ) continue;
      iZ0 = iZ;
      jZ0 = jZ;

   /****  Continue until find first cell with count of 1 */

      for( iZ=iZ0+1, jZ=jZ0+nXYCells; iZ<nZCells; ++iZ, jZ+=nXYCells ) {
	 if( cellsBool[ jZ ] == 1 ) break;
      }

      if( iZ >= nZCells ) continue;

   /****  This is the start of the desired span in Z */

      iZ0 = iZ;
      jZ0 = jZ;

   /****  From end of Z axis search backwards for first null cell */

      for( iZ=eZCell, jZ=jX+eZCell*nXYCells; iZ>iZ0; --iZ, jZ-=nXYCells ) {
	 if( !cellsBool[ jZ ] ) break;
      }

      if( iZ <= iZ0 ) continue;
      iZ1 = iZ;
      jZ1 = jZ;

   /****  Continue backwards until find first cell at 1 */
   /****  This search cannot fail but it may end when iZ=iZ0 */

      for( iZ=iZ1, jZ=jZ1; iZ>iZ0; --iZ, jZ-=nXYCells ) {
	 if( cellsBool[ jZ ] == 1 ) break;
      }

   /****  This is the end of the desired Z span */

      iZ1 = iZ;

   /****  For each Z in (iZ0,iZ1) increase longevity of cells at 1 */
   /****  and strings of their similar X neighbors */

      for( iZ=iZ0, jZ=jZ0; iZ<=iZ1; ++iZ, jZ+=nXYCells ) {
	 if( cellsBool[ jZ ] != 1 ) continue;

      /****  From [jZ] work towards +X end */

	 for( kX=iX, lX=jZ; kX<nXCells; ++kX, ++lX ) {
	    if( !cellsBool[ lX ] ) break;
	    ++cellsBool[ lX ];
	 }

      /****  From [jZ] work towards -X end */

	 for( kX=iX-1, lX=jZ-1; kX>=0; --kX, --lX ) {
	    if( !cellsBool[ lX ] ) break;
	    ++cellsBool[ lX ];
	 }
      }
   }

   kOut;
   return( kStatus );
}

/*******************
********************
*******************/

static int applyRange( int t0, int t1 )
{
   int i, *c, *cb;
   kIn( pickCellsRange );

   nCells = 0;
   c = cells;
   cb = cellsBool;
   memset( cb, 0, nXYZCells * sizeof(int) );

   for( i=0; i<nXYZCells; ++i, ++c, ++cb ) {
      if( ( *c < t0 ) || ( *c > t1 ) ) continue;
      *cb = 1;
      ++nCells;
   }

   kOut;
   return( kStatus );
}

/*******************
********************
*******************/

static int pickCells( DrawRule rule, int t0, int t1 )
{
   int iX, jX, iY, jY, iZ, jZ;
   /* char s[80]; */
   kIn( pickCells );

   applyRange( t0, t1 );

/****  Exclude cells at outer walls and corners if requested */

   if( rule == DR_ExclWallCorner ) {
      for( iY=jY=0; iY<nYCells; ++iY, jY+=nXCells ) {
         /*
	 sprintf( s, "orig %d by %d at iY=%d", nXCells, nZCells, iY );
	 plot( s, iY );
         */

      /****  Add at least 1 to boolean flags for lines of cells in X and Z */

	 keepXLines( jY );
	 keepZLines( jY );

      /****  Remove anything not a line */

	 for( iZ=0, jZ=jY; iZ<nZCells; ++iZ, jZ+=nXYCells ) {
	    for( iX=0, jX=jZ; iX<nXCells; ++iX, ++jX ) {
	       if( cellsBool[ jX ] < 2 ) cellsBool[ jX ] = 0;
	    }
	 }

         /* plot( "done", iY ); */
      }
   }

   kOut;
   return( kStatus );
}

/******************************************************************************/
/*******************                                        *******************/
/*******************        Build list of drawings          *******************/
/*******************                                        *******************/
/******************************************************************************/

static int addDrawing( char *name, char *qName )
{
   Drawing *drawing;
   kIn( addDrawing );

   ++nDrawings;
   drawing = drawings;
   drawings = (Drawing *)malloc( nDrawings * sizeof(Drawing) );

   if( nDrawings > 1 ) {
      memcpy( drawings, drawing, ( nDrawings - 1 ) * sizeof(Drawing) );
      free( drawing );
   }

   drawing = drawings + nDrawings - 1;
   drawing->name = strdup( name );
   drawing->nCells = nCells;
   drawing->eName = strdup( qName );
   drawing->nE = nElements;
   drawing->e = elements;

   kOut;
   return( kStatus );
}

/*******************
********************
*******************/

static int addQuad( int i0, int i1, int i2, int i3 )
{
   int *e2;
   kIn( addQuad );

/****  Make room if needed */

   if( nElements >= nElementsAlloc ) {
      nElementsAlloc += 100;
      e2 = elements;
      elements = (int *)malloc( 4 * nElementsAlloc * sizeof(int) );

      if( nElements > 0 ) {
	 memcpy( elements, e2, 4 * nElements * sizeof(int) );
	 free( e2 );
      }
   }

/****  Add new quad to end of list */

   e2 = elements + 4 * nElements;
   e2[0] = i0;
   e2[1] = i1;
   e2[2] = i2;
   e2[3] = i3;
   ++nElements;

   kOut;
   return( kStatus );
}

/*******************
********************
*******************/

static int makeOneDrawing( char *name, DrawRule rule, int t0, int t1 )
{
   int i, nYZ, nZ, *cb, iX, iY, iZ;
   kIn( makeOneDrawing );
   elements = 0;
   nElements = nElementsAlloc = 0;
   pickCells( rule, t0, t1 );
   cb = cellsBool;

   for( iZ=0; iZ<nZCells; ++iZ ) {
      nZ = 1 + iZ * nXYNodes;

      for( iY=0; iY<nYCells; ++iY ) {
	 nYZ = iY * nXNodes + nZ;

	 for( iX=0; iX<nXCells; ++iX, ++cb ) {
	    if( !(*cb) ) continue;  /* Skip cells that are not part of the drawing */
	    i = iX + nYZ;

	 /****  2D Case: Always put a quad at the location of the cell */

	    if( nZCells == 1 ) {
	       addQuad( i, i+1, i+1+nXNodes, i+nXNodes );
	    }

	 /****  3D Case: Put a quad on any side that is at the edge or where the adjacent */
          /****  cell is not part of the drawing */

	    else {
               /****  Add +X side of cell */
	       if( ( iX == eXCells ) || !cb[1] )
		  addQuad( i+1, i+1+nXNodes, i+1+nXNodes+nXYNodes,i+1+nXYNodes);

               /****  Add -X side of cell */
	       if( ( iX == 0 ) || !cb[-1] )
		  addQuad( i, i+nXYNodes, i+nXNodes+nXYNodes, i+nXNodes );

               /****  Add +Y side of cell */
	       if( ( iY == eYCells ) || !cb[ nXCells ] )
		  addQuad( i+nXNodes, i+nXNodes+nXYNodes, i+1+nXNodes+nXYNodes,
		     i+1+nXNodes );

               /****  Add -Y side of cell */
	       if( ( iY == 0 ) || !cb[ -nXCells ] )
		  addQuad( i, i+1, i+1+nXYNodes, i+nXYNodes );

               /****  Add +Z side of cell */
	       if( ( iZ == eZCell ) || !cb[ nXYCells ] )
		  addQuad( i+nXYNodes, i+1+nXYNodes, i+1+nXNodes+nXYNodes,
		     i+nXNodes+nXYNodes );

               /****  Add -Z side of cell */
	       if( ( iZ == 0 ) || !cb[ -nXYCells ] )
		  addQuad( i, i+nXNodes, i+1+nXNodes, i+1 );
	    }
	 }
      }
   }

   if( nElements > 0 ) addDrawing( name, "quad4" );
   kOut;
   return( kStatus );
}

/******************************************************************************/
/*******************                                        *******************/
/*******************           External routines            *******************/
/*******************                                        *******************/
/******************************************************************************/

int cell2geomFree( Drawing *d )
{
   int i;
   kIn( cell2geomFree );

   if( d == 0 ) goto done;
   nDrawings = d->nDrawings;

   for( i=0; i<nDrawings; ++i ) {
      free( (d[i]).name );
      free( (d[i]).eName );
      free( (d[i]).e );
   }

   free(d);

done:
   kOut;
   return( kStatus );
}

/*******************
********************
*******************/

int cell2geom(
   int dim[3],     /* Number of cells in X, Y, Z */
   int *cells2,     /* Cell types at (0,0,0), (1,0,0), ... */
   int nT,         /* Number of cell types to be drawn */
   DrawType *t,    /* Cell types to be drawn */
   Drawing **d     /* Return: drawings, or null if none found */
)
{
   int i;
   kIn( cell2geom );

   if( cell2geom_Id ) nXCells = 1; /* Done only to force var into object code */
   nDrawings = 0;
   drawings = 0;

/****  Move inputs to local variables */

   cells = cells2;

   nXCells = dim[0];
   nYCells = dim[1];
   nZCells = dim[2];
   eXCells = nXCells - 1;
   eYCells = nYCells - 1;
   eZCell = nZCells - 1;
   nXYCells = nXCells * nYCells;
   nXYZCells = nXYCells * nZCells;

   nXNodes = nXCells + 1;
   nYNodes = nYCells + 1;
   nZNodes = nZCells + 1;
   nXYNodes = nXNodes * nYNodes;

   cellsBool = (int *)malloc( nXYZCells * sizeof(int) );

/****  Create drawing for each cell type range */

   for( i=0; i<nT; ++i, ++t ) {
      if( !makeOneDrawing( t->name, t->rule, t->t0, t->t1 ) ) break;
   }

/*****  Return possibly empty list of drawings */

   if( drawings ) {
      for( i=0; i<nDrawings; ++i ) drawings[i].nDrawings = nDrawings;
   }

   *d = drawings;
   if( cellsBool != 0 ) free( cellsBool );
   kOut;
   return( kStatus );
}

/****************************  End of cell2geom.c  ****************************/
