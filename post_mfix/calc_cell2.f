CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: CALC_CELL2 (REACTOR_LOC,D_DIR,N_DIR,CELL_LOC)          C
C  Purpose: calculate the cell index for a reactor location            C
C                                                                      C
C  Author: P. Nicoletti                               Date: 02-DEC-91  C
C  Reviewer: M.SYAMLAL, W.ROGERS, P.NICOLETTI         Date: 24-JAN-92  C
C                                                                      C
C  Revision Number:                                                    C
C  Purpose:                                                            C
C  Author:                                            Date: dd-mmm-yy  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Literature/Document References:                                     C
C                                                                      C
C  Variables referenced:  None                                         C
C  Variables modified: None                                            C
C                                                                      C
C  Local variables: LC, CELL_START, CELL_END                           C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
C
      SUBROUTINE CALC_CELL2 (REACTOR_LOC,D_DIR,N_DIR,CELL_LOC)
C
      IMPLICIT NONE
      INCLUDE 'param.inc'
C
C passed arguments
C
C    REACTOR_LOC - location along one of the axis for which the cell
C                  index is to be found
C    D_DIR       - cell lengths (DX,DY,DZ)
C    N_DIR       - number of cells in this direction (IMAX2,JMAX2,KMAX2)
C    CELL_LOC    - cell index corresponding to REACTOR_LOC
C
C local variables
C    LC         -  loop counter
C    CELL_START -  start coordinate for cell
C    CELL_END   -  end   coordinate for cell
C
      INTEGER           N_DIR , CELL_LOC , LC
      DOUBLE PRECISION  REACTOR_LOC , D_DIR(*) , CELL_START , CELL_END
C
      CELL_LOC = - 1
      CELL_START = 0.0
      DO 100 LC = 2,N_DIR+1
         CELL_END   = CELL_START + D_DIR(LC)
         IF (REACTOR_LOC .GE.CELL_START .AND. 
     &                           REACTOR_LOC .LE. CELL_END)  THEN
            CELL_LOC = LC
            RETURN
         END IF
      CELL_START = CELL_END
100   CONTINUE
C
      RETURN
      END
