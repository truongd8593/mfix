CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: CALC_QUANTITIES                                        C
C  Purpose: Calculate various quantities                               C
C                                                                      C
C  Author: P. Nicoletti                               Date: dd-mmm-yy  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Revision Number:                                                    C
C  Purpose:                                                            C
C  Author:                                            Date: dd-mmm-yy  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Literature/Document References:                                     C
C                                                                      C
C  Variables referenced:                                               C
C  Variables modified:                                                 C
C                                                                      C
C  Local variables:                                                    C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
C
      SUBROUTINE CALC_QUANTITIES
C
      IMPLICIT NONE
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE 'post3d.inc'
C
C
10    WRITE (*,*) 
     &  '  0   - Return to main menu'
      WRITE (*,*) 
     &  '  1   - Print (ASCII) MU_s and THETA calculated from RES file'
      WRITE (*,*) 
     &  '  2   - Write MU_s and THETA in .SP1 format'
      WRITE (*,*) 
     &  '  3   - Calculate Gas Flux : AVG(EP_g * V_g)'
      WRITE (*,*) 
     &  '  4   - Calculate Solids Flux AVG(EP_sm * V_sm)'
      WRITE (*,*) 
     &  '  5   - Calculate Correlations'
C
      CALL GET_SELECTION (SELECTION)
      IF (SELECTION .EQ. 0) THEN
        RETURN
      ELSEIF (SELECTION .EQ. 1) THEN
        CALL GET_MU_s
      ELSEIF (SELECTION.EQ.2) THEN
        CALL GET_GRANULAR_QTY
      ELSEIF (SELECTION.EQ.3) THEN
         CALL GAS_FLUX
      ELSEIF (SELECTION.EQ.4) THEN
         CALL SOL_FLUX
      ELSEIF (SELECTION.EQ.5) THEN
         CALL CALC_CORR_TYPE_1
      END IF
      GOTO 10 
C
      END
