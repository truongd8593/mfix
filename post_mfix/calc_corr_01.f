CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: CALC_CORR_01(FINISH)                                   C
C  Purpose: Calculate correlations between different variables in a    C
C           fluidized bed hydrodynamic model                           C
C                                                                      C
C  Author: M. Syamlal                                 Date: 01-AUG-92  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Revision Number:                                                    C
C  Purpose:                                                            C
C  Author:                                            Date: dd-mmm-yy  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Literature/Document References:                                     C
C                                                                      C
C  Variables referenced: KMAX2, JMAX2, IMAX2, EP_g, V_g                C
C  Variables modified: STARTED, NSUM, I, J, K, IJK, SUM_EP_g, SUM_V_g  C
C                      SUM_EPxEP_g, SUM_VxV_g, AVG_EP_g, AVG_V_g       C
C                      SDV_EP_g, SDV_V_g                               C
C                                                                      C
C  Local variables: SDV_2                                              C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
C
      SUBROUTINE CALC_CORR_01(FINISH,INTER)
C
      IMPLICIT NONE
C
C                      A flag to check whether the subroutine is called for
C                      the first time.
      INTEGER          RAND_NO
      PARAMETER (RAND_NO=946235187)
C
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE 'fldvar.inc'
      INCLUDE 'physprop.inc'
      INCLUDE 'geometry.inc'
      INCLUDE 'indices.inc'
      INCLUDE 'correl.inc'
      INCLUDE 'xforms.inc'
C
C  passed variables
C
C                      Logical variable to tell the subroutine to finish
C                      calculations. i.e. To compute averages and variances
      LOGICAL          FINISH,INTER
C
C
C  local variables
C
C                       temporary stograge for SDV**2
      DOUBLE PRECISION  SDV_2
C
C                       indices
      INTEGER           I, J, K, IJK
C
      INCLUDE 'function.inc'
C
C  If finish command is not given, do summations for correlations.
C
      IF(.NOT.FINISH)THEN
C
C  First check whether the subroutine is called the first time, by checking
C  whether the flag STARTED is correctly set to a specific number.  If the
C  subroutine is called the first time, do initializations.
C
        IF(STARTED .NE. RAND_NO)THEN
          STARTED = RAND_NO
          NSUM = 0
        DO 100 K = 1, KMAX2
        DO 100 J = 1, JMAX2
        DO 100 I = 1, IMAX2
          IJK = FUNIJK(I, J, K)
          SUM_EP_g(IJK) = ZERO
          SUM_EPxEP_g(IJK) = ZERO
          SUM_V_g(IJK) = ZERO
          SUM_VxV_g(IJK) = ZERO
100     CONTINUE
        ENDIF
C
C       Sum field variable values in all fluid cells
C
        NSUM = NSUM + 1
        DO 1000 K = 1, KMAX2
        DO 1000 J = 1, JMAX2
        DO 1000 I = 1, IMAX2
          IJK = FUNIJK(I, J, K)
          IF( FLUID_AT(IJK) ) THEN
            SUM_EP_g(IJK) = SUM_EP_g(IJK) + EP_g(IJK)
            SUM_EPxEP_g(IJK) = SUM_EPxEP_g(IJK) + EP_g(IJK)**2
            SUM_V_g(IJK) = SUM_V_g(IJK) + V_g(IJK)
            SUM_VxV_g(IJK) = SUM_VxV_g(IJK) + V_g(IJK)**2
          ENDIF
1000    CONTINUE
C
C  If finish command is given, calculate averages, variances,
C  and correlations
C
      ELSE
        IF(NSUM .LE. 0)THEN
          WRITE(*,5000)
          STOP
        ENDIF
        DO 2000 K = 1, KMAX2
        DO 2000 J = 1, JMAX2
        IF (DO_XFORMS) THEN
           CALL CHECK_INTER(INTER)
           IF (INTER) RETURN
        END IF
        DO 2000 I = 1, IMAX2
          IJK = FUNIJK(I, J, K)
          IF( FLUID_AT(IJK) ) THEN
            AVG_EP_g(IJK) = SUM_EP_g(IJK) / NSUM
            SDV_2 = SUM_EPxEP_g(IJK) / NSUM - AVG_EP_g(IJK)**2 
            SDV_2 = MAX(ZERO,SDV_2)
            SDV_EP_g(IJK) = SQRT( SDV_2 )
            AVG_V_g(IJK) = SUM_V_g(IJK) / NSUM
            SDV_2 = SUM_VxV_g(IJK) / NSUM - AVG_V_g(IJK)**2
            SDV_2 = MAX(ZERO,SDV_2)
            SDV_V_g(IJK) = SQRT( SDV_2 )
          ENDIF
2000    CONTINUE
      ENDIF
      RETURN
5000  FORMAT(/1X,70('*')//' From: CALC_CORR',
     &/' Message: Attempting averaging before doing summations',
     &/1X, 70('*')/)
      END
