!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CHECK_CONVERGENCE(NIT, errorpercent, MUSTIT, IER)      C
!  Purpose: Monitor convergence                                        C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 8-JUL-96   C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CHECK_CONVERGENCE(NIT, errorpercent, MUSTIT, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE geometry
      USE indices
      USE physprop
      USE run
      USE residual
      USE toleranc 
      USE mpi_utility
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!Maximum % error allowed in fluid continuity
!      DOUBLE PRECISION, PARAMETER :: MaxErrorPercent = 1.0E-6  
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Iteration number
      INTEGER          NIT
!
!                      whether to iterate (1) or not (0).
      INTEGER          MUSTIT
!
!                      Error index
      INTEGER          IER
!
!                      %error in fluid mass balance
      DOUBLE PRECISION errorpercent
!!
!                      sum of residuals
      DOUBLE PRECISION SUM, SUM_T, SUM_X, SUM_Th
!
!                      max of residuals
      DOUBLE PRECISION maxres
!
!                      index
      INTEGER          L, M, maxL, maxM, N, maxN
!
!                      to indicate undefined residual in species eq at the
!                      begining of iterations
      LOGICAL          NO_RESID
!-----------------------------------------------
!
!//SP
!      if(abs(errorpercent) > MaxErrorPercent)then
        SUM = RESID(RESID_P,0) 
!      else
!        SUM = zero
!      endif

      if(MMAX > 0) SUM = SUM + RESID(RESID_P,1) 
!
      DO M = 0, MMAX 
         SUM = SUM + RESID(RESID_RO,M) 
      END DO 
      DO M = 0, MMAX 
         SUM = SUM + RESID(RESID_U,M) 
      END DO 
      DO M = 0, MMAX 
         SUM = SUM + RESID(RESID_V,M) 
      END DO 
      IF (DO_K) THEN 
         DO M = 0, MMAX 
            SUM = SUM + RESID(RESID_W,M) 
         END DO 
      ENDIF 
!
      SUM_Th = zero
      IF (GRANULAR_ENERGY) THEN 
         DO M = 1, MMAX 
            SUM_Th = SUM_Th + RESID(RESID_TH,M) 
         END DO 
      ENDIF 
      
!//SP
!    call global_all_sum(SUM)
!
      SUM_T = ZERO 
      IF (ENERGY_EQ) THEN 
         DO M = 0, MMAX 
            SUM_T = SUM_T + RESID(RESID_T,M) 
         END DO 
      ENDIF 
!//SP
!    call global_all_sum(SUM_T)
      SUM_X = ZERO 
      NO_RESID = .FALSE. 
      DO M = 0, MMAX 
         IF (SPECIES_EQ(M)) THEN 
            DO N = 1, NMAX(M) 
               IF (RESID(RESID_X+(N-1),M) == UNDEFINED) NO_RESID = .TRUE. 
               SUM_X = SUM_X + RESID(RESID_X+(N-1),M) 
            END DO 
         ENDIF 
      END DO 
!//SP
!    call global_all_sum(SUM_X)
      IF (NO_RESID) SUM_X = TOL_RESID_X + ONE 
      
!
!  Find the variable with maximum residual
!
      IF (RESID_INDEX(MAX_RESID_INDEX,1) == UNDEFINED_I) THEN 
         MAXRES = ZERO 
         DO L = 1, NRESID 
            DO M = 0, MMAX 
               IF (RESID(L,M) >= MAXRES) THEN 
                  MAXRES = RESID(L,M) 
                  MAXL = L 
                  MAXM = M 
                  IF (L >= RESID_X) THEN 
                     MAXN = L - RESID_X + 1 
                  ELSE 
                     MAXN = UNDEFINED_I 
                  ENDIF 
               ENDIF 
            END DO 
         END DO 
         IF (MAXN == UNDEFINED_I) THEN 
            WRITE (RESID_STRING(MAX_RESID_INDEX), '(A1,I1)') RESID_PREFIX(MAXL)&
               , MAXM 
         ELSE 
            WRITE (RESID_STRING(MAX_RESID_INDEX), '(A1,I1,I2.0)') 'X', MAXM, &
               MAXN 
         ENDIF 
      ENDIF 
!
!  Detect whether the run is stalled
!
      IF (DETECT_STALL) THEN 
         IF (MOD(NIT,5) == 0) THEN 
            IF (NIT > 10) THEN 
               IF (SUM5_RESID <= SUM) THEN 
                  MUSTIT = 2                     !stalled 
                  RETURN  
               ENDIF 
            ENDIF 
            SUM5_RESID = SUM 
         ENDIF 
      ENDIF 
!
!    total residual
!
      IF(NIT == 1) THEN
        MUSTIT = 1
	RETURN
      ENDIF
      
      IF(SUM<=TOL_RESID .AND. SUM_T<=TOL_RESID_T .AND. &
         RESID(RESID_sc,0)<=TOL_RESID_Scalar .AND. SUM_X<=TOL_RESID_X &
	  .AND. RESID(RESID_ke,0)<=TOL_RESID_K_Epsilon &
          .AND. SUM_Th <=TOL_RESID_Th)THEN 
         MUSTIT = 0                              !converged 
      ELSE IF (SUM>=TOL_DIVERGE .OR. SUM_T>=TOL_DIVERGE .OR.&
                RESID(RESID_sc,0)>= TOL_DIVERGE .OR. SUM_X>=TOL_DIVERGE&
            .OR. RESID(RESID_ke,0)>= TOL_DIVERGE &
            .OR. SUM_Th >= TOL_DIVERGE ) THEN 
         IF (NIT /= 1) THEN 
            MUSTIT = 2                           !diverged 
         ELSE 
            MUSTIT = 1                           !not converged 
         ENDIF 
      ELSE 
         MUSTIT = 1                              !not converged 
      ENDIF 
!
      RETURN  
      END SUBROUTINE CHECK_CONVERGENCE 
      
!// Comments on the modifications for DMP version implementation            
!// 400 Added mpi_utility module and other global reduction (bcast) calls
