!//? either we bcast all necessary info to root PE and then let it proceed with
!//? convergence check or let all PEs execute following conv. check provided 
!//? that all have the same residuals (which is the global residuals for each
!//? variable thru out the domain?) what do you think Sreekanth?


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CHECK_CONVERGENCE(NIT, MUSTIT, IER)                    C
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
      SUBROUTINE CHECK_CONVERGENCE(NIT, MUSTIT, IER) 
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
      USE mpi_utility !//SP
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
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
!                      sum of residuals
      DOUBLE PRECISION SUM, SUM_T, SUM_X
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
      SUM = RESID(RESID_P,0) + RESID(RESID_P,1) 
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
      IF (GRANULAR_ENERGY) THEN 
         DO M = 1, MMAX 
            SUM = SUM + RESID(RESID_TH,M) 
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
      IF(SUM<=TOL_RESID.AND.SUM_T<=TOL_RESID_T.AND.SUM_X<=TOL_RESID_X)THEN 
         MUSTIT = 0                              !converged 
      ELSE IF (SUM>=TOL_DIVERGE .OR. SUM_T>=TOL_DIVERGE .OR. SUM_X>=TOL_DIVERGE&
            ) THEN 
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
