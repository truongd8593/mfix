!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: RRATES0(IER)                                           C
!  Purpose: Calculate reaction rates for various reactions in cell ijk C
!           using information from the data file                       C
!                                                                      C
!  Author: M. Syamlal                                 Date: 3-10-98    C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: MMAX, IJK, T_g, T_s1, D_p, X_g, X_s, EP_g,    C
!            P_g, HOR_g, HOR_s                                         C
!                                                                      C
!                                                                      C
!  Variables modified: M, N, R_gp, R_sp, RoX_gc, RoX_sc, SUM_R_g,      C
!                      SUM_R_s                                         C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
!
      SUBROUTINE RRATES0(IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE fldvar
      USE rxns
      USE energy
      USE geometry
      USE run
      USE indices
      USE physprop
      USE constant
      USE funits 
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Error index
      INTEGER          IER
!
!                      Local phase and species indices
      INTEGER          L, LM, M, N, LR, ID

!                      cell index
      INTEGER          IJK

      DOUBLE PRECISION stmw, ex, Tr, EP, RATE(DIMENSION_RXN) 
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'ep_s2.inc'
!
!
!     Initialize arrays to zero
      CALL ZERO_ARRAY (SUM_R_G, IJKMAX2, IER) 
      CALL ZERO_ARRAY (HOR_G, IJKMAX2, IER) 
      DO N = 1, NMAX(0) 
         CALL ZERO_ARRAY (R_GP(1,N), IJKMAX2, IER) 
         CALL ZERO_ARRAY (ROX_GC(1,N), IJKMAX2, IER) 
      END DO 
      DO M = 1, MMAX 
         CALL ZERO_ARRAY (SUM_R_S(1,M), IJKMAX2, IER) 
         CALL ZERO_ARRAY (HOR_S(1,M), IJKMAX2, IER) 
         DO N = 1, NMAX(M) 
            CALL ZERO_ARRAY (R_SP(1,M,N), IJKMAX2, IER) 
            CALL ZERO_ARRAY (ROX_SC(1,M,N), IJKMAX2, IER) 
         END DO 
      END DO 
      
      
!!$omp  parallel do private( IJK, L, LM, M, N, LR, ID, stmw, ex, Tr, EP, RATE )    
      DO IJK = 1, IJKMAX2 
         IF (FLUID_AT(IJK)) THEN 
!
!
!  Define the reaction rates and enthalpy changes due to reaction
!  (not given in the data file) at the end of Section 1.
!
!1111111111111111111111111111111111111111111111111111111111111111111111111111111
!
! 1. Write the rates of various reactions:
!
            DO LR = 1, NO_OF_RXNS 
!
               IF (GOT_RATE(LR)) THEN 
!
                  IF (RATE_M4T(LR) == 0) THEN 
                     TR = T_G(IJK) 
                     EP = EP_G(IJK) 
                  ELSE 
                     TR = T_S(IJK,RATE_M4T(LR)) 
                     EP = EP_S(IJK,RATE_M4T(LR)) 
                  ENDIF 
!
                  IF (EP > ZERO) THEN 
                     RATE(LR) = RATE_FAC(LR,1)*TR**RATE_FAC(LR,2)*EXP((-&
                        RATE_FAC(LR,3)/TR))*EP 
!
                     DO ID = 1, N_ALL 
                        EX = RATE_EXP(LR,ID) 
                        IF (EX /= ZERO) THEN 
!
                           M = SPECIES_ID2N(ID,1) 
                           N = SPECIES_ID2N(ID,2) 
!
                           IF (M == 0) THEN 
                              RATE(LR) = RATE(LR)*(RO_G(IJK)*X_G(IJK,N)/MW_G(N)&
                                 )**EX 
                           ELSE 
                              RATE(LR) = RATE(LR)*(RO_S(M)*X_S(IJK,M,N)/MW_S(M,&
                                 N))**EX 
                           ENDIF 
!
                        ENDIF 
                     END DO 
                  ELSE 
                     RATE(LR) = ZERO 
                  ENDIF 
!
               ELSE 
                  RATE(LR) = UNDEFINED 
               ENDIF 
            END DO 
            DO LR = 1, NO_OF_RXNS 
!
               IF (RATE(LR) /= UNDEFINED) THEN 
!
                  DO ID = 1, N_ALL 
                     STMW = STOICHXMW(LR,ID) 
                     IF (STMW /= ZERO) THEN 
!
                        M = SPECIES_ID2N(ID,1) 
                        N = SPECIES_ID2N(ID,2) 
!
                        IF (M == 0) THEN 
                           IF (STMW > ZERO) THEN !production 
                              R_GP(IJK,N) = R_GP(IJK,N) + RATE(LR)*STMW 
!
                           ELSE 
!                                            !consumption
                              EX = RATE_EXP(LR,ID) 
                              IF (EX <= ONE) THEN 
                                 IF (X_G(IJK,N) > SMALL_NUMBER) ROX_GC(IJK,N)&
                                     = ROX_GC(IJK,N) - RATE(LR)*STMW/X_G(IJK,N) 
                              ELSE 
                                 R_GP(IJK,N) = R_GP(IJK,N) - (EX - ONE)*RATE(LR&
                                    )*STMW 
!
                                 IF (X_G(IJK,N) > SMALL_NUMBER) ROX_GC(IJK,N)&
                                     = ROX_GC(IJK,N) - EX*RATE(LR)*STMW/X_G(IJK&
                                    ,N) 
                              ENDIF 
!
                           ENDIF 
                        ELSE 
                           IF (STMW > ZERO) THEN !production 
                              R_SP(IJK,M,N) = R_SP(IJK,M,N) + RATE(LR)*STMW 
!
                           ELSE 
!                                            !consumption
!
                              EX = RATE_EXP(LR,ID) 
                              IF (EX <= ONE) THEN 
                                 IF (X_S(IJK,M,N) > SMALL_NUMBER) ROX_SC(IJK,M,&
                                    N) = ROX_SC(IJK,M,N) - RATE(LR)*STMW/X_S(&
                                    IJK,M,N) 
                              ELSE 
                                 R_SP(IJK,M,N) = R_SP(IJK,M,N) - (EX - ONE)*&
                                    RATE(LR)*STMW 
!
                                 IF (X_S(IJK,M,N) > SMALL_NUMBER) ROX_SC(IJK,M,&
                                    N) = ROX_SC(IJK,M,N) - EX*RATE(LR)*STMW/X_S&
                                    (IJK,M,N) 
                              ENDIF 
!
                           ENDIF 
                        ENDIF 
!
                     ENDIF 
                  END DO 
               ELSE 
                  WRITE (UNIT_LOG, 1000) LR, RXN_NAME(LR) 
                  STOP  
               ENDIF 
            END DO 
            DO LR = 1, NO_OF_RXNS 
!
               M = RATE_M4T(LR) 
               IF (M == 0) THEN 
                  HOR_G(IJK) = HOR_G(IJK) + RATE(LR)*DELTA_H(LR) 
               ELSE 
                  HOR_S(IJK,M) = HOR_S(IJK,M) + RATE(LR)*DELTA_H(LR) 
               ENDIF 
            END DO 
            IF (SPECIES_EQ(0)) THEN 
               N = 1 
               IF (NMAX(0) > 0) THEN 
                  SUM_R_G(IJK) = SUM_R_G(IJK) + SUM(R_GP(IJK,:NMAX(0))-ROX_GC(&
                     IJK,:NMAX(0))*X_G(IJK,:NMAX(0))) 
                  N = NMAX(0) + 1 
               ENDIF 
            ENDIF 
!
            DO M = 1, MMAX 
               IF (SPECIES_EQ(M)) THEN 
                  N = 1 
                  IF (NMAX(M) > 0) THEN 
                     SUM_R_S(IJK,M) = SUM_R_S(IJK,M) + SUM(R_SP(IJK,M,:NMAX(M))&
                        -ROX_SC(IJK,M,:NMAX(M))*X_S(IJK,M,:NMAX(M))) 
                     N = NMAX(M) + 1 
                  ENDIF 
               ENDIF 
            END DO 
            DO L = 0, MMAX 
               DO M = L + 1, MMAX 
                  LM = L + 1 + (M - 1)*M/2 
                  R_PHASE(IJK,LM) = ZERO 
               END DO 
            END DO 
            DO LR = 1, NO_OF_RXNS 
               DO L = 0, MMAX 
                  DO M = L + 1, MMAX 
                     LM = L + 1 + (M - 1)*M/2 
                     IF (R_TEMP(LR,L,M) /= UNDEFINED) THEN 
                        R_PHASE(IJK,LM) = R_PHASE(IJK,LM) + R_TEMP(LR,L,M)*RATE&
                           (LR) 
                     ELSE IF (R_TEMP(LR,M,L) /= UNDEFINED) THEN 
                        R_PHASE(IJK,LM) = R_PHASE(IJK,LM) - R_TEMP(LR,M,L)*RATE&
                           (LR) 
                     ELSE 
                        CALL START_LOG 
                        WRITE (UNIT_LOG, 1010) L, M 
                        CALL END_LOG 
                        STOP  
                     ENDIF 
                  END DO 
               END DO 
            END DO 
         ENDIF 
      END DO 
      RETURN  
 1000 FORMAT(/1X,70('*')//' From: RRATES0',/' Error: ',&
         'Reaction rate for reaction ',I2,' (',A,') not specified',/1X,70('*')/&
         ) 
 1010 FORMAT(/1X,70('*')//' From: RRATES0',/&
         ' Error: Mass transfer between phases ',I2,' and ',I2,&
         ' (R_temp) not specified',/1X,70('*')/) 
!
      END SUBROUTINE RRATES0 
