!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CHECK_DATA_20                                          C
!  Purpose: Check whether the sum of volume fractions is 1.0           C
!           and EP_g >= EP_star. Set miscellaneous constants           C
!                                                                      C
!  Author: M. Syamlal                                 Date: 30-JAN-92  C
!  Reviewer: P. Nicoletti, W. Rogers, S. Venkatesan   Date: 31-JAN-92  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: IMAX2, JMAX2, KMAX2, MMAX, EP_g, ROP_s, RO_s  C
!  Variables modified: I, J, K, M, MAX_DELP                            C
!                                                                      C
!  Local variables: ABORT, DIF                                         C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CHECK_DATA_20 
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
      USE toleranc 
      USE fldvar
      USE run
      USE geometry
      USE constant
      USE physprop
      USE indices
      USE funits 
      USE visc_g
      USE rxns
      USE compar   !//d
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!
!                      Indices
      INTEGER          I, J, K, IJK, IMJK, IPJK, IJMK, IJPK, IJKM, IJKP, &
                      IJKW, IJKE, IJKS, IJKN, IJKB, IJKT, &
                      IM, JM, KM, Lr
!
!                      Solids phase
      INTEGER          M
!
!                      Species index
      INTEGER          N
!
!                      Logical variable to set, if there is an error
      LOGICAL          ABORT
!
!                      Whether L_scale is nonzero
      LOGICAL          NONZERO
!
!                      Local index
      INTEGER          L
!
!                      1.0 - sum of all volume fractions
      DOUBLE PRECISION DIF
!-----------------------------------------------
      INCLUDE 'function.inc'
!
      CALL START_LOG 
      ABORT = .FALSE. 
      NONZERO = .FALSE. 
!
!  Check whether all field variables are initialized
!
      DO K = 1, KMAX2 
         DO J = 1, JMAX2 
            DO I = 1, IMAX2 
               IJK = FUNIJK(I,J,K) 
               IF (.NOT.WALL_AT(IJK)) THEN 
                  CALL SET_INDEX1 (IJK, I, J, K, IMJK, IPJK, IJMK, IJPK, IJKM, &
                     IJKP, IJKW, IJKE, IJKS, IJKN, IJKB, IJKT, IM, JM, KM) 
                  IF (EP_G(IJK) == UNDEFINED) THEN 
                     IF (.NOT.ABORT) THEN 
                        WRITE (UNIT_LOG, 1000) 
                        ABORT = .TRUE. 
                     ENDIF 
                     WRITE (UNIT_LOG, 1010) I, J, K, 'EP_g' 
                  ENDIF 
                  IF (P_G(IJK) == UNDEFINED) THEN 
                     IF (.NOT.ABORT) THEN 
                        WRITE (UNIT_LOG, 1000) 
                        ABORT = .TRUE. 
                     ENDIF 
                     WRITE (UNIT_LOG, 1010) I, J, K, 'P_g' 
                  ENDIF 
                  IF (P_STAR(IJK) == UNDEFINED) THEN 
                     IF (.NOT.ABORT) THEN 
                        WRITE (UNIT_LOG, 1000) 
                        ABORT = .TRUE. 
                     ENDIF 
                     WRITE (UNIT_LOG, 1010) I, J, K, 'P_star' 
                  ENDIF 
                  IF (RO_G(IJK) == UNDEFINED) THEN 
                     IF (.NOT.ABORT) THEN 
                        WRITE (UNIT_LOG, 1000) 
                        ABORT = .TRUE. 
                     ENDIF 
                     WRITE (UNIT_LOG, 1010) I, J, K, 'RO_g' 
                  ENDIF 
                  IF (ROP_G(IJK) == UNDEFINED) THEN 
                     IF (.NOT.ABORT) THEN 
                        WRITE (UNIT_LOG, 1000) 
                        ABORT = .TRUE. 
                     ENDIF 
                     WRITE (UNIT_LOG, 1010) I, J, K, 'ROP_g' 
                  ENDIF 
                  IF (T_G(IJK) == UNDEFINED) THEN 
                     IF (ENERGY_EQ .OR. RO_G0==UNDEFINED .OR. MU_G0==UNDEFINED&
                        ) THEN 
                        IF (.NOT.ABORT) THEN 
                           WRITE (UNIT_LOG, 1000) 
                           ABORT = .TRUE. 
                        ENDIF 
                        WRITE (UNIT_LOG, 1010) I, J, K, 'T_g' 
                     ENDIF 
                  ENDIF 
                  IF (U_G(IMJK) == UNDEFINED) THEN 
                     IF (.NOT.ABORT) THEN 
                        WRITE (UNIT_LOG, 1000) 
                        ABORT = .TRUE. 
                     ENDIF 
                     WRITE (UNIT_LOG, 1010) I - 1, J, K, 'U_g' 
                  ENDIF 
                  IF (U_G(IJK)==UNDEFINED .AND. I/=IMAX2) THEN 
                     IF (.NOT.ABORT) THEN 
                        WRITE (UNIT_LOG, 1000) 
                        ABORT = .TRUE. 
                     ENDIF 
                     WRITE (UNIT_LOG, 1010) I, J, K, 'U_g' 
                  ENDIF 
                  IF (V_G(IJMK) == UNDEFINED) THEN 
                     IF (.NOT.ABORT) THEN 
                        WRITE (UNIT_LOG, 1000) 
                        ABORT = .TRUE. 
                     ENDIF 
                     WRITE (UNIT_LOG, 1010) I, J - 1, K, 'V_g' 
                  ENDIF 
                  IF (V_G(IJK)==UNDEFINED .AND. J/=JMAX2) THEN 
                     IF (.NOT.ABORT) THEN 
                        WRITE (UNIT_LOG, 1000) 
                        ABORT = .TRUE. 
                     ENDIF 
                     WRITE (UNIT_LOG, 1010) I, J, K, 'V_g' 
                  ENDIF 
                  IF (W_G(IJKM) == UNDEFINED) THEN 
                     IF (.NOT.ABORT) THEN 
                        WRITE (UNIT_LOG, 1000) 
                        ABORT = .TRUE. 
                     ENDIF 
                     WRITE (UNIT_LOG, 1010) I, J, K - 1, 'W_g' 
                  ENDIF 
                  IF (W_G(IJK)==UNDEFINED .AND. K/=KMAX2) THEN 
                     IF (.NOT.ABORT) THEN 
                        WRITE (UNIT_LOG, 1000) 
                        ABORT = .TRUE. 
                     ENDIF 
                     WRITE (UNIT_LOG, 1010) I, J, K, 'W_g' 
                  ENDIF 
                  IF (SPECIES_EQ(0) .OR. RO_G0==UNDEFINED .AND. MW_AVG==&
                     UNDEFINED) THEN 
                     DO N = 1, NMAX(0) 
                        IF (X_G(IJK,N) == UNDEFINED) THEN 
                           IF (.NOT.ABORT) THEN 
                              WRITE (UNIT_LOG, 1000) 
                              ABORT = .TRUE. 
                           ENDIF 
                           WRITE (UNIT_LOG, 1012) I, J, K, N, 'X_g' 
                        ENDIF 
                     END DO 
                  ENDIF 
                  DO M = 1, MMAX 
                     IF (ROP_S(IJK,M) == UNDEFINED) THEN 
                        IF (.NOT.ABORT) THEN 
                           WRITE (UNIT_LOG, 1000) 
                           ABORT = .TRUE. 
                        ENDIF 
                        WRITE (UNIT_LOG, 1011) I, J, K, M, 'ROP_s' 
                     ENDIF 
                     IF (T_S(IJK,M) == UNDEFINED) THEN 
                        IF (ENERGY_EQ) THEN 
                           IF (.NOT.ABORT) THEN 
                              WRITE (UNIT_LOG, 1000) 
                              ABORT = .TRUE. 
                           ENDIF 
                           WRITE (UNIT_LOG, 1011) I, J, K, M, 'T_s' 
                        ENDIF 
                     ENDIF 
                     IF (U_S(IMJK,M) == UNDEFINED) THEN 
                        IF (.NOT.ABORT) THEN 
                           WRITE (UNIT_LOG, 1000) 
                           ABORT = .TRUE. 
                        ENDIF 
                        WRITE (UNIT_LOG, 1011) I - 1, J, K, M, 'U_s' 
                     ENDIF 
                     IF (U_S(IJK,M)==UNDEFINED .AND. I/=IMAX2) THEN 
                        IF (.NOT.ABORT) THEN 
                           WRITE (UNIT_LOG, 1000) 
                           ABORT = .TRUE. 
                        ENDIF 
                        WRITE (UNIT_LOG, 1011) I, J, K, M, 'U_s' 
                     ENDIF 
                     IF (V_S(IJMK,M) == UNDEFINED) THEN 
                        IF (.NOT.ABORT) THEN 
                           WRITE (UNIT_LOG, 1000) 
                           ABORT = .TRUE. 
                        ENDIF 
                        WRITE (UNIT_LOG, 1011) I, J - 1, K, M, 'V_s' 
                     ENDIF 
                     IF (V_S(IJK,M)==UNDEFINED .AND. J/=JMAX2) THEN 
                        IF (.NOT.ABORT) THEN 
                           WRITE (UNIT_LOG, 1000) 
                           ABORT = .TRUE. 
                        ENDIF 
                        WRITE (UNIT_LOG, 1011) I, J, K, M, 'V_s' 
                     ENDIF 
                     IF (W_S(IJKM,M) == UNDEFINED) THEN 
                        IF (.NOT.ABORT) THEN 
                           WRITE (UNIT_LOG, 1000) 
                           ABORT = .TRUE. 
                        ENDIF 
                        WRITE (UNIT_LOG, 1011) I, J, K - 1, M, 'W_s' 
                     ENDIF 
                     IF (W_S(IJK,M)==UNDEFINED .AND. K/=KMAX2) THEN 
                        IF (.NOT.ABORT) THEN 
                           WRITE (UNIT_LOG, 1000) 
                           ABORT = .TRUE. 
                        ENDIF 
                        WRITE (UNIT_LOG, 1011) I, J, K, M, 'W_s' 
                     ENDIF 
                     IF (SPECIES_EQ(M)) THEN 
                        DO N = 1, NMAX(M) 
                           IF (X_S(IJK,M,N) == UNDEFINED) THEN 
                              IF (.NOT.ABORT) THEN 
                                 WRITE (UNIT_LOG, 1000) 
                                 ABORT = .TRUE. 
                              ENDIF 
                              WRITE (UNIT_LOG, 1013) I, J, K, M, N, 'X_s' 
                           ENDIF 
                        END DO 
                     ENDIF 
                  END DO 
               ENDIF 
            END DO 
         END DO 
      END DO 
      IF (ABORT) THEN 
         WRITE (UNIT_LOG, 1300) 
         CALL MFIX_EXIT 
      ENDIF 
!
      DO K = 1, KMAX2 
         DO J = 1, JMAX2 
            DO I = 1, IMAX2 
               IJK = FUNIJK(I,J,K) 
               IF (FLAG(IJK)==1 .OR. FLAG(IJK)==20) THEN 
!
!  Check the sum of volume fractions
!
                  DIF = ONE - EP_G(IJK) 
                  M = 1 
                  IF (MMAX > 0) THEN 
                     DIF = DIF - SUM(ROP_S(IJK,:MMAX)/RO_S(:MMAX)) 
                     M = MMAX + 1 
                  ENDIF 
                  IF (ABS(DIF) > ZERO_EP_S) THEN 
                     IF (.NOT.ABORT) THEN 
                        WRITE (UNIT_LOG, 1050) 
                        ABORT = .TRUE. 
                     ENDIF 
                     WRITE (UNIT_LOG, 1100) I, J, K 
                  ENDIF 
!
!  Check whether L_scale is non-zero anywhere
!
                  IF (L_SCALE(IJK) /= ZERO) NONZERO = .TRUE. 
               ENDIF 
            END DO 
         END DO 
      END DO 
      IF (ABORT) THEN 
         WRITE (UNIT_LOG, 1300) 
         CALL MFIX_EXIT 
      ENDIF 
!
!  Check whether MU_gmax is specified
!
      IF (NONZERO .AND. MU_GMAX==UNDEFINED) THEN 
         WRITE (UNIT_LOG, 1350) 
         CALL MFIX_EXIT 
      ENDIF 

      CALL END_LOG 
      RETURN  
 1000 FORMAT(/1X,70('*')//' From: CHECK_DATA_20',/&
         ' Message: The following field variables are undefined') 
 1010 FORMAT(1X,'I = ',I4,' J = ',I4,' K = ',I4,5X,A) 
 1011 FORMAT(1X,'I = ',I4,' J = ',I4,' K = ',I4,' M = ',I4,5X,A) 
 1012 FORMAT(1X,'I = ',I4,' J = ',I4,' K = ',I4,' N = ',I4,5X,A) 
 1013 FORMAT(1X,'I = ',I4,' J = ',I4,' K = ',I4,' M = ',I4,' N = ',I4,5X,A) 
 1050 FORMAT(/1X,70('*')//' From: CHECK_DATA_20',/&
         ' Message: The sum of volume fractions is not equal to 1',/&
         '          in the following cells:',/4X,'I',T14,'J',T24,'K') 
 1100 FORMAT(1X,I4,T11,I4,T21,I4,'  Sum of EP .NE. 1') 
 1300 FORMAT(/1X,70('*')/) 
 1350 FORMAT(/1X,70('*')//' From: CHECK_DATA_20',/&
         ' Message: Turbulent length scale is nonzero. Specify MU_gmax.',/1X,70&
         ('*')/) 
      END SUBROUTINE CHECK_DATA_20 
