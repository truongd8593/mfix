!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: GET_WALLS_BC                                           C
!  Purpose: Find and validate i, j, k locations for walls BC's         C
!                                                                      C
!  Author: P. Nicoletti                               Date: 10-DEC-91  C
!  Reviewer: M.SYAMLAL, W.ROGERS, P.NICOLETTI         Date: 27-JAN-92  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: BC_TYPE, BC_DEFINED, BC_X_w, BC_X_e, BC_Y_s   C
!                        BC_Y_n, BC_Z_b, BC_Z_t, DX, DY, DZ, IMAX      C
!                        JMAX, KMAX, IMAX2, JAMX2, KMAX2               C
!  Variables modified: BC_I_w, BC_I_e, BC_J_s, BC_J_n, BC_K_b, BC_K_t  C
!                      ICBC_FLAG                                       C
!                                                                      C
!  Local variables: I_w, I_e, J_s, J_n, K_b, K_t, BC, L1, L2, L3, L4   C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE GET_WALLS_BC 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE geometry
      USE fldvar
      USE physprop
      USE bc
      USE indices
      USE funits 
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
!             loop/variable indices
      INTEGER I , J , K , IJK
!
!             calculated indices of the wall boundary
      INTEGER I_w , I_e , J_s , J_n , K_b , K_t
!
!             loop index
      INTEGER BCV
!
!             Last twodigits of BC
      INTEGER BC2
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      LOGICAL , EXTERNAL :: COMPARE 
!-----------------------------------------------
      INCLUDE 'function.inc'
!
! FIND THE WALLS
!
      DO BCV = 1, DIMENSION_BC 
         IF (BC_DEFINED(BCV)) THEN 
            IF (BC_TYPE(BCV)=='FREE_SLIP_WALL' .OR. BC_TYPE(BCV)=='NO_SLIP_WALL'&
                .OR. BC_TYPE(BCV)=='PAR_SLIP_WALL') THEN 
!
               IF (BC_X_W(BCV)/=UNDEFINED .AND. BC_X_E(BCV)/=UNDEFINED) THEN 
                  IF (NO_I) THEN 
                     I_W = 1 
                     I_E = 1 
                  ELSE 
                     CALL CALC_CELL (BC_X_W(BCV), DX, IMAX, I_W) 
                     I_W = I_W + 1 
                     CALL CALC_CELL (BC_X_E(BCV), DX, IMAX, I_E) 
                     IF (BC_X_W(BCV) == BC_X_E(BCV)) THEN 
                        IF (COMPARE(BC_X_W(BCV),XMIN)) THEN 
                           I_W = 1 
                           I_E = 1 
                        ELSE IF (COMPARE(BC_X_W(BCV),XLENGTH)) THEN 
                           I_W = IMAX2 
                           I_E = IMAX2 
                        ENDIF 
                     ENDIF 
                  ENDIF 
                  IF (BC_I_W(BCV)/=UNDEFINED_I .OR. BC_I_E(BCV)/=UNDEFINED_I) &
                     THEN 
                     CALL LOCATION_CHECK (BC_I_W(BCV), I_W, BCV, 'BC - west') 
                     CALL LOCATION_CHECK (BC_I_E(BCV), I_E, BCV, 'BC - east') 
                  ELSE 
                     BC_I_W(BCV) = I_W 
                     BC_I_E(BCV) = I_E 
                  ENDIF 
               ENDIF 
!
               IF (BC_Y_S(BCV)/=UNDEFINED .AND. BC_Y_N(BCV)/=UNDEFINED) THEN 
                  IF (NO_J) THEN 
                     J_S = 1 
                     J_N = 1 
                  ELSE 
                     CALL CALC_CELL (BC_Y_S(BCV), DY, JMAX, J_S) 
                     J_S = J_S + 1 
                     CALL CALC_CELL (BC_Y_N(BCV), DY, JMAX, J_N) 
                     IF (BC_Y_S(BCV) == BC_Y_N(BCV)) THEN 
                        IF (COMPARE(BC_Y_S(BCV),ZERO)) THEN 
                           J_S = 1 
                           J_N = 1 
                        ELSE IF (COMPARE(BC_Y_S(BCV),YLENGTH)) THEN 
                           J_S = JMAX2 
                           J_N = JMAX2 
                        ENDIF 
                     ENDIF 
                  ENDIF 
                  IF (BC_J_S(BCV)/=UNDEFINED_I .OR. BC_J_N(BCV)/=UNDEFINED_I) &
                     THEN 
                     CALL LOCATION_CHECK (BC_J_S(BCV), J_S, BCV, 'BC - south') 
                     CALL LOCATION_CHECK (BC_J_N(BCV), J_N, BCV, 'BC - north') 
                  ELSE 
                     BC_J_S(BCV) = J_S 
                     BC_J_N(BCV) = J_N 
                  ENDIF 
               ENDIF 
!
               IF (BC_Z_B(BCV)/=UNDEFINED .AND. BC_Z_T(BCV)/=UNDEFINED) THEN 
                  IF (NO_K) THEN 
                     K_B = 1 
                     K_T = 1 
                  ELSE 
                     CALL CALC_CELL (BC_Z_B(BCV), DZ, KMAX, K_B) 
                     K_B = K_B + 1 
                     CALL CALC_CELL (BC_Z_T(BCV), DZ, KMAX, K_T) 
                     IF (BC_Z_B(BCV) == BC_Z_T(BCV)) THEN 
                        IF (COMPARE(BC_Z_B(BCV),ZERO)) THEN 
                           K_B = 1 
                           K_T = 1 
                        ELSE IF (COMPARE(BC_Z_B(BCV),ZLENGTH)) THEN 
                           K_B = KMAX2 
                           K_T = KMAX2 
                        ENDIF 
                     ENDIF 
                  ENDIF 
                  IF (BC_K_B(BCV)/=UNDEFINED_I .OR. BC_K_T(BCV)/=UNDEFINED_I) &
                     THEN 
                     CALL LOCATION_CHECK (BC_K_B(BCV), K_B, BCV, 'BC - bottom') 
                     CALL LOCATION_CHECK (BC_K_T(BCV), K_T, BCV, 'BC - top') 
                  ELSE 
                     BC_K_B(BCV) = K_B 
                     BC_K_T(BCV) = K_T 
                  ENDIF 
               ENDIF 
!
! CHECK FOR VALID VALUES
!
               IF (BC_K_B(BCV)<1 .OR. BC_K_B(BCV)>KMAX2) GO TO 900 
               IF (BC_J_S(BCV)<1 .OR. BC_J_S(BCV)>JMAX2) GO TO 900 
               IF (BC_I_W(BCV)<1 .OR. BC_I_W(BCV)>IMAX2) GO TO 900 
               IF (BC_K_T(BCV)<1 .OR. BC_K_T(BCV)>KMAX2) GO TO 900 
               IF (BC_J_N(BCV)<1 .OR. BC_J_N(BCV)>JMAX2) GO TO 900 
               IF (BC_I_E(BCV)<1 .OR. BC_I_E(BCV)>IMAX2) GO TO 900 
               IF (BC_K_B(BCV) > BC_K_T(BCV)) GO TO 900 
               IF (BC_J_S(BCV) > BC_J_N(BCV)) GO TO 900 
               IF (BC_I_W(BCV) > BC_I_E(BCV)) GO TO 900 
!
               DO K = BC_K_B(BCV), BC_K_T(BCV) 
                  DO J = BC_J_S(BCV), BC_J_N(BCV) 
                     DO I = BC_I_W(BCV), BC_I_E(BCV) 
                        IJK = FUNIJK(I,J,K) 
                        SELECT CASE (BC_TYPE(BCV))  
                        CASE ('FREE_SLIP_WALL')  
                           ICBC_FLAG(IJK)(1:1) = 'S' 
                        CASE ('NO_SLIP_WALL')  
                           ICBC_FLAG(IJK)(1:1) = 'W' 
                        CASE ('PAR_SLIP_WALL')  
                           ICBC_FLAG(IJK)(1:1) = 's' 
                        END SELECT 
                        BC2 = MOD(BCV,100) 
                        WRITE (ICBC_FLAG(IJK)(2:3), 1000) BC2 
                     END DO 
                  END DO 
               END DO 
            ENDIF 
!
         ENDIF 
      END DO 
      RETURN  
 1000 FORMAT(I2.2) 
!
! HERE IF AN ERROR OCCURRED
!
  900 CONTINUE 
      CALL ERROR_ROUTINE ('GET_WALLS_BC', 'Invalid BC location specified', 0, 2&
         ) 
      WRITE (UNIT_LOG, *) ' BC number = ', BCV 
      WRITE (UNIT_LOG, *) ' BC_I_w(BCV) = ', BC_I_W(BCV) 
      WRITE (UNIT_LOG, *) ' BC_I_e(BCV) = ', BC_I_E(BCV) 
      WRITE (UNIT_LOG, *) ' BC_J_s(BCV) = ', BC_J_S(BCV) 
      WRITE (UNIT_LOG, *) ' BC_J_n(BCV) = ', BC_J_N(BCV) 
      WRITE (UNIT_LOG, *) ' BC_K_b(BCV) = ', BC_K_B(BCV) 
      WRITE (UNIT_LOG, *) ' BC_K_t(BCV) = ', BC_K_T(BCV) 
      CALL ERROR_ROUTINE (' ', ' ', 1, 3) 
      RETURN  
      END SUBROUTINE GET_WALLS_BC 
