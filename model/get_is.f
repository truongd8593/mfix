!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: GET_IS                                                 C
!  Purpose: Find and validate i, j, k locations for IS's               C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-OCT-92  C
!  Reviewer: M. Syamlal                               Date: 11-DEC-92  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: IS_DEFINED, IS_X_w, IS_X_e, IS_Y_s, IS_Y_n    C
!                        IS_Z_b, IS_Z_t, DX, DY, DZ, IMAX, JMAX, KMAX  C
!  Variables modified: IS_I_w, IS_I_e, IS_J_s, IS_J_n, IS_K_b, IS_K_t  C
!                                                                      C
!  Local variables: IS, I_w, I_e, J_s, J_n, K_b, K_t                   C
!                   X_CONSTANT, Y_CONSTANT, Z_CONSTANT                 C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE GET_IS 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE geometry
      USE is
      USE indices
      USE funits 
      USE compar
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
! note : this routine will not work if indices are specified ... 
!        x_constant etc.
!
! local variables
!
!     loop/variable indices
      INTEGER ISV
!
!             calculated indices of the wall boundary
      INTEGER I_w , I_e , J_s , J_n , K_b , K_t
!
!             surface indictors
      LOGICAL X_CONSTANT, Y_CONSTANT, Z_CONSTANT
!-----------------------------------------------
!
! FIND THE FLOW SURFACES
!
      DO ISV = 1, DIMENSION_IS 
         IF (.NOT.IS_DEFINED(ISV)) CYCLE  
!
         X_CONSTANT = .TRUE. 
         Y_CONSTANT = .TRUE. 
         Z_CONSTANT = .TRUE. 
!
         IF (IS_X_W(ISV)/=UNDEFINED .AND. IS_X_E(ISV)/=UNDEFINED) THEN 
            CALL CALC_CELL (IS_X_W(ISV), DX, IMAX, I_W) 
            CALL CALC_CELL (IS_X_E(ISV), DX, IMAX, I_E) 
            IF (IS_X_W(ISV) /= IS_X_E(ISV)) THEN 
               X_CONSTANT = .FALSE. 
               I_W = I_W + 1 
               IF (IS_I_W(ISV)/=UNDEFINED_I .OR. IS_I_E(ISV)/=UNDEFINED_I) THEN 
                  CALL LOCATION_CHECK (IS_I_W(ISV), I_W, ISV, 'IS - west') 
                  CALL LOCATION_CHECK (IS_I_E(ISV), I_E, ISV, 'IS - east') 
               ENDIF 
            ENDIF 
            IS_I_W(ISV) = I_W 
            IS_I_E(ISV) = I_E 
         ELSE 
            IF (IS_I_W(ISV) /= UNDEFINED_I) CALL CALC_LOC (XMIN, DX, IS_I_W(ISV)&
               , IS_X_W(ISV)) 
!
            IF (IS_I_E(ISV) /= UNDEFINED_I) CALL CALC_LOC (XMIN, DX, IS_I_E(ISV)&
               , IS_X_E(ISV)) 
!
            IF (IS_X_W(ISV) /= IS_X_E(ISV)) X_CONSTANT = .FALSE. 
         ENDIF 
!
!  If there is no variation in the I direction set indices to 1
!
         IF (NO_I) THEN 
            IS_I_W(ISV) = 1 
            IS_I_E(ISV) = 1 
         ENDIF 
!
         IF (IS_Y_S(ISV)/=UNDEFINED .AND. IS_Y_N(ISV)/=UNDEFINED) THEN 
            CALL CALC_CELL (IS_Y_S(ISV), DY, JMAX, J_S) 
            CALL CALC_CELL (IS_Y_N(ISV), DY, JMAX, J_N) 
            IF (IS_Y_S(ISV) /= IS_Y_N(ISV)) THEN 
               Y_CONSTANT = .FALSE. 
               J_S = J_S + 1 
               IF (IS_J_S(ISV)/=UNDEFINED_I .OR. IS_J_N(ISV)/=UNDEFINED_I) THEN 
                  CALL LOCATION_CHECK (IS_J_S(ISV), J_S, ISV, 'IS - south') 
                  CALL LOCATION_CHECK (IS_J_N(ISV), J_N, ISV, 'IS - north') 
               ENDIF 
            ENDIF 
            IS_J_S(ISV) = J_S 
            IS_J_N(ISV) = J_N 
         ELSE 
            IF (IS_J_S(ISV) /= UNDEFINED_I) CALL CALC_LOC (ZERO, DY, IS_J_S(ISV)&
               , IS_Y_S(ISV)) 
!
            IF (IS_J_N(ISV) /= UNDEFINED_I) CALL CALC_LOC (ZERO, DY, IS_J_N(ISV)&
               , IS_Y_N(ISV)) 
!
            IF (IS_Y_S(ISV) /= IS_Y_N(ISV)) Y_CONSTANT = .FALSE. 
         ENDIF 
!
!  If there is no variation in the J direction set indices to 1
!
         IF (NO_J) THEN 
            IS_J_S(ISV) = 1 
            IS_J_N(ISV) = 1 
         ENDIF 
!
         IF (IS_Z_B(ISV)/=UNDEFINED .AND. IS_Z_T(ISV)/=UNDEFINED) THEN 
            CALL CALC_CELL (IS_Z_B(ISV), DZ, KMAX, K_B) 
            CALL CALC_CELL (IS_Z_T(ISV), DZ, KMAX, K_T) 
            IF (IS_Z_B(ISV) /= IS_Z_T(ISV)) THEN 
               Z_CONSTANT = .FALSE. 
               K_B = K_B + 1 
               IF (IS_K_B(ISV)/=UNDEFINED_I .OR. IS_K_T(ISV)/=UNDEFINED_I) THEN 
                  CALL LOCATION_CHECK (IS_K_B(ISV), K_B, ISV, 'IS - bottom') 
                  CALL LOCATION_CHECK (IS_K_T(ISV), K_T, ISV, 'IS - top') 
               ENDIF 
            ENDIF 
            IS_K_B(ISV) = K_B 
            IS_K_T(ISV) = K_T 
         ELSE 
            IF (IS_K_B(ISV) /= UNDEFINED_I) CALL CALC_LOC (ZERO, DZ, IS_K_B(ISV)&
               , IS_Z_B(ISV)) 
!
            IF (IS_K_T(ISV) /= UNDEFINED_I) CALL CALC_LOC (ZERO, DZ, IS_K_T(ISV)&
               , IS_Z_T(ISV)) 
!
            IF (IS_Z_B(ISV) /= IS_Z_T(ISV)) Z_CONSTANT = .FALSE. 
         ENDIF 
!
!  If there is no variation in the K direction set indices to 1
!
         IF (NO_K) THEN 
            IS_K_B(ISV) = 1 
            IS_K_T(ISV) = 1 
         ENDIF 
!
!  Check whether the boundary is a plane parallel to one of the three
!  coordinate planes, else check whether a direction is specified by IS_TYPE
!
         IF (X_CONSTANT .OR. Y_CONSTANT .OR. Z_CONSTANT) THEN 
            IF (IS_X_W(ISV)/=UNDEFINED .AND. IS_Y_S(ISV)/=UNDEFINED .AND. IS_Z_B(&
               ISV)/=UNDEFINED) CALL CHECK_PLANE (X_CONSTANT, Y_CONSTANT, &
               Z_CONSTANT, ISV, 'IS') 
         ELSE 
            IF (IS_TYPE(ISV)(1:1)/='X' .AND. IS_TYPE(ISV)(1:1)/='Y' .AND. IS_TYPE&
               (ISV)(1:1)/='Z') THEN 
               WRITE (UNIT_LOG, 1000) ISV 
               CALL MFIX_EXIT(myPE) 
            ENDIF 
         ENDIF 
!
! CHECK FOR VALID VALUES
!
         IF (IS_I_W(ISV)<1 .OR. IS_I_W(ISV)>IMAX2) GO TO 900 
         IF (IS_I_E(ISV)<1 .OR. IS_I_E(ISV)>IMAX2) GO TO 900 
         IF (IS_J_S(ISV)<1 .OR. IS_J_S(ISV)>JMAX2) GO TO 900 
         IF (IS_J_N(ISV)<1 .OR. IS_J_N(ISV)>JMAX2) GO TO 900 
         IF (IS_K_B(ISV)<1 .OR. IS_K_B(ISV)>KMAX2) GO TO 900 
         IF (IS_K_T(ISV)<1 .OR. IS_K_T(ISV)>KMAX2) GO TO 900 
         IF (IS_K_B(ISV) > IS_K_T(ISV)) GO TO 900 
         IF (IS_J_S(ISV) > IS_J_N(ISV)) GO TO 900 
         IF (IS_I_W(ISV) > IS_I_E(ISV)) GO TO 900 
!
      END DO 
      RETURN  
!
  900 CONTINUE 
      CALL ERROR_ROUTINE ('GET_IS', 'Invalid IS location specified', 0, 2) 
      WRITE (UNIT_LOG, *) ' IS number = ', ISV 
      WRITE (UNIT_LOG, *) ' IS_I_w(ISV) = ', IS_I_W(ISV) 
      WRITE (UNIT_LOG, *) ' IS_I_e(ISV) = ', IS_I_E(ISV) 
      WRITE (UNIT_LOG, *) ' IS_J_s(ISV) = ', IS_J_S(ISV) 
      WRITE (UNIT_LOG, *) ' IS_J_n(ISV) = ', IS_J_N(ISV) 
      WRITE (UNIT_LOG, *) ' IS_K_b(ISV) = ', IS_K_B(ISV) 
      WRITE (UNIT_LOG, *) ' IS_K_t(ISV) = ', IS_K_T(ISV) 
      CALL ERROR_ROUTINE (' ', ' ', 1, 3) 
      RETURN  
 1000 FORMAT(/1X,70('*')//' From: GET_IS','     IS # = ',I2,/&
         ' Message: Since a volume is specified IS_TYPE should have a',/&
         ' prefix (X_, Y_, or Z_) to show the IS orientation') 
      END SUBROUTINE GET_IS 
