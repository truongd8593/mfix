!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: GET_FLOW_BC                                            C
!  Purpose: Find and validate i, j, k locations for flow BC's          C
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
!  Variables referenced: BC_DEFINED, BC_X_w, BC_X_e, BC_Y_s, BC_Y_n    C
!                        BC_Z_b, BC_Z_t, DX, DY, DZ, IMAX, JMAX, KMAX  C
!  Variables modified: BC_I_w, BC_I_e, BC_J_s, BC_J_n, BC_K_b, BC_K_t  C
!                      ICBC_FLAG, BC_PLANE                             C
!                                                                      C
!  Local variables: BC, I, J, K, IJK, I_w, I_e, J_s, J_n, K_b, K_t   C
!                   ERROR, X_CONSTANT, Y_CONSTANT, Z_CONSTANT          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE GET_FLOW_BC 
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
      USE compar        !//d
      USE sendrecv      !//SP
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
      INTEGER BCV , I, J, K, IJK
!
!     Last two digits of BC
      INTEGER BC2
!
!             calculated indices of the wall boundary
      INTEGER I_w , I_e , J_s , J_n , K_b , K_t
!
!             Indices for error checking
      INTEGER I_wall, I_fluid, J_wall, J_fluid, K_wall, K_fluid, & 
             IJK_wall, IJK_fluid
!
!             error indicator
      LOGICAL ERROR
!
!             surface indictors
      LOGICAL X_CONSTANT, Y_CONSTANT, Z_CONSTANT
!-----------------------------------------------
      INCLUDE 'function.inc'

!//SP
       call send_recv(icbc_flag,2)

!//AIKEPARDBG
!      write(*,"('(PE ',I2,'): entered get_flow_bc')") myPE	!//AIKEPARDBG
!//AIKEPARDBG dump the ICBC_FLAG in matrix form to verify with serial version
!      DO K = Kstart3, Kend3                               !//AIKEPARDBG
!         write(UNIT_LOG,"('K = ',I5)") K                !//AIKEPARDBG 
!	 write(UNIT_LOG,"(7X,14(I3,2X))") (I,i=IMIN3,IMAX3)  !//AIKEPARDBG
!         DO J = Jstart3, Jend3                            !//AIKEPARDBG
!           write(UNIT_LOG,"(I5,')',$)") J               !//AIKEPARDBG	
!           DO I = Istart3, Iend3                          !//AIKEPARDBG
!             IJK = FUNIJK(I,J,K)                     !//AIKEPARDBG
!             write(UNIT_LOG,"(2X,A3,$)") ICBC_FLAG(IJK) !//AIKEPARDBG
!           END DO                                       !//AIKEPARDBG
!           write(UNIT_LOG,"(/)")                        !//AIKEPARDBG
!         END DO                                         !//AIKEPARDBG
!      END DO                                            !//AIKEPARDBG
!      call mfix_exit(myPE)	!//AIKEPARDBG


!
! FIND THE FLOW SURFACES
!
      ERROR = .FALSE. 
!
      DO BCV = 1, DIMENSION_BC 
!

         IF (BC_DEFINED(BCV)) THEN 
            IF (.NOT.(BC_TYPE(BCV)=='FREE_SLIP_WALL' .OR. BC_TYPE(BCV)==&
               'NO_SLIP_WALL' .OR. BC_TYPE(BCV)=='PAR_SLIP_WALL')) THEN 
!
               X_CONSTANT = .TRUE. 
               Y_CONSTANT = .TRUE. 
               Z_CONSTANT = .TRUE. 
!
               IF (BC_X_W(BCV)/=UNDEFINED .AND. BC_X_E(BCV)/=UNDEFINED) THEN 
                  CALL CALC_CELL (BC_X_W(BCV), DX, IMAX, I_W) 
                  CALL CALC_CELL (BC_X_E(BCV), DX, IMAX, I_E) 
                  IF (BC_X_W(BCV) /= BC_X_E(BCV)) THEN 
                     X_CONSTANT = .FALSE. 
                     I_W = I_W + 1 
                     IF (BC_I_W(BCV)/=UNDEFINED_I .OR. BC_I_E(BCV)/=UNDEFINED_I) &
                        THEN 
                        CALL LOCATION_CHECK (BC_I_W(BCV), I_W, BCV, 'BC - west') 
                        CALL LOCATION_CHECK (BC_I_E(BCV), I_E, BCV, 'BC - east') 
                     ENDIF 
                  ENDIF 
                  BC_I_W(BCV) = I_W 
                  BC_I_E(BCV) = I_E 
               ELSE 
                  IF (BC_I_W(BCV) /= UNDEFINED_I) CALL CALC_LOC (XMIN, DX, &
                     BC_I_W(BCV), BC_X_W(BCV)) 
!
                  IF (BC_I_E(BCV) /= UNDEFINED_I) CALL CALC_LOC (XMIN, DX, &
                     BC_I_E(BCV), BC_X_E(BCV)) 
!
                  IF (BC_X_W(BCV) /= BC_X_E(BCV)) X_CONSTANT = .FALSE. 
               ENDIF 
!
!  If there is no variation in the I direction set indices to 1
!
               IF (NO_I) THEN 
                  BC_I_W(BCV) = 1 
                  BC_I_E(BCV) = 1 
               ENDIF 
!
               IF (BC_Y_S(BCV)/=UNDEFINED .AND. BC_Y_N(BCV)/=UNDEFINED) THEN 
                  CALL CALC_CELL (BC_Y_S(BCV), DY, JMAX, J_S) 
                  CALL CALC_CELL (BC_Y_N(BCV), DY, JMAX, J_N) 
                  IF (BC_Y_S(BCV) /= BC_Y_N(BCV)) THEN 
                     Y_CONSTANT = .FALSE. 
                     J_S = J_S + 1 
                     IF (BC_J_S(BCV)/=UNDEFINED_I .OR. BC_J_N(BCV)/=UNDEFINED_I) &
                        THEN 
                        CALL LOCATION_CHECK (BC_J_S(BCV), J_S, BCV, 'BC - south') 
                        CALL LOCATION_CHECK (BC_J_N(BCV), J_N, BCV, 'BC - north') 
                     ENDIF 
                  ENDIF 
                  BC_J_S(BCV) = J_S 
                  BC_J_N(BCV) = J_N 
               ELSE 
                  IF (BC_J_S(BCV) /= UNDEFINED_I) CALL CALC_LOC (ZERO, DY, &
                     BC_J_S(BCV), BC_Y_S(BCV)) 
!
                  IF (BC_J_N(BCV) /= UNDEFINED_I) CALL CALC_LOC (ZERO, DY, &
                     BC_J_N(BCV), BC_Y_N(BCV)) 
!
                  IF (BC_Y_S(BCV) /= BC_Y_N(BCV)) Y_CONSTANT = .FALSE. 
               ENDIF 
!
!  If there is no variation in the J direction set indices to 1
!
               IF (NO_J) THEN 
                  BC_J_S(BCV) = 1 
                  BC_J_N(BCV) = 1 
               ENDIF 
!
               IF (BC_Z_B(BCV)/=UNDEFINED .AND. BC_Z_T(BCV)/=UNDEFINED) THEN 
                  CALL CALC_CELL (BC_Z_B(BCV), DZ, KMAX, K_B) 
                  CALL CALC_CELL (BC_Z_T(BCV), DZ, KMAX, K_T) 
                  IF (BC_Z_B(BCV) /= BC_Z_T(BCV)) THEN 
                     Z_CONSTANT = .FALSE. 
                     K_B = K_B + 1 
                     IF (BC_K_B(BCV)/=UNDEFINED_I .OR. BC_K_T(BCV)/=UNDEFINED_I) &
                        THEN 
                        CALL LOCATION_CHECK (BC_K_B(BCV), K_B, BCV, 'BC - bottom'&
                           ) 
                        CALL LOCATION_CHECK (BC_K_T(BCV), K_T, BCV, 'BC - top') 
                     ENDIF 
                  ENDIF 
                  BC_K_B(BCV) = K_B 
                  BC_K_T(BCV) = K_T 
               ELSE 
                  IF (BC_K_B(BCV) /= UNDEFINED_I) CALL CALC_LOC (ZERO, DZ, &
                     BC_K_B(BCV), BC_Z_B(BCV)) 
!
                  IF (BC_K_T(BCV) /= UNDEFINED_I) CALL CALC_LOC (ZERO, DZ, &
                     BC_K_T(BCV), BC_Z_T(BCV)) 
!
                  IF (BC_Z_B(BCV) /= BC_Z_T(BCV)) Z_CONSTANT = .FALSE. 
               ENDIF 
!
!  If there is no variation in the K direction set indices to 1
!
               IF (NO_K) THEN 
                  BC_K_B(BCV) = 1 
                  BC_K_T(BCV) = 1 
               ENDIF 
!
!  Check whether the boundary is a plane parallel to one of the three
!  coordinate planes
!
               IF (BC_X_W(BCV)/=UNDEFINED .AND. BC_Y_S(BCV)/=UNDEFINED .AND. &
                  BC_Z_B(BCV)/=UNDEFINED) CALL CHECK_PLANE (X_CONSTANT, &
                  Y_CONSTANT, Z_CONSTANT, BCV, 'BC') 
!
!  Find the value of index for the boundary cell in the direction
!  normal to the boundary plane based on one of the two values of
!  indices in the parallel direction.  For example, find whether
!  I_w or I_e should be the index for a boundary plane parallel to the
!  y-z plane based on J_s and K_b.  Then ensure that it is true for J_s
!  and K_b - K_t
!
               IF (X_CONSTANT .AND. BC_X_W(BCV)/=UNDEFINED) THEN 
                  CALL MOD_BC_I (BCV, BC_I_W(BCV), BC_I_E(BCV), BC_J_S(BCV), BC_K_B&
                     (BCV), BC_PLANE(BCV)) 
                  I_WALL = BC_I_W(BCV) 
                  IF (BC_PLANE(BCV) == 'W') THEN 
                     I_FLUID = I_WALL - 1 
                  ELSE 
                     I_FLUID = I_WALL + 1 
                  ENDIF 
                  DO K = BC_K_B(BCV), BC_K_T(BCV) 
                     DO J = BC_J_S(BCV), BC_J_N(BCV) 
!// 360 1025 Check if current i,j,k resides on this PE		     
   		       IF (.NOT.IS_ON_myPE_plus2layers(I_FLUID,J,K)) CYCLE
   		       IF (.NOT.IS_ON_myPE_plus2layers(I_WALL,J,K)) CYCLE
		       
!// 220 1004 Need to use local FUNIJK		     
                        IJK_WALL = FUNIJK(I_WALL,J,K) 
                        IJK_FLUID = FUNIJK(I_FLUID,J,K) 
                        IF (.NOT.(WALL_ICBC_FLAG(IJK_WALL) .AND. ICBC_FLAG(&
                           IJK_FLUID)(1:1)=='.')) THEN 
                           WRITE (UNIT_LOG, 1100) BCV, I_WALL, I_FLUID, J, K, &
                              IJK_WALL, ICBC_FLAG(IJK_WALL), IJK_FLUID, &
                              ICBC_FLAG(IJK_FLUID) 
                           CALL MFIX_EXIT 
                        ENDIF 
                     END DO 
                  END DO 
               ENDIF 
!
               IF (Y_CONSTANT .AND. BC_Y_S(BCV)/=UNDEFINED) THEN 
                  CALL MOD_BC_J (BCV, BC_I_W(BCV), BC_J_S(BCV), BC_J_N(BCV), BC_K_B&
                     (BCV), BC_PLANE(BCV)) 
                  J_WALL = BC_J_S(BCV) 
                  IF (BC_PLANE(BCV) == 'S') THEN 
                     J_FLUID = J_WALL - 1 
                  ELSE 
                     J_FLUID = J_WALL + 1 
                  ENDIF 
                  DO K = BC_K_B(BCV), BC_K_T(BCV) 
                     DO I = BC_I_W(BCV), BC_I_E(BCV) 
!//? Add filter
!// 360 1025 Check if current i,j,k resides on this PE		     
    		       IF (.NOT.IS_ON_myPE_plus2layers(I,J_FLUID,K)) CYCLE
    		       IF (.NOT.IS_ON_myPE_plus2layers(I,J_WALL,K)) CYCLE
!//SP
!     write(*,*) 'pass1', myPE, BCV, K, I
!// 360 Check if current k resides on this PE
!//SP
!		       if(k .ge. kstart3_all(myPE) .AND. k .le. kend3_all(myPE)) then		     
!// 220 1004 Need to use local FUNIJK		     
                        IJK_WALL = FUNIJK(I,J_WALL,K) 
                        IJK_FLUID = FUNIJK(I,J_FLUID,K) 
                        IF (.NOT.(WALL_ICBC_FLAG(IJK_WALL) .AND. ICBC_FLAG(&
                           IJK_FLUID)(1:1)=='.')) THEN 
                           WRITE (UNIT_LOG, 1200) BCV, I, J_WALL, J_FLUID, K, &
                              IJK_WALL, ICBC_FLAG(IJK_WALL), IJK_FLUID, &
                              ICBC_FLAG(IJK_FLUID) 
                           WRITE (*, 1200) BCV, I, J_WALL, J_FLUID, K, &
                              IJK_WALL, ICBC_FLAG(IJK_WALL), IJK_FLUID, &
                              ICBC_FLAG(IJK_FLUID) 
                           CALL MFIX_EXIT 
                        ENDIF 
!		       endif
                     END DO 
                  END DO 
               ENDIF 
!
               IF (Z_CONSTANT .AND. BC_Z_B(BCV)/=UNDEFINED) THEN 
                  CALL MOD_BC_K (BCV, BC_I_W(BCV), BC_J_S(BCV), BC_K_B(BCV), BC_K_T&
                     (BCV), BC_PLANE(BCV)) 
                  K_WALL = BC_K_B(BCV) 
                  IF (BC_PLANE(BCV) == 'B') THEN 
                     K_FLUID = K_WALL - 1 
                  ELSE 
                     K_FLUID = K_WALL + 1 
                  ENDIF 
                  DO J = BC_J_S(BCV), BC_J_N(BCV) 
                     DO I = BC_I_W(BCV), BC_I_E(BCV) 
!//? Add filter
!// 360 1025 Check if current i,j,k resides on this PE		     
    		       IF (.NOT.IS_ON_myPE_plus2layers(I,J,K_FLUID)) CYCLE
    		       IF (.NOT.IS_ON_myPE_plus2layers(I,J,K_WALL)) CYCLE		     
!// 220 1004 Need to use local FUNIJK		     		     
                        IJK_WALL = FUNIJK(I,J,K_WALL) 
                        IJK_FLUID = FUNIJK(I,J,K_FLUID) 
                        IF (.NOT.(WALL_ICBC_FLAG(IJK_WALL) .AND. ICBC_FLAG(&
                           IJK_FLUID)(1:1)=='.')) THEN 
                           WRITE (UNIT_LOG, 1300) BCV, I, J, K_WALL, K_FLUID, &
                              IJK_WALL, ICBC_FLAG(IJK_WALL), IJK_FLUID, &
                              ICBC_FLAG(IJK_FLUID) 
                           CALL MFIX_EXIT 
                        ENDIF 
                     END DO 
                  END DO 
               ENDIF 
!
! CHECK FOR VALID VALUES
!
               IF (BC_I_W(BCV)<1 .OR. BC_I_W(BCV)>IMAX2) GO TO 900 
               IF (BC_I_E(BCV)<1 .OR. BC_I_E(BCV)>IMAX2) GO TO 900 
               IF (BC_J_S(BCV)<1 .OR. BC_J_S(BCV)>JMAX2) GO TO 900 
               IF (BC_J_N(BCV)<1 .OR. BC_J_N(BCV)>JMAX2) GO TO 900 
               IF (BC_K_B(BCV)<1 .OR. BC_K_B(BCV)>KMAX2) GO TO 900 
               IF (BC_K_T(BCV)<1 .OR. BC_K_T(BCV)>KMAX2) GO TO 900 
               IF (BC_K_B(BCV) > BC_K_T(BCV)) GO TO 900 
               IF (BC_J_S(BCV) > BC_J_N(BCV)) GO TO 900 
               IF (BC_I_W(BCV) > BC_I_E(BCV)) GO TO 900 
!
               DO K = BC_K_B(BCV), BC_K_T(BCV) 
                  DO J = BC_J_S(BCV), BC_J_N(BCV) 
                     DO I = BC_I_W(BCV), BC_I_E(BCV) 
!//? Add filter
!// 360 1025 Check if current i,j,k resides on this PE		     
    		       IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE		     		     
!// 220 1004 Need to use local FUNIJK as dimension of ICBC_FLAG is DIM_3L		     		     
                        IJK = FUNIJK(I,J,K) 
!// 360 1104 Check if current i,j,k resides on this PE before printing error msg			
!                        IF (.NOT.WALL_ICBC_FLAG(IJK)) THEN 
                        IF (.NOT.WALL_ICBC_FLAG(IJK)) THEN 			
                           WRITE (UNIT_LOG, 1500) BCV, ICBC_FLAG(IJK), I, J, K 
                           ERROR = .TRUE. 
                        ENDIF 
                        SELECT CASE (BC_TYPE(BCV))  
                        CASE ('P_OUTFLOW')  
                           ICBC_FLAG(IJK)(1:1) = 'P' 
                        CASE ('MASS_INFLOW')  
                           ICBC_FLAG(IJK)(1:1) = 'I' 
                        CASE ('MASS_OUTFLOW')  
                           ICBC_FLAG(IJK)(1:1) = 'O' 
                        CASE ('OUTFLOW')  
                           ICBC_FLAG(IJK)(1:1) = 'o' 
                        CASE ('P_INFLOW')  
                           ICBC_FLAG(IJK)(1:1) = 'p' 
                        END SELECT 
                        BC2 = MOD(BCV,100) 
                        WRITE (ICBC_FLAG(IJK)(2:3), 950) BC2 
                     END DO 
                  END DO 
               END DO 
            ENDIF 
         ENDIF 
      END DO 


      IF (ERROR) call mfix_exit(myPE)  
!
!//AIKEPARDBG dump the ICBC_FLAG in matrix form to verify with serial version
!      DO K = Kstart3, Kend3                               !//AIKEPARDBG
!         write(UNIT_LOG,"('K = ',I5)") K                !//AIKEPARDBG 
!	 write(UNIT_LOG,"(7X,14(I3,2X))") (I,i=IMIN3,IMAX3)  !//AIKEPARDBG
!         DO J = Jstart3, Jend3                            !//AIKEPARDBG
!           write(UNIT_LOG,"(I5,')',$)") J               !//AIKEPARDBG	
!           DO I = Istart3, Iend3                          !//AIKEPARDBG
!             IJK = FUNIJK(I,J,K)                     !//AIKEPARDBG
!             write(UNIT_LOG,"(2X,A3,$)") ICBC_FLAG(IJK) !//AIKEPARDBG
!           END DO                                       !//AIKEPARDBG
!           write(UNIT_LOG,"(/)")                        !//AIKEPARDBG
!         END DO                                         !//AIKEPARDBG
!      END DO                                            !//AIKEPARDBG

      RETURN  
!
  900 CONTINUE 
  950 FORMAT(I2.2) 
!
! NOTE: WANT TO PRINT OUT ALL THE VALUES HERE ...
!
      CALL ERROR_ROUTINE ('GET_FLOW_BC', 'Invalid BC location specified', 0, 2) 
      WRITE (UNIT_LOG, *) ' BC number = ', BCV 
      WRITE (UNIT_LOG, *) ' BC_I_w(BCV) = ', BC_I_W(BCV) 
      WRITE (UNIT_LOG, *) ' BC_I_e(BCV) = ', BC_I_E(BCV) 
      WRITE (UNIT_LOG, *) ' BC_J_s(BCV) = ', BC_J_S(BCV) 
      WRITE (UNIT_LOG, *) ' BC_J_n(BCV) = ', BC_J_N(BCV) 
      WRITE (UNIT_LOG, *) ' BC_K_b(BCV) = ', BC_K_B(BCV) 
      WRITE (UNIT_LOG, *) ' BC_K_t(BCV) = ', BC_K_T(BCV) 
      CALL ERROR_ROUTINE (' ', ' ', 1, 3) 
      RETURN  
 1100 FORMAT(/70('*')//' From: GET_FLOW_BC'/'Message: Boundary ','condition ',&
         I3,' has illegal geometry',/,'I_wall  = ',I3,/,'I_fluid = ',I3,/,&
         'J       = ',I3,/,'K       = ',I3,/,'ICBC_FLAG(',I4,') = ',A3,&
         'This should be a wall cell',/,'ICBC_FLAG(',I4,') = ',A3,&
         'This should be a fluid cell',/70('*')/) 
 1200 FORMAT(/70('*')//' From: GET_FLOW_BC'/'Message: Boundary ','condition ',&
         I3,' has illegal geometry',/,'I       = ',I3,/,'J_wall  = ',I3,/,&
         'J_fluid = ',I3,/,'K       = ',I3,/,'ICBC_FLAG(',I4,') = ',A3,&
         'This should be a wall cell',/,'ICBC_FLAG(',I4,') = ',A3,&
         'This should be a fluid cell',/70('*')/) 
 1300 FORMAT(/70('*')//' From: GET_FLOW_BC'/'Message: Boundary ','condition ',&
         I3,' has illegal geometry',/,'I       = ',I3,/,'J       = ',I3,/,&
         'K_wall  = ',I3,/,'K_fluid = ',I3,/,'ICBC_FLAG(',I4,') = ',A3,&
         'This should be a wall cell',/,'ICBC_FLAG(',I4,') = ',A3,&
         'This should be a fluid cell',/70('*')/) 
 1500 FORMAT(/70('*')//' From: GET_FLOW_BC'/'Message: Boundary ','condition ',&
         I3,' overlaps boundary condition ',A3,' at'/' I = ',I3,'  J = ',I3,&
         '  K = ',I3/70('*')/) 
      END SUBROUTINE GET_FLOW_BC 
