!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: MOD_BC_K (BC, I_w, J_s, K_b, K_t, PLANE)               C
!  Purpose: modify the "K" values for the b.c. plane                   C
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
!  Variables referenced: None                                          C
!  Variables modified: ICBC_FLAG                                       C
!                                                                      C
!  Local variables: IJK1, IJK                                         C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE MOD_BC_K(BC, I_W, J_S, K_B, K_T, PLANE) 
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
      USE indices
      USE funits 
      USE compar        !//d
      USE mpi_utility   !//d
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!              boundary condition index
      INTEGER  BC
!
!              calculated cell indices in I,J,K directions
      INTEGER  I_w, J_s, K_b, K_t
!
!               the flow surface plane
      CHARACTER PLANE
!
! local variables
!
!              'IJK' indices
      INTEGER  IJK1 , IJK, bcast_root
!-----------------------------------------------
      INCLUDE 'function.inc'
!
!// 220 1004 Replaced with global FUNIJK
!// SP
    IF(IS_ON_myPE_owns(I_W,J_S,K_B)) then
      bcast_root = myPE
      call global_all_sum(bcast_root,bcast_root)
    ELSE
      bcast_root = 0
      call global_all_sum(bcast_root,bcast_root)
    ENDIF

    IF(IS_ON_myPE_owns(I_W,J_S,K_B)) then
      IJK1 = FUNIJK(I_W,J_S,K_B) 
      IJK = FUNIJK(I_W,J_S,K_B + 1) 
      
      IF (WALL_ICBC_FLAG(IJK1) .AND. ICBC_FLAG(IJK)(1:1)=='.') THEN 
         K_B = K_B 
         K_T = K_T 
         PLANE = 'T' 
      ELSE IF (WALL_ICBC_FLAG(IJK) .AND. ICBC_FLAG(IJK1)(1:1)=='.') THEN 
         K_B = K_B + 1 
         K_T = K_T + 1 
         PLANE = 'B' 
      ELSE 
         WRITE (UNIT_LOG, 1000) BC, I_W, J_S, K_B, K_T, ICBC_FLAG(IJK1), &
            ICBC_FLAG(IJK) 
         call mfix_exit(myPE)  
      ENDIF 
    ENDIF

!// SP
      CALL bcast(K_B,bcast_root)
      CALL bcast(K_T,bcast_root)
      CALL bcast(PLANE,bcast_root)
!
      RETURN  
 1000 FORMAT(/70('*')//'From: MOD_BC_K'/'Message: Cannot locate the ',&
         'flow plane for boundary condition ',I3,/' I West   = ',I3,/&
         ' J South  = ',I3,/' K Bottom = ',I3,/' K Top    = ',I3,/&
         ' One of the following should be a wall cell and the other a',&
         ' fluid cell:',/5X,A3,5X,A3,/&
         ' May be no IC was specified for the fluid cell.',/70('*')/) 
      END SUBROUTINE MOD_BC_K 
