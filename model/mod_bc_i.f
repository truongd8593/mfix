!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: MOD_BC_I(BC, I_w, I_e, J_s, K_b, PLANE)                C
!  Purpose: modify the "I" values for the b.c. plane                   C
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
!  Variables referenced: ICBC_FLAG                                     C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables: IJK1, IJK                                         C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE MOD_BC_I(BC, I_W, I_E, J_S, K_B, PLANE) 
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
      INTEGER  I_w, I_e, J_s, K_b
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
      IJK = FUNIJK(I_W + 1,J_S,K_B) 
      IF (WALL_ICBC_FLAG(IJK1) .AND. ICBC_FLAG(IJK)(1:1)=='.') THEN 
         I_W = I_W 
         I_E = I_E 
         PLANE = 'E' 
      ELSE IF (WALL_ICBC_FLAG(IJK) .AND. ICBC_FLAG(IJK1)(1:1)=='.') THEN 
         I_W = I_W + 1 
         I_E = I_E + 1 
         PLANE = 'W' 
      ELSE 
         WRITE (UNIT_LOG, 1000) BC, I_W, I_E, J_S, K_B, ICBC_FLAG(IJK1), &
            ICBC_FLAG(IJK) 
         call mfix_exit(myPE)  
      ENDIF 
    ENDIF

!/SP
      CALL bcast(I_W,bcast_root)
      CALL bcast(I_E,bcast_root)
      CALL bcast(PLANE,bcast_root)
!
      RETURN  
 1000 FORMAT(/70('*')//'From: MOD_BC_I'/'Message: Cannot locate the ',&
         'flow plane for boundary condition ',I3,/' I West   = ',I3,/&
         ' I East   = ',I3,/' J South  = ',I3,/' K Bottom = ',I3,/&
         ' One of the following should be a wall cell and the other a',&
         ' fluid cell:',/5X,A3,5X,A3,/&
         ' May be no IC was specified for the fluid cell.',/70('*')/) 
      END SUBROUTINE MOD_BC_I 
