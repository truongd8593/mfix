!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: MOD_BC_I                                                C
!  Purpose: modify the "I" values for the b.c. plane                   C
!     This is a yz plane                                               C
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
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE MOD_BC_I(BC, I_W, I_E, J_S, K_B, PLANE) 

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param 
      USE param1 
      USE geometry
      USE fldvar
      USE physprop
      USE indices
      USE funits 
      USE compar 
      USE mpi_utility 
      IMPLICIT NONE

!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! boundary condition index
      INTEGER, INTENT(IN) :: BC
! i cell indices defining location of yz plane
      INTEGER, INTENT(INOUT) :: I_w, I_e
! south/bottom j,k cell indices of yz plane
      INTEGER, INTENT(IN) :: J_s, K_b
! the flow surface plane
      CHARACTER, INTENT(OUT) :: PLANE

!-----------------------------------------------
! Local variables
!-----------------------------------------------
!     'IJK' indices
      INTEGER  IJK1 , IJK, bcast_root

!-----------------------------------------------
! Include statement functions
!-----------------------------------------------
    INCLUDE 'function.inc'
!-----------------------------------------------

      IF(IS_ON_myPE_owns(I_W,J_S,K_B)) then
         bcast_root = myPE
         call global_all_sum(bcast_root)
      ELSE
         bcast_root = 0
         call global_all_sum(bcast_root)
      ENDIF

      IF(IS_ON_myPE_owns(I_W,J_S,K_B)) THEN
         IJK1 = FUNIJK(I_W,J_S,K_B)
         IJK = FUNIJK(I_W + 1,J_S,K_B) 
         IF (WALL_ICBC_FLAG(IJK1) .AND. ICBC_FLAG(IJK)(1:1)=='.') THEN 
! flow in/out on west boundary (fluid cell on east)
            I_W = I_W 
            I_E = I_E 
            PLANE = 'E' 
         ELSEIF (WALL_ICBC_FLAG(IJK) .AND. ICBC_FLAG(IJK1)(1:1)=='.') THEN 
! flow in/out on east boundary (fluid cell on west)
            I_W = I_W + 1 
            I_E = I_E + 1 
            PLANE = 'W' 
         ELSE 
            IF(DMP_LOG)WRITE (UNIT_LOG, 1000) BC, I_W, I_E, J_S, K_B, ICBC_FLAG(IJK1), &
               ICBC_FLAG(IJK) 
            call mfix_exit(myPE)  
         ENDIF 
      ENDIF

      CALL bcast(I_W,bcast_root)
      CALL bcast(I_E,bcast_root)
      CALL bcast(PLANE,bcast_root)

      RETURN  
 1000 FORMAT(/70('*')//'From: MOD_BC_I'/'Message: Cannot locate the ',&
         'flow plane for boundary condition ',I3,/' I West   = ',I3,/&
         ' I East   = ',I3,/' J South  = ',I3,/' K Bottom = ',I3,/&
         ' One of the following should be a wall cell and the other a',&
         ' fluid cell:',/5X,A3,5X,A3,/&
         ' May be no IC was specified for the fluid cell.',/70('*')/) 
      END SUBROUTINE MOD_BC_I 
      
