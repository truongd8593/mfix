!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SET_RO_g                                               C
!  Purpose: Initialize the gas densities                               C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-JAN-92  C
!  Reviewer:M. Syamlal, S. Venkatesan, P. Nicoletti,  Date: 29-JAN-92  C
!           W. Rogers                                                  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: IJKMAX2, MW_AVG, P_g, T_g, EP_g, RO_g         C
!                                                                      C
!  Variables modified: RO_g, ROP_g, IJK                                C
!                                                                      C
!  Local variables: NONE                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE SET_RO_G 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE physprop
      USE geometry
      USE fldvar
      USE constant
      USE indices
      USE compar        !//d
      USE funits        !//AIKEPARDBG
      USE sendrecv      !// 400
!      USE dbg_util      !//AIKEPARDBG
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
      INTEGER :: IJK
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      DOUBLE PRECISION , EXTERNAL :: EOSG 
!-----------------------------------------------
      INCLUDE 'function.inc'
!
      IF (RO_G0 == UNDEFINED) THEN 

!!$omp parallel do private(IJK)  

!//? 1024 make sure the values of the dependent variables are consistent at the
!//? overlapping subdomain boundaries, need debug print during testing.

!//? make sure the values in ghosts are set appropriately to avoid the division
!//? by zero problem in function EOSG

!// 350 1025 change do loop limits: 1,ijkmax2-> ijkstart3, ijkend3    
         DO IJK = ijkstart3, ijkend3 
            IF (.NOT.WALL_AT(IJK)) THEN 
               RO_G(IJK) = EOSG(MW_MIX_G(IJK),P_G(IJK),T_G(IJK)) 
               ROP_G(IJK) = EP_G(IJK)*RO_G(IJK) 
            ENDIF 
         END DO 
      ELSE 

!!$omp   parallel do private(ijk)  
!// 350 1025 change do loop limits: 1,ijkmax2-> ijkstart3, ijkend3    
         DO IJK = ijkstart3, ijkend3 
            IF (.NOT.WALL_AT(IJK)) THEN 
               RO_G(IJK) = RO_G0 
               ROP_G(IJK) = EP_G(IJK)*RO_G(IJK) 
            ENDIF 
         END DO 
      ENDIF 


!       call prnfield(RO_G,'RO_G','BEF')   !//AIKEPARDBG

!
!//S No communication is necessary as the RO_G and ROP_G are calculated based
!//  on already available vars and these portions are already overlapping.
!//  A good check for the execution is to see whether RO_G and ROP_G have
!//  same values on both PEs although they were calculated independently.
!//? Mike's implementation but check necessity and also syntax as it give compilation error
!	CALL SEND_RECV(RO_G, 2)
!	CALL SEND_RECV(ROP_G, 2)
	

!       call prnfield(RO_G,'RO_G','BEF')    !//AIKEPARDBG

      RETURN  
      END SUBROUTINE SET_RO_G 
