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

         DO IJK = 1, IJKMAX2 
            IF (.NOT.WALL_AT(IJK)) THEN 
               RO_G(IJK) = EOSG(MW_MIX_G(IJK),P_G(IJK),T_G(IJK)) 
               ROP_G(IJK) = EP_G(IJK)*RO_G(IJK) 
            ENDIF 
         END DO 
      ELSE 

!!$omp   parallel do private(ijk)  
         DO IJK = 1, IJKMAX2 
            IF (.NOT.WALL_AT(IJK)) THEN 
               RO_G(IJK) = RO_G0 
               ROP_G(IJK) = EP_G(IJK)*RO_G(IJK) 
            ENDIF 
         END DO 
      ENDIF 
!
      RETURN  
      END SUBROUTINE SET_RO_G 
