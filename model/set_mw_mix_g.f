!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SET_MW_MIX_g                                           C
!  Purpose: calculate gas mixture molecular weights                    C
!                                                                      C
!  Author: M. Syamlal                                 Date: 19-OCT-92  C
!  Reviewer: S. Venkatesan                            Date: 11-DEC-92  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: IJKMAX2, X_g                                  C
!                                                                      C
!  Variables modified: MW_MIX_g                                        C
!                                                                      C
!  Local variables: NONE                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE SET_MW_MIX_G 
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
      DOUBLE PRECISION , EXTERNAL :: CALC_MW 
!-----------------------------------------------
      INCLUDE 'function.inc'
!
      IF (MW_AVG /= UNDEFINED) RETURN  

!!$omp parallel do private(ijk) &
!!$omp schedule(dynamic,chunk_size)

      DO IJK = IJKMIN1, IJKMAX1 
         IF (FLUID_AT(IJK)) MW_MIX_G(IJK) = CALC_MW(X_G,DIMENSION_3,IJK,NMAX(0)&
            ,MW_G) 
      END DO 
      RETURN  
      END SUBROUTINE SET_MW_MIX_G 
