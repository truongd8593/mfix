!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: TRANSPORT_PROP(VISC, COND, DIFF, IER)                  C
!  Purpose: Calculate transport properties that vary with time         C
!                                                                      C
!  Author: M. Syamlal                                 Date: 17-JUL-92  C
!  Reviewer: P. Nicoletti                             Date: 11-DEC-92  C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: Mods for MFIX 2.0 (old name CALC_PHYSPROP)                 C
!  Author: M. Syamlal                                 Date: 23-APR-96  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: None                                          C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables: None                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE TRANSPORT_PROP(VISC, COND, DIFF, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE fldvar
      USE physprop
      USE geometry
      USE indices
      USE run
      USE toleranc 
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                       Error index
      INTEGER          IER
!
!                      Phase index
      INTEGER          M
!
!                      Flags to tell whether to calculate or not
      LOGICAL          VISC(0:DIMENSION_M), COND(0:DIMENSION_M), &
                       DIFF(0:DIMENSION_M)
!                 
!-----------------------------------------------
!
! 1.2.1 Fluid viscosity
      IF (VISC(0)) CALL CALC_MU_G (IER) 
!
! 1.2.2 Solids viscosity
      DO M = 1, MMAX 
         IF (VISC(M)) CALL CALC_MU_S (M, IER) 
      END DO 
!
! 1.2.3 Fluid conductivity
      IF (COND(0)) CALL CALC_K_G (IER) 
!
! 1.2.4 Solids conductivity
      DO M = 1, MMAX 
         IF (COND(M)) CALL CALC_K_S (M, IER) 
      END DO 
!
! 1.2.5 Fluid diffusivity
      IF (DIFF(0)) CALL CALC_DIF_G (IER) 
!
! 1.2.6 Solids diffusivity
      DO M = 1, MMAX 
         IF (DIFF(M)) CALL CALC_DIF_S (M, IER) 
      END DO 
      RETURN  
      END SUBROUTINE TRANSPORT_PROP 
