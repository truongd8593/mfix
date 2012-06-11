!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  SUBROUTINE: TRANSPORT_PROP                                          C
!  Purpose: Calculate the indicated transport properties that vary     C
!           with time if directed to do so by the corresponding flag   C
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
!  Local variables: None                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE TRANSPORT_PROP(VISC, COND, DIFF, GRAN_DISS, IER) 

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param 
      USE param1 
      USE fldvar
      USE physprop
      USE geometry
      USE indices
      USE run
      USE toleranc 
      USE compar
      USE discretelement
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Flags to tell whether to calculate or not
      LOGICAL, INTENT(IN) ::  VISC(0:DIMENSION_M), &
                              COND(0:DIMENSION_M), &
                              DIFF(0:DIMENSION_M)
! Flags to tell whether to calculate or not
      LOGICAL, INTENT(IN) :: GRAN_DISS(0:DIMENSION_M)
!  Error index
      INTEGER, INTENT(INOUT) :: IER
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Phase index
      INTEGER :: M                
!-----------------------------------------------

! Fluid viscosity
      IF (VISC(0)) CALL CALC_MU_G (IER) 

! Fluid conductivity
      IF (COND(0)) CALL CALC_K_G (IER) 

! Fluid diffusivity
      IF (DIFF(0)) CALL CALC_DIF_G (IER) 

! note that the code could be changed so that the logical flags
! directing this routine are appropriately set for the discrete_element
! case
      IF (.NOT.DISCRETE_ELEMENT .OR. DES_CONTINUUM_HYBRID) THEN
         DO M = 1, MMAX 
! Solids conductivity
            IF (COND(M)) CALL CALC_K_S (M, IER) 

! Solids viscosity
            IF (VISC(M)) CALL CALC_MU_S (M, IER) 

! Solids diffusivity
            IF (DIFF(M)) CALL CALC_DIF_S (M, IER) 

! Particle-Particle Energy Dissipation
            IF (GRAN_DISS(M)) THEN
               IF (TRIM(KT_TYPE) .EQ. 'IA_NONEP') THEN
                  CALL CALC_IA_NONEP_ENERGY_DISSIPATION_SS(M, IER)
               ELSEIF (TRIM(KT_TYPE) .EQ. 'GD_99') THEN
                  CALL CALC_GD_99_ENERGY_DISSIPATION_SS(M, IER)
               ENDIF
            ENDIF
         ENDDO   ! end do (m=1,mmax)
      ENDIF   ! end if (.not.discrete_element .or des_continuum_hybrid)
      
      RETURN  
      END SUBROUTINE TRANSPORT_PROP 

