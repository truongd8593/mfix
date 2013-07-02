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
      SUBROUTINE TRANSPORT_PROP(IER)

! Global Variables:
!----------------------------------------------------------------------
! Number of solids phases.
      use physprop, only: MMAX
! Flags for calculating viscosity.
      use coeff, only: VISC
! Flags for calculating conductivity.
      use coeff, only: COND
! Flags for calculating diffusivity.
      use coeff, only: DIFF
! Flags for calculating particle-particle energy dissipation.
      use coeff, only: GRAN_DISS
! Kinetic theory model.
      use run, only: KT_TYPE

      implicit none

! Dummy arguments
!-----------------------------------------------------------------------
!  Error index
      INTEGER, intent(inout) :: IER

! Local variables
!-----------------------------------------------------------------------
! Loop counter
      INTEGER :: M ! Solids phase


      IF (VISC(0)) CALL CALC_MU_G (IER)    ! Fluid viscosity
      IF (COND(0)) CALL CALC_K_G (IER)     ! Fluid conductivity
      IF (DIFF(0)) CALL CALC_DIF_G (IER)   ! Fluid diffusivity

      DO M = 1, MMAX
         IF (COND(M)) CALL CALC_K_S (M, IER)   ! Solids conductivity
         IF (VISC(M)) CALL CALC_MU_S (M, IER)  ! Solids viscosity
         IF (DIFF(M)) CALL CALC_DIF_S (M, IER) ! Solids diffusivity

! Particle-Particle Energy Dissipation
         IF (GRAN_DISS(M)) THEN
            IF (TRIM(KT_TYPE) .EQ. 'IA_NONEP') THEN
               CALL CALC_IA_NONEP_ENERGY_DISSIPATION_SS(M, IER)
            ELSEIF (TRIM(KT_TYPE) .EQ. 'GD_99') THEN
               CALL CALC_GD_99_ENERGY_DISSIPATION_SS(M, IER)
            ENDIF
         ENDIF
      ENDDO

      RETURN  
      END SUBROUTINE TRANSPORT_PROP
