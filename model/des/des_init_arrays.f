!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                         C
!     Subrourtine: DES_INIT_ARRAYS                                        C
!     Purpose: initialize arrays from des_allocate arrays                 C
!                                                                         C
!     Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!     Reviewer:                                          Date:            C
!                                                                         C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE DES_INIT_ARRAYS

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE discretelement
      USE indices
      USE geometry
      USE compar
      USE physprop
      USE des_bc
      USE run
      use desgrid
      use desmpi
      USE des_thermo
      USE des_rxns

      IMPLICIT NONE

      INTEGER :: II
!-----------------------------------------------

! Pradeep: parallel processing
      iglobal_id = 0

! particle properties
      DES_RADIUS(:) = ZERO
      PMASS(:) = ZERO
      PVOL(:) = ZERO
      OMOI(:) = ZERO
      RO_Sol(:) = ZERO

! particle position, velocity, etc

      DES_POS_NEW(:,:) = ZERO
      DES_VEL_NEW(:,:) = ZERO
      OMEGA_NEW(:,:) = ZERO

      IF (DO_OLD) THEN
         DES_POS_OLD(:,:) = ZERO
         DES_VEL_OLD(:,:) = ZERO
         DES_ACC_OLD(:,:) = ZERO
         OMEGA_OLD(:,:) = ZERO
         ROT_ACC_OLD(:,:) = ZERO
      ENDIF

! Initializing user defined array (Surya Dec 10, 2014)
      DES_USR_VAR(:,:) = ZERO   

      FC(:,:) = ZERO
      TOW(:,:) = ZERO

      PPOS(:,:) = ZERO

      PINC(:) = ZERO
      PIJK(:,:) = ZERO

      DES_U_s(:,:) = ZERO
      DES_V_s(:,:) = ZERO
      DES_W_s(:,:) = ZERO
      DES_ROP_S(:,:) = ZERO
      DES_ROP_SO(:,:) = ZERO

      P_FORCE(:,:) = ZERO

      IF(allocated(DRAG_AM)) DRAG_AM = ZERO
      IF(allocated(DRAG_BM)) DRAG_BM = ZERO

      F_GDS = ZERO
      VXF_GDS = ZERO

      IF (DES_CONTINUUM_HYBRID) THEN
         F_SDS = ZERO
         VXF_SDS = ZERO
         SDRAG_AM = ZERO
         SDRAG_BM = ZERO
      ENDIF

      GRAV(:) = ZERO

      DO II = 1, SIZE(particle_wall_collisions)
         nullify(particle_wall_collisions(II)%pp)
      ENDDO

! Cohesion VDW forces
      IF(USE_COHESION) THEN
         PostCohesive (:) = ZERO
      ENDIF

! J.Musser: DEM particle tracking quantity
      PEA(:,:) = .FALSE.

! J.Musser: Energy and Species Equation Arrays
      IF(ENERGY_EQ)THEN
         DES_T_s_OLD(:) = ZERO  !Changed from Undefined to Zero to avoid mpi problems (Surya Oct 29, 2014)
         DES_T_s_NEW(:) = ZERO  !Changed from Undefined to Zero to avoid mpi problems (Surya OCt 29, 2014)
!        DES_C_PS(:) = UNDEFINED
         DES_C_PS(:) = ZERO
         DES_X_s(:,:) = ZERO
         Q_Source(:) = ZERO
         avgDES_T_s(:) = ZERO
         DES_ENERGY_SOURCE(:) = ZERO
         IF (INTG_ADAMS_BASHFORTH) &
            Q_Source0(:) = ZERO
      ENDIF

      IF(ANY_SPECIES_EQ)THEN
         DES_R_sp(:,:) = ZERO
         DES_R_sc(:,:) = ZERO
         Qint(:) = ZERO
         IF (INTG_ADAMS_BASHFORTH) THEN
            dMdt_OLD(:) = ZERO
            dXdt_OLD(:,:) = ZERO
         ENDIF
      ENDIF

      RETURN
      END SUBROUTINE DES_INIT_ARRAYS


