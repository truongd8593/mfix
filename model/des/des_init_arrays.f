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
!-----------------------------------------------

! Pradeep: parallel processing
      iglobal_id = 0

! T.Li: Hertzian collision model
      g_mod(:) = zero
      hert_kn(:,:) = zero
      hert_kwn(:) = zero
      hert_kt(:,:) = zero
      hert_kwt(:) = zero

! particle properties
      DES_RADIUS(:) = ZERO
      PMASS(:) = ZERO
      PVOL(:) = ZERO
      OMOI(:) = ZERO
      RO_Sol(:) = ZERO

! particle position, velocity, etc
      DES_POS_OLD(:,:) = ZERO
      DES_POS_NEW(:,:) = ZERO
      DES_VEL_OLD(:,:) = ZERO
      DES_VEL_NEW(:,:) = ZERO

      DES_ACC_OLD(:,:) = ZERO

      OMEGA_OLD(:,:) = ZERO
      OMEGA_NEW(:,:) = ZERO
      ROT_ACC_OLD(:,:) = ZERO

      FC(:,:) = ZERO
      TOW(:,:) = ZERO

      PN(:,:) = -1
      PN_WALL(:,:) = -1
      PN(1,:) = 0
      PN_WALL(1,:) = 0
      PV(:,:) = .TRUE.
      PV_WALL(:,:) = .TRUE.
      PFT(:,:,:) = ZERO
      PFT_WALL(:,:,:) = ZERO
      PPOS(:,:) = ZERO

      PINC(:) = ZERO
      PIJK(:,:) = ZERO

      DES_U_s(:,:) = ZERO
      DES_V_s(:,:) = ZERO
      DES_W_s(:,:) = ZERO
      DES_ROP_S(:,:) = ZERO
      DES_ROP_SO(:,:) = ZERO

      P_FORCE(:,:) = ZERO

      IF (DES_INTERP_ON) THEN
         DRAG_AM(:,:) = ZERO
         DRAG_BM(:,:,:) = ZERO
      ENDIF

      IF (DES_CONTINUUM_HYBRID) THEN
         F_GDS(:,:) = ZERO
         F_SDS(:,:,:) = ZERO
         VXF_GDS(:,:) = ZERO
         VXF_SDS(:,:,:) = ZERO
      ENDIF



      GRAV(:) = ZERO

      NEIGHBOURS(:,:) = -1
      NEIGHBOURS(:,1) = 0


! Cohesion VDW forces
      IF(USE_COHESION) THEN
         Fcohesive(:,:) = ZERO
         PostCohesive (:) = ZERO
      ENDIF

      IF(DES_CALC_CLUSTER) THEN
         InACluster(:) = .FALSE.
         PostCluster(:) = ZERO
      ENDIF

! J.Musser: DEM particle tracking quantity
      PEA(:,:) = .FALSE.

! J.Musser: Energy and Species Equation Arrays
      IF(ENERGY_EQ)THEN
         DES_T_s_OLD(:) = UNDEFINED
         DES_T_s_NEW(:) = UNDEFINED
         DES_C_PS(:) = UNDEFINED
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


