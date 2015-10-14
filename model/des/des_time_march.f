!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!     Subroutine: DES_TIME_MARCH                                       !
!     Author: Jay Boyalakuntla                        Date: 21-Jun-04  !
!                                                                      !
!     Purpose: Main DEM driver routine                                 !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DES_TIME_MARCH

      use des_bc, only: DEM_BCMI, DEM_BCMO
      use des_thermo, only: DES_ENERGY_SOURCE
      use desgrid, only: desgrid_pic
      use discretelement
      use error_manager
      use fldvar, only: EP_g, ROP_g, ROP_s
      use fldvar, only: EP_g, ROP_g, ROP_s
      use functions
      use machine
      use mpi_funs_des, only: DES_PAR_EXCHANGE
      use mpi_utility
      use run, only: ANY_SPECIES_EQ
      use run, only: ANY_SPECIES_EQ
      use run, only: CALL_USR
      use run, only: ENERGY_EQ
      use run, only: NSTEP
      use run, only: TIME, TSTOP, DT
      use sendrecv
      use sort
        USE des_rxns
        USE des_thermo
        USE mfix_pic
        USE discretelement
        USE particle_filter
        USE run


      IMPLICIT NONE
!------------------------------------------------
! Local variables
!------------------------------------------------
! Total number of particles
      INTEGER, SAVE :: NP=0

! time step loop counter index
      INTEGER :: NN

! loop counter index for any initial particle settling incoupled cases
      INTEGER :: FACTOR,ii
      INTEGER, dimension(:), allocatable :: tmp_int

! Temporary variables when des_continuum_coupled is T to track
! changes in solid time step
      DOUBLE PRECISION :: TMP_DTS, DTSOLID_TMP

! Numbers to calculate wall time spent in DEM calculations.
      DOUBLE PRECISION :: TMP_WALL

           DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: tmp_des_radius
           DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: tmp_RO_Sol
           DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: tmp_PVOL
           DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: tmp_PMASS
           DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: tmp_OMOI
           DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: tmp_DES_POS_NEW
           DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: tmp_DES_VEL_NEW
           DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: tmp_OMEGA_NEW
           DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: tmp_PPOS
           !call byte_grow(PARTICLE_STATE
           INTEGER, ALLOCATABLE, DIMENSION(:) :: tmp_PARTICLE_STATE
           INTEGER, ALLOCATABLE, DIMENSION(:) :: tmp_orig_index
           INTEGER, ALLOCATABLE, DIMENSION(:) :: tmp_iglobal_id
           INTEGER, ALLOCATABLE, DIMENSION(:,:) :: tmp_pijk
           INTEGER, ALLOCATABLE, DIMENSION(:) :: tmp_dg_pijk
           INTEGER, ALLOCATABLE, DIMENSION(:) :: tmp_dg_pijkprv
           LOGICAL, ALLOCATABLE, DIMENSION(:) :: tmp_ighost_updated
           DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: tmp_FC
           DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: tmp_TOW
           DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: tmp_F_GP
           INTEGER, ALLOCATABLE, DIMENSION(:,:) :: tmp_WALL_COLLISION_FACET_ID
           DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: tmp_WALL_COLLISION_PFT
           DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: tmp_DRAG_FC

           INTEGER, ALLOCATABLE, DIMENSION(:,:) :: tmp_NEIGHBOR_INDEX
           INTEGER, ALLOCATABLE, DIMENSION(:,:) :: tmp_NEIGHBOR_INDEX_OLD

           DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: tmp_ORIENTATION

              INTEGER, ALLOCATABLE, DIMENSION(:,:) :: tmp_FILTER_CELL
              DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: tmp_FILTER_WEIGHT

              DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: tmp_DES_STAT_WT
              DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: tmp_PS_GRAD
              DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: tmp_AVGSOLVEL_P
              DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: tmp_EPG_P

              DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: tmp_PostCohesive

              DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: tmp_DES_POS_OLD
              DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: tmp_DES_VEL_OLD
              DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: tmp_DES_ACC_OLD
              DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: tmp_OMEGA_OLD
              DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: tmp_ROT_ACC_OLD

              DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: tmp_DES_T_s_OLD
              DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: tmp_DES_T_s_NEW
              DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: tmp_DES_C_PS
              DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: tmp_DES_X_s
              DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: tmp_Q_Source

                   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: tmp_Q_Source0

              DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: tmp_DES_R_sp
              DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: tmp_DES_R_sc

                 DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: tmp_dMdt_OLD
                 DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: tmp_dXdt_OLD
              DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: tmp_Qint
              DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: tmp_DES_USR_VAR

! In case of restarts assign S_TIME from MFIX TIME
      S_TIME = TIME
      TMP_DTS = ZERO
      DTSOLID_TMP = ZERO
      TMP_WALL = WALL_TIME()

! Initialize time stepping variables for coupled gas/solids simulations.
      IF(DES_CONTINUUM_COUPLED) THEN
         IF(DT.GE.DTSOLID) THEN
            FACTOR = CEILING(real(DT/DTSOLID))
         ELSE
            FACTOR = 1
            DTSOLID_TMP = DTSOLID
            DTSOLID = DT
         ENDIF

! Initialize time stepping variable for pure granular simulations.
      ELSE
         FACTOR = CEILING(real((TSTOP-TIME)/DTSOLID))
         DT = DTSOLID
         CALL OUTPUT_MANAGER(.FALSE., .FALSE.)
      ENDIF   ! end if/else (des_continuum_coupled)

      NP = PIP - IGHOST_CNT
      CALL GLOBAL_ALL_SUM(NP)

      WRITE(ERR_MSG, 1000) trim(iVal(factor)), trim(iVAL(NP))
      CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE., LOG=.FALSE.)

 1000 FORMAT(/'DEM NITs: ',A,3x,'Total PIP: ', A)

      IF(CALL_USR) CALL USR0_DES

      IF(DES_CONTINUUM_COUPLED) THEN
         IF(DES_EXPLICITLY_COUPLED) THEN
            CALL DRAG_GS_DES1
         ELSE
            IF(ANY_SPECIES_EQ) CALL ZERO_RRATE_DES
            IF(ENERGY_EQ) CALL ZERO_ENERGY_SOURCE
         ENDIF
         CALL CALC_PG_GRAD
      ENDIF


! Main DEM time loop
!----------------------------------------------------------------->>>
      DO NN = 1, FACTOR

         IF(DES_CONTINUUM_COUPLED) THEN
! If the current time in the discrete loop exceeds the current time in
! the continuum simulation, exit the discrete loop
            IF(S_TIME.GE.(TIME+DT)) EXIT
! If next time step in the discrete loop will exceed the current time
! in the continuum simulation, modify the discrete time step so final
! time will match
            IF((S_TIME+DTSOLID).GT.(TIME+DT)) THEN
               TMP_DTS = DTSOLID
               DTSOLID = TIME + DT - S_TIME
            ENDIF
         ENDIF

! Calculate forces acting on particles (collisions, drag, etc).
         CALL CALC_FORCE_DEM
! Calculate or distribute fluid-particle drag force.
         CALL CALC_DRAG_DES

! Update the old values of particle position and velocity with the new
! values computed
         IF (DO_OLD) CALL CFUPDATEOLD
! Calculate thermochemical sources (energy and  rates of formation).
         CALL CALC_THERMO_DES
! Call user functions.
         IF(CALL_USR) CALL USR1_DES
! Update position and velocities
         CALL CFNEWVALUES
! Update particle temperatures
         CALL DES_THERMO_NEWVALUES
! Update particle from reactive chemistry process.
         CALL DES_REACTION_MODEL

! Set DO_NSEARCH before calling DES_PAR_EXCHANGE.
         DO_NSEARCH = (NN == 1 .OR. MOD(NN,NEIGHBOR_SEARCH_N) == 0)

! Add/Remove particles to the system via flow BCs.
         IF(DEM_BCMI > 0) CALL MASS_INFLOW_DEM
         IF(DEM_BCMO > 0) CALL MASS_OUTFLOW_DEM(DO_NSEARCH)

! Call exchange particles - this will exchange particle crossing
! boundaries as well as updates ghost particles information
         IF (DO_NSEARCH .OR. (numPEs>1) .OR. DES_PERIODIC_WALLS) THEN
            CALL DESGRID_PIC(.TRUE.)
            CALL DES_PAR_EXCHANGE
         ENDIF

         IF(DO_NSEARCH) THEN

            do ii=1,max_pip
               orig_index(ii) = ii
            enddo

            allocate(tmp_des_radius(size(des_radius)))
            allocate(tmp_RO_Sol(size(RO_Sol)))
            allocate(tmp_PVOL(size(PVOL)))
            allocate(tmp_PMASS(size(PMASS)))
            allocate(tmp_OMOI(size(OMOI)))
            allocate(tmp_DES_POS_NEW(size(DES_POS_NEW),3))
            allocate(tmp_DES_VEL_NEW(size(DES_VEL_NEW),3))
            allocate(tmp_OMEGA_NEW(size(OMEGA_NEW),3))
            allocate(tmp_PPOS(size(PPOS),3))
            allocate(tmp_PARTICLE_STATE(size(PARTICLE_STATE)))
            allocate(tmp_orig_index(size(orig_index)))
            allocate(tmp_iglobal_id(size(iglobal_id)))
            allocate(tmp_pijk(size(pijk),5))
            allocate(tmp_dg_pijk(size(dg_pijk)))
            allocate(tmp_dg_pijkprv(size(dg_pijkprv)))
            allocate(tmp_ighost_updated(size(ighost_updated)))
            allocate(tmp_FC(size(FC),3))
            allocate(tmp_TOW(size(TOW),3))
            allocate(tmp_F_GP(size(F_GP)))
            allocate(tmp_WALL_COLLISION_FACET_ID(Collision_Array_Max,size(WALL_COLLISION_FACET_ID)))
            allocate(tmp_WALL_COLLISION_PFT(DIMN,COLLISION_ARRAY_MAX,size(WALL_COLLISION_PFT)))
            allocate(tmp_DRAG_FC(size(DRAG_FC),3))
            allocate(tmp_NEIGHBOR_INDEX(2,size(NEIGHBOR_INDEX)))
            allocate(tmp_NEIGHBOR_INDEX_OLD(2,size(NEIGHBOR_INDEX_OLD)))

           IF(PARTICLE_ORIENTATION) allocate(tmp_ORIENTATION(DIMN,size(ORIENTATION)))

           IF(FILTER_SIZE > 0) THEN
             allocate(tmp_FILTER_CELL(FILTER_SIZE,size(FILTER_CELL)))
            allocate(tmp_FILTER_WEIGHT(FILTER_SIZE,size(FILTER_WEIGHT)))
           ENDIF

           IF(MPPIC) THEN
             allocate(tmp_DES_STAT_WT(size(DES_STAT_WT)))
            allocate(tmp_PS_GRAD(size(PS_GRAD),DIMN))
           allocate(tmp_AVGSOLVEL_P(DIMN,size(AVGSOLVEL_P)))
          allocate(tmp_EPG_P(size(EPG_P)))
           ENDIF

           IF(USE_COHESION) THEN
             allocate(tmp_PostCohesive(size(PostCohesive)))
           ENDIF

           IF (DO_OLD) THEN
             allocate(tmp_DES_POS_OLD(size(DES_POS_OLD),3))
            allocate(tmp_DES_VEL_OLD(size(DES_VEL_OLD),3))
           allocate(tmp_DES_ACC_OLD(3,size(DES_ACC_OLD)))
          allocate(tmp_OMEGA_OLD(size(OMEGA_OLD),3))
         allocate(tmp_ROT_ACC_OLD(DIMN,size(ROT_ACC_OLD)))
           ENDIF

           IF(ENERGY_EQ)THEN
             allocate(tmp_DES_T_s_OLD(size(DES_T_s_OLD)))
            allocate(tmp_DES_T_s_NEW(size(DES_T_s_NEW)))
           allocate(tmp_DES_C_PS(size(DES_C_PS)))
          allocate(tmp_DES_X_s(size(DES_X_s),DIMENSION_N_S))
         allocate(tmp_Q_Source(size(Q_Source)))

              IF (INTG_ADAMS_BASHFORTH) &
                  allocate(tmp_Q_Source0(size(Q_Source0)))
           ENDIF

           IF(ANY_SPECIES_EQ)THEN
             allocate(tmp_DES_R_sp(size(DES_R_sp),DIMENSION_N_S))
            allocate(tmp_DES_R_sc(size(DES_R_sc), DIMENSION_N_S))

              IF (INTG_ADAMS_BASHFORTH) THEN
                allocate(tmp_dMdt_OLD(size( dMdt_OLD)))
               allocate(tmp_dXdt_OLD(size( dXdt_OLD),DIMENSION_N_S))
              ENDIF

             allocate(tmp_Qint(size(Qint)))
           ENDIF

           IF(DES_USR_VAR_SIZE > 0) &
               allocate(tmp_DES_USR_VAR(DES_USR_VAR_SIZE,size(DES_USR_VAR)))






            tmp_des_radius = des_radius
            tmp_RO_Sol = RO_Sol
            tmp_PVOL = PVOL
            tmp_PMASS = PMASS
            tmp_OMOI = OMOI
            tmp_DES_POS_NEW = DES_POS_NEW
            tmp_DES_VEL_NEW = DES_VEL_NEW
            tmp_OMEGA_NEW = OMEGA_NEW
            tmp_PPOS = PPOS
            tmp_PARTICLE_STATE = PARTICLE_STATE
            tmp_orig_index = orig_index
            tmp_iglobal_id = iglobal_id
            tmp_pijk = pijk
            tmp_dg_pijk = dg_pijk
            tmp_dg_pijkprv = dg_pijkprv
            tmp_ighost_updated = ighost_updated
            tmp_FC = FC
            tmp_TOW = TOW
            tmp_F_GP = F_GP
            tmp_WALL_COLLISION_FACET_ID = WALL_COLLISION_FACET_ID
            tmp_WALL_COLLISION_PFT = WALL_COLLISION_PFT
            tmp_DRAG_FC = DRAG_FC

            tmp_NEIGHBOR_INDEX = NEIGHBOR_INDEX
            tmp_NEIGHBOR_INDEX_OLD = NEIGHBOR_INDEX_OLD

           IF(PARTICLE_ORIENTATION) tmp_ORIENTATION = ORIENTATION

           IF(FILTER_SIZE > 0) THEN
             tmp_FILTER_CELL = FILTER_CELL
            tmp_FILTER_WEIGHT = FILTER_WEIGHT
           ENDIF

           IF(MPPIC) THEN
             tmp_DES_STAT_WT = DES_STAT_WT
            tmp_PS_GRAD = PS_GRAD
           tmp_AVGSOLVEL_P = AVGSOLVEL_P
          tmp_EPG_P = EPG_P
           ENDIF

           IF(USE_COHESION) THEN
             tmp_PostCohesive = PostCohesive
           ENDIF

           IF (DO_OLD) THEN
             tmp_DES_POS_OLD = DES_POS_OLD
            tmp_DES_VEL_OLD = DES_VEL_OLD
           tmp_DES_ACC_OLD = DES_ACC_OLD
          tmp_OMEGA_OLD = OMEGA_OLD
         tmp_ROT_ACC_OLD = ROT_ACC_OLD
           ENDIF

           IF(ENERGY_EQ)THEN
             tmp_DES_T_s_OLD = DES_T_s_OLD
            tmp_DES_T_s_NEW = DES_T_s_NEW
           tmp_DES_C_PS = DES_C_PS
          tmp_DES_X_s = DES_X_s
         tmp_Q_Source = Q_Source

              IF (INTG_ADAMS_BASHFORTH) &
                  tmp_Q_Source0 = Q_Source0
           ENDIF

           IF(ANY_SPECIES_EQ)THEN
             tmp_DES_R_sp = DES_R_sp
            tmp_DES_R_sc = DES_R_sc

              IF (INTG_ADAMS_BASHFORTH) THEN
                tmp_dMdt_OLD =  dMdt_OLD
               tmp_dXdt_OLD =  dXdt_OLD
              ENDIF

             tmp_Qint = Qint
           ENDIF

           IF(DES_USR_VAR_SIZE > 0) &
               tmp_DES_USR_VAR = DES_USR_VAR











            tmp_des_radius = des_radius
            tmp_RO_Sol = RO_Sol
            tmp_PVOL = PVOL
            tmp_PMASS = PMASS
            tmp_OMOI = OMOI
            tmp_DES_POS_NEW = DES_POS_NEW
            tmp_DES_VEL_NEW = DES_VEL_NEW
            tmp_OMEGA_NEW = OMEGA_NEW
            tmp_PPOS = PPOS
            tmp_PARTICLE_STATE = PARTICLE_STATE
            tmp_orig_index = orig_index
            tmp_iglobal_id = iglobal_id
            tmp_pijk = pijk
            tmp_dg_pijk = dg_pijk
            tmp_dg_pijkprv = dg_pijkprv
            tmp_ighost_updated = ighost_updated
            tmp_FC = FC
            tmp_TOW = TOW
            tmp_F_GP = F_GP
            tmp_WALL_COLLISION_FACET_ID = WALL_COLLISION_FACET_ID
            tmp_WALL_COLLISION_PFT = WALL_COLLISION_PFT
            tmp_DRAG_FC = DRAG_FC

            tmp_NEIGHBOR_INDEX = NEIGHBOR_INDEX
            tmp_NEIGHBOR_INDEX_OLD = NEIGHBOR_INDEX_OLD

           IF(PARTICLE_ORIENTATION) tmp_ORIENTATION = ORIENTATION

           IF(FILTER_SIZE > 0) THEN
             tmp_FILTER_CELL = FILTER_CELL
            tmp_FILTER_WEIGHT = FILTER_WEIGHT
           ENDIF

           IF(MPPIC) THEN
             tmp_DES_STAT_WT = DES_STAT_WT
            tmp_PS_GRAD = PS_GRAD
           tmp_AVGSOLVEL_P = AVGSOLVEL_P
          tmp_EPG_P = EPG_P
           ENDIF

           IF(USE_COHESION) THEN
             tmp_PostCohesive = PostCohesive
           ENDIF

           IF (DO_OLD) THEN
             tmp_DES_POS_OLD = DES_POS_OLD
            tmp_DES_VEL_OLD = DES_VEL_OLD
           tmp_DES_ACC_OLD = DES_ACC_OLD
          tmp_OMEGA_OLD = OMEGA_OLD
         tmp_ROT_ACC_OLD = ROT_ACC_OLD
           ENDIF

           IF(ENERGY_EQ)THEN
             tmp_DES_T_s_OLD = DES_T_s_OLD
            tmp_DES_T_s_NEW = DES_T_s_NEW
           tmp_DES_C_PS = DES_C_PS
          tmp_DES_X_s = DES_X_s
         tmp_Q_Source = Q_Source

              IF (INTG_ADAMS_BASHFORTH) &
                  tmp_Q_Source0 = Q_Source0
           ENDIF

           IF(ANY_SPECIES_EQ)THEN
             tmp_DES_R_sp = DES_R_sp
            tmp_DES_R_sc = DES_R_sc

              IF (INTG_ADAMS_BASHFORTH) THEN
                tmp_dMdt_OLD =  dMdt_OLD
               tmp_dXdt_OLD =  dXdt_OLD
              ENDIF

             tmp_Qint = Qint
           ENDIF

           IF(DES_USR_VAR_SIZE > 0) &
               tmp_DES_USR_VAR = DES_USR_VAR


            allocate(tmp_int(max_pip))
            !tmp_int(:,:) = pijk(:,:)
            tmp_int(:) = particle_state(:)
            if (mype.eq.1) print *,"PARTICLE 1122 IS ",pijk(1122,1),orig_index(1122),particle_state(1122)

            CALL SORT_PARTICLES(1,size(PARTICLE_STATE),.false.)

            do ii=1,max_pip
               if (orig_index(ii).eq.1122) then
                  if (mype.eq.1) print *,"particle 1122 is       ",ii,pijk(ii,1),orig_index(ii),particle_state(ii)
               endif
            enddo

            CALL SORT_PARTICLES(1,size(PARTICLE_STATE),.true.)

            if (mype.eq.1) print *,"PARTICLE 1122 IS  ",pijk(1122,1),orig_index(1122),particle_state(1122)





            if (any(tmp_des_radius .ne. des_radius)) stop 12399
            if (any(tmp_RO_Sol .ne. RO_Sol)) stop 1234
            if (any(tmp_PVOL .ne. PVOL)) stop 1235
            if (any(tmp_PMASS .ne. PMASS)) stop 1236
            if (any(tmp_OMOI .ne. OMOI)) stop 1237
            if (any(tmp_DES_POS_NEW .ne. DES_POS_NEW)) stop 1238
            if (any(tmp_DES_VEL_NEW .ne. DES_VEL_NEW)) stop 123
            if (any(tmp_OMEGA_NEW .ne. OMEGA_NEW)) stop 123
            if (any(tmp_PPOS .ne. PPOS)) stop 123
            if (any(tmp_PARTICLE_STATE .ne. PARTICLE_STATE)) stop 123
            if (any(tmp_orig_index .ne. orig_index)) stop 123
            if (any(tmp_iglobal_id .ne. iglobal_id)) stop 123
            if (any(tmp_pijk .ne. pijk)) stop 123
            if (any(tmp_dg_pijk .ne. dg_pijk)) stop 123
            if (any(tmp_dg_pijkprv .ne. dg_pijkprv)) stop 123
            if (any(tmp_ighost_updated .neqv. ighost_updated)) stop 123
            if (any(tmp_FC .ne. FC)) stop 123
            if (any(tmp_TOW .ne. TOW)) stop 123
            if (any(tmp_F_GP .ne. F_GP)) stop 123
            if (any(tmp_WALL_COLLISION_FACET_ID .ne. WALL_COLLISION_FACET_ID)) stop 123
            if (any(tmp_WALL_COLLISION_PFT .ne. WALL_COLLISION_PFT)) stop 123
            if (any(tmp_DRAG_FC .ne. DRAG_FC)) stop 123

            if (any(tmp_NEIGHBOR_INDEX .ne. NEIGHBOR_INDEX)) stop 123
            if (any(tmp_NEIGHBOR_INDEX_OLD .ne. NEIGHBOR_INDEX_OLD)) stop 123

            IF(PARTICLE_ORIENTATION) tmp_ORIENTATION = ORIENTATION

            IF(FILTER_SIZE > 0) THEN
               if (any(tmp_FILTER_CELL .ne. FILTER_CELL)) stop 123
               if (any(tmp_FILTER_WEIGHT .ne. FILTER_WEIGHT)) stop 123
            ENDIF

            IF(MPPIC) THEN
               if (any(tmp_DES_STAT_WT .ne. DES_STAT_WT)) stop 123
               if (any(tmp_PS_GRAD .ne. PS_GRAD)) stop 123
               if (any(tmp_AVGSOLVEL_P .ne. AVGSOLVEL_P)) stop 123
               if (any(tmp_EPG_P .ne. EPG_P)) stop 123
            ENDIF

            IF(USE_COHESION) THEN
               if (any(tmp_PostCohesive .ne. PostCohesive)) stop 123
            ENDIF

            IF (DO_OLD) THEN
               if (any(tmp_DES_POS_OLD .ne. DES_POS_OLD)) stop 123
               if (any(tmp_DES_VEL_OLD .ne. DES_VEL_OLD)) stop 123
               if (any(tmp_DES_ACC_OLD .ne. DES_ACC_OLD)) stop 123
               if (any(tmp_OMEGA_OLD .ne. OMEGA_OLD)) stop 123
               if (any(tmp_ROT_ACC_OLD .ne. ROT_ACC_OLD)) stop 123
            ENDIF

            IF(ENERGY_EQ)THEN
               if (any(tmp_DES_T_s_OLD .ne. DES_T_s_OLD)) stop 123
               if (any(tmp_DES_T_s_NEW .ne. DES_T_s_NEW)) stop 123
               if (any(tmp_DES_C_PS .ne. DES_C_PS)) stop 123
               if (any(tmp_DES_X_s .ne. DES_X_s)) stop 123
               if (any(tmp_Q_Source .ne. Q_Source)) stop 123
               
               IF (INTG_ADAMS_BASHFORTH) then
                  if (any(tmp_Q_Source0 .ne. Q_Source0)) stop 123
               ENDIF
            ENDIF
               
            IF(ANY_SPECIES_EQ)THEN
               if (any(tmp_DES_R_sp .ne. DES_R_sp)) stop 123
               if (any(tmp_DES_R_sc .ne. DES_R_sc)) stop 123
               
               IF (INTG_ADAMS_BASHFORTH) THEN
                  if (any(tmp_dMdt_OLD .ne.  dMdt_OLD)) stop 123
                  if (any(tmp_dXdt_OLD .ne.  dXdt_OLD)) stop 123
               ENDIF
               
               if (any(tmp_Qint .ne. Qint)) stop 123
            ENDIF
            
            IF(DES_USR_VAR_SIZE > 0) then
               if (any(tmp_DES_USR_VAR .ne. DES_USR_VAR)) stop 123
            end IF



            deallocate(tmp_des_radius)
            deallocate(tmp_RO_Sol)
            deallocate(tmp_PVOL)
            deallocate(tmp_PMASS)
            deallocate(tmp_OMOI)
            deallocate(tmp_DES_POS_NEW)
            deallocate(tmp_DES_VEL_NEW)
            deallocate(tmp_OMEGA_NEW)
            deallocate(tmp_PPOS)
            deallocate(tmp_PARTICLE_STATE)
            deallocate(tmp_orig_index)
            deallocate(tmp_iglobal_id)
            deallocate(tmp_pijk)
            deallocate(tmp_dg_pijk)
            deallocate(tmp_dg_pijkprv)
            deallocate(tmp_ighost_updated)
            deallocate(tmp_FC)
            deallocate(tmp_TOW)
            deallocate(tmp_F_GP)
            deallocate(tmp_WALL_COLLISION_FACET_ID)
            deallocate(tmp_WALL_COLLISION_PFT)
            deallocate(tmp_DRAG_FC)
            deallocate(tmp_NEIGHBOR_INDEX)
            deallocate(tmp_NEIGHBOR_INDEX_OLD)

           IF(PARTICLE_ORIENTATION) deallocate(tmp_ORIENTATION)

           IF(FILTER_SIZE > 0) THEN
             deallocate(tmp_FILTER_CELL)
            deallocate(tmp_FILTER_WEIGHT)
           ENDIF

           IF(MPPIC) THEN
             deallocate(tmp_DES_STAT_WT)
            deallocate(tmp_PS_GRAD)
           deallocate(tmp_AVGSOLVEL_P)
          deallocate(tmp_EPG_P)
           ENDIF

           IF(USE_COHESION) THEN
             deallocate(tmp_PostCohesive)
           ENDIF

           IF (DO_OLD) THEN
             deallocate(tmp_DES_POS_OLD)
            deallocate(tmp_DES_VEL_OLD)
           deallocate(tmp_DES_ACC_OLD)
          deallocate(tmp_OMEGA_OLD)
         deallocate(tmp_ROT_ACC_OLD)
           ENDIF

           IF(ENERGY_EQ)THEN
             deallocate(tmp_DES_T_s_OLD)
            deallocate(tmp_DES_T_s_NEW)
           deallocate(tmp_DES_C_PS)
          deallocate(tmp_DES_X_s)
         deallocate(tmp_Q_Source)

              IF (INTG_ADAMS_BASHFORTH) &
                  deallocate(tmp_Q_Source0)
           ENDIF

           IF(ANY_SPECIES_EQ)THEN
             deallocate(tmp_DES_R_sp)
            deallocate(tmp_DES_R_sc)

              IF (INTG_ADAMS_BASHFORTH) THEN
                deallocate(tmp_dMdt_OLD)
               deallocate(tmp_dXdt_OLD)
              ENDIF

             deallocate(tmp_Qint)
           ENDIF

           IF(DES_USR_VAR_SIZE > 0) &
               deallocate(tmp_DES_USR_VAR)











            do ii=1,max_pip
                if (tmp_int(ii) .ne. particle_state(ii)) then
                   print *,"PIJK   ",mype,"||",ii,"||",tmp_int(ii),"||",particle_state(ii)
                   stop 99239
                endif

                ! if (tmp_int(ii,1) .ne. pijk(ii,1)) then
                !    print *,"PIJK   ",mype,"||",ii,"||",tmp_int(ii,1),"||",pijk(ii,1)
                !    stop 99239
                ! endif

                ! if (tmp_int(ii,2) .ne. pijk(ii,2)) then
                !    print *,"PIJK   ",mype,"||",ii,"||",tmp_int(ii,2),"||",pijk(ii,2)
                !    stop 99239
                ! endif

                ! if (tmp_int(ii,3) .ne. pijk(ii,3)) then
                !    print *,"PIJK   ",mype,"||",ii,"||",tmp_int(ii,3),"||",pijk(ii,3)
                !    stop 99239
                ! endif                

                ! if (tmp_int(ii,4) .ne. pijk(ii,4)) then
                !    print *,"PIJK   ",mype,"||",ii,"||",tmp_int(ii,4),"||",pijk(ii,4)
                !    stop 99239
                ! endif

                ! if (tmp_int(ii,5) .ne. pijk(ii,5)) then
                !    print *,"PIJK   ",mype,"||",ii,"||",tmp_int(ii,5),"||",pijk(ii,5)
                !    stop 99239
                ! endif
            enddo
            deallocate(tmp_int)

            CALL FIND_STATE_BOUNDS
            CALL NEIGHBOUR
         ENDIF

! Explicitly coupled simulations do not need to rebin particles to
! the fluid grid every time step. However, this implies that the
! fluid cell information and interpolation weights become stale.
         IF(DES_CONTINUUM_COUPLED .AND. &
            .NOT.DES_EXPLICITLY_COUPLED) THEN
! Bin particles to fluid grid.
            CALL PARTICLES_IN_CELL
! Calculate interpolation weights
            CALL CALC_INTERP_WEIGHTS
! Calculate mean fields (EPg).
            CALL COMP_MEAN_FIELDS
         ENDIF

! Update time to reflect changes
         S_TIME = S_TIME + DTSOLID

! The following section targets data writes for DEM only cases:
         IF(.NOT.DES_CONTINUUM_COUPLED) THEN
! Keep track of TIME and number of steps for DEM simulations
            TIME = S_TIME
            NSTEP = NSTEP + 1
! Call the output manager to write RES and SPx data.
            CALL OUTPUT_MANAGER(.FALSE., .FALSE.)
         ENDIF  ! end if (.not.des_continuum_coupled)

         IF(CALL_USR) CALL USR2_DES

      ENDDO ! end do NN = 1, FACTOR

! END DEM time loop
!-----------------------------------------------------------------<<<

      IF(CALL_USR) CALL USR3_DES

!      CALL DIFFUSE_MEAN_FIELDS
!      CALL CALC_EPG_DES

! When coupled, and if needed, reset the discrete time step accordingly
      IF(DT.LT.DTSOLID_TMP) THEN
         DTSOLID = DTSOLID_TMP
      ENDIF

      IF(TMP_DTS.NE.ZERO) THEN
         DTSOLID = TMP_DTS
         TMP_DTS = ZERO
      ENDIF

      IF(.NOT.DES_CONTINUUM_COUPLED)THEN
         WRITE(ERR_MSG,"('<---------- END DES_TIME_MARCH ----------')")
         CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
      ELSE
         call send_recv(ep_g,2)
         call send_recv(rop_g,2)
         call send_recv(des_u_s,2)
         call send_recv(des_v_s,2)
         if(do_K) call send_recv(des_w_s,2)
         call send_recv(rop_s,2)
         if(ENERGY_EQ) call send_recv(des_energy_source,2)

         TMP_WALL = WALL_TIME() - TMP_WALL
         IF(TMP_WALL > 1.0d-10) THEN
            WRITE(ERR_MSG, 9000) trim(iVal(dble(FACTOR)/TMP_WALL))
         ELSE
            WRITE(ERR_MSG, 9000) '+Inf'
         ENDIF
         CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE., LOG=.FALSE.)

 9000 FORMAT('    NITs/SEC = ',A)

      ENDIF

      RETURN
      END SUBROUTINE DES_TIME_MARCH
