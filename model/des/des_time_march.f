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
      INTEGER :: FACTOR,ii,jj
      INTEGER, dimension(:), allocatable :: tmp_int

! Temporary variables when des_continuum_coupled is T to track
! changes in solid time step
      DOUBLE PRECISION :: TMP_DTS, DTSOLID_TMP

! Numbers to calculate wall time spent in DEM calculations.
      DOUBLE PRECISION :: TMP_WALL

      LOGICAL, SAVE :: FIRST_PASS = .TRUE.

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

      IF(FIRST_PASS) THEN
         CALL SORT_PARTICLES(1,size(PARTICLE_STATE),.false.)
         FIRST_PASS = .false.
      ENDIF

         if (.not.sorted) then
            print *,"NOT SORTED?!?!"
            stop 100100100
         endif

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

         if (.not.sorted) then
            print *,"NOT SORTED?!?!"
            stop 323
         endif

if (.false. .and. mype.eq.1) then
         print *,mype,":  BEFORE_CALC_FD neighbors of 208(",iglobal_id(208),") are:"
do ii=NEIGHBOR_INDEX(1,208), NEIGHBOR_INDEX(2,208)-1
         print *,mype,":  BEFORE_CALC_FD ",neighbors(ii),"(",iglobal_id(neighbors(ii)),")"
enddo
endif

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

         if (.not.sorted) then
            print *,"NOT SORTED?!?!"
            stop 345543111
         endif

! Call exchange particles - this will exchange particle crossing
! boundaries as well as updates ghost particles information
         IF (DO_NSEARCH .OR. (numPEs>1) .OR. DES_PERIODIC_WALLS) THEN
            CALL DESGRID_PIC(.TRUE.)

         if (.not.sorted) then
            print *,"NOT SORTED?!?!"
            stop 77661
         endif


if (.false. .and. mype.eq.1) then
         print *,mype,":  BEFORE_PAR_EX neighbors of 208(",iglobal_id(208),") are:"
do ii=NEIGHBOR_INDEX(1,208), NEIGHBOR_INDEX(2,208)-1
         print *,mype,":  BEFORE_PAR_EX ",neighbors(ii),"(",iglobal_id(neighbors(ii)),")"
enddo
endif



            CALL DES_PAR_EXCHANGE

if (.false. .and. mype.eq.1) then
         print *,mype,":  AFTER_PAR_EX neighbors of 208(",iglobal_id(208),") are:"
do ii=NEIGHBOR_INDEX(1,208), NEIGHBOR_INDEX(2,208)-1
         print *,mype,":  AFTER_PAR_EX ",neighbors(ii),"(",iglobal_id(neighbors(ii)),")"
enddo
endif

            do ii=1,max_pip
               orig_index(ii) = ii
            enddo

!            sorted = .true.
            CALL SORT_PARTICLES(1,size(PARTICLE_STATE),.false.)
            !CALL SORT_PARTICLES(1,size(PARTICLE_STATE),.true.)


if (.false. .and. mype.eq.1) then
         print *,mype,":  AFTER_SORT neighbors of 208(",iglobal_id(208),") are:"
do ii=NEIGHBOR_INDEX(1,208), NEIGHBOR_INDEX(2,208)-1
         print *,mype,":  AFTER_SORT ",neighbors(ii),"(",iglobal_id(neighbors(ii)),")"
      enddo
   endif

         ENDIF

         IF(DO_NSEARCH) THEN

            CALL FIND_STATE_BOUNDS

         if (.not.sorted) then
            print *,"NOT SORTED?!?!"
            stop 21212
         endif

            CALL NEIGHBOUR

         if (.not.sorted) then
            print *,"NOT SORTED?!?!"
            stop 1111
         endif

         ENDIF

         if (.not.sorted) then
            print *,"NOT SORTED?!?!"
            stop 88888
         endif

if (.false.) then

do jj=1, max_pip
if (1317.eq.iglobal_id(jj)) then
!if (.false. .and. mype.eq.1) then
         print *,mype,":  AFTER_NEIGH neighbors of 208(",iglobal_id(jj),") are:"
do ii=NEIGHBOR_INDEX(1,jj), NEIGHBOR_INDEX(2,jj)-1
         print *,mype,":  AFTER_NEIGH ",neighbors(ii),"(",iglobal_id(neighbors(ii)),")"
      enddo
   endif
enddo
endif

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

         if (.not.sorted) then
            print *,"NOT SORTED?!?!"
            stop 555
         endif

if (.false. .and. mype.eq.1) then
         print *,mype,":  AFTER_COMP neighbors of 208(",iglobal_id(208),") are:"
do ii=NEIGHBOR_INDEX(1,208), NEIGHBOR_INDEX(2,208)-1
         print *,mype,":  AFTER_COMP ",neighbors(ii),"(",iglobal_id(neighbors(ii)),")"
      enddo
   endif

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

         if (.not.sorted) then
            print *,"NOT SORTED?!?!"
            stop 2828
         endif

if (.false. .and. mype.eq.1) then
         print *,mype,":  END_LOOP neighbors of 208(",iglobal_id(208),") are:"
do ii=NEIGHBOR_INDEX(1,208), NEIGHBOR_INDEX(2,208)-1
         print *,mype,":  END_LOOP ",neighbors(ii),"(",iglobal_id(neighbors(ii)),")"
      enddo
   endif

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
