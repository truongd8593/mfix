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
      use fldvar, only: EP_g, ROP_g, ROP_s
      use mfix_pic, only: MPPIC
      use output, only: SPX_DT
      use run, only: ANY_SPECIES_EQ
      use run, only: CALL_USR
      use run, only: ENERGY_EQ
      use run, only: RUN_TYPE
      use run, only: ANY_SPECIES_EQ

      use output, only: SPX_DT
      USE fldvar, only: EP_g, ROP_g, ROP_s
      use des_thermo, only: DES_ENERGY_SOURCE !Added by Surya Oct 22, 2014    
      use run, only: TIME, TSTOP, DT

      use desmpi, only: DES_PAR_EXCHANGE

      use discretelement
      use error_manager
      use functions
      use mpi_utility
      use sendrecv
      use stl

      IMPLICIT NONE
!------------------------------------------------
! Local variables
!------------------------------------------------
! Total number of particles
      INTEGER, SAVE :: NP=0

! time step loop counter index
      INTEGER :: NN
! loop counter index for any initial particle settling incoupled cases
      INTEGER :: FACTOR

      INTEGER :: II, JJ, KK, IJK, CELL_ID, I_CELL, J_CELL, K_CELL, COUNT, NF
      INTEGER :: IMINUS1, IPLUS1, JMINUS1, JPLUS1, KMINUS1, KPLUS1, PHASELL, LOC_MIN_PIP

! Local variables to keep track of time when dem restart and des
! write data need to be written when des_continuum_coupled is F
      DOUBLE PRECISION :: DES_RES_TIME, DES_SPX_TIME, NTIME

! Temporary variables when des_continuum_coupled is T to track
! changes in solid time step
      DOUBLE PRECISION :: TMP_DTS, DTSOLID_TMP

! Numbers to calculate wall time spent in DEM calculations.
      DOUBLE PRECISION :: WALL_TIME, TMP_WALL


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

! Set the DES_SPX and RES variables.
         IF(RUN_TYPE .EQ. 'NEW') THEN
            DES_SPX_TIME = S_TIME
            DES_RES_TIME = S_TIME
         ELSE
            NTIME = (S_TIME+0.1d0*DTSOLID)
            DES_SPX_TIME = (INT(NTIME/DES_SPX_DT) + 1)*DES_SPX_DT
            DES_RES_TIME = (INT(NTIME/DES_RES_DT) + 1)*DES_RES_DT
         ENDIF
      ENDIF   ! end if/else (des_continuum_coupled)

      NP = PIP - IGHOST_CNT
      CALL GLOBAL_ALL_SUM(NP)

      WRITE(ERR_MSG, 1000) trim(iVal(factor)), trim(iVAL(NP))
      CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)

 1000 FORMAT(/'DEM NITs: ',A,3x,'Total PIP: ', A)


      IF(DES_CONTINUUM_COUPLED) THEN
         DES_SPX_DT = SPX_DT(1)
         CALL CALC_PG_GRAD
         IF(.NOT.EXPLICITLY_COUPLED) THEN
            IF(ANY_SPECIES_EQ) CALL ZERO_RRATE_DES
            IF(ENERGY_EQ) CALL ZERO_ENERGY_SOURCE
         ENDIF
      ENDIF


      IF(CALL_USR) CALL USR0_DES

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
! Call exchange particles - this will exchange particle crossing
! boundaries as well as updates ghost particles information

         IF(DO_NSEARCH) THEN
            CALL DES_PAR_EXCHANGE
            CALL NEIGHBOUR
         ELSEIF (1 < numPEs) THEN
            CALL DES_PAR_EXCHANGE
         ENDIF

         IF(.NOT.EXPLICITLY_COUPLED) THEN
! Bin particles to fluid grid.
            CALL PARTICLES_IN_CELL
! Calculate interpolation weights
            CALL CALC_INTERP_WEIGHTS
! Calculate mean fields (EPg).
            CALL COMP_MEAN_FIELDS
         ENDIF

! Update time to reflect changes
         S_TIME = S_TIME + DTSOLID

! When coupled, all write calls are made in time_march (the continuum
! portion) according to user settings for spx_time and res_time.
! The following section targets data writes for DEM only cases:
         IF(.NOT.DES_CONTINUUM_COUPLED) THEN
! Keep track of TIME for DEM simulations
            TIME = S_TIME
! Calculate the 'next time' output is written.
            NTIME = S_TIME + 0.1d0*DTSOLID

! Write data using des_spx_time and des_res_time; note the time will
! reflect current position of particles
            IF(PRINT_DES_DATA) THEN
               IF((NTIME >= DES_SPX_TIME) .OR. (NTIME >= TSTOP)&
                  .OR. (NN == FACTOR) ) THEN
                  DES_SPX_TIME =(INT(NTIME/DES_SPX_DT)+1)*DES_SPX_DT

! Granular temperature subroutine should be called/calculated when
! writing DES data
                  CALL DES_GRANULAR_TEMPERATURE
                  IF (DES_CALC_BEDHEIGHT) CALL CALC_DES_BEDHEIGHT
                  CALL WRITE_DES_DATA
               ENDIF
            ENDIF

! Write out the restart infomration here. Note that there is a call
! to write out the continuum variables.
            IF((NTIME >= DES_RES_TIME) .OR. (NTIME >= TSTOP) .OR.      &
               (NN == FACTOR)) THEN
               DES_RES_TIME = (INT(NTIME/DES_RES_DT)+1)*DES_RES_DT
               CALL WRITE_RES0_DES
               CALL WRITE_RES1
            ENDIF
         ENDIF  ! end if (.not.des_continuum_coupled)

! Seed new particles entering the system.
         IF(DEM_BCMI > 0) CALL MASS_INFLOW_DEM
         IF(DEM_BCMO > 0) CALL MASS_OUTFLOW_DEM

         IF(CALL_USR) CALL USR2_DES

      ENDDO     ! end do NN = 1, FACTOR
! END DEM time loop
!-----------------------------------------------------------------<<<

      IF(CALL_USR) CALL USR3_DES


!      CALL DIFFUSE_MEAN_FIELDS
!      CALL CALC_EPG_DES

! When coupled the granular temperature subroutine is only calculated at end
! of the current DEM simulation
      IF(DES_CONTINUUM_COUPLED) THEN
         CALL DES_GRANULAR_TEMPERATURE()
         IF (DES_CALC_BEDHEIGHT) CALL CALC_DES_BEDHEIGHT()
! the call to identify clusters is now done in time_march, uncomment
! line below to compute clusters each fluid time step.
      ENDIF

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
         if(ENERGY_EQ) call send_recv(des_energy_source,2) !Added by Surya Oct 22, 2014

         TMP_WALL = WALL_TIME() - TMP_WALL
         IF(TMP_WALL > 1.0d-10) THEN
            WRITE(ERR_MSG, 9000) trim(iVal(dble(FACTOR)/TMP_WALL))
         ELSE
            WRITE(ERR_MSG, 9000) '+Inf'
         ENDIF
         CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)

 9000 FORMAT('    NITs/SEC = ',A)

      ENDIF

      RETURN
      END SUBROUTINE DES_TIME_MARCH
