!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Subroutine: PIC_TIME_MARCH                                          !
!  Author: R. Garg                                                     !
!                                                                      !
!  Purpose: Main PIC driver routine.                                   !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE PIC_TIME_MARCH

! Global variables
!---------------------------------------------------------------------//
! Fluid time, simulation end time, time step size, number of time steps
      use run, only: TIME, TSTOP, DT, NSTEP
! Discrete particle time, time step size
      use discretelement, only: S_TIME, DTSOLID
! MPPIC model step-size bounds
      use mfix_pic, only: DTPIC_MAX, DTPIC_CFL, DTPIC_TAUP
! Local particle count
      use discretelement, only: PIP
! Flag: Coupled fluid-solids simulation
      use discretelement, only: DES_CONTINUUM_COUPLED
! Flag: Store _OLD arrays
      use discretelement, only: DO_OLD


! Module procedures
!---------------------------------------------------------------------//
      use mpi_utility, only: GLOBAL_ALL_SUM
      use mpi_funs_des, only: DES_PAR_EXCHANGE

      use error_manager

      IMPLICIT NONE

! Local variables
!---------------------------------------------------------------------//
! time till which the PIC loop will be run
      double precision :: TEND_PIC_LOOP
! number of PIC time steps
      Integer :: PIC_ITERS
! Global number of parcels.
      INTEGER :: gPIP
!......................................................................!


! Set solids time to fluid time.
      S_TIME = TIME

      IF(DES_CONTINUUM_COUPLED) THEN
         TEND_PIC_LOOP = TIME+DT
      ELSE
         TEND_PIC_LOOP = TSTOP
         DTSOLID = DT
      ENDIF
      PIC_ITERS = 0

! Compute the gas-phase pressure gradient
      IF(DES_CONTINUUM_COUPLED) CALL CALC_PG_GRAD

! If the current time in the discrete loop exceeds the current time in
! the continuum simulation, exit the lagrangian loop
      DO WHILE(S_TIME.LT.TEND_PIC_LOOP)

         PIC_ITERS  = PIC_ITERS + 1

! Set the solids time step
!         DTSOLID = MERGE(MIN(DTPIC_MAX, DT), DTPIC_MAX,                &
!            DES_CONTINUUM_COUPLED)

! If next time step in the discrete loop will exceed the current time
! in the continuum simulation, modify the discrete time step so final
! time will match
         IF(S_TIME + DTSOLID > TEND_PIC_LOOP) &
            DTSOLID = TEND_PIC_LOOP - S_TIME

! Exchange particle crossing processor boundaries
         CALL DES_PAR_EXCHANGE
! Bin particles to the fluid grid
         CALL PARTICLES_IN_CELL
! Calculate mean fields
         CALL COMP_MEAN_FIELDS

! This was moved from particles in cell and the passed variables should
! be added to particles in cell or made global.
         CALL REPORT_STATS_PIC

! Calculate the solids pressure
         CALL CALC_PS_PIC
         CALL CALC_PS_GRAD_PIC
         CALL INTERPOLATE_PIC

         IF(DES_CONTINUUM_COUPLED) CALL CALC_DRAG_DES

         IF (DO_OLD) CALL CFUPDATEOLD

         CALL INTEGRATE_TIME_PIC 

! Impose the wall-particle boundary condition for pic case
! The same routine also applies the mass inflow/outflow BCs as well
         CALL PIC_APPLY_WALLBC_STL

! Update time to reflect changes
         S_TIME = S_TIME + DTSOLID

         DTPIC_MAX = MIN(DTPIC_CFL, DTPIC_TAUP)

! When coupled, all write calls are made in time_march (the continuum
! portion) according to user settings for spx_time and res_time.
! The following section targets data writes for DEM only cases:
         IF(.NOT.DES_CONTINUUM_COUPLED) THEN
! Keep track of TIME for DEM simulations
            TIME = S_TIME
            NSTEP = NSTEP + 1
! Call the output manager to write RES and SPx data.
            CALL OUTPUT_MANAGER(.FALSE., .FALSE.)
         ENDIF  ! end if (.not.des_continuum_coupled)

      ENDDO

      CALL GLOBAL_ALL_SUM(PIP, gPIP)
      WRITE(ERR_MSG, 3000) trim(iVal(PIC_ITERS)), trim(iVal(gPIP))
      CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)

 3000 FORMAT(/'PIC NITs: ',A,3x,'Total PIP: ', A)

      RETURN
      END SUBROUTINE PIC_TIME_MARCH




!         !DTPIC_MAX = MIN( 1e-04, DTPIC_MAX)
!         IF(MOD(PIC_ITERS, 10).eq.0) then
!            IF(DES_CONTINUUM_COUPLED) then
!               WRITE(ERR_MSG, 2000) DTSOLID, DTPIC_CFL, DTPIC_TAUP, DT
!            ELSE
!               WRITE(ERR_MSG, 2001) S_TIME, DTSOLID, DTPIC_CFL, DTPIC_TAUP, DT
!            ENDIF
!            CALL FLUSH_ERR_MSG(HEADER = .FALSE., FOOTER = .FALSE.)
!         ENDIF
!
! 2000 FORMAT(/5x,'DTSOLID CURRENT  = ',g17.8,/5x,'DTPIC_CFL',8x,'= ',  &
!         g17.8, /5x,'DTPIC TAUP',7x,'= ',g17.8,/5x,'DT FLOW',10x,'= ', &
!         g17.8)
!
! 2001 FORMAT(/5x,'TIME',13X,'= ',g17.8,/5x,'DTSOLID CURRENT  = ',g17.8,&
!         /5x,'DTPIC_CFL',8X,'= ', g17.8,/5x,'DTPIC TAUP',7x,'= ',g17.8,&
!         /5x,'DT FLOW',10X,'= ', g17.8)
