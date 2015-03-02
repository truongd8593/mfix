!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  SUBROUTINE: TIME_MARCH                                              !
!  Author: M. Syamlal                                 Date: 21-JAN-92  !
!                                                                      !
!  Purpose: Controlling module for time marching and finding the       !
!           solution of equations from TIME to TSTOP at intervals of   !
!           DT, updating the b.c.'s, and creating output.              !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE TIME_MARCH

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE run
      USE physprop
      USE fldvar
      USE geometry
      USE pgcor
      USE pscor
      USE cont
      USE tau_g
      USE tau_s
      USE visc_g
      USE visc_s
      USE funits
      USE vshear
      USE scalars
      USE toleranc
      USE drag
      USE rxns, only: nRR
      USE compar
      USE time_cpu
      USE discretelement
      USE leqsol, only: SOLVER_STATISTICS, REPORT_SOLVER_STATS
      use mpi_utility
      USE cdist
      USE MFIX_netcdf
      USE cutcell
      USE vtk
      USE qmom_kinetic_equation
      USE dashboard
      USE indices
      USE bc
      USE coeff
      USE stiff_chem, only : STIFF_CHEMISTRY, STIFF_CHEM_SOLVER
      USE vtp

      use output, only: RES_DT, NLOG


      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Flag to indicate one pass through iterate for steady
! state conditions.
      LOGICAL :: FINISH

! Loop indices
      INTEGER :: L, M , I, IJK, N
! Error index
      INTEGER :: IER
! Number of iterations
      INTEGER :: NIT, NIT_TOTAL
! used for activating check_data_30
      INTEGER :: NCHECK, DNCHECK
! dummy logical variable for initializing adjust_dt
      LOGICAL :: dummy

! AEOLUS : stop trigger mechanism to terminate MFIX normally before
! batch queue terminates
      DOUBLE PRECISION :: CPU_STOP
      LOGICAL :: eofBATCHQ

!-----------------------------------------------
! External functions
!-----------------------------------------------
! use function MAX_VEL_INLET to compute max. velocity at inlet
      DOUBLE PRECISION :: MAX_VEL_INLET
! use function vavg_v_g to catch NaN's
      DOUBLE PRECISION :: VAVG_U_G, VAVG_V_G, VAVG_W_G, X_vavg

      LOGICAL , EXTERNAL :: ADJUST_DT
!-----------------------------------------------

      IF(AUTOMATIC_RESTART) RETURN


      FINISH  = .FALSE.
      NCHECK  = NSTEP
      DNCHECK = 1
      CPU_IO  = ZERO
      NIT_TOTAL = 0
      eofBATCHQ = .FALSE.


      CALL INIT_OUTPUT_VARS

! Parse residual strings
      CALL PARSE_RESID_STRING (IER)

! Call user-defined subroutine to set constants, check data, etc.
      IF (CALL_USR) CALL USR0

      CALL RRATES_INIT(IER)

! Calculate all the coefficients once before entering the time loop
      CALL INIT_COEFF(IER)

      DO M=1, MMAX
         CALL ZERO_ARRAY (F_gs(1,M), IER)
      ENDDO

! Remove undefined values at wall cells for scalars
      CALL UNDEF_2_0 (ROP_G, IER)
      DO M = 1, MMAX
         CALL UNDEF_2_0 (ROP_S(1,M), IER)
      ENDDO

! Initialize d's and e's to zero
      DO M = 0, MMAX
         CALL ZERO_ARRAY (D_E(1,M), IER)
         CALL ZERO_ARRAY (D_N(1,M), IER)
         CALL ZERO_ARRAY (D_T(1,M), IER)
      ENDDO
      CALL ZERO_ARRAY (E_E, IER)
      CALL ZERO_ARRAY (E_N, IER)
      CALL ZERO_ARRAY (E_T, IER)

! Initialize adjust_ur
      dummy = ADJUST_DT(100, 0)

! calculate shear velocities if periodic shear BCs are used
      IF(SHEAR) CALL CAL_D(V_sh)

! Initialize check_mass_balance.  This routine is not active by default.
! Specify a reporting interval (hard-wired in the routine) to activate
! the routine.
      Call check_mass_balance (0)

! sof modification: now it's only needed to do this once before time-loop
! Mark the phase whose continuity will be solved and used to correct
! void/volume fraction in calc_vol_fr (see subroutine for details)
      CALL MARK_PHASE_4_COR (PHASE_4_P_G, PHASE_4_P_S, DO_CONT, MCP,&
          DO_P_S, SWITCH_4_P_G, SWITCH_4_P_S, IER)

! uncoupled discrete element simulations do not need to be within
! the two fluid model time-loop
      IF(DISCRETE_ELEMENT.AND.(.NOT.DES_CONTINUUM_COUPLED))  THEN
         IF(WRITE_VTK_FILES) THEN
            DO L = 1, DIMENSION_VTK
               ! CALL WRITE_VTP_FILE(L)
            ENDDO
         ENDIF
         CALL DES_TIME_MARCH
         CALL CPU_TIME(CPU_STOP)
         CPU_STOP = CPU_STOP - CPU00
         IF(myPE.EQ.PE_IO) &
            write(*,"('Elapsed CPU time = ',E15.6,' sec')") CPU_STOP
         CALL PARALLEL_FIN
         STOP
      ENDIF


! The TIME loop begins here.............................................
 100  CONTINUE

! AEOLUS: stop trigger mechanism to terminate MFIX normally before batch
! queue terminates
      IF (CHK_BATCHQ_END) CALL CHECK_BATCH_QUEUE_END

      IF (CALL_USR) CALL USR1

! Remove solids from cells containing very small quantities of solids
      IF(.NOT.(DISCRETE_ELEMENT .OR. QMOMK) .OR. &
         DES_CONTINUUM_HYBRID) THEN
         IF(KT_TYPE_ENUM == GHD_2007) THEN
            CALL ADJUST_EPS_GHD
         ELSE
            CALL ADJUST_EPS
         ENDIF
      ENDIF

! sof modification: uncomment code below and modify MARK_PHASE_4_COR to
! use previous MFIX algorithm. Nov 22 2010.
! Mark the phase whose continuity will be solved and used to correct
! void/volume fraction in calc_vol_fr (see subroutine for details)
!      CALL MARK_PHASE_4_COR (PHASE_4_P_G, PHASE_4_P_S, DO_CONT, MCP,&
!           DO_P_S, SWITCH_4_P_G, SWITCH_4_P_S, IER)

! Set wall boundary conditions and transient flow b.c.'s
      CALL SET_BC1

      CALL OUTPUT_MANAGER(eofBATCHQ, FINISH)

      IF (DT == UNDEFINED) THEN
         IF (FINISH) THEN
            RETURN
         ELSE
            FINISH = .TRUE.
         ENDIF

! Mechanism to terminate MFIX normally before batch queue terminates.
      ELSEIF (TIME + 0.1d0*DT >= TSTOP .OR. eofBATCHQ) THEN
         IF(SOLVER_STATISTICS) &
            CALL REPORT_SOLVER_STATS(NIT_TOTAL, NSTEP)
         RETURN
      ENDIF

! Update previous-time-step values of field variables
      CALL UPDATE_OLD

! Calculate coefficients
      CALL CALC_COEFF_ALL (0, IER)

! Calculate the stress tensor trace and cross terms for all phases.
      CALL CALC_TRD_AND_TAU(IER)

! Calculate additional solid phase momentum source terms
! that arise from kinetic theory constitutive relations
      IF (.NOT.DISCRETE_ELEMENT .OR. DES_CONTINUUM_HYBRID) THEN
         CALL CALC_KTMOMSOURCE_U_S (IER)
         CALL CALC_KTMOMSOURCE_V_S (IER)
         CALL CALC_KTMOMSOURCE_W_S (IER)
      ENDIF

! Check rates and sums of mass fractions every NLOG time steps
      IF (NSTEP == NCHECK) THEN
         IF (DNCHECK < 256) DNCHECK = DNCHECK*2
         NCHECK = NCHECK + DNCHECK
! Upate the reaction rates for checking
         CALL CALC_RRATE(IER)
         CALL CHECK_DATA_30
      ENDIF

! Double the timestep for 2nd order accurate time implementation
      IF ((CN_ON.AND.NSTEP>1.AND.RUN_TYPE == 'NEW') .OR. &
          (CN_ON.AND.RUN_TYPE /= 'NEW' .AND. NSTEP >= (NSTEPRST+1))) THEN
         DT = 0.5d0*DT
         ODT = ODT * 2.0d0
      ENDIF

! Check for maximum velocity at inlet to avoid convergence problems
      MAX_INLET_VEL = 100.0d0*MAX_VEL_INLET()
! if no inlet velocity is specified, use an upper limit defined in
! toleranc_mod.f
      IF(MAX_INLET_VEL <= SMALL_NUMBER) THEN
         MAX_INLET_VEL = MAX_ALLOWED_VEL
         IF (UNITS == 'SI') MAX_INLET_VEL = 1D-2 * MAX_ALLOWED_VEL
      ENDIF
! Scale the value using a user defined scale factor
      MAX_INLET_VEL = MAX_INLET_VEL * MAX_INLET_VEL_FAC

! Advance the solution in time by iteratively solving the equations
 150  CALL ITERATE (IER, NIT)
      IF(AUTOMATIC_RESTART) RETURN

! Just to Check for NaN's, Uncomment the following lines and also lines
! of code in  VAVG_U_G, VAVG_V_G, VAVG_W_G to use.
!      X_vavg = VAVG_U_G ()
!      X_vavg = VAVG_V_G ()
!      X_vavg = VAVG_W_G ()
!      IF(AUTOMATIC_RESTART) RETURN

      DO WHILE (ADJUST_DT(IER,NIT))
         CALL ITERATE (IER, NIT)
      ENDDO


      IF(DT.LT.DT_MIN) THEN
         IF(TIME.LE.RES_DT .AND. AUTO_RESTART) THEN
            WRITE(ERR_MSG,1000)
            CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
            AUTO_RESTART = .FALSE.
         ENDIF

 1000 FORMAT('Automatic restart not possible as Total Time < RES_DT')

         IF(AUTO_RESTART) THEN
            AUTOMATIC_RESTART = .TRUE.
            RETURN
         ELSE
            IF(WRITE_DASHBOARD) THEN
               RUN_STATUS = 'DT < DT_MIN.  Recovery not possible!'
               CALL UPDATE_DASHBOARD(NIT,0.0d0,'    ')
            ENDIF
            CALL MFIX_EXIT(MyPE)
         ENDIF
      ENDIF

! Stiff Chemistry Solver.
      IF(STIFF_CHEMISTRY) THEN
         CALL STIFF_CHEM_SOLVER(DT, IER)
         IF(IER /= 0) THEN
            dummy = ADJUST_DT(IER, NIT)
            GOTO 150
         ENDIF
      ENDIF

! Check over mass and elemental balances.  This routine is not active by default.
! Edit the routine and specify a reporting interval to activate it.
      CALL CHECK_MASS_BALANCE (1)

! Othe solids model implementations
      IF (DEM_SOLIDS) CALL DES_TIME_MARCH
      IF (PIC_SOLIDS) CALL PIC_TIME_MARCH
      IF (QMOMK) CALL QMOMK_TIME_MARCH
      IF (CALL_DQMOM) CALL USR_DQMOM

! Advance the time step and continue
      IF((CN_ON.AND.NSTEP>1.AND.RUN_TYPE == 'NEW') .OR. &
         (CN_ON.AND.RUN_TYPE /= 'NEW' .AND. NSTEP >= (NSTEPRST+1))) THEN
! Double the timestep for 2nd order accurate time implementation
         DT = 2.d0*DT
         ODT = ODT * 0.5d0
! Perform the explicit extrapolation for CN implementation
         CALL CN_EXTRAPOL
      ENDIF

      IF (DT /= UNDEFINED) THEN
         IF(use_DT_prev) THEN
            TIME = TIME + DT_PREV
         ELSE
            TIME = TIME + DT
         ENDIF
         use_DT_prev = .FALSE.
         NSTEP = NSTEP + 1
      ENDIF

      NIT_TOTAL = NIT_TOTAL+NIT

! write (*,"('Compute the Courant number')")
! call get_stats(IER)

      FLUSH (6)
! The TIME loop ends here....................................................
      GOTO 100

      IF(SOLVER_STATISTICS) CALL REPORT_SOLVER_STATS(NIT_TOTAL, NSTEP)

      contains

!----------------------------------------------------------------------!
!                                                                      !
!  Subroutine: CHECK_BATCH_QUEUE_END                                   !
!  Author: A.Gel                                      Date:            !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE CHECK_BATCH_QUEUE_END

      use time_cpu, only: WALL_START

      use error_manager

! Logical flags for hault cases.
      LOGICAL :: USER_HAULT, WALL_HAULT
! Elapsed wall time, and fancy formatted buffer/batch queue times.
      DOUBLE PRECISION :: WALL_STOP, FANCY_BUFF, FANCY_BATCH
! Time units for formatted output.
      CHARACTER(LEN=4) :: WT_UNIT, BF_UNIT, BC_UNIT
! External function
      DOUBLE PRECISION :: WALL_TIME

! Calculate the current elapsed wall time.
      WALL_STOP = WALL_TIME()
      WALL_STOP = WALL_STOP - WALL_START

! Set flags for wall time exceeded or user specified hault.
      WALL_HAULT = ((WALL_STOP+TERM_BUFFER) >= BATCH_WALLCLOCK)
      INQUIRE(file="MFIX.STOP", exist=USER_HAULT)

! Report that the max user wall time was reached and exit.
      IF(WALL_HAULT) THEN
         CALL GET_TUNIT(WALL_STOP,WT_UNIT)
         FANCY_BUFF = TERM_BUFFER
         CALL GET_TUNIT(FANCY_BUFF, BF_UNIT)
         FANCY_BATCH = BATCH_WALLCLOCK
         CALL GET_TUNIT(FANCY_BATCH, BC_UNIT)
         WRITE(ERR_MSG, 1100) WALL_STOP, WT_UNIT, FANCY_BUFF, BF_UNIT, &
            FANCY_BATCH, BC_UNIT
         CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
      ENDIF

 1100 FORMAT(2/,15('='),' REQUESTED CPU TIME LIMIT REACHED ',('='),/   &
         'Batch Wall Time:',3X,F9.2,1X,A,/'Elapsed Wall Time: ',F9.2,  &
         1X,A,/'Term Buffer:',7X,F9.2,A,/15('='),' REQUESTED CPU ',    &
         'TIME LIMIT REACHED ',('='))

! Report that the hault signal was detected.
      IF(USER_HAULT) THEN
         WRITE(ERR_MSG, 1200)
         CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
      ENDIF

 1200 FORMAT(2/,19('='),' MFIX STOP SIGNAL DETECTED ',19('='),/'MFIX.',&
         'STOP file detected in run directory. Terminating MFIX.',/    &
         'Please DO NOT FORGET to erase the MFIX.STOP file before ',   &
         'restarting',/19('='),'MFIX STOP SIGNAL DETECTED ',19('='))

! This routine was restructured so all MPI ranks to the same action. As
! a result, broadcasting the BATCHQ flag may not be needed.
      eofBATCHQ = (WALL_HAULT .OR. USER_HAULT)
      call bcast (eofBATCHQ,PE_IO)

      END SUBROUTINE CHECK_BATCH_QUEUE_END

      END SUBROUTINE TIME_MARCH

