!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  SUBROUTINE: TIME_MARCH                                              C
!  Purpose: Controlling module for time marching and finding the       C
!           solution of equations from TIME to TSTOP at intervals of   C
!           DT, updating the b.c.'s, and creating output.              C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-JAN-92  C
!  Reviewer:M. Syamlal, S. Venkatesan, P. Nicoletti,  Date: 29-JAN-92  C
!           W. Rogers                                                  C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: Change subroutine name from SET_GASFLUX to SET_FLUXES      C
!           Add a CALC_THETA call for granular stresses                C
!  Author: M. Syamlal                                 Date: 20-FEB-92  C
!  Reviewer: S. Venkatesan                            Date: 11-DEC-92  C
!  Revision Number: 2                                                  C
!  Purpose: Changes for MFIX 2.0                                       C
!  Author: M. Syamlal                                 Date: 12-APR-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number: 2                                                  C
!  Purpose: To call DES related routines when doing DES                C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!                                                                      C
!  Revision Number: 3                                                  C
!  Purpose: To call ISAT to calculate chemical rxns                    C
!  Author: Nan Xie                                    Date: 02-Aug-04  C
!                                                                      C
!  Revision Number: 4                                                  C
!  Purpose: To call Cartesian grid subroutines, and update dasboard    C
!  Author: Jeff Dietiker                              Date: 01-Jul-09  C
!                                                                      C
!  Revision Number: 5                                                  C
!  Purpose: Incorporation of QMOM for the solution of the particle     C
!           kinetic equation                                           C
!  Author: Alberto Passalacqua - Fox Research Group   Date: 02-Dec-09  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE TIME_MARCH

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE run
      USE output
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
      USE rxns
      USE compar
      USE time_cpu
      USE discretelement
      USE leqsol
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

      IMPLICIT NONE
!-----------------------------------------------
! Local parameters
!-----------------------------------------------
      DOUBLE PRECISION, PARAMETER :: ONEMEG = 1048576
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Flag to indicate one pass through iterate for steady
! state conditions.
      LOGICAL :: FINISH
! Time at which standard output is to be written
      DOUBLE PRECISION :: OUT_TIME
! Time at which restart file is to be written
      DOUBLE PRECISION :: RES_TIME
! Time at which REAL restart file is to be written
      DOUBLE PRECISION :: SPX_TIME(N_SPX)
! Disk space needed for one variable and each SPX file
      DOUBLE PRECISION :: DISK_ONE, DISK(N_SPX)
! Total Disk space
      DOUBLE PRECISION :: DISK_TOT
! number SPX writes
      INTEGER :: ISPX

      LOGICAL :: RES_MSG, SPX_MSG
! Time at which special output is to be written
      DOUBLE PRECISION :: USR_TIME (DIMENSION_USR)
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

      CHARACTER(LEN=35) ::  EXT_END
! AEOLUS : stop trigger mechanism to terminate MFIX normally before
! batch queue terminates
      DOUBLE PRECISION :: CPU_STOP
      LOGICAL :: eofBATCHQ
! not used remove after verification
      LOGICAL :: bWrite_netCDF_files

      DOUBLE PRECISION :: TIME_START
      DOUBLE PRECISION :: WALL_START ! wall time at the beginning
      DOUBLE PRECISION :: WALL_NOW   ! wall time at the end of each timestep
      DOUBLE PRECISION :: WALL_LEFT, WALL_TOTAL
      CHARACTER(LEN=4) TUNIT, TOT_UNIT

!-----------------------------------------------
! External functions
!-----------------------------------------------
! use function MAX_VEL_INLET to compute max. velocity at inlet
      DOUBLE PRECISION :: MAX_VEL_INLET
! use function vavg_v_g to catch NaN's
      DOUBLE PRECISION :: VAVG_U_G, VAVG_V_G, VAVG_W_G, X_vavg

      LOGICAL , EXTERNAL :: ADJUST_DT
      DOUBLE PRECISION :: WALL_TIME
!-----------------------------------------------

      IF(AUTOMATIC_RESTART) RETURN

      EXT_END = '123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'

      FINISH  = .FALSE.
      NCHECK  = NSTEP
      DNCHECK = 1
      CPU_IO  = ZERO
      NIT_TOTAL = 0
      eofBATCHQ = .FALSE.


! Initialize times for writing outputs
      OUT_TIME = TIME
      LOG_HEADER = .TRUE.
      RES_MSG = .TRUE.
      SPX_MSG = .TRUE.

! Initialize disk space calculations
      DISK_TOT = ZERO
      DISK_ONE = 4.*IJKMAX2/ONEMEG
      DISK(1) = DISK_ONE
      DISK(2) = 2.*DISK_ONE
      DISK(3) = 3.*DISK_ONE
      DISK(4) = MMAX*3.*DISK_ONE
      DISK(5) = MMAX*DISK_ONE
      DISK(6) = 3.*DISK_ONE
      DISK(7) = NMAX(0)*DISK_ONE
      M = 1
      IF (MMAX > 0) THEN
         DISK(7) = DISK(7) + SUM(NMAX(1:MMAX)*DISK_ONE)
         M = MMAX + 1
      ENDIF
      DISK(8) = MMAX*DISK_ONE
      DISK(9) = NScalar*DISK_ONE
      DISK(10) = nRR*DISK_ONE
      IF(K_Epsilon) THEN
         DISK(11) = 2.*DISK_ONE
      ELSE
         DISK(11) = ZERO
      ENDIF


      IF (RUN_TYPE == 'NEW') THEN
         RES_TIME = TIME
         SPX_TIME(:N_SPX) = TIME
         L = N_SPX + 1
      ELSE
         IF (DT /= UNDEFINED) THEN
            RES_TIME = (INT((TIME + 0.1d0*DT)/RES_DT) + 1)*RES_DT
            SPX_TIME(:N_SPX) = (INT((TIME+0.1d0*DT)/SPX_DT(:N_SPX))+1)*&
               SPX_DT(:N_SPX)
            L = N_SPX + 1
         ENDIF
      ENDIF

      DO L = 1, DIMENSION_USR
         USR_TIME(L) = UNDEFINED
         IF (USR_DT(L) /= UNDEFINED) THEN
            IF (RUN_TYPE == 'NEW') THEN
               USR_TIME(L) = TIME
            ELSE
               USR_TIME(L) = (INT((TIME+0.1d0*DT)/USR_DT(L))+1)*USR_DT(L)
            ENDIF
         ENDIF
      ENDDO

! JFD modification: cartesian grid implementation
! Initialize VTK_TIME
      IF(WRITE_VTK_FILES) THEN
         DO L = 1, DIMENSION_VTK
            VTK_TIME(L) = UNDEFINED
            IF (VTK_DT(L) /= UNDEFINED) THEN
               IF (RUN_TYPE == 'NEW'.OR.RUN_TYPE=='RESTART_2') THEN
                  VTK_TIME(L) = TIME
               ELSE
                  VTK_TIME(L) = (INT((TIME + 0.1d0*DT)/VTK_DT(L))+1)*VTK_DT(L)
               ENDIF
            ENDIF
         ENDDO
      ENDIF

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

      WALL_START = WALL_TIME()
      TIME_START = TIME

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

      CALL CPU_TIME(CPU0_IO)

      CALL OUTPUT_MANAGEMENT


      IF (DT == UNDEFINED) THEN
         IF (FINISH) THEN
            RETURN
         ELSE
            FINISH = .TRUE.
         ENDIF

! AEOLUS: stop trigger mechanism to terminate MFIX normally before batch
! queue terminates
      ELSEIF (TIME + 0.1d0*DT >= TSTOP.OR.eofBATCHQ) THEN
         IF(solver_statistics) then
            WRITE(*,*) 'Total number of non-linear iterations', &
               NIT_TOTAL
            WRITE(*,*) 'Average number per time-step', NIT_TOTAL/NSTEP
            WRITE(*,*) 'Equation number', '-----', &
               'Number of linear solves'
            DO I = 1, 10
               WRITE(*,*) I, '---------',  iter_tot(I)
            ENDDO
            WRITE(*,*) 'Equation number', '-----', &
               'Avg. number of linear solves for NIT'
            DO I = 1, 10
               WRITE(*,*) I, '---------',  iter_tot(I)/NIT_TOTAL
            ENDDO
         ENDIF
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

      WALL_NOW = WALL_TIME()
      WALL_LEFT = (WALL_NOW-WALL_START)*(TSTOP-TIME)/                  &
         max(TIME-TIME_START,1.0d-6)
      CALL GET_TUNIT(WALL_LEFT,TUNIT)
      WALL_TOTAL = (WALL_NOW-WALL_START)*(TSTOP-TIME_START)/           &
         max(TIME-TIME_START,1.0d-6)
      CALL GET_TUNIT(WALL_TOTAL,TOT_UNIT)
      IF(DMP_LOG) WRITE (*,2000) 'Remaining', WALL_LEFT, trim(TUNIT),  &
         'Total',WALL_TOTAL, trim(TOT_UNIT)

 2000 FORMAT('Estimated Wall Time: ',2(3x,A,F9.3,1X,A))

      IF(DT.LT.DT_MIN) THEN
         IF(TIME.LE.RES_DT .AND. AUTO_RESTART) THEN
            IF (DMP_LOG)WRITE(UNIT_LOG,*) &
                 'Automatic restart not possible as Total Time < RES_DT'
            AUTO_RESTART = .FALSE.
         ENDIF

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
      Call check_mass_balance (1)

! DES
      IF (DEM_SOLIDS) CALL DES_TIME_MARCH
      IF (PIC_SOLIDS) CALL PIC_TIME_MARCH


! Alberto Passalacqua: QMOMK
      IF (QMOMK) CALL QMOMK_TIME_MARCH

      IF (CALL_DQMOM) CALL USR_DQMOM

! Advance the time step and continue
      IF ((CN_ON.AND.NSTEP>1.AND.RUN_TYPE == 'NEW') .OR. &
          (CN_ON.AND.RUN_TYPE /= 'NEW' .AND. NSTEP >= (NSTEPRST+1))) THEN
! Double the timestep for 2nd order accurate time implementation
         DT = 2.d0*DT
         ODT = ODT * 0.5d0
! Perform the explicit extrapolation for CN implementation
         CALL CN_EXTRAPOL
      ENDIF

      IF (DT /= UNDEFINED) THEN
         TIME = TIME + DT
         NSTEP = NSTEP + 1
      ENDIF

      NIT_TOTAL = NIT_TOTAL+NIT

! AIKEDEBUG 081101
! write (*,"('Compute the Courant number')")
! call get_stats(IER)

      FLUSH (6)
      GOTO 100

      IF(solver_statistics) then
         WRITE(*,*) 'Total number of non-linear iterations', NIT_TOTAL
         WRITE(*,*) 'Average number per time-step', NIT_TOTAL/NSTEP
         WRITE(*,*) 'Equation number', '-----', &
            'Number of linear solves'
         DO I = 1, 10
            Write(*,*) I, '---------',  iter_tot(I)
         ENDDO
         WRITE(*,*) 'Equation number', '-----',&
            'Avg. number of linear solves for NIT'
         DO I = 1, 10
            Write(*,*) I, '---------',  iter_tot(I)/NIT_TOTAL
         ENDDO
      ENDIF

! The TIME loop ends here....................................................





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

      use error_manager

! Logical flags for hault cases.
      LOGICAL :: USER_HAULT, WALL_HAULT
! Elapsed wall time, and fancy formatted buffer/batch queue times.
      DOUBLE PRECISION :: WALL_STOP, FANCY_BUFF, FANCY_BATCH
! Time units for formatted output.
      CHARACTER(LEN=4) :: WT_UNIT, BF_UNIT, BC_UNIT

! Calculate the current elapsed wall time.
      WALL_STOP = WALL_TIME()
      WALL_STOP = WALL_STOP - WALL_START

! Set flags for wall time exceeded or user specified hault.
      WALL_HAULT = ((WALL_STOP+TERM_BUFFER) >= BATCH_WALLCLOCK)
      INQUIRE(file="MFIX.STOP", exist=USER_HAULT)

! Fancy format the elapsed wall time and get its unit.
      CALL GET_TUNIT(WALL_STOP,WT_UNIT)

! Write out the amount of elapsed wall time.
      WRITE(ERR_MSG,1000) WALL_STOP, WT_UNIT
      CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)

 1000 FORMAT('Elapsed Wall time: ',F9.2,1X,A)

! Report that the max user wall time was reached and exit.
      IF(WALL_HAULT) THEN
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


!----------------------------------------------------------------------!
!                                                                      !
!  Subroutine: OUTPUT_MANAGEMENT                                       !
!  Author: J.Musser                                   Date:            !
!                                                                      !
!  Purpose: Relocate calls to write output files (RES, SPx, VTP). This !
!  was done to simplify the time_march code.                           !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE OUTPUT_MANAGEMENT

! Write standard output, if needed
      IF (OUT_DT /= UNDEFINED) THEN
         IF (DT == UNDEFINED) THEN
            CALL WRITE_OUT1
         ELSE IF (TIME + 0.1d0*DT>=OUT_TIME .OR. TIME+0.1d0*DT>=TSTOP) THEN
            OUT_TIME = (INT((TIME + 0.1d0*DT)/OUT_DT) + 1)*OUT_DT
            CALL WRITE_OUT1
         ENDIF
      ENDIF

! Write SPx files, if needed
      bWrite_netCDF_files = .false.
      ISPX = 0
      DO L = 1, N_SPX
         IF (DT == UNDEFINED) THEN
            IF (FINISH) THEN
               bWrite_netCDF_files = .true.
               CALL WRITE_SPX1 (L, 0)
               DISK_TOT = DISK_TOT + DISK(L)
               ISPX = ISPX + 1

               IF (SPX_MSG) THEN
                  IF (RES_MSG) THEN
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1001,  ADVANCE='NO') TIME
                     IF (FULL_LOG .and. myPE.eq.PE_IO) WRITE (*, 1001,  ADVANCE='NO') TIME
                  ELSE
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1002,  ADVANCE='NO')
                     IF (FULL_LOG .and. myPE.eq.PE_IO) WRITE (*, 1002,  ADVANCE='NO')
                  ENDIF
                  SPX_MSG = .FALSE.
               ENDIF
               IF(DMP_LOG)WRITE (UNIT_LOG, 1011,  ADVANCE='NO') EXT_END(L:L)
               IF (FULL_LOG .and. myPE.eq.PE_IO) WRITE (*, 1011,  ADVANCE='NO') EXT_END(L:L)
            ENDIF
! AEOLUS: stop trigger mechanism to terminate MFIX normally before batch
! queue terminates
         ELSEIF (TIME + 0.1d0*DT>=SPX_TIME(L) .OR. &
                 TIME+0.1d0*DT>=TSTOP.OR.eofBATCHQ) THEN
            SPX_TIME(L) = (INT((TIME + 0.1d0*DT)/SPX_DT(L))+1)*SPX_DT(L)
            CALL WRITE_SPX1 (L, 0)
            bWrite_netCDF_files = .true.
            DISK_TOT = DISK_TOT + DISK(L)
            ISPX = ISPX + 1
! remove this redundant call here to write_des_data in case of new
! coupled runs that contain at least a particle and involve some initial
! settling of the system.
! the call made in des_time_march is a better call for capturing the
! initial state of such a des continuum coupled system
            IF(DISCRETE_ELEMENT.AND.PRINT_DES_DATA .AND. L.EQ.1 .AND. &
               .NOT.(TRIM(RUN_TYPE)=='NEW' .AND. PARTICLES /=0 .AND. &
                     NFACTOR >0 .AND. TIME == ZERO)) THEN
                  CALL WRITE_DES_DATA
            ENDIF

            IF (SPX_MSG) THEN
               IF (RES_MSG) THEN
                  IF(DMP_LOG)WRITE (UNIT_LOG, 1001,  ADVANCE='NO') TIME
                  IF (FULL_LOG .and. myPE.eq.PE_IO) WRITE (*, 1001,  ADVANCE='NO') TIME
               ELSE
                  IF(DMP_LOG)WRITE (UNIT_LOG, 1002,  ADVANCE='NO')
                  IF (FULL_LOG .and. myPE.eq.PE_IO) WRITE (*, 1002,  ADVANCE='NO')
               ENDIF
               SPX_MSG = .FALSE.
            ENDIF
            IF(DMP_LOG)WRITE (UNIT_LOG, 1011,  ADVANCE='NO') EXT_END(L:L)
            IF (FULL_LOG .and. myPE.eq.PE_IO) &
               WRITE (*, 1011,  ADVANCE='NO') EXT_END(L:L)
         ENDIF
      ENDDO

      if (bWrite_netCDF_files) call write_netcdf(0,0,time)

      IF (.NOT.SPX_MSG) THEN
         DO L = 1, N_SPX - ISPX
            IF(DMP_LOG)WRITE (UNIT_LOG, '(A)', ADVANCE='NO') '   '
            IF (FULL_LOG .and. myPE.eq.PE_IO) &
               WRITE (*, '(A)', ADVANCE='NO') '   ' !//
         ENDDO
         IF(DMP_LOG)WRITE (UNIT_LOG, 1015) DISK_TOT
         IF (FULL_LOG.and.myPE.eq.PE_IO) WRITE (*, 1015) DISK_TOT !//
      ELSEIF (.NOT.RES_MSG) THEN
         IF(DMP_LOG)WRITE (UNIT_LOG, *)
         IF (FULL_LOG .and. myPE.eq.PE_IO) WRITE (*, *) !//
      ENDIF


! Write restart file, if needed
      CALL START_LOG
      IF (DT == UNDEFINED) THEN
         IF (FINISH) THEN
            CALL WRITE_RES1
            RES_MSG = .FALSE.
            IF(DMP_LOG)WRITE (UNIT_LOG, '(" t=",F10.4, "  Wrote RES;")', ADVANCE='NO') TIME
            IF (FULL_LOG .and. myPE.eq.PE_IO) THEN
               WRITE (*, 1000,  ADVANCE="NO") TIME
            ENDIF
         ENDIF
! AEOLUS: stop trigger mechanism to terminate MFIX normally before batch
! queue terminates
      ELSEIF (TIME + 0.1d0*DT>=RES_TIME .OR. TIME+0.1d0*DT>=TSTOP &
              .OR. eofBATCHQ) THEN
         RES_TIME = (INT((TIME + 0.1d0*DT)/RES_DT) + 1)*RES_DT
         CALL WRITE_RES1
         IF(DISCRETE_ELEMENT) CALL WRITE_RES0_DES
         IF(QMOMK) CALL QMOMK_WRITE_RESTART
         RES_MSG = .FALSE.
         IF(DMP_LOG)WRITE (UNIT_LOG, 1000,  ADVANCE='NO') TIME
         IF (FULL_LOG .and. myPE.eq.PE_IO) WRITE (*, 1000,  ADVANCE='NO') TIME
      ENDIF

      RES_MSG = .TRUE.
      SPX_MSG = .TRUE.
      CALL END_LOG

      CALL CPU_TIME(CPU1_IO)
      CPU_IO = CPU_IO + (CPU1_IO-CPU0_IO)

! Write special output, if needed
      DO L = 1, DIMENSION_USR
         IF (DT == UNDEFINED) THEN
            IF (FINISH) CALL WRITE_USR1 (L)
! AEOLUS: stop trigger mechanism to terminate MFIX normally before batch
! queue terminates
         ELSEIF (USR_TIME(L)/=UNDEFINED .AND. &
                 TIME+0.1d0*DT>=USR_TIME(L) .OR. eofBATCHQ) THEN
            USR_TIME(L) = (INT((TIME + 0.1d0*DT)/USR_DT(L))+1)*USR_DT(L)
            CALL WRITE_USR1 (L)
         ENDIF
      ENDDO

! JFD modification: cartesian grid implementation
! Write vtk file, if needed
      IF(WRITE_VTK_FILES) THEN
         DO L = 1, DIMENSION_VTK
            IF (DT == UNDEFINED) THEN
               IF (FINISH) THEN
                  CALL WRITE_VTU_FILE(L)
                  IF(DISCRETE_ELEMENT.AND.(DES_CONTINUUM_COUPLED))   CALL WRITE_VTP_FILE(L)
               ENDIF
            ELSEIF (VTK_TIME(L)/=UNDEFINED .AND.  TIME+0.1d0*DT>=VTK_TIME(L)) THEN
               VTK_TIME(L) = (INT((TIME + 0.1d0*DT)/VTK_DT(L))+1)*VTK_DT(L)
               CALL WRITE_VTU_FILE(L)
               IF(DISCRETE_ELEMENT.AND.(DES_CONTINUUM_COUPLED))   CALL WRITE_VTP_FILE(L)
            ENDIF
         ENDDO
      ENDIF

 1000 FORMAT(' t=',F10.4,'  Wrote RES;')
 1001 FORMAT(' t=',F10.4,'  Wrote      SPx:')
 1002 FORMAT(' SPx:')
 1010 FORMAT(I2,',')
 1011 FORMAT(A2,',')
 1015 FORMAT(14X,'Disk=',F7.2,' Mb')

      END SUBROUTINE OUTPUT_MANAGEMENT


      END SUBROUTINE TIME_MARCH

