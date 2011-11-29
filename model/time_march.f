!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: TIME_MARCH                                             C
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
!  kinetic equation                                                    C
!  Author: Alberto Passalacqua - Fox Research Group   Date: 02-Dec-09  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: RUN_TYPE, TIME, DT, NSTEP, TSTOP, OUT_DT,     C
!                        RES_DT, USR_DT, SPX_DT                        C
!                                                                      C
!  Variables modified: TIME, NSTEP                                     C
!                                                                      C
!  Local variables: OUT_TIME, RES_TIME, USR_TIME, SPX_TIME, USR_TIME,  C
!                   L, FINISH                                          C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE TIME_MARCH 
!...  Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...  Switches: -xf
!-----------------------------------------------
!     M o d u l e s 
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
      USE coeff
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
      USE mchem
      USE leqsol
! netcdf
      USE cdist
      USE MFIX_netcdf
! AEOLUS: stop trigger mechanism to terminate MFIX normally before
! batch queue terminates
      use mpi_utility
! JFD modification: cartesian grid implementation 
      USE cutcell
      USE vtk
      USE dashboard
! QMOMK - Alberto Passalacqua
      USE qmom_kinetic_equation
! QMOMK - End
      IMPLICIT NONE

!-----------------------------------------------
!     G l o b a l   P a r a m e t e r s
!-----------------------------------------------

!-----------------------------------------------
!     L o c a l   P a r a m e t e r s
!-----------------------------------------------
      DOUBLE PRECISION, PARAMETER :: ONEMEG = 1048576
  
!-----------------------------------------------
!     L o c a l   V a r i a b l e s
!-----------------------------------------------
     
! Flag to indicate one pass through iterate for steady 
! state conditions. 
      LOGICAL          FINISH 
     
! Time at which standard output is to be written 
      DOUBLE PRECISION OUT_TIME 
     
! Time at which restart file is to be written 
      DOUBLE PRECISION RES_TIME 
     
! Time at which REAL restart file is to be written 
      DOUBLE PRECISION SPX_TIME(N_SPX) 
     
! Disk space needed for one variable and each SPX file
      DOUBLE PRECISION DISK_ONE, DISK(N_SPX) 

! Total Disk space 
      DOUBLE PRECISION DISK_TOT 
     
! number SPX writes 
      INTEGER          ISPX 
      
      LOGICAL          RES_MSG, SPX_MSG 
     
! Time at which special output is to be written 
      DOUBLE PRECISION USR_TIME (DIMENSION_USR) 
     
! Loop indices 
      INTEGER          L, M , I
     
! Error index 
      INTEGER          IER 

! Number of iterations 
      INTEGER          NIT, NIT_TOTAL 
     
! used for activating check_data_30 
      INTEGER          NCHECK, DNCHECK,ijk 

! dummy logical variable for initializing adjust_dt 
      LOGICAL          dummy 

      CHARACTER        EXT_END*35 

! use function vavg_v_g to catch NaN's 
      DOUBLE PRECISION VAVG_U_G, VAVG_V_G, VAVG_W_G, X_vavg     

! use function MAX_VEL_INLET to compute max. velocity at inlet
      DOUBLE PRECISION MAX_VEL_INLET

! AEOLUS : stop trigger mechanism to terminate MFIX normally before
! batch queue terminates
      DOUBLE PRECISION CPU_STOP 
      LOGICAL AlreadyThere
      LOGICAL eofBATCHQ
! not used remove after verification
      INTEGER CHKBATCHQ_FLAG
      logical :: bWrite_netCDF_files
     
!-----------------------------------------------
!     E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      LOGICAL , EXTERNAL :: ADJUST_DT 
!-----------------------------------------------
     
      IF(AUTOMATIC_RESTART) RETURN
     
      EXT_END = '123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'
     
      FINISH  = .FALSE. 
      NCHECK  = NSTEP 
      DNCHECK = 1 
      CPU_IO  = ZERO 
      NIT_TOTAL = 0
! AEOLUS: stop before batch queue terminates
      AlreadyThere = .FALSE.
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
         VTK_TIME = UNDEFINED 
         IF (VTK_DT /= UNDEFINED) THEN 
            IF (RUN_TYPE == 'NEW') THEN 
               VTK_TIME = TIME 
            ELSE 
               VTK_TIME = (INT((TIME + 0.1d0*DT)/VTK_DT)+1)*VTK_DT
            ENDIF 
         ENDIF 
      ENDIF
 
! Parse residual strings
      CALL PARSE_RESID_STRING (IER) 

     
! Call user-defined subroutine to set constants, check data, etc.
      IF (CALL_USR) CALL USR0 
     
! CHEM & ISAT (Nan Xie)
      IF (CALL_DI) THEN 
         CALL MCHEM_INIT
         CALL MCHEM_ODEPACK_INIT
      ENDIF
      IF (CALL_ISAT) THEN 
         CALL MCHEM_INIT
         CALL MCHEM_ODEPACK_INIT
         CALL MISAT_TABLE_INIT
      ENDIF


      CALL RRATES_INIT(IER)

      DO M=1, MMAX 
         CALL ZERO_ARRAY (F_gs(1,M), IER)
      ENDDO

! DES 
! This call to make_arrays_des has now been moved ahead of calc_coeff_all 
! so that on the call to des/drag_fgs.f, the particle in cell info and also
! Ep_s are known.  Rahul Garg 
!      IF(DISCRETE_ELEMENT) THEN
!         CALL MAKE_ARRAYS_DES
!      END IF
! rahul: make_arrays_des (along with check_des_data) is now
! called from mfix.f. This consolidated the various DEM
! intialization related calls to one place.  
      IF (QMOMK) THEN
          CALL QMOMK_MAKE_ARRAYS
      END IF
! QMOMK - End

! Calculate all the coefficients once before entering the time loop
      IF(TRIM(KT_TYPE) == UNDEFINED_C) CALL CALC_COEFF_ALL (0, IER) 

     
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
      IF (SHEAR) THEN
         call CAL_D(V_sh)
      ENDIF

! Initialize check_mass_balance.  This routine is not active by default.
! Specify a reporting interval (hard-wired in the routine) to activate 
! the routine.
      Call check_mass_balance (0)

! sof modification: now it's only needed to do this once before time-loop     
! Mark the phase whose continuity will be used for forming Pp_g and Pp_s eqs.
      CALL MARK_PHASE_4_COR (PHASE_4_P_G, PHASE_4_P_S, DO_CONT, MCP,&
          DO_P_S, SWITCH_4_P_G, SWITCH_4_P_S, IER) 

 

! The TIME loop begins here.............................................
 100  CONTINUE
      
      IF(DISCRETE_ELEMENT.AND.(.NOT.DES_CONTINUUM_COUPLED))  THEN 
         IF(WRITE_VTK_FILES) THEN
! rahul: in order to write vtk files for cut-cell
            CALL WRITE_VTK_FILE
         ENDIF 
         CALL DES_TIME_MARCH
         call CPU_TIME(CPU_STOP)
! pradeep added  finalize
         CPU_STOP = CPU_STOP - CPU00
         if(mype.eq.pe_io)write(*,"('Elapsed CPU time = ',E15.6,' sec')") CPU_STOP
         call parallel_fin  
         STOP
      ENDIF 

! AEOLUS: stop trigger mechanism to terminate MFIX normally before batch 
! queue terminates
      IF (CHK_BATCHQ_END) THEN
         CHKBATCHQ_FLAG = 0 
         IF (myPE.eq.PE_IO) THEN
            CALL CPU_TIME(CPU_STOP)
! need to use CPU00, a timestamp from first line of mfix.f to take 
! account the time spent in I/O
            CPU_STOP = CPU_STOP - CPU00
            write(*,"('Elapsed CPU time = ',E15.6,' sec')") CPU_STOP
        
            IF ((CPU_STOP+TERM_BUFFER) .ge. BATCH_WALLCLOCK) THEN
               write(*,'(/,A,A)') '=============== REQUESTED CPU ',&
                  'TIME LIMIT REACHED ==========='
               write(*,*) 'Elapsed CPU time                        = ',&
                  CPU_STOP,' sec'
               write(*,*) 'Buffer CPU time before triggering abort = ',&
                  TERM_BUFFER,' sec'
               write(*,*) 'Elapsed+Buffer CPU time = ',&
                  (CPU_STOP+TERM_BUFFER), &
                  ' sec >= Allocated Wallclock ',&
                  BATCH_WALLCLOCK, ' sec'
!               write(*,'(A,A,E10.2,A)') 'Buffer CPU time before ',&
!                  'triggering abort = ', TERM_BUFFER, ' sec'
!               write(*,'(A,E10.2,A,E10.2,A)') &
!                  'Elapsed+Buffer CPU time = ', (CPU_STOP+TERM_BUFFER), &
!                  ' sec >= Allocated Wallclock ', &
!                  BATCH_WALLCLOCK, ' sec'
               write(*,'(A,A,/)') '=============== REQUESTED CPU ',&
                  'TIME LIMIT REACHED ==========='
               eofBATCHQ = .TRUE.
               CHKBATCHQ_FLAG = 1
!              FINISH = .TRUE.	  
!              DT = UNDEFINED
            ENDIF 
        
            INQUIRE(file="MFIX.STOP",exist=AlreadyThere)
            IF (AlreadyThere) THEN
               write(*,'(/,A,A)') '=============== MFIX STOP SIGNAL ',&
                  'DETECTED ==========='
               write(*,'(A,A)') '  MFIX.STOP file detected in ',&
                  'working directory, terminating MFIX run'
               write(*,'(A,A)') '  Please DO NOT FORGET to erase ',&
                  'MFIX.STOP file before next run'
!               write(*,'(A,E10.2,A)') '  Elapsed CPU time = ',&
!                 CPU_STOP, ' sec'
               write(*,*) ' Elapsed CPU time = ',CPU_STOP,' sec'
               write(*,'(A,A,/)') '=============== MFIX STOP ',&
                  'SIGNAL DETECTED ==========='
               eofBATCHQ = .TRUE.
               CHKBATCHQ_FLAG = 1
!              FINISH = .TRUE.
!              DT = UNDEFINED          
               AlreadyThere = .FALSE.
            ENDIF
         ENDIF     ! myPE = PE_IO
! Try to move this bcast call to another location where there is barrier or bcast 
         call bcast (eofBATCHQ,PE_IO)
      ENDIF   

      IF (CALL_USR) CALL USR1 

! Remove solids from cells containing very small quantities of solids
      IF(.NOT.(DISCRETE_ELEMENT .OR. QMOMK)) THEN
         IF(TRIM(KT_TYPE) == 'GHD') THEN
            CALL ADJUST_EPS_GHD
         ELSE 
            CALL ADJUST_EPS
         ENDIF
      ENDIF

     
! sof modification: uncomment code below and modify MARK_PHASE_4_COR to use
! previous MFIX algorithm. Oct. 18 2010.
! Mark the phase whose continuity will be used for forming Pp_g and Pp_s eqs.
!      CALL MARK_PHASE_4_COR (PHASE_4_P_G, PHASE_4_P_S, DO_CONT, MCP, DO_P_S, &
!          SWITCH_4_P_G, SWITCH_4_P_S, IER) 

! Set wall boundary conditions and transient flow b.c.'s
      CALL SET_BC1 

      CALL CPU_TIME(CPU0_IO)
     
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
            IF(DISCRETE_ELEMENT.AND.L.EQ.1.AND.PRINT_DES_DATA) CALL WRITE_DES_DATA
     
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
            IF(DMP_LOG)WRITE (UNIT_LOG, '(A,$)') '   ' 
            IF (FULL_LOG .and. myPE.eq.PE_IO) WRITE (*, '(A,$)') '   ' !//
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
! Pradeep Changing the routine name 
         IF(DISCRETE_ELEMENT) call des_write_restart !CALL WRITE_DES_RESTART
	 ! QMOMK - Alberto Passsalacqua
	 ! Storing QMOM restart data for the solution of the Boltzmann equation
         IF(QMOMK) CALL QMOMK_WRITE_RESTART         
         ! QMOMK - End

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
            IF (FINISH.and.myPE.eq.PE_IO) CALL WRITE_USR1 (L) 
! AEOLUS: stop trigger mechanism to terminate MFIX normally before batch
! queue terminates
         ELSEIF (USR_TIME(L)/=UNDEFINED .AND. &
                 TIME+0.1d0*DT>=USR_TIME(L) .OR. eofBATCHQ) THEN 
            USR_TIME(L) = (INT((TIME + 0.1d0*DT)/USR_DT(L))+1)*USR_DT(L)
            if (myPE.eq.PE_IO) CALL WRITE_USR1 (L) 
         ENDIF 
      ENDDO 

! JFD modification: cartesian grid implementation
! Write vtk file, if needed
      IF(WRITE_VTK_FILES) THEN
         IF (DT == UNDEFINED) THEN 
            IF (FINISH) CALL WRITE_VTK_FILE
         ELSEIF (VTK_TIME/=UNDEFINED .AND. TIME+0.1d0*DT>=VTK_TIME) THEN
            VTK_TIME = (INT((TIME + 0.1d0*DT)/VTK_DT)+1)*VTK_DT 
            CALL WRITE_VTK_FILE
         ENDIF 
      ENDIF

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
     
! Calculate the trace of the stress tensor
      CALL CALC_TRD_G (TRD_G, IER) 
      IF (.NOT.DISCRETE_ELEMENT) CALL CALC_TRD_S (TRD_S, IER)

! Calculate the cross terms of the stress tensor
      CALL CALC_TAU_U_G (TAU_U_G, IER) 
      CALL CALC_TAU_V_G (TAU_V_G, IER) 
      CALL CALC_TAU_W_G (TAU_W_G, IER) 
      IF (.NOT.DISCRETE_ELEMENT) THEN
         CALL CALC_TAU_U_S (TAU_U_S, IER) 
         CALL CALC_TAU_V_S (TAU_V_S, IER) 
         CALL CALC_TAU_W_S (TAU_W_S, IER) 
      ENDIF

! Calculate additional solid phase momentum source terms 
! that arise from kinetic theory constitutive relations
      IF (.NOT.DISCRETE_ELEMENT) THEN
         CALL CALC_KTMOMSOURCE_U_S (IER)
         CALL CALC_KTMOMSOURCE_V_S (IER)
         CALL CALC_KTMOMSOURCE_W_S (IER)
      ENDIF

! Check rates and sums of mass fractions every NLOG time steps
      IF (NSTEP == NCHECK) THEN 
         IF (DNCHECK < 256) DNCHECK = DNCHECK*2 
         NCHECK = NCHECK + DNCHECK 
! Upate the reaction rates for checking
         IF (ANY_SPECIES_EQ) RRATE = .TRUE. 
         CALL CALC_RRATE(RRATE)
         CALL CHECK_DATA_30 
      ENDIF 

! AE TIME 041601 Double the timestep for 2nd order accurate time implementation
! IF ((CN_ON.AND.NSTEP>1)) THEN
      IF ((CN_ON.AND.NSTEP>1.AND.RUN_TYPE == 'NEW') .OR. & 
      (CN_ON.AND.RUN_TYPE /= 'NEW' .AND. NSTEP >= (NSTEPRST+1))) THEN
         DT = 0.5d0*DT
         ODT = ODT * 2.0d0
      ENDIF
      
! Check for maximum velocity at inlet to avoid convergence problems 
      MAX_INLET_VEL = 100.0d0*MAX_VEL_INLET ()
! if no inlet velocity is specified, use an upper limit defined in 
! toleranc_mod.f
      IF(MAX_INLET_VEL == ZERO) THEN
         MAX_INLET_VEL = MAX_ALLOWED_VEL
         IF (UNITS == 'SI') MAX_INLET_VEL = 1D-2 * MAX_ALLOWED_VEL
      ENDIF
! Scale the value using a user defined scale factor      
      MAX_INLET_VEL = MAX_INLET_VEL * MAX_INLET_VEL_FAC
      
     
! Advance the solution in time by iteratively solving the equations 
      CALL ITERATE (IER, NIT)
      IF(AUTOMATIC_RESTART) RETURN
      
! Just to Check for NaN's, Uncomment the following lines and also lines  
! of code in  VAVG_U_G, VAVG_V_G, VAVG_W_G to use.
!      X_vavg = VAVG_U_G ()
!      X_vavg = VAVG_V_G ()
!      X_vavg = VAVG_W_G ()
!      IF(AUTOMATIC_RESTART) RETURN

     
! Adjust time step and reiterate if necessary
      DO WHILE (ADJUST_DT(IER,NIT))
         CALL ITERATE (IER, NIT) 
      ENDDO

      IF(DT.LT.DT_MIN) THEN
         IF(TIME.LE.RES_DT) THEN
            IF (AUTO_RESTART .AND. DMP_LOG)WRITE(UNIT_LOG,*) &
               'Automatic restart not possible as Total Time < RES_DT'

! JFD modification: cartesian grid implementation
            IF(WRITE_DASHBOARD) THEN
               RUN_STATUS = 'DT < DT_MIN.  Recovery not possible!'
               CALL UPDATE_DASHBOARD(NIT,0.0d0,'    ')
            ENDIF

            CALL MFIX_EXIT(MyPE)
         ENDIF
         IF(AUTO_RESTART) AUTOMATIC_RESTART = .TRUE.
         RETURN
      ENDIF

      
! Check over mass and elemental balances.  This routine is not active by default.
! Edit the routine and specify a reporting interval to activate it.
      Call check_mass_balance (1)

! DES
      IF (DISCRETE_ELEMENT.AND.DES_CONTINUUM_COUPLED) CALL DES_TIME_MARCH

! QMOMK - Alberto Passalacqua
      IF (QMOMK) CALL QMOMK_TIME_MARCH
!     END QMOMB

! CHEM & ISAT (Nan Xie)
! Advance the time step and continue
      IF (CALL_ISAT .OR. CALL_DI) THEN 
         CALL MCHEM_TIME_MARCH
      ENDIF
      IF (CALL_DQMOM) CALL USR_DQMOM


! Advance the time step and continue
! AE TIME 041601 Double the timestep for 2nd order accurate time implementation
! IF (CN_ON.AND.NSTEP>1) then 
      IF ((CN_ON.AND.NSTEP>1.AND.RUN_TYPE == 'NEW') .OR. & 
      (CN_ON.AND.RUN_TYPE /= 'NEW' .AND. NSTEP >= (NSTEPRST+1))) THEN
         DT = 2.d0*DT      
         ODT = ODT * 0.5d0
      ENDIF      

! AE TIME 043001 Perform the explicit extrapolation for CN implementation
! IF (CN_ON.AND.NSTEP>1) then
      IF ((CN_ON.AND.NSTEP>1.AND.RUN_TYPE == 'NEW') .OR. & 
      (CN_ON.AND.RUN_TYPE /= 'NEW' .AND. NSTEP >= (NSTEPRST+1))) THEN
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
      
      CALL FLUSH (6) 
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


!     1000  FORMAT(' t=',F10.4, '  Wrote RES;',$)
!     1001  FORMAT(' t=',F10.4, '  Wrote      SPx:',$)
!     1002  FORMAT(' SPx:',$)
!     1010  FORMAT(I2,',',$)
!     1015  FORMAT(14X, 'Disk=', F7.2,' Mb')
 1000 FORMAT(' t=',F10.4,'  Wrote RES;') 
 1001 FORMAT(' t=',F10.4,'  Wrote      SPx:') 
 1002 FORMAT(' SPx:') 
 1010 FORMAT(I2,',') 
 1011 FORMAT(A2,',') 
 1015 FORMAT(14X,'Disk=',F7.2,' Mb') 

      END SUBROUTINE TIME_MARCH 
     
