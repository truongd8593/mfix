! -*- f90 -*-
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  MODULE: MAIN                                                        !
!                                                                      !
!  Purpose: Main module for top level mfix subroutines.                !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      MODULE MAIN

      use exit, only: mfix_exit

!-----------------------------------------------
! Module variables
!-----------------------------------------------
! Final value of CPU time.
      DOUBLE PRECISION :: CPU1
! time used for computations.
      DOUBLE PRECISION :: CPUTIME_USED, WALLTIME_USED
! DISTIO variable for specifying the mfix version
      CHARACTER(LEN=512) :: version
! environment variable
!$ CHARACTER(LEN=512) :: omp_num_threads
!$ INTEGER :: length
!$ INTEGER :: status

! Number of iterations
      INTEGER :: NIT_TOTAL
! used for activating check_data_30
      INTEGER :: NCHECK, DNCHECK

! Flag to save results and cleanly exit.
      LOGICAL :: EXIT_SIGNAL = .FALSE.

      CHARACTER(LEN=80), DIMENSION(100) :: CMD_LINE_ARGS
      INTEGER :: CMD_LINE_ARGS_COUNT = 0
      CHARACTER(LEN=80) :: MFIX_DAT = "mfix.dat"

      CONTAINS

      SUBROUTINE INITIALIZE

#ifdef MPI
      USE mpi, only: mpi_comm_world, mpi_barrier
#endif
      USE cdist, only: bdoing_postmfix
      USE cdist, only: bglobalnetcdf, bstart_with_one_res, bdist_io, bwrite_netcdf
      USE check, only: check_mass_balance
      USE check_data_cg, only: check_bc_flags, report_best_processor_size
      USE coeff, only: init_coeff
      USE compar, only: mpierr, mype, pe_io
      USE compar, only: mype, pe_io
      USE cont, only: do_cont
      USE cutcell, only: cartesian_grid, re_indexing, set_corner_cells
      USE dbg, only: debug_write_layout, write_parallel_info
      USE discretelement, only: discrete_element
      USE drag, only: f_gs
      USE error_manager, only: err_msg, flush_err_msg
      USE error_manager, only: init_err_msg, finl_err_msg
      USE fldvar, only: rop_g, rop_s
      USE funits, only: dmp_log, unit_log, unit_res
      USE machine, only: start_log, end_log
      USE machine, only: wall_time, pc_quickwin, machine_cons, get_run_id, start_log, end_log
      USE mfix_netcdf, only: mfix_usingnetcdf
      USE output, only: dbgprn_layout
      USE output_man, only: init_output_vars, output_manager
      USE parallel_mpi, only: parallel_init, parallel_fin
      USE param1, only: n_spx, undefined, zero
      USE pgcor, only: d_e, d_n, d_t, phase_4_p_g, switch_4_p_g
      USE physprop, only: mmax
      USE pscor, only: e_e, e_n, e_t, do_p_s, phase_4_p_s, mcp, switch_4_p_s
      USE qmom_kinetic_equation, only: qmomk
      USE run, only: call_usr, dem_solids, dt_max, dt_min
      USE run, only: id_version
      USE run, only: ier
      USE run, only: nstep, pic_solids, run_type, dt, shear, time, v_sh
      USE run, only: time
      USE time_cpu, only: CPU00, wall0
      USE time_cpu, only: cpu_io, cpu_nlog, cpu0, cpuos, time_nlog
      USE vtk, only: write_vtk_files

      IMPLICIT NONE

      ! Temporary storage for DT
     DOUBLE PRECISION :: DT_tmp
      ! Save TIME in input file for RESTART_2
     DOUBLE PRECISION :: TIME_SAVE

     INTEGER :: LL, MM

!$    INTEGER num_threads, threads_specified, omp_id
!$    INTEGER omp_get_num_threads
!$    INTEGER omp_get_thread_num

! DISTIO
! If you change the value below in this subroutine, you must also
! change it in write_res0.f and the value should also be consistent
! with the check in read_res0
      version = 'RES = 01.6'

      bDoing_postmfix = .false.

! Invoke MPI initialization routines and get rank info.
      CALL PARALLEL_INIT
      CALL GEN_LOG_BASENAME

! we want only PE_IO to write out common error messages
      DMP_LOG = (myPE == PE_IO)

! set the version.release of the software
      ID_VERSION = '2016-1'

! set automatic restart flag to false
!      AUTOMATIC_RESTART = .FALSE.
!      ITER_RESTART      = 1

! specify the number of processors to be used
!$        call get_environment_variable("OMP_NUM_THREADS",omp_num_threads,length,status, .true.)
!$      if (status.eq.0 .and. length.ne.0) then
!$        read(omp_num_threads,*) threads_specified
!$      else
!$        WRITE(*,'(A)') 'Enter the number of threads to be used for SMP: '
!$        READ(*,*) threads_specified
!$      endif

!$      call omp_set_num_threads(threads_specified)

! Find the number of processors used
!$omp  parallel
!$      num_threads = omp_get_num_threads()
!$      omp_id = omp_get_thread_num()
!$      if(omp_id.eq.0) Write(*,*)' Number of threads used for SMP = ',  num_threads
!$omp  end parallel

! Set machine dependent constants
      CALL MACHINE_CONS

! Get the date and time. They give the unique run_id in binary output
! files
      CALL GET_RUN_ID

! AEOLUS: stop trigger mechanism to terminate MFIX normally before batch
! queue terminates. timestep at the beginning of execution
      CALL CPU_TIME (CPU00)
      WALL0 = WALL_TIME()

! Read input data, check data, do computations for IC and BC locations
! and flows, and set geometry parameters such as X, X_E, DToDX, etc.
      CALL GET_DATA

! Check data, do computations for IC and BC locations
! and flows, and set geometry parameters such as X, X_E, DToDX, etc.
      CALL CHECK_DATA

! Write the initial part of the standard output file
      CALL WRITE_OUT0
      IF(.NOT.CARTESIAN_GRID)  CALL WRITE_FLAGS

! Write the initial part of the special output file(s)
      CALL WRITE_USR0

!$    CALL START_LOG
!$    IF(DMP_LOG)WRITE (UNIT_LOG, *) ' '
!$    IF(DMP_LOG)WRITE (UNIT_LOG, *) ' Number of processors used = ', threads_specified
!$    IF(DMP_LOG)WRITE (UNIT_LOG, *) ' '
!$    CALL END_LOG

!  setup for PC quickwin application
      CALL PC_QUICKWIN

      CALL INIT_ERR_MSG('MFIX')


! if not netcdf writes asked for ... globally turn off netcdf
      if(MFIX_usingNETCDF()) then
         bGlobalNetcdf = .false.
         do LL = 1,20
            if (bWrite_netcdf(LL)) bGlobalNetcdf = .true.
         enddo
      endif

      DT_TMP = DT
      SELECT CASE (TRIM(RUN_TYPE))

      CASE ('NEW')
! Write the initial part of the restart files
         CALL WRITE_RES0
         DO LL = 1, N_SPX
            CALL WRITE_SPX0 (LL, 0)
         ENDDO

      CASE ('RESTART_1')
! Read the time-dependent part of the restart file
         CALL READ_RES1
         WRITE(ERR_MSG, 1010) TIME, NSTEP
         CALL FLUSH_ERR_MSG()

      CASE ('RESTART_2')
         TIME_SAVE = TIME
! DISTIO
         if (myPE .ne. PE_IO .and. bDist_IO .and. bStart_with_one_res) then
            write (unit_res,rec=1) version
            write (unit_res,rec=2) 4
            write (unit_res,rec=3) 4
         endif

         CALL READ_RES1
         TIME = TIME_SAVE

1010     FORMAT('Message 1010: Read in data from .RES file for TIME = ',&
              G12.5,/'Time step number (NSTEP) =',I7)

         WRITE(ERR_MSG, 1010) TIME, NSTEP
         CALL FLUSH_ERR_MSG()

         CALL WRITE_RES0

! Writing the RES1 and SPX1 can only be done here when re-indexing is turned off
! This will be done after the cell re-indexing is done later in this file.
! This allows restarting independently of the re-indexing setting between
! the previous and current run.
         IF(.NOT.RE_INDEXING) THEN
            CALL WRITE_RES1
            DO LL = 1, N_SPX
               CALL WRITE_SPX0 (LL, 0)
               CALL WRITE_SPX1 (LL, 0)
            END DO
            call write_netcdf(0,0,time)
         ENDIF

      CASE DEFAULT
         CALL START_LOG
         IF(DMP_LOG)WRITE (UNIT_LOG, *) &
              ' MFIX: Do not know how to process'
         IF(DMP_LOG)WRITE (UNIT_LOG, *) ' RUN_TYPE in data file'
         CALL END_LOG
         call mfix_exit(myPE)

      END SELECT

#ifdef MPI
      call MPI_Barrier(MPI_COMM_WORLD,mpierr)
#endif

      IF (DT_TMP /= UNDEFINED) THEN
         DT = MAX(DT_MIN,MIN(DT_MAX,DT))
      ELSE
         DT = DT_TMP
      ENDIF

! Set arrays for computing indices. A secondary call is made
! after cut cell-preprocessing to update array indices.
      IF(CARTESIAN_GRID) THEN
         CALL SET_INCREMENTS
         CALL SET_INCREMENTS3
      ENDIF

!      IF(.NOT.RE_INDEXING) CALL WRITE_IJK_VALUES

! Set the flags for wall surfaces impermeable and identify flow
! boundaries using FLAG_E, FLAG_N, and FLAG_T
      CALL SET_FLAGS1

!  Update flags for Cartesian_GRID.
      IF(CARTESIAN_GRID) CALL CHECK_BC_FLAGS

! Calculate cell volumes and face areas
      IF(.NOT.CARTESIAN_GRID) CALL SET_GEOMETRY1

! Find corner cells and set their face areas to zero
      IF(.NOT.CARTESIAN_GRID)  THEN
         CALL GET_CORNER_CELLS()
      ELSE
         IF (SET_CORNER_CELLS)  CALL GET_CORNER_CELLS ()
      ENDIF

! Set constant physical properties
      CALL SET_CONSTPROP

! Set initial conditions
      CALL SET_IC

! Set point sources.
      CALL SET_PS

! Set boundary conditions
      CALL ZERO_NORM_VEL
      CALL SET_BC0

! Cartesian grid implementation
      IF(CARTESIAN_GRID) CALL CG_SET_BC0

! Set gas mixture molecular weight
      CALL SET_MW_MIX_G

! Set the pressure field for a fluidized bed
      IF (RUN_TYPE == 'NEW') CALL SET_FLUIDBED_P

! Initialize densities.
      IF (RUN_TYPE == 'NEW') CALL SET_RO_G
      IF (RUN_TYPE == 'NEW') CALL SET_RO_S

! Initialize time dependent boundary conditions
      CALL SET_BC1

! Check the field variable data and report errors.
      IF(.NOT.CARTESIAN_GRID)  CALL CHECK_DATA_20

!=======================================================================
! JFD: START MODIFICATION FOR RE-INDEXING CELLS
!=======================================================================
      IF(CARTESIAN_GRID.AND.RE_INDEXING) THEN

         IF(myPE == PE_IO) THEN
            WRITE(*,"(72('='))")
            WRITE(*,*)' RE-INDEXING CELLS FOR CARTESIAN GRID...'
         ENDIF
         CALL RE_INDEX_ARRAYS


         !IF(myPE == PE_IO)print*,'Calling REPORT_BEST_IJK_SIZE:'
         !CALL REPORT_BEST_IJK_SIZE
         CALL REPORT_BEST_PROCESSOR_SIZE
         !IF(myPE == PE_IO)print*,'Exiting MFIX after REPORT_BEST_IJK_SIZE.'

         IF(myPE == PE_IO) WRITE(*,"(72('='))")

! In case of a RESTART_2, write the RES1 and SPX1 files here
! This was commented out earlier in this file.
         IF(RUN_TYPE == 'RESTART_2') THEN
            CALL WRITE_RES1
            DO LL = 1, N_SPX
               CALL WRITE_SPX0 (LL, 0)
               CALL WRITE_SPX1 (LL, 0)
            END DO
            call write_netcdf(0,0,time)
         ENDIF
      ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR RE-INDEXING CELLS
!=======================================================================

! Setup VTK data for regular (no cut cells) grid
      IF(.NOT.CARTESIAN_GRID.AND.WRITE_VTK_FILES) CALL SETUP_VTK_NO_CUTCELL

      IF(DISCRETE_ELEMENT) CALL MAKE_ARRAYS_DES
      IF(QMOMK) CALL QMOMK_MAKE_ARRAYS

! Set the inflow/outflow BCs for DEM solids
      IF(DEM_SOLIDS) CALL SET_BC_DEM
! Set the inflow/outflow BC for PIC solids
      IF(PIC_SOLIDS) CALL SET_BC_PIC

! Set the inital properties of each particle.
      IF(DEM_SOLIDS) CALL SET_IC_DEM

! AEOLUS: debug prints
      if (DBGPRN_LAYOUT .or. bdist_io) then
         !write (*,*) myPE , ' E.4 ... version = ' , version(1:33)
         call debug_write_layout()
         call write_parallel_info()
      endif

! Initializations for CPU time calculations in iterate
      CPUOS = 0.
      CALL CPU_TIME (CPU1)
      CPU_NLOG = CPU1
      TIME_NLOG = TIME - DT

! Get the initial value of CPU time
      CALL CPU_TIME (CPU0)

! Find the solution of the equations from TIME to TSTOP at
! intervals of DT

!-----------------------------------------------

      NCHECK  = NSTEP
      DNCHECK = 1
      CPU_IO  = ZERO
      NIT_TOTAL = 0

      CALL INIT_OUTPUT_VARS

! Parse residual strings
      CALL PARSE_RESID_STRING ()

! Call user-defined subroutine to set constants, check data, etc.
      IF (CALL_USR) CALL USR0

      CALL RRATES_INIT()

! Calculate all the coefficients once before entering the time loop
      CALL INIT_COEFF(IER)

      DO MM=1, MMAX
         F_gs(:,MM) = ZERO
      ENDDO

! Remove undefined values at wall cells for scalars
      CALL UNDEF_2_0 (ROP_G)
      DO MM = 1, MMAX
         CALL UNDEF_2_0 (ROP_S(:,MM))
      ENDDO

! Initialize d's and e's to zero
      DO MM = 0, MMAX
         D_E(:,MM) = ZERO
         D_N(:,MM) = ZERO
         D_T(:,MM) = ZERO
      ENDDO
      E_E(:) = ZERO
      E_N(:) = ZERO
      E_T(:) = ZERO

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
           DO_P_S, SWITCH_4_P_G, SWITCH_4_P_S)

      END SUBROUTINE INITIALIZE

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  SUBROUTINE: GET_DATA                                                C
!  Purpose: read and verify input data, open files                     C
!                                                                      C
!  Author: P. Nicoletti                               Date: 04-DEC-91  C
!  Reviewer: M.SYAMLAL, W.ROGERS, P.NICOLETTI         Date: 24-JAN-92  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE GET_DATA

!-----------------------------------------------
! Modules
!-----------------------------------------------

      USE compar, only: adjust_partition, mype, nodesi, nodesj, nodesk, pe_io
      USE mpi_utility, only: bcast
      USE run, only: run_type

      IMPLICIT NONE

      LOGICAL :: PRESENT

! This module call routines to initialize the namelist variables.
      CALL INIT_NAMELIST
! Read in the namelist variables from the ascii input file.
      CALL READ_NAMELIST(0,MFIX_DAT)
! Set RUN_TYPE to RESTART_1 when adjusting partition
! and read partition layout in gridmap.dat if it exists
      IF(ADJUST_PARTITION) THEN
         RUN_TYPE = 'RESTART_1'

         INQUIRE(FILE='gridmap.dat',EXIST=PRESENT)
         IF(PRESENT) THEN
            IF(MyPE == PE_IO) THEN
               WRITE(*,*)'Reading partition layout from grimap.dat...'
               OPEN(UNIT=777, FILE='gridmap.dat', STATUS='OLD')

                READ (777, *) NODESI,NODESJ,NODESK

                CLOSE(777)
            ENDIF

            CALL BCAST(NODESI)
            CALL BCAST(NODESJ)
            CALL BCAST(NODESK)
         ENDIF

      ENDIF

      RETURN

    END SUBROUTINE GET_DATA

      SUBROUTINE CHECK_DATA
         USE check_data_cg, only: adjust_ijk_size, check_data_cartesian
         USE constant, only: set_constants
         USE cut_cell_preproc, only: cut_cell_preprocessing
         USE cutcell, ONLY: cartesian_grid
         USE dashboard, ONLY: dtmax, dtmin, init_time, n_dashboard, nit_max, nit_min, smmax, smmin
         USE des_allocate
         USE desgrid, only: DESGRID_INIT
         USE discretelement, ONLY: discrete_element
         USE error_manager
         USE gridmap
         USE iterate, ONLY: max_nit
         USE mpi_init_des, only: DESMPI_INIT
         USE param1, ONLY: large_number, n_spx
         USE qmom_kinetic_equation, ONLY: qmomk
         USE run, ONLY: run_name, run_type, time
         USE stl_preproc_des, only: DES_STL_PREPROCESSING
         IMPLICIT NONE

         !-----------------------------------------------
         ! Local variables
         !-----------------------------------------------
         ! shift DX, DY and DZ values
         LOGICAL, PARAMETER :: SHIFT = .TRUE.

! Initialize the error manager. This call occurs after the mfix.dat
! is read so that message verbosity can be set and the .LOG file
! can be opened.
      CALL INIT_ERROR_MANAGER

! Write header in the .LOG file and to screen.
! Not sure if the values are correct or useful
      CALL WRITE_HEADER

! Open files
      CALL OPEN_FILES(RUN_NAME, RUN_TYPE, N_SPX)

! These checks verify that sufficient information was provided
! to setup the domain indices and DMP gridmap.
      CALL CHECK_GEOMETRY_PREREQS
      CALL CHECK_DMP_PREREQS

! Set up the physical domain indicies (cell index max/min values).
      CALL SET_MAX2

! Set constants
      CALL SET_CONSTANTS

! Adjust partition for better load balance (done when RE_INDEXING is .TRUE.)
      CALL ADJUST_IJK_SIZE

! Partition the domain and set indices
      CALL GRIDMAP_INIT

! Check the minimum solids phase requirements.
      CALL CHECK_SOLIDS_MODEL_PREREQS

      CALL CHECK_RUN_CONTROL
      CALL CHECK_NUMERICS
      CALL CHECK_OUTPUT_CONTROL

      CALL CHECK_GAS_PHASE
      CALL CHECK_SOLIDS_PHASES
      CALL SET_PARAMETERS

! Basic geometry checks.
      CALL CHECK_GEOMETRY(SHIFT)
      IF(DISCRETE_ELEMENT) CALL CHECK_GEOMETRY_DES

! Set grid spacing variables.
      CALL SET_GEOMETRY
      IF(DISCRETE_ELEMENT) CALL SET_GEOMETRY_DES

      CALL CHECK_INITIAL_CONDITIONS
      CALL CHECK_BOUNDARY_CONDITIONS
      CALL CHECK_INTERNAL_SURFACES
      CALL CHECK_POINT_SOURCES

      CALL CHECK_CHEMICAL_RXNS
      CALL CHECK_ODEPACK_STIFF_CHEM



!----------------------  DOMAIN SPECIFIC CHECKS  --------------------!


! This call needs to occur before any of the IC/BC checks.
      CALL SET_ICBC_FLAG

! Compute area of boundary surfaces.
      CALL GET_BC_AREA

! Convert (mass, volume) flows to velocities.
      CALL SET_BC_FLOW

! Set the flags for identifying computational cells
      CALL SET_FLAGS
! Set arrays for computing indices
      CALL SET_INCREMENTS
      CALL SET_INCREMENTS3

! Cartesian grid implementation
      CALL CHECK_DATA_CARTESIAN
      IF(CARTESIAN_GRID) THEN
         CALL CUT_CELL_PREPROCESSING
      ELSE
         CALL ALLOCATE_DUMMY_CUT_CELL_ARRAYS
      ENDIF

      IF(DISCRETE_ELEMENT) THEN
         CALL DESGRID_INIT
         CALL DESMPI_INIT
         CALL DES_STL_PREPROCESSING
      ENDIF

!--------------------------  ARRAY ALLOCATION -----------------------!

! Allocate array storage.
      CALL ALLOCATE_ARRAYS
      IF(DISCRETE_ELEMENT) CALL DES_ALLOCATE_ARRAYS
      IF(QMOMK) CALL QMOMK_ALLOCATE_ARRAYS

! Initialize arrays.
      CALL INIT_FVARS
      IF(DISCRETE_ELEMENT) CALL DES_INIT_ARRAYS

!======================================================================
! Data initialization for Dashboard
!======================================================================
      INIT_TIME = TIME
      SMMIN =  LARGE_NUMBER
      SMMAX = -LARGE_NUMBER

      DTMIN =  LARGE_NUMBER
      DTMAX = -LARGE_NUMBER

      NIT_MIN = MAX_NIT
      NIT_MAX = 0

      N_DASHBOARD = 0


      END SUBROUTINE CHECK_DATA

      SUBROUTINE FINALIZE

      USE cutcell, only: cartesian_grid
      USE dashboard
      USE cut_cell_preproc, only: close_cut_cell_files
      USE error_manager, only: finl_err_msg
      USE machine, only: wall_time
      USE parallel_mpi, only: parallel_fin
      USE run, only: dt, call_usr, dt_min, get_tunit, tunit,time,tstop
      USE time_cpu
      USE compar, only:ADJUST_PARTITION
      USE output_man
      IMPLICIT NONE

! Call user-defined subroutine after time-loop.
      IF (CALL_USR) CALL USR3

! Get the final value of CPU time.  The difference gives the
! CPU time used for the computations.
      CALL CPU_TIME (CPU1)

! Compute the CPU time and write it out in the .OUT file.
      CPUTIME_USED = CPU1 - CPU0 - CPU_IO
      WALLTIME_USED = WALL_TIME() - WALL0
      CALL WRITE_OUT3 (CPUTIME_USED, WALLTIME_USED, CPU_IO)

! JFD: cartesian grid implementation
      IF(WRITE_DASHBOARD) THEN
         IF(DT>=DT_MIN) THEN
            RUN_STATUS = 'Complete.'
         ELSE
            RUN_STATUS = 'DT < DT_MIN.  Recovery not possible!'
         ENDIF
         CALL GET_TUNIT(CPUTIME_USED,TUNIT)
         CALL UPDATE_DASHBOARD(0,CPUTIME_USED,TUNIT)
      ENDIF
      IF(CARTESIAN_GRID)  CALL CLOSE_CUT_CELL_FILES

! JFD: Dynamic load balance
      IF(TIME+0.1d0*DT>=TSTOP)  ADJUST_PARTITION=.FALSE.
      IF(ADJUST_PARTITION) THEN
         CALL OUTPUT_MANAGER(.TRUE., .FALSE.)
         CALL DEALLOCATE_ARRAYS
      ELSE
! Finalize and terminate MPI
         call parallel_fin
      ENDIF


      CALL FINL_ERR_MSG

      END SUBROUTINE FINALIZE

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: GEN_LOG_BASENAME                                        !
!  Author: Aytekin Gel                                Date: 19-SEP-03  !
!                                                                      !
!  Purpose: Generate the file base for DMP logs.                       !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE GEN_LOG_BASENAME

      use compar, only: myPE
      use compar, only: fbname

      implicit none

! Variables for generating file basename with processor id
      INTEGER :: i1, i10, i100, i1000, i10000

! PAR_I/O Generate file basename for LOG files
      i10000 = int(myPE/10000)
      i1000  = int((myPE-i10000*10000)/1000)
      i100   = int((myPE-i10000*10000-i1000*1000)/100)
      i10    = int((myPE-i10000*10000-i1000*1000-i100*100)/10)
      i1     = int((myPE-i10000*10000-i1000*1000-i100*100-i10*10)/1)

      i10000 = i10000 + 48
      i1000  = i1000  + 48
      i100   = i100   + 48
      i10    = i10    + 48
      i1     = i1     + 48

      fbname=char(i10000)//char(i1000)//char(i100)//char(i10)//char(i1)

      RETURN
      END SUBROUTINE GEN_LOG_BASENAME




   !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
   !                                                                      !
   !  SUBROUTINE: DO_MPI_BCAST                                            !
   !                                                                      !
   !  Purpose: Used by pymfix for broadcasting commands.                  !
   !                                                                      !
   !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   function do_mpi_bcast(str)
#ifdef MPI
      use mpi
#endif
      implicit none

      ! TODO: return a dynamically allocated string, instead of fixed size of 100,000 bytes
      character(len=100000),intent(in) :: str
      character :: aa(100000)
      character :: do_mpi_bcast(100000)
      integer :: ii
      integer :: ierr

      do ii = 1,len(str)
         aa(ii) = str(ii:ii)
      end do

#ifdef MPI
      call mpi_bcast(aa,100000,mpi_character,0,mpi_comm_world,ierr)
#endif

      do ii = 1,100000
         do_mpi_bcast(ii:ii) = aa(ii)
      end do

   end function do_mpi_bcast

   subroutine do_write_dbg_vtu_and_vtp_files
      implicit none
      call write_dbg_vtu_and_vtp_files
   end subroutine do_write_dbg_vtu_and_vtp_files

   subroutine do_backupres
      use output_man, only: backup_res
      implicit none
      call backup_res
   end subroutine do_backupres

   subroutine do_reinit(filename)
      use reinit, only: reinitialize
      implicit none
      ! filename of uploaded mfix.dat file
      character(len=*), intent(in) :: filename
      call reinitialize(filename)
   end subroutine do_reinit

   subroutine do_abort
      use compar, only: mype
      implicit none
      call mfix_exit(mype)
   end subroutine do_abort

   !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
   !  Subroutine: ADD_COMMAND_LINE_KEYWORD                                !
   !  Author: M.Meredith                                 Date: 03-FEB-16  !
   !                                                                      !
   !  Purpose: Save command line arguments in CMD_LINE_ARGS array.        !
   !           Used by both mfix.f and pymfix.                            !
   !                                                                      !
   !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   SUBROUTINE ADD_COMMAND_LINE_KEYWORD(ARG)
      implicit none
      CHARACTER(LEN=80), INTENT(IN) :: ARG

      CMD_LINE_ARGS_COUNT = CMD_LINE_ARGS_COUNT + 1

      IF (CMD_LINE_ARGS_COUNT > SIZE(CMD_LINE_ARGS)) THEN
         print *, "Too many command line arguments: ", COMMAND_ARGUMENT_COUNT()
         print *, "Only ",SIZE( CMD_LINE_ARGS )," or fewer command line arguments are supported."
         STOP
      ENDIF

      CMD_LINE_ARGS(CMD_LINE_ARGS_COUNT) = arg

   END SUBROUTINE ADD_COMMAND_LINE_KEYWORD

   !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
   !  Subroutine: PRINT_FLAGS                                             !
   !  Author: M.Meredith                                 Date: 27-APR-16  !
   !                                                                      !
   !  Purpose: Print the configure flags MFIX was built with.             !
   !                                                                      !
   !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   SUBROUTINE PRINT_FLAGS
      implicit none

#ifdef MPI
      write(*,"(A)",advance="no") "dmp "
#endif

#ifdef MKL
      write(*,"(A)",advance="no") "mkl "
#endif

#ifdef NETCDF
      write(*,"(A)",advance="no") "netcdf "
#endif

#ifdef PYMFIX
      write(*,"(A)",advance="no") "python "
#endif

#ifdef _OPENMP
      write(*,"(A)",advance="no") "smp "
#endif

      write(*,"(A)",advance="yes") ""

   END SUBROUTINE PRINT_FLAGS

END MODULE MAIN
