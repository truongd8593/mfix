! -*- f90 -*-
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: MFIX                                                    !
!  Author: M. Syamlal                                 Date: 29-JAN-92  !
!                                                                      !
!  Purpose: The main module in the MFIX program                        !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!
!> \mainpage Multiphase Flow with Interphase eXchanges
!!
!! MFIX is a general-purpose computer code developed at the National
!! Energy Technology Laboratory, NETL, for describing the hydrodynamics,
!! heat transfer, and chemical reactions in fluid-solid systems.
!!
!! It has been used for describing bubbling and circulating fluidized
!! beds and spouted beds. MFiX calculations give transient data on the
!! three-dimensional distribution of pressure, velocity, temperature,
!! and species mass fractions. MFiX code is based on a generally
!! accepted set of multiphase flow equations. The code is used as a
!! "test-stand" for testing and developing multiphase flow constitutive
!!  equations.
!!
!! \section Notice
!! Neither the United States Government nor any agency thereof, nor any
!! of their employees, makes any warranty, expressed or implied, or
!! assumes any legal liability or responsibility for the accuracy,
!! completeness, or usefulness of any information, apparatus, product,
!! or process disclosed or represents that its use would not infringe
!! privately owned rights.
!!
!! * MFIX is provided without any user support for applications in the
!!   user's immediate organization. It should not be redistributed in
!!   whole or in part.
!!
!! * The use of MFIX is to be acknowledged in any published paper based
!!   on computations using this software by citing the MFIX theory
!!   manual. Some of the submodels are being developed by researchers
!!   outside of NETL. The use of such submodels is to be acknowledged
!!   by citing the appropriate papers of the developers of the submodels.
!!
!! * The authors would appreciate receiving any reports of bugs or other
!!   difficulties with the software, enhancements to the software, and
!!   accounts of practical applications of this software.
!!
!! \section Disclaimer
!! This report was prepared as an account of work sponsored by an agency
!! of the United States Government. Neither the United States Government
!! nor any agency thereof, nor any of their employees, makes any
!! warranty, express or implied, or assumes any legal liability or
!! responsibility for the accuracy, completeness, or usefulness of any
!! information, apparatus, product, or process disclosed, or represents
!! that its use would not infringe privately owned rights. Reference
!! herein to any specific commercial product, process, or service by
!! trade name, trademark, manufacturer, or otherwise does not
!! necessarily constitute or imply its endorsement, recommendation, or
!! favoring by the United States Government or any agency thereof. The
!! views and opinions of authors expressed herein do not necessarily
!! state or reflect those of the United States Government or any
!! agency thereof.

      PROGRAM MFIX

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE COMPAR, only:ADJUST_PARTITION
      USE DES_TIME_MARCH, ONLY: DES_TIME_INIT, DES_TIME_STEP, DES_TIME_END, FACTOR, EXIT_LOOP
      USE DISCRETELEMENT, ONLY: DES_CONTINUUM_COUPLED, DISCRETE_ELEMENT
      USE ERROR_MANAGER, ONLY: INIT_ERROR_MANAGER
      USE ITERATE, ONLY: CONVERGED, DIVERGED, ADJUSTDT
      USE ITERATE, ONLY: ITERATE_INIT, DO_ITERATION, POST_ITERATE
      USE ITERATE, ONLY: LOG_CONVERGED, LOG_DIVERGED, NIT, MAX_NIT
      USE MAIN, ONLY: ADD_COMMAND_LINE_KEYWORD
      USE MAIN, ONLY: GET_DATA, CHECK_DATA, PRE_INIT, INITIALIZE, FINALIZE
      USE MAIN, ONLY: EXIT_SIGNAL, MFIX_DAT, PRINT_FLAGS
      USE MPI_UTILITY
      USE RUN, ONLY:  DT, DEM_SOLIDS, PIC_SOLIDS, STEADY_STATE, TIME, TSTOP, RUN_NAME, RUN_TYPE
      USE STEP, ONLY: CHECK_LOW_DT, CHEM_MASS, TIME_STEP_INIT, TIME_STEP_END
      USE param1, only: n_spx
      USE parallel_mpi, only: parallel_init
      IMPLICIT NONE

!-----------------------------------------------
! Local variables
!-----------------------------------------------

      INTEGER :: II
      CHARACTER(LEN=80) :: tmp
      LOGICAL :: read_dat

      read_dat = .false.
      DO II=1, COMMAND_ARGUMENT_COUNT()

         IF (read_dat) THEN
            CALL GET_COMMAND_ARGUMENT(II,mfix_dat)
            read_dat = .false.
            CYCLE
         ELSE
            CALL GET_COMMAND_ARGUMENT(II,tmp)
         ENDIF

         IF (tmp=='-h'.or.tmp=='--help') THEN
            print *, "Usage: mfix [-h,--help] [-p,--print-flags] &
               &[-f,--file <filename>] [<KEYWORD>=<VALUE> ...]"
            print *, "       -h,--help: display this help message"
            print *, "       -p,--print-flags: print flags MFIX was &
               &built with (if any): dmp mkl netcdf python smp"
            print *, "       -f,--file <filename>: &
               &specify filename of input file (Default: mfix.dat)"
            print *, "       <KEYWORD>=<VALUE>: specify keyword on &
               &command line, overrides values in mfix.dat"
            STOP
         ELSE IF (tmp=='-p'.or.tmp=='--print-flags') THEN
            CALL PRINT_FLAGS
            STOP
         ELSE IF (tmp=='-f'.or.tmp=='--file') THEN
            read_dat = .true.
            CYCLE
         ELSE
            CALL ADD_COMMAND_LINE_KEYWORD(tmp)
         ENDIF
      ENDDO

! Invoke MPI/SMP initialization and get rank info.
      CALL PARALLEL_INIT

! Dynamic load balance loop
      DO

! Read input file and prepare to do check input
      CALL PRE_INIT

! Check data, do computations for IC and BC locations
! and flows, and set geometry parameters such as X, X_E, DToDX, etc.
      CALL CHECK_DATA

! Initialize the simulation
      CALL INITIALIZE

! Time march loop.

      IF(DISCRETE_ELEMENT .AND. .NOT.DES_CONTINUUM_COUPLED) THEN

! Uncoupled discrete element simulations
         IF (DEM_SOLIDS) THEN
            CALL DES_TIME_INIT
            DO II = 1, FACTOR
               CALL DES_TIME_STEP(II)
               IF ( EXIT_LOOP ) EXIT
            ENDDO
            CALL DES_TIME_END
         ENDIF
         IF (PIC_SOLIDS) CALL PIC_TIME_MARCH

      ELSE

! Transient or steady state simulation
         DO WHILE (TIME + 0.1d0*DT < TSTOP .AND. .NOT. EXIT_SIGNAL)
            CALL TIME_STEP_INIT
            DO
               CALL ITERATE_INIT
               DO WHILE (NIT<MAX_NIT .AND. .NOT.(CONVERGED.OR.DIVERGED))
                  NIT = NIT + 1
                  CALL DO_ITERATION
               ENDDO

               CALL POST_ITERATE

               IF(STEADY_STATE) EXIT
               IF(.NOT.ADJUSTDT()) EXIT
            ENDDO

! Exit if DT < DT_MIN
            CALL CHECK_LOW_DT

! Stiff Chemistry Solver.
            CALL CHEM_MASS

! DEM time loop
            IF (DEM_SOLIDS) THEN
               CALL DES_TIME_INIT
               DO II = 1, FACTOR
                  CALL DES_TIME_STEP(II)
                  IF ( EXIT_LOOP ) EXIT
               ENDDO
               CALL DES_TIME_END
            ENDIF

            CALL TIME_STEP_END
            IF (STEADY_STATE .OR. ADJUST_PARTITION) EXIT
         ENDDO

      ENDIF


      CALL FINALIZE
      IF(.NOT.ADJUST_PARTITION) EXIT
      ENDDO


      STOP
      END PROGRAM MFIX
