!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: MFIX                                                   C
!  Purpose: The main module in the MFIX program                        C
!                                                                      C
!  Author: M. Syamlal                                 Date: 29-JAN-92  C
!  Reviewer: P. Nicoletti, W. Rogers, S. Venkatesan   Date: 29-JAN-92  C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: Add version                                                C
!  Author: P.Nicoletti                                Date: 07-FEB-92  C
!  Revision Number: 2                                                  C
!  Purpose: Add a routine to set velocity components normal to a wall  C
!           to zero                                                    C
!  Author: M. Syamlal                                 Date: 14-MAY-92  C
!  Reviewer: W. Rogers                                Date: 11-DEC-92  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:  Modifications for initializing & finalizing MPI through   C
!            parallel_mpi module developed by ORNL and also starts     C
!            the MPI Timer (MPI_WTIME)                                 C
!
!  Author:   Aeolus Res. Inc.                         Date: 01-AUG-99  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: RUN_TYPE                                      C
!  Variables modified: CPU0, ID_VERSION                                C
!                                                                      C
!  Local variables:       CPU1, CPUTIME_USED,L                         C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      PROGRAM MFIX 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE run
      USE time_cpu 
      USE funits 
      USE output
      USE compar      !// 001 Include MPI header file
      USE mpi_utility !//     added for exitMPI calls for debugging     
      USE parallel_mpi
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------

!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!
!                       Final value of CPU time.
       DOUBLE PRECISION CPU1
!
!                       CPU time used for the computations.
       DOUBLE PRECISION CPUTIME_USED
!
!                       Save TIME in input file for RESTART_2
       DOUBLE PRECISION TIME_SAVE
!
!                       Temporary storage for DT
       DOUBLE PRECISION DT_tmp
!
!                       loop counter
      INTEGER           L
!
!                      Error index
      INTEGER          IER
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
!-----------------------------------------------
!
!$      INTEGER num_threads, threads_specified
!$      INTEGER mp_numthreads, omp_get_num_threads
      INTEGER IJK
      INCLUDE 'function.inc'

!// 010 Initialize MPI & get ranks & total PEs employed
    call parallel_init
!//PARDBG 0801 Printout total PEs and also signal each alive PE
!   if(myPE == 0) then
!       write(*,"('Total number of PEs  :',I3)") numPEs
!   endif
!   write(*,"('Hello, I am PE :',I2)") myPE



!
! set the version.release of the software
!
      ID_VERSION = '(V03.10)' 
!
!   Specify the number of processors to be used
!
!$      WRITE(*,'(A,$)') 'Enter the number of processors to be used: '
!$      READ(*,*) threads_specified
!$      call omp_set_num_threads(threads_specified)
!
!       Find the number of processors used
!
!$      num_threads = omp_get_num_threads()
!$      Write(*,*)' Number of processors used = ',  num_threads
!
!  Set machine dependent constants
!
      CALL MACHINE_CONS 
!
!  Get the date and time.
!      They give the unique run_id in binary output files
!
      CALL GET_RUN_ID 
!
!  Get the initial value of CPU time
!

!// 888 Start MPI's own clock to measure wall time
!      CALL CPU_TIME (CPU0) 
      CPU0 = MPI_WTIME()

!
!   Read input data, check data, do computations for IC and BC locations
!   and flows, and set geometry parameters such as X, X_E, DToDX, etc.
!
      CALL GET_DATA 

!//AIKEPARDBG
!      write(*,"('(PE ',I2,'): aft get_data in mfix')") myPE    !//AIKEPARDBG
!      call exitMPI(myPE)   !//AIKEPARDBGSTOP

!
!  Initialize all field variables as undefined
!
      CALL INIT_FVARS 

!
!  Set the flags for identifying computational cells
!
      CALL SET_FLAGS 

!//AIKEPARDBG
!      write(*,"('(PE ',I2,'): aft set_flags in mfix')") myPE    !//AIKEPARDBG
!      call exitMPI(myPE)   !//AIKEPARDBGSTOP

!
!  Set constant physical properties
!
      CALL SET_CONSTPROP 
      
!//AIKEPARDBG
!      write(*,"('(PE ',I2,'): aft set_constprop in mfix')") myPE    !//AIKEPARDBG
!      call exitMPI(myPE)   !//AIKEPARDBGSTOP

!
!
!  Write the initial part of the standard output file
!
      CALL WRITE_OUT0 
!      call MPI_Barrier(MPI_COMM_WORLD,mpierr)
!//AIKEPARDBG
!      write(*,"('(PE ',I2,'): after write_out0 in mfix')") myPE    !//AIKEPARDBG
!      call exitMPI(myPE)   !//AIKEPARDBGSTOP

!
!  Write the initial part of the special output file(s)
!
      CALL WRITE_USR0 

!//AIKEPARDBG
!      write(*,"('(PE ',I2,'): after write_usr0 in mfix')") myPE    !//AIKEPARDBG
!      call exitMPI(myPE)   !//AIKEPARDBGSTOP
      
!$
!$    CALL START_LOG 
!$    WRITE (UNIT_LOG, *) ' '  
!$    WRITE (UNIT_LOG, *) ' Number of processors used = ', threads_specified  
!$    WRITE (UNIT_LOG, *) ' '  
!$    CALL END_LOG 

!
!  setup for PC quickwin application
!
      CALL PC_QUICKWIN 
!
      DT_TMP = DT 
      SELECT CASE (RUN_TYPE)  
      CASE ('NEW')  
!
!  Write the initial part of the restart files
!
         CALL WRITE_RES0 
!      call MPI_Barrier(MPI_COMM_WORLD,mpierr)
!      write(*,"('(PE ',I2,'): after write_res0 in mfix')") myPE    !//AIKEPARDBG
!      call exitMPI(myPE)   !//AIKEPARDBGSTOP

         DO L = 1, N_SPX 
            CALL WRITE_SPX0 (L) 
         END DO 
      CASE ('RESTART_1')  
!
! Read the time-dependant part of the restart file
!
         CALL READ_RES1 
!
         CALL START_LOG 
         WRITE (UNIT_LOG, 1010) TIME, NSTEP 
         CALL END_LOG 
         WRITE (UNIT_OUT, 1010) TIME, NSTEP 
         IF (FULL_LOG) WRITE (*, 1010) TIME, NSTEP 
!
      CASE ('RESTART_2')  
!
         TIME_SAVE = TIME 
         CALL READ_RES1 
         TIME = TIME_SAVE 
!
         CALL START_LOG 
         WRITE (UNIT_LOG, 1010) TIME, NSTEP 
         CALL END_LOG 
         WRITE (UNIT_OUT, 1010) TIME, NSTEP 
         IF (FULL_LOG) WRITE (*, 1010) TIME, NSTEP 
!
         CALL WRITE_RES0 
         CALL WRITE_RES1 
         DO L = 1, N_SPX 
            CALL WRITE_SPX0 (L) 
            CALL WRITE_SPX1 (L) 
         END DO 
      CASE DEFAULT 
!
         CALL START_LOG 
         WRITE (UNIT_LOG, *) ' MFIX: Do not know how to process' 
         WRITE (UNIT_LOG, *) ' RUN_TYPE in data file' 
         CALL END_LOG 
         STOP  
!
      END SELECT 
!
!//AIKEPARDBG
       call MPI_Barrier(MPI_COMM_WORLD,mpierr)   !//AIKEPARDBG
!      write(*,"('(PE ',I2,'): after write_XXXX in mfix')") myPE    !//AIKEPARDBG
!      call exitMPI(myPE)   !//AIKEPARDBGSTOP

      IF (DT_TMP /= UNDEFINED) THEN 
         DT = MAX(DT_MIN,MIN(DT_MAX,DT)) 
!
      ELSE 
         DT = DT_TMP 
!
      ENDIF 
!
!  Set arrays for computing indices
!
      CALL SET_INCREMENTS 
!
!      write(*,"('(PE ',I2,'): after set_increments in mfix')") myPE    !//AIKEPARDBG
!
!  Set the flags for wall surfaces impermeable and identify flow boundaries
!  using FLAG_E, FLAG_N, and FLAG_T
!
      CALL SET_FLAGS1 

!      write(*,"('(PE ',I2,'): after set_flags1 in mfix')") myPE    !//AIKEPARDBG
!      call exitMPI(myPE)   !//AIKEPARDBGSTOP

!
!  Calculate cell volumes and face areas
!
      CALL SET_GEOMETRY1 

!      write(*,"('(PE ',I2,'): after set_geometry1 in mfix')") myPE    !//AIKEPARDBG
!      call exitMPI(myPE)   !//AIKEPARDBGSTOP

!
!  Find corner cells and set their face areas to zero
!
      CALL GET_CORNER_CELLS (IER) 

!      write(*,"('(PE ',I2,'): after get_corner_cells in mfix')") myPE    !//AIKEPARDBG
!
!  Set initial conditions
!
      CALL SET_IC 

!      write(*,"('(PE ',I2,'): after set_ic in mfix')") myPE    !//AIKEPARDBG
!      call exitMPI(myPE)   !//AIKEPARDBGSTOP

!
!  Set boundary conditions
!
      CALL ZERO_NORM_VEL 
      CALL SET_BC0 

!      write(*,"('(PE ',I2,'): after set_bc0 in mfix')") myPE    !//AIKEPARDBG
!      call exitMPI(myPE)   !//AIKEPARDBGSTOP

!
!  Set gas mixture molecular weight
!
      CALL SET_MW_MIX_G 

!      write(*,"('(PE ',I2,'): after set_mw_mix_g in mfix')") myPE    !//AIKEPARDBG
!
!  Set the pressure field for a fluidized bed
!
      IF (RUN_TYPE == 'NEW') CALL SET_FLUIDBED_P 

!      write(*,"('(PE ',I2,'): after set_fluidbed_p in mfix')") myPE    !//AIKEPARDBG
!
!  Initialize gas densities
!
      CALL SET_RO_G 

!      write(*,"('(PE ',I2,'): after set_ro_g in mfix')") myPE    !//AIKEPARDBG
!
!  Initialize time dependent boundary conditions
!
      CALL SET_BC1 

!      write(*,"('(PE ',I2,'): after set_bc1 in mfix')") myPE    !//AIKEPARDBG
!      call exitMPI(myPE)   !//AIKEPARDBGSTOP

!
!  Check the field variable data and report errors.
!
      CALL CHECK_DATA_20 
!
!  Initializations for CPU time calculations in iterate
!
      CPUOS = 0. 
      CALL CPU_TIME (CPU1) 
      CPU_NLOG = CPU1 
      TIME_NLOG = TIME - DT 

!     write(*,"('(PE ',I2,'): reached beginning of time march in mfix')") myPE    !//AIKEPARDBG
!     call exitMPI(myPE)   !//AIKEPARDBGSTOP

!
!  Find the solution of the equations from TIME to TSTOP at
!  intervals of DT
!
      CALL TIME_MARCH 
!
!  Get the final value of CPU time.  The difference gives the
!  CPU time used for the computations.
!

!// 888 Read time from MPI's own clock to measure wall time
!      CALL CPU_TIME (CPU1) 
      CPU1 = MPI_WTIME()
!
!  Compute the CPU time and write it out in the .OUT file.
!
      CPUTIME_USED = CPU1 - CPU0 
      CALL WRITE_OUT3 (CPUTIME_USED) 

!// 010 Finalize and terminate MPI
      call parallel_fin
            
      STOP  
 1000 FORMAT(/1X,'MFIX ',A,' Simulation:'/) 
 1010 FORMAT(/1X,70('*')//' From: MFIX',/&
         ' Message: Read in data from .RES file for TIME = ',G12.5,/&
         ' Time step number (NSTEP) =',I7,/1X,70('*')/) 
      END PROGRAM MFIX 
      
