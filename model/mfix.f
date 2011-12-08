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
!  Revision Number: 3                                                  C
!  Purpose:  Modifications for initializing & finalizing MPI through   C
!            parallel_mpi module developed by ORNL and also starts     C
!            the MPI Timer (MPI_WTIME)                                 C
!                                                                      C
!  Author:   Aeolus Res. Inc.                         Date: 01-AUG-99  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number: 4                                                  C
!  Purpose: Modifications to call Cartesian grid subroutines           C
!           and dashboard                                              C
!                                                                      C
!  Author:   Jeff Dietiker                            Date: 01-JUL-09  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
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
      USE compar      
      USE mpi_utility     
      USE parallel_mpi
      USE discretelement
      USE mfix_pic

!DISTIO      
      USE cdist      
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
      USE fldvar
      USE cutcell
      USE quadric
      USE dashboard
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================


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
! JFD
!                       CPU time unit.
      CHARACTER(LEN=4)  :: TUNIT
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
!DISTIO
      character :: version*512      
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      external cpu_time
!-----------------------------------------------
!
!$      INTEGER num_threads, threads_specified, omp_id
!$      INTEGER mp_numthreads, omp_get_num_threads
!$      INTEGER omp_get_thread_num      
      INTEGER IJK
      INCLUDE 'function.inc'

!DISTIO
      ! note: the value below for version must be the same as
      !       the value of version in mfix.f
      !       If you change it in this subroutine, you must
      !       change it in write_res0.f also
      !   
      !       The value should also be consistent with the check
      !       in read_res0
      
	version = 'RES = 01.6'

	bDoing_postmfix = .false.
	
!// 010 Initialize MPI & get ranks & total PEs employed
    call parallel_init
    
    ! we want only PE_IO to write out common error messages
    if(enable_dmp_log)then
      dmp_log = .true.
    elseif(mype == pe_io) then
      dmp_log = .true.
    else
      dmp_log = .false.
    endif
!
! set the version.release of the software
!
      ID_VERSION = '2010-1'

!
! set automatic restart flag to false
!      AUTOMATIC_RESTART = .FALSE.
!      ITER_RESTART      = 1
!
!!   Specify the number of processors to be used
!
!$      WRITE(*,'(A,$)') 'Enter the number of threads to be used for SMP: '
!$      READ(*,*) threads_specified
!$      call omp_set_num_threads(threads_specified)
!
!!       Find the number of processors used
!
!$omp  parallel
!$      num_threads = omp_get_num_threads()
!$      omp_id = omp_get_thread_num()
!$      if(omp_id.eq.0) Write(*,*)' Number of threads used for SMP = ',  num_threads
!$omp  end parallel
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
!AEOLUS STOP Trigger mechanism to terminate MFIX normally before batch queue terminates
!            timestamp at the beginning of execution
      CALL CPU_TIME (CPU00)
!
!   Read input data, check data, do computations for IC and BC locations
!   and flows, and set geometry parameters such as X, X_E, DToDX, etc.
!
      CALL GET_DATA 



!
!  Initialize all field variables as undefined
!
      CALL INIT_FVARS 

!
!  Set the flags for identifying computational cells
!
      CALL SET_FLAGS 

!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
      CALL CHECK_DATA_CARTESIAN
      IF(CARTESIAN_GRID) THEN
         CALL CUT_CELL_PREPROCESSING
      ELSE
         CALL ALLOCATE_DUMMY_CUT_CELL_ARRAYS
      ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
!
!  Set constant physical properties
!
      CALL SET_CONSTPROP 
      
!
!
!  Write the initial part of the standard output file
!
      CALL WRITE_OUT0
      CALL WRITE_FLAGS 

!
!  Write the initial part of the special output file(s)
!
      CALL WRITE_USR0 
      
!$
!$    CALL START_LOG 
!$    IF(DMP_LOG)WRITE (UNIT_LOG, *) ' '  
!$    IF(DMP_LOG)WRITE (UNIT_LOG, *) ' Number of processors used = ', threads_specified  
!$    IF(DMP_LOG)WRITE (UNIT_LOG, *) ' '  
!$    CALL END_LOG 

!
!  setup for PC quickwin application
!
      CALL PC_QUICKWIN 
!
  101 CONTINUE
      IF(AUTOMATIC_RESTART) THEN
         RUN_TYPE = 'RESTART_1'
         AUTOMATIC_RESTART = .FALSE.
         ITER_RESTART = ITER_RESTART + 1
         CALL CHECK_DATA_06
         CALL CHECK_DATA_07
         CALL CHECK_DATA_08
         CALL CHECK_DATA_09
         CALL SET_FLAGS
         CALL SET_CONSTPROP
      ENDIF
!
      DT_TMP = DT 
      SELECT CASE (TRIM(RUN_TYPE))  
      CASE ('NEW')  
!
!  Write the initial part of the restart files
!
         CALL WRITE_RES0 

         DO L = 1, N_SPX 
            CALL WRITE_SPX0 (L, 0) 
         END DO 
      CASE ('RESTART_1')  
!
! Read the time-dependant part of the restart file
!
         CALL READ_RES1 
!
         CALL START_LOG 
         IF(DMP_LOG)WRITE (UNIT_LOG, 1010) TIME, NSTEP 
         CALL END_LOG 
         if (myPE == PE_IO) WRITE (UNIT_OUT, 1010) TIME, NSTEP 
         IF (FULL_LOG) WRITE (*, 1010) TIME, NSTEP 
!
      CASE ('RESTART_2')  
!
         TIME_SAVE = TIME 
!DISTIO
         if (myPE .ne. PE_IO .and. bDist_IO .and. bStart_with_one_res) then
             write (unit_res,rec=1) version
             write (unit_res,rec=2) 4
             write (unit_res,rec=3) 4
         end if	 
	 
         CALL READ_RES1 
         TIME = TIME_SAVE 
!
         CALL START_LOG 
         IF(DMP_LOG)WRITE (UNIT_LOG, 1010) TIME, NSTEP 
         CALL END_LOG 
         if (myPE == PE_IO) WRITE (UNIT_OUT, 1010) TIME, NSTEP 
         IF (FULL_LOG) WRITE (*, 1010) TIME, NSTEP 
!
         CALL WRITE_RES0 
         CALL WRITE_RES1 
         DO L = 1, N_SPX 
            CALL WRITE_SPX0 (L, 0) 
            CALL WRITE_SPX1 (L, 0) 
         END DO 
      CASE DEFAULT 
!
         CALL START_LOG 
         IF(DMP_LOG)WRITE (UNIT_LOG, *) ' MFIX: Do not know how to process' 
         IF(DMP_LOG)WRITE (UNIT_LOG, *) ' RUN_TYPE in data file' 
         CALL END_LOG 
         call mfix_exit(myPE)  
!
      END SELECT 
!
       call MPI_Barrier(MPI_COMM_WORLD,mpierr)   !//AIKEPARDBG

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
!
!  Set arrays for computing indices for higher order scheme
!
      CALL SET_INCREMENTS3 
!
!  Set the flags for wall surfaces impermeable and identify flow boundaries
!  using FLAG_E, FLAG_N, and FLAG_T
!
      CALL SET_FLAGS1 

!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
!
!  Calculate cell volumes and face areas
!
      IF(.NOT.CARTESIAN_GRID)  THEN
         CALL SET_GEOMETRY1 
       ELSE
         CALL SET_GEOMETRY
       ENDIF
!
!  Find corner cells and set their face areas to zero
!
      IF(.NOT.CARTESIAN_GRID)  THEN
         CALL GET_CORNER_CELLS (IER) 
      ELSE
         IF (SET_CORNER_CELLS)  CALL GET_CORNER_CELLS (IER)
      ENDIF

!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
!
!  Set initial conditions
!
      CALL SET_IC 

!
!  Set boundary conditions
!
      CALL ZERO_NORM_VEL 
      CALL SET_BC0 

!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
      IF(CARTESIAN_GRID) CALL CG_SET_BC0 
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
!
!  Set gas mixture molecular weight
!
      CALL SET_MW_MIX_G 


!
!  Set the pressure field for a fluidized bed
!
      IF (RUN_TYPE == 'NEW') CALL SET_FLUIDBED_P 

!
!  Initialize gas densities
!
      CALL SET_RO_G 

!
!  Initialize time dependent boundary conditions
!
      CALL SET_BC1 

!
!  Check the field variable data and report errors.
!
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
      IF(.NOT.CARTESIAN_GRID)      CALL CHECK_DATA_20
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================



! Set constants and allocate/initialize DEM variables
      IF(DISCRETE_ELEMENT) THEN
! moving all of DES related initializations from get_data.f to here.
! This is because it is best to initialize DES once all the fluid and
! geometry information has been obtained and initial fields set. 
! RG 2/15/2011         
         CALL CHECK_DES_DATA
         CALL CHECK_DES_RXNS
         CALL CHECK_DES_THERMO
         CALL CHECK_DES_IC
         CALL CHECK_DES_BC
         CALL MAKE_ARRAYS_DES
         !STOP
      ELSE
! If discrete_element is .false. then overwrite the following user DES
! logicals which may be set to true in the input file 
         DES_CONTINUUM_COUPLED = .FALSE.
         DES_INTERP_ON = .FALSE.
         TSUJI_DRAG = .FALSE.
         PRINT_DES_DATA = .FALSE.
         MPPIC = .FALSE. 
         DES_ONEWAY_COUPLED = .false. 
      ENDIF

!DISTIO
! for creating files needed by post_mfix with distributed IO

!AEOLUS DEBUG PRINT 
    if (DBGPRN_LAYOUT .or. bdist_io) then
!     write (*,*) myPE , ' E.4 ... version = ' , version(1:33)
      call debug_write_layout(1,ier)
      call write_parallel_info(1,ier)
    endif
!
!  Initializations for CPU time calculations in iterate
!
      CPUOS = 0. 
      CALL CPU_TIME (CPU1) 
      CPU_NLOG = CPU1 
      TIME_NLOG = TIME - DT 

!  Get the initial value of CPU time
!
      CALL CPU_TIME (CPU0) 
!
!  Find the solution of the equations from TIME to TSTOP at
!  intervals of DT
!
      CALL TIME_MARCH
      IF(AUTO_RESTART.AND.AUTOMATIC_RESTART.AND.ITER_RESTART.LE.10) GOTO 101

!  Call user-defined subroutine after time-loop.
      IF (CALL_USR) CALL USR3
!
!  Get the final value of CPU time.  The difference gives the
!  CPU time used for the computations.
!
      CALL CPU_TIME (CPU1) 
!
!  Compute the CPU time and write it out in the .OUT file.
!
      CPUTIME_USED = CPU1 - CPU0 - CPU_IO
      if(myPE.eq.root) then
      WRITE(*,*) '************** CPU TIME for IO **********************',CPU_IO
      endif
      CALL WRITE_OUT3 (CPUTIME_USED) 

!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
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
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================

!// Finalize and terminate MPI
      call parallel_fin
            
      STOP  
 1000 FORMAT(/1X,'MFIX ',A,' Simulation:'/) 
 1010 FORMAT(/1X,70('*')//' From: MFIX',/&
         ' Message: Read in data from .RES file for TIME = ',G12.5,/&
         ' Time step number (NSTEP) =',I7,/1X,70('*')/) 
      END PROGRAM MFIX 
      
!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: debug_write(debuglvl, ier )                            C
!  Purpose: Write out full geometry index setup information for the case C
!                                                                      C
!  Author: Aytekin Gel                                Date: 19-SEP-03  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE debug_write_layout(debuglvl, IER)

!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE matrix
      USE geometry
      USE compar
      USE mpi_utility
      USE sendrecv
      USE sendrecv3      
      USE indices
      USE leqsol
!DISTIO      
      USE cdist      
      USE funits
      USE run       ! added for AMG input parameters entered in mfix.dat
!AE AMG 091503
      USE time_cpu
!      USE hypre
      IMPLICIT NONE
      
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Error indicator
      INTEGER ::          IER
!                      phase index
      INTEGER ::          M
!                      debug level
      INTEGER ::          debuglvl
      
      INTEGER ::      i,j,k,ijk,ijk_GL,ijk_PROC,ijk_IO
!
      integer :: i_of_g,j_of_g,k_of_g
      integer :: indxA,indxA_gl,indxB,indxB_gl,indxC,indxC_gl
      integer :: indxD,indxD_gl,indxE,indxE_gl,indxF,indxF_gl
      integer :: indxG,indxG_gl,indxH,indxH_gl      
!            
      logical :: amgdbg = .TRUE.

      
      character fname*80
!
      INCLUDE 'function.inc'
      
  
       k_of_g(ijk) = int( (ijk-1)/( (imax3-imin3+1)*(jmax3-jmin3+1) ) ) + kmin3
       i_of_g(ijk) = int( ( (ijk-  (k_of_g(ijk)-kmin3)*((imax3-imin3+1)*(jmax3-jmin3+1))) &
                     - 1)/(jmax3-jmin3+1)) + imin3
       j_of_g(ijk) = ijk - (i_of_g(ijk)-imin3)*(jmax3-jmin3+1) - &
                     (k_of_g(ijk)-kmin3)*((imax3-imin3+1)*(jmax3-jmin3+1)) - 1 + jmin3  
                      

!DISTIO
!      fname = "layout_xxxx.txt"
!      write (fname(8:11),'(i4.4)') myPE
      fname = "layout_xxxxx.txt"
      write (fname(8:12),'(i5.5)') myPE      
      open (unit=11,file=fname,status='unknown')
      
      write (11,*) ' ********************************************'
      write (11,*) ' ********************************************'
      write (11,*) ' ********************************************'
      write (11,*) ' ********************************************'
      write (11,*) ' '
      write (11,*) ' '
      write (11,*) ' myPE =           ' , myPE
      write (11,*) ' '
      write (11,*) ' '
 

      if (AMGDBG .or. bDist_IO) then
      write(11,"('BLK1: Running from istart3,iend3 .AND. jstart3, jend3 .AND. kstart3, kend3')")
      write(11,"(' (   i ,    j,     k) =>    ijk      ijk_GL     ijk_PROC    ijk_IO')") 
      write(11,"(' ====================      =====     =======    ========    ======')")
      DO k = kstart3, kend3 
        DO i = istart3,iend3
          DO j = jstart3, jend3
             ijk = FUNIJK(i,j,k)
             ijk_GL = FUNIJK_GL(i,j,k)
             ijk_PROC = FUNIJK_PROC(i,j,k,myPE)
             ijk_IO = FUNIJK_IO(i,j,k)            
             write(11,"(' (',I4,' , ',I4,' , ',I4,') => ',4(I8,' , '))") &
                                         i,j,k,ijk,ijk_GL,ijk_PROC,ijk_IO
          END DO
        END DO
      END DO

      write(11,"(/,/,'BLK2: Print out Bottom, South, West, East, North, Top neighbors')")
      write(11,"(' (   i ,    j,     k) =>    ijk    ijk_GL    B_of    S_of    W_of    E_of    N_of    T_of')") 
      write(11,"(' ====================      =====   =======  ======  ======  ======  ======  ======  ======')")
      DO k = kstart3, kend3 
        DO i = istart3,iend3
          DO j = jstart3, jend3
             ijk = FUNIJK(i,j,k)
             ijk_GL = FUNIJK_GL(i,j,k)
             write(11,"(' (',I4,' , ',I4,' , ',I4,') => ',2(I7,' , '),6(I7,2X))") &
                                         i,j,k,ijk,ijk_GL,bottom_of(ijk),south_of(ijk),west_of(ijk),&
                                         east_of(ijk),north_of(ijk),top_of(ijk)
          END DO
        END DO
      END DO

      write(11,"(/,/,'BLK3: Print out km, jm, im, ip, jp, kp neighbors')")      
      write(11,"(' (   i ,    j,     k) =>    ijk    ijk_GL    km_of   jm_of   im_of   ip_of   jp_of   kp_of')") 
      write(11,"(' ====================      =====   =======  ======  ======  ======  ======  ======  ======')")
      DO k = kstart3, kend3 
        DO i = istart3,iend3
          DO j = jstart3, jend3
             ijk = FUNIJK(i,j,k)
             ijk_GL = FUNIJK_GL(i,j,k)
             write(11,"(' (',I4,' , ',I4,' , ',I4,') => ',2(I7,' , '),6(I7,2X))") &
                                         i,j,k,ijk,ijk_GL,km_of(ijk),jm_of(ijk),im_of(ijk),&
                                         ip_of(ijk),jp_of(ijk),kp_of(ijk)
          END DO
        END DO
      END DO

      write(11,"(/,'BLK4a: Active Fluid Cells:FLUID_AT(ijk)=.T.',/,&
     &           ' (   i ,    j,     k) =>    ijk  [   x ,     ,     z]')") 
      write(11,"(' ====================      =====  ====================')")
       DO ijk = ijkstart3, ijkend3
         I = I_OF(IJK) 
         J = J_OF(IJK) 
         K = K_OF(IJK)        
        
!         IF (FLOW_AT_E(IJK)) THEN 
         IF (FLUID_AT(IJK)) THEN         
!          write(11,"(' (',I4,' , ',I4,' , ',I4,') => ',I8)") I,J,K,ijk
           write(11,"(' (',I4,' , ',I4,' , ',I4,') => ',I8,' [',E12.5,',',E12.5,' ]')") I,J,K,ijk,X(i),Z(k)
         ENDIF
      END DO

      write(11,"(/,'BLK4b: Cells that are (.NOT.WALL_AT(IJK)) = .T.',/,&
     &           ' (   i ,    j,     k) =>    ijk  [   x ,     ,     z]')")
      write(11,"(' ====================      =====  ====================')")
       DO ijk = ijkstart3, ijkend3
         I = I_OF(IJK)
         J = J_OF(IJK)
         K = K_OF(IJK)

         IF (.NOT.WALL_AT(IJK)) THEN
!          write(11,"(' (',I4,' , ',I4,' , ',I4,') => ',I8)") I,J,K,ijk
           write(11,"(' (',I4,' , ',I4,' , ',I4,') => ',I8,' [',E12.5,',',E12.5,' ]')") I,J,K,ijk,X(i),Z(k)
         ENDIF
      END DO

      DO k = kstart3, kend3 
        DO i = istart3,iend3
          DO j = jstart3, jend3
             ijk = FUNIJK(i,j,k)
             ijk_GL = FUNIJK_GL(i,j,k)

             if (i == istart2 .AND. j == jstart2) then 
                 indxA = ijk 
                 indxA_gl = ijk_GL
             endif       
             if (i == istart1 .AND. j == jstart1) then 
                 indxE = ijk
                 indxE_gl = ijk_GL
             endif                       
             if (i == istart2 .AND. j == jend2) then 
                 indxB = ijk
                 indxB_gl = ijk_GL
             endif                       
             if (i == istart1 .AND. j == jend1) then 
                 indxF = ijk
                 indxF_gl = ijk_GL
             endif                       
             if (i == iend1 .AND. j == jstart1) then 
                 indxH = ijk
                 indxH_gl = ijk_GL
             endif                       
             if (i == iend2 .AND. j == jstart2) then 
                 indxD = ijk
                 indxD_gl = ijk_GL
             endif                       
             if (i == iend1 .AND. j == jend1) then 
                 indxG = ijk
                 indxG_gl = ijk_GL
             endif                       
             if (i == iend2 .AND. j == jend2) then 
                 indxC = ijk
                 indxC_gl = ijk_GL
             endif                       

          END DO
        END DO
        write(11,"('BLK5:')")
        write(11,"(57('='))")
        write(11,"('k= ',I5,/,57('='))") k
        write(11,"('B= ',I5,' (',I7,')',20X,'C= ',I5,' (',I7,')',/)") indxB, indxB_gl, &
                        indxC, indxC_gl
!        write(UNIT_LOG,"(' \',34X,'/')")
!        write(UNIT_LOG,"(2X,'\',32X,'/')")
        write(11,"(3X,'F= ',I5,' (',I7,')',12X,'G= ',I5,' (',I7,')')") indxF, indxF_gl, &
                        indxG, indxG_gl
        write(11,"(4(9X,'|',29X,'|',/))")
        write(11,"(3X,'E= ',I5,' (',I7,')',12X,'H= ',I5,' (',I7,')',/)") indxE, indxE_gl, &
                        indxH, indxH_gl
!        write(UNIT_LOG,"(2X,'/',32X,'\')")
!        write(UNIT_LOG,"('/',34X,'\')")
        write(11,"('A= ',I5,' (',I7,')',20X,'D= ',I5,' (',I7,')',/,/)") indxA, indxA_gl, &
                        indxD, indxD_gl

!        write(UNIT_LOG,"(' (',I4,' , ',I4,' , ',I4,') => ',2(I7,' , '),6(I7,2X))") &
!                                         i,j,k,ijk,ijk_GL,bottom_of(ijk),south_of(ijk),west_of(ijk),&
!                                        east_of(ijk),north_of(ijk),top_of(ijk)

      END DO
  
!      write(UNIT_LOG,"(/,' (   i ,    j,     k) =>    ijk (Active Fluid)')") 
!      write(UNIT_LOG,"(' ====================      =====')")
!       DO ijk = ijkstart3, ijkend3
!         I = I_OF(IJK) 
!         J = J_OF(IJK) 
!         K = K_OF(IJK)               
        
!         IF (FLOW_AT_E(IJK)) THEN 
!         IF (FLUID_AT(IJK)) THEN        
!           write(UNIT_LOG,"(' (',I4,' , ',I4,' , ',I4,') => ',I8)") I,J,K,ijk
!         ENDIF
!      END DO


      endif     ! if (AMGDBG) branch
            
      M = 0 
!      CALL WRITE_AB_M (A_M, B_M, IJKMAX2, M, IER)      
      


      if (AMGDBG .or. bDist_IO) then
      write(11,"(/,/,'BLK6: ========= ORIGINAL MFIX VARIABLES ===========')")
      write(11,"('PE ',I5,': imin1  = ',I6,3X,'imax1= ',I6,/,'PE ',I5,': jmin1  = ',I6,3X,'jmax1= ',I6)") & 
             myPE,imin1,imax1,myPE,jmin1,jmax1
      write(11,"('PE ',I5,': kmin1  = ',I6,3X,'kmax1= ',I6)") myPE,kmin1,kmax1
      write(11,"('-----')")      
      write(11,"('PE ',I5,': imin2  = ',I6,3X,'imax2= ',I6,/,'PE ',I5,': jmin2  = ',I6,3X,'jmax2= ',I6)") & 
             myPE,imin2,imax2,myPE,jmin2,jmax2
      write(11,"('PE ',I5,': kmin2  = ',I6,3X,'kmax2= ',I6)") myPE,kmin2,kmax2       
      write(11,"('----- Below xxx3 set is DMP extension ------------')")
      write(11,"('PE ',I5,': imin3  = ',I6,3X,'imax3= ',I6,/,'PE ',I5,': jmin3  = ',I6,3X,'jmax3= ',I6)") & 
             myPE,imin3,imax3,myPE,jmin3,jmax3
      write(11,"('PE ',I5,': kmin3  = ',I6,3X,'kmax3= ',I6)") myPE,kmin3,kmax3
      write(11,"('----- End of Below xxx3 set is DMP extension -----')")
!      write(11,"('PE ',I5,': ijkmax2= ',I6)") myPE,ijkmax2
      write(11,"('PE ',I5,': ijmax2 = ',I6)") myPE,ijmax2
      write(11,"('PE ',I5,': ijkmin1= ',I6,' ijkmax1= ',I12)") myPE,ijkmin1, ijkmax1
      write(11,"('PE ',I5,':          ',6X,' ijkmax2= ',I12)") myPE,ijkmax2           
      write(11,"('PE ',I5,':          ',6X,' ijkmax3= ',I12)") myPE,ijkmax3      
      write(11,"('PE ',I5,': ijkmin4= ',I6,' ijkmax4= ',I12)") myPE,ijkmin4, ijkmax4      


      write(11,"(/,/,' ========= DMP EXTENSION VARIABLES ===========')")
!      write(UNIT_LOG,"('PE ',I5,': ijksize  = ',I6)") myPE,ijksize
      write(11,"('PE ',I5,': ijksize3 = ',I6,3X,'ijksize3_all = ',I6)") myPE,ijksize3,ijksize3_all(myPE)
      write(11,"('PE ',I5,': ijksize4 = ',I6,3X,'ijksize4_all = ',I6)") myPE,ijksize4,ijksize4_all(myPE)      
      write(11,"('PE ',I5,': ijkstart3  = ',I6,3X,'ijkend3  = ',I6)") myPE,ijkstart3, ijkend3
      write(11,"('PE ',I5,': ijkstart3_all = ',I6,3X,'ijkstart4_all = ',I6)") myPE,ijkstart3_all(myPE),ijkstart4_all(myPE)
      write(11,"('PE ',I5,': istart_all = ',I6,3X,'iend_all = ',I6,/,'PE ',I5,': jstart_all = ',I6,3X,'jend_all = ',I6)") & 
             myPE,istart_all(myPE),iend_all(myPE),myPE,jstart_all(myPE),jend_all(myPE)
      write(11,"('PE ',I5,': kstart_all = ',I6,3X,'kend_all = ',I6,/,'----------------------')") & 
             myPE,kstart_all(myPE),kend_all(myPE)

      write(11,"('PE ',I5,': istart1_all= ',I6,3X,'iend1_all= ',I6,/,'PE ',I5,': jstart1_all= ',I6,3X,'jend3_all= ',I6)") & 
             myPE,istart1_all(myPE),iend1_all(myPE),myPE,jstart1_all(myPE),jend1_all(myPE)
      write(11,"('PE ',I5,': kstart1_all= ',I6,3X,'kend1_all= ',I6,/,'----------------------')") & 
             myPE,kstart1_all(myPE),kend1_all(myPE)

      write(11,"('PE ',I5,': istart2_all= ',I6,3X,'iend2_all= ',I6,/,'PE ',I5,': jstart2_all= ',I6,3X,'jend3_all= ',I6)") & 
             myPE,istart2_all(myPE),iend2_all(myPE),myPE,jstart2_all(myPE),jend2_all(myPE)
      write(11,"('PE ',I5,': kstart2_all= ',I6,3X,'kend2_all= ',I6,/,'----------------------')") & 
             myPE,kstart2_all(myPE),kend2_all(myPE)

      write(11,"('PE ',I5,': istart3_all= ',I6,3X,'iend3_all= ',I6,/,'PE ',I5,': jstart3_all= ',I6,3X,'jend3_all= ',I6)") & 
             myPE,istart3_all(myPE),iend3_all(myPE),myPE,jstart3_all(myPE),jend3_all(myPE)
      write(11,"('PE ',I5,': kstart3_all= ',I6,3X,'kend3_all= ',I6,/,'----------------------')") & 
             myPE,kstart3_all(myPE),kend3_all(myPE)

      write(11,"('PE ',I5,': istart1= ',I6,3X,'iend1= ',I6,/,'PE ',I5,': jstart1= ',I6,3X,'jend1= ',I6)") & 
             myPE,istart1,iend1,myPE,jstart1,jend1
      write(11,"('PE ',I5,': kstart1= ',I6,3X,'kend1= ',I6,/,'----------------------')") & 
             myPE,kstart1,kend1
      write(11,"('PE ',I5,': istart2= ',I6,3X,'iend2= ',I6,/,'PE ',I5,': jstart2= ',I6,3X,'jend2= ',I6)") & 
             myPE,istart2,iend2,myPE,jstart2,jend2
      write(11,"('PE ',I5,': kstart2= ',I6,3X,'kend2= ',I6,/,'----------------------')") & 
             myPE,kstart2,kend2
      write(11,"('PE ',I5,': istart3= ',I6,3X,'iend3= ',I6,/,'PE ',I5,': jstart3= ',I6,3X,'jend3= ',I6)") & 
             myPE,istart3,iend3,myPE,jstart3,jend3
      write(11,"('PE ',I5,': kstart3= ',I6,3X,'kend3= ',I6,/,'----------------------')") & 
             myPE,kstart3,kend3
      endif     ! if (AMGDBG) branch
      
      close(unit=11)


      return
      
      
      END SUBROUTINE DEBUG_WRITE_LAYOUT


      SUBROUTINE write_parallel_info(debuglvl, IER)

!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE matrix
      USE geometry
      USE compar
      USE mpi_utility
      USE sendrecv
      USE sendrecv3      
      USE indices
      USE leqsol
      USE funits
      USE run       ! added for AMG input parameters entered in mfix.dat
!AE AMG 091503
      USE time_cpu
!      USE hypre
      IMPLICIT NONE
      
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Error indicator
      INTEGER ::          IER
!                      phase index
      INTEGER ::          M
!                      debug level
      INTEGER ::          debuglvl
      
      INTEGER ::      i,j,k,ijk,ijk_GL,ijk_PROC,ijk_IO
!
      integer :: i_of_g,j_of_g,k_of_g
      integer :: indxA,indxA_gl,indxB,indxB_gl,indxC,indxC_gl
      integer :: indxD,indxD_gl,indxE,indxE_gl,indxF,indxF_gl
      integer :: indxG,indxG_gl,indxH,indxH_gl      
!            
      logical :: amgdbg = .TRUE.
      
      character fname*80
!
      INCLUDE 'function.inc'
      
  
       k_of_g(ijk) = int( (ijk-1)/( (imax3-imin3+1)*(jmax3-jmin3+1) ) ) + kmin3
       i_of_g(ijk) = int( ( (ijk-  (k_of_g(ijk)-kmin3)*((imax3-imin3+1)*(jmax3-jmin3+1))) &
                     - 1)/(jmax3-jmin3+1)) + imin3
       j_of_g(ijk) = ijk - (i_of_g(ijk)-imin3)*(jmax3-jmin3+1) - &
                     (k_of_g(ijk)-kmin3)*((imax3-imin3+1)*(jmax3-jmin3+1)) - 1 + jmin3  
                      
!DISTIO
!      fname = "p_info_xxxx.txt"
!      write (fname(8:11),'(i4.4)') myPE
      fname = "p_info_xxxxx.txt"
      write (fname(8:12),'(i5.5)') myPE      
      open (unit=11,file=fname,status='unknown')
      
      write (11,*) myPe , ' = myPE'
 
  
       write (11,*) myPE , istart3,iend3
       write (11,*) myPE , jstart3,jend3
       write (11,*) myPE , kstart3,kend3

      
      write(11,"('BLK1: Running from istart3,iend3 .AND. jstart3, jend3 .AND. kstart3, kend3')")
      write(11,"(' (   i ,    j,     k)       ijk      ijk_GL     ijk_PROC    ijk_IO')") 
      write(11,"(' ====================      =====     =======    ========    ======')")
      DO k = kstart3, kend3 
        DO i = istart3,iend3
          DO j = jstart3, jend3
             ijk = FUNIJK(i,j,k)
             ijk_GL = FUNIJK_GL(i,j,k)
             ijk_PROC = FUNIJK_PROC(i,j,k,myPE)
             ijk_IO = FUNIJK_IO(i,j,k)            
             write(11,"('  ',I4,'   ',I4,'   ',I4,'     ',4(I8,'   '))" ) &
                                         i,j,k,ijk,ijk_GL,ijk_PROC,ijk_IO
          END DO
        END DO
      END DO

            
      M = 0 
!      CALL WRITE_AB_M (A_M, B_M, IJKMAX2, M, IER)      
      


      close(unit=11)


      return
      END SUBROUTINE write_parallel_info

