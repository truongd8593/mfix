!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: STIFF_CHEM_DEBUG                                       !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!  Author: J.Musser                                   Date:            !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      MODULE STIFF_CHEM_DBG

      PRIVATE

! Variable Access:
!---------------------------------------------------------------------//
      PUBLIC :: ODE_DEBUG_LEVEL

! Subroutine Access:
!---------------------------------------------------------------------//
      PUBLIC :: ODE_UPDATE_OLD,   &
                ODE_RESET,        &
                INIT_ODE_STATS0,  &
                INIT_ODE_STATS,   &
                UPDATE_ODE_STATS, &
                WRITE_ODE_STATS,  &
                CHECK_ODE_DATA,   &
                WRITE_ODE_LOG


! Routine used to compare to values.
      LOGICAL, external :: COMPARE


! Static variables/parameters.
!---------------------------------------------------------------------//
! Debug Level: (Messages)
! 0 - None
! 1 - Limited
! 2 - Aggressive
      INTEGER :: ODE_DEBUG_LEVEL = 2

! Frequency to report the number of steps distribution.
      INTEGER, parameter :: reportNST_Freq = 10

! File unit for ODE Error Log.
      INTEGER, parameter :: OEL_Unit = 6589


! Variables updated once each call to the stiff solver.
!---------------------------------------------------------------------//
! Frequency to report the number of steps distribution.
      INTEGER :: reportNST
      INTEGER :: failedCount_total

! Variables updated every IJK loop cycle.
!---------------------------------------------------------------------//

! Original field variables. Used to reset after failed integration.
      DOUBLE PRECISION, allocatable :: ODE_VARS_o(:)

! The minimum number of integrations needed (over all IJK)
      INTEGER, allocatable :: minNST(:)                     ! local
      INTEGER, allocatable :: minNST_all(:)                 ! global

! The maximum number of integrations needed (over all IJK)
      INTEGER, allocatable :: maxNST(:)                     ! local
      INTEGER, allocatable :: maxNST_all(:)                 ! global

! An array that stores the distrubtion of the number of steps needed
! to integrate ODES.
      INTEGER, allocatable :: countNST(:)                   ! local
      INTEGER, allocatable :: countNST_all(:)               ! global

! Number of cells that only have homogeneous chemical reactions.
      INTEGER, allocatable :: lHomogeneous(:)               ! local
      INTEGER, allocatable :: lHomogeneous_all(:)           ! global

! Number of cells that only have homogeneous and/or heterogeneous
! chemical reactions.
      INTEGER, allocatable :: lHeterogeneous(:)             ! local
      INTEGER, allocatable :: lHeterogeneous_all(:)         ! global

! Number of cells that failed to successfully integration ODEs.
      INTEGER, allocatable :: failedCount(:)                ! local
      INTEGER, allocatable :: failedCount_all(:)            ! global

! Maximum number of attempts to integrate.
      INTEGER, allocatable :: maxAttempts(:)                ! local
      INTEGER, allocatable :: maxAttempts_all(:)            ! global



      DOUBLE PRECISION :: ODE_StartTime
      LOGICAL :: PRINT_ERR_HEADER

      contains


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Module name: ODE_UPDATE_OLD                                         !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!  Author: J.Musser                                   Date:            !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE ODE_UPDATE_OLD(lODE_DIMN, lODE_VARS)

      implicit none

! The number ODEs (maximum).
      INTEGER, intent(in) :: lODE_DIMN

! MFIX field variables mapped into the ODE format.
      DOUBLE PRECISION, intent(in) :: lODE_VARS(lODE_DIMN)

      ODE_VARS_o = lODE_VARS

      RETURN
      END SUBROUTINE ODE_UPDATE_OLD



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Module name: ODE_RESET                                              !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!  Author: J.Musser                                   Date:            !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE ODE_RESET(lDIMN, lVARS, lAtps)

      use compar,     only : myPE

      implicit none

! The number ODEs (maximum).
      INTEGER, intent(in) :: lDIMN

! MFIX field variables mapped into the ODE format.
      DOUBLE PRECISION, intent(inout) :: lVARS(lDIMN)

      INTEGER, intent(in) :: lAtps

      lVARS = ODE_VARS_o

      IF(lAtps == 3) failedCount(myPE) = failedCount(myPE) + 1

      RETURN
      END SUBROUTINE ODE_RESET

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: INIT_ODE_STATS0                                        !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!  Author: J.Musser                                   Date:            !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE INIT_ODE_STATS0(lODE_DIMN)

      use compar,     only : numPEs

      implicit none

! The number ODEs (maximum).
      INTEGER, intent(in) :: lODE_DIMN

      reportNST = 1
      failedCount_total = 0

      allocate( ODE_VARS_o(lODE_DIMN) )


! The minimum number of integrations needed (over all IJK)
      allocate( minNST(0:numPEs-1))      ! local
      allocate( minNST_all(0:numPEs-1) )  ! global

! The maximum number of integrations needed (over all IJK)
      allocate( maxNST(0:numPEs-1) )     ! local
      allocate( maxNST_all(0:numPEs-1) ) ! global

! An array that stores the distrubtion of the number of steps needed
! to integrate ODES.
      allocate( countNST(7) )     ! local
      allocate( countNST_all(7) ) ! global

! Number of cells that only have homogeneous chemical reactions.
      allocate( lHomogeneous(0:numPEs-1) )     ! local
      allocate( lHomogeneous_all(0:numPEs-1) ) ! global

! Number of cells that only have homogeneous and/or heterogeneous
! chemical reactions.
      allocate( lHeterogeneous(0:numPEs-1) )     ! local
      allocate( lHeterogeneous_all(0:numPEs-1) ) ! global

! Number of cells that failed to successfully integration ODEs.
      allocate( failedCount(0:numPEs-1) )     ! local
      allocate( failedCount_all(0:numPEs-1) ) ! global

! Maximum number of attempts to integrate.
      allocate( maxAttempts(0:numPEs-1) )     ! local
      allocate( maxAttempts_all(0:numPEs-1) ) ! global


      RETURN
      END SUBROUTINE INIT_ODE_STATS0


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: INIT_ODE_STATS                                         !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!  Author: J.Musser                                   Date:            !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE INIT_ODE_STATS

      use compar,   only : myPE, PE_IO

      implicit none

      CHARACTER*64 :: lFile

! use the subroutine from machine.f
      external CPU_TIME   


      if(myPE == PE_IO) &
         write(*,"(/3x,'Entering the stiff chemistry solver...',$)")

      PRINT_ERR_HEADER = .TRUE.

      maxNST = 0
      minNST = 0
      minNST(myPE) = 5000
      maxAttempts = 0

      IF(reportNST == 1) then
         countNST = 0
         countNST_all = 0
      ENDIF

      failedCount = 0
      lHomogeneous = 0
      lHeterogeneous = 0


      CALL CPU_TIME(ODE_StartTime)


      IF(ODE_DEBUG_LEVEL >= 1) THEN
         lFile = ''; write(lFile,"('StiffChem_',I2.2,'.log')")myPE
         open(OEL_Unit, file=lFile, status='unknown', position='append')
      ENDIF


      END SUBROUTINE INIT_ODE_STATS



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: CALC_ODE_STATS                                         !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!  Author: J.Musser                                   Date:            !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE UPDATE_ODE_STATS(lNEQ, lNEQ_DIMN, lNST, lODE_DIMN, &
         lAtps)

      use compar,   only : myPE, PE_IO

      implicit none

! The number of steps needed to integrate.
      INTEGER, intent(in) :: lNEQ_DIMN

! (1) :: Number of ODEs
! (2) :: Fluid cell index (IJK) passed into ODEPACK
      INTEGER, dimension(lNEQ_DIMN), intent(in) :: lNEQ
! The number of steps needed to integrate.
      INTEGER, intent(in) :: lNST
! The number ODEs (maximum).
      INTEGER, intent(in) :: lODE_DIMN

! The number of attempts.
      INTEGER, intent(in) :: lAtps


      IF(lNEQ(1) == lODE_DIMN) THEN
         lHeterogeneous(myPE) = lHeterogeneous(myPE) + 1
      ELSE
         lHomogeneous(myPE) = lHomogeneous(myPE) + 1
      ENDIF

      maxAttempts(myPE) = max(lAtps, maxAttempts(myPE))

      minNST(myPE) = min(minNST(myPE), lNST)
      maxNST(myPE) = max(maxNST(myPE), lNST)

      IF (lNST <           10) THEN
         countNST(1) = countNST(1) + 1 
      ELSE IF (lNST <     100) THEN
         countNST(2) = countNST(2) + 1 
      ELSE IF (lNST <    1000) THEN
         countNST(3) = countNST(3) + 1 
      ELSE IF (lNST <   10000) THEN
         countNST(4) = countNST(4) + 1 
      ELSE IF (lNST <  100000) THEN
         countNST(5) = countNST(5) + 1 
      ELSE IF (lNST < 1000000) THEN
         countNST(6) = countNST(6) + 1 
      ELSE
         countNST(7) = countNST(7) + 1 
      ENDIF


      RETURN
      END SUBROUTINE UPDATE_ODE_STATS



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Module name: WRITE_ODE_STATS                                        !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!  Author: J.Musser                                   Date:            !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE WRITE_ODE_STATS

      use compar,   only : myPE, PE_IO

      use mpi_utility

      implicit none

! Message buffer.
      CHARACTER*64 :: lMsg0, lMsg1

      DOUBLE PRECISION :: lODE_EndTime, lODE_RunTime

! use the subroutine from machine.f
      external CPU_TIME   

! Update screen message.
      IF(myPE == PE_IO) WRITE(*,"(2x,'DONE.',/)")


      CALL CPU_TIME(lODE_EndTime)
      lODE_RunTime = lODE_EndTime - ODE_StartTime


! Collect stats on min/max number of steps.
      minNST_all = 0; CALL global_sum(minNST, minNST_all)
      maxNST_all = 0; CALL global_sum(maxNST, maxNST_all)

! Collect stats on the number of cells with pure homogeneous reactions.
      lHomogeneous_all = 0;
      CALL global_sum(lHomogeneous, lHomogeneous_all)

! Collect stats on the number of cells with heterogeneous and 
! homogeneous reactions.
      lHeterogeneous_all = 0
      CALL global_sum(lHeterogeneous, lHeterogeneous_all)

! Collect stats on the maximum number of integration attempts.
      maxAttempts_all = 0
      CALL global_sum(maxAttempts, maxAttempts_all)

! Collect stats on the number of failed integrations.
      failedCount_all = 0
      CALL global_sum(failedCount, failedCount_all)


! Display stiff solver summary.
      IF(myPE == PE_IO) THEN

! Report Min/Max steps:
         lMsg0=''; write(lMsg0,*) minval(minNST_all)
         lMsg1=''; write(lMsg1,*) maxval(maxNST_all)
         write(*,1000)  trim(adjustl(lMsg0)), trim(adjustl(lMsg1))

! Report Homogeneous/Heterogeneous reactions:
         lMsg0=''; write(lMsg0,*) sum(lHomogeneous_all)
         lMsg1=''; write(lMsg1,*) sum(lHeterogeneous_all)
         write(*,1001) trim(adjustl(lMsg0)), trim(adjustl(lMsg1))

! Report Max attempts:
         lMsg0=''; write(lMsg0,*) maxval(maxAttempts_all)
         write(*,1004)  trim(adjustl(lMsg0))


! Report failed integrations:
         failedCount_total = failedCount_total + sum(failedCount_all)

         lMsg0=''; write(lMsg0,*) sum(failedCount_all)
         lMsg1=''; write(lMsg1,*) failedCount_total
         write(*,1002) trim(adjustl(lMsg0)), trim(adjustl(lMsg1))


         IF(lODE_RunTime > 3.6d3) THEN
            lMsg0=''; write(lMsg0,"(f8.4)") lODE_RunTime/3.6d3
            lMsg1='hrs'
         ELSEIF(lODE_RunTime > 6.0d1) THEN
            lMsg0=''; write(lMsg0,"(f8.4)") lODE_RunTime/6.0d1
            lMsg1='min'
         ELSE
            lMsg0=''; write(lMsg0,"(f8.4)") lODE_RunTime
            lMsg1='sec'
         ENDIF
         write(*,1003) trim(adjustl(lMsg0)), trim(adjustl(lMsg1))

      ENDIF


      if(reportNST == reportNST_Freq) then
! Collect the number of steps distrubutions.
         countNST_all = 0;
         CALL global_sum(countNST, countNST_all)

         countNST_all = int(countNST_all/reportNST_Freq)

         if(myPE == PE_IO) then
            write(*,"(/5x,'Average Integration Distrubution:')")
            write(*,"(7x,'NST <      10: ', I6)")countNST_all(1)
            write(*,"(7x,'NST <     100: ', I6)")countNST_all(2)
            write(*,"(7x,'NST <    1000: ', I6)")countNST_all(3)
            write(*,"(7x,'NST <   10000: ', I6)")countNST_all(4)
            write(*,"(7x,'NST <  100000: ', I6)")countNST_all(5)
            write(*,"(7x,'NST < 1000000: ', I6)")countNST_all(6)
            write(*,"(7x,'NST > 1000000: ', I6)")countNST_all(7)
         endif
! Reset the reporting counter.
         reportNST = 1
      else
! Increment the reporting counter.
         reportNST = reportNST + 1
      endif

      if(myPE == PE_IO)write(*,"(/' ')")


      close(OEL_Unit)

      RETURN

 1000 Format(5x,'Minimum/Maximum number of steps over all cells: ',A,'/',A)
 1001 Format(5x,'Number of cells with Homogeneous/Heterogeneous reactions: ',A,'/',A)
 1002 Format(5x,'Number of Current/Cumulative failed integrations: ',A,'/',A)
 1003 Format(5x,'CPU Time Used: ',A,' ',A)
 1004 Format(5x,'Maximum number of integration attempts: ',A)

      END SUBROUTINE WRITE_ODE_STATS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: CHECK_ODE_DATA                                         !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!  Author: J.Musser                                   Date:            !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_ODE_DATA(lnD, lNEQ, loD, VARS, lState, lErr)

      use param1,   only : ZERO, SMALL_NUMBER, ONE
      use physprop, only : NMAX, MMAX
      use toleranc, only : TMax, TMin

      implicit none

      INTEGER, intent(in) :: lnD

      INTEGER, intent(in) :: lNEQ(lnD)

! The number ODEs (maximum).
      INTEGER, intent(in) :: loD

! MFIX field variables mapped into the ODE format.
      DOUBLE PRECISION, intent(in) :: VARS(loD)
! Sate value returned from ODEPACK
      INTEGER, intent(in)  :: lState
! Error Flag
      INTEGER, intent(out) :: lErr


! Loop indicies:
      INTEGER :: M  ! phase
      INTEGER :: N  ! species
      INTEGER :: Node  ! other


      lErr = 0

! Check ODE_VAR for NaNs.
!---------------------------------------------------------------------//
      DO Node=1,loD
         IF(isNaN(VARS(node))) THEN
            lErr = -8
            return
         ENDIF
      ENDDO


! Check ODEPACK State
!---------------------------------------------------------------------//
      IF(lState == -1) THEN
         lErr = lState
         return
      ELSEIF(lState == 2) THEN
         lErr = 0
      ELSE
         lErr = lState
      ENDIF


! Initialize.
      Node = 1

! Check variables for physical values.
!---------------------------------------------------------------------//

! Gas phase density.
      IF(VARS(Node) <= ZERO) THEN
         lErr = -(100 + Node)
         return
      ENDIF
      Node = Node + 1


! Gas phase temperature.
      IF(VARS(Node) > Tmax .OR. VARS(Node) < Tmin) THEN
         lErr = -(100 + Node)
         return
      ENDIF
      Node = Node + 1


! Gas phase species mass fractions.
      DO N=1,NMAX(0)
         IF(VARS(Node) > 1.009d0) THEN
            lErr = -(100 + Node)
            return
         ENDIF
         IF(VARS(Node) < 0.0d0) THEN
            IF(abs(VARS(Node)) > 1.0d-5) THEN
               lErr = -(100 + Node)
               return
            ENDIF
         ENDIF
         Node = Node + 1
      ENDDO


! Solids temperature.
      DO M = 1, MMAX
         IF(VARS(Node) > Tmax .OR. VARS(Node) < Tmin) THEN
            lErr = -(100 + Node)
            return
         ENDIF
         Node = Node + 1
      ENDDO



      DO M = 1, MMAX
         IF(lNEQ(2+M) == 1) THEN

! Solids volume fraction.
            IF(VARS(Node) <= ZERO) THEN
               lErr = -(100 + Node)
               return
            ENDIF
            Node = Node + 1

! Solids phase species mass fractions.
            DO N=1,NMAX(M)

               IF(VARS(Node) > 1.009d0) THEN
                  lErr = -(100 + Node)
                  return
               ENDIF
               IF(VARS(Node) < 0.0d0) THEN
                  IF(abs(VARS(Node)) > 1.0d-5) THEN
                     lErr = -(100 + Node)
                     return
                  ENDIF
               ENDIF
               Node = Node + 1

            ENDDO
         ENDIF
      ENDDO   


      
      RETURN
      END SUBROUTINE CHECK_ODE_DATA


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: ODE_ErrorLog                                           !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!  Author: J.Musser                                   Date:            !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE WRITE_ODE_LOG(lErr, lnD, lNEQ, loD, lVARS)

      use compar,   only : myPE
      use param1,   only : SMALL_NUMBER
      use run,      only : TIME
      use rxns,     only : NO_OF_RXNS
      use physprop, only : MMAX

      use indices

      implicit none

      INTEGER, intent(in) :: lErr

      INTEGER, intent(in) :: lnD
      INTEGER, intent(in) :: lNEQ(lnD)

! The number ODEs (maximum).
      INTEGER, intent(in) :: loD
      DOUBLE PRECISION, intent(in) :: lVARS(loD)

      INTEGER :: lc1
      INTEGER :: IJK

      DOUBLE PRECISION :: ddt_lVARS(loD)

      IJK = lNEQ(2)

      if(PRINT_ERR_HEADER) then
         WRITE(OEL_Unit,9000) Time
         PRINT_ERR_HEADER = .FALSE.
      endif


      WRITE(OEL_Unit,9001) IJK, I_OF(IJK), J_OF(IJK), myPE
      WRITE(OEL_Unit,9002) lErr
      WRITE(OEL_Unit,9003) lNEQ(1)

      DO lc1=1,MMAX
         IF(lNEQ(lc1+2) == 1) THEN
            WRITE(OEL_Unit,9004) lc1
         ELSEIF(lNEQ(lc1+2) == 0) THEN
            WRITE(OEL_Unit,9005) lc1
         ELSE
            WRITE(OEL_Unit,9006) lc1, lNEQ(lc1+2)
         ENDIF
      ENDDO

      WRITE(OEL_Unit,"(/5x,'Field Variables:')")
      WRITE(OEL_Unit,9007)
      DO lc1=1, loD
         WRITE(OEL_Unit,9008) lc1, ODE_VARS_o(lc1), lVARS(lc1), &
            (lVARS(lc1) - ODE_VARS_o(lc1))
      ENDDO

      CALL STIFF_CHEM_RRATES(lNEQ, Time, lVARS, ddt_lVARS)

      WRITE(OEL_Unit,"(/5x,'Derivatives:')")
      WRITE(OEL_Unit,9009)
      DO lc1=1, loD
         WRITE(OEL_Unit,9010) lc1, ddt_lVARS(lc1)
      ENDDO


      RETURN

 9000 FORMAT(//3x,'Time: ',g15.8)
 9001 FORMAT(//5x,'Fluid Cell: ',I4,4x,'(',I3,',',I3,')',10x,&
         'Process: ',I4)
 9002 FORMAT(  5x,'Error Flag/ODEPACK Status: ',I4,/)

 9003 FORMAT(  5x,'Number of ODEs Solved: ', I4,/)
 9004 FORMAT(  5x,'Solving solids phase ',I2,&
         ' bulk density and species equations.')

 9005 FORMAT(  5x,'NOT solving solids phase ',I2,&
         ' bulk density and species equations.')

 9006 FORMAT(  5x,'Unknown flag for solids phase ',I2,'  Flag = ',I4)


 9007 FORMAT(  5x,'Var',6x,'Incoming',13x,'Returned',13x,'Difference')
 9008 FORMAT(  5x,I3,3(3x,G18.8))


 9009 FORMAT(  5x,'Var',6x,'  YDOT  ')
 9010 FORMAT(  5x,I3,3x,G18.8)

      END SUBROUTINE WRITE_ODE_LOG


      END MODULE STIFF_CHEM_DBG
