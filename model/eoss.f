!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: EOSS0                                                  !
!                                                                      !
!  Author: J.Musser                                   Date: 02-Dec-13  !
!  Reviewer:                                                           !
!                                                                      !
!  Purpose: Calculate the initial solids density. This calculation is  !
!  only valid at time zero. Thus, this routine should only be invoked  !
!  by the initialization routines.                                     !
!                                                                      !
!  Literature/Document References:                                     !
!                                                                      !
!  Local variables: None                                               !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      DOUBLE PRECISION FUNCTION EOSS0(M)

! Global Variables:
!---------------------------------------------------------------------//
! Baseline/initial solids density
      use physprop, only: RO_Xs0
! Baseline/initial solids mass fractions.
      use physprop, only: X_s0
! Number of species comprising each phase.
      use physprop, only: NMAX
! Process Rank
      use compar, only: myPE
! Unit number for RUN_NAME.LOG file
      use funits, only: UNIT_LOG
! Logical for who writes error messages
      use funits, only: DMP_LOG

! Global parameters
!---------------------------------------------------------------------/
      use param1, only: ONE
      use param1, only: ZERO
      use param1, only: UNDEFINED
      use param1, only: SMALL_NUMBER

      implicit none

! Passed Arguments:
!---------------------------------------------------------------------/
! Solids phase index.
      INTEGER, intent(in) :: M

! Local Variables:
!---------------------------------------------------------------------/
! Alias for inert species index.
      DOUBLE PRECISION :: OoRO_s0
! Character string for error messages
      CHARACTER(len=64) :: MSG

! Evaluate the first part of the calculation.
      OoRO_s0 = sum(X_S0(M,:NMAX(M))/RO_Xs0(M,:NMAX(M)))
! If the value is physical (positive) finish the calculation and return.
      IF(OoRO_s0 > ZERO) THEN
         EOSS0 = ONE/OoRO_s0
         return
      ENDIF

! This is an extra sanity check that should be caught be one ore more
! of the data checks.
      MSG=''
      IF(abs(OoRO_S0) <= SMALL_NUMBER) THEN
         WRITE(MSG,"('Infinity')")
      ELSE
         WRITE(MSG,*) ONE/OoRO_s0
      ENDIF

      IF(DMP_LOG) THEN
         WRITE(*,1000) M, trim(MSG)
         WRITE(UNIT_LOG,1000) M, trim(MSG)
      ENDIF

      CALL MFIX_EXIT(myPE)

 1000 FORMAT(//1X,70('*')/' From: EOSS',/,' Error 1300:'               &
         ' Unphysical baseline density calculated:',/' RO_s(',I2,') = '&
         ,A,/' Please refer to the Readme file on the required input', &
         ' and make',/' the necessary corrections to the data file.',  &
         /1X,70('*')//)

      END FUNCTION EOSS0



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: EOSS                                                   !
!                                                                      !
!  Author: J.Musser                                   Date: 09-Oct-13  !
!  Reviewer:                                                           !
!                                                                      !
!  Purpose: Calculate solid density - runtime.                         !
!                                                                      !
!  Literature/Document References:                                     !
!                                                                      !
!  Local variables: None                                               !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      DOUBLE PRECISION FUNCTION EOSS(pBase, Xs0_INERT, Xs_INERT)

      implicit none

! Passed Arguments:
!---------------------------------------------------------------------/
! Baseline phase density (unreacted)
      DOUBLE PRECISION, intent(in) :: pBase
! Baseline inert mass fraction
      DOUBLE PRECISION, intent(in) :: Xs0_INERT
! Current mass fraction of inert
      DOUBLE PRECISION, intent(in) :: Xs_INERT

! Evaluate the solids EOS.
      EOSS = pBase * Xs0_INERT / max(Xs_INERT, 1.0d-8)

      RETURN  
      END FUNCTION EOSS
