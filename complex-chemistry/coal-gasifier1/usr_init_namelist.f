!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: USR_INIT_NAMELIST                                      C
!  Purpose: initialize user_defined NAMELIST variables                 C
!                                                                      C
!  Author: S. Venkatesan, M. Syamlal                  Date:6-21-93     C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: None                                          C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE USR_INIT_NAMELIST
!
      Use param
      Use param1
      Use usr
      IMPLICIT NONE
      INCLUDE 'usrnlst.inc'
!
!
!     logical DEVOL set to FALSE as default
      USE_DEVOL       = .FALSE.
! INITIALIZE THE SOLID PHASE DEVOLOATILIZATION SECTION
!
      COAL = UNDEFINED_I
      PAFC = UNDEFINED
      PAVM = UNDEFINED
      PAA  = UNDEFINED
      PAM  = UNDEFINED
      UAC  = UNDEFINED
      UAH  = UNDEFINED
      UAO  = UNDEFINED
      UAN  = UNDEFINED
      UAS  = UNDEFINED
      HHVC = UNDEFINED
      HHVT = UNDEFINED
      RETURN
      END SUBROUTINE USR_INIT_NAMELIST
