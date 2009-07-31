!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFNOCONTACT(L)                                         C
!>
!!  Purpose: DES - Zeroing values when particles are not in contact     
!<
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CFNOCONTACT(L)

      USE param1      
      USE discretelement
      IMPLICIT NONE

!-----------------------------------------------
! Local Variables
!-----------------------------------------------
      INTEGER L    
!-----------------------------------------------

      
         FN(L,:) = ZERO
         FT(L,:) = ZERO
         FC(L,:) = ZERO
         TOW(L,:) = ZERO
         OMEGA_NEW(L,:) = ZERO

      RETURN
      END SUBROUTINE CFNOCONTACT


