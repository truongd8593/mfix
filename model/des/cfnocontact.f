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

      INTEGER K, L    
!-----------------------------------------------------------------------------

      
         FN(L,:) = ZERO
         FNS1(:) = ZERO
         FT(L,:) = ZERO
         FTS1(:) = ZERO
         FC(L,:) = ZERO
         TOW(L,:) = ZERO
         OMEGA_NEW(L,:) = ZERO

      RETURN
      END SUBROUTINE CFNOCONTACT


