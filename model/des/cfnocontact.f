!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFNOCONTACT(L)                                         C
!  Purpose: DES - Zeroing values when particles are not in contact     C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CFNOCONTACT(L)
      
      USE discretelement
      IMPLICIT NONE

      INTEGER K, L    
!-----------------------------------------------------------------------------

      
      DO K = 1, DIMN
         FN(K,L) = 0.0
         FNS1(K) = 0.0
         FT(K,L) = 0.0
         FTS1(K) = 0.0
         FC(K,L) = 0.0
         TOW(K,L) = 0.0
         TOW(3,L) = 0.0
         OMEGA_NEW(K,L) = 0.0
         OMEGA_NEW(3,L) = 0.0
      END DO

      RETURN
      END SUBROUTINE CFNOCONTACT


