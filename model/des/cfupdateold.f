!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFUPDATEOLD(PARTS)                                     C
!  Purpose: DES - Update arrays                                        C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CFUPDATEOLD(PARTS)
      
      USE discretelement
      IMPLICIT NONE

      INTEGER LL, K, PARTS
!     
!---------------------------------------------------------------------

      DO LL = 1, PARTICLES
         DO K = 1, DIMN
            DES_POS_OLD(K,LL) = DES_POS_NEW(K,LL)
            DES_VEL_OLD(K,LL) = DES_VEL_NEW(K,LL)
            OMEGA_OLD(K,LL) = OMEGA_NEW(K,LL)
         END DO
         OMEGA_OLD(3,LL) = OMEGA_NEW(3,LL)
      END DO

      RETURN
      END SUBROUTINE CFUPDATEOLD


