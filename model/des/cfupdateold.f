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
            DES_POS_OLD(LL,:) = DES_POS_NEW(LL,:)
            DES_VEL_OLD(LL,:) = DES_VEL_NEW(LL,:)
            OMEGA_OLD(LL,:) = OMEGA_NEW(LL,:)
      END DO

      RETURN
      END SUBROUTINE CFUPDATEOLD


