!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFFT(L, Vtan, OVERLP_T, TANGNT                         C
!  Purpose: DES - Calculate tangential force on a particle             C
!           due to inter-particle collision                            C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CFFT(L, Vtan, OVERLP_T, TANGNT)
      
      USE discretelement
      IMPLICIT NONE
      
      INTEGER L, K
      DOUBLE PRECISION OVERLP_T, TANGNT(NDIM), Vtan
!---------------------------------------------------------------------


      DO K = 1, DIMN
         FTS1(K) = 0 - KT*OVERLP_T*TANGNT(K)
         FT(K,L) = FTS1(K) - ETA_DES_T*Vtan*TANGNT(K)
      END DO

!     PRINT *, FT(1,L), FT(2,L)

      RETURN
      END SUBROUTINE CFFT


