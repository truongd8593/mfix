!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFFCTOW(L, II, NORM)                                       C
!  Purpose: DES - Calculate the total force and torque on a particle   C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CFFCTOW(L, II,  NORM)
     
      USE param1 
      USE discretelement
      IMPLICIT NONE
      
      INTEGER L, II, K
      DOUBLE PRECISION NORM(DIMN), CHECK
!---------------------------------------------------------------------
!     

         FC(L,:) = FC(L,:) + FN(L,:) + FT(L,:) 
         FC(II,:) = FC(II,:) - FN(L,:) - FT(L,:)

      IF(DIMN.EQ.3) THEN
         TOW(L,1) = TOW(L,1) + DES_RADIUS(L)*(NORM(2)*FT(L,3) - NORM(3)*&
                    FT(L,2))
         TOW(L,2) = TOW(L,2) + DES_RADIUS(L)*(NORM(3)*FT(L,1) - NORM(1)*&
                    FT(L,3))
         TOW(L,3) = TOW(L,3) + DES_RADIUS(L)*(NORM(1)*FT(L,2) - NORM(2)*&
                    FT(L,1))

         TOW(II,1) = TOW(II,1) - DES_RADIUS(II)*(NORM(2)*FT(L,3) - NORM(3)*&
                    FT(L,2))
         TOW(II,2) = TOW(II,2) - DES_RADIUS(II)*(NORM(3)*FT(L,1) - NORM(1)*&
                    FT(L,3))
         TOW(II,3) = TOW(II,3) - DES_RADIUS(II)*(NORM(1)*FT(L,2) - NORM(2)*&
                    FT(L,1))
      ELSE
         TOW(L,:) = TOW(L,:) + DES_RADIUS(L)*(NORM(1)*FT(L,2) - NORM(2)*&
                    FT(L,1))
         TOW(II,:) = TOW(L,:) - DES_RADIUS(II)*(NORM(1)*FT(L,2) - NORM(2)*&
                    FT(L,1))
      END IF

      RETURN
      END SUBROUTINE CFFCTOW


