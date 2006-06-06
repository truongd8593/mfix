!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFSLIDEWALL(L, TANGNT)                                 C
!  Purpose: DES - Calculate slide between particles and walls          C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CFSLIDEWALL(L, TANGNT)
      
      USE discretelement
      USE param1
      IMPLICIT NONE

      DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT

      INTEGER L, K
      DOUBLE PRECISION FTMD, FNMD, TANGNT(DIMN), FT1(DIMN)
!     
!---------------------------------------------------------------------


      FTMD = SQRT(DES_DOTPRDCT(FT,FT))
      FNMD = SQRT(DES_DOTPRDCT(FN,FN))

      IF (FTMD.GT.(MEW_W*FNMD)) THEN
         PARTICLE_SLIDE = .TRUE.
         FT(L,:) = -MEW_W*FNMD*TANGNT(:)
      END IF

      RETURN
      END SUBROUTINE CFSLIDEWALL


