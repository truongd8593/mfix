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
      DOUBLE PRECISION FTMD, FNMD, TANGNT(DIMN)
      DOUBLE PRECISION TEMP_FT(DIMN), TEMP_FN(DIMN)
!     
!---------------------------------------------------------------------


      TEMP_FT(:) = FT(L,:)
      TEMP_FN(:) = FN(L,:)

      FTMD = SQRT(DES_DOTPRDCT(TEMP_FT,TEMP_FT))
      FNMD = SQRT(DES_DOTPRDCT(TEMP_FN,TEMP_FN))

      IF (FTMD.GT.(MEW_W*FNMD)) THEN
         PARTICLE_SLIDE = .TRUE.
         FT(L,:) = -MEW_W*FNMD*TANGNT(:)
      END IF

      RETURN
      END SUBROUTINE CFSLIDEWALL


