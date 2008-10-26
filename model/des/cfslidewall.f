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

      SUBROUTINE CFSLIDEWALL(L, TANGNT, TEMP_FT)
      
      USE discretelement
      USE param1
      IMPLICIT NONE

      DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT

      INTEGER L, K
      DOUBLE PRECISION FTMD, FNMD, TANGNT(DIMN)
      DOUBLE PRECISION TEMP_FT(DIMN), TEMP_FN(DIMN)
!     
!---------------------------------------------------------------------


      TEMP_FN(:) = FN(L, :)

      FTMD = SQRT(DES_DOTPRDCT(TEMP_FT,TEMP_FT))
      FNMD = SQRT(DES_DOTPRDCT(TEMP_FN,TEMP_FN))

      IF (FTMD.GT.(MEW_W*FNMD)) THEN
         IF(DEBUG_DES) PRINT*,'From cfslidewall.f'
         IF(DEBUG_DES) PRINT*,'SLIDE, FTMD, mu*FNMD = ', FTMD, MEW*FNMD
         PARTICLE_SLIDE = .TRUE.
         FT(L,:) = - MEW_W*FNMD*TANGNT(:)
      ELSE
         FT(L,:) = TEMP_FT(:)
      END IF
      

      RETURN
      END SUBROUTINE CFSLIDEWALL


