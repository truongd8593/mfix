!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFSLIDE(L, TANGNT_VREL, FT1)                           C
!  Purpose: DES - calculate sliding between particles                  C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Comments: Implements Eqns 9 & 10 from the following paper           C
!  Tsuji Y., Kawaguchi T., and Tanak T., "Lagrangian numerical         C
!  simulation of plug glow of cohesionless particles in a              C
!  horizontal pipe", Powder technology, 71, 239-250, 1992              C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CFSLIDE(L, TANGNT, TEMP_FT)

      USE param1
      USE discretelement
      IMPLICIT NONE

      DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT

      INTEGER L, K
      DOUBLE PRECISION FTMD, FNMD, TANGNT(DIMN)
      DOUBLE PRECISION TEMP_FT(DIMN), TEMP_FN(DIMN)
!     
!---------------------------------------------------------------------

      TEMP_FN(:) = FN(L,:)

      FTMD = SQRT(DES_DOTPRDCT(TEMP_FT,TEMP_FT))
      FNMD = SQRT(DES_DOTPRDCT(TEMP_FN,TEMP_FN))

      IF (FTMD.GT.(MEW*FNMD)) THEN
          PARTICLE_SLIDE = .TRUE.
          FT(L,:) = - MEW*FNMD*TANGNT(:)
      ELSE
         FT(L,:) = TEMP_FT(:)
      END IF

      RETURN
      END SUBROUTINE CFSLIDE


