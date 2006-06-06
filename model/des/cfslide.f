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

      SUBROUTINE CFSLIDE(L, TANGNT, FT1)

      USE param1
      USE discretelement
      IMPLICIT NONE

      DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT

      INTEGER L, K
      DOUBLE PRECISION FTMD, FNMD,  FT1(DIMN), TANGNT(DIMN)
!     
!---------------------------------------------------------------------

      FTMD = SQRT(DES_DOTPRDCT(FT1,FT1))
      FNMD = SQRT(DES_DOTPRDCT(FN,FN))

      IF (FTMD.GT.(MEW*FNMD)) THEN
          PARTICLE_SLIDE = .TRUE.
          FT(L,:) = - MEW*FNMD*TANGNT(:)
      ELSE
         FT(L,:) = FT1(:)
      END IF

      RETURN
      END SUBROUTINE CFSLIDE


