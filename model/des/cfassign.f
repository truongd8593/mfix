!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFASSIGN(PARTS)                                        C
!  Purpose: DES - To assign values for particles                       C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CFASSIGN(PARTS)

      USE discretelement
      USE param
      USE param1
      USE parallel
      USE fldvar
      USE run
      USE geometry
      USE matrix
      USE indices
      USE physprop
      USE drag
      USE constant
      USE compar
      USE sendrecv

      IMPLICIT NONE

      INTEGER L, PARTS, IJK, M, I,J, K
      DOUBLE PRECISION FOUR_BY_THREE, RAD2, MINMASS
!     
!---------------------------------------------------------------------
!     Assignments
!---------------------------------------------------------------------
!     

      INCLUDE 'b_force1.inc'
      INCLUDE 'b_force2.inc'
      
      FOUR_BY_THREE = 4.0/3.0
      MINMASS = LARGE_NUMBER
      
         DO L = 1, PARTS
            RAD2 = DES_RADIUS(L)**2
            PVOL(L) = FOUR_BY_THREE*Pi*RAD2*DES_RADIUS(L)
            PMASS(L) = PVOL(L)*RO_Sol(L) 
            OMOI(L) = 2.5/(PMASS(L)*RAD2) !one over MOI
            IF(PMASS(L).LT.MINMASS) MINMASS = PMASS(L) 
         END DO

         RADIUS_EQ = DES_RADIUS(1)*1.05D0
         NEIGHBOR_SEARCH_DIST = NEIGHBOR_SEARCH_RAD_RATIO * RADIUS_EQ

         DTSOLID = DTSOLID_FACTOR*2D0*PI*SQRT(MINMASS/(15*KN)) ! DTs - Rotational Constraint
!         DTSOLID = DTSOLID_FACTOR*2D0*PI*SQRT(MINMASS/(6*KN)) ! DTs - Translational Constraint

         WX1 = ZERO 
         EX2 = XLENGTH 
         BY1 = ZERO
         TY2 = YLENGTH 
         SZ1 = ZERO
         NZ2 = ZLENGTH
         SZ1 = ZERO 
         NZ2 = 2*RADIUS_EQ

         IF((DIMN.EQ.2).AND.(COORDINATES == 'CARTESIAN')) THEN
            DZ(:) = 2D0*RADIUS_EQ
         END IF

         GRAV(1) = BFX_s(1,1)
         GRAV(2) = BFY_s(1,1)
         IF(DIMN.EQ.3) GRAV(3) = BFZ_s(1,1)

      RETURN
      END SUBROUTINE CFASSIGN



