!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFASSIGN                                               C
!  Purpose: DES - To assign values for particles                       C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer: Sreekanth Pannala                        Date: 09-Nov-06  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CFASSIGN

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

      INTEGER L, IJK, M, I,J, K
      DOUBLE PRECISION FOUR_BY_THREE, RAD2, MINMASS
!     
!---------------------------------------------------------------------
!     Assignments
!---------------------------------------------------------------------
!     

      INCLUDE 'b_force1.inc'
      INCLUDE 'b_force2.inc'
      
      FOUR_BY_THREE = 4.0d0/3.0d0
      MINMASS = LARGE_NUMBER
      MAX_RADIUS = ZERO
      MIN_RADIUS = LARGE_NUMBER
      DO L = 1, PARTICLES
         RAD2 = DES_RADIUS(L)**2
         PVOL(L) = FOUR_BY_THREE*Pi*RAD2*DES_RADIUS(L)
         PMASS(L) = PVOL(L)*RO_Sol(L) 
         OMOI(L) = 2.5d0/(PMASS(L)*RAD2) !one over MOI
         MAX_RADIUS = MAX(MAX_RADIUS, DES_RADIUS(L))
         MIN_RADIUS = MIN(MIN_RADIUS, DES_RADIUS(L))
         IF(PMASS(L).LT.MINMASS) MINMASS = PMASS(L) 
      END DO
      
      
      RADIUS_EQ = DES_RADIUS(1)*1.05D0
      NEIGHBOR_SEARCH_DIST = NEIGHBOR_SEARCH_RAD_RATIO * RADIUS_EQ

      DTSOLID = DTSOLID_FACTOR*2.0D0*PI*SQRT(MINMASS/(15*KN)) ! DTs - Rotational Constraint
!     DTSOLID = DTSOLID_FACTOR*2D0*PI*SQRT(MINMASS/(6*KN)) ! DTs - Translational Constraint

      !Print*,'DTSOLID = ', dtsolid
      !read(*,*)
      WX1 = ZERO 
      EX2 = XLENGTH 
      BY1 = ZERO
      TY2 = YLENGTH 
      IF((DIMN.EQ.2).AND.(COORDINATES == 'CARTESIAN')) THEN
         SZ1 = ZERO 
         NZ2 = 2*RADIUS_EQ
      ELSE
         SZ1 = ZERO
         NZ2 = ZLENGTH
        ! WRITE(*,*) 'XLENGHT =', XLENGTH, YLENGTH, ZLENGTH
      ENDIF

      !IF((DIMN.EQ.2).AND.(COORDINATES == 'CARTESIAN')) THEN
      !   DZ(:) = 2D0*RADIUS_EQ
      !END IF

      GRAV(1) = BFX_s(1,1)
      GRAV(2) = BFY_s(1,1)
      IF(DIMN.EQ.3) GRAV(3) = BFZ_s(1,1)

      RETURN
      END SUBROUTINE CFASSIGN
      
