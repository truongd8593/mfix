!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFNEWVALUES(L)                                         C
!  Purpose: DES - Calculate the new values of particle velocity,       C
!           position, angular velocity etc                             C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CFNEWVALUES(L)

      USE param
      USE param1
      USE parallel
      USE matrix
      USE scales
      USE constant
      USE compar
      USE sendrecv
      USE physprop
      USE fldvar
      USE visc_g
      USE rxns
      USE run
      USE geometry
      USE indices
      USE drag
      USE discretelement

      IMPLICIT NONE
      
      INTEGER L, KK, K, NSPLIT, CHECK
      INTEGER IJK, I, J, KKK 
      DOUBLE PRECISION TEMPTIME, PASSTIME, DIST, V
!     
!---------------------------------------------------------------------
!     Calculate new values
!---------------------------------------------------------------------
!     
      CHECK = 0
      TEMPTIME = DTSOLID

      DO KK = 1, DIMN

!     PRINT *,KK, L,':', FC(KK,L), PMASS(L), FC(KK,L)/PMASS(L), GRAV(KK) 

         DES_VEL_NEW(KK,L) = FC(KK,L)/PMASS(L) + GRAV(KK)

         DES_VEL_NEW(KK,L) = DES_VEL_OLD(KK,L)+DES_VEL_NEW(KK,L)*DTSOLID

	 DES_POS_NEW(KK,L)=DES_VEL_NEW(KK,L)*DTSOLID+DES_POS_OLD(KK,L)

	 OMEGA_NEW(KK,L) = OMEGA_OLD(KK,L) + (TOW(KK,L)/MOI(L))*DTSOLID
      END DO

      OMEGA_NEW(3,L) = OMEGA_OLD(3,L) + (TOW(3,L)/MOI(L))*DTSOLID


!     IF(L.EQ.1) THEN
!     PRINT *,CALLED, L,':',FC(2,L), DES_VEL_OLD(2,L), DES_VEL_NEW(2,L), DES_POS_NEW(2,L), DTSOLID, gr
!     END IF

      DIST = 0.0
      DO KK=1, DIMN
         DIST = DIST + (DES_POS_NEW(KK,L)-DES_POS_OLD(KK,L))**2
      END DO
      DIST = SQRT(DIST)

      IF(DIST.GE.2*DES_RADIUS(L)) THEN
         PRINT *,'MOVEMENT UNDESIRED', CALLED, L, DTSOLID, DIST, DES_RADIUS(L)
         PRINT *, (DES_POS_OLD(KK,L),KK=1,DIMN), (DES_POS_NEW(KK,L), KK= 1,DIMN)
         PRINT *, (DES_VEL_OLD(KK,L),KK=1,DIMN), (DES_VEL_NEW(KK,L),KK=1,DIMN)
         PRINT *, (FC(KK,L),KK=1,DIMN)
         PRINT *,(NEIGHBOURS(KK,L), KK=1, MAXNEIGHBORS)
         STOP
      END IF

      IF(DES_PERIODIC_WALLS) THEN
        IF(DES_PERIODIC_WALLS_X) THEN
         IF(DES_POS_NEW(1,L).GT.EX2) THEN
            DES_POS_NEW(1,L) = DES_POS_NEW(1,L) - (EX2 - WX1)
         ELSE IF(DES_POS_NEW(1,L).LT.WX1) THEN
            DES_POS_NEW(1,L) = DES_POS_NEW(1,L) + (EX2 - WX1)
         END IF
        END IF
      END IF

      IF(INLET_OUTLET) THEN
         IF(DES_POS_NEW(1,L).GT.(EX2+DES_RADIUS(L))) THEN
            DES_POS_NEW(1,L) = -1000
         END IF
      END IF

      V = 0.0
      DO K = 1, DIMN
         V = V + DES_VEL_NEW(K,L)**2
      END DO 

      DES_KE = DES_KE + 0.5*V*PMASS(L)
      DES_PE = DES_PE + PMASS(L)*SQRT(GRAV(1)**2 +&
              GRAV(2)**2 + GRAV(3)**2)*DES_POS_NEW(2,L)

      DTSOLID = TEMPTIME

      RETURN
      END SUBROUTINE CFNEWVALUES


