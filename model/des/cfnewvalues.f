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
     

         DES_VEL_NEW(L,:) = FC(L,:)/PMASS(L) + GRAV(:)
         DES_VEL_NEW(L,:) = DES_VEL_OLD(L,:) + DES_VEL_NEW(L,:)*DTSOLID

         DES_POS_NEW(L,:) = DES_POS_OLD(L,:) + DES_VEL_NEW(L,:)*DTSOLID 

	 OMEGA_NEW(L,:) = OMEGA_OLD(L,:) + TOW(L,:)*OMOI(L)*DTSOLID


         
      DIST = ZERO 
      DO KK=1, DIMN
         DIST = DIST + (DES_POS_NEW(L,KK)-DES_POS_OLD(L,KK))**2
      END DO
      DIST = SQRT(DIST)

      IF(DIST.GE.2*DES_RADIUS(L)) THEN
         PRINT *,'MOVEMENT UNDESIRED'
         PRINT *, 'T, DT', S_TIME, DTSOLID 
         PRINT *, 'L, DIST, RADIUS', L, DIST, DES_RADIUS(L)
         PRINT *, 'POS: OLD, NEW', (DES_POS_OLD(L,KK),KK=1,DIMN), (DES_POS_NEW(L,KK), KK= 1,DIMN)
         PRINT *, 'VEL: OLD, NEW', (DES_VEL_OLD(L,KK),KK=1,DIMN), (DES_VEL_NEW(L,KK),KK=1,DIMN)
         PRINT *, 'FC', (FC(L,KK),KK=1,DIMN)
         PRINT *, 'NEIGHBORS', (NEIGHBOURS(L,KK), KK=1, MAXNEIGHBORS)
         STOP
      END IF

      IF(DES_PERIODIC_WALLS) THEN
        IF(DES_PERIODIC_WALLS_X) THEN
         IF(DES_POS_NEW(L,1).GT.EX2) THEN
            DES_POS_NEW(L,1) = DES_POS_NEW(L,1) - (EX2 - WX1)
         ELSE IF(DES_POS_NEW(L,1).LT.WX1) THEN
            DES_POS_NEW(L,1) = DES_POS_NEW(L,1) + (EX2 - WX1)
         END IF
        END IF
        IF(DES_PERIODIC_WALLS_Y) THEN
         IF(DES_POS_NEW(L,2).GT.TY2) THEN
            DES_POS_NEW(L,2) = DES_POS_NEW(L,2) - (TY2 - BY1)
         ELSE IF(DES_POS_NEW(L,2).LT.BY1) THEN
            DES_POS_NEW(L,2) = DES_POS_NEW(L,2) + (TY2 - BY1)
         END IF
        END IF
        IF(DES_PERIODIC_WALLS_Z) THEN
         IF(DES_POS_NEW(L,3).GT.NZ2) THEN
            DES_POS_NEW(L,3) = DES_POS_NEW(L,3) - (NZ2 - SZ1)
         ELSE IF(DES_POS_NEW(L,3).LT.SZ1) THEN
            DES_POS_NEW(L,3) = DES_POS_NEW(L,3) + (NZ2 - SZ1)
         END IF
        END IF
      END IF

      IF(INLET_OUTLET) THEN
         IF(DES_POS_NEW(L,1).GT.(EX2+DES_RADIUS(L))) THEN
            DES_POS_NEW(L,1) = -1000
         END IF
      END IF

      V = ZERO 
      DO K = 1, DIMN
         V = V + DES_VEL_NEW(L,K)**2
      END DO 

      DES_KE = DES_KE + HALF*V*PMASS(L) 

      IF(DIMN.EQ.3) THEN
         DES_PE = DES_PE + PMASS(L)*SQRT(GRAV(1)**2 +&
                  GRAV(2)**2 + GRAV(3)**2)*DES_POS_NEW(L,2)
      ELSE
         DES_PE = DES_PE + PMASS(L)*SQRT(GRAV(1)**2 +&
                  GRAV(2)**2)*DES_POS_NEW(L,2)
      END IF

      DTSOLID = TEMPTIME

      FC(L,:) = ZERO
      TOW(L,:) = ZERO

      RETURN
      END SUBROUTINE CFNEWVALUES


