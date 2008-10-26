!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFNEWVALUES(L)                                         C
!  Purpose: DES - Calculate the new values of particle velocity,       C
!           position, angular velocity etc                             C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C 
!  Comments: Implements Eqns 1, 2, 3, 4 & 5  from the following paper  C
!  Tsuji Y., Kawaguchi T., and Tanak T., "Lagrangian numerical         C
!  simulation of plug glow of cohesionless particles in a              C
!  horizontal pipe", Powder technology, 71, 239-250, 1992              C
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

      DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT
      
      INTEGER L, KK, K, NSPLIT, CHECK
      INTEGER IJK, I, J, KKK 
      DOUBLE PRECISION TEMPTIME, PASSTIME, D(DIMN), DIST,  V, rhat(dimn), rhat2, force_tmpy
!     
!---------------------------------------------------------------------
!     Calculate new values
!---------------------------------------------------------------------
!     
      IF(L.EQ.1) THEN 
         DES_KE = ZERO
         DES_PE = ZERO 
         DES_VEL_AVG = ZERO
      ENDIF
      CHECK = 0 

      TEMPTIME = DTSOLID
      FC(L, :) = FC(L,:)/PMASS(L) + GRAV(:)

      DES_VEL_NEW(L,:) = FC(L,:)
      
      DES_VEL_NEW(L,:) = DES_VEL_OLD(L,:) + DES_VEL_NEW(L,:)*DTSOLID
      DES_POS_NEW(L,:) = DES_POS_OLD(L,:) + DES_VEL_NEW(L,:)*DTSOLID 
      OMEGA_NEW(L,:) = OMEGA_OLD(L,:) + TOW(L,:)*OMOI(L)*DTSOLID
      IF(DES_POS_NEW(L,2).LT.ZERO.AND..NOT.DES_PERIODIC_WALLS) THEN
         PRINT*,'POSITION LE ZERO FOR L = ', L, DES_VEL_NEW(L,:), DES_POS_NEW(L,:)
      ENDIF 
!PRINT*,'grav = ', GRAV(2), pgrad(2)*PVOL(L)/PMASS(L), RO_sol(L)
      IF(.NOT.DO_NSEARCH) THEN
         D(:) = DES_POS_NEW(L,:) - PPOS(L,:)
         DIST = SQRT(DES_DOTPRDCT(D,D))
         
         NEIGHBOR_SEARCH_DIST = NEIGHBOR_SEARCH_RAD_RATIO*DES_RADIUS(L)
         IF(DIST.GE.NEIGHBOR_SEARCH_DIST) DO_NSEARCH = .TRUE.
      END IF

!     Chacking if the partcile moved more than a dia in a solid time step   
      D(:) = DES_POS_NEW(L,:) - DES_POS_OLD(L,:)
      DIST = SQRT(DES_DOTPRDCT(D,D))
      IF(DIST.GE.DES_RADIUS(L)) THEN
         PRINT *,'MOVEMENT UNDESIRED: PARTICLE', L
         STOP
      END IF

!     Periodic treatment
      IF(DES_PERIODIC_WALLS) THEN
         IF(DES_PERIODIC_WALLS_X) THEN
            IF(DES_POS_NEW(L,1).GT.EX2) THEN
               DES_POS_NEW(L,1) = DES_POS_NEW(L,1) - (EX2 - WX1)
               PIJK(L,1) = 2
            ELSE IF(DES_POS_NEW(L,1).LT.WX1) THEN
               DES_POS_NEW(L,1) = DES_POS_NEW(L,1) + (EX2 - WX1)
               PIJK(L,1) = IMAX1
            END IF
         END IF
         IF(DES_PERIODIC_WALLS_Y) THEN
            IF(DES_POS_NEW(L,2).GT.TY2) THEN
               DES_POS_NEW(L,2) = DES_POS_NEW(L,2) - (TY2 - BY1)
               PIJK(L,2) = 2
            ELSE IF(DES_POS_NEW(L,2).LT.BY1) THEN
               DES_POS_NEW(L,2) = DES_POS_NEW(L,2) + (TY2 - BY1)
               PIJK(L,2) = JMAX1
            END IF
         END IF
         IF(DES_PERIODIC_WALLS_Z) THEN
            IF(DES_POS_NEW(L,3).GT.NZ2) THEN
               DES_POS_NEW(L,3) = DES_POS_NEW(L,3) - (NZ2 - SZ1)
               PIJK(L,3) = 2
            ELSE IF(DES_POS_NEW(L,3).LT.SZ1) THEN
               DES_POS_NEW(L,3) = DES_POS_NEW(L,3) + (NZ2 - SZ1)
               PIJK(L,3) = KMAX1
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
      
      DES_VEL_AVG(:) =  DES_VEL_AVG(:) + DES_VEL_NEW(L,:)

      DES_KE = DES_KE + HALF*V*PMASS(L) 
      DES_PE = DES_PE + PMASS(L)*DES_POS_NEW(L,2)*SQRT(DES_DOTPRDCT(GRAV,GRAV))

      DTSOLID = TEMPTIME

      FC(L,:) = ZERO
      TOW(L,:) = ZERO

      RETURN
      END SUBROUTINE CFNEWVALUES


