!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: PERIODIC_WALL_CALC_FORCE_DES(C)                        C
!  Purpose: DES - DES calculaton for a periodic boundary flow          C
!                                                                      C
!                                                                C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!  REVISED aug 24 2005                                                                    C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE PERIODIC_WALL_CALC_FORCE_DES

      USE param1
      USE run
      USE discretelement
      USE geometry
      IMPLICIT NONE
      
      INTEGER LL, I, II, K, LC, C, CO, IW, KK, TEMP, TEMPN
      DOUBLE PRECISION OVERLAP_N, OVERLAP_T, TEMPFN(DIMN), TEMPFT(DIMN)
      DOUBLE PRECISION NORMAL(DIMN), VRE(DIMN)
      DOUBLE PRECISION TANGENT(DIMN)
      DOUBLE PRECISION Vn, Vt
      INTEGER NI, NWS, IJ, WALLCHECK
      INTEGER WWALL(1000), EWALL(1000)
      INTEGER BWALL(1000), TWALL(1000)
      INTEGER SWALL(1000), NWALL(1000)
!     
!---------------------------------------------------------------------
!     Calculate new values
!---------------------------------------------------------------------
!     

      IF (TIME.LE.DTSOLID) THEN
         DO K = 1, DIMN
            VRE(K) = ZERO
            TANGENT(K) = ZERO
            NORMAL(K) = ZERO
         END DO
         DO LC = 1, PARTICLES
            DO K = 1, DIMN
               FN(LC,K) = ZERO
               FT(LC,K) = ZERO
            END DO
            DO KK = 1, MN
               PN(LC,KK) = -1
               PV(LC,KK) = 1
                    DO K = 1, DIMN
                       PFN(LC,KK,K) = ZERO
                       PFT(LC,KK,K) = ZERO
                    END DO
                 END DO
                 PN(LC,1) = 0
              END DO
           END IF
           
           IF(DES_CONTINUUM_COUPLED) THEN
              CALL PARTICLES_IN_CELL(PARTICLES)
           END IF

!     
!---------------------------------------------------------------------
!     Calculate contact force and torque
!---------------------------------------------------------------------
!     
           
           DO K = 1, 1000
              WWALL(K) = 0
              EWALL(K) = 0
              BWALL(K) = 0
              TWALL(K) = 0
              SWALL(K) = 0
              NWALL(K) = 0
           END DO

              IF(DES_PERIODIC_WALLS_X) THEN

                CALL CF_PERIODIC_WALL_X(WWALL, EWALL)
                      
                IF (WWALL(1).GT.0) THEN
                DO K = 2, WWALL(1) + 1
                   LL = WWALL(K)
                   CALL CF_PERIODIC_WALL_NEIGHBOR_X(LL, WWALL, EWALL)
                END DO
                END IF
                
                IF(EWALL(1).GT.0) THEN
                DO K = 2, EWALL(1) + 1
                   LL = EWALL(K)
                   CALL CF_PERIODIC_WALL_NEIGHBOR_X(LL, WWALL, EWALL)
                END DO
                END IF
              END IF
              
              IF(DES_PERIODIC_WALLS_Y) THEN

                CALL CF_PERIODIC_WALL_Y(BWALL,TWALL)
 
                IF(BWALL(1).GT.0) THEN
                DO K = 2, BWALL(1) + 1
                   LL = BWALL(K)
                   CALL CF_PERIODIC_WALL_NEIGHBOR_Y(LL, BWALL, TWALL)
                END DO
                END IF
                
                IF(TWALL(1).GT.0) THEN
                DO K = 2, TWALL(1) + 1
                   LL = TWALL(K)
                   CALL CF_PERIODIC_WALL_NEIGHBOR_Y(LL, BWALL, TWALL)
                END DO
                END IF
              END IF
              
              IF(DIMN.EQ.3) THEN
              IF(DES_PERIODIC_WALLS_z) THEN

                CALL CF_PERIODIC_WALL_Z(SWALL, NWALL)
                      
                IF(SWALL(1).GT.0) THEN
                DO K = 2, SWALL(1) + 1
                   LL = SWALL(K)
                   CALL CF_PERIODIC_WALL_NEIGHBOR_Z(LL, SWALL, NWALL)
                END DO
                END IF
                
                IF(NWALL(1).GT.0) THEN
                DO K = 2, NWALL(1) + 1
                   LL = NWALL(K)
                   CALL CF_PERIODIC_WALL_NEIGHBOR_Z(LL, SWALL, NWALL)
                END DO
                END IF
              END IF
              END IF
           
           DO LL = 1, PARTICLES
            
              NI = 0
              TEMPFN(:) = ZERO
              TEMPFT(:) = ZERO
              K = 0
              KK = 0

              IF(PN(LL,1).GT.0) THEN
                 DO K = 2, PN(LL,1)+1
                    IF(PV(LL,K).EQ.0) THEN
                       PN(LL,K) = UNDEFINED_I 
                       PFN(LL,K,:) = ZERO  
                       PFT(LL,K,:) = ZERO
                    END IF
                 END DO
              END IF

              IF(PN(LL,1).GT.1) THEN
                 DO CO = 3, PN(LL,1)+1
                    TEMPN = PN(LL,CO)
                    TEMPFN(:) = PFN(LL,CO,:)
                    TEMPFT(:) = PFT(LL,CO,:)
                    NI = CO - 1
                    DO WHILE((NI.GT.1).AND.(PN(LL,NI).GT.TEMPN))
                       PN(LL,NI+1) = PN(LL,NI)
                       PFN(LL,NI+1,:) = PFN(LL,NI,:)
                       PFT(Ll,NI+1,:) = PFT(LL,NI,:)
                       NI = NI-1
                    END DO
                    PN(LL,NI+1) = TEMPN
                    PFN(LL,NI+1,:) = TEMPFN(:)
                    PFT(LL,NI+1,:) = TEMPFT(:)
                 END DO
              END IF

              IF(PN(LL,1).GT.0) THEN
                 IJ = 0
                 DO NI = 2, PN(LL,1)+1 
                    IF(PN(LL,NI).GE.UNDEFINED_I) THEN
                       PN(LL,NI) = -1
                       IJ = IJ + 1
                    END IF
                 END DO
                 PN(LL,1) = PN(LL,1) - IJ 
              END IF

              IF (PN(LL,1).EQ.0) THEN
                 DO NI = 1, MAXNEIGHBORS 
                    PFN(LL,NI,:) = ZERO
                    PFT(LL,NI,:) = ZERO
                 END DO
              END IF

              DO K = 2, MAXNEIGHBORS
                 PV(LL,K) = 0
              END DO

              IF(WALLDTSPLIT) THEN
                 WALLCHECK = 0
                 NWS = 2*DIMN
                 DO IW = 1, NWS
                    IF(DES_PERIODIC_WALLS_X.AND.(IW.NE.1).AND.(IW.NE.2)) THEN
                       WALLCONTACT = 0
                       CALL CFWALLCONTACT(IW, LL, WALLCONTACT)
                       IF(WALLCONTACT.EQ.1) THEN
                          WALLCHECK = 1
                          CALL CFWALLPOSVEL(LL, IW)
                          I = PARTICLES + IW
                             DES_POS_NEW(I,:) = DES_WALL_POS(IW,:)
                             DES_VEL_NEW(I,:) = DES_WALL_VEL(IW,:)
                             OMEGA_NEW(I,:) = ZERO
                          DES_RADIUS(I) = DES_RADIUS(LL)
                          CALL CFNORMAL(LL, I, NORMAL)
                          CALL CFTANGENT(TANGENT, NORMAL, VRE)
                          CALL CFRELVEL(LL, I, VRE, TANGENT)
!                          CALL CFTANGENT(TANGENT, NORMAL, VRE)
                          CALL CFVRN(Vn, VRE, NORMAL)
                          CALL CFVRT(Vt, VRE, TANGENT)
                          CALL CFTOTALOVERLAPS(LL, I, Vt, OVERLAP_N, OVERLAP_T)
                          CALL CFFNWALL(LL, Vn, OVERLAP_N, NORMAL)
                          CALL CFFTWALL(LL, Vt, OVERLAP_T, TANGENT)
                          CALL CFSLIDEWALL(LL, TANGENT)
                          CALL CFFCTOWALL(LL, NORMAL)
                       END IF
                    END IF

                    IF(DES_PERIODIC_WALLS_Y.AND.(IW.NE.3).AND.(IW.NE.4)) THEN
                       WALLCONTACT = 0
                       CALL CFWALLCONTACT(IW, LL, WALLCONTACT)
                       IF(WALLCONTACT.EQ.1) THEN
                          WALLCHECK = 1
                          CALL CFWALLPOSVEL(LL, IW)
                          I = PARTICLES + IW
                             DES_POS_NEW(I,:) = DES_WALL_POS(IW,:)
                             DES_VEL_NEW(I,:) = DES_WALL_VEL(IW,:)
                             OMEGA_NEW(I,:) = ZERO
                          DES_RADIUS(I) = DES_RADIUS(LL)
                          CALL CFNORMAL(LL, I, NORMAL)
                          CALL CFTANGENT(TANGENT, NORMAL, VRE)
                          CALL CFRELVEL(LL, I, VRE, TANGENT)
!                          CALL CFTANGENT(TANGENT, NORMAL, VRE)
                          CALL CFVRN(Vn, VRE, NORMAL)
                          CALL CFVRT(Vt, VRE, TANGENT)
                          CALL CFTOTALOVERLAPS(LL, I, Vt, OVERLAP_N, OVERLAP_T)
                          CALL CFFNWALL(LL, Vn, OVERLAP_N, NORMAL)
                          CALL CFFTWALL(LL, Vt, OVERLAP_T, TANGENT)
                          CALL CFSLIDEWALL(LL, TANGENT)
                          CALL CFFCTOWALL(LL, NORMAL)
                       END IF
                    END IF

                    IF(DES_PERIODIC_WALLS_Z.AND.(IW.NE.5).AND.(IW.NE.6)) THEN
                       WALLCONTACT = 0
                       CALL CFWALLCONTACT(IW, LL, WALLCONTACT)
                       IF(WALLCONTACT.EQ.1) THEN
                          WALLCHECK = 1
                          CALL CFWALLPOSVEL(LL, IW)
                          I = PARTICLES + IW
                             DES_POS_NEW(I,:) = DES_WALL_POS(IW,:)
                             DES_VEL_NEW(I,:) = DES_WALL_VEL(IW,:)
                             OMEGA_NEW(I,:) = ZERO
                          DES_RADIUS(I) = DES_RADIUS(LL)
                          CALL CFNORMAL(LL, I, NORMAL)
                          CALL CFTANGENT(TANGENT, NORMAL, VRE)
                          CALL CFRELVEL(LL, I, VRE, TANGENT)
!                          CALL CFTANGENT(TANGENT, NORMAL, VRE)
                          CALL CFVRN(Vn, VRE, NORMAL)
                          CALL CFVRT(Vt, VRE, TANGENT)
                          CALL CFTOTALOVERLAPS(LL, I, Vt, OVERLAP_N, OVERLAP_T)
                          CALL CFFNWALL(LL, Vn, OVERLAP_N, NORMAL)
                          CALL CFFTWALL(LL, Vt, OVERLAP_T, TANGENT)
                          CALL CFSLIDEWALL(LL, TANGENT)
                          CALL CFFCTOWALL(LL, NORMAL)
                       END IF
                    END IF
                    
                 END DO
              END IF

              IF (NEIGHBOURS(LL,1).GT.0) THEN
                 DO II = 2, NEIGHBOURS(LL,1)+1
                   I = NEIGHBOURS(LL,II)
                   IF(I.GT.LL) THEN
	      	    CO = 0
		    NI = 2
                    IF(PN(1,LL).GT.0) THEN
                       DO WHILE((CO.EQ.0).AND.(NI.LE.(PN(LL,1)+1)))
                          IF(I.EQ.PN(LL,NI)) THEN
                             CO = 1
                             PV(LL,NI) = 1
                             CALL CFNORMAL(LL, I, NORMAL)
                             CALL CFTANGENT(TANGENT, NORMAL, VRE)
                             CALL CFRELVEL(LL, I, VRE, TANGENT)
!                             CALL CFTANGENT(TANGENT, NORMAL, VRE)
                             CALL CFVRN(Vn, VRE, NORMAL)
                             CALL CFVRT(Vt, VRE, TANGENT)
                             CALL CFINCREMENTALOVERLAPS(Vn, Vt, OVERLAP_N, OVERLAP_T)
                             CALL CFFN(LL, Vn, OVERLAP_N, NORMAL)
                             CALL CFFT(LL, Vt, OVERLAP_T, TANGENT)
                             FN(LL,:) = FN(LL,:) + PFN(LL,NI,:)
                             TEMPFT(:) = FT(LL,:) + PFT(LL,NI,:)
                             CALL CFSLIDE(LL, TANGENT, TEMPFT)
                             CALL CFFCTOW(LL, I, NORMAL)
                             PFN(LL,NI,:) = PFN(LL,NI,:) + FNS1(:)
                             PFT(LL,NI,:) = PFT(LL,NI,:) + FTS1(:)
                          ELSE
                             NI = NI + 1
                          END IF

                       END DO
                    END IF
		    IF(CO.EQ.0) THEN
                       PN(LL,1) = PN(LL,1) + 1
                       NI = PN(LL,1) + 1
                       PN(LL,NI) = I
                       PV(LL,NI) = 1
                       CALL CFNORMAL(LL, I, NORMAL)
                       CALL CFTANGENT(TANGENT, NORMAL, VRE)
                       CALL CFRELVEL(LL, I, VRE, TANGENT)
!                       CALL CFTANGENT(TANGENT, NORMAL, VRE)
                       CALL CFVRN(Vn, VRE, NORMAL)
                       CALL CFVRT(Vt, VRE, TANGENT)
                       CALL CFTOTALOVERLAPS(LL, I, Vt, OVERLAP_N, OVERLAP_T)
                       CALL CFFN(LL, Vn, OVERLAP_N, NORMAL)
                       CALL CFFT(LL, Vt, OVERLAP_T, TANGENT)
                       TEMPFT(:) = FT(LL,:)
                       CALL CFSLIDE(LL, TANGENT, TEMPFT)
                       CALL CFFCTOW(LL, I, NORMAL)
                       PFN(LL,NI,:) = FNS1(:)
                       PFT(LL,NI,:) = FTS1(:)
		    END IF
                  END IF
                 END DO
              END IF

              IF((NEIGHBOURS(LL,1).EQ.0).AND.(WALLCHECK.EQ.0)) THEN
                 CALL CFNOCONTACT(LL)
              END IF

           END DO

!!COHESION
        IF(USE_COHESION)THEN
          CALL CALC_COHESIVE_FORCES
        END IF
!!COHESION

           IF(DES_CONTINUUM_COUPLED) THEN
              CALL DRAG_FGS(PARTICLES)
           END IF


!-------------------------------------------------------------------
!     Update old values with new values
!-------------------------------------------------------------------
           CALL CFUPDATEOLD(PARTICLES)

           RETURN
           END SUBROUTINE PERIODIC_WALL_CALC_FORCE_DES

