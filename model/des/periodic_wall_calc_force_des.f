!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: PERIODIC_WALL_CALC_FORCE_DES(C)                        C
!  Purpose: DES - DES calculaton for a periodic boundary flow          C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE PERIODIC_WALL_CALC_FORCE_DES(C)

      USE discretelement
      USE geometry
      IMPLICIT NONE
      
      INTEGER LL, I, II, K, LC, C, CO, IW, KK, TEMP, TEMPN
      DOUBLE PRECISION OVERLAP_N, OVERLAP_T, TEMPFN(NDIM), TEMPFT(NDIM)
      DOUBLE PRECISION NORMAL(NDIM), VRE(NDIM)
      DOUBLE PRECISION TANGENT(NDIM)
      DOUBLE PRECISION Vn, Vt, S_TIME
      INTEGER IJK, NWS, IJ, WALLCHECK
      INTEGER WWALL(1000), EWALL(1000)
      INTEGER BWALL(1000), TWALL(1000)
      INTEGER SWALL(1000), NWALL(1000)
!     
!---------------------------------------------------------------------
!     Assignments
!---------------------------------------------------------------------
      IF(CALLED.EQ.0) THEN
         DO LL= 1, PARTICLES
            CALL CFASSIGN(LL)
         END DO
      END IF
!---------------------------------------------------------------------
!     Calculate new values
!---------------------------------------------------------------------
!     

      IF (CALLED.EQ.0) THEN
         PRINT *,'**** PERIODIC_WALL_CALC_FORCE_DES ****'
         DO LC = 1, PARTICLES
            DO K = 1, DIMN
               FN(K,LC) = 0.0
               FT(K,LC) = 0.0
            END DO
            DO KK = 1, MN
               PN(KK,LC) = -1
               PV(KK,LC) = 1
                    DO K = 1, DIMN
                       PFN(K,KK,LC) = 0.0
                       PFT(K,KK,LC) = 0.0
                    END DO
                 END DO
                 PN(1,LC) = 0
              END DO
           END IF

           DO LC = 1, PARTICLES
              CALL CFNEWVALUES(LC)
           END DO
           
           IF(DES_CONTINUUM_COUPLED) THEN
              CALL PARTICLES_IN_CELL(PARTICLES)
           END IF

!     
!---------------------------------------------------------------------
!     Calculate contact force and torque
!---------------------------------------------------------------------
!     

           S_TIME = CALLED*DTSOLID
           CALL NEIGHBOUR(PARTICLES)

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
           ELSE IF(DES_PERIODIC_WALLS_Y) THEN
              CALL CF_PERIODIC_WALL_Y(BWALL,TWALL)
           ELSE IF(DES_PERIODIC_WALLS_Z) THEN
              CALL CF_PERIODIC_WALL_Z(SWALL, NWALL)
           END IF

           DO LL = 1, PARTICLES

              IF(DES_PERIODIC_WALLS_X) THEN
                IF(DES_POS_NEW(1,LL).LE.(WX1 + 2*RADIUS_EQ)) THEN
                   CALL CF_PERIODIC_WALL_NEIGHBOUR_X(LL, WWALL, EWALL)
                ELSE IF(DES_POS_NEW(1,LL).GE.(EX2 - 2*RADIUS_EQ))THEN
                   CALL CF_PERIODIC_WALL_NEIGHBOUR_X(LL, WWALL, EWALL)
                END IF
              END IF
              IF(DES_PERIODIC_WALLS_Y) THEN
                IF(DES_POS_NEW(2,LL).LE.(BY1 + 2*RADIUS_EQ)) THEN
                   CALL CF_PERIODIC_WALL_NEIGHBOUR_Y(LL, BWALL, TWALL)
                ELSE IF(DES_POS_NEW(2,LL).GE.(TY2 - 2*RADIUS_EQ))THEN
                   CALL CF_PERIODIC_WALL_NEIGHBOUR_Y(LL, BWALL, TWALL)
                END IF
              END IF
              IF(DES_PERIODIC_WALLS_z) THEN
                IF(DES_POS_NEW(3,LL).LE.(SZ1 + 2*RADIUS_EQ)) THEN
                   CALL CF_PERIODIC_WALL_NEIGHBOUR_Z(LL, SWALL, NWALL)
                ELSE IF(DES_POS_NEW(3,LL).GE.(NZ2 - 2*RADIUS_EQ))THEN
                   CALL CF_PERIODIC_WALL_NEIGHBOUR_Z(LL, SWALL, NWALL)
                END IF
              END IF

              IJK = 0
              DO K = 1, DIMN
                 TEMPFN(K) = 0.0
                 TEMPFT(K) = 0.0
              END DO
              K = 0
              KK = 0

              IF(NEIGHBOURS(1,LL).GT.1) THEN
                 DO CO = 3, NEIGHBOURS(1,LL)+1
                    TEMP = NEIGHBOURS(CO,LL)
                    IJK = CO - 1
                    DO WHILE((IJK.GT.1).AND.(NEIGHBOURS(IJK,LL).GT.TEMP))
                       NEIGHBOURS(IJK+1,LL) = NEIGHBOURS(IJK,LL)
                       IJK = IJK-1
                    END DO
                    NEIGHBOURS(IJK+1,LL) = TEMP
                 END DO
              END IF

              IF(PN(1,LL).GT.0) THEN
                 DO K = 2, PN(1,LL)+1
                    IF(PV(K,LL).EQ.0) THEN
                       PN(K,LL) = K*PARTICLES
                       DO KK = 1, DIMN
                          PFN(KK,K,LL) = 0.0
                          PFT(KK,K,LL) = 0.0
                       END DO
                    END IF
                 END DO
              END IF

              IF(PN(1,LL).GT.1) THEN
                 DO CO = 3, PN(1,LL)+1
                    TEMPN = PN(CO,LL)
                    DO K = 1, DIMN
                       TEMPFN(K) = PFN(K,CO,LL)
                       TEMPFT(K) = PFT(K,CO,LL)
                    END DO
                    IJK = CO - 1
                    DO WHILE((IJK.GT.1).AND.(PN(IJK,LL).GT.TEMPN))
                       PN(IJK+1,LL) = PN(IJK,LL)
                       DO K = 1, DIMN
                          PFN(K,IJK+1,LL) = PFN(K,IJK,LL)
                          PFT(K,IJK+1,LL) = PFT(K,IJK,LL)
                       END DO
                       IJK = IJK-1
                    END DO
                    PN(IJK+1,LL) = TEMPN
                    DO K = 1, DIMN
                       PFN(K,IJK+1,LL) = TEMPFN(K)
                       PFT(K,IJK+1,LL) = TEMPFT(K)
                    END DO
                 END DO
              END IF

              IF(PN(1,LL).GT.0) THEN
                 IJ = PN(1,LL)+1
                 DO K = 2, IJ
                    IF(PN(K,LL).GT.PARTICLES) THEN
                       PN(1,LL) = PN(1,LL) - 1
                       PN(K,LL) = -1
                    END IF
                 END DO
              END IF

              IF (PN(1,LL).EQ.0) THEN
                 DO IJK = 1, MN
		    DO K = 1, DIMN
                       PFN(K,IJK,LL) = 0.0
                       PFT(K,IJK,LL) = 0.0
		    END DO
                 END DO
              END IF

              DO K = 2, MN
                 PV(K,LL) = 0
              END DO

              DO K = 1, DIMN
                 FC(K,LL) = 0.0
                 TOW(K,LL) = 0.0
                 TOW(3,LL) = 0.0
              END DO

              IF(WALLDTSPLIT) THEN
                 WALLCHECK = 0
                 NWS = 2*DIMN
                 DO IW = 1, NWS
                    IF(IW.GE.3) THEN
                       WALLCONTACT = 0
                       CALL CFWALLCONTACT(IW, LL, S_TIME, WALLCONTACT)
                       IF(WALLCONTACT.EQ.1) THEN
                          WALLCHECK = 1
                          CALL CFWALLPOSVEL(LL, S_TIME, IW)
                          I = PARTICLES + IW
                          DO K = 1, DIMN
                             DES_POS_NEW(K,I) = DES_WALL_POS(K,IW)
                             DES_VEL_NEW(K,I) = DES_WALL_VEL(K,IW)
                             OMEGA_NEW(K,I) = 0.0
                             OMEGA_NEW(3,I) = 0.0
                          END DO
                          DES_RADIUS(I) = DES_RADIUS(LL)
                          CALL CFNORMAL(LL, I, NORMAL)
!                          CALL CFTANGENT(TANGENT, NORMAL)
                          CALL CFRELVEL(LL, I, VRE, TANGENT)
                          CALL CFTANGENT(TANGENT, NORMAL, VRE)
                          CALL CFVRN(Vn, VRE, NORMAL)
                          CALL CFVRT(Vt, VRE, TANGENT)
                          CALL CFTOTALOVERLAPS(LL, I, Vt, OVERLAP_N, OVERLAP_T)
                          CALL CFFNWALL(LL, Vn, OVERLAP_N, NORMAL)
                          CALL CFFTWALL(LL, Vt, OVERLAP_T, TANGENT)
                          CALL CFSLIDEWALL(LL, TANGENT)
                          CALL CFFCTOW(LL, NORMAL)
                       END IF
                    END IF
                 END DO
              END IF

              NEIGHBOR = 0
              NEIGHBOR = NEIGHBOURS(1,LL)
              IF (NEIGHBOR.GT.0) THEN
                 DO II = 2, NEIGHBOR+1
                    I = NEIGHBOURS(II,LL)
	      	    CO = 0
		    IJK = 2
                    IF(PN(1,LL).GT.0) THEN
                       DO WHILE((CO.EQ.0).AND.(IJK.LE.(PN(1,LL)+1)))
                          IF(I.EQ.PN(IJK,LL)) THEN
                             CO = 1
                             PV(IJK,LL) = 1
                             CALL CFNORMAL(LL, I, NORMAL)
!                             CALL CFTANGENT(TANGENT, NORMAL)
                             CALL CFRELVEL(LL, I, VRE, TANGENT)
                             CALL CFTANGENT(TANGENT, NORMAL, VRE)
                             CALL CFVRN(Vn, VRE, NORMAL)
                             CALL CFVRT(Vt, VRE, TANGENT)
                             CALL CFINCREMENTALOVERLAPS(Vn, Vt, OVERLAP_N, OVERLAP_T)
                             CALL CFFN(LL, Vn, OVERLAP_N, NORMAL)
                             CALL CFFT(LL, Vt, OVERLAP_T, TANGENT)
                             DO K = 1, DIMN
                                FN(K,LL) = FN(K,LL) + PFN(K,IJK,LL)
                                TEMPFT(K) = FT(K,LL) + PFT(K,IJK,LL)
                             END DO
                             CALL CFSLIDE(LL, TANGENT, TEMPFT)
                             CALL CFFCTOW(LL, NORMAL)
                             DO K = 1, DIMN
                                PFN(K,IJK,LL) = PFN(K,IJK,LL) + FNS1(K)
                                PFT(K,IJK,LL) = PFT(K,IJK,LL) + FTS1(K)
                             END DO
                          ELSE
                             IJK = IJK + 1
                          END IF

                       END DO
                    END IF
		    IF(CO.EQ.0) THEN
                       PN(1,LL) = PN(1,LL) + 1
                       IJK = PN(1,LL) + 1
                       PN(IJK,LL) = I
                       PV(IJK,LL) = 1
                       CALL CFNORMAL(LL, I, NORMAL)
!                       CALL CFTANGENT(TANGENT, NORMAL)
                       CALL CFRELVEL(LL, I, VRE, TANGENT)
                       CALL CFTANGENT(TANGENT, NORMAL, VRE)
                       CALL CFVRN(Vn, VRE, NORMAL)
                       CALL CFVRT(Vt, VRE, TANGENT)
                       CALL CFTOTALOVERLAPS(LL, I, Vt, OVERLAP_N, OVERLAP_T)
                       CALL CFFN(LL, Vn, OVERLAP_N, NORMAL)
                       CALL CFFT(LL, Vt, OVERLAP_T, TANGENT)
                       DO K = 1, DIMN
                          TEMPFT(K) = FT(K,LL)
                       END DO
                       CALL CFSLIDE(LL, TANGENT, TEMPFT)
                       CALL CFFCTOW(LL, NORMAL)
                       DO K = 1, DIMN
                          PFN(K,IJK,LL) = FNS1(K)
                          PFT(K,IJK,LL) = FTS1(K)
                       END DO
		    END IF
                 END DO
              END IF

              IF((NEIGHBOR.EQ.0).AND.(WALLCHECK.EQ.0)) THEN
                 CALL CFNOCONTACT(LL)
              END IF

           END DO

           IF(DES_CONTINUUM_COUPLED) THEN
              CALL DRAG_FGS(PARTICLES)
           END IF


!-------------------------------------------------------------------
!     Update old values with new values
!-------------------------------------------------------------------
           CALL CFUPDATEOLD(PARTICLES)

           RETURN
           END SUBROUTINE PERIODIC_WALL_CALC_FORCE_DES

