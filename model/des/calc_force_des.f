!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_FORCE_DES(C)                                      C
!  Purpose: DES calculations of force acting on a particle,            C
!           its velocity and its position                              C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_FORCE_DES

      USE run
      USE param1
      USE discretelement
      USE geometry

      IMPLICIT NONE
      
      INTEGER LL, I, II, K, LC, C, CO, IW, KK, TEMP, TEMPN
      DOUBLE PRECISION OVERLAP_N, OVERLAP_T, TEMPFN(DIMN), TEMPFT(DIMN)
      DOUBLE PRECISION NORMAL(DIMN), VRE(DIMN)
      DOUBLE PRECISION TANGENT(DIMN)
      DOUBLE PRECISION Vn, Vt
      INTEGER NI, IJ, WALLCHECK
      LOGICAL ALREADY_EXISTS
!     
!---------------------------------------------------------------------
!     Calculate new values
!---------------------------------------------------------------------
!     

      IF (S_TIME.LE.DTSOLID) THEN
            VRE(:) = ZERO
            TANGENT(:) = ZERO
            NORMAL(:) = ZERO
         DO LC = 1, PARTICLES
               FC(LC,:) = ZERO
               FN(LC,:) = ZERO
               FT(LC,:) = ZERO
            DO KK = 1, MN
               PN(LC,KK) = -1
               PV(LC,KK) = 1
               PFN(LC,KK,:) = ZERO
               PFT(LC,KK,:) = ZERO
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
                       PFT(LL,NI+1,:) = PFT(LL,NI,:)
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
                   IF(PN(LL,NI).GE.2*PARTICLES) THEN
                       PN(LL,NI) = - 1
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
                 DO IW = 1, NWALLS
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
                       CALL CFVRN(Vn, VRE, NORMAL)
                       CALL CFVRT(Vt, VRE, TANGENT)
                       CALL CFTOTALOVERLAPS(LL, I, Vt, OVERLAP_N, OVERLAP_T)
                       CALL CFFNWALL(LL, Vn, OVERLAP_N, NORMAL)
                       CALL CFFTWALL(LL, Vt, OVERLAP_T, TANGENT)
                       CALL CFSLIDEWALL(LL, TANGENT)
                       CALL CFFCTOWALL(LL, NORMAL)
!--DEBUGGING
      !!impulse is effectively doubled for wall interactions
      IF(LL.eq.FOCUS_PARTICLE)THEN
      INQUIRE(FILE='debug_file',EXIST=ALREADY_EXISTS)
      IF(ALREADY_EXISTS)THEN
        OPEN(UNIT=1,FILE='debug_file',STATUS='OLD',POSITION='APPEND')
        WRITE(1,'(A)')'CALC_FORCE-WALL'
      ELSE
         OPEN(UNIT=1,FILE='debug_file',STATUS='NEW')
        WRITE(1,'(A)')'CALC_FORCE-WALL'
      END IF
      CLOSE (1)
      END IF
!--END DEBUGGING
                    END IF
                 END DO
              END IF

              IF (NEIGHBOURS(LL,1).GT.0) THEN
                 DO II = 2, NEIGHBOURS(LL,1)+1
                   I = NEIGHBOURS(LL,II)
                   IF(I.GT.LL) THEN
	      	    CO = 0
		    NI = 2
                    IF(PN(LL,1).GT.0) THEN
                       DO WHILE((CO.EQ.0).AND.(NI.LE.(PN(LL,1)+1)))
                          IF(I.EQ.PN(LL,NI)) THEN
                             CO = 1
                             PV(LL,NI) = 1
                             CALL CFNORMAL(LL, I, NORMAL)
                             CALL CFTANGENT(TANGENT, NORMAL, VRE)
                             CALL CFRELVEL(LL, I, VRE, TANGENT)
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

!--DEBUGGING
      !!impulse is effectively doubled for wall interactions
      IF(LL.eq.FOCUS_PARTICLE)THEN
      INQUIRE(FILE='debug_file',EXIST=ALREADY_EXISTS)
      IF(ALREADY_EXISTS)THEN
        OPEN(UNIT=1,FILE='debug_file',STATUS='OLD',POSITION='APPEND')
        WRITE(1,'(A,I5)')'CALC FORCE -- NEIGHBOR',II
        WRITE(1,'(1X,L5)')DES_CONTINUUM_COUPLED
        WRITE(1,'(2(1x,A,E12.5))')&
            'FNx=',FN(LL,1)&
            ,'FNy=',FN(LL,2)
      ELSE
        OPEN(UNIT=1,FILE='debug_file',STATUS='NEW')
        WRITE(1,'(A,I5)')'CALC FORCE -- NEIGHBOR',II
        WRITE(1,'(1X,L5)')DES_CONTINUUM_COUPLED
        WRITE(1,'(2(1x,A,E12.5))')&
            'FNx=',FN(LL,1)&
            ,'FNy=',FN(LL,2)
      END IF
      CLOSE (1)
      END IF
!--END DEBUGGING
                  END IF
                 END DO
              END IF
              
              IF((NEIGHBOURS(LL,1).EQ.0).AND.(WALLCHECK.EQ.0)) THEN
                 CALL CFNOCONTACT(LL)
              END IF

           END DO
           
           IF(DES_CONTINUUM_COUPLED) THEN
              CALL DRAG_FGS(PARTICLES)
           END IF

!!COHESION
        IF(USE_COHESION)THEN
          CALL CALC_COHESIVE_FORCES
        END IF
!!COHESION

!-------------------------------------------------------------------
!     Update old values with new values
!-------------------------------------------------------------------
        CALL CFUPDATEOLD(PARTICLES)

           RETURN
           END SUBROUTINE CALC_FORCE_DES

