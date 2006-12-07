!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_FORCE_DES(C)                                      C
!  Purpose: DES calculations of force acting on a particle,            C
!           its velocity and its position                              C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer: Sreekanth Pannala                        Date: 06-Dec-06  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_FORCE_DES

      USE run
      USE param1
      USE discretelement
      USE geometry
      USE compar

      IMPLICIT NONE
      
      INTEGER LL, I, II, K, LC, C, CO, IW, KK, TEMP, TEMPN
      INTEGER NI, IJ, WALLCHECK, NLIM, N_NOCON
      DOUBLE PRECISION OVERLAP_N, OVERLAP_T
      DOUBLE PRECISION V_REL_TRANS(DIMN), V_SLIP(DIMN)
      DOUBLE PRECISION NORMAL(DIMN), TANGENT(DIMN)
      DOUBLE PRECISION V_REL_TRANS_NORM, V_REL_TRANS_TANG
      DOUBLE PRECISION TEMPFN(DIMN), TEMPFT(DIMN)
      LOGICAL ALREADY_EXISTS
!     
!---------------------------------------------------------------------
!     Calculate new values
!---------------------------------------------------------------------
!     

      IF (S_TIME.LE.DTSOLID) THEN
         V_REL_TRANS(:) = ZERO
         TANGENT(:) = ZERO
         NORMAL(:) = ZERO
         FC(:,:) = ZERO
         FN(:,:) = ZERO
         FT(:,:) = ZERO
!        PN(:,1) = 0
!        PN(:,2:MAXNEIGHBORS) = -1
!        PV(:,:) = 1
         PFN(:,:,:) = ZERO
         PFT(:,:,:) = ZERO
      END IF

      IF(DES_CONTINUUM_COUPLED) THEN
         CALL PARTICLES_IN_CELL
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

!     Compact the neighboring particles by eliminating the particles out of contact (PV(LL,K)=0).

         IF(PN(LL,1).GE.1) THEN
            NLIM = PN(LL,1)+1
            N_NOCON = 0
            DO NI = 2, PN(LL,1)+1
               IF(PV(LL,NI-N_NOCON).EQ.0.AND.NI-N_NOCON.GE.2) THEN
                  PN(LL,NI-N_NOCON:NLIM-1) = PN(LL,NI-N_NOCON+1:NLIM) 
                  PV(LL,NI-N_NOCON:NLIM-1) = PV(LL,NI-N_NOCON+1:NLIM) 
                  PFN(LL,NI-N_NOCON:NLIM-1,:) = PFN(LL,NI-N_NOCON+1:NLIM,:) 
                  PFT(LL,NI-N_NOCON:NLIM-1,:) = PFN(LL,NI-N_NOCON+1:NLIM,:) 
                  N_NOCON = N_NOCON + 1
                  PN(LL,1) = PN(LL,1) - 1
                  NLIM = PN(LL,1) + 1
               ENDIF
            END DO
         ENDIF

!     Initializing rest of the neighbor list which is not in contact

         NLIM = MAX(2,PN(LL,1) + 1) 
         PN(LL,NLIM:MAXNEIGHBORS) = -1
         PFN(LL,NLIM:MAXNEIGHBORS,:) = ZERO
         PFT(LL,NLIM:MAXNEIGHBORS,:) = ZERO

!     Initializing the neighbor list contact information when particles are not in contact

         IF (PN(LL,1).EQ.0) THEN
           PFN(LL,:,:) = ZERO
           PFT(LL,:,:) = ZERO
         END IF

!     Initializing the particle
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
                  CALL CFNORMALWALL(LL, I, NORMAL)
                  CALL CFRELVEL(LL, I, V_REL_TRANS)
                  CALL CFVRN(V_REL_TRANS_NORM, V_REL_TRANS, NORMAL)
                  CALL CFSLIPVEL(LL, I, V_SLIP, V_REL_TRANS, V_REL_TRANS_NORM, NORMAL)
                  CALL CFTANGENT(V_SLIP, TANGENT, NORMAL)
                  CALL CFVRT(V_REL_TRANS_TANG, V_REL_TRANS, TANGENT)
                  CALL CFTOTALOVERLAPSWALL(LL, I, V_REL_TRANS_TANG, OVERLAP_N, OVERLAP_T)
                  CALL CFFNWALL(LL, V_REL_TRANS_NORM, OVERLAP_N, NORMAL)
                  CALL CFFTWALL(LL, V_REL_TRANS_TANG, OVERLAP_T, TANGENT)
                  CALL CFSLIDEWALL(LL, TANGENT)
                  CALL CFFCTOWALL(LL, NORMAL)
!--   DEBUGGING
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
!--   END DEBUGGING
               END IF
            END DO
         END IF

         IF (NEIGHBOURS(LL,1).GT.0) THEN
            DO II = 2, NEIGHBOURS(LL,1)+1
               I = NEIGHBOURS(LL,II)
               IF(I.GT.LL) THEN
                  CO = 0 ! already in contact - 1, new contact - 0
                  NI = 2

                  IF(PN(LL,1).GT.0) THEN ! check to see whether already in contact
                     IF((CO.EQ.0).AND.(NI.LE.(PN(LL,1)+1))) THEN
 20                     CONTINUE
                        IF(I.EQ.PN(LL,NI)) THEN
                           CO = 1
                           PV(LL,NI) = 1
                           CALL CFNORMAL(LL, I, II, NORMAL)
                           CALL CFRELVEL(LL, I, V_REL_TRANS)
                           CALL CFVRN(V_REL_TRANS_NORM, V_REL_TRANS, NORMAL)
                           CALL CFSLIPVEL(LL, I, V_SLIP, V_REL_TRANS, V_REL_TRANS_NORM, NORMAL)
                           CALL CFTANGENT(V_SLIP, TANGENT, NORMAL)
                           CALL CFVRT(V_REL_TRANS_TANG, V_REL_TRANS, TANGENT)
                           CALL CFINCREMENTALOVERLAPS(V_REL_TRANS_NORM, V_REL_TRANS_TANG, OVERLAP_N, OVERLAP_T)
                           CALL CFFN(LL, V_REL_TRANS_NORM, OVERLAP_N, NORMAL)
                           CALL CFFT(LL, V_REL_TRANS_TANG, OVERLAP_T, TANGENT)
                           FN(LL,:) = FN(LL,:) + PFN(LL,NI,:)
                           TEMPFT(:) = FT(LL,:) + PFT(LL,NI,:)
                           CALL CFSLIDE(LL, TANGENT, TEMPFT)
                           CALL CFFCTOW(LL, I, NORMAL)
                           PFN(LL,NI,:) = PFN(LL,NI,:) + FNS1(:)
                           IF(.NOT.PARTICLE_SLIDE) THEN
                              PFT(LL,NI,:) = PFT(LL,NI,:) + FTS1(:)
                           ELSE
                              PFT(LL,NI,:) = FT(LL,:)
                              PARTICLE_SLIDE = .FALSE.
                           END IF
                        ELSE
                           NI = NI + 1
                        END IF
                        IF((CO.EQ.0).AND.(NI.LE.(PN(LL,1)+1))) GO TO 20
                     END IF
                  END IF

                  IF(CO.EQ.0) THEN ! New contact
                     PN(LL,1) = PN(LL,1) + 1
                     NI = PN(LL,1) + 1
                     PN(LL,NI) = I
                     PV(LL,NI) = 1
                     CALL CFNORMAL(LL, I, II, NORMAL)
                     CALL CFRELVEL(LL, I, V_REL_TRANS)
                     CALL CFVRN(V_REL_TRANS_NORM, V_REL_TRANS, NORMAL)
                     CALL CFSLIPVEL(LL, I, V_SLIP, V_REL_TRANS, V_REL_TRANS_NORM, NORMAL)
                     CALL CFTANGENT(V_SLIP, TANGENT, NORMAL)
                     CALL CFVRT(V_REL_TRANS_TANG, V_REL_TRANS, TANGENT)
                     CALL CFTOTALOVERLAPS(LL, I, II, V_REL_TRANS_TANG, OVERLAP_N, OVERLAP_T)
                     CALL CFFN(LL, V_REL_TRANS_NORM, OVERLAP_N, NORMAL)
                     CALL CFFT(LL, V_REL_TRANS_TANG, OVERLAP_T, TANGENT)
                     TEMPFT(:) = FT(LL,:)
                     CALL CFSLIDE(LL, TANGENT, TEMPFT)
                     CALL CFFCTOW(LL, I, NORMAL)
                     PFN(LL,NI,:) = FNS1(:)
                     IF(.NOT.PARTICLE_SLIDE) THEN
                        PFT(LL,NI,:) =  FTS1(:)
                     ELSE
                        PFT(LL,NI,:) = FT(LL,:)
                        PARTICLE_SLIDE = .FALSE.
                     END IF
                  END IF

!--   DEBUGGING
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
!--   END DEBUGGING
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

!     !COHESION
      IF(USE_COHESION)THEN
         CALL CALC_COHESIVE_FORCES
      END IF
!     !COHESION

!-------------------------------------------------------------------
!     Update old values with new values
!-------------------------------------------------------------------
      CALL CFUPDATEOLD(PARTICLES)

      RETURN
      END SUBROUTINE CALC_FORCE_DES

