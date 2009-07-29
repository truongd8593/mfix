!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_FORCE_DES(C)                                      C
!>
!!  Purpose: DES calculations of force acting on a particle, 
!!           its velocity and its position                  
!<
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer: Sreekanth Pannala                        Date: 06-Dec-06  C
!  Reviewer: Rahul Garg                               Date: 02-Aug-07  C
!  Comments: Now includes particle-wall interaction history.           C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_FORCE_DES

      USE run
      USE param1
      USE discretelement
      USE geometry
      USE compar
      USE constant
      IMPLICIT NONE
      
      INTEGER LL, I, K, J, II, IW, FOCUS_PART2, NEIGH_L
      INTEGER NI, NLIM, N_NOCON 
      INTEGER OVERLAP_MAXP
      INTEGER WALLCONTACT, WALLCHECK

      DOUBLE PRECISION OVERLAP_N, OVERLAP_T, TEMPX, TEMPY, TEMPZ, TEMPD
      DOUBLE PRECISION V_REL_TRANS_NORM, V_REL_TRANS_TANG 
      DOUBLE PRECISION FTMD, FNMD
      DOUBLE PRECISION TEMPFN(DIMN), TEMPFT(DIMN), DIST(DIMN), R_LM
      DOUBLE PRECISION NORMAL(DIMN), TANGENT(DIMN), DISTMOD

      DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT 
      DOUBLE PRECISION, SAVE :: FTHIST(2)
      INTEGER, SAVE ::  CONTACT_COUNT = 1

      LOGICAL ALREADY_EXISTS
      LOGICAL CHECK_CON, ALREADY_NEIGHBOURS, OVERLAP_MAX_WALL

!     
!---------------------------------------------------------------------
!     Calculate new values
!---------------------------------------------------------------------
      OVERLAP_MAXP = UNDEFINED_I
      OVERLAP_MAX_WALL = .FALSE.
     
      OVERLAP_MAX = ZERO
      NEIGH_MAX = -1
      FOCUS_PARTICLE = 0
      FOCUS_PART2 = 200000000
      IF (S_TIME.LE.DTSOLID) THEN
         TANGENT(:) = ZERO
         NORMAL(:) = ZERO
         FC(:,:) = ZERO
         FN(:,:) = ZERO
         FT(:,:) = ZERO
      ENDIF

!---------------------------------------------------------------------
!     Calculate contact force and torque
!---------------------------------------------------------------------
     
      DO LL = 1, PARTICLES
         
         IF(DEBUG_DES.AND.LL.EQ.FOCUS_PARTICLE) THEN 
            PRINT*, DES_POS_NEW(LL,1), DES_POS_NEW(LL,2)
            PRINT*, DES_VEL_NEW(LL,1), DES_VEL_NEW(LL,2)
            PRINT*,'PFN', PFN(LL,NI,1:2)
         ENDIF

         NI = 0
         TEMPFN(:) = ZERO
         TEMPFT(:) = ZERO
         K = 0
         IF(PN(LL,1).GE.1) THEN
            NLIM = PN(LL,1)+1
            N_NOCON = 0
            DO NI = 2, PN(LL,1)+1
               IF(PV(LL,NI-N_NOCON).EQ.0.AND.NI-N_NOCON.GE.2) THEN
                  PN(LL,NI-N_NOCON:NLIM-1) = PN(LL,NI-N_NOCON+1:NLIM) 
                  PV(LL,NI-N_NOCON:NLIM-1) = PV(LL,NI-N_NOCON+1:NLIM) 
                  PFN(LL,NI-N_NOCON:NLIM-1,:) = PFN(LL,NI-N_NOCON+1:NLIM,:) 
                  PFT(LL,NI-N_NOCON:NLIM-1,:) = PFT(LL,NI-N_NOCON+1:NLIM,:) 
                  N_NOCON = N_NOCON + 1
                  PN(LL,1) = PN(LL,1) - 1
                  NLIM = PN(LL,1) + 1
               ENDIF
            END DO
         ENDIF
         
!     Initializing rest of the neighbor list which is not in contact

         NLIM = MAX(2,PN(LL,1) + 2) 
         PN(LL,NLIM:MAXNEIGHBORS) = -1
         PFN(LL,NLIM:MAXNEIGHBORS,:) = ZERO
         PFT(LL,NLIM:MAXNEIGHBORS,:) = ZERO
         
         NEIGH_MAX = MAX(NEIGH_MAX, PN(LL,1)+1)
!     Initializing the neighbor list contact information when particles are not in contact

         IF (PN(LL,1).EQ.0) THEN
            PFN(LL,:,:) = ZERO
            PFT(LL,:,:) = ZERO
         END IF
         
!     Initializing the particle
         DO K = 2, MAXNEIGHBORS
            PV(LL,K) = 0
         END DO
         
!     Treats wall interaction also as a two-particle interaction but accounting for the wall properties
         IF(WALLDTSPLIT) THEN
            WALLCHECK = 0
            DO IW = 1, NWALLS
               WALLCONTACT = 0
!     Check to see if a particle is in contact with any of the walls
               CALL CFWALLCONTACT(IW, LL, WALLCONTACT)
               IF(WALLCONTACT.EQ.1) THEN
                  WALLCHECK = 1
                  IF(DEBUG_DES) PRINT*, 'WAll = ', IW, LL
                  I = PARTICLES + IW
                  
                  ALREADY_NEIGHBOURS=.FALSE.
                  
                  IF(PN(LL,1).GT.0) THEN
                     
                     DO NEIGH_L = 2, PN(LL,1)+1
                        IF(I.EQ. PN(LL,NEIGH_L)) THEN 
                           ALREADY_NEIGHBOURS=.TRUE.
                           NI = NEIGH_L
                           EXIT
                        end IF
                        
                     ENDDO
                  end IF
                  
!     Assign the wall particle a position and velocity
                  CALL CFWALLPOSVEL(LL, IW)
                  DES_POS_NEW(I,:) = DES_WALL_POS(IW,:)
                  DES_VEL_NEW(I,:) = DES_WALL_VEL(IW,:)
                  OMEGA_NEW(I,:) = ZERO
                  DES_RADIUS(I) = DES_RADIUS(LL)
                  R_LM = DES_RADIUS(LL) + DES_RADIUS(I)
                  DIST(:) = DES_POS_NEW(I,:) - DES_POS_NEW(LL,:)
                  DISTMOD = SQRT(DES_DOTPRDCT(DIST,DIST))

                  IF(R_LM - DISTMOD.gt.SMALL_NUMBER) then 
                     
                     IF((((R_LM-DISTMOD)/R_LM)*100.d0).GT.OVERLAP_MAX) THen
                        OVERLAP_MAX = (((R_LM-DISTMOD)/R_LM)*100.d0)
                        OVERLAP_MAXP = LL
                        OVERLAP_MAX_WALL = .TRUE.
                     ENDIF
                     CHECk_CON = .TRUE.
                     IF(DISTMOD.NE.ZERO) THEN
                        NORMAL(:)= DIST(:)/DISTMOD
                     ELSE 
                        PRINT *,'DISTMOD IS ZERO', I,LL
                        STOP
                     ENDIF

!     Calculate the translational relative velocity for a contacting particle pair
                     CALL CFRELVEL(LL, I, V_REL_TRANS_NORM, V_REL_TRANS_TANG, TANGENT, NORMAL)
                     
                     IF(ALREADY_NEIGHBOURS) THEN 
                        PV(LL,NI) = 1
                        OVERLAP_N = V_REL_TRANS_NORM*DTSOLID
                        OVERLAP_T = V_REL_TRANS_TANG*DTSOLID
                        contact_count = contact_count + 1
                     ELSE 
                        PN(LL,1) = PN(LL,1) + 1
                        NI = PN(LL,1) + 1
                        PN(LL,NI) = I
                        PV(LL,NI) = 1
                        OVERLAP_N = V_REL_TRANS_NORM*DTSOLID ! R_LM-DISTMOD !
                        OVERLAP_T = V_REL_TRANS_TANG*DTSOLID
                        contact_count = 1
                        FTHIST = ZERO
                     ENDIF
                  ELSE
                     CHECk_CON = .FALSE.
                     GOTO 200
                  ENDIF
                  
                  
                  ETA_N_W = DES_ETAN_WALL(PIJK(LL,5))
                  ETA_T_W = DES_ETAT_WALL(PIJK(LL,5))
                  
                  
                  FNS1(:) = -KN_W*((OVERLAP_N))*NORMAL(:)
                  FNS2(:) = -ETA_N_W*V_REL_TRANS_NORM*NORMAL(:)
                  
                  FTS1(:) = -KT_W*((OVERLAP_T)) *TANGENT(:)
                  FTS2(:) = -ETA_T_W*V_REL_TRANS_TANG*TANGENT(:)
                  
                  
                  FT(LL,:) = FTS1(:) + FTS2(:) 
                  FN(LL,:) = FNS1(:) + FNS2(:) 
                  
                  FN(LL, :) = FN(LL,:) + PFN(LL,NI,:)
                  TEMPFT(:) = FT(LL, :) + PFT(LL,NI,:)
                  
!     Check for Coulomb’s friction law and limit the maximum value of the tangential force on a particle in contact with a wall
                  CALL CFSLIDEWALL(LL, TANGENT, TEMPFT)
                  
!     Calculate the total force Fc and Tow on a particle in a particle-wall collision
                  CALL CFFCTOWALL(LL, NORMAL)

                  PFN(LL,NI,:) =  PFN(LL,NI,:) + FNS1(:)
                  IF(DEBUG_DES.AND.LL.EQ.FOCUS_PARTICLE) THEN 
                     PRINT*,'WALL CONTACT ON', NI
                     PRINT*, 'ALREADY_NEIGHBOURS? = ', ALREADY_NEIGHBOURS   
                     PRINT*, 'STIME, DTSOLID = ',S_TIME, DTSOLID
                     PRINT*, 'DES_VEL = ',DES_VEL_NEW(LL,1:DIMN), des_radius(LL)*OMEGA_NEW(LL,1)
                     PRINT*,'MAA'
                     WRITE(*,*) 'V-OMEGA R = ',DES_VEL_NEW(LL,1)+ des_radius(LL)* OMEGA_NEW(LL,1),   &
                     (DES_VEL_NEW(LL,1)+ des_radius(LL)*OMEGA_NEW(LL,1))*DTSOLID

                     PRINT*,'Mg = ', PMASS(LL)*gravity
                     PRINT*,'K_N, ETA_N = ', KN_W, ETA_N_W, KT_W, ETA_T_W
                     PRINT*,'TANGENT= ', TANGENT
                     PRINT*,'HIST = ', PFT(LL,NI,1:2), PFN(LL,NI,1:2)
                     PRINT*, 'PARTICLE_SLIDE ? ', PARTICLE_SLIDE
                     PRINT*,' FT and FN= ', FT( LL,:), FN(LL,:)
                     PRINT*,'KW*OT*TAN = ', KT_W*((OVERLAP_T)) *TANGENT(:)
                     PRINT*,'OVERLAP_T = ', OVERLAP_T, TANGENT
                     FTMD = SQRT(DES_DOTPRDCT(TEMPFT,TEMPFT))
                     FNMD = SQRT(DES_DOTPRDCT(FN(LL,1:DIMN),FN(LL,1:DIMN)))
                     PRINT*,'FTMD, mu FNMD = ', FTMD, MEW_W*FNMD
                     FTHIST(:) = FTHIST(:) + FT(LL,:)
                     PRINT*,'FT AVG = ', FTHIST/contact_count, mew_w*PMASS(LL)*gravity, contact_count
                     READ(*,*)
                  ENDIF

                  
                  IF(.NOT.PARTICLE_SLIDE) THEN
                     PFT(LL,NI,:) = PFT(LL,NI,:) + FTS1(:)
                  ELSE
                     PFT(LL,NI,:) = PFT(LL,NI,:) + FTS1(:)
                     PARTICLE_SLIDE = .FALSE.
                  ENDIF
                  
                  

               END IF           !Wall Contact
 200           CONTINUE
            ENDDO
         ENDIF                 !if(walldtsplit)
         
         
         IF (NEIGHBOURS(LL,1).GT.0) THEN
            DO II = 2, NEIGHBOURS(LL,1)+1
               I = NEIGHBOURS(LL,II)

               IF(I.GT.LL) THEN
                  
                  ALREADY_NEIGHBOURS=.FALSE.
                  
                  IF(PN(LL,1).GT.0) THEN
                     
                     DO NEIGH_L = 2, PN(LL,1)+1
                        IF(I.EQ. PN(LL,NEIGH_L)) THEN 
                           ALREADY_NEIGHBOURS=.TRUE.
                           NI = NEIGH_L
                           EXIT
                        ENDIF
                     ENDDO
                  ENDIF
                  
                  IF(DES_PERIODIC_WALLS) THEN
                     TEMPX = DES_POS_NEW(I,1)
                     TEMPY = DES_POS_NEW(I,2)
                     IF(DIMN.EQ.3) TEMPZ = DES_POS_NEW(I,3)        
                     TEMPD = ABS(DES_POS_NEW(LL,1) - DES_POS_NEW(I,1))
                     IF(TEMPD.GT.4.d0*MAX_RADIUS.AND.DES_PERIODIC_WALLS_X) THEN 
                        IF(DEBUG_DES) THEN
                           PRINT*,'II, L = ', I,LL
                           PRINT*,'OLD POS, I', DES_POS_NEW(I,:)
                           PRINT*,'OLD POS, LL = ', DES_POS_NEW(LL,:)
                        ENDIF
                        IF(TEMPX.GT.DES_POS_NEW(LL,1)) THEN 
                           DES_POS_NEW(I,1) = DES_POS_NEW(I,1) - (EX2-WX1)
                           IF(DEBUG_DES) PRINT*,'NEW POS WEST= ', DES_POS_NEW(I,1)
                        ELSE
                           DES_POS_NEW(I,1) = DES_POS_NEW(I,1) + EX2 - WX1
                           
                           IF(DEBUG_DES) PRINT*,'NEW POS EAST = ', DES_POS_NEW(I,1)
                        ENDIF
                     ENDIF
                     
                     TEMPD = ABS(DES_POS_NEW(LL,2) - DES_POS_NEW(I,2))
                     IF(TEMPD.GT.4.d0*MAX_RADIUS.AND.DES_PERIODIC_WALLS_Y) THEN
                        IF(TEMPY.GT.DES_POS_NEW(LL,2)) THEN 
                           DES_POS_NEW(I,2) = DES_POS_NEW(I,2) - (TY2-BY1)
                        ELSE
                           DES_POS_NEW(I,2) = DES_POS_NEW(I,2) + (TY2-BY1)
                        ENDIF
                     ENDIF
                     
                     IF(DIMN.EQ.3) THEN
                        TEMPD = ABS(DES_POS_NEW(LL,3) - DES_POS_NEW(I,3))
                        IF(TEMPD.GT.4.d0*MAX_RADIUS.AND.DES_PERIODIC_WALLS_Z) THEN 
                           IF(TEMPZ.GT.DES_POS_NEW(LL,3)) THEN 
                              DES_POS_NEW(I,3) = DES_POS_NEW(I,3) -(NZ2 - SZ1)
                           ELSE
                              DES_POS_NEW(I,3) = DES_POS_NEW(I,3) + (NZ2-SZ1)
                           ENDIF
                        ENDIF
                     ENDIF
                  ENDIF
                  
                  R_LM = DES_RADIUS(LL) + DES_RADIUS(I)
                  DIST(:) = DES_POS_NEW(I,:) - DES_POS_NEW(LL,:)
                  DISTMOD = SQRT(DES_DOTPRDCT(DIST,DIST))
                  
                  IF(DES_PERIODIC_WALLS) THEN
                     DES_POS_NEW(I,1) = TEMPX
                     DES_POS_NEW(I,2) = TEMPY
                     IF (DIMN.EQ.3) DES_POS_NEW(I,3) = TEMPZ  
                  ENDIF

                  IF(R_LM - DISTMOD.gt.SMALL_NUMBER) then 
                     IF(LL.EQ.FOCUS_PARTICLE) Print*, 'NEIGHBORS', NEIGHBOURS(LL,:)

                     
                     IF((((R_LM-DISTMOD)/R_LM)*100.d0).GT.OVERLAP_MAX) THEN
                        OVERLAP_MAX = (((R_LM-DISTMOD)/R_LM)*100.d0)
                        OVERLAP_MAXP = LL
                        OVERLAP_MAX_WALL = .FALSE.
                     ENDIF

                     CHECk_CON = .TRUE.
                     IF(DISTMOD.NE.ZERO) THEN
                        NORMAL(:)= DIST(:)/DISTMOD
                     ELSE 
                        PRINT *,'DISTMOD IS ZERO FOR PART. PAIR', I,LL
                        STOP
                     ENDIF
!     Calculate the translational relative velocity for a contacting particle pair.
                     CALL CFRELVEL(LL, I, V_REL_TRANS_NORM, V_REL_TRANS_TANG, TANGENT, NORMAL)
                     IF(ALREADY_NEIGHBOURS) THEN 
                        PV(LL,NI) = 1
                        OVERLAP_N = V_REL_TRANS_NORM*DTSOLID
                        OVERLAP_T = V_REL_TRANS_TANG*DTSOLID
                     ELSE 
                        IF(DEBUG_DES) THEN
                           PRINT*,'INITIATING CONTACT'
                           PRINT*,'OVERLAP = ', ((R_LM-DISTMOD)/R_LM)*100.d0
                        ENDIF
                        PN(LL,1) = PN(LL,1) + 1
                        NI = PN(LL,1) + 1
                        PN(LL,NI) = I
                        PV(LL,NI) = 1
                        OVERLAP_N = R_LM-DISTMOD !V_REL_TRANS_NORM*DTSOLID! 
                        OVERLAP_T = V_REL_TRANS_TANG*DTSOLID
                     END IF
                  ELSE
                     GOTO 300
                  ENDIF
                  
                  ETA_DES_N = DES_ETAN(PIJK(LL,5), PIJK(I,5))
                  ETA_DES_T = DES_ETAT(PIJK(LL,5), PIJK(I,5))
                  FNS1(:) = -KN*((OVERLAP_N))*NORMAL(:)
                  FNS2(:) = -ETA_DES_N*V_REL_TRANS_NORM*NORMAL(:)
                  
                  FTS1(:) = -KT*((OVERLAP_T)) *TANGENT(:)
                  FTS2(:) = -ETA_DES_T*V_REL_TRANS_TANG*TANGENT(:)
                  
                  FT(LL,:) = FTS1(:) + FTS2(:) 
                  FN(LL,:) = FNS1(:) + FNS2(:) 
                  
                  FN(LL, :) = FN(LL,:) + PFN(LL,NI,:)
                  TEMPFT(:) = FT(LL, :) + PFT(LL,NI,:)
                  
                  IF(DEBUG_DES.AND.LL.EQ.FOCUS_PARTICLE) THEN 
                     PRINT*, 'I = ', I
                     PRINT*,'ETAs =  ', ETA_DES_N, ETA_DES_T
                     PRINT*,'overlap = ', overlap_n, (R_LM - DISTMOD)*100.d0/R_LM
                     
                     PRINT*,'rad ratio = ', DES_RADIUS(LL)/DES_RADIUS(I)
                     PRINT*, 'FNS1 and FNS2 = ', FNS1(:), FNS2(:)
                     PRINT*, 'PFN = ', PFN(LL,NI,:)
                     PRINT*, 'PFT = ', PFT(LL,NI,:)
                     PRINT*, 'FORCEST = ', FT(LL,:)
                     PRINT*, 'FORCESN = ', FN(LL,:)
                     PRINT*, 'FORCEST = ', FT(LL,:)
                  ENDIF

!     Check for Coulomb’s friction law and limit the maximum value of the tangential force on a particle in contact with another particle
                  CALL CFSLIDE(LL, TANGENT, TEMPFT)
                  
                  PFN(LL,NI,:) = PFN(LL,NI,:) +  FNS1(:)

!     Calculate the total force Fc and Tow on a particle in a particle-particle collision
                  CALL CFFCTOW(LL, I, NORMAL)

                  IF(.NOT.PARTICLE_SLIDE) THEN
                     PFT(LL,NI,:) = PFT(LL,NI,:) + FTS1(:)
                  ELSE
                     PFT(LL,NI,:) = PFT(LL,NI,:) + FTS1(:)
                     PARTICLE_SLIDE = .FALSE.
                  ENDIF
                  
!     !impulse is effectively doubled for wall interactions

                  IF(DEBUG_DES.AND.LL.eq.FOCUS_PARTICLE)THEN
                     INQUIRE(FILE='debug_file',EXIST=ALREADY_EXISTS)
                     IF(ALREADY_EXISTS)THEN
                        OPEN(UNIT=1,FILE='debug_file',STATUS='OLD',POSITION='APPEND')
                        WRITE(1,'(A,I5)')'CALC FORCE -- NEIGHBOR',II
                        WRITE(1,'(2(1x,A,E12.5))')&
                        'FNx=',FN(LL,1)&
                        ,'FNy=',FN(LL,2)
                     ELSE
                        OPEN(UNIT=1,FILE='debug_file',STATUS='NEW')
                        WRITE(1,'(A,I5)')'CALC FORCE -- NEIGHBOR',II
                        WRITE(1,'(2(1x,A,E12.5))')&
                        'FNx=',FN(LL,1)&
                        ,'FNy=',FN(LL,2)
                     ENDIF
                     CLOSE (1)
                     Print*, 'PN', PN(LL,:)
                  ENDIF
!--   END DEBUGGING
                  
               ENDIF

               
 300           CONTINUE
               
               
            ENDDO              !II = 2, NEIGHBOURS(LL,1)+I
         ENDIF                 !(NEIGHBOURS(LL,1).GT.0)
         
         IF((NEIGHBOURS(LL,1).EQ.0).AND.(WALLCHECK.EQ.0)) THEN
!     The subroutine sets all forces on a particle to zero is the particle is found to have no neighbors (either particles or walls)
            CALL CFNOCONTACT(LL)
         ENDIF

      ENDDO

      IF(DES_CONTINUUM_COUPLED) THEN
!     Treats wall interaction also as a two-particle interaction but accounting for the wall properties
         CALL DRAG_FGS
      ENDIF
!     COHESION
      IF(USE_COHESION)THEN
         CALL CALC_COHESIVE_FORCES
      ENDIF
!     !COHESION
!-------------------------------------------------------------------
!     Update old values with new values
!-------------------------------------------------------------------
      
!     Update the old values of particle position and velocity with the new values computed
      CALL CFUPDATEOLD
      RETURN
      END SUBROUTINE CALC_FORCE_DES

