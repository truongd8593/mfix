!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_FORCE_DES(C)                                      C
!  Purpose: DES calculations of force acting on a particle,            C
!           its velocity and its position                              C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer: Sreekanth Pannala                        Date: 06-Dec-06  C
!  Reviewer: Rahul Garg                               Date: 02-Aug-07  C
!  Comments: Now includes particle-wall interaction history.  
!  Reviewer: Tingwen Li                               Date: 18-Jan-08  C
!  Comments: Non-rectangular wall boundary conditions are implemented. C
!            Only particles in the cells neighboring to walls or       C
!            internal surfaces are checked                             C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_FORCE_DES

      USE run
      USE param1
      USE discretelement
      USE geometry
      USE compar
      
      IMPLICIT NONE
      
      INTEGER LL, I, II, K, LC, C, CO, IW, KK, TEMP, TEMPN, JJ, J, FOCUS_PART2 
      INTEGER NI, IJ, WALLCHECK, NLIM, N_NOCON, IP2
      INTEGER KM1, KP1, IM1, IP1, JM1, JP1, PNO, NPG, PC(3)
      DOUBLE PRECISION OVERLAP_N, OVERLAP_T
      DOUBLE PRECISION V_REL_TRANS(DIMN), V_SLIP(DIMN)
      DOUBLE PRECISION NORMAL(DIMN), TANGENT(DIMN)
      DOUBLE PRECISION V_REL_TRANS_NORM, V_REL_TRANS_TANG, VSLIPMAG
      DOUBLE PRECISION TEMPFN(DIMN), TEMPFT(DIMN), ALPHA, DIST(DIMN), DISTMAG, R_LM
      LOGICAL ALREADY_EXISTS, calc_fc, CALLFROMDES
      DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT 
                  LOGICAL CHECK_CON
      integer tempijk  ! temporary variable 

!     
                  !---------------------------------------------------------------------
!     Calculate new values
!---------------------------------------------------------------------
!     
      FOCUS_PARTICLE = 2000000
      FOCUS_PART2 = 200000000
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
!         PFN(:,:,:) = ZERO
!         PFT(:,:,:) = ZERO
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
!			Compact the neighboring particles by eliminating particles out of contact
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
         

!     Initializing the neighbor list contact information when particles are not in contact

         IF (PN(LL,1).EQ.0) THEN
           PFN(LL,:,:) = ZERO
           PFT(LL,:,:) = ZERO
         END IF
         
!     Initializing the particle
         DO K = 2, MAXNEIGHBORS
            PV(LL,K) = 0
         END DO
         
!     by Tingwen 18/01/2008 10:57:57 AM         
         tempijk=pijk(LL,4)     ! locate the cell that particle is in
!     IF(WALLDTSPLIT) THEN  
         
         IF(WALLDTSPLIT .and. c_near_w(tempijk,7) .ge. 0) THEN
            WALLCHECK = 0
!-------------------------------
            if(NON_RECT_BC)then ! not simple domain shape, rectangular boundary
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
               if(c_near_w(tempijk,7).eq. 0) then ! the cell contact with wall through face
                  do iw=1,nwalls
                     wallcontact=0
                     if(c_near_w(tempijk,iw) .eq. 1)then
                        call wallfacecontact(iw,LL,wallcontact)

                        IF(WALLCONTACT .EQ. 1) THEN ! determine the wall
                           WALLCHECK = 1    
                           NI = 2
                           CO = 0 
                           I = PARTICLES + IW

                           IF(PN(LL,1).GT.0) THEN ! check to see whether already in contact
                              IF((CO.EQ.0).AND.(NI.LE.(PN(LL,1)+1))) THEN
 101                             CONTINUE
                                 IF(I.EQ.PN(LL,NI)) THEN
                                    CO = 1
                                    PV(LL,NI) = 1
                                    
                                    DES_POS_NEW(I,:) = DES_WALL_POS(IW,:)
                                    DES_VEL_NEW(I,:) = DES_WALL_VEL(IW,:)
                                    OMEGA_NEW(I,:) = ZERO
                                    DES_RADIUS(I) = DES_RADIUS(LL)
                                    CALL CFNORMALWALL(LL, I, NORMAL)
                                    CALL CFRELVEL(LL, I, V_REL_TRANS, V_REL_TRANS_NORM, V_REL_TRANS_TANG, TANGENT, NORMAL)
                                    CALL CFINCREMENTALOVERLAPS(LL,I, V_REL_TRANS_NORM,  V_REL_TRANS_TANG, OVERLAP_N, OVERLAP_T, CHECK_CON)
                                    
                                    IF(.not.CHECK_CON) THEN 
                                       PV(LL,NI) = 0 !No More Contact 
                                       GOTO 201
                                    ENDIF
                                    
                                    FNS1(:) = -KN_W*((OVERLAP_N))*NORMAL(:)
                                    FNS2(:) = -ETA_N_W*V_REL_TRANS_NORM*NORMAL(:)
                                    
                                    FTS1(:) = -KT_W*((OVERLAP_T)) *TANGENT(:)
                                    FTS2(:) = - ETA_T_W*V_REL_TRANS_TANG*TANGENT(:)                            
                                    
                                    FT(LL,:) = FTS1(:) + FTS2(:) 
                                    FN(LL,:) = FNS1(:) + FNS2(:) 
                                    
                                    FN(LL, :) = FN(LL,:) + PFN(LL,NI,:)
                                    TEMPFT(:) = FT(LL, :) + PFT(LL,NI,:)
                                    
                                    CALL CFSLIDEWALL(LL, TANGENT, TEMPFT)
                                    
                                    CALL CFFCTOWALL(LL, NORMAL)
                                    
                                    PFN(LL,NI,:) =  PFN(LL,NI,:) + FNS1(:)
                                    
                                    IF(.NOT.PARTICLE_SLIDE) THEN
                                       PFT(LL,NI,:) = PFT(LL,NI,:) + FTS1(:)
                                    ELSE
                                       PFT(LL,NI,:) =  FT(LL,:)
                                       PARTICLE_SLIDE = .FALSE.
                                    END IF                      
                                 ELSE
                                    NI = NI + 1
                                 END IF
                                 
                                 IF((CO.EQ.0).AND.(NI.LE.(PN(LL,1)+1))) GO TO 101
                                 
                              END IF                        
                           END IF ! PN(LL,1).GT.0
                           IF(CO.EQ.0) THEN ! New contact
                              DES_POS_NEW(I,:) = DES_WALL_POS(IW,:)
                              DES_VEL_NEW(I,:) = DES_WALL_VEL(IW,:)
                              OMEGA_NEW(I,:) = ZERO
                              DES_RADIUS(I) = DES_RADIUS(LL)
                              CALL CFNORMALWALL(LL, I, NORMAL)
                              CALL CFRELVEL(LL, I, V_REL_TRANS, V_REL_TRANS_NORM, V_REL_TRANS_TANG, TANGENT, NORMAL)
                              CALL CFTOTALOVERLAPS(LL, I, II, V_REL_TRANS_NORM, V_REL_TRANS_TANG, OVERLAP_N, OVERLAP_T, CHECK_CON)
                              
                              IF(.not.CHECK_CON) THEN 
                                 GOTO 201
                              ENDIF
                              
                              PN(LL,1) = PN(LL,1) + 1
                              NI = PN(LL,1) + 1
                              PN(LL,NI) = I
                              PV(LL,NI) = 1
                              
                              FNS1(:) = -KN_W*((OVERLAP_N))*NORMAL(:)
                              FNS2(:) = -ETA_N_W*V_REL_TRANS_NORM*NORMAL(:)
                              
                              FTS1(:) = -KT_W*((OVERLAP_T)) *TANGENT(:)
                              FTS2(:) = - ETA_T_W*V_REL_TRANS_TANG*TANGENT(:)                        
                              
                              FT(LL,:) = FTS1(:) + FTS2(:) 
                              FN(LL,:) = FNS1(:) + FNS2(:)                         
                              
                              TEMPFT(:) = FT(LL, :)
                              
                              CALL CFSLIDEWALL(LL, TANGENT, TEMPFT)
                              CALL CFFCTOWALL(LL, NORMAL)
                              
                              PFN(LL,NI,:) = FNS1(:)
                              
                              IF(.NOT.PARTICLE_SLIDE) THEN
                                 PFT(LL,NI,:) =  FTS1(:)
                              ELSE
                                 PFT(LL,NI,:) = FT(LL,:)
                                 PARTICLE_SLIDE = .FALSE.
                              END IF                         
                           END IF
                        END IF  !wall contact                     
                     end if     !c_near_w(tmpijk,iw).eq.1
 201                 continue            		  
                  end do        ! iw=1,nwalls    
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
               elseif(c_near_w(tempijk,7) .eq. 1) then ! the cell contact with wall through edge
                  wallcontact=0
                  call walledgecontact(tempijk,LL,wallcontact)

                  IF(WALLCONTACT.EQ.1) THEN
                     WALLCHECK = 1
                     NI = 2
                     CO = 0 
                     I = PARTICLES + nwalls + 1

                     IF(PN(LL,1).GT.0) THEN ! check to see whether already in contact
                        IF((CO.EQ.0).AND.(NI.LE.(PN(LL,1)+1))) THEN
 102                       CONTINUE
                           IF(I.EQ.PN(LL,NI)) THEN
                              CO = 1
                              PV(LL,NI) = 1
                              
                              DES_POS_NEW(I,:) = DES_WALL_POS(1,:)
                              DES_VEL_NEW(I,:) = DES_WALL_VEL(1,:)
                              OMEGA_NEW(I,:) = ZERO
                              DES_RADIUS(I) = DES_RADIUS(LL)
                              CALL CFNORMALWALL(LL, I, NORMAL)
                              CALL CFRELVEL(LL, I, V_REL_TRANS, V_REL_TRANS_NORM, V_REL_TRANS_TANG, TANGENT, NORMAL)
                              CALL CFINCREMENTALOVERLAPS(LL,I, V_REL_TRANS_NORM,  V_REL_TRANS_TANG, OVERLAP_N, OVERLAP_T, CHECK_CON)
                              
                              IF(.not.CHECK_CON) THEN 
                                 PV(LL,NI) = 0 !No More Contact 
                                 GOTO 202
                              ENDIF
                              
                              FNS1(:) = -KN_W*((OVERLAP_N))*NORMAL(:)
                              FNS2(:) = -ETA_N_W*V_REL_TRANS_NORM*NORMAL(:)

                              FTS1(:) = -KT_W*((OVERLAP_T)) *TANGENT(:)
                              FTS2(:) = - ETA_T_W*V_REL_TRANS_TANG*TANGENT(:)                         
                              
                              FT(LL,:) = FTS1(:) + FTS2(:) 
                              FN(LL,:) = FNS1(:) + FNS2(:) 
                              
                              FN(LL, :) = FN(LL,:) + PFN(LL,NI,:)
                              TEMPFT(:) = FT(LL, :) + PFT(LL,NI,:)
                              
                              CALL CFSLIDEWALL(LL, TANGENT, TEMPFT)
                              
                              CALL CFFCTOWALL(LL, NORMAL)

                              PFN(LL,NI,:) =  PFN(LL,NI,:) + FNS1(:)
                              
                              IF(.NOT.PARTICLE_SLIDE) THEN
                                 PFT(LL,NI,:) = PFT(LL,NI,:) + FTS1(:)
                              ELSE
                                 PFT(LL,NI,:) =  FT(LL,:)
                                 PARTICLE_SLIDE = .FALSE.
                              END IF
                           ELSE
                              NI = NI + 1
                           END IF
                           
                           IF((CO.EQ.0).AND.(NI.LE.(PN(LL,1)+1))) GO TO 102
                           
                        END IF  !(CO.EQ.0).AND.(NI.LE.(PN(LL,1)+1))                   
                     END IF     !PN(LL,1).GT.0
                     IF(CO.EQ.0) THEN ! New contact
                        DES_POS_NEW(I,:) = DES_WALL_POS(1,:)
                        DES_VEL_NEW(I,:) = DES_WALL_VEL(1,:)
                        OMEGA_NEW(I,:) = ZERO
                        DES_RADIUS(I) = DES_RADIUS(LL)
                        CALL CFNORMALWALL(LL, I, NORMAL)
                        CALL CFRELVEL(LL, I, V_REL_TRANS, V_REL_TRANS_NORM, V_REL_TRANS_TANG, TANGENT, NORMAL)
                        CALL CFTOTALOVERLAPS(LL, I, II, V_REL_TRANS_NORM, V_REL_TRANS_TANG, OVERLAP_N, OVERLAP_T, CHECK_CON)

                        IF(.not.CHECK_CON) THEN 
                           GOTO 202
                        ENDIF
                        
                        PN(LL,1) = PN(LL,1) + 1
                        NI = PN(LL,1) + 1
                        PN(LL,NI) = I
                        PV(LL,NI) = 1
                        
                        FNS1(:) = -KN_W*((OVERLAP_N))*NORMAL(:)
                        FNS2(:) = -ETA_N_W*V_REL_TRANS_NORM*NORMAL(:)
                        
                        FTS1(:) = -KT_W*((OVERLAP_T)) *TANGENT(:)
                        FTS2(:) = - ETA_T_W*V_REL_TRANS_TANG*TANGENT(:)                    
                        
                        FT(LL,:) = FTS1(:) + FTS2(:) 
                        FN(LL,:) = FNS1(:) + FNS2(:)                     

                        TEMPFT(:) = FT(LL, :)

                        CALL CFSLIDEWALL(LL, TANGENT, TEMPFT)
                        CALL CFFCTOWALL(LL, NORMAL)
                        
                        PFN(LL,NI,:) = FNS1(:)
                        
                        IF(.NOT.PARTICLE_SLIDE) THEN
                           PFT(LL,NI,:) =  FTS1(:)
                        ELSE
                           PFT(LL,NI,:) = FT(LL,:)
                           PARTICLE_SLIDE = .FALSE.
                        END IF                     
                     END IF
                  END IF        !wall contact 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
               elseif(c_near_w(tempijk,7) .eq. 2)then ! the cell contact with wall through node
                  wallcontact=0
                  call wallnodecontact(tempijk,LL,wallcontact)
                  
                  IF(WALLCONTACT.EQ.1) THEN
                     WALLCHECK = 1
                     
                     NI = 2
                     CO = 0 
                     I = PARTICLES + nwalls + 2
                     
                     IF(PN(LL,1).GT.0) THEN ! check to see whether already in contact
                        IF((CO.EQ.0).AND.(NI.LE.(PN(LL,1)+1))) THEN
 103                       CONTINUE
                           IF(I.EQ.PN(LL,NI)) THEN
                              CO = 1
                              PV(LL,NI) = 1
                              
                              DES_POS_NEW(I,:) = DES_WALL_POS(1,:)
                              DES_VEL_NEW(I,:) = DES_WALL_VEL(1,:)
                              OMEGA_NEW(I,:) = ZERO
                              DES_RADIUS(I) = DES_RADIUS(LL)
                              CALL CFNORMALWALL(LL, I, NORMAL)
                              CALL CFRELVEL(LL, I, V_REL_TRANS, V_REL_TRANS_NORM, V_REL_TRANS_TANG, TANGENT, NORMAL)
                              CALL CFINCREMENTALOVERLAPS(LL,I, V_REL_TRANS_NORM,  V_REL_TRANS_TANG, OVERLAP_N, OVERLAP_T, CHECK_CON)
                              
                              IF(.not.CHECK_CON) THEN 
                                 PV(LL,NI) = 0 !No More Contact 
                                 GOTO 202
                              ENDIF
                              
                              FNS1(:) = -KN_W*((OVERLAP_N))*NORMAL(:)
                              FNS2(:) = -ETA_N_W*V_REL_TRANS_NORM*NORMAL(:)
                              
                              FTS1(:) = -KT_W*((OVERLAP_T)) *TANGENT(:)
                              FTS2(:) = - ETA_T_W*V_REL_TRANS_TANG*TANGENT(:)                        
                              
                              FT(LL,:) = FTS1(:) + FTS2(:) 
                              FN(LL,:) = FNS1(:) + FNS2(:) 
                              
                              FN(LL, :) = FN(LL,:) + PFN(LL,NI,:)
                              TEMPFT(:) = FT(LL, :) + PFT(LL,NI,:)
                              
                              CALL CFSLIDEWALL(LL, TANGENT, TEMPFT)
                              
                              CALL CFFCTOWALL(LL, NORMAL)
                              
                              PFN(LL,NI,:) =  PFN(LL,NI,:) + FNS1(:)
                              
                              IF(.NOT.PARTICLE_SLIDE) THEN
                                 PFT(LL,NI,:) = PFT(LL,NI,:) + FTS1(:)
                              ELSE
                                 PFT(LL,NI,:) =  FT(LL,:)
                                 PARTICLE_SLIDE = .FALSE.
                              END IF                        
                           ELSE
                              NI = NI + 1
                           END IF
                           
                           IF((CO.EQ.0).AND.(NI.LE.(PN(LL,1)+1))) GO TO 103
                           
                        END IF  !(CO.EQ.0).AND.(NI.LE.(PN(LL,1)+1))                      
                     END IF     !PN(LL,1).GT.0
                     IF(CO.EQ.0) THEN ! New contact
                        DES_POS_NEW(I,:) = DES_WALL_POS(1,:)
                        DES_VEL_NEW(I,:) = DES_WALL_VEL(1,:)
                        OMEGA_NEW(I,:) = ZERO
                        DES_RADIUS(I) = DES_RADIUS(LL)
                        CALL CFNORMALWALL(LL, I, NORMAL)
                        CALL CFRELVEL(LL, I, V_REL_TRANS, V_REL_TRANS_NORM, V_REL_TRANS_TANG, TANGENT, NORMAL)
                        CALL CFTOTALOVERLAPS(LL, I, II, V_REL_TRANS_NORM, V_REL_TRANS_TANG, OVERLAP_N, OVERLAP_T, CHECK_CON)
                        
                        IF(.not.CHECK_CON) THEN 
                           GOTO 202
                        ENDIF
                        
                        PN(LL,1) = PN(LL,1) + 1
                        NI = PN(LL,1) + 1
                        PN(LL,NI) = I
                        PV(LL,NI) = 1
                        
                        FNS1(:) = -KN_W*((OVERLAP_N))*NORMAL(:)
                        FNS2(:) = -ETA_N_W*V_REL_TRANS_NORM*NORMAL(:)
                        
                        FTS1(:) = -KT_W*((OVERLAP_T)) *TANGENT(:)
                        FTS2(:) = - ETA_T_W*V_REL_TRANS_TANG*TANGENT(:)                    
                        
                        FT(LL,:) = FTS1(:) + FTS2(:) 
                        FN(LL,:) = FNS1(:) + FNS2(:)                     
                        
                        TEMPFT(:) = FT(LL, :)
                        
                        CALL CFSLIDEWALL(LL, TANGENT, TEMPFT)
                        CALL CFFCTOWALL(LL, NORMAL)
                        
                        PFN(LL,NI,:) = FNS1(:)
                        
                        IF(.NOT.PARTICLE_SLIDE) THEN
                           PFT(LL,NI,:) =  FTS1(:)
                        ELSE
                           PFT(LL,NI,:) = FT(LL,:)
                           PARTICLE_SLIDE = .FALSE.
                        END IF                    
                     END IF     !CO.EQ.0
                  END IF        !wall contact 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
               elseif(c_near_w(tempijk,7) .eq. 3)then ! the cell contact with internal surface through face
                  do iw=1,nwalls
                     wallcontact=0
                     if(c_near_w(tempijk,iw) .eq. 1)then
                        call wallfacecontact(iw,LL,wallcontact)

                        IF(WALLCONTACT .EQ. 1) THEN
                           WALLCHECK = 1    
                           NI = 2
                           CO = 0 
                           I = PARTICLES + nwalls + 2 + IW

                           IF(PN(LL,1).GT.0) THEN ! check to see whether already in contact
                              IF((CO.EQ.0).AND.(NI.LE.(PN(LL,1)+1))) THEN
 104                             CONTINUE
                                 IF(I.EQ.PN(LL,NI)) THEN
                                    CO = 1
                                    PV(LL,NI) = 1
                                    
                                    DES_POS_NEW(I,:) = DES_WALL_POS(IW,:)
                                    DES_VEL_NEW(I,:) = DES_WALL_VEL(IW,:)
                                    OMEGA_NEW(I,:) = ZERO
                                    DES_RADIUS(I) = DES_RADIUS(LL)
                                    CALL CFNORMALWALL(LL, I, NORMAL)
                                    CALL CFRELVEL(LL, I, V_REL_TRANS, V_REL_TRANS_NORM, V_REL_TRANS_TANG, TANGENT, NORMAL)
                                    CALL CFINCREMENTALOVERLAPS(LL,I, V_REL_TRANS_NORM,  V_REL_TRANS_TANG, OVERLAP_N, OVERLAP_T, CHECK_CON)
                                    
                                    IF(.not.CHECK_CON) THEN 
                                       PV(LL,NI) = 0 !No More Contact 
                                       GOTO 203
                                    ENDIF
                                    
                                    FNS1(:) = -KN_W*((OVERLAP_N))*NORMAL(:)
                                    FNS2(:) = -ETA_N_W*V_REL_TRANS_NORM*NORMAL(:)
                                    
                                    FTS1(:) = -KT_W*((OVERLAP_T)) *TANGENT(:)
                                    FTS2(:) = - ETA_T_W*V_REL_TRANS_TANG*TANGENT(:)                            
                                    
                                    FT(LL,:) = FTS1(:) + FTS2(:) 
                                    FN(LL,:) = FNS1(:) + FNS2(:) 
                                    
                                    FN(LL, :) = FN(LL,:) + PFN(LL,NI,:)
                                    TEMPFT(:) = FT(LL, :) + PFT(LL,NI,:)
                                    
                                    CALL CFSLIDEWALL(LL, TANGENT, TEMPFT)
                                    
                                    CALL CFFCTOWALL(LL, NORMAL)
                                    
                                    PFN(LL,NI,:) =  PFN(LL,NI,:) + FNS1(:)
                                    
                                    IF(.NOT.PARTICLE_SLIDE) THEN
                                       PFT(LL,NI,:) = PFT(LL,NI,:) + FTS1(:)
                                    ELSE
                                       PFT(LL,NI,:) =  FT(LL,:)
                                       PARTICLE_SLIDE = .FALSE.
                                    END IF                      
                                 ELSE
                                    NI = NI + 1
                                 END IF
                                 
                                 IF((CO.EQ.0).AND.(NI.LE.(PN(LL,1)+1))) GO TO 104
                                 
                              END IF                        
                           END IF ! PN(LL,1).GT.0
                           IF(CO.EQ.0) THEN ! New contact
                              DES_POS_NEW(I,:) = DES_WALL_POS(IW,:)
                              DES_VEL_NEW(I,:) = DES_WALL_VEL(IW,:)
                              OMEGA_NEW(I,:) = ZERO
                              DES_RADIUS(I) = DES_RADIUS(LL)
                              CALL CFNORMALWALL(LL, I, NORMAL)
                              CALL CFRELVEL(LL, I, V_REL_TRANS, V_REL_TRANS_NORM, V_REL_TRANS_TANG, TANGENT, NORMAL)
                              CALL CFTOTALOVERLAPS(LL, I, II, V_REL_TRANS_NORM, V_REL_TRANS_TANG, OVERLAP_N, OVERLAP_T, CHECK_CON)
                              
                              IF(.not.CHECK_CON) THEN 
                                 GOTO 203
                              ENDIF
                              
                              PN(LL,1) = PN(LL,1) + 1
                              NI = PN(LL,1) + 1
                              PN(LL,NI) = I
                              PV(LL,NI) = 1
                              
                              FNS1(:) = -KN_W*((OVERLAP_N))*NORMAL(:)
                              FNS2(:) = -ETA_N_W*V_REL_TRANS_NORM*NORMAL(:)
                              
                              FTS1(:) = -KT_W*((OVERLAP_T)) *TANGENT(:)
                              FTS2(:) = - ETA_T_W*V_REL_TRANS_TANG*TANGENT(:)                        
                              
                              FT(LL,:) = FTS1(:) + FTS2(:) 
                              FN(LL,:) = FNS1(:) + FNS2(:)                         
                              
                              TEMPFT(:) = FT(LL, :)
                              
                              CALL CFSLIDEWALL(LL, TANGENT, TEMPFT)
                              CALL CFFCTOWALL(LL, NORMAL)
                              
                              PFN(LL,NI,:) = FNS1(:)
                              
                              IF(.NOT.PARTICLE_SLIDE) THEN
                                 PFT(LL,NI,:) =  FTS1(:)
                              ELSE
                                 PFT(LL,NI,:) = FT(LL,:)
                                 PARTICLE_SLIDE = .FALSE.
                              END IF                         
                           END IF
                        END IF  !wall contact                     
                     end if     !c_near_w(tmpijk,iw).eq.1
 203                 continue            		  
                  end do        ! iw=1,nwalls              
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
               elseif(c_near_w(tempijk,7) .eq. 4)then ! the cell contact with internal surface through edge
                  wallcontact=0
                  call walledgecontact(tempijk,LL,wallcontact)

                  IF(WALLCONTACT.EQ.1) THEN
                     WALLCHECK = 1
                     NI = 2
                     CO = 0 
                     I = PARTICLES + nwalls + 2 + nwalls + 1

                     IF(PN(LL,1).GT.0) THEN ! check to see whether already in contact
                        IF((CO.EQ.0).AND.(NI.LE.(PN(LL,1)+1))) THEN
 105                       CONTINUE
                           IF(I.EQ.PN(LL,NI)) THEN
                              CO = 1
                              PV(LL,NI) = 1
                              
                              DES_POS_NEW(I,:) = DES_WALL_POS(1,:)
                              DES_VEL_NEW(I,:) = DES_WALL_VEL(1,:)
                              OMEGA_NEW(I,:) = ZERO
                              DES_RADIUS(I) = DES_RADIUS(LL)
                              CALL CFNORMALWALL(LL, I, NORMAL)
                              CALL CFRELVEL(LL, I, V_REL_TRANS, V_REL_TRANS_NORM, V_REL_TRANS_TANG, TANGENT, NORMAL)
                              CALL CFINCREMENTALOVERLAPS(LL,I, V_REL_TRANS_NORM,  V_REL_TRANS_TANG, OVERLAP_N, OVERLAP_T, CHECK_CON)
                              
                              IF(.not.CHECK_CON) THEN 
                                 PV(LL,NI) = 0 !No More Contact 
                                 GOTO 202
                              ENDIF
                              
                              FNS1(:) = -KN_W*((OVERLAP_N))*NORMAL(:)
                              FNS2(:) = -ETA_N_W*V_REL_TRANS_NORM*NORMAL(:)

                              FTS1(:) = -KT_W*((OVERLAP_T)) *TANGENT(:)
                              FTS2(:) = - ETA_T_W*V_REL_TRANS_TANG*TANGENT(:)                         
                              
                              FT(LL,:) = FTS1(:) + FTS2(:) 
                              FN(LL,:) = FNS1(:) + FNS2(:) 
                              
                              FN(LL, :) = FN(LL,:) + PFN(LL,NI,:)
                              TEMPFT(:) = FT(LL, :) + PFT(LL,NI,:)
                              
                              CALL CFSLIDEWALL(LL, TANGENT, TEMPFT)
                              
                              CALL CFFCTOWALL(LL, NORMAL)

                              PFN(LL,NI,:) =  PFN(LL,NI,:) + FNS1(:)
                              
                              IF(.NOT.PARTICLE_SLIDE) THEN
                                 PFT(LL,NI,:) = PFT(LL,NI,:) + FTS1(:)
                              ELSE
                                 PFT(LL,NI,:) =  FT(LL,:)
                                 PARTICLE_SLIDE = .FALSE.
                              END IF
                           ELSE
                              NI = NI + 1
                           END IF
                           
                           IF((CO.EQ.0).AND.(NI.LE.(PN(LL,1)+1))) GO TO 105
                           
                        END IF  !(CO.EQ.0).AND.(NI.LE.(PN(LL,1)+1))                   
                     END IF     !PN(LL,1).GT.0
                     IF(CO.EQ.0) THEN ! New contact
                        DES_POS_NEW(I,:) = DES_WALL_POS(1,:)
                        DES_VEL_NEW(I,:) = DES_WALL_VEL(1,:)
                        OMEGA_NEW(I,:) = ZERO
                        DES_RADIUS(I) = DES_RADIUS(LL)
                        CALL CFNORMALWALL(LL, I, NORMAL)
                        CALL CFRELVEL(LL, I, V_REL_TRANS, V_REL_TRANS_NORM, V_REL_TRANS_TANG, TANGENT, NORMAL)
                        CALL CFTOTALOVERLAPS(LL, I, II, V_REL_TRANS_NORM, V_REL_TRANS_TANG, OVERLAP_N, OVERLAP_T, CHECK_CON)

                        IF(.not.CHECK_CON) THEN 
                           GOTO 202
                        ENDIF
                        
                        PN(LL,1) = PN(LL,1) + 1
                        NI = PN(LL,1) + 1
                        PN(LL,NI) = I
                        PV(LL,NI) = 1
                        
                        FNS1(:) = -KN_W*((OVERLAP_N))*NORMAL(:)
                        FNS2(:) = -ETA_N_W*V_REL_TRANS_NORM*NORMAL(:)
                        
                        FTS1(:) = -KT_W*((OVERLAP_T)) *TANGENT(:)
                        FTS2(:) = - ETA_T_W*V_REL_TRANS_TANG*TANGENT(:)                    
                        
                        FT(LL,:) = FTS1(:) + FTS2(:) 
                        FN(LL,:) = FNS1(:) + FNS2(:)                     

                        TEMPFT(:) = FT(LL, :)

                        CALL CFSLIDEWALL(LL, TANGENT, TEMPFT)
                        CALL CFFCTOWALL(LL, NORMAL)
                        
                        PFN(LL,NI,:) = FNS1(:)
                        
                        IF(.NOT.PARTICLE_SLIDE) THEN
                           PFT(LL,NI,:) =  FTS1(:)
                        ELSE
                           PFT(LL,NI,:) = FT(LL,:)
                           PARTICLE_SLIDE = .FALSE.
                        END IF                     
                     END IF
                  END IF        !wall contact              
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
               else
                  write(*,*) 'Error from calc_force_des'
                  write(*,*) 'c_near_w(ijk,7) is incorrect', c_near_w(tempijk,7)
                  write(*,*) 'It should be -2, -1, 0, 1, 2, 3 or 4'
                  stop
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
               end if
!-------------------------------
            else                ! NON_RECT_BC 
                                ! simple domain shape, rectangular boundary, should be used if the periodic boundary condition is enabled
!-------------------------------
               DO IW = 1, NWALLS
                  WALLCONTACT = 0
                  CALL CFWALLCONTACT(IW, LL, WALLCONTACT)
                  IF(WALLCONTACT.EQ.1) THEN
                     WALLCHECK = 1

                     NI = 2
                     CO = 0 
                     I = PARTICLES + IW

                     IF(PN(LL,1).GT.0) THEN ! check to see whether already in contact
                        IF((CO.EQ.0).AND.(NI.LE.(PN(LL,1)+1))) THEN
 100                       CONTINUE
                           IF(I.EQ.PN(LL,NI)) THEN
                              CO = 1
                              PV(LL,NI) = 1
                              
                              CALL CFWALLPOSVEL(LL, IW)
                              DES_POS_NEW(I,:) = DES_WALL_POS(IW,:)
                              DES_VEL_NEW(I,:) = DES_WALL_VEL(IW,:)
                              OMEGA_NEW(I,:) = ZERO
                              DES_RADIUS(I) = DES_RADIUS(LL)
                              CALL CFNORMALWALL(LL, I, NORMAL)
                              CALL CFRELVEL(LL, I, V_REL_TRANS, V_REL_TRANS_NORM, V_REL_TRANS_TANG, TANGENT, NORMAL)
                              CALL CFINCREMENTALOVERLAPS(LL,I, V_REL_TRANS_NORM,  V_REL_TRANS_TANG, OVERLAP_N, OVERLAP_T, CHECK_CON)
                              
                              IF(.not.CHECK_CON) THEN 
                                 PV(LL,NI) = 0 !No More Contact 
                                 GOTO 200
                              ENDIF
                              
                              FNS1(:) = -KN_W*((OVERLAP_N))*NORMAL(:)
                              FNS2(:) = -ETA_N_W*V_REL_TRANS_NORM*NORMAL(:)

                              FTS1(:) = -KT_W*((OVERLAP_T)) *TANGENT(:)
                              FTS2(:) = - ETA_T_W*V_REL_TRANS_TANG*TANGENT(:)
                              
                              
                              FT(LL,:) = FTS1(:) + FTS2(:) 
                              FN(LL,:) = FNS1(:) + FNS2(:) 
                              
                              FN(LL, :) = FN(LL,:) + PFN(LL,NI,:)
                              TEMPFT(:) = FT(LL, :) + PFT(LL,NI,:)
                              
                              CALL CFSLIDEWALL(LL, TANGENT, TEMPFT)
                              
                              CALL CFFCTOWALL(LL, NORMAL)

                              PFN(LL,NI,:) =  PFN(LL,NI,:) + FNS1(:)
                              
                              IF(.NOT.PARTICLE_SLIDE) THEN
                                 PFT(LL,NI,:) = PFT(LL,NI,:) + FTS1(:)
                              ELSE
                                 PFT(LL,NI,:) =  FT(LL,:)
                                 PARTICLE_SLIDE = .FALSE.
                              END IF

                                !--   DEBUGGING
                                !!impulse is effectively doubled for wall interactions
                              
                           ELSE
                              NI = NI + 1
                           END IF
                           IF((CO.EQ.0).AND.(NI.LE.(PN(LL,1)+1))) GO TO 100
                        END IF
                        
                     END IF
                     IF(CO.EQ.0) THEN ! New contact
                        CALL CFWALLPOSVEL(LL, IW)
                        DES_POS_NEW(I,:) = DES_WALL_POS(IW,:)
                        DES_VEL_NEW(I,:) = DES_WALL_VEL(IW,:)
                        OMEGA_NEW(I,:) = ZERO
                        DES_RADIUS(I) = DES_RADIUS(LL)
                        CALL CFNORMALWALL(LL, I, NORMAL)
                        CALL CFRELVEL(LL, I, V_REL_TRANS, V_REL_TRANS_NORM, V_REL_TRANS_TANG, TANGENT, NORMAL)
                        CALL CFTOTALOVERLAPS(LL, I, II, V_REL_TRANS_NORM, V_REL_TRANS_TANG, OVERLAP_N, OVERLAP_T, CHECK_CON)

                        IF(.not.CHECK_CON) THEN 
                           GOTO 200
                        ENDIF
                        
                        PN(LL,1) = PN(LL,1) + 1
                        NI = PN(LL,1) + 1
                        PN(LL,NI) = I
                        PV(LL,NI) = 1
                        
                        FNS1(:) = -KN_W*((OVERLAP_N))*NORMAL(:)
                        FNS2(:) = -ETA_N_W*V_REL_TRANS_NORM*NORMAL(:)
                        
                        FTS1(:) = -KT_W*((OVERLAP_T)) *TANGENT(:)
                        FTS2(:) = - ETA_T_W*V_REL_TRANS_TANG*TANGENT(:)
                        
                        
                        FT(LL,:) = FTS1(:) + FTS2(:) 
                        FN(LL,:) = FNS1(:) + FNS2(:) 
                        

                        TEMPFT(:) = FT(LL, :)

                        CALL CFSLIDEWALL(LL, TANGENT, TEMPFT)
                        CALL CFFCTOWALL(LL, NORMAL)
                        
                        PFN(LL,NI,:) = FNS1(:)
                        
                        IF(.NOT.PARTICLE_SLIDE) THEN
                           PFT(LL,NI,:) =  FTS1(:)
                        ELSE
                           PFT(LL,NI,:) = FT(LL,:)
                           PARTICLE_SLIDE = .FALSE.
                        END IF
                        
                     end IF
                  end IF        !wall contact 
 200              continue
               end DO
            end if              !NON_RECT_BC
                                !domain shape, rectangular boundary
!-------------------------------
         end IF                 !if(walldtsplit
 202     continue      
         
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
                           CALL CFRELVEL(LL, I, V_REL_TRANS, V_REL_TRANS_NORM, V_REL_TRANS_TANG, TANGENT, NORMAL)

                           CALL CFINCREMENTALOVERLAPS(LL,I, V_REL_TRANS_NORM,  V_REL_TRANS_TANG, OVERLAP_N, OVERLAP_T, CHECK_CON)
                           
                           IF(.not.CHECK_CON) THEN 
                              PV(LL,NI) = 0 !No More Contact 
                              GOTO 300 
                           ENDIF
                           IF(LL.EQ.FOCUS_PART2) Print*, 'CHECK CON =', CHECK_CON, I
                           FNS1(:) = -KN*((OVERLAP_N))*NORMAL(:)
                           FNS2(:) = -ETA_DES_N*V_REL_TRANS_NORM*NORMAL(:)
                           
                           FTS1(:) = -KT*((OVERLAP_T)) *TANGENT(:)
                           FTS2(:) = - ETA_DES_T*V_REL_TRANS_TANG*TANGENT(:)
                           
                                                      
                           FT(LL,:) = FTS1(:) + FTS2(:) 
                           FN(LL,:) = FNS1(:) + FNS2(:) 
                           
                           FN(LL, :) = FN(LL,:) + PFN(LL,NI,:)
                           TEMPFT(:) = FT(LL, :) + PFT(LL,NI,:)
                           
                           CALL CFSLIDE(LL, TANGENT, TEMPFT)
                           
                           PFN(LL,NI,:) = PFN(LL,NI,:) +  FNS1(:)

                           CALL CFFCTOW(LL, I, NORMAL)

                           IF(.NOT.PARTICLE_SLIDE) THEN
                              PFT(LL,NI,:) = PFT(LL,NI,:) + FTS1(:)
                           ELSE
                              PFT(LL,NI,:) =  FT(LL,:)
                              PARTICLE_SLIDE = .FALSE.
                           END IF

                        ELSE
                           NI = NI + 1
                        END IF
                        IF((CO.EQ.0).AND.(NI.LE.(PN(LL,1)+1))) GO TO 20
                     END IF
                     
                  END IF

                  IF(CO.EQ.0) THEN ! New contact
                     
                     CALL CFNORMAL(LL, I, II, NORMAL)
                     
                     CALL CFRELVEL(LL, I, V_REL_TRANS, V_REL_TRANS_NORM, V_REL_TRANS_TANG, TANGENT, NORMAL)
                     CALL CFTOTALOVERLAPS(LL, I, II, V_REL_TRANS_NORM, V_REL_TRANS_TANG, OVERLAP_N, OVERLAP_T, CHECK_CON)
                     
                     IF(.not.CHECK_CON) THEN 
			goto 300
                     end IF
                     PN(LL,1) = PN(LL,1) + 1
                     NI = PN(LL,1) + 1
                     PN(LL,NI) = I
                     PV(LL,NI) = 1

                     FNS1(:) = -KN*((OVERLAP_N))*NORMAL(:)
                     FNS2(:) = -ETA_DES_N*V_REL_TRANS_NORM*NORMAL(:)
                     
                     FTS1(:) = -KT*((OVERLAP_T)) *TANGENT(:)
                     FTS2(:) = - ETA_DES_T*V_REL_TRANS_TANG*TANGENT(:)
                                                                           
                     FT(LL,:) = FTS1(:) + FTS2(:) 
                     FN(LL,:) = FNS1(:) + FNS2(:) 

                     TEMPFT(:) = FT(LL, :)
                     CALL CFSLIDE(LL, TANGENT, TEMPFT)
                     
                     CALL CFFCTOW(LL, I, NORMAL)
                     PFN(LL,NI,:) = FNS1(:)
                     IF(.NOT.PARTICLE_SLIDE) THEN
                        PFT(LL,NI,:) =  FTS1(:)
                     ELSE
                        PFT(LL,NI,:) = FT(LL,:)
                        PARTICLE_SLIDE = .FALSE.
                     END IF
!400                  continue
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
300            CONTINUE
            END DO              !II = 2, NEIGHBOURS(LL,1)+I
         END IF
         
         IF((NEIGHBOURS(LL,1).EQ.0).AND.(WALLCHECK.EQ.0)) THEN
            CALL CFNOCONTACT(LL)
         END IF
         
      END DO
      calc_fc = .true.
      IF(DES_CONTINUUM_COUPLED) THEN
        CALLFROMDES = .TRUE.
        CALL DRAG_FGS(calc_fc,CALLFROMDES)
        CALLFROMDES = .FALSE.
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

