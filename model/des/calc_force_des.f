!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_FORCE_DES                                         C
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

!-----------------------------------------------
! Local variables
!-----------------------------------------------      
      INTEGER I, J, LL, II, IW
      INTEGER NI, NLIM, N_NOCON, NEIGH_L
      INTEGER OVERLAP_MAXP
      INTEGER WALLCONTACT, WALLCHECK

      DOUBLE PRECISION OVERLAP_N, OVERLAP_T
      DOUBLE PRECISION TEMPX, TEMPY, TEMPZ, TEMPD
      DOUBLE PRECISION V_REL_TRANS_NORM, V_REL_TRANS_TANG 
      DOUBLE PRECISION FTMD, FNMD
      DOUBLE PRECISION TEMPFT(DIMN)
      DOUBLE PRECISION FNS1(DIMN), FNS2(DIMN)
      DOUBLE PRECISION FTS1(DIMN), FTS2(DIMN)
      DOUBLE PRECISION NORMAL(DIMN), TANGENT(DIMN)
      DOUBLE PRECISION DIST(DIMN), DISTMOD, R_LM

      DOUBLE PRECISION, SAVE :: FTHIST(2)
      INTEGER, SAVE ::  CONTACT_COUNT = 1

      LOGICAL ALREADY_EXISTS
      LOGICAL ALREADY_NEIGHBOURS, OVERLAP_MAX_WALL

! Logical for local debug warnings
      LOGICAL DES_LOC_DEBUG
! index to track accounted for particles      
      INTEGER PC
! local values used damping coefficients
      DOUBLE PRECISION ETA_DES_N, ETA_DES_NW, ETA_DES_T, ETA_DES_TW
!-----------------------------------------------      
! Functions
!-----------------------------------------------      
      DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT 

!-----------------------------------------------      



!---------------------------------------------------------------------
! Calculate new values
!---------------------------------------------------------------------
      OVERLAP_MAXP = UNDEFINED_I
      OVERLAP_MAX_WALL = .FALSE.
      DES_LOC_DEBUG = .FALSE.

      OVERLAP_MAX = ZERO
      NEIGH_MAX = -1
      FOCUS_PARTICLE = 0

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
     
      PC = 1      
      DO LL = 1, MAX_PIS

         IF(PC .GT. PIS) EXIT
         IF(.NOT.PEA(LL,1)) CYCLE
         
         IF(DEBUG_DES.AND.LL.EQ.FOCUS_PARTICLE) THEN 
            IF (.NOT.DES_LOC_DEBUG) THEN
               DES_LOC_DEBUG = .TRUE.
               WRITE(*,1000)
            ENDIF
            WRITE(*,*) '     X,Y POS: ', DES_POS_NEW(LL,1), DES_POS_NEW(LL,2)
            WRITE(*,*) '     X,Y VEL: ', DES_VEL_NEW(LL,1), DES_VEL_NEW(LL,2)
         ENDIF

         FTS1(:) = ZERO
         FTS2(:) = ZERO
         FNS1(:) = ZERO
         FNS2(:) = ZERO

         TEMPFT(:) = ZERO

         IF(PN(LL,1).GE.1) THEN
            NLIM = PN(LL,1)+1
            N_NOCON = 0
! For each particle listed as in contact with particle LL in array PN,
! check the flag array PV to determine if particles remained in contact 
! after the previous call of CALC_FORCE_DES. 
            DO NI = 2, NLIM
               IF(PV(LL,NI-N_NOCON).EQ.0) THEN
! For each particle in PN(LL,2:MAXNEIGHBORS) that is no longer in 
! contact, shift the remaining particle contact information PN, PN,
! PFT left by one and reduce PN(LL,1) by one.
                  PN(LL,(NI-N_NOCON):MAXNEIGHBORS-1) = &
                     PN(LL,(NI-N_NOCON+1):MAXNEIGHBORS) 
                  PV(LL,(NI-N_NOCON):(MAXNEIGHBORS-1)) = &
                     PV(LL,(NI-N_NOCON+1):MAXNEIGHBORS)
                  PFT(LL,(NI-N_NOCON):(MAXNEIGHBORS-1),:) = &
                     PFT(LL,(NI-N_NOCON+1):MAXNEIGHBORS,:)
                  N_NOCON = N_NOCON + 1
                  PN(LL,1) = PN(LL,1) - 1
               ENDIF
            ENDDO
         ENDIF

! Initializing rest of the neighbor list which is not in contact and
! clean up after the above array left shifts
         NLIM = MAX(2,PN(LL,1) + 2) 
         PN(LL,NLIM:MAXNEIGHBORS) = -1
         PFT(LL,NLIM:MAXNEIGHBORS,:) = ZERO
         IF (PN(LL,1) .GT. NEIGH_MAX) NEIGH_MAX = PN(LL,1)


! Initializing the neighbor list contact information when particles are
! not in contact; i.e. when particle LL has no neighbors         
         IF (PN(LL,1).EQ.0) THEN
            PFT(LL,:,:) = ZERO
         ENDIF
         
! Reset the flag array PV; during each call to calc_force_des this
! variable tracks whether particle LL has any current neighbors
! the array is used in the next call to calc_force_des to update
! particle LL neighbor history above
         PV(LL,2:MAXNEIGHBORS) = 0


! Check particle LL for wall contacts
!---------------------------------------------------------------------
! Treats wall interaction also as a two-particle interaction but accounting
! for the wall properties
         IF(WALLDTSPLIT .AND. .NOT.PEA(LL,2) .AND. .NOT.PEA(LL,3)) THEN
            WALLCHECK = 0
            DO IW = 1, NWALLS
               WALLCONTACT = 0
! Check to see if a particle is in contact with any of the walls
               CALL CFWALLCONTACT(IW, LL, WALLCONTACT)
               IF(WALLCONTACT.EQ.1) THEN
                  WALLCHECK = 1

                  IF(DEBUG_DES) THEN
                     IF (.NOT.DES_LOC_DEBUG) THEN
                        DES_LOC_DEBUG = .TRUE.
                        WRITE(*,1000)
                     ENDIF
                     WRITE(*,*) '     Wall No. & Particle No.: ', IW, LL
                  ENDIF

! J.Musser : changed particles to max_pis                  
                  I = MAX_PIS + IW
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
                  
! Assign the wall particle a position and velocity
                  CALL CFWALLPOSVEL(LL, IW)

                  R_LM = DES_RADIUS(LL) + DES_RADIUS(LL)
                  DIST(:) = DES_WALL_POS(IW,:) - DES_POS_NEW(LL,:)
                  DISTMOD = SQRT(DES_DOTPRDCT(DIST,DIST))

                  IF(R_LM - DISTMOD.GT.SMALL_NUMBER) THEN 
                     
                     IF((((R_LM-DISTMOD)/R_LM)*100.d0).GT.OVERLAP_MAX) THEN
                        OVERLAP_MAX = (((R_LM-DISTMOD)/R_LM)*100.d0)
                        OVERLAP_MAXP = LL
                        OVERLAP_MAX_WALL = .TRUE.
                     ENDIF

                     IF(DISTMOD.NE.ZERO) THEN
                        NORMAL(:)= DIST(:)/DISTMOD
                     ELSE 
                        IF (.NOT.DES_LOC_DEBUG) THEN
                           DES_LOC_DEBUG = .TRUE.
                           WRITE(*,1000)
                        ENDIF                             
                        WRITE(*,*) '     Particle ', LL, ' completely',&
                           ' crossed through wall ', I, &
                           ' so that DISTMOD is zero'
                        STOP
                     ENDIF

! Calculate the translational relative velocity for a contacting particle pair
                     CALL CFRELVEL(LL, IW, V_REL_TRANS_NORM, &
                        V_REL_TRANS_TANG, TANGENT, NORMAL, 1)

! The normal overlap calculation was changed so that it no longer
! depends on the contact history (i.e., integration of incremental
! overlap found by velocity*dtsolid).  Now overlap is based purely on
! position of neighbors.  The change was made because the former
! method was found not to conserve energy 
                     OVERLAP_N =  R_LM-DISTMOD 

                     IF(ALREADY_NEIGHBOURS) THEN 
                        PV(LL,NI) = 1
                        contact_count = contact_count + 1
                        OVERLAP_T = V_REL_TRANS_TANG*DTSOLID
                     ELSE 
                        PN(LL,1) = PN(LL,1) + 1
                        NI = PN(LL,1) + 1
                        PN(LL,NI) = I
                        PV(LL,NI) = 1
                        contact_count = 1
                        FTHIST = ZERO
                        OVERLAP_T = V_REL_TRANS_TANG*DTSOLID
                     ENDIF

                  ELSE
                     GOTO 200
                  ENDIF
                  
                  ETA_DES_NW = DES_ETAN_WALL(PIJK(LL,5))
                  ETA_DES_TW = DES_ETAT_WALL(PIJK(LL,5))
                  
                  FNS1(:) = -KN_W*((OVERLAP_N))*NORMAL(:)
                  FNS2(:) = -ETA_DES_NW*V_REL_TRANS_NORM*NORMAL(:)
                  
                  FTS1(:) = -KT_W*((OVERLAP_T)) *TANGENT(:)
                  FTS2(:) = -ETA_DES_TW*V_REL_TRANS_TANG*TANGENT(:)
                  
                  FT(LL,:) = FTS1(:) + FTS2(:) 
                  FN(LL,:) = FNS1(:) + FNS2(:) 
                  
                  TEMPFT(:) = FT(LL, :) + PFT(LL,NI,:)
                  
! Check for Coulombs friction law and limit the maximum value of the
! tangential force on a particle in contact with a wall
                  CALL CFSLIDEWALL(LL, TANGENT, TEMPFT)
                  
! Calculate the total force Fc and Tow on a particle in a particle-wall
! collision
                  CALL CFFCTOWALL(LL, NORMAL)
                  
                  PFT(LL,NI,:) = PFT(LL,NI,:) + FTS1(:)

                  IF(DEBUG_DES.AND.LL.EQ.FOCUS_PARTICLE) THEN
                     IF (.NOT.DES_LOC_DEBUG) THEN
                        DES_LOC_DEBUG = .TRUE.
                        WRITE(*,1000)
                     ENDIF                          
                     WRITE(*,*) '     WALL CONTACT ON ', NI
                     WRITE(*,*) '     ALREADY_NEIGHBOURS? = ',&
                        ALREADY_NEIGHBOURS
                     WRITE(*,*) '     STIME, DTSOLID = ', S_TIME, DTSOLID
                     WRITE(*,*) '     DES_VEL = ', DES_VEL_NEW(LL,1:DIMN),&
                        des_radius(LL)*OMEGA_NEW(LL,1)
                     WRITE(*,*) '     MAA'
                     WRITE(*,*) '     V-OMEGA R = ', &
                        DES_VEL_NEW(LL,1)+des_radius(LL)* OMEGA_NEW(LL,1),&
                        (DES_VEL_NEW(LL,1)+des_radius(LL)*OMEGA_NEW(LL,1))*DTSOLID
                     WRITE(*,*) '     Mg = ', PMASS(LL)*gravity
                     WRITE(*,*) '     KN_W, ETA_DES_NW, KT_W, ETA_DES_TW = ',&
                        KN_W, ETA_DES_NW, KT_W, ETA_DES_TW
                     WRITE(*,*) '     TANGENT= ', TANGENT
                     WRITE(*,*) '     HIST = ', PFT(LL,NI,1:2)
                     WRITE(*,*) '     PARTICLE_SLIDE ? ', PARTICLE_SLIDE
                     WRITE(*,*) '     FT and FN= ', FT( LL,:), FN(LL,:)
                     WRITE(*,*) '     KW*OT*TAN = ', &
                        KT_W*((OVERLAP_T)) *TANGENT(:)
                     WRITE(*,*) '     OVERLAP_T = ', OVERLAP_T, TANGENT
                     FTMD = SQRT(DES_DOTPRDCT(TEMPFT,TEMPFT))
                     FNMD = SQRT(DES_DOTPRDCT(FN(LL,1:DIMN),FN(LL,1:DIMN)))
                     WRITE(*,*) '     FTMD, mu FNMD = ', FTMD, MEW_W*FNMD
                     FTHIST(:) = FTHIST(:) + FT(LL,:)
                     PRINT*, '     FT AVG = ', FTHIST/contact_count, mew_w*PMASS(LL)*gravity, contact_count
                     READ(*,*)
                  ENDIF

                  PARTICLE_SLIDE = .FALSE.

               ENDIF   !wall contact
 200           CONTINUE
            ENDDO
         ENDIF   !if(walldtsplit .and. .not.pea(LL,2))
!---------------------------------------------------------------------
! End check particle LL for wall contacts         


! Check particle LL neighbour contacts         
!---------------------------------------------------------------------
         IF (NEIGHBOURS(LL,1).GT.0) THEN
            DO II = 2, NEIGHBOURS(LL,1)+1
               I = NEIGHBOURS(LL,II)

               IF(I.GT.LL .AND. PEA(I,1)) THEN
                  
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
                           IF (.NOT.DES_LOC_DEBUG) THEN
                              DES_LOC_DEBUG = .TRUE.
                              WRITE(*,1000)
                           ENDIF                                
                           WRITE(*,*) '     PARTICLE & NEIGHBOR ',&
                              'LL, I: ', LL, I
                           WRITE(*,*) '     OLD POS, I', DES_POS_NEW(I,:)
                           WRITE(*,*),'     OLD POS, LL = ', DES_POS_NEW(LL,:)
                        ENDIF

                        IF(TEMPX.GT.DES_POS_NEW(LL,1)) THEN 
                           DES_POS_NEW(I,1) = DES_POS_NEW(I,1) - (EX2-WX1)

                           IF(DEBUG_DES) THEN
                              IF (.NOT.DES_LOC_DEBUG) THEN
                                 DES_LOC_DEBUG = .TRUE.
                                 WRITE(*,1000) 
                              ENDIF                                
                              WRITE(*,*) '     NEW POS WEST= ', DES_POS_NEW(I,1)
                           ENDIF
                        ELSE
                           DES_POS_NEW(I,1) = DES_POS_NEW(I,1) + EX2 - WX1

                           IF(DEBUG_DES) THEN
                              IF (.NOT.DES_LOC_DEBUG) THEN
                                 DES_LOC_DEBUG = .TRUE.
                                 WRITE(*,1000)
                              ENDIF
                              WRITE(*,*) '     NEW POS EAST = ', DES_POS_NEW(I,1)
                           ENDIF
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

                  IF(R_LM - DISTMOD.GT.SMALL_NUMBER) THEN

                     IF(DEBUG_DES .AND. LL.EQ.FOCUS_PARTICLE) THEN
                        IF (.NOT.DES_LOC_DEBUG) THEN
                           DES_LOC_DEBUG = .TRUE.
                           WRITE(*,1000) 
                        ENDIF
                        WRITE(*,*) '     NEIGHBORS: ', NEIGHBOURS(LL,:)
                     ENDIF

                     IF((((R_LM-DISTMOD)/R_LM)*100.d0).GT.OVERLAP_MAX) THEN
                        OVERLAP_MAX = (((R_LM-DISTMOD)/R_LM)*100.d0)
                        OVERLAP_MAXP = LL
                        OVERLAP_MAX_WALL = .FALSE.
                     ENDIF

                     IF(DISTMOD.NE.ZERO) THEN
                        NORMAL(:)= DIST(:)/DISTMOD
                     ELSE 
                        IF (.NOT.DES_LOC_DEBUG) THEN
                           DES_LOC_DEBUG = .TRUE.
                           WRITE(*,1000)
                        ENDIF
                        WRITE(*,*) '     Particle pair ', LL, I, &
                           'crossed through one another ', &
                           'so that DISTMOD is zero'
                        STOP
                     ENDIF

! Calculate the translational relative velocity for a contacting particle pair.
                     CALL CFRELVEL(LL, I, V_REL_TRANS_NORM, &
                        V_REL_TRANS_TANG, TANGENT, NORMAL, 0)

! Overlap calculation changed from history based to current position 
                     OVERLAP_N = R_LM-DISTMOD

                     IF(ALREADY_NEIGHBOURS) THEN 
                        PV(LL,NI) = 1
                        OVERLAP_T = V_REL_TRANS_TANG*DTSOLID
                     ELSE 
                        IF(DEBUG_DES) THEN
                           IF (.NOT.DES_LOC_DEBUG) THEN
                              DES_LOC_DEBUG = .TRUE.
                              WRITE(*,1000)
                           ENDIF
                           WRITE(*,*) '     Overlap for ',&
                              ' particle pair: ', LL, I, &
                              ' is: ', ((R_LM-DISTMOD)/R_LM)*100.d0
                        ENDIF
                        PN(LL,1) = PN(LL,1) + 1
                        NI = PN(LL,1) + 1
                        PN(LL,NI) = I
                        PV(LL,NI) = 1
                        OVERLAP_T = V_REL_TRANS_TANG*DTSOLID
                     ENDIF
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
                  
                  TEMPFT(:) = FT(LL, :) + PFT(LL,NI,:)
                  
                  IF(DEBUG_DES.AND.LL.EQ.FOCUS_PARTICLE) THEN 
                     IF (.NOT.DES_LOC_DEBUG) THEN
                        DES_LOC_DEBUG = .TRUE.
                        WRITE(*,1000) 
                     ENDIF

                     PRINT*, '     I = ', I
                     PRINT*, '     EtaN, EtaT =  ', ETA_DES_N, ETA_DES_T
                     PRINT*, '     Overlap = ', overlap_n, (R_LM - DISTMOD)*100.d0/R_LM
                     PRINT*, '     rad ratio = ', DES_RADIUS(LL)/DES_RADIUS(I)
                     PRINT*, '     FNS1 and FNS2 = ', FNS1(:), FNS2(:)
                     PRINT*, '     PFT = ', PFT(LL,NI,:)
                     PRINT*, '     FORCEST = ', FT(LL,:)
                     PRINT*, '     FORCESN = ', FN(LL,:)
                     PRINT*, '     FORCEST = ', FT(LL,:)
                  ENDIF

! Check for Coulombs friction law and limit the maximum value of the tangential
! force on a particle in contact with another particle
                  CALL CFSLIDE(LL, TANGENT, TEMPFT)
                  
! Calculate the total force Fc and Tow on a particle in a particle-particle collision
                  CALL CFFCTOW(LL, I, NORMAL)

                  PFT(LL,NI,:) = PFT(LL,NI,:) + FTS1(:)
                  PARTICLE_SLIDE = .FALSE.
                  
                  IF(DEBUG_DES.AND.LL.eq.FOCUS_PARTICLE)THEN
                     INQUIRE(FILE='debug_file',EXIST=ALREADY_EXISTS)
                     IF(ALREADY_EXISTS)THEN
                        OPEN(UNIT=1,FILE='debug_file',STATUS='OLD',POSITION='APPEND')
                        WRITE(1,'(A,I5)')'CALC FORCE -- NEIGHBOR',II
                        WRITE(1,'(2(1x,A,E12.5))')&
                        'FNx=',FN(LL,1), 'FNy=',FN(LL,2)
                     ELSE
                        OPEN(UNIT=1,FILE='debug_file',STATUS='NEW')
                        WRITE(1,'(A,I5)')'CALC FORCE -- NEIGHBOR',II
                        WRITE(1,'(2(1x,A,E12.5))')&
                        'FNx=',FN(LL,1),'FNy=',FN(LL,2)
                     ENDIF
                     CLOSE (1)
                     Print*, 'PN', PN(LL,:)
                  ENDIF
                  
               ENDIF         ! IF (I>LL .AND. PEA(I,1))
 300           CONTINUE
            ENDDO            ! DO II = 2, NEIGHBOURS(LL,1)+I
         ENDIF               ! IF(NEIGHBOURS(LL,1).GT.0)

!---------------------------------------------------------------------
! End check particle LL neighbour contacts         


         IF((NEIGHBOURS(LL,1).EQ.0).AND.(WALLCHECK.EQ.0)) THEN
! The subroutine sets all forces on a particle to zero is the particle is 
! found to have no neighbors (either particles or walls)
            CALL CFNOCONTACT(LL)
         ENDIF
      
      PC = PC + 1
      ENDDO   ! end loop over paticles LL


! Treats wall interaction also as a two-particle interaction but accounting
! for the wall properties
      IF(DES_CONTINUUM_COUPLED) THEN
         CALL DRAG_FGS
      ENDIF

! COHESION
      IF(USE_COHESION)THEN
         CALL CALC_COHESIVE_FORCES
      ENDIF
      
! Update the old values of particle position and velocity with the new values computed
      CALL CFUPDATEOLD

      IF (DES_LOC_DEBUG) WRITE(*,1001)

 1000 FORMAT('---------- START CALC_FORCE_DES ---------->')
 1001 FORMAT('<---------- END CALC_FORCE_DES ----------') 

      RETURN
      END SUBROUTINE CALC_FORCE_DES

