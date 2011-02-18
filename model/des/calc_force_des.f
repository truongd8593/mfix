!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_FORCE_DES                                         C
!
!  Purpose: DES calculations of force acting on a particle, 
!           its velocity and its position                  
!
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

      DOUBLE PRECISION OVERLAP_N, OVERLAP_T, SQRT_OVERLAP
      DOUBLE PRECISION FRAC_OVERLAP1, FRAC_OVERLAP2
! percent of particle radius when excess overlap will be flagged
      DOUBLE PRECISION, PARAMETER :: flag_overlap = 0.20d0

      DOUBLE PRECISION V_REL_TRANS_NORM, V_REL_TRANS_TANG, &
                       V_REL_NORM_OLD, VRN_OLD(DIMN)
! temporary storage of tangential force
      DOUBLE PRECISION FT_TMP(DIMN)
! temporary storage of tangential displacement
      DOUBLE PRECISION PFT_TMP(DIMN)
! local normal and tangential forces      
      DOUBLE PRECISION FNS1(DIMN), FNS2(DIMN)
      DOUBLE PRECISION FTS1(DIMN), FTS2(DIMN)
! tmp variables to calculate magnitude of forces
      DOUBLE PRECISION FTMD, FNMD
      DOUBLE PRECISION NORMAL(DIMN), NORMAL_OLD(DIMN), TANGENT(DIMN)
      DOUBLE PRECISION DIST(DIMN), DIST_OLD(DIMN), DISTMOD, &
                       DISTMOD_OLD, R_LM
      DOUBLE PRECISION DTSOLID_TMP

! temporary variables for periodic/LE boundaries      
      DOUBLE PRECISION TMP_PART_POS(DIMN), TMP_PART_VEL(DIMN)

! logic flag telling whether contact pair is old      
      LOGICAL ALREADY_NEIGHBOURS
! logic flag for local debug warnings
      LOGICAL DES_LOC_DEBUG
      LOGICAL ALREADY_EXISTS
! index to track accounted for particles      
      INTEGER PC
! store solids phase index of particle (i.e. pijk(np,5))
      INTEGER PHASEI, PHASELL
! local values used spring constants and damping coefficients
      DOUBLE PRECISION ETAN_DES, ETAN_DES_W, ETAT_DES, ETAT_DES_W,&
                       KN_DES, KN_DES_W, KT_DES, KT_DES_W

!-----------------------------------------------      
! Functions
!-----------------------------------------------      
      DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT 

!-----------------------------------------------      


! Initialization
      OVERLAP_MAXP = UNDEFINED_I
      DES_LOC_DEBUG = .FALSE.
      OVERLAP_MAX = ZERO
      NEIGH_MAX = -1



! Calculate contact force and torque
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
            WRITE(*,'(5X,A,I10)') 'On Particle ', LL
            WRITE(*,'(5X,A,2(ES15.7))') &
               'X,Y POS: ', DES_POS_NEW(LL,1), DES_POS_NEW(LL,2)
            WRITE(*,'(5X,A,2(ES15.7))') &
               'X,Y VEL: ', DES_VEL_NEW(LL,1), DES_VEL_NEW(LL,2)
         ENDIF

         TANGENT(:) = ZERO
         NORMAL(:) = ZERO
         FTS1(:) = ZERO
         FTS2(:) = ZERO
         FNS1(:) = ZERO
         FNS2(:) = ZERO
         PFT_TMP(:) = ZERO

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
! for the wall properties; make sure the particle is not classified as
! a new 'entering' particle or is already marked as a potential exiting particle
         IF(WALLDTSPLIT .AND. .NOT.PEA(LL,2) .AND. .NOT.PEA(LL,3)) THEN
            WALLCHECK = 0

            DO IW = 1, NWALLS
               WALLCONTACT = 0

! Check to see if a particle is in contact with any of the walls
               CALL CFWALLCONTACT(IW, LL, WALLCONTACT)

               IF(WALLCONTACT.EQ.1) THEN
                  WALLCHECK = 1

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
! Particle overlap detected (i.e. resolve collision)
!---------------------------------- 

! Error reporting                     
                     FRAC_OVERLAP1 = (R_LM-DISTMOD)/DES_RADIUS(LL)
                     IF (FRAC_OVERLAP1 > flag_overlap) THEN
                        WRITE(*,'(5X,A,A,ES15.7)') &
                           'WARNING: excessive overlap detected ', &
                           'at time ', S_TIME
                        WRITE(*,'(7X,A,I10,2X,A,I5,2X,A)') &
                           'between particle ', LL, 'and wall ',&
                           IW, 'with'
                        WRITE(*,'(7X,A,ES15.7,2X,A,ES15.7)') &
                          'overlap = ', (R_LM-DISTMOD), &
                           ' radius = ', DES_RADIUS(LL)
                     ENDIF

! Des_time_march periodically reports max overlap info in the screen log                     
                     IF((((R_LM-DISTMOD)/R_LM)*100.d0).GT.OVERLAP_MAX) THEN
                        OVERLAP_MAX = (((R_LM-DISTMOD)/R_LM)*100.d0)
                        OVERLAP_MAXP = LL
                     ENDIF

                     IF(DISTMOD.NE.ZERO) THEN
                        NORMAL(:)= DIST(:)/DISTMOD
                     ELSE 
                        IF (.NOT.DES_LOC_DEBUG) THEN
                           DES_LOC_DEBUG = .TRUE.
                           WRITE(*,1000)
                        ENDIF                             
                        WRITE(*,'(5X,A,I10,A,I5)') &
                           'DISTMOD is zero between particle ', LL, &
                           ' and wall ', IW
                        STOP
                     ENDIF

! Calculate the translational relative velocity for a contacting particle pair
                     CALL CFRELVEL(LL, IW, V_REL_TRANS_NORM, &
                        V_REL_TRANS_TANG, TANGENT, NORMAL, DISTMOD, 1)

! The normal overlap calculation was changed so that it no longer
! depends on the contact history (i.e., integration of incremental
! overlap found by velocity*dtsolid).  Now overlap is based purely on
! position of neighbors.  The change was made because the former
! method was found not to conserve energy 
                     OVERLAP_N =  R_LM-DISTMOD 

                     IF(ALREADY_NEIGHBOURS) THEN 
                        PV(LL,NI) = 1
                        OVERLAP_T = V_REL_TRANS_TANG*DTSOLID
                     ELSE 
                        PN(LL,1) = PN(LL,1) + 1
                        NI = PN(LL,1) + 1
                        PN(LL,NI) = I
                        PV(LL,NI) = 1
                        IF (V_REL_TRANS_NORM .GT. ZERO) THEN
                           DTSOLID_TMP = OVERLAP_N/(V_REL_TRANS_NORM)
                        ELSEIF (V_REL_TRANS_NORM .LT. ZERO) THEN
                          DTSOLID_TMP = DTSOLID
                          WRITE(*,'(5X,A,A,ES15.7)') &
                             'WARNING: normal relative velocity less ',&
                             'than zero at time ', S_TIME
                          WRITE(*,'(7X,A,I10,2X,A,I5,2X,A)') &
                             'for first contact between particle', LL, &
                             'and wall ', IW, 'with'
                          WRITE(*,'(7X,A,ES15.7)') &
                             'V_REL_NORM = ', V_REL_TRANS_NORM
                        ELSE
                           DTSOLID_TMP = OVERLAP_N/&
                              (V_REL_TRANS_NORM+SMALL_NUMBER)
                        ENDIF
                        OVERLAP_T = V_REL_TRANS_TANG*&
                           MIN(DTSOLID,DTSOLID_TMP)
                     ENDIF

                     phaseLL = PIJK(LL,5) 

! T.Li : Hertz vs linear spring-dashpot contact model
                     IF (TRIM(DES_COLL_MODEL) .EQ. 'HERTZIAN') THEN
                        sqrt_overlap = SQRT(OVERLAP_N)
                        KN_DES_W = hert_kwn(phaseLL)*sqrt_overlap
                        KT_DES_W = hert_kwt(phaseLL)*sqrt_overlap
                        sqrt_overlap = SQRT(sqrt_overlap)
                        ETAN_DES_W = DES_ETAN_WALL(phaseLL)*sqrt_overlap
                        ETAT_DES_W = DES_ETAT_WALL(phaseLL)*sqrt_overlap
                     ELSE
                        KN_DES_W = KN_W
                        KT_DES_W = KT_W
                        ETAN_DES_W = DES_ETAN_WALL(phaseLL)
                        ETAT_DES_W = DES_ETAT_WALL(phaseLL)
                     ENDIF

! Calculate the normal contact force
                     FNS1(:) = -KN_DES_W * OVERLAP_N * NORMAL(:)
                     FNS2(:) = -ETAN_DES_W * V_REL_TRANS_NORM*NORMAL(:)
                     FN(LL,:) = FNS1(:) + FNS2(:) 

! Calculate the tangential displacement which is integration of tangential 
! relative velocity with respect to contact time. Correction in the tangential 
! direction is imposed                   
                     PFT(LL,NI,:) = PFT(LL,NI,:)+OVERLAP_T*TANGENT(:)
                     PFT_TMP(:) = PFT(LL,NI,:)   ! update pft_tmp before it is used
                     PFT_TMP(:) = PFT(LL,NI,:) - &
                        DES_DOTPRDCT(PFT_TMP,NORMAL)*NORMAL(:)
! Calculate the tangential contact force
                     FTS1(:) = -KT_DES_W * PFT_TMP(:)
                     FTS2(:) = -ETAT_DES_W * V_REL_TRANS_TANG*TANGENT(:)
                     FT(LL,:) = FTS1(:) + FTS2(:) 

! Temporary storage of tangential contact force for reporting
                     FT_TMP(:) = FT(LL,:)

! Check for Coulombs friction law and limit the maximum value of the
! tangential force on a particle in contact with a wall
                     CALL CFSLIDEWALL(LL, TANGENT)
                  
! Calculate the total force FC and torque TOW on a particle in a
! particle-wall collision
                     CALL CFFCTOWALL(LL, NORMAL, DISTMOD)
                  
! Save the tangential displacement history with the correction of Coulomb's law
                     IF (PARTICLE_SLIDE) THEN
! Since FT might be corrected during the call to cfslide, the tangental
! displacement history needs to be changed accordingly                          
                        PFT(LL,NI,:) = -( FT(LL,:) - FTS2(:) ) / KT_DES_W
                     ELSE
                        PFT(LL,NI,:) = PFT_TMP(:)
                     ENDIF 

! Reporting info                     
                     IF(DEBUG_DES.AND.LL.EQ.FOCUS_PARTICLE) THEN
                        IF (.NOT.DES_LOC_DEBUG) THEN
                           DES_LOC_DEBUG = .TRUE.
                           WRITE(*,1000)
                        ENDIF                          
                        WRITE(*,*) '     STIME, DTSOLID = ', S_TIME, DTSOLID
                        WRITE(*,*) '     WALL CONTACT ON WALL =', IW
                        WRITE(*,*) '     ALREADY_NEIGHBOURS = ',&
                           ALREADY_NEIGHBOURS
                        WRITE(*,*) '     DES_VEL = ', DES_VEL_NEW(LL,1:DIMN),&
                           des_radius(LL)*OMEGA_NEW(LL,1)
                        WRITE(*,*) '     V-OMEGA R = ', &
                           DES_VEL_NEW(LL,1)+des_radius(LL)* OMEGA_NEW(LL,1),&
                           (DES_VEL_NEW(LL,1)+des_radius(LL)*OMEGA_NEW(LL,1))*DTSOLID
                        WRITE(*,*) '     M*g = ', PMASS(LL)*gravity
                        WRITE(*,*) '     KN_W, ETAN_W, KT_W, ETAT_W = ',&
                           KN_DES_W, ETAN_DES_W, KT_DES_W, ETAT_DES_W
                        WRITE(*,*) '     TANGENT= ', TANGENT
                        WRITE(*,*) '     HIST = ', PFT(LL,NI,1:2)
                        WRITE(*,*) '     PARTICLE_SLIDE ? ', PARTICLE_SLIDE
                        WRITE(*,*) '     FT and FN= ', FT( LL,:), FN(LL,:)
                        WRITE(*,*) '     KT_W*OT*TAN = ', &
                           KT_DES_W*((OVERLAP_T)) *TANGENT(:)
                        WRITE(*,*) '     OVERLAP_T = ', OVERLAP_T, TANGENT
                        FTMD = SQRT(DES_DOTPRDCT(FT_TMP,FT_TMP))
                        FNMD = SQRT(DES_DOTPRDCT(FN(LL,1:DIMN),FN(LL,1:DIMN)))
                        WRITE(*,*) '     FTMD, mu FNMD = ', FTMD, MEW_W*FNMD
                     ENDIF

                     PARTICLE_SLIDE = .FALSE.

                  ENDIF   ! IF(R_LM - DISTMOD.GT.SMALL_NUMBER) -> resolve collision


              ENDIF   ! IF(WALLCONTACT.EQ.1) 
            ENDDO   ! DO IW = 1, NWALLS
         ENDIF   ! IF(WALLDTSPLIT .AND. .NOT.PEA(LL,2) .AND. .NOT.PEA(LL,3))
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
                  
! If necessary shift particle LL position/velocity accordingly to
! boundary conditions 
                  IF (DES_LE_BC) THEN
                     CALL DES_LEBC_NEIGHBOR_CHECK(LL,I,TMP_PART_POS,&
                        TMP_PART_VEL)
                  ELSEIF(DES_PERIODIC_WALLS) THEN
                     CALL DES_PERIODIC_NEIGHBOR_CHECK(LL,I,&
                        TMP_PART_POS)
                  ENDIF

                  R_LM = DES_RADIUS(LL) + DES_RADIUS(I)
                  DIST(:) = DES_POS_NEW(I,:) - DES_POS_NEW(LL,:)
                  DISTMOD = SQRT(DES_DOTPRDCT(DIST,DIST))
                  

                  IF(R_LM - DISTMOD.GT.SMALL_NUMBER) THEN
! Particle overlap detected (i.e. resolve collision)
!---------------------------------- 

! Error reporting
                     IF(DEBUG_DES .AND. LL.EQ.FOCUS_PARTICLE) THEN
                        IF (.NOT.DES_LOC_DEBUG) THEN
                           DES_LOC_DEBUG = .TRUE.
                           WRITE(*,1000) 
                        ENDIF
                        WRITE(*,'(5X,A,I10)') 'NEIGHBORS: ', NEIGHBOURS(LL,:)
                     ENDIF
                     
                     FRAC_OVERLAP1 = (R_LM-DISTMOD)/DES_RADIUS(LL)
                     FRAC_OVERLAP2 = (R_LM-DISTMOD)/DES_RADIUS(I)
                     IF (FRAC_OVERLAP1 > flag_overlap .OR. &
                         FRAC_OVERLAP2 > flag_overlap) THEN
                        WRITE(*,'(5X,A,A,ES15.7)') &
                           'WARNING: excessive overlap detected ', &
                           'at time ', S_TIME
                        WRITE(*,'(7X,A,I10,2X,A,I5,2X,A)') &
                           'between particles ', LL, 'and ',&
                           I, 'with'
                        WRITE(*,'(7X,A,ES15.7,2X,A,ES15.7,2X,ES15.7)') &
                           'overlap = ', (R_LM-DISTMOD), &
                           ' radii = ', DES_RADIUS(LL), DES_RADIUS(I)
                     ENDIF

! Des_time_march periodically reports max overlap info in the screen log
                     IF((((R_LM-DISTMOD)/R_LM)*100.d0).GT.OVERLAP_MAX) THEN
                        OVERLAP_MAX = (((R_LM-DISTMOD)/R_LM)*100.d0)
                        OVERLAP_MAXP = LL
                     ENDIF

                     IF(DISTMOD.NE.ZERO) THEN
                        NORMAL(:)= DIST(:)/DISTMOD
                     ELSE 
                        IF (.NOT.DES_LOC_DEBUG) THEN
                           DES_LOC_DEBUG = .TRUE.
                           WRITE(*,1000)
                        ENDIF
                        WRITE(*,'(5X,A,I10,I10)') &
                           'DISTMOD is zero between particle-pair ',&
                           LL, I
                        STOP
                     ENDIF

! Calculate the translational relative velocity for a contacting particle pair.
                     CALL CFRELVEL(LL, I, V_REL_TRANS_NORM, &
                        V_REL_TRANS_TANG, TANGENT, NORMAL, DISTMOD, 0)

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
                           WRITE(*,'(5X,A,2(I10,X),A,ES15.7)') &
                              'Normal overlap for particle pair ',&
                              LL, I, ' : ', OVERLAP_N 
                        ENDIF
                        PN(LL,1) = PN(LL,1) + 1
                        NI = PN(LL,1) + 1
                        PN(LL,NI) = I
                        PV(LL,NI) = 1
                        IF (V_REL_TRANS_NORM .GT. ZERO) THEN
                           DTSOLID_TMP = OVERLAP_N/(V_REL_TRANS_NORM)
                        ELSEIF (V_REL_TRANS_NORM .LT. ZERO) THEN
                          DTSOLID_TMP = DTSOLID
                          DIST_OLD(:)=DES_POS_OLD(I,:)-DES_POS_OLD(LL,:)
                          DISTMOD_OLD=SQRT(DES_DOTPRDCT(DIST_OLD,DIST_OLD))
                          NORMAL_OLD(:)=DIST_OLD(:)/DISTMOD_OLD
                          VRN_OLD(:)=DES_VEL_OLD(LL,:)-DES_VEL_OLD(I,:)
                          V_REL_NORM_OLD=DES_DOTPRDCT(VRN_OLD,NORMAL_OLD)
                          WRITE(*,'(5X,A,A,ES15.7)') &
                             'WARNING: normal relative velocity less ',&
                             'than zero at time ', S_TIME
                          WRITE(*,'(7X,A,I10,2X,A,I10,2X,A)') &
                             'for first contact between particles', LL, &
                             'and ', I, 'with'
                          WRITE(*,'(7X,A,ES15.7,2X,A,ES15.7)') &
                             'V_REL_NORM = ', V_REL_TRANS_NORM, &
                             'and V_REL_NORM_OLD = ', V_REL_NORM_OLD
                        ELSE
                           DTSOLID_TMP = OVERLAP_N/&
                              (V_REL_TRANS_NORM+SMALL_NUMBER)
                        ENDIF
                        OVERLAP_T = V_REL_TRANS_TANG*&
                           MIN(DTSOLID,DTSOLID_TMP)
                     ENDIF

                     phaseLL = PIJK(LL,5)                  
                     phaseI = PIJK(I,5)

! T.Li : Hertz vs linear spring-dashpot contact model
                     IF (TRIM(DES_COLL_MODEL) .EQ. 'HERTZIAN') THEN
                        sqrt_overlap = SQRT(OVERLAP_N)
                        KN_DES = hert_kn(phaseLL,phaseI)*sqrt_overlap
                        KT_DES = hert_kt(phaseLL,phaseI)*sqrt_overlap
                        sqrt_overlap = SQRT(sqrt_overlap)
                        ETAN_DES = DES_ETAN(phaseLL,phaseI)*sqrt_overlap
                        ETAT_DES = DES_ETAT(phaseLL,phaseI)*sqrt_overlap
                     ELSE
                        KN_DES = KN
                        KT_DES = KT
                        ETAN_DES = DES_ETAN(phaseLL,phaseI)
                        ETAT_DES = DES_ETAT(phaseLL,phaseI)
                     ENDIF
   
! Calculate the normal contact force                  
                     FNS1(:) = -KN_DES * OVERLAP_N*NORMAL(:)
                     FNS2(:) = -ETAN_DES * V_REL_TRANS_NORM*NORMAL(:)
                     FN(LL,:) = FNS1(:) + FNS2(:)       

! Calculate the tangential displacement which is integration of tangential 
! relative velocity with respect to contact time. Correction in the tangential 
! direction is imposed                   
                     PFT(LL,NI,:) = PFT(LL,NI,:) + OVERLAP_T*TANGENT(:)
                     PFT_TMP(:) = PFT(LL,NI,:)   ! update pft_tmp before it used 
                     PFT_TMP(:) = PFT(LL,NI,:) - &
                        DES_DOTPRDCT(PFT_TMP,NORMAL)*NORMAL(:)
! Calculate the tangential contact force                     
                     FTS1(:) = -KT_DES * PFT_TMP(:)
                     FTS2(:) = -ETAT_DES * V_REL_TRANS_TANG*TANGENT(:)
                     FT(LL,:) = FTS1(:) + FTS2(:) 
   
! Temporary storage of tangential force for reporting                  
                     FT_TMP(:) = FT(LL,:)

! Check for Coulombs friction law and limit the maximum value of the
! tangential force on a particle in contact with another particle
                     CALL CFSLIDE(LL, TANGENT)
                  
! Calculate the total force FC and torque TOW on a particle in a
! particle-particle collision
                     CALL CFFCTOW(LL, I, NORMAL, DISTMOD)

! Save the tangential displacement history with the correction of Coulomb's law
                     IF (PARTICLE_SLIDE) THEN
! Since FT might be corrected during the call to cfslide, the tangental
! displacement history needs to be changed accordingly                  
                        PFT(LL,NI,:) = -( FT(LL,:) - FTS2(:) ) / KT_DES
                     ELSE
                        PFT(LL,NI,:) = PFT_TMP(:)
                     ENDIF

! Reporting info                   
                     IF(DEBUG_DES.AND.LL.EQ.FOCUS_PARTICLE) THEN 
                        IF (.NOT.DES_LOC_DEBUG) THEN
                           DES_LOC_DEBUG = .TRUE.
                           WRITE(*,1000) 
                        ENDIF

                        PRINT*, '     EtaN, EtaT =  ', ETAN_DES, ETAT_DES
                        PRINT*, '     Percent overlap = ', (R_LM - DISTMOD)*100.d0/R_LM
                        PRINT*, '     rad ratio = ', DES_RADIUS(LL)/DES_RADIUS(I)
                        PRINT*, '     FNS1 and FNS2 = ', FNS1(:), FNS2(:)
                        PRINT*, '     PFT = ', PFT(LL,NI,:)
                        PRINT*, '     FORCEST = ', FT(LL,:)
                        PRINT*, '     FORCESN = ', FN(LL,:)
                        PRINT*, '     FORCEST = ', FT(LL,:)
                     ENDIF
   
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
                        PRINT*, 'PN', PN(LL,:)
                     ENDIF
   
                     PARTICLE_SLIDE = .FALSE.

                  ENDIF   ! IF(R_LM - DISTMOD.GT.SMALL_NUMBER) -> resolve collision 

! Return particle to its original position/velocity (if moved/changed)
! since they are not needed for any further calculations
                  IF(DES_LE_BC) THEN
                     DES_POS_NEW(LL,:) = TMP_PART_POS(:)
                     DES_VEL_NEW(LL,:) = TMP_PART_VEL(:)
                  ELSEIF(DES_PERIODIC_WALLS) THEN
                     DES_POS_NEW(LL,:) = TMP_PART_POS(:)
!                     DES_POS_NEW(I,:) = TMP_PART_POS(:)
                  ENDIF

               ENDIF   ! IF (I>LL .AND. PEA(I,1))
            ENDDO   ! DO II = 2, NEIGHBOURS(LL,1)+I
         ENDIF   ! IF(NEIGHBOURS(LL,1).GT.0)
!---------------------------------------------------------------------
! End check particle LL neighbour contacts         

      
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

 1000 FORMAT(5X,'---------- START CALC_FORCE_DES ---------->')
 1001 FORMAT(5X,'<---------- END CALC_FORCE_DES ----------') 

      RETURN
      END SUBROUTINE CALC_FORCE_DES

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!
!  Module name: DES_PERIODIC_NEIGHBOR_CHECK
!
!  Purpose: Calculates the distance between particle LL and a potential
!     neighbor knowing the system contains periodic boundaries 
!
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE DES_PERIODIC_NEIGHBOR_CHECK(LL, I, TMP_PART_POS)


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
! Given particle no id. of particle LL (current target particle) and
! potential neighbor particle I
      INTEGER LL, I

! Stores the current position of particle LL (target particle) before
! any modification required from periodic boundaries
      DOUBLE PRECISION TMP_PART_POS(DIMN)

! Various x, y, z distances between particle LL and I 
      DOUBLE PRECISION DELTA_X, DIST_X, DELTA_Y, DIST_Y, &
                       DELTA_Z, DIST_Z

! System dimensions
      DOUBLE PRECISION LX, LY, LZ 

! local variables for x, y, z position of the particle LL and I
      DOUBLE PRECISION XPOS_LL, YPOS_LL, ZPOS_LL, &
                       XPOS_I, YPOS_I, ZPOS_I

! Check whether to print local debug messages
      INTEGER FLAG_PERIODIC_MOVE
!-----------------------------------------------      
! Functions
!-----------------------------------------------      


!-----------------------------------------------

! store the current position of particle LL 
      TMP_PART_POS(:) = DES_POS_NEW(LL,:)
!      TMP_PART_POS(:) = DES_POS_NEW(I,:)

! assign temporary local variables for quick reference
      LX = EX2 - WX1 !=XLENGTH =XE(IMAX1) - XE(1)
      LY = TY2 - BY1 !=YLENGTH =YN(JMAX1) - YN(1)
      LZ = NZ2 - SZ1 !=ZLENGTH =ZT(KMAX1) - ZT(1)

! assign temporary local variables for manipulation/use
      XPOS_LL = DES_POS_NEW(LL,1)
      YPOS_LL = DES_POS_NEW(LL,2)
      XPOS_I = DES_POS_NEW(I,1)
      YPOS_I = DES_POS_NEW(I,2)
      IF (DIMN.EQ.3) THEN
         ZPOS_LL = DES_POS_NEW(LL,3)
         ZPOS_I = DES_POS_NEW(I,3)
      ENDIF

! initialize      
      FLAG_PERIODIC_MOVE = 0


! check x-position
! --------------------
! the x-distance between particle LL and potential neighbor I
      DELTA_X = XPOS_I - XPOS_LL
      DIST_X = DABS(DELTA_X)
! if x-separation distance greater than 1/2 x-domain length, reposition.
      IF(DIST_X>=0.5d0*LX .AND. DES_PERIODIC_WALLS_X) THEN 
! If particle I is east of particle LL, then shift particle LL by +LX
! (i.e. east). Else if particle I is west of particle LL, then shift
! particle LL by -LX (i.e. west).
         XPOS_LL = XPOS_LL + (DELTA_X/DIST_X)*LX
         XPOS_I = XPOS_I - (DELTA_X/DIST_X)*LX
         FLAG_PERIODIC_MOVE = FLAG_PERIODIC_MOVE + 1
      ENDIF

! check y-position
! --------------------
      DELTA_Y = YPOS_I - YPOS_LL
      DIST_Y = DABS(DELTA_Y)
      IF(DIST_Y>=0.5d0*LY .AND. DES_PERIODIC_WALLS_Y) THEN
         YPOS_LL = YPOS_LL + (DELTA_Y/DIST_Y)*LY
         YPOS_I = YPOS_I - (DELTA_Y/DIST_Y)*LY
         FLAG_PERIODIC_MOVE = FLAG_PERIODIC_MOVE + 1         
      ENDIF

! check z-position
! --------------------
      IF (DIMN .EQ. 3) THEN
         DELTA_Z = ZPOS_I - ZPOS_LL
         DIST_Z = DABS(DELTA_Z)
         IF(DIST_Z>=0.50d0*LZ .AND. DES_PERIODIC_WALLS_Z) THEN
            ZPOS_LL = ZPOS_LL + (DELTA_Z/DIST_Z)*LZ
            ZPOS_I = ZPOS_I - (DELTA_Z/DIST_Z)*LZ
            FLAG_PERIODIC_MOVE = FLAG_PERIODIC_MOVE + 1         
         ENDIF      
      ENDIF

! adjust position of particle LL for any periodic movement
      DES_POS_NEW(LL,1) = XPOS_LL
      DES_POS_NEW(LL,2) = YPOS_LL
      IF (DIMN .EQ. 3) DES_POS_NEW(LL,3) = ZPOS_LL
!       DES_POS_NEW(I,1) = XPOS_I
!       DES_POS_NEW(I,2) = YPOS_I
!       IF (DIMN .EQ. 3) DES_POS_NEW(I,3) = ZPOS_I

! Error reporting
      IF(DEBUG_DES .AND. FLAG_PERIODIC_MOVE >0) THEN
         WRITE(*,'(7X,A,I10,X,I10)') &
            'PARTICLE LL & NEIGHBOR I: ', LL, I
         WRITE(*,'(7X,A,(ES15.7))') &
            'I DES_POS = ', DES_POS_NEW(I,:)
         WRITE(*,'(7X,A,(ES15.7))') &
            'ORIGINAL LL DES_POS = ', TMP_PART_POS(:)
         WRITE(*,'(7X,A,(ES15.7))') &
            'NEW LL DES_POS = ', DES_POS_NEW(LL,:)
      ENDIF


      RETURN

      END SUBROUTINE DES_PERIODIC_NEIGHBOR_CHECK



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!
!  Module name: DES_LEBC_NEIGHBOR_CHECK
!
!  Purpose: Calculates the distance between particle LL and a potential
!     neighbor (particle I) knowing the system contains periodic
!     boundaries 
!
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      SUBROUTINE DES_LEBC_NEIGHBOR_CHECK(LL, I, TMP_PART_POS, &
         TMP_PART_VEL)

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
! Given particle no id. of particle LL (current target particle) and
! potential neighbor particle I
      INTEGER LL, I
! Stores the current position and velocity of particle LL (target 
! particle) before any modification required from LE boundaries
      DOUBLE PRECISION TMP_PART_POS(DIMN), TMP_PART_VEL(DIMN)

! local variables for x, y, z position of the particle LL and I
      DOUBLE PRECISION XPOS_LL, YPOS_LL, ZPOS_LL, &
                       XPOS_I, YPOS_I, ZPOS_I

! local variables for x, y, z velocity of the particle LL
      DOUBLE PRECISION XVEL_LL, YVEL_LL, ZVEL_LL

! Various x, y, z distances between particle LL and I 
      DOUBLE PRECISION DELTA_X, DIST_X, DELTA_Y, DIST_Y, DELTA_Z, DIST_Z

! System dimensions
      DOUBLE PRECISION LX, LY, LZ 

! local variable for shear direction
      CHARACTER*4 SHEAR_DIR      

! local variable for relative velocity of shear
      DOUBLE PRECISION REL_VEL

! the distance between the periodic boundaries corresponding to the
! direction the shear is acting. for du/dy shear this corresponds to the
! x domain length
      DOUBLE PRECISION DOMAIN_SIZE

! determined by first calcaulating the distance the LE boundary (cell) 
! that was originally aligned with the center cell traveled in a given
! time step.  then integer multiples of the domain size are subtracted
! from this quantity until a distance less than the domain size remains
      DOUBLE PRECISION OFFSET_DISTANCE

! Check whether to print local debug messages
      INTEGER FLAG_PERIODIC_MOVE
!-----------------------------------------------      
! Functions
!-----------------------------------------------      


!-----------------------------------------------

! store the current position & velocity of particle LL 
      TMP_PART_POS(:) = DES_POS_NEW(LL,:)
      TMP_PART_VEL(:) = DES_VEL_NEW(LL,:)

! assign temporary local variables for quick reference
      LX = EX2 - WX1 !=XLENGTH =XE(IMAX1) - XE(1)
      LY = TY2 - BY1 !=YLENGTH =YN(JMAX1) - YN(1)
      LZ = NZ2 - SZ1 !=ZLENGTH =ZT(KMAX1) - ZT(1)

      REL_VEL = DES_LE_REL_VEL
      SHEAR_DIR = TRIM(DES_LE_SHEAR_DIR)

! assign temporary local variables for manipulation/use
      XPOS_LL = DES_POS_NEW(LL,1)
      YPOS_LL = DES_POS_NEW(LL,2)
      XPOS_I  = DES_POS_NEW(I,1)
      YPOS_I  = DES_POS_NEW(I,2)
      XVEL_LL = DES_VEL_NEW(LL,1)
      YVEL_LL = DES_VEL_NEW(LL,2)
      IF (DIMN.EQ.3) THEN
         ZPOS_LL = DES_POS_NEW(LL,3)
         ZPOS_I  = DES_POS_NEW(I,3)         
         ZVEL_LL = DES_VEL_NEW(LL,3)
      ENDIF

! initialize
      FLAG_PERIODIC_MOVE = 0
      OFFSET_DISTANCE = ZERO

      IF (DIMN .EQ. 2) THEN

! the x, y distances between particle LL and potential neighbor I
         DELTA_X = XPOS_I - XPOS_LL
         DELTA_Y = YPOS_I - YPOS_LL
         DIST_X = DABS(DELTA_X)      
         DIST_Y = DABS(DELTA_Y)

! 2D shear : du/dy
! ----------------------------------------               
         IF(TRIM(SHEAR_DIR).EQ.'DUDY') THEN
            DOMAIN_SIZE = LX                 
            IF (REL_VEL .NE. ZERO) THEN
               OFFSET_DISTANCE = REL_VEL*S_TIME - DOMAIN_SIZE*&
                  DBLE( FLOOR( REAL( (REL_VEL*S_TIME/DOMAIN_SIZE) )) )
            ENDIF

! if y-separation distance greater than 1/2 y-domain length, reposition.
            IF (DIST_Y >= 0.5d0*LY) THEN
! adjust y-position by moving it one y-domain length.  if particle I is
! above particle LL, then shift particle LL by +Ly (north).  Else if
! particle I is below particle LL, then shift particle LL by -Ly (south).
               YPOS_LL = YPOS_LL + (DELTA_Y/DIST_Y)*LY

! adjust x-position for Lees & Edwards Boundaries
               XPOS_LL = XPOS_LL + (DELTA_Y/DIST_Y)*OFFSET_DISTANCE

! adjust x-velocity for Lees & Edwards Boundaries                  
               XVEL_LL = XVEL_LL + (DELTA_Y/DIST_Y)*REL_VEL

! calculate distance terms with corrected position
               DELTA_X = XPOS_I - XPOS_LL
               DIST_X = DABS(DELTA_X)      
               FLAG_PERIODIC_MOVE = FLAG_PERIODIC_MOVE + 1
            ENDIF

! if x-separation distance greater than 1/2 x-domain length, reposition.
            IF (DIST_X >= 0.5d0*DOMAIN_SIZE) THEN
! adjust x-position by moving it one x-domain length.  
               XPOS_LL = XPOS_LL + (DELTA_X/DIST_X)*DOMAIN_SIZE

! calculate distance terms with corrected position
               DELTA_X = XPOS_I - XPOS_LL
               DIST_X = DABS(DELTA_X)      
               FLAG_PERIODIC_MOVE = FLAG_PERIODIC_MOVE + 1
            ENDIF

! if x-separation distance greater than 1/2 x-domain length, reposition
! for the second and last time (2nd time possibly required due to LE BC).
            IF (DIST_X >= 0.5d0*DOMAIN_SIZE) THEN
! adjust x-position by moving it one x-domain length.  
               XPOS_LL = XPOS_LL + (DELTA_X/DIST_X)*DOMAIN_SIZE
               FLAG_PERIODIC_MOVE = FLAG_PERIODIC_MOVE + 1
            ENDIF

! 2D shear : dv/dx
! ----------------------------------------            
         ELSEIF(TRIM(SHEAR_DIR).EQ.'DVDX') THEN
            DOMAIN_SIZE = LY                 
            IF (REL_VEL .NE. ZERO) THEN
               OFFSET_DISTANCE = REL_VEL*S_TIME - DOMAIN_SIZE*&
                  DBLE( FLOOR( REAL( (REL_VEL*S_TIME/DOMAIN_SIZE) )) )
            ENDIF

! if x-separation distance greater than 1/2 x-domain length, reposition.
            IF (DIST_X >= 0.5d0*LX) THEN
! adjust x-position by moving it one x-domain length.  if particle I is
! above particle LL, then shift particle LL by +Lx (east).  Else if
! particle I is below particle LL, then shift particle LL by -Lx (west).
               XPOS_LL = XPOS_LL + (DELTA_X/DIST_X)*LX

! adjust y-position for Lees & Edwards Boundaries
               YPOS_LL = YPOS_LL + (DELTA_X/DIST_X)*OFFSET_DISTANCE

! adjust y-velocity for Lees & Edwards Boundaries                  
               YVEL_LL = YVEL_LL + (DELTA_X/DIST_X)*REL_VEL

! calculate distance terms with corrected position
               DELTA_Y = YPOS_I - YPOS_LL
               DIST_Y = DABS(DELTA_Y)
               FLAG_PERIODIC_MOVE = FLAG_PERIODIC_MOVE + 1
            ENDIF

! if y-separation distance greater than 1/2 y-domain length, reposition.
            IF (DIST_Y >= 0.5d0*DOMAIN_SIZE) THEN
! adjust y-position by moving it one y-domain length.  
               YPOS_LL = YPOS_LL + (DELTA_Y/DIST_Y)*DOMAIN_SIZE

! calculate distance terms with corrected position
               DELTA_Y = YPOS_I - YPOS_LL
               DIST_Y = DABS(DELTA_Y)
               FLAG_PERIODIC_MOVE = FLAG_PERIODIC_MOVE + 1
            ENDIF

! if y-separation distance greater than 1/2 y-domain length, reposition
! for the second and last time (2nd time possibly required due to LE BC).
            IF (DIST_Y >= 0.5d0*DOMAIN_SIZE) THEN
! adjust x-position by moving it one x-domain length.  
               YPOS_LL = YPOS_LL + (DELTA_Y/DIST_Y)*DOMAIN_SIZE
               FLAG_PERIODIC_MOVE = FLAG_PERIODIC_MOVE + 1
            ENDIF

         ENDIF   ! endif shear_dir == dudy or dvdx
      ENDIF  ! if dimn == 2


      IF (DIMN .EQ. 3) THEN

! the x, y, z distances between particle LL and potential neighbor I
         DELTA_X = XPOS_I - XPOS_LL
         DELTA_Y = YPOS_I - YPOS_LL
         DELTA_Z = ZPOS_I - ZPOS_LL         
         DIST_X = DABS(DELTA_X)      
         DIST_Y = DABS(DELTA_Y)
         DIST_Z = DABS(DELTA_Z)

! 3D shear : du/dy
! ----------------------------------------               
         IF(TRIM(SHEAR_DIR).EQ.'DUDY') THEN
            DOMAIN_SIZE = LX                 
            IF (REL_VEL .NE. ZERO) THEN
               OFFSET_DISTANCE = REL_VEL*S_TIME - DOMAIN_SIZE*&
                  DBLE( FLOOR( REAL( (REL_VEL*S_TIME/DOMAIN_SIZE) )) )
            ENDIF

! if y-separation distance greater than 1/2 y-domain length, reposition.
            IF (DIST_Y >= 0.5d0*LY) THEN
! adjust y-position by moving it one y-domain length.  if particle I is
! above particle LL, then shift particle LL by +Ly (north).  Else if
! particle I is below particle LL, then shift particle LL by -Ly (south).
               YPOS_LL = YPOS_LL + (DELTA_Y/DIST_Y)*LY

! adjust x-position for Lees & Edwards Boundaries
               XPOS_LL = XPOS_LL + (DELTA_Y/DIST_Y)*OFFSET_DISTANCE

! adjust x-velocity for Lees & Edwards Boundaries                  
               XVEL_LL = XVEL_LL + (DELTA_Y/DIST_Y)*REL_VEL

! calculate distance terms with corrected position
               DELTA_X = XPOS_I - XPOS_LL
               DIST_X = DABS(DELTA_X)      
               FLAG_PERIODIC_MOVE = FLAG_PERIODIC_MOVE + 1
            ENDIF

! if x-separation distance greater than 1/2 x-domain length, reposition.
            IF (DIST_X >= 0.5d0*DOMAIN_SIZE) THEN
! adjust x-position by moving it one x-domain length.  
               XPOS_LL = XPOS_LL + (DELTA_X/DIST_X)*DOMAIN_SIZE
! calculate distance terms with corrected position
               DELTA_X = XPOS_I - XPOS_LL
               DIST_X = DABS(DELTA_X)      
               FLAG_PERIODIC_MOVE = FLAG_PERIODIC_MOVE + 1
            ENDIF

! if x-separation distance greater than 1/2 x-domain length, reposition
! for the second and last time (2nd time possibly required due to LE BC).
            IF (DIST_X >= 0.5d0*DOMAIN_SIZE) THEN
! adjust x-position by moving it one x-domain length.  
               XPOS_LL = XPOS_LL + (DELTA_X/DIST_X)*DOMAIN_SIZE
               FLAG_PERIODIC_MOVE = FLAG_PERIODIC_MOVE + 1
            ENDIF

! if z-separation distance greater than 1/2 z-domain length, reposition.
            IF (DIST_Z >= 0.5d0*LZ) THEN
               ZPOS_LL = ZPOS_LL + (DELTA_Z/DIST_Z)*LZ
               FLAG_PERIODIC_MOVE = FLAG_PERIODIC_MOVE + 1
            ENDIF

         ENDIF   ! endif shear_dir == dudy
      ENDIF  ! if dimn == 3


! adjust position/velocity of particle LL for any BC movement
      DES_POS_NEW(LL,1) = XPOS_LL
      DES_POS_NEW(LL,2) = YPOS_LL
      DES_VEL_NEW(LL,1) = XVEL_LL
      DES_VEL_NEW(LL,2) = YVEL_LL
      IF (DIMN .EQ. 3) THEN
         DES_POS_NEW(LL,3) = ZPOS_LL
         DES_VEL_NEW(LL,3) = ZVEL_LL
      ENDIF

! Error reporting
      IF(DEBUG_DES .AND. FLAG_PERIODIC_MOVE >0) THEN
         WRITE(*,'(7X,A,I10,X,I10)') &
            'PARTICLE LL & NEIGHBOR I: ', LL, I
         WRITE(*,'(7X,A,(ES15.7))') &
            'I DES_POS = ', DES_POS_NEW(I,:)
         WRITE(*,'(7X,A,(ES15.7))') &
            'ORIGINAL LL DES_POS = ', TMP_PART_POS(:)
         WRITE(*,'(7X,A,(ES15.7))') &
            'NEW LL DES_POS = ', DES_POS_NEW(LL,:)
         WRITE(*,'(7X,A,(ES15.7))') &
            'ORIGINAL LL DES_VEL = ', TMP_PART_VEL(:)
         WRITE(*,'(7X,A,(ES15.7))') &
            'NEW LL DES_VEL = ', DES_VEL_NEW(LL,:)

      ENDIF

      RETURN

      END SUBROUTINE DES_LEBC_NEIGHBOR_CHECK            
