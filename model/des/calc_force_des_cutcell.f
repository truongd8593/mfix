!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_FORCE_DES_CUTCELL                                 C
!
!  Purpose: DES calculations of force acting on a particle for cutcell C
!           treatment of boundaries. 
!           This routine is still under development. Works well        C
!           for reasonably difficult geometries.                       C
!
!                                                                      C
!  Reviewer: Rahul Garg                             Date: 01-Jan-2012  C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_FORCE_DES_CUTCELL

      USE run
      USE param1
      USE discretelement
      USE geometry
      USE compar
      USE constant
      USE cutcell 
      USE funits
      USE indices
      USE physprop
      USE parallel
      USE softspring_funcs_cutcell
      IMPLICIT NONE

!-----------------------------------------------
! Local variables
!-----------------------------------------------      
      INTEGER I, J, LL, II, IW, IDIM, WALL_COUNT
      INTEGER NI, NLIM, N_NOCON, NEIGH_L
      INTEGER OVERLAP_MAXP

      DOUBLE PRECISION OVERLAP_N, OVERLAP_T, SQRT_OVERLAP
      DOUBLE PRECISION FRAC_OVERLAP1, FRAC_OVERLAP2
! percent of particle radius when excess overlap will be flagged
      DOUBLE PRECISION, PARAMETER :: flag_overlap = 0.20d0

      DOUBLE PRECISION V_REL_TRANS_NORM, V_REL_TRANS_TANG, &
                       V_REL_NORM_OLD, VRN_OLD(DIMN)
      DOUBLE PRECISION PFT_TMP(DIMN), FT_TMP(DIMN)
! local normal and tangential forces      
      DOUBLE PRECISION FNS1(DIMN), FNS2(DIMN)
      DOUBLE PRECISION FTS1(DIMN), FTS2(DIMN)
! tmp variables to calculate magnitude of forces
      DOUBLE PRECISION FTMD, FNMD
      DOUBLE PRECISION NORMAL(DIMN), NORMAL_OLD(DIMN), TANGENT(DIMN)
      DOUBLE PRECISION DIST(DIMN), DIST_OLD(DIMN), DISTMOD, &
                       DISTMOD_OLD, R_LM
      DOUBLE PRECISION DTSOLID_TMP

! logic flag telling whether contact pair is old      
      LOGICAL ALREADY_NEIGHBOURS
! logic flag for local debug warnings
      LOGICAL DES_LOC_DEBUG
      LOGICAL ALREADY_EXISTS
! index to track accounted for particles      
      INTEGER PC
! store solids phase index of particle (i.e. pijk(np,5))
      INTEGER PHASEI, PHASELL, IJK
! local values used spring constants and damping coefficients
      DOUBLE PRECISION ETAN_DES,  ETAT_DES, KN_DES, KT_DES

      Double precision :: OVERLAP_FRAC, OVERLAP_FRAC_MAX_ALLOWED
      
      LOGICAL PARTICLE_SLIDE  
!-----------------------------------------------      
! Functions
!-----------------------------------------------      
      DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT 

!-----------------------------------------------      
      
      INCLUDE 'function.inc'

! Calculate new values
!---------------------------------------------------------------------
      OVERLAP_FRAC_MAX_ALLOWED = 0.05
      OVERLAP_MAXP = UNDEFINED_I
      OVERLAP_MAX = ZERO
      NEIGH_MAX = -1
      !DES_LOC_DEBUG = .true. ;    DEBUG_DES = .true. 
      DES_LOC_DEBUG = .false. ;      DEBUG_DES = .false. 

      FOCUS_PARTICLE = 87
      IF (S_TIME.LE.DTSOLID) THEN
         TANGENT(:) = ZERO
         NORMAL(:) = ZERO
         FC(:,:) = ZERO
      ENDIF

      IF(USE_STL) then
         CALL CALC_FORCE_WITH_WALL_CUTFACE_STL
      ELSE 
         CALL CALC_FORCE_WITH_WALL_CUTFACE
      ENDIF
!     Calculate contact force and torque
!---------------------------------------------------------------------
      
      PC = 1      
      DO LL = 1, MAX_PIP
         
         IF(LL.eq.focus_particle) then 
            !debug_des = .true. 
         else 
            debug_des = .false. 
         endif
         IF (PEA(LL,1)) PC = PC + 1 
         IF(.NOT.PEA(LL,1) .OR. PEA(LL,4)) CYCLE

         

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



! Check particle LL neighbour contacts         
!---------------------------------------------------------------------
         IF (NEIGHBOURS(LL,1).GT.0) THEN
            DO II = 2, NEIGHBOURS(LL,1)+1
               I = NEIGHBOURS(LL,II)

               !IF(I.GT.LL .AND. PEA(I,1)) THEN
               
               IF(PEA(I,1)) THEN                  
   
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
                  
                  
                  R_LM = DES_RADIUS(LL) + DES_RADIUS(I)
                  DIST(:) = DES_POS_NEW(I,:) - DES_POS_NEW(LL,:)
                  DISTMOD = SQRT(DES_DOTPRDCT(DIST,DIST))
                  
                  
                  IF(R_LM - DISTMOD.GT.SMALL_NUMBER) THEN

                     IF(DEBUG_DES .AND. LL.EQ.FOCUS_PARTICLE) THEN
                        IF (.NOT.DES_LOC_DEBUG) THEN
                           DES_LOC_DEBUG = .TRUE.
                           WRITE(*,1000) 
                        ENDIF
                        WRITE(*,'(5X,A,I10)') 'NEIGHBORS: ', NEIGHBOURS(LL,:)
                     ENDIF

                     IF((((R_LM-DISTMOD)/R_LM)*100.d0).GT.OVERLAP_MAX) THEN
                        OVERLAP_MAX = (((R_LM-DISTMOD)/R_LM)*100.d0)
                        OVERLAP_MAXP = LL
                     ENDIF
                     
                     FRAC_OVERLAP1 = (R_LM-DISTMOD)/DES_RADIUS(LL)
                     FRAC_OVERLAP2 = (R_LM-DISTMOD)/DES_RADIUS(I)
                     IF (FRAC_OVERLAP1 > flag_overlap .OR. &
                         FRAC_OVERLAP2 > flag_overlap) THEN
!                        WRITE(*,'(5X,A,A,ES15.7)') &
!                           'WARNING: excessive overlap detected ', &
!                           'at time ', S_TIME
!                        WRITE(*,'(7X,A,I10,2X,A,I5,2X,A)') &
!                           'between particles ', LL, 'and ',&
!                           I, 'with'
!                        WRITE(*,'(7X,A,ES15.7,2X,A,ES15.7,2X,ES15.7)') &
!                           'overlap = ', (R_LM-DISTMOD), &
!                           ' radii = ', DES_RADIUS(LL), DES_RADIUS(I)
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
                     CALL CFRELVEL2(LL, I, V_REL_TRANS_NORM, &
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
                  ELSE
                     GOTO 300
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

                  FNS1(:) = -KN_DES * OVERLAP_N * NORMAL(:)
                  FNS2(:) = -ETAN_DES * V_REL_TRANS_NORM*NORMAL(:)
                  FNORM(:) = FNS1(:) + FNS2(:)       

! Calculate the tangential displacement which is integration of tangential 
! relative velocity with respect to contact time. Correction in the tangential 
! direction is imposed                   
                  PFT(LL,NI,:) = PFT(LL,NI,:) + OVERLAP_T * TANGENT(:)
                  PFT_TMP(:) = PFT(LL,NI,:)   ! for an easy pass to des_dotprdct
                  PFT_TMP(:) = PFT(LL,NI,:) - &
                     DES_DOTPRDCT(PFT_TMP,NORMAL)*NORMAL(:)
                  FTS1(:) = -KT_DES * PFT_TMP(:)
                  FTS2(:) = -ETAT_DES * V_REL_TRANS_TANG * TANGENT(:)
                  FTAN(:) = FTS1(:) + FTS2(:) 
                  
! Check for Coulombs friction law and limit the maximum value of the tangential
! force on a particle in contact with another particle
                  CALL CFSLIDE2(TANGENT, PARTICLE_SLIDE)
                  
! Calculate the total force FC and TOW on a particle in a particle-particle collision
                  CALL CFFCTOW2(LL, I, NORMAL, DISTMOD)

! Save the tangential displacement history with the correction of Coulomb's law
                  IF (PARTICLE_SLIDE) THEN
                     PFT(LL,NI,:) = -( FTAN(:) - FTS2(:) ) / KT_DES
                  ELSE
                     PFT(LL,NI,:) = PFT_TMP(:)
                  ENDIF
                  
                  IF(DEBUG_DES.AND.LL.EQ.FOCUS_PARTICLE) THEN 
                     IF (.NOT.DES_LOC_DEBUG) THEN
                        DES_LOC_DEBUG = .TRUE.
                        WRITE(*,1000) 
                     ENDIF
                     
                     WRITE(*,*)'NORMAL =         ;', NORMAL
                     WRITE(*,*)'TANGNT =         ;', TANGENT
                     PRINT*, '     EtaN, EtaT =  ', ETAN_DES, ETAT_DES
                     PRINT*, '     Percent overlap = ', (R_LM - DISTMOD)*100.d0/R_LM
                     PRINT*, '     rad ratio = ', DES_RADIUS(LL)/DES_RADIUS(I)
                     PRINT*, '     FTS1 and FTS2 = ', FTS1(:), FTS2(:)
                     PRINT*, '     FNS1 and FNS2 = ', FNS1(:), FNS2(:)
                     PRINT*, '     PFT = ', PFT(LL,NI,:)
                     PRINT*, '     FORCEST = ', FTAN(:)
                     PRINT*, '     FORCESN = ', FNORM(:)
                     READ(*,*)
                  ENDIF

                  IF(DEBUG_DES.AND.LL.eq.FOCUS_PARTICLE)THEN
                     INQUIRE(FILE='debug_file',EXIST=ALREADY_EXISTS)
                     IF(ALREADY_EXISTS)THEN
                        OPEN(UNIT=1,FILE='debug_file',STATUS='OLD',POSITION='APPEND')
                        WRITE(1,'(A,I5)')'CALC FORCE -- NEIGHBOR',II
                        WRITE(1,'(2(1x,A,E12.5))')&
                        'FNx=',FNORM(1), 'FNy=',FNORM(2)
                     ELSE
                        OPEN(UNIT=1,FILE='debug_file',STATUS='NEW')
                        WRITE(1,'(A,I5)')'CALC FORCE -- NEIGHBOR',II
                        WRITE(1,'(2(1x,A,E12.5))')&
                        'FNx=',FNORM(1), 'FNy=',FNORM(2)
                     ENDIF
                     CLOSE (1)
                     PRINT*, 'PN', PN(LL,:)
                  ENDIF

                  PARTICLE_SLIDE = .FALSE.

               ENDIF         ! IF (I>LL .AND. PEA(I,1))

 300           CONTINUE
            ENDDO            ! DO II = 2, NEIGHBOURS(LL,1)+I
         ENDIF               ! IF(NEIGHBOURS(LL,1).GT.0)

!---------------------------------------------------------------------
! End check particle LL neighbour contacts         
      ENDDO   ! end loop over paticles LL
      
      
      !IF(ABS(FC(1,1)).GT.ZERO)
! Treats wall interaction also as a two-particle interaction but accounting
! for the wall properties

! Calculate gas-solids drag force on particle 
      IF(DES_CONTINUUM_COUPLED) THEN
         CALL CALC_DES_DRAG_GS
      ENDIF

! Calculate solids-solids drag force on particle 
      IF (DES_CONTINUUM_HYBRID) THEN
         CALL CALC_DES_DRAG_SS
      ENDIF

! Update the old values of particle position and velocity with the new
! values computed
      CALL CFUPDATEOLD

      IF (DES_LOC_DEBUG) WRITE(*,1001)

 1000 FORMAT(5X,'---------- START CALC_FORCE_DES ---------->')
 1001 FORMAT(5X,'<---------- END CALC_FORCE_DES ----------') 

      RETURN
      END SUBROUTINE CALC_FORCE_DES_CUTCELL

      
      SUBROUTINE CALC_FORCE_WITH_WALL_CUTFACE
      USE run
      USE param1
      USE discretelement
      USE geometry
      USE compar
      USE constant
      USE cutcell 
      USE funits
      USE indices
      USE physprop
      USE parallel
      USE softspring_funcs_cutcell
      Implicit none 
      INTEGER I, J, LL, II, IW, IDIM, WALL_COUNT, IJK, COUNT, COUNT_BC, IJK_WALL, PC
      DOUBLE PRECISION OVERLAP_N, OVERLAP_T, SQRT_OVERLAP, FRAC_OVERLAP1, &
      OVERLAP_FRAC

      DOUBLE PRECISION V_REL_TRANS_NORM, V_REL_TRANS_TANG, FT_TMP(DIMN)
! local normal and tangential forces      
      DOUBLE PRECISION FNS1(DIMN), FNS2(DIMN)
      DOUBLE PRECISION FTS1(DIMN), FTS2(DIMN)
      DOUBLE PRECISION NORMAL(DIMN), TANGENT(DIMN), DIST(DIMN), DISTMOD, R_LM, &
      WALL_COOR(DIMN)
      
      LOGICAL DES_LOC_DEBUG, consider_bc, consider_bc_temp, particle_slide 
      
      DOUBLE PRECISION DIR_X, DIR_Y, DIR_Z, DIR_X_ARR(DIMN), XREF, YREF, ZREF
      
      DOUBLE PRECISION X1MINX0(3), X2MINX1(3), TEMP_CROSSP(3), X1(3)

      integer RESTRICT_COUNT, RESTRICT_ARR(10), COUNT_REST, OVERLAP_MAXP, phaseLL

! normal vector components for sending to get_del_h (compute perp distance from
! particle to cut-cell face)
      DOUBLE PRECISION  NORM1, NORM2, NORM3, NORMAL_MAG
      CHARACTER*100 :: WALL_TYPE

! temporary variables for periodic boundaries      
      DOUBLE PRECISION TEMPX, TEMPY, TEMPZ, TEMPD      
! local values used spring constants and damping coefficients
      DOUBLE PRECISION ETAN_DES_W, ETAT_DES_W, KN_DES_W, KT_DES_W
      INCLUDE 'function.inc'

      DES_LOC_DEBUG = .false. ;      DEBUG_DES = .false. 
      OVERLAP_MAXP = UNDEFINED_I

      FOCUS_PARTICLE = -1 

      PC = 1      
      DO LL = 1, MAX_PIP
         
         IF(LL.eq.focus_particle) then 
            !debug_des = .true. 
         else 
            debug_des = .false. 
         endif
         IF (PEA(LL,1)) PC = PC + 1 
         IF(.NOT.PEA(LL,1) .OR. PEA(LL,4)) CYCLE

         
         
         FTS1(:) = ZERO
         FTS2(:) = ZERO
         FNS1(:) = ZERO
         FNS2(:) = ZERO


         ! Check particle LL for wall contacts
!---------------------------------------------------------------------
! Treats wall interaction also as a two-particle interaction but accounting
! for the wall properties; make sure the particle is not classified as
! a new 'entering' particle or is already marked as a potential exiting particle
         IF( .NOT.PEA(LL,2) .AND. .NOT.PEA(LL,3)) THEN
            WALL_COUNT = 0 
            IJK = PIJK(LL,4)
            COUNT_BC = DES_CELLWISE_BCDATA(IJK)%COUNT_DES_BC
            RESTRICT_COUNT = 0
            R_LM = ZERO 
            DISTMOD = ZERO 
            
            
            DO COUNT = 1, COUNT_BC 
               IJK_WALL = DES_CELLWISE_BCDATA(IJK)%BDRY_LIST(COUNT)%IJK_SCAL
               WALL_TYPE = DES_CELLWISE_BCDATA(IJK)%BDRY_LIST(COUNT)%DES_BC_TYPE               
               SELECT CASE (TRIM(WALL_TYPE)) 
               
               CASE('NORMAL_WALL')
                  NORMAL(1:DIMN) = DES_CELLWISE_BCDATA(IJK)%BDRY_LIST(COUNT)%NORMAL(1:DIMN)
                  !perpendicular distance from center of Lth particle
                  WALL_COOR(1) = (NORMAL(1)+ONE)*XLENGTH*HALF
                  !For west wall, normal(1) = -1, therefore, WALL_COOR(1) = ZERO and other coordinates will be zero 
                  !For east wall, normal(1) = 1, thetefore, wall_coor(1) = XLENGTH
                  
                  
                  WALL_COOR(2) = (NORMAL(2)+ONE)*YLENGTH*HALF
                  IF(DIMN.EQ.3) WALL_COOR(3) = (NORMAL(3)+ONE)*ZLENGTH*HALF
                  
                  !perpendicular distance from center of LLth partcle to WALL 
                  
                  DO IDIM  = 1, DIMN
                     DIST(IDIM) = WALL_COOR(IDIM) - DES_POS_NEW(LL,IDIM)
                     DIST(IDIM) = DIST(IDIM)*ABS(NORMAL(IDIM))
                   !this is because the distance between the wall_coor and particle
                     !matters only in the direction of the normal
                  ENDDO
                  
                  DISTMOD = SQRT(DOT_PRODUCT(DIST,DIST))
                  R_LM = DES_RADIUS(LL)
                  
!                  WRITE(*,*) 'WALL_COOR2 = ', WALL_COOR(2), R_LM, DISTMOD
                  !Note that R_LM does not include another radius as was the practice in
                  !the older version of the code. This is because now the wall coordinate has 
                  !not been reduced by the particle radius
                  
                  IW = ABS(NORMAL(1))*1 + ABS(NORMAL(2))*2
                  IF(DIMN.EQ.3) IW = IW +  ABS(NORMAL(3))*3

                  IF(DEBUG_DES.AND.LL.EQ.FOCUS_PARTICLE) THEN 
                     
                     WRITE(*,'(A, 2x, i4, 2x, g17.8)'      ) 'WALL ID, DISTMOD ', IW, DISTMOD
                     WRITE(*,'(A, 3(2x, g17.8))') 'WALL_COR = ', WALL_COOR(:)
                     WRITE(*,'(A, 3(2x, g17.8))') 'DIST   = ', DIST(:)
                     WRITE(*,'(A, 3(2x, g17.8))') 'NORMAL   = ', NORMAL(:)
                  ENDIF
                  
               CASE('CUT_FACE')
                  
                  CONSIDER_BC = .TRUE.
                  !first check if this wall shud be considered or not 
                  !this will only happen if cut_face_line was also 
                  !encountered as a bc earlier. This will happen only for 
                  !particles belonging to cut-cells. 
                  DO COUNT_REST  = 1, RESTRICT_COUNT
                     IF(RESTRICT_ARR(COUNT_REST).EQ.IJK_WALL) THEN
                        WRITE(*,*) 'RESTRICTING A CUT-FACE BECAUSE A CUT-FACE LINE WAS ALREADY ACCOUNTED AS WALL'
                        
                        WRITE(*,*) 'S_TIME = ', S_TIME
                        WRITE(*,*) 'PARTICLE ID = ', LL
                        WRITE(*,*) 'RESTRICT_COUNT = ', RESTRICT_COUNT
                        WRITE(*,*) 'BC ID and #of BCs ', COUNT, COUNT_BC
                        
                        WRITE(*,'(A,3(2x,i4))') 'I, J, K NAT = ', I_OF(IJK), J_OF(IJK), K_OF(IJK)
                        WRITE(*,'(A,3(2x,i4))') 'I, J, K WALL = ', I_OF(IJK_WALL), J_OF(IJK_WALL), K_OF(IJK_WALL)
                        !read(*,*) 
                        CONSIDER_BC = .FALSE.
                        EXIT !THE COUNT_REST LOOP
                     endif
                  ENDDO
                  !Rahul: Jan 24, 2012.
                  !if the particle belongs to non-cutcell then check if
                  !particle is close close enough to this cut-face.
                  !very vague comment above. hopefully, documentation will
                  !do some justice
                  
                  IF(.NOT.CUT_CELL_AT(IJK)) THEN 
                     DIR_X = I_OF(IJK_WALL) - I_OF(IJK)
                     DIR_Y = J_OF(IJK_WALL) - J_OF(IJK)
                     DIR_Z = K_OF(IJK_WALL) - K_OF(IJK)
                     DIR_X_ARR(1) = DIR_X
                     DIR_X_ARR(2) = DIR_Y
                     IF(DIMN.eq.3) DIR_X_ARR(3) = DIR_Z

                     XREF = XG_E(I_OF(IJK)) - 0.5d0*DX(I_OF(IJK)) + DIR_X*0.5d0*DX(I_OF(IJK)) 
                     YREF = YG_N(J_OF(IJK)) - 0.5d0*DY(J_OF(IJK)) + DIR_Y*0.5d0*DY(J_OF(IJK)) 
                     ZREF = ZG_T(K_OF(IJK)) - 0.5d0*DZ(K_OF(IJK)) + DIR_Z*0.5d0*DZ(K_OF(IJK)) 
                     DIST(1) = (XREF - DES_POS_NEW(LL,1))*DIR_X
                     DIST(2) = (YREF - DES_POS_NEW(LL,2))*DIR_Y
                     IF(DIMN.eq.3) DIST(3) = (ZREF - DES_POS_NEW(LL,3))*DIR_Z
                     
                     
                     
                     CONSIDER_BC_TEMP = .TRUE.
                     DO IDIM = 1, DIMN
                        IF(ABS(DIR_X_ARR(IDIM)).GT.ZERO) THEN 
                           
                           CONSIDER_BC_TEMP = CONSIDER_BC_TEMP.and.(DIST(IDIM).LT.DES_RADIUS(LL))
                                                   
                           
                           IF( CONSIDER_BC_TEMP.and.DEBUG_DES) THEN 
                           !IF( CONSIDER_BC_TEMP) THEN 
                              WRITE(*,*) 'IDIM, DIST, DP  = ', IDIM, DIST(IDIM), DES_RADIUS(LL)
                              
                              WRITE(*,*) 'consider_bc_temp = ', consider_bc_temp
                           ENDIF
                        ENDIF
                        
                     ENDDO

                     IF( CONSIDER_BC_TEMP.and.DEBUG_DES) THEN 
                     !IF( CONSIDER_BC_TEMP) THEN 
                   
                     WRITE(*,*) 'IJK NAT = ', IJK, I_OF(IJK), J_OF(IJK)
                     WRITE(*,*) 'IJK BC  = ', IJK_WALL, I_OF(IJK_WALL), J_OF(IJK_WALL)
                     WRITE(*, '(A10, 3(2x, g17.8))')'DX, DY, DZ = ', DX(I_OF(IJK)), DY(J_OF(IJK)), DZ(K_OF(IJK))
                     
                     
                     WRITE(*, '(A10, 3(2x, g17.8))')'XE, YN, ZT = ', XG_E(I_OF(IJK)), YG_N(J_OF(IJK)), ZG_T(K_OF(IJK))
                     
                     WRITE(*,'(A10, 3(2x, g17.8))') 'DIR_X = ', dir_x, dir_y, dir_z
                     WRITE(*, '(A10, 3(2x, g17.8))') 'XREF = ', Xref, yref, zref 
                     WRITE(*, '(A10, 3(2x, g17.8))') 'XPOST = ', (DES_POS_NEW(LL,IDIM), IDIM = 1, DIMN) 
                     
                     WRITE(*, '(A10, 3(2x, g17.8))') 'DIST = ', (DIST(IDIM), IDIM = 1, DIMN) 
                     
                   !  read(*,*)
                     ENDIF

                     CONSIDER_BC = CONSIDER_BC_TEMP
                  ENDIF
                  IF(.NOT.CONSIDER_BC) CYCLE !THE COUNT = 1, COUNT_BC LOOP 
                  
                    
                  TEMPX = DES_POS_NEW(LL,1)
                  TEMPY = DES_POS_NEW(LL,2)
                  TEMPZ = 0.d0 
                  IF(DIMN.EQ.3) TEMPZ = DES_POS_NEW(LL,3)
                  
                  
                  CALL GET_DEL_H_DES(IJK_WALL,'SCALAR',TEMPX, TEMPY, TEMPZ, &
                       &  DISTMOD, NORM1, NORM2, NORM3, .true.)
                  
                  !Remember the normal from get_del_h is from wall to particle. But we need 
                  !the normal pointing from particle toward wall. 
                  NORMAL(1) = -norm1
                  NORMAL(2) = -norm2 
                  IF(DIMN.EQ.3) NORMAL(3) = -norm3 
                  NORMAL_MAG = DOT_PRODUCT(NORMAL, NORMAL)
                  
                  IF(DISTMOD.LT.ZERO) THEN 
                     WRITE(*,1010) mype, IJK, I_OF(IJK), & 
                          & J_OF(IJK), K_OF(IJK), CUT_CELL_AT(IJK), IJK_WALL, I_OF(IJK_WALL), J_OF(IJK_WALL), K_OF(IJK_WALL),  &
                          & TEMPX, TEMPY, TEMPZ
                     WRITE(*,*) 'NORMAL OF CUTCELL = ', NORMAL(1:3)
                     IF(DMP_LOG) WRITE(UNIT_LOG,1010) mype, IJK, I_OF(IJK), & 
                          & J_OF(IJK), K_OF(IJK), CUT_CELL_AT(IJK), IJK_WALL, I_OF(IJK_WALL), J_OF(IJK_WALL), K_OF(IJK_WALL),  &
                          & TEMPX, TEMPY, TEMPZ

            
1010                 FORMAT( 10x,'FROM PROC = ', 2x, i10, /10x, &
                          & 'ERROR IN CALC_FORCE_DES_CUT_CELL:' , /10x, &
                          & 'NEGATIVE DISTANCE:', /10x, &
                          & 'PARTICLE CELL ATTRIBUTES (IJK, I, J, K)', 4(2x,i10) , /10x, &
                          & 'CUT_CELL_AT_PARTICLE CELL ?', 2x, L5 , /10x, &
                          & 'WALL CELL ATTRIBUTES     (IJK, I, J, K)', 4(2x,i10) , /10x, &
                          & 'PARTICLE POSITION CORDINATES', 3(2x,g17.8) , /10x, &
                          & 'TERMINAL ERROR: STOPPING' , /)
                     
                     CALL write_des_data
                     CALL write_VTU_FILE
                     CALL mfix_exit(mype)
                  ENDIF
                  
                  DIST(:) = DISTMOD*NORMAL(:) 
                  R_LM = DES_RADIUS(LL)

                  IW = 4

                  
               CASE('CUT_FACE_LINE')
                  TEMPX = DES_POS_NEW(LL,1)
                  TEMPY = DES_POS_NEW(LL,2)
                  TEMPZ = 0.d0 
                  IF(DIMN.EQ.3) TEMPZ = DES_POS_NEW(LL,3)
                  !in 3-d, minimum distance from point x0 to line is given by 
                  !d = mag({x2-x1} CROSS {x1-x0})/mag({x2-x1})
                  !reference 
                  !http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html

                  !the equation for line in vector form is assumed as 
                  !r = x1+t(x2_x1)
                  !however, in our terminology, we have stored line as 
                  !r = cnot+t vec
                  !therefore, cnot = x1, and vec = x2-x1
                  !therefore, the distance is d = mag{vec CROSS {cnot-X0}}/mag{vec}
                  !and also note that we already normalized vec to a unit vector 
                  X1(:) = DES_CELLWISE_BCDATA(IJK)%BDRY_LIST(COUNT)%CNOT(1:3)
                  X2MINX1(:) = DES_CELLWISE_BCDATA(IJK)%BDRY_LIST(COUNT)%VEC(1:3)
                  X1MINX0(1) = X1(1) - TEMPX
                  X1MINX0(2) = X1(2) - TEMPY
                  X1MINX0(3) = ZERO 
                  IF(DIMN.EQ.3)  X1MINX0(3) = X1(3) - TEMPZ

                  CALL DES_CROSSPRDCT_3D(TEMP_CROSSP(1:3), X2MINX1(1:3), X1MINX0(1:3))
                  
                  DISTMOD = SQRT(DOT_PRODUCT(TEMP_CROSSP(1:3), TEMP_CROSSP(1:3)))
                  
                  !Remember the normal is from line to particle. But we need 
                  !the normal pointing from particle toward wall. 
                  NORMAL(1:DIMN) = - DES_CELLWISE_BCDATA(IJK)%BDRY_LIST(COUNT)%NORMAL(1:DIMN)

                  
                  DIST(:) = DISTMOD*NORMAL(:) 
                  R_LM = DES_RADIUS(LL)

                  IW = 5

               CASE DEFAULT
                  WRITE(*,*)'EROR IN SUBROUTINE CALC_FORCE_DES:'
                  WRITE(*,*)'UNKNOWN TYPE OF DES_WALL_TYPE:',WALL_TYPE
                  WRITE(*,*)'ACCEPTABLE TYPES ARE:' 
                  WRITE(*,*)'NORMAL_WALL' 
                  WRITE(*,*)'CUT_FACE' 
                  WRITE(*,*)'CUT_FACE_LINE' 
                  WRITE(*,*) 'TERMINAL ERROR: STOPPING'
                  CALL mfix_exit(mype)
               END SELECT
               
               I = MAX_PIS + IW
               
               IF(R_LM - DISTMOD.GT.SMALL_NUMBER) THEN 
                  OVERLAP_FRAC = (R_LM - DISTMOD)/R_LM
                     
                  WALL_COUNT = WALL_COUNT+1
!                  WRITE(*,*) 'DETECTED PARTICLE-WALL OVERLAP, IW = ', IW

                  IF(IW.eq.5) THen 
                     !this can only happen if cutcell method is on. So use cutcell structures without any concern for seg errors 
                     IF(CUT_CELL_AT(IJK)) THEN
                        RESTRICT_COUNT = RESTRICT_COUNT +1
                        RESTRICT_ARR(RESTRICT_COUNT )  = IJK_WALL 
                        !for the cut-cell, ijk and ijk_wall shud be same
                     ENDIF
                  ENDIF
                  IF((((R_LM-DISTMOD)/R_LM)*100.d0).GT.OVERLAP_MAX) THEN
                     OVERLAP_MAX = (((R_LM-DISTMOD)/R_LM)*100.d0)
                     OVERLAP_MAXP = LL
                  ENDIF
                                    
                  !Calculate the translational relative velocity for a contacting particle pair
                  CALL CFRELVEL_WALL2(LL, IW, V_REL_TRANS_NORM, &
                  & V_REL_TRANS_TANG, TANGENT, NORMAL, DISTMOD, 1)

! The normal overlap calculation was changed so that it no longer
! depends on the contact history (i.e., integration of incremental
! overlap found by velocity*dtsolid).  Now overlap is based purely on
! position of neighbors.  The change was made because the former
! method was found not to conserve energy 
                  OVERLAP_N =  R_LM-DISTMOD 
                  OVERLAP_T = V_REL_TRANS_TANG*DTSOLID
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
                  
                  FNS1(:) = -KN_DES_W * OVERLAP_N * NORMAL(:)
                  FNS2(:) = -ETAN_DES_W * V_REL_TRANS_NORM * NORMAL(:)
                  FNORM(:) = FNS1(:) + FNS2(:) 
                  
! Calculate the tangential displacement which is integration of tangential 
! relative velocity with respect to contact time. Correction in the tangential 
! direction is imposed                   
                  FTS1(:) = -KT_DES_W * OVERLAP_T*TANGENT(:)
                  FTS2(:) = -ETAT_DES_W * V_REL_TRANS_TANG * TANGENT(:)
                  FTAN(:) = FTS1(:) + FTS2(:) 
                  
                  
                  FT_TMP(:) = FTAN(:)

! Check for Coulombs friction law and limit the maximum value of the
! tangential force on a particle in contact with a wall
                  CALL CFSLIDEWALL2(TANGENT, PARTICLE_SLIDE)
                  
! Calculate the total force FC and TOW on a particle in a particle-wall
! collision
                  CALL CFFCTOWALL2(LL, NORMAL, DISTMOD)
                  
        
                  
                  PARTICLE_SLIDE = .FALSE.
                  
               ENDIF            !wall contact
 200           CONTINUE
            ENDDO
         ENDIF                  !if(walldtsplit .and. .not.pea(LL,2))

      ENDDO
!---------------------------------------------------------------------
! End check particle LL for wall contacts         

      END SUBROUTINE CALC_FORCE_WITH_WALL_CUTFACE
      

      SUBROUTINE CALC_FORCE_WITH_WALL_CUTFACE_STL(PART_ID, OVERLAP_EXISTS)
      
      USE run
      USE param1
      USE discretelement
      USE geometry
      USE compar
      USE constant
      USE cutcell 
      USE funits
      USE indices
      USE physprop
      USE parallel
      USE softspring_funcs_cutcell
      USE stl 
      USE des_stl_functions 
      Implicit none 
      INTEGER, INTENT(IN), OPTIONAL :: PART_ID
      LOGICAL, INTENT(OUT), OPTIONAL :: OVERLAP_EXISTS
      
      INTEGER I, J,K, LL, II, IW, IDIM, IJK, PC, NF, wall_count 
      DOUBLE PRECISION OVERLAP_N, OVERLAP_T, SQRT_OVERLAP, OVERLAP_PERCENT

      DOUBLE PRECISION V_REL_TRANS_NORM, V_REL_TRANS_TANG, &
      DISTSQ, RADSQ, CLOSEST_PT(DIMN) , FT_TMP(DIMN)
! local normal and tangential forces      
      DOUBLE PRECISION FNS1(DIMN), FNS2(DIMN)
      DOUBLE PRECISION FTS1(DIMN), FTS2(DIMN)
      DOUBLE PRECISION NORMAL(DIMN), TANGENT(DIMN), DIST(DIMN), DISTMOD
      
      LOGICAL :: CONTACT_ALREADY_FACET(DIM_STL),DES_LOC_DEBUG, PARTICLE_SLIDE, &
      test_overlap_and_exit
      INTEGER :: COUNT_FAC, COUNT, COUNT2, list_of_cont_facets(20), &
      contact_facet_count, NEIGH_CELLS, NEIGH_CELLS_NONNAT, &
      LIST_OF_CELLS(27), CELL_ID, I_CELL, J_CELL, K_CELL, cell_count , & 
      IMINUS1, IPLUS1, JMINUS1, JPLUS1, KMINUS1, KPLUS1, PHASELL, LOC_MIN_PIP, &
      LOC_MAX_PIP 

! local values used spring constants and damping coefficients
      DOUBLE PRECISION ETAN_DES_W, ETAT_DES_W, KN_DES_W, KT_DES_W      
      INCLUDE 'function.inc'
      
      CONTACT_ALREADY_FACET = .false. 
      DES_LOC_DEBUG = .false. ;      DEBUG_DES = .false. 
      FOCUS_PARTICLE = -1 
      
      !When sent from main routine the loop shud go from 1, max_pip
      !adding the capability to test for a specified particle 
      LOC_MIN_PIP = 1
      LOC_MAX_PIP = MAX_PIP
      if(present(part_id)) then
         LOC_MIN_PIP = part_id
         LOC_MAX_PIP = part_id 
      endif
      test_overlap_and_exit = .false.
      IF(present(OVERLAP_EXISTS)) then 
         test_overlap_and_exit = .true.
      endif

      PC = 1      
      
      DO LL = LOC_MIN_PIP, LOC_MAX_PIP 
         
         IF (NO_NEIGHBORING_FACET_DES(PIJK(LL,4))) CYCLE 

         IF(LL.EQ.FOCUS_PARTICLE) then 
                                !debug_des = .true. 
         else 
            DEBUG_DES = .FALSE. 
         endif
         IF (PEA(LL,1)) PC = PC + 1 
         IF(.NOT.PEA(LL,1) .OR. PEA(LL,4)) CYCLE
         
            
         IF(DEBUG_DES.AND.LL.EQ.FOCUS_PARTICLE) THEN 
            IJK = PIJK(LL,4)
            COUNT_FAC = LIST_FACET_AT_DES(IJK)%COUNT_FACETS
            
            WRITE(*,*) 'NUMBER OF FACETS = ', COUNT_FAC, I_OF(IJK), J_OF(IJK), K_OF(IJK)
            
            WRITE(*,'(A, 3(2x, g17.8))') 'POS = ', DES_POS_NEW(LL, :)
         ENDIF

         FTS1(:) = ZERO
         FTS2(:) = ZERO
         FNS1(:) = ZERO
         FNS2(:) = ZERO


! Check particle LL for wall contacts
!---------------------------------------------------------------------
! Treats wall interaction also as a two-particle interaction but accounting
! for the wall properties; make sure the particle is not classified as
! a new 'entering' particle or is already marked as a potential exiting particle
         IF( .NOT.PEA(LL,2) .AND. .NOT.PEA(LL,3)) THEN
            
            LIST_OF_CELLS(:) = -1
            NEIGH_CELLS = 0
            NEIGH_CELLS_NONNAT  = 0
            WALL_COUNT = 0 
            CELL_ID = PIJK(LL,4)
            COUNT_FAC = LIST_FACET_AT_DES(CELL_ID)%COUNT_FACETS
            RADSQ = DES_RADIUS(LL)*DES_RADIUS(LL)
            
            IF (COUNT_FAC.gt.0)   then 
               NEIGH_CELLS = NEIGH_CELLS + 1
               LIST_OF_CELLS(NEIGH_CELLS) = CELL_ID 
            ENDIF
            
            I_CELL = I_OF(CELL_ID)
            J_CELL = J_OF(CELL_ID)
            K_CELL = K_OF(CELL_ID) 

            IPLUS1   =  MIN( I_CELL + 1, IEND2)
            IMINUS1 =  MAX(I_CELL - 1, ISTART2)

            JPLUS1   =  MIN (J_CELL + 1, JEND2)
            JMINUS1 =  MAX(J_CELL - 1, JSTART2)

            KPLUS1   =  MIN (K_CELL + 1, KEND2)
            KMINUS1 =  MAX(K_CELL - 1, KSTART2)
           ! WRITE(*,*) '---------------------------------------------------'
           ! WRITE(*,'(A10, 4(2x,i5))') 'PCELL  = ', CELL_ID, I_CELL, J_CELL, K_CELL
           ! WRITE(*,'(A10, (2x,i5), 2(2x,g17.8))') '# of FACETS = ', COUNT_FAC, des_pos_new(LL,2), des_vel_new(ll,2)

            DO K = KMINUS1, KPLUS1
               DO J = JMINUS1, JPLUS1
                  DO I = IMINUS1, IPLUS1
                     IJK = FUNIJK(I,J,K)
                     COUNT_FAC = LIST_FACET_AT_DES(IJK)%COUNT_FACETS
                     IF(COUNT_FAC.EQ.0) CYCLE 
                     distsq = zero 
                     IF(DES_POS_NEW( LL , 1) > XG_E(I)) DISTSQ = DISTSQ &
                     + (DES_POS_NEW(LL,1)-XG_E(I))*(DES_POS_NEW(LL,1)-XG_E(I))

                     IF(DES_POS_NEW( LL , 1) < XG_E(I) - DX(I)) DISTSQ = DISTSQ &
                     + (XG_E(I) - DX(I) - DES_POS_NEW(LL,1))*(XG_E(I) - DX(I) - DES_POS_NEW(LL,1))

                     IF(DES_POS_NEW( LL , 2) > YG_N(J)) DISTSQ = DISTSQ &
                     + (DES_POS_NEW(LL,2)-YG_N(J))* (DES_POS_NEW(LL,2)-YG_N(J))

                     IF(DES_POS_NEW( LL , 2) < YG_N(J) - DY(J)) DISTSQ = DISTSQ &
                     + (YG_N(J) - DY(J) - DES_POS_NEW(LL,2))* (YG_N(J) - DY(J) - DES_POS_NEW(LL,2))

                     IF(DES_POS_NEW( LL , 3) > ZG_T(K)) DISTSQ = DISTSQ &
                     + (DES_POS_NEW(LL,3)-ZG_T(K))*(DES_POS_NEW(LL,3)-ZG_T(K))

                     IF(DES_POS_NEW( LL , 3) < ZG_T(K) - DZ(K)) DISTSQ = DISTSQ &
                     + (ZG_T(K) - DZ(K) - DES_POS_NEW(LL,3))*(ZG_T(K) - DZ(K) - DES_POS_NEW(LL,3))
                     IF (DISTSQ < RADSQ) then
                        NEIGH_CELLS_NONNAT = NEIGH_CELLS_NONNAT + 1 
                        NEIGH_CELLS = NEIGH_CELLS + 1
                        LIST_OF_CELLS(NEIGH_CELLS) = IJK
                        !WRITE(*,'(A10, 4(2x,i5))') 'WCELL  = ', IJK, I,J,K
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
            
            CONTACT_FACET_COUNT = 0 

            DO CELL_COUNT = 1, NEIGH_CELLS
               IJK = LIST_OF_CELLS(CELL_COUNT)
               
               DO COUNT = 1, LIST_FACET_AT_DES(IJK)%COUNT_FACETS 
                  NF = LIST_FACET_AT_DES(IJK)%FACET_LIST(COUNT)
                  IF(CONTACT_ALREADY_FACET(NF)) CYCLE 
                  
                  CALL ClosestPtPointTriangle(DES_POS_NEW(LL,:), VERTEX(NF, 1,:), &
                       VERTEX(NF, 2,:), VERTEX(NF, 3,:), CLOSEST_PT(:))
                                   
                  DIST(:) = CLOSEST_PT(:) - DES_POS_NEW(LL,:)
                  DISTSQ = DOT_PRODUCT(DIST, DIST)
                  OVERLAP_N = ZERO
                  OVERLAP_PERCENT = ZERO 
                  IF (DISTSQ < RADSQ) THEN 
                     !Overlap detected 
                     DISTMOD = SQRT(DISTSQ)
                     OVERLAP_N = DES_RADIUS(LL) - DISTMOD
                     OVERLAP_PERCENT = (OVERLAP_N/DES_RADIUS(LL))*100.D0
                     CONTACT_ALREADY_FACET(NF) = .TRUE.
                     CONTACT_FACET_COUNT = CONTACT_FACET_COUNT + 1 
                     LIST_OF_CONT_FACETS(CONTACT_FACET_COUNT) = NF
                     !WRITE(*, '(A10, 2x,i5, 5(2x,g17.8))') 'overlap with NF',NF, overlap_n, overlap_percent 

                     NORMAL(:) = -NORM_FACE(NF,:)

                     !Calculate the translational relative velocity for a contacting particle pair
                     CALL CFRELVEL_WALL2(LL, IW, V_REL_TRANS_NORM, &
                     & V_REL_TRANS_TANG, TANGENT, NORMAL, DISTMOD, 1)

! The normal overlap calculation was changed so that it no longer
! depends on the contact history (i.e., integration of incremental
! overlap found by velocity*dtsolid).  Now overlap is based purely on
! position of neighbors.  The change was made because the former
! method was found not to conserve energy 
                     OVERLAP_T = V_REL_TRANS_TANG*DTSOLID
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
                     
                     FNS1(:) = -KN_DES_W * OVERLAP_N * NORMAL(:)
                     FNS2(:) = -ETAN_DES_W * V_REL_TRANS_NORM * NORMAL(:)
                     FNORM(:) = FNS1(:) + FNS2(:) 
                     
! Calculate the tangential displacement which is integration of tangential 
! relative velocity with respect to contact time. Correction in the tangential 
! direction is imposed                   
                     FTS1(:) = -KT_DES_W * OVERLAP_T*TANGENT(:)
                     FTS2(:) = -ETAT_DES_W * V_REL_TRANS_TANG * TANGENT(:)
                     FTAN(:) = FTS1(:) + FTS2(:) 
                     
                  
                     FT_TMP(:) = FTAN(:)

! Check for Coulombs friction law and limit the maximum value of the
! tangential force on a particle in contact with a wall
                     CALL CFSLIDEWALL2(TANGENT, PARTICLE_SLIDE)
                  
! Calculate the total force FC and TOW on a particle in a particle-wall
! collision
                     CALL CFFCTOWALL2(LL, NORMAL, DISTMOD)
                     
                  END IF
                  
               ENDDO
              
            end DO

            !!if (CONTACT_FACET_COUNT.gt.1) read(*,*) 
            
            DO COUNT2 = 1, CONTACT_FACET_COUNT
               CONTACT_ALREADY_FACET(LIST_OF_CONT_FACETS(COUNT2)) = .FALSE. 
            ENDDO
         END IF

      ENDDO
!---------------------------------------------------------------------
! End check particle LL for wall contacts         

         
      END SUBROUTINE CALC_FORCE_WITH_WALL_CUTFACE_STL
         
               
