module softspring_funcs_cutcell
  
  CONTAINS
    
      
      
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
      USE stl 
      USE des_stl_functions 
      Implicit none 
      INTEGER, INTENT(IN), OPTIONAL :: PART_ID
      LOGICAL, INTENT(OUT), OPTIONAL :: OVERLAP_EXISTS
      
      INTEGER I, J,K, LL, II, IW, IDIM, IJK, PC, NF, wall_count 
      DOUBLE PRECISION OVERLAP_N, OVERLAP_T, SQRT_OVERLAP, OVERLAP_PERCENT

      DOUBLE PRECISION V_REL_TRANS_NORM, V_REL_TRANS_TANG, &
      DISTSQ, RADSQ, CLOSEST_PT(DIMN) , FT_TMP(DIMN), DISTSIGN 
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

      IF(test_overlap_and_exit) write(101,'(A, 2(2x,i6), 2x,L1)') 'min, max pip = ', LOC_MIN_PIP, LOC_MAX_PIP,NO_NEIGHBORING_FACET_DES(PIJK(part_id,4))
      
      DO LL = LOC_MIN_PIP, LOC_MAX_PIP 
         IF(LL.EQ.FOCUS_PARTICLE) then 
            DEBUG_DES = .TRUE. 
         else 
            DEBUG_DES = .FALSE. 
         endif

         ! skipping non-existent particles or ghost particles      
         IF(.NOT.PEA(LL,1) .OR. PEA(LL,4)) CYCLE

         IF (NO_NEIGHBORING_FACET_DES(PIJK(LL,4))) CYCLE 
         
         PC = PC + 1 
         
         IF(DEBUG_DES.AND.LL.EQ.FOCUS_PARTICLE) THEN 
            IJK = PIJK(LL,4)
            COUNT_FAC = LIST_FACET_AT_DES(IJK)%COUNT_FACETS
            
            WRITE(*,*) 'NUMBER OF FACETS = ', I_OF(IJK), J_OF(IJK), K_OF(IJK), IJK
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

            IPLUS1  =  MIN( I_CELL + 1, IEND2)
            IMINUS1 =  MAX( I_CELL - 1, ISTART2)

            JPLUS1  =  MIN (J_CELL + 1, JEND2)
            JMINUS1 =  MAX( J_CELL - 1, JSTART2)

            KPLUS1  =  MIN (K_CELL + 1, KEND2)
            KMINUS1 =  MAX( K_CELL - 1, KSTART2)
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
                  
                  CALL ClosestPtPointTriangle(DES_POS_NEW(LL,:), &
                       VERTEX(NF, 1,:), VERTEX(NF, 2,:), VERTEX(NF, 3,:), &
                       CLOSEST_PT(:))
                  
                  DIST(:) = DES_POS_NEW(LL,:) - CLOSEST_PT(:) 
                  DISTSQ = DOT_PRODUCT(DIST, DIST)
                  OVERLAP_N = ZERO
                  OVERLAP_PERCENT = ZERO 
                  
                  IF(DISTSQ .GE. RADSQ) CYCLE !No overlap exists 

                  !Overlap detected 
                  IF(test_overlap_and_exit) then 
                  !for the special case where we only 
                  !want to test if a particle overlaps any of the walls
                     OVERLAP_EXISTS = .true.
                     return 
                  endif
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
1001  FORMAT(/1X,70('*'),/5x, & 
           'From: CALC_FORCE_WITH_WALL_CUTFACE_STL: Error 1001',/5x, & 
           'Particle center on the non fluid side of the facet', /5x, & 
           'Particle and Facets IDs ', 2(2x, i10), /5x, & 
           'Facet info: Vert 1 ', 3(2x, g17.8), /5x, & 
           'Facet info: Vert 2 ', 3(2x, g17.8), /5x, & 
           'Facet info: Vert 3 ', 3(2x, g17.8), /5x, & 
           'Facet info: Norm   ', 3(2x, g17.8), /5x, & 
           'Particle Center    ', 3(2x, g17.8), /5x, & 
           'Closest pt on facet', 3(2x, g17.8), /5x, & 
           'Distance vector    ', 3(2x, g17.8), /5x, & 
           'Distance sign      ', (2x, g17.8), /5x, & 
           'See the stl and vtp files for this offending pair', /5x, & 
           'Exiting.',/, 1X,70('*')/)

      !WRITE(UNIT_LOG, 1001) LL, NF, VERTEX(NF,1,:), &
      !     VERTEX(NF,2,:),VERTEX(NF,3,:),NORM_FACE(NF,:), &
      !     DES_POS_NEW(LL,:), CLOSEST_PT(:), DIST(:), DISTSIGN
      
         
    END SUBROUTINE CALC_FORCE_WITH_WALL_CUTFACE_STL
         
    SUBROUTINE write_this_facet_and_part(FID, PID)                
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
      USE stl 
      USE des_stl_functions 
      Implicit none 
      !facet id and particle id 
      Integer, intent(in) :: fid, pid 
      Integer :: stl_unit, vtp_unit , k
      CHARACTER*100 :: stl_fname, vtp_fname 
      real :: temp_array(3)

      stl_unit = 1001
      vtp_unit = 1002      
      
      WRITE(vtp_fname,'(A,"_OFFENDING_PARTICLE",".vtp")') TRIM(RUN_NAME)
      WRITE(stl_fname,'(A,"_STL_FACE",".stl")') TRIM(RUN_NAME)
      
      open(vtp_unit, file = vtp_fname, form='formatted')
      open(stl_unit, file = stl_fname, form='formatted')

      write(vtp_unit,"(a)") '<?xml version="1.0"?>'
      write(vtp_unit,"(a,es24.16,a)") '<!-- time =',s_time,'s -->'
      write(vtp_unit,"(a,a)") '<VTKFile type="PolyData"',&
           ' version="0.1" byte_order="LittleEndian">'
      write(vtp_unit,"(3x,a)") '<PolyData>'
      write(vtp_unit,"(6x,a,i10.10,a,a)")&
           '<Piece NumberOfPoints="',1,'" NumberOfVerts="0" ',&
           'NumberOfLines="0" NumberOfStrips="0" NumberOfPolys="0">'
      write(vtp_unit,"(9x,a)")&
           '<PointData Scalars="Diameter" Vectors="Velocity">'
      write(vtp_unit,"(12x,a)")&
           '<DataArray type="Float32" Name="Diameter" format="ascii">'
      write (vtp_unit,"(15x,es12.6)") (2*des_radius(pid))
      write(vtp_unit,"(12x,a)") '</DataArray>'

      temp_array = zero 
      temp_array(1:DIMN) = des_vel_new(pid, 1:dimn)
      write(vtp_unit,"(12x,a,a)") '<DataArray type="Float32" ',&
           'Name="Velocity" NumberOfComponents="3" format="ascii">'
      write (vtp_unit,"(15x,3(es13.6,3x))")&
           ((temp_array(k)),k=1,3)
      write(vtp_unit,"(12x,a,/9x,a)") '</DataArray>','</PointData>'
      ! skip cell data
      write(vtp_unit,"(9x,a)") '<CellData></CellData>'

      temp_array = zero 
      temp_array(1:dimn) = des_pos_new(pid, 1:dimn)
      write(vtp_unit,"(9x,a)") '<Points>'
      write(vtp_unit,"(12x,a,a)") '<DataArray type="Float32" ',&
           'Name="Position" NumberOfComponents="3" format="ascii">'
      write (vtp_unit,"(15x,3(es13.6,3x))")&
           ((temp_array(k)),k=1,3)
      write(vtp_unit,"(12x,a,/9x,a)")'</DataArray>','</Points>'
      ! Write tags for data not included (vtp format style)
      write(vtp_unit,"(9x,a,/9x,a,/9x,a,/9x,a)")'<Verts></Verts>',&
           '<Lines></Lines>','<Strips></Strips>','<Polys></Polys>'
      write(vtp_unit,"(6x,a,/3x,a,/a)")&
           '</Piece>','</PolyData>','</VTKFile>'
      
      !Now write the facet info 
      
      write(stl_unit,*)'solid vcg'      
      
      write(stl_unit,*) '   facet normal ', NORM_FACE(FID,1:3)
      write(stl_unit,*) '      outer loop' 
      write(stl_unit,*) '         vertex ', VERTEX(FID,1,1:3)
      write(stl_unit,*) '         vertex ', VERTEX(FID,2,1:3)
      write(stl_unit,*) '         vertex ', VERTEX(FID,3,1:3)
      write(stl_unit,*) '      endloop'
      write(stl_unit,*) '   endfacet'
               
      write(stl_unit,*)'endsolid vcg'

      close(vtp_unit, status = 'keep')
      close(stl_unit, status = 'keep')      
      
    end SUBROUTINE write_this_facet_and_part
    
      SUBROUTINE CFSLIDEWALL2(TANGNT,PARTICLE_SLIDE)
      
      USE discretelement
      USE param1
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(DIMN), INTENT(IN) :: TANGNT
      
      logical PARTICLE_SLIDE
!-----------------------------------------------
! Local Variables
!-----------------------------------------------

      INTEGER :: K
      DOUBLE PRECISION FTMD, FNMD

!-----------------------------------------------      
! Functions
!-----------------------------------------------         
      DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT
      
!----------------------------------------------- 

      PARTICLE_SLIDE = .FALSE.
      FTMD = SQRT(DES_DOTPRDCT(FTAN, FTAN))
      FNMD = SQRT(DES_DOTPRDCT(FNORM,FNORM))

      IF (FTMD.GT.(MEW_W*FNMD)) THEN
         PARTICLE_SLIDE = .TRUE.
         IF(DES_DOTPRDCT(TANGNT,TANGNT).EQ.0) THEN
            FTAN(:) =  MEW_W * FNMD * FTAN(:)/FTMD
         ELSE
            FTAN(:) = -MEW_W * FNMD * TANGNT(:)
         ENDIF
      ENDIF

      IF(DEBUG_DES .AND. PARTICLE_SLIDE) THEN
         WRITE(*,'(7X,A)') &
            'FROM CFSLIDEWALL2.F ---------->'
         WRITE(*,'(9X,A)') 'PARTICLE_SLIDE = T'
         WRITE(*,'(9X,A,2(ES15.7,X))') &
         'FTMD, mu_w*FNMD = ', FTMD, MEW_W*FNMD
         WRITE(*,'(7X,A)') '<----------END CFSLIDEWALL2.F'
      ENDIF

      RETURN
      END SUBROUTINE CFSLIDEWALL2



      SUBROUTINE CFFCTOWALL2(L, NORM, DIST_LI)
     
      USE param1 
      USE discretelement
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  L
      DOUBLE PRECISION, INTENT(IN), DIMENSION(DIMN) ::  NORM

      
! distance between centers of particles
      DOUBLE PRECISION, INTENT(IN) :: DIST_LI

!-----------------------------------------------
! Local variables
!-----------------------------------------------

      DOUBLE PRECISION :: CROSSP(DIMN)

! distance from the contact point to the particle centers 
      DOUBLE PRECISION :: DIST_CL, DIST_CI      

      FC(L,:) = FC(L,:) + FNORM(:) + FTAN(:) 

      
! calculate the distance from the particle center to the wall
      DIST_CL = DIST_LI - DES_RADIUS(L)
! remember that the dist_CL for the particle-wall case is actually the distance 
! from the center of L particle to wall and not the substituted particle as was 
! done previously

      IF(DIMN.EQ.3) THEN 
         CALL DES_CROSSPRDCT(CROSSP, NORM, FTAN)
         TOW(L,:) = TOW(L,:) + DIST_CL*CROSSP(:)
      ELSE 
         CROSSP(1) = NORM(1)*FTAN(2) - NORM(2)*FTAN(1)
         TOW(L,1) =  TOW(L,1) + DIST_CL*CROSSP(1)
      ENDIF 

      RETURN
      END SUBROUTINE CFFCTOWALL2

      SUBROUTINE CFRELVEL2(L, II, VRN, VRT, TANGNT, NORM, DIST_LI, &
                          WALLCONTACT)
      
      USE discretelement
      USE param1
      IMPLICIT NONE     

!-----------------------------------------------
! Local variables
!----------------------------------------------- 
      INTEGER L, II

! marker for particle/wall contact (1=p/w)
      INTEGER WALLCONTACT

      DOUBLE PRECISION TANGNT(DIMN), NORM(DIMN)
      DOUBLE PRECISION TANMOD, VRN, VRT
      DOUBLE PRECISION VRELTRANS(DIMN)
      DOUBLE PRECISION VSLIP(DIMN), &
                       V_ROT(DIMN), OMEGA_SUM(DIMN)

! distance between particles
      DOUBLE PRECISION DIST_LI
! distance from the contact point to the particle centers 
      DOUBLE PRECISION DIST_CL, DIST_CI      

!----------------------------------------------- 
! Functions
!----------------------------------------------- 
      DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT 
!----------------------------------------------- 
      
      VRELTRANS(:) = (DES_VEL_NEW(L,:) - DES_VEL_NEW(II,:))

! calculate the distance from the particle center to the contact point,
! which is taken as the radical line
! dist_ci+dist_cl=dist_li; dist_ci^2+a^2=ri^2;  dist_cl^2+a^2=rl^2       
      DIST_CL = (DIST_LI**2 + DES_RADIUS(L)**2 - DES_RADIUS(II)**2)/&
      (2.d0*DIST_LI)
      DIST_CI = DIST_LI - DIST_CL
      IF(DIMN.EQ.3) THEN
         OMEGA_SUM(:) = OMEGA_NEW(L,:)*DIST_CL + &
         OMEGA_NEW(II,:)*DIST_CI
      ELSE
         OMEGA_SUM(1) = OMEGA_NEW(L,1)*DIST_CL + &
         OMEGA_NEW(II,1)*DIST_CI
         OMEGA_SUM(2) = ZERO
      ENDIF

      CALL DES_CROSSPRDCT(V_ROT, OMEGA_SUM, NORM)

! total relative velocity 
      VRELTRANS(:) =  VRELTRANS(:) + V_ROT(:)

! normal component of relative velocity (scalar)
      VRN = DES_DOTPRDCT(VRELTRANS,NORM)

! slip velocity of the contact point 
! Equation (8) in Tsuji et al. 1992      
      VSLIP(:) =  VRELTRANS(:) - VRN*NORM(:)
      
! the magnitude of the tangential vector      
      TANMOD = SQRT(DES_DOTPRDCT(VSLIP,VSLIP))     
      IF(TANMOD.NE.ZERO) THEN
! the unit vector in the tangential direction  
         TANGNT(:) = VSLIP(:)/TANMOD
      ELSE
         TANGNT(:) = ZERO
      ENDIF

! tangential component of relative surface velocity (scalar)
      VRT  = DES_DOTPRDCT(VRELTRANS,TANGNT)
      
      
      IF(DEBUG_DES) THEN 
         WRITE(*,*) 'IN CFRELVEL2 ---------------------------------'
         
         WRITE(*,'(3(2x,g17.8))') 'VEL LL = ', DES_VEL_NEW(L,:)
         WRITE(*,'(3(2x,g17.8))') 'VEL I = ', DES_VEL_NEW(II,:)
         WRITE(*,'(3(2x,g17.8))') 'OMEGA LL = ',OMEGA_NEW(L,:)
         WRITE(*,'(3(2x,g17.8))') 'OMEGA I = ', OMEGA_NEW(II,:)
         WRITE(*,'(3(2x,g17.8))') 'NORMAL = ', NORM(:)
         WRITE(*,'(3(2x,g17.8))') 'TANGENT = ', TANGNT(:)
         
         WRITE(*,*) 'DIST_CL, DIST_CI = ', DIST_CL, DIST_CI
         
         WRITE(*,'(3(2x,g17.8))') 'VRN, VRT = ', VRN, VRT
         WRITE(*,*) 'OUT OF CFRELVEL2 ---------------------------------'
      ENDIF
      RETURN
      END SUBROUTINE CFRELVEL2
      
      SUBROUTINE CFSLIDE2(TANGNT,PARTICLE_SLIDE)    
      USE discretelement
      USE param1
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(DIMN), INTENT(IN) :: TANGNT
!-----------------------------------------------
! Local Variables
!-----------------------------------------------

      INTEGER :: K
      DOUBLE PRECISION FTMD, FNMD
      LOGICAL PARTICLE_SLIDE

!-----------------------------------------------      
! Functions
!-----------------------------------------------         
      DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT
      
!----------------------------------------------- 


      FTMD = SQRT(DES_DOTPRDCT(FTAN, FTAN))
      FNMD = SQRT(DES_DOTPRDCT(FNORM,FNORM))

      IF (FTMD.GT.(MEW*FNMD)) THEN
         PARTICLE_SLIDE = .TRUE.
         IF(DES_DOTPRDCT(TANGNT,TANGNT).EQ.0) THEN
            FTAN(:) =  MEW * FNMD * FTAN(:)/FTMD
         ELSE
            FTAN(:) = -MEW * FNMD * TANGNT(:)
         ENDIF
      ENDIF

      IF(DEBUG_DES .AND. PARTICLE_SLIDE) THEN
         WRITE(*,'(7X,A)') &
            'FROM CFSLIDE2.F ---------->'
         WRITE(*,'(9X,A)') 'PARTICLE_SLIDE = T'
         WRITE(*,'(9X,A,2(ES15.7,X))')&
         'FTMD, mu*FNMD = ', FTMD, MEW*FNMD
         WRITE(*,'(7X,A)') '<----------END CFSLIDE2.F'
      ENDIF

      RETURN
      END SUBROUTINE CFSLIDE2

      SUBROUTINE CFFCTOW2(L, II,  NORM, DIST_LI)
     
      USE param1 
      USE discretelement
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  L, II
      DOUBLE PRECISION, INTENT(IN), DIMENSION(DIMN) ::  NORM(DIMN)

      
! distance between particles
      DOUBLE PRECISION, INTENT(IN) :: DIST_LI

!-----------------------------------------------
! Local variables
!-----------------------------------------------

      DOUBLE PRECISION :: CROSSP(DIMN)

! distance from the contact point to the particle centers 
      DOUBLE PRECISION :: DIST_CL, DIST_CI      

!-----------------------------------------------

      FC(L,:) = FC(L,:) + FNORM(:) + FTAN(:) 
!      FC(II,:) = FC(II,:) - FNORM(:) - FTAN(:)


! calculate the distance from the particle center to the contact point,
! which is taken as the radical line
! dist_ci+dist_cl=dist_li; dist_ci^2+a^2=ri^2;  dist_cl^2+a^2=rl^2       
      DIST_CL = (DIST_LI**2 + DES_RADIUS(L)**2 - DES_RADIUS(II)**2)/&
         (2.d0*DIST_LI)
      DIST_CI = DIST_LI - DIST_CL

      IF(DIMN.EQ.3) THEN 
         CALL DES_CROSSPRDCT(CROSSP, NORM, FTAN)
         TOW(L,:)  = TOW(L,:)  + DIST_CL*CROSSP(:)
!         TOW(II,:) = TOW(II,:) + DIST_CI*CROSSP(:)
! remember torque is R cross FT, which, compared to I particle, are
! both negative for the J particle.  Therefore, the toqrue, unlike tangential
! and normal contact forces, will be in the same direction for both the 
! particles making the pair 
      ELSE 
         CROSSP(1) = NORM(1)*FTAN(2) - NORM(2)*FTAN(1)
         TOW(L,1)  = TOW(L,1)  + DIST_CL*CROSSP(1)
!        TOW(II,1) = TOW(II,1) + DIST_CI*CROSSP(1)
      ENDIF 


      RETURN
      END SUBROUTINE CFFCTOW2

      SUBROUTINE CFRELVEL_WALL2(L, II, VRN, VRT, TANGNT, NORM, DIST_LI, &
                          WALLCONTACT)
      
      USE discretelement
      USE param1
      IMPLICIT NONE     

!-----------------------------------------------
! Local variables
!----------------------------------------------- 
      INTEGER L, II

! marker for particle/wall contact (1=p/w)
      INTEGER WALLCONTACT

      DOUBLE PRECISION TANGNT(DIMN), NORM(DIMN)
      DOUBLE PRECISION TANMOD, VRN, VRT
      DOUBLE PRECISION VRELTRANS(DIMN)
      DOUBLE PRECISION VSLIP(DIMN), &
                       V_ROT(DIMN), OMEGA_SUM(DIMN)

! distance between particles
      DOUBLE PRECISION DIST_LI
! distance from the contact point to the particle centers 
      DOUBLE PRECISION DIST_CL, DIST_CI      

!----------------------------------------------- 
! Functions
!----------------------------------------------- 
      DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT 
!----------------------------------------------- 

! translational relative velocity 
      VRELTRANS(:) = DES_VEL_NEW(L,:)

! rotational contribution  : v_rot
! calculate the distance from the particle center to the wall
      DIST_CL = DIST_LI         !- DES_RADIUS(L)
      IF(DIMN.EQ.3) THEN
         OMEGA_SUM(:) = OMEGA_NEW(L,:)*DIST_CL
      ELSE
         OMEGA_SUM(1) = OMEGA_NEW(L,1)*DIST_CL
         OMEGA_SUM(2) = ZERO
      ENDIF
      
      CALL DES_CROSSPRDCT(V_ROT, OMEGA_SUM, NORM)

! total relative velocity 
      VRELTRANS(:) =  VRELTRANS(:) + V_ROT(:)

! normal component of relative velocity (scalar)
      VRN = DES_DOTPRDCT(VRELTRANS,NORM)

! slip velocity of the contact point 
! Equation (8) in Tsuji et al. 1992      
      VSLIP(:) =  VRELTRANS(:) - VRN*NORM(:)
      
! the magnitude of the tangential vector      
      TANMOD = SQRT(DES_DOTPRDCT(VSLIP,VSLIP))     
      IF(TANMOD.NE.ZERO) THEN
! the unit vector in the tangential direction  
         TANGNT(:) = VSLIP(:)/TANMOD
      ELSE
         TANGNT(:) = ZERO
      ENDIF

! tangential component of relative surface velocity (scalar)
      VRT  = DES_DOTPRDCT(VRELTRANS,TANGNT)
      
      
      IF(DEBUG_DES) THEN 
         WRITE(*,*) 'IN CFRELVEL_WALL2------------------------------'
         
         WRITE(*,'(3(2x,g17.8))') 'VEL LL = ', DES_VEL_NEW(L,:)
         WRITE(*,'(3(2x,g17.8))') 'OMEGA LL = ',OMEGA_NEW(L,:)
         WRITE(*,'(3(2x,g17.8))') 'NORMAL = ', NORM(:)
         WRITE(*,'(3(2x,g17.8))') 'TANGENT = ', TANGNT(:)
         
         WRITE(*,*) 'DIST_CL, DIST_CI = ', DIST_CL, DIST_CI
         
         WRITE(*,'(3(2x,g17.8))') 'VRN, VRT = ', VRN, VRT
         WRITE(*,*) 'OUT OF CFRELVEL_WALL2 ---------------------------------'
      ENDIF
      RETURN
      END SUBROUTINE CFRELVEL_WALL2

 end module softspring_funcs_cutcell
    
