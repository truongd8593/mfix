
      SUBROUTINE DES_WALLBC_PREPROCSSING 
      
      USE param1
      USE funits
      USE run
      USE compar      
      USE discretelement
      USE mfix_pic
      USE cutcell
      
      USE indices
      USE physprop
      USE parallel
      USE geometry
      USE bc
      USE funits 
      
      IMPLICIT NONE 
      INTEGER :: I, II, IPLUS1, IMINUS1 ! X-coordinate loop indices
      INTEGER :: J, JJ, JPLUS1, JMINUS1   ! Y-coordinate loop indices
      INTEGER :: K, KK, KPLUS1, KMINUS1   ! Z-coordinate loop indices

      INTEGER :: IJK, IJK_WALL, IJK2, COUNT_BC, L, INEW, JNEW, KNEW, IJK_OLD
      INTEGER :: I1, I2, J1, J2, K1, K2
      INTEGER :: FLUID_IND, COUNT, INTCELL_IND, CUTCELL_IND, BCELL_IND
      LOGICAL :: ADD_GHOST_CELL_WALL, LOOP_OVER
      DOUBLE PRECISION :: WALL_NORM(DIMN)
      DOUBLE PRECISION :: XCOR(4), YCOR(4), ZCOR(4)
      DOUBLE PRECISION :: NORM_DIST(4), NORM_VEC(4, 3),  NORMAL_TEMP(3)
      CHARACTER*100 :: WALL_TYPE,FILENAME
      LOGICAL :: REPLACED_BCTYPE,  ADD_THIS_CUT_FACE

      INTEGER :: CELL_ID, IJK_TEMP
      CHARACTER*100 :: WALL_TYPE_TEMP

      INCLUDE 'function.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'fun_avg2.inc'


      ALLOCATE(DES_CELLWISE_BCDATA(DIMENSION_3))
      
      DO IJK = IJKSTART3, IJKEND3
         ALLOCATE(DES_CELLWISE_BCDATA(IJK)%BDRY_LIST(MAX_DES_BC_CELL ))
         
         DES_CELLWISE_BCDATA(IJK)%COUNT_DES_BC = 0
      ENDDO

      !First do it for x=0 plane 

      !Check if this processor contains the first cell
      IF(ISTART1.EQ.IMIN1) THEN 
         I = ISTART1
         DO K = KSTART1, KEND1 
            DO J = JSTART1, JEND1
               ADD_GHOST_CELL_WALL = .TRUE.
               IJK  = FUNIJK(I,J,K)
               IJK_WALL = FUNIJK(IMIN2, J, K)
               WALL_NORM = ZERO
               IF(FLUID_AT(IJK)) THEN 
                  IF(CARTESIAN_GRID) THEN 
                     IF(CUT_CELL_AT(IJK))  THEN 

! If cartesian grid and cut_cell, set it to false for now
! and then later check if this ghost cell wall should also 
! be classfied as a wall. This is done in routine CHECK_IF_GHOST_CELL_WALL_NEEDED
! where all the vertices of the ghost cell wall are tested for lying inside 
! of the cut-face. If any of the vertices lie inside of the cut-face, then
! the ghost cell wall should also be included as a wall 

! Also remember that a cut-cell will always take 
! precedence in this case and the normal-wall condition 
! will only be included in list of BC's if the cut-face 
! leaves a part of this wall in the fluid domain 
                        
                     
                        ADD_GHOST_CELL_WALL = .FALSE.
                        XCOR(1:4) = XG_E(I) - DX(I)
                        YCOR(1) = YG_N(J) - DY(J)
                        ZCOR(1) = ZG_T(K)
                        
                        YCOR(2) = YG_N(J)
                        ZCOR(2) = ZG_T(K)
                        
                        YCOR(3) = YG_N(J) - DY(J)
                        ZCOR(3) = ZG_T(K) - DZ(K) 
                        
                        YCOR(4) = YG_N(J)
                        ZCOR(4) = ZG_T(K) - DZ(K)
                        
                        CALL CHECK_IF_GHOST_CELL_WALL_NEEDED(IJK, IJK_WALL, & 
                        XCOR(1:4), YCOR(1:4), ZCOR(1:4), NORM_DIST(1:4), &
                        NORM_VEC(1:4,1:3), ADD_GHOST_CELL_WALL)
                        
                        WALL_TYPE = 'CUT_FACE'
                        CALL DES_ADD_WALLBC_TO_CELL(IJK, IJK, WALL_NORM, WALL_TYPE)
                     ENDIF
                  ENDIF

                  IF(ADD_GHOST_CELL_WALL) THEN                     
                     WALL_TYPE = 'NORMAL_WALL'
                     WALL_NORM(1) = -ONE
                     IF(DES_PERIODIC_WALLS_X) WALL_TYPE = 'PERIODIC_WALL'
                     CALL DES_ADD_WALLBC_TO_CELL(IJK, IJK_WALL, WALL_NORM, WALL_TYPE)
                  ENDIF
               ENDIF
               
            ENDDO
         ENDDO
      ENDIF
      
      !Now do it for x=Lx plane 

      !Check if this processor contains the last physical cell in the x- direction
      IF(IEND1.EQ.IMAX1) THEN 
         I = IMAX1
         DO K = KSTART1, KEND1 
            DO J = JSTART1, JEND1
               
               ADD_GHOST_CELL_WALL = .TRUE.
               IJK  = FUNIJK(I,J,K)
               IJK_WALL = FUNIJK(IMAX2, J, K)
               WALL_NORM = ZERO
               IF(FLUID_AT(IJK)) THEN 
                  IF(CARTESIAN_GRID) THEN 
                     
                     IF(CUT_CELL_AT(IJK))  THEN 
                        ADD_GHOST_CELL_WALL = .FALSE.
                        XCOR(1:4) = XG_E(I)
                        
                        YCOR(1) = YG_N(J) - DY(J)
                        ZCOR(1) = ZG_T(K)
                        
                        YCOR(2) = YG_N(J)
                        ZCOR(2) = ZG_T(K)
                        
                        YCOR(3) = YG_N(J) - DY(J)
                        ZCOR(3) = ZG_T(K) - DZ(K) 
                        
                        YCOR(4) = YG_N(J)
                        ZCOR(4) = ZG_T(K) - DZ(K)
                        
                        CALL CHECK_IF_GHOST_CELL_WALL_NEEDED(IJK, IJK_WALL, & 
                        XCOR(1:4), YCOR(1:4), ZCOR(1:4), NORM_DIST(1:4), &
                        NORM_VEC(1:4,1:3), ADD_GHOST_CELL_WALL)
                        
                        WALL_TYPE = 'CUT_FACE'
                        
                        CALL DES_ADD_WALLBC_TO_CELL(IJK, IJK, WALL_NORM, WALL_TYPE)
                     ENDIF
                  ENDIF

                  IF(ADD_GHOST_CELL_WALL) THEN
                     WALL_TYPE = 'NORMAL_WALL'
                     WALL_NORM(1) = ONE
                     
                     IF(DES_PERIODIC_WALLS_X) WALL_TYPE = 'PERIODIC_WALL'
                                !Wall normal points away from the fluid 
                     
                     CALL DES_ADD_WALLBC_TO_CELL(IJK, IJK_WALL, WALL_NORM, WALL_TYPE)
                  ENDIF
               ENDIF
               
            ENDDO
         ENDDO
      ENDIF

      
     !Now do it for y = 0 plane 

      !Check if this processor contains the first cell
      IF(JSTART1.EQ.JMIN1) THEN 
         J = JMIN1
         DO K = KSTART1, KEND1 
            DO I = ISTART1, IEND1
               ADD_GHOST_CELL_WALL = .TRUE.
               IJK  = FUNIJK(I,J,K)
               IJK_WALL = FUNIJK(I, JMIN2, K)
               WALL_NORM = ZERO
               IF(FLUID_AT(IJK)) THEN 
                  IF(CARTESIAN_GRID) THEN 

                     IF(CUT_CELL_AT(IJK))  THEN 
                        
                        ADD_GHOST_CELL_WALL = .FAlSE.
                        YCOR(1:4) = YG_N(J) - DY(J)
                        
                        XCOR(1) = XG_E(I) - DX(I)
                        ZCOR(1) = ZG_T(K)
                        
                        XCOR(2) = XG_E(I) 
                        ZCOR(2) = ZG_T(K)
                        
                        XCOR(3) = XG_E(I) - DX(I)
                        ZCOR(3) = ZG_T(K) - DZ(K) 
                        
                        XCOR(4) = XG_E(I) 
                        ZCOR(4) = ZG_T(K) - DZ(K)
                        CALL CHECK_IF_GHOST_CELL_WALL_NEEDED(IJK, IJK_WALL, & 
                        XCOR(1:4), YCOR(1:4), ZCOR(1:4), NORM_DIST(1:4), &
                        NORM_VEC(1:4,1:3), ADD_GHOST_CELL_WALL)
                        
                        WALL_TYPE = 'CUT_FACE'
                        
                        CALL DES_ADD_WALLBC_TO_CELL(IJK, IJK, WALL_NORM, WALL_TYPE)
                     ENDIF

                  ENDIF
                  !Wall normal points away from the fluid 
                  IF(ADD_GHOST_CELL_WALL) THEN
                     WALL_TYPE = 'NORMAL_WALL'
                     WALL_NORM(2) = -ONE
                     IF(DES_PERIODIC_WALLS_Y) WALL_TYPE = 'PERIODIC_WALL'
                     CALL DES_ADD_WALLBC_TO_CELL(IJK, IJK_WALL, WALL_NORM, WALL_TYPE)
                  ENDIF
               ENDIF
               
            ENDDO
         ENDDO
      ENDIF


      !Now do it for y = Ly plane 
      !Check if this processor contains the last physical cell in the y- direction

      IF(JEND1.EQ.JMAX1) THEN 
         J = JMAX1
         DO K = KSTART1, KEND1 
            DO I = ISTART1, IEND1
               
               ADD_GHOST_CELL_WALL = .TRUE.
               IJK  = FUNIJK(I,J,K)
               IJK_WALL = FUNIJK(I, JMAX2, K)
               WALL_NORM = ZERO
               IF(FLUID_AT(IJK)) THEN 
                  IF(CARTESIAN_GRID) THEN 
                     IF(CUT_CELL_AT(IJK)) THEN 

                        ADD_GHOST_CELL_WALL = .FAlSE.
                        YCOR(1:4) = YG_N(J)
                        
                        XCOR(1) = XG_E(I) - DX(I)
                        ZCOR(1) = ZG_T(K)
                        
                        XCOR(2) = XG_E(I) 
                        ZCOR(2) = ZG_T(K)
                        
                        XCOR(3) = XG_E(I) - DX(I)
                        ZCOR(3) = ZG_T(K) - DZ(K) 
                        
                        XCOR(4) = XG_E(I) 
                        ZCOR(4) = ZG_T(K) - DZ(K)
                        
                        CALL CHECK_IF_GHOST_CELL_WALL_NEEDED(IJK, IJK_WALL, XCOR(1:4), YCOR(1:4), ZCOR(1:4), NORM_DIST(1:4), NORM_VEC(1:4,1:3), ADD_GHOST_CELL_WALL)                     
                        WALL_TYPE = 'CUT_FACE'
                                                
                        CALL DES_ADD_WALLBC_TO_CELL(IJK, IJK, WALL_NORM, WALL_TYPE)

                     ENDIF
                  ENDIF
                  !Wall normal points away from the fluid 
                  IF(ADD_GHOST_CELL_WALL) THEN
                     WALL_TYPE = 'NORMAL_WALL'
                     WALL_NORM(2) = ONE
                     
                     IF(DES_PERIODIC_WALLS_Y) WALL_TYPE = 'PERIODIC_WALL'
                     CALL DES_ADD_WALLBC_TO_CELL(IJK, IJK_WALL, WALL_NORM, WALL_TYPE)
                  ENDIF
               ENDIF
               
            ENDDO
         ENDDO
      ENDIF

     !Now do it for z = 0 plane 

      !Check if this processor contains the first cell
      IF(KSTART1.EQ.KMIN1.AND.DIMN.EQ.3) THEN  
         K = KMIN1 
         DO J = JSTART1, JEND1 
            DO I = ISTART1, IEND1
               
               ADD_GHOST_CELL_WALL = .TRUE.
               IJK  = FUNIJK(I,J,K)
               IJK_WALL = FUNIJK(I, J, KMIN2)
               WALL_NORM = ZERO
               IF(FLUID_AT(IJK)) THEN 
                  IF(CARTESIAN_GRID) THEN 
                     IF(CUT_CELL_AT(IJK)) THEN 
                        
                        ADD_GHOST_CELL_WALL = .FALSE.
                        
                        ZCOR(1:4) = ZG_T(K) - DZ(K)
                                                
                        YCOR(1) = YG_N(J) - DY(J)
                        XCOR(1) = XG_E(I) 

                        YCOR(2) = YG_N(J)
                        XCOR(2) = XG_E(I) 
                        
                        YCOR(3) = YG_N(J) - DY(J)
                        XCOR(3) = XG_E(I) - DX(I)

                        YCOR(4) = YG_N(J)
                        XCOR(4) = XG_E(I) - DX(I)
                        
                        CALL CHECK_IF_GHOST_CELL_WALL_NEEDED(IJK, IJK_WALL, XCOR(1:4), YCOR(1:4), ZCOR(1:4), NORM_DIST(1:4), NORM_VEC(1:4,1:3), ADD_GHOST_CELL_WALL)                     
                        WALL_TYPE = 'CUT_FACE'
                                                
                        CALL DES_ADD_WALLBC_TO_CELL(IJK, IJK, WALL_NORM, WALL_TYPE)


                     ENDIF
                                !need to have a special treatment here... 
                  ENDIF
                  WALL_TYPE = 'NORMAL_WALL'
                  WALL_NORM(3) = -ONE
                  !Wall normal points away from the fluid 
                  
                  IF(DES_PERIODIC_WALLS_Z) WALL_TYPE = 'PERIODIC_WALL'
                  IF(ADD_GHOST_CELL_WALL) THEN
                     CALL DES_ADD_WALLBC_TO_CELL(IJK, IJK_WALL, WALL_NORM, WALL_TYPE)
                  ENDIF
               ENDIF
               
            ENDDO
         ENDDO
      ENDIF


     !Now do it for z = Lz plane 

      !Check if this processor contains the last cell in the z- direction
      IF(KEND1.EQ.KMAX1.AND.DIMN.EQ.3) THEN  
         K = KMAX1
         DO J = JSTART1, JEND1 
            DO I = ISTART1, IEND1
               
               ADD_GHOST_CELL_WALL = .TRUE.
               IJK  = FUNIJK(I,J,K)
               IJK_WALL = FUNIJK(I, J, KMAX2)
               WALL_NORM = ZERO
               IF(FLUID_AT(IJK)) THEN 
                  IF(CARTESIAN_GRID) THEN 
                     IF(CUT_CELL_AT(IJK))  THEN 
                        
                        ADD_GHOST_CELL_WALL = .FALSE.
                        ZCOR(1:4) = ZG_T(K) 
                                                
                        YCOR(1) = YG_N(J) - DY(J)
                        XCOR(1) = XG_E(I) 

                        YCOR(2) = YG_N(J)
                        XCOR(2) = XG_E(I) 
                        
                        YCOR(3) = YG_N(J) - DY(J)
                        XCOR(3) = XG_E(I) - DX(I)

                        YCOR(4) = YG_N(J)
                        XCOR(4) = XG_E(I) - DX(I)
                        
                        CALL CHECK_IF_GHOST_CELL_WALL_NEEDED(IJK, IJK_WALL, XCOR(1:4), YCOR(1:4), ZCOR(1:4), NORM_DIST(1:4), NORM_VEC(1:4,1:3), ADD_GHOST_CELL_WALL)                     
                        WALL_TYPE = 'CUT_FACE'
                         CALL DES_ADD_WALLBC_TO_CELL(IJK, IJK, WALL_NORM, WALL_TYPE)

                     ENDIF
                                !need to have a special treatment here... 
                  ENDIF
                  WALL_TYPE = 'NORMAL_WALL'
                  WALL_NORM(3) = ONE
                  IF(DES_PERIODIC_WALLS_Z) WALL_TYPE = 'PERIODIC_WALL'
                  !Wall normal points away from the fluid 
                  IF(ADD_GHOST_CELL_WALL) THEN
                     CALL DES_ADD_WALLBC_TO_CELL(IJK, IJK_WALL, WALL_NORM, WALL_TYPE)
                  ENDIF
               ENDIF
               
            ENDDO
         ENDDO
      ENDIF

!      RETURN 
!     Now check for the cells not adjacent to the ghost cells. This will only 
!     happen for cartesian grid 

      
      CG: IF(CARTESIAN_GRID) THEN 
         DO K = KSTART2, KEND2
            DO J = JSTART2, JEND2
               DO I = ISTART2, IEND2
                  IJK = funijk(i,j,k)
                  !IF(SMALL_CELL_AT(IJK).AND.(.NOT.CUT_CELL_AT(IJK))) THEN
                  !   WRITE(*, '(A, 3(2x,L1))') 'SMALL NE CUT', SMALL_CELL_AT(IJK), CUT_CELL_AT(IJK), FLUID_AT(IJK)
                  !   WRITE(*,'(A, 5(2x,i4))' ) 'IJK =' , IJK, I, J, K
                  !   READ(*,*) 
                  !ENDIF
                  !do not add the small cell for the soft spring case 
                  IF(SMALL_CELL_AT(IJK).and.MPPIC) CALL ADD_SMALL_CELL_BC(IJK)

                  IF(.NOT.FLUID_AT(IJK)) CYCLE 
                  LOOP_OVER = .FALSE.

                  !no need to again check the cells adjacent to ghost cells
                  IF(I.EQ.IMIN1.OR.I.EQ.IMAX1.OR.J.EQ.JMIN1.OR.J.EQ.JMAX1) LOOP_OVER=.TRUE.
                  IF((K.EQ.KMIN1.OR.K.EQ.KMAX1).AND.(DIMN.EQ.3)) LOOP_OVER=.TRUE.
                  !but if IJK is a cut cell then check again for the surrounding cut cells
                  LOOP_OVER = LOOP_OVER.AND..NOT.CUT_CELL_AT(IJK)
                  
                  IF(LOOP_OVER) CYCLE
                  
                  !FLUID: IF(FLUID_AT(IJK)) THEN 
                  !Now check for surrounding 9 scalar cells for 2-D or 27 scalar cells for 3-d to determine which one of them is a cut-cell. Add any surrounding cut-cell to the list of walls for this cell. 
                  IF(CUT_CELL_AT(IJK)) THEN 
                     WALL_TYPE = 'CUT_FACE'
                     !First check the cell itself 
                     !If the cell is a cut-cell then add it first 
                     CALL DES_ADD_WALLBC_TO_CELL(IJK, IJK, WALL_NORM, WALL_TYPE)
                  ENDIF

                  
                  IF(DIMN.EQ.2) THEN 
                     KPLUS1 = K
                     KMINUS1 = K
                  ELSE
                     KPLUS1 =  MIN(K+1, KEND3) !in serial kend3 = kend2. therefore, the min and max ops
                     KMINUS1 = MAX(K-1, KSTART3)
                  ENDIF
                  
                  JPLUS1 =  MIN(J+1, JEND3)
                  JMINUS1 = MAX(J-1, JSTART3)

                  IPLUS1 =  MIN(I+1, IEND3)
                  IMINUS1 = MAX(I-1, ISTART3)

                  DO KK = KMINUS1, KPLUS1
                     DO JJ = JMINUS1, JPLUS1
                        DO II = IMINUS1, IPLUS1 
                           IJK2 = FUNIJK(II,JJ,KK)
                           IF(CUT_CELL_AT(IJK2).and.FLUID_AT(IJK2)) THEN 
                              WALL_TYPE = 'CUT_FACE'
                              CALL CHECK_IF_TOADD_THIS_CUT_FACE(IJK, IJK2, WALL_NORM, WALL_TYPE, ADD_THIS_CUT_FACE)
                              IF(ADD_THIS_CUT_FACE) CALL DES_ADD_WALLBC_TO_CELL(IJK, IJK2, WALL_NORM, WALL_TYPE)
                           ENDIF
                        ENDDO
                     ENDDO
                  ENDDO

                  COUNT_BC  = DES_CELLWISE_BCDATA(IJK)%COUNT_DES_BC
                  IF(COUNT_BC.GE.2.AND.CUT_CELL_AT(IJK).AND.(.NOT.MPPIC)) CALL REORDER_BCS_FOR_CUTCELL(IJK)
                  !IF(COUNT_BC.GE.2.and.(.NOT.CUT_CELL_AT(IJK))) CALL CHECK_NON_CUTCELLBC(IJK)
                  !IF(COUNT_BC.GE.2) WRITE(*,*) 'BC GE 2', COUNT_BC, CUT_CELL_AT(IJK)
                  
               ENDDO
            ENDDO
         ENDDO
      ENDIF CG

      DO L = 1, DIMENSION_BC 
         IF (BC_DEFINED(L)) THEN 
            IF (BC_TYPE(L)=='NO_SLIP_WALL' .OR. BC_TYPE(L)=='FREE_SLIP_WALL'&
                 .OR. BC_TYPE(L)=='PAR_SLIP_WALL') THEN 
               !do nothing 
            ELSEIF (BC_TYPE(L)=='P_OUTFLOW' ) THEN 
               !CYCLE 
               IF(.not.BC_APPLY_TO_MPPIC(L)) THEN 
                  IF(DMP_LOG) write(unit_log,'(2x,A,/2x,A,2x,i4,2x,A,/2x,A)') & 
                  'FROM des_wallbc_preprocessing', 'Pressure outflow  BC # ', L, &
                  ' will not be applied to the discrete phase', &
                  'i.e., particles will not exit the domain from this boundary plane'

                  IF(PRINT_DES_SCREEN) & 
                  write(*,'(2x,A,/2x,A,2x,i4,2x,A,/2x,A)') & 
                  'FROM des_wallbc_preprocessing', 'Pressure outflow  BC # ', L, &
                  ' will not be applied to the discrete phase', & 
                  'i.e., particles will not exit the domain from this boundary plane'
                  
                  !read(*,*) 
                  CYCLE
               ENDIF
               I1 = BC_I_W(L) 
               I2 = BC_I_E(L) 
               J1 = BC_J_S(L) 
               J2 = BC_J_N(L) 
               K1 = BC_K_B(L) 
               K2 = BC_K_T(L) 
               DO K = K1, K2 
                  DO J = J1, J2 
                     DO I = I1, I2  
                        IF (.NOT.IS_ON_myPE_owns(I, J, K)) CYCLE
                        INEW = I
                        JNEW = J 
                        KNEW = K 
                        SELECT CASE (TRIM(BC_PLANE(L)))  
                        CASE ('E') 
                           INEW = I+1
                        CASE ('W')  
                           INEW = I-1 
                        CASE ('N')  
                           JNEW = J+1
                        CASE ('S')  
                           JNEW = J-1
                        CASE ('T')  
                           KNEW = K+1
                        CASE ('B')  
                           KNEW = K-1
                        END SELECT 

! Rahul:
! confusing ? I know! A little explanation on the logic here. 
! If a pressure outflow BC is specified, then it is important 
! to flag fluid cells from where particles can exit. Recall that 
! earlier all fluid (emphasis on fluid) cells next to ghost cells
! were set up as normal walls for particle reflection. Now need to
! look at the specified exits and RE-flag (emphasis on RE) the 
! fluid cells appropriately. 
! in the mfix convention, if the northmost face (Y=ylength) is 
! specified 
! as PO, then J1 = J2 = IMAX+2 or IMAX2, which is the ghost cell.
! Since the BC plane is at the bottom of this ghost cell, BC_PLane
! is set as 'S' in the native mfix pre-processing. 

                        IJK = FUNIJK(INEW,JNEW,KNEW) 
! This is the cell where the new bc will go 

                        IF (.NOT.IS_ON_myPE_wobnd(INEW, JNEW, KNEW)) CYCLE
! For periodic BC's, the east, bottom or southmost cells (i.e., i=1 
! or j = 1, or k = 1, will be flagged as fluid
! and the following fluid_at(IJK) test alone will not 
! be enough to avoid replacement error below.

! However, this IJK could also fall in a non-fluid cell for cut-cell and there
! will be an error in finding the replacement. 
! Therefore, if IJK is not a fluid cell then do not bother and cycle.

                        IF(.not.FLUID_AT(IJK)) CYCLE 
                        !WRITE(*,*) 'BC_PLANE  = ', TRIM(BC_PLANE(L)), J1, J2
                        
                        IJK_WALL = FUNIJK(I,J,K)  !this is the id of the new bc 
                        
                        REPLACED_BCTYPE = .FALSE.
                        BCLOOP : DO COUNT  = 1, DES_CELLWISE_BCDATA(IJK)%COUNT_DES_BC
                           IJK_OLD = DES_CELLWISE_BCDATA(IJK)%BDRY_LIST(COUNT)%IJK_SCAL
                           IF(IJK_OLD.EQ.IJK_WALL) THEN 
                              !then overwrite the des_bc_type as outflow 
                              DES_CELLWISE_BCDATA(IJK)%BDRY_LIST(COUNT)%DES_BC_TYPE = 'OUTFLOW'
                              REPLACED_BCTYPE = .TRUE.
                              EXIT BCLOOP
                           end IF
                        end DO BCLOOP
                        
                        IF(.NOT.REPLACED_BCTYPE) THEN
                                          
                           IF(DMP_LOG) THEN 
                              WRITE(UNIT_LOG, *) 'ERROR IN DES_CUT_CELL_PRE_PROCESSOR'
                              WRITE(UNIT_LOG, *) 'COULD NOT FIND AN OLDER WALL TO REPLACE FOR OUTFLOW BC CASE'
                              WRITE(UNIT_LOG, *) 'LOOKING FOR ', IJK_WALL, ' AS AN EXISTING BC FOR CELL ID', IJK
                              
                              WRITE(*,*) 'ERROR WITH FINDING THE EXISTING WALL TO REPLACE FOR OUTFLOW CASE'
                              WRITE(*,*) 'SEE THE LOG FILE'
                  
                              
                              call mfix_exit(myPE)

                           end IF
                           
                        end IF
                        
                        !WRITE(*,*) 'I, J, K, WALL = ', I, J, K, TRIM(BC_PLANE(L))
                        !WRITE(*,*) 'I, J, K, CELL = ', INEW, JNEW, KNEW
                     end DO
                  end DO
               end DO
               
            ELSEIF (BC_TYPE(L)=='MASS_INFLOW' ) THEN 
               CYCLE 
               !dont add the MI to des_bc information 
               !this is implemented directly in mppic_walbc
               I1 = BC_I_W(L) 
               I2 = BC_I_E(L) 
               J1 = BC_J_S(L) 
               J2 = BC_J_N(L) 
               K1 = BC_K_B(L) 
               K2 = BC_K_T(L) 
               DO K = K1, K2 
                  DO J = J1, J2 
                     DO I = I1, I2  
                        IF (.NOT.IS_ON_myPE_owns(I, J, K)) CYCLE
                        !IS_ON_myPE_owns is true for istart<i<iend  and therefore includes 
                        !the ghost cells as well. 
                        INEW = I
                        JNEW = J 
                        KNEW = K 
                        WALL_NORM(:) = ZERO 
                        SELECT CASE (TRIM(BC_PLANE(L)))  
                        CASE ('E') 
                           INEW = I+1
                           WALL_NORM(1) = -ONE
                        CASE ('W')  
                           INEW = I-1 
                           WALL_NORM(1) = ONE
                        CASE ('N')  
                           JNEW = J+1
                           WALL_NORM(2) = -ONE
                        CASE ('S')  
                           JNEW = J-1
                           WALL_NORM(2) = ONE
                        CASE ('T')  
                           KNEW = K+1
                           WALL_NORM(3) = -ONE
                        CASE ('B')  
                           KNEW = K-1
                           WALL_NORM(3) = ONE
                        END SELECT 
                        IJK = FUNIJK(INEW,JNEW,KNEW) !this is the cell where the new bc will go 
                        IF(.NOT.FLUID_AT(IJK)) CYCLE 
                        IJK_WALL = FUNIJK(I,J,K)  !this is the id of the new bc 
                        COUNT_BC = DES_CELLWISE_BCDATA(IJK)%COUNT_DES_BC 
                        COUNT_BC = COUNT_BC + 1
                        DES_CELLWISE_BCDATA(IJK)%COUNT_DES_BC = COUNT_BC 
                        
                        IF(COUNT_BC.GT.MAX_DES_BC_CELL) THEN
                           IF(DMP_LOG) THEN 
                              WRITE(UNIT_LOG, *) 'ERROR IN DES_CUT_CELL_PRE_PROCESSOR'
                              
                              WRITE(UNIT_LOG,*) 'INCREASE MAX_DES_BC_CELL from the current value of', MAX_DES_BC_CELL
                              WRITE(*,*) 'SEE THE LOG FILE'
                  
                              call mfix_exit(myPE)
                              
                           ENDIF
                        end IF
                        
                        DES_CELLWISE_BCDATA(IJK)%BDRY_LIST(COUNT_BC)%IJK_SCAL = IJK_WALL
                        
                        ALLOCATE(DES_CELLWISE_BCDATA(IJK)%BDRY_LIST(COUNT_BC)%NORMAL(DIMN))
                        
                        DES_CELLWISE_BCDATA(IJK)%BDRY_LIST(COUNT_BC)%NORMAL(:) = WALL_NORM(:)

                        DES_CELLWISE_BCDATA(IJK)%BDRY_LIST(COUNT_BC)%DES_BC_TYPE = 'INFLOW'

                        ALLOCATE(DES_CELLWISE_BCDATA(IJK)%BDRY_LIST(COUNT_BC)%MI_BCID(1))
                        DES_CELLWISE_BCDATA(IJK)%BDRY_LIST(COUNT_BC)%MI_BCID(1) = L 

                        !WRITE(*,*) 'I, J, K, WALL = ', I, J, K, TRIM(BC_PLANE(L)), L
                        !WRITE(*,*) 'I, J, K, CELL = ', INEW, JNEW, KNEW
                     end DO
                  end DO
               end DO
            ELSE
               
               IF(DMP_LOG) THEN 
                  WRITE(UNIT_LOG, *) 'ERROR IN DES_CUT_CELL_PRE_PROCESSOR'
                  WRITE(UNIT_LOG, *) 'THIS WALL BC', BC_TYPE(L),'  IS NOT SUPPORTED FOR MP-PIC' 
                  WRITE(*,*) 'ERROR WITH WALL BC SPECIFICATION, SEE THE LOG FILE'
                  
                  call mfix_exit(myPE)
               end IF
               
            end IF
         end IF
      end DO
      
      !
      
      write(filename,'(A,"_",A,"_",I5.5,".dat")') TRIM(RUN_NAME), 'DES_WALL_BC',myPE
      OPEN(1001, file = TRIM(filename), form ='formatted')
      
      DO K=KSTART2, KEND2
         DO J=JSTART2, JEND2
            DO I=ISTART2, IEND2
               CELL_ID = FUNIJK(I,J,K)
               COUNT_BC = DES_CELLWISE_BCDATA(CELL_ID)%COUNT_DES_BC 

               IF(COUNT_BC.eq.0) CYCLE 

               WRITE(1001, '("**************************************************")')

               write(1001, '(2X, "CELL IJK, I, J, K =        = ", i20, 2x, 4(2x,i10))') CELL_ID, I_OF(CELL_ID), J_OF(CELL_ID), K_OF(CELL_ID)
               IF(CARTESIAN_GRID) write(1001, '(2x, "FLUID, CUTCELL, SMALL CELL ? ", 3(2x, L4))') FLUID_AT(CELL_ID), CUT_CELL_AT(CELL_ID), SMALL_CELL_AT(CELL_ID)
               WRITE(1001, '(2x, "TOTAL BCs                  = ", 3(2x, i10))') COUNT_BC
               DO COUNT = 1, COUNT_BC 
         
                  IJK_TEMP = DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT)%IJK_SCAL 
                  WALL_TYPE_TEMP = DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT)%DES_BC_TYPE
                  NORMAL_TEMP = ZERO 
                  IF(TRIM(WALL_TYPE_TEMP).eq.'CUT_FACE') NORMAL_TEMP(1:3) = NORMAL_S(IJK_TEMP,1:3)
                  IF(TRIM(WALL_TYPE_TEMP).eq.'NORMAL_WALL') NORMAL_TEMP(1:DIMN) = DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT)%NORMAL(1:DIMN) 
                  IF(TRIM(WALL_TYPE_TEMP).eq.'CUT_FACE_LINE') NORMAL_TEMP(1:DIMN) = DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT)%NORMAL(1:DIMN) 
                  
                  SELECT CASE (TRIM(WALL_TYPE_TEMP))
                  CASE('NORMAL_WALL')
                     WRITE(1001, '("----------------------------------------------------")')
                     WRITE(1001, '(2x, "BC #                 = ", 3(2x, i10))') COUNT
                     WRITE(1001, '(2x, "BC TYPE              = ", (2x, A20))') TRIM(WALL_TYPE_TEMP)
                     WRITE(1001, '(2x, "BC IJK, I, J, K      = ", i20, 2x, 4(2x,i10))') IJK_TEMP, I_OF(IJK_TEMP), J_OF(IJK_TEMP), K_OF(IJK_TEMP)
                     WRITE(1001, '(2x, "BC NORMAL            = ", 3(2x, g17.8))') NORMAL_TEMP(1:DIMN)

                     !WRITE(1001, '("----------------------------------------------------")')
                     
                  CASE('CUT_FACE') 
                     
                     WRITE(1001, '("----------------------------------------------------")')
                     WRITE(1001, '(2x, "BC #                 = ", 3(2x, i10))') COUNT
                     WRITE(1001, '(2x, "BC TYPE              = ", (2x, A20))') TRIM(WALL_TYPE_TEMP)
                     WRITE(1001, '(2x, "BC IJK, I, J, K      = ", i20, 2x, 4(2x,i10))') IJK_TEMP, I_OF(IJK_TEMP), J_OF(IJK_TEMP), K_OF(IJK_TEMP)
                     WRITE(1001, '(2x, "FLUID, CUTCELL, SMALL CELL ? ", 3(2x, L4))') FLUID_AT(IJK_TEMP), CUT_CELL_AT(IJK_TEMP), SMALL_CELL_AT(IJK_TEMP)
                     WRITE(1001, '(2x, "BC NORMAL            = ", 3(2x, g17.8))') NORMAL_TEMP(1:DIMN)
                     
                     !WRITE(1001, '("----------------------------------------------------")')
                     
                     
                  CASE('CUT_FACE_LINE') 
                     
                     WRITE(1001, '("----------------------------------------------------")')
                     WRITE(1001, '(2x, "BC #                 = ", 3(2x, i10))') COUNT
                     WRITE(1001, '(2x, "BC TYPE              = ", (2x, A20))') TRIM(WALL_TYPE_TEMP)
                     WRITE(1001, '(2x, "BC IJK, I, J, K      = ", i20, 2x, 4(2x,i10))') IJK_TEMP, I_OF(IJK_TEMP), J_OF(IJK_TEMP), K_OF(IJK_TEMP)
                     WRITE(1001, '(2x, "FLUID, CUTCELL, SMALL CELL ? ", 3(2x, L4))') FLUID_AT(IJK_TEMP), CUT_CELL_AT(IJK_TEMP), SMALL_CELL_AT(IJK_TEMP)
                     WRITE(1001, '(2x, "BC NORMAL            = ", 3(2x, g17.8))') NORMAL_TEMP(1:DIMN)
                     WRITE(1001, '(2x, "CNOT                 = ", 3(2x, g17.8))') DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT)%CNOT(1:3)
                     WRITE(1001, '(2x, "VEC                  = ", 3(2x, g17.8))') DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT)%VEC(1:3)

                     !WRITE(1001, '("----------------------------------------------------")')

                  CASE DEFAULT
                  end SELECT
                  
               ENDDO
            end DO
         end DO
      end DO
      
      close(1001, status = 'keep')

      
!      RETURN 
      write(filename,'(A,"_",I5.5,".dat")') 'CUTCELL_IDEN',myPE
      OPEN(1000, file = TRIM(filename), form ='formatted')
      write(1000,*)'VARIABLES= ',' "X" ',' "Y" ',' "Z" ',' "CUCELL" ', ' "INTCELL" ', ' "FLUID" ', ' "BC" ', ' "VOL" ' 
      write(1000,*)'ZONE F=POINT, I=', (IEND2-ISTART2)+1,  ', J=', JEND2-JSTART2+1, ', K=', KEND2-KSTART2 + 1

      IF(CARTESIAN_GRID) THEN 
      
         DO K=KSTART2, KEND2
            DO J=JSTART2, JEND2
               DO I=ISTART2, IEND2
                  IJK  = FUNIJK(I,J,K)
                  IF(FLUID_AT(IJK)) THEN 
                     FLUID_IND = 1
                  ELSE 
                     FLUID_IND = 0
                  END IF
                  
                  IF(INTERIOR_CELL_AT(IJK)) THEN 
                     INTCELL_IND = 1
                  ELSE 
                     INTCELL_IND = 0
                  ENDIF
                  IF(CUT_CELL_AT(IJK)) THEN 
                     CUTCELL_IND = 1
                  ELSE 
                     CUTCELL_IND = 0
                  ENDIF
                  
                  IF(BLOCKED_CELL_AT(IJK)) THEN 
                     BCELL_IND = 1
                  ELSE
                     BCELL_IND = 0 
                  ENDIF 
                  
!                  write(1000,'(3(2x,g17.8),4(2x,i4), 2x, g17.8)') XE(I),YN(J),ZT(K), CUTCELL_IND, INTCELL_IND ,FLUID_IND, BCELL_IND, VOL(IJK) ! DES_CELLWISE_BCDATA(IJK)%COUNT_DES_BC!,  (DES_CELLWISE_BCDATA(IJK)%COUNT_DES_BC)*1
                  write(1000,'(3(2x,i10),4(2x,i4), 2x, g17.8)') (I),(J),(K), CUTCELL_IND, INTCELL_IND ,FLUID_IND, BCELL_IND, VOL(IJK) ! DES_CELLWISE_BCDATA(IJK)%COUNT_DES_BC!,  (DES_CELLWISE_BCDATA(IJK)%COUNT_DES_BC)*1
               enddo
            enddo
         enddo
         
      ELSE 
            
            DO K=KSTART2, KEND2
               DO J=JSTART2, JEND2
                  DO I=ISTART2, IEND2
                  IJK  = FUNIJK(I,J,K)
                  IF(FLUID_AT(IJK)) THEN 
                     FLUID_IND = 1
                  ELSE 
                     FLUID_IND = 0
                  END IF
                  
                  write(1000,'(3(2x,g17.8), 3(2x,i4), 2x, g17.8)') XE(I-1),YN(J-1),ZT(K-1), 0,FLUID_IND, DES_CELLWISE_BCDATA(IJK)%COUNT_DES_BC,  VOL(IJK)
               ENDDO
            ENDDO
         ENDDO
      ENDIF
      close(1000, status = 'keep')
 200  continue

      return 
      
      CELL_ID = funijk(6,189,7)
      COUNT_BC = DES_CELLWISE_BCDATA(CELL_ID)%COUNT_DES_BC 

      
      DO COUNT = 1, COUNT_BC 
         
         IJK_TEMP = DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT)%IJK_SCAL 
         WALL_TYPE_TEMP = DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT)%DES_BC_TYPE
         NORMAL_TEMP = ZERO 
         IF(TRIM(WALL_TYPE_TEMP).eq.'CUT_FACE') NORMAL_TEMP(1:3) = NORMAL_S(IJK_TEMP,1:3)
         IF(TRIM(WALL_TYPE_TEMP).eq.'NORMAL_WALL') NORMAL_TEMP(1:DIMN) = DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT)%NORMAL(1:DIMN) 

         IF(myPE.eq.pe_IO) WRITE(*,1037) CELL_ID, I_OF(CELL_ID), J_OF(CELL_ID), K_OF(CELL_ID), CUT_CELL_AT(CELL_ID), COUNT, TRIM(WALL_TYPE_TEMP), IJK_TEMP, I_OF(IJK_TEMP), J_OF(IJK_TEMP), K_OF(IJK_TEMP), NORMAL_TEMP(1:3)
      ENDDO

1037   FORMAT(  & 
            & 'CELL ID, CUT_CELL ?: ', 4(2x,i8), 2x, L1, /10x, & 
            & 'COUNT_BC: ', i3, /10x, & 
            & 'BC_TYPE : ', A, /10x, &
            & 'BC_IJK  : ', i10, /10X, &
            & 'BC_I    : ', i4, /10X, &
            & 'BC_J    : ', i4, /10X, &
            & 'BC_K    : ', i4, /10X, &
            & 'NORMAL  : ', 3(2x,g17.8))
      
      CALL mfix_exit(mype)
      RETURN
      END SUBROUTINE DES_WALLBC_PREPROCSSING

      SUBROUTINE ADD_SMALL_CELL_BC(CELL_ID)
       USE param1
       USE funits
       USE run
       USE compar      
       USE discretelement
       USE cutcell
       USE constant 
       USE indices
       USE physprop
       USE parallel
       USE geometry
       
       Implicit none 
       
       INTEGER, INTENT(IN) :: CELL_ID

       double precision :: NORM_CF(3) !normal of the cut_face 
       double precision :: NORM_FACE(6, 3) 
       integer :: count_face_loop, count
       double precision :: TOL_NORM_ANG, MAX_ANG_PLANE, NORM_SQ
       
       CHARACTER*100 :: WALL_TYPE

       NORM_CF (1:DIMN) = NORMAL_S(CELL_ID,1:DIMN)

       MAX_ANG_PLANE = 90.d0 
       TOL_NORM_ANG =  ASIN(1.d0)*MAX_ANG_PLANE/90.d0

       TOL_NORM_ANG  =  COS(TOL_NORM_ANG)

       IF(DIMN.EQ.2) THEN 
          COUNT_FACE_LOOP = 4
       ELSE
          COUNT_FACE_LOOP = 6
       ENDIF
       
       !store the normal information for the faces 
       !think of the whole cell as a obstruction and assume fluid all
       !around. 
       !And, keep the convention for the face normal same 
       !as cut_face normal, i.e., pointing toward the fluid. 
       !later on change the sign in and store this bc as normal_wall 
       !that will make it consistent with normal definition in normal_wall:
       !pointing away from the fluid 

       NORM_FACE  = ZERO
       !1 = East
       !2 = West
       !3 = North
       !4 = South
       !5 = Top
       !6 = Bottom

       
       NORM_FACE(1, 1) =  one  !from east pointing toward fluid will be +1 
       NORM_FACE(2, 1) = -one

       NORM_FACE(3, 2) =  one
       NORM_FACE(4, 2) = -one

       NORM_FACE(5, 3) =  one
       NORM_FACE(6, 3) = -one
       
       !WRITE(*,2001)  
       !WRITE(*,'(A,3(2x,i5))') 'SMALL CELL:', CELL_ID, I_OF(CELL_id), j_of(cell_id), k_of(cell_id)
       !WRITE(*,'(A,3(2x,g17.8))') 'NORMAL OF CUT FACE  =', NORM_CF( 1:dimn)

       DO count  = 1, COUNT_FACE_LOOP 
          NORM_SQ = DOT_PRODUCT(NORM_FACE(count, 1:dimn), NORM_CF(1:dimn))
          WALL_TYPE = 'SMALL_CELL'
          IF(NORM_SQ.GT.TOL_NORM_ANG) THEN 
             
             CALL DES_ADD_WALLBC_TO_CELL(CELL_ID, CELL_ID, NORM_FACE(count, 1:dimn), WALL_TYPE)
             !WRITE(*,'(A,3(2x,g17.8))') 'NORMAL OF WALL NORMAL  =', -norm_face(count, 1:dimn)
             
             !WRITE(*,*) 'count_bc = ',        DES_CELLWISE_BCDATA(CELL_ID)%COUNT_DES_BC
          ENDIF
       end DO
       
       
       !WRITE(*,2001)  
 2001 FORMAT(/1X,70('*'))
      end SUBROUTINE ADD_SMALL_CELL_BC
      
      SUBROUTINE CHECK_IF_TOADD_THIS_CUT_FACE(CELL_ID, WALL_CELL_ID, WALL_NORM, DES_WALL_TYPE, ADD_THIS_CUT_FACE)
       USE param1
       USE funits
       USE run
       USE compar      
       USE discretelement
       USE cutcell
       USE constant 
       USE indices
       USE physprop
       USE parallel
       USE geometry
       USE mfix_pic
       IMPLICIT NONE 
       !Cell_id is the IJK of the cell to which the new bc is being considered for addition
       !wall_cell_id is the IJK of the cell where the new bc is present 
       INTEGER, INTENT(IN) :: CELL_ID
       INTEGER, INTENT(INOUT) :: WALL_CELL_ID
       DOUBLE PRECISION, DIMENSION(DIMN), INTENT(INOUT) :: WALL_NORM
       CHARACTER*100, INTENT(INOUT) :: DES_WALL_TYPE 
       LOGICAL, INTENT(OUT) :: ADD_THIS_CUT_FACE
       !Local variables


       INTEGER :: COUNT_BC, IJK_TEMP, IJK_TEMP2, COUNT, I, J, K
       CHARACTER*100 :: WALL_TYPE_TEMP, WALL_TYPE_TEMP2, DES_WALL_TYPE_FINAL
       LOGICAL :: CHECK_NORMAL_AND_DIST, INT_OUTSIDECELL, REJECT_NEW, ATLEAST_ONEINSIDE, FORCE_INTPOINT_INSIDE, INT_LINE_COR_OUTSIDECELL
       DOUBLE PRECISION :: DIAGONAL, TOL_DELH_DES, WALL_NORM_FINAL(DIMN)
       DOUBLE PRECISION :: NORM1_NEW, NORM2_NEW, NORM3_NEW
       DOUBLE PRECISION :: NORM1_OLD, NORM2_OLD, NORM3_OLD
       DOUBLE PRECISION :: NORM_SQ, DIST_OLD, DIST_NEW, DIST_NEW_INT, NORMAL_OLD(DIMN)
       DOUBLE precision :: xcor(8), ycor(8), zcor(8)!store the vertices of the scalar cell with id cell_id
       DOUBLE PRECISION :: XINT_OLD, YINT_OLD, ZINT_OLD, XMINCOR, XMAXCOR, YMINCOR, YMAXCOR, ZMINCOR, ZMAXCOR 
       INTEGER ::  COUNT_VERT_LOOP, COUNT2, COUNT3, COUNT_FACE_LOOP, COUNTBCIN
       INTEGER :: DIST_CORDS_OLD, DIST_CORDS_NEW, WALL_CELL_ID_FINAL
       
       DOUBLE PRECISION :: coeff_a1, coeff_b1, coeff_c1, coeff_a2, coeff_b2, coeff_c2, rootx, rooty, det, NORM_AVG(DIMN), NORM_MAG, TOL_CONVEX_ANGLE_DES, TOL_CONVEX_ANGLE, TMP_DOT, XINT_COR_LINE(8, 3)

       DOUBLE PRECISION ::  CNOT_ARR(3), VEC_ARR(3), DISTFROMLINE, PTXYZ_ON_LINE(3)
       DOUBLE PRECISION :: TOL_LINEPLANE_DIST_DES, TOL_LINEINT, TOL_NORM_ANG, TOL_NORM_ANG_DES

       DOUBLE PRECISION :: XREF_FACE(6), YREF_FACE(6), ZREF_FACE(6), NORMAL_FACE(6,3), NORM_ANG, MIDWAY_PT(3)

       DOUBLE PRECISION :: TOL_DELH_DIAG, TOL_DELH_DIAG_DIST
       
       !distance in grid units
       DOUBLE PRECISION :: DIST_GU_OLD, DIST_GU_NEW

       LOGICAL :: DEBUG_LOCAL 
       
       DEBUG_LOCAL = .false. 

       ADD_THIS_CUT_FACE = .FALSE. 
       TOL_LINEPLANE_DIST_DES = 0.01
       
       TOL_CONVEX_ANGLE_DES  = 5.d0
       TOL_CONVEX_ANGLE = ASIN(1.d0)*TOL_CONVEX_ANGLE_DES/90.d0
       TOL_CONVEX_ANGLE  =  COS(TOL_CONVEX_ANGLE)
       
       TOL_NORM_ANG_DES = 10.d0
       TOL_NORM_ANG =  ASIN(1.d0)*TOL_NORM_ANG_DES/90.d0
       TOL_NORM_ANG  =  COS(TOL_NORM_ANG)
       

       TOL_DELH_DES = 0.05
       
       TOL_DELH_DIAG = 0.01
       
       WALL_CELL_ID_FINAL = WALL_CELL_ID
       DES_WALL_TYPE_FINAL = DES_WALL_TYPE
       WALL_NORM_FINAL(:) = WALL_NORM(:)

       

       IF(DIMN.EQ.2) THEN 
          COUNT_VERT_LOOP = 4   
          COUNT_FACE_LOOP = 4
       ELSE
          COUNT_VERT_LOOP = 8
          COUNT_FACE_LOOP = 6
       ENDIF

       CHECK_NORMAL_AND_DIST = .FALSE.
       COUNT_BC = DES_CELLWISE_BCDATA(CELL_ID)%COUNT_DES_BC

       !First check if the same wall_cell_id with the identical wall_type does not already exist 

       DO COUNT = 1, COUNT_BC 
          IJK_TEMP = DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT)%IJK_SCAL
          WALL_TYPE_TEMP = DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT)%DES_BC_TYPE
          IF((IJK_TEMP.EQ.WALL_CELL_ID).AND.(WALL_TYPE_TEMP.EQ.DES_WALL_TYPE))  RETURN
       ENDDO

       ! if the new wall bc in question is cut-face then test it against the existing cut-face BC's for this cell
       
       !store the normal information of the new cut-face                                
       NORM1_NEW = NORMAL_S(WALL_CELL_ID, 1)
       NORM2_NEW = NORMAL_S(WALL_CELL_ID, 2)
       NORM3_NEW = NORMAL_S(WALL_CELL_ID, 3)
       I = I_OF(CELL_ID)
       J = J_OF(CELL_ID)
       K = K_OF(CELL_ID)

       Diagonal = (DX(I)**2 + DY(J)**2)
       IF(DIMN.eq.3)  Diagonal = diagonal + DZ(K)**2
       Diagonal = dsqrt(diagonal)
       TOL_LINEINT = - DIAGONAL*TOL_LINEPLANE_DIST_DES  
       TOL_DELH_DIAG_DIST = TOL_DELH_DIAG*Diagonal
       !store the coordinate information of the scalar cell for the cell to which new bc is being added 
       !this will be used to find if the new cut-face is between any of the eight vertices and an existing
       !cut-face.              
       ZCOR(1:4) = ZG_T(K)
       
       XCOR(1) = XG_E(I)
       YCOR(1) = YG_N(J)
       
       XCOR(2) = XG_E(I) - DX(I)
       YCOR(2) = YG_N(J)
       
       XCOR(3) = XG_E(I)
       YCOR(3) = YG_N(J) - DY(J)
       
       XCOR(4) = XG_E(I) - DX(I)
       YCOR(4) = YG_N(J) - DY(J)
       
       ZCOR(5:8) = ZG_T(K) - DZ(K)
       
       XCOR(5) = XG_E(I)
       YCOR(5) = YG_N(J)
       
       XCOR(6) = XG_E(I) - DX(I)
       YCOR(6) = YG_N(J)
       
       XCOR(7) = XG_E(I)
       YCOR(7) = YG_N(J) - DY(J)
       
       XCOR(8) = XG_E(I) - DX(I)
       YCOR(8) = YG_N(J) - DY(J)
       
       XMINCOR = XG_E(I) - DX(I)
       XMAXCOR = XG_E(I) 
       YMINCOR = YG_N(J) - DY(J)
       YMAXCOR = YG_N(J)
       
       ZMINCOR = ZG_T(K) - DZ(K)
       ZMAXCOR = ZG_T(K)

       IF(.NOT.CUT_CELL_AT(CELL_ID)) THEN 
          !first check if the new cut-face bc is not cutting the non-cut cell if extended
          !this test is ok for non cut cells. For cutcells, the relative
          !distance of each vertex from cutfaces has to be compared 
          VERTLOOP1: DO COUNT2 = 1, COUNT_VERT_LOOP
             CALL GET_DEL_H_DES(WALL_CELL_ID,'SCALAR', XCOR(COUNT2), YCOR(COUNT2), ZCOR(COUNT2), & 
                  & DIST_NEW, NORM1_NEW, NORM2_NEW, NORM3_NEW, .true.)

             !If the new cut-face has a negative distance from any one of the vertices of a normal cell, then reject this wall bc 

             !!IF(DIST_NEW.LT.ZERO) RETURN
             !Rahul: change on Jan 23, 2012.
             !the non-cut cells will only exclude cut-faces that have normal's 
             !within tolerance 
             
          ENDDO VERTLOOP1
       ENDIF

       !check the normal of this cut-face against the normal of 
       !existing cut-faces. If the normals are nearly aligned then
       !reject this new cut-face 
       DO COUNT = 1, COUNT_BC 
          IJK_TEMP = DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT)%IJK_SCAL
          WALL_TYPE_TEMP = DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT)%DES_BC_TYPE
          CUTFACE: IF(TRIM(WALL_TYPE_TEMP).EQ.'CUT_FACE') THEN 
             NORM1_OLD = NORMAL_S(IJK_TEMP, 1)
             NORM2_OLD = NORMAL_S(IJK_TEMP, 2)
             NORM3_OLD = NORMAL_S(IJK_TEMP, 3)

             NORM_SQ = NORM1_OLD*NORM1_NEW + NORM2_OLD*NORM2_NEW + NORM3_OLD*NORM3_NEW

             ! now check if there already exists a cut-face wall with the same normal as the 
             ! one being added 
             
             !NORM_ANG = 180.d0*(ACOS(NORM_SQ))/PI
             IF(NORM_SQ.GT.TOL_NORM_ANG.and.CUT_CELL_AT(CELL_ID).and..true.) THEN 
!!$                WRITE(*,*) '***************************************************************'
!!$                WRITE(*,*) 'THE CURRENT CUT-CELL WITH THE SAME NORMAL ALREADY EXISTS'
!!$                ! WRITE(*,*) 'WALL_TYPE OLD AND NEW= ', TRIM(WALL_TYPE_TEMP), TRIM(DES_WALL_TYPE)
!!$                WRITE(*,*) 'BC # and TOTAL BCs', COUNT, COUNT_BC
!!$                WRITE(*,*) 'WALL_CELL ID OLD AND NEW =', IJK_TEMP, WALL_CELL_ID
!!$                
!!$                WRITE(*,'(A,3(2x,g17.8))') 'NORMAL OLD =', NORM1_OLD, NORM2_OLD, NORM3_OLD
!!$                WRITE(*,'(A,3(2x,g17.8))') 'NORMAL NEW =', NORM1_NEW, NORM2_NEW, NORM3_NEW
!!$
!!$                WRITE(*,*) 'NORM_ANG = ', NORM_ANG, NORM_SQ
!!$                WRITE(*,'(A,3(2x,i4))') 'I, J, K NAT = ', I, J, K
!!$                WRITE(*,'(A,3(2x,i4))') 'I, J, K OLD = ', I_OF(IJK_TEMP), J_OF(IJK_TEMP), K_OF(IJK_TEMP)
!!$                WRITE(*,'(A,3(2x,i4))') 'I, J, K NEW = ', I_OF(WALL_CELL_ID), J_OF(WALL_CELL_ID), K_OF(WALL_CELL_ID)                      
!!$                WRITE(*,*) '*********************************************************************'
                IF(CUT_CELL_AT(CELL_ID)) THEN 
                   RETURN
                ELSE
                   !compare the distance in grid units. 
                   !keep the closer one
                   DIST_GU_OLD = ABS(I_OF(CELL_ID)-I_OF(IJK_TEMP)) + ABS(J_OF(CELL_ID)-J_OF(IJK_TEMP))
                   
                   IF(DIMN.Eq.3) DIST_GU_OLD = DIST_GU_OLD + ABS(K_OF(CELL_ID)-K_OF(IJK_TEMP))
                   DIST_GU_NEW = ABS(I_OF(CELL_ID)-I_OF(WALL_CELL_ID)) + ABS(J_OF(CELL_ID)-J_OF(WALL_CELL_ID))
                   
                   IF(DIMN.Eq.3) DIST_GU_NEW = DIST_GU_NEW + ABS(K_OF(CELL_ID)-K_OF(WALL_CELL_ID))
                   
                   IF(DIST_GU_NEW.LT.DIST_GU_OLD) then
                      
                      !replace ijk_scal information for the existing
                      !BC with the newer bc ijk
                      DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT)%IJK_SCAL = WALL_CELL_ID
                      DEBUG_LOCAL = .false.
                      IF(DEBUG_LOCAL.and.dmp_log) THEN
                         WRITE(*,*) 'REPLACING NEARLY SAME NORMALS WITH CLOSER CUT-FACE'
                         WRITE(*,*) 'CUT_FACE CELL ? = ', CUT_CELL_AT(CELL_ID)
                         WRITE(*,'(A,3(2x,g17.8))') 'NORMAL OLD =', NORM1_OLD, NORM2_OLD, NORM3_OLD
                         WRITE(*,'(A,3(2x,g17.8))') 'NORMAL NEW =', NORM1_NEW, NORM2_NEW, NORM3_NEW

                         WRITE(*,'(A,3(2x,i4))') 'I, J, K NAT = ', I, J, K
                         WRITE(*,'(A,3(2x,i4))') 'I, J, K OLD = ', I_OF(IJK_TEMP), J_OF(IJK_TEMP), K_OF(IJK_TEMP)
                         WRITE(*,'(A,3(2x,i4))') 'I, J, K NEW = ', I_OF(WALL_CELL_ID), J_OF(WALL_CELL_ID), K_OF(WALL_CELL_ID)
                         READ(*,*)
                      ENDIF
                      DEBUG_LOCAL = .false.
                   ELSE
                      !nothing more to do. Exit this routine 
                      RETURN 
                   ENDIF
                ENDIF
             end IF
          END IF CUTFACE
       end DO


       CUT_CELL_SPECIAL: IF(COUNT_BC.GE.1.AND.CUT_CELL_AT(CELL_ID).and..true.) THEN 
          !Do the remanining checks only if the current cell is a cut-cell
          COUNT = 1 !assuming that the first bc for a cut-cell cell will be the cut-face
          !that is native to this cell. So far in the development, this will be case.
          !adding another sanity check that will stop the code if otherwise
          
          IJK_TEMP = DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT)%IJK_SCAL
          WALL_TYPE_TEMP = DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT)%DES_BC_TYPE
          
          IF(CELL_ID.NE.IJK_TEMP.AND.TRIM(WALL_TYPE_TEMP).NE.'CUT_FACE')THEN
             IF(DMP_LOG) WRITE(UNIT_LOG, 2001) COUNT_BC, I, J, K, & 
             & I_OF(IJK_TEMP), J_OF(IJK_TEMP), K_OF(IJK_TEMP)

             IF(PRINT_DES_SCREEN) WRITE(*, 2001) COUNT_BC, I, J, K, & 
             & I_OF(IJK_TEMP), J_OF(IJK_TEMP), K_OF(IJK_TEMP)
             
             CALL mfix_exit(myPE)
          ENDIF

          
          !if the code has not returned so far, compare the perpendicular distance of all vertices of this cell from the new and old cut-faces. If the distance from the new cut-face is less than the distance from the old cut-face from any of the vertices, then do not add this new cut-face.
          
          ATLEAST_ONEINSIDE = .false.
          FORCE_INTPOINT_INSIDE = .TRUE. 
          
          ! 300      CONTINUE 
          
          VERTLOOP:  DO COUNT2 = 1, COUNT_VERT_LOOP
             !First find the normal distance from the vertex to the existing cut-face 
             CALL GET_DEL_H_DES(IJK_TEMP,'SCALAR',XCOR(COUNT2), YCOR(COUNT2), ZCOR(COUNT2), & 
                  & DIST_OLD, NORM1_OLD, NORM2_OLD, NORM3_OLD, .true.)
             
             !IF(DIST_OLD.LT.TOL_DELH_DES*DIAGONAL) CYCLE !if the vertex in question is not in the fluid domain, then cycle to the next vertex 
             IF(DIST_OLD.LE.ZERO) CYCLE 
             !here => dist_old is postitive. Now find the point of intersection for the normal from this vertex to the existing cut-face 
             XINT_OLD  = XCOR(COUNT2) - DIST_OLD*NORM1_OLD
             YINT_OLD  = YCOR(COUNT2) - DIST_OLD*NORM2_OLD
             ZINT_OLD  = ZCOR(COUNT2) - DIST_OLD*NORM3_OLD
             
             !at sharp edges the intersection point might lie outside of the cell containing the cut-face... and this point as a result mmight lie on the wrong side of the cut-face being considered
             !if the new cut-face is indeed before the existing cut-face then it will show up this way from other vertices where the intersection point will also lie within the cell 
             INT_OUTSIDECELL = .false.
             
             IF(XINT_OLD.LT.XMINCOR.OR.XINT_OLD.GT.XMAXCOR.OR.YINT_OLD.LT.YMINCOR.OR.YINT_OLD.GT.YMAXCOR) INT_OUTSIDECELL =.true.
             IF(DIMN.EQ.3.AND.(ZINT_OLD.LT.ZMINCOR.OR.ZINT_OLD.GT.ZMAXCOR)) INT_OUTSIDECELL =.true.
             !now find the distance from the vertex and the intersection point to the new cut-face in question 
             
             CALL GET_DEL_H_DES(WALL_CELL_ID,'SCALAR', XCOR(COUNT2), YCOR(COUNT2), ZCOR(COUNT2), & 
                  & DIST_NEW, NORM1_NEW, NORM2_NEW, NORM3_NEW, .true.)
             
             CALL GET_DEL_H_DES(WALL_CELL_ID,'SCALAR', XINT_OLD, YINT_OLD, ZINT_OLD, & 
                  & DIST_NEW_INT, NORM1_NEW, NORM2_NEW, NORM3_NEW, .true.)
             
             NORM_SQ = NORM1_OLD*NORM1_NEW + NORM2_OLD*NORM2_NEW + NORM3_OLD*NORM3_NEW
             IF(DIST_OLD.EQ.DIST_NEW.AND.DIST_OLD.EQ.DIST_NEW_INT) CYCLE
             REJECT_NEW = .FALSE.
             IF(FORCE_INTPOINT_INSIDE) THEN 
                IF(INT_OUTSIDECELL ) THEN 
                   IF(DIST_NEW.LT.ZERO) REJECT_NEW = .true.
                ELSE
                   ATLEAST_ONEINSIDE = .TRUE.
                   IF(DIST_NEW_INT.LT.ZERO) REJECT_NEW = .true.
                ENDIF
             ELSE
                IF(DIST_NEW.LT.ZERO.OR.DIST_NEW_INT.LT.ZERO)  REJECT_NEW = .true.
             ENDIF

             IF(FORCE_INTPOINT_INSIDE) THEN 
                REJECT_NEW = REJECT_NEW.AND.ATLEAST_ONEINSIDE 
                !reject will be true only if the intersection point is found 
                !to be inside the cell 
             ELSE
                !If this does not happen then this loop will be repeated with 
                !FORCE_INTPOINT_INSIDE = .false. 
                !Jan 23, 2012. 
                !not anymore, this loop will not be retested by setting 
                !force_point_inside = .false.
                !rather the distances will be compared from a point 
                !midway on line joining the reference points on the two
                !cut-faces
             ENDIF
             
             
             IF(REJECT_NEW) THEN
                DEBUG_LOCAL = .false.
                IF(DEBUG_LOCAL.and.dmp_log) THEN
                   WRITE(*,*) 
                   WRITE(*,*) '*********************************************************************'
                   WRITE(*,*) 'THE NEW CUT-FACE IS CLOSER TO ONE OF THE VERTICES OR THE INTERSECTION POINT'
                   WRITE(*,*) 'FORCING THE INTERSECTION POINT TO INSIDE CELL ?', FORCE_INTPOINT_INSIDE
                   WRITE(*,*) 'NOT ADDING THE NEW CUT-FACE '
                   WRITE(*,*) 'WALL_TYPE OLD AND NEW= ', TRIM(WALL_TYPE_TEMP),  TRIM(DES_WALL_TYPE), -TOL_DELH_DES*DIAGONAL
                   
                   WRITE(*,'(A,3(2x,i4))') 'CELL ID, CELL ID OLD, WALL CELL ID =', CELL_ID, IJK_TEMP, WALL_CELL_ID
                   WRITE(*,'(A,3(2x,i4))') 'I, J, K NAT = ', I, J, K
                   WRITE(*,'(A,3(2x,i4))') 'I, J, K OLD = ', I_OF(IJK_TEMP), J_OF(IJK_TEMP), K_OF(IJK_TEMP)
                   WRITE(*,'(A,3(2x,i4))') 'I, J, K NEW = ', I_OF(WALL_CELL_ID), J_OF(WALL_CELL_ID), K_OF(WALL_CELL_ID)
                   WRITE(*,'(A,3(2x,g17.8))') 'NORMAL OLD =', NORM1_OLD, NORM2_OLD, NORM3_OLD
                   WRITE(*,'(A,3(2x,g17.8))') 'NORMAL NEW =', NORM1_NEW, NORM2_NEW, NORM3_NEW
                   WRITE(*,'(A,3(2x,g17.8))') 'DIST OLD AND NEW=', DIST_OLD, DIST_NEW, DIST_NEW_INT
                   WRITE(*,*) 'INTERSECTION POINT WITHIN THE CELL ?', .NOT.INT_OUTSIDECELL
                   WRITE(*,'(A,3(2x,g17.8))') 'INTERSECTION POINT WITH THE OLD =', XINT_OLD, YINT_OLD, ZINT_OLD
                   
                   WRITE(*,'(A,3(2x,g17.8))') 'EAST NORTH AND TOP OF THE CELL =', XG_E(I), YG_N(J), ZG_T(K)
                   WRITE(*,'(A,3(2x,g17.8))') 'WEST SOUTH AND BOTTOM OF THE CELL =', XG_E(I)-DX(I), YG_N(J)-DY(J), ZG_T(K) - DZ(K)
                   WRITE(*,*) 'VERTEX NUMBER = ', COUNT2, MPPIC
                   !IF(I.EQ.2) read(*,*) !.and.J.eq.25) read(*,*)
                   !IF(.NOT.FORCE_INTPOINT_INSIDE) READ(*,*)
                   
                   DEBUG_LOCAL = .false.
                ENDIF
                IF(DIST_OLD.EQ.UNDEFINED.OR.DIST_NEW.EQ.UNDEFINED) THEN 
                   
                   !IF(DEBUG_LOCAL.and.dmp_log) THEN
                   IF(DMP_LOG) THEN 
                      WRITE(UNIT_LOG,*) 'ERROR MESSAGE FROM DES_ADD_WALLBC_TO_CELL'
                      WRITE(UNIT_LOG,*) 'IN COMPARING DISTANCE FROM VERTICES, ONE OF THE DISTANCE IS UNDEFINED'
                      WRITE(UNIT_LOG,'(A,3(2x,g17.8))') 'NORMAL OLD =', NORM1_OLD, NORM2_OLD, NORM3_OLD
                      WRITE(UNIT_LOG,'(A,3(2x,g17.8))') 'NORMAL NEW =', NORM1_NEW, NORM2_NEW, NORM3_NEW
                      WRITE(UNIT_LOG,'(A,3(2x,g17.8))') 'DIST OLD AND NEW=', DIST_OLD, DIST_NEW
                   
                      WRITE(UNIT_LOG,*) 'TERMINAL ERROR: STOPPING'
                   ENDIF
                   CALL mfix_exit(mype)
                ENDIF
                
                !IF(MPPIC) RETURN
                RETURN
                !what follows below is more relevant for soft spring like
                !collions where convex edge has to be replaced by a line
                !which is the intersection of planes 

                write(*,*) 'now will try to add cut-face line for this case'
                NORM_SQ = NORM1_OLD*NORM1_NEW + NORM2_OLD*NORM2_NEW + NORM3_OLD*NORM3_NEW
                
                IF(NORM_SQ.GT.TOL_CONVEX_ANGLE) THEN
                   write(*,*) 'not adding as the convex angle threshold not reached' 
                                !do not add cut-face-line as a new bc 
                   RETURN 
                ELSE
                   IJK_TEMP2  = WALL_CELL_ID
                   CALL DES_CROSSPRDCT_3D(VEC_ARR(:), NORMAL_S(IJK_TEMP,:), NORMAL_S(IJK_TEMP2,:))
                   VEC_ARR(:) = VEC_ARR(:)/SQRT(DOT_PRODUCT(VEC_ARR(:), VEC_ARR(:)))
                   coeff_a1 = NORMAL_S(IJK_TEMP,1) 
                   coeff_b1 = NORMAL_S(IJK_TEMP,2) 
                   coeff_a2 = NORMAL_S(IJK_TEMP2,1) 
                   coeff_b2 = NORMAL_S(IJK_TEMP2,2) 
                   
                   coeff_c1 = DOT_PRODUCT(REFP_S(IJK_TEMP,:), NORMAL_S(IJK_TEMP,:))
                   coeff_c2 = DOT_PRODUCT(REFP_S(IJK_TEMP2,:), NORMAL_S(IJK_TEMP2,:))
                   
                   det = coeff_a1*coeff_b2 - coeff_a2*coeff_b1
                   IF(ABS(DET).LT.0.00001) THEN
                      WRITE(*,*) 'ERROR MESSAGE FROM CHECK_IF_TOADD_THIS_CUT_FACE'
                      WRITE(*,*) 'PLANES DO NOT INTERSECT AS DETERMINANT NEARLY EQUAL TO ZERO', DET
                      WRITE(*,'(A35, 3(2x,g17.8))') 'COeff for first plane  = ', coeff_a1,coeff_b1, coeff_c1
                      WRITE(*,'(A35, 3(2x,g17.8))') 'COeff for second plane  = ', coeff_a2,coeff_b2, coeff_c2
                      
                      WRITE(*,'(A35, 3(2x,i4))') 'I, J, K of first plane', I_OF(IJK_TEMP), J_OF(IJK_TEMP), K_OF(IJK_TEMP)
                      WRITE(*,'(A35, 3(2x,i4))') 'I, J, K of second plane', I_OF(IJK_TEMP2), J_OF(IJK_TEMP2), K_OF(IJK_TEMP2)
                      WRITE(*,'(A35, 3(2x,g17.8))') 'NORMAL 1 =', NORMAL_S(IJK_TEMP,: )
                      WRITE(*,'(A35, 3(2x,g17.8))') 'NORMAL 2 =', NORMAL_S(IJK_TEMP2,: )
                      
                      WRITE(*,*) 'TERMINAL ERROR: STOPPING'
                      call mfix_exit(myPE)
                   ENDIF
                   
                   ROOTX = (COEFF_C1*COEFF_B2 - COEFF_C2*COEFF_B1)/DET
                   ROOTY = (COEFF_A1*COEFF_C2 - COEFF_A2*COEFF_C1)/DET
                   CNOT_ARR(1) = ROOTX
                   CNOT_ARR(2) = ROOTY
                   CNOT_ARR(3) = ZERO
                   
                   
                   NORM_AVG(:) = HALF*(NORMAL_S(IJK_TEMP,:)+NORMAL_S(IJK_TEMP2,:))
                   
                   NORM_MAG = SQRT(DOT_PRODUCT(NORM_AVG, NORM_AVG))
                   NORM_AVG(:) = NORM_AVG(:)/NORM_MAG
                   
                   NORM_MAG = SQRT(DOT_PRODUCT(NORM_AVG, NORM_AVG))
                   !change the wall information
                   WALL_NORM_FINAL(:) = NORM_AVG(:)
                   DES_WALL_TYPE_FINAL = 'CUT_FACE_LINE'
                   WALL_CELL_ID_FINAL = IJK_TEMP
                   IF(NORM_MAG.LT.0.99) THEN
                      WRITE(*,*) 'ERROR MESSAGE FROM CHECK_IF_TOADD_THIS_CUT_FACE'
                      WRITE(*,*) 'NEW NORMAL MAGNITUDE NE 1', NORM_MAG
                      WRITE(*,*) 'NOMRAL_AVG = ', NORM_AVG(:)
                      WRITE(*,*) 'TERMINAL ERROR: STOPPING'
                      STOP
                   ENDIF
                   
                   WRITE(*,*) 'CHANGED THE BC TO CUT_FACE_LINE IN THE PARAMETERIC FORM vec(CNOT) + t*vec(v)'
                   WRITE(*,'(A35, 3(2x,g17.8))') 'LINE DIRECTION COSINE OR VEC(V)= ', VEC_ARR(:)
                   WRITE(*,'(A35, 3(2x,g17.8))') 'LINE INTERCEPT OR VEC(CNOT) = ', CNOT_ARR(:)
                   
                   WRITE(*,'(A35, 3(2x,g17.8))') 'OLD NORMAL 1 =', NORMAL_S(IJK_TEMP,: )
                   WRITE(*,'(A35, 3(2x,g17.8))') 'OLD NORMAL 2 =', NORMAL_S(IJK_TEMP2,: )
                   WRITE(*,'(A35, 3(2x,g17.8))') 'NEW NORMAL  = ', WALL_NORM_FINAL(:)
                   
                   
                   WRITE(*,'(A35, 3(2x,g17.8))') 'REF PT 1 =', REFP_S(IJK_TEMP,:), IJK_TEMP          
                   WRITE(*,'(A35, 3(2x,g17.8))') 'REF PT 2 =', REFP_S(IJK_TEMP2,:), IJK_TEMP2          
                   WRITE(*,'(A35, 3(2x,g17.8))') 'COeff for first plane  = ', coeff_a1,coeff_b1, coeff_c1
                   WRITE(*,'(A35, 3(2x,g17.8))') 'COeff for second plane  = ', coeff_a2,coeff_b2, coeff_c2
                   
                   
                   WRITE(*,'(80(A))') '***********************************************************************'
                   !goto 500
                   NORMAL_FACE(6,:) = 0.d0
                   !East face 
                   XREF_FACE(1) = XG_E(I)
                   YREF_FACE(1) = YG_N(J)
                   ZREF_FACE(1) = ZG_T(K)
                   NORMAL_FACE(1,1) = -1.d0
                   
                   !west face 
                   XREF_FACE(2) = XG_E(I) - DX(I)
                   YREF_FACE(2) = YG_N(J)
                   ZREF_FACE(2) = ZG_T(K)
                   NORMAL_FACE(2,1) = 1.d0
                   
                   !North face 
                   XREF_FACE(3) = XG_E(I)
                   YREF_FACE(3) = YG_N(J)
                   ZREF_FACE(3) = ZG_T(K)
                   NORMAL_FACE(3,2) = -1.d0
                   
                   
                   !South face 
                   XREF_FACE(4) = XG_E(I)
                   YREF_FACE(4) = YG_N(J) - DY(J)
                   ZREF_FACE(4) = ZG_T(K)
                   NORMAL_FACE(4,2) = 1.d0
                   
                   !Top face 
                   XREF_FACE(5) = XG_E(I)
                   YREF_FACE(5) = YG_N(J) 
                   ZREF_FACE(5) = ZG_T(K)
                   NORMAL_FACE(5,3) = -1.d0
                   
                   !Bottom face 
                   XREF_FACE(6) = XG_E(I)
                   YREF_FACE(6) = YG_N(J) 
                   ZREF_FACE(6) = ZG_T(K) - DZ(K)
                   NORMAL_FACE(6,3) = 1.d0
                      

                   FACELOOP:  DO COUNT3 = 1, COUNT_FACE_LOOP
                      !first check if the line and plane do not intersect 
                      TMP_DOT = DOT_PRODUCT(VEC_ARR(1:3), NORMAL_FACE(COUNT3,1:3))
                      IF(ABS(TMP_DOT).GT.0.00001)  CYCLE !if the dot product is finite, that implies the line and the 
                      !plane intersect and computing normal distance is useless 
                      PTXYZ_ON_LINE(1) = CNOT_ARR(1) + 1.d0*VEC_ARR(1)
                      PTXYZ_ON_LINE(2) = CNOT_ARR(2) + 1.d0*VEC_ARR(2)
                      PTXYZ_ON_LINE(3) = CNOT_ARR(3) + 1.d0*VEC_ARR(3)
                      DISTFROMLINE = NORMAL_FACE(COUNT3,1) * (PTXYZ_ON_LINE(1) - XREF_FACE(COUNT3)) + & 
                           &  NORMAL_FACE(COUNT3,2) * (PTXYZ_ON_LINE(2) - YREF_FACE(COUNT3)) + &
                           &  NORMAL_FACE(COUNT3,3) * (PTXYZ_ON_LINE(3) - ZREF_FACE(COUNT3)) 
                      
                      !IF(DISTFROMLINE.LT.TOL_LINEINT) THEN 
                         
                      IF(DISTFROMLINE.LT.ZERO) THEN 
                         WRITE(*,*) 'INTERSECTION POINT OUTSIDE THE CELL '
                         WRITE(*,*) 'FACE NUMBER  = ', COUNT3, TOL_LINEINT
                         
                         WRITE(*,'(A,3(2x,g17.8))') 'FACE COORDINATES =', XREF_FACE(COUNT3), yref_face(count3), zref_face(count3) 
                         
                         WRITE(*,'(A,3(2x,g17.8))') 'POINT ON LINE = ', PTXYZ_ON_LINE(1:3)
                         WRITE(*,'(A,6(2x,g17.8))') 'DISTANCE =',  DISTFROMLINE
                         
                         !READ(*,*)
                         !don't add this line in this case 
                         RETURN
                      end IF
                   END DO FACELOOP
                   
                   goto 500 
                   !also check if the same bc type with the same attributes does not already exist
                   DO COUNTBCIN = 1, COUNT_BC 
                      IJK_TEMP2 = DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNTBCIN)%IJK_SCAL
                      WALL_TYPE_TEMP2 = DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNTBCIN)%DES_BC_TYPE
                      CUTFACELINE: IF(TRIM(WALL_TYPE_TEMP2).EQ.'CUT_FACE_LINE') THEN 
                         NORMAL_OLD(:) = DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNTBCIN)%NORMAL(:)
                         NORM_SQ = DOT_PRODUCT(WALL_NORM_FINAL(:), NORMAL_OLD(:))
                         
                         IF(NORM_SQ.GT.0.99d0) THEN 
                            WRITE(*,*) 'WAS GOING TO ADD A NEW CUT_FACE_LINE, BUT FOUND A CUT_FACE_LINE BC WITH THE SAME NORMAL'
                            WRITE(*,*) 'BC ID and TOTAL BCs =', COUNTBCIN, COUNT_BC
                            WRITE(*,'(A,3(2x,g17.8))') 'VEC = ', DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNTBCIN)%VEC(1:3)
                            WRITE(*,'(A,3(2x,g17.8))') 'CNOT = ', DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNTBCIN)%CNOT(1:3)
                            WRITE(*,'(A,3(2x,g17.8))') 'NORMAL OLD = ', DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNTBCIN)%NORMAL(:)
                            WRITE(*,'(A,3(2x,g17.8))') 'NORMAL NEW = ', WALL_NORM_FINAL(:) 
                                  
                            WRITE(*,'(A,3(2x,i4))') 'I, J, K NAT = ', I, J, K
                            WRITE(*,'(A,3(2x,i4))') 'I, J, K OLD = ', I_OF(IJK_TEMP2), J_OF(IJK_TEMP2), K_OF(IJK_TEMP2)
                            WRITE(*,'(A,3(2x,i4))') 'I, J, K NEW = ', I_OF(WALL_CELL_ID_FINAL), J_OF(WALL_CELL_ID_FINAL), K_OF(WALL_CELL_ID_FINAL)                      
                            !read(*,*)
                            RETURN
                         end IF

500                      continue 
                      end IF CUTFACELINE
                   end DO
                      
                   !get out of the vertloop
                   EXIT 
                   !do not return from this subroutine. Instead, add this new bc 
                end IF
                
             END IF
          end DO VERTLOOP
             
          IF(.NOT.REJECT_NEW.AND.(FORCE_INTPOINT_INSIDE.AND.(.NOT.ATLEAST_ONEINSIDE))) THEN
             !redo the checking but this time do not require the 
             !intersection point to lie inside the cell 
             !DEBUG_LOCAL  = .false.
             
             IF(DMP_LOG.and.debug_local) WRITE(UNIT_LOG, 2002) I, J, K, & 
             & I_OF(IJK_TEMP), J_OF(IJK_TEMP), K_OF(IJK_TEMP), & 
             & I_OF(WALL_CELL_ID), J_OF(WALL_CELL_ID), K_OF(WALL_CELL_ID), &
             & REFP_S(IJK_TEMP, 1:3), REFP_S(WALL_CELL_ID, 1:3) 

             IF(PRINT_DES_SCREEN.and.debug_local) WRITE(*, 2002) I, J, K, & 
             & I_OF(IJK_TEMP), J_OF(IJK_TEMP), K_OF(IJK_TEMP), & 
             & I_OF(WALL_CELL_ID), J_OF(WALL_CELL_ID), K_OF(WALL_CELL_ID), &
             & REFP_S(IJK_TEMP, 1:3), REFP_S(WALL_CELL_ID, 1:3) 
             
             !Jan 23, 2012. 
             !Rahul: Don't rerun with this condition. This can lead to 
             !wrong conclusions in regions of sharp curvature.  

             !FORCE_INTPOINT_INSIDE = .FALSE. 
             
             !GOTO 300   
             
             !midway point lies mid-way on the line joining the reference points 
             !on the two cut-faces in question. 

             MIDWAY_PT(1:3) = 0.5d0*(REFP_S(IJK_TEMP, 1:3) + REFP_S(WALL_CELL_ID, 1:3))
          
             CALL GET_DEL_H_DES(IJK_TEMP,'SCALAR',MIDWAY_PT(1), MIDWAY_PT(2), MIDWAY_PT(3), & 
             & DIST_OLD, NORM1_OLD, NORM2_OLD, NORM3_OLD, .true.)
             
             CALL GET_DEL_H_DES(WALL_CELL_ID,'SCALAR',MIDWAY_PT(1), MIDWAY_PT(2), MIDWAY_PT(3), & 
             & DIST_NEW, NORM1_NEW, NORM2_NEW, NORM3_NEW, .true.)

             IF(DIST_OLD*DIST_NEW.LT.ZERO.and.debug_local)THEN
                IF(DMP_LOG)  WRITE(UNIT_LOG,2000) DIST_OLD, DIST_NEW,I,J,K, &
                & I_OF(IJK_TEMP), J_OF(IJK_TEMP), K_OF(IJK_TEMP), & 
                & I_OF(WALL_CELL_ID), J_OF(WALL_CELL_ID), K_OF(WALL_CELL_ID), & 
                & REFP_S(IJK_TEMP, :), REFP_S(WALL_CELL_ID, :),  & 
                & NORM1_OLD, NORM2_OLD, NORM3_OLD, & 
                NORM1_NEW, NORM2_NEW, NORM3_NEW, MIDWAY_PT(:)
                
                IF(PRINT_DES_SCREEN) WRITE(*,2000) DIST_OLD, DIST_NEW,I,J,K, &
                & I_OF(IJK_TEMP), J_OF(IJK_TEMP), K_OF(IJK_TEMP), & 
                & I_OF(WALL_CELL_ID), J_OF(WALL_CELL_ID), K_OF(WALL_CELL_ID), & 
                & REFP_S(IJK_TEMP, :), REFP_S(WALL_CELL_ID, :),  & 
                & NORM1_OLD, NORM2_OLD, NORM3_OLD, & 
                NORM1_NEW, NORM2_NEW, NORM3_NEW, MIDWAY_PT(:) 
                !CALL mfix_exit(myPE)
             ENDIF
          
             IF(DIST_OLD.LT.ZERO.OR.DIST_NEW.LT.ZERO) RETURN 
           
          ENDIF
       end IF CUT_CELL_SPECIAL
       
       !gave him no chance, but made it here...
       !.. deserves to be added :-)

       ADD_THIS_CUT_FACE  = .TRUE. 
       IF(TRIM(DES_WALL_TYPE_FINAL).eq.'CUT_FACE_LINE') ADD_THIS_CUT_FACE  = .false. 
       RETURN
       !as this case will be added below in this routine only. 

       
       !reset the wall variables just in case if they were 
       !modified during the call to this routine (like, for example, 
       !the cut-face-line BC changed the attributes. 
       WALL_NORM(:) = WALL_NORM_FINAL(:) 
       WALL_CELL_ID = WALL_CELL_ID_FINAL 
       DES_WALL_TYPE = DES_WALL_TYPE_FINAL
       

       
       SELECT CASE (TRIM(DES_WALL_TYPE_FINAL)) 
          !Rahul:
          !for the case of cut_face_line, add this BC 
          !in this routine only as cnot is not a global
          !array. This can be fixed later but right not 
          !i do not care as this case will not pop up
          !for the MPPIC case. For the DEM, this can be
          !improved upon later. 
       CASE('CUT_FACE_LINE')
          ADD_THIS_CUT_FACE  = .FALSE. 

          COUNT_BC = COUNT_BC + 1
          DES_CELLWISE_BCDATA(CELL_ID)%COUNT_DES_BC = COUNT_BC 
          
          IF(COUNT_BC.GT.MAX_DES_BC_CELL) THEN
             IF(DMP_LOG) THEN 
                WRITE(UNIT_LOG,*) 'ERROR MESSAGE FROM CHECK_IF_TOADD_THIS_CUT_FACE'
                WRITE(UNIT_LOG,*) 'INCREASE MAX_DES_BC_CELL from the current value of', MAX_DES_BC_CELL
                WRITE(UNIT_LOG,*) 'TERMINAL ERROR: STOPPING'
             ENDIF
             CALL mfix_exit(myPE)
          ENDIF

          DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT_BC)%IJK_SCAL = WALL_CELL_ID_FINAL
          
          ALLOCATE(DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT_BC)%NORMAL(DIMN))
          
          DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT_BC)%NORMAL(:) = WALL_NORM_FINAL(:)

          DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT_BC)%DES_BC_TYPE = DES_WALL_TYPE_FINAL 
          
          !for this special case, allocate the line attributes in 3-d
          ALLOCATE(DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT_BC)%CNOT(3))
                
          ALLOCATE(DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT_BC)%VEC(3))

          DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT_BC)%CNOT(:) = CNOT_ARR(:)
          DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT_BC)%VEC(:) = VEC_ARR(:)

       CASE DEFAULT
          RETURN 
       END SELECT



 2000  FORMAT(/1X,70('*'),//,1X,  & 
       & 'ERROR IN CHECK_IF_TOADD_THIS_CUT_FACE', /1x, & 
       & 'MID WAY POINT BETWEEN TWO ADJACENT CUT FACES NOT GIVING CONSISTENT DISTANCE FROM THE PLANES', /1x, &
       & 'DIST_OLD AND DIST_NEW = ', 2(2x, g17.8), /1X, &
       & 'I, J, K NAT = ', 3(2x, i5), /1X, &
       & 'I, J, K OLD = ', 3(2x, i5), /1X, &
       & 'I, J, K NEW = ', 3(2x, i5), /1X, &
       & 'REF PT OLD  = ', 3(2x, g17.8), /1X, &
       & 'REF PT NEW  = ', 3(2x, g17.8), /1X, &
       & 'NORM OLD  = ', 3(2x, g17.8), /1X, &
       & 'NORM NEW  = ', 3(2x, g17.8), /1X, &
       & 'MID WAY POINT  = ', 3(2x, g17.8), /1X, &
       & 'BOTH DISTANCES SHOULD HAVE BEEN ETIHER POSITIVE OR NEGATIVE  = ',  /1X, &
       & 'BUT ONE POSITIVE AND ONE NEGATIVE IS NOT RIGHT  = ',  /1X, &
       & 'NOT STOPPING BUT MONITOR THIS WARNING', &
       & /1X,70('*')/)
     
 2001  FORMAT(/1X,70('*'),//,1X,  & 
       & 'ERROR IN CHECK_IF_TOADD_THIS_CUT_FACE', /1x, & 
       & 'COMPARING AGAINST A CUT-FACE NOT NATIVE TO THIS CUT-CELL', /1x, &
       & 'OR THE FIRST WALL BC IN THIS CELL NOT A CUT-FACE',  /1X, &
       & 'TOTAL BCs = ', i5, /1X, &
       & 'I, J, K NAT = ', 3(2x, i5), /1X, &
       & 'I, J, K OLD = ', 3(2x, i5), /1X, &
       & 'I, J, K NEW = ', 3(2x, i5), /1X, &
       & 'TERMINAL ERROR: STOPPING', &
       & /1X,70('*')/)
       

 2002  FORMAT(/1X,70('*'),//,1X,  & 
       & 'NOT REJECTED BUT NOT A SINGLE INTERSECTION POINT IN THE CELL', /1x, & 
       & 'CHECKING BASED ON THE POINT MID-WAY BETWEEN CUT FACES IN QUESTION', /1x, &
       & 'I, J, K NAT = ', 3(2x, i5), /1X, &
       & 'I, J, K OLD = ', 3(2x, i5), /1X, &
       & 'I, J, K NEW = ', 3(2x, i5), /1X, &
       & 'REF PT OLD  = ', 3(2x, g17.8), /1X, &
       & 'REF PT NEW  = ', 3(2x, g17.8), /1X, & 
       & /1X,70('*')/)       

     END SUBROUTINE CHECK_IF_TOADD_THIS_CUT_FACE
     
       SUBROUTINE DES_ADD_WALLBC_TO_CELL(CELL_ID, WALL_CELL_ID, WALL_NORM, DES_WALL_TYPE)
       USE param1
       USE funits
       USE run
       USE compar      
       USE discretelement
       USE cutcell

       USE indices
       USE physprop
       USE parallel
       USE geometry

       IMPLICIT NONE 
       !Cell_id is the IJK of the cell to which the new bc is being considered for addition
       !wall_cell_id is the IJK of the cell where the new bc is present 
       INTEGER, INTENT(IN) :: CELL_ID
       INTEGER, INTENT(IN) :: WALL_CELL_ID
       DOUBLE PRECISION, DIMENSION(DIMN), INTENT(IN) :: WALL_NORM
       CHARACTER*100, INTENT(IN) :: DES_WALL_TYPE 
       INTEGER :: COUNT_BC, IJK_TEMP, COUNT , IJK
       CHARACTER*100 :: WALL_TYPE_TEMP
       DOUBLE PRECISION :: NORMAL_TEMP(3)
       LOGICAL :: SMALL_CELL 
       INCLUDE 'function.inc'
       COUNT_BC = DES_CELLWISE_BCDATA(CELL_ID)%COUNT_DES_BC
       
       SMALL_CELL = .false. 

       IF(CARTESIAN_GRID) THEN 
          IF(SMALL_CELL_AT(CELL_ID)) SMALL_CELL = .true.
       ENDIF
       DO COUNT = 1, COUNT_BC 
          IJK_TEMP = DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT)%IJK_SCAL
          WALL_TYPE_TEMP = DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT)%DES_BC_TYPE
          IF((IJK_TEMP.EQ.WALL_CELL_ID).AND.(WALL_TYPE_TEMP.EQ.DES_WALL_TYPE).and.(.not.SMALL_CELL)) THEN 
             !for small cells, the various bc will all have
             !the same wall_cell_id = cell_id
             RETURN
          ENDIF
       ENDDO


       COUNT_BC = COUNT_BC + 1
       DES_CELLWISE_BCDATA(CELL_ID)%COUNT_DES_BC = COUNT_BC 


       IF(COUNT_BC.GT.MAX_DES_BC_CELL) THEN
          IF(DMP_LOG) THEN 
             WRITE(UNIT_LOG, 1036) MAX_DES_BC_CELL, CELL_ID, I_OF(CELL_ID), J_OF(CELL_ID), K_OF(CELL_ID), IS_ON_myPE_owns(I_OF(CELL_ID), J_OF(CELL_ID), K_OF(CELL_ID))
             IF(myPE.eq.pe_IO) WRITE(*, 1036) MAX_DES_BC_CELL, CELL_ID, I_OF(CELL_ID), J_OF(CELL_ID), K_OF(CELL_ID), IS_ON_myPE_owns(I_OF(CELL_ID), J_OF(CELL_ID), K_OF(CELL_ID))

             WRITE(UNIT_LOG, *) 'EXISTING BCs FOR THIS CELL', '   CELL CUT_CELL ?', CUT_CELL_AT(CELL_ID)
             IF(myPE.eq.pe_IO) WRITE(*, *) 'EXISTING BCs FOR THIS CELL', '   CELL CUT_CELL ?', CUT_CELL_AT(CELL_ID)
             DO COUNT = 1, COUNT_BC - 1 
                IJK_TEMP = DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT)%IJK_SCAL 
                WALL_TYPE_TEMP = DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT)%DES_BC_TYPE
                NORMAL_TEMP = ZERO 
                IF(TRIM(WALL_TYPE_TEMP).eq.'CUT_FACE') NORMAL_TEMP(1:3) = NORMAL_S(IJK_TEMP,1:3)

                WRITE(UNIT_LOG,1037) COUNT, TRIM(WALL_TYPE_TEMP), IJK_TEMP, I_OF(IJK_TEMP), J_OF(IJK_TEMP), K_OF(IJK_TEMP), NORMAL_TEMP(1:3)
                IF(myPE.eq.pe_IO) WRITE(*,1037) COUNT, TRIM(WALL_TYPE_TEMP), IJK_TEMP, I_OF(IJK_TEMP), J_OF(IJK_TEMP), K_OF(IJK_TEMP), NORMAL_TEMP(1:3)
             ENDDO
          ENDIF
          CALL mfix_exit(myPE)
       ENDIF


1036   FORMAT(/1X,70('*'),//,10X,  & 
            & 'ERROR MESSAGE FROM DES_ADD_WALLBC_TO_CELL', /10x, & 
            & 'INCREASE MAX_DES_BC_CELL from the current value of', i3, /10x, &
            & 'HAPPENING FOR CELL IJK, I, J, K = ', 4(2x, i5), /10X, &
            & 'TERMINAL ERROR: STOPPING', /10X,  &
            & 'IS ON MY PROCESSOR ? = ', L5,  &
            & /1X,70('*')/)

1037   FORMAT(  & 
            & 'COUNT_BC: ', i3, /10x, & 
            & 'BC_TYPE : ', A, /10x, &
            & 'BC_IJK  : ', i10, /10X, &
            & 'BC_I    : ', i4, /10X, &
            & 'BC_J    : ', i4, /10X, &
            & 'BC_K    : ', i4, /10X, &
            & 'NORMAL  : ', 3(2x,g17.8))

       SELECT CASE (TRIM(DES_WALL_TYPE)) 
       CASE('NORMAL_WALL')
          
          DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT_BC)%IJK_SCAL = WALL_CELL_ID 
          
          ALLOCATE(DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT_BC)%NORMAL(DIMN))

          DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT_BC)%NORMAL(:) = WALL_NORM (:)
          
          DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT_BC)%DES_BC_TYPE = DES_WALL_TYPE  
       CASE('SMALL_CELL')
                    
          DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT_BC)%IJK_SCAL = WALL_CELL_ID 
          
          ALLOCATE(DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT_BC)%NORMAL(DIMN))

          DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT_BC)%NORMAL(:) = WALL_NORM (:)
          
          DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT_BC)%DES_BC_TYPE = DES_WALL_TYPE  

       CASE('PERIODIC_WALL')

          DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT_BC)%IJK_SCAL = WALL_CELL_ID 

          ALLOCATE(DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT_BC)%NORMAL(DIMN))

          DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT_BC)%NORMAL(:) = WALL_NORM (:)

          DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT_BC)%DES_BC_TYPE = DES_WALL_TYPE  

       CASE('CUT_FACE')
          DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT_BC)%IJK_SCAL = WALL_CELL_ID 

          DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT_BC)%DES_BC_TYPE = DES_WALL_TYPE 
          !no need to store the normal, as depening on the particle position, the 
          !new normal from cut-face to particle center will be calculated by 
          !get_del_h_des
          
!!$       CASE('CUT_FACE_LINE')
!!$          
!!$          DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT_BC)%IJK_SCAL = WALL_CELL_ID 
!!$          
!!$          ALLOCATE(DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT_BC)%NORMAL(DIMN))
!!$          
!!$          DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT_BC)%NORMAL(:) = WALL_NORM (:)
!!$
!!$          DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT_BC)%DES_BC_TYPE = DES_WALL_TYPE  
!!$          
!!$          !for this special case, allocate the line attributes in 3-d
!!$          ALLOCATE(DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT_BC)%CNOT(3))
!!$                
!!$          ALLOCATE(DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT_BC)%VEC(3))
!!$
!!$          DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT_BC)%CNOT(:) = CNOT_ARR(:)
!!$          DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT_BC)%VEC(:) = VEC_ARR(:)
       CASE DEFAULT
          IF(DMP_LOG.OR.myPE.eq.pe_IO) THEN 
             WRITE(UNIT_LOG, 1035) 'DES_ADD_WALLBC_TO_CELL', TRIM(DES_WALL_TYPE)
             WRITE(*, 1035) 'DES_ADD_WALLBC_TO_CELL', TRIM(DES_WALL_TYPE)

          ENDIF
          CALL mfix_exit(myPE)
       END SELECT

 1035    FORMAT(/1X,70('*'),//,10X,  & 
         & 'ERROR IN SUBROUTINE....', A, /,10X, &
         & 'THIS BC....(', A,')....FOR PARTICLE-WALL CANT BE ADDED HERE   ', /, 10X, & 
         & 'CHECK THIS CODE AGAIN' , /, 10X, & 
         & 'ACCEPTABLE TYPES ARE:', /10x, & 
         & 'NORMAL_WALL' , /10X, &
         & 'SMALL_FACE', /10X, & 
         & 'PERIODIC_WALL', /10X, & 
         & 'CUT_FACE', /10X, & 
         & 'TERMINAL ERROR: STOPPING', &
         & /1X,70('*')/)


      END SUBROUTINE DES_ADD_WALLBC_TO_CELL

      SUBROUTINE CHECK_IF_GHOST_CELL_WALL_NEEDED(CELL_ID, WALL_CID, XCORD, YCORD, ZCORD, PERP_DIST, NORM_VEC, ADD_GHOST_CELL_WALL)

        USE discretelement
        USE indices
        USE physprop
        USE parallel
        USE geometry
        USE compar      

       IMPLICIT NONE 
       
       DOUBLE PRECISION, INTENT(IN) :: XCORD(4), YCORD(4), ZCORD(4)
       DOUBLE PRECISION, INTENT(OUT) :: PERP_DIST(4), NORM_VEC(4, 3)
       LOGICAL, INTENT(OUT) :: ADD_GHOST_CELL_WALL 
       INTEGER, INTENT(IN) :: CELL_ID, WALL_CID !Wall cell id 
       !local variables 
       INTEGER ::  COUNT, COUNT_WALL_LOOP, IJK
       DOUBLE PRECISION :: NORMAL_MAG, diagonal 
       
       INCLUDE 'function.inc'
       
       IF(DIMN.EQ.2) THEN 
          COUNT_WALL_LOOP = 2   !This is used to check the vertices of the cell defined 
                                !by ghost cell to determine if is this wall is inside or outside the domain
       ELSE
          COUNT_WALL_LOOP = 4
       ENDIF

       Diagonal = DX(I_OF(CELL_ID))**2 + DY(J_OF(CELL_ID))**2 
       IF(DIMN.eq.3) Diagonal = Diagonal + DZ(K_OF(CELL_ID))**2
       Diagonal = sqrt(diagonal) 
       
       DO COUNT = 1, COUNT_WALL_LOOP
          CALL GET_DEL_H_DES(CELL_ID,'SCALAR',XCORD(COUNT),YCORD(COUNT),ZCORD(COUNT), & 
          & PERP_DIST(COUNT), NORM_VEC(COUNT, 1), NORM_VEC(COUNT, 2), NORM_VEC(COUNT, 3), .true.)
       
          IF(PERP_DIST(COUNT).GT.0.01*diagonal) THEN !ZERO) THEN 
             IF(CELL_ID.EQ.-1) THEN 
                WRITE(*,*) 'GHOST CELL VERTEX FOUND INSIDE OF THE CUT CELL FACE'
                WRITE(*,*) 'WALL CELL = ', WALL_CID, I_OF(WALL_CID), J_OF(WALL_CID), K_OF(WALL_CID)

                WRITE(*,'(A10,3(2x,g17.8))') 'XYZ COR = ', XCORD(COUNT), YCORD(COUNT), ZCORD(COUNT)
                
                WRITE(*,'(A10,3(2x,g17.8))') 'DIST = ', PERP_DIST(COUNT)

                WRITE(*,'(A10,3(2x,g17.8))') 'NORMAL VEC = ', NORM_VEC(COUNT, 1:3)
             
             ENDIF
             ADD_GHOST_CELL_WALL = .TRUE. 
             EXIT !no need to check the remaining vertices 
             
          ENDIF
       ENDDO
       
       RETURN 
       END SUBROUTINE CHECK_IF_GHOST_CELL_WALL_NEEDED
       
       
      SUBROUTINE REORDER_BCS_FOR_CUTCELL(CELL_ID)
        !the purpose of this routine is to re-order the wall bc's for cut-cells that include the cut_face_line.
        !bring the cut_face_line to the top of wall bcs. if the particle is found in contact with the cut_face_line,
        !then the cut_face will be restricted to be used as a bc. This will prevent double push to the particle from 
        !both cut-face and cut-face-line 
        USE param1
        USE funits
        USE run
        USE compar      
        USE discretelement
        USE cutcell
        
        USE indices
        USE physprop
        USE parallel
        USE geometry
        IMPLICIT NONE 
        !Cell_id is the IJK of the cell to which the new bc is being considered for addition
        !wall_cell_id is the IJK of the cell where the new bc is present 
        INTEGER, INTENT(IN) :: CELL_ID
        !Local variables
        INTEGER, PARAMETER :: MAX_BC = 10
        INTEGER :: COUNT_BC, IJK_TEMP, COUNT, COUNT_BC_OLD
        INTEGER ::  IJK_ARR_TEMP(MAX_BC), COUNT_BC_NEW, MAX_DIST_LOC(1), MAX_DIST2_LOC(1),  COUNT_CUT_FACE_LINE
        
        CHARACTER*100 :: WALL_TYPE_ARR_TEMP(MAX_BC), WALL_TYPE_ARR_TEMP2(MAX_BC), WALL_TYPE_TEMP 
        DOUBLE PRECISION :: NORMAL_ARR_TEMP(MAX_BC, DIMN), CNOT_ARR_TEMP(MAX_BC, 3), VEC_ARR_TEMP(MAX_BC,3)
        LOGICAL :: REORDER, BC_REORDER_ARR(MAX_BC)


        COUNT_BC = DES_CELLWISE_BCDATA(CELL_ID)%COUNT_DES_BC
        NORMAL_ARR_TEMP = ZERO
        CNOT_ARR_TEMP = ZERO
        VEC_ARR_TEMP = ZERO
        REORDER = .FALSE.

        !first check if re-ordering is even needed
        DO COUNT = 1, COUNT_BC 
           IJK_TEMP = DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT)%IJK_SCAL
           WALL_TYPE_TEMP = DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT)%DES_BC_TYPE
           IF(TRIM(WALL_TYPE_TEMP).EQ.'CUT_FACE_LINE') THEN
              REORDER = .TRUE.
              EXIT !THE COUNT LOOP
           end IF
        end DO
        
        IF(.NOT.REORDER) RETURN 

        BC_REORDER_ARR = .FALSE.
        
        IF(COUNT_BC.GT.MAX_BC) THEN 
           WRITE(UNIT_LOG, '(A, i5, A, i5, A,/, A)') 'SIZE OF TEMP ARRAYS', MAX_BC, ' IN REORDER_BCS_FOR_CUTCELL SMALLER THAN COUNT_BC ( = ', COUNT_BC, ' ) FOR THIS CELL', 'TERMINAL ERROR: STOPPING'
           IF(PRINT_DES_SCREEN) WRITE(*, '(A, i5, A, i5, A,/, A)') 'SIZE OF TEMP ARRAYS', MAX_BC, ' IN REORDER_BCS_FOR_CUTCELL SMALLER THAN COUNT_BC ( = ', COUNT_BC, ' ) FOR THIS CELL', 'TERMINAL ERROR: STOPPING'
           CALL mfix_exit(myPE)
        ENDIF
        
        DO COUNT = 1, COUNT_BC 
           IJK_TEMP = DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT)%IJK_SCAL
           WALL_TYPE_TEMP = DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT)%DES_BC_TYPE
           IJK_ARR_TEMP(COUNT) = IJK_TEMP
           WALL_TYPE_ARR_TEMP(COUNT) = WALL_TYPE_TEMP
           
           SELECT CASE (TRIM(WALL_TYPE_TEMP))
           CASE('NORMAL_WALL')
              NORMAL_ARR_TEMP(COUNT,:) = DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT)%NORMAL(:)
              DEALLOCATE(DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT)%NORMAL)
              
           CASE('CUT_FACE')
              !nothing more to do in this case. Normal is obtained from the cut-cell data structures

           CASE('CUT_FACE_LINE')
              NORMAL_ARR_TEMP(COUNT,:) = DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT)%NORMAL(:)
              DEALLOCATE(DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT)%NORMAL)

              CNOT_ARR_TEMP(COUNT,:) = DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT)%CNOT(:)
              DEALLOCATE(DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT)%CNOT)
              
              VEC_ARR_TEMP(COUNT,:) = DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT)%VEC(:)
              DEALLOCATE(DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT)%VEC)
              
              !mark this to be reordered 
              BC_REORDER_ARR(COUNT) = .TRUE.

              
           CASE DEFAULT
              WRITE(*,*)' EROR IN SUBROUTINE REORDER_BCS_FOR_CUTCELL'
              WRITE(*,*)'UNKNOWN TYPE OF DES_WALL_TYPE:',TRIM(WALL_TYPE_ARR_TEMP(COUNT))
              WRITE(*,*)'ACCEPTABLE TYPES ARE:' 
              WRITE(*,*)'NORMAL_WALL' 
              WRITE(*,*)'CUT_FACE' 
           END SELECT
        END DO
        
        COUNT_BC_OLD = COUNT_BC
        COUNT_BC = 0
        DO COUNT = 1, COUNT_BC_OLD
           
           IF(BC_REORDER_ARR(COUNT)) THEN
              COUNT_BC = COUNT_BC+1
              IF(TRIM(WALL_TYPE_ARR_TEMP(COUNT)).NE.'CUT_FACE_LINE') THEN
                 WRITE(*,*) 'ERROR IN REORDER_BCS_FOR_CUTCELL'
                 WRITE(*,*) 'WAS EXPECting cut_face_line bc'
                 WRITE(*,*) 'BUT FOUND            ', TRIM(WALL_TYPE_ARR_TEMP(COUNT))
                 
                 WRITE(*,*) 'TERMINAL ERROR: STOPPING'
                 STOP
              end IF
              
              DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT_BC)%IJK_SCAL = IJK_ARR_TEMP(COUNT)
              
              ALLOCATE(DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT_BC)%NORMAL(DIMN))
              
              DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT_BC)%NORMAL(:) = NORMAL_ARR_TEMP(COUNT,:)
              
              DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT_BC)%DES_BC_TYPE = WALL_TYPE_ARR_TEMP(COUNT)
              
              !for this special case, allocate the line attributes in 3-d
              ALLOCATE(DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT_BC)%CNOT(3))
              
              ALLOCATE(DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT_BC)%VEC(3))
              
              DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT_BC)%CNOT(:) = CNOT_ARR_TEMP(COUNT,:)
              DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT_BC)%VEC(:) = VEC_ARR_TEMP(COUNT,:)
              
           end IF
        end DO
        
        DO COUNT = 1, COUNT_BC_OLD
           
           IF(.NOT.BC_REORDER_ARR(COUNT)) THEN
              COUNT_BC = COUNT_BC+1
              
              SELECT CASE (TRIM(WALL_TYPE_ARR_TEMP(COUNT)))
              CASE('NORMAL_WALL')

                 DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT_BC)%IJK_SCAL = IJK_ARR_TEMP(COUNT)
                 
                 ALLOCATE(DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT_BC)%NORMAL(DIMN))
                 
                 DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT_BC)%NORMAL(:) = NORMAL_ARR_TEMP(COUNT, :)

                 DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT_BC)%DES_BC_TYPE = WALL_TYPE_ARR_TEMP(COUNT)
          
              CASE('CUT_FACE')
                 DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT_BC)%IJK_SCAL =  IJK_ARR_TEMP(COUNT)

                 DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT_BC)%DES_BC_TYPE =  WALL_TYPE_ARR_TEMP(COUNT)
                 
              CASE DEFAULT
                 
                 WRITE(*,*)' ERROR IN SUBROUTINE REORDER_BCS_FOR_CUTCELL:'
                 WRITE(*,*)'UNKNOWN TYPE OF DES_WALL_TYPE:  ', TRIM(WALL_TYPE_ARR_TEMP(COUNT))
                 WRITE(*,*)'ACCEPTABLE TYPES ARE:' 
                 WRITE(*,*)'NORMAL_WALL' 
                 WRITE(*,*)'CUT_FACE' 
                 WRITE(*,*) 'TERMINAL ERROR: STOPPING'
                 STOP
              END SELECT
           end IF
        end DO
        
        
        IF(COUNT_BC.NE.COUNT_BC_OLD)  THEN 
           
           WRITE(*,*)' ERROR IN SUBROUTINE REORDER_BCS_FOR_CUTCELL:'
           WRITE(*,*)'NUMBER OF NEW AND OLD BCs DIFFERENT AFTER REORDERING ', COUNT_BC, COUNT_BC_OLD
           WRITE(*,*) 'TERMINAL ERROR: STOPPING'
           STOP
        end IF
        
        return 
      END SUBROUTINE REORDER_BCS_FOR_CUTCELL
      
      SUBROUTINE CHECK_NON_CUTCELLBC(CELL_ID) 
      USE param1
      USE funits
      USE run
      USE compar      
      USE discretelement
      USE cutcell
 
      USE indices
      USE physprop
      USE parallel
      USE geometry
      
       IMPLICIT NONE 
       !Cell_id is the IJK of the cell to which the new bc is being considered for addition
       !wall_cell_id is the IJK of the cell where the new bc is present 
       INTEGER, INTENT(IN) :: CELL_ID
       !Local variables
       
       INTEGER :: COUNT_BC, IJK_TEMP, IJK_TEMP2, IJK_TEMP3, COUNT, I, J, K, I1, J1, K1, COUNT_VERT, DEL_COUNT,  COUNT_FACE_LOOP, COUNTBCIN, COUNT3
       CHARACTER*100 :: WALL_TYPE_TEMP, WALL_TYPE_TEMP2, WALL_TYPE_TEMP3
       LOGICAL :: INT_OUTSIDECELL, REJECT_NEW, ATLEAST_ONEINSIDE, FORCE_INTPOINT_INSIDE
              
       DOUBLE PRECISION :: NORM1_NEW, NORM2_NEW, NORM3_NEW
       DOUBLE PRECISION :: NORM1_OLD, NORM2_OLD, NORM3_OLD
       DOUBLE PRECISION :: NORM_SQ, DIST_OLD, DIST_NEW, DIST_NEW_INT
       DOUBLE PRECISION :: XCOR(8), YCOR(8), ZCOR(8)!store the vertices of the scalar cell with id cell_id

       DOUBLE PRECISION :: XINT_OLD, YINT_OLD, ZINT_OLD, XMINCOR, XMAXCOR, YMINCOR, YMAXCOR, ZMINCOR, ZMAXCOR, DIAGONAL, DIAGONAL2, TMP_DOT
       INTEGER ::  COUNT_VERT_LOOP, COUNT2, IJK_ARR_TEMP(6), IJK_ARR_TEMP2(6), COUNT_BC_NEW, COUNT_CUT_FACE_LINE

       CHARACTER*100 :: WALL_TYPE_ARR_TEMP(6), WALL_TYPE_ARR_TEMP2(6)
       DOUBLE PRECISION :: NORMAL_ARR_TEMP(6, DIMN), NORMAL_ARR_TEMP2(6, DIMN), CNOT_ARR_TEMP(6, 3), VEC_ARR_TEMP(6,3), DIFF_DIST, NORMAL_OLD(DIMN)

       DOUBLE PRECISION :: COEFF_A1, COEFF_B1, COEFF_C1, COEFF_A2, COEFF_B2, COEFF_C2, ROOTX, ROOTY, DET, NORM_AVG(DIMN), NORM_MAG, TOL_CONVEX_ANGLE_DES, TOL_CONVEX_ANGLE
       LOGICAL :: BC_DEL_ARR(6), ADD_CUT_FACE_LINE, OUT_OF_PLANE1, OUT_OF_PLANE2
       
       DOUBLE PRECISION :: XREF_FACE(6), YREF_FACE(6), ZREF_FACE(6), NORMAL_FACE(6,3)
       DOUBLE PRECISION :: XREF_FACE2(6), YREF_FACE2(6), ZREF_FACE2(6), NORMAL_FACE2(6,3)
       
       DOUBLE PRECISION :: TOL_LINEPLANE_DIST_DES, TOL_LINEINT, TOL_LINEINT2, PTXYZ_ON_LINE(3), DISTFROMLINE, DISTFROMLINE2
       
       INTEGER :: DIST_CORDS_1, DIST_CORDS_2, DIST12X, DIST12Y, DIST12Z
              
       TOL_LINEPLANE_DIST_DES = 0.01

       DEL_COUNT = 0
       COUNT_BC = DES_CELLWISE_BCDATA(CELL_ID)%COUNT_DES_BC
       BC_DEL_ARR = .false.
       TOL_CONVEX_ANGLE_DES  = 10.d0
       TOL_CONVEX_ANGLE = ASIN(1.d0)*TOL_CONVEX_ANGLE_DES/90.d0
       TOL_CONVEX_ANGLE  =  COS(TOL_CONVEX_ANGLE)
       
       WRITE(*,*) 'PLANES WHOSE NORMALS ENCOLSE MORE THAN', TOL_CONVEX_ANGLE,' WILL BE REMOVED AND A CUT_FACE_LINE CREATED'

       
       DO COUNT = 1, COUNT_BC 
          IJK_TEMP = DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT)%IJK_SCAL
          WALL_TYPE_TEMP = DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT)%DES_BC_TYPE
          IJK_ARR_TEMP(COUNT) = IJK_TEMP
          WALL_TYPE_ARR_TEMP(COUNT) = WALL_TYPE_TEMP
          IF(TRIM(WALL_TYPE_TEMP).EQ.'CUT_FACE') THEN 
             !no normal 
          ELSE
             NORMAL_ARR_TEMP(COUNT,:) = DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT)%NORMAL(:)
             DEALLOCATE(DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT)%NORMAL)
          ENDIF
       ENDDO

       WRITE(*,*) 'CHECKING BCS FOR NON_CUT_CELL = ', CELL_ID
       IF(DIMN.EQ.2) THEN 
          COUNT_VERT_LOOP = 4   
          COUNT_FACE_LOOP = 4
       ELSE
          COUNT_VERT_LOOP = 8
          COUNT_FACE_LOOP = 6
       ENDIF
       
       COUNT_CUT_FACE_LINE = 0

       DO COUNT = 1, COUNT_BC 
          IJK_TEMP = DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT)%IJK_SCAL
          WALL_TYPE_TEMP = DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT)%DES_BC_TYPE
          IF(TRIM(WALL_TYPE_TEMP).NE.'CUT_FACE') THEN 
             !A non cut cell can only have a cut-face as one of the bc's
             WRITE(*,*) 'ERROR MESSAGE FROM CHECK_NON_CUTCELLBC'
             WRITE(*,*) 'A NON CUT-CELL HAS A NON CUT_FACE BC AS ONE OF THE BCs'
             WRITE(*,*) 'TERMINAL ERROR: STOPPING'
             STOP
          ENDIF
             
          I = I_OF(IJK_TEMP)
          J = J_OF(IJK_TEMP)
          K = K_OF(IJK_TEMP)
          
          !store the coordinate information of the scalar cell for this cell
          !the distance of all the cut-faces in cell with id cell_id will be computed from the 
          !vertices of this cell. If any of the cut-faces are found between this cut-face and vertices
          !of the cell this cut-face lies in than that cut-face will be rejected as a bc
          ZCOR(1:4) = ZG_T(K)
          
          XCOR(1) = XG_E(I)
          YCOR(1) = YG_N(J)
          
          XCOR(2) = XG_E(I) - DX(I)
          YCOR(2) = YG_N(J)

          XCOR(3) = XG_E(I)
          YCOR(3) = YG_N(J) - DY(J)
          
          XCOR(4) = XG_E(I) - DX(I)
          YCOR(4) = YG_N(J) - DY(J)
          
          ZCOR(5:8) = ZG_T(K) - DZ(K)
          
          XCOR(5) = XG_E(I)
          YCOR(5) = YG_N(J)
          
          XCOR(6) = XG_E(I) - DX(I)
          YCOR(6) = YG_N(J)

          XCOR(7) = XG_E(I)
          YCOR(7) = YG_N(J) - DY(J)

          XCOR(8) = XG_E(I) - DX(I)
          YCOR(8) = YG_N(J) - DY(J)
          
          XMINCOR = XG_E(I) - DX(I)
          XMAXCOR = XG_E(I) 
          YMINCOR = YG_N(J) - DY(J)
          YMAXCOR = YG_N(J)
          
          ZMINCOR = ZG_T(K) - DZ(K)
          ZMAXCOR = ZG_T(K)

          NORMAL_FACE(6,:) = 0.d0
          !East face 
          XREF_FACE(1) = XG_E(I)
          YREF_FACE(1) = YG_N(J)
          ZREF_FACE(1) = ZG_T(K)
          NORMAL_FACE(1,1) = -1.d0
          
          !west face 
          XREF_FACE(2) = XG_E(I) - DX(I)
          YREF_FACE(2) = YG_N(J)
          ZREF_FACE(2) = ZG_T(K)
          NORMAL_FACE(2,1) = 1.d0
          
          !North face 
          XREF_FACE(3) = XG_E(I)
          YREF_FACE(3) = YG_N(J)
          ZREF_FACE(3) = ZG_T(K)
          NORMAL_FACE(3,2) = -1.d0
          
          
          !South face 
          XREF_FACE(4) = XG_E(I)
          YREF_FACE(4) = YG_N(J) - DY(J)
          ZREF_FACE(4) = ZG_T(K)
          NORMAL_FACE(4,2) = 1.d0
          
          !Top face 
          XREF_FACE(5) = XG_E(I)
          YREF_FACE(5) = YG_N(J) 
          ZREF_FACE(5) = ZG_T(K)
          NORMAL_FACE(5,3) = -1.d0
          
          !Bottom face 
          XREF_FACE(6) = XG_E(I)
          YREF_FACE(6) = YG_N(J) 
          ZREF_FACE(6) = ZG_T(K) - DZ(K)
          NORMAL_FACE(6,3) = 1.d0

          
          Diagonal = (DX(I)**2 + DY(J)**2)
          IF(DIMN.eq.3)  Diagonal = diagonal + DZ(K)**2
          Diagonal = dsqrt(diagonal)
          TOL_LINEINT = - DIAGONAL*TOL_LINEPLANE_DIST_DES  

          DIST_CORDS_1  = ABS(I_OF(IJK_TEMP)-I_OF(CELL_ID)) + ABS(J_OF(IJK_TEMP)-J_OF(CELL_ID))          
          IF(DIMN.Eq.3) DIST_CORDS_1  = DIST_CORDS_1 + ABS(K_OF(IJK_TEMP)-K_OF(CELL_ID))

          DO COUNT2 = COUNT+1, COUNT_BC
          !DO COUNT2 = 1, COUNT_BC
             IF(COUNT2.EQ.COUNT) CYCLE 
             IJK_TEMP2 = DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT2)%IJK_SCAL
             WALL_TYPE_TEMP2 = DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT2)%DES_BC_TYPE

             
             
             DIST_CORDS_2  = ABS(I_OF(IJK_TEMP2)-I_OF(CELL_ID)) + ABS(J_OF(IJK_TEMP2)-J_OF(CELL_ID))          
             IF(DIMN.Eq.3) DIST_CORDS_2  = DIST_CORDS_2 + ABS(K_OF(IJK_TEMP2)-K_OF(CELL_ID))

             DIST12X = ABS(I_OF(IJK_TEMP2)-I_OF(IJK_TEMP))
             DIST12Y = ABS(J_OF(IJK_TEMP2)-J_OF(IJK_TEMP))
             DIST12Z = ABS(K_OF(IJK_TEMP2)-K_OF(IJK_TEMP))

             !only check those cut-faces that are in the adjacent or diagonal cells 
             IF(DIST12X.GT.1.OR.DIST12Y.GT.1.OR.DIST12Z.GT.1) CYCLE
             IF(TRIM(WALL_TYPE_TEMP2).NE.'CUT_FACE') THEN 
             !A non cut cell can only have a cut-face as one of the bc's
                WRITE(*,*) 'ERROR MESSAGE FROM CHECK_NON_CUTCELLBC'
                WRITE(*,*) 'A NON CUT-CELL HAS A NON CUT_FACE BC AS ONE OF THE BCs'
                WRITE(*,*) 'TERMINAL ERROR: STOPPING'
                STOP
             ENDIF

             ATLEAST_ONEINSIDE = .false.
             FORCE_INTPOINT_INSIDE = .TRUE. 

600          CONTINUE 

             VERTLOOP: DO COUNT_VERT = 1, COUNT_VERT_LOOP
                CALL GET_DEL_H_DES(IJK_TEMP,'SCALAR',XCOR(COUNT_VERT), YCOR(COUNT_VERT), &
                & ZCOR(COUNT_VERT), DIST_OLD, NORM1_OLD, NORM2_OLD, NORM3_OLD, .true.)

                IF(DIST_OLD.LT.ZERO) CYCLE !if the vertex in question is not in the fluid domain, then cycle to the next vertex 
                
                XINT_OLD  = XCOR(COUNT_VERT) - DIST_OLD*NORM1_OLD
                YINT_OLD  = YCOR(COUNT_VERT) - DIST_OLD*NORM2_OLD
                ZINT_OLD  = ZCOR(COUNT_VERT) - DIST_OLD*NORM3_OLD
                   
                !at sharp edges the intersection point my lie outside of the cell containing the cut-face... and this point as a result mmight lie on the wrong side of the cut-face being considered
                !if the new cut-face is indeed before the existing cut-face then it will show up this way from other vertices where the intersection point will also lie within the cell 
                INT_OUTSIDECELL = .false.

                IF(XINT_OLD.LT.XMINCOR.OR.XINT_OLD.GT.XMAXCOR.OR.YINT_OLD.LT.YMINCOR.OR.YINT_OLD.GT.YMAXCOR) INT_OUTSIDECELL =.true.
                IF(DIMN.EQ.3.AND.(ZINT_OLD.LT.ZMINCOR.OR.ZINT_OLD.GT.ZMAXCOR)) INT_OUTSIDECELL =.true.
                
                CALL GET_DEL_H_DES(IJK_TEMP2,'SCALAR',XCOR(COUNT_VERT), YCOR(COUNT_VERT), & 
                & ZCOR(COUNT_VERT), DIST_NEW, NORM1_NEW, NORM2_NEW, NORM3_NEW, .true.)


                CALL GET_DEL_H_DES(IJK_TEMP2,'SCALAR', XINT_OLD, YINT_OLD, ZINT_OLD, & 
                & DIST_NEW_INT, NORM1_NEW, NORM2_NEW, NORM3_NEW, .true.)

                IF(DIST_OLD.EQ.DIST_NEW.AND.DIST_OLD.EQ.DIST_NEW_INT) CYCLE


                REJECT_NEW = .FALSE.
                IF(FORCE_INTPOINT_INSIDE) THEN 
                   IF(INT_OUTSIDECELL ) THEN 
                      IF(DIST_NEW.LT.ZERO) REJECT_NEW = .true.
                   ELSE
                      ATLEAST_ONEINSIDE = .TRUE.
                      IF(DIST_NEW_INT.LT.ZERO) REJECT_NEW = .true.
                   ENDIF
                ELSE
                   IF(DIST_NEW.LT.ZERO.OR.DIST_NEW_INT.LT.ZERO)  REJECT_NEW = .true.
                ENDIF

                REJECT_NEW = REJECT_NEW.AND.ATLEAST_ONEINSIDE 
                IF(I_OF(CELL_ID).EQ.-1.and.J_OF(CELL_ID).EQ.-1) THEN 
                   WRITE(*,*) 'TRouble cell'
                   WRITE(*,*) 'I,J,K NAT ', I_OF(CELL_id), j_of(cell_id), CUT_CELL_AT(CELL_ID)
                   WRITE(*,*) 'I,J,K OLD ', I_OF(IJK_TEMP), j_of(IJK_TEMP), CUT_CELL_AT(IJK_TEMP)
                   WRITE(*,*) 'I,J,K NEW ', I_OF(IJK_TEMP2), j_of(IJK_TEMP2), CUT_CELL_AT(IJK_TEMP2)
                   WRITE(*,*) 'VERTEX NUMBER = ', COUNT_VERT
                                         
                   WRITE(*,'(A,3(2x,g17.8))') 'DIST OLD AND NEW=', DIST_OLD, DIST_NEW, DIST_NEW_INT

                   WRITE(*,'(A,3(2x,g17.8))') 'INTERSECTION POINT WITH THE OLD =', XINT_OLD, YINT_OLD, ZINT_OLD
                   WRITE(*,'(A,3(2x,g17.8))') 'NORMAL OLD =', NORM1_OLD, NORM2_OLD, NORM3_OLD
                   WRITE(*,'(A,3(2x,g17.8))') 'NORMAL NEW =', NORM1_NEW, NORM2_NEW, NORM3_NEW

                   WRITE(*,'(A,3(2x,g17.8))') 'EAST NORTH AND TOP OF THE CELL =', XG_E(I), YG_N(J), ZG_T(K)
                   WRITE(*,'(A,3(2x,g17.8))') 'WEST SOUTH AND BOTTOM OF THE CELL =', XG_E(I)-DX(I), YG_N(J)-DY(J), ZG_T(K) - DZ(K)

                   READ(*,*)
                ENDIF
                
                IF(REJECT_NEW) THEN 

                   NORM_SQ = NORM1_OLD*NORM1_NEW + NORM2_OLD*NORM2_NEW + NORM3_OLD*NORM3_NEW
                   IF(NORM_SQ.GT.TOL_CONVEX_ANGLE) THEN

                      !keep any one plane 
                      IF(.NOT.BC_DEL_ARR(COUNT2)) THEN 
                         DEL_COUNT = DEL_COUNT+1
                         BC_DEL_ARR(COUNT2) = .TRUE.
                      ENDIF
                      
                      !exit the vertex do loop 
                      EXIT 
                      
                   ELSE
                      
                         
                      !REJECT both count and count2 and add a new cut-face-line as a bc which is the interestion of these two planes
                      
                      IF(.NOT.BC_DEL_ARR(COUNT)) THEN 
                         DEL_COUNT = DEL_COUNT+1
                         BC_DEL_ARR(COUNT) = .TRUE.
                      ENDIF
                      
                      IF(.NOT.BC_DEL_ARR(COUNT2)) THEN 
                         DEL_COUNT = DEL_COUNT+1
                         BC_DEL_ARR(COUNT2) = .TRUE.
                      ENDIF
                   END IF
                   
                   !here => both planes have been deactivated. 
                   !Add a new bc (cut_face_line) which is the intersection of these two planes 
                   !before deciding to add, make sure the line falls on one of the faces of cell given by count 
                   COUNT_CUT_FACE_LINE =  COUNT_CUT_FACE_LINE+1

                   CALL DES_CROSSPRDCT_3D(VEC_ARR_TEMP(COUNT_CUT_FACE_LINE,:), NORMAL_S(IJK_TEMP,:), NORMAL_S(IJK_TEMP2,:))
                   VEC_ARR_TEMP(COUNT_CUT_FACE_LINE,:) = VEC_ARR_TEMP(COUNT_CUT_FACE_LINE,:)/SQRT(DOT_PRODUCT(VEC_ARR_TEMP(COUNT_CUT_FACE_LINE,:), VEC_ARR_TEMP(COUNT_CUT_FACE_LINE,:)))
                   coeff_a1 = NORMAL_S(IJK_TEMP,1) 
                   coeff_b1 = NORMAL_S(IJK_TEMP,2) 
                   coeff_a2 = NORMAL_S(IJK_TEMP2,1) 
                   coeff_b2 = NORMAL_S(IJK_TEMP2,2) 
          
                   coeff_c1 = DOT_PRODUCT(REFP_S(IJK_TEMP,:), NORMAL_S(IJK_TEMP,:))
                   coeff_c2 = DOT_PRODUCT(REFP_S(IJK_TEMP2,:), NORMAL_S(IJK_TEMP2,:))

                   det = coeff_a1*coeff_b2 - coeff_a2*coeff_b1
                   IF(ABS(DET).LT.0.00001) THEN
                      WRITE(*,*) 'ERROR MESSAGE FROM CHECK_NON_CUTCELLBC'
                      WRITE(*,*) 'PLANES DO NOT INTERSECT AS DETERMINANT NEARLY EQUAL TO ZERO', DET
                      WRITE(*,'(A35, 3(2x,g17.8))') 'COeff for first plane  = ', coeff_a1,coeff_b1, coeff_c1
                      WRITE(*,'(A35, 3(2x,g17.8))') 'COeff for second plane  = ', coeff_a2,coeff_b2, coeff_c2
                      
                      WRITE(*,'(A35, 3(2x,i4))') 'I, J, K of first plane', I_OF(IJK_TEMP), J_OF(IJK_TEMP), K_OF(IJK_TEMP)
                      WRITE(*,'(A35, 3(2x,i4))') 'I, J, K of second plane', I_OF(IJK_TEMP2), J_OF(IJK_TEMP2), K_OF(IJK_TEMP2)
                      WRITE(*,'(A35, 3(2x,g17.8))') 'NORMAL 1 =', NORMAL_S(IJK_TEMP,: )
                      WRITE(*,'(A35, 3(2x,g17.8))') 'NORMAL 2 =', NORMAL_S(IJK_TEMP2,: )
                      
                      WRITE(*,*) 'TERMINAL ERROR: STOPPING'
                      STOP
                   ENDIF

                   ROOTX = (COEFF_C1*COEFF_B2 - COEFF_C2*COEFF_B1)/DET
                   ROOTY = (COEFF_A1*COEFF_C2 - COEFF_A2*COEFF_C1)/DET
                   
                   CNOT_ARR_TEMP(COUNT_CUT_FACE_LINE, 1) = ROOTX
                   CNOT_ARR_TEMP(COUNT_CUT_FACE_LINE, 2) = ROOTY
                   CNOT_ARR_TEMP(COUNT_CUT_FACE_LINE, 3) = ZERO
                   
          
                   NORM_AVG(:) = HALF*(NORMAL_S(IJK_TEMP,:)+NORMAL_S(IJK_TEMP2,:))
                   
                   NORM_MAG = SQRT(DOT_PRODUCT(NORM_AVG, NORM_AVG))
                   NORM_AVG(:) = NORM_AVG(:)/NORM_MAG
                   
                   NORM_MAG = SQRT(DOT_PRODUCT(NORM_AVG, NORM_AVG))
                   IF(NORM_MAG.LT.0.99) THEN
                      WRITE(*,*) 'ERROR MESSAGE FROM CHECK_NON_CUTCELLBC'
                      WRITE(*,*) 'NEW NORMAL MAGNITUDE NE 1', NORM_MAG
                      WRITE(*,*) 'NOMRAL_AVG = ', NORM_AVG(:)
                      WRITE(*,*) 'TERMINAL ERROR: STOPPING'
                      STOP
                   ENDIF
                   ADD_CUT_FACE_LINE = .TRUE.
                   NORMAL_ARR_TEMP2(COUNT_CUT_FACE_LINE,:) = NORM_AVG(:)
                   IJK_ARR_TEMP2(COUNT_CUT_FACE_LINE) = IJK_TEMP
                   WALL_TYPE_ARR_TEMP2(COUNT_CUT_FACE_LINE) = 'CUT_FACE_LINE'
                   
                   WRITE(*,*) 'CHANGED THE BC TO CUT_FACE_LINE IN THE PARAMETERIC FORM vec(CNOT) + t*vec(v)'
                   WRITE(*,'(A35, 3(2x,g17.8))') 'LINE DIRECTION COSINE OR VEC(V)= ', VEC_ARR_TEMP(COUNT_CUT_FACE_LINE,:)
                   WRITE(*,'(A35, 3(2x,g17.8))') 'LINE INTERCEPT OR VEC(CNOT) = ', CNOT_ARR_TEMP(COUNT_CUT_FACE_LINE,:)
                   
                   WRITE(*,'(A35, 3(2x,g17.8))') 'OLD NORMAL 1 =', NORMAL_S(IJK_TEMP,: )
                   WRITE(*,'(A35, 3(2x,g17.8))') 'OLD NORMAL 2 =', NORMAL_S(IJK_TEMP2,: )
                   WRITE(*,'(A35, 3(2x,g17.8))') 'NEW NORMAL  = ', NORMAL_ARR_TEMP2(COUNT_CUT_FACE_LINE,:)
          
                   
                   WRITE(*,'(A35, 3(2x,g17.8))') 'REF PT 1 =', REFP_S(IJK_TEMP,:), IJK_TEMP          
                   WRITE(*,'(A35, 3(2x,g17.8))') 'REF PT 2 =', REFP_S(IJK_TEMP2,:), IJK_TEMP2          
                   WRITE(*,'(A35, 3(2x,g17.8))') 'COeff for first plane  = ', coeff_a1,coeff_b1, coeff_c1
                   WRITE(*,'(A35, 3(2x,g17.8))') 'COeff for second plane  = ', coeff_a2,coeff_b2, coeff_c2
                   
                   !Now do ur best to not add this new bc 

                   I1 = I_OF(IJK_TEMP2)
                   J1 = J_OF(IJK_TEMP2)
                   K1 = K_OF(IJK_TEMP2)

                   NORMAL_FACE2(6,:) = 0.d0
                                !East face 
                   XREF_FACE2(1) = XG_E(I1)
                   YREF_FACE2(1) = YG_N(J1)
                   ZREF_FACE2(1) = ZG_T(K1)
                   NORMAL_FACE2(1,1) = -1.d0
                   
                                !west face 
                   XREF_FACE2(2) = XG_E(I1) - DX(I1)
                   YREF_FACE2(2) = YG_N(J1)
                   ZREF_FACE2(2) = ZG_T(K1)
                   NORMAL_FACE2(2,1) = 1.d0
                   
                                !North face 
                   XREF_FACE2(3) = XG_E(I1)
                   YREF_FACE2(3) = YG_N(J1)
                   ZREF_FACE2(3) = ZG_T(K1)
                   NORMAL_FACE2(3,2) = -1.d0
                   
                   
                                !South face 
                   XREF_FACE2(4) = XG_E(I1)
                   YREF_FACE2(4) = YG_N(J1) - DY(J1)
                   ZREF_FACE2(4) = ZG_T(K1)
                   NORMAL_FACE2(4,2) = 1.d0
                   
                                !Top face 
                   XREF_FACE2(5) = XG_E(I1)
                   YREF_FACE2(5) = YG_N(J1) 
                   ZREF_FACE2(5) = ZG_T(K1)
                   NORMAL_FACE2(5,3) = -1.d0
                   
                                !Bottom face 
                   XREF_FACE2(6) = XG_E(I1)
                   YREF_FACE2(6) = YG_N(J1) 
                   ZREF_FACE2(6) = ZG_T(K1) - DZ(K1)
                   NORMAL_FACE2(6,3) = 1.d0
                   
                   
                   Diagonal2 = (DX(I1)**2 + DY(J1)**2)
                   IF(DIMN.eq.3)  Diagonal2 = diagonal2 + DZ(K1)**2
                   Diagonal2 = dsqrt(diagonal2)
                   TOL_LINEINT2 = - DIAGONAL2*TOL_LINEPLANE_DIST_DES  
          
                   OUT_OF_PLANE1 = .FALSE.
                   OUT_OF_PLANE2 = .FALSE.
                   
                   PTXYZ_ON_LINE(:) = CNOT_ARR_TEMP(COUNT_CUT_FACE_LINE, :) 
                      
                   FACELOOP:  DO COUNT3 = 1, COUNT_FACE_LOOP
                      
                      !first check if the line and plane do not intersect 
                      TMP_DOT = DOT_PRODUCT(VEC_ARR_TEMP(COUNT_CUT_FACE_LINE,1:3), NORMAL_FACE(COUNT3,1:3))
                      
                      IF(ABS(TMP_DOT).GT.0.00001)  CYCLE !if the dot product is finite, that implies the line and the 
                      !plane intersect and computing normal distance is useless 
                      
                      DISTFROMLINE = NORMAL_FACE(COUNT3,1) * (PTXYZ_ON_LINE(1) - XREF_FACE(COUNT3)) + & 
                           &  NORMAL_FACE(COUNT3,2) * (PTXYZ_ON_LINE(2) - YREF_FACE(COUNT3)) + &
                           &  NORMAL_FACE(COUNT3,3) * (PTXYZ_ON_LINE(3) - ZREF_FACE(COUNT3)) 
                      
                      IF(DISTFROMLINE.LT.TOL_LINEINT)  THEN 
                         OUT_OF_PLANE1 = .TRUE.
                         EXIT !the faceloop
                      ENDIF
                      
                   end DO FACELOOP
                   
                   
                   
                   FACELOOP2:  DO COUNT3 = 1, COUNT_FACE_LOOP
                      !first check if the line and plane do not intersect 
                      TMP_DOT = DOT_PRODUCT(VEC_ARR_TEMP(COUNT_CUT_FACE_LINE,1:3), NORMAL_FACE2(COUNT3,1:3))
                      
                      IF(ABS(TMP_DOT).GT.0.00001)  CYCLE !if the dot product is finite, that implies the line and the 
                      !plane intersect and computing normal distance is useless 
                      
                      DISTFROMLINE2 = NORMAL_FACE2(COUNT3,1) * (PTXYZ_ON_LINE(1) - XREF_FACE2(COUNT3)) + & 
                           &  NORMAL_FACE2(COUNT3,2) * (PTXYZ_ON_LINE(2) - YREF_FACE2(COUNT3)) + &
                           &  NORMAL_FACE2(COUNT3,3) * (PTXYZ_ON_LINE(3) - ZREF_FACE2(COUNT3)) 
                      
                      IF(DISTFROMLINE2.LT.TOL_LINEINT2)  THEN 
                         OUT_OF_PLANE2 = .TRUE.
                         EXIT !the faceloop2
                      ENDIF
                      
                   end DO FACELOOP2
                   IF(OUT_OF_PLANE1.AND.OUT_OF_PLANE2) THEN 
                      WRITE(*,*) 'INTERSECTION POINT OUTSIDE THE CELL: CUT_CELL ? ', CUT_CELL_AT(CELL_ID)
                      WRITE(*,*) 'FACE NUMBER  = ', COUNT3, TOL_LINEINT, TOL_LINEINT2
                      
                      WRITE(*,'(A,3(2x,g17.8))') 'FACE COORDINATES =', XREF_FACE(COUNT3), yref_face(count3), zref_face(count3) 
                      
                      WRITE(*,'(A,3(2x,g17.8))') 'POINT ON LINE = ', PTXYZ_ON_LINE(1:3)
                      WRITE(*,'(A,6(2x,g17.8))') 'DISTANCE =',  DISTFROMLINE,  DISTFROMLINE2
                      
                      WRITE(*,'(A,3(2x,g17.8))') 'EAST NORTH AND TOP OF THE CELL 1=', XG_E(I), YG_N(J), ZG_T(K)
                      WRITE(*,'(A,3(2x,g17.8))') 'WEST SOUTH AND BOTTOM OF THE CELL 1=', XG_E(I)-DX(I), YG_N(J)-DY(J), ZG_T(K) - DZ(K)
                      
                      WRITE(*,'(A,3(2x,g17.8))') 'EAST NORTH AND TOP OF THE CELL 2=', XG_E(I1), YG_N(J1), ZG_T(K1)
                      WRITE(*,'(A,3(2x,g17.8))') 'WEST SOUTH AND BOTTOM OF THE CELL 2=', XG_E(I1)-DX(I1), YG_N(J1)-DY(J1), ZG_T(K1) - DZ(K1)
                      
                      !READ(*,*)
                      ADD_CUT_FACE_LINE = .FALSE.
                      !Exit the face loop. no need to test any further 
                      EXIT 
                   end IF
                   
                   !if the add_cut_face_line is still true, then check against the newly created cut_face_line bc's to prevent redundancy

                   IF(ADD_CUT_FACE_LINE) THEN 
                      !NOW CHECK IF A SIMILAR HAS NOT NOT ALREADY BEEN ADDED 
                         
                      DO COUNTBCIN = 1, COUNT_CUT_FACE_LINE - 1 !dont check against the one that was just added 
                         IJK_TEMP3 = IJK_ARR_TEMP2(COUNTBCIN)
                         WALL_TYPE_TEMP3 = WALL_TYPE_ARR_TEMP2(COUNTBCIN)
                         CUTFACELINE: IF(TRIM(WALL_TYPE_TEMP3).EQ.'CUT_FACE_LINE') THEN 
                            NORMAL_OLD(:) =  NORMAL_ARR_TEMP2(COUNTBCIN,:)
                            NORM_SQ = DOT_PRODUCT(NORM_AVG(:), NORMAL_OLD(:))
                            
                            IF(NORM_SQ.GT.0.99d0) THEN 
                               WRITE(*,*) 'WAS GOING TO ADD A NEW CUT_FACE_LINE TO NON CUT-CELL, BUT FOUND A WITH THE SAME NORMAL'
                               WRITE(*,*) 'BC ID and TOTAL BCs =', COUNTBCIN, COUNT_CUT_FACE_LINE
                               WRITE(*,'(A,3(2x,g17.8))') 'VEC = ',  VEC_ARR_TEMP(COUNTBCIN, :)
                               WRITE(*,'(A,3(2x,g17.8))') 'CNOT = ', CNOT_ARR_TEMP(COUNTBCIN, :)
                               WRITE(*,'(A,3(2x,g17.8))') 'NORMAL OLD = ', NORMAL_OLD(:)
                               WRITE(*,'(A,3(2x,g17.8))') 'NORMAL NEW = ', NORM_AVG(:) 
                                  
                               WRITE(*,'(A,3(2x,i4))') 'I, J, K NAT = ',  I_OF(CELL_ID), J_OF(CELL_ID), K_OF(CELL_ID)               
                               WRITE(*,'(A,3(2x,i4))') 'I, J, K OLD = ', I_OF(IJK_TEMP3), J_OF(IJK_TEMP3), K_OF(IJK_TEMP3)
                               WRITE(*,'(A,3(2x,i4))') 'I, J, K NEW = ', I_OF(IJK_TEMP), J_OF(IJK_TEMP), K_OF(IJK_TEMP)               
                               READ(*,*)
                               ADD_CUT_FACE_LINE = .TRUE.
                            END IF
                         END IF CUTFACELINE
                      END DO
                   END IF
                      
                   IF(.NOT.ADD_CUT_FACE_LINE) COUNT_CUT_FACE_LINE = COUNT_CUT_FACE_LINE - 1
                      
                      
                   WRITE(*,'(80(A))') '***********************************************************************'
                   IF(COUNT_CUT_FACE_LINE.LT.0) THEN 
                      WRITE(*,*) 'COUNT_CUT_FACE_LINE = ', COUNT_CUT_FACE_LINE
                      WRITE(*,'(A35, 3(2x,i4))') 'I, J, K of first plane', I_OF(IJK_TEMP), J_OF(IJK_TEMP), K_OF(IJK_TEMP)
                      WRITE(*,'(A35, 3(2x,i4))') 'I, J, K of second plane', I_OF(IJK_TEMP2), J_OF(IJK_TEMP2), K_OF(IJK_TEMP2)
                      read(*,*) 
                   ENDIF
                   
                   !exit out of the vertloop. 
                   EXIT
                   
                ENDIF
                
                
             ENDDO VERTLOOP

             IF(.NOT.ATLEAST_ONEINSIDE) THEN

                WRITE(*,*)'NOT A SINGLE INTERSECTION POINT IN THE CELL', ATLEAST_ONEINSIDE
                WRITE(*,*) 'RECHECKING FOR THIS CELL. THIS TIME WILL BE CHECKING THE DISTANCE TO THE NEW CUT-FACE SIMULTANEOUSLY FROM THE VERTEX AND THE INTERSECTION POINT'
                ATLEAST_ONEINSIDE = .TRUE.
                FORCE_INTPOINT_INSIDE = .FALSE. 

                WRITE(*,'(A,3(2x,i4))') 'CELL ID, CELL ID OLD, WALL CELL ID =', CELL_ID, IJK_TEMP, IJK_TEMP2

                WRITE(*,'(A,3(2x,i4))') 'I, J, K NAT = ', I, J, K
                WRITE(*,'(A,3(2x,i4))') 'I, J, K OLD = ', I_OF(IJK_TEMP), J_OF(IJK_TEMP), K_OF(IJK_TEMP)
                WRITE(*,'(A,3(2x,i4))') 'I, J, K NEW = ', I_OF(IJK_TEMP2), J_OF(IJK_TEMP2), K_OF(IJK_TEMP2)

                !                   READ(*,*)
                GOTO 600
             ENDIF
             
             
          ENDDO
       ENDDO

       DES_CELLWISE_BCDATA(CELL_ID)%COUNT_DES_BC = COUNT_BC - DEL_COUNT + COUNT_CUT_FACE_LINE
       COUNT_BC_NEW = 0

       DO COUNT = 1, COUNT_BC
          IF(.NOT.BC_DEL_ARR(COUNT)) THEN 
             COUNT_BC_NEW = COUNT_BC_NEW+1
             SELECT CASE (TRIM(WALL_TYPE_ARR_TEMP(COUNT)))
             CASE('NORMAL_WALL')

                DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT_BC_NEW)%IJK_SCAL = IJK_ARR_TEMP(COUNT)
                
                ALLOCATE(DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT_BC_NEW)%NORMAL(DIMN))
                
                DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT_BC_NEW)%NORMAL(:) = NORMAL_ARR_TEMP(COUNT,:)
                
                DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT_BC_NEW)%DES_BC_TYPE = WALL_TYPE_ARR_TEMP(COUNT)
       
             CASE('CUT_FACE')
          
                DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT_BC_NEW)%IJK_SCAL = IJK_ARR_TEMP(COUNT)
                DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT_BC_NEW)%DES_BC_TYPE = WALL_TYPE_ARR_TEMP(COUNT)
                
                !no need to store the normal, as depening on the particle position, the 
                !new normal from cut-face to particle center will be calculated by 
                !get_del_h_des
             CASE DEFAULT
                WRITE(*,*)' EROR IN SUBROUTINE CHECK_NON_CUTCELLBC'
                WRITE(*,*)'UNKNOWN TYPE OF DES_WALL_TYPE:',TRIM(WALL_TYPE_ARR_TEMP(COUNT))
                WRITE(*,*)'ACCEPTABLE TYPES ARE:' 
                WRITE(*,*)'NORMAL_WALL' 
                WRITE(*,*)'CUT_FACE' 
             
                
                STOP
             END SELECT
          ENDIF
       ENDDO
       
       
       DO COUNT = 1, COUNT_CUT_FACE_LINE
          COUNT_BC_NEW = COUNT_BC_NEW+1
          SELECT CASE (TRIM(WALL_TYPE_ARR_TEMP2(COUNT)))
          CASE('CUT_FACE_LINE')
             DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT_BC_NEW)%IJK_SCAL = IJK_ARR_TEMP2(COUNT)!this will really not matter 
             ALLOCATE(DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT_BC_NEW)%NORMAL(DIMN))
                
             DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT_BC_NEW)%NORMAL(:) = NORMAL_ARR_TEMP2(COUNT,:)
                
             DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT_BC_NEW)%DES_BC_TYPE = WALL_TYPE_ARR_TEMP2(COUNT)
                
                !for this special case, allocate the line attributes in 3-d
             ALLOCATE(DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT_BC_NEW)%CNOT(3))
                
             ALLOCATE(DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT_BC_NEW)%VEC(3))
             DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT_BC_NEW)%CNOT(:) = CNOT_ARR_TEMP(COUNT,:)
             DES_CELLWISE_BCDATA(CELL_ID)%BDRY_LIST(COUNT_BC_NEW)%VEC(:) = VEC_ARR_TEMP(COUNT,:)
          CASE DEFAULT
             WRITE(*,*)' EROR IN SUBROUTINE CHECK_NON_CUTCELLBC'
             WRITE(*,*)'UNKNOWN TYPE OF DES_WALL_TYPE:',TRIM(WALL_TYPE_ARR_TEMP2(COUNT))
             WRITE(*,*)'ONLY LOOKING FOR CUT_FACE_TYPE HERE:' 
             STOP
          END SELECT
       ENDDO

       IF(COUNT_BC_NEW.NE.DES_CELLWISE_BCDATA(CELL_ID)%COUNT_DES_BC) THEN 
          WRITE(*,*) 'problem in check_non_cutcellbc'
          write(*,*) 'count_bc_new not equal to new number of bc', COUNT_BC_NEW, DES_CELLWISE_BCDATA(CELL_ID)%COUNT_DES_BC
          WRITE(*,*) 'TERMINAL ERROR'
          WRITE(*,*) 'STOPPING'
       ENDIF
       END SUBROUTINE CHECK_NON_CUTCELLBC

       SUBROUTINE CG_DEL_OUTOFDOMAIN_PARTS
       USE discretelement
       USE cutcell 
       USE compar 
       USE fldvar
       use mpi_utility
       USE sendrecv
       USE mfix_pic
       IMPLICIT NONE 
       
       INTEGER :: LCURPAR,LPIP_ALL(0:NUMPES-1),LGLOBAL_ID, LGHOST_CNT_ALL(0:NUMPES-1), LL, PIP_INIT, PIP_FINAL
       INTEGER :: IJK, K
       INCLUDE 'function.inc'
      
       IF(mype.eq.pe_IO) THEN 
          WRITE(*,*) 'DELETING PARTICLES OUTSIDE THE CUT-FACE BCs'
       ENDIF

       !if the particles in cell was called then the ghost particles
       !were also added to the particle araay on each processor. 
       !Although this ghost particles shud not be there for 
       !MPPIC, in the current implementation that is the case
       !First reset the pip on each processor to reflect only 
       !the particles present on that proc. 
       DO LL = 1, MAX_PIP
          IF(.NOT.PEA(LL,1)) CYCLE  
          IF(PEA(LL,4)) CYCLE  
          PIJK(LL, 1) = DES_GETINDEXFROMPOS(ISTART1,IEND1,DES_POS_NEW(LL,1),XE(1:size(XE,1)),'X','I')
          PIJK(LL, 2) = DES_GETINDEXFROMPOS(JSTART1,JEND1,DES_POS_NEW(LL,2),YN(1:size(YN,1)),'Y','J')
          PIJK(LL, 3) = 1
          IF(DIMN.eq.3) PIJK(LL, 3) = DES_GETINDEXFROMPOS(KSTART1,KEND1,DES_POS_NEW(LL,3),ZT(1:size(ZT,1)),'Z','K')
          PIJK(LL, 4) = FUNIJK(PIJK(LL, 1), PIJK(LL, 2), PIJK(LL, 3))         
       ENDDO
      !**********DELETE START 

      !OPEN(UNIT=24,FILE='PARTICLES_BEFORE_DEL.dat',&
      !STATUS='REPLACE')
      !WRITE(24, *) '"X"', '"Y"', '"Z"', '"Dp"', '"I"', '"J"', '"K"', '"IJK"'
      !DO LL = 1, MAX_PIP
      !   IF(.NOT.PEA(LL,1)) CYCLE  
      !   IF(PEA(LL,4)) CYCLE  

      ! WRITE(24,'(4(2x, g17.8), 4(2x, i20))')&
      !   (DES_POS_NEW(LL,K),K=1,DIMN), 2.d0*DES_RADIUS(LL), PIJK(LL, 1:4)
      !ENDDO
      !CLOSE(24)
      !call mfix_exit(mype)
      
      !**********DELETE END 

       DO LL = PIP - IGHOST_CNT + 1, PIP, 1 
          !deactivate the particle and also reset its ghost 
          !particle identifier flag pea(ll,4) 
          PEA(LL, 1) = .false. 
          PEA(LL, 4) = .false. 
       ENDDO
       
       PIP = PIP - IGHOST_CNT 
       IGHOST_CNT  = 0 

!!$       LGHOST_CNT_ALL = 0 
!!$       LGHOST_CNT_ALL(MYPE) = IGHOST_CNT
!!$       
!!$       CALL GLOBAL_ALL_SUM(LGHOST_CNT_ALL)
!!$
       LPIP_ALL = 0 
       LPIP_ALL(MYPE) = PIP - IGHOST_CNT 
       CALL GLOBAL_ALL_SUM(LPIP_ALL)
       PIP_INIT = SUM(LPIP_ALL(:))
       
!!$
!!$       IF(mype.eq.pe_IO) THEN 
!!$          WRITE(*,*) 'PIP = ', LPIP_ALL(:)
!!$
!!$          WRITE(*,*) 'GHOST COUNT = ', LGHOST_CNT_All(0:numpes-1)
!!$          WRITE(*,*) 'PIP ADJUSTED = ', LPIP_ALL(:) - LGHOST_CNT_all(:)
!!$          WRITE(*,*) 'ACTIVE LAST = ', PEA(PIP, 1)
!!$          
!!$       end IF


       IF(MPPIC) THEN 
          CALL LIST_PARTS_TOBE_DEL_MPPIC
       ELSE 
          CALL LIST_PARTS_TOBE_DEL_DEM
       END IF

       !in the above call, particle arrays are also rebuilt
       !so that all the active particles are in a continuous manner
       
       
       !reset the global id's of particles after this deletion
       !and array compaction steps
       
       LPIP_ALL = 0 
       LPIP_ALL(MYPE) = PIP 
       CALL GLOBAL_ALL_SUM(LPIP_ALL)
       PIP_FINAL = SUM(LPIP_ALL(:))
       LGLOBAL_ID = SUM(LPIP_ALL(0:MYPE-1))
       DO LCURPAR  = 1,PIP 
          LGLOBAL_ID = LGLOBAL_ID + 1
          IGLOBAL_ID(LCURPAR) = LGLOBAL_ID 
       END DO

       

!!$       LGHOST_CNT_ALL = 0 
!!$       LGHOST_CNT_ALL(MYPE) = IGHOST_CNT
!!$       
!!$       CALL GLOBAL_ALL_SUM(LGHOST_CNT_ALL)
!!$
!!$       LPIP_ALL = 0 
!!$       LPIP_ALL(MYPE) = PIP 
!!$       CALL GLOBAL_ALL_SUM(LPIP_ALL)
!!$
!!$       IF(mype.eq.pe_IO) THEN 
!!$          WRITE(*,*) 'AFTER PARTILCES OUTSIDE THE CUT-FACE BCs'
!!$          WRITE(*,*) 'PIP = ', LPIP_ALL(:)
!!$
!!$          WRITE(*,*) 'GHOST COUNT = ', LGHOST_CNT_All(0:numpes-1)
!!$          WRITE(*,*) 'PIP ADJUSTED = ', LPIP_ALL(:) - LGHOST_CNT_all(:)
!!$          WRITE(*,*) 'ACTIVE LAST = ', PEA(PIP, 1)
!!$          
!!$       end IF

       if(myPE.EQ.pe_IO) THEN 
          WRITE(*,'(A, i10)') 'INITIAL TOTAL NUMBER OF PARTICLES IN THE SYSTEM = ', PIP_INIT
          WRITE(*,'(A, i10)') 'NUMBER OF PARTICLES DELETED = ', PIP_INIT -  PIP_FINAL
          WRITE(*,'(A, i10)') 'FINAL TOTAL NUMBER OF PARTICLES IN THE SYSTEM = ', PIP_FINAL
       ENDIF
          
       END SUBROUTINE CG_DEL_OUTOFDOMAIN_PARTS

       SUBROUTINE LIST_PARTS_TOBE_DEL_MPPIC
       USE param1
       USE funits
       USE run
       USE compar      
       USE discretelement
       USE cutcell
       USE indices
       USE physprop
       USE parallel
       USE geometry
      
       IMPLICIT NONE 
       
       INTEGER :: I, J, K, IJK, LL, COUNT_BC, COUNT, IJK_WALL, PC, IDIM
       DOUBLE PRECISION :: NORMAL(DIMN), WALL_COOR(DIMN), NORM1, NORM2, NORM3, DIST(DIMN), DISTMOD, NORM_LINE_TO_PART(DIMN)
         
       DOUBLE PRECISION TEMPX, TEMPY, TEMPZ
       DOUBLE PRECISION :: NORMAL_MAG, NORM_DOT, TEMP_DOT
       CHARACTER*100 :: WALL_TYPE 
       INTEGER :: COUNT_DELETE, PIP_OLD, TMP2, TMP1
       LOGICAL :: DELETE_PART
       CHARACTER*100 :: FILENAME
       DOUBLE PRECISION X1MINX0(3), X2MINX1(3), TEMP_CROSSP(3), X1(3)
       INCLUDE 'function.inc'
       
       IF(mype.eq.pe_IO) WRITE(*,*) 'IN LIST_PARTS_TOBE_DEL_MPPIC: FLAGGING PARTICLES OUTSIDE THE DOMAIN'

       COUNT_DELETE  = 0
       PIP_OLD = PIP

       ALLOCATE(TOBE_DELETED(MAX_PIP))
       TOBE_DELETED(:) = .false. 
       
       PC = 1      
       PARTLOOP: DO LL = 1, MAX_PIP

          IF(PC .GT. PIP) EXIT PARTLOOP
          PC = PC + 1 

          IJK = PIJK(LL,4)
          
          IF(PIJK(LL,4).EQ.0) THEN 
             WRITE(*,*) 'PIJK = ', PIJK(LL,4), LL, PIP, PC
             READ(*,*) 
          ENDIF
          COUNT_BC = DES_CELLWISE_BCDATA(IJK)%COUNT_DES_BC
          FLUID: IF(FLUID_AT(IJK)) THEN
             !Check if any of the particles are overlapping with the walls 
             !this will include the normal ghost cell walls and also the cut cells 
             DELETE_PART = .FALSE.
         
             TMP1 = COUNT_DELETE
             IF(.NOT.CUT_CELL_AT(IJK)) CYCLE PARTLOOP

                  
             TEMPX = DES_POS_NEW(LL,1)
             TEMPY = DES_POS_NEW(LL,2)
             TEMPZ = 0.d0 
             IF(DIMN.EQ.3) TEMPZ = DES_POS_NEW(LL,3)
                     
             CALL GET_DEL_H_DES(IJK,'SCALAR',TEMPX, TEMPY, TEMPZ, &
                  &  DISTMOD, NORM1, NORM2, NORM3, .TRUE.)
             
             IF(DISTMOD.LE.ZERO) DELETE_PART = .TRUE.
      
             IF(DELETE_PART) THEN 
                TOBE_DELETED(LL) = .TRUE.
                COUNT_DELETE = COUNT_DELETE + 1
             END IF
             TMP2 = COUNT_DELETE
             IF(TMP2-TMP1.GT.1) WRITE(*,*) 'COUNT_DELETE JUMPED BY MORE THAN ONE = ', TMP2, TMP1
          ELSE 
             TOBE_DELETED(LL) = .TRUE.
             COUNT_DELETE = COUNT_DELETE + 1
               
          END IF FLUID
       END DO PARTLOOP
         
         
       !Now remove the particles just marked as tobe deleted 
       CALL REMOVE_PARTICLES(COUNT_DELETE)

       !Call particles in cell again to get things back in right shape

       IF(mype.eq.pe_IO) WRITE(*,*) 'END IN LIST_PARTICLES_TOBE_DELETED'
         
       DEALLOCATE(TOBE_DELETED)
         
       IF(PIP_OLD.NE.(PIP+COUNT_DELETE)) THEN 
          IF(DMP_LOG) THEN 
            WRITE(UNIT_LOG,*) 'LIST_PARTS_TOBE_DEL_MPPIC'
            WRITE(UNIT_LOG,*) 'PIP_OLD NE PIP_NEW +  DELETED', PIP_OLD, PIP, COUNT_DELETE
            WRITE(UNIT_LOG,*) 'LOGIC IN LIST_PARTILCES_TOBE_DELETED SCREWED UP FOR THE ABOVE MESSAGE TO SHOW UP'
            
            WRITE(UNIT_LOG,*) 'TERMINAL ERROR: STOPPING'
         ENDIF
         CALL mfix_exit(myPE)

      END IF

      
      WRITE(FILENAME,'(A,"_",I5.5,".DAT")') 'PARTS_AFTERDEL',MYPE
      OPEN(UNIT=24,FILE=TRIM(FILENAME),&
           STATUS='UNKNOWN')
      
      PC = 1
      
      DO LL = 1, MAX_PIP
         IF(PC .GT. PIP) EXIT
         WRITE(24,'(10(X,ES15.5))')&
              (DES_POS_NEW(LL,K),K=1,DIMN), DES_RADIUS(LL),&
              RO_Sol(LL), (DES_VEL_NEW(LL,K),K=1,DIMN) 
         PC = PC+1
      ENDDO
      CLOSE(24, STATUS='keep')
         
    END SUBROUTINE LIST_PARTS_TOBE_DEL_MPPIC

    SUBROUTINE LIST_PARTS_TOBE_DEL_DEM
       USE param1
       USE funits
       USE run
       USE compar      
       USE discretelement
       USE cutcell
 
       USE indices
       USE physprop
       USE parallel
       USE geometry
       
       IMPLICIT NONE 
       
       INTEGER :: I, J, K, IJK, LL, COUNT_BC, COUNT, IJK_WALL, PC, IDIM
       DOUBLE PRECISION :: NORMAL(DIMN), WALL_COOR(DIMN), NORM1, NORM2, NORM3, DIST(DIMN), DISTMOD, NORM_LINE_TO_PART(DIMN)
         
       DOUBLE PRECISION TEMPX, TEMPY, TEMPZ
       DOUBLE PRECISION :: NORMAL_MAG, NORM_DOT, TEMP_DOT
       CHARACTER*100 :: WALL_TYPE 
       INTEGER :: COUNT_DELETE, PIP_OLD, TMP2, TMP1
       LOGICAL :: DELETE_PART, PRINT_PART
       
       DOUBLE PRECISION X1MINX0(3), X2MINX1(3), TEMP_CROSSP(3), X1(3)
       CHARACTER*100 :: FILENAME
       INCLUDE 'function.inc'
       
       WRITE(*,*) 'IN LIST_PARTICLES_TOBE_DELETED: FLAGGING PARTICLES OUTSIDE THE DOMAIN'

       COUNT_DELETE  = 0
       PIP_OLD = PIP
       
       ALLOCATE(TOBE_DELETED(MAX_PIP))
       
       TOBE_DELETED(:) = .false. 

      PC = 1      
      DO LL = 1, MAX_PIP

         IF(PC .GT. PIP) EXIT
         
         IJK = PIJK(LL,4)
         COUNT_BC = DES_CELLWISE_BCDATA(IJK)%COUNT_DES_BC
         I = I_OF(IJK)
         J = J_OF(IJK)
         K = K_OF(IJK)
         print_part = .false.
         IF(K.eq.kend1) print_part = .true.
         !IF(PRINT_PART) THEN 
            
         !   WRITE(*,'(A,3(2x,i4))') 'I, J, K NAT = ', I, J, K
         !   DO COUNT = 1, COUNT_BC 
         !      IJK_WALL = DES_CELLWISE_BCDATA(IJK)%BDRY_LIST(COUNT)%IJK_SCAL
         !      WALL_TYPE = DES_CELLWISE_BCDATA(IJK)%BDRY_LIST(COUNT)%DES_BC_TYPE
               
         !      WRITE(*,*) 'WALL TYPE = ', TRIM(WALL_TYPE)
         !      WRITE(*,'(A,3(2x,i4))') 'I, J, K WALL = ', I_OF(IJK_WALL), J_OF(IJK_WALL), K_OF(IJK_WALL)                      
         !   ENDDO
         !   READ(*,*) 
         !ENDIF
         FLUID: IF(FLUID_AT(IJK)) THEN
            !Check if any of the particles are overlapping with the walls 
            !this will include the normal ghost cell walls and also the cut cells 
         DELETE_PART = .FALSE.

         TMP1 = COUNT_DELETE
            CHECK_BC: DO COUNT = 1, COUNT_BC 
               IJK_WALL = DES_CELLWISE_BCDATA(IJK)%BDRY_LIST(COUNT)%IJK_SCAL
               WALL_TYPE = DES_CELLWISE_BCDATA(IJK)%BDRY_LIST(COUNT)%DES_BC_TYPE
               
               IF(TRIM(WALL_TYPE).EQ.'NORMAL_WALL'.or.TRIM(WALL_TYPE).EQ.'OUTFLOW') WALL_TYPE = 'NORMAL_WALL'

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
                  DISTMOD = DISTMOD - DES_RADIUS(LL) 
                  
                  IF(DISTMOD.LT.ZERO) THEN 
                     DELETE_PART = .TRUE.
!                     TOBE_DELETED(LL) = .TRUE.
!                     COUNT_DELETE = COUNT_DELETE + 1
                  END IF
                  
               CASE('CUT_FACE')
                  TEMPX = DES_POS_NEW(LL,1)
                  TEMPY = DES_POS_NEW(LL,2)
                  TEMPZ = 0.d0 
                  IF(DIMN.EQ.3) TEMPZ = DES_POS_NEW(LL,3)
                     
                  CALL GET_DEL_H_DES(IJK_WALL,'SCALAR',TEMPX, TEMPY, TEMPZ, &
                       &  DISTMOD, NORM1, NORM2, NORM3, .TRUE.)

                  IF(DISTMOD.LE.ZERO) DELETE_PART = .TRUE.

                  DISTMOD = DISTMOD - DES_RADIUS(LL)
                  
                  IF(DISTMOD.LT.ZERO) DELETE_PART = .TRUE.

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
                  NORMAL(1:DIMN) = DES_CELLWISE_BCDATA(IJK)%BDRY_LIST(COUNT)%NORMAL(1:DIMN)
                  !normal for cut-faces or cut-face-lines are from line/face into the fluid 

                  X1MINX0(1) = X1(1) - TEMPX
                  X1MINX0(2) = X1(2) - TEMPY
                  X1MINX0(3) = ZERO 
                  IF(DIMN.EQ.3)  X1MINX0(3) = X1(3) - TEMPZ

                  CALL DES_CROSSPRDCT_3D(TEMP_CROSSP(1:3), X2MINX1(1:3), X1MINX0(1:3))
                  
                  DISTMOD = SQRT(DOT_PRODUCT(TEMP_CROSSP(1:3), TEMP_CROSSP(1:3)))
                  !need to find the normal from des_pos_new to the point where the perpendicular 
                  !line from des_pos_new meets the cut_face_line
                  !this normal, pointing from cut-face-line toward the particle, can be derived as 
                  !norm = (des_pos - cnot) - {(des_pos - cnot) \cdot vec } vec 
                  !where cnot and vec are from the definition of line 
                  TEMP_DOT = (TEMPX-X1(1))*X2MINX1(1) + (TEMPY-X1(2))*X2MINX1(2)
                  
                  IF(DIMN.EQ.3) TEMP_DOT = TEMP_DOT + (TEMPZ-X1(3))*X2MINX1(3)
                  
                  IF(ABS(TEMP_DOT).GT.ZERO.AND.DIMN.EQ.2) THEN 
                     WRITE(*,*) 'TEMP_DOT = ', TEMP_DOT 
                     READ(*,*)
                  end IF
                  
                  NORM_LINE_TO_PART(1:DIMN) = DES_POS_NEW(LL,1:DIMN)-X1(1:DIMN) - TEMP_DOT*X2MINX1(1:DIMN)
                  
                  NORM_DOT = DOT_PRODUCT(NORMAL(1:DIMN), NORM_LINE_TO_PART(1:DIMN))
                  
                  IF(NORM_DOT.LT.ZERO) DELETE_PART = .TRUE.

                  DISTMOD = DISTMOD - DES_RADIUS(LL)
                  
                  IF(DISTMOD.LT.ZERO) DELETE_PART = .TRUE.
               CASE DEFAULT
                  TEMPX = DES_POS_NEW(LL,1)
                  TEMPY = DES_POS_NEW(LL,2)
                  TEMPZ = 0.d0 
                  IF(DIMN.EQ.3) TEMPZ = DES_POS_NEW(LL,3)

                  WRITE(*, 1001) mype, TRIM(WALL_TYPE), IJK, I_OF(IJK), & 
                  & J_OF(IJK), K_OF(IJK), IJK_WALL, I_OF(IJK_WALL), J_OF(IJK_WALL), K_OF(IJK_WALL),  &
                  & TEMPX, TEMPY, TEMPZ

                  if(DMP_LOG) WRITE(UNIT_LOG, 1001) mype, TRIM(WALL_TYPE), IJK, I_OF(IJK), & 
                  & J_OF(IJK), K_OF(IJK), IJK_WALL, I_OF(IJK_WALL), J_OF(IJK_WALL), K_OF(IJK_WALL),  &
                  & TEMPX, TEMPY, TEMPZ

                  
 1001             FORMAT( 10x,'FROM PROC = ', 2x, i10, /10x, &
                  & 'EROR IN LIST_PARTILCES_TOBE_DELETED:' , /10x, &
                  & 'UNKNOWN TYPE OF DES_WALL_TYPE:', 2x, A, /10x, &
                  & 'PARTICLE CELL ATTRIBUTES (IJK, I, J, K)', 4(2x,i10) , /10x, &
                  & 'WALL CELL ATTRIBUTES     (IJK, I, J, K)', 4(2x,i10) , /10x, &
                  & 'PARTICLE POSITION CORDINATES', 3(2x,g17.8) , /10x, &
                  & 'ACCEPTABLE TYPES ARE:' , /10x, &
                  & 'NORMAL_WALL' , /10x, &
                  & 'CUT_FACE' , /10x, &
                  & 'CUT_FACE_LINE' , /10x, &                  
                  & 'TERMINAL ERROR: STOPPING' , /)
      
                  CALL mfix_exit(mype) 
                  END SELECT
                  
               END DO CHECK_BC
               
               IF(DELETE_PART) THEN 
                  TOBE_DELETED(LL) = .TRUE.
                  COUNT_DELETE = COUNT_DELETE + 1
               END IF
               TMP2 = COUNT_DELETE
               IF(TMP2-TMP1.GT.1) WRITE(*,*) 'COUNT_DELETE JUMPED BY MORE THAN ONE = ', TMP2, TMP1
            ELSE 
               TOBE_DELETED(LL) = .TRUE.
               COUNT_DELETE = COUNT_DELETE + 1
               
            END IF FLUID
            PC = PC+1
         END DO
         
         
         !Now remove the particles just marked as tobe deleted 
         CALL REMOVE_PARTICLES(COUNT_DELETE)

         WRITE(*,*) 'END IN LIST_PARTICLES_TOBE_DELETED'
         
         DEALLOCATE(TOBE_DELETED)
         
         IF(PIP_OLD.NE.(PIP+COUNT_DELETE)) THEN 
            WRITE(*,*) 'LIST_PARTILCES_TOBE_DELETED'
            WRITE(*,*) 'PIP_OLD NE PIP_NEW +  DELETED', PIP_OLD, PIP, COUNT_DELETE
            WRITE(*,*) 'LOGIC IN LIST_PARTILCES_TOBE_DELETED SCREWED UP FOR THE ABOVE MESSAGE TO SHOW UP'
            
          WRITE(*,*) 'TERMINAL ERROR: STOPPING'
          STOP

         END IF
         
         WRITE(FILENAME,'(A,"_",I5.5,".DAT")') 'PARTS_AFTERDEL',MYPE
         OPEN(UNIT=24,FILE=TRIM(FILENAME),&
           STATUS='UNKNOWN')

         PC = 1
         
         DO LL = 1, MAX_PIP
            IF(PC .GT. PIP) EXIT
            WRITE(24,'(10(X,ES15.5))')&
               (DES_POS_NEW(LL,K),K=1,DIMN), DES_RADIUS(LL),&
               RO_Sol(LL), (DES_VEL_NEW(LL,K),K=1,DIMN) 
            PC = PC+1
         ENDDO
         CLOSE(24)
         
       END SUBROUTINE LIST_PARTS_TOBE_DEL_DEM 
       
       
       SUBROUTINE REMOVE_PARTICLES(NUMBER_DELETIONS)
       USE DISCRETELEMENT
       USE funits 
       USE compar
       USE mfix_pic
       IMPLICIT NONE

       INTEGER, INTENT(IN) :: NUMBER_DELETIONS
       
       !LOCAL VARIABLES

       INTEGER :: PC, LL, COUNT_ACTIVE, COUNT_ACTIVE2, II, COUNT_DELETE
       LOGICAL :: FOUND_REPL_PART
       
       WRITE(*,*) 'IN REMOVE_PARTICLES: REMOVING FLAGGED PARTICLES'
                
       PC = 1
       COUNT_ACTIVE  = 0
       COUNT_DELETE = 0
       DO LL = 1, MAX_PIP
          
          IF(PC .GT. PIP) EXIT

          ACTI: IF(PEA(LL,1)) THEN  !if active 
             DELETE: IF(TOBE_DELETED(LL)) THEN 
                !Traverse a reverse loop from the bottom and find a 
                !particle to swap this particle 
                FOUND_REPL_PART = .FALSE. 
                DO II = MAX_PIP, LL+1, -1 
                   IF(FOUND_REPL_PART) EXIT 
                   !The first particle that is active and not marked as 
                   !to be deleted will be chosen to replace LL
                   IF(PEA(II,1).and..NOT.TOBE_DELETED(II)) THEN
                      DES_POS_NEW(LL, :) = DES_POS_NEW(II,:)
                      DES_POS_OLD(LL, :) = DES_POS_OLD(II,:)
                      
                      DES_VEL_NEW(LL, :) = DES_VEL_NEW(II,:)
                      DES_VEL_OLD(LL, :) = DES_VEL_OLD(II,:)

                      OMEGA_NEW(LL, :) = OMEGA_NEW(II, :)
                      OMEGA_OLD(LL, :) = OMEGA_OLD(II, :)
                                            
                      DES_RADIUS(LL) = DES_RADIUS(II)
                      RO_Sol(LL) = RO_Sol(II)
                      PVOL(LL) = PVOL(II)
                      PMASS(LL) = PMASS(II)
                      OMOI(LL) = OMOI(II)

                      PIJK(LL,:) = PIJK(II,:)

                      IF(MPPIC) DES_STAT_WT(LL) = DES_STAT_WT(II)
                      !mark II as inactive 
                      PEA(II,1) = .FALSE.
                      
                      FOUND_REPL_PART = .TRUE.
                      
                      COUNT_DELETE = COUNT_DELETE + 1
                   END IF
                END DO
                
                IF(.NOT.FOUND_REPL_PART) THEN 
                   
                   PEA(LL,1) = .FALSE.
                   
                   COUNT_DELETE = COUNT_DELETE + 1

!                   WRITE(*,*) 'MESSAGE FROM REMOVE_PARTICLES'
!                   WRITE(*,*) 'REPLACEMENT PARTICLE NOT FOUND FOR PARTICLE', LL
!                   WRITE(*,*) 'PIP AND MAX_PIP = ', PIP, MAX_PIP
!                   WRITE(*,*) 'DELETED SO FAR AND TO DELETE', COUNT_DELETE, NUMBER_DELETIONS
                   
!                   WRITE(*,*) 'THIS CAN HAPPEN IF LL IS THE LAST PARTICLE'
!                   WRITE(*,*) 'OR IF THE PARTICLES FOLLOWING LL ARE ALSO MARKED AS TO BE DELETED'

                   
!                   WRITE(*,*) 'NOT AN ERROR'
                   
!                   WRITE(*,*) 'END OF MESSAGE FROM REMOVE_PARTICLES'
                   !in this case set the LL particle to inactive 
                   
                ENDIF
                   
             ENDIF DELETE
             
          END IF ACTI
          
          !if the particle is active then increment the count_active 
          !remember active can be false also. See the above logic for found_repl_part 
          IF(PEA(LL,1)) COUNT_ACTIVE = COUNT_ACTIVE + 1
             
          PC = PC + 1
       ENDDO                    ! end loop over paticles LL
      
       !PIP now equals COUNT_ACTIVE 
       PIP = COUNT_ACTIVE 
       
       IF(DMP_LOG) THEN 
          WRITE(UNIT_LOG,*) 'TOTAL NUMBER OF PARTICLES DELETED = ', COUNT_DELETE
          WRITE(UNIT_LOG,*) 'NUMBER OF PARTICLES LEFT = ', PIP
       ENDIF
       
       !Sanity check 
       COUNT_ACTIVE2 = 0 
       DO LL = 1, MAX_PIP
          IF(PEA(LL,1)) COUNT_ACTIVE2 = COUNT_ACTIVE2 + 1
       ENDDO

       IF(COUNT_ACTIVE2.NE.COUNT_ACTIVE) THEN
          IF(DMP_LOG) THEN 
             WRITE(UNIT_LOG,*) 'ERROR MESSAGE FROM REMOVE_PARTICLES'
             WRITE(UNIT_LOG,*) 'COUNT_ACTIVE2 NE COUNT_ACTIVE',COUNT_ACTIVE2, COUNT_ACTIVE
             WRITE(UNIT_LOG,*) 'LOGIC IN REMOVE PARTICLES SCREWED UP FOR THE ABOVE MESSAGE TO SHOW UP'
          end IF
          CALL mfix_exit(myPE)
       ENDIF
       
       
       if(mype.eq.pe_IO) WRITE(*,*) 'END REMOVE_PARTICLES'
      END SUBROUTINE REMOVE_PARTICLES
     

       
!       SUBROUTINE  GET_DEL_H_DES
!       END SUBROUTINE  GET_DEL_H_DES
       
