      MODULE mppic_wallbc
      USE param 
      USE param1 
      USE discretelement   
      USE bc 
      USE geometry
      USE compar 
      USE indices 
      USE funits
      USE mpi_utility
      USE constant 
      USE physprop
      USE randomno
      USE cutcell
      USE fldvar 
      USE mfix_pic
      IMPLICIT NONE
      PRIVATE 
      
      PUBLIC:: mppic_apply_wallbc, MPPIC_FIND_NEW_CELL
      
      LOGICAL :: INSIDE_DOMAIN, INSIDE_SMALL_CELL, REFLECT_FROM_ORIG_CELL
      INTEGER :: PIP_DEL_COUNT, PIP_ADD_COUNT, REFLECTING_CELL 
      
      INTEGER :: PIJK_OLD(5), PIJK_INT(5), REFLECT_COUNT 
      DOUBLE PRECISION :: DIST_OFFSET 
      CONTAINS 

      SUBROUTINE MPPIC_APPLY_WALLBC
      
      INTEGER :: I, J, K, IJK_CELL, IDIM, NPG, IDIR, IP, NP, IJK
      INTEGER :: COUNT, COUNT_BC , PC, IPROC
      

      INTEGER :: LPIP_DEL_COUNT_ALL(0:numPEs-1), LPIP_ADD_COUNT_ALL(0:numPEs-1)
      LOGICAL :: DEBUG_DES_LOCAL, tmp_logical, reflected 
      
      LOGICAL, SAVE :: FIRST_PASS = .TRUE.
      DOUBLE PRECISION :: POS_INIT(DIMN)
      
      DOUBLE PRECISION :: NORM_CF(3), XPOS, YPOS, ZPOS, DIST
      
      INCLUDE 'function.inc'
      IF(FIRST_PASS) THEN 
         open(1000, file='parts_out_mppic.dat', form="formatted", status="unknown")
         write(1000,*) 'partices outside domain'
         FIRST_PASS = .false. 
      ELSE
         open(1000, file='parts_out_mppic.dat', form="formatted", status="old", position="append")
      ENDIF
      !write(1000, *) 's_time =', s_time 

      PIP_DEL_COUNT = ZERO 
      DEBUG_DES_LOCAL = .FALSE.
      PC = 1
      PART_LOOP: DO NP = 1, MAX_PIP
      REFLECT_COUNT = 0 
! pradeep: skip ghost particles
         IF(PC.GT.PIP) EXIT
         IF(.NOT.PEA(NP,1)) CYCLE PART_LOOP
         PC = PC+1
         IF(PEA(NP,4)) CYCLE PART_LOOP
         IJK = PIJK(NP, 4) 

         COUNT_BC = DES_CELLWISE_BCDATA(IJK)%COUNT_DES_BC 

         IF(COUNT_BC.EQ.0) CYCLE PART_LOOP
         !this implies that this particle can't be in a cell 
         !that is close to boundary 
         POS_INIT(1:DIMN) = DES_POS_NEW(NP, 1:DIMN)
         I = PIJK(NP, 1) 
         J = PIJK(NP, 2) 
         K = PIJK(NP, 3) 
         
         PIJK_OLD(1:5) = PIJK(NP, 1:5) 

         CALL MPPIC_FIND_NEW_CELL(NP)
         
         !IF(SMALL_CELL_AT(PIJK_OLD(4)).or.SMALL_CELL_AT(PIJK(NP, 4))) WRITE(1000,*) NP, SMALL_CELL_AT(PIJK_OLD(4)), SMALL_CELL_AT(PIJK(NP, 4)), DES_CELLWISE_BCDATA(PIJK_OLD(4))%COUNT_DES_BC, DES_CELLWISE_BCDATA(PIJK(NP, 4))%COUNT_DES_BC  

         CALL MPPIC_CHECK_IF_INSIDE_DOMAIN(NP)
            
         tmp_logical = reflect_from_orig_cell
         !WRITE(*,*) 'NP = ', NP, INSIDE_DOMAIN, REFLECT_FROM_CUTFACE, REFLECT_FROM_ORIG_CELL, PIJK(NP, 1:3)
         !WRITE(*,*) 'NP = ', NP, DES_POS_NEW(NP,:)
         
         REFLECTED = .false.
         IF(.NOT.INSIDE_DOMAIN) THEN
            IJK_CELL = PIJK(NP, 4) 
            DIST_OFFSET = ZERO 
            IF(REFLECT_FROM_ORIG_CELL)  THEN 
               IJK_CELL = PIJK_OLD(4) 
            ENDIF
            REFLECTING_CELL = IJK_CELL 
            CALL MPPIC_APPLY_BC2PART(NP, IJK_CELL)
            reflected  = .true. 
         ENDIF

         IF(.NOT.PEA(NP, 1)) CYCLE PART_LOOP
         PIJK_INT(1:5) = PIJK(NP, 1:5) 

         CALL MPPIC_FIND_NEW_CELL(NP)

         CALL MPPIC_CHECK_IF_INSIDE_DOMAIN(NP)

         IF(.NOT.INSIDE_DOMAIN)THEN 
            WRITE(1000,'(i10,6(2x,L1),9(2x,i10),6(2x, g17.8))') NP, &
            CUT_CELL_AT(PIJK_OLD(4)), CUT_CELL_AT(PIJK_INT(4)),     &
            CUT_CELL_AT(PIJK(NP,4)), FLUID_AT(PIJK_OLD(4)),         &
            FLUID_AT(PIJK_INT(4)), FLUID_AT(PIJK(NP,4)), PIJK_OLD(1:4),&
            PIJK_INT(4), PIJK(NP,1:4), (POS_INIT(IDIM), IDIM = 1, DIMN) ,&
            (DES_POS_NEW(NP,IDIM), IDIM = 1, DIMN)
            !IF(NP.Eq.2614) 
            write(1000,'(A, i10,2x,i3, 2(2x, L1))') 'IJK_CELL, REF COUNT, REF FROM ORIG?   =', IJK_CELL, REFLECT_COUNT, tmp_logical 
            
            !IF(tmp_logical.and.CUT_CELL_AT(IJK_CELL)) THEN 
            IJK_CELL = PIJK(NP, 4) 
            IF(CUT_CELL_AT(IJK_CELL)) THEN 
            
               XPOS = DES_POS_NEW(NP,1) 
               YPOS = DES_POS_NEW(NP,2)
               ZPOS = ZERO 
               IF (DIMN .EQ. 3) THEN
                  ZPOS = DES_POS_NEW(NP,3)
               ENDIF
            
            
               CALL GET_DEL_H_DES(IJK_CELL,'SCALAR',XPOS , YPOS, ZPOS,& 
               & DIST, NORM_CF(1), NORM_CF(2), NORM_CF(3), .true.)
               
               write(1000,'(A, 2(2x, g17.8))') 'DIST FROM CUTFACE  =', DIST, DIST-DIST_OFFSET
               write(1000,'(A, 3(2x, g17.8))') 'CUT FACE NORMAL  =', NORM_CF(1:3) 
               write(1000,'(A, 3(2x, g17.8))') 'CUT FACE REFPNT  =', REFP_S(IJK_CELL,1:3)
            ENDIF
            
            WRITE(1000, '(A/)') '--------------------------------------------------------'
            
             IF(PRINT_DES_SCREEN)  WRITE(*,'(i10,6(2x,L1), 3(2x,i10))') NP, CUT_CELL_AT(PIJK_OLD(4)), CUT_CELL_AT(PIJK_INT(4)), CUT_CELL_AT(PIJK(NP,4)), FLUID_AT(PIJK_OLD(4)), FLUID_AT(PIJK_INT(4)), FLUID_AT(PIJK(NP,4)), PIJK_OLD(4), PIJK_INT(4), PIJK(NP,4)
         ENDIF
      END DO PART_LOOP
      
      CLOSE(1000, status='keep')
      !CALL mfix_exit(mype)
      PIP = PIP - PIP_DEL_COUNT
      
      LPIP_DEL_COUNT_ALL(:) = 0
      LPIP_DEL_COUNT_ALL(mype) = PIP_DEL_COUNT 
      CALL GLOBAL_ALL_SUM(LPIP_DEL_COUNT_ALL) 

      IF((DMP_LOG).AND.SUM(LPIP_DEL_COUNT_ALL(:)).GT.0) THEN 
         IF(PRINT_DES_SCREEN) WRITE(*,'(/,2x,A,2x,i10)') 'TOTAL NUMBER OF PARTICLES DELETED GLOBALLY = ', SUM(LPIP_DEL_COUNT_ALL(:))
         WRITE(UNIT_LOG,'(/,2x,A,2x,i10)') 'TOTAL NUMBER OF PARTICLES DELETED GLOBALLY = ', SUM(LPIP_DEL_COUNT_ALL(:))
         DO IPROC = 0, NUMPES-1 
            WRITE(UNIT_LOG, '(/,A,i4,2x,A,2x,i5)') 'PARTICLES DELETED ON PROC:', IPROC,' EQUAL TO', LPIP_DEL_COUNT_ALL(IPROC)
         ENDDO
      ENDIF
      
      !now apply the mass inflow bc 
      CALL MPPIC_MI_BC 
     
! The following lines could be moved to MPPIC_MI_BC
      LPIP_ADD_COUNT_ALL(:) = 0
      LPIP_ADD_COUNT_ALL(mype) = PIP_ADD_COUNT 
      CALL GLOBAL_ALL_SUM(LPIP_ADD_COUNT_ALL) 

      !WRITE(*,*) 'PIP2 = ', PIP, PIP_ADD_COUNT,  LPIP_ADD_COUNT_ALL
      !WRITE(*,*) 'PIP3 = ', LPIP_ADD_COUNT_ALL
      IF((DMP_LOG).AND.SUM(LPIP_ADD_COUNT_ALL(:)).GT.0) THEN 
         IF(PRINT_DES_SCREEN) WRITE(*,'(/,2x,A,2x,i10)') 'TOTAL NUMBER OF PARTICLES ADDED GLOBALLY = ', SUM(LPIP_ADD_COUNT_ALL(:))
         WRITE(UNIT_LOG,'(/,2x,A,2x,i10)') 'TOTAL NUMBER OF PARTICLES ADDED GLOBALLY = ', SUM(LPIP_ADD_COUNT_ALL(:))
         DO IPROC = 0, NUMPES-1 
            WRITE(UNIT_LOG, '(/,A,i4,2x,A,2x,i5)') 'PARTICLES ADDED ON PROC:', IPROC,' EQUAL TO', LPIP_ADD_COUNT_ALL(IPROC)
         ENDDO
         !READ(*,*)
      ENDIF

      
      end SUBROUTINE MPPIC_APPLY_WALLBC
      
      SUBROUTINE  MPPIC_MI_BC 
      IMPLICIT NONE 
      
      INTEGER :: L, I, J, K, IDIR, IDIM, IPCOUNT, LAST_EMPTY_SPOT, NEW_SPOT
      INTEGER :: I1, I2, J1, J2, K1, K2, INEW, JNEW, KNEW, IJK, IJK_WALL,  M
      DOUBLE PRECISION :: WALL_NORM(DIMN), CORD_START(DIMN), DOML(DIMN)
      DOUBLE PRECISION :: EPS_INFLOW(MMAX), AREA_INFLOW, VEL_INFLOW(DIMN), STAT_WT
      
      DOUBLE PRECISION REAL_PARTS(DIM_M), COMP_PARTS(DIM_M), VEL_NORM_MAG, VOL_INFLOW, VOL_IJK
       
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: RANDPOS
      INTEGER :: CNP_CELL_COUNT
      integer :: lglobal_id
      integer, dimension(0:numpes-1) :: add_count_all
      type :: ty_spotlist
         integer spot
         type(ty_spotlist),pointer :: next
      end type ty_spotlist
      type(ty_spotlist),pointer :: root_spotlist,cur_spotlist,prev_spotlist
       
      INCLUDE 'function.inc'
       
      PIP_ADD_COUNT = 0 
      LAST_EMPTY_SPOT = 0
      allocate(root_spotlist); nullify(root_spotlist%next)
      cur_spotlist => root_spotlist

      DO L = 1, DIMENSION_BC 
         IF (BC_DEFINED(L)) THEN 
            IF (BC_TYPE(L)=='MASS_INFLOW' ) THEN 
               !dont add the MI to des_bc information 
               I1 = BC_I_W(L) 
               I2 = BC_I_E(L) 
               J1 = BC_J_S(L) 
               J2 = BC_J_N(L) 
               K1 = BC_K_B(L) 
               K2 = BC_K_T(L) 
               DO K = K1, K2 
                  DO J = J1, J2 
                     DO I = I1, I2  
!                        IF (.NOT.IS_ON_myPE_owns(I, J, K)) CYCLE
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
                        !wall_norms are pointing away from the fluid.


                        IJK = FUNIJK(INEW,JNEW,KNEW) !this is the cell where the new particles will be intialized  
                        IJK_WALL = FUNIJK(I,J,K)  !this is the id of the new bc 
                        
                        IF(.NOT.FLUID_AT(IJK)) CYCLE 
                        
                        !Find the direction of the normal 
                        IDIR = 0
                        DO IDIM = 1, DIMN
                           IDIR = IDIR + ABS(WALL_NORM(IDIM))*IDIM 
                        end DO
                        
                        CORD_START(1) = XE(INEW-1)
                        CORD_START(2) = YN(JNEW-1)
                        IF(DIMN.EQ.3) CORD_START(3) = ZT(KNEW-1)
                        DOML(1) = DX(INEW)
                        DOML(2) = DY(JNEW)
                        IF(DIMN.EQ.3) DOML(3) = DZ(KNEW)

                        IF(DIMN.eq.2) THEN 
                           !compute the normal area 
                           AREA_INFLOW = DOML(1)*DOML(2)*ZLENGTH/DOML(IDIR)
                           VOL_IJK = DOML(1)*DOML(2)*ZLENGTH
                        ELSE
                           AREA_INFLOW = DOML(1)*DOML(2)*DOML(3)/DOML(IDIR)
                           
                           VOL_IJK = DOML(1)*DOML(2)*DOML(3)
                        ENDIF

                        IF(WALL_NORM(IDIR).GT.0) THEN 
                           CORD_START(IDIR) = CORD_START(IDIR) + DOML(IDIR)
                        ENDIF
                        DOML(IDIR) = ZERO 
                        !set this to zero as the particles will 
                        !be seeded only on the BC plane 

                        !WRITE(*,*) 'BC_P =', BC_PLANE(L), CORD_START(2), WALL_NORM(2)

!let's say if the wall is a east wall, then the wall cordinate 
!will be XE(INEW) which corresponds to XE(INEW-1) + DX(INEW). 
                        
                        DO M = 1, MMAX 
                           EPS_INFLOW(M) = BC_ROP_S(L, M)/RO_S(M)
                           VEL_INFLOW(1) = BC_U_S(L, M)
                           VEL_INFLOW(2) = BC_V_S(L, M)
                           IF(DIMN.eq.3) VEL_INFLOW(3) = BC_W_S(L, M)
                           VEL_NORM_MAG = ABS(DOT_PRODUCT(VEL_INFLOW(1:DIMN), WALL_NORM(1:DIMN)))
                           VOL_INFLOW = AREA_INFLOW*VEL_NORM_MAG*DTSOLID
                           
                           REAL_PARTS(M) = 6.d0*EPS_INFLOW(M)*VOL_INFLOW/(PI*(D_p0(M)**3.d0))
                           COMP_PARTS(M) = zero
                           IF(CONSTANTNPC) THEN 
                              IF(EPS_INFLOW(M).GT.ZERO) COMP_PARTS(M) = REAL(NPC_PIC(M))*VOL_INFLOW/VOL_IJK
                           ELSEIF(CONSTANTWT) THEN
                              COMP_PARTS(M) = REAL_PARTS(M)/REAL(STATWT_PIC(M))
                           ENDIF
                           CNP_ARRAY(IJK,0) = CNP_ARRAY(IJK,0) + REAL_PARTS(M)
                           CNP_ARRAY(IJK,M) = CNP_ARRAY(IJK,M) + COMP_PARTS(M)
                           

                           IF(CNP_ARRAY(IJK,M).GE.1.d0) THEN 
                              CNP_CELL_COUNT = INT(CNP_ARRAY(IJK,M))
                              CNP_ARRAY(IJK,M) =  CNP_ARRAY(IJK,M) - CNP_CELL_COUNT
                              REAL_PARTS(M) = CNP_ARRAY(IJK, 0)
                              !set cnp_array(ijk,0) to zero to reflect that all real particles have been seeded 
                              CNP_ARRAY(IJK,0) = 0.d0 
                              
                              ALLOCATE(RANDPOS(CNP_CELL_COUNT, DIMN))
                              RANDPOS = ZERO 
                              DO IDIM = 1, DIMN
                                 CALL UNI_RNO(RANDPOS(1:CNP_CELL_COUNT, IDIM))
                              ENDDO

                              IF(CONSTANTNPC) THEN 
                                 !calculate the statistical weight
                                 !for CP's belonging to this 
                                 !solid phase 

                                 STAT_WT = REAL_PARTS(M)/REAL(CNP_CELL_COUNT)
                              ELSEIF(CONSTANTWT) THEN
                                 STAT_WT = STATWT_PIC(M) 
                              ENDIF
                              
                              
                              DO IPCOUNT = 1, CNP_CELL_COUNT
                                 
                                 CALL MPPIC_FIND_EMPTY_SPOT(LAST_EMPTY_SPOT, NEW_SPOT)
                                 
                                 DES_POS_OLD(NEW_SPOT, :) =  CORD_START(:) + RANDPOS(IPCOUNT, :)*DOML(:)
                                 DES_POS_NEW(NEW_SPOT, :) = DES_POS_OLD(NEW_SPOT, :)
                                 DES_VEL_OLD(NEW_SPOT, :) = VEL_INFLOW(:)
                                 
                                 DES_VEL_NEW(NEW_SPOT, :) = DES_VEL_OLD(NEW_SPOT, :)
                                 
                                 DES_RADIUS(NEW_SPOT) = D_p0(M)*HALF
                                 
                                 RO_Sol(NEW_SPOT) =  RO_S(M)
                                 
                                 DES_STAT_WT(NEW_SPOT) = STAT_WT 

                                 MARK_PART(NEW_SPOT) = myPE
                                 
                                 PIJK(NEW_SPOT, 1) = INEW
                                 PIJK(NEW_SPOT, 2) = JNEW
                                 PIJK(NEW_SPOT, 3) = KNEW
                                 PIJK(NEW_SPOT, 4) = IJK
                                 PIJK(NEW_SPOT, 5) = M 
                                 
                                 PVOL(NEW_SPOT) = (4.0d0/3.0d0)*Pi*DES_RADIUS(L)**3
                                 PMASS(NEW_SPOT) = PVOL(NEW_SPOT)*RO_SOL(NEW_SPOT) 
                                 CALL MPPIC_CHECK_IF_INSIDE_DOMAIN(NEW_SPOT)
                                 
                                 IF(INSIDE_DOMAIN) THEN 
                                    PIP = PIP+1
                                    PIP_ADD_COUNT = PIP_ADD_COUNT + 1 
                                    PEA(NEW_SPOT, 1) = .true. 
                                    PEA(NEW_SPOT, 2:4) = .false. 
                                 ! addd to the list 
                                    cur_spotlist%spot = new_spot
                                    allocate(cur_spotlist%next)
                                    cur_spotlist => cur_spotlist%next
                                    nullify(cur_spotlist%next)
                                 ELSE
                                    PEA(NEW_SPOT, 1) = .false. 
                                    LAST_EMPTY_SPOT = NEW_SPOT - 1 
                                 ENDIF


                                 
                                 
                                 !WRITE(*,'(A,2(2x,i5), 2x, A,2x, 3(2x,i2),2x, A, 3(2x,g17.8))') 'NEW PART AT ', NEW_SPOT, MAX_PIP, 'I, J, K = ', INEW, JNEW, KNEW, 'POS =', DES_POS_NEW(NEW_SPOT,:)
                                 !IF(DMP_LOG) WRITE(UNIT_LOG,'(A,2x,i5, 2x, A,2x, 3(2x,i2),2x, A, 3(2x,g17.8))') 'NEW PART AT ', NEW_SPOT, 'I, J, K = ', INEW, JNEW, KNEW, 'POS =', DES_POS_NEW(NEW_SPOT,:) 

                                 !WRITE(*,*) 'IDIR, DOML = ', IDIR, DOML(:)
                              ENDDO
                              DEALLOCATE(RANDPOS)
                              
                           end IF
                        end DO
                     end DO
                  end DO
               end DO
                 
            ELSE 
               !do nothing 
            END IF
         END IF
      END DO
! Pradeep : Assign global id to new particles added 
      add_count_all(:) = 0 
      add_count_all(mype) = pip_add_count 
      call global_all_sum(add_count_all(0:numpes-1))
      lglobal_id = imax_global_id + sum(add_count_all(0:mype-1))
      cur_spotlist=>root_spotlist
      do l = 1,pip_add_count 
         lglobal_id = lglobal_id + 1
         iglobal_id(cur_spotlist%spot)= lglobal_id  
         prev_spotlist=> cur_spotlist
         cur_spotlist => cur_spotlist%next
         deallocate(prev_spotlist)
      end do  
      deallocate(cur_spotlist)
      imax_global_id = imax_global_id+sum(add_count_all(0:numpes-1))
 
      END SUBROUTINE MPPIC_MI_BC
    
      SUBROUTINE MPPIC_FIND_EMPTY_SPOT(LAST_INDEX, EMPTY_SPOT)
      
      IMPLICIT NONE 
      INTEGER, INTENT(INOUT) :: LAST_INDEX 
      INTEGER, INTENT(OUT) :: EMPTY_SPOT 
      LOGICAL :: SPOT_FOUND
      INTEGER :: LL
      
      IF(LAST_INDEX.EQ.MAX_PIP) THEN
         IF(DMP_LOG) then 
            WRITE(UNIT_LOG,2001) 
            IF(mype.eq.pe_IO) write(*,2001) 
         ENDIF
         call mfix_exit(mype)
      ENDIF
      SPOT_FOUND = .false. 

      DO LL = LAST_INDEX+1, MAX_PIP 
         !Rahul:
         !for Pradeep: Im not clear how to proceed here as I havent
         !yet scene or gone over the particle exchange information 
         !for inflow/outflow BC 
         if(.NOT.PEA(LL,1)) THEN 
            EMPTY_SPOT = LL
            LAST_INDEX = LL 
            SPOT_FOUND = .true. 
            EXIT 
         ENDIF
      ENDDO

      IF(.NOT.SPOT_FOUND) THEN 
         IF(DMP_LOG) then 
            WRITE(UNIT_LOG,2002) 
            IF(mype.eq.pe_IO) write(*,2002) 
         ENDIF
         call mfix_exit(mype)
      ENDIF
      
 2001 FORMAT(/1X,70('*'),//,10X,  & 
      & 'ERROR IN MPPIC_FIND_EMPTY_SPOT', /10X, &
      & 'NO MORE EMPTY SPOT IN THE PARTICLE ARRAY TO ADD A NEW PARTICLE',/10X &
      & 'TERMINAL ERROR: STOPPING (CALL RESTART NOT YET IMPLEMENTED BEFORE THIS EXIT)', &
      & /1X,70('*')/)
      
 2002 FORMAT(/1X,70('*'),//,10X,  & 
      & 'ERROR IN MPPIC_FIND_EMPTY_SPOT', /10X, &
      & 'COULD NOT FIND A SPOT FOR ADDING NEW PARTICLE',/10X &
      & 'INCREASE THE SIZE OF THE INITIAL ARRAYS', 10X, &
      & 'TERMINAL ERROR: STOPPING (CALL RESTART NOT YET IMPLEMENTED BEFORE THIS EXIT)', &
      & /1X,70('*')/)
      END SUBROUTINE MPPIC_FIND_EMPTY_SPOT

      SUBROUTINE  MPPIC_APPLY_BC2PART(LL, IJK_CELL)
      IMPLICIT NONE 
      
      INTEGER, INTENT(IN) :: LL, IJK_CELL
      INTEGER :: I, J, K, IDIM, IDIR, IJK
      INTEGER :: COUNT, COUNT_BC, IJK_WALL
      
      DOUBLE PRECISION :: NORMAL(DIMN), WALL_COOR(DIMN), DIST, WALLCOR_MIN(DIMN), WALLCOR_MAX(DIMN) 
      
      DOUBLE PRECISION :: NORM_CF(3), XPOS, YPOS, ZPOS
      DOUBLE PRECISION :: XPOS_OLD, YPOS_OLD, ZPOS_OLD, DIST_OLD
      CHARACTER*100 :: WALL_TYPE


      INCLUDE 'function.inc'


      COUNT_BC = DES_CELLWISE_BCDATA(IJK_CELL)%COUNT_DES_BC 

      IF(COUNT_BC.EQ.0) THEN 
         
         
         IF(CARTESIAN_GRID) THEN 
            IF(DMP_LOG) WRITE(UNIT_LOG,1002) 'APPLY_WALLBC_PIC', INSIDE_DOMAIN, INSIDE_SMALL_CELL, REFLECT_FROM_ORIG_CELL, IJK_CELL,& 
                 I_OF(IJK_CELL), J_OF(IJK_CELL), K_OF(IJK_CELL),  FLUID_AT(IJK) , CUT_CELL_AT(IJK), SMALL_CELL_AT(IJK), &
                 PIJK_OLD(4),PIJK_OLD(1:3), FLUID_AT(PIJK_OLD(4)) , CUT_CELL_AT(PIJK_OLD(4)) , SMALL_CELL_AT(PIJK_OLD(4)) 
 
            IF(PRINT_DES_SCREEN) WRITE(*,1002) 'APPLY_WALLBC_PIC',  INSIDE_DOMAIN, INSIDE_SMALL_CELL, REFLECT_FROM_ORIG_CELL, IJK_CELL, & 
                 I_OF(IJK_CELL), J_OF(IJK_CELL), K_OF(IJK_CELL),  FLUID_AT(IJK) , CUT_CELL_AT(IJK), SMALL_CELL_AT(IJK), &
                 PIJK_OLD(4),PIJK_OLD(1:3), FLUID_AT(PIJK_OLD(4)) , CUT_CELL_AT(PIJK_OLD(4)) , SMALL_CELL_AT(PIJK_OLD(4)) 

            
         ELSE
            IF(DMP_LOG) WRITE(UNIT_LOG,1001) 'APPLY_WALLBC_PIC', INSIDE_DOMAIN, INSIDE_SMALL_CELL, REFLECT_FROM_ORIG_CELL, IJK_CELL, &
                 I_OF(IJK_CELL), J_OF(IJK_CELL), K_OF(IJK_CELL), FLUID_AT(IJK), &
                 PIJK_OLD(4),PIJK_OLD(1:3), FLUID_AT(PIJK_OLD(4)) 
            IF(PRINT_DES_SCREEN) WRITE(*,1001) 'APPLY_WALLBC_PIC', INSIDE_DOMAIN, INSIDE_SMALL_CELL, REFLECT_FROM_ORIG_CELL,  IJK_CELL,& 
                 I_OF(IJK_CELL), J_OF(IJK_CELL), K_OF(IJK_CELL), FLUID_AT(IJK), &
                 PIJK_OLD(4),PIJK_OLD(1:3), FLUID_AT(PIJK_OLD(4)) 
         ENDIF
         CALL mfix_exit(myPE)
      ENDIF

      XPOS_OLD = DES_POS_OLD(LL, 1) 
      YPOS_OLD = DES_POS_OLD(LL, 2) 
      ZPOS_OLD = ZERO
      IF(DIMN.eq.3) ZPOS_OLD = DES_POS_OLD(LL, 3)

      DO COUNT = 1, COUNT_BC
         IF(.NOT.PEA(LL,1)) CYCLE
         !a particle cud be marked as inactive during this loop
         !as it may encounter an outflow BC 

         IJK_WALL  = DES_CELLWISE_BCDATA(IJK_CELL)%BDRY_LIST(COUNT)%IJK_SCAL 
         WALL_TYPE = DES_CELLWISE_BCDATA(IJK_CELL)%BDRY_LIST(COUNT)%DES_BC_TYPE
         
         DIST_OFFSET = ZERO 

         SELECT CASE (TRIM(WALL_TYPE)) 
         
         CASE('NORMAL_WALL')
            !Use IJK_CELL below because IJK_WALL will point to a 
            !ghost cell in this case 
            I = I_OF(IJK_CELL)
            J = J_OF(IJK_CELL)
            K = K_OF(IJK_CELL)
            
            WALLCOR_MIN(1) = XE(I-1)
            WALLCOR_MAX(1) = XE(I)
            WALLCOR_MIN(2) = YN(J-1)
            WALLCOR_MAX(2) = YN(J)
            
            IF(DIMN.EQ.3) THEN 
               WALLCOR_MIN(3) = ZT(K-1)
               WALLCOR_MAX(3) = ZT(K)
            END IF
            NORMAL(1:DIMN) = DES_CELLWISE_BCDATA(IJK_CELL)%BDRY_LIST(COUNT)%NORMAL(1:DIMN)
            !Find the direction of the normal 
            IDIR = 0
            DO IDIM = 1, DIMN
               IDIR = IDIR + ABS(NORMAL(IDIM))*IDIM 
            end DO

            WALL_COOR(1:DIMN)  = WALLCOR_MIN(1:DIMN)
            
            IF(NORMAL(IDIR).GT.0) WALL_COOR(IDIR) = WALLCOR_MAX(IDIR)
!let's say if the wall is the east wall of this scalar cell, then the wall cordinate 
!will be XE(I) which corresponds to WALLCOR_MAX 
               
            DIST =  NORMAL(IDIR)*(DES_POS_NEW(LL, IDIR) - WALL_COOR(IDIR))
            
            
            !WRITE(*,*) 'LL  =', LL, WALL_COOR(:), NORMAL(:), DIST
!according to the above convention for normal for 'normal_walls', 
!distance will be negative for particles inside the domain and 
!positive for particles outside the domain. This is becuase
!the wall normal points away from the physical domain, and any 
!point inside the domain will have negative distance. 
            IF(DIST.GT.ZERO) CALL MPPIC_REFLECT_PART(LL, (DIST), -NORMAL(1:DIMN))
!in the routine mppic_reflect_part, it is assumed that the normal 
!points into the system (i.e., toward the fluid, hence the negative
!sign in front of normal)

         CASE('SMALL_CELL')
            !Use either IJK_CELL or IJK_WALL coz both are same 
            !ghost cell in this case 
            I = I_OF(IJK_CELL)
            J = J_OF(IJK_CELL)
            K = K_OF(IJK_CELL)
            
            WALLCOR_MIN(1) = XE(I-1)
            WALLCOR_MAX(1) = XE(I)
            WALLCOR_MIN(2) = YN(J-1)
            WALLCOR_MAX(2) = YN(J)
            
            IF(DIMN.EQ.3) THEN 
               WALLCOR_MIN(3) = ZT(K-1)
               WALLCOR_MAX(3) = ZT(K)
            END IF
            NORMAL(1:DIMN) = DES_CELLWISE_BCDATA(IJK_CELL)%BDRY_LIST(COUNT)%NORMAL(1:DIMN)
            !Find the direction of the normal 
            IDIR = 0
            DO IDIM = 1, DIMN
               IDIR = IDIR + ABS(NORMAL(IDIM))*IDIM 
            end DO

            WALL_COOR(1:DIMN)  = WALLCOR_MIN(1:DIMN)
            
            IF(NORMAL(IDIR).GT.0) WALL_COOR(IDIR) = WALLCOR_MAX(IDIR)
!let's say if the wall is the east wall of this scalar cell, then the wall cordinate 
!will be XE(I) which corresponds to WALLCOR_MAX 
               
            DIST =  NORMAL(IDIR)*(DES_POS_NEW(LL, IDIR) - WALL_COOR(IDIR))
            
!!$            IF(IJK_CELL.eq.3413) THEN 
!!$               WRITE(*,*) 'SMALL CELL ? ', SMALL_CELL_AT(IJK_CELL), LL
!!$               WRITE(*,*) 'BC COUNT = ', COUNT
!!$               WRITE(*,'(A, 3(2x,g17.8))') 'NORMAL = ', NORMAL(:)
!!$               WRITE(*,'(A, 3(2x,g17.8))') 'POS = ', DES_POS_NEW(LL, :)
!!$               WRITE(*,'(A, 3(2x,g17.8))') 'VEL = ', DES_VEL_NEW(LL, :)
!!$               WRITE(*,'(A, 3(2x,g17.8))') 'DIST = ', DIST
!!$            ENDIF

            IF(DIST.LT.ZERO) CALL MPPIC_REFLECT_PART(LL, ABS(DIST), NORMAL(1:DIMN))
!recall the normal for small cell points toward the fluid. 
!!$            IF(IJK_CELL.eq.3413) THEN 
!!$               WRITE(*,'(A, 3(2x,g17.8))') 'POS PC = ', DES_POS_NEW(LL, :)
!!$               WRITE(*,'(A, 3(2x,g17.8))') 'VEL PC= ', DES_VEL_NEW(LL, :)
!!$               READ(*,*) 
!!$            ENDIF

            
            
         CASE('CUT_FACE')
            XPOS = DES_POS_NEW(LL,1) 
            YPOS = DES_POS_NEW(LL,2)
            ZPOS = ZERO 
            IF (DIMN .EQ. 3) THEN
               ZPOS = DES_POS_NEW(LL,3)
            ENDIF
            
            
            CALL GET_DEL_H_DES(IJK_WALL,'SCALAR',XPOS , YPOS, ZPOS,& 
            & DIST, NORM_CF(1), NORM_CF(2), NORM_CF(3), .true.)
            
            
            DIST_OFFSET  = 0.02*DES_RADIUS(LL)
            !this is done so as to 
            !prevent particles that are very close to the cut-face
            !but yet out of the domain from returning positive distance
            !and escaping reflection. 
            !Due to the tolerances used in constructing cut-faces, this 
            !is possible and can lead to particles staying in 
            !non-fluid cells and drifting until they magically appear
            !in other cells.
            DIST = DIST - DIST_OFFSET !this is done so as to 
            !prevent particles that are very close to the cut-face
            !but yet out of the domain from returning positive distance.
            !Due to the tolerances used in constructing cut-faces, this 
            !is possible. 
            NORMAL(1:DIMN)  = NORM_CF(1:DIMN) 
!cut-face normal points into the system. Therefore, negative 
!distance will imply a particle out of domain 
            IF(DIST.LT.ZERO) THEN 
               !CALL GET_DEL_H_DES(IJK_WALL,'SCALAR',XPOS_OLD , YPOS_OLD, ZPOS_OLD,& 
               !& DIST_OLD, NORM_CF(1), NORM_CF(2), NORM_CF(3), .true.)

               !IF(DIST_OLD.GE.0.d0)
               CALL MPPIC_REFLECT_PART(LL, ABS(DIST), NORMAL(1:DIMN))
            ENDIF

         CASE ('OUTFLOW')
            !This case is similar to WALL_NORMAL in terms of the 
            !checks for out of bounds 
            
            !Use IJK_CELL below because IJK_WALL will point to a 
            !ghost cell in this case 
            !WRITE(*,*) 'DETECTED AN OUTFLOW CELL'
            !READ(*,*)
            I = I_OF(IJK_CELL)
            J = J_OF(IJK_CELL)
            K = K_OF(IJK_CELL)
            
            WALLCOR_MIN(1) = XE(I-1)
            WALLCOR_MAX(1) = XE(I)
            WALLCOR_MIN(2) = YN(J-1)
            WALLCOR_MAX(2) = YN(J)
            
            IF(DIMN.EQ.3) THEN 
               WALLCOR_MIN(3) = ZT(K-1)
               WALLCOR_MAX(3) = ZT(K)
            END IF
            NORMAL(1:DIMN) = DES_CELLWISE_BCDATA(IJK_CELL)%BDRY_LIST(COUNT)%NORMAL(1:DIMN)
            !Find the direction of the normal 
            IDIR = 0
            DO IDIM = 1, DIMN
               IDIR = IDIR + ABS(NORMAL(IDIM))*IDIM 
            end DO

            WALL_COOR(1:DIMN)  = WALLCOR_MIN(1:DIMN)
            
            IF(NORMAL(IDIR).GT.0) WALL_COOR(IDIR) = WALLCOR_MAX(IDIR)
!let's say if the wall is a east wall, then the wall cordinate 
!will be XE(I) which corresponds to WALLCOR_MAX 
            
            DIST =  NORMAL(IDIR)*(DES_POS_NEW(LL, IDIR) - WALL_COOR(IDIR))
            !same convention for normal as explained in normal walls
            IF(DIST.GT.ZERO) THEN 
               !currently i'm just flagging a particle as inactive
               !if it is out of domain. In order to maintain the 
               !contiguity of arrays, it is better to move the 
               !last particle to this location 
               ! 
               PEA(LL, 1)  = .false. 
               PIP_DEL_COUNT = PIP_DEL_COUNT+1
            ENDIF
         CASE DEFAULT
            IF(DMP_LOG.OR.mype.eq.pe_IO) THEN 
               WRITE(UNIT_LOG, 1035) 'MPPIC_APPLY_BC2PART', TRIM(WALL_TYPE)
               WRITE(*, 1035) 'MPPIC_APPLY_BC2PART', TRIM(WALL_TYPE)
            ENDIF
            CALL mfix_exit(myPE)
         END SELECT
      END DO

 1001    FORMAT(/1X,70('*'),//,10X,  & 
         & 'ERROR IN SUBROUTINE....', A, /,10X, &
         & 'INSIDE_DOMAIN, INSIDE_SMALL_CELL, REFLECT_FROM_ORIG_CELL ?', 3(2x, L1) , /10X, &
         & 'COUNT_BC SHUD HAVE BEEN GE 1 FOR THIS CELL', /, 10X, & 
         & 'IJK, I, J, K  = ', 4(2x,i10), /10x, & 
         & 'FLUID_AT IJK  ? ', 1(2x, L1) , /10X, &
         & 'OLD CELL IJK, I, J, K = ',  4(2x,i10), /10x, & 
         & 'FLUID_AT OLD IJK  ? ', 1(2x, L1) , /10X, &
         & 'TERMINAL ERROR: STOPPING', &
         & /1X,70('*')/)
         
 1002    FORMAT(/1X,70('*'),//,10X,  & 
         & 'ERROR IN SUBROUTINE....', A, /,10X, &
         & 'INSIDE_DOMAIN, INSIDE_SMALL_CELL, REFLECT_FROM_ORIG_CELL ?', 3(2x, L1) , /10X, &
         & 'COUNT_BC SHUD HAVE BEEN GE 1 FOR THIS CELL', /, 10X, & 
         & 'IJK, I, J, K  = ', 4(2x,i10), /10x, & 
         & 'FLUID_AT, CUT_CELL_AT, SMALL_CELL_AT IJK  ? ', 3(2x, L1) , /10X,  &
         & 'OLD CELL IJK, I, J, K = ', 4(2x,i10), /10x, & 
         & 'FLUID_AT, CUT_CELL_AT, SMALL_CELL_AT IJK  ? ', 3(2x, L1) , /10X,  &
         & 'TERMINAL ERROR: STOPPING', &
         & /1X,70('*')/)
         
 1035 FORMAT(/1X,70('*'),//,10X,  & 
         & 'ERROR IN SUBROUTINE....', A, /,10X, &
         & 'THIS BC....(', A,')....FOR PARTICLE-WALL NOT SUPPORTED YET   ', /, 10X, & 
         & 'ACCEPTABLE TYPES ARE:', /10x, & 
         & 'NORMAL_WALL' , /10X, &
         & 'CUT_FACE', /10X, &
         & 'OUTFLOW', /10X, & 
         & 'TERMINAL ERROR: STOPPING', &
         & /1X,70('*')/)



      END SUBROUTINE  MPPIC_APPLY_BC2PART

      SUBROUTINE MPPIC_REMOVE_PARTICLE(LL) 
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: LL
        
        INTEGER :: II 

        LOGICAL ::  REPLACE_PART_FOUND 
        REPLACE_PART_FOUND = .TRUE.
        !Don't call this routine in its current form. 
        !this routine will bring the IIth particle to LL's position
        !But doing so, this particle will not be tested for
        !BC treatment, since the PART_LOOP in main routine is
        !already past LL 

        IILOOP: DO II = MAX_PIP, LL+1, -1 
           IF(PEA(II,1).AND..NOT.PEA(II,4)) THEN  
              !II exists and also belongs to this processor,
              !i.e., not marked as a ghost particle 
              REPLACE_PART_FOUND = .TRUE.
              EXIT  IILOOP
           ENDIF
        END DO IILOOP 
        IF (REPLACE_PART_FOUND) THEN 
           !CALL MPPIC_SWAP_PARTS(LL, II)
        ELSE
           PEA(LL,1) = .FALSE.
        ENDIF
        
        PIP_DEL_COUNT = PIP_DEL_COUNT + 1
       
        
      END SUBROUTINE MPPIC_REMOVE_PARTICLE
      
      SUBROUTINE MPPIC_REFLECT_PART(LL, DISTFROMWALL, WALL_NORM)
              
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: LL
      DOUBLE PRECISION, INTENT(IN) :: DISTFROMWALL, WALL_NORM(DIMN)
      !magnitude of pre-collisional normal and tangential velocity components 
      DOUBLE PRECISION :: VEL_NORMMAG_APP, VEL_TANGMAG_APP, TANGENT(DIMN)
      
      !pre collisional normal and tangential velocity components in vector format
      !APP ==> approach 
      DOUBLE PRECISION :: VEL_NORM_APP(DIMN), VEL_TANG_APP(DIMN)

      
      !post collisional normal and tangential velocity components in vector format
      !SEP ==> separation 
      DOUBLE PRECISION :: VEL_NORM_SEP(DIMN), VEL_TANG_SEP(DIMN)

      DOUBLE PRECISION :: COEFF_REST_EN, COEFF_REST_ET

      DOUBLE PRECISION :: VELMOD, NORMMOD, VELDOTNORM 
!bring the particle inside the domain 
      REFLECT_COUNT = REFLECT_COUNT + 1 
      DES_POS_NEW(LL, 1:DIMN) = DES_POS_NEW(LL, 1:DIMN) + &
           & 1.001d0*DISTFROMWALL*WALL_NORM(1:DIMN) 

      !VELMOD = SQRT(DOT_PRODUCT(DES_VEL_NEW(LL,1:DIMN), DES_VEL_NEW(LL, 1:DIMN))    
      !NORMMOD = SQRT(DOT_PRODUCT(WALL_NORM(1:DIMN), WALL_NORM( 1:DIMN))
      !sometime the magnitude of the normal from cut-face may not be
      !exactly equal to one 
      
      VELDOTNORM = DOT_PRODUCT(DES_VEL_NEW(LL,1:DIMN), WALL_NORM(1:DIMN))
      !VELDOTNORM = VELDOTNORM/(NORMMOD*VELMOD)

      
      IF(VELDOTNORM.GT.ZERO) THEN
         WRITE(*,*) 'NO NEED TO REFLECT VELOCITY: L = ', LL
         if(dmp_log) WRITE(unit_log,*) 'NO NEED TO REFLECT VELOCITY: L = ', LL
         
         RETURN    
                                !in this case only correct the position 
      ENDIF
      !and do not reflect the velocity 
      !this kind of situation shud happen very rarely. 

      VEL_NORMMAG_APP = DOT_PRODUCT(WALL_NORM(1:DIMN), DES_VEL_NEW(LL, 1:DIMN))

!currently assuming that wall is at rest. Needs improvement for moving wall 
      
      VEL_NORM_APP(1:DIMN) = VEL_NORMMAG_APP*WALL_NORM(1:DIMN)

      VEL_TANG_APP(:) = DES_VEL_NEW(LL, 1:DIMN) - VEL_NORM_APP(1:DIMN)

      VEL_TANGMAG_APP = SQRT(DOT_PRODUCT(VEL_TANG_APP(1:DIMN) , VEL_TANG_APP(1:DIMN)))
                  

!post collisional velocities  
                  
      COEFF_REST_EN = REAL_EN_WALL(PIJK(LL,5))
      !if(ep_g(PIJK(LL,4)).lt.0.42) coeff_rest_en = 1.05
      VEL_NORM_SEP(1:DIMN) = -COEFF_REST_EN*VEL_NORM_APP(1:DIMN) 

      VEL_TANG_SEP(1:DIMN) = VEL_TANG_APP(1:DIMN) 

      !IF(MEW_W.GT.ZERO) THEN 
      COEFF_REST_ET = 1.d0!REAL_ET_WALL(PIJK(LL,5))
      VEL_TANG_SEP(1:DIMN) = COEFF_REST_ET*VEL_TANG_APP(1:DIMN) 
      !ENDIF

      DES_VEL_NEW(LL, 1:DIMN) = VEL_NORM_SEP(1:DIMN) + VEL_TANG_SEP(1:DIMN)
                  
      END SUBROUTINE MPPIC_REFLECT_PART

      SUBROUTINE MPPIC_CHECK_IF_INSIDE_DOMAIN(LL) 
      
      INTEGER, INTENT(IN) :: LL 
      
      DOUBLE PRECISION :: XPOS, YPOS, ZPOS, WNORM(3), DISTMOD 
      INTEGER :: I, J, K, IJK
      
      
      INCLUDE 'function.inc'
      IJK = PIJK(LL, 4)

      INSIDE_DOMAIN = .TRUE.
      REFLECT_FROM_ORIG_CELL = .FALSE.
      INSIDE_SMALL_CELL = .false.

      IF(FLUID_AT(IJK)) THEN
         CG: IF(CARTESIAN_GRID) THEN 
            IF(CUT_CELL_AT(IJK)) THEN 
               XPOS = DES_POS_NEW(LL,1) 
               YPOS = DES_POS_NEW(LL,2)
               ZPOS = ZERO 
               IF (DIMN .EQ. 3) THEN
                  ZPOS = DES_POS_NEW(LL,3)
               ENDIF
               
               CALL GET_DEL_H_DES(IJK,'SCALAR', XPOS, YPOS, ZPOS,& 
               & DISTMOD, WNORM(1), WNORM(2), WNORM(3), .true.)
               
               IF(DISTMOD.LT.ZERO) THEN 
                  INSIDE_DOMAIN  = .FALSE. 
               ENDIF
               
            ENDIF 
         ENDIF CG 
      ELSE
         INSIDE_DOMAIN  = .FALSE. 
         REFLECT_FROM_ORIG_CELL = .TRUE. 
         IF(CARTESIAN_GRID) THEN 
            IF(SMALL_CELL_AT(IJK)) THEN 
               INSIDE_SMALL_CELL = .TRUE. 
               REFLECT_FROM_ORIG_CELL = .FALSE. 
            ENDIF
            !In this case reflect using the small cell bc's
         ENDIF
      ENDIF
      END SUBROUTINE MPPIC_CHECK_IF_INSIDE_DOMAIN
      SUBROUTINE MPPIC_FIND_NEW_CELL(LL)
      
      IMPLICIT NONE 
      INTEGER, INTENT(IN) :: LL 
      DOUBLE PRECISION :: XPOS, YPOS, ZPOS 
      INTEGER :: I, J, K, IJK
      
      INCLUDE 'function.inc'
      I = PIJK(LL,1)
      J = PIJK(LL,2)
      K = PIJK(LL,3)
      XPOS = DES_POS_NEW(LL,1) 
      YPOS = DES_POS_NEW(LL,2)
      IF (DIMN .EQ. 3) THEN
         ZPOS = DES_POS_NEW(LL,3)
      ENDIF
      
      IF(XPOS >= XE(I-1) .AND. XPOS < XE(I)) THEN 
         PIJK(LL,1) = I
      ELSEIF(XPOS >= XE(I)) THEN 
         PIJK(LL,1) = I+1
      ELSE 
         PIJK(LL,1) = I-1 
      END IF 

      IF(YPOS >= YN(J-1) .AND. YPOS < YN(J)) THEN 
         PIJK(LL,2) = J
      ELSEIF(YPOS >= YN(J))THEN 
         PIJK(LL,2) = J+1
      ELSE
         PIJK(LL,2) = J-1
      END IF 
      
      IF(DIMN.EQ.2) THEN
         PIJK(LL,3) = 1
      ELSE
         IF(ZPOS >= ZT(K-1) .AND. ZPOS < ZT(K)) THEN
            PIJK(LL,3) = K
         ELSEIF(ZPOS >= ZT(K)) THEN 
            PIJK(LL,3) = K+1
         ELSE
            PIJK(LL,3) = K-1
         END IF 
      ENDIF          
    
      I = PIJK(LL,1)
      J = PIJK(LL,2)
      K = PIJK(LL,3)

      IF(I.EQ.IEND1+1) then
         IF(XPOS >= XE(IEND1-1) .AND. XPOS <= XE(IEND1)) PIJK(LL,1) = IEND1
      ENDIF
      
      IF(J.EQ.JEND1+1) then 
         IF(YPOS >= YN(JEND1-1) .AND. YPOS <= YN(JEND1)) PIJK(LL,2) = JEND1
      ENDIF
      
      IF(DIMN.EQ.3.AND.K.EQ.KEND1+1) THEN 
         IF(ZPOS >= ZT(KEND1-1) .AND. ZPOS <= ZT(KEND1)) PIJK(LL,3) = KEND1
      ENDIF

      PIJK(LL,4) = FUNIJK(PIJK(LL,1), PIJK(LL,2),PIJK(LL,3))

      END SUBROUTINE MPPIC_FIND_NEW_CELL
      END MODULE mppic_wallbc
