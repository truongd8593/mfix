!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: PIC_APPLY_WALLBC_STL                                    !
!  Author: R. Garg                                                     !
!                                                                      !
!  Purpose: Detect collisions between PIC particles and STLs.          !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE PIC_APPLY_WALLBC_STL

! Global Variables:
!---------------------------------------------------------------------//
! Particle postions: Current/Previous
      use discretelement, only: DES_POS_NEW, DES_POS_OLD
! Paricle velocities
      use discretelement, only: DES_VEL_NEW
! Particle radius
      use discretelement, only: DES_RADIUS
! Max number of particles on this process
      use discretelement, only: MAX_PIP
! The number of neighbor facets for each DES grid cell
      use discretelement, only: cellneighbor_facet_num
! Facet informatin binned to the DES grid
      use discretelement, only: cellneighbor_facet
! Map from particle to DES grid cell.
      use discretelement, only: DG_PIJK, DTSOLID
! Flag indicating that the index cell contains no STLs
      use stl, only: NO_NEIGHBORING_FACET_DES
! STL Vertices
      use stl, only: VERTEX
! STL Facet normals
      use stl, only: NORM_FACE
! The number of mass inflow/outflow BCs
      use pic_bc, only: PIC_BCMI, PIC_BCMO
! Minimum velocity to offset gravitational forces
      use pic_bc, only: minVEL, minVEL_MAG, OoMinVEL_MAG
! Solids time step size
      use discretelement, only: DTSOLID
! Gravitational force vector
      use discretelement, only: GRAV

      use discretelement, only: MAX_RADIUS

! Global Parameters:
!---------------------------------------------------------------------//
      use param1, only: ZERO, SMALL_NUMBER, ONE, UNDEFINED
      use functions, only: IS_NONEXISTENT

! Module Procedures:
!---------------------------------------------------------------------//
      use des_stl_functions

      use error_manager

      IMPLICIT NONE

! Local Variables:
!---------------------------------------------------------------------//
! Loop counters.
      INTEGER NF, NP, LC1, LC2
! STL normal vector
      DOUBLE PRECISION :: NORM_PLANE(3)
! line is parameterized as p = p_ref + t * dir_line, t is line_param
      DOUBLE PRECISION :: LINE_T

! distance parcel travelled out of domain along the facet normal
      DOUBLE PRECISION :: DIST(3)

      INTEGER :: AXIS
      DOUBLE PRECISION :: PARTICLE_MIN(3), PARTICLE_MAX(3)

      DOUBLE PRECISION :: DT_RBND, RADSQ, OFFSET

      CALL INIT_ERR_MSG("PIC_APPLY_WALLBC_STL")

! Minimum velocity needed to offset gravity.
      minVEL = -DTSOLID*GRAV(:)
      minVEL_MAG = dot_product(minVEL,minVEL)
      OoMinVEL_MAG = 1.0d0/minVEL_MAG

      DT_RBND = 1.01d0*DTSOLID

      OFFSET = 2.0*MAX_RADIUS

      DO NP = 1, MAX_PIP

! Skip non-existent particles
         IF(IS_NONEXISTENT(NP)) CYCLE

! If no neighboring facet in the surrounding 27 cells, then exit
         IF (NO_NEIGHBORING_FACET_DES(DG_PIJK(NP))) CYCLE

         RADSQ = DES_RADIUS(NP)*DES_RADIUS(NP)

! Loop through the STLs associated with the DES grid cell.
         CELLLOOP: DO LC1=1, CELLNEIGHBOR_FACET_NUM(DG_PIJK(NP))

! Get the index of the neighbor facet
            NF = CELLNEIGHBOR_FACET(DG_PIJK(NP))%P(LC1)

            NORM_PLANE = NORM_FACE(:,NF)

            IF(dot_product(NORM_PLANE, DES_VEL_NEW(:,NP)) > 0.d0)      &
               CYCLE CELLLOOP

            CALL INTERSECTLNPLANE(DES_POS_NEW(:,NP), DES_VEL_NEW(:,NP),&
                VERTEX(1,:,NF), NORM_FACE(:,NF), LINE_T)

            IF(abs(LINE_T) <= DT_RBND) THEN

! Avoid direct collisions with STLs that are not next to the parcel
               DO LC2=1,3
                  IF(CELLNEIGHBOR_FACET(DG_PIJK(NP))%extentMIN(LC1) >  &
                     DES_POS_NEW(LC2,NP)+OFFSET) CYCLE CELLLOOP
                  IF(CELLNEIGHBOR_FACET(DG_PIJK(NP))%extentMAX(LC1) <  &
                     DES_POS_NEW(LC2,NP)-OFFSET) CYCLE CELLLOOP
               ENDDO

! Correct the position of particles found too close to the STL.
               DIST = LINE_T*DES_VEL_NEW(:,NP)
               IF(dot_product(DIST,DIST) <= RADSQ .OR. LINE_T <= ZERO) &
                  DES_POS_NEW(:,NP) = DES_POS_NEW(:,NP) + DIST(:) +    &
                  DES_RADIUS(NP)*NORM_PLANE

               CALL PIC_REFLECT_PART(NP, NORM_PLANE(:))
            ENDIF

         ENDDO CELLLOOP
      END DO

! Seed new parcels entering the system.
      IF(PIC_BCMI > 0) CALL PIC_MI_BC
      IF(PIC_BCMO > 0) CALL PIC_MO_BC

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE PIC_APPLY_WALLBC_STL


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: PIC_REFLECT_PARTICLE                                  C
!  Purpose:                                                            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE PIC_REFLECT_PART(LL, WALL_NORM)

      USE discretelement, only : dimn, DES_VEL_NEW
      USE mfix_pic, only : MPPIC_COEFF_EN_WALL, MPPIC_COEFF_ET_WALL

! Minimum velocity to offset gravitational forces
      use pic_bc, only: minVEL, minVEL_MAG, OoMinVEL_MAG


      use param1, only: SMALL_NUMBER

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
      INTEGER, INTENT(IN) :: LL
      DOUBLE PRECISION, INTENT(IN) ::  WALL_NORM(DIMN)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      !magnitude of pre-collisional normal and tangential velocity components
      DOUBLE PRECISION :: VEL_NORMMAG_APP

      !pre collisional normal and tangential velocity components in vector format
      !APP ==> approach
      DOUBLE PRECISION :: VEL_NORM_APP(DIMN), VEL_TANG_APP(DIMN)


      !post collisional normal and tangential velocity components in vector format
      !SEP ==> separation
      DOUBLE PRECISION :: VEL_NORM_SEP(DIMN), VEL_TANG_SEP(DIMN)

      DOUBLE PRECISION :: COEFF_REST_EN, COEFF_REST_ET

! Minimum particle velocity with respect to gravity. This ensures that
! walls "push back" to offset gravitational forces.
      DOUBLE PRECISION :: projGRAV(3)

      INTEGER :: LC
 
!-----------------------------------------------





      VEL_NORMMAG_APP = DOT_PRODUCT(WALL_NORM(:), DES_VEL_NEW(:,LL))

! Currently assuming that wall is at rest. Needs improvement for moving wall

      VEL_NORM_APP(:) = VEL_NORMMAG_APP*WALL_NORM(:)
      VEL_TANG_APP(:) = DES_VEL_NEW(:,LL) - VEL_NORM_APP(:)


!post collisional velocities

      COEFF_REST_EN = MPPIC_COEFF_EN_WALL
      COEFF_REST_ET = MPPIC_COEFF_ET_WALL

      VEL_NORM_SEP(:) = -COEFF_REST_EN*VEL_NORM_APP(:)
      VEL_TANG_SEP(:) =  COEFF_REST_ET*VEL_TANG_APP(:)

      DES_VEL_NEW(:,LL) = VEL_NORM_SEP(:) + VEL_TANG_SEP(:)

      IF(dot_product(WALL_NORM, minVEL) > 0.0d0)THEN

! Projection of rebound velocity onto mininum velocity.
         projGRAV = dot_product(minVEL, DES_VEL_NEW(:,LL))*OoMinVEL_MAG
         projGRAV = minVEL*projGRAV

         IF(dot_product(projGRAV, projGRAV) < minVEL_MAG) THEN
            DO LC=1,3
               IF(minVEL(LC) > SMALL_NUMBER) THEN
                  DES_VEL_NEW(LC,LL)=max(DES_VEL_NEW(LC,LL),minVEL(LC))
               ELSEIF(minVEL(LC) < -SMALL_NUMBER) THEN
                  DES_VEL_NEW(LC,LL)=min(DES_VEL_NEW(LC,LL),minVEL(LC))
               ENDIF
            ENDDO

         ENDIF
      ENDIF

      RETURN
      END SUBROUTINE PIC_REFLECT_PART

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: Mass_OUTFLOW_PIC                                        !
!  Author: R. Garg                                   Date: 23-Jun-14   !
!                                                                      !
!  Purpose:  Routine to delete out of domain parcels for PIC           !
!  implementation                                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE PIC_MO_BC

      USE error_manager
      USE mpi_utility
      use bc
      use discretelement
      use functions
      use pic_bc

      implicit none

      INTEGER :: IJK
      INTEGER :: LC, LP, NP
      INTEGER :: BCV, BCV_I

      DOUBLE PRECISION :: DIST

      INTEGER :: PIP_DEL_COUNT
      INTEGER, dimension(0:numpes-1) :: del_count_all
      INTEGER :: IPROC


      CALL INIT_ERR_MSG("PIC_MO_BC")

      PIP_DEL_COUNT = 0
      DO BCV_I = 1, PIC_BCMO

         BCV = PIC_BCMO_MAP(BCV_I)

         SELECT CASE (BC_PLANE(BCV))
         END SELECT

         DO LC=PIC_BCMO_IJKSTART(BCV_I), PIC_BCMO_IJKEND(BCV_I)
            IJK = PIC_BCMO_IJK(LC)
            DO LP= 1,PINC(IJK)

               NP = PIC(IJK)%p(LP)

               SELECT CASE (BC_PLANE(BCV))
               CASE('S'); DIST = YN(BC_J_s(BCV)-1) - DES_POS_NEW(2,NP)
               CASE('N'); DIST = DES_POS_NEW(2,NP) - YN(BC_J_s(BCV))
               CASE('W'); DIST = XE(BC_I_w(BCV)-1) - DES_POS_NEW(1,NP)
               CASE('E'); DIST = DES_POS_NEW(1,NP) - XE(BC_I_w(BCV))
               CASE('B'); DIST = ZT(BC_K_b(BCV)-1) - DES_POS_NEW(3,NP)
               CASE('T'); DIST = DES_POS_NEW(3,NP) - ZT(BC_K_b(BCV))
               END SELECT

               IF(DIST .lt. zero ) THEN

! The parcel has left this bc_plane

                  CALL DELETE_PARCEL(NP)
                  PIP_DEL_COUNT = PIP_DEL_COUNT+1
               ENDIF

            ENDDO
         ENDDO
      ENDDO


      IF(PIC_REPORT_DELETION_STATS) then

         del_count_all(:) = 0
         del_count_all(mype) = pip_del_count
         call global_all_sum(del_count_all(0:numpes-1))

         IF(SUM(del_COUNT_ALL(:)).GT.0) THEN
            WRITE(err_msg,'(/,2x,A,2x,i10)')&
            'TOTAL NUMBER OF PARCELS DELETED GLOBALLY = ', SUM(DEL_COUNT_ALL(:))

            call flush_err_msg(header = .false., footer = .false.)

            DO IPROC = 0, NUMPES-1
               IF(DMP_LOG) WRITE(UNIT_LOG, '(/,A,i4,2x,A,2x,i5)') &
                    'PARCELS ADDED ON PROC:', IPROC,' EQUAL TO', pip_del_count
            ENDDO
         ENDIF

      ENDIF


      CALL FINL_ERR_MSG
      RETURN
      END SUBROUTINE PIC_MO_BC




!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: DELETE_PARCEL                                           !
!  Author: R. Garg                                    Date: 23-Jun-14  !
!                                                                      !
!  Purpose:  Routine to delete parcel                                  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DELETE_PARCEL(NP)

      USE compar
      USE constant
      USE des_bc
      USE discretelement
      USE funits
      USE geometry
      USE indices
      USE param1
      USE physprop
      USE mfix_pic
      USE functions

      IMPLICIT NONE


      INTEGER, INTENT(IN) :: NP

      CALL SET_NORMAL(NP)

      DES_POS_OLD(:,NP) = ZERO
      DES_POS_NEW(:,NP) = ZERO
      DES_VEL_OLD(:,NP) = ZERO
      DES_VEL_NEW(:,NP) = ZERO
      OMEGA_OLD(:,NP) = ZERO
      OMEGA_NEW(:,NP) = ZERO
      DES_RADIUS(NP) = ZERO
      PMASS(NP) = ZERO
      PVOL(NP) = ZERO
      RO_Sol(NP) = ZERO
      OMOI(NP) = ZERO

      DES_STAT_WT(NP) = ZERO

      FC(:,NP) = ZERO

      PIP = PIP - 1

      RETURN
      END SUBROUTINE DELETE_PARCEL

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: MPPIC_MI_BC                                             C
!  Purpose:                                                            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE  PIC_MI_BC

      USE run
      USE discretelement
      USE geometry
      USE constant
      USE cutcell
      USE mfix_pic
      USE pic_bc
      USE bc
      USE mpi_utility
      USE randomno
      use error_manager
      USE functions

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER :: WDIR, IDIM, IPCOUNT, LAST_EMPTY_SPOT, NEW_SPOT
      INTEGER :: BCV, BCV_I, L, LC, PIP_ADD_COUNT, IPROC
      INTEGER ::  IFLUID, JFLUID, KFLUID, IJK, M
      DOUBLE PRECISION :: CORD_START(3), DOML(3),WALL_NORM(3)
      DOUBLE PRECISION :: AREA_INFLOW, VEL_INFLOW(DIMN), STAT_WT

      LOGICAL :: DELETE_PART
      DOUBLE PRECISION :: EPS_INFLOW(DES_MMAX)
      DOUBLE PRECISION REAL_PARTS(DES_MMAX), COMP_PARTS(DES_MMAX), VEL_NORM_MAG, VOL_INFLOW, VOL_IJK


! Temp logical variables for checking constant npc and statwt specification
      LOGICAL :: CONST_NPC, CONST_STATWT

      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: RANDPOS
      INTEGER :: CNP_CELL_COUNT
      DOUBLE PRECISION :: RNP_CELL_COUNT
      integer :: lglobal_id
      integer, dimension(0:numpes-1) :: add_count_all

      LOGICAL, parameter :: setDBG = .true.
      LOGICAL :: dFlag

      type :: ty_spotlist
         integer spot
         type(ty_spotlist),pointer :: next => NULL()
      end type ty_spotlist

      type(ty_spotlist),pointer :: &
           cur_spotlist => NULL(), &
           prev_spotlist => NULL(), &
           temp_spotlist => NULL()

!-----------------------------------------------

      CALL INIT_ERR_MSG("PIC_MI_BC")

      dFlag = (DMP_LOG .AND. setDBG)
      PIP_ADD_COUNT = 0
      LAST_EMPTY_SPOT = 0

      ALLOCATE(cur_spotlist); cur_spotlist%spot = -1

      DO BCV_I = 1, PIC_BCMI
         BCV = PIC_BCMI_MAP(BCV_I)

         WALL_NORM(1:3) =  PIC_BCMI_NORMDIR(BCV_I,1:3)

         !Find the direction of the normal for this wall
         WDIR = 0
         DO IDIM = 1, DIMN
            WDIR = WDIR + ABS(WALL_NORM(IDIM))*IDIM
         end DO

         DO LC=PIC_BCMI_IJKSTART(BCV_I), PIC_BCMI_IJKEND(BCV_I)
            IJK = PIC_BCMI_IJK(LC)

            IF(.NOT.FLUID_AT(IJK)) CYCLE

            IFLUID = I_OF(IJK)
            JFLUID = J_OF(IJK)
            KFLUID = K_OF(IJK)

            CORD_START(1) = XE(IFLUID) - PIC_BCMI_OFFSET (BCV_I,1)*DX(IFLUID)

            CORD_START(2) = YN(JFLUID) - PIC_BCMI_OFFSET (BCV_I,2)*DY(JFLUID)


            CORD_START(3) = merge(zero, ZT(KFLUID) - PIC_BCMI_OFFSET (BCV_I,3)*DZ(KFLUID), no_k)

            DOML(1) = DX(IFLUID)
            DOML(2) = DY(JFLUID)
            DOML(3) = MERGE(DZ(1), DZ(KFLUID), NO_K)

            AREA_INFLOW = DOML(1)*DOML(2)*DOML(3)/DOML(WDIR)

            VOL_IJK = DOML(1)*DOML(2)*DOML(3)

            DOML(WDIR) = ZERO
            !set this to zero as the particles will
            !be seeded only on the BC plane

            DO M = 1, DES_MMAX

               IF(SOLIDS_MODEL(M) /= 'PIC') CYCLE
               EPS_INFLOW(M) = BC_ROP_S(BCV, M)/DES_RO_S(M)
               VEL_INFLOW(1) = BC_U_S(BCV, M)
               VEL_INFLOW(2) = BC_V_S(BCV, M)
               VEL_INFLOW(3) = BC_W_S(BCV, M)

               VEL_NORM_MAG = ABS(DOT_PRODUCT(VEL_INFLOW(1:DIMN), WALL_NORM(1:DIMN)))
               VOL_INFLOW   = AREA_INFLOW*VEL_NORM_MAG*DTSOLID

               REAL_PARTS(M) = 6.d0*EPS_INFLOW(M)*VOL_INFLOW/(PI*(DES_D_p0(M)**3.d0))
               COMP_PARTS(M) = zero

               CONST_NPC    = (BC_PIC_MI_CONST_NPC   (BCV, M) .ne. 0)
               CONST_STATWT = (BC_PIC_MI_CONST_STATWT(BCV, M) .ne. ZERO)
               IF(CONST_NPC) THEN
                  IF(EPS_INFLOW(M).GT.ZERO) &
                  COMP_PARTS(M) = REAL(BC_PIC_MI_CONST_NPC(BCV, M))* &
                  VOL_INFLOW/VOL_IJK
               ELSEIF(CONST_STATWT) THEN
                  COMP_PARTS(M) = REAL_PARTS(M)/ &
                  BC_PIC_MI_CONST_STATWT(BCV, M)
               ENDIF

               pic_bcmi_rnp(LC,M) = pic_bcmi_rnp(LC,M) + REAL_PARTS(M)

               pic_bcmi_cnp(LC,M) = pic_bcmi_cnp(LC,M) + COMP_PARTS(M)


               IF(pic_bcmi_cnp(LC,M).GE.1.d0) THEN
                  CNP_CELL_COUNT = INT(pic_bcmi_cnp(LC,M))
                  pic_bcmi_cnp(LC,M) = pic_bcmi_cnp(LC,M) - CNP_CELL_COUNT

                  RNP_CELL_COUNT = pic_bcmi_rnp(LC,M)
                  pic_bcmi_rnp(LC,M)  = zero

                  !set pic_bcmi_rnp to zero to reflect that all real particles have been seeded
                  STAT_WT = RNP_CELL_COUNT/REAL(CNP_CELL_COUNT)
                  ALLOCATE(RANDPOS(CNP_CELL_COUNT*DIMN))
                  CALL UNI_RNO(RANDPOS(:))

                  DO IPCOUNT = 1, CNP_CELL_COUNT


                     CALL PIC_FIND_EMPTY_SPOT(LAST_EMPTY_SPOT, NEW_SPOT)

                     DES_POS_OLD(1:DIMN, NEW_SPOT) =  CORD_START(1:DIMN) &
                          & + RANDPOS((IPCOUNT-1)*DIMN+1: &
                          & (IPCOUNT-1)*DIMN+DIMN)*DOML(1:DIMN)
                     DES_POS_NEW(:, NEW_SPOT) = DES_POS_OLD(:,NEW_SPOT)
                     DES_VEL_OLD(1:DIMN,NEW_SPOT) = VEL_INFLOW(1:DIMN)

                     DES_VEL_NEW(:, NEW_SPOT) = DES_VEL_OLD(:, NEW_SPOT)

                     DES_RADIUS(NEW_SPOT) = DES_D_p0(M)*HALF

                     RO_Sol(NEW_SPOT) =  DES_RO_S(M)

                     DES_STAT_WT(NEW_SPOT) = STAT_WT

                     PIJK(NEW_SPOT, 1) = IFLUID
                     PIJK(NEW_SPOT, 2) = JFLUID
                     PIJK(NEW_SPOT, 3) = KFLUID
                     PIJK(NEW_SPOT, 4) = IJK
                     PIJK(NEW_SPOT, 5) = M

                     PVOL(NEW_SPOT) = (4.0d0/3.0d0)*Pi*DES_RADIUS(NEW_SPOT)**3
                     PMASS(NEW_SPOT) = PVOL(NEW_SPOT)*RO_SOL(NEW_SPOT)

                     FC(:, NEW_SPOT) = zero
                     DELETE_PART = .false.
                     IF(PIC_BCMI_INCL_CUTCELL(BCV_I)) &
                          CALL CHECK_IF_PARCEL_OVERLAPS_STL &
                          (des_pos_new(1:dimn, NEW_SPOT), &
                          DELETE_PART)

                     IF(.NOT.DELETE_PART) THEN

                        PIP = PIP+1
                        PIP_ADD_COUNT = PIP_ADD_COUNT + 1
                        CALL SET_NORMAL(NEW_SPOT)
                        ! add to the list
                        ALLOCATE(temp_spotlist)
                        temp_spotlist%spot = new_spot
                        temp_spotlist%next => cur_spotlist
                        cur_spotlist => temp_spotlist
                        nullify(temp_spotlist)

                     ELSE
                        CALL SET_NONEXISTENT(NEW_SPOT)
                        LAST_EMPTY_SPOT = NEW_SPOT - 1
                     ENDIF

                     !WRITE(*,'(A,2(2x,i5), 2x, A,2x, 3(2x,i2),2x, A, 3(2x,g17.8))') &
                     !   'NEW PART AT ', NEW_SPOT, MAX_PIP, 'I, J, K = ', IFLUID, JFLUID, KFLUID, 'POS =', DES_POS_NEW(:,NEW_SPOT)
                     !IF(DMP_LOG) WRITE(UNIT_LOG,'(A,2x,i5, 2x, A,2x, 3(2x,i2),2x, A, 3(2x,g17.8))') &
                     !    'NEW PART AT ', NEW_SPOT, 'I, J, K = ', IFLUID, JFLUID, KFLUID, 'POS =', DES_POS_NEW(:,NEW_SPOT)

                     !WRITE(*,*) 'WDIR, DOML = ', WDIR, DOML(:)
                  END DO
                  DEALLOCATE(RANDPOS)

               end IF
            end DO
         end DO
      end DO

!Now assign global id to new particles added
      add_count_all(:) = 0
      add_count_all(mype) = pip_add_count
      call global_all_sum(add_count_all(0:numpes-1))
      lglobal_id = imax_global_id + sum(add_count_all(0:mype-1))

      do l = 1,pip_add_count
         lglobal_id = lglobal_id + 1
         iglobal_id(cur_spotlist%spot)= lglobal_id
         prev_spotlist=> cur_spotlist
         cur_spotlist => cur_spotlist%next
         deallocate(prev_spotlist)
      end do
      deallocate(cur_spotlist)
      imax_global_id = imax_global_id + sum(add_count_all(0:numpes-1))

      IF(PIC_REPORT_SEEDING_STATS) then

         IF(SUM(ADD_COUNT_ALL(:)).GT.0) THEN
            WRITE(err_msg,'(/,2x,A,2x,i10)') &
               'TOTAL NUMBER OF PARCELS ADDED GLOBALLY = ',&
                SUM(ADD_COUNT_ALL(:))

            call flush_err_msg(header = .false., footer = .false.)

            DO IPROC = 0, NUMPES-1
               IF(DMP_LOG) WRITE(UNIT_LOG, '(/,A,i4,2x,A,2x,i5)') &
                  'PARCELS ADDED ON PROC:', IPROC,&
                  ' EQUAL TO', ADD_COUNT_ALL(IPROC)
            ENDDO
         ENDIF

      ENDIF

      CALL FINL_ERR_MSG
      RETURN
      END SUBROUTINE PIC_MI_BC



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: PIC_FIND_EMPTY_SPOT                                   C
!  Purpose:                                                            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE PIC_FIND_EMPTY_SPOT(LAST_INDEX, EMPTY_SPOT)
      USE funits
      USE mpi_utility
      USE error_manager
      USE discretelement, only: max_pip
      USE functions, only: is_nonexistent
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
      INTEGER, INTENT(INOUT) :: LAST_INDEX
      INTEGER, INTENT(OUT) :: EMPTY_SPOT
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      LOGICAL :: SPOT_FOUND
      INTEGER :: LL
!-----------------------------------------------

      CALL INIT_ERR_MSG("PIC_FIND_EMPTY_SPOT")

      IF(LAST_INDEX.EQ.MAX_PIP) THEN
         WRITE(ERR_MSG,2001)

         CALL FLUSH_ERR_MSG(abort = .true.)
         call mfix_exit(mype)
      ENDIF
      SPOT_FOUND = .false.

      DO LL = LAST_INDEX+1, MAX_PIP

         if(IS_NONEXISTENT(LL)) THEN
            EMPTY_SPOT = LL
            LAST_INDEX = LL
            SPOT_FOUND = .true.
            EXIT
         ENDIF
      ENDDO

      IF(.NOT.SPOT_FOUND) THEN
         WRITE(ERR_MSG,2002)
         CALL FLUSH_ERR_MSG(abort = .true.)
      ENDIF

 2001 FORMAT(/,5X,  &
      & 'ERROR IN PIC_FIND_EMPTY_SPOT', /5X, &
      & 'NO MORE EMPTY SPOT IN THE PARTICLE ARRAY TO ADD A NEW PARTICLE',/5X, &
      & 'TERMINAL ERROR: STOPPING')

 2002 FORMAT(/,5X,  &
      & 'ERROR IN PIC_FIND_EMPTY_SPOT', /5X, &
      & 'COULD NOT FIND A SPOT FOR ADDING NEW PARTICLE',/5X, &
      & 'INCREASE THE SIZE OF THE INITIAL ARRAYS', 5X, &
      & 'TERMINAL ERROR: STOPPING')

      CALL FINL_ERR_MSG
      END SUBROUTINE PIC_FIND_EMPTY_SPOT





!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: PIC_FIND_NEW_CELL                                     C
!  Purpose:                                                            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE PIC_FIND_NEW_CELL(LL)
      USE discretelement
      use mpi_utility
      USE cutcell
      USE functions

      IMPLICIT NONE

!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
      INTEGER, INTENT(IN) :: LL
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      DOUBLE PRECISION :: XPOS, YPOS, ZPOS
      INTEGER :: I, J, K
!-----------------------------------------------

      I = PIJK(LL,1)
      J = PIJK(LL,2)
      K = PIJK(LL,3)
      XPOS = DES_POS_NEW(1,LL)
      YPOS = DES_POS_NEW(2,LL)
      IF (DIMN .EQ. 3) THEN
         ZPOS = DES_POS_NEW(3,LL)
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

      END SUBROUTINE PIC_FIND_NEW_CELL

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  SUBROUTINE: CHECK_IF_PARCEL_OVERLAPS_STL                            C
!  Authors: Rahul Garg                               Date: 21-Mar-2014 C
!                                                                      C
!  Purpose: This subroutine is special written to check if a particle  C
!          overlaps any of the STL faces. The routine exits on         C
!          detecting an overlap. It is called after initial            C
!          generation of lattice configuration to remove out of domain C
!          particles                                                   C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CHECK_IF_PARCEL_OVERLAPS_STL(POS, OVERLAP_EXISTS)

      USE discretelement, only: MAX_RADIUS
      USE geometry, only: do_k

      USE des_stl_functions
      USE desgrid
      USE stl

      Implicit none

      DOUBLE PRECISION, INTENT(IN) :: POS(3)
      LOGICAL, INTENT(OUT) :: OVERLAP_EXISTS

      INTEGER I, J, K, IJK, NF, LC

! line is parameterized as p = p_ref + t * dir_line, t is line_param
      DOUBLE PRECISION :: LINE_T
      DOUBLE PRECISION :: DIST(3), RADSQ


      K = 1
      IF(DO_K) K = min(DG_KEND2,max(DG_KSTART2,KOFPOS(POS(3))))
      J = min(DG_JEND2,max(DG_JSTART2,JOFPOS(POS(2))))
      I = min(DG_IEND2,max(DG_ISTART2,IOFPOS(POS(1))))

      IJK = DG_FUNIJK(I,J,K)
      IF (NO_NEIGHBORING_FACET_DES(IJK)) RETURN

      OVERLAP_EXISTS = .TRUE.

! Pad the max radius to keep parcels sufficiently away from walls.
      RADSQ = (1.1d0*MAX_RADIUS)**2

! Parametrize a line as p = p_0 + t normal and intersect with the
! triangular plane.

! If t>0, then point is on the non-fluid side of the plane, if the
! plane normal is assumed to point toward the fluid side.
      DO LC = 1, LIST_FACET_AT_DES(IJK)%COUNT_FACETS
         NF = LIST_FACET_AT_DES(IJK)%FACET_LIST(LC)

         CALL INTERSECTLNPLANE(POS, NORM_FACE(:,NF), &
            VERTEX(1,:,NF), NORM_FACE(:,NF), LINE_T)

! Orthogonal projection puts the point outside of the domain or less than
! one particle radius to the facet.
         DIST = LINE_T*NORM_FACE(:,NF)
         IF(LINE_T > ZERO .OR. dot_product(DIST,DIST) <= RADSQ) RETURN

      ENDDO

      OVERLAP_EXISTS = .FALSE.

      RETURN
      END SUBROUTINE CHECK_IF_PARCEL_OVERLAPS_STL


!----------------------------------------------------------------------!
!                                                                      !
!                                                                      !
!                                                                      !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE write_this_facet_and_parcel(FID, position, velocity)
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
      double precision, intent(in), dimension(dimn) :: position, velocity
      Integer, intent(in) :: fid
      Integer :: stl_unit, vtp_unit , k
      CHARACTER(LEN=100) :: stl_fname, vtp_fname
      real :: temp_array(3)

      stl_unit = 1001
      vtp_unit = 1002

      WRITE(vtp_fname,'(A,"_OFFENDING_PARTICLE",".vtp")') TRIM(RUN_NAME)
      WRITE(stl_fname,'(A,"_STL_FACE",".stl")') TRIM(RUN_NAME)

      open(vtp_unit, file = vtp_fname, form='formatted',convert='big_endian')
      open(stl_unit, file = stl_fname, form='formatted',convert='big_endian')

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
      write (vtp_unit,"(15x,es13.6)") (1.d0)
      write(vtp_unit,"(12x,a)") '</DataArray>'

      temp_array = zero
      temp_array(1:DIMN) = velocity(1:dimn)
      write(vtp_unit,"(12x,a,a)") '<DataArray type="Float32" ',&
           'Name="Velocity" NumberOfComponents="3" format="ascii">'
      write (vtp_unit,"(15x,3(es13.6,3x))")&
           ((temp_array(k)),k=1,3)
      write(vtp_unit,"(12x,a,/9x,a)") '</DataArray>','</PointData>'
      ! skip cell data
      write(vtp_unit,"(9x,a)") '<CellData></CellData>'

      temp_array = zero
      temp_array(1:dimn) = position(1:dimn)
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

      write(stl_unit,*) '   facet normal ', NORM_FACE(1:3,FID)
      write(stl_unit,*) '      outer loop'
      write(stl_unit,*) '         vertex ', VERTEX(1,1:3,FID)
      write(stl_unit,*) '         vertex ', VERTEX(2,1:3,FID)
      write(stl_unit,*) '         vertex ', VERTEX(3,1:3,FID)
      write(stl_unit,*) '      endloop'
      write(stl_unit,*) '   endfacet'

      write(stl_unit,*)'endsolid vcg'

      close(vtp_unit, status = 'keep')
      close(stl_unit, status = 'keep')
      write(*,*) 'wrote a facet and a parcel. now waiting'
      read(*,*)
      end SUBROUTINE write_this_facet_and_parcel

