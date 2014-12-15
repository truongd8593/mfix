!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  SUBROUTINE: CHECK_IF_PARCEL_OVELAPS_STL                             C
!                                                                      C
!  Purpose: This subroutine is special written to check if a particle  C
!          overlaps any of the STL faces. The routine exits on         C
!          detecting an overlap. It is called after initial            C
!          generation of lattice configuration to remove out of domain C
!          particles                                                   C
!                                                                      C
!  Authors: Rahul Garg                               Date: 21-Mar-2014 C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CHECK_IF_PARCEL_OVELAPS_STL(POSITION, PCELL, &
      OVERLAP_EXISTS)

      USE run
      USE param1
      USE discretelement, only: dimn
      USE geometry
      USE constant
      USE cutcell
      USE indices
      USE stl
      USE des_stl_functions
      USE mpi_utility
      USE functions
      Implicit none

      DOUBLE PRECISION, INTENT(IN) :: POSITION(DIMN)
      INTEGER , INTENT(IN) :: PCELL(4)
      LOGICAL, INTENT(OUT) :: OVERLAP_EXISTS

      INTEGER I, J, K, IJK, NF, FOCUS_PARTICLE

      DOUBLE PRECISION :: RADSQ, DISTSQ, DIST(DIMN), CLOSEST_PT(DIMN)
      INTEGER :: COUNT_FAC, COUNT, contact_facet_count, NEIGH_CELLS, &
      NEIGH_CELLS_NONNAT, &
      LIST_OF_CELLS(27), CELL_ID, I_CELL, J_CELL, K_CELL, cell_count , &
      IMINUS1, IPLUS1, JMINUS1, JPLUS1, KMINUS1, KPLUS1, PHASELL, LOC_MIN_PIP, &
      LOC_MAX_PIP

      double precision :: velocity(dimn)
      !reference point and direction of the line
      double precision, dimension(dimn) :: ref_line,  dir_line
      !reference point and normal of the plane
      double precision, dimension(dimn) :: ref_plane, norm_plane
      !line is parameterized as p = p_ref + t * dir_line, t is line_param
      double precision :: line_t
      !logical for determining if a point is on the triangle or not
      logical :: ontriangle
      !The line and plane intersection point
      double precision, dimension(dimn) :: point_onplane

      FOCUS_PARTICLE = -1

      OVERLAP_EXISTS = .false.

      IF (NO_NEIGHBORING_FACET_DES(PCELL(4))) RETURN

      LIST_OF_CELLS(:) = -1
      NEIGH_CELLS = 0
      NEIGH_CELLS_NONNAT  = 0
      CELL_ID = PCELL(4)
      COUNT_FAC = LIST_FACET_AT_DES(CELL_ID)%COUNT_FACETS


      IF (COUNT_FAC.gt.0)   then
         !first add the facets in the cell the particle currently resides in
         NEIGH_CELLS = NEIGH_CELLS + 1
         LIST_OF_CELLS(NEIGH_CELLS) = CELL_ID
      ENDIF

      I_CELL = PCELL(1)
      J_CELL = PCELL(2)
      K_CELL = PCELL(3)

      IPLUS1  =  MIN( I_CELL + 1, IEND2)
      IMINUS1 =  MAX( I_CELL - 1, ISTART2)

      JPLUS1  =  MIN (J_CELL + 1, JEND2)
      JMINUS1 =  MAX( J_CELL - 1, JSTART2)

      KPLUS1  =  MIN (K_CELL + 1, KEND2)
      KMINUS1 =  MAX( K_CELL - 1, KSTART2)

      DO K = KMINUS1, KPLUS1
         DO J = JMINUS1, JPLUS1
            DO I = IMINUS1, IPLUS1
               IJK = FUNIJK(I,J,K)
               COUNT_FAC = LIST_FACET_AT_DES(IJK)%COUNT_FACETS
               IF(COUNT_FAC.EQ.0) CYCLE
               NEIGH_CELLS_NONNAT = NEIGH_CELLS_NONNAT + 1
               NEIGH_CELLS = NEIGH_CELLS + 1
               LIST_OF_CELLS(NEIGH_CELLS) = IJK
               !WRITE(*,'(A10, 4(2x,i5))') 'WCELL  = ', IJK, I,J,K
            ENDDO
         ENDDO
      ENDDO

      CONTACT_FACET_COUNT = 0

      DO CELL_COUNT = 1, NEIGH_CELLS
         IJK = LIST_OF_CELLS(CELL_COUNT)

         DO COUNT = 1, LIST_FACET_AT_DES(IJK)%COUNT_FACETS
            NF = LIST_FACET_AT_DES(IJK)%FACET_LIST(COUNT)
            line_t  = -Undefined
            !-undefined, because non zero values will imply intersection
            !with the plane
            ontriangle = .false.

            !parametrize a line as p = p_0 + t normal
            !and intersect with the triangular plane.
            !if t>0, then point is on the
            !non-fluid side of the plane, if the plane normal
            !is assumed to point toward the fluid side
            ref_line(1:dimn) = position(1:dimn)
            dir_line(1:dimn) = NORM_FACE(1:dimn,NF)
            !Since this is for checking static config, line's direction
            !is the same as plane's normal. For moving particles,
            !the line's normal will be along the point joining new
            !and old positions.

            norm_plane(1:dimn) = NORM_FACE(1:dimn,NF)
            ref_plane(1:dimn)  = VERTEX(1, 1:dimn,NF)
            CALL intersectLnPlane(ref_line, dir_line, ref_plane, &
            norm_plane, line_t)
            if(line_t.gt.zero) then
               !this implies by orthogonal projection
               !that the point is on non-fluid side of the
               !facet
               point_onplane(1:dimn) = ref_line(1:dimn) + &
               line_t*dir_line(1:dimn)
               !Now check through barycentric coordinates if
               !this orthogonal projection lies on the facet or not
               !If it does, then this point is deemed to be on the
               !non-fluid side of the facet
               CALL checkPTonTriangle(point_onplane(1:dimn), &
               VERTEX(1,:,NF), VERTEX(2,:,NF), VERTEX(3,:,NF), &
               ontriangle)

               if(ontriangle) then
                  OVERLAP_EXISTS = .true.
                  !velocity = zero
                  !write(*,*) 'over lap detected', line_t
                  !call write_this_facet_and_parcel(NF, position, velocity)
                  RETURN
               endif
            endif

         ENDDO

      end DO

      RETURN

      END SUBROUTINE CHECK_IF_PARCEL_OVELAPS_STL


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

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: PIC_APPLY_WALLBC_STL                                  C
!  Purpose:                                                            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE PIC_APPLY_WALLBC_STL

      USE run
      USE param1
      USE discretelement
      USE geometry
      USE constant
      USE cutcell
      USE indices
      USE stl
      USE des_stl_functions
      USE mpi_utility
      use error_manager
      USE pic_bc
      USE functions
      use desgrid 

      Implicit none

      INTEGER I, J, K, IJK, NF, LL

      DOUBLE PRECISION :: RADSQ, DISTSQ, DIST(DIMN)
      INTEGER :: COUNT_FAC, COUNT, contact_facet_count, NEIGH_CELLS, &
      NEIGH_CELLS_NONNAT, LIST_OF_CELLS(27), IPROC, COUNT2, &
      CELL_ID, I_CELL, J_CELL, K_CELL, cell_count , &
      IMINUS1, IPLUS1, JMINUS1, JPLUS1, KMINUS1, KPLUS1

      !reference point and direction of the line
      double precision, dimension(dimn) :: ref_line,  dir_line
      !reference point and normal of the plane
      double precision, dimension(dimn) :: ref_plane, norm_plane
      !line is parameterized as p = p_ref + t * dir_line, t is line_param
      double precision :: line_t
      !logical for determining if a point is on the triangle or not
      logical :: ontriangle
      !The line and plane intersection point
      double precision, dimension(dimn) :: point_onplane
      !Old and new position of the parcel
      double precision, dimension(dimn) :: pt_old, pt_new

      !distance parcel travelled out of domain along the facet normal
      double precision :: dist_from_facet
      INTEGER, Parameter :: MAX_FACET_CONTS = 500
      INTEGER :: list_of_checked_facets(max_facet_conts)
!      INTEGER :: FOCUS_CELLID

      double precision :: veldotnorm

      CALL INIT_ERR_MSG("PIC_APPLY_WALLBC_STL")

!      FOCUS_CELLID = funijk(14,24,2)

      DO LL = 1, MAX_PIP

         ! skipping non-existent particles or ghost particles
         IF(.NOT.PEA(LL,1) .OR. PEA(LL,4)) CYCLE

         IF(LL.EQ.FOCUS_PARTICLE) DEBUG_DES = .TRUE.


         ! If no neighboring facet in the surrounding 27 cells, then exit
         IF (NO_NEIGHBORING_FACET_DES(DG_PIJK(LL))) cycle

         LIST_OF_CELLS(:) = -1
         NEIGH_CELLS = 0
         NEIGH_CELLS_NONNAT  = 0

         CELL_ID = DG_PIJK(LL)

         COUNT_FAC = LIST_FACET_AT_DES(CELL_ID)%COUNT_FACETS
         RADSQ = DES_RADIUS(LL)*DES_RADIUS(LL)

         PT_OLD = DES_POS_OLD(:,LL)
         PT_NEW = DES_POS_NEW(:,LL)

! First add the facets in the cell the particle currently resides in
         IF (COUNT_FAC.gt.01)   then
            NEIGH_CELLS = NEIGH_CELLS + 1
            LIST_OF_CELLS(NEIGH_CELLS) = CELL_ID
         ENDIF

         I_CELL = IOFPOS(PT_NEW(1))
         J_CELL = JOFPOS(PT_NEW(2))
         K_CELL = KOFPOS(PT_NEW(3))

         IPLUS1  =  MIN (I_CELL+1, DG_IEND2)
         IMINUS1 =  MAX (I_CELL-1, DG_ISTART2)

         JPLUS1  =  MIN (J_CELL+1, DG_JEND2)
         JMINUS1 =  MAX (J_CELL-1, DG_JSTART2)

         KPLUS1  =  MIN (K_CELL+1, DG_KEND2)
         KMINUS1 =  MAX (K_CELL-1, DG_KSTART2)

         DO K = KMINUS1, KPLUS1
         DO J = JMINUS1, JPLUS1
         DO I = IMINUS1, IPLUS1
            IJK = DG_FUNIJK(I,J,K)
            COUNT_FAC = LIST_FACET_AT_DES(IJK)%COUNT_FACETS
            IF(COUNT_FAC.EQ.0.or.ijk.eq.cell_id) CYCLE
            NEIGH_CELLS_NONNAT = NEIGH_CELLS_NONNAT + 1
            NEIGH_CELLS = NEIGH_CELLS + 1
            LIST_OF_CELLS(NEIGH_CELLS) = IJK
         ENDDO
         ENDDO
         ENDDO

         CONTACT_FACET_COUNT = 0

         CELLLOOP: DO CELL_COUNT = 1, NEIGH_CELLS

            IJK = LIST_OF_CELLS(CELL_COUNT)

            DO COUNT = 1, LIST_FACET_AT_DES(IJK)%COUNT_FACETS
               NF = LIST_FACET_AT_DES(IJK)%FACET_LIST(COUNT)

!only check for facets flagged as normal
               IF(STL_FACET_TYPE(NF).ne.facet_type_normal) cycle

! Neighboring cells will share facets with same facet ID so it is
! important to make sure a facet is checked (for speed) and accounted
! (for accuracy) only once.
               IF (ANY(NF.EQ.LIST_OF_CHECKED_FACETS                    &
                  (1:CONTACT_FACET_COUNT))) CYCLE

               CONTACT_FACET_COUNT = CONTACT_FACET_COUNT + 1
               LIST_OF_CHECKED_FACETS(CONTACT_FACET_COUNT) = NF

! Set to -UNDEFINED, because non zero values imply intersection.
               !with the plane
               LINE_T  = -UNDEFINED

               ONTRIANGLE = .FALSE.

! Parametrize a line as p = p_0 + t (p_1 - p_0) and intersect with the 
! triangular plane. If 0<t<1, then this ray connect old and new
! positions intersects with the facet.
               REF_LINE(1:DIMN) = PT_OLD(1:DIMN)
               DIR_LINE(1:DIMN) = PT_NEW(1:DIMN) - PT_OLD(1:DIMN)

               NORM_PLANE(1:DIMN) = NORM_FACE(1:DIMN,NF)
               REF_PLANE(1:DIMN)  = VERTEX(1, 1:DIMN,NF)

               VELDOTNORM = DOT_PRODUCT(NORM_PLANE(1:DIMN),            &
                    DES_VEL_NEW(:,LL))

! If the normal velocity of parcel is in the same direction as facet
! normal, then the parcel could not have crossed this facet.
               IF(VELDOTNORM.GT.0.d0) CYCLE

               CALL INTERSECTLNPLANE(REF_LINE, DIR_LINE, REF_PLANE,    &
                    NORM_PLANE, LINE_T)

! line_t = one implies that the particle is sitting on the stl face.
               IF(LINE_T.GE.ZERO.AND.LINE_T.LT.ONE) THEN

                  POINT_ONPLANE = REF_LINE + LINE_T*DIR_LINE

! Check through barycentric coordinates to determine if the orthogonal
! projection lies on the facet or not. If it does, then this point is
! deemed to be on the non-fluid side of the facet.
                  CALL CHECKPTONTRIANGLE(POINT_ONPLANE(:),             &
                       VERTEX(1,:,NF), VERTEX(2,:,NF), VERTEX(3,:,NF), &
                       ONTRIANGLE)

                  IF(ONTRIANGLE) THEN

                     DIST_FROM_FACET = ABS(DOT_PRODUCT(PT_NEW(:) -     &
                        POINT_ONPLANE(:), NORM_PLANE(:)))

                     DES_POS_NEW(:,LL) = POINT_ONPLANE(1:DIMN) + &
                        DIST_FROM_FACET*(NORM_PLANE(1:DIMN))

! In the reflect routine, it is assumed that the normal points into the
! system which is consitent with the definition of normal for facets.
                     CALL PIC_REFLECT_PART(LL, NORM_PLANE(:))

                     EXIT CELLLOOP
                  ENDIF
               ENDIF

            ENDDO

         END DO CELLLOOP

      END DO


      ! Seed new parcels entering the system.
      IF(PIC_BCMI > 0) CALL PIC_MI_BC
      IF(PIC_BCMO > 0) CALL PIC_MO_BC

      CALL FINL_ERR_MSG
      RETURN

      END SUBROUTINE PIC_APPLY_WALLBC_STL



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

      use discretelement
      use pic_bc
      use bc
      USE error_manager
      USE mpi_utility


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

               IF(PEA(NP,4)) CYCLE

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
            WRITE(err_msg,'(/,2x,A,2x,i10)') 'TOTAL NUMBER OF PARCELS DELETED GLOBALLY = ', SUM(DEL_COUNT_ALL(:))

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

!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER I,J,K,IJK       ! Unused indices at this time
      INTEGER PC              ! Loop counter
      INTEGER NPT             ! Temp Particle Index
! Tmp holder of particle position
      DOUBLE PRECISION XPOS, YPOS, ZPOS

! Loop indices used for clearing forces associated with exiting particles
      INTEGER NLIMNP, NLIM, NEIGHNP, NLIM_NEIGHNP

! Logical for local debug warnings
      LOGICAL DES_LOC_DEBUG
!-----------------------------------------------

      PEA(NP,:) = .FALSE.

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
                     MARK_PART(NEW_SPOT) = myPE

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
                          CALL CHECK_IF_PARCEL_OVELAPS_STL &
                          (des_pos_new(1:dimn, NEW_SPOT), &
                          PIJK(NEW_SPOT, 1:4), DELETE_PART)

                     IF(.NOT.DELETE_PART) THEN

                        PIP = PIP+1
                        PIP_ADD_COUNT = PIP_ADD_COUNT + 1
                        PEA(NEW_SPOT, 1) = .true.
                        PEA(NEW_SPOT, 2:4) = .false.
                        ! add to the list
                        ALLOCATE(temp_spotlist)
                        temp_spotlist%spot = new_spot
                        temp_spotlist%next => cur_spotlist
                        cur_spotlist => temp_spotlist
                        nullify(temp_spotlist)

                     ELSE
                        PEA(NEW_SPOT, 1) = .false.
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
            WRITE(err_msg,'(/,2x,A,2x,i10)') 'TOTAL NUMBER OF PARCELS ADDED GLOBALLY = ', SUM(ADD_COUNT_ALL(:))

            call flush_err_msg(header = .false., footer = .false.)

            DO IPROC = 0, NUMPES-1
               IF(DMP_LOG) WRITE(UNIT_LOG, '(/,A,i4,2x,A,2x,i5)') 'PARCELS ADDED ON PROC:', IPROC,' EQUAL TO', ADD_COUNT_ALL(IPROC)
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
      USE discretelement, only: max_pip, pea
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

         if(.NOT.PEA(LL,1)) THEN
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
!  Subroutine: PIC_REFLECT_PARTICLE                                  C
!  Purpose:                                                            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE PIC_REFLECT_PART(LL, WALL_NORM)

      USE discretelement, only : dimn, DES_VEL_NEW
      USE mfix_pic, only : MPPIC_COEFF_EN_WALL, &
      & MPPIC_COEFF_ET_WALL
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
      DOUBLE PRECISION :: VEL_NORMMAG_APP, VEL_TANGMAG_APP, TANGENT(DIMN)

      !pre collisional normal and tangential velocity components in vector format
      !APP ==> approach
      DOUBLE PRECISION :: VEL_NORM_APP(DIMN), VEL_TANG_APP(DIMN)


      !post collisional normal and tangential velocity components in vector format
      !SEP ==> separation
      DOUBLE PRECISION :: VEL_NORM_SEP(DIMN), VEL_TANG_SEP(DIMN)

      DOUBLE PRECISION :: COEFF_REST_EN, COEFF_REST_ET

      DOUBLE PRECISION :: VELMOD, NORMMOD, VELDOTNORM
!-----------------------------------------------


      VEL_NORMMAG_APP = DOT_PRODUCT(WALL_NORM(1:DIMN), DES_VEL_NEW(:, LL))


!currently assuming that wall is at rest. Needs improvement for moving wall

      VEL_NORM_APP(1:DIMN) = VEL_NORMMAG_APP*WALL_NORM(1:DIMN)

      VEL_TANG_APP(:) = DES_VEL_NEW(:, LL) - VEL_NORM_APP(1:DIMN)



!post collisional velocities

      !COEFF_REST_EN = REAL_EN_WALL(PIJK(LL,5))
      COEFF_REST_EN = MPPIC_COEFF_EN_WALL
      COEFF_REST_ET = MPPIC_COEFF_ET_WALL

      !if(ep_g(PIJK(LL,4)).lt.0.42) coeff_rest_en = 1.05
      VEL_NORM_SEP(1:DIMN) = -COEFF_REST_EN*VEL_NORM_APP(1:DIMN)

      VEL_TANG_SEP(1:DIMN) =  COEFF_REST_ET*VEL_TANG_APP(1:DIMN)

      DES_VEL_NEW(:, LL) = VEL_NORM_SEP(1:DIMN) + VEL_TANG_SEP(1:DIMN)

      END SUBROUTINE PIC_REFLECT_PART




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
      INTEGER :: I, J, K, IJK
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
