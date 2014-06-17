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
      Implicit none

      DOUBLE PRECISION, INTENT(IN) :: POSITION(DIMN)
      INTEGER , INTENT(IN) :: PCELL(4)
      LOGICAL, INTENT(OUT) :: OVERLAP_EXISTS
      LOGICAL :: INSIDE_DOMAIN, INSIDE_SMALL_CELL, REFLECT_FROM_ORIG_CELL
      INTEGER :: PIP_DEL_COUNT, PIP_ADD_COUNT, REFLECTING_CELL

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

      INCLUDE 'function.inc'

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
      write (vtp_unit,"(15x,es12.6)") (1.d0)
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
      Implicit none


      INTEGER I, J, K, IJK, NF, LL

      DOUBLE PRECISION :: RADSQ, DISTSQ, DIST(DIMN)
      INTEGER :: COUNT_FAC, COUNT, contact_facet_count, NEIGH_CELLS, &
      NEIGH_CELLS_NONNAT, IPROC, COUNT2, &
      LIST_OF_CELLS(27), CELL_ID, I_CELL, J_CELL, K_CELL, cell_count , &
      IMINUS1, IPLUS1, JMINUS1, JPLUS1, KMINUS1, KPLUS1, PHASELL, LOC_MIN_PIP, &
      LOC_MAX_PIP
      INTEGER :: LPIP_DEL_COUNT_ALL(0:numPEs-1), LPIP_ADD_COUNT_ALL(0:numPEs-1)

      
      INTEGER :: PIP_DEL_COUNT, PIP_ADD_COUNT, REFLECTING_CELL
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
      !Old and new position of the parcel 
      double precision, dimension(dimn) :: pt_old, pt_new 
      
      Logical :: checked_facet_already
      !distance parcel travelled out of domain along the facet normal
      double precision :: dist_from_facet
      INTEGER, Parameter :: MAX_FACET_CONTS = 500
      INTEGER :: list_of_checked_facets(max_facet_conts)

      INCLUDE 'function.inc'

      PIP_DEL_COUNT = ZERO

      DO LL = 1, MAX_PIP


         IF(LL.EQ.FOCUS_PARTICLE) DEBUG_DES = .TRUE.

         ! skipping non-existent particles or ghost particles
         IF(.NOT.PEA(LL,1) .OR. PEA(LL,4)) CYCLE


         ! If no neighboring facet in the surrounding 27 cells, then exit
         IF (NO_NEIGHBORING_FACET_DES(PIJK(LL,4))) cycle

         LIST_OF_CELLS(:) = -1
         NEIGH_CELLS = 0
         NEIGH_CELLS_NONNAT  = 0

         CELL_ID = PIJK(LL,4)
         COUNT_FAC = LIST_FACET_AT_DES(CELL_ID)%COUNT_FACETS
         RADSQ = DES_RADIUS(LL)*DES_RADIUS(LL)

         pt_old(1:dimn) = des_pos_old(LL, 1:dimn)
         pt_new(1:dimn) = des_pos_new(LL, 1:dimn)

         IF (COUNT_FAC.gt.0)   then
            !first add the facets in the cell the particle currently resides in
            NEIGH_CELLS = NEIGH_CELLS + 1
            LIST_OF_CELLS(NEIGH_CELLS) = CELL_ID
         ENDIF

         I_CELL = I_OF(CELL_ID)
         J_CELL = J_OF(CELL_ID)
         K_CELL = K_OF(CELL_ID)


         IPLUS1  =  MIN (I_CELL + 1, IEND2)
         IMINUS1 =  MAX (I_CELL - 1, ISTART2)

         JPLUS1  =  MIN (J_CELL + 1, JEND2)
         JMINUS1 =  MAX (J_CELL - 1, JSTART2)

         KPLUS1  =  MIN (K_CELL + 1, KEND2)
         KMINUS1 =  MAX (K_CELL - 1, KSTART2)

         DO K = KMINUS1, KPLUS1
            DO J = JMINUS1, JPLUS1
               DO I = IMINUS1, IPLUS1
                  IJK = FUNIJK(I,J,K)
                  COUNT_FAC = LIST_FACET_AT_DES(IJK)%COUNT_FACETS
                  IF(COUNT_FAC.EQ.0.or.ijk.eq.cell_id) CYCLE
                  NEIGH_CELLS_NONNAT = NEIGH_CELLS_NONNAT + 1
                  NEIGH_CELLS = NEIGH_CELLS + 1
                  LIST_OF_CELLS(NEIGH_CELLS) = IJK
                  !WRITE(*,'(A10, 4(2x,i5))') 'WCELL  = ', IJK, I,J,K
               ENDDO
            ENDDO
         ENDDO

         CONTACT_FACET_COUNT = 0

         CELLLOOP: DO CELL_COUNT = 1, NEIGH_CELLS
            IJK = LIST_OF_CELLS(CELL_COUNT)

            DO COUNT = 1, LIST_FACET_AT_DES(IJK)%COUNT_FACETS
               NF = LIST_FACET_AT_DES(IJK)%FACET_LIST(COUNT)
               ! Neighboring cells will share facets with same facet ID
               ! So it is important to make sure a facet is checked (for speed)
               ! and accounted (for accuracy) only once
               checked_facet_already = .false.
               DO COUNT2 = 1, CONTACT_FACET_COUNT
                  checked_facet_already = (NF.eq.LIST_OF_CHECKED_FACETS(count2))
                  IF(checked_facet_already) exit
               enddo

               IF(checked_facet_already) CYCLE

               CONTACT_FACET_COUNT = CONTACT_FACET_COUNT + 1
               LIST_OF_CHECKED_FACETS(CONTACT_FACET_COUNT) = NF

               line_t  = -Undefined
               !-undefined, because non zero values will imply intersection
               !with the plane
               ontriangle = .false.

               !parametrize a line as p = p_0 + t (p_1 - p_0)
               !and intersect with the triangular plane.
               !if 0<t<1, then this ray connect old and new positions 
               !intersects with the facet. 
               ref_line(1:dimn) = pt_old(1:dimn)
               dir_line(1:dimn) = pt_new(1:dimn) - pt_old(1:dimn)

               norm_plane(1:dimn) = NORM_FACE(1:dimn,NF)
               ref_plane(1:dimn)  = VERTEX(1, 1:dimn,NF)
               CALL intersectLnPlane(ref_line, dir_line, ref_plane, &
                    norm_plane, line_t)
               
                  
                  
               if(line_t.ge.zero.and.line_t.lt.one) then
                  !line_t = one would imply particle sitting on the stl 
                  !face, so don't reflect it yet

                  point_onplane(1:dimn) = ref_line(1:dimn) + &
                       line_t*dir_line(1:dimn)
                  !Now check through barycentric coordinates if
                  !this orthogonal projection lies on the facet or not
                  !If it does, then this point is deemed to be on the
                  !non-fluid side of the facet
                  CALL checkPTonTriangle(point_onplane(1:dimn), &
                       VERTEX(1,:,NF), VERTEX(2,:,NF), VERTEX(3,:,NF), &
                       ontriangle)
                  IF(LL.eq.-11.and.pt_new(3).gt.zlength)then  
                     !IF(LL.eq.1) then 
                     
                     write(*,*) 'POS OLD = ', pt_old(1:dimn)
                     write(*,*) 'POS NEW = ', pt_new(1:dimn)
                     write(*,*) 'plane ref= ', ref_plane(1:dimn)
                     write(*,*) 'norm plane= ', norm_plane(1:dimn)
                     write(*,*) 'dir_line = ', dir_line(1:dimn)
                     
                     WRITE(*,*) 'line_t = ',  line_t,ontriangle
                     write(*,*) 'vel OLD = ', des_vel_old(LL,1:dimn)
                     write(*,*) 'vel NEW = ', des_vel_new(LL,1:dimn)
                  
                     read(*,*)
               endif

                  if(ontriangle) then
                     dist_from_facet = abs(dot_product(pt_new(1:dimn) & 
                          & - point_onplane(1:dimn), norm_plane(1:dimn)))

                     DES_POS_NEW(LL, 1:DIMN) = point_onplane(1:dimn) + & 
                          & dist_from_facet*(norm_plane(1:dimn)) 

                     IF(STL_FACET_TYPE(NF).eq.facet_type_normal.or. &
                          STL_FACET_TYPE(NF).eq.facet_type_mi) then 
                        CALL PIC_REFLECT_PART(LL, norm_plane(1:DIMN))
                        !In the reflect routine, it is assumed that the normal
                        !points into the system which is consitent with the 
                        !definition of normal for facets 

                     ELSEIF(STL_FACET_TYPE(NF).eq.FACET_TYPE_PO) then 
                        PEA(LL, 1)  = .false.
                        PIP_DEL_COUNT = PIP_DEL_COUNT+1
                     ENDIF


                     Exit CELLLOOP
                  endif
               endif

            ENDDO

         end DO CELLLOOP

      end DO

      RETURN 
      
! The following lines could be moved to PIC_MI_BC
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

      END SUBROUTINE PIC_APPLY_WALLBC_STL





!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: PIC_FIND_EMPTY_SPOT                                   C
!  Purpose:                                                            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE PIC_FIND_EMPTY_SPOT(LAST_INDEX, EMPTY_SPOT)
      USE funits
      USE mpi_utility

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
      & 'ERROR IN PIC_FIND_EMPTY_SPOT', /10X, &
      & 'NO MORE EMPTY SPOT IN THE PARTICLE ARRAY TO ADD A NEW PARTICLE',/10X &
      & 'TERMINAL ERROR: STOPPING (CALL RESTART NOT YET IMPLEMENTED BEFORE THIS EXIT)', &
      & /1X,70('*')/)

 2002 FORMAT(/1X,70('*'),//,10X,  &
      & 'ERROR IN PIC_FIND_EMPTY_SPOT', /10X, &
      & 'COULD NOT FIND A SPOT FOR ADDING NEW PARTICLE',/10X &
      & 'INCREASE THE SIZE OF THE INITIAL ARRAYS', 10X, &
      & 'TERMINAL ERROR: STOPPING (CALL RESTART NOT YET IMPLEMENTED BEFORE THIS EXIT)', &
      & /1X,70('*')/)
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


      VEL_NORMMAG_APP = DOT_PRODUCT(WALL_NORM(1:DIMN), DES_VEL_NEW(LL, 1:DIMN))

!currently assuming that wall is at rest. Needs improvement for moving wall

      VEL_NORM_APP(1:DIMN) = VEL_NORMMAG_APP*WALL_NORM(1:DIMN)

      VEL_TANG_APP(:) = DES_VEL_NEW(LL, 1:DIMN) - VEL_NORM_APP(1:DIMN)

      

!post collisional velocities

      !COEFF_REST_EN = REAL_EN_WALL(PIJK(LL,5))
      COEFF_REST_EN = MPPIC_COEFF_EN_WALL
      COEFF_REST_ET = MPPIC_COEFF_ET_WALL

      !if(ep_g(PIJK(LL,4)).lt.0.42) coeff_rest_en = 1.05
      VEL_NORM_SEP(1:DIMN) = -COEFF_REST_EN*VEL_NORM_APP(1:DIMN)

      VEL_TANG_SEP(1:DIMN) =  COEFF_REST_ET*VEL_TANG_APP(1:DIMN)

      DES_VEL_NEW(LL, 1:DIMN) = VEL_NORM_SEP(1:DIMN) + VEL_TANG_SEP(1:DIMN)

      END SUBROUTINE PIC_REFLECT_PART



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: PIC_CHECK_IF_INSIDE_DOMAIN                            C
!  Purpose:                                                            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE PIC_CHECK_IF_INSIDE_DOMAIN(LL, INSIDE_DOMAIN, & 
      REFLECT_FROM_ORIG_CELL,INSIDE_SMALL_CELL )
      
      USE discretelement
      use mpi_utility
      USE cutcell 
      
      IMPLICIT NONE 
      
      LOGICAL, intent(out) ::  INSIDE_DOMAIN, & 
      REFLECT_FROM_ORIG_CELL,INSIDE_SMALL_CELL
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
      INTEGER, INTENT(IN) :: LL
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      DOUBLE PRECISION :: XPOS, YPOS, ZPOS, WNORM(3), DISTMOD
      INTEGER :: I, J, K, IJK
!-----------------------------------------------
! Include statement functions
!-----------------------------------------------
      INCLUDE 'function.inc'
!-----------------------------------------------

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
      END SUBROUTINE PIC_CHECK_IF_INSIDE_DOMAIN



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
! Include statement functions
!-----------------------------------------------
      INCLUDE 'function.inc'
!-----------------------------------------------

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

      END SUBROUTINE PIC_FIND_NEW_CELL
