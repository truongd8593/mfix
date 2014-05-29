!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: des_stl_functions                                      !
!  Author: Rahul Garg                                 Date: 24-Oct-13  !
!                                                                      !
!  Purpose: This module containd routines for geometric interaction    ! 
!  required for STL files.                                             !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      MODULE des_stl_functions

      IMPLICIT NONE

! Use this module only to define functions and subroutines.
      CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: des_stl_preprocessing                                   !
!  Author: Rahul Garg                                 Date: 24-Oct-13  !
!                                                                      !
!  Purpose:                                                            !
!                                                                      ! 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DES_STL_PREPROCESSING
      USE stl
      USE param
      USE error_manager
      USE discretelement, only: DES_CONVERT_BOX_TO_FACETS
      USE cutcell, only: use_stl
      implicit none
      integer :: ijk, count
      integer :: max_des_facets

      CALL INIT_ERR_MSG("DES_STL_PREPROCESSING")

      !max_des_facets  = 2*(2*(imax*jmax) + 2*(jmax*kmax) + 2*(kmax*imax))
      !allocate(vertex_des(max_des_facets))

      CALL ALLOCATE_DES_STL_ARRAYS 

! Pre-procssing for the des in order to assign facets to grid cells.
      WRITE(ERR_MSG,"('Pre-Processing geometry for DES.')")
      CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
      N_FACETS_DES = 0

! Set N_facets_des to add any more facets needed by dem and not to.
! contaminate the Eulerian-Eulerian CG stuff 
      IF(USE_STL) N_FACETS_DES = N_FACETS

! Triangulate default walls (bounding box of the simulation)
      CALL CG_DES_CONVERT_TO_FACETS

      CALL BIN_FACETS_TO_GRID_DES

      DO IJK = 1, DIMENSION_3
         COUNT = LIST_FACET_AT_DES(IJK)%COUNT_FACETS
         IF(COUNT.eq.0) DEALLOCATE(LIST_FACET_AT_DES(IJK)%FACET_LIST)
      ENDDO

      !CALL DEBUG_WRITE_GRID_FACEINFO
      CALL DEBUG_write_stl_from_grid_facet(WRITE_FACETS_EACH_CELL=.false.)

!      call DEBUG_write_all_readin_facets

      WRITE(ERR_MSG,"('Done pre-Processing geometry for DES.')")
      CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE DES_STL_PREPROCESSING

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: ALLOCATE_DES_STL_ARRAYS                                 !
!  Author: Rahul Garg                                 Date: 24-Oct-13  !
!                                                                      !
!  Purpose:                                                            !
!                                                                      ! 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE ALLOCATE_DES_STL_ARRAYS
      USE param
      USE stl

      IMPLICIT NONE

      integer :: ijk

      ALLOCATE(LIST_FACET_AT_DES(DIMENSION_3))

      DO IJK = 1, DIMENSION_3
         LIST_FACET_AT_DES(IJK)%COUNT_FACETS = 0
         ALLOCATE(LIST_FACET_AT_DES(IJK)%FACET_LIST(MAX_FACETS_PER_CELL_DES))
      ENDDO

      ALLOCATE(NO_NEIGHBORING_FACET_DES(DIMENSION_3))
      NO_NEIGHBORING_FACET_DES = .false.

      END SUBROUTINE ALLOCATE_DES_STL_ARRAYS


      SUBROUTINE TestTriangleAABB(vert0, vert1, vert2, tri_norm, &
      box_origin, box_extents, sa_exist, sa, i,j,k)
      USE stl
      Implicit none

      !Separating axis test algorithm from Realtimecollision book by
      !Christer Ericson
      double precision, intent(in), dimension(3) :: vert0, vert1, vert2, tri_norm
      double precision, intent(in), dimension(3) :: box_origin, box_extents
      logical, intent(out) :: sa_exist
!      double precision, intent(out) :: proj_box, proj_tri(3)
      Integer, intent(out) :: sa
      Integer, intent(in) :: i,j,k

      double precision :: p0, p1, p2, r, e0, e1, e2

      double precision :: c_orig(3), c(3)

      !box face normals
      double precision, dimension(3) :: u0, u1, u2

      double precision, dimension(3) :: v0, v1, v2 , vert0_orig
      double precision, dimension(3) :: f0, f1, f2

      !13 possible Separating axes
      double precision, dimension(9, 3) :: sep_axis
      !triangle plane n and d
      double precision :: n(3), d, tol_tri_aabb_proj

      double precision :: axis_magsq, s
      integer :: count
      sa_exist = .false.
      sa = 0
      tol_tri_aabb_proj = 1e-06

      !Initially set e0, e1, and e2 to box extents

      !Compute box center and extents (if not already given in that format)
      c_orig(:) = box_origin(:) + box_extents(:)* 0.5d0

      !Translate triangle as conceptually moving AABB's center to origin

      v0 = vert0 - c_orig
      v1 = vert1 - c_orig
      v2 = vert2 - c_orig

      !Scale everything by box_extents
      v0(1) = v0(1)/box_extents(1)
      v1(1) = v1(1)/box_extents(1)
      v2(1) = v2(1)/box_extents(1)

      v0(2) = v0(2)/box_extents(2)
      v1(2) = v1(2)/box_extents(2)
      v2(2) = v2(2)/box_extents(2)

      v0(3) = v0(3)/box_extents(3)
      v1(3) = v1(3)/box_extents(3)
      v2(3) = v2(3)/box_extents(3)

      !set the box origin to (0,0,0)
      c = 0.d0
      !e0, e1, and e2 are half box extents which equals 0.5d0 in the scaled co-ordinates
      e0 = 0.5d0
      e1 = 0.5d0
      e2 = 0.5d0
      u0 = (/ 1.d0, 0.d0, 0.d0/)
      u1 = (/ 0.d0, 1.d0, 0.d0/)
      u2 = (/ 0.d0, 0.d0, 1.d0/)

!!$  write(*,'(A20,5(2x,g17.8))') 'box center orgi = ',c_orig(:)
!!$  write(*,'(A20,5(2x,g17.8))') 'box center new = ',c(:)
!!$  write(*,'(A20,5(2x,g17.8))') 'box extents = ',e0,e1,e2
!!$  write(*,'(A20,5(2x,g17.8))') 'tri vert1 = ',v0(:)
!!$  write(*,'(A20,5(2x,g17.8))') 'tri vert2 = ',v1(:)
!!$  write(*,'(A20,5(2x,g17.8))') 'tri vert3 = ',v2(:)

      !Compute edge vectors for triangle
      f0 = v1 - v0
      f1 = v2 - v1
      f2 = v0 - v2
      !Category 3 (edge-edge cross products)...
      !total of 9 (3 from triangle cross 3 independent from box)
      sep_axis(1, 1:3) = (/0.d0, -f0(3), f0(2)/)
      sep_axis(2, 1:3) = (/0.d0, -f1(3), f1(2)/)
      sep_axis(3, 1:3) = (/0.d0, -f2(3), f2(2)/)
      sep_axis(4, 1:3) = (/f0(3),  0.d0, -f0(1)/)
      sep_axis(5, 1:3) = (/f1(3),  0.d0, -f1(1)/)
      sep_axis(6, 1:3) = (/f2(3),  0.d0, -f2(1)/)
      sep_axis(7, 1:3) = (/-f0(2), f0(1), 0.d0/)
      sep_axis(8, 1:3) = (/-f1(2), f1(1), 0.d0/)
      sep_axis(9, 1:3) = (/-f2(2), f2(1), 0.d0/)

      do count = 1, 9
         axis_magsq = (sep_axis(count,1)**2+sep_axis(count,2)**2 &
         +sep_axis(count,3)**2)
         if(axis_magsq < TOL_STL*TOL_STL) THEN
            !nearly a zero vector. This will happen if the cross-products
            !from category 3 are formed from co-planar vectors
            !WRITE(*,'(/5x,A20, i2,/5x, A, 3(2x,g17.8))') 'Axis number:', count, &
            !     'Ignored due to zero vec:', sep_axis(count,:)
            CYCLE
         endif

         r =    e0*abs(dot_product(u0(:), sep_axis(count,:))) &
         + e1*abs(dot_product(u1(:), sep_axis(count,:))) &
         + e2*abs(dot_product(u2(:), sep_axis(count,:)))
         p0 = dot_product(v0(:), sep_axis(count, :))
         p1 = dot_product(v1(:), sep_axis(count, :))
         p2 = dot_product(v2(:), sep_axis(count, :))

         IF (I .eq. 2.and.j.eq.3.and.k.eq.2.and..false.) then
            !if(.false.) then

            write(*,'(A6,5(2x,g17.8))')'c:',c_orig
            write(*,'(A6,5(2x,g17.8))')'vorig:',vert0
            write(*,'(A6,5(2x,g17.8))')'v0:',v0
            write(*,'(A6,5(2x,g17.8))')'v1:',v1
            write(*,'(A6,5(2x,g17.8))')'extents:',e0,e1,e2
            write(*,'(A6,5(2x,g17.8))')'sa:', count, sep_axis(count,:)
            write(*,'(A6,5(2x,g17.8))')'data:',  r, p0, p1, p2
            write(*,'(A6,(2x,g17.8),"<",2(2x,g17.8), ">", 2x,g17.8)')'data2:', max(p0, p1, p2)*1e8, -r*1e8, min(p0, p1, p2)*1e8,r*1e8

            write(*,'(A6,5(2x,L2))') 'data2:', max(p0, p1, p2)< -r, min(p0, p1, p2) > r

             write(*,'(A6,(2x,g17.8),"<",2(2x,g17.8), ">", 2x,g17.8)') 'data2:', (max(p0, p1, p2)+tol_tri_aabb_proj)*1e8, -r*1e8, min(p0, p1, p2)*1e8, (r+tol_tri_aabb_proj)*1e8
            write(*,'(A6,5(2x,L2))') 'data2:', max(p0, p1, p2)+tol_tri_aabb_proj< -r, min(p0, p1, p2) > r+tol_tri_aabb_proj
         endif

         if (max(p0, p1, p2) +tol_tri_aabb_proj < -r .or. min(p0,p1,p2) > r+tol_tri_aabb_proj) then
            sa = count
            sa_exist = .true.
            return
         endif
      enddo
      ! Test the three axes corresponding to the face normals of AABB b (category 1).
      !Could have been through dot products like above, but explotiing the
      !the simple definition of box face normals to optimize the code.

      ! Exit if...
      ! ... [-e0, e0] and [min(v0.x,v1.x,v2.x), max(v0.x,v1.x,v2.x)] do not overlap
      if (Max(v0(1), v1(1), v2(1)) +tol_tri_aabb_proj< -e0 .or. Min(v0(1), v1(1), v2(1)) > e0+tol_tri_aabb_proj) then
         sa = 10
         sa_exist = .true.
         return
      endif
      !... [-e1, e1] and [min(v0.y,v1.y,v2.y), max(v0.y,v1.y,v2.y)] do not overlap
      if (Max(v0(2), v1(2), v2(2)) +tol_tri_aabb_proj< -e1 .or. Min(v0(2), v1(2), v2(2)) > e1+tol_tri_aabb_proj) then
         sa = 11
         sa_exist = .true.
         return
      endif

      ! ... [-e2, e2] and [min(v0.z,v1.z,v2.z), max(v0.z,v1.z,v2.z)] do not overlap
      if (Max(v0(3), v1(3), v2(3)) +tol_tri_aabb_proj< -e2 .or. Min(v0(3), v1(3), v2(3)) > e2+tol_tri_aabb_proj) then
         sa = 12
         sa_exist = .true.
         return
      endif

      r =    e0*abs(dot_product(u0(:), tri_norm(:))) &
      + e1*abs(dot_product(u1(:), tri_norm(:))) &
      + e2*abs(dot_product(u2(:), tri_norm(:)))

      !Compute distance of box center from plane
      !s = Dot_product(tri_norm(:), c(:) - v0(:))
      !c = zero in the shfted co-ordinates
      s = Dot_product(tri_norm(:), v0(:))
      !write(*,*) r,s
      !Intersection occurs when distance s falls within [-r,+r] interval
      if(Abs(s) > r+ tol_tri_aabb_proj) then
         !projections do not intersect, so the separating axis exists
         sa_exist = .true.
         sa = 13
      endif
      end subroutine TestTriangleAABB


      Subroutine ClosestPtPointTriangle(pointp, pointa, pointb, pointc, closest_point)
      USE param1, only: zero, one
      USE discretelement, only: dimn
      !USE funits
      !USE run
      !USE compar
      !USE mfix_pic
      !USE cutcell
      !Use stl
      !USE indices
      !USE geometry
      !USE bc
      !USE funits
      !USE mpi_utility

      IMPLICIT NONE
      !point a, pointb, and pointc are the three nodes of the triangle
      !point p is the sphere center
      double precision, intent(in), dimension(3) :: pointa, pointb, pointc
      double precision, intent(in), dimension(dimn) :: pointp
      double precision, intent(out), dimension(dimn) ::  closest_point
      !Local variables
      double precision, dimension(dimn) :: ab, ac, ap, bp,cp
      double precision :: d1, d2, d3, d4, vc, v, d5, d6, vb, w, va, denom
      ab = pointb - pointa
      ac = pointc - pointa
      ap = pointp - pointa
      d1 = DOT_PRODUCT(ab, ap)
      d2 = DOT_PRODUCT(ac, ap)

      IF(d1 <= Zero .AND. d2 <= zero) then
         closest_point = pointa ! barycentric coordinates (1,0,0)
         return
      end if

      !Check if P in vertex region outside B
      bp = pointp - pointb
      d3 = DOT_PRODUCT(ab, bp);
      d4 = DOT_PRODUCT(ac, bp);
      if (d3 >= zero .and. d4 <= d3) then
         closest_point = pointb !barycentric coordinates (0,1,0)
         return
      endif

      ! Check if P in edge region of AB, if so return projection of P onto AB
      vc = d1*d4 - d3*d2;
      if (vc <= zero .and. d1 >= zero .and. d3 <= zero) then
         v = d1 / (d1 - d3);
         closest_point =  pointa + v * ab; ! barycentric coordinates (1-v,v,0)
         return
      end if

      !Check if P in vertex region outside C
      cp = pointp - pointc
      d5 = DOT_PRODUCT(ab, cp)
      d6 = DOT_PRODUCT(ac, cp)
      if (d6 >= zero .and. d5 <= d6) then
         closest_point  = pointc ! barycentric coordinates (0,0,1)
         return
      endif

      !Check if P in edge region of AC, if so return projection of P onto AC
      vb = d5*d2 - d1*d6

      if (vb <= zero .and. d2 >= zero .and. d6 <= zero) then
         w = d2 / (d2 - d6)
         closest_point = pointa + w * ac ! barycentric coordinates (1-w,0,w)
         return
      end if

      !Check if P in edge region of BC, if so return projection of P onto BC
      va = d3*d6 - d5*d4
      if (va <= zero .and.(d4 - d3) >= zero .and. (d5 - d6) >= zero) then
         w = (d4 - d3) / ((d4 - d3) + (d5 - d6))
         closest_point = pointb + w * (pointc - pointb) !barycentric coordinates (0,1-w,w)
         return
      end if


      !P inside face region. Compute Q through its barycentric coordinates (u,v,w)
      denom = one / (va + vb + vc)
      v = vb * denom
      w = vc * denom
      closest_point = pointa + ab * v + ac * w; ! = u*a + v*b + w*c, u = va * denom = 1.0f - v - w
      return
      end Subroutine ClosestPtPointTriangle

      subroutine checkPTonTriangle(pointp, pointa, pointb, pointc, on_trian)
      Use param1, only : one
      USE discretelement, only: dimn
      USE stl, only: tol_stl
      IMPLICIT NONE
      !point a, pointb, and pointc are the three nodes of the triangle
      !point p is the test point
      double precision, intent(in), dimension(3) :: pointa, pointb, pointc
      double precision, intent(in), dimension(dimn) :: pointp
      logical, intent(out)::  on_trian
      !Local variables
      !triangle edges
      double precision, dimension(dimn) :: v0, v1, v2
      double precision :: d00, d01, d11, d20, d21, denom
      !barcycentric coordinates
      double precision :: v, w

      logical :: v_positive, w_positive, VplusW_positive
      v0 = pointb - pointa
      v1 = pointc - pointa
      v2 = pointp - pointa;
      d00 = DOT_PRODUCT(v0, v0)
      d01 = DOT_PRODUCT(v0, v1);
      d11 = DOT_PRODUCT(v1, v1);
      d20 = DOT_PRODUCT(v2, v0);
      d21 = DOT_PRODUCT(v2, v1);
      denom = d00 * d11 - d01 * d01;
      v = (d11 * d20 - d01 * d21) / denom;
      w = (d00 * d21 - d01 * d20) / denom;
      !u = 1.0f - v - w;

      V_POSITIVE = (v>=-TOL_STL)
      W_POSITIVE = (w>=-TOL_STL)
      VplusW_positive = ((v+w)<=ONE+TOL_STL)

      ON_TRIAN = (V_POSITIVE.AND.W_POSITIVE.AND.VplusW_positive)

      RETURN
      end subroutine checkPTonTriangle


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: intersectLnPlane                                        !
!  Author: Rahul Garg                                 Date: 24-Oct-13  !
!                                                                      !
!  Purpose:                                                            !
!                                                                      ! 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      Subroutine intersectLnPlane(ref_line, dir_line, ref_plane, norm_plane, line_param)
      USE discretelement, only: dimn
      USE param1, only: zero
      IMPLICIT NONE
      !reference point and direction of the line
      double precision, intent(in), dimension(dimn) :: ref_line,  dir_line
      !reference point and normal of the plane
      double precision, intent(in), dimension(dimn) :: ref_plane, norm_plane

      !line is parameterized as p = p_ref + t * dir_line, t is line_param
      double precision, intent(out) :: line_param

      !local vars
      double precision :: denom


      denom = DOT_PRODUCT(dir_line, norm_plane)
      if(denom.gt.zero) then
         line_param = DOT_PRODUCT(ref_plane(:) - ref_line(:), norm_plane(:))
         line_param = line_param/denom
      endif
      return
      end subroutine intersectLnPlane

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: cg_des_convert_to_facets                                !
!  Author: Rahul Garg                                 Date: 24-Oct-13  !
!                                                                      !
!  Purpose:                                                            !
!                                                                      ! 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      Subroutine CG_DES_CONVERT_TO_FACETS
      USE param1
      USE run
      USE compar
      Use stl
      USE indices
      USE geometry
      USE compar
      USE error_manager
      Implicit None
      INTEGER :: I, J, K, IJK, NF

      INCLUDE 'function.inc'

      CALL INIT_ERR_MSG("CG_DES_CONVERT_TO_FACETS")

      I = IMIN1 !West Face
      DO K = KMIN1, KMAX1
         DO J = JMIN1, JMAX1
            IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
            !IF(.not.fluid_at(FUNIJK(I+1,J,K))) cycle

            IJK  = FUNIJK(I,J,K)
            N_FACETS_DES = N_FACETS_DES + 1
            NF = N_FACETS_DES
            NORM_FACE(:,NF) = (/one, zero, zero/)
            !NORM_FACE(:,NF) = (/zero, zero, zero/)
            !For stl, the vertices are stored in CCW order looking from outside.
            !In the stl convention, the normal points outwards.
            !However, in MFIX, cutcell and facets normals point into the fluid
            !So, the normal's are written as pointing to the fluid but when the stl
            !is viewed in paraview, it will surely calculate its own normals based on the
            !order ot the vertices specified here and may look opposite.
            VERTEX(1,:,NF) = (/get_nodes(i,j,k, 'w'), get_nodes(i,j,k, 's'), get_nodes(i,j,k, 'b')/)
            VERTEX(2,:,NF) = (/get_nodes(i,j,k, 'w'), get_nodes(i,j,k, 'n'), get_nodes(i,j,k, 't')/)
            VERTEX(3,:,NF) = (/get_nodes(i,j,k, 'w'), get_nodes(i,j,k, 's'), get_nodes(i,j,k, 't')/)

            N_FACETS_DES = N_FACETS_DES + 1
            NF = N_FACETS_DES
            NORM_FACE(:,NF) = (/one, zero, zero/)
            !NORM_FACE(:,NF) = (/zero, zero, zero/)
            VERTEX(1,:,NF) = (/get_nodes(i,j,k, 'w'), get_nodes(i,j,k, 's'), get_nodes(i,j,k, 'b')/)
            VERTEX(2,:,NF) = (/get_nodes(i,j,k, 'w'), get_nodes(i,j,k, 'n'), get_nodes(i,j,k, 'b')/)
            VERTEX(3,:,NF) = (/get_nodes(i,j,k, 'w'), get_nodes(i,j,k, 'n'), get_nodes(i,j,k, 't')/)
         enddo
      enddo

      I = IMAX1 !East Face
      DO K = KMIN1, KMAX1
         DO J = JMIN1, JMAX1
            IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
            !IF(.not.fluid_at(FUNIJK(I-1,J,K))) cycle

            IJK  = FUNIJK(I,J,K)
            N_FACETS_DES = N_FACETS_DES + 1
            NF = N_FACETS_DES
            NORM_FACE(:,NF) = (/-one, zero, zero/)
            !NORM_FACE(:,NF) = (/zero, zero, zero/)
            VERTEX(1,:,NF) = (/get_nodes(i,j,k, 'e'), get_nodes(i,j,k, 's'), get_nodes(i,j,k, 'b')/)
            VERTEX(2,:,NF) = (/get_nodes(i,j,k, 'e'), get_nodes(i,j,k, 's'), get_nodes(i,j,k, 't')/)
            VERTEX(3,:,NF) = (/get_nodes(i,j,k, 'e'), get_nodes(i,j,k, 'n'), get_nodes(i,j,k, 't')/)

            N_FACETS_DES = N_FACETS_DES + 1
            NF = N_FACETS_DES
            NORM_FACE(:,NF) = (/-one, zero, zero/)
            !NORM_FACE(:,NF) = (/zero, zero, zero/)
            VERTEX(1,:,NF) = (/get_nodes(i,j,k, 'e'), get_nodes(i,j,k, 's'), get_nodes(i,j,k, 'b')/)
            VERTEX(2,:,NF) = (/get_nodes(i,j,k, 'e'), get_nodes(i,j,k, 'n'), get_nodes(i,j,k, 't')/)
            VERTEX(3,:,NF) = (/get_nodes(i,j,k, 'e'), get_nodes(i,j,k, 'n'), get_nodes(i,j,k, 'b')/)
         enddo
      enddo

      J = JMIN1 !south face
      DO K = KMIN1, KMAX1
         DO I = IMIN1, IMAX1
            IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
            !IF(.not.fluid_at(FUNIJK(I,J+1,K))) cycle

            IJK  = FUNIJK(I,J,K)
            N_FACETS_DES = N_FACETS_DES + 1
            NF = N_FACETS_DES
            NORM_FACE(:,NF) = (/zero, one, zero/)
            VERTEX(1,:,NF) = (/get_nodes(i,j,k, 'w'), get_nodes(i,j,k, 's'), get_nodes(i,j,k, 'b')/)
            VERTEX(3,:,NF) = (/get_nodes(i,j,k, 'e'), get_nodes(i,j,k, 's'), get_nodes(i,j,k, 'b')/)
            VERTEX(2,:,NF) = (/get_nodes(i,j,k, 'e'), get_nodes(i,j,k, 's'), get_nodes(i,j,k, 't')/)

            N_FACETS_DES = N_FACETS_DES + 1
            NF = N_FACETS_DES
            NORM_FACE(:,NF) = (/zero, one, zero/)
            VERTEX(1,:,NF) = (/get_nodes(i,j,k, 'w'), get_nodes(i,j,k, 's'), get_nodes(i,j,k, 'b')/)
            VERTEX(2,:,NF) = (/get_nodes(i,j,k, 'w'), get_nodes(i,j,k, 's'), get_nodes(i,j,k, 't')/)
            VERTEX(3,:,NF) = (/get_nodes(i,j,k, 'e'), get_nodes(i,j,k, 's'), get_nodes(i,j,k, 't')/)
         enddo
      enddo

      J = JMAX1 !north  face
      DO K = KMIN1, KMAX1
         DO I = IMIN1, IMAX1
            IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE

            !IF(.not.fluid_at(FUNIJK(I,J-1,K))) cycle
            IJK  = FUNIJK(I,J,K)
            N_FACETS_DES = N_FACETS_DES + 1
            NF = N_FACETS_DES
            NORM_FACE(:,NF) = (/zero, -one, zero/)
            VERTEX(1,:,NF) = (/get_nodes(i,j,k, 'w'), get_nodes(i,j,k, 'n'), get_nodes(i,j,k, 'b')/)
            VERTEX(2,:,NF) = (/get_nodes(i,j,k, 'e'), get_nodes(i,j,k, 'n'), get_nodes(i,j,k, 'b')/)
            VERTEX(3,:,NF) = (/get_nodes(i,j,k, 'e'), get_nodes(i,j,k, 'n'), get_nodes(i,j,k, 't')/)

            N_FACETS_DES = N_FACETS_DES + 1
            NF = N_FACETS_DES
            NORM_FACE(:,NF) = (/zero, -one, zero/)
            VERTEX(1,:,NF) = (/get_nodes(i,j,k, 'w'), get_nodes(i,j,k, 'n'), get_nodes(i,j,k, 'b')/)
            VERTEX(2,:,NF) = (/get_nodes(i,j,k, 'e'), get_nodes(i,j,k, 'n'), get_nodes(i,j,k, 't')/)
            VERTEX(3,:,NF) = (/get_nodes(i,j,k, 'w'), get_nodes(i,j,k, 'n'), get_nodes(i,j,k, 't')/)
         enddo
      enddo

      IF(DO_K) THEN      
         K = KMIN1 !bottom face
         DO J = JMIN1, JMAX1
         DO I = IMIN1, IMAX1
            IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE

            !IF(.not.fluid_at(FUNIJK(I,J,K+1))) cycle
            IJK  = FUNIJK(I,J,K)
            N_FACETS_DES = N_FACETS_DES + 1
            NF = N_FACETS_DES
            NORM_FACE(:,NF) = (/zero, zero, one/)
            VERTEX(1,:,NF) = (/get_nodes(i,j,k, 'w'), get_nodes(i,j,k, 's'), get_nodes(i,j,k, 'b')/)
            VERTEX(2,:,NF) = (/get_nodes(i,j,k, 'e'), get_nodes(i,j,k, 'n'), get_nodes(i,j,k, 'b')/)
            VERTEX(3,:,NF) = (/get_nodes(i,j,k, 'w'), get_nodes(i,j,k, 'n'), get_nodes(i,j,k, 'b')/)
            N_FACETS_DES = N_FACETS_DES + 1
            NF = N_FACETS_DES

            NORM_FACE(:,NF) = (/zero, zero, one/)
            VERTEX(1,:,NF) = (/get_nodes(i,j,k, 'w'), get_nodes(i,j,k, 's'), get_nodes(i,j,k, 'b')/)
            VERTEX(2,:,NF) = (/get_nodes(i,j,k, 'e'), get_nodes(i,j,k, 's'), get_nodes(i,j,k, 'b')/)
            VERTEX(3,:,NF) = (/get_nodes(i,j,k, 'e'), get_nodes(i,j,k, 'n'), get_nodes(i,j,k, 'b')/)
         enddo
         enddo

         K = KMAX1 !top face
         DO J = JMIN1, JMAX1
         DO I = IMIN1, IMAX1
            IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE

            !IF(.not.fluid_at(FUNIJK(I,J,K-1))) cycle
            IJK  = FUNIJK(I,J,K)
            N_FACETS_DES = N_FACETS_DES + 1
            NF = N_FACETS_DES
            NORM_FACE(:,NF) = (/zero, zero, -one/)
            VERTEX(1,:,NF) = (/get_nodes(i,j,k, 'w'), get_nodes(i,j,k, 's'), get_nodes(i,j,k, 't')/)
            VERTEX(2,:,NF) = (/get_nodes(i,j,k, 'w'), get_nodes(i,j,k, 'n'), get_nodes(i,j,k, 't')/)
            VERTEX(3,:,NF) = (/get_nodes(i,j,k, 'e'), get_nodes(i,j,k, 'n'), get_nodes(i,j,k, 't')/)

            N_FACETS_DES = N_FACETS_DES + 1
            NF = N_FACETS_DES
            NORM_FACE(:,NF) = (/zero, zero, -one/)
            VERTEX(1,:,NF) = (/get_nodes(i,j,k, 'w'), get_nodes(i,j,k, 's'), get_nodes(i,j,k, 't')/)
            VERTEX(2,:,NF) = (/get_nodes(i,j,k, 'e'), get_nodes(i,j,k, 'n'), get_nodes(i,j,k, 't')/)
            VERTEX(3,:,NF) = (/get_nodes(i,j,k, 'e'), get_nodes(i,j,k, 's'), get_nodes(i,j,k, 't')/)
         enddo
         enddo
      endif

      CALL FINL_ERR_MSG
      end Subroutine cg_des_convert_to_facets

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: bin_facets_to_grid_des                                  !
!  Author: Rahul Garg                                 Date: 24-Oct-13  !
!                                                                      !
!  Purpose:                                                            !
!                                                                      ! 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      Subroutine BIN_FACETS_TO_GRID_DES

      USE param1, only: zero
      USE discretelement, only: dimn, xe, yn, zt
      Use stl, only: vertex, NO_NEIGHBORING_FACET_DES, LIST_FACET_AT_DES
      Use stl, only: tol_stl, n_facets_des
      USE indices
      USE geometry
      USE mpi_utility
      implicit none
      INTEGER :: IJK,I,J,K, I1, I2, J1, J2, K1, K2, N, II, JJ, KK, count_fac
      INTEGER :: IM,IP,JM,JP,KM,KP,IMJK,IPJK,IJMK,IJPK,IJKM,IJKP
      INTEGER :: IJPKP,IPJKP,IPJPK

      DOUBLE PRECISION:: x1,y1,z1,x2,y2,z2,x3,y3,z3

      INTEGER :: IJK2,CURRENT_I,CURRENT_J,CURRENT_K

      include "function.inc"    

!      CHARACTER (LEN=3) :: CAD_PROPAGATE_ORDER

      DO N = 1,N_FACETS_DES

         X1 = MINVAL(VERTEX(1:3,1,N))
         X2 = MAXVAL(VERTEX(1:3,1,N))
         Y1 = MINVAL(VERTEX(1:3,2,N))
         Y2 = MAXVAL(VERTEX(1:3,2,N))
         Z1 = MINVAL(VERTEX(1:3,3,N))
         Z2 = MAXVAL(VERTEX(1:3,3,N))

         I1 = IEND3
         I2 = ISTART3

         IF(X2>=ZERO.AND.X1<=XLENGTH+TOL_STL ) THEN
            DO I = ISTART3, IEND3
               IP = I+1
               IF(XE(I)>=X1-TOL_STL) THEN
                  I1=I
                  EXIT
               ENDIF
            ENDDO

            DO I = IEND3, ISTART3,-1
               IP = I+1
               IF(XE(I)-DX(I)<=X2+TOL_STL) THEN
                  I2=I
                  EXIT
               ENDIF
            ENDDO
         ENDIF

         J1 = JEND3
         J2 = JSTART3

         IF(Y2>=ZERO.AND.Y1<=YLENGTH+TOL_STL) THEN
            DO J = JSTART3, JEND3
               JP = J+1
               IF(YN(J)>=Y1-TOL_STL) THEN
                  J1=J
                  EXIT
               ENDIF
            ENDDO

            DO J = JEND3, JSTART3,-1
               JP=J+1
               IF(YN(J)-DY(J)<=Y2+TOL_STL) THEN
                  J2=J
                  EXIT
               ENDIF
            ENDDO
         ENDIF

         K1 = KEND3
         K2 = KSTART3

         IF(DO_K) THEN
            IF(Z2>=ZERO.AND.Z1<=ZLENGTH+TOL_STL) THEN
               DO K = KSTART3, KEND3
                  KP=K+1

                  IF(ZT(K)>=Z1-TOL_STL) THEN
                     K1=K
                     EXIT
                  ENDIF
               ENDDO

               DO K = KEND3, KSTART3,-1
                  KP = K+1
                  IF(ZT(K)-DZ(K)<=Z2+TOL_STL) THEN
                     K2=K
                     EXIT
                  ENDIF
               ENDDO
            ENDIF
         ENDIF

         IF(N.eq.92.and..false.) then
            write(*,*) 'vertex x ', VERTEX(1:3,1,N)
            write(*,*) 'vertex y ', VERTEX(1:3,2,N)
            write(*,*) 'vertex z ', VERTEX(1:3,3,N)
            write(*,*) 'I1, I2', I1, I2, J1, J2, K1, K2
            read(*,*)
         endif

         DO K=K1,K2
            DO J=J1,J2
               DO I=I1,I2

                  IJK = FUNIJK(I,J,K)
                  CALL ADD_FACET_FOR_DES(I,J,K,IJK,N)
               enddo
            enddo
         enddo
      enddo

      DO K = KSTART1, KEND1
         DO J = JSTART1, JEND1
            DO I = ISTART1, IEND1
               IJK = FUNIJK(I, J, K)

               I1 =  MAX( I - 1, ISTART2)
               I2 =  MIN( I + 1, IEND2)

               J1 =  MAX( J - 1, JSTART2)
               J2 =  MIN (J + 1, JEND2)

               K1 =  MAX( K - 1, KSTART2)
               K2 =  MIN (K + 1, KEND2)

               count_fac = 0
               NO_NEIGHBORING_FACET_DES(IJK)  = .false.
               DO KK = K1, k2
                  DO JJ = J1, j2
                     DO II = I1, i2
                        IJK2 = FUNIJK(II, JJ, KK)
                        count_fac = count_fac +  LIST_FACET_AT_DES(IJK2)%COUNT_FACETS
                        ENDDO
                     ENDDO
                  ENDDO

                  if(count_fac.eq.0)  then
                     NO_NEIGHBORING_FACET_DES(IJK)  = .true.
                  endif
               ENDDO
            ENDDO
         ENDDO

      End Subroutine bin_facets_to_grid_des

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: ADD_FACET_for_des                                      C
!  Purpose: Add facet to list in IJK scalar cell for the               C
!           discrete modules.                                          C
!                                                                      C
!  Author: Rahul Garg                              Date: 24-Oct-13     C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE ADD_FACET_FOR_DES(I,J,K,IJK,N)
      USE param1, only: zero, one
      USE discretelement, only: dimn, xe, yn, zt
      Use stl
      USE indices
      USE geometry
      USE mpi_utility
      USE error_manager
      USE run
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: I,j,k,IJK, N
      !Local variables
      INTEGER ::  CURRENT_COUNT, COUNT

      double precision ::   box_origin(3), box_extents(3), tri_norm(3)
      Logical :: sa_exist
      integer :: sep_axis

      CHARACTER*100 :: FNAME
      integer :: stl_unit, fid
      INCLUDE 'function.inc'

      stl_unit = 1001
      IF (I.lt.IMIN1.OR.I.gt.IMAX1) RETURN
      IF (J.lt.JMIN1.OR.J.gt.JMAX1) RETURN
      IF (K.lt.KMIN1.OR.K.gt.KMAX1) RETURN


      box_origin(1) = xe(I) - dx(I)
      box_origin(2) = yn(J) - dy(J)
      box_origin(3) = merge(0.0d0, zt(K) - dz(K), NO_K)

      box_extents(1) = dx(I)
      box_extents(2) = dy(J)
      box_extents(3) = merge(ZLENGTH, dz(K), NO_K)

      !Do the separating axis test to check if a separating axis exist. If the separating
      !axis exsit then the cell and facet cannot intersect, so return without adding.
      CALL TestTriangleAABB(vertex(1,:,N), vertex(2,:,N), vertex(3,:,N), norm_face(:,N), &
      box_origin(:), box_extents(:), sa_exist, sep_axis,i,j,k )

      IF (I .eq. 2.and.j.eq.3.and.k.eq.2.and..not.sa_exist.and..false.) then

         write(*, *) box_origin(:)
         write(*, '(5x,A10, 2x, i10) ') 'Facet number: ', N
         write(*,'(5x,A10, 3(2x, g17.8))') 'vert1: ', vertex(1,:,N)
         write(*,'(5x,A10,3(2x, g17.8))') 'vert1: ', vertex(2,:,N)
         write(*,'(5x,A10,3(2x, g17.8))') 'vert1: ', vertex(3,:,N)
         write(*,'(5x,A10,3(2x, g17.8))') 'norm: ', norm_face(:,N)
         write(*,'(5x,A25, L2, 2x, i2)') 'sep_axis exist, axis #:', &
         sa_exist, sep_axis
         read(*,*)
      endif

      IF (sa_exist) then
         !write(*,'(5x,A25, L2, 2x, i2)') 'sep_axis exist, axis #:', &
         !sa_exist, sep_axis
         return
      ENDIF
      !separating axis does not exist ==> cell and the triangle do intersect.

      CURRENT_COUNT = LIST_FACET_AT_DES(IJK)%COUNT_FACETS

      IF(CURRENT_COUNT .LT. MAX_FACETS_PER_CELL_DES) THEN
         LIST_FACET_AT_DES(IJK)%COUNT_FACETS = CURRENT_COUNT+1
         LIST_FACET_AT_DES(IJK)%FACET_LIST(CURRENT_COUNT+1) = N

      ELSE
         CALL INIT_ERR_MSG("add_facets_for_des  under des_stl_functions_mod")
         WRITE(err_msg, 200) MAX_FACETS_PER_CELL_DES, IJK, &
         I, J, K, mype, &
         IS_ON_myPE_owns(I, J, K)
         CALL flush_err_msg(footer = .false.)


         write(err_msg, *) "current_list for this cell is"
         CALL flush_err_msg(header = .false., footer = .false.)


         IF(nodesI*nodesJ*nodesK.gt.1) then

            WRITE(fname,'(A,"_TROUBLE_CELL",A, I5.5, 3(A,I5.5), ".stl")') &
            TRIM(RUN_NAME), '_pid', mype, '_I', I, '_J', J, '_K', K
         else

            WRITE(fname,'(A,"_TROUBLE_CELL", 3(A,I5.5), ".stl")') &
            TRIM(RUN_NAME), '_I', I, '_J', J, '_K', K
         endif

         open(stl_unit, file = fname, form='formatted')
         write(stl_unit,*)'solid vcg'

         DO COUNT  = 1, CURRENT_COUNT

             !write(*, '(/,I10)') LIST_FACET_AT(IJK)%FACET_LIST(COUNT)
            FID = LIST_FACET_AT_DES(IJK)%FACET_LIST(COUNT)
            write(err_msg, '(I10)') FID

            CALL flush_err_msg(header = .false., footer = .false.)


             write(stl_unit,*) '   facet normal ', NORM_FACE(1:3,FID)
             write(stl_unit,*) '      outer loop'
             write(stl_unit,*) '         vertex ', VERTEX(1,1:3,FID)
             write(stl_unit,*) '         vertex ', VERTEX(2,1:3,FID)
             write(stl_unit,*) '         vertex ', VERTEX(3,1:3,FID)
             write(stl_unit,*) '      endloop'
             write(stl_unit,*) '   endfacet'

          ENDDO

          write(stl_unit,*)'endsolid vcg'
          close(stl_unit, status = 'keep')


          write(err_msg, *) 'Stopping'
          CALL flush_err_msg(header = .false., abort = .true.)
      ENDIF



 200  FORMAT(&
      & 'ERROR MESSAGE FROM CUT_CELL_PREPROCESSING', /10x, &
      & 'INCREASE MAX_FACETS_PER_CELL_DES from the current value of', i3, /10x, &
      & 'Happening for cell IJK, I, J, K = ', 4(2x, i5), /10X, &
      & 'mype, Is on myPe? ', I6, L2, /10X, &
      & 'see the file TROUBLE_CELL for all the current facets in this cell')

      END SUBROUTINE ADD_FACET_FOR_DES



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: DEBUG_WRITE_GRID_FACEINFO                               !
!  Author: Rahul Garg                                 Date: 24-Oct-13  !
!                                                                      !
!  Purpose:                                                            !
!                                                                      ! 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DEBUG_WRITE_GRID_FACEINFO
      USE param
      USE param1
      USE run
      USE geometry
      USE indices
      USE compar
      USE stl

      IMPLICIT NONE

      INTEGER :: CELL_ID, I, J, K, COUNT, COUNT_FACETS, IJK

      CHARACTER*100 :: FILENAME

      INCLUDE 'function.inc'


      IF(nodesI*nodesJ*nodesK.gt.1) then
         write(filename,'(A,"_FACETS_GRID_CELL_",I5.5,".dat")')  trim(run_name), myPE
      else
         write(filename,'(A,"_FACETS_GRID_CELL",".dat")')  trim(run_name)
      endif
      OPEN(1001, file = TRIM(filename), form ='formatted')
      DO K=KSTART2, KEND2
         DO J=JSTART2, JEND2
            DO I=ISTART2, IEND2
               CELL_ID = FUNIJK(I,J,K)
               COUNT_FACETS =  LIST_FACET_AT_DES(CELL_ID)%COUNT_FACETS
               IF(COUNT_FACETS.eq.0) cycle
               WRITE(1001, '("**************************************************")')

               WRITE(1001, '(2X, "CELL IJK, I, J, K =        = ", i20, 2x, 4(2x,i10))') CELL_ID, I, J, K

               WRITE(1001, '(2x, "TOTAL FACETS                  = ", 3(2x, i10))') COUNT_FACETS

               DO COUNT = 1, COUNT_FACETS
                  WRITE(1001, '(2x, i20)')  LIST_FACET_AT_DES(CELL_ID)%FACET_LIST(COUNT)
               ENDDO

            ENDDO
         ENDDO
      ENDDO

      CLOSE(1001, STATUS = "keep")
      END    SUBROUTINE DEBUG_WRITE_GRID_FACEINFO

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: DEBUG_write_all_readin_facets                           !
!  Author: Rahul Garg                                 Date: 24-Oct-13  !
!                                                                      !
!  Purpose:                                                            !
!                                                                      ! 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      Subroutine  DEBUG_write_all_readin_facets
      USE stl
      USE compar
      IMPLICIT NONE

      INTEGER ::  CELL_ID, N, I, J, K, COUNT, COUNT_FACETS, IJK
      CHARACTER*100 :: FILENAME


      OPEN(UNIT=444, FILE='geometry_from_readin_facets.stl')
      write(444,*)'solid vcg'
      write(*,*) 'NFACETS, NFACETS DES =', N_FACETS, N_FACETS_DES
    !  DO N = N_FACETS+1, N_FACETS_DES
      DO N = 1, N_FACETS_DES
         write(444,*) '   facet normal ', NORM_FACE(1:3,N)
         write(444,*) '      outer loop'
         write(444,*) '         vertex ', VERTEX(1,1:3,N)
         write(444,*) '         vertex ', VERTEX(2,1:3,N)
         write(444,*) '         vertex ', VERTEX(3,1:3,N)
         write(444,*) '      endloop'
         write(444,*) '   endfacet'

      ENDDO

      write(444,*)'endsolid vcg'

      close(444)

      IF(MyPE == PE_IO) THEN
         WRITE(*,*) ' The file geometry_from_readin_facets.stl was sucessfully written.'
         WRITE(*,*) ' This is the based on readin and newly generated facets'
      ENDIF
      END Subroutine DEBUG_write_all_readin_facets

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: DEBUG_write_stl_from_grid_facet                         !
!  Author: Rahul Garg                                 Date: 24-Oct-13  !
!                                                                      !
!  Purpose:                                                            !
!                                                                      ! 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      Subroutine  DEBUG_write_stl_from_grid_facet(WRITE_FACETS_EACH_CELL)
      use run
      USE stl

      USE geometry
      USE indices
      USE compar

      IMPLICIT NONE

      LOGICAL, INTENT(IN),optional  :: WRITE_FACETS_EACH_CELL
      INTEGER ::  CELL_ID, N, I, J, K, COUNT, COUNT_FACETS, IJK
      CHARACTER*100 :: FILENAME
      LOGICAL :: write_each_cell
      LOGICAL, DIMENSION(:), allocatable :: FACET_WRITTEN

      INCLUDE 'function.inc'

      ALLOCATE (FACET_WRITTEN(DIM_STL))

      write_each_cell = .false.
      if(present(WRITE_FACETS_EACH_CELL)) then
         write_each_cell = WRITE_FACETS_EACH_CELL
      endif

      FACET_WRITTEN = .false.

      IF(nodesI*nodesJ*nodesK.gt.1) then
         write(filename,'(A,"_GEOMETRY_FROM_GRID_FACETS_",I5.5,".stl")')  &
         trim(run_name), myPE
      else
         write(filename,'(A,"_GEOMETRY_FROM_GRID_FACETS",".stl")')  &
         trim(run_name)
      endif
      OPEN(UNIT=444, FILE=trim(filename))
      write(444,*)'solid vcg'

      DO K=KSTART2, KEND2
         DO J=JSTART2, JEND2
            DO I=ISTART2, IEND2
               CELL_ID = FUNIJK(I,J,K)
               COUNT_FACETS =  LIST_FACET_AT_DES(CELL_ID)%COUNT_FACETS
               IF(COUNT_FACETS.eq.0) cycle

               if(write_each_cell) then
                  write(filename, '(A,"_geometry_", i3.3, "_", i3.3, "_", i3.3,"_" , i8.8, ".stl")') trim(run_name) , I,J,K,CELL_ID
                  OPEN(UNIT=445, FILE=filename)
                  write(445,*)'solid vcg'
               endif

               DO COUNT = 1, COUNT_FACETS
                  N = LIST_FACET_AT_DES(CELL_ID)%FACET_LIST(COUNT)

                  if(write_each_cell) then
                     write(445,*) '   facet normal ', NORM_FACE(:,N)
                     write(445,*) '      outer loop'
                     write(445,*) '         vertex ', VERTEX(1,1:3,N)
                     write(445,*) '         vertex ', VERTEX(2,1:3,N)
                     write(445,*) '         vertex ', VERTEX(3,1:3,N)
                     write(445,*) '      endloop'
                     write(445,*) '   endfacet'
                  endif

                  if (facet_written(n)) cycle
                  write(444,*) '   facet normal ', NORM_FACE(:,N)
                  write(444,*) '      outer loop'
                  write(444,*) '         vertex ', VERTEX(1,1:3,N)
                  write(444,*) '         vertex ', VERTEX(2,1:3,N)
                  write(444,*) '         vertex ', VERTEX(3,1:3,N)
                  write(444,*) '      endloop'
                  write(444,*) '   endfacet'
                  facet_written(N) = .true.
               ENDDO

               if(write_each_cell) then
                  write(445,*)'endsolid vcg'
                  close(445)
               endif

            ENDDO
         ENDDO
      ENDDO
      write(444,*)'endsolid vcg'

      close(444)

      DEALLOCATE (FACET_WRITTEN)

      IF(MyPE == PE_IO) THEN
         WRITE(*,*) ' The file geometry_from_grid_facets.stl was sucessfully written.'
         WRITE(*,*) ' This is the based on facets stored grid wise for DES modules'
         WRITE(*,*) ' and is provided for convenience and debugging (it is not used).'
      ENDIF
      END Subroutine  DEBUG_write_stl_from_grid_facet


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Function: get_nodes                                                 !
!  Author: Rahul Garg                                 Date: 24-Oct-13  !
!                                                                      !
!  Purpose: Return the coordinate of the I/J/K for direction IDIR.     !
!                                                                      ! 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      DOUBLE PRECISION FUNCTION GET_NODES(I, J, K, IDIR) 

! Scalar grid face locations.
      use discretelement, only: XE, YN, ZT
! Scalar grid dimensions.
      use geometry, only: DX, DY, DZ
! Flag: Use Kth direction.
      use geometry, only: DO_K, ZLENGTH

      use error_manager 

      IMPLICIT NONE

      INTEGER, INTENT(in) :: I, J, K
      CHARACTER*1, INTENT(in) :: Idir

      SELECT CASE((TRIM(IDIR)))


      CASE('w'); GET_NODES = XE(I) - DX(I)
      CASE('e'); GET_NODES = XE(I)

      CASE('n'); GET_NODES = YN(J)
      CASE('s'); GET_NODES = YN(J) - DY(J)

      CASE('t'); GET_NODES = merge(ZT(K), ZLENGTH, DO_K)
      CASE('b'); GET_NODES = merge(ZT(K) - DZ(K), 0.0d0, DO_K)

! Catch illegal usage.
      CASE DEFAULT
         CALL INIT_ERR_MSG("DES_STL_FUNCTIONS --> GET_NODES")
         WRITE(ERR_MSG, 1100) trim(idir)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      END SELECT

 1100 FORMAT('Error 1100: Illegal usage. UNKNOWN direction: ',A,/      &
         'Valid directions are E, W, N, S, T, B.')

      RETURN
      END FUNCTION GET_NODES

      END MODULE DES_STL_FUNCTIONS      


