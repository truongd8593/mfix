
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                                                                                C
!  Module name: des_stl_functions                                                                          C
!  Purpose: This module will contain routines for geometric interacti                       C 
!  required for STL files                                                                                             C
!                                                                                                                                C
!  Author: Rahul Garg                                                                    Date: 24-Oct-13  C
!  Reviewer:                                                                                          Date:            C
!                                                                                                                                C
!  Revision Number #                                                                    Date: ##-###-##  C
!  Author: #                                                                                                                C
!  Purpose: #                                                                                                             C
!                                                                                                                                C 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      MODULE des_stl_functions
      USE param1
      USE funits
      USE run
      USE compar      
      USE discretelement
      USE mfix_pic
      USE cutcell
      Use stl 
      USE indices 
      USE geometry
      USE bc
      USE funits 
      USE mpi_utility
      IMPLICIT NONE      
      !Dont declare any global variables here as it could lead to cyclic 
      !dependency issues 
      !Use this module only to define functions and subroutines 
      CONTAINS 
     
      SUBROUTINE TestTriangleAABB(vert0, vert1, vert2, tri_norm, &
      box_origin, box_extents, sa_exist, sa, i,j,k)
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

         !IF (I .eq. 21.and.j.eq.21.and.k.eq.2.and..false.) then 
         if(.false.) then
            
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
      IMPLICIT NONE
      !point a, pointb, and pointc are the three nodes of the triangle
      !point p is the sphere center 
      double precision, intent(in), dimension(dimn) :: pointp, pointa, pointb, pointc
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

      Subroutine cg_des_convert_to_facets
      USE param1
      USE funits
      USE run
      USE compar      
      USE discretelement
      USE mfix_pic
      USE cutcell
      Use stl
      USE indices 
      USE geometry
      USE bc
      USE funits 

      Implicit None 
      INTEGER :: I, J, K, IJK, NF

      INCLUDE 'function.inc'

      IF(mype.eq.pe_IO) write(*,*) 'des_convert_box_to_facets specified as true', &
      'converting outside box into facets for particle-wall interaction'
      IF(DMP_LOG) write(Unit_log,*) 'des_convert_box_to_facets specified as true', &
      'converting outside box into facets for particle-wall interaction'
      I = IMIN1 !West Face 
      DO K = KMIN1, KMAX1
         DO J = JMIN1, JMAX1               
            IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
            IJK  = FUNIJK(I,J,K)                        
            N_FACETS_DES = N_FACETS_DES + 1 
            NF = N_FACETS_DES
            NORM_FACE(NF,:) = (/one, zero, zero/)
            !NORM_FACE(NF,:) = (/zero, zero, zero/)
            !For stl, the vertices are stored in CCW order looking from outside. 
            !In the stl convention, the normal points outwards. 
            !However, in MFIX, cutcell and facets normals point into the fluid
            !So, the normal's are written as pointing to the fluid but when the stl 
            !is viewed in paraview, it will surely calculate its own normals based on the 
            !order ot the vertices specified here and may look opposite. 
            VERTEX(NF, 1,:) = (/get_nodes(i,j,k, 'w'), get_nodes(i,j,k, 's'), get_nodes(i,j,k, 'b')/)
            VERTEX(NF, 2,:) = (/get_nodes(i,j,k, 'w'), get_nodes(i,j,k, 'n'), get_nodes(i,j,k, 't')/)
            VERTEX(NF, 3,:) = (/get_nodes(i,j,k, 'w'), get_nodes(i,j,k, 's'), get_nodes(i,j,k, 't')/)
            
            N_FACETS_DES = N_FACETS_DES + 1 
            NF = N_FACETS_DES
            NORM_FACE(NF,:) = (/one, zero, zero/)
            !NORM_FACE(NF,:) = (/zero, zero, zero/)
            VERTEX(NF, 1,:) = (/get_nodes(i,j,k, 'w'), get_nodes(i,j,k, 's'), get_nodes(i,j,k, 'b')/)
            VERTEX(NF, 2,:) = (/get_nodes(i,j,k, 'w'), get_nodes(i,j,k, 'n'), get_nodes(i,j,k, 'b')/)
            VERTEX(NF, 3,:) = (/get_nodes(i,j,k, 'w'), get_nodes(i,j,k, 'n'), get_nodes(i,j,k, 't')/)
         enddo
      enddo
      
      I = IMAX1 !East Face 
      DO K = KMIN1, KMAX1
         DO J = JMIN1, JMAX1 
            IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
            IJK  = FUNIJK(I,J,K)            
            N_FACETS_DES = N_FACETS_DES + 1 
            NF = N_FACETS_DES
            NORM_FACE(NF,:) = (/-one, zero, zero/)
            !NORM_FACE(NF,:) = (/zero, zero, zero/)
            VERTEX(NF, 1,:) = (/get_nodes(i,j,k, 'e'), get_nodes(i,j,k, 's'), get_nodes(i,j,k, 'b')/)
            VERTEX(NF, 2,:) = (/get_nodes(i,j,k, 'e'), get_nodes(i,j,k, 's'), get_nodes(i,j,k, 't')/)
            VERTEX(NF, 3,:) = (/get_nodes(i,j,k, 'e'), get_nodes(i,j,k, 'n'), get_nodes(i,j,k, 't')/)
            
            N_FACETS_DES = N_FACETS_DES + 1 
            NF = N_FACETS_DES
            NORM_FACE(NF,:) = (/-one, zero, zero/)
            !NORM_FACE(NF,:) = (/zero, zero, zero/)
            VERTEX(NF, 1,:) = (/get_nodes(i,j,k, 'e'), get_nodes(i,j,k, 's'), get_nodes(i,j,k, 'b')/)
            VERTEX(NF, 2,:) = (/get_nodes(i,j,k, 'e'), get_nodes(i,j,k, 'n'), get_nodes(i,j,k, 't')/)
            VERTEX(NF, 3,:) = (/get_nodes(i,j,k, 'e'), get_nodes(i,j,k, 'n'), get_nodes(i,j,k, 'b')/)
         enddo
      enddo

      J = JMIN1 !south face
      DO K = KMIN1, KMAX1
         DO I = IMIN1, IMAX1 
            IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
            IJK  = FUNIJK(I,J,K)            
            N_FACETS_DES = N_FACETS_DES + 1 
            NF = N_FACETS_DES
            NORM_FACE(NF,:) = (/zero, one, zero/)
            VERTEX(NF, 1,:) = (/get_nodes(i,j,k, 'w'), get_nodes(i,j,k, 's'), get_nodes(i,j,k, 'b')/)
            VERTEX(NF, 3,:) = (/get_nodes(i,j,k, 'e'), get_nodes(i,j,k, 's'), get_nodes(i,j,k, 'b')/)
            VERTEX(NF, 2,:) = (/get_nodes(i,j,k, 'e'), get_nodes(i,j,k, 's'), get_nodes(i,j,k, 't')/)
            
            N_FACETS_DES = N_FACETS_DES + 1 
            NF = N_FACETS_DES
            NORM_FACE(NF,:) = (/zero, one, zero/)
            VERTEX(NF, 1,:) = (/get_nodes(i,j,k, 'w'), get_nodes(i,j,k, 's'), get_nodes(i,j,k, 'b')/)
            VERTEX(NF, 2,:) = (/get_nodes(i,j,k, 'w'), get_nodes(i,j,k, 's'), get_nodes(i,j,k, 't')/)
            VERTEX(NF, 3,:) = (/get_nodes(i,j,k, 'e'), get_nodes(i,j,k, 's'), get_nodes(i,j,k, 't')/)
 
           
         enddo
      enddo

      J = JMAX1 !north  face
      DO K = KMIN1, KMAX1
         DO I = IMIN1, IMAX1 
            IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
            IJK  = FUNIJK(I,J,K)            
            N_FACETS_DES = N_FACETS_DES + 1 
            NF = N_FACETS_DES
            NORM_FACE(NF,:) = (/zero, -one, zero/)
            VERTEX(NF, 1,:) = (/get_nodes(i,j,k, 'w'), get_nodes(i,j,k, 'n'), get_nodes(i,j,k, 'b')/)
            VERTEX(NF, 2,:) = (/get_nodes(i,j,k, 'e'), get_nodes(i,j,k, 'n'), get_nodes(i,j,k, 'b')/)
            VERTEX(NF, 3,:) = (/get_nodes(i,j,k, 'e'), get_nodes(i,j,k, 'n'), get_nodes(i,j,k, 't')/)
            
            N_FACETS_DES = N_FACETS_DES + 1 
            NF = N_FACETS_DES
            NORM_FACE(NF,:) = (/zero, -one, zero/)
            VERTEX(NF, 1,:) = (/get_nodes(i,j,k, 'w'), get_nodes(i,j,k, 'n'), get_nodes(i,j,k, 'b')/)
            VERTEX(NF, 2,:) = (/get_nodes(i,j,k, 'e'), get_nodes(i,j,k, 'n'), get_nodes(i,j,k, 't')/)
            VERTEX(NF, 3,:) = (/get_nodes(i,j,k, 'w'), get_nodes(i,j,k, 'n'), get_nodes(i,j,k, 't')/)
          
         enddo
      enddo

      K = KMIN1 !bottom face 
      DO J = JMIN1, JMAX1
         DO I = IMIN1, IMAX1 
            IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
            IJK  = FUNIJK(I,J,K)            
            N_FACETS_DES = N_FACETS_DES + 1 
            NF = N_FACETS_DES
            NORM_FACE(NF,:) = (/zero, zero, one/)
            VERTEX(NF, 1,:) = (/get_nodes(i,j,k, 'w'), get_nodes(i,j,k, 's'), get_nodes(i,j,k, 'b')/)
            VERTEX(NF, 2,:) = (/get_nodes(i,j,k, 'e'), get_nodes(i,j,k, 'n'), get_nodes(i,j,k, 'b')/)
            VERTEX(NF, 3,:) = (/get_nodes(i,j,k, 'w'), get_nodes(i,j,k, 'n'), get_nodes(i,j,k, 'b')/)
            N_FACETS_DES = N_FACETS_DES + 1 
            NF = N_FACETS_DES

            NORM_FACE(NF,:) = (/zero, zero, one/)           
            VERTEX(NF, 1,:) = (/get_nodes(i,j,k, 'w'), get_nodes(i,j,k, 's'), get_nodes(i,j,k, 'b')/)
            VERTEX(NF, 2,:) = (/get_nodes(i,j,k, 'e'), get_nodes(i,j,k, 's'), get_nodes(i,j,k, 'b')/)
            VERTEX(NF, 3,:) = (/get_nodes(i,j,k, 'e'), get_nodes(i,j,k, 'n'), get_nodes(i,j,k, 'b')/)
         enddo
      enddo
      
      K = KMAX1 !top face 
      DO J = JMIN1, JMAX1
         DO I = IMIN1, IMAX1 
            IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
            IJK  = FUNIJK(I,J,K)            
            N_FACETS_DES = N_FACETS_DES + 1 
            NF = N_FACETS_DES
            NORM_FACE(NF,:) = (/zero, zero, -one/)
            VERTEX(NF, 1,:) = (/get_nodes(i,j,k, 'w'), get_nodes(i,j,k, 's'), get_nodes(i,j,k, 't')/)
            VERTEX(NF, 2,:) = (/get_nodes(i,j,k, 'w'), get_nodes(i,j,k, 'n'), get_nodes(i,j,k, 't')/)
            VERTEX(NF, 3,:) = (/get_nodes(i,j,k, 'e'), get_nodes(i,j,k, 'n'), get_nodes(i,j,k, 't')/)
            
            N_FACETS_DES = N_FACETS_DES + 1 
            NF = N_FACETS_DES
            NORM_FACE(NF,:) = (/zero, zero, -one/)          
            VERTEX(NF, 1,:) = (/get_nodes(i,j,k, 'w'), get_nodes(i,j,k, 's'), get_nodes(i,j,k, 't')/)
            VERTEX(NF, 2,:) = (/get_nodes(i,j,k, 'e'), get_nodes(i,j,k, 'n'), get_nodes(i,j,k, 't')/)
            VERTEX(NF, 3,:) = (/get_nodes(i,j,k, 'e'), get_nodes(i,j,k, 's'), get_nodes(i,j,k, 't')/)
         enddo
      enddo
      !CALL  DEBUG_write_all_readin_facets
      !call mfix_exit(mype)
      end Subroutine cg_des_convert_to_facets

      Subroutine bin_facets_to_grid_des
      implicit none     
      INTEGER :: IJK,I,J,K, I1, I2, J1, J2, K1, K2, N, II, JJ, KK, count_fac
      INTEGER :: IM,IP,JM,JP,KM,KP,IMJK,IPJK,IJMK,IJPK,IJKM,IJKP
      INTEGER :: IJPKP,IPJKP,IPJPK

      DOUBLE PRECISION:: x1,y1,z1,x2,y2,z2,x3,y3,z3

      INTEGER :: IJK2,CURRENT_I,CURRENT_J,CURRENT_K


!      CHARACTER (LEN=3) :: CAD_PROPAGATE_ORDER

      include "function.inc"    

      DO N = 1,N_FACETS_DES


         X1 = MINVAL(VERTEX(N,1:3,1))
         X2 = MAXVAL(VERTEX(N,1:3,1))
         Y1 = MINVAL(VERTEX(N,1:3,2))
         Y2 = MAXVAL(VERTEX(N,1:3,2))
         Z1 = MINVAL(VERTEX(N,1:3,3))
         Z2 = MAXVAL(VERTEX(N,1:3,3))
         
         I1 = IEND3
         I2 = ISTART3

         IF(X2>=ZERO.AND.X1<=XLENGTH+TOL_STL ) THEN
            DO I = ISTART3, IEND3
               IP = I+1
               IF(XG_E(I)>=X1-TOL_STL) THEN
                  I1=I
                  EXIT
               ENDIF
            ENDDO

            DO I = IEND3, ISTART3,-1
               IP = I+1
               IF(XG_E(I)-DX(I)<=X2+TOL_STL) THEN
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
               IF(YG_N(J)>=Y1-TOL_STL) THEN
                  J1=J
                  EXIT
               ENDIF
            ENDDO

            DO J = JEND3, JSTART3,-1
               JP=J+1
               IF(YG_N(J)-DY(J)<=Y2+TOL_STL) THEN
                  J2=J
                  EXIT
               ENDIF
            ENDDO
         ENDIF

         K1 = KEND3
         K2 = KSTART3

         IF(Z2>=ZERO.AND.Z1<=ZLENGTH+TOL_STL) THEN
            DO K = KSTART3, KEND3
               KP=K+1

               IF(ZG_T(K)>=Z1-TOL_STL) THEN
                  K1=K
                  EXIT
               ENDIF
            ENDDO

            DO K = KEND3, KSTART3,-1
               KP = K+1
               IF(ZG_T(K)-DZ(K)<=Z2+TOL_STL) THEN
                  K2=K
                  EXIT
               ENDIF
            ENDDO
         ENDIF

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
!  Module name: ADD_FACET_for_des                                              C
!  Purpose: Add facet to list in IJK scalar cell for the Lagrangian modules                       C
!                                                                      C
!  Author: Rahul Garg                              Date: 24-Oct-13  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE ADD_FACET_FOR_DES(I,J,K,IJK,N)
    
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: I,j,k,IJK, N
      !Local variables
      INTEGER ::  CURRENT_COUNT, COUNT
      CHARACTER*10000 :: err_message 
      CHARACTER*80 :: temp_message 
      
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
      
      box_origin(1) = xg_e(I) - dx(I)
      box_origin(2) = yg_n(J) - dy(J)
      box_origin(3) = zg_t(K) - dz(K)
      
      box_extents(1) = dx(I)
      box_extents(2) = dy(J)
      box_extents(3) = dz(K)
      
      !Do the separating axis test to check if a separating axis exist. If the separating 
      !axis exsit then the cell and facet cannot intersect, so return without adding. 
      CALL TestTriangleAABB(vertex(N, 1,:), vertex(N,2,:), vertex(N,3,:), norm_face(N,:), &
      box_origin(:), box_extents(:), sa_exist, sep_axis,i,j,k )

      IF (I .eq. -1.and.j.eq.21.and.k.eq.2) then 
         write(*, '(5x,A10, 2x, i10) ') 'Facet number: ', N
         write(*,'(5x,A10, 3(2x, g17.8))') 'vert1: ', vertex(N,1,:)
         write(*,'(5x,A10,3(2x, g17.8))') 'vert1: ', vertex(N,2,:)
         write(*,'(5x,A10,3(2x, g17.8))') 'vert1: ', vertex(N,3,:)
         write(*,'(5x,A10,3(2x, g17.8))') 'norm: ', norm_face(N,:)
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
         IF(PRINT_DES_SCREEN) WRITE(*, 200) MAX_FACETS_PER_CELL_DES, IJK, & 
         I_OF(IJK), J_OF(IJK), K_OF(IJK), mype, & 
         IS_ON_myPE_owns(I_OF(IJK), J_OF(IJK), K_OF(IJK))
         IF(DMP_LOG) WRITE(UNIT_LOG, 200) MAX_FACETS_PER_CELL_DES, IJK, & 
         I_OF(IJK), J_OF(IJK), K_OF(IJK), mype, & 
         IS_ON_myPE_owns(I_OF(IJK), J_OF(IJK), K_OF(IJK))
         
         err_message  = "current_list for this cell is"
!         write(err_message, '(/,10X,"current_list")')         
         
         IF(nodesI*nodesJ*nodesK.gt.1) then 
            
            WRITE(fname,'(A,"_TROUBLE_CELL",A, I5.5, 3(A,I5.5), ".stl")') & 
            TRIM(RUN_NAME), '_pid', mype, '_I', I_OF(IJK), '_J', J_OF(IJK), '_K', K_OF(IJK)         
         else
            
            WRITE(fname,'(A,"_TROUBLE_CELL", 3(A,I5.5), ".stl")') & 
            TRIM(RUN_NAME), '_I', I_OF(IJK), '_J', J_OF(IJK), '_K', K_OF(IJK) 
         endif

         open(stl_unit, file = fname, form='formatted')
         write(stl_unit,*)'solid vcg'      

         DO COUNT  = 1, CURRENT_COUNT 
            
             !write(*, '(/,I10)') LIST_FACET_AT(IJK)%FACET_LIST(COUNT)
             write(err_message, '(A, I10)') TRIM(err_message), LIST_FACET_AT_DES(IJK)%FACET_LIST(COUNT)
             
             FID = LIST_FACET_AT_DES(IJK)%FACET_LIST(COUNT)
             
             write(stl_unit,*) '   facet normal ', NORM_FACE(FID,1:3)
             write(stl_unit,*) '      outer loop' 
             write(stl_unit,*) '         vertex ', VERTEX(FID,1,1:3)
             write(stl_unit,*) '         vertex ', VERTEX(FID,2,1:3)
             write(stl_unit,*) '         vertex ', VERTEX(FID,3,1:3)
             write(stl_unit,*) '      endloop'
             write(stl_unit,*) '   endfacet'
             
          ENDDO
          
          write(stl_unit,*)'endsolid vcg'
          close(stl_unit, status = 'keep')
          IF(myPE.eq.pe_IO) WRITE(*, *) trim(err_message)
          IF(DMP_LOG) WRITE(UNIT_LOG, *) trim(err_message)

         CALL MFIX_EXIT(myPE)
      ENDIF
      
      
 200  FORMAT(/1X,70('*'),//,10X,  & 
      & 'ERROR MESSAGE FROM CUT_CELL_PREPROCESSING', /10x, & 
      & 'INCREASE MAX_FACETS_PER_CELL_DES from the current value of', i3, /10x, &
      & 'Happening for cell IJK, I, J, K = ', 4(2x, i5), /10X, &
      & 'mype, Is on myPe? ', I6, L2, /10X, &
      & 'see the file TROUBLE_CELL for all the current facets in this cell', /10X, &
      & 'Terminal Error: STOPPING', /10X, &
      & /1X,70('*')/)
      
      END SUBROUTINE ADD_FACET_FOR_DES

     

      SUBROUTINE DEBUG_WRITE_GRID_FACEINFO
      USE param
      USE param1
      USE parallel
      USE constant
      USE run
      USE toleranc
      USE geometry
      USE indices  
      USE compar
      USE sendrecv
      USE quadric
      USE cutcell
      USE polygon
      USE stl
      USE mpi_utility

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

               WRITE(1001, '(2X, "CELL IJK, I, J, K =        = ", i20, 2x, 4(2x,i10))') CELL_ID, I_OF(CELL_ID), J_OF(CELL_ID), K_OF(CELL_ID)

               WRITE(1001, '(2x, "TOTAL FACETS                  = ", 3(2x, i10))') COUNT_FACETS
               
               DO COUNT = 1, COUNT_FACETS
                  WRITE(1001, '(2x, i20)')  LIST_FACET_AT_DES(CELL_ID)%FACET_LIST(COUNT)
               ENDDO
               
            ENDDO
         ENDDO
      ENDDO
      
      CLOSE(1001, STATUS = "keep")
      END    SUBROUTINE DEBUG_WRITE_GRID_FACEINFO
      
      Subroutine  DEBUG_write_all_readin_facets
      IMPLICIT NONE 
      
      INTEGER ::  CELL_ID, N, I, J, K, COUNT, COUNT_FACETS, IJK
      CHARACTER*100 :: FILENAME

      INCLUDE 'function.inc'

      OPEN(UNIT=444, FILE='geometry_from_readin_facets.stl') 
      write(444,*)'solid vcg'      
      write(*,*) 'NFACETS, NFACETS DES =', N_FACETS, N_FACETS_DES
    !  DO N = N_FACETS+1, N_FACETS_DES
      DO N = 1, N_FACETS_DES
         write(444,*) '   facet normal ', NORM_FACE(N,1:3)
         write(444,*) '      outer loop' 
         write(444,*) '         vertex ', VERTEX(N,1,1:3)
         write(444,*) '         vertex ', VERTEX(N,2,1:3)
         write(444,*) '         vertex ', VERTEX(N,3,1:3)
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

      Subroutine  DEBUG_write_stl_from_grid_facet(WRITE_FACETS_EACH_CELL)
      IMPLICIT NONE 
      
      LOGICAL, INTENT(IN),optional  :: WRITE_FACETS_EACH_CELL
      INTEGER ::  CELL_ID, N, I, J, K, COUNT, COUNT_FACETS, IJK
      CHARACTER*100 :: FILENAME
      LOGICAL :: FACET_WRITTEN(DIM_STL), write_each_cell 

      INCLUDE 'function.inc'

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
                     write(445,*) '   facet normal ', NORM_FACE(N,:)
                     write(445,*) '      outer loop'
                     write(445,*) '         vertex ', VERTEX(N,1,1:3)
                     write(445,*) '         vertex ', VERTEX(N,2,1:3)
                     write(445,*) '         vertex ', VERTEX(N,3,1:3)
                     write(445,*) '      endloop'
                     write(445,*) '   endfacet'
                  endif
                  
                  if (facet_written(n)) cycle 
                  write(444,*) '   facet normal ', NORM_FACE(N,:)
                  write(444,*) '      outer loop'
                  write(444,*) '         vertex ', VERTEX(N,1,1:3)
                  write(444,*) '         vertex ', VERTEX(N,2,1:3)
                  write(444,*) '         vertex ', VERTEX(N,3,1:3)
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

      IF(MyPE == PE_IO) THEN
         WRITE(*,*) ' The file geometry_from_grid_facets.stl was sucessfully written.'
         WRITE(*,*) ' This is the based on facets stored grid wise for DES modules'
         WRITE(*,*) ' and is provided for convenience and debugging (it is not used).'
      ENDIF
      END Subroutine  DEBUG_write_stl_from_grid_facet


      double precision function get_nodes(I, J,K, IDIR) 
      implicit none     
      Integer, intent(in) :: I, J, K
      Character*1, intent(in) :: Idir
      Integer :: IJK 
      INCLUDE 'function.inc'
      
      SELECT CASE((TRIM(IDIR)))
      CASE('w' )
         GET_NODES = XG_E(I) - DX(I)
      CASE('e')
         GET_NODES = XG_E(I)
      CASE('n')
         GET_NODES = YG_N(J)  
      CASE('s')
         GET_NODES = YG_N(J) - DY(J) 
      CASE('t')
         GET_NODES = ZG_T(K)
      CASE('b')
         GET_NODES = ZG_T(K) - DZ(K) 
      CASE DEFAULT
         
         WRITE(*,*)'ERROR IN get_nodes under des_stl_functions.'
         WRITE(*,*)'UNKNOWN direction:  ', trim(idir)
         WRITE(*,*)'ACCEPTABLE directions are E,W,N,S,T,B:'
                                !            WRITE(*,*)'STL'
         CALL MFIX_EXIT(myPE)

      END SELECT
         
      return 
      end function get_nodes



      end       MODULE des_stl_functions
      
      
      
