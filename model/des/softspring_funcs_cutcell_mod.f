!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: softspring_funcs_cutcell                               C
!                                                                      C
!  Purpose: subroutines for particle-wall collisions when cutcell is   C
!           used. Also contains rehack of routines for cfslide and     C
!           cffctow which might be different from the stand alone      C
!           routines. Eventually all the DEM routines will be          C
!           consolidated.                                              C
!                                                                      C
!                                                                      C
! CALC_FORCE_WITH_WALL_CUTFACE_STL: Particle-wall driver routine when  C
!    stl files are used                                                C
!                                                                      C
!  Author: Rahul Garg                               Date: 1-Dec-2013   C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      module softspring_funcs_cutcell

      PRIVATE 
      PUBLIC:: CHECK_IF_PARTICLE_OVELAPS_STL, CALC_DEM_FORCE_WITH_WALL_STL, &
      cffctow2, CFRELVEL2, CFSLIDE2  
        
      
      CONTAINS
      

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  SUBROUTINE: CHECK_IF_PARTICLE_OVELAPS_STL                           C
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
      SUBROUTINE CHECK_IF_PARTICLE_OVELAPS_STL(POSITION, PCELL, RADIUS, & 
      OVERLAP_EXISTS)

      USE run
      USE param1
      USE discretelement, only: dimn, xe, yn, zt
      USE geometry
      USE constant
      USE cutcell
      USE indices
      USE stl
      USE compar 
      USE des_stl_functions
      Implicit none

      DOUBLE PRECISION, INTENT(IN) :: POSITION(DIMN), RADIUS 
      INTEGER , INTENT(IN) :: PCELL(4)
      LOGICAL, INTENT(OUT) :: OVERLAP_EXISTS

      INTEGER I, J, K, IJK, NF

      DOUBLE PRECISION :: RADSQ, DISTSQ, DIST(DIMN), CLOSEST_PT(DIMN)
      INTEGER :: COUNT_FAC, COUNT, contact_facet_count, NEIGH_CELLS, &
      NEIGH_CELLS_NONNAT, &
      LIST_OF_CELLS(27), CELL_ID, I_CELL, J_CELL, K_CELL, cell_count , &
      IMINUS1, IPLUS1, JMINUS1, JPLUS1, KMINUS1, KPLUS1, PHASELL, LOC_MIN_PIP, &
      LOC_MAX_PIP, focus_particle 
      
      
      INCLUDE 'function.inc'

      FOCUS_PARTICLE = -1

      OVERLAP_EXISTS = .false. 
           
      IF (NO_NEIGHBORING_FACET_DES(PCELL(4))) RETURN 

      LIST_OF_CELLS(:) = -1
      NEIGH_CELLS = 0
      NEIGH_CELLS_NONNAT  = 0
      CELL_ID = PCELL(4)
      COUNT_FAC = LIST_FACET_AT_DES(CELL_ID)%COUNT_FACETS

      RADSQ = RADIUS*RADIUS

      IF (COUNT_FAC.gt.0)   then
         !first add the facets in the cell the particle currently resides in
         NEIGH_CELLS = NEIGH_CELLS + 1
         LIST_OF_CELLS(NEIGH_CELLS) = CELL_ID
      ENDIF

      I_CELL = PCELL(1)
      J_CELL = PCELL(2) 
      K_CELL = PCELL(3)

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
               IF(COUNT_FAC.EQ.0) CYCLE
               distsq = zero
               IF(POSITION(1) > XE(I)) DISTSQ = DISTSQ &
               + (POSITION(1)-XE(I))*(POSITION(1)-XE(I))
               
               IF(POSITION(1) < XE(I) - DX(I)) DISTSQ = DISTSQ &
               + (XE(I) - DX(I) - POSITION(1))*(XE(I) - DX(I) - POSITION(1))
               
               IF(POSITION(2) > YN(J)) DISTSQ = DISTSQ &
               + (POSITION(2)-YN(J))* (POSITION(2)-YN(J))

               IF(POSITION(2) < YN(J) - DY(J)) DISTSQ = DISTSQ &
               + (YN(J) - DY(J) - POSITION(2))* (YN(J) - DY(J) - POSITION(2))
               
               IF(POSITION( 3) > ZT(K)) DISTSQ = DISTSQ &
               + (POSITION(3)-ZT(K))*(POSITION(3)-ZT(K))
               
               IF(POSITION(3) < ZT(K) - DZ(K)) DISTSQ = DISTSQ &
               + (ZT(K) - DZ(K) - POSITION(3))*(ZT(K) - DZ(K) - POSITION(3))
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
            
            CALL ClosestPtPointTriangle(POSITION(:), &
            VERTEX(NF, 1,:), VERTEX(NF, 2,:), VERTEX(NF, 3,:), &
            CLOSEST_PT(:))

            DIST(:) = POSITION(:) - CLOSEST_PT(:)
            DISTSQ = DOT_PRODUCT(DIST, DIST)

            IF(DISTSQ .GE. RADSQ) CYCLE !No overlap exists, move on to the next facet
                       
            !Overlap detected
            !Set overlap_exists to true and exit 
            OVERLAP_EXISTS = .true.
            RETURN

         ENDDO
         
      end DO
      
      RETURN 
      
      END SUBROUTINE CHECK_IF_PARTICLE_OVELAPS_STL


      SUBROUTINE CALC_DEM_FORCE_WITH_WALL_STL

      USE run
      USE param1
      USE discretelement
      USE geometry
      USE compar
      USE constant
      USE indices
      USE stl
      USE des_stl_functions
      Implicit none
      
      INTEGER :: LL 
      INTEGER I, J,K, II, IW, IDIM, IJK, NF, wall_count
      DOUBLE PRECISION OVERLAP_N, OVERLAP_T, SQRT_OVERLAP, OVERLAP_PERCENT

      DOUBLE PRECISION V_REL_TRANS_NORM, V_REL_TRANS_TANG, &
      DISTSQ, RADSQ, CLOSEST_PT(DIMN) 
! local normal and tangential forces
      DOUBLE PRECISION FNS1(DIMN), FNS2(DIMN)
      DOUBLE PRECISION FTS1(DIMN), FTS2(DIMN)
      DOUBLE PRECISION NORMAL(DIMN), TANGENT(DIMN), DIST(DIMN), DISTMOD

      LOGICAL :: checked_facet_already,DES_LOC_DEBUG, PARTICLE_SLIDE, &
      test_overlap_and_exit
      INTEGER :: COUNT_FAC, COUNT, COUNT2, & 
      contact_facet_count, NEIGH_CELLS, NEIGH_CELLS_NONNAT, &
      LIST_OF_CELLS(27), CELL_ID, I_CELL, J_CELL, K_CELL, cell_count
      INTEGER :: IMINUS1, IPLUS1, JMINUS1, JPLUS1, KMINUS1, KPLUS1, PHASELL

      DOUBLE PRECISION :: CROSSP(DIMN)
      DOUBLE PRECISION :: FTMD, FNMD
! local values used spring constants and damping coefficients
      DOUBLE PRECISION ETAN_DES_W, ETAT_DES_W, KN_DES_W, KT_DES_W

      !reference point and direction of the line
      double precision, dimension(dimn) :: ref_line,  dir_line 
      !reference point and normal of the plane
      double precision, dimension(dimn) :: ref_plane, norm_plane
      !line is parameterized as p = p_ref + t * dir_line, t is line_param
      double precision :: line_t
      !flag to tell if the orthogonal projection of sphere center to 
      !extended plane detects an overlap 
      LOGICAL :: ortho_proj_cut 
      INTEGER, Parameter :: MAX_FACET_CONTS = 200 
      INTEGER :: list_of_checked_facets(max_facet_conts)

      INCLUDE 'function.inc'

      DES_LOC_DEBUG = .false. ;      DEBUG_DES = .false.
      FOCUS_PARTICLE = -1
      
      DO LL = 1, MAX_PIP 
         
         IF(LL.EQ.FOCUS_PARTICLE) DEBUG_DES = .TRUE.
        
! skipping non-existent particles or ghost particles
         IF(.NOT.PEA(LL,1) .OR. PEA(LL,4)) CYCLE
         
!---------------------------------------------------------------------
! make sure the particle is not classified as a new 'entering' particle 
! or is already marked as a potential exiting particle

         IF( .NOT.PEA(LL,2) .AND. .NOT.PEA(LL,3)) CYCLE

! If no neighboring facet in the surrounding 27 cells, then exit 
         IF (NO_NEIGHBORING_FACET_DES(PIJK(LL,4))) cycle  
              

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
                  IF(COUNT_FAC.EQ.0) CYCLE
                  distsq = zero
                  IF(DES_POS_NEW( LL , 1) > XE(I)) DISTSQ = DISTSQ &
                  + (DES_POS_NEW(LL,1)-XE(I))*(DES_POS_NEW(LL,1)-XE(I))
                  
                  IF(DES_POS_NEW( LL , 1) < XE(I) - DX(I)) DISTSQ = DISTSQ &
                  + (XE(I) - DX(I) - DES_POS_NEW(LL,1))*(XE(I) - DX(I) - DES_POS_NEW(LL,1))
                  
                  IF(DES_POS_NEW( LL , 2) > YN(J)) DISTSQ = DISTSQ &
                  + (DES_POS_NEW(LL,2)-YN(J))* (DES_POS_NEW(LL,2)-YN(J))
                  
                  IF(DES_POS_NEW( LL , 2) < YN(J) - DY(J)) DISTSQ = DISTSQ &
                  + (YN(J) - DY(J) - DES_POS_NEW(LL,2))* (YN(J) - DY(J) - DES_POS_NEW(LL,2))
                    
                  IF(DES_POS_NEW( LL , 3) > ZT(K)) DISTSQ = DISTSQ &
                  + (DES_POS_NEW(LL,3)-ZT(K))*(DES_POS_NEW(LL,3)-ZT(K))
                  
                  IF(DES_POS_NEW( LL , 3) < ZT(K) - DZ(K)) DISTSQ = DISTSQ &
                  + (ZT(K) - DZ(K) - DES_POS_NEW(LL,3))*(ZT(K) - DZ(K) - DES_POS_NEW(LL,3))
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

               !Checking all the facets is time consuming due to the 
               !expensive separating axis test. Remove this facet from 
               !contention based on a simple orthogonal projection test. 
               
               !parametrize a line as p = p_0 + t normal 
               !and intersect with the triangular plane. 
               !if t>0, then point is on the 
               !non-fluid side of the plane, if the plane normal
               !is assumed to point toward the fluid side 
               

               line_t  = Undefined 
               !-undefined, because non zero values will imply the sphere center
               !is on the non-fluid side of the plane. Since the testing 
               !is with extended plane, this could very well happen even 
               !when the particle is well inside the domain (assuming the plane
               !normal points toward the fluid). See the pic below. So check
               !only when line_t is negative
    
!                            \   Solid  /
!                             \  Side  /
!                              \      /        
!                               \    /
!            Wall 1, fluid side  \  /  Wall 2, fluid side
!                                 \/
!                                   o particle 
!                  line_t will be positive for wall 1 (incorrectly indicating center
!                  is outside the domain) 
!                  line_t will be negative for wall 2 
!              
!                Therefore, only stick with this test when line_t is negative and let the
!                separating axis test take care of the other cases. 

               !Assume the orthogonal projection detects an overlap 
               ortho_proj_cut = .true. 
               
               ref_line(1:dimn) = des_pos_new(LL, 1:dimn)
               dir_line(1:dimn) = NORM_FACE(NF,1:dimn)
               !Since this is for checking static config, line's direction
               !is the same as plane's normal. For moving particles, 
               !the line's normal will be along the point joining new 
               !and old positions. 
               
               norm_plane(1:dimn) = NORM_FACE(NF,1:dimn)
               ref_plane(1:dimn)  = VERTEX(NF, 1, 1:dimn)
               CALL intersectLnPlane(ref_line, dir_line, ref_plane, & 
                    norm_plane, line_t) 
               !k - rad >= tol_orth, where k = -line_t, then orthogonal 
               !projection is false. Substituting for k
            !=> line_t + rad <= -tol_orth
               !choosing tol_orth = 0.01% of des_radius = 0.0001*des_radius 
               if(line_t.le.zero.and. & 
                    (line_t+des_radius(LL).le.-0.0001d0*des_radius(LL))) ortho_proj_cut = .false. 
               !Orthogonal projection will detect false postitives even
               !when the particle does not overlap the triangle.
               !However, if the orthgonal projection shows no overlap, then 
               !that is a big fat nagative and overlaps are not possible. 
               if(.not.ortho_proj_cut) cycle 
               
               CALL ClosestPtPointTriangle(DES_POS_NEW(LL,:), &
                    VERTEX(NF, 1,:), VERTEX(NF, 2,:), VERTEX(NF, 3,:), &
                    CLOSEST_PT(:))
            
               DIST(:) = DES_POS_NEW(LL,:) - CLOSEST_PT(:)
               DISTSQ = DOT_PRODUCT(DIST, DIST)
               OVERLAP_N = ZERO
               OVERLAP_PERCENT = ZERO

               IF(DISTSQ .GE. RADSQ) CYCLE !No overlap exists
            
               phaseLL = PIJK(LL,5)
               ! Hertz vs linear spring-dashpot contact model
               IF (DES_COLL_MODEL_ENUM .EQ. HERTZIAN) THEN
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
               
               DISTMOD = SQRT(DISTSQ)
               OVERLAP_N = DES_RADIUS(LL) - DISTMOD
               OVERLAP_PERCENT = (OVERLAP_N/DES_RADIUS(LL))*100.D0
               
               NORMAL(:) = -NORM_FACE(NF,:)
               
               !Calculate the translational relative velocity for a contacting particle pair
               CALL CFRELVEL_WALL2(LL, V_REL_TRANS_NORM, &
                    & V_REL_TRANS_TANG, TANGENT, NORMAL, DISTMOD)
               
               FNS1(:) = -KN_DES_W * OVERLAP_N * NORMAL(:)
               FNS2(:) = -ETAN_DES_W * V_REL_TRANS_NORM * NORMAL(:)
               FNORM(:) = FNS1(:) + FNS2(:)

! Calculate the tangential displacement. The tangential displacement for facet
! treatment is not based on time integration of forces. It is given just as 
! overlap_t = vrel_tang*dtsolid 
            
               OVERLAP_T = V_REL_TRANS_TANG*DTSOLID
               
               FTS1(:) = -KT_DES_W * OVERLAP_T*TANGENT(:)
               FTS2(:) = -ETAT_DES_W * V_REL_TRANS_TANG * TANGENT(:)
               FTAN(:) =  FTS1(:) + FTS2(:)

               
! Check for Coulombs friction law and limit the maximum value of the
! tangential force on a particle in contact with a wall
            
               FTMD = sqrt(dot_product(FTAN, FTAN))   
               FNMD = sqrt(dot_product(FNORM,FNORM))  
               !the Square roots could be removed to further optimize this code

               IF (FTMD.GT.(MEW_W*FNMD)) THEN
                  IF(dot_product(TANGENT,TANGENT).EQ.zero) THEN
                     FTAN(:) =  MEW_W * FNMD * FTAN(:)/FTMD
                  ELSE
                     FTAN(:) = -MEW_W * FNMD * TANGENT(:)
                  ENDIF
               ENDIF
               
               !add the force 
               FC(:,LL) = FC(:,LL) + FNORM(:) + FTAN(:)
               
               !now add the torque             
               !Using particle radius as the moment arm for computing the 
               !torque 
               IF(DIMN.EQ.3) THEN
                  CALL DES_CROSSPRDCT(CROSSP, NORMAL, FTAN)
                  TOW(:,LL) = TOW(:,LL) + DES_RADIUS(LL)*CROSSP(:)
               ELSE
                  CROSSP(1) = NORMAL(1)*FTAN(2) - NORMAL(2)*FTAN(1)
                  TOW(1,LL) = TOW(1,LL) + DES_RADIUS(LL)*CROSSP(1)
               ENDIF
            
            ENDDO

         end DO
         
      end DO
      
      RETURN
    END SUBROUTINE CALC_DEM_FORCE_WITH_WALL_STL

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
      VRN = DOT_PRODUCT(VRELTRANS,NORM)

! slip velocity of the contact point
! Equation (8) in Tsuji et al. 1992
      VSLIP(:) =  VRELTRANS(:) - VRN*NORM(:)

! the magnitude of the tangential vector
      TANMOD = SQRT(DOT_PRODUCT(VSLIP,VSLIP))
      IF(TANMOD.NE.ZERO) THEN
! the unit vector in the tangential direction
         TANGNT(:) = VSLIP(:)/TANMOD
      ELSE
         TANGNT(:) = ZERO
      ENDIF

! tangential component of relative surface velocity (scalar)
      VRT  = DOT_PRODUCT(VRELTRANS,TANGNT)


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


      FTMD = SQRT(DOT_PRODUCT(FTAN, FTAN))
      FNMD = SQRT(DOT_PRODUCT(FNORM,FNORM))

      IF (FTMD.GT.(MEW*FNMD)) THEN
         PARTICLE_SLIDE = .TRUE.
         IF(DOT_PRODUCT(TANGNT,TANGNT).EQ.0) THEN
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
      DOUBLE PRECISION :: DIST_CL

!-----------------------------------------------

      FC(:,L) = FC(:,L) + FNORM(:) + FTAN(:)

! calculate the distance from the particle center to the contact point,
! which is taken as the radical line
! dist_ci+dist_cl=dist_li; dist_ci^2+a^2=ri^2;  dist_cl^2+a^2=rl^2
      DIST_CL = (DIST_LI**2 + DES_RADIUS(L)**2 - DES_RADIUS(II)**2)/&
         (2.d0*DIST_LI)


      IF(DIMN.EQ.3) THEN
         CALL DES_CROSSPRDCT(CROSSP, NORM, FTAN)
         TOW(:,L)  = TOW(:,L)  + DIST_CL*CROSSP(:)
      ELSE
         CROSSP(1) = NORM(1)*FTAN(2) - NORM(2)*FTAN(1)
         TOW(1,L)  = TOW(1,L)  + DIST_CL*CROSSP(1)
      ENDIF


      RETURN
      END SUBROUTINE CFFCTOW2

      SUBROUTINE CFRELVEL_WALL2(L,  VRN, VRT, TANGNT, NORM, DIST_LI)

      USE discretelement
      USE param1
      IMPLICIT NONE

!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER L


      DOUBLE PRECISION TANGNT(DIMN), NORM(DIMN)
      DOUBLE PRECISION TANMOD, VRN, VRT
      DOUBLE PRECISION VRELTRANS(DIMN)
      DOUBLE PRECISION VSLIP(DIMN), &
                       V_ROT(DIMN), OMEGA_SUM(DIMN)

! distance between particles
      DOUBLE PRECISION DIST_LI
! distance from the contact point to the particle centers
      DOUBLE PRECISION DIST_CL, DIST_CI


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
      VRN = dot_product(VRELTRANS,NORM)

! slip velocity of the contact point
! Equation (8) in Tsuji et al. 1992
      VSLIP(:) =  VRELTRANS(:) - VRN*NORM(:)

! the magnitude of the tangential vector
      TANMOD = SQRT(DOT_PRODUCT(VSLIP,VSLIP))
      IF(TANMOD.NE.ZERO) THEN
! the unit vector in the tangential direction
         TANGNT(:) = VSLIP(:)/TANMOD
      ELSE
         TANGNT(:) = ZERO
      ENDIF

! tangential component of relative surface velocity (scalar)
      VRT  = dot_product(VRELTRANS,TANGNT)


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

