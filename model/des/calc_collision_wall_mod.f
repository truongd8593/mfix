!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_COLLISION_WALL                                    C
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

      module CALC_COLLISION_WALL

      PRIVATE
      PUBLIC:: CHECK_IF_PARTICLE_OVELAPS_STL, CALC_DEM_FORCE_WITH_WALL_STL, ADD_FACET

      CONTAINS

        subroutine add_facet(cell_id, facet_id)
          Use discretelement
          USE stl
          implicit none
          INTEGER, INTENT(IN) :: cell_id, facet_id

          INTEGER, DIMENSION(:), ALLOCATABLE :: int_tmp
          DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: real_tmp

          INTEGER :: lSIZE1, lSIZE2, ii
          DOUBLE PRECISION :: smallest_extent, min_temp, max_temp

          IF(STL_FACET_TYPE(facet_id).ne.FACET_TYPE_NORMAL) return !Skip this facet

          DO ii = 1, cellneighbor_facet_num(cell_id)
             IF(facet_id .eq. cellneighbor_facet(cell_id)%p(ii)) return
          ENDDO

          cellneighbor_facet_num(cell_id) = cellneighbor_facet_num(cell_id) + 1

          IF(cellneighbor_facet_num(cell_id) > cellneighbor_facet_max(cell_id)) THEN
             cellneighbor_facet_max(cell_id) = 2*cellneighbor_facet_max(cell_id)

             lSIZE2 = size(cellneighbor_facet(cell_id)%p)
             allocate(int_tmp(cellneighbor_facet_max(cell_id)))
             int_tmp(1:lSIZE2) = cellneighbor_facet(cell_id)%p(1:lSIZE2)
             call move_alloc(int_tmp,cellneighbor_facet(cell_id)%p)

             lSIZE2 = size(cellneighbor_facet(cell_id)%extentdir)
             allocate(int_tmp(cellneighbor_facet_max(cell_id)))
             int_tmp(1:lSIZE2) = cellneighbor_facet(cell_id)%extentdir(1:lSIZE2)
             call move_alloc(int_tmp,cellneighbor_facet(cell_id)%extentdir)

             lSIZE2 = size(cellneighbor_facet(cell_id)%extentmin)
             allocate(real_tmp(cellneighbor_facet_max(cell_id)))
             real_tmp(1:lSIZE2) = cellneighbor_facet(cell_id)%extentmin(1:lSIZE2)
             call move_alloc(real_tmp,cellneighbor_facet(cell_id)%extentmin)

             lSIZE2 = size(cellneighbor_facet(cell_id)%extentmax)
             allocate(real_tmp(cellneighbor_facet_max(cell_id)))
             real_tmp(1:lSIZE2) = cellneighbor_facet(cell_id)%extentmax(1:lSIZE2)
             call move_alloc(real_tmp,cellneighbor_facet(cell_id)%extentmax)

          ENDIF

          cellneighbor_facet(cell_id)%p(cellneighbor_facet_num(cell_id)) = facet_id
          smallest_extent = huge(0.0)
          do ii=1,3
             min_temp = minval(VERTEX(:,ii,facet_id))
             max_temp = maxval(VERTEX(:,ii,facet_id))
             if ( abs(max_temp - min_temp) < smallest_extent ) then
                cellneighbor_facet(cell_id)%extentdir(cellneighbor_facet_num(cell_id)) = ii
                cellneighbor_facet(cell_id)%extentmin(cellneighbor_facet_num(cell_id)) = min_temp
                cellneighbor_facet(cell_id)%extentmax(cellneighbor_facet_num(cell_id)) = max_temp
                smallest_extent = abs(max_temp - min_temp)
             endif
          enddo

        end subroutine add_facet

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
      USE functions
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

            CALL ClosestPtPointTriangle(POSITION(:), VERTEX(:,:,NF), CLOSEST_PT(:))

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



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!                                                                      !
!                                                                      !
!                                                                      !
!                                                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_DEM_FORCE_WITH_WALL_STL

      USE run
      USE param1
      USE desgrid
      USE discretelement
      USE geometry
      USE compar
      USE constant
      USE indices
      USE stl
      USE des_stl_functions
      USE functions
      Implicit none

      INTEGER :: LL
      INTEGER I, J,K, II, IW, IDIM, IJK, NF, wall_count
      DOUBLE PRECISION OVERLAP_N, SQRT_OVERLAP

      DOUBLE PRECISION V_REL_TRANS_NORM, DISTSQ, RADSQ, CLOSEST_PT(DIMN)
! local normal and tangential forces
      DOUBLE PRECISION FNS1(DIMN), FNS2(DIMN)
      DOUBLE PRECISION FTS1(DIMN), FTS2(DIMN)
      DOUBLE PRECISION NORMAL(DIMN), TANGENT(DIMN), DIST(DIMN), DISTMOD
      DOUBLE PRECISION, DIMENSION(DIMN) :: FTAN, FNORM, OVERLAP_T

      LOGICAL :: checked_facet_already,DES_LOC_DEBUG, PARTICLE_SLIDE, &
      test_overlap_and_exit
      INTEGER :: COUNT_FAC, COUNT, COUNT2, &
      contact_facet_count, &
      LIST_OF_CELLS(27), CELL_ID, I_CELL, J_CELL, K_CELL, cell_count
      INTEGER :: IMINUS1, IPLUS1, JMINUS1, JPLUS1, KMINUS1, KPLUS1, PHASELL

      DOUBLE PRECISION :: CROSSP(DIMN)
      DOUBLE PRECISION :: FTMD, FNMD
! local values used spring constants and damping coefficients
      DOUBLE PRECISION ETAN_DES_W, ETAT_DES_W, KN_DES_W, KT_DES_W

      double precision :: line_t
      !flag to tell if the orthogonal projection of sphere center to
      !extended plane detects an overlap
      LOGICAL :: ortho_proj_cut
      INTEGER, Parameter :: MAX_FACET_CONTS = 200
      INTEGER :: list_of_checked_facets(max_facet_conts)

      DOUBLE PRECISION :: FORCE_HISTORY(DIMN), DTSOLID_TMP


      DOUBLE PRECISION :: MAX_DISTSQ, DISTAPART, FORCE_COH, R_LM
      INTEGER :: MAX_NF, axis
      DOUBLE PRECISION, DIMENSION(3) :: PARTICLE_MIN, PARTICLE_MAX

      type(facet_linked_list), POINTER :: current_p, previous_p

      DES_LOC_DEBUG = .false. ;      DEBUG_DES = .false.
      FOCUS_PARTICLE = -1

      DO LL = 1, MAX_PIP

         IF(LL.EQ.FOCUS_PARTICLE) DEBUG_DES = .TRUE.

! skipping non-existent particles or ghost particles
         IF(.NOT.PEA(LL,1) .OR. PEA(LL,4)) CYCLE

!---------------------------------------------------------------------
! make sure the particle is not classified as a new 'entering' particle
! or is already marked as a potential exiting particle

         IF( PEA(LL,2) .OR. PEA(LL,3)) CYCLE

! If no neighboring facet in the surrounding 27 cells, then exit
         IF (NO_NEIGHBORING_FACET_DES(DG_PIJK(LL))) cycle
!         IF (NO_NEIGHBORING_FACET_DES(PIJK(LL,4))) cycle

        IF(DEBUG_DES.AND.LL.EQ.FOCUS_PARTICLE) THEN
            IJK = PIJK(LL,4)
            COUNT_FAC = LIST_FACET_AT_DES(IJK)%COUNT_FACETS

            WRITE(*,*) 'NUMBER OF FACETS = ', I_OF(IJK), J_OF(IJK), K_OF(IJK), IJK
            WRITE(*,*) 'NUMBER OF FACETS = ', COUNT_FAC, I_OF(IJK), J_OF(IJK), K_OF(IJK)

            WRITE(*,'(A, 3(2x, g17.8))') 'POS = ', DES_POS_NEW(:, LL)
         ENDIF

         FTS1(:) = ZERO
         FTS2(:) = ZERO
         FNS1(:) = ZERO
         FNS2(:) = ZERO

! Check particle LL for wall contacts

         LIST_OF_CELLS(:) = -1
         CELL_ID = DG_PIJK(LL)
         COUNT_FAC = LIST_FACET_AT_DES(CELL_ID)%COUNT_FACETS
         RADSQ = DES_RADIUS(LL)*DES_RADIUS(LL)

         particle_max(:) = des_pos_new(:, LL) + des_radius(LL)
         particle_min(:) = des_pos_new(:, LL) - des_radius(LL)

         DO CELL_COUNT = 1, cellneighbor_facet_num(cell_id)

            axis = cellneighbor_facet(cell_id)%extentdir(cell_count)

            NF = cellneighbor_facet(cell_id)%p(cell_count)

! Compute particle-particle VDW cohesive short-range forces
            IF(USE_COHESION .AND. VAN_DER_WAALS) THEN

               CALL ClosestPtPointTriangle(DES_POS_NEW(:,LL), VERTEX(:,:,NF), CLOSEST_PT(:))
               DIST(:) = CLOSEST_PT(:) - DES_POS_NEW(:,LL)
               DISTSQ = DOT_PRODUCT(DIST, DIST)
               R_LM = 2*DES_RADIUS(LL)

               IF(DISTSQ < (R_LM+VDW_OUTER_CUTOFF)**2) THEN
                  IF(DISTSQ > (VDW_INNER_CUTOFF+R_LM)**2) THEN
                     DistApart = (SQRT(DISTSQ)-R_LM)
                     FORCE_COH = HAMAKER_CONSTANT * DES_RADIUS(LL) /           &
                          (12d0*DistApart**2) * (Asperities/(Asperities+    &
                          DES_RADIUS(LL)) + ONE/(ONE+Asperities/DistApart)**2 )
                  ELSE
                     FORCE_COH = 2d0 * PI * SURFACE_ENERGY * DES_RADIUS(LL) *  &
                          (Asperities/(Asperities+DES_RADIUS(LL)) + ONE/          &
                          (ONE+Asperities/VDW_INNER_CUTOFF)**2 )
                  ENDIF
                  FC(:,LL) = FC(:,LL) + DIST(:)*FORCE_COH/SQRT(DISTSQ)
               ENDIF
            ENDIF

            if (cellneighbor_facet(cell_id)%extentmin(cell_count) > particle_max(axis)) then
               call remove_collision()
               cycle
            endif

            if (cellneighbor_facet(cell_id)%extentmax(cell_count) < particle_min(axis)) then
               call remove_collision()
               cycle
            endif

               !Recall the facets on the MI plane are re-classified
               !as MI only when the user specifies BC_MI_AS_WALL_FOR_DES
               !as false. The default is to account for MI BC plane
               !as a wall and the facets making this plane are by
               !default classified as normal.

               !Checking all the facets is time consuming due to the
               !expensive separating axis test. Remove this facet from
               !contention based on a simple orthogonal projection test.

               !parametrize a line as p = p_0 + t normal
               !and intersect with the triangular plane.
               !if t>0, then point is on the
               !non-fluid side of the plane, if the plane normal
               !is assumed to point toward the fluid side

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

            !Since this is for checking static config, line's direction
            !is the same as plane's normal. For moving particles,
            !the line's normal will be along the point joining new
            !and old positions.

            line_t = DOT_PRODUCT(VERTEX(1, 1:dimn,NF) - des_pos_new(1:dimn, LL), NORM_FACE(1:dimn,NF))
            !k - rad >= tol_orth, where k = -line_t, then orthogonal
            !projection is false. Substituting for k
            !=> line_t + rad <= -tol_orth
            !choosing tol_orth = 0.01% of des_radius = 0.0001*des_radius

            !Orthogonal projection will detect false positives even
            !when the particle does not overlap the triangle.
            !However, if the orthogonal projection shows no overlap, then
            !that is a big fat negative and overlaps are not possible.
            if((line_t.le.-1.0001d0*des_radius(LL))) then  ! no overlap
               call remove_collision()
               CYCLE
            ENDIF

            CALL ClosestPtPointTriangle(DES_POS_NEW(:,LL), VERTEX(:,:,NF), CLOSEST_PT(:))

            DIST(:) = CLOSEST_PT(:) - DES_POS_NEW(:,LL)
            DISTSQ = DOT_PRODUCT(DIST, DIST)

            IF(DISTSQ .GE. RADSQ) THEN !No overlap exists
               call remove_collision()
               CYCLE
            ENDIF

!               IF(DISTSQ < MAX_DISTSQ)THEN
                  MAX_DISTSQ = DISTSQ
                  MAX_NF = NF
!               ENDIF

!            ENDDO
!         ENDDO

!         IF(MAX_DISTSQ /= UNDEFINED) THEN
! Assign the collision normal based on the facet with the
! largest overlap.
!                  NORMAL(:) = DIST(:)/sqrt(DISTSQ)
                  NORMAL(:) = DIST(:)/max(sqrt(DISTSQ),0.0000001)

                  !NORMAL(:) = -NORM_FACE(:,MAX_NF)
               !facet's normal is correct normal only when the
               !intersection is with the face. When the intersection
               !is with edge or vertex, then the normal is
               !based on closest pt and sphere center. The
               !definition above of the normal is generic enough to
               !account for differences between vertex, edge, and facet.

! Calculate the particle/wall overlap.
               DISTMOD = SQRT(MAX_DISTSQ)
               OVERLAP_N = DES_RADIUS(LL) - DISTMOD

! Calculate the translational relative velocity for a contacting particle pair
               CALL CFRELVEL_WALL2(LL, V_REL_TRANS_NORM, &
                  TANGENT, NORMAL, DISTMOD)

! Calculate the spring model parameters.
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

! Calculate the normal contact force
               FNS1(:) = -KN_DES_W * OVERLAP_N * NORMAL(:)
               FNS2(:) = -ETAN_DES_W * V_REL_TRANS_NORM * NORMAL(:)
               FNORM(:) = FNS1(:) + FNS2(:)


! Calculate the tangential displacement. Note that only the maximum
! wall collision is considered. Therefore, enduring contact can exist
! from facet-to-facet as long as the particle remains in contact with
! one or more facets.

               if (.not. associated(particle_wall_collisions(LL)%pp)) then
                  allocate(particle_wall_collisions(LL)%pp)
                  current_p => particle_wall_collisions(LL)%pp
                  current_p%PFT(:) = ZERO
                  current_p%facet_id = nf
               else
                  current_p => particle_wall_collisions(LL)%pp
                  do while (associated(current_p))
                     if (current_p%facet_id .eq. nf) exit
                     previous_p => current_p
                     current_p => current_p%next
                  enddo

                  if (.not. associated(current_p)) then
                     allocate(current_p)
                     previous_p%next => current_p
                     current_p%PFT(:) = ZERO
                     current_p%facet_id = nf
                     exit
                  endif
               endif

               IF(abs(sum(FORCE_HISTORY)) .gt. small_number) THEN
                  OVERLAP_T(:) = TANGENT(:) * DTSOLID
               ELSE
                  IF(V_REL_TRANS_NORM .GT. ZERO) THEN
                     DTSOLID_TMP = OVERLAP_N/(V_REL_TRANS_NORM)
                  ELSEIF(V_REL_TRANS_NORM .LT. ZERO) THEN
                     DTSOLID_TMP = DTSOLID
                  ELSE
                     DTSOLID_TMP = OVERLAP_N /                         &
                        (V_REL_TRANS_NORM+SMALL_NUMBER)
                  ENDIF
                  OVERLAP_T(:) = TANGENT(:) * MIN(DTSOLID,DTSOLID_TMP)
               ENDIF

! Update the tangential history.
               current_p%PFT(:) = current_p%PFT(:) + OVERLAP_T(:)

               FORCE_HISTORY(:) = current_p%PFT(:) - &
                  DOT_PRODUCT(current_p%PFT(:),NORMAL)*NORMAL(:)

! Calculate the tangential collision force.
               FTS1(:) = -KT_DES_W * FORCE_HISTORY(:)
               FTS2(:) = -ETAT_DES_W * TANGENT(:)
               FTAN(:) =  FTS1(:) + FTS2(:)

               PARTICLE_SLIDE = .FALSE.

! Check for Coulombs friction law and limit the maximum value of the
! tangential force on a particle in contact with a wall.
               FTMD = DOT_PRODUCT(FTAN, FTAN)
               FNMD = DOT_PRODUCT(FNORM,FNORM)
               IF (FTMD.GT.(MEW_W*MEW_W*FNMD)) THEN
                  PARTICLE_SLIDE = .TRUE.
                  IF(all(TANGENT.EQ.zero)) THEN
                     FTAN(:) =  MEW_W * FTAN(:) * SQRT(FNMD/FTMD)
                  ELSE
                     FTAN(:) = -MEW_W * TANGENT(:) * SQRT(FNMD/dot_product(TANGENT,TANGENT))
                  ENDIF
               ENDIF

! Add the collision force to the total forces acting on the particle.
               FC(:,LL) = FC(:,LL) + FNORM(:) + FTAN(:)

! Add the torque: The particle radius is used as the moment arm
               IF(DO_K) THEN
                  CALL DES_CROSSPRDCT(CROSSP, NORMAL, FTAN)
                  TOW(:,LL) = TOW(:,LL) + DES_RADIUS(LL)*CROSSP(:)
               ELSE
                  CROSSP(1) = NORMAL(1)*FTAN(2) - NORMAL(2)*FTAN(1)
                  TOW(1,LL) = TOW(1,LL) + DISTMOD*CROSSP(1)
               ENDIF

! Save the tangential displacement history with the correction of Coulomb's law
               IF (PARTICLE_SLIDE) THEN
! Since FT might be corrected during the call to cfslide, the tangential
! displacement history needs to be changed accordingly
                  current_p%PFT(:) = -( FTAN(:) - FTS2(:) ) / KT_DES_W
               ELSE
                  current_p%PFT(:) = FORCE_HISTORY(:)
               ENDIF

         ENDDO

!         ENDIF
      ENDDO

      RETURN

       contains

         subroutine remove_collision
           implicit none

               if (associated(particle_wall_collisions(LL)%pp)) then

                  current_p => particle_wall_collisions(LL)%pp

                  if (current_p%facet_id .eq. nf) then
                     particle_wall_collisions(LL)%pp => current_p%next
                     deallocate(current_p)
                  else

                     previous_p => current_p
                     current_p => current_p%next

                     do while (associated(current_p))
                        if (current_p%facet_id .eq. nf) exit
                        previous_p => current_p
                        current_p => current_p%next
                     enddo

                     if (associated(current_p)) then
                        previous_p%next => current_p%next
                        deallocate(current_p)
                     endif

                  endif
               endif

         end subroutine remove_collision

    END SUBROUTINE CALC_DEM_FORCE_WITH_WALL_STL

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!                                                                      !
!                                                                      !
!                                                                      !
!                                                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
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
      write (vtp_unit,"(15x,es13.6)") (2*des_radius(pid))
      write(vtp_unit,"(12x,a)") '</DataArray>'

      temp_array = zero
      temp_array(:) = des_vel_new(:,pid)
      write(vtp_unit,"(12x,a,a)") '<DataArray type="Float32" ',&
           'Name="Velocity" NumberOfComponents="3" format="ascii">'
      write (vtp_unit,"(15x,3(es13.6,3x))")&
           ((temp_array(k)),k=1,3)
      write(vtp_unit,"(12x,a,/9x,a)") '</DataArray>','</PointData>'
      ! skip cell data
      write(vtp_unit,"(9x,a)") '<CellData></CellData>'

      temp_array = zero
      temp_array(1:dimn) = des_pos_new(1:dimn, pid)
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

    end SUBROUTINE write_this_facet_and_part

!----------------------------------------------------------------------!
!                                                                      !
!                                                                      !
!                                                                      !
!                                                                      !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE CFRELVEL_WALL2(LL,  VRN, VSLIP, NORM, DIST_LI)

      USE discretelement
      USE param1

      use geometry, only: DO_K
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: LL
      DOUBLE PRECISION, DIMENSION(DIMN), INTENT(IN) :: NORM
      DOUBLE PRECISION, DIMENSION(DIMN), INTENT(OUT):: VSLIP
      DOUBLE PRECISION, INTENT(OUT):: VRN
! distance between particles
      DOUBLE PRECISION, INTENT(IN) :: DIST_LI

!-----------------------------------------------
! Local variables
!-----------------------------------------------

      DOUBLE PRECISION, DIMENSION(DIMN) :: V_ROT, OMEGA_SUM, VRELTRANS

! distance from the contact point to the particle centers
      DOUBLE PRECISION DIST_CL, DIST_CI

! translational relative velocity
      VRELTRANS(:) = DES_VEL_NEW(:,LL)

! rotational contribution  : v_rot
! calculate the distance from the particle center to the wall
      DIST_CL = DIST_LI         !- DES_RADIUS(L)
      IF(DO_K) THEN
         OMEGA_SUM(:) = OMEGA_NEW(:,LL)*DIST_CL
      ELSE
         OMEGA_SUM(1) = OMEGA_NEW(1,LL)*DIST_CL
         OMEGA_SUM(2) = ZERO
         OMEGA_SUM(3) = ZERO
      ENDIF

      CALL DES_CROSSPRDCT(V_ROT, OMEGA_SUM, NORM)

! total relative velocity
      VRELTRANS(:) =  DES_VEL_NEW(:,LL) + V_ROT(:)

! normal component of relative velocity (scalar)
      VRN = DOT_PRODUCT(VRELTRANS,NORM)

! slip velocity of the contact point
! Equation (8) in Tsuji et al. 1992
      VSLIP(:) =  VRELTRANS(:) - VRN*NORM(:)

      RETURN
      END SUBROUTINE CFRELVEL_WALL2

    end module CALC_COLLISION_WALL

