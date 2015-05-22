!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Module name: MAKE_ARRAYS_DES                                        !
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  !
!                                                                      !
!  Purpose: DES - allocating DES arrays
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE MAKE_ARRAYS_DES

      USE calc_collision_wall
      USE compar
      USE cutcell
      USE des_rxns
      USE des_stl_functions
      USE des_thermo
      USE desgrid
      USE discretelement
      USE error_manager
      USE functions
      USE funits
      USE geometry
      USE mpi_utility
      USE param1
      USE run
      USE stl
      use mfix_pic, only:  MPPIC
      use mpi_funs_des, only: DES_PAR_EXCHANGE
      use constant, only: PI
      use stl_preproc_des, only: add_facet


      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER :: I, J, K, L, IJK, PC, SM_CELL
      INTEGER :: I1, I2, J1, J2, K1, K2, II, JJ, KK, IJK2
      INTEGER :: lface, lcurpar, lpip_all(0:numpes-1), lglobal_id  , lparcnt
      INTEGER :: CELL_ID, I_CELL, J_CELL, K_CELL, COUNT, NF
      INTEGER :: IMINUS1, IPLUS1, JMINUS1, JPLUS1, KMINUS1, KPLUS1

      INTEGER :: FACTOR

! MPPIC related quantities
      DOUBLE PRECISION :: DTPIC_TMPX, DTPIC_TMPY, DTPIC_TMPZ
      CALL INIT_ERR_MSG("MAKE_ARRAYS_DES")

! Check interpolation input.
      CALL SET_FILTER_DES

! cfassign and des_init_bc called before reading the particle info
      CALL CFASSIGN

      VOL_SURR(:) = ZERO

      ! initialize VOL_SURR array
      DO K = KSTART2, KEND1
         DO J = JSTART2, JEND1
            DO I = ISTART2, IEND1
               IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
               IJK = funijk(I,J,K)
               I1 = I
               I2 = I+1
               J1 = J
               J2 = J+1
               K1 = K
               K2 = merge(K, K+1, NO_K)

! looping over stencil points (node values)
               DO KK = K1, K2
                  DO JJ = J1, J2
                     DO II = I1, I2
                        IF (DEAD_CELL_AT(II,JJ,KK)) CYCLE  ! skip dead cells
                        IJK2 = funijk_map_c(II, JJ, KK)
                        IF(FLUID_AT(IJK2)) VOL_SURR(IJK) = &
                           VOL_SURR(IJK)+VOL(IJK2)
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      VERTEX(1,:,WEST_FACEID) = (/zero, zero, zero/)
      VERTEX(2,:,WEST_FACEID) = (/zero, 2*YLENGTH, zero/)
      VERTEX(3,:,WEST_FACEID) = (/zero, zero, 2*ZLENGTH/)

      VERTEX(1,:,EAST_FACEID) = (/XLENGTH, zero, zero/)
      VERTEX(2,:,EAST_FACEID) = (/XLENGTH, 2*YLENGTH, zero/)
      VERTEX(3,:,EAST_FACEID) = (/XLENGTH, zero, 2*ZLENGTH/)

      VERTEX(1,:,SOUTH_FACEID) = (/zero, zero, zero/)
      VERTEX(2,:,SOUTH_FACEID) = (/2*XLENGTH, zero, zero/)
      VERTEX(3,:,SOUTH_FACEID) = (/zero, zero, 2*ZLENGTH/)

      VERTEX(1,:,NORTH_FACEID) = (/zero, YLENGTH, zero/)
      VERTEX(2,:,NORTH_FACEID) = (/2*XLENGTH, YLENGTH, zero/)
      VERTEX(3,:,NORTH_FACEID) = (/zero, YLENGTH, 2*ZLENGTH/)

      VERTEX(1,:,BOTTOM_FACEID) = (/zero, zero, zero/)
      VERTEX(2,:,BOTTOM_FACEID) = (/2*XLENGTH, zero, zero/)
      VERTEX(3,:,BOTTOM_FACEID) = (/zero, 2*YLENGTH, zero/)

      VERTEX(1,:,TOP_FACEID) = (/zero, zero, ZLENGTH/)
      VERTEX(2,:,TOP_FACEID) = (/2*XLENGTH, zero, ZLENGTH/)
      VERTEX(3,:,TOP_FACEID) = (/zero, 2*YLENGTH, ZLENGTH/)


      NORM_FACE(:,WEST_FACEID) = (/one, zero, zero/)
      NORM_FACE(:,EAST_FACEID) = (/-one, zero, zero/)
      NORM_FACE(:,SOUTH_FACEID) = (/zero, one, zero/)
      NORM_FACE(:,NORTH_FACEID) = (/zero, -one, zero/)
      NORM_FACE(:,BOTTOM_FACEID) = (/zero, zero, one/)
      NORM_FACE(:,TOP_FACEID) = (/zero, zero, -one/)


      STL_FACET_TYPE(WEST_FACEID) = FACET_TYPE_NORMAL
      STL_FACET_TYPE(EAST_FACEID) = FACET_TYPE_NORMAL
      STL_FACET_TYPE(NORTH_FACEID) = FACET_TYPE_NORMAL
      STL_FACET_TYPE(SOUTH_FACEID) = FACET_TYPE_NORMAL
      STL_FACET_TYPE(TOP_FACEID) = FACET_TYPE_NORMAL
      STL_FACET_TYPE(BOTTOM_FACEID) = FACET_TYPE_NORMAL

! initialize CELLNEIGHBOR_FACET array
      DO K_CELL = DG_KSTART2, DG_KEND2
      DO J_CELL = DG_JSTART2, DG_JEND2
      DO I_CELL = DG_ISTART2, DG_IEND2

         IPLUS1  =  MIN (I_CELL + 1, DG_IEND2)
         IMINUS1 =  MAX (I_CELL - 1, DG_ISTART2)

         JPLUS1  =  MIN (J_CELL + 1, DG_JEND2)
         JMINUS1 =  MAX (J_CELL - 1, DG_JSTART2)

         KPLUS1  =  MIN (K_CELL + 1, DG_KEND2)
         KMINUS1 =  MAX (K_CELL - 1, DG_KSTART2)

         CELL_ID = DG_FUNIJK(I_CELL,J_CELL,K_CELL)

         IF(.NOT.DES_PERIODIC_WALLS_X)THEN
            IF(I_CELL == DG_IMIN1) CALL ADD_FACET(CELL_ID,WEST_FACEID)
            IF(I_CELL == DG_IMAX1) CALL ADD_FACET(CELL_ID,EAST_FACEID)
         ENDIF

         IF(.NOT.DES_PERIODIC_WALLS_Y)THEN
            IF(J_CELL == DG_JMIN1) CALL ADD_FACET(CELL_ID,SOUTH_FACEID)
            IF(J_CELL == DG_JMAX1) CALL ADD_FACET(CELL_ID,NORTH_FACEID)
         ENDIF

         IF (DO_K .AND. .NOT.DES_PERIODIC_WALLS_Z)THEN
            IF(K_CELL == DG_KMIN1) CALL ADD_FACET(CELL_ID,BOTTOM_FACEID)
            IF(K_CELL == DG_KMAX1) CALL ADD_FACET(CELL_ID,TOP_FACEID)
         ENDIF

         DO KK = KMINUS1, KPLUS1
         DO JJ = JMINUS1, JPLUS1
         DO II = IMINUS1, IPLUS1
            IJK = DG_FUNIJK(II,JJ,KK)

            DO COUNT = 1, LIST_FACET_AT_DES(IJK)%COUNT_FACETS
               NF = LIST_FACET_AT_DES(IJK)%FACET_LIST(COUNT)
               CALL ADD_FACET(CELL_ID,NF)
            ENDDO
         ENDDO
         ENDDO
         ENDDO
      ENDDO
      ENDDO
      ENDDO

! Set the initial particle data.
      IF(RUN_TYPE == 'NEW') THEN

         IF(PARTICLES /= 0) THEN
            IF(GENER_PART_CONFIG) THEN
               CALL COPY_PARTICLE_CONFIG_FROMLISTS
            ELSE
               CALL READ_PAR_INPUT
            ENDIF
         ENDIF

! Set the global ID for the particles and set the ghost cnt
         ighost_cnt = 0
         lpip_all = 0
         lpip_all(mype) = pip
         call global_all_sum(lpip_all)
         lglobal_id = sum(lpip_all(0:mype-1))
         imax_global_id = 0
         do lcurpar  = 1,pip
            lglobal_id = lglobal_id + 1
            iglobal_id(lcurpar) = lglobal_id
            imax_global_id = iglobal_id(pip)
         end do
         call global_all_max(imax_global_id)

! Initialize old values
         omega_new(:,:)   = zero

! Particle orientation
         IF(PARTICLE_ORIENTATION) THEN
            ORIENTATION(1,:) = INIT_ORIENTATION(1)
            ORIENTATION(2,:) = INIT_ORIENTATION(2)
            ORIENTATION(3,:) = INIT_ORIENTATION(3)
         ENDIF


         IF (DO_OLD) THEN
            omega_old(:,:)   = zero
            des_pos_old(:,:) = des_pos_new(:,:)
            des_vel_old(:,:) = des_vel_new(:,:)
         ENDIF

! Read the restart file.
      ELSEIF(RUN_TYPE == 'RESTART_1' .OR. RUN_TYPE == 'RESTART_2') THEN

         CALL READ_RES0_DES
         imax_global_id = maxval(iglobal_id(1:pip))
         call global_all_max(imax_global_id)

! Initizlie the old values.
         IF (DO_OLD) THEN
            omega_old(:,:)   = omega_new(:,:)
            des_pos_old(:,:) = des_pos_new(:,:)
            des_vel_old(:,:) = des_vel_new(:,:)
         ENDIF
         IF(ENERGY_EQ) DES_T_s_OLD(:) = DES_T_s_NEW(:)

      ELSE

         WRITE(ERR_MSG, 1100)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
 1100 FORMAT('Error 1100: Unsupported RUN_TYPE for DES.')

      ENDIF

      IF(RUN_TYPE == 'RESTART_2') VTP_FINDEX=0

! setting additional particle properties now that the particles
! have been identified
      DO L = 1, MAX_PIP
! Skip 'empty' locations when populating the particle property arrays.
         IF(.NOT.PEA(L,1)) CYCLE
         IF(PEA(L,4)) CYCLE
         PVOL(L) = (4.0D0/3.0D0)*PI*DES_RADIUS(L)**3
         PMASS(L) = PVOL(L)*RO_SOL(L)
         OMOI(L) = 2.5D0/(PMASS(L)*DES_RADIUS(L)**2) !ONE OVER MOI
      ENDDO

      CALL SET_PHASE_INDEX
      CALL INIT_PARTICLES_IN_CELL

! do_nsearch should be set before calling particle in cell
      DO_NSEARCH =.TRUE.
      CALL DES_PAR_EXCHANGE
      CALL PARTICLES_IN_CELL

      IF(DEM_SOLIDS) THEN
         CALL NEIGHBOUR
         CALL INIT_SETTLING_DEM
      ENDIF

! Calculate interpolation weights
      CALL CALC_INTERP_WEIGHTS
! Calculate mean fields using either interpolation or cell averaging.
      CALL COMP_MEAN_FIELDS

      IF(MPPIC) CALL CALC_DTPIC

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE MAKE_ARRAYS_DES
