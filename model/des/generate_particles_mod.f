MODULE GENERATE_PARTICLES

      TYPE PARTICLE
      INTEGER :: CELL(4)
      LOGICAL :: INDOMAIN
! Solid phase
      INTEGER :: M
      double precision :: RAD, DENS, STATWT
      double precision :: VELOCITY(3), POSITION(3)

      !could be made allocatable later to reduce memory usage for 2-D runs
      TYPE (PARTICLE), POINTER :: NEXT=>NULL() ! NEXT PARTICLE ADDRESS
      TYPE (PARTICLE), POINTER :: PREV=>NULL() ! PREVIOUS PARTICLE ADDRESS

      END TYPE PARTICLE

! This is the linked list of first set of particles
      TYPE (PARTICLE), POINTER :: ORIG_PART_LIST => NULL()

CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  SUBROUTINE: GENERATE_PARTICLE_CONFIG                                C
!                                                                      C
!  Purpose: Generate particle configuration based on maximum particle  C
!           radius and filling from top to bottom within specified     C
!           bounds                                                     C
!                                                                      C
!                                                                      C
!  Authors: Rahul Garg                                Date: 19-Mar-14  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE GENERATE_PARTICLE_CONFIG

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE mfix_pic, only: MPPIC
      USE discretelement, only: particles

      USE mpi_utility

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager

      IMPLICIT NONE
      INTEGER :: IPE

      INTEGER :: LORIG_ALL(0:Numpes-1)
      INTEGER :: LDEL_ALL(0:Numpes-1), LREM_ALL(0:Numpes-1)

! Initialize the error manager.
      CALL INIT_ERR_MSG("Generate_Particle_Config")

      IF(MPPIC) THEN
         CALL GENERATE_PARTICLE_CONFIG_MPPIC(LORIG_ALL, LREM_ALL, LDEL_ALL)
      ELSE
         CALL GENERATE_PARTICLE_CONFIG_DEM(LORIG_ALL, LREM_ALL, LDEL_ALL)
      ENDIF

      write(err_msg, '(/,70("-"),/,A)') 'Particle statistics reporting'

      call flush_err_msg(header = .false., footer = .false.)

      !total number of particles used to allocate the arrays is set equal to sum
      !of remaining particles on all processors
      PARTICLES  = SUM(LREM_ALL(0:numPEs-1))

      Do ipe = 0, numPEs - 1
         write(err_msg, 1002) ipe, Lorig_all(ipe), Ldel_all(ipe), Lrem_all(ipe)
 1002    FORMAT(/, &
         'For proc', 2x, i5, / &
         'Total number of particles originally seeded       : ', 2x, I15, / &
         'Total number of particles found outside the domain: ', 2x, I15, / &
         'Total number of particles remaining               : ', 2x, I15, /)

         CALL FLUSH_ERR_MSG (header = .false., Footer = .false.)

         if( Lrem_all(ipe) + Ldel_all(ipe) .ne. Lorig_all(ipe)) then

            write(err_msg, 1003) ipe, Lrem_all(ipe) + Ldel_all(ipe), Lorig_all(ipe)
            CALL FLUSH_ERR_MSG (ABORT = .true.)
         endif

 1003    Format('Sanity check failed for particle seeding on proc:', 2x, i5, / &
         '# of particles deleted + # of particles remaining ',  2x, I15, / &
         'not equal to original number of particles',  2x, I15, / &
         'MFIX will now exit')

      ENDDO

      write(err_msg, 1004) PARTICLES
 1004 FORMAT(/, &
      'Total number of particles in the system: ', 2x, i15)

      CALL FLUSH_ERR_MSG (header = .false.)

      CALL FINL_ERR_MSG

      END SUBROUTINE GENERATE_PARTICLE_CONFIG

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  SUBROUTINE: GENERATE_PARTICLE_CONFIG                                C
!                                                                      C
!  Purpose: Generate particle configuration for DEM solids for each IC C
!           region. Now using the particle linked lists for initial    C
!           build                                                      C
!           This routine will ultimately supersede the older rouine    C
!           that has not been deleted yet                              C
!                                                                      C
!  Authors: Rahul Garg                               Date: 21-Mar-2014 C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE GENERATE_PARTICLE_CONFIG_DEM(LORIG_ALL, LREM_ALL, LDEL_ALL)


! Global Variables:
! particle radius and density
      USE discretelement, only: DES_RADIUS, RO_Sol
! particle position new and old
      USE discretelement, only: des_pos_new, des_pos_old
! particle velocity new and old
      USE discretelement, only: des_vel_new, des_vel_old
! Simulation dimension (2D/3D)
      USE discretelement, only: DIMN
! X , Y, and Z position of cell faces of computational fluid grid
      USE discretelement, only: XE, YN, ZT
! Number of particles in the system (current)
      USE discretelement, only:  PIP
! Number of DEM solids phases.
      USE discretelement, only: DES_MMAX
! DEM solid phase diameters and densities.
      USE discretelement, only: DES_D_p0, DES_RO_s
! Number of particles estimated to be seeded, per phase in each IC region
      USE discretelement, only: PART_MPHASE_BYIC
! Number of particles actually seeded, per phase in each IC region
      USE discretelement, only: REALPART_MPHASE_BYIC

      USE discretelement, only: DO_OLD, OMEGA_OLD, OMEGA_NEW, PIJK

! Constant: 3.14159...

      USE discretelement, only:DES_GETINDEXFROMPOS
      USE constant, only: PI
! Flag indicating that the IC region is defined.
      USE ic, only: IC_DEFINED
! IC Region solids volume fraction.
      USE ic, only: IC_EP_S
! min and max physical co-ordinates of IC regions in each direction
      USE ic, only: IC_X_w, IC_X_e, IC_Y_s, IC_Y_n, IC_Z_b, IC_Z_t
! initally specified velocity field and granular temperature
      USE ic, only: IC_U_s, IC_V_s, IC_W_s, IC_Theta_M
! Flag to extend the lattice distribution in a given IC to available area
      Use ic, only: IC_DES_FIT_TO_REGION
! Parameter for detecting unspecified values, zero, and one
      USE param1, only: UNDEFINED, UNDEFINED_I, ZERO, ONE, Half
! Parameter for small and large numbers
      USE param1, only: SMALL_NUMBER, LARGE_NUMBER
! Maximum number of initial conditions
      USE param, only: DIMENSION_IC
! first and last index of the physical cells in regular MFIX grid
      USE compar, only: istart1, iend1, jstart1, jend1, kstart1, kend1
! to access random number generator subroutines
      USE randomno
      USE mpi_utility
      USE discretelement, only: SET_NORMAL

      USE error_manager

! direction wise spans of the domain and grid spacing in each direction
      USE geometry, only: xlength, ylength, zlength

      USE functions
      USE cutcell, only : CARTESIAN_GRID, CUT_CELL_AT
      USE STL_PREPROC_DES, only: CHECK_IF_PARTICLE_OVERLAPS_STL
      USE run, only: solids_model
      USE mfix_pic, only: des_stat_wt
      use des_allocate, only: PARTICLE_GROW
      IMPLICIT NONE

      INTEGER, INTENT(OUT) :: LORIG_ALL(0:Numpes-1)
      INTEGER, INTENT(OUT) :: LDEL_ALL(0:Numpes-1), LREM_ALL(0:Numpes-1)

!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER :: M, ICV, I, J, K, idim, IJK, II, JJ, KK
      INTEGER :: lproc_parcount, pcount_byic_byphase(dimension_ic, DES_MMAX)
      INTEGER :: seed_x, seed_y, seed_z
      INTEGER :: TOTAL_PARTS_IC, TMP_PART_COUNT_INTHIS_IC
      double precision lmax_dia,lfac,xp,yp,zp
      double precision :: XSTART_IC, YSTART_IC, ZSTART_IC, adj_dia, ep_sm
      double precision :: XEND_IC, YEND_IC, ZEND_IC
      double precision :: xinit, yinit, zinit, ymax , doml(3)
      !max point where the particle center cud be placed
      double precision ::  max_ic_pt(3)
      !factor to re-fit the configuration to the IC region size
      double precision :: refit_fac(3), volijk
      !particles min and max positions
      double precision :: PART_CEN_MIN(DIMN), PART_CEN_MAX(DIMN)

      !Particle mean velocity and standard deviation
      double precision :: vel_mean(3), vel_sig(3)

      double precision, dimension(:,:), allocatable :: pvel_temp

      type(particle) :: part

      INTEGER :: count_part, phase
      LOGICAL :: DELETE_PART

      CALL INIT_ERR_MSG("GENERATE_PARTICLE_CONFIG_DEM")

! initializing particle count
      lproc_parcount = 0

      PCOUNT_BYIC_BYPHASE(:,:)  = 0
! setting particle seed spacing grid to be slightly greater than
! the maximum particle diameter. seed at ~particle radii
      lfac = 1.05d0

      LORIG_ALL = 0
      LDEL_ALL = 0
      LREM_ALL = 0
      count_part = 0

      if(.not.allocated(part_mphase_byic)) &
           ALLOCATE(PART_MPHASE_BYIC(DIMENSION_IC, DES_MMAX))
      IF(.not.allocated(REALPART_MPHASE_BYIC)) &
           ALLOCATE(REALPART_MPHASE_BYIC(DIMENSION_IC, DES_MMAX))

      write(err_msg, 2015)
      CALL FLUSH_ERR_MSG(header = .false., FOOTER = .false.)

 2015 FORMAT('IC region wise report on particle initialization')

      IC_LOOP : DO ICV = 1, DIMENSION_IC

         PART_MPHASE_BYIC(ICV,:) = 0
         REALPART_MPHASE_BYIC(ICV,:) = 0
         IF (.not.IC_DEFINED(ICV)) cycle

         PART_CEN_MIN(:) = LARGE_NUMBER
         PART_CEN_MAX(:) = SMALL_NUMBER

         TOTAL_PARTS_IC = 0

         LMAX_DIA = SMALL_NUMBER
         DO M = 1, DES_MMAX
            DOML(1) = IC_X_E(ICV) - IC_X_W(ICV)
            DOML(2) = IC_Y_N(ICV) - IC_Y_S(ICV)
            DOML(3) = IC_Z_T(ICV) - IC_Z_B(ICV)

            IF(NO_K) DOML(3) = ZLENGTH

            VOLIJK = DOML(1)*DOML(2)*DOML(3)

            PART_MPHASE_BYIC(ICV, M) = FLOOR((6.D0*IC_EP_S(ICV, M)*VOLIJK)/&
                 (PI*(DES_D_P0(M)**3)))
            ! setting a local value of maximum diameter for each IC region
            if(PART_MPHASE_BYIC(ICV,M).gt.0) then
               total_parts_ic = total_parts_ic + PART_MPHASE_BYIC(ICV,M)
               LMAX_DIA = MAX(LMAX_DIA, DES_D_P0(M))
            endif
         ENDDO

         IF(total_parts_ic.eq.0) cycle IC_LOOP

         WRITE(ERR_MSG,2016) ICV

2016     FORMAT(/1X,70('-')/, 5x, &
              'IC number         = ', I4)

         CALL FLUSH_ERR_MSG(HEADER = .false., Footer = .false.)


         XSTART_IC = IC_X_W(ICV)
         YSTART_IC = IC_Y_S(ICV)
         ZSTART_IC = IC_Z_B(ICV)
         XEND_IC   = IC_X_E(ICV)
         YEND_IC   = IC_Y_N(ICV)
         ZEND_IC   = IC_Z_T(ICV)

         !Max location of any particle in this IC
         MAX_IC_PT(1) = XEND_IC - 0.5d0*LMAX_DIA*LFAC
         MAX_IC_PT(2) = YEND_IC - 0.5d0*LMAX_DIA*LFAC
         MAX_IC_PT(3) = ZEND_IC - 0.5d0*LMAX_DIA*LFAC
         XINIT = XSTART_IC
         YINIT = YSTART_IC
         ZINIT = ZSTART_IC

         DO M = 1, DES_MMAX
            if(PART_MPHASE_BYIC(ICV,M).eq.0) cycle

            VEL_MEAN(1) = IC_U_S(ICV, M)
            VEL_MEAN(2) = IC_V_S(ICV, M)
            IF(DO_K) VEL_MEAN(3) = IC_W_S(ICV, M)
            !granular temp is defined as (Variance uprime + Variance vprime + Variance wprime)/3
            !assuming equal energy in each direction
            !Variance uprime  = IC_Theta
            !Stdev (or sigma) = sqrt(Variance)

            VEL_SIG(:) = sqrt(IC_Theta_M(ICV, M))
            write(ERR_MSG,2022) M,  &
                 vel_mean(:), IC_theta_m(ICV, M), vel_sig(:)
            CALL FLUSH_ERR_MSG(HEADER = .false., Footer = .false.)


 2022       FORMAT(1X,70('.'),/5x, &
            'PHASE INDEX, M                              =  ', I5,2X, /5x, &
            'INITIALIZING SOLIDS VELOCITY FIELD', /5x, &
            'Mean velocity direction wise                =  ', 3(G15.8,2X), /5x, &
            'Use specified initial granular temperature  =  ', (G15.8,2X), /5x, &
            'Velocity standard deviation direction wise  =  ', 3(G15.8,2X))

            allocate(pvel_temp( PART_MPHASE_BYIC(ICV,M) , DIMN))

            do IDIM = 1, DIMN
               pvel_temp(:, idim) = vel_mean(idim)

               IF(vel_sig(idim).gt.zero) then
                  CALL nor_rno(pvel_temp(1:PART_MPHASE_BYIC(ICV,M),IDIM), &
                       vel_mean(idim),vel_sig(idim))
               ENDIF
            ENDDO

            SEED_X = 1
            SEED_Y = 1
            SEED_Z = 1

            ADJ_DIA = LFAC*DES_D_P0(M)

            SEED_X = FLOOR((XEND_IC - XINIT)/ADJ_DIA)
            SEED_Y = FLOOR((YEND_IC - YINIT)/ADJ_DIA)
            SEED_Z = FLOOR((ZEND_IC - ZINIT)/ADJ_DIA)
            !write(*,*) 'adj_dia = ', adj_dia, lfac, lmax_dia
            !write(*,*) 'seedx  = ', seed_x, seed_y, seed_z
            if(NO_K) seed_z = 1

            TMP_PART_COUNT_INTHIS_IC = (PCOUNT_BYIC_BYPHASE(ICV,M)) + 1

            DO  J = 1, SEED_Y
               DO  K = 1, SEED_Z
                  DO  I = 1, SEED_X
                     XP = XINIT + I*ADJ_DIA - DES_D_P0(M)*0.5D0
                     YP = YINIT + J*ADJ_DIA - DES_D_P0(M)*0.5D0
                     ZP = ZINIT + K*ADJ_DIA - DES_D_P0(M)*0.5D0

                     TMP_PART_COUNT_INTHIS_IC = (TMP_PART_COUNT_INTHIS_IC) + 1

                     IF(TMP_PART_COUNT_INTHIS_IC.GT.PART_MPHASE_BYIC(ICV,M)) EXIT

                     PART_CEN_MIN(1)  = MIN(XP , PART_CEN_MIN(1))
                     PART_CEN_MIN(2)  = MIN(YP , PART_CEN_MIN(2))

                     PART_CEN_MAX(1)  = MAX(XP , PART_CEN_MAX(1))
                     PART_CEN_MAX(2)  = MAX(YP , PART_CEN_MAX(2))
                     IF(DO_K) THEN
                        PART_CEN_MIN(3)  = MIN(ZP, PART_CEN_MIN(3))
                        PART_CEN_MAX(3)  = MAX(ZP, PART_CEN_MAX(3))
                     ENDIF

                  end DO
               end DO
            end DO

         refit_fac = 1.d0
         IF(IC_DES_FIT_TO_REGION(ICV)) THEN

            DO IDIM = 1, DIMN
               IF((PART_CEN_MAX(IDIM)-PART_CEN_MIN(IDIM).GT.LMAX_DIA).AND. &
                    (MAX_IC_PT(IDIM) - PART_CEN_MAX(IDIM).GT.LMAX_DIA)) THEN

                  REFIT_FAC(IDIM)  = MAX_IC_PT(IDIM)/PART_CEN_MAX(IDIM)
                  !write(*,*) ' REFI, IDIM =', IDIM, REFIT_FAC(IDIM)
               END IF
            END DO

            write(err_msg, 2020) ICV, refit_fac(1:dimn)
            CALL FLUSH_ERR_MSG(HEADER = .false., Footer = .false.)

 2020       Format(/5x, 'Refitting to box for IC', I4, /5x,   &
            'Refitting factors (1-DIMN): ', 3(2x,g17.8))

         END IF              !DES_IC_FIT_TO_REGION

            outerloop :  DO  JJ = 1, SEED_Y
               DO  KK = 1, SEED_Z
                  DO  II = 1, SEED_X
                     XP = XINIT + II*ADJ_DIA - DES_D_P0(M)*0.5D0
                     YP = YINIT + JJ*ADJ_DIA - DES_D_P0(M)*0.5D0
                     ZP = ZINIT + KK*ADJ_DIA - DES_D_P0(M)*0.5D0

                     TMP_PART_COUNT_INTHIS_IC = (PCOUNT_BYIC_BYPHASE(ICV,M)) + 1

                     IF(TMP_PART_COUNT_INTHIS_IC.GT.PART_MPHASE_BYIC(ICV,M)) EXIT outerloop

                     !Associate this new particle to the solid phase id based on the map_to_proc defined
                     !earlier

                     PCOUNT_BYIC_BYPHASE(ICV,M)  = PCOUNT_BYIC_BYPHASE(ICV,M) + 1

                     PART%RAD = DES_D_P0(M)*HALF
                     PART%DENS = DES_RO_S(M)
                     PART%POSITION(1) = XP
                     PART%POSITION(2) = YP
                     IF(DO_K) PART%POSITION(3) = ZP
                     PART%VELOCITY(1:DIMN) = PVEL_TEMP(PCOUNT_BYIC_BYPHASE(ICV,M),1:DIMN)

                     do idim = 1, dimn
                        PART%position(idim) = &
                             PART%position(idim)*REFIT_FAC(IDIM)
                     ENDDO

                     I = DES_GETINDEXFROMPOS(ISTART1, IEND1, &
                          PART%position(1), XE(1:size(XE,1)-1),'X','I')
                     !-1 above since xe ranges from 0:IMAX3, so size is imax3 + 1.
                     !therefore, (1:size(xe,1)) will give 1:imax3 + 1, resulting in a seg error.

                     J = DES_GETINDEXFROMPOS(JSTART1,JEND1, &
                          PART%position(2),YN(1:size(YN,1)-1),'Y','J')

                     K = 1
                     IF(DO_K) K = DES_GETINDEXFROMPOS(KSTART1, &
                          KEND1,PART%position(3),ZT(1:size(ZT,1)-1),'Z','K')

                     IF(.not.IS_ON_myPE_wobnd(I,J,K)) cycle

                     IF(DEAD_CELL_AT(I,J,K)) cycle

                     IJK = FUNIJK(I, J, K)
                     IF(.not.FLUID_AT(IJK)) cycle

                     PART%M = M
                     PART%INDOMAIN = .true.

                     PART%STATWT = 1.d0

                     IF(CARTESIAN_GRID) then
! Only look for particles that are in domain. Particles in dead cells
! have already been flagged as outside of the domain in
!  BIN_PARTICLES_TO_CELL.
! Since checking if a particle is on other side of a triangle is tricky,
! for safety, initially remove all the particles in the cut-cell.
! now cut-cell and actualy geometry could be off, so this might not work
! very robustly.
                        IF(.not.CUT_CELL_AT(IJK)) THEN
! Check if any of the particles are overlapping with the walls. This
! will include the normal ghost cell walls and also the cut cells.
                           DELETE_PART = .false.
                           CALL CHECK_IF_PARTICLE_OVERLAPS_STL(PART%POSITION(:), &
                                PART%RAD, DELETE_PART)
                           IF(DELETE_PART) PART%INDOMAIN = .FALSE.
                        ELSE
                           PART%INDOMAIN = .FALSE.
                        ENDIF
                     ENDIF

                     LORIG_ALL(mype) = LORIG_ALL(mype) + 1

                     IF(part%INDOMAIN) then
                        count_part = count_part + 1
                        CALL PARTICLE_GROW(count_part)
                        des_pos_new(1:dimn, count_part) = part%position(1:dimn)
                        des_vel_new(1:dimn, count_part) = part%velocity(1:dimn)

                        DES_RADIUS(count_part) = part%rad
                        RO_SOL(count_part)  = part%dens

                        PIJK(count_part,5) = part%M
                        phase  = part%M
                        IF (DO_OLD) THEN
                           des_vel_old(:, count_part) = des_vel_new(:, count_part)
                           des_pos_old(:, count_part) = des_pos_new(:, count_part)
                        ENDIF

                        CALL SET_NORMAL(count_part)

                        IF(SOLIDS_MODEL(phase).eq.'PIC') then
                           DES_STAT_WT(count_part) = part%STATWT
                        ELSEIF(SOLIDS_MODEL(phase).eq.'DEM') then
                           IF (DO_OLD) OMEGA_OLD(:, count_part) = zero
                           OMEGA_NEW(:, count_part) = zero
                        ENDIF

                        LREM_ALL(mype) = LREM_ALL(mype) + 1
                     ELSE
                        LDEL_ALL(mype) = LDEL_ALL(mype) + 1
                     ENDIF

                     REALPART_MPHASE_BYIC(ICV, part%M) =  &
                          REALPART_MPHASE_BYIC(ICV, part%M) + 1

                     YMAX = YP + DES_D_P0(M)*0.5D0

                  end DO
               end DO
            end DO outerloop

            YINIT = YMAX

            DEALLOCATE(pvel_temp)

            EP_SM = IC_EP_S(ICV,M)
            WRITE(ERR_MSG,2017) EP_SM
            CALL FLUSH_ERR_MSG(HEADER = .false., Footer = .false.)
 2017       FORMAT(5x, 'Solids vol fraction for M   =  ', (G15.8,2X))

         end DO                 ! DO M = 1, DES_MMAX

         CALL global_all_sum(REALPART_MPHASE_BYIC(ICV, 1:DES_MMAX))

         DO M = 1, DES_MMAX
            WRITE(ERR_MSG,'(//, 70("-"),/5x, A, I5,/5x,A, I15, /5x, A, I15)') &
                 'For Phase M:              ', M,  &
                 '# of particles estimated      :',  PART_MPHASE_BYIC(ICV, M), &
                 '# of particles acutally seeded:', INT(REALPART_MPHASE_BYIC(ICV, M))

            CALL FLUSH_ERR_MSG(header = .false. ,footer = .false.)
         END DO

      end DO IC_LOOP

      CALL global_all_sum(LORIG_ALL)
      CALL global_all_sum(LREM_ALL)
      CALL global_all_sum(LDEL_ALL)

      PIP = LREM_ALL(mype) !Setting pip on each processor equal to remaining paticles

      if(count_part.ne.PIP) then
         write(err_msg, 1000) count_part, PIP
         call flush_err_msg(abort = .true.)
      endif

 1000 format('Particles copied from linked lists:', 2x, i15, /, &
      'not equal to particles earlier calculated to be inside domain', 2x, i15, / &
      'This should not have happened, exiting')
      CALL FINL_ERR_MSG

      RETURN
    END SUBROUTINE GENERATE_PARTICLE_CONFIG_DEM

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: GENERATE_PARTICLE_CONFIG_MMPPIC                         C
!  Purpose: generates particle position distribution for MPPIC         C
!  Author: Rahul Garg                                 Date: 3-May-2011 C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE GENERATE_PARTICLE_CONFIG_MPPIC(LORIG_ALL, LREM_ALL, LDEL_ALL)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE cutcell, only : CARTESIAN_GRID
      USE discretelement, only: dimn, xe, yn, zt

      USE discretelement, only:DES_GETINDEXFROMPOS
! Number of DES solids phases.
      USE discretelement, only: DES_MMAX

! DEM solid phase diameters and densities.
      USE discretelement, only: DES_D_p0, DES_RO_s
! Number of particles seeded, per phase in each IC region
      USE discretelement, only: PART_MPHASE_BYIC
! Number of implied real particles by phase by ic region
      USE discretelement, only: REALPART_MPHASE_BYIC

! Global Number of particles seeded per phase
      USE discretelement, only: PART_MPHASE
! Global Number of implied real particles
      USE discretelement, only: REALPART_MPHASE
      USE discretelement

! Implied total number of physical particles
      USE mfix_pic, only: rnp_pic
! Total number of computational particles/parcels
      USE mfix_pic, only: cnp_pic

      USE mfix_pic, only: des_stat_wt

! Flag indicating that the IC region is defined.
      USE ic, only: IC_DEFINED
! IC Region bulk density (RO_s * EP_s)
      USE ic, only: IC_ROP_s
! IC Region gas volume fraction.
      USE ic, only: IC_EP_G

      USE ic, only: IC_X_w, IC_X_e, IC_Y_s, IC_Y_n, IC_Z_b, IC_Z_t
      USE ic, only: IC_I_w, IC_I_e, IC_J_s, IC_J_n, IC_K_b, IC_K_t

! initally specified velocity field and granular temperature
      USE ic, only: IC_U_s, IC_V_s, IC_W_s, IC_Theta_M

! MPPIC specific IC region specification.
      USE ic, only: IC_PIC_CONST_NPC, IC_PIC_CONST_STATWT

! Cut_cell identifier array
      USE cutcell, only: cut_cell_at

! Maximum number of IC regions and solids phases
      USE param, only: DIMENSION_IC
      USE param, only: DIM_M
      USE param1, only: UNDEFINED, UNDEFINED_I, ZERO, ONE, HALF

! Constant: 3.14159...
      USE constant, only: PI

      USE mpi_utility

      USE cutcell, only : CARTESIAN_GRID
      USE randomno
      USE error_manager
      USE functions
      USE STL_PREPROC_DES, only: CHECK_IF_PARTICLE_OVERLAPS_STL
      USE run, only: solids_model
      USE des_allocate, only: PARTICLE_GROW
      IMPLICIT NONE

      INTEGER, INTENT(OUT) :: LORIG_ALL(0:Numpes-1)
      INTEGER, INTENT(OUT) :: LDEL_ALL(0:Numpes-1), LREM_ALL(0:Numpes-1)

!-----------------------------------------------
! Local variables
!-----------------------------------------------
      DOUBLE PRECISION :: EP_SM
! Temp logical variables for checking constant npc and statwt specification
      LOGICAL :: CONST_NPC, CONST_STATWT

      INTEGER :: I, J, K, IJK, ICV, M, COUNT_IC, IDIM

! Number of real and comp. particles in a cell.
      DOUBLE PRECISION ::  REAL_PARTS(DIM_M)
      INTEGER :: COMP_PARTS(DIM_M), IPCOUNT
      DOUBLE PRECISION :: DOML(3), CORD_START(3) , CORD_END(3)
      DOUBLE PRECISION :: STAT_WT

      DOUBLE PRECISION :: RAD, DENS, POSITION(DIMN), VELOCITY(DIMN)
      DOUBLE PRECISION :: VOLIJK, VOLIJK_UNCUT
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: RANDPOS
      double precision :: vel_mean(3), vel_sig(3)

      LOGICAL :: DELETE_PART

      INTEGER :: count_part, phase

      double precision, dimension(:,:), allocatable :: pvel_temp

      type(particle), pointer :: part_list_byic
      type(particle), pointer :: part => null()

      CALL INIT_ERR_MSG("GENERATE_PARTICLE_CONFIG_MPPIC")
      !ALLOCATE(CNP_ARRAY(DIM3_TEMP, DES_MMAX))

      ALLOCATE(RNP_PIC(DES_MMAX))
      ALLOCATE(CNP_PIC(DES_MMAX))

      if(.not.allocated(part_mphase_byic)) &
           ALLOCATE(PART_MPHASE_BYIC(DIMENSION_IC, DES_MMAX))
      IF(.not.allocated(REALPART_MPHASE_BYIC)) &
           ALLOCATE(REALPART_MPHASE_BYIC(DIMENSION_IC, DES_MMAX))


      RNP_PIC(:) = ZERO
      CNP_PIC(:) = ZERO

      COUNT_IC = 0

      write(ERR_MSG,'(/,A,/,A,/)') 'Generating particle configuration for MPPIC model', &
           'Per IC, per Solid phase reporting follows'
      CALL FLUSH_ERR_MSG(header = .false. ,footer = .false.)

      DO ICV = 1, DIMENSION_IC
         PART_MPHASE_BYIC(ICV,:) = 0
         REALPART_MPHASE_BYIC(ICV,:) = 0

         IF(associated(part_list_byic)) Nullify(part_list_byic)

         IF(.not.ic_defined(icv)) cycle

         IF (IC_EP_G(ICV).eq.ONE) cycle

         COUNT_IC = COUNT_IC + 1


         MLOOP: DO M = 1, DES_MMAX
            COMP_PARTS(M) = 0
            REAL_PARTS(M) = zero

            EP_SM = IC_ROP_S(ICV,M)/DES_RO_s(M)
            IF(EP_SM.eq.zero) cycle
            write(ERR_MSG,'(/,70("-"), /, 5x, A, i4, 2x, A, i2)') 'Generating for IC#', ICV, ' and phase', M
            CALL FLUSH_ERR_MSG(header = .false., footer = .false.)

            WRITE(ERR_MSG,2017) DES_D_P0(M), EP_SM
            CALL FLUSH_ERR_MSG(HEADER = .false., FOOTER = .false.)

2017        FORMAT(/, 5x, &
                 'Diameter                    =  ', (ES15.7,2X), /5x, &
                 'Solids vol fraction for M   =  ', (ES15.7,2X), /5x)

            VEL_MEAN = ZERO
            VEL_SIG = ZERO
            VEL_MEAN(1) = IC_U_S(ICV, M)
            VEL_MEAN(2) = IC_V_S(ICV, M)
            IF(DO_K) VEL_MEAN(3) = IC_W_S(ICV, M)
            !granular temp is defined as (Variance uprime + Variance vprime + Variance wprime)/3
            !assuming equal energy in each direction
            !Variance uprime  = IC_Theta
            !Stdev (or sigma) = sqrt(Variance)

            VEL_SIG(:) = sqrt(IC_Theta_M(ICV, M))

            write(ERR_MSG,2022) &
                 vel_mean(:), IC_theta_m(ICV, M), vel_sig(:)
            CALL FLUSH_ERR_MSG(HEADER = .false., Footer = .false.)


2022        FORMAT(/,10x, 'Initializing solids velocity field', /10x, &
                 'Mean velocity direction wise                =  ', 3(G15.8,2X), /10x, &
                 'Use specified initial granular temperature  =  ', (G15.8,2X), /10x, &
                 'Velocity standard deviation direction wise  =  ', 3(G15.8,2X))

            CONST_NPC    = (IC_PIC_CONST_NPC   (ICV, M) .ne. 0) &
                 .AND. (EP_SM.gt.Zero)

            CONST_STATWT = (IC_PIC_CONST_STATWT(ICV, M) .ne. ZERO  ) &
                 .AND. (EP_SM.gt.Zero)

            IF(CONST_NPC) then
               !Estimate the number of parcels
               COMP_PARTS(M) = IC_PIC_CONST_NPC(ICV,M)*(IC_K_T(ICV) - IC_K_B(ICV) +1)* &
                    (IC_J_N(ICV) - IC_J_S(ICV) + 1)* (IC_I_E(ICV) - IC_I_W(ICV) + 1)

               WRITE(ERR_MSG, '(/, 10x, A, 2x, i5,/10x, A, I20)')  &
                    'Constant Number of particles per cell specified         = ', IC_PIC_CONST_NPC(ICV, M), &
                    'Number of estimated parcels  =                          = ', comp_parts(m)

               CALL FLUSH_ERR_MSG(header = .false., FOOTER = .false.)
               comp_parts(m) = 0

               DO K = IC_K_B(ICV), IC_K_T(ICV)
                  DO J = IC_J_S(ICV), IC_J_N(ICV)
                     DO I = IC_I_W(ICV), IC_I_E(ICV)

                        IF(.not.IS_ON_myPE_wobnd(I,J,K)) cycle
                        IF(DEAD_CELL_AT(I,J,K)) cycle
                        IJK = FUNIJK(I,J,K)

                        IF(.not.FLUID_AT(IJK)) cycle

                        CORD_START(1) = XE(I-1)
                        CORD_START(2) = YN(J-1)
                        CORD_START(3) = ZT(K-1)
                        DOML(1) = DX(I)
                        DOML(2) = DY(J)
                        IF(DO_K) DOML(3) = DZ(K)

                        !VOL(IJK) is calculated by set_geometry1, which is
                        !called in mfix.f after calling get_data. At this
                        !stage, VOL(IJK) = zero
                        !VOLIJK = VOL(IJK)
                        VOLIJK       = DX(I)*DY(J)*DZ(K)
                        VOLIJK_UNCUT = DX(I)*DY(J)*DZ(K)

                        REAL_PARTS(M) = 0.
                        COMP_PARTS(M) = 0

                        REAL_PARTS(M) = 6.d0*EP_SM*VOLIJK/ &
                             (PI*(Des_D_P0(M)**3.d0))


                        COMP_PARTS(M) = IC_PIC_CONST_NPC(ICV,M)
                        IF(CUT_CELL_AT(IJK)) COMP_PARTS(M) = &
                             MAX(1,INT(VOLIJK*real(COMP_PARTS(M))/VOLIJK_UNCUT))

                        STAT_WT = REAL_PARTS(M)/REAL(COMP_PARTS(M))


                        ALLOCATE(RANDPOS    (COMP_PARTS(M), DIMN))
                        ALLOCATE(PVEL_TEMP  (COMP_PARTS(M), DIMN))


                        DO IDIM = 1, DIMN
                           CALL UNI_RNO(RANDPOS(1:COMP_PARTS(M), IDIM))
                           PVEL_TEMP(:, IDIM) = VEL_MEAN(IDIM)
                           IF(VEL_SIG(IDIM).GT.ZERO) THEN
                              CALL NOR_RNO(PVEL_TEMP(1:COMP_PARTS(M),IDIM), &
                                   VEL_MEAN(IDIM),VEL_SIG(IDIM))
                           ENDIF
                        ENDDO

                        DO IPCOUNT = 1, COMP_PARTS(M)

                           RNP_PIC(M) = RNP_PIC(M) + STAT_WT
                           CNP_PIC(M) = CNP_PIC(M) + 1
                           PART_MPHASE_BYIC(ICV, M) =  PART_MPHASE_BYIC(ICV, M) &
                                + 1
                           REALPART_MPHASE_BYIC(ICV, M) =  REALPART_MPHASE_BYIC(ICV, M) &
                                + STAT_WT

                           RAD = DES_D_P0(M)*HALF
                           DENS  =  DES_RO_S(M)
                           POSITION(:) = CORD_START(:) + RANDPOS(IPCOUNT, :)*DOML(:)
                           VELOCITY(1:DIMN) = PVEL_TEMP(IPCOUNT,1:DIMN)

                           CALL GEN_AND_ADD_TO_PART_LIST(PART_LIST_BYIC, M, POSITION(1:DIMN), &
                                VELOCITY(1:DIMN), RAD, DENS, STAT_WT)
                           PART_LIST_BYIC%cell(1) = I
                           PART_LIST_BYIC%cell(2) = J
                           PART_LIST_BYIC%cell(3) = K
                           PART_LIST_BYIC%cell(4) = IJK
                        ENDDO
                        deallocate(pvel_temp)
                        deallocate(randpos)
                     ENDDO
                  ENDDO
               ENDDO

            ELSEIF(CONST_STATWT) THEN

               CORD_START(1) = IC_X_W(ICV)
               CORD_START(2) = IC_Y_S(ICV)
               CORD_START(3) = IC_Z_B(ICV)
               CORD_END(1)   = IC_X_E(ICV)
               CORD_END(2)   = IC_Y_N(ICV)
               CORD_END(3)   = IC_Z_T(ICV)

               DOML(:) = CORD_END(:) - CORD_START(:)

               IF(NO_K) DOML(3) = DZ(1)

               VOLIJK = DOML(1)*DOML(2)*DOML(3)

               REAL_PARTS(M) = 0.
               COMP_PARTS(M) = 0

               REAL_PARTS(M) = 6.d0*EP_SM*VOLIJK/ &
                    (PI*(Des_D_P0(M)**3.d0))

               COMP_PARTS(M) = MAX(1, &
                    INT(REAL_PARTS(M)/REAL(IC_PIC_CONST_STATWT(ICV,M))))

               ! although the weight was specified in the input file, the integer
               ! number of CP's requires recalculating statistical weight. This
               ! slightly different statistical weight will ensure that the initial
               ! volume fraction is as inputted. If the input statwt_pic is used,
               ! then the initial volume fraction might be slightly less than the
               ! input initial volume fraction
               STAT_WT = REAL_PARTS(M)/REAL(COMP_PARTS(M))

               WRITE(ERR_MSG, '(/10x, A, 2x, ES15.7, /10x, A, ES15.7, /10x, A, i20)') &
                    'Constant statistical weight specified                  = ', IC_PIC_CONST_STATWT(ICV, M), &
                    'Statistical weight computed (could be a little different from specified value)  = ', STAT_WT, &
                    'Number of estimated parcels                             = ', comp_parts(m)

               CALL FLUSH_ERR_MSG(header = .false., FOOTER = .false.)

               ALLOCATE(RANDPOS    (COMP_PARTS(M), DIMN))
               ALLOCATE(PVEL_TEMP  (COMP_PARTS(M), DIMN))


               DO IDIM = 1, DIMN
                  CALL UNI_RNO(RANDPOS(1:COMP_PARTS(M), IDIM))
                  PVEL_TEMP(:, IDIM) = VEL_MEAN(IDIM)
                  IF(VEL_SIG(IDIM).GT.ZERO) THEN
                     CALL NOR_RNO(PVEL_TEMP(1:COMP_PARTS(M),IDIM), &
                          VEL_MEAN(IDIM),VEL_SIG(IDIM))
                  ENDIF
               ENDDO

               DO IPCOUNT = 1, COMP_PARTS(M)

                  POSITION(:) = CORD_START(:) + RANDPOS(IPCOUNT, :)*DOML(:)

                  I = DES_GETINDEXFROMPOS(ISTART1, IEND1, &
                       position(1), XE(1:size(XE,1)-1),'X','I')
                  !-1 above since xe ranges from 0:IMAX3, so size is imax3 + 1.
                  !therefore, (1:size(xe,1)) will give 1:imax3 + 1, resulting in a seg error.

                  J = DES_GETINDEXFROMPOS(JSTART1,JEND1,position(2), &
                       YN(1:size(YN,1)-1),'Y','J')

                  K = 1
                  IF(DO_K) K = DES_GETINDEXFROMPOS(KSTART1, &
                       KEND1,position(3),ZT(1:size(ZT,1)-1),'Z','K')

                  IF(.not.IS_ON_myPE_wobnd(I,J,K)) cycle

                  IF(DEAD_CELL_AT(I,J,K)) cycle

                  IJK = FUNIJK(I, J, K)
                  IF(.not.FLUID_AT(IJK)) cycle

                  RNP_PIC(M) = RNP_PIC(M) + STAT_WT
                  CNP_PIC(M) = CNP_PIC(M) + 1
                  PART_MPHASE_BYIC(ICV, M) =  PART_MPHASE_BYIC(ICV, M) &
                       + 1

                  REALPART_MPHASE_BYIC(ICV, M) =  REALPART_MPHASE_BYIC(ICV, M) &
                       + STAT_WT

                  RAD = DES_D_P0(M)*HALF
                  DENS  =  DES_RO_S(M)
                  VELOCITY(1:DIMN) = PVEL_TEMP(IPCOUNT,1:DIMN)

                  CALL GEN_AND_ADD_TO_PART_LIST(PART_LIST_BYIC, M, POSITION(1:DIMN), &
                       VELOCITY(1:DIMN), RAD, DENS, STAT_WT)
                  PART_LIST_BYIC%cell(1) = I
                  PART_LIST_BYIC%cell(2) = J
                  PART_LIST_BYIC%cell(3) = K
                  PART_LIST_BYIC%cell(4) = IJK
               ENDDO
               deallocate(pvel_temp)
               deallocate(randpos)

            ENDIF

            CALL global_all_sum(PART_MPHASE_BYIC(ICV, M))
            CALL global_all_sum(REALPART_MPHASE_BYIC(ICV, M))


            WRITE(ERR_MSG, '(10x, A, I15, /, 10x, A, ES15.7)') &
                 '# of computational particles or parcels acutally seeded = ', PART_MPHASE_BYIC(ICV, M), &
                 '# of implied real particles for above parcel count          = ', REALPART_MPHASE_BYIC(ICV, M)

            CALL FLUSH_ERR_MSG(header = .false., FOOTER = .false.)

         ENDDO MLOOP
         ! Now merge this IC specific list to the master list
         if(associated(part_list_byic)) then
            If(.not.associated(orig_part_list)) then
               !just point the orig_part_list to part_list_byic
               orig_part_list => part_list_byic
            else
               nullify(part)
               part => orig_part_list
               DO WHILE (ASSOCIATED(part%next))
                  part => part%next
               ENDDO

               part%next => part_list_byic
               part_list_byic%prev => part
            endif

            nullify(part_list_byic)
         endif
      ENDDO                     !ICV = 1, DIMENSION_IC


      CALL global_all_sum(RNP_PIC)
      CALL global_all_sum(CNP_PIC)

      PART_MPHASE(1:DES_MMAX)     = CNP_PIC(1:DES_MMAX)
      REALPART_MPHASE(1:DES_MMAX) = RNP_PIC(1:DES_MMAX)
      !particles will be set in gener_part_config routine

! Local reporting for the user
      WRITE(ERR_MSG, 2015)  COUNT_IC, DES_MMAX
      CALL FLUSH_ERR_MSG(HEADER = .false., FOOTER = .false.)

2015  FORMAT(/, 'Done generating the initial configuration for PIC model', /, &
      'Total IC region count with non zero solids = ', I5,2X, /, &
      'Total discrete solid phases                = ', I5,2X, /, &
      'Global reporting per solid phase follows')

      DO M = 1, DES_MMAX
         WRITE(err_msg, 2023) M, PART_MPHASE(M), REALPART_MPHASE(M)
         CALL FLUSH_ERR_MSG(HEADER = .false., FOOTER = .false.)
      enddo


2023  FORMAT(/5x, &
           'For solid phase M:                          ',  I5, /5x, &
           'Total number of parcels                    = ', I15, /5x, &
           'Total number of implied physical particles = ', ES15.7)

      WRITE(ERR_MSG,*) ''
      CALL FLUSH_ERR_MSG(HEADER = .false.)

      IF(CARTESIAN_GRID) then
         write(ERR_MSG, 3000)
         CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
 3000 FORMAT('Cartesian grid detected with particle configuration,'/&
         'deleting particles located outsidte the domain.')

         write(err_msg, '(A)') 'Now Marking particles to be deleted'
         call flush_err_msg(header = .false., footer = .false.)

      part => orig_part_list
      DO WHILE (ASSOCIATED(part))
         DELETE_PART = .false.
         IJK = part%cell(4)

! Only look for particles that are in domain. Particles in dead cells
! have already been flagged as outside of the domain in
!  BIN_PARTICLES_TO_CELL.
         IF(part%indomain) THEN
               IF(FLUID_AT(IJK)) then
                  CALL CHECK_IF_PARCEL_OVERLAPS_STL(part%position(1:dimn), &
                       DELETE_PART)
               ELSE
                  DELETE_PART = .true.
               ENDIF
         ENDIF

         IF(DELETE_PART) PART%INDOMAIN = .FALSE.
         PART => PART%NEXT
       ENDDO

         write(ERR_MSG,"('Finished marking particles to delete.')")
         CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)

      ENDIF

      LORIG_ALL = 0
      LDEL_ALL = 0
      LREM_ALL = 0

      part => orig_part_list
      DO WHILE (ASSOCIATED(part))
         LORIG_ALL(mype) = LORIG_ALL(mype) + 1

         IF(part%INDOMAIN) then
            LREM_ALL(mype) = LREM_ALL(mype) + 1
         else
            LDEL_ALL(mype) = LDEL_ALL(mype) + 1
         ENDIF

         part => part%next
      ENDDO

      CALL global_all_sum(LORIG_ALL)
      CALL global_all_sum(LREM_ALL)
      CALL global_all_sum(LDEL_ALL)

      PIP = LREM_ALL(mype) !Setting pip on each processor equal to remaining paticles

      CALL PARTICLE_GROW(PIP)

      part => orig_part_list
      count_part = 0
      DO WHILE (ASSOCIATED(part))
         IF(part%INDOMAIN) then
            count_part = count_part + 1

            des_pos_new(1:dimn, count_part) = part%position(1:dimn)
            des_vel_new(1:dimn, count_part) = part%velocity(1:dimn)

            DES_RADIUS(count_part) = part%rad
            RO_SOL(count_part)  = part%dens

            PIJK(count_part,5) = part%M
            phase  = part%M
            IF (DO_OLD) THEN
               des_vel_old(:, count_part) = des_vel_new(:, count_part)
               des_pos_old(:, count_part) = des_pos_new(:, count_part)
            ENDIF

            CALL SET_NORMAL(count_part)

            IF(SOLIDS_MODEL(phase).eq.'PIC') then
               DES_STAT_WT(count_part) = part%STATWT
            ELSEIF(SOLIDS_MODEL(phase).eq.'DEM') then
               IF (DO_OLD) OMEGA_OLD(:, count_part) = zero
               OMEGA_NEW(:, count_part) = zero
            ENDIF

         ENDIF
         part => part%next

      ENDDO
      if(count_part.ne.pip) then
         write(err_msg, 1000) count_part, pip
         call flush_err_msg(abort = .true.)
      endif

 1000 format('Particles copied from linked lists:', 2x, i15, /, &
      'not equal to particles earlier calculated to be inside domain', 2x, i15, / &
      'This should not have happened, exitting')
      CALL FINL_ERR_MSG

      CALL DEALLOC_PART_LIST(orig_part_list)

      END SUBROUTINE GENERATE_PARTICLE_CONFIG_MPPIC

      SUBROUTINE DEALLOC_PART_LIST(LIST_NAME)

      IMPLICIT NONE

      TYPE(PARTICLE), POINTER  :: LIST_NAME, part => NULL(), part_old => NULL()

      IF(.not.associated(list_name)) return !do nothing

      part =>  list_name

      DO WHILE (ASSOCIATED(part))
         part_old => part
         part => part%next
         DEALLOCATE(part_old)
      enddo
      RETURN
      END SUBROUTINE DEALLOC_PART_LIST
      SUBROUTINE GEN_AND_ADD_TO_PART_LIST(LIST_NAME, PHASE, POS, VEL,RAD, DENS, STATWT)
      USE discretelement, only : dimn

      IMPLICIT NONE

      TYPE(PARTICLE), POINTER  :: LIST_NAME
      INTEGER, INTENT(IN) ::  PHASE
      double precision, INTENT(IN), DIMENSION(DIMN) :: POS, VEL
      double precision, INTENT(IN) :: RAD, DENS, STATWT

      TYPE(PARTICLE), POINTER  :: NEW_PART => NULL()

      ALLOCATE(NEW_PART)

      NEW_PART%M = PHASE
      NEW_PART%INDOMAIN = .true.

      NEW_PART%POSITION(1:DIMN) = POS(1:DIMN)
      NEW_PART%VELOCITY(1:DIMN) = VEL(1:DIMN)
      NEW_PART%STATWT = STATWT

      NEW_PART%RAD = RAD
      NEW_PART%DENS = DENS
      NEW_PART%STATWT = STATWT

      IF(.NOT.ASSOCIATED(LIST_NAME)) THEN !FIRST ENTRY
         LIST_NAME => NEW_PART  !make the new entry as head
      ELSE
         NEW_PART%NEXT => LIST_NAME
         LIST_NAME%PREV => NEW_PART
         LIST_NAME =>  NEW_PART
      END IF
      END SUBROUTINE GEN_AND_ADD_TO_PART_LIST

    END MODULE GENERATE_PARTICLES
