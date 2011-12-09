 !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CHECK_DES_DATA                                         C
!  Purpose: DES - Check user input Data   
!                                                                      C
!                                                                      C
!  Author: S. Benyahia                                Date: 10-Nov-08  C
!  Reviewer:                                                           C
!  Comments: All user's input data are checked here                    C
!  Revision : Some of the checks made in des_allocate_arrays are moved 
!             here. In addition write statments are now made to log file
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CHECK_DES_DATA
      
      USE param1
      USE geometry
      USE funits
      USE compar      
      USE discretelement
      USE run
      USE constant
      USE physprop
      use desgrid 
      USE indices
      USE fldvar
      USE toleranc 
      USE mpi_utility
      USE output 
      USE mfix_pic
      USE cutcell
      IMPLICIT NONE
!-----------------------------------------------
! Local Variables
!-----------------------------------------------      
      INTEGER CHECK_MPI, M, IJK
      LOGICAL FLAG_WARN
      INTEGER I, J, K
! quantities used to check/specify mesh size if des_neighbor_search = 4      
      DOUBLE PRECISION DL_TMP, TMP_FACTOR     
! domain volume      
      DOUBLE PRECISION :: VOL_DOMAIN
! the maximum and minimum specified particle diameter 
      DOUBLE PRECISION MAX_DIAM, MIN_DIAM
! for gener_part_config, the total solids volume fraction
      DOUBLE PRECISION TOT_VOL_FRAC
! number of real and comp. particles in a cell 
! for MP-PIC case 

      DOUBLE PRECISION REAL_PARTS(DIM_M), COMP_PARTS(DIM_M)
      
! volume of the cell 
      DOUBLE PRECISION VOLIJK, VOLIJK_UNCUT
      
!----------------------------------------------- 
      INCLUDE 'function.inc'

      INCLUDE 'ep_s1.inc'
      INCLUDE 'ep_s2.inc'

      if(dmp_log.and.debug_des) WRITE(unit_log,'(1X,A)')&
         '---------- START CHECK_DES_DATA ---------->'

      IF (.NOT.DES_CONTINUUM_HYBRID) THEN
! MMAX, D_p0 and RO_s are to be strictly associated with the continuum
! model and are no longer to be used in the DEM. However, when not using
! DES_CONTINUUM_HYBRID the current DEM code still assumes that the user
! specifies MMAX, D_P0 and RO_S.  These are then linked to their
! respective DES variables.  Note that valid values of MMAX, D_P0 and 
! RO_s are ensured by check_data_04.
         DES_MMAX = MMAX
         DES_D_p0(:) = D_p0(:)
         DES_RO_s(:) = RO_s(:)
      ENDIF

      IF(COORDINATES == 'CYLINDRICAL') THEN
         if(dmp_log) WRITE (UNIT_LOG, 1000)
         CALL MFIX_EXIT(myPE)
      ENDIF
      
      PRINT_DES_SCREEN = FULL_LOG.and.myPE.eq.pe_IO
      
      PRINT_DES_SCREEN = .true.
      IF(DMP_LOG) WRITE(UNIT_LOG,*) 'FULL_LOG  = ', PRINT_DES_SCREEN, ' FOR DES ON IO PROC'
      IF(mype.eq.pe_IO) WRITE(*,*) 'FULL_LOG  = ', PRINT_DES_SCREEN, ' FOR DES ON IO PROC'
      !read(*,*)
! remove the number of processor check 
!      CHECK_MPI = NODESI * NODESJ * NODESK
!      IF(CHECK_MPI.NE.1) THEN
!         WRITE (UNIT_LOG, 1001)
!         CALL MFIX_EXIT(myPE)
!      ENDIF

! Check dimension
      IF(NO_I.OR.NO_J) THEN ! redundant, this check is made in check_data_03
         IF(DMP_LOG) WRITE(UNIT_LOG,1039)
         CALL MFIX_EXIT(myPE)
      ENDIF
      IF(DIMN == UNDEFINED_I) THEN
         IF(NO_K) THEN
            DIMN = 2 
            IF(DMP_LOG)  WRITE(UNIT_LOG, 1040)
            WRITE(*,1040) 
         ELSE              
            IF(DMP_LOG) WRITE(UNIT_LOG, 1041)
            CALL MFIX_EXIT(myPE)
         ENDIF
      ENDIF
      IF(DIMN > 3) THEN
         IF(DMP_LOG) WRITE(UNIT_LOG, 1042)
         CALL MFIX_EXIT(myPE)
      ENDIF

      IF(IMAX.GT.1.AND.JMAX.GT.1.AND.KMAX.GT.1) THEN 
         IF(DIMN.NE.3) THEN 
            IF(DMP_LOG)  WRITE(*,2001) DIMN 
            CALL MFIX_EXIT(myPE)
         ENDIF
      ENDIF

!!!! code commented below to let users use cohesion
! if you want to run DEM with the existing cohesion implementation then
! comment out the following 4 lines; also note that this feature is
! unlikely to function correctly with the current version of DEM so
! proceed at your own risk
!      IF(USE_COHESION) THEN
!         if(dmp_log) WRITE(UNIT_LOG, 1035)
!         CALL MFIX_EXIT(myPE)
!      ENDIF

      IF ( (DES_PERIODIC_WALLS_X .OR. DES_PERIODIC_WALLS_Z .OR. &
            DES_PERIODIC_WALLS_Z) .AND. .NOT.DES_PERIODIC_WALLS) THEN
         DES_PERIODIC_WALLS = .TRUE.
         if(dmp_log) WRITE(UNIT_LOG, 1017)
      ENDIF


      IF(DES_NEIGHBOR_SEARCH.EQ.1) THEN
         if(dmp_log) WRITE(*, 1045) 1, 'N-SQUARE'
      ELSEIF(DES_NEIGHBOR_SEARCH.EQ.2) THEN
         if(dmp_log) WRITE(*, 1045) 2, 'QUADTREE'
      ELSEIF(DES_NEIGHBOR_SEARCH.EQ.3) THEN
         if(dmp_log) WRITE(*, 1045) 3, 'OCTREE'
      ELSEIF(DES_NEIGHBOR_SEARCH.EQ.4) THEN
         if(dmp_log) WRITE(*, 1045) 4, 'GRID BASED'
      ELSE
         if(dmp_log) WRITE(UNIT_LOG, 1003)
         CALL MFIX_EXIT(myPE)         
      ENDIF


! Check the output file format 
      IF (DES_OUTPUT_TYPE /= UNDEFINED_C) THEN
         IF(TRIM(DES_OUTPUT_TYPE) /= 'TECPLOT') THEN
            if(dmp_log) WRITE (*, 1038)
            CALL MFIX_EXIT(myPE)   
         ENDIF
      ENDIF

      IF(DES_PERIODIC_WALLS) THEN
         IF(.NOT.DES_PERIODIC_WALLS_X .AND. .NOT.DES_PERIODIC_WALLS_Y .AND. &
            .NOT.DES_PERIODIC_WALLS_Z) THEN
            if(dmp_log) WRITE (UNIT_LOG, 1002)
            CALL MFIX_EXIT(myPE)
         ENDIF
         IF (DES_NEIGHBOR_SEARCH .EQ. 1) THEN
            if(dmp_log) WRITE (UNIT_LOG, 1018)
         ENDIF
         IF (DES_NEIGHBOR_SEARCH .EQ. 2 .OR. &
          DES_NEIGHBOR_SEARCH .EQ. 3) THEN
            if(dmp_log) WRITE (UNIT_LOG, 1019)
            CALL MFIX_EXIT(myPE)
         ENDIF
      ENDIF

      IF(DES_LE_BC) THEN
         IF (DES_CONTINUUM_COUPLED) THEN
            WRITE(UNIT_LOG, 1064)
             CALL MFIX_EXIT(myPE)
         ENDIF 

         IF (DES_NEIGHBOR_SEARCH .NE. 4) THEN
            WRITE(UNIT_LOG, 1060)
            CALL MFIX_EXIT(myPE)
         ENDIF

! not all possible shear directions are fully coded         
         IF (DIMN .EQ. 2) THEN
            IF(TRIM(DES_LE_SHEAR_DIR) .NE. 'DUDY' .AND. &
               TRIM(DES_LE_SHEAR_DIR) .NE. 'DVDX') THEN
               WRITE(UNIT_LOG, 1061)
               CALL MFIX_EXIT(myPE)
            ENDIF
         ELSEIF(DIMN.EQ.3) THEN
            IF(TRIM(DES_LE_SHEAR_DIR) .NE. 'DUDY') THEN ! .AND. & 
!               TRIM(DES_LE_SHEAR_DIR) .NE. 'DUDZ' .AND. &
!               TRIM(DES_LE_SHEAR_DIR) .NE. 'DVDX' .AND. &
!               TRIM(DES_LE_SHEAR_DIR) .NE. 'DVDZ' .AND. &
!               TRIM(DES_LE_SHEAR_DIR) .NE. 'DWDX' .AND. &
!               TRIM(DES_LE_SHEAR_DIR) .NE. 'DWDY') THEN
               WRITE(UNIT_LOG, 1062)
               CALL MFIX_EXIT(myPE)
            ENDIF
         ENDIF

         IF (DES_PERIODIC_WALLS) THEN
            DES_PERIODIC_WALLS = .FALSE.
            DES_PERIODIC_WALLS_X = .FALSE.
            DES_PERIODIC_WALLS_Y = .FALSE.
            DES_PERIODIC_WALLS_Z = .FALSE.            
            WRITE(UNIT_LOG, 1063)
            WRITE(*,1063)
         ENDIF

      ENDIF

      IF (TRIM(DES_INTG_METHOD) /= 'ADAMS_BASHFORTH' .AND. &
          TRIM(DES_INTG_METHOD) /= 'EULER') THEN
! stop if the specified integration method is unavailable
         if(dmp_log) WRITE (UNIT_LOG, 1034)
         CALL MFIX_EXIT(myPE)
      ENDIF

      IF(.NOT.MPPIC) THEN 
      IF (DES_COLL_MODEL /= UNDEFINED_C) THEN      
         IF (TRIM(DES_COLL_MODEL) .NE. 'HERTZIAN') THEN
            if(dmp_log) WRITE(UNIT_LOG, 1006)
            CALL MFIX_EXIT(myPE)
         ENDIF
      ENDIF

      IF (TRIM(DES_COLL_MODEL) .EQ. 'HERTZIAN') THEN
! check young's modulus and poisson ratio
         IF(EW_YOUNG == UNDEFINED ) THEN
            if(dmp_log) WRITE (UNIT_LOG, 1020)
            CALL MFIX_EXIT(myPE)
         ENDIF
         IF(VW_POISSON == UNDEFINED) THEN
            if(dmp_log) WRITE (UNIT_LOG, 1021)
            CALL MFIX_EXIT(myPE)
         ELSE
           IF (VW_POISSON > 0.5d0 .OR. VW_POISSON <= -ONE) THEN
              if(dmp_log) WRITE (UNIT_LOG, 1033)
           ENDIF
         ENDIF
         DO M = 1, MMAX
            IF(E_YOUNG(M) == UNDEFINED) THEN
               if(dmp_log) WRITE (UNIT_LOG, 1022)
               CALL MFIX_EXIT(myPE)
            ENDIF
            IF(V_POISSON(M) == UNDEFINED) THEN
               if(dmp_log) WRITE (UNIT_LOG, 1023)
               CALL MFIX_EXIT(myPE)
            ELSE
              IF (V_POISSON(M) > 0.5d0 .OR. V_POISSON(M) <= -ONE) THEN
                 if(dmp_log) WRITE (UNIT_LOG, 1033)
              ENDIF
            ENDIF
         ENDDO
! check particle-particle tangential restitution coefficient      
         DO M = 1, MMAX+MMAX*(MMAX-1)/2
            IF(DES_ET_INPUT(M) == UNDEFINED) THEN
               if(dmp_log) WRITE (UNIT_LOG, 1024) MMAX+MMAX*(MMAX-1)/2
               CALL MFIX_EXIT(myPE)
            ENDIF
         ENDDO
         DO M = 1, MMAX+MMAX*(MMAX-1)/2
            IF(DES_ET_INPUT(M) > ONE .OR. DES_ET_INPUT(M) < ZERO) THEN
               if(dmp_log) WRITE (UNIT_LOG, 1025)
               CALL MFIX_EXIT(myPE)
            ENDIF
         ENDDO
! check particle-wall tangential restitution coefficient
         DO M = 1, MMAX
            IF(DES_ET_WALL_INPUT(M) == UNDEFINED) THEN
               if(dmp_log) WRITE (UNIT_LOG, 1026)
               CALL MFIX_EXIT(myPE)
            ENDIF
         ENDDO
         DO M = 1, MMAX
            IF(DES_ET_WALL_INPUT(M) > ONE .OR. DES_ET_WALL_INPUT(M) < ZERO) THEN
               if(dmp_log) WRITE (UNIT_LOG, 1027)
               CALL MFIX_EXIT(myPE)
            ENDIF
         ENDDO
! if following are assigned warn user they are discarded
         IF(KN .NE. UNDEFINED .OR. KN_W .NE. UNDEFINED) then 
            if(dmp_log) WRITE (UNIT_LOG, 1028)
         end if 
         IF(KT_FAC .NE. UNDEFINED .OR. KT_W_FAC .NE. UNDEFINED) then 
            if(dmp_log) WRITE (UNIT_LOG, 1029)
         end if 
         IF(DES_ETAT_FAC .NE. UNDEFINED .OR. &
            DES_ETAT_W_FAC .NE. UNDEFINED) then 
            if(dmp_log) WRITE (UNIT_LOG, 1030)
         end if 

      ELSE  ! default linear spring-dashpot model

! check normal spring constants 
         IF(KN == UNDEFINED .OR. KN_W == UNDEFINED) THEN
            if(dmp_log) WRITE (UNIT_LOG, 1004)
            CALL MFIX_EXIT(myPE)
         ENDIF
! check for tangential spring constant factors         
         IF(KT_FAC == UNDEFINED .OR. KT_W_FAC == UNDEFINED) THEN
            if(dmp_log) WRITE (UNIT_LOG, 1005)
         ENDIF
         IF(KT_FAC .NE. UNDEFINED) THEN
            IF(KT_FAC > ONE .OR. KT_FAC < ZERO) THEN
               if(dmp_log) WRITE (UNIT_LOG, 1016)
               CALL MFIX_EXIT(myPE)
            ENDIF
         ENDIF
         IF(KT_W_FAC .NE. UNDEFINED) THEN
            IF(KT_W_FAC > ONE .OR. KT_W_FAC < ZERO) THEN
               if(dmp_log) WRITE (UNIT_LOG, 1016)
               CALL MFIX_EXIT(myPE)
            ENDIF
         ENDIF
! check for tangential damping factor
         IF(DES_ETAT_FAC == UNDEFINED) THEN
            if(dmp_log) WRITE (UNIT_LOG, 1010)
         ELSEIF(DES_ETAT_FAC > ONE .OR. DES_ETAT_FAC < ZERO) THEN
            if(dmp_log) WRITE (UNIT_LOG, 1015)
            CALL MFIX_EXIT(myPE)
         ENDIF
         IF(DES_ETAT_W_FAC == UNDEFINED) THEN
            if(dmp_log) WRITE (UNIT_LOG, 1011)
         ELSEIF(DES_ETAT_W_FAC > ONE .OR. DES_ETAT_W_FAC < ZERO) THEN
            if(dmp_log) WRITE (UNIT_LOG, 1015)
            CALL MFIX_EXIT(myPE)
         ENDIF
! if following are assigned warn user they are discarded
         FLAG_WARN = .FALSE.
         DO M = 1, MMAX+MMAX*(MMAX-1)/2
            IF(DES_ET_INPUT(M) .NE. UNDEFINED) FLAG_WARN = .TRUE.
         ENDDO
         IF (FLAG_WARN) WRITE(UNIT_LOG,1031)
         FLAG_WARN = .FALSE.
         DO M = 1, MMAX
            IF(DES_ET_WALL_INPUT(M) .NE. UNDEFINED) FLAG_WARN = .TRUE.
         ENDDO
         IF (FLAG_WARN) WRITE(UNIT_LOG,1032)
         FLAG_WARN = .FALSE.

      ENDIF   ! endif des_coll_model .eq. hertzian


! check particle-particle normal restitution coefficient      
      DO M = 1, MMAX+MMAX*(MMAX-1)/2
         IF(DES_EN_INPUT(M) == UNDEFINED) THEN
            if(dmp_log) WRITE (UNIT_LOG, 1008) MMAX+MMAX*(MMAX-1)/2
            CALL MFIX_EXIT(myPE)
         END IF
      ENDDO
      DO M = 1, MMAX+MMAX*(MMAX-1)/2
         IF(DES_EN_INPUT(M) > ONE .OR. DES_EN_INPUT(M) < ZERO) THEN
            if(dmp_log) WRITE (UNIT_LOG, 1012)
            CALL MFIX_EXIT(myPE)
         ENDIF
      ENDDO

! check particle-wall normal restitution coefficient
      DO M = 1, MMAX
         IF(DES_EN_WALL_INPUT(M) == UNDEFINED) THEN
            if(dmp_log) WRITE (UNIT_LOG, 1009)
            CALL MFIX_EXIT(myPE)
         ENDIF
      ENDDO
      DO M = 1, MMAX
         IF(DES_EN_WALL_INPUT(M) > ONE .OR. DES_EN_WALL_INPUT(M) < ZERO) THEN
            if(dmp_log) WRITE (UNIT_LOG, 1013)
            CALL MFIX_EXIT(myPE)
         ENDIF
      ENDDO

! check coefficient friction 
      IF(MEW > ONE .OR. MEW_W > ONE .OR. MEW < ZERO .OR. MEW_W < ZERO) THEN
         if(dmp_log) WRITE (UNIT_LOG, 1014)
         CALL MFIX_EXIT(myPE)
      END IF

      ENDIF

! Overwrite user's input in case of DEM (no fluid)
      IF(.NOT.DES_CONTINUUM_COUPLED) DES_INTERP_ON = .FALSE.
      

! TSUJI drag option only available with syam_obrien drag model
      IF(DES_CONTINUUM_COUPLED .AND. TSUJI_DRAG .AND. &
         TRIM(DRAG_TYPE) /= 'SYAM_OBRIEN') THEN
         IF(DMP_LOG) WRITE(UNIT_LOG, 1056)
      ENDIF


      IF(MPPIC) THEN 
         !DES_INTERP_ON = .TRUE.
         !IF(DMP_LOG) WRITE(UNIT_LOG,*) 'INTERPOLATION IS ALWAYS ON FOR MPPIC'
         !IF(PRINT_DES_SCREEN) WRITE(*,*) 'INTERPOLATION IS ALWAYS ON FOR MPPIC'

         if(MPPIC_SOLID_STRESS_SNIDER) WRITE(UNIT_LOG,*) 'USING THE SNIDER MODEL FOR SOLID STRESS AND THE INTEGRATION APPROACH'
         IF(MPPIC_SOLID_STRESS_SNIDER.and.PRINT_DES_SCREEN)  WRITE(*,*) 'USING THE SNIDER MODEL FOR SOLID STRESS AND THE INTEGRATION APPROACH'

         
! check particle-wall normal restitution coefficient
         DO M = 1, MMAX
            IF(DES_EN_WALL_INPUT(M) == UNDEFINED) THEN
               if(dmp_log) WRITE (UNIT_LOG, 1009)
               CALL MFIX_EXIT(myPE)
            ENDIF
         ENDDO
         DO M = 1, MMAX
            IF(DES_EN_WALL_INPUT(M) > ONE .OR. DES_EN_WALL_INPUT(M) < ZERO) THEN
               if(dmp_log) WRITE (UNIT_LOG, 1013)
               CALL MFIX_EXIT(myPE)
            ENDIF
         ENDDO
         
      ENDIF
!***********************************************
! Pradeep moved following code from des_allocate_arrays 
! Valid D_p0(M) are needed here if gener_part_config.  A second check 
! for realistic d_p0 values is made in check_data_04 but this routine 
! is called after des_allocate_arrays (i.e. would be too late).
! Valid D_P0(M) are also needed to identify which solids phase each
! particle belongs in (although the sorting is performed after 
! check_data_04 is called), and to determine the maximum particle size
! in the system (MAX_RADIUS), which in turn is used for various tasks
       MAX_DIAM = ZERO
       MIN_DIAM = LARGE_NUMBER
       DO M = 1,MMAX   
          IF (D_P0(M)<ZERO .OR. D_P0(M)==UNDEFINED) THEN    
             if(dmp_log)write(unit_log,'(3X,A,A)') &   
                'D_P0 must be defined and >0 in mfix.dat ',&   
                'for M = 1,MMAX'   
                CALL MFIX_EXIT(myPE)   
          ENDIF   
          MAX_DIAM = MAX(MAX_DIAM, D_P0(M))
          MIN_DIAM = MIN(MIN_DIAM, D_P0(M))
       ENDDO   
       DO M = MMAX+1, DIMENSION_M   
          IF (D_P0(M) /= UNDEFINED) THEN   
             if(dmp_log)write(unit_log,'(3X,A,A)') &   
                'Too many D_P0 are defined for given MMAX'   
             CALL MFIX_EXIT(myPE)   
          ENDIF   
       ENDDO
       MAX_RADIUS = 0.5d0*MAX_DIAM
       MIN_RADIUS = 0.5d0*MIN_DIAM

! Check that the depth of the simulation in 2D exceeds the largest 
! particle size to ensure correct calculation of volume fraction.  This
! is important for coupled simulations (not essential for pure granular
! simulations)
      IF (DES_CONTINUUM_COUPLED .AND. DIMN == 2) THEN
         IF (2.0d0*MAX_RADIUS > ZLENGTH) THEN
            WRITE(UNIT_LOG, 1036)
            CALL MFIX_EXIT(myPE)
         ENDIF
      ENDIF


!Rahul: 
!If run_type is not 'NEW' then force the gener_part_config to .false. 
       IF(RUN_TYPE .ne. 'NEW'.and.GENER_PART_CONFIG ) THEN 
          GENER_PART_CONFIG = .false. 
          IF(DMP_LOG) WRITE(UNIT_LOG, '(3(/2X, A))') &
          & 'IN CHECK_DES_DATA:', &
          & 'SETTING GENER_PART_CONFIG to FALSE', & 
          & 'BECAUSE A RESTART DETECTED'
          IF(mype.eq.pe_IO)  WRITE(*, '(3(/2X, A))') &
          & 'IN CHECK_DES_DATA:', &
          & 'SETTING GENER_PART_CONFIG to FALSE', & 
          & 'BECAUSE A RESTART DETECTED'
       ENDIF
       
       IF(MPPIC) THEN 
          !allocate the cnp_array here only as it will be used
          !if gener_part_config is true. 
          ALLOCATE(CNP_ARRAY(DIMENSION_3, 0:DIMENSION_M))
          
          
          CNP_ARRAY(:, :) = 0

         !rahul:
         !cnp_array(ijk, 0) will contain the cumulative
         !number of real particles later in the handling
         !of inlfow BC for MPPIC. See the mppic_mi_bc in 
         !mppic_wallbc_mod.f 
       ENDIF
!If gener_part_config ensure various quantities
!are defined and valid
      IF(GENER_PART_CONFIG.AND..NOT.MPPIC) THEN 
         TOT_VOL_FRAC = ZERO

! Check if des_eps_xstart, ystart, zstart are set
         IF(DES_EPS_XSTART.EQ.UNDEFINED.OR.&
            DES_EPS_YSTART.EQ.UNDEFINED) THEN
             IF(DMP_LOG) WRITE(UNIT_LOG, 1047)
            CALL MFIX_EXIT(myPE)
         ENDIF
         IF (DIMN.EQ.3 .AND. DES_EPS_ZSTART.EQ.UNDEFINED) THEN
            IF(DMP_LOG) WRITE(UNIT_LOG, 1048)
            CALL MFIX_EXIT(myPE)
         ENDIF


! Check for quantities needed for gener_part_config option and make sure
! they are physically realistic
! Check for quantities needed for gener_part_config option and make sure
! they are physically realistic
         DO M = 1, MMAX
            IF(VOL_FRAC(M) == UNDEFINED) THEN
               IF(DMP_LOG) WRITE(UNIT_LOG, 1049) 
               CALL MFIX_EXIT(myPE)
            ENDIF
            IF (DES_CONTINUUM_COUPLED) THEN
               IF(VOL_FRAC(M) < ZERO .OR. VOL_FRAC(M) > (ONE-EP_STAR)) THEN
                  IF(DMP_LOG) WRITE(UNIT_LOG, 1050) 
                  CALL MFIX_EXIT(myPE)
               ENDIF
            ELSE
               IF(VOL_FRAC(M) < ZERO .OR. VOL_FRAC(M) > ONE) THEN
                  IF(DMP_LOG)  WRITE(UNIT_LOG, 1051)
                  CALL MFIX_EXIT(myPE)
               ENDIF
            ENDIF 
            TOT_VOL_FRAC = TOT_VOL_FRAC + VOL_FRAC(M)
         ENDDO
         
         IF(DES_CONTINUUM_COUPLED) THEN
            IF (TOT_VOL_FRAC > (ONE-EP_STAR)) THEN
               IF(DMP_LOG) WRITE(UNIT_LOG, 1052) 
               CALL MFIX_EXIT(myPE)
            ENDIF
         ELSE
            IF (TOT_VOL_FRAC > ONE) THEN
               IF(DMP_LOG) WRITE(UNIT_LOG, 1053)
               CALL MFIX_EXIT(myPE)
            ENDIF
         ENDIF


         IF(DIMN.EQ.2) THEN 
! Values of DZ(1) or zlength are not guaranteed at this point; however,
! some value is needed to calculate the number of particles
            IF (DZ(1) == UNDEFINED .AND. ZLENGTH == UNDEFINED) THEN
               IF (MAX_DIAM .EQ. ZERO) MAX_DIAM = ONE ! for calculations
               if(dmp_log)write(unit_log,'(5X,A,A,/11X,A,A,ES15.7,/11X,A,A,/11X,A,/)') &
                  'NOTE: neither zlength or dz(1) were specified ',&
                  'so ZLENGTH is being set','to the maximum particle ',&
                  'diameter of ', MAX_DIAM, 'to provide a basis ',&
                  'for calculating the number of particles',&
                  'in the system'
! set zlength to ensure consistency with calculations later on especially 
! when conducting coupled simulations when zlength/dz(1) are not set in 
! the mfix.dat file
               ZLENGTH = MAX_DIAM 
               VOL_DOMAIN  = DES_EPS_XSTART*DES_EPS_YSTART*MAX_DIAM
            ELSEIF (DZ(1) == UNDEFINED) THEN
               if(dmp_log)write(unit_log,'(5X,A,G15.8,A,/11X,A,A,/11X,A,/)') &
                  'NOTE: specified zlength of ', ZLENGTH,&
                  ' is used','to provide a basis for calculating ',&
                  'the number of','particles in the system'
               VOL_DOMAIN  = DES_EPS_XSTART*DES_EPS_YSTART*ZLENGTH
            ELSE
               if(dmp_log)write(unit_log,'(5X,A,G15.8,A,/11X,A,A,/11X,A,/)') &
                  'NOTE: specified dz(1) of ', DZ(1),&
                  ' is used','to provide a basis for calculating ',&
                  'the number of','particles in the system'
               VOL_DOMAIN  = DES_EPS_XSTART*DES_EPS_YSTART*DZ(1)
            ENDIF
         ELSE 
            VOL_DOMAIN  = DES_EPS_XSTART*DES_EPS_YSTART*DES_EPS_ZSTART
         ENDIF

         DO M = 1, MMAX
            PART_MPHASE(M) = FLOOR((6.D0*VOL_FRAC(M)*VOL_DOMAIN)/&
               (PI*(D_P0(M)**3)))
         ENDDO
         PARTICLES = 0
         PARTICLES = SUM(PART_MPHASE(1:MMAX))


         if(dmp_log) WRITE(*,'(/3X,A,A,/)') 'Particle configuration will ', &
         'automatically be generated'
         
         if(dmp_log)write(unit_log,'(5X,A,I5,2X,A,G15.7)') 'MMAX = ', MMAX, &
           ' VOL_DOMAIN = ', VOL_DOMAIN
         if(dmp_log)write(unit_log,'(5X,A,/7X,(ES15.7,2X,$))') 'D_P0(M) = ', &
            D_P0(1:MMAX)
         if(dmp_log)write(unit_log,*)
         if(dmp_log)write(unit_log,'(5X,A,/7X,(G15.8,2X,$))') &
            'VOL_FRAC(M) (solids volume fraction of phase M) = ', &
            VOL_FRAC(1:MMAX)
         if(dmp_log)write(unit_log,*)
         if(dmp_log)write(unit_log,'(5X,A,/7X,(I10,2X,$))') &
            'PART_MPHASE(M) (number particles in phase M) = ', &
            PART_MPHASE(1:MMAX)
         if(dmp_log)write(unit_log,*)
      ENDIF 

      IF(GENER_PART_CONFIG.AND.MPPIC) THEN 
         ALLOCATE(RNP_PIC(MMAX))
         ALLOCATE(CNP_PIC(MMAX))
         RNP_PIC = ZERO
         CNP_PIC = ZERO 
         IF(DIMN.EQ.2) THEN 
!Earlier DZ(1) or ZLENGTH were tested if they were inputted or not.
!this is confusign and especially more so for the cutcells.
!now require that DZ(1) to be specified for 2-dimensional case. 
            IF (DZ(1) == ONE) then 
               WRITE(*,'(5X,A,/5x,A)') &
               'For DIMN = 2, Specify a value for DZ(1) in mfix.dat which is not equal to one', 'If you really want it to be one, then set it very close to one but not exactly one!'
               if(mype.eq.pe_IO) WRITE(UNIT_LOG,'(5X,A,/5x,A)') &
               'For DIMN = 2, Specify a value for DZ(1) in mfix.dat which is not equal to one', 'If you really want it to be one, then set it very close to one but not exactly one!'
               CALL mfix_exit(mype)
            ENDIF

            IF (DZ(1) .ne. ZLENGTH) then 
               !this condition will probably never occur. Redundancy doesn't hurt, however!
               WRITE(*,'(5x, 2(A20,2x, g17.8))') 'DZ(1) = ', DZ(1), 'ZLENGTH = ', ZLENGTH
               WRITE(*,'(5X,A,/5x,A)') &
               'For DIMN = 2, DZ(1) and ZLENGTH are used interchangeably', ' Specify same values for DZ(1) and Zlength'  
               if(myPE.eq.pe_IO) then 
                  WRITE(UNIT_LOG,'(5X,A,/5x,A)') 'For DIMN = 2, DZ(1) and ZLENGTH are used interchangeably', ' Specify same values for DZ(1) and Zlength'  
                  WRITE(UNIT_LOG,'(5x, 2(A20,2x, g17.8))') 'DZ(1) = ', DZ(1), 'ZLENGTH = ', ZLENGTH
               ENDIF
            ENDIF
         ENDIF

         RNP_PIC = ZERO
         CNP_PIC = ZERO 

         DO K = KSTART1, KEND1 
            DO J = JSTART1, JEND1
               DO I = ISTART1, IEND1 
                  IJK  = FUNIJK(I,J,K)
                  IF(.NOT.FLUID_AT(IJK)) CYCLE 
                  IF(EP_G(IJK).GE.1.d0-DIL_EP_s) CYCLE 
                  VOLIJK = VOL(IJK)
                  VOLIJK_UNCUT = DX(I)*DY(J)*DZ(K) 
                  DO M = 1, MMAX
                     REAL_PARTS(M) = 6.d0*EP_S(IJK,M)*VOLIJK/(PI*(D_p0(M)**3.d0))
                     !WRITE(*,*) 'NUMBER OF REAL PARTS  = ', REAL_PARTS(M)
                     IF(CONSTANTNPC) THEN 
                        COMP_PARTS(M) = NPC_PIC(M)
                        IF(CUT_CELL_AT(IJK)) COMP_PARTS(M) = INT(VOLIJK*real(COMP_PARTS(M))/VOLIJK_UNCUT)
                     ELSEIF(CONSTANTWT) THEN
                        COMP_PARTS(M) = MAX(1, INT(REAL_PARTS(M)/REAL(STATWT_PIC(M))))
                     ENDIF
                     
                     RNP_PIC(M) = RNP_PIC(M) + REAL_PARTS(M)

                     CNP_PIC(M) = CNP_PIC(M) + COMP_PARTS(M)
                     CNP_ARRAY(IJK,M) = COMP_PARTS(M)
                  ENDDO
                  
                  
               ENDDO
            ENDDO
         ENDDO
         PART_MPHASE(1:MMAX) = CNP_PIC(1:MMAX)
         PARTICLES = SUM(PART_MPHASE(1:MMAX))
         !WRITE(*,*) 'FROM pe =', mype, 'particles = ', particles 
         CALL global_all_sum(PARTICLES)
      ENDIF !  end if gener_part_config


! pradeep: parallel processing 
! the entire checking and setting up indices for desgridsearch
! moved to desgrod_mod
  
      CALL desgrid_check()

!! Ensure settings for grid based neighbor search method. This section
!! was placed here since the domain lengths (xlength, ylength, zlength)
!! are needed and their value is not gauranteed until after 
!! check_data_03 is called
!! ------------------------------------------------------------
!      IF (DES_NEIGHBOR_SEARCH .EQ. 4) THEN
!         MAX_DIAM = 2.0d0*MAX_RADIUS
!
!         TMP_FACTOR = 3.0d0*(MAX_DIAM)
!
!! If the search grid is undefined then set it to approximately 3 times
!! the maximum particle. Otherwise, check to see that the user set search
!! grid is at least greater than or equal to the maximum particle
!! diameter and warn the user if not         
!         IF (DESGRIDSEARCH_IMAX == UNDEFINED_I) THEN
!            DL_TMP = XLENGTH/TMP_FACTOR
!            DESGRIDSEARCH_IMAX = INT(DL_TMP)
!            IF (DESGRIDSEARCH_IMAX <= 0) DESGRIDSEARCH_IMAX = 1
!            WRITE(*,'(3X,A,I8)') &
!               'DESGRIDSEARCH_IMAX was set to ', DESGRIDSEARCH_IMAX
!         ELSE
!            DL_TMP = XLENGTH/DBLE(DESGRIDSEARCH_IMAX)
!            IF (DL_TMP < MAX_DIAM) THEN
!               WRITE(*,1037) 'x', 'x', 'i', 'i'                    
!               CALL MFIX_EXIT(myPE)
!            ENDIF
!         ENDIF
!         IF (DESGRIDSEARCH_JMAX == UNDEFINED_I) THEN
!            DL_TMP = YLENGTH/TMP_FACTOR
!            DESGRIDSEARCH_JMAX = INT(DL_TMP)
!            IF (DESGRIDSEARCH_JMAX <= 0) DESGRIDSEARCH_JMAX = 1
!            WRITE(*,'(3X,A,I8)') &
!               'DESGRIDSEARCH_JMAX was set to ', DESGRIDSEARCH_JMAX
!         ELSE
!            DL_TMP = YLENGTH/DBLE(DESGRIDSEARCH_JMAX)
!            IF (DL_TMP < MAX_DIAM) THEN
!               WRITE(*,1037) 'y', 'y', 'j', 'j'
!               CALL MFIX_EXIT(myPE)
!            ENDIF
!         ENDIF
!         IF (DIMN .EQ. 2) THEN
!            IF (DESGRIDSEARCH_KMAX == UNDEFINED_I) THEN
!               DESGRIDSEARCH_KMAX = 1
!            ELSEIF(DESGRIDSEARCH_KMAX /= 1) THEN
!               DESGRIDSEARCH_KMAX = 1            
!               WRITE(*,'(3X,A,I8)') &
!                  'DESGRIDSEARCH_KMAX was set to ', DESGRIDSEARCH_KMAX
!            ENDIF            
!         ELSE
!            IF (DESGRIDSEARCH_KMAX == UNDEFINED_I) THEN
!                DL_TMP = ZLENGTH/TMP_FACTOR
!                DESGRIDSEARCH_KMAX = INT(DL_TMP)
!                IF (DESGRIDSEARCH_KMAX <= 0) DESGRIDSEARCH_KMAX = 1
!                WRITE(*,'(3X,A,I8)') &
!                  'DESGRIDSEARCH_KMAX was set to ', DESGRIDSEARCH_KMAX
!            ELSE
!                DL_TMP = ZLENGTH/DBLE(DESGRIDSEARCH_KMAX)
!                IF (DL_TMP < MAX_DIAM) THEN
!                   WRITE(*,1037) 'z', 'z', 'k', 'k'
!                   CALL MFIX_EXIT(myPE)
!               ENDIF
!            ENDIF
!         ENDIF   ! end if/else dimn == 2
!
!         DESGS_IMAX2= DESGRIDSEARCH_IMAX+2
!         DESGS_JMAX2 = DESGRIDSEARCH_JMAX+2
!         IF (DIMN .EQ. 2) THEN
!            DESGS_KMAX2 = DESGRIDSEARCH_KMAX
!         ELSE
!            DESGS_KMAX2 = DESGRIDSEARCH_KMAX+2
!         ENDIF         
!
!
!! Variable that stores the particle in cell information (ID) on the
!! computational grid defined by cell/grid based search.  Similar to the
!! variable PIC but tailored for the grid based neighbor search option
!         ALLOCATE(DESGRIDSEARCH_PIC(DESGS_IMAX2,DESGS_JMAX2,DESGS_KMAX2))
!         DO K = 1,DESGS_KMAX2
!            DO J = 1,DESGS_JMAX2
!               DO I = 1,DESGS_IMAX2
!                  NULLIFY(DESGRIDSEARCH_PIC(I,J,K)%p) 
!               ENDDO 
!             ENDDO 
!          ENDDO                      
!
!      ENDIF   ! end if des_neighbor_search == 4
!! End checks if grid based neighbor search         
!! ------------------------------------------------------------

      
      if(dmp_log.and.debug_des) WRITE(unit_log,'(1X,A)')&
         '<---------- END CHECK_DES_DATA ----------'


      RETURN
 1000 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'DES should only be run using CARTESIAN coordinates',&
         /1X,70('*')/)
 1001 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'DES being run on multiple processors. Only serial runs',/10X,&
         'allowed',/1X,70('*')/)
 1002 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'Direction of periodicity not defined in mfix.dat',&
         /1X,70('*')/)

 1003 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'Invalid option for DES_NEIGHBOR_SEARCH in mfix.dat',/10X,&
         'Must be > 0 or < 5',/1X,70('*')/)

 1004 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'Spring constants KN or KN_W not specified in mfix.dat',&
         /1X,70('*')/)
 1005 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'WARNING: Tangential spring factors KT_FAC or KT_W_FAC not',/10X,&
         'specified in mfix.dat.  These factors will be defined in ',&
         'cfassign.f',/10X,'as 2/7.  See subroutine for references.',&
         /1X,70('*')/)
 1006 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'Change the value DES_COLL_MODEL in mfix.dat. Options are',/10X,&
         'leave it undefined to use the (default) linear spring-',/10X,&
         'dashpot model or set it to HERTZIAN for hertz model.',&
          /1X,70('*')/)
 1007 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'Friction coefficients MEW or MEW_W not specified in mfix.dat',&
         /1X,70('*')/)
 1008 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'Particle-particle restitution coefficient DES_EN_INPUT(M)',/10X,&
         'not specified in mfix.dat for interactions M = 1 to ',I5,&
          /1X,70('*')/)
 1009 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'Particle-wall restitution coefficients DES_EN_WALL_INPUT(M)',&
         /10X,'not specified in mfix.dat for interactions M= 1 to MMAX',&
         /1X,70('*')/)
 1010 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'WARNING: Tangential damping factors DES_ETAT_FAC not ',&
         'specified in',/10X, 'mfix.dat. This factor will be set in ',&
         'cfassign.f as 1/2. See subroutine for references.',&
         /1X,70('*')/)
 1011 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'WARNING: Tangential damping factors DES_ETAT_W_FAC not ',&
         'specified in mfix.dat',/10X,'This factor will be set in ',&
         'cfassign.f as 1/2. See subroutine for references.',&
         /1X,70('*')/)
 1012 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'Unphysical ( > 1 or < 0) values of DES_EN_INPUT(M)',&
         /1X,70('*')/)
 1013 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'Unphysical ( > 1 or < 0) values of DES_EN_WALL_INPUT(M)',&
         /1X,70('*')/)
 1014 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'Unphysical ( > 1 or < 0) values of friction coefficients',&
         /1X,70('*')/)
 1015 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'Values of DES_ETAT_FAC or DES_ETAT_W_FAC unphysical ',/10X,&
         '(< 0 or > 1) defined in mfix.dat',/1X,70('*')/)
 1016 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'Values of KT_FAC or KT_W_FAC unphysical (< 0 or > 1) ',/10X,&
         'defined in mfix.dat',/1X,70('*')/)

 1017 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'WARNING: A direction of periodicity was defined (i.e. ',/10X,&
         'DES_PERIODIC_WALLS_X, _Y or _Z=T) but DES_PERIODIC_WALLS ',&
         'was set to F.',/10X,&
         'The latter was forced to T for consistency.',/1X,70('*')/)
 1018 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'WARNING: nsquare neighbor search may be slow with periodic',/10X,&
         'boundaries.  Grid based search is recommended.',/1X,70('*')/)
 1019 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'octree and quadtree neighbor search methods do not',/10X,&
         'currently work with periodic boundaries. Change ',&
         'value of',/10X,'DES_NEIGHBOR_SEARCH in mfix.dat',/1X,70('*')/)

 1020 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'Wall value for Youngs modulus (EW_YOUNG) must be,'/10X,&
         'specified in mfix.dat',/1X,70('*')/)
 1021 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'Wall value for Poissons ratio (VW_POISSON) must be',/10X,&
         'specified in mfix.dat',/1X,70('*')/)
 1022 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'Youngs modulus (E_YOUNG) must be specified in mfix.dat',/10X,&
         'for all particle types M = 1 to MMAX',/1X,70('*')/)
 1023 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'Poissons ratio (V_POISSON) must be specified in mfix.dat',/10X,&
         'for all particle types M = 1 to MMAX',/1X,70('*')/)
 1024 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'Particle-particle tangential restitution coefficient',/10X,&
         'DES_ET_INPUT(M) must be specified in mfix.dat for all',/10X,&
         'interactions M = 1 to ',I5,/1X,70('*')/)
 1025 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'Unphysical ( > 1 or < 0) values of DES_ET_INPUT(M)',/1X,70('*')/)
 1026 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'Particle-wall tangential restitution coefficient',/10X,&
         'DES_ET_WALL_INPUT(M) must be specified in mfix.dat for all',/10X,&
         'wall interactions M = 1 to MMAX',/1X,70('*')/)
 1027 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'Unphysical ( > 1 or < 0) values of DES_ET_WALL_INPUT(M)',/1X,70('*')/)
 1028 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'WARNING: specified values of KN and KN_W are not used with',/10X,&
         'Hertz collision model',/1X,70('*')/)
 1029 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'WARNING: specified values of KT_FAC and KT_W_WALL are not',/10X,&
         'used with Hertz collision model',/1X,70('*')/)
 1030 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'WARNING: specified values of DES_ETAT_FAC and',/10X,&
         'DES_ETAT_W_FAC are not used with Hertz collision model',/1X,70('*')/)
 1031 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'WARNING: specified values of DES_ET_INPUT(M) are not used',/10X,&
         'with default collision model',/1X,70('*')/)
 1032 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'WARNING: specified values of DES_ET_WALL_INPUT(M) are not',/10X,&
         'used with default collision model',/1X,70('*')/)
 1033 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'WARNING: Values of VW_POISSON OR V_POISSON unphysical',/10X,&
         '(> 0.50 or =< -1.d0) defined in mfix.dat',/1X,70('*')/)

 1034 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'Invalid option for DES_INTG_METHOD in mfix.dat.  Must be',&
         /10X,'EULER (default/undefined) or ADAMS_BASHFORTH',&
         /1X,70('*')/)

 1035 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'USE_COHESION and corresponding code has been disabled as',/10X,&
         'this feature has not verified with the current DEM code. To',&
         /10X,'activate this feature comment out this check and ',&
         'proceed at',/10X,'your own risk',/1X,70('*')/)

 1036 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         '2D coupled simulation with a particle diameter > ZLENGTH.',/10X,&
         'This will create problems for calculations of void ',&
         'fraction. Check',/10X, 'mfix.dat file.',/1X,70('*')/)

 1037 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
          'The neighbor search grid is too fine in the ',A, &
          '-direction',/10X,'with a particle diameter > ',A, &
          'length/desgridsearch_',A,'max. This will',/10X,'create ',&
          'problems for the search method and detecting neighbors',/10X,&
          'Decrease desgridsearch_',A,'max in mfix.dat to coarsen ',&
          'grid.',/1X,70('*')/)
  
 1038 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
          'The only option for DES_OUTPUT_DATA is TECPLOT',&
          /1X,70('*')/)

 1039 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'DES can only be run in XY plane in 2D.',/1X,70('*')/)
 1040 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'WARNING: DIMN left unspecified.  Since NO_K is T ',/10X,&
         'DIMN was set to 2. 2D DES running in XY plane.',/1X,70('*')/)
 1041 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'Geometry dimension DIMN not specified.',/1X,70('*')/)
 1042 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'Physical dimension DIMN cannot be > 3.',/1X,70('*')/)

         
 1043 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'D_P0 must be defined and >0 in mfix.dat for M=1,',&
         'DES_MMAX.',/1X,70('*')/)
 1044 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'Too many D_P0 are defined in mfix.dat for given ',&
         'DES_MMAX.',/1X,70('*')/)

 1045 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'DES_NEIGHBOR_SEARCH set to ', I2, ' ',A,/1X,70('*')/)


 1047 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'Values of DES_EPS_XSTART and/or DES_EPS_YSTART must be ',&
         'specified',/10X,'when GENER_PART_CONFIG.',/1X,70('*')/)

 1048 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'Value of DES_EPS_ZSTART must be specified when DIMN = 3',&
         /10X,'and GENER_PART_CONFIG.',/1X,70('*')/)

 1049 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'VOL_FRAC(M) must be defined in mfix.dat for M=1, MMAX ',&
         'when',/10X,'GENER_PART_CONFIG',/1X,70('*')/)

 1050  FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
          'Unphysical ( > 1-EP_STAR or < 0) values of VOL_FRAC(M) ',&
          'set in',/10X,'mfix.dat for GENER_PART_CONFIG option',&
          /1X,70('*')/)
 1051  FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
          'Unphysical ( > 1 or < 0) values of VOL_FRAC(M) ',&
          'set in',/10X,'mfix.dat for GENER_PART_CONFIG option',&
          /1X,70('*')/)
 1052  FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
          'Total solids volume fraction should not exceed 1-EP_STAR ',&
          'for',/10X,'GENER_PART_CONFIG option where EP_STAR = ',&
          G12.5,/10X,'and the total solids volume fraction = ',G12.5, &
          /1X,70('*')/)
 1053  FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
          'Total solids volume fraction should not exceed 1 for',/10X,&
          'GENER_PART_CONFIG option where the total solids volume ',&
          'fraction = ',G12.5,/1X,70('*')/)

 1054 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
          'WARNING: zlength or dz(1) is used to calculate the ',&
          'number',/10X,'of particles in the 2D simulation when ',&
          'GENER_PART_CONFIG is T and DIMN = 2.',/10X,'This depth ',&
          'does not equal D_P0(1).',/1X,70('*'))

 1055 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
          'WARNING: zlength or dz(1) is used to calculate the ',&
          'number',/10X,'of particles in the 2D simulation when ',&
          'GENER_PART_CONFIG is T and DIMN = 2.',/1X,70('*'))

 1056 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
          'WARNING: TSUJI_DRAG option only available with ',&
          'SYAM_OBRIEN',/10X,'drag model.',/1X,70('*'))

     
 1060 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'Only the grid based search option is allowed when using',&
         'using',/10X,'Lees & Edwards BC.',/1X,70('*')/)

 1061 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'Invalid option for DES_LE_SHEAR_DIR in mfix.dat. When',/10X,&
         'DIMN=2 shear options are DUDY or DVDX',/1X,70('*')/)

 1062 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'Invalid option for DES_LE_SHEAR_DIR in mfix.dat. When',/10X,&
         'DIMN=3 shear options are DUDY, DUDZ, DVDX, DVDZ, DWDX or',&
         'DWDY.',/1X,70('*')/)
     
 1063 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'WARNING: DES_PERIODIC_WALLS set to false when DES_LE_BC.',&
         /10X,'DES_LE_BC implies periodic walls, however, the ',&
         'periodicity is handled',/10X, 'independently of ',&
         'DES_PERIODIC_WALLS option and so it is shut off.',&
         /1X,70('*')/)

 1064 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'DES_CONTINUUM_COUPLED cannot be true when using ',&
         'DES_LE_BC.',/1X,70('*')/)
 
 1070 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'MMAX not specified or unphysical in mfix.dat.' ,&
         /1X,70('*')/)

     
 2001    FORMAT(/1X,70('*')//' From: mfix.f',/' Message: ',&
         'Looks like a 3-D case but dimn equal to ',i1,/1X,70('*')/)

         END SUBROUTINE CHECK_DES_DATA
