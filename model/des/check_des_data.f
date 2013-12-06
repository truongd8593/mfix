!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CHECK_DES_DATA                                          C
!  Purpose: Check user input data                                      C
!                                                                      C
!                                                                      C
!  Author: S. Benyahia                                Date: 10-Nov-08  C
!  Reviewer:                                                           C
!  Comments: Most all user's input data are checked here               C
!  Revision: Some of the checks made in des_allocate_arrays are        C
!            moved here. In addition write statments are now made to   C
!            log file                                                  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CHECK_DES_DATA

!-----------------------------------------------
! Modules 
!-----------------------------------------------      
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
      USE qmom_kinetic_equation
      USE ic 
      IMPLICIT NONE
!-----------------------------------------------
! Local Variables
!-----------------------------------------------      
      INTEGER :: M
      INTEGER :: I, J, K, IJK
      LOGICAL :: FLAG_WARN
! domain volume      
      DOUBLE PRECISION :: VOL_DOMAIN
! the maximum and minimum specified particle diameter 
      DOUBLE PRECISION :: MAX_DIAM, MIN_DIAM
! for gener_part_config, the total solids volume fraction
      DOUBLE PRECISION :: TOT_VOL_FRAC
! number of real and comp. particles in a cell 
! for MP-PIC case 
      DOUBLE PRECISION REAL_PARTS(DIM_M), COMP_PARTS(DIM_M)
! volume of the cell 
      DOUBLE PRECISION :: VOLIJK, VOLIJK_UNCUT
! temp variable for solids volume fraction and particle count 
      DOUBLE PRECISION :: EP_SM, PARTS_TEMP 
!IC region count with non-zero solids 
      INTEGER :: IC_COUNT, ICV

!-----------------------------------------------
! Functions
!-----------------------------------------------
      LOGICAL, EXTERNAL :: COMPARE
!-----------------------------------------------
! Include statement functions
!----------------------------------------------- 
      INCLUDE 'function.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'ep_s2.inc'
!----------------------------------------------- 


      IF(DMP_LOG.AND.DEBUG_DES) WRITE(UNIT_LOG,'(1X,A)')&
         '---------- START CHECK_DES_DATA ---------->'


! want separate flag to toggle between screen output for the discrete modules
      PRINT_DES_SCREEN = (myPE.EQ.pe_IO)
      IF(DMP_LOG) WRITE(UNIT_LOG,1001) PRINT_DES_SCREEN
      IF(DMP_LOG) WRITE(*,1001) PRINT_DES_SCREEN

      IF (.NOT.DES_CONTINUUM_HYBRID) THEN
! MMAX, D_p0 and RO_s0 are to be strictly associated with the continuum
! model and are no longer to be used in the DEM. However, when not using
! DES_CONTINUUM_HYBRID the current DEM code still assumes that the user
! specifies MMAX, D_P0 and RO_S0.  These are then linked to their
! respective DES variables.  Note that valid values of MMAX, D_P0 and 
! RO_s0 are ensured by check_data_04.
         DES_MMAX = MMAX
         DO M = 1, DES_MMAX
            DES_D_p0(M) = D_p0(M)
            DES_RO_s(M) = RO_s0(M)
         ENDDO
      ENDIF

      IF (DES_MMAX == UNDEFINED_I) THEN
         IF (DMP_LOG) WRITE(UNIT_LOG, 1090)
         CALL MFIX_EXIT(myPE)
      ELSE
! The following checks on DES_mmax, density and diameter are similar to
! those in check_data_04.  If .not.des_continuum_hybrid then these
! checks are basically redundant, however, they will facilitate a move
! to make DEM operate separately from MFIX's check_data_04.  
! ---------------------------------------------------------------->>>    
! Check mmax: 
         IF (DES_MMAX<0 .OR. DES_MMAX>DIM_M) THEN 
            IF(DMP_LOG) WRITE(UNIT_LOG, 1070)
            CALL MFIX_EXIT(myPE)
         ENDIF
! Check for valid diameter values. 
! Valid diameters are needed: 1) when using gener_part_config the
! diameters are used to identify the number of particles in a given
! solids phase and 2) to identify which solids phase each particle
! belongs in (this sorting is conducted in particles_in_cell)
         DO M = 1,DES_MMAX
            IF (DES_D_P0(M)<ZERO .OR. DES_D_P0(M)==UNDEFINED) THEN
               IF(DMP_LOG) WRITE(UNIT_LOG, 1043)
               CALL MFIX_EXIT(myPE)
            ENDIF
         ENDDO
         DO M = DES_MMAX+1, DIM_M
            IF (DES_D_P0(M) /= UNDEFINED) THEN
               IF(DMP_LOG) WRITE(UNIT_LOG, 1044)
               CALL MFIX_EXIT(myPE)
            ENDIF
         ENDDO
! Check for valid density values
         DO M = 1,DES_MMAX
            IF (DES_RO_S(M)<ZERO .OR. DES_RO_S(M)==UNDEFINED) THEN
               IF(DMP_LOG) WRITE(UNIT_LOG, 1071)
               CALL MFIX_EXIT(myPE)
            ENDIF
         ENDDO
         DO M = DES_MMAX+1, DIM_M
            IF (DES_RO_S(M) /= UNDEFINED) THEN
               IF(DMP_LOG) WRITE(UNIT_LOG, 1072)
               CALL MFIX_EXIT(myPE)
            ENDIF
         ENDDO

      ENDIF   ! end if/else(des_mmax==undefined_i)
! End checks on DES_MMAX, density and diameter      
! ----------------------------------------------------------------<<<

! Determine the maximum particle size in the system (MAX_RADIUS), which
! in turn is used for various tasks        
      MAX_DIAM = ZERO
      MIN_DIAM = LARGE_NUMBER
      DO M = 1,DES_MMAX
         MAX_DIAM = MAX(MAX_DIAM, DES_D_P0(M))
         MIN_DIAM = MIN(MIN_DIAM, DES_D_P0(M))
       ENDDO
      MAX_RADIUS = 0.5d0*MAX_DIAM
      MIN_RADIUS = 0.5d0*MIN_DIAM


! Set flags for energy equations
         CALL CHECK_DES_ENERGY
! Check thermodynamic properties of discrete solids.
         CALL CHECK_DES_THERMO


      IF(COORDINATES == 'CYLINDRICAL') THEN
         IF(DMP_LOG) WRITE (UNIT_LOG, 1000)
         CALL MFIX_EXIT(myPE)
      ENDIF

! Verify that there are no internal obstacles.
      IF(.NOT.CARTESIAN_GRID) THEN
         DO K = KSTART1, KEND1 
            DO J = JSTART1, JEND1
               DO I = ISTART1, IEND1 
                  IJK  = FUNIJK(I,J,K)
                  IF(.NOT.FLUID_AT(IJK)) THEN
                     IF(DMP_LOG) WRITE(*,2010)
                     CALL MFiX_EXIT(myPE)
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDIF


! Check dimension
      IF(NO_I.OR.NO_J) THEN ! redundant, this check is made in check_data_03
         IF(DMP_LOG) WRITE(UNIT_LOG,1039)
         CALL MFIX_EXIT(myPE)
      ENDIF
      IF(DIMN == UNDEFINED_I) THEN
         IF(NO_K) THEN
            DIMN = 2 
            IF(DMP_LOG) WRITE(UNIT_LOG, 1040)
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
! this check will make sure DIMN is correctly set to 3
! for a truly 3-D case. 
         IF(DIMN.NE.3) THEN 
            IF(DMP_LOG)  WRITE(*,2001) DIMN 
            CALL MFIX_EXIT(myPE)
         ENDIF
      ENDIF


! Check settings on cohesion model
! ---------------------------------------------------------------->>>
      IF(USE_COHESION) THEN
         IF (SQUARE_WELL .AND. VAN_DER_WAALS) THEN
            IF(DMP_LOG) WRITE(UNIT_LOG, 1081)
            CALL MFIX_EXIT(myPE)            
         ELSEIF(.NOT.SQUARE_WELL .AND. .NOT.VAN_DER_WAALS) THEN
            IF(DMP_LOG) WRITE(UNIT_LOG,1082)
            CALL MFIX_EXIT(myPE)
         ELSEIF (SQUARE_WELL) THEN     
! note that square_well cohesion code may not function correctly with
! all aspects of the current DEM code so proceed at your own risk. 
            IF(DMP_LOG) WRITE(UNIT_LOG, 1080)
            IF (MASTER_WELL_DEPTH .EQ. UNDEFINED .AND. &
                MASTER_WALL_WELL_DEPTH .EQ. UNDEFINED .AND. &
                RADIUS_RATIO .EQ. UNDEFINED .AND. &
                WALL_RADIUS_RATIO .EQ. UNDEFINED) THEN
               IF(DMP_LOG) WRITE(UNIT_LOG,1083)
               CALL MFIX_EXIT(myPE)
            ENDIF
         ELSEIF (VAN_DER_WAALS) THEN
            IF (VDW_INNER_CUTOFF .EQ. UNDEFINED .AND. &
                VDW_OUTER_CUTOFF .EQ. UNDEFINED .AND. &
                HAMAKER_CONSTANT .EQ. UNDEFINED) THEN
               IF(DMP_LOG) WRITE(UNIT_LOG,1084)
               CALL MFIX_EXIT(myPE)         
            ENDIF
            IF (WALL_VDW_INNER_CUTOFF .EQ. UNDEFINED .AND. &
                WALL_HAMAKER_CONSTANT .EQ. UNDEFINED) THEN
               IF(DMP_LOG) WRITE(UNIT_LOG,1085)
               CALL MFIX_EXIT(myPE)
            ENDIF
            IF (FACTOR_RLM < &
                1.d0+VDW_OUTER_CUTOFF/(2.d0*MAX_RADIUS)) THEN
               IF(DMP_LOG) WRITE(UNIT_LOG,1100)
               IF(DMP_LOG) WRITE(*,1100)
            ENDIF
         ENDIF                
      ELSE
! override any of the following settings if cohesion not used              
          WALL_VDW_OUTER_CUTOFF = ZERO
          SQUARE_WELL = .FALSE.
          VAN_DER_WAALS = .FALSE.
      ENDIF
! end check settings on cohesion model
! ----------------------------------------------------------------<<<


! Check settings on cluster identification
! ---------------------------------------------------------------->>>

      IF(DES_CALC_CLUSTER) THEN
         IF(CLUSTER_LENGTH_CUTOFF .EQ. UNDEFINED) THEN
            IF(DMP_LOG) WRITE(UNIT_LOG,1101)
            IF(DMP_LOG) WRITE(*,1101)
            CALL MFIX_EXIT(myPE)
         ENDIF
         IF(FACTOR_RLM < &
            1.d0+CLUSTER_LENGTH_CUTOFF/(2.d0*MAX_RADIUS)) THEN
            IF(DMP_LOG) WRITE(UNIT_LOG,1102)
            IF(DMP_LOG) WRITE(*,1102)
            CALL MFIX_EXIT(myPE)
         ENDIF
      ENDIF
! ----------------------------------------------------------------<<<


! des_periodic_walls must be true if any wall is periodic      
      IF ( (DES_PERIODIC_WALLS_X .OR. DES_PERIODIC_WALLS_Z .OR. &
            DES_PERIODIC_WALLS_Z) .AND. .NOT.DES_PERIODIC_WALLS) THEN
         DES_PERIODIC_WALLS = .TRUE.
         IF(DMP_LOG) WRITE(UNIT_LOG, 1017)
      ENDIF

! Check for valid neighbor search option      
      IF(DES_NEIGHBOR_SEARCH.EQ.1) THEN
         IF(DMP_LOG) WRITE(unit_log, 1045) 1, 'N-SQUARE'
      ELSEIF(DES_NEIGHBOR_SEARCH.EQ.2) THEN
         IF(DMP_LOG) WRITE(unit_log, 1045) 2, 'QUADTREE'
      ELSEIF(DES_NEIGHBOR_SEARCH.EQ.3) THEN
         IF(DMP_LOG) WRITE(unit_log, 1045) 3, 'OCTREE'
      ELSEIF(DES_NEIGHBOR_SEARCH.EQ.4) THEN
         IF(DMP_LOG) WRITE(unit_log, 1045) 4, 'GRID BASED'
      ELSE
         IF(DMP_LOG) WRITE(UNIT_LOG, 1003)
         CALL MFIX_EXIT(myPE)         
      ENDIF

! Periodicity does not work with most neighbor search options
      IF(DES_PERIODIC_WALLS) THEN
         IF(.NOT.DES_PERIODIC_WALLS_X .AND. .NOT.DES_PERIODIC_WALLS_Y .AND. &
            .NOT.DES_PERIODIC_WALLS_Z) THEN
            IF(DMP_LOG) WRITE (UNIT_LOG, 1002)
            CALL MFIX_EXIT(myPE)
         ENDIF
         IF (DES_NEIGHBOR_SEARCH .EQ. 1) THEN
            IF(DMP_LOG) WRITE (UNIT_LOG, 1018)
         ENDIF
         IF (DES_NEIGHBOR_SEARCH .EQ. 2 .OR. &
          DES_NEIGHBOR_SEARCH .EQ. 3) THEN
            IF(DMP_LOG) WRITE (UNIT_LOG, 1019)
            CALL MFIX_EXIT(myPE)
         ENDIF
      ENDIF


! Lees Edwards BC functionality has been lost in current DEM code
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

! Check the output file format 
      IF (DES_OUTPUT_TYPE /= UNDEFINED_C) THEN
         IF(TRIM(DES_OUTPUT_TYPE) /= 'TECPLOT') THEN
            IF(DMP_LOG) WRITE (*, 1038)
            CALL MFIX_EXIT(myPE)   
         ENDIF
      ENDIF

! Check for valid integration method
      IF (TRIM(DES_INTG_METHOD) == 'ADAMS_BASHFORTH') THEN
         INTG_ADAMS_BASHFORTH = .TRUE.
         INTG_EULER = .FALSE.
      ELSEIF(TRIM(DES_INTG_METHOD) == 'EULER') THEN
         INTG_ADAMS_BASHFORTH = .FALSE.
         INTG_EULER = .TRUE.
! stop if the specified integration method is unavailable
      ELSE
         IF(DMP_LOG) WRITE (UNIT_LOG, 1034)
         CALL MFIX_EXIT(myPE)
      ENDIF

! Check for valid integration method
      IF (TRIM(DES_INTG_METHOD) /= 'ADAMS_BASHFORTH' .AND. &
          TRIM(DES_INTG_METHOD) /= 'EULER') THEN
! stop if the specified integration method is unavailable
         IF(DMP_LOG) WRITE (UNIT_LOG, 1034)
         CALL MFIX_EXIT(myPE)
      ENDIF



! Check settings for collision models
! ---------------------------------------------------------------->>>      
      IF(.NOT.MPPIC) THEN 
         IF (DES_COLL_MODEL /= UNDEFINED_C) THEN      
            IF (TRIM(DES_COLL_MODEL) .NE. 'HERTZIAN') THEN
               IF(DMP_LOG) WRITE(UNIT_LOG, 1006)
               CALL MFIX_EXIT(myPE)
            ENDIF
         ENDIF

         IF (TRIM(DES_COLL_MODEL) .EQ. 'HERTZIAN') THEN
! check young's modulus and poisson ratio
            IF(EW_YOUNG == UNDEFINED ) THEN
               IF(DMP_LOG) WRITE (UNIT_LOG, 1020)
               CALL MFIX_EXIT(myPE)
            ENDIF
            IF(VW_POISSON == UNDEFINED) THEN
               IF(DMP_LOG) WRITE (UNIT_LOG, 1021)
               CALL MFIX_EXIT(myPE)
            ELSE
              IF (VW_POISSON > 0.5d0 .OR. VW_POISSON <= -ONE) THEN
                 IF(DMP_LOG) WRITE (UNIT_LOG, 1033)
              ENDIF
            ENDIF
            DO M = 1, DES_MMAX
               IF(E_YOUNG(M) == UNDEFINED) THEN
                  IF(DMP_LOG) WRITE (UNIT_LOG, 1022)
                  CALL MFIX_EXIT(myPE)
               ENDIF
               IF(V_POISSON(M) == UNDEFINED) THEN
                  IF(DMP_LOG) WRITE (UNIT_LOG, 1023)
                  CALL MFIX_EXIT(myPE)
               ELSE
                 IF(V_POISSON(M) > 0.5d0 .OR. V_POISSON(M) <= -ONE) THEN
                    IF(DMP_LOG) WRITE (UNIT_LOG, 1033)
                 ENDIF
               ENDIF
            ENDDO
! check particle-particle tangential restitution coefficient      
            DO M = 1, DES_MMAX+DES_MMAX*(DES_MMAX-1)/2
               IF(DES_ET_INPUT(M) == UNDEFINED) THEN
                  IF(DMP_LOG) WRITE (UNIT_LOG, 1024) DES_MMAX+DES_MMAX*(DES_MMAX-1)/2
                  CALL MFIX_EXIT(myPE)
               ENDIF
            ENDDO
            DO M = 1, DES_MMAX+DES_MMAX*(DES_MMAX-1)/2
               IF(DES_ET_INPUT(M) > ONE .OR. DES_ET_INPUT(M) < ZERO) THEN
                  IF(DMP_LOG) WRITE (UNIT_LOG, 1025)
                  CALL MFIX_EXIT(myPE)
               ENDIF
            ENDDO
! check particle-wall tangential restitution coefficient
            DO M = 1, DES_MMAX
               IF(DES_ET_WALL_INPUT(M) == UNDEFINED) THEN
                  IF(DMP_LOG) WRITE (UNIT_LOG, 1026)
                  CALL MFIX_EXIT(myPE)
               ENDIF
            ENDDO
            DO M = 1, DES_MMAX
               IF(DES_ET_WALL_INPUT(M) > ONE .OR. DES_ET_WALL_INPUT(M) < ZERO) THEN
                  IF(DMP_LOG) WRITE (UNIT_LOG, 1027)
                  CALL MFIX_EXIT(myPE)
               ENDIF
            ENDDO
! if following are assigned warn user they are discarded
            IF(KN .NE. UNDEFINED .OR. KN_W .NE. UNDEFINED) THEN 
               IF(DMP_LOG) WRITE (UNIT_LOG, 1028)
            ENDIF
            IF(KT_FAC .NE. UNDEFINED .OR. KT_W_FAC .NE. UNDEFINED) then 
               IF(DMP_LOG) WRITE (UNIT_LOG, 1029)
            ENDIF
            IF(DES_ETAT_FAC .NE. UNDEFINED .OR. &
               DES_ETAT_W_FAC .NE. UNDEFINED) then 
               IF(DMP_LOG) WRITE (UNIT_LOG, 1030)
            ENDIF


         ELSE  ! default linear spring-dashpot model

! check normal spring constants 
            IF(KN == UNDEFINED .OR. KN_W == UNDEFINED) THEN
               IF(DMP_LOG) WRITE (UNIT_LOG, 1004)
               CALL MFIX_EXIT(myPE)
            ENDIF
! check for tangential spring constant factors         
            IF(KT_FAC == UNDEFINED .OR. KT_W_FAC == UNDEFINED) THEN
               IF(DMP_LOG) WRITE (UNIT_LOG, 1005)
            ENDIF
            IF(KT_FAC .NE. UNDEFINED) THEN
               IF(KT_FAC > ONE .OR. KT_FAC < ZERO) THEN
                  IF(DMP_LOG) WRITE (UNIT_LOG, 1016)
                  CALL MFIX_EXIT(myPE)
               ENDIF
            ENDIF
            IF(KT_W_FAC .NE. UNDEFINED) THEN
               IF(KT_W_FAC > ONE .OR. KT_W_FAC < ZERO) THEN
                  IF(DMP_LOG) WRITE (UNIT_LOG, 1016)
                  CALL MFIX_EXIT(myPE)
               ENDIF
            ENDIF
! check for tangential damping factor
            IF(DES_ETAT_FAC == UNDEFINED) THEN
               IF(DMP_LOG) WRITE (UNIT_LOG, 1010)
            ELSEIF(DES_ETAT_FAC > ONE .OR. DES_ETAT_FAC < ZERO) THEN
               IF(DMP_LOG) WRITE (UNIT_LOG, 1015)
               CALL MFIX_EXIT(myPE)
            ENDIF
            IF(DES_ETAT_W_FAC == UNDEFINED) THEN
               IF(DMP_LOG) WRITE (UNIT_LOG, 1011)
            ELSEIF(DES_ETAT_W_FAC > ONE .OR. DES_ETAT_W_FAC < ZERO) THEN
               IF(DMP_LOG) WRITE (UNIT_LOG, 1015)
               CALL MFIX_EXIT(myPE)
            ENDIF
! if following are assigned warn user they are discarded
            FLAG_WARN = .FALSE.
            DO M = 1, DES_MMAX+DES_MMAX*(DES_MMAX-1)/2
               IF(DES_ET_INPUT(M) .NE. UNDEFINED) FLAG_WARN = .TRUE.
            ENDDO
            IF (FLAG_WARN) WRITE(UNIT_LOG,1031)
            FLAG_WARN = .FALSE.
            DO M = 1, DES_MMAX
               IF(DES_ET_WALL_INPUT(M) .NE. UNDEFINED) FLAG_WARN = .TRUE.
            ENDDO
            IF (FLAG_WARN) WRITE(UNIT_LOG,1032)
            FLAG_WARN = .FALSE.

         ENDIF   ! end if/else(des_coll_model.eq.hertzian)

! check particle-particle normal restitution coefficient      
         DO M = 1, DES_MMAX+DES_MMAX*(DES_MMAX-1)/2
            IF(DES_EN_INPUT(M) == UNDEFINED) THEN
               IF(DMP_LOG) WRITE (UNIT_LOG, 1008) DES_MMAX+DES_MMAX*(DES_MMAX-1)/2
               CALL MFIX_EXIT(myPE)
            ENDIF
         ENDDO
         DO M = 1, DES_MMAX+DES_MMAX*(DES_MMAX-1)/2
            IF(DES_EN_INPUT(M) > ONE .OR. DES_EN_INPUT(M) < ZERO) THEN
               IF(DMP_LOG) WRITE (UNIT_LOG, 1012)
               CALL MFIX_EXIT(myPE)
            ENDIF
         ENDDO

! check particle-wall normal restitution coefficient (needed by all
! current collision models)
! MPPIC has separate key words now for the particle-wall interaction
         DO M = 1, DES_MMAX
            IF(DES_EN_WALL_INPUT(M) == UNDEFINED) THEN
               IF(DMP_LOG) WRITE (UNIT_LOG, 1009)
               CALL MFIX_EXIT(myPE)
            ENDIF
         ENDDO
         DO M = 1, DES_MMAX
            IF(DES_EN_WALL_INPUT(M) > ONE .OR. DES_EN_WALL_INPUT(M) < ZERO) THEN
               IF(DMP_LOG) WRITE (UNIT_LOG, 1013)
               CALL MFIX_EXIT(myPE)
            ENDIF
         ENDDO

! check coefficient friction 
         IF(MEW /= UNDEFINED .AND. MEW_W /= UNDEFINED) THEN
            IF (MEW> ONE .OR. MEW_W > ONE .OR. &
                MEW < ZERO .OR. MEW_W < ZERO) THEN
               IF(DMP_LOG) WRITE (UNIT_LOG, 1014)
               CALL MFIX_EXIT(myPE)
            ENDIF
         ELSE
            WRITE(UNIT_LOG, 1035)
            CALL MFIX_EXIT(myPE)
         ENDIF

      ENDIF   ! end if (.NOT.MPPIC)
! End check settings for collision models
! ----------------------------------------------------------------<<<



      IF (DES_CONTINUUM_HYBRID) THEN
! des_continuum_coupled must be true if des_continuum_hybrid              
         IF(.NOT.DES_CONTINUUM_COUPLED) THEN
            WRITE(UNIT_LOG, 1094)
            CALL MFIX_EXIT(myPE)
         ENDIF
! DES_CONTINUUM_HYBRID does not work with GHD or QMOMK
         IF (TRIM(KT_TYPE)=='GHD') THEN
            IF(DMP_LOG) WRITE(UNIT_LOG, 1091)
            CALL MFIX_EXIT(myPE)
         ENDIF
         IF (QMOMK) THEN
            IF(DMP_LOG) WRITE(UNIT_LOG, 1092)
            CALL MFIX_EXIT(myPE)
         ENDIF
         IF (C_F==UNDEFINED .AND. MMAX>0 .AND. DES_MMAX >0) THEN
            IF(DMP_LOG) WRITE(UNIT_LOG, 1095)
            CALL MFIX_EXIT(myPE)
         ENDIF                 
      ENDIF

      IF (MPPIC .AND. DES_CONTINUUM_HYBRID) THEN
         IF(DMP_LOG) WRITE(UNIT_LOG, 1093)
         CALL MFIX_EXIT(myPE)
      ENDIF         
      

      IF(.NOT.DES_CONTINUUM_COUPLED) THEN
! Overwrite user's input in case of DEM (no fluid)
         DES_INTERP_ON = .FALSE. 
      ENDIF


! Check that the depth of the simulation in 2D exceeds the largest 
! particle size to ensure correct calculation of volume fraction.  This
! is important for coupled simulations (not essential for pure granular
! simulations)
      IF (DES_CONTINUUM_COUPLED .AND. DIMN == 2) THEN
         IF (2.0d0*MAX_RADIUS > ZLENGTH) THEN
            IF(DMP_LOG) WRITE(UNIT_LOG, 1036)
            CALL MFIX_EXIT(myPE)
         ENDIF
      ENDIF


! If run_type is not 'NEW' then force the gener_part_config to .false. 
! This will prevent overwriting over the number of particles which could
! have potentially changed depending on inlet/outlet      
      IF(TRIM(RUN_TYPE) .NE. 'NEW'.AND.GENER_PART_CONFIG) THEN 
         GENER_PART_CONFIG = .FALSE. 
         IF(DMP_LOG) WRITE(UNIT_LOG, 1037)
         IF(PRINT_DES_SCREEN) WRITE(*, 1037)
      ENDIF
     
       
      IF(gener_part_config) then 
!Check if any of the deprecated flags are not present in the input file 
         if(des_eps_xstart.ne.undefined.or.des_eps_ystart.ne.undefined.or.des_eps_xstart.ne.undefined) then 
            if (PRINT_DES_SCREEN) write(*, 1056) 
            if (dmp_log) write(unit_log, 1056) 
            call mfix_exit(mype) 
         endif
         
         DO M = 1, DES_MMAX
            if(vol_frac(M).ne.undefined) then 
               if (print_des_screen) write(*, 1057) 
               if (dmp_log) write(unit_log, 1057) 
               call mfix_exit(mype) 
            endif
         enddo
      endif

! If gener_part_config ensure various quantities are defined and valid
! ---------------------------------------------------------------->>>
      IF(TRIM(RUN_TYPE).EQ. 'NEW' .AND. GENER_PART_CONFIG .AND.&
        .NOT.MPPIC) THEN 

! Determine the domain volume which is used to calculate the total
! number of particles and the number of particles in each phase. 
! Values of DZ(1) or zlength are guaranteed at this point due to
! check_data_03. If the user left both undefined and NO_K = .T., then
! they are set to ONE. If dz(1) is undefined but zlength is defined,
! then dz(1) is set to zlength (and vice versa).  If both are defined
! they must be equal.
         IF(DIMN.EQ.2) THEN 
            IF (DES_MMAX.EQ.1) THEN
! Warn the user if the domain depth is not equal to the particle
! diameter as it may cause problems for coupled simulations.
! The user should also be aware of this when interpreting
! volume/void fraction calculations (including bulk density).
               IF(.NOT.COMPARE(ZLENGTH,DES_D_P0(1))) THEN
                  IF(DMP_LOG) WRITE(UNIT_LOG, 1054)
                  WRITE(*,1054)
               ENDIF
            ELSE
! Let the user know basis of depth dimension for calculating number of
! particles. this will also be important when considering volume/void
! fraction calculations.
               IF(DMP_LOG) WRITE(UNIT_LOG, 1055)
               IF(print_des_screen) WRITE(*,1055)
            ENDIF
         ENDIF

         CALL CHECK_DES_IC_EPS 

         
         PARTICLES = 0
         IC_COUNT = 0 
         DO ICV = 1, DIMENSION_IC 
            IF(.not.ic_defined(icv)) cycle 

            PART_MPHASE_BYIC(ICV,:) = 0

            IF (IC_EP_G(ICV).lt.ONE) THEN 
               IC_COUNT = IC_COUNT + 1 
               DO M = 1, DES_MMAX
                  EP_SM = IC_ROP_S(ICV,M)/DES_RO_s(M)
                  PARTS_TEMP = FLOOR((6.D0*EP_SM*VOL_IC_REGION(ICV))/&
                  (PI*(DES_D_P0(M)**3)))
                  PART_MPHASE_BYIC(ICV, M) = PARTS_TEMP
                  PARTICLES = PARTICLES + PARTS_TEMP
               ENDDO
            ENDIF
         ENDDO
    

! Local reporting for the user 

         IF(DMP_LOG)        WRITE(UNIT_LOG,2015) IC_COUNT, DES_MMAX
         IF(PRINT_DES_SCREEN)  WRITE(*,2015) IC_COUNT, DES_MMAX

         DO ICV = 1, DIMENSION_IC 
            IF(.not.ic_defined(icv)) cycle 

            IF (IC_EP_G(ICV).lt.ONE) THEN 

               IF(DMP_LOG) WRITE(UNIT_LOG,2016) ICV, VOL_IC_REGION(ICV)
               IF(MYPE.EQ.PE_IO)  WRITE(*,2016) ICV, VOL_IC_REGION(ICV)
               DO M = 1, DES_MMAX
                  EP_SM = IC_ROP_S(ICV,M)/DES_RO_s(M)
                  IF(DMP_LOG)  WRITE(UNIT_LOG,2017) M,  &
                  DES_D_P0(M), EP_SM, PART_MPHASE_BYIC(ICV, M) 
                  IF(PRINT_DES_SCREEN)  WRITE(*,2017) M,  &
                  DES_D_P0(M), EP_SM, PART_MPHASE_BYIC(ICV, M) 
               ENDDO

            ENDIF

         ENDDO

         IF(DMP_LOG) WRITE(UNIT_LOG,2018)
         IF(PRINT_DES_SCREEN)  WRITE(*,2018)

            
      ENDIF 
! ----------------------------------------------------------------<<< 


! Check quantities related to MPPIC      
! ---------------------------------------------------------------->>>
      IF(MPPIC) THEN 

! For MPPIC, the periodicity is based on the continuum variables such 
! as cyclic_x, cyclic_x_pd, etc. 
         
         IF(CYCLIC_X.OR.CYCLIC_X_PD.OR.CYCLIC_Y.OR.CYCLIC_Y_PD.OR.&
            CYCLIC_Z.OR.CYCLIC_Z_PD ) THEN 
         
            IF(myPE.eq.pe_IO) WRITE(*,2005) 
            IF(DMP_LOG) WRITE(UNIT_LOG,2005) 
         
            DES_PERIODIC_WALLS = .TRUE. 
            DES_PERIODIC_WALLS_X = CYCLIC_X.OR.CYCLIC_X_PD
            DES_PERIODIC_WALLS_Y = CYCLIC_Y.OR.CYCLIC_Y_PD 
            DES_PERIODIC_WALLS_Z = CYCLIC_Z.OR.CYCLIC_Z_PD 
            
            IF(myPE.eq.pe_IO) WRITE(*,2006)  DES_PERIODIC_WALLS, &
            DES_PERIODIC_WALLS_X, DES_PERIODIC_WALLS_Y, DES_PERIODIC_WALLS_Z
            IF(DMP_LOG) WRITE(UNIT_LOG,2006)  DES_PERIODIC_WALLS, &
            DES_PERIODIC_WALLS_X, DES_PERIODIC_WALLS_Y, DES_PERIODIC_WALLS_Z     

         ENDIF
! check particle-wall normal restitution coefficient
         IF(MPPIC_COEFF_EN_WALL == UNDEFINED) THEN
            if(dmp_log) WRITE (UNIT_LOG, 2002) 'MPPIC_COEFF_EN_WALL'
            CALL MFIX_EXIT(myPE)
         ENDIF
         
         !IF(MPPIC_COEFF_ET_WALL == UNDEFINED) THEN
         !   if(dmp_log) WRITE (UNIT_LOG, 2002) 'MPPIC_COEFF_ET_WALL'
         !   CALL MFIX_EXIT(myPE)
         !ENDIF

         IF(MPPIC_COEFF_EN1 == UNDEFINED) THEN
            if(dmp_log) WRITE (UNIT_LOG, 2004) 'MPPIC_COEFF_EN1'
            CALL MFIX_EXIT(myPE)
         ENDIF

         IF(MPPIC_COEFF_EN2 == UNDEFINED) THEN
            if(dmp_log) WRITE (UNIT_LOG, 2004) 'MPPIC_COEFF_EN2'
            CALL MFIX_EXIT(myPE)
         ENDIF

         IF(MPPIC_COEFF_EN_WALL > ONE .OR. MPPIC_COEFF_EN_WALL < ZERO) THEN
            if(dmp_log) WRITE (UNIT_LOG, 2003) 'MPPIC_COEFF_EN_WALL'
            CALL MFIX_EXIT(myPE)
         ENDIF
         
         IF(MPPIC_COEFF_ET_WALL > ONE .OR. MPPIC_COEFF_ET_WALL < ZERO) THEN
            if(dmp_log) WRITE (UNIT_LOG, 2003) 'MPPIC_COEFF_ET_WALL'
            CALL MFIX_EXIT(myPE)
         ENDIF
        
         IF(MPPIC_COEFF_EN1 > ONE .OR. MPPIC_COEFF_EN1 < ZERO) THEN
            if(dmp_log) WRITE (UNIT_LOG, 2003) 'MPPIC_COEFF_EN1'
            CALL MFIX_EXIT(myPE)
         ENDIF

         IF(MPPIC_COEFF_EN2 > ONE .OR. MPPIC_COEFF_EN2 < ZERO) THEN
            if(dmp_log) WRITE (UNIT_LOG, 2003) 'MPPIC_COEFF_EN2'
            CALL MFIX_EXIT(myPE)
         ENDIF


         IF(MPPIC_SOLID_STRESS_SNIDER) WRITE(UNIT_LOG,*) &
         'USING THE SNIDER MODEL FOR SOLID STRESS AND THE ',&
         'INTEGRATION APPROACH'
         IF(MPPIC_SOLID_STRESS_SNIDER.AND.PRINT_DES_SCREEN) &
         WRITE(*,*) 'USING THE SNIDER MODEL FOR SOLID STRESS AND ',&
         'THE INTEGRATION APPROACH'
         

         IF(MPPIC_CONSTANTNPC.AND.MPPIC_CONSTANTWT) then 
            WRITE(UNIT_LOG, 2011)
            if(PRINT_DES_SCREEN) WRITE(*, 2011)
            call  mfix_exit(mype)
         ENDIF

         IF(.not.MPPIC_CONSTANTNPC.AND.(.not.MPPIC_CONSTANTWT)) then 
            WRITE(UNIT_LOG, 2012)
            if(PRINT_DES_SCREEN) WRITE(*, 2012)
            call  mfix_exit(mype)
         ENDIF

         IF(MPPIC_CONSTANTNPC) then 
            DO M = 1, DES_MMAX
               IF(NPC_PIC(M).Eq.UNDEFINED_I) then 
                  WRITE(UNIT_LOG, 2013) M
                  if(PRINT_DES_SCREEN) WRITE(*, 2013) M
                  call  mfix_exit(mype)
               ENDIF
            ENDDO
         ENDIF
         
         IF(MPPIC_CONSTANTWT) then 
            DO M = 1, DES_MMAX
               IF(STATWT_PIC(M).Eq.UNDEFINED) then 
                  WRITE(UNIT_LOG, 2014) M
                  if(PRINT_DES_SCREEN) WRITE(*, 2014) M
                  call  mfix_exit(mype)
               ENDIF
            ENDDO
         ENDIF

      
! cnp_array(ijk, 0) will contain the cumulative number of real 
! particles later in the handling of inflow BC for MPPIC. See 
! the mppic_mi_bc in mppic_wallbc_mod.f 
         ALLOCATE(CNP_ARRAY(DIMENSION_3, 0:DIMENSION_M))
         CNP_ARRAY(:, :) = 0

         IF(GENER_PART_CONFIG) THEN
            ALLOCATE(RNP_PIC(DES_MMAX))
            ALLOCATE(CNP_PIC(DES_MMAX))
            RNP_PIC = ZERO
            CNP_PIC = ZERO 
            IF(DIMN.EQ.2) THEN 
! require that DZ(1)/ZLENGTH be specified for 2-dimensional case.  
! unclear why this cannot be one - other than the user may be unaware 
! that a depth has been set (a value of one implies default setting) 
               IF (DZ(1) == ONE) THEN
                  IF(mype.eq.pe_IO) WRITE(*,'(5X,A,A,/5X,A,A)') &
                     'For DIMN = 2, specify a value for DZ(1) or ',&
                     'ZLENGTH in mfix.dat which is not',&
                     'equal to one. If you want it to be one then ',&
                     'set it close to one but not exactly one'

                  IF(DMP_LOG) WRITE(UNIT_LOG,'(5X,A,A,/5X,A,A)') &
                     'For DIMN = 2, specify a value for DZ(1) or ',&
                     'ZLENGTH in mfix.dat which is not',&
                     'equal to one. If you want it to be one then ',&
                     'set it close to one but not exactly one'
                  CALL mfix_exit(mype)
               ENDIF

               IF (DZ(1) .NE. ZLENGTH) THEN 
!this condition will probably never occur. Redundancy doesn't hurt, however!
                  
                  IF(myPE.eq.pe_IO) THEN
                     WRITE(*,'(5X,A,/5x,A)') &
                        'For DIMN = 2, DZ(1) and ZLENGTH are used ',&
                        'interchangeably', ' Specify same values for ',&
                        'DZ(1) and ZLENGTH'  
                     WRITE(*,'(5X,2(A20,2X,G17.8))') &
                        'DZ(1) = ', DZ(1), 'ZLENGTH = ', ZLENGTH
                  ENDIF
                  IF(DMP_LOG) THEN
                     WRITE(UNIT_LOG,'(5X,A,/5x,A)') &
                        'For DIMN = 2, DZ(1) and ZLENGTH are used ',&
                        'interchangeably', ' Specify same values for ',&
                        'DZ(1) and ZLENGTH'  
                     WRITE(UNIT_LOG,'(5X,2(A20,2X,G17.8))') &
                        'DZ(1) = ', DZ(1), 'ZLENGTH = ', ZLENGTH
                  ENDIF
               ENDIF
            ENDIF

            DO K = KSTART1, KEND1 
               DO J = JSTART1, JEND1
                  DO I = ISTART1, IEND1 
                     IJK  = FUNIJK(I,J,K)
                     IF(.NOT.FLUID_AT(IJK)) CYCLE 
                     IF(EP_G(IJK).GE.1.d0-DIL_EP_s) CYCLE 
                     VOLIJK = VOL(IJK)
                     VOLIJK_UNCUT = DX(I)*DY(J)*DZ(K) 
                     DO M = 1, DES_MMAX
                        REAL_PARTS(M) = 6.d0*EP_S(IJK,M)*VOLIJK/&
                           (PI*(D_p0(M)**3.d0))
                        IF(MPPIC_CONSTANTNPC) THEN 
                           COMP_PARTS(M) = NPC_PIC(M)
                           IF(CUT_CELL_AT(IJK)) COMP_PARTS(M) = &
                            INT(VOLIJK*real(COMP_PARTS(M))/VOLIJK_UNCUT)
                        ELSEIF(MPPIC_CONSTANTWT) THEN
                           COMP_PARTS(M) = MAX(1,INT(REAL_PARTS(M)/REAL(STATWT_PIC(M))))
                        ENDIF
                        
                        RNP_PIC(M) = RNP_PIC(M) + REAL_PARTS(M)
                        CNP_PIC(M) = CNP_PIC(M) + COMP_PARTS(M)
                        CNP_ARRAY(IJK,M) = COMP_PARTS(M)
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
            PART_MPHASE(1:DES_MMAX) = CNP_PIC(1:DES_MMAX)
            PARTICLES = SUM(PART_MPHASE(1:DES_MMAX))
            CALL global_all_sum(PARTICLES)
         ENDIF !  end if(gener_part_config)

      ENDIF   ! end if(mppic)
! ----------------------------------------------------------------<<<


      IF(MPPIC.OR.CARTESIAN_GRID) then
         IF(.NOT.DES_INTERP_MEAN_FIELDS) THEN
            DES_INTERP_MEAN_FIELDS = .TRUE.
            IF(myPe.EQ.pe_IO) WRITE(*,'(2(10X,A,/))') &
               'MPPIC OR CG CASE DETECTED:', 'Mean DES fields will be ',&
               'obtained by backward interpolation '
            IF(DMP_LOG) WRITE(UNIT_LOG,'(2(10X,A,/))') &
               'MPPIC OR CG CASE DETECTED:', 'Mean DES fields will be ',&
               'obtained by backward interpolation '
         ENDIF

      ELSE IF(DES_CONTINUUM_COUPLED.and.DES_INTERP_ON) then 
! in the previous versions, ep_s was obtained by interpolation in des_interp_on was 
! true. 
! For backward compatibility, des_interp_mean_fieldsis forced to true
! if des_interp_on is true 

         IF(.not.DES_INTERP_MEAN_FIELDS) then 
            IF(myPe.eq.pe_IO) WRITE(*,'(2(10X,A,/))') &
            'DES_INTERP_ON TRUE FOR DEM CASE:', &
            'Mean DES fields will be obtained by backward interpolation'
            IF(DMP_LOG) WRITE(UNIT_LOG,'(2(10X,A,/))') & 
            'DES_INTERP_ON TRUE  FOR DEM CASE:', & 
            'Mean DES fields will be obtained by backward interpolation'

            DES_INTERP_MEAN_FIELDS= .TRUE.
         ENDIF
         
      ELSE
         IF(DES_INTERP_MEAN_FIELDS) THEN 
            IF(myPe.EQ.pe_IO) WRITE(*,'(/,10X,A,L4,(/,10X,A))') & 
               'For DES_INTERP_MEAN_FIELDS = ', &
               DES_INTERP_MEAN_FIELDS, 'Mean DES fields will be ',&
               'obtained by backward interpolation'
            IF(DMP_LOG) WRITE(UNIT_LOG,'(/,10X, A, L4, (/,10X,A))') & 
               'For DES_INTERP_MEAN_FIELDS = ',& 
               DES_INTERP_MEAN_FIELDS, 'Mean DES fields will be ',&
               'obtained by backward interpolation'
         ELSE
            IF(myPe.eq.pe_IO) WRITE(*,'(/,10X, A, L4, (/,10X,A))') & 
               'For DES_INTERP_MEAN_FIELDS = ',&
               DES_INTERP_MEAN_FIELDS, 'Mean DES fields will be ',&
               'obtained by simple cell averaging'
            IF(DMP_LOG) WRITE(UNIT_LOG,'(/,10X, A, L4, (/,10X,A))') & 
               'For DES_INTERP_MEAN_FIELDS = ', &
               DES_INTERP_MEAN_FIELDS, 'Mean DES fields will be ',&
               'obtained by simple cell averaging'
         ENDIF
      ENDIF
            
! the entire checking and setting up indices for desgridsearch
! moved to desgrid_mod to accomodate parallelization
! this is now conducted regardless of neighbor search option      
      CALL DESGRID_CHECK()

      
      IF(DMP_LOG.AND.DEBUG_DES) WRITE(UNIT_LOG,'(1X,A)')&
         '<---------- END CHECK_DES_DATA ----------'

      RETURN

 1000 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'DES should only be run using CARTESIAN coordinates',&
         /1X,70('*')/)

 1001 FORMAT(/2X,'FULL_LOG = ',L3,' FOR DES ON IO PROC') 

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
         'Friction coefficients (MEW and MEW_w) must be ',&
         'specified in mfix.dat',/1X,70('*')/)

 1036 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         '2D coupled simulation with a particle diameter > ZLENGTH.',/10X,&
         'This will create problems for calculations of void ',&
         'fraction. Check',/10X, 'mfix.dat file.',/1X,70('*')/)

 1037 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
          'GENER_PAR_CONFIG set to false because a restart detected',&
          /1X,70('*')/)
  
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
         'MMAX.',/1X,70('*')/)
 1044 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'Too many D_P0 are defined in mfix.dat for given ',&
         'MMAX.',/1X,70('*')/)

 1045 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'DES_NEIGHBOR_SEARCH set to ', I2, ' ',A,/1X,70('*')/)

 1054 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
          'WARNING: zlength or dz(1) is used to calculate the ',&
          'number',/10X,'of particles in the 2D simulation when ',&
          'GENER_PART_CONFIG is T and DIMN = 2.',/10X,'This depth ',&
          'does not equal D_P0(1).',/1X,70('*'))

 1055 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
          'WARNING: zlength or dz(1) is used to calculate the ',&
          'number',/10X,'of particles in the 2D simulation when ',&
          'GENER_PART_CONFIG is T and DIMN = 2.',/1X,70('*'))

 1056     FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
          'Gener_part_config is true but using deprecated flags', /, &
          'DES_EPS_XSTART, DES_EPS_YSTART, and DES_EPS_ZSTART are obsolete flags.', /, &
          'IC region is now based on usual IC_flags', /, & 
          'Delete these flags from mfix.dat and restart', /, &
          'Exitting',/1X,70('*'))

 1057     FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
          'Gener_part_config is true but using deprecated flags', /, &
          'VOL_FRAC is an obsolete flag', /, &
          'IC region is now based on usual IC_flags', /, & 
          'Delete these flags from mfix.dat and restart', /, &
          'Exitting',/1X,70('*'))

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

 1071 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'RO_s must be defined and >0 in mfix.dat for M=1,',&
         'MMAX.',/1X,70('*')/)
 1072 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'Too many RO_s are defined in mfix.dat for given ',&
         'MMAX.',/1X,70('*')/)


 1080 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'USE_COHESION has not been verified with the current',/10X,&
         'DEM code so proceed at your own risk.',/1X,70('*')/)
 1081 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'Cannot use both SQUARE_WELL and VAN_DER_WAALS',/10X,&
         'cohesion models'/1X,70('*')/)            
 1082 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'A cohesion model has not been selected'/1X,70('*')/)
 1083 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'MASTER_WELL_DEPTH, MASTER_WALL_WELL_DEPTH',/10X,&
         'RADIUS_RATIO and/or WALL_RADIUS_RATIO are undefined ',&
         'but must be',/10X,'defined when SQUARE_WELL=T',/1X,70('*')/)
 1084 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'VDW_INNER_CUTOFF, VDW_OUTER_CUTOFF, HAMAKER_CONSTANT',/10X,&
         'are undefined but must be defined when ',&
         'VAN_DER_WAALS=T',/1X,70('*')/)
 1085 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'WALL_VDW_INNER_CUTOFF, WALL_HAMAKER_CONSTANT are ',&
         'undefined',/10X, 'but must be defined when ',&
         'VAN_DER_WAALS=T',/1X,70('*')/)

 1090 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'DES_MMAX not defined in mfix.dat. It must be defined',/10X,&
         'when using DES_CONTINUUM_HYBRID',/1X,70('*')/)

 1091 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'KT_TYPE cannot be set to GHD when using ',&
         'DES_CONTINUUM_HYBRID.',/1X,70('*')/)
 1092 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'QMOMK cannot be TRUE when using ',&
         'DES_CONTINUUM_HYBRID.',/1X,70('*')/)
 1093 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'MPPIC and DES_CONTINUUM_HYBRID cannot both be TRUE.',&
         /1X,70('*')/)
 1094 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'DES_CONTINUUM_COUPLED must be to true when using ',/10X,&
         'DES_CONTINUUM_HYBRID.',&
         /1X,70('*')/)
 1095 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'C_F must be defined when DES_CONTINUUM_HYBRID TRUE ',/10X,&
         'and both continuum and discrete solids phases are ',&
         'present'/10X,'(MAX>=1 and DES_MMAX>=1).',/1X,70('*')/)

 1100 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'WARNING: VDW_OUTER_CUTOFF outside of neighbor ',/10X,&
         'search distance. Increase FACTOR_RLM.',/1X,70('*')/) 
 1101 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'ERROR: CLUSTER_LENGTH_CUTOFF must be defined when ',/10X,&
         'using function DES_CALC_CLUSTER.',/1X,70('*')/)          
 1102 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'WARNING: CLUSTER_LENGTH_CUTOFF outside of neighbor ',/10X,&
         'search distance. Increase FACTOR_RLM.',/1X,70('*')/) 


 2001 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'Looks like a 3-D case (IMAX, JMAX & KMAX all >1) ',&
         'but DIMN equal to ',I1,/1X,70('*')/)

 2002    FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'Particle-wall restitution coefficient:  ',  A, &
         /10X,'not specified in mfix.dat for the MPPIC model',&
         /1X,70('*')/)

 2003    FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'Unphysical ( > 1 or < 0) values of:  ', A, &
         /1X,70('*')/)

 2004    FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'Friction restitution coefficient:  ',  A, &
         /10X,'not specified in mfix.dat for the MPPIC model',&
         /1X,70('*')/)
 2005    FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ', &
         'MPPIC Case: One of the cyclic_flags found to be true', /10X, & 
         'from  continuum flags, so forcing des_periodic_walls to true')

 2006    FORMAT(/10X,'Updated Periodicity related flags:', /10X,&
         'DES_PERIODIC_WALLS   :', L4, /10X, &
         'DES_PERIODIC_WALLS_X :', L4, /10X, &
         'DES_PERIODIC_WALLS_Y :', L4, /10X, &
         'DES_PERIODIC_WALLS_Z :', L4, /10X, &
         /1X,70('*')/)

 2010 FORMAT(/1X,70('*')/' From: CHECK_DES_DATA',/'Error 2010:',       &
         ' DES simulations cannot have defined internal obstacles.',/  &
         ' Please correct the data file.', 1X,70('*')/)
         
 2011    FORMAT(/1X,70('*')/' From: CHECK_DES_DATA',/' Error 2011:',       &
         ' In MPPIC model, both MPPIC_CONSTANTNPC and MPPIC_CONSTANTWT', /,  & 
         ' set to TRUE. Set at least one to false. See MFIX readme', /,  &
         ' Please correct the data file. Exiting.', /,1X,70('*')/)

 2012    FORMAT(/1X,70('*')/' From: CHECK_DES_DATA',/' Error 2012:',       &
         ' In MPPIC model, both MPPIC_CONSTANTNPC and MPPIC_CONSTANTWT', /,  & 
         ' set to false. Set at least one to true. See MFIX readme', /,  &
         ' Please correct the data file. Exiting.', /,1X,70('*')/)

 2013    FORMAT(/1X,70('*')/' From: CHECK_DES_DATA',/' Error 2013:',       &
         ' In MPPIC model, MPPIC_CONSTANTNPC specified', &
         ' but NPC_PIC undefined', /,&
         ' for phase', I4, /, &
         ' See MFIX readme', /,  &
         ' Please correct the data file. Exiting.',/, 1X,70('*')/)

 2014    FORMAT(/1X,70('*')/' From: CHECK_DES_DATA',/' Error 2014:',       &
         ' In MPPIC model, MPPIC_CONSTANTWT specified', &
         ' but STATWT_PIC undefined', /,&
         ' for phase', I4, /, &
         ' See MFIX readme.', /,  &
         ' Please correct the data file. Exiting.',/, 1X,70('*')/)

 2015    FORMAT( 1X,70('*')/ & 
         'From: CHECK_DES_DATA    ', /5x, &
         'gener_part_config specified as true', /,5x, &
         'Total IC region count with non zero solids = ', I5,2X, /5x, & 
         'Total discrete solid phases = ', I5,2X, /5x, & 
         'Particle configuration will be automatically generated: ')

 2016    FORMAT(/1X,70('-')/, 5x, &
         'IC number and Volume        = ', I4, 2x, (ES15.7,2X) )
         
 2017    FORMAT(1X,70('.'),/,5x, &
         'PHASE Index, M              =  ', I5,2X, /5x, & 
         'D_P0(M)                     =  ', (ES15.7,2X), /5x, &
         'Solids vol fraction for M   =  ', (G15.8,2X), /5x, & 
         '# of particles for phase M  = ',  (I10,2X))

 2018    FORMAT( 1X,70('*')/)

         END SUBROUTINE CHECK_DES_DATA
