!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CHECK_DES_DATA                                         C
!  Purpose: DES - Check user input Data   
!                                                                      C
!                                                                      C
!  Author: S. Benyahia                                Date: 10-Nov-08  C
!  Reviewer:                                                           C
!  Comments: All user's input data are checked here                    C
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

      IMPLICIT NONE
!-----------------------------------------------
! Local Variables
!-----------------------------------------------      
      INTEGER CHECK_MPI, M
      LOGICAL FLAG_WARN
      INTEGER I, J, K
! quantities used to check/specify mesh size if des_neighbor_search = 4      
      DOUBLE PRECISION DL_TMP, TMP_FACTOR     
! needed if des_neighbor_search=4, the maximum specified particle diameter 
      DOUBLE PRECISION MAX_DIAM      
!----------------------------------------------- 

      WRITE(*,'(1X,A)')&
         '---------- START CHECK_DES_DATA ---------->'


      IF(COORDINATES == 'CYLINDRICAL') THEN
         WRITE (UNIT_LOG, 1000)
         CALL MFIX_EXIT(myPE)
      ENDIF

      CHECK_MPI = NODESI * NODESJ * NODESK
      IF(CHECK_MPI.NE.1) THEN
         WRITE (UNIT_LOG, 1001)
         CALL MFIX_EXIT(myPE)
      ENDIF

! if you want to run DEM with the existing cohesion implementation then
! comment out the following 4 lines; also note that this feature is
! unlikely to function correctly with the current version of DEM so
! proceed at your own risk
      IF(USE_COHESION) THEN
         WRITE(UNIT_LOG, 1035)
         CALL MFIX_EXIT(myPE)
      ENDIF

      IF ( (DES_PERIODIC_WALLS_X .OR. DES_PERIODIC_WALLS_Z .OR. &
            DES_PERIODIC_WALLS_Z) .AND. .NOT.DES_PERIODIC_WALLS) THEN
         DES_PERIODIC_WALLS = .TRUE.
         WRITE(UNIT_LOG, 1017)
      ENDIF

      IF(DES_NEIGHBOR_SEARCH.EQ.1) THEN
         WRITE(*,'(3X,A)') &
            'DES_NEIGHBOR_SEARCH set to 1 (N-SQUARE)'
      ELSEIF(DES_NEIGHBOR_SEARCH.EQ.2) THEN
         WRITE(*,'(3X,A)') &
            'DES_NEIGHBOR_SEARCH set to 2 (QUADTREE)'
      ELSEIF(DES_NEIGHBOR_SEARCH.EQ.3) THEN
         WRITE(*,'(3X,A)') &
            'DES_NEIGHBOR_SEARCH set to 3 (OCTREE)'
      ELSEIF(DES_NEIGHBOR_SEARCH.EQ.4) THEN
         WRITE(*,'(3X,A)') &
            'DES_NEIGHBOR_SEARCH set to 4 (GRID BASED)'
      ELSE
         WRITE (UNIT_LOG, 1003)
         WRITE (*, 1003)
         CALL MFIX_EXIT(myPE)         
      ENDIF

      IF(DES_PERIODIC_WALLS) THEN
         IF(.NOT.DES_PERIODIC_WALLS_X .AND. .NOT.DES_PERIODIC_WALLS_Y .AND. &
            .NOT.DES_PERIODIC_WALLS_Z) THEN
            WRITE (UNIT_LOG, 1002)
            CALL MFIX_EXIT(myPE)
         ENDIF
         IF (DES_NEIGHBOR_SEARCH .EQ. 1) THEN
            WRITE (UNIT_LOG, 1018)
            WRITE (*,1018)
         ENDIF
         IF (DES_NEIGHBOR_SEARCH .EQ. 2 .OR. &
          DES_NEIGHBOR_SEARCH .EQ. 3) THEN
            WRITE (UNIT_LOG, 1019)
            CALL MFIX_EXIT(myPE)
         ENDIF
      ENDIF


      IF (TRIM(DES_INTG_METHOD) /= 'ADAMS_BASHFORTH' .AND. &
          TRIM(DES_INTG_METHOD) /= 'EULER') THEN
! stop if the specified integration method is unavailable
         WRITE (UNIT_LOG, 1034)
         CALL MFIX_EXIT(myPE)
      ENDIF


      IF (DES_COLL_MODEL /= UNDEFINED_C) THEN      
         IF (TRIM(DES_COLL_MODEL) .NE. 'HERTZIAN') THEN
            WRITE(UNIT_LOG, 1006)
            CALL MFIX_EXIT(myPE)
         ENDIF
      ENDIF

      IF (TRIM(DES_COLL_MODEL) .EQ. 'HERTZIAN') THEN
! check young's modulus and poisson ratio
         IF(EW_YOUNG == UNDEFINED ) THEN
            WRITE (UNIT_LOG, 1020)
            CALL MFIX_EXIT(myPE)
         ENDIF
         IF(VW_POISSON == UNDEFINED) THEN
            WRITE (UNIT_LOG, 1021)
            CALL MFIX_EXIT(myPE)
         ELSE
           IF (VW_POISSON > 0.5d0 .OR. VW_POISSON <= -ONE) THEN
              WRITE (UNIT_LOG, 1033)
           ENDIF
         ENDIF
         DO M = 1, MMAX
            IF(E_YOUNG(M) == UNDEFINED) THEN
               WRITE (UNIT_LOG, 1022)
               CALL MFIX_EXIT(myPE)
            ENDIF
            IF(V_POISSON(M) == UNDEFINED) THEN
               WRITE (UNIT_LOG, 1023)
               CALL MFIX_EXIT(myPE)
            ELSE
              IF (V_POISSON(M) > 0.5d0 .OR. V_POISSON(M) <= -ONE) THEN
                 WRITE (UNIT_LOG, 1033)
              ENDIF
            ENDIF
         ENDDO
! check particle-particle tangential restitution coefficient      
         DO M = 1, MMAX+MMAX*(MMAX-1)/2
            IF(DES_ET_INPUT(M) == UNDEFINED) THEN
               WRITE (UNIT_LOG, 1024) MMAX+MMAX*(MMAX-1)/2
               CALL MFIX_EXIT(myPE)
            ENDIF
         ENDDO
         DO M = 1, MMAX+MMAX*(MMAX-1)/2
            IF(DES_ET_INPUT(M) > ONE .OR. DES_ET_INPUT(M) < ZERO) THEN
               WRITE (UNIT_LOG, 1025)
               CALL MFIX_EXIT(myPE)
            ENDIF
         ENDDO
! check particle-wall tangential restitution coefficient
         DO M = 1, MMAX
            IF(DES_ET_WALL_INPUT(M) == UNDEFINED) THEN
               WRITE (UNIT_LOG, 1026)
               CALL MFIX_EXIT(myPE)
            ENDIF
         ENDDO
         DO M = 1, MMAX
            IF(DES_ET_WALL_INPUT(M) > ONE .OR. DES_ET_WALL_INPUT(M) < ZERO) THEN
               WRITE (UNIT_LOG, 1027)
               CALL MFIX_EXIT(myPE)
            ENDIF
         ENDDO
! if following are assigned warn user they are discarded
         IF(KN .NE. UNDEFINED .OR. KN_W .NE. UNDEFINED) &
            WRITE (UNIT_LOG, 1028)
         IF(KT_FAC .NE. UNDEFINED .OR. KT_W_FAC .NE. UNDEFINED) &
            WRITE (UNIT_LOG, 1029)
         IF(DES_ETAT_FAC .NE. UNDEFINED .OR. &
            DES_ETAT_W_FAC .NE. UNDEFINED) WRITE (UNIT_LOG, 1030)

      ELSE  ! default linear spring-dashpot model

! check normal spring constants 
         IF(KN == UNDEFINED .OR. KN_W == UNDEFINED) THEN
            WRITE (UNIT_LOG, 1004)
            CALL MFIX_EXIT(myPE)
         ENDIF
! check for tangential spring constant factors         
         IF(KT_FAC == UNDEFINED .OR. KT_W_FAC == UNDEFINED) THEN
            WRITE (UNIT_LOG, 1005)
         ENDIF
         IF(KT_FAC .NE. UNDEFINED) THEN
            IF(KT_FAC > ONE .OR. KT_FAC < ZERO) THEN
               WRITE (UNIT_LOG, 1016)
               CALL MFIX_EXIT(myPE)
            ENDIF
         ENDIF
         IF(KT_W_FAC .NE. UNDEFINED) THEN
            IF(KT_W_FAC > ONE .OR. KT_W_FAC < ZERO) THEN
               WRITE (UNIT_LOG, 1016)
               CALL MFIX_EXIT(myPE)
            ENDIF
         ENDIF
! check for tangential damping factor
         IF(DES_ETAT_FAC == UNDEFINED) THEN
            WRITE (UNIT_LOG, 1010)
         ELSEIF(DES_ETAT_FAC > ONE .OR. DES_ETAT_FAC < ZERO) THEN
            WRITE (UNIT_LOG, 1015)
            CALL MFIX_EXIT(myPE)
         ENDIF
         IF(DES_ETAT_W_FAC == UNDEFINED) THEN
            WRITE (UNIT_LOG, 1011)
         ELSEIF(DES_ETAT_W_FAC > ONE .OR. DES_ETAT_W_FAC < ZERO) THEN
            WRITE (UNIT_LOG, 1015)
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
            WRITE (UNIT_LOG, 1008) MMAX+MMAX*(MMAX-1)/2
            CALL MFIX_EXIT(myPE)
         END IF
      ENDDO
      DO M = 1, MMAX+MMAX*(MMAX-1)/2
         IF(DES_EN_INPUT(M) > ONE .OR. DES_EN_INPUT(M) < ZERO) THEN
            WRITE (UNIT_LOG, 1012)
            CALL MFIX_EXIT(myPE)
         ENDIF
      ENDDO

! check particle-wall normal restitution coefficient
      DO M = 1, MMAX
         IF(DES_EN_WALL_INPUT(M) == UNDEFINED) THEN
            WRITE (UNIT_LOG, 1009)
            CALL MFIX_EXIT(myPE)
         ENDIF
      ENDDO
      DO M = 1, MMAX
         IF(DES_EN_WALL_INPUT(M) > ONE .OR. DES_EN_WALL_INPUT(M) < ZERO) THEN
            WRITE (UNIT_LOG, 1013)
            CALL MFIX_EXIT(myPE)
         ENDIF
      ENDDO

! check coefficient friction 
      IF(MEW > ONE .OR. MEW_W > ONE .OR. MEW < ZERO .OR. MEW_W < ZERO) THEN
         WRITE (UNIT_LOG, 1014)
         CALL MFIX_EXIT(myPE)
      END IF

! Overwrite user's input in case of DEM (no fluid)
      IF(.NOT.DES_CONTINUUM_COUPLED) DES_INTERP_ON = .FALSE.

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

! Ensure settings for grid based neighbor search method. This section
! was placed here since the domain lengths (xlength, ylength, zlength)
! are needed and their value is not gauranteed until after 
! check_data_03 is called
! ------------------------------------------------------------
      IF (DES_NEIGHBOR_SEARCH .EQ. 4) THEN
         MAX_DIAM = 2.0d0*MAX_RADIUS

         TMP_FACTOR = 3.0d0*(MAX_DIAM)

! If the search grid is undefined then set it to approximately 3 times
! the maximum particle. Otherwise, check to see that the user set search
! grid is at least greater than or equal to the maximum particle
! diameter and warn the user if not         
         IF (DESGRIDSEARCH_IMAX == UNDEFINED_I) THEN
            DL_TMP = XLENGTH/TMP_FACTOR
            DESGRIDSEARCH_IMAX = INT(DL_TMP)
            IF (DESGRIDSEARCH_IMAX <= 0) DESGRIDSEARCH_IMAX = 1
            WRITE(*,'(3X,A,I8)') &
               'DESGRIDSEARCH_IMAX was set to ', DESGRIDSEARCH_IMAX
         ELSE
            DL_TMP = XLENGTH/DBLE(DESGRIDSEARCH_IMAX)
            IF (DL_TMP < MAX_DIAM) THEN
               WRITE(*,1037) 'x', 'x', 'i', 'i'                    
               CALL MFIX_EXIT(myPE)
            ENDIF
         ENDIF
         IF (DESGRIDSEARCH_JMAX == UNDEFINED_I) THEN
            DL_TMP = YLENGTH/TMP_FACTOR
            DESGRIDSEARCH_JMAX = INT(DL_TMP)
            IF (DESGRIDSEARCH_JMAX <= 0) DESGRIDSEARCH_JMAX = 1
            WRITE(*,'(3X,A,I8)') &
               'DESGRIDSEARCH_JMAX was set to ', DESGRIDSEARCH_JMAX
         ELSE
            DL_TMP = YLENGTH/DBLE(DESGRIDSEARCH_JMAX)
            IF (DL_TMP < MAX_DIAM) THEN
               WRITE(*,1037) 'y', 'y', 'j', 'j'
               CALL MFIX_EXIT(myPE)
            ENDIF
         ENDIF
         IF (DIMN .EQ. 2) THEN
            IF (DESGRIDSEARCH_KMAX == UNDEFINED_I) THEN
               DESGRIDSEARCH_KMAX = 1
            ELSEIF(DESGRIDSEARCH_KMAX /= 1) THEN
               DESGRIDSEARCH_KMAX = 1            
               WRITE(*,'(3X,A,I8)') &
                  'DESGRIDSEARCH_KMAX was set to ', DESGRIDSEARCH_KMAX
            ENDIF            
         ELSE
            IF (DESGRIDSEARCH_KMAX == UNDEFINED_I) THEN
                DL_TMP = ZLENGTH/TMP_FACTOR
                DESGRIDSEARCH_KMAX = INT(DL_TMP)
                IF (DESGRIDSEARCH_KMAX <= 0) DESGRIDSEARCH_JMAX = 1
                WRITE(*,'(3X,A,I8)') &
                  'DESGRIDSEARCH_KMAX was set to ', DESGRIDSEARCH_KMAX
            ELSE
                DL_TMP = ZLENGTH/DBLE(DESGRIDSEARCH_KMAX)
                IF (DL_TMP < MAX_DIAM) THEN
                   WRITE(*,1037) 'z', 'z', 'k', 'k'
                   CALL MFIX_EXIT(myPE)
               ENDIF
            ENDIF
         ENDIF   ! end if/else dimn == 2

         DESGS_IMAX2= DESGRIDSEARCH_IMAX+2
         DESGS_JMAX2 = DESGRIDSEARCH_JMAX+2
         IF (DIMN .EQ. 2) THEN
            DESGS_KMAX2 = DESGRIDSEARCH_KMAX
         ELSE
            DESGS_KMAX2 = DESGRIDSEARCH_KMAX+2
         ENDIF         

! Variable that stores the particle in cell information (ID) on the
! computational grid defined by cell/grid based search.  Similar to the
! variable PIC but tailored for the grid based neighbor search option
         ALLOCATE(DESGRIDSEARCH_PIC(DESGS_IMAX2,DESGS_JMAX2,DESGS_KMAX2))
         DO K = 1,DESGS_KMAX2
            DO J = 1,DESGS_JMAX2
               DO I = 1,DESGS_IMAX2
                  NULLIFY(DESGRIDSEARCH_PIC(I,J,K)%p) 
               ENDDO 
             ENDDO 
          ENDDO                      

      ENDIF   ! end if des_neighbor_search == 4
! End checks if grid based neighbor search         
! ------------------------------------------------------------

      
      WRITE(*,'(1X,A)')&
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
          'length/dessearchgrid_',A,'max. This will',/10X,'create ',&
          'problems for the search method and detecting neighbors',/10X,&
          'Decrease desgridsearch_',A,'max in mfix.dat to coarsen ',&
          'grid.',/1X,70('*')/)

         END SUBROUTINE CHECK_DES_DATA
