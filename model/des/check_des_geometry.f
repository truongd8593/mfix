!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_DES_GEOMETRY                                      !
!  Purpose: Check user input data                                      !
!                                                                      !
!  Author: J.Musser                                   Date: 11-DEC-13  !
!                                                                      !
!  Comments: Geometry checks were moved here from CHECK_DES_DATA.      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_DES_GEOMETRY

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


! DEM is restriced to CARTESIAN coordinates.
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


! Check dimension. This is redundant with check_data_03.
      IF(NO_I.OR.NO_J) THEN
         IF(DMP_LOG) WRITE(UNIT_LOG,1039)
         CALL MFIX_EXIT(myPE)
      ENDIF

! Verify that DIMN was specified.
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


! des_periodic_walls must be true if any wall is periodic      
      IF ( (DES_PERIODIC_WALLS_X .OR. DES_PERIODIC_WALLS_Z .OR. &
            DES_PERIODIC_WALLS_Z) .AND. .NOT.DES_PERIODIC_WALLS) THEN
         DES_PERIODIC_WALLS = .TRUE.
         IF(DMP_LOG) WRITE(UNIT_LOG, 1017)
      ENDIF


! Periodicity does not work with all neighbor search options
      IF(DES_PERIODIC_WALLS) THEN
         IF(.NOT.DES_PERIODIC_WALLS_X .AND. &
            .NOT.DES_PERIODIC_WALLS_Y .AND. &
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



      RETURN

 1000 FORMAT(/1X,70('*')//' From: CHECK_DES_GEOMETRY',/' Message: ',&
         'DES should only be run using CARTESIAN coordinates',&
         /1X,70('*')/)

 1002 FORMAT(/1X,70('*')//' From: CHECK_DES_GEOMETRY',/' Message: ',&
         'Direction of periodicity not defined in mfix.dat',&
         /1X,70('*')/)

 1017 FORMAT(/1X,70('*')//' From: CHECK_DES_GEOMETRY',/' Message: ',&
         'WARNING: A direction of periodicity was defined (i.e. ',/10X,&
         'DES_PERIODIC_WALLS_X, _Y or _Z=T) but DES_PERIODIC_WALLS ',&
         'was set to F.',/10X,&
         'The latter was forced to T for consistency.',/1X,70('*')/)

 1018 FORMAT(/1X,70('*')//' From: CHECK_DES_GEOMETRY',/' Message: ',&
         'WARNING: nsquare neighbor search may be slow with periodic',/10X,&
         'boundaries.  Grid based search is recommended.',/1X,70('*')/)

 1019 FORMAT(/1X,70('*')//' From: CHECK_DES_GEOMETRY',/' Message: ',&
         'octree and quadtree neighbor search methods do not',/10X,&
         'currently work with periodic boundaries. Change ',&
         'value of',/10X,'DES_NEIGHBOR_SEARCH in mfix.dat',/1X,70('*')/)

 1036 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         '2D coupled simulation with a particle diameter > ZLENGTH.',/10X,&
         'This will create problems for calculations of void ',&
         'fraction. Check',/10X, 'mfix.dat file.',/1X,70('*')/)

 1039 FORMAT(/1X,70('*')//' From: CHECK_DES_GEOMETRY',/' Message: ',&
         'DES can only be run in XY plane in 2D.',/1X,70('*')/)

 1040 FORMAT(/1X,70('*')//' From: CHECK_DES_GEOMETRY',/' Message: ',&
         'WARNING: DIMN left unspecified.  Since NO_K is T ',/10X,&
         'DIMN was set to 2. 2D DES running in XY plane.',/1X,70('*')/)

 1041 FORMAT(/1X,70('*')//' From: CHECK_DES_GEOMETRY',/' Message: ',&
         'Geometry dimension DIMN not specified.',/1X,70('*')/)

 1042 FORMAT(/1X,70('*')//' From: CHECK_DES_GEOMETRY',/' Message: ',&
         'Physical dimension DIMN cannot be > 3.',/1X,70('*')/)
         
 2001 FORMAT(/1X,70('*')//' From: CHECK_DES_GEOMETRY',/' Message: ',&
         'Looks like a 3-D case (IMAX, JMAX & KMAX all >1) ',&
         'but DIMN equal to ',I1,/1X,70('*')/)

 2010 FORMAT(/1X,70('*')/' From: CHECK_DES_GEOMETRY',/'Error 2010:',       &
         ' DES simulations cannot have defined internal obstacles.',/  &
         ' Please correct the data file.', 1X,70('*')/)
         
         END SUBROUTINE CHECK_DES_GEOMETRY
