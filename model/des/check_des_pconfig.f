!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_DES_PCONFIG                                       !
!  Author: J.Musser                                   Date: 11-Dec-13  !
!                                                                      !
!  Purpose: Check user input data for generating initial particle      !
!  configurations.                                                     !
!                                                                      !
!  Comments: Original code moded here from CHECK_DES_DATA.             !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_DES_PCONFIG


! Global Variables:
!---------------------------------------------------------------------//
! Runtime Flag: Invoke MPPIC model.
      USE mfix_pic, only: MPPIC
! Runtime Flag: Generate initial particle configuation.
      USE discretelement, only: GENER_PART_CONFIG
! Runtime Flag: Printing screen messages. PE_IO Only.
      USE discretelement, only: PRINT_DES_SCREEN

! Simulation dimension (2D/3D)
      USE discretelement, only: DIMN
! Number of DEM solids phases.
      USE discretelement, only: DES_MMAX
! DEM solid phase diameters and densities.
      USE discretelement, only: DES_D_p0, DES_RO_s
! Computed volume of IC region for seeding
      USE discretelement, only: VOL_IC_REGION
! Number of particles seeded, per phase in each IC region
      USE discretelement, only: PART_MPHASE_BYIC
! Number of particles to read from input file.
      USE discretelement, only: PARTICLES
! Type of run (New, Restart)
      USE run, only: RUN_TYPE
! Constant: 3.14159...
      USE constant, only: PI
! Flag indicating that the IC region is defined.
      USE ic, only: IC_DEFINED
! IC Region bulk density (RO_s * EP_s)
      USE ic, only: IC_ROP_s
! IC Region gas volume fraction.
      USE ic, only: IC_EP_G
! Parameter for detecting unspecified values.
      USE param1, only: UNDEFINED
! File unit for LOG messages.
      USE funits, only: UNIT_LOG

! Deprecated varaibles for specifying particle generation.
      USE discretelement, only: VOL_FRAC
      USE discretelement, only: DES_EPS_XSTART
      USE discretelement, only: DES_EPS_YSTART
      USE discretelement, only: DES_EPS_ZSTART

      USE mpi_utility

      IMPLICIT NONE

! Local Variables:
!---------------------------------------------------------------------//
! Loop index.
      INTEGER :: M   ! Phase
      INTEGER :: ICV ! IC Region
! Solids phase volume fraction.
      DOUBLE PRECISION :: EP_SM
! Temp variable for storing number of particles.
      DOUBLE PRECISION :: PARTS_TEMP 
! Number of IC regions with solids.
      INTEGER :: IC_COUNT

! External Functions:
!---------------------------------------------------------------------//
      LOGICAL, EXTERNAL :: COMPARE

! If run_type is not 'NEW' then force the gener_part_config to .false. 
! This will prevent overwriting over the number of particles which could
! have potentially changed depending on inlet/outlet      
      IF(TRIM(RUN_TYPE) .NE. 'NEW'.AND.GENER_PART_CONFIG) THEN 
         GENER_PART_CONFIG = .FALSE. 
         IF(DMP_LOG) WRITE(UNIT_LOG, 1037)
         IF(PRINT_DES_SCREEN) WRITE(*, 1037)
      ENDIF


! Check if any of the deprecated flags are in the input file.
      if(des_eps_xstart.ne.undefined .or. &
         des_eps_ystart.ne.undefined .or. &
         des_eps_xstart.ne.undefined) then 
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

! Only new DEM runs proceed with the following checks.
      IF(TRIM(RUN_TYPE) .NE. 'NEW') RETURN
      IF(MPPIC) RETURN

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

      IF(DMP_LOG) WRITE(UNIT_LOG,2015) IC_COUNT, DES_MMAX
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


      RETURN


 1037 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
          'GENER_PAR_CONFIG set to false because a restart detected',&
          /1X,70('*')/)

 1054 FORMAT(/1X,70('*')//' From: CHECK_DES_PCONFIG',/' Message: ',&
          'WARNING: zlength or dz(1) is used to calculate the ',&
          'number',/10X,'of particles in the 2D simulation when ',&
          'GENER_PART_CONFIG is T and DIMN = 2.',/10X,'This depth ',&
          'does not equal D_P0(1).',/1X,70('*'))

 1055 FORMAT(/1X,70('*')//' From: CHECK_DES_PCONFIG',/' Message: ',&
          'WARNING: zlength or dz(1) is used to calculate the ',&
          'number',/10X,'of particles in the 2D simulation when ',&
          'GENER_PART_CONFIG is T and DIMN = 2.',/1X,70('*'))

 1056 FORMAT(/1X,70('*')//' From: CHECK_DES_PCONFIG',/' Message: ',&
          'Gener_part_config is true but using deprecated flags', /, &
          'DES_EPS_XSTART, DES_EPS_YSTART, and DES_EPS_ZSTART are obsolete flags.', /, &
          'IC region is now based on usual IC_flags', /, & 
          'Delete these flags from mfix.dat and restart', /, &
          'Exitting',/1X,70('*'))

 1057 FORMAT(/1X,70('*')//' From: CHECK_DES_PCONFIG',/' Message: ',&
          'Gener_part_config is true but using deprecated flags', /, &
          'VOL_FRAC is an obsolete flag', /, &
          'IC region is now based on usual IC_flags', /, & 
          'Delete these flags from mfix.dat and restart', /, &
          'Exitting',/1X,70('*'))


 2015    FORMAT( 1X,70('*')/ & 
         'From: CHECK_DES_PCONFIG    ', /5x, &
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

         END SUBROUTINE CHECK_DES_PCONFIG
