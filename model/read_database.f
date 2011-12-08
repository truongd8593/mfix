!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: read_database(Ier)                                     C
!  Purpose: read thermochemical database                               C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-OCT-05  C
!                                                                      C
!  Modification 1: J. Musser                          Date: 02-May-11  C
!  Purpose: Provided support for DEM access to database.               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE READ_DATABASE(Ier) 

      USE param 
      USE param1 
      USE physprop
      USE constant
      USE compar
      USE rxns
      USE funits 
      USE discretelement
      USE des_rxns

      IMPLICIT NONE

! External function
      DOUBLE PRECISION, EXTERNAL ::calc_ICpoR
! Molecular weight obtained from database. Only employed if the user
! has not specified a molecular weight in the mfix.dat file.
      DOUBLE PRECISION :: MW

! Error message returned from Read_Therm and sent to calling routine
      INTEGER IER
! Loop indicies for mass phase and species
      INTEGER M, N
! Loop counter for continuum and discrete species
      INTEGER Nsp, DES_Nsp
! File unit of databases searched
      INTEGER FUNIT
! Loop counter for checking all three locations for data
      INTEGER FILE
! Input/Output error status ( 0 is no error)
      INTEGER IOS

! Identifies if an error was detected.
      LOGICAL ErrorFlag, NOT_WARNED

! File name of Burcat and Ruscic database
      CHARACTER(len=10) :: THERM = 'BURCAT.THR'
! Full path to model directory
      INCLUDE 'mfix_directory_path.inc'
! Full path to Burcat and Ruscic database
      CHARACTER(len=147) FILENAME

! Identifies if a particular continuum/ discrete species has been found
! in one of the databases.
      LOGICAL SPECIES_READ(DIM_N_ALL)
      LOGICAL DES_SPECIES_READ(DIM_N_ALL)

! Return to the calling routine if the database has already been called.
      IF(database_read)RETURN

! Set the flag identifying that the database has been read.
      database_read = .TRUE.

! Initialize logicals identifying that the data for a species has been
! read from the database.
      SPECIES_READ(:) = .FALSE.
      DES_SPECIES_READ(:) = .FALSE.

! Initialize the file unit to be used.
      FUNIT = UNIT_DAT  ! .dat file unit
! Read data from mfix.dat or from BURCAT.THR in run directory or
! mfix directory.
      FILE = 0
      DO
         FILE = FILE + 1
! Check for thermochemical data in the mfix.dat file.
         IF(FILE == 1) then
           OPEN(UNIT=FUNIT, FILE='mfix.dat', STATUS='OLD', IOSTAT= IOS)
	          IF(IOS /= 0) CYCLE       ! Cycle on opening error
           WRITE(*,1000)'mfix.dat'  ! Identify read file. (screen)
           IF(DMP_LOG) WRITE(UNIT_LOG,1000)'mfix.dat' ! (log file)
! Read thermochemical data from the BURCAT.THR database in the local
! run directory.
	        ELSEIF(FILE == 2) THEN
            OPEN(UNIT=FUNIT,FILE=TRIM(THERM), STATUS='OLD', IOSTAT= IOS)
	           IF(IOS /= 0) CYCLE         ! Cycle on opening error
            WRITE(*,1000) TRIM(THERM)  ! Identify read file.
            IF(DMP_LOG) WRITE(UNIT_LOG,1000) TRIM(THERM)  ! (log file)
! Read thermochemical data from the BURCAT.THR database in the model
! directory (model/thermochemical/BURCAT.THR).
     	   ELSEIF(file == 3) then
            FILENAME = trim(MFIX_PATH)//'/thermochemical/'//TRIM(THERM)
            OPEN(UNIT=FUNIT,FILE=TRIM(FILENAME), STATUS='OLD',IOSTAT= IOS)
	           IF(IOS /= 0) CYCLE            ! Cycle on opening error
            WRITE(*,1000) TRIM(FILENAME)  ! Identify read file.
            IF(DMP_LOG) WRITE(UNIT_LOG,1000) TRIM(FILENAME)  ! (log file)
        	ELSE
            EXIT 
        	ENDIF
! Initialize counters
         Nsp = 0
         DES_Nsp = 0

! Read species data for the gas phase.
!-----------------------------------------------------------------------
         DO N = 1, NMAX(0)
        	   Nsp = Nsp + 1
! Skip species that have already been read from a previous source.
	           IF(SPECIES_READ(Nsp)) CYCLE
	           REWIND(UNIT=funit)
! If a species name was not specified in mfix.dat, flag error and exit.
	           IF(SPECIES_NAME(Nsp) == UNDEFINED_C) THEN
               WRITE(*,1010) N         ! screen
               IF(DMP_LOG) WRITE(UNIT_LOG,1010) N  ! log file
	              CALL MFIX_EXIT(mypE)
	           ENDIF
! Read the database.
        	   CALL Read_Therm(FUNIT, SPECIES_NAME(Nsp), Thigh_g(N),      &
               Tlow_g(N), Tcom_g(N), MW, Ahigh_g(1,N), Alow_g(1,N),    &
               HfrefoR_g(N), IER)
         	  IF(IER == 0)THEN
! If the user did not supply a value for the gas phase molecular weight
! in the mfix.dat file, use the value from the database.
	              IF(MW_g(N) == UNDEFINED) MW_g(N) = MW
! There are a number of species with Tlow as 300, for which the
! following calculation will produce an error because T_ref = 298.  So
! slightly extend validity of the correaltion.
               IF(ABS(Tlow_g(N)-T_ref)<=2.0D0 .AND. Tlow_g(N)>T_ref)   &
                  Tlow_g(N) = T_ref
! Calculate IC data.
	              IC_PGrefoR(N) = calc_ICpoR(T_ref, Thigh_g(N),           &
                  Tlow_g(N), Tcom_g(N), Ahigh_g(1,N), Alow_g(1,N))
! Identify that the species data was successfully read from the file.
               WRITE(*,1001) SPECIES_NAME(Nsp)        ! screen
               IF(DMP_LOG) WRITE(UNIT_LOG,1001) SPECIES_NAME(Nsp) ! log file
	              SPECIES_READ(Nsp) = .TRUE.
        	   ENDIF
	        ENDDO

! Read species data for the continuum solids phases.
!-----------------------------------------------------------------------
! Skip reading the database for the continuum solids phase if the 
! simulation is only employing discrete solids.
         IF(.NOT.DISCRETE_ELEMENT .OR. DES_CONTINUUM_HYBRID)THEN
          	 DO M = 1, MMAX
               DO N = 1, NMAX(M)
	                 Nsp = Nsp + 1
! Skip species that have already been read from a different source.
	                 IF(SPECIES_READ(Nsp)) CYCLE
	                 REWIND(UNIT=funit)

! If a species name was not specified in mfix.dat, flag error and exit.
	                 IF(SPECIES_NAME(NSP) == UNDEFINED_C)THEN
                     WRITE(*,1011)'continuum', M, N ! screen
                     IF(DMP_LOG) WRITE(UNIT_LOG,1011)'continuum', M, N
	                    CALL MFIX_EXIT(mypE)
	                 ENDIF

! Read the database.
	                 CALL Read_Therm(FUNIT, SPECIES_NAME(Nsp),            &
                     Thigh_s(M,N), Tlow_s(M,N), Tcom_s(M,N), MW,       &
                     Ahigh_s(1,M,N), Alow_s(1,M,N), HfrefoR_s(M,N), IER)

   	              IF(IER == 0)THEN
! If the user did not supply a value for the solids phase molecular
! weight of species in the mfix.dat file, use the value from the database.
	                    IF(MW_s(M,N)== UNDEFINED) MW_s(M,N) = MW
	       
! There are a number of species with Tlow as 300, for which the
! following calculation will produce an error because T_ref = 298.  So
! slightly extend validity of the correaltion.
                     IF( abs(Tlow_s(M,N)-T_ref) <= 2.0D0 .AND.         &
                        Tlow_s(M,N) > T_ref ) Tlow_s(M,N) = T_ref
! Calculate IC_PsrefoR.
   	                 IC_PsrefoR(M,N) = calc_ICpoR(T_ref, Thigh_s(M,N), &
                        Tlow_s(M,N), Tcom_s(M,N), Ahigh_s(1,M,N),      &
                        Alow_s(1,M,N))
! Identify the species as being successfully read from the database.
                     WRITE(*,1002)'continuum', M, SPECIES_NAME(Nsp)
                     IF(DMP_LOG) WRITE(UNIT_LOG,1002)'continuum', M,   &
                        SPECIES_NAME(Nsp)
	                    SPECIES_READ(Nsp) = .TRUE.
	                 ENDIF
	              ENDDO   ! N=1, NMAX(M)
	           ENDDO   ! M=1, MMAX
         ENDIF

! Read species data for discrete solids phases.
!-----------------------------------------------------------------------
         IF(DISCRETE_ELEMENT .OR. DES_CONTINUUM_HYBRID)THEN
          	 DO M = 1, MMAX
               DO N = 1, NMAX(M)
	                 DES_Nsp = DES_Nsp + 1
! Skip species that have already been read from a different source.
	                 IF(DES_SPECIES_READ(DES_Nsp)) CYCLE
	                 REWIND(UNIT=funit)
! If a species name was not specified in mfix.dat, flag error and exit.
	                 IF(DES_SPECIES_NAME(DES_Nsp) == UNDEFINED_C)THEN
                     WRITE(*,1011)'discrete', M, N ! screen
                     IF(DMP_LOG) WRITE(UNIT_LOG,1011)'discrete', M, N
	                    CALL MFIX_EXIT(mypE)
	                 ENDIF
! Read the database.
	                 CALL Read_Therm(FUNIT, DES_SPECIES_NAME(DES_Nsp),    &
                     DES_Thigh_s(M,N), DES_Tlow_s(M,N),                &
                     DES_Tcom_s(M,N),  MW, DES_Ahigh_s(1,M,N),         &
                     DES_Alow_s(1,M,N), DES_HfrefoR_s(M,N), IER)
   	              IF(IER == 0)THEN
! If the user did not supply a value for the solids phase molecular
! weight of species in the mfix.dat file, use the value from the database.
	                    IF(DES_MW_s(M,N)== UNDEFINED) DES_MW_s(M,N) = MW
! There are a number of species with Tlow as 300, for which the
! following calculation will produce an error because T_ref = 298.  So
! slightly extend validity of the correaltion.
                     IF( abs(DES_Tlow_s(M,N)-T_ref) <= 2.0D0 .AND.     &
                        DES_Tlow_s(M,N) > T_ref )                      &
                        DES_Tlow_s(M,N) = T_ref
! Calculate and store DES_C_PSoR at ref tempature
                     DES_IC_PSrefoR(M,N) = calc_ICpoR(T_ref, &
                        DES_Thigh_s(M,N), DES_Tlow_s(M,N), &
                        DES_Tcom_s(M,N), DES_Ahigh_s(1,M,N), &
                        DES_Alow_s(1,M,N))
! Identify the species as being successfully read from the database.
                     WRITE(*,1002)'discrete', M,                       &
                        DES_SPECIES_NAME(DES_Nsp)
                     IF(DMP_LOG) WRITE(UNIT_LOG,1002)'discrete', M,    &
                        DES_SPECIES_NAME(DES_Nsp)
	                    DES_SPECIES_READ(DES_Nsp) = .TRUE.
	                 ENDIF
	              ENDDO   ! N=1, NMAX(M)
	           ENDDO   ! M=1, MMAX
         ENDIF

! Close the file
	        CLOSE(UNIT=funit)
      ENDDO

! Error message control.
!-----------------------------------------------------------------------
! Check whether all continuum species data have been successfully read.
      ErrorFlag = .FALSE.
      NOT_WARNED = .TRUE.
      DO N = 1, Nsp
         IF(.NOT.SPECIES_READ(N))THEN
            IF(NOT_WARNED)THEN
               WRITE(*,1012)
               IF(DMP_LOG) WRITE(UNIT_LOG,1012)
               NOT_WARNED = .FALSE.
            ENDIF
            WRITE(*,1013) SPECIES_NAME(N)
            IF(DMP_LOG) WRITE(UNIT_LOG,1013) SPECIES_NAME(N)
	           ErrorFlag = .TRUE.
	        ENDIF
      ENDDO
! Check whether all continuum species data have been successfully read.
      DO N = 1, DES_Nsp
         IF(.NOT.DES_SPECIES_READ(N))THEN
            IF(NOT_WARNED)THEN
               WRITE(*,1012)
               IF(DMP_LOG) WRITE(UNIT_LOG,1012)
               NOT_WARNED = .FALSE.
            ENDIF
            WRITE(*,1013) DES_SPECIES_NAME(N)
            IF(DMP_LOG) WRITE(UNIT_LOG,1013) DES_SPECIES_NAME(N)
	           ErrorFlag = .TRUE.
	        ENDIF
      ENDDO
! Write remaining error message if needed.
      IF(ErrorFlag) THEN
        WRITE(*,1014) TRIM(FILENAME)
        IF(DMP_LOG) WRITE(UNIT_LOG,1014) TRIM(FILENAME)
        call mfix_exit(mype)
      ENDIF

      RETURN  

! Messages
!-----------------------------------------------------------------------
 1000 FORMAT(1X,'Checking for thermochemical data in ',A)
 1001 FORMAT(4X,'Species data read for gas phase species ',A)
 1002 FORMAT(4X,'Species data read for ',A,' solids phase ',I2,        &
         ' species ',A)

! Error Flags
!-----------------------------------------------------------------------
 1010 FORMAT(/1X,70('*')/, ' From: READ_DATABASE',/, ' Message: ',     &
         'No SPECIES_NAME provided for gas phase species ',I3,'.',/' ',&
         'Check mfix.dat.',/1X,70('*')/)
 1011 FORMAT(/1X,70('*')/, ' From: READ_DATABASE',/, ' Message: ',     &
         'No SPECIES_NAME provided for ',A,' solids phase ',I2,', ',/  &
         ' species ',I3,'.',/' Check mfix.dat.',/1X,70('*')/)
 1012 FORMAT(/1X,70('*')/, ' From: READ_DATABASE',/, ' Message: ',     &
         'No thermochemical data found in database for the following:')
 1013 FORMAT(4X,A)
 1014 FORMAT(' SUGGESTION: Search the database for the exact species', &
         ' name. The',/' species names are case sensitive and should', &
         ' mathch the names in',/' BURCAT.THR exactly excluding',      &
         ' trailing blanks and tabs. Also verify',/' that the data',   &
         ' section in the mfix.dat file (if any) is below a line',/    &
         ' that starts with THERMO DATA.',//' Database location:',     &
         /1X,A,/1X,70('*')/)


      END SUBROUTINE read_database
