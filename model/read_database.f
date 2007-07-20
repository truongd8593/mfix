!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: read_database(Ier)                                     C
!  Purpose: read thermochemical database                               C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-OCT-05  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE read_database(Ier) 

      USE param 
      USE param1 
      USE physprop
      USE constant
      USE compar
      USE rxns
      USE funits 

      IMPLICIT NONE
      
      CHARACTER(len=10) :: THERM = 'BURCAT.THR'
      INCLUDE 'mfix_directory_path.inc'
      CHARACTER(len=147) FILENAME
      
      DOUBLE PRECISION ::MW
      DOUBLE PRECISION, EXTERNAL ::calc_ICpoR
      INTEGER IER, M, N, Nsp
      integer funit, file, ios
      LOGICAL ErrorFlag
      LOGICAL SPECIES_READ(DIM_N_ALL)
      
      if(database_read)return
      database_read = .TRUE.
      SPECIES_READ(:) = .FALSE.

      funit = UNIT_DAT
      
      !Read data from mfix.dat or from BURCAT.THR in run directory or mfix directory
!      DO file = 1, 3
      file = 0
      DO
        file = file + 1
        ! Open files
        IF(file == 1) then
          OPEN(UNIT=funit, FILE='mfix.dat', STATUS='OLD', IOSTAT= ios) 
	  if(ios /= 0) cycle
	  PRINT *, 'Reading thermochemical data from mfix.dat'
          
	ELSEif(file == 2) then
          OPEN(UNIT=funit,FILE=TRIM(THERM), STATUS='OLD', IOSTAT= ios)
	  if(ios /= 0) cycle
          PRINT *, 'Reading thermochemical data from ', TRIM(THERM)

	ELSEif(file == 3) then
          FILENAME = trim(MFIX_PATH) // '/thermochemical/' // TRIM(THERM)
          OPEN(UNIT=funit,FILE=TRIM(FILENAME), STATUS='OLD', IOSTAT= ios)
	  if(ios /= 0) cycle
          PRINT *, 'Reading thermochemical data from ', TRIM(FILENAME)

	ELSE
          EXIT 
	ENDIF
      
        ! Read species data
        Nsp = 0      
        DO N = 1, NMAX(0)
	  Nsp = Nsp + 1
	  IF(SPECIES_READ(Nsp)) CYCLE
	  REWIND(UNIT=funit)
	  if(SPECIES_NAME(Nsp) == UNDEFINED_C)then
             PRINT *, 'read_database: Specify gas species ', N, ' species_name.'
	     CALL MFIX_EXIT(mypE)
	  endif
	  Call Read_Therm(funit, SPECIES_NAME(Nsp), Thigh_g(N), Tlow_g(N), Tcom_g(N), MW,&
	    Ahigh_g(1,N), Alow_g(1,N), HfrefoR_g(N), IER)
	    
	  IF(IER == 0)THEN
	    if(MW_g(N) == UNDEFINED) MW_g(N) = MW
            PRINT *, '  Read data for ',SPECIES_NAME(Nsp)
!           There are a number of species with Tlow as 300, for which the following calculation will 
!           produce an error because T_ref = 298.  So slightly extend validity of the correaltion
            if( abs(Tlow_g(N)-T_ref) <= 2.0D0 .and. Tlow_g(N) > T_ref ) Tlow_g(N) = T_ref
	    IC_PGrefoR(N) = calc_ICpoR(T_ref, Thigh_g(N), Tlow_g(N), Tcom_g(N), Ahigh_g(1,N), Alow_g(1,N))
	    SPECIES_READ(Nsp) = .TRUE.
	  ENDIF

	ENDDO
	DO M = 1, MMAX
          DO N = 1, NMAX(M)
	    Nsp = Nsp + 1
	    IF(SPECIES_READ(Nsp)) CYCLE
	    REWIND(UNIT=funit)
	    if(SPECIES_NAME(Nsp) == UNDEFINED_C)then
              PRINT *, 'read_database: Specify solids ', M, ' species ', N, ' species_name.'
	      CALL MFIX_EXIT(mypE)
	    endif
	    Call Read_Therm(funit, SPECIES_NAME(Nsp), Thigh_s(M,N), Tlow_s(M,N), Tcom_s(M,N),&
	      MW, Ahigh_s(1,M,N), Alow_s(1,M,N), HfrefoR_s(M,N), IER)
	      
	    IF(IER == 0)THEN
	      if(MW_s(M,N)== UNDEFINED) MW_s(M,N) = MW
	      
              PRINT *, '  Read data for ',SPECIES_NAME(Nsp)
!             There are a number of species with Tlow as 300, for which the following calculation will 
!             produce an error because T_ref = 298.  So slightly extend validity of the correaltion
              if( abs(Tlow_s(M,N)-T_ref) <= 2.0D0 .and. Tlow_s(M,N) > T_ref ) Tlow_s(M,N) = T_ref
	      IC_PsrefoR(M,N) = calc_ICpoR(T_ref, Thigh_s(M,N), Tlow_s(M,N), Tcom_s(M,N),&
	                         Ahigh_s(1,M,N), Alow_s(1,M,N))
	      SPECIES_READ(Nsp) = .TRUE.
	    ENDIF

	  ENDDO
	ENDDO
	CLOSE(UNIT=funit)
      ENDDO

!     Check whether all species data have been read
      ErrorFlag = .FALSE.
      DO N = 1, Nsp
        IF(.NOT.SPECIES_READ(N))THEN
            PRINT *, 'Thermochemical data for ',SPECIES_NAME(N), ' not found in Database!'
	    ErrorFlag = .TRUE.
	ENDIF
      ENDDO
      IF(ErrorFlag) THEN
        PRINT *, 'READ_DATABASE: Search for exact species name in the database file ', TRIM(FILENAME)
        PRINT *, 'Also ensure that the data section is below a line that starts with THERMO DATA '
        call mfix_exit(mype)
      ENDIF

      RETURN  
      END SUBROUTINE read_database
