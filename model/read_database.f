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

      IMPLICIT NONE
      DOUBLE PRECISION ::MW
      DOUBLE PRECISION, EXTERNAL ::calc_ICpoR
      INTEGER IER, M, N, Nsp
      CHARACTER PATH*147
      INCLUDE 'mfix_directory_path.inc'

      if(database_read)return
      database_read = .TRUE.
      
      PATH = trim(MFIX_PATH)//'/thermochemical'

        Nsp = 0      
        DO N = 1, NMAX(0)
	  Nsp = Nsp + 1
	  if(SPECIES_NAME(Nsp) == UNDEFINED_C)then
             PRINT *, 'read_database: Specify gas species ', N, ' species_name.'
	     CALL MFIX_EXIT(mypE)
	  endif
	  Call Read_Therm(PATH, SPECIES_NAME(Nsp), Thigh_g(N), Tlow_g(N), Tcom_g(N), MW,&
	    Ahigh_g(1,N), Alow_g(1,N), HfrefoR_g(N))
	  if(MW_g(N) == UNDEFINED) MW_g(N) = MW
          PRINT *, 'Read data for ',SPECIES_NAME(Nsp)
!         There are a number of species with Tlow as 300, for which the following calculation will 
!         produce an error because T_ref = 298.  So slightly extend validity of the correaltion
          if( Tlow_g(N)-T_ref <= 2.0D0) Tlow_g(N) = 298.D0
	  IC_PGrefoR(N) = calc_ICpoR(T_ref, Thigh_g(N), Tlow_g(N), Tcom_g(N), Ahigh_g(1,N), Alow_g(1,N))

	ENDDO
	DO M = 1, MMAX
          DO N = 1, NMAX(M)
	    Nsp = Nsp + 1
	    if(SPECIES_NAME(Nsp) == UNDEFINED_C)then
               PRINT *, 'read_database: Specify solids ', M, ' species ', N, ' species_name.'
	       CALL MFIX_EXIT(mypE)
	    endif
	    Call Read_Therm(PATH, SPECIES_NAME(Nsp), Thigh_s(M,N), Tlow_s(M,N), Tcom_s(M,N),&
	      MW, Ahigh_s(1,M,N), Alow_s(1,M,N), HfrefoR_s(M,N))
	    if(MW_s(M,N)== UNDEFINED) MW_s(M,N) = MW
	      
            PRINT *, 'Read data for ',SPECIES_NAME(Nsp)
!         There are a number of species with Tlow as 300, for which the following calculation will 
!         produce an error because T_ref = 298.  So slightly extend validity of the correaltion
          if( Tlow_s(M,N)-T_ref <= 2.0D0) Tlow_s(M,N) = 298.D0
	    IC_PsrefoR(M,N) = calc_ICpoR(T_ref, Thigh_s(M,N), Tlow_s(M,N), Tcom_s(M,N),&
	                         Ahigh_s(1,M,N), Alow_s(1,M,N))

	  ENDDO
	ENDDO

      RETURN  
      END SUBROUTINE read_database
