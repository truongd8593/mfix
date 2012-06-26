!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_DES_RXNS                                          !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!                                                                      !
!  Author: J.Musser                                   Date: 17-Jun-11  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_DES_RXNS

!-----------------------------------------------
! Modules
!-----------------------------------------------      
      Use compar
      USE des_rxns
      Use discretelement
      Use funits  
      Use run
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Dummy loop indices for solids phase and species
      INTEGER M, N
! Number of processors used. (DES reactive chemistry is currently limited
! to serial runs!)
      INTEGER CHECK_MPI
! Logical indicating that an error message has been sent to the user.
      LOGICAL WARNED
!-----------------------------------------------

      FIRST_CALL = .TRUE.

      IF(DMP_LOG) WRITE(*,'(1X,A)') &
         '---------- START CHECK_DES_RXNS ---------->'

! Initialize MAX_DES_NMAX. This is the maximum number of species found
! over all solids phases. This is used for allocating DES_X_s
      MAX_DES_NMAX = 1

! Check DES_NMAX - this value must be defined if DES_SPECIES_EQ(M) is 
! true. 
      DO M = 1, DES_MMAX 
         IF (DES_NMAX(M) == UNDEFINED_I) THEN 
! If the species equation is being solved, then the user must specify
! the number of species in the discrete solids phase.
            IF (DES_SPECIES_EQ(M)) THEN 
               IF(DMP_LOG) THEN 
                  WRITE(*,1000)M
                  WRITE(UNIT_LOG,1000)M
               ENDIF
               CALL MFIX_EXIT(myPE)
            ELSE 
! Set the species count to be 1 if the speices count is not being 
! solved.
               DES_NMAX(M) = 1 
            ENDIF
         ELSE
            MAX_DES_NMAX = MAX(MAX_DES_NMAX, DES_NMAX(M))
! Ensure that the number of species specified by the user is below the
! max limit.
            IF (DES_NMAX(M) > DIM_N_S) THEN
               IF(DMP_LOG) THEN 
                  WRITE(*,1001)M
                  WRITE(UNIT_LOG,1001)M
               ENDIF
               CALL MFIX_EXIT(myPE)
            ENDIF 
         ENDIF

      ENDDO 

! Initialize error message flag.
      WARNED = .FALSE.
! Check if the molecular weights are provided by the user.
      DO M = 1, DES_MMAX
         IF (DES_SPECIES_EQ(M)) THEN
            DO N = 1, DES_NMAX(M)
               IF (DES_MW_s(M,N) == UNDEFINED) THEN

! If the molecular weight has not been provided in the mfix.dat file,
! notify the user that it will be obtained from the thermodynamic
! database.
                  IF(.NOT.WARNED)THEN
                     IF(DMP_LOG) THEN
                        WRITE(*,1005)
                        WRITE(UNIT_LOG, 1005)
                     ENDIF
                     WARNED = .TRUE.
                  ENDIF
                  IF(DMP_LOG) WRITE(*, 1006) M, N
                  IF(DMP_LOG) WRITE(UNIT_LOG) M, N
               ENDIF 
            ENDDO 
            IF(WARNED)THEN
               IF(DMP_LOG) THEN
                  WRITE(*,1007)
                  WRITE(UNIT_LOG,1007)
               ENDIF
            ENDIF

! If the molecular weight is specified for species with an index greater
! than DES_NMAX, flag the error and exit.
            DO N = DES_NMAX(M) + 1, DIM_N_S
               IF (DES_MW_S(M,N) /= UNDEFINED) THEN 
                  IF(DMP_LOG) THEN
                     WRITE(*,1008) M
                     WRITE(UNIT_LOG,1008) M
                  ENDIF
               ENDIF 
            ENDDO 
         ENDIF 
      ENDDO 

! Set the flag identifying that at least one of the species equations
! are being solved.
      ANY_DES_SPECIES_EQ = ANY(DES_SPECIES_EQ(:DES_MMAX))

! Check the number of processors. DES reactive chemistry is currently 
! limited to serial runs.
      IF(ANY_DES_SPECIES_EQ) THEN
         CHECK_MPI = NODESI * NODESJ * NODESK
         IF(CHECK_MPI.NE.1) THEN
            WRITE (*, 1009)
            WRITE (UNIT_LOG, 1009)
            CALL MFIX_EXIT(myPE)
         ENDIF
      ENDIF

      IF(DMP_LOG) WRITE(*,'(1X,A)') &
         '<---------- END CHECK_DES_RXNS ----------'

      RETURN

! Error Messages
!-----------------------------------------------------------------------
 1000 FORMAT(/1X,70('*')/, ' From: CHECK_DES_RXNS',/, ' Message: ',    &
         'Number of species for discrete phase ',I2,' not specified.', &
         /1X,70('*')/)

 1001 FORMAT(/1X,70('*')/, ' From: CHECK_DES_RXNS',/, ' Message: ',    &
         'DES_NMAX is too large for discrete phase ',I2,'.',/1X,70('*'))

 1002 FORMAT(/1X,70('*')/, ' From: CHECK_DES_RXNS',/, ' Message:',     &
         ' Species molecular weight is undefinded for phase ',I2,',',/ &
         ' species ',I3,' in mfix.dat file. These values will be',     &
         ' obtained from the',/'thermodynamic database.',/1X,70('*')/)

 1005 FORMAT(/1X,70('*')/, ' From: CHECK_DES_RXNS',/, ' Message:',     &
         ' Molecular weight not defined for the following:')

 1006 FORMAT(4X,'Solids phase: ',I2,' Species: ',I3)

 1007 FORMAT(' The thermodynamic databases will be consulted.',/       &
         1X,70('*')/)

 1008 FORMAT(/1X,70('*')/, ' From: CHECK_DES_RXNS',/, ' Message:',     &
         ' DES_MW_s is defined for an index value greater than',       &
         ' DES_NMAX for',/' solids phase ',I2,'. Check mfixdat file.',/&
         1X,70('*')/)

 1009 FORMAT(/1X,70('*')/ ' From: CHECK_DES_RXNS',/' Message: ',&
         'DES reacive chemisty modules are currently limited to serial',/&
         '  runs. Set nodesi, modesj, and nodesk to 1 in the mfix.dat',&
         ' file.',/1X,70('*')/)


      END SUBROUTINE CHECK_DES_RXNS
