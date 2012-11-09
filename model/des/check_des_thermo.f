!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_DES_THERMO                                        !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!                                                                      !
!  Author: J.Musser                                   Date: 17-Jun-10  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE CHECK_DES_THERMO

!-----------------------------------------------
! Modules
!-----------------------------------------------
      Use compar
      Use des_thermo
      Use des_rxns
      Use discretelement
      Use funits  
      Use interpolation
      Use physprop
      Use run
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------      

! Loop index
      INTEGER M, N ! Phase, Species

! Flag that the user was already warned why a call to the thermo-
! chemical database is being made.
      LOGICAL WARNED_USR

! Flag that the energy equations are solved and specified solids phase
! specific heat is undefined.
! If true, a call to the thermochemical database is made.
      LOGICAL EEQ_CPS

! Flag that the solids phase species equations are solved and the 
! molecular weight for a species are not given in the data file.
! If true, a call to the thermochemical database is made.
      LOGICAL SEQ_MWs
!-----------------------------------------------      

      IF(DMP_LOG) WRITE(*,'(1X,A)') &
         '---------- START CHECK_DES_THERMO ---------->'

! Set the flag identifying that at least one of the species equations
! are being solved.
      ANY_DES_SPECIES_EQ = ANY(DES_SPECIES_EQ(:DES_MMAX))

! Initialize MAX_DES_NMAX. This is the maximum number of species found
! over all solids phases. This is used for allocating DES_X_s
      MAX_DES_NMAX = 1

! Loop over all solids phases
      DO M=1,DES_MMAX
! Verify that the thermal conductivity values are physical.
         IF(DES_C_ps0(M) .LT. ZERO)THEN
            IF(DMP_LOG) THEN
               WRITE(UNIT_LOG,1001)'DES_C_ps0','unphysical',M
               WRITE(*,1001)'DES_C_ps0','unphysical',M
            ENDIF
            CALL MFIX_EXIT(myPE)
         ENDIF
! Verify that the thermal conductivity values are physical and defined.
         IF(DES_K_s0(M) .LT. ZERO)THEN
            IF(DMP_LOG) THEN
               WRITE(UNIT_LOG,1001)'DES_K_s0','unphysical',M
               WRITE(*,1001)'DES_K_s0','unphysical',M
            ENDIF
            CALL MFIX_EXIT(myPE)
         ENDIF
! Verify that a emmisivity value is specified for each solids pase
         IF(DES_Em(M) .LT. ZERO .OR. &
            DES_Em(M) .GT. ONE) THEN
            IF(DMP_LOG) THEN
               WRITE(UNIT_LOG,1001)'DES_Em','unphysical',M
               WRITE(*,1001)'DES_Em','unphysical',M
            ENDIF
            CALL MFIX_EXIT(myPE)
         ENDIF

! Check required values.
         IF(DES_COND_EQ .AND. DES_K_s0(M) == UNDEFINED) THEN
            IF(DMP_LOG) THEN
               WRITE(UNIT_LOG,1001)'DES_K_s0','undefined',M
               WRITE(*,1001)'DES_K_s0','undefined',M
            ENDIF
            CALL MFIX_EXIT(myPE)
         ENDIF

         IF(DES_RADI_EQ .AND. DES_Em(M) == UNDEFINED) THEN
            IF(DMP_LOG) THEN
               WRITE(UNIT_LOG,1001)'DES_Em','undefined',M
               WRITE(*,1001)'DES_Em','undefined',M
            ENDIF
            CALL MFIX_EXIT(myPE)
         ENDIF

! Flag indicating if the user was already warned.
         WARNED_USR = .FALSE.
! Flag that the energy equations are solved and specified solids phase
! specific heat is undefined.
         EEQ_CPS = .FALSE.
         IF(DES_ENERGY_EQ .AND. DES_C_PS0(M) == UNDEFINED)             &
            EEQ_CPS = .TRUE.

         DES_rDatabase(M,:) = .FALSE.


! Check DES_NMAX_s - This value must be defined if solving species eqs
         IF (DES_NMAX_s(M) == UNDEFINED_I) THEN 
! If the species equation is being solved, then the user must specify
! the number of species in the discrete solids phase.
            IF (DES_SPECIES_EQ(M)) THEN 
               IF(DMP_LOG) THEN 
                  WRITE(*,1006) M
                  WRITE(UNIT_LOG,1006) M
               ENDIF
               CALL MFIX_EXIT(myPE)
            ELSE 
! Set the species count to be 1 if the speices count is not being 
! solved.
               DES_NMAX_s(M) = 1 
            ENDIF
         ELSE
! Ensure that the number of species specified by the user is below the
! max limit.
            IF (DES_NMAX_s(M) > DIM_N_S) THEN
               IF(DMP_LOG) THEN 
                  WRITE(*,1005) M, DIM_N_S
                  WRITE(UNIT_LOG,1005) M, DIM_N_S
               ENDIF
               CALL MFIX_EXIT(myPE)
            ENDIF 
            MAX_DES_NMAX = MAX(MAX_DES_NMAX, DES_NMAX_s(M))
         ENDIF

         DO N = 1, DES_NMAX_s(M)
            SEQ_MWs = .FALSE.
            IF(DES_SPECIES_EQ(M) .AND. DES_MW_S(M,N) == UNDEFINED)     &
               SEQ_MWs = .TRUE.

! The thermodynamic data base provides the specific heat, molecular
! weights, and heat of formation.
! Check the thermodynamic database if:
!   1) the energy equation is solved and constant solids phase specific
!      heat isn't given. (EEQ_CPS = .TRUE.)
!   2) the species energy equation is solved and the molecular weights
!      for the solids phase species are not given. (SEQ_MWs = .TRUE.)
! A final thermochemical check is preformed in check_des_rxns. If
! neither of the above conditions result in species data being read from
! the database AND a particular species is referenced by a chemical
! equation then a call to read_database is forced.
            IF(EEQ_CPS .OR. SEQ_MWs) THEN
! Notify the user of the reason the thermochemical database is used.
               IF(.NOT.WARNED_USR) THEN
                  IF(EEQ_CPS .AND. DMP_LOG) THEN
                     WRITE(*,1002)
                     WRITE(UNIT_LOG,1002)
                  ENDIF
                  IF(SEQ_MWs .AND. DMP_LOG) THEN
                     WRITE(*,1003) M
                     WRITE(UNIT_LOG,1003) M
                  ENDIF
! Set a flag to prevent the user from getting the same message over
! and over again.
                  WARNED_USR = .TRUE.
               ENDIF
! Flag that the species name is not provided.
               IF(DES_SPECIES_s(M,N) == UNDEFINED_C) THEN
                  IF(DMP_LOG) THEN
                     WRITE(*,1004) M, N
                     WRITE(UNIT_LOG,1004) M, N
                  ENDIF
                  CALL MFIX_EXIT(myPE)
               ENDIF

! Read the database.
               IF(DMP_LOG) THEN
                  WRITE(*,1100) M, N
                  WRITE(UNIT_LOG,1100) M, N
               ENDIF
               CALL READ_DATABASE('DEM', M, N, DES_SPECIES_s(M,N),     &
                  DES_MW_S(M,N))
! Flag variable to stating that the database was read.
               DES_rDatabase(M,N) = .TRUE.
            ENDIF

         ENDDO ! DES_NMAX_s
      ENDDO ! DES_MMAX

      IF(DMP_LOG) WRITE(*,'(1X,A)')&
         '<---------- END CHECK_DES_THERMO ----------'

      RETURN

 1001 FORMAT(/1X,70('*')/, ' From: CHECK_DES_THERMO',/,' Error 1001: ',&
         A,' is ',A,' for discrete solids phase ',I2,/' Please',       &
         ' correct the data file.',/1X,70('*')/)

 1002 FORMAT(/1X,70('*')/' From: CHECK_DES_THERMO',/' Message 1002:',  &
         ' The discrete phase energy equations are being solved',/     &
         ' (DES_ENERGY_EQ) and the specified constant solids specific',&
         ' heat is',/' undefined (DES_C_PS0). Thus, the',              &
         ' thermochemical database will be used',/' to gather',        &
         ' specific heat data on the individual soids phase species.',/&
         1X,70('*')/)

 1003 FORMAT(/1X,70('*')/' From: CHECK_DES_THERMO',/' Message 1003:',  &
         ' Discrete solids phase ',I2,' species equations are being',/ &
         ' solved, and one or more species molecular weights are',     &
         ' undefined. Thus,',/' the thermochemical database will be',  &
         ' used to gather molecular weight',/' data on the solids',    &
         ' phase speicies.',/1X,70('*')/)

 1004 FORMAT(/1X,70('*')/' From: CHECK_DES_THERMO',/' Error 1004:',    &
         ' Discrete solids phase ',I2,' species ',I3,' name',          &
         ' (DES_SPECIES_s)',/' is undefined. Please correct the data', &
         ' file. ',/1X,70('*')/)

 1005 FORMAT(/1X,70('*')/, ' From: CHECK_DES_RXNS',/, ' Error 1005:',  &
         ' DES_NMAX_s is too large for discrete phase ',I2,'.',        &
         ' The maximum',/' number of species is ',I3,'.',/1X,70('*')/)

 1006 FORMAT(/1X,70('*')/, ' From: CHECK_DES_RXNS',/, ' Error 1006:',  &
         ' Number of species for discrete phase ',I2,' is not',        &
         ' specified.',/1X,70('*')/)


 1100 FORMAT(/'  Searching thermochemical databases for discrete',     &
         ' solids phase ',I2,', species ',I2)

      END SUBROUTINE CHECK_DES_THERMO
