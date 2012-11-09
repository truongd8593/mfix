
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  SUBROUTINE Name: DES_SET_IC                                         !
!                                                                      !
!  Purpose: Assign initial conditions to particles basded upon thier   !
!  location within the domain.                                         !
!                                                                      !
!  Author: J.Musser                                   Date: 15-Feb-11  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DES_SET_IC

      USE compar
      USE des_thermo
      USE discretelement
      USE des_ic
      USE des_rxns
      USE funits

      IMPLICIT NONE

! Dummy indices
      INTEGER PC, I, M, ICV

! X,Y, and Z values of particle position
      DOUBLE PRECISION X_POS, Y_POS, Z_POS

! Logical indicating that no IC was define in the region containing
! the specified particle.
      LOGICAL IC_FOUND

      LOGICAL , EXTERNAL :: COMPARE 


! Loop through the particles in the system, skipping entring, exiting
! and "non-existing" particles.
      PC = 1
      DO I = 1, MAX_PIP
         IF(PC .GT. PIP) EXIT
         IF(.NOT.PEA(I,1)) CYCLE
! Set the solids phase of the particle
         M = PIJK(I,5)
         X_POS = DES_POS_NEW(I,1)
         Y_POS = DES_POS_NEW(I,2)
         IF(DIMN == 3) Z_POS = DES_POS_NEW(I,3)
         IC_FOUND = .FALSE.

! Loop through the initial conditions
         DO ICV = 1, DIMENSION_IC 
            IF (DES_IC_DEFINED(ICV)) THEN 
! Check to see if the position of the particle is located inside the
! initial condition region.
               IF(X_POS .GE. DES_IC_X_W(ICV) .AND. &
                  X_POS .LE. DES_IC_X_E(ICV) .AND. &
                  Y_POS .GE. DES_IC_Y_S(ICV) .AND. &
                  Y_POS .LE. DES_IC_Y_N(ICV) ) THEN
                  IF(DIMN == 2 .OR. &
                    (Z_POS .GE. DES_IC_Z_B(ICV) .AND. &
                     Z_POS .LE. DES_IC_Z_T(ICV)) ) THEN
                  ENDIF
! Change the flag to indicate that an IC regioon has been found.
                  IC_FOUND = .TRUE.
! Exit the loop and assign the initial conditions to particle I.
                  EXIT
               ENDIF
            ENDIF
         ENDDO

         IF(IC_FOUND)THEN
! Assign the initial values for the specified models.
            CALL SET_DES_TEMPERATURE(I,ICV,M)
            CALL SET_DES_MASS_FRACTION(I,ICV,M)
         ELSE
            IF(DIMN == 2) Z_POS = 0.0d0
            WRITE(*,1000)I, X_POS, Y_POS, Z_POS
            WRITE(UNIT_LOG,1000)I, X_POS, Y_POS, Z_POS
            CALL MFIX_EXIT(myPE)
         ENDIF
! icrement the particle counter
         PC = PC + 1
      ENDDO

 1000 FORMAT(/1X,70('*')/,' From: SET_DES_IC',/, ' Message:',          &
         ' There is no initial condition associated with the region',/ &
         ' containing particle ',I6,'.',/,' X Position: ',F9.4, /      &
         ' Y Position: ',F9.4,/' Z Position: ',F9.4,/' Please check',  &
         ' the mfix.dat file.',/                                       &
         1X,70('*')/)

      CONTAINS

!......................................................................!
! Subroutine: SET_DES_TEMPERATURE                                      !
!                                                                      !
! Purpose: Set the initial temperature of a partile.                   !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE SET_DES_TEMPERATURE(I,ICV,M)
      INTEGER, INTENT(IN) :: I, ICV, M

! Initial temperatures must be specified for particles when solving
! the internal energy equation and species equation for a particle.
      IF(.NOT.ANY_DES_SPECIES_EQ  .AND. .NOT.DES_ENERGY_EQ ) RETURN

      IF(DES_IC_T_s(ICV,M) == UNDEFINED) THEN
         WRITE(*,1000)ICV,M,I
         WRITE(UNIT_LOG,1000)ICV,M,I
         CALL MFIX_EXIT(myPE)
      ELSE
         DES_T_s_NEW(I) = DES_IC_T_s(ICV,M)
      ENDIF

! Set the value of the old temperature to the initial value.
      DES_T_s_OLD(I) = DES_T_s_NEW(I)

      RETURN

 1000 FORMAT(/1X,70('*')/,' From: SET_DES_IC',/, ' Message: The',      &
         ' initial temperature is UNDEFINED for particle ',I6,'.',/    &
         ' This is required when solving the internal energy equation',&
         ' and species',/' equation for a particle. Check mfix.dat',   &
         ' file.',/1X,70('*')/)

      END SUBROUTINE SET_DES_TEMPERATURE

!......................................................................!
! Subroutine: SET_DES_MASS_FRAC                                        !
!                                                                      !
! Purpose: Set the initial species mass fraction of a partile.         !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE SET_DES_MASS_FRACTION(I,ICV,M)
      INTEGER, INTENT(IN) :: I, ICV, M

      INTEGER N ! species loop index

      DOUBLE PRECISION SUM

! The initial species mass fractions for a particle must be supplied if
! the species equation is being solved or if DES_C_PS0 has not been
! specified in conjunction with the energy equations being solved.
      IF(.NOT.ANY_DES_SPECIES_EQ  .AND. (.NOT.DES_ENERGY_EQ .OR. &
         (DES_ENERGY_EQ .AND. DES_C_PS0(M) /= UNDEFINED) ) ) RETURN

! Verify that the species mass fraction sums to one.
      SUM = ZERO
      DO N=1,DES_NMAX_s(M)
         IF(DES_IC_X_s(ICV,M,N) /= UNDEFINED) THEN
            DES_X_s(I,N) = DES_IC_X_s(ICV,M,N)
            SUM = SUM + DES_IC_X_s(ICV,M,N)
         ELSE
            DES_X_s(I,N) = ZERO
         ENDIF
      ENDDO
      IF(.NOT.COMPARE(ONE,SUM))THEN
         WRITE(*,1000)I
         WRITE(UNIT_LOG,1000)I
         CALL MFIX_EXIT(myPE)
      ENDIF

      RETURN

 1000 FORMAT(/1X,70('*')/,' From: SET_DES_IC',/, ' Message: The',      &
         ' species mass fraction of particle ',I6,' does not sum',/    &
         ' to one. This is required if sovling the discrete particle', &
         ' species',/' equation or if DES_C_PS0 is not specified when',&
         ' solving the internal',/' energy equation for a particle.',  &
         ' Check mfix.dat file.',/1X,70('*')/)

      END SUBROUTINE SET_DES_MASS_FRACTION


      END SUBROUTINE DES_SET_IC
