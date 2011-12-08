!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Subroutine name: PHYSICAL_PROP                                      !
!                                                                      !
!  Purpose: Calculate physical properties that vary with time.         !
!                                                                      !
!  Author: J.Musser                                   Date: 09-May-11  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DES_PHYSICAL_PROP(NP, FOCUS)

      Use des_rxns
      Use des_thermo
      Use discretelement
      Use funits 
      Use param
      Use param1 
      Use physprop
      Use run

      IMPLICIT NONE

! Passed Variables
!-----------------------------------------------------------------------
! Index of particle
      INTEGER, INTENT(IN) :: NP
! Logical indicating to write additional data about the particle for
! debugging purposes.
      LOGICAL, INTENT(IN) :: FOCUS

! Local Variables
!-----------------------------------------------------------------------
! Gas constant
      DOUBLE PRECISION, PARAMETER :: RGAS = 1.987207D0  !cal/mol.K
! Dummy indices
      INTEGER M, N
! error indicator
      INTEGER IER

! External fuctions
!-----------------------------------------------------------------------
! Calculate the specific heat from polynomical data obtained from the
! thermodynamic databases.
      DOUBLE PRECISION , EXTERNAL :: calc_CpoR

! Get the solids phase of the particle
      M = PIJK(NP,5)

! Specific heat
!-----------------------------------------------------------------------
! This only needs calculated when solving the energy equations.
      IF(DES_ENERGY_EQ) THEN
! If a constant value specific heat has not been defined in the mfix.dat
! file, calculate the temperature dependent value based on data from the
! thermodynamic databases.
         IF (DES_C_PS0(M) == UNDEFINED) THEN
! Read the thermodynamic database if it has not already been done.
            IF(.NOT.DATABASE_READ) CALL READ_DATABASE(IER)
            DES_C_PS(NP) = ZERO
! Calculate the specific heat based on the species composition of the
! particle and the data from the thermodynamic databases.
            DO N = 1, DES_NMAX(M)
    	          DES_C_PS(NP) = DES_C_PS(NP) + DES_X_s(NP,N) *           &
                  calc_CpoR(DES_T_s_NEW(NP), DES_Thigh_s(M, N),        &
 		               DES_Tlow_s(M, N), DES_Tcom_s(M, N),                  &
                  DES_Ahigh_s(1,M,N), DES_Alow_s(1,M,N))*              &
                  RGAS/DES_MW_s(M,N)
            ENDDO
! Convert to SI units if needed.
            IF (UNITS == 'SI') DES_C_PS(NP) = 4183.925D0*DES_C_PS(NP)
         ELSE
! If a constant value specific heat has been assigned to the particle
! in the mfix.dat file, use this value.
            DES_C_PS(NP) = DES_C_PS0(M)
         ENDIF
      ENDIF

!      WRITE(*,"(3X,A,I2,A,F10.6)")'DES_C_PS(',NP,'): ',DES_C_PS(NP)

      RETURN

      END SUBROUTINE DES_PHYSICAL_PROP
