!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: DES_REACTION_MODEL                                     !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!                                                                      !
!  Author: J.Musser                                   Date: 16-Jun-10  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DES_REACTION_MODEL(NP, FOCUS)

      Use constant
      Use des_rxns
      Use des_thermo
      Use discretelement
      Use param1

      IMPLICIT NONE

! Passed variables
!-----------------------------------------------
! Index value of particle
      INTEGER, INTENT(IN) :: NP
! Logical indicating that the specified particle is of special interest
      LOGICAL, INTENT(IN) :: FOCUS

! Local variables
!-----------------------------------------------  
! index of solids phase and species
      INTEGER M, N, NN
! total rate of consumption/production of species (g/sec)
      DOUBLE PRECISION SUM_DES_R_sc, SUM_DES_Rs
! masses of species comprising the particle (g)
      DOUBLE PRECISION S_MASS( DES_NMAX(PIJK(NP,5)) )

      DOUBLE PRECISION PERCENT_REDUCE
      LOGICAL REDUCED
      LOGICAL ALL_GONE( DES_NMAX(PIJK(NP,5)) )

      DOUBLE PRECISION, PARAMETER :: P43 = 4.0d0/3.0d0
      DOUBLE PRECISION DELTA_CORE
      DOUBLE PRECISION DOM, K1, K2, K3, K4

      DOUBLE PRECISION dMdt, dXdt, dRdt

! Logical for Adams-Bashfort integration.
      LOGICAL,SAVE:: FIRST_PASS = .TRUE.

! Initialize local variables
      M = PIJK(NP,5)

! Return to the calling routine if the species equation for this mass
! phase is not being solved.
      IF(.NOT.DES_SPECIES_EQ(M))RETURN

! Check to see if any of the reaction rate will consume more than the
! available species in the particle.
      REDUCED = .TRUE.
      ALL_GONE(:) = .FALSE.
      DO WHILE(REDUCED)
! Initialize local variables
         SUM_DES_Rs = ZERO
         SUM_DES_R_sc = ZERO
! Reset the flag
         REDUCED = .FALSE.
         DO N=1,DES_NMAX(M)
! Calculate the current mass of species N in the particle
            S_MASS(N) = DES_X_s(NP,N) * PMASS(NP)
            IF((S_MASS(N) - DTSOLID*DES_R_sc(NP,N)) < ZERO) THEN
! Flag that the consumption and production rates need reduced
               REDUCED = .TRUE.
! Indicate that the last of the species is gone.
               ALL_GONE(N) = .TRUE.
! Calculate the amount of reduction needed
               PERCENT_REDUCE =  S_MASS(N)/(DTSOLID*DES_R_sc(NP,N))
! Reduce the production/consumption rates for each species
               DO NN=1,DES_NMAX(M)
                  DES_R_sc(NP,NN) = DES_R_sc(NP,NN)*PERCENT_REDUCE
                  DES_R_sp(NP,NN) = DES_R_sp(NP,NN)*PERCENT_REDUCE
               ENDDO
            ELSE
! Calculate the total rate of change in particle mass
               SUM_DES_Rs = SUM_DES_Rs + (DES_R_sp(NP,N)-DES_R_sc(NP,N))
! Calculate the total rate of consumption
               SUM_DES_R_sc = SUM_DES_R_sc + DES_R_sc(NP,N)
            ENDIF
         ENDDO
      ENDDO

! Update the particle's mass.
! The mass of the particle is updated first so that it can be used in 
! updating the species mass percent of the particle.
!-----------------------------------------------------------------------
      dMdt = SUM_DES_Rs
      IF (TRIM(DES_INTG_METHOD) .EQ. 'EULER') THEN 
! First-order method: Euler
         PMASS(NP) = PMASS(NP) + DTSOLID*dMdt
      ELSEIF (TRIM(DES_INTG_METHOD) .EQ. 'ADAMS_BASHFORTH')THEN
! Second-order Adams-Bashforth scheme
         IF(FIRST_PASS)THEN
            PMASS(NP) = PMASS(NP) + DTSOLID*dMdt
         ELSE
            PMASS(NP) = PMASS(NP) + DTSOLID * &
               (1.50d0*dMdt-0.5d0*dMdt_OLD(NP))
         ENDIF
! Store the current value as the old value.
         dMdt_OLD(NP) = dMdt
      ENDIF


! Update the species mass percent
!-----------------------------------------------------------------------
      DO N=1,DES_NMAX(M)
         IF(.NOT.ALL_GONE(N))THEN
            dXdt = ( (DES_R_sp(NP,N)-DES_R_sc(NP,N)) - DES_X_s(NP,N) * &
               SUM_DES_Rs)/PMASS(NP)
            IF (TRIM(DES_INTG_METHOD) .EQ. 'EULER') THEN 
! First-order method: Euler
               DES_X_s(NP,N) = DES_X_s(NP,N) + DTSOLID*dXdt
            ELSEIF (TRIM(DES_INTG_METHOD) .EQ. 'ADAMS_BASHFORTH')THEN
! Second-order Adams-Bashforth scheme
               IF(FIRST_PASS)THEN
                  DES_X_s(NP,N) = DES_X_s(NP,N) + DTSOLID*dXdt
               ELSE
                  DES_X_s(NP,N) = DES_X_s(NP,N) + DTSOLID * &
                     (1.50d0*dXdt - 0.5d0*dXdt_OLD(NP,N))
               ENDIF
! Store the current value as the old value.
               dXdt_OLD(NP,N) = dXdt
            ENDIF
         ELSE
            DES_X_s(NP,N) = ZERO
         ENDIF
      ENDDO




! Shrinking core model.
!-----------------------------------------------------------------------
      IF(TRIM(REACTION_MODEL) == 'SHRINKING_CORE')THEN
         IF(SUM_DES_R_sc /= ZERO)THEN
            DOM = 4.0d0*Pi*CORE_Rho(NP)*(CORE_RAD(NP))**2
! Restrict the the size of the denominator to prevent floating point
! issues.
            IF( DOM > 1.0d-6)THEN
! Advance the solution to the core's radius
               dRdt = (-SUM_DES_R_sc)/DOM
               IF (TRIM(DES_INTG_METHOD) .EQ. 'EULER') THEN 
! First-order method: Euler
                  DELTA_CORE = DTSOLID*dRdt
               ELSEIF (TRIM(DES_INTG_METHOD) .EQ. 'ADAMS_BASHFORTH')THEN
! Second-order Adams-Bashforth scheme
                  IF(FIRST_PASS)THEN
                     FIRST_PASS = .FALSE.
                     DELTA_CORE = DTSOLID*dRdt
                  ELSE
                     DELTA_CORE = DTSOLID * &
                        (1.5d0*dRdt - 0.5d0*dRdt_OLD(NP))
                  ENDIF
! Store the current value as the old value.
                  dRdt_OLD(NP) = dRdt
               ENDIF
! If the radius will remain postive, update it's value
               IF((CORE_RAD(NP) + DELTA_CORE) > ZERO)THEN
                  CORE_RAD(NP) = CORE_RAD(NP) + DELTA_CORE
               ELSE
! To prevent the radius from being less than zero, limit value to zero.
                  CORE_RAD(NP) = ZERO
               ENDIF
            ELSE
! To prevent the radius from being less than zero, limit value to zero.
               CORE_RAD(NP) = ZERO
            ENDIF
         ENDIF
! update variables
         RO_Sol(NP) = PMASS(NP)/PVOL(NP)

! Shrinking particle model.
!-----------------------------------------------------------------------
      ELSEIF(TRIM(REACTION_MODEL) == 'SHRINKING_PARTICLE')THEN
         DES_RADIUS(NP) = (PMASS(NP)/(P43*Pi*Ro_Sol(NP)))**(1.0d0/3.0d0)
! update variables
         PVOL(NP) = PMASS(NP)/Ro_Sol(NP)

! Progressive conversion model
!-----------------------------------------------------------------------
      ELSEIF(TRIM(REACTION_MODEL) == 'PROGRESSIVE_CONVERSION')THEN
! update variables
         RO_Sol(NP) = PMASS(NP)/PVOL(NP)
      ENDIF

! Update one over the particle's moment of inertia
      OMOI(NP) = 5.0d0 / (2.0d0 * PMASS(NP) * DES_RADIUS(NP)**2) 

! Clear the necessary variables.
      DES_R_sp(NP,:) = ZERO
      DES_R_sc(NP,:) = ZERO

! Flag that the first pass is over
      IF(FIRST_PASS) FIRST_PASS = .FALSE.

      RETURN
      END SUBROUTINE DES_REACTION_MODEL
