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

!----------------------------------------------------------------------->>> 
      IF(FIRST_PASS) THEN
         CALL WRITE_DES_MASS
      ENDIF
!-----------------------------------------------------------------------<<<


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


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Module name: WRITE_DES_MASS                                         !
!                                                                      !
! 	THIS ROUTINE IS NOT DESIGNED TO BE CHECKED INTO CVS. IT IS STRICTLY !
!  HERE FOR PRINTING MESSAGES USED FOR DEBUGGING AND V&V WORK.         !
!                                                                      !
!  Author: J.Musser                                                    !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE WRITE_DES_MASS

      Use constant
      Use des_rxns
      Use des_thermo
      Use discretelement
      Use param1
      USE run

      IMPLICIT NONE

!-----------------------------------------------
! Local variables
!-----------------------------------------------
! file name
      CHARACTER*64 :: FNAME
! logical used for testing is the data file already exists
      LOGICAL :: F_EXISTS
! file unit for heat transfer data
      INTEGER, PARAMETER :: MASS_UNIT = 2034
      INTEGER, PARAMETER :: CORE_UNIT = 2038
      INTEGER, PARAMETER :: XB_UNIT = 2040
      INTEGER, PARAMETER :: XD_UNIT = 2041
      INTEGER, PARAMETER :: XI_UNIT = 2042

      INTEGER M

      DOUBLE PRECISION m_p0, r_B, MW_B, MW_D, M_p, ABS_ERR
      DOUBLE PRECISION X_B0, X_D0, X_I0
      DOUBLE PRECISION X_B, X_D, X_I
      DOUBLE PRECISION C_RAD
      DOUBLE PRECISION T_COMP


      r_B  = 0.00117810000665485d0

      M = PIJK(2,5)

      m_p0 = 0.0353429173528852d0
      MW_B = 15.0d0
      MW_D = 25.0d0

      X_B0 = 0.5d0
      X_D0 = 0.0d0
      X_I0 = 0.5d0

! Freeze the analytic solution at 14.99996484 seconds (the approximate time
! when the reaction should stop due to a lack of species B reactant).
      IF(S_TIME < 14.99996484d0) THEN
         T_COMP  = S_TIME
      ELSE
         T_COMP = 14.99996484d0
      ENDIF


      FNAME = TRIM(RUN_NAME)//'_MASS.dat'
      INQUIRE(FILE=FNAME,EXIST=F_EXISTS)
      IF (.NOT.F_EXISTS) THEN
         OPEN(UNIT=MASS_UNIT,FILE=FNAME,STATUS='NEW')

         write(MASS_UNIT,"(3X,A)")'Initial Conditions:'
         write(MASS_UNIT,"(15X,A,7X,A)")'Analytic','MFIX-DEM'
         write(MASS_UNIT,"(6X,A,3X,F12.8,3X,F12.8)")'Mass',m_p0,PMASS(2)

         write(MASS_UNIT,"(6X,A,3X,F12.8,3X,F12.8)")'MW_B',MW_B,DES_MW_s(M,1)
         write(MASS_UNIT,"(6X,A,3X,F12.8,3X,F12.8)")'MW_D',MW_D,DES_MW_s(M,2)
         write(MASS_UNIT,"(6X,A,4X,F12.8)")'r_B',r_B

         write(MASS_UNIT,"(//7X,A,10X,A,11X,A,8X,A)")&
            'S_TIME','Mp','Mp-DEM','ABS ERR'
      ELSE
         OPEN(UNIT=MASS_UNIT,FILE=FNAME,&
            POSITION="APPEND",STATUS='OLD')
      ENDIF

! Calculate the analytic solution for the particle's mass
      M_p = m_p0 + T_COMP*((MW_D/MW_B)/2.0d0 - 1.0d0)*r_B

! Calculate the absolute error
      ABS_ERR = abs(M_p - PMASS(2))
      WRITE(MASS_UNIT,"(4(3X,F12.8))")S_TIME, M_p, PMASS(2),ABS_ERR
      CLOSE(MASS_UNIT)



      FNAME = TRIM(RUN_NAME)//'_Xs_B.dat'
      INQUIRE(FILE=FNAME,EXIST=F_EXISTS)
      IF (.NOT.F_EXISTS) THEN
         OPEN(UNIT=XB_UNIT,FILE=FNAME,STATUS='NEW')

         write(XB_UNIT,"(7X,A,10X,A,11X,A,8X,A)")&
            'S_TIME','XB','XB-DEM','ABS ERR'
      ELSE
         OPEN(UNIT=XB_UNIT,FILE=FNAME,&
            POSITION="APPEND",STATUS='OLD')
      ENDIF
! Calculate the analytic solution for species B mass fraction
      X_B = (M_p0*X_B0 - r_B*T_COMP)/(M_p)
! Calculate the abolute error
      ABS_ERR = ABS(X_B - DES_X_s(2,1))
      WRITE(XB_UNIT,"(4(3X,F12.8))")S_TIME, X_B, DES_X_s(2,1), ABS_ERR
      CLOSE(XB_UNIT)



      FNAME = TRIM(RUN_NAME)//'_Xs_D.dat'
      INQUIRE(FILE=FNAME,EXIST=F_EXISTS)
      IF (.NOT.F_EXISTS) THEN
         OPEN(UNIT=XD_UNIT,FILE=FNAME,STATUS='NEW')

         write(XD_UNIT,"(7X,A,10X,A,11X,A,8X,A)")&
            'S_TIME','XD','XD-DEM','ABS ERR'
      ELSE
         OPEN(UNIT=XD_UNIT,FILE=FNAME,&
            POSITION="APPEND",STATUS='OLD')
      ENDIF
! Calculate the analytic solution for species B mass fraction
      X_D = (M_p0*X_D0 + T_COMP*((MW_D/MW_B)/2.0d0)*r_B)/M_p
! Calculate the abolute error
      ABS_ERR = ABS(X_D - DES_X_s(2,2))
      WRITE(XD_UNIT,"(4(3X,F12.8))")S_TIME, X_D, DES_X_s(2,2), ABS_ERR
      CLOSE(XD_UNIT)


      FNAME = TRIM(RUN_NAME)//'_Xs_I.dat'
      INQUIRE(FILE=FNAME,EXIST=F_EXISTS)
      IF (.NOT.F_EXISTS) THEN
         OPEN(UNIT=XI_UNIT,FILE=FNAME,STATUS='NEW')

         write(XI_UNIT,"(7X,A,10X,A,11X,A,8X,A)")&
            'S_TIME','XI','XI-DEM','ABS ERR'
      ELSE
         OPEN(UNIT=XI_UNIT,FILE=FNAME,&
            POSITION="APPEND",STATUS='OLD')
      ENDIF
! Calculate the analytic solution for species B mass fraction
      X_I = (M_p0*X_I0)/M_p
! Calculate the abolute error
      ABS_ERR = ABS(X_I - DES_X_s(2,3))
      WRITE(XI_UNIT,"(4(3X,F12.8))")S_TIME, X_I, DES_X_s(2,3), ABS_ERR
      CLOSE(XI_UNIT)


      FNAME = TRIM(RUN_NAME)//'_CORE.dat'
      INQUIRE(FILE=FNAME,EXIST=F_EXISTS)
      IF (.NOT.F_EXISTS) THEN
         OPEN(UNIT=CORE_UNIT,FILE=FNAME,STATUS='NEW')
         write(CORE_UNIT,"(7X,A,7X,A,5X,A,6X,A)")&
            'S_TIME','CORE RAD','CORE RAD-DEM','ABS ERR'

      ELSE
         OPEN(UNIT=CORE_UNIT,FILE=FNAME,&
            POSITION="APPEND",STATUS='OLD')
      ENDIF

      C_RAD = (DES_RADIUS(2)**3 - (T_COMP*r_B)/ &
         ((4.0d0/3.0d0)*PI*CORE_Rho(2)))**(1.0d0/3.0d0)

      ABS_ERR = ABS(C_RAD - CORE_RAD(2))

      WRITE(CORE_UNIT,"(4(3X,F12.8))")S_TIME, C_RAD, CORE_RAD(2), ABS_ERR
      CLOSE(CORE_UNIT)


      RETURN
      END SUBROUTINE WRITE_DES_MASS
