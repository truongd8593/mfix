!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: DES_THERMO_NEWVALUES                                   !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!                                                                      !
!  Author: J.Musser                                   Date: 16-Jun-10  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DES_THERMO_NEWVALUES(NP, FOCUS)

      Use des_thermo
      Use des_rxns
      Use discretelement
      Use param1
      Use physprop

      IMPLICIT NONE

! Passed variables
!-----------------------------------------------
! Index value of particle
      INTEGER, INTENT(IN) :: NP
! Logical indicating that the particle is of special interst
      LOGICAL, INTENT(IN) :: FOCUS

! Local variables
!-----------------------------------------------  

! index of solids phase
      INTEGER M
! Total amount of energy transfer TO the particle NP. If this value is
! positive, then the particle is heating, if it is negative, then the 
! particle is cooling.
      DOUBLE PRECISION  Qtotal

      DOUBLE PRECISION Qtotal_MAX
      DOUBLE PRECISION Qpp_MAX, Qpfp_MAX, Qcv_MAX, Qrd_MAX, Qint_MAX

! Logical for Adams-Bashfort integration.
      LOGICAL,SAVE:: FIRST_PASS = .TRUE.
      LOGICAL NOISY

! Initialize local variables
      Qcv_MAX = ZERO
      Qpp_MAX = ZERO
      Qpfp_MAX = ZERO
      Qrd_MAX = ZERO
      Qint_MAX = ZERO
      Qtotal_MAX = ZERO

! Initialize variables
      Qtotal = ZERO
      M = PIJK(NP,5)

! Particle-fluid convection
      IF(DES_CONV_EQ)THEN 
         Qtotal = Qtotal + Qcv(NP)
         Qcv_MAX = MAX( ABS(Qcv(NP)), Qcv_MAX)
! Clear storage variable
         Qcv(NP) = ZERO
      ENDIF

! Particle-particle conduction
      IF(DES_COND_EQ_PP) THEN
         Qtotal = Qtotal + Qpp(NP)
         Qpp_MAX = MAX( ABS(Qpp(NP)), Qpp_MAX)
! Clear storage variable
         Qpp(NP) = ZERO
      ENDIF

! Particle-fluid-particle conduction
      IF(DES_COND_EQ_PFP)THEN
         Qtotal = Qtotal + Qpfp(NP)
         Qpfp_MAX = MAX( ABS(Qpfp(NP)), Qpfp_MAX)
! Clear storage variable
         Qpfp(NP) = ZERO
      ENDIF

! Energy from radiation
      IF(DES_RADI_EQ)THEN
         Qtotal = Qtotal + Qrd(NP)
         Qrd_MAX = MAX( ABS(Qrd(NP)), Qrd_MAX)
! Clear storage variable
         Qrd(NP) = ZERO
      ENDIF

! Energy from reaction
      IF(DES_SPECIES_EQ(M))THEN
         Qtotal = Qtotal - Qint(NP)
         Qint_MAX = MAX( ABS(Qint(NP)), Qint_MAX)
! Clear storage variable
         Qint(NP) = ZERO
      ENDIF

      Qtotal_MAX = MAX( ABS(Qtotal), Qtotal_MAX)

! Advance particle position, velocity
      IF (TRIM(DES_INTG_METHOD) .EQ. 'EULER') THEN 
! First-order method              

!----------------------------------------------------------------------->>> REMOVE JM
         IF(DES_T_s_OLD(NP) .NE. DES_T_s_OLD(NP)) THEN
            WRITE(*,*)'   DES_T_s_OLD NAN FOR NP: ',NP
            STOP
         ENDIF
         IF(Qtotal .NE. Qtotal) THEN
            WRITE(*,*)'   Qtotal NAN FOR NP: ',NP
            WRITE(*,*)'      Qcv  : ',Qcv(NP)
            WRITE(*,*)'      Qpp  : ',Qpp(NP)
            WRITE(*,*)'      Qpfp : ',Qpfp(NP)
            WRITE(*,*)'      Qrd  : ',Qrd(NP)
            WRITE(*,*)'      Qint : ',Qint(NP)
            STOP
         ENDIF            
         
         IF(PMASS(NP) .NE. PMASS(NP))THEN
            WRITE(*,*)'   PMASS NAN FOR NP: ',NP
            STOP
         ENDIF

         IF(DES_C_ps(NP) .NE. DES_C_ps(NP))THEN
            WRITE(*,*)'   DES_C_ps NAN FOR NP: ',NP
            STOP
         ENDIF
!-----------------------------------------------------------------------<<<
         DES_T_s_NEW(NP) = DES_T_s_OLD(NP) + &
            DTSOLID*(Qtotal / ( PMASS(NP) * DES_C_ps(NP) ))
      ELSEIF (TRIM(DES_INTG_METHOD) .EQ. 'ADAMS_BASHFORTH') THEN
! Second-order Adams-Bashforth scheme
         IF(FIRST_PASS)THEN
            FIRST_PASS = .FALSE.
            DES_T_s_NEW(NP) = DES_T_s_OLD(NP) + &
               DTSOLID*(Qtotal / ( PMASS(NP) * DES_C_ps(NP) ))
            Qtotal_OLD(NP) = Qtotal / (PMASS(NP)*DES_C_ps(NP))
         ELSE
            DES_T_s_NEW(NP) = DES_T_s_OLD(NP) + &
              ( 1.5d0 * (Qtotal/(PMASS(NP)*DES_C_ps(NP))) - &
                0.5d0 * Qtotal_OLD(NP) ) * DTSOLID
            Qtotal_OLD(NP) = Qtotal / (PMASS(NP)*DES_C_ps(NP))
         ENDIF
      ENDIF

! Update the old temperature value
      DES_T_s_OLD(NP) = DES_T_s_NEW(NP)


!------------------------------------------------------------------------->>> 
      IF(OUTPUT_DATA_TIME .LE. S_TIME)THEN
         IF(NP==1) CALL WRITE_DES_TP
         OUTPUT_DATA_TIME = OUTPUT_DATA_TIME + 0.10d0
      ENDIF
!-------------------------------------------------------------------------<<<


      RETURN

      CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Module name: WRITE_DES_TP                                           !
!                                                                      !
! 	THIS ROUTINE IS NOT DESIGNED TO BE CHECKED INTO CVS. IT IS STRICTLY !
!  HERE FOR PRINTING MESSAGES USED FOR DEBUGGING AND V&V WORK.         !
!                                                                      !
!  Author: J.Musser                                                    !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE WRITE_DES_TP

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
      INTEGER, PARAMETER :: TP_UNIT = 2030

      DOUBLE PRECISION Tg, Tp0, Coeff, Tp,REL_ERR

      FNAME = TRIM(RUN_NAME)//'_DES_TP.dat'
      INQUIRE(FILE=FNAME,EXIST=F_EXISTS)
      IF (.NOT.F_EXISTS) THEN
         OPEN(UNIT=TP_UNIT,FILE=FNAME,STATUS='NEW')
         WRITE(TP_UNIT,"(7X,'Time',12X,'Tp',10X,'Tp_MFIX',8X,'REL ERR')")
      ELSE
         OPEN(UNIT=TP_UNIT,FILE=FNAME,&
            POSITION="APPEND",STATUS='OLD')
      ENDIF

! Calculate the value for the analytic solution.
      Tg = 298.15d0
      Tp0 = 373.14994357d0
      Coeff =  -0.878629815377096

      Tp = Tg - (Tg - Tp0)*exp(coeff*S_TIME)

! Calculate the relative error.
      REL_ERR = (ABS(Tp - DES_T_S_NEW(1))/ABS(Tp))*100

! Write the data to a file.
      WRITE(TP_UNIT,"(4(3X,F12.8))")S_TIME,Tp,DES_T_S_NEW(1),REL_ERR
      CLOSE(TP_UNIT)

      END SUBROUTINE WRITE_DES_TP

      END SUBROUTINE DES_THERMO_NEWVALUES
