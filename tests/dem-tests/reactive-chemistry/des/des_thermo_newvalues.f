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

!------------------------------------------------------------------------->>>>
! Prohibit any modes of heat transfer except HOR.
         Qcv(NP) = ZERO
!-------------------------------------------------------------------------<<<<
         Qtotal = Qtotal + Qcv(NP)
         Qcv_MAX = MAX( ABS(Qcv(NP)), Qcv_MAX)
! Clear storage variable
         Qcv(NP) = ZERO
      ENDIF

! Particle-particle conduction
      IF(DES_COND_EQ_PP) THEN
!------------------------------------------------------------------------->>>>
! Prohibit any modes of heat transfer except HOR.
         Qpp(NP) = ZERO
!-------------------------------------------------------------------------<<<<
         Qtotal = Qtotal + Qpp(NP)
         Qpp_MAX = MAX( ABS(Qpp(NP)), Qpp_MAX)
! Clear storage variable
         Qpp(NP) = ZERO
      ENDIF

! Particle-fluid-particle conduction
      IF(DES_COND_EQ_PFP)THEN
!------------------------------------------------------------------------->>>>
! Prohibit any modes of heat transfer except HOR.
         Qpfp(NP) = ZERO
!-------------------------------------------------------------------------<<<<
         Qtotal = Qtotal + Qpfp(NP)
         Qpfp_MAX = MAX( ABS(Qpfp(NP)), Qpfp_MAX)
! Clear storage variable
         Qpfp(NP) = ZERO
      ENDIF

! Energy from radiation
      IF(DES_RADI_EQ)THEN
!------------------------------------------------------------------------->>>>
! Prohibit any modes of heat transfer except HOR.
         Qrd(NP) = ZERO
!-------------------------------------------------------------------------<<<<
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

!----------------------------------------------------------------------->>> 
      IF(OUTPUT_DATA_TIME .LE. S_TIME)THEN
         IF(NP==2) CALL WRITE_DES_TP
      ENDIF
!-----------------------------------------------------------------------<<<

      RETURN

      END SUBROUTINE DES_THERMO_NEWVALUES


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

      Use des_thermo
      Use des_rxns
      Use discretelement
      Use param1
      Use physprop
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


      DOUBLE PRECISION, PARAMETER :: a0_A = 3.5d0
      DOUBLE PRECISION, PARAMETER :: a0_B = 1.4d0
      DOUBLE PRECISION, PARAMETER :: a0_C = 4.6d0
      DOUBLE PRECISION, PARAMETER :: a0_D = 1.5d0
      DOUBLE PRECISION, PARAMETER :: a0_I = 3.0d0
      DOUBLE PRECISION, PARAMETER :: MW_A = 25.0d0
      DOUBLE PRECISION, PARAMETER :: MW_B = 15.0d0
      DOUBLE PRECISION, PARAMETER :: MW_C = 30.0d0
      DOUBLE PRECISION, PARAMETER :: MW_D = 25.0d0
      DOUBLE PRECISION, PARAMETER :: MW_I = 20.0d0
      DOUBLE PRECISION, PARAMETER :: HoR_C = -1666.67d0
      DOUBLE PRECISION, PARAMETER :: r_B  = 0.00117810000665485d0
      DOUBLE PRECISION, PARAMETER :: UGC = 1.987207d0
      DOUBLE PRECISION, PARAMETER :: m_p0 = 0.0353429173528852d0
      DOUBLE PRECISION, PARAMETER :: X_B0 = 0.5d0
      DOUBLE PRECISION, PARAMETER :: X_D0 = 0.0d0
      DOUBLE PRECISION, PARAMETER :: X_I0 = 0.5d0
      DOUBLE PRECISION, PARAMETER :: TREF = 298.0d0


      DOUBLE PRECISION, PARAMETER :: RK4_DT_DEFAULT = 0.00001

      DOUBLE PRECISION, SAVE :: RK4_TIME = 0.0d0
      DOUBLE PRECISION, SAVE :: Tp_RK4 = 500.0d0

      DOUBLE PRECISION RK4_DT, RK4_DT_LAST
      DOUBLE PRECISION TIME_INTERVAL
      DOUBLE PRECISION K1, K2, K3, K4

      INTEGER I, RK4_STEPS

      DOUBLE PRECISION a1, a2, a3, X_B, X_D, t

      DOUBLE PRECISION REL_ERR


      X_B(t) = (m_p0*X_B0 - r_B*t)/(m_p0 + t*(0.5d0*(MW_D/MW_B)-1.0d0)*r_B)

      X_D(t) = (m_p0*X_D0 + 0.5d0*(MW_D/MW_B)*r_B*t) / &
         (m_p0 + t*(0.5d0*(MW_D/MW_B)-1.0d0)*r_B)

      a1(t) = (X_D(t) * 0.5d0*(MW_D/MW_B)*r_B * a0_D *UGC)/MW_D - &
           (X_B(t)*r_B*a0_B* UGC)/MW_B

      a2(t) = (X_B(t)*(         r_B         )*a0_B*TREF*UGC)/MW_B - &
              (X_D(t)*(0.5d0*(MW_D/MW_B)*r_B)*a0_D*TREF*UGC)/MW_D + &
              (r_B - 0.5d0*(MW_D/MW_B)*r_B)*HoR_C*(UGC/MW_C)

      a3(t) = ((m_p0*X_B0 - (         r_B*t         ))*a0_B*UGC)/MW_B + &
              ((m_p0*X_D0 + (0.5d0*(MW_D/MW_B)*r_B*t))*a0_D*UGC)/MW_D + &
              (m_p0*X_I0*a0_I*UGC)/MW_I


      FNAME = TRIM(RUN_NAME)//'_DES_TP.dat'
      INQUIRE(FILE=FNAME,EXIST=F_EXISTS)
      IF (.NOT.F_EXISTS) THEN
         OPEN(UNIT=TP_UNIT,FILE=FNAME,STATUS='NEW')
         write(TP_UNIT,"(7X,A,10X,A,11X,A,8X,A)")&
            'S_TIME','Tp','Tp-DEM','REL ERR'
      ELSE
         OPEN(UNIT=TP_UNIT,FILE=FNAME,&
            POSITION="APPEND",STATUS='OLD')
      ENDIF


! Freeze the analytic solution at 14.99996484 seconds (the approximate time
! when the reaction should stop due to a lack of species B reactant).
      IF(S_TIME < 14.99996484d0) THEN

         TIME_INTERVAL = S_TIME - RK4_TIME
         IF(TIME_INTERVAL .LE. RK4_DT) THEN
            RK4_STEPS = 1
            RK4_DT = TIME_INTERVAL
            RK4_DT_LAST = UNDEFINED
         ELSE
            RK4_STEPS = floor(real(TIME_INTERVAL/RK4_DT_DEFAULT))
            RK4_DT = RK4_DT_DEFAULT
            RK4_DT_LAST = S_TIME - (RK4_TIME + RK4_STEPS*RK4_DT)
         ENDIF


         DO I=1, RK4_STEPS

            K1 = RK4_DT * (-1.0d0/a3(RK4_TIME)) * &
               (a1(RK4_TIME)*(TP_RK4) + &
                a2(RK4_TIME))

            K2 = RK4_DT * (-1.0d0/a3(RK4_TIME + 0.5d0*RK4_DT)) * &
               (a1(RK4_TIME + 0.5d0*RK4_DT)*(TP_RK4 + 0.5d0*K1) + &
                a2(RK4_TIME + 0.5d0*RK4_DT))

            K3 = RK4_DT * (-1.0d0/a3(RK4_TIME + 0.5d0*RK4_DT)) * &
               (a1(RK4_TIME + 0.5d0*RK4_DT)*(TP_RK4 + 0.5d0*K2) + &
                a2(RK4_TIME + 0.5d0*RK4_DT))


            K4 = RK4_DT * (-1.0d0/a3(RK4_TIME + RK4_DT)) * &
               (a1(RK4_TIME + RK4_DT)*(TP_RK4 + K3) + &
                a2(RK4_TIME + RK4_DT))

            TP_RK4 = TP_RK4 + (K1 + 2.0d0*K2 + 2.0d0*K3 + K4)/6.0d0

            RK4_TIME = RK4_TIME + RK4_DT

         ENDDO

         

         IF(RK4_DT_LAST .NE. UNDEFINED) THEN

            RK4_DT = RK4_DT_LAST

            K1 = RK4_DT * (-1.0d0/a3(TP_RK4)) * &
               (a1(RK4_TIME)*(TP_RK4) + a2(RK4_TIME))

            K2 = RK4_DT * (-1.0d0/a3(TP_RK4 + 0.5d0*RK4_DT)) * &
               (a1(RK4_TIME + 0.5d0*RK4_DT)*(TP_RK4 + 0.5d0*K1) + &
                a2(RK4_TIME + 0.5d0*RK4_DT))

            K3 = RK4_DT * (-1.0d0/a3(TP_RK4 + 0.5d0*RK4_DT)) * &
               (a1(RK4_TIME + 0.5d0*RK4_DT)*(TP_RK4 + 0.5d0*K2) + &
                a2(RK4_TIME + 0.5d0*RK4_DT))


            K4 = RK4_DT * (-1.0d0/a3(TP_RK4 + RK4_DT)) * &
               (a1(RK4_TIME + RK4_DT)*(TP_RK4 + K3) + &
                a2(RK4_TIME + RK4_DT))


            TP_RK4 = TP_RK4 + (K1 + 2.0d0*K2 + 2.0d0*K3 + K4)/6.0d0

            RK4_TIME = RK4_TIME + RK4_DT

         ENDIF
      ENDIF

      REL_ERR = ABS(TP_RK4 - DES_T_S_NEW(2))/TP_RK4

      WRITE(TP_UNIT,"(4(3X,F12.8))")S_TIME, TP_RK4, DES_T_S_NEW(2), REL_ERR
      CLOSE(TP_UNIT)

      RETURN
      END SUBROUTINE WRITE_DES_TP
