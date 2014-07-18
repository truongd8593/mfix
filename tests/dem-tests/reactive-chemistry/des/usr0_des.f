!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: URS0_DES                                               !
!                                                                      !
!  Purpose: This routine is called before the discrete phase time loop !
!  and is user-definable. The user may insert code in this routine or  !
!  call appropriate user defined subroutines.                          !
!                                                                      !
!  This routien is not called from a loop, hence all indicies are      !
!  undefined.                                                          !
!                                                                      !
!  Author: J.Musser                                   Date: 06-Nov-12  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE USR0_DES

      Use usr

      IMPLICIT NONE

      LOGICAL, SAVE :: FIRST_PASS = .TRUE.

      IF(FIRST_PASS) OUTPUT_DATA_TIME = ZERO

      FIRST_PASS = .FALSE.

      RETURN
      END SUBROUTINE USR0_DES




!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Module name: WRITE_DES_MASS                                         !
!                                                                      !
!  THIS ROUTINE IS NOT DESIGNED TO BE CHECKED INTO CVS. IT IS STRICTLY !
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
      use physprop, only: MW_s

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

         write(MASS_UNIT,"(6X,A,3X,F12.8,3X,F12.8)")'MW_B',MW_B,MW_s(M,1)
         write(MASS_UNIT,"(6X,A,3X,F12.8,3X,F12.8)")'MW_D',MW_D,MW_s(M,2)
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

      C_RAD = UNDEFINED
      ABS_ERR = UNDEFINED
!      C_RAD = (DES_RADIUS(2)**3 - (T_COMP*r_B)/ &
!         ((4.0d0/3.0d0)*PI*CORE_Rho(2)))**(1.0d0/3.0d0)

!      ABS_ERR = ABS(C_RAD - CORE_RAD(2))

!      WRITE(CORE_UNIT,"(4(3X,F12.8))")S_TIME, C_RAD, CORE_RAD(2), ABS_ERR
      CLOSE(CORE_UNIT)


      RETURN
      END SUBROUTINE WRITE_DES_MASS



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Module name: WRITE_DES_TP                                           !
!                                                                      !
!  THIS ROUTINE IS NOT DESIGNED TO BE CHECKED INTO CVS. IT IS STRICTLY !
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


