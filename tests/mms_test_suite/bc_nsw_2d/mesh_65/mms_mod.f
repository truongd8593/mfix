!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: MMS                                                    !
!  Purpose: Global storage containers for MMS variables.               !
!                                                                      !
!  Author: Aniruddha Choudhary                        Date: 17-Oct-11  !
!  email: anirudd@vt.edu                                               !
!                                                                      !
!  Reviewer: J.Musser                                 Date: 04-Dec-13  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      MODULE mms

! By default the MMS functions are unavailable.
      LOGICAL, parameter :: USE_MMS = .TRUE.

!! Method of Manufactured Solutions (MMS) and Tecplot variables :

! Gas volume fraction
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MMS_Ep_g
! Gas pressure
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MMS_P_g

! Gas velocity components
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MMS_U_g
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MMS_V_g
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MMS_W_g

! Gas temperature
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MMS_T_g

! Solid bulk density
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MMS_ROP_s
! Solids velocity components
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MMS_U_s
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MMS_V_s
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MMS_W_s
! Solids temperature
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MMS_T_s
! Granular temperature
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MMS_Theta_m

! Gas continuity MMS source term
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MMS_ROP_g_Src

! Gas Momentum source terms
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MMS_U_g_Src
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MMS_V_g_Src
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MMS_W_g_Src

! Gas energy equation MMS source term
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MMS_T_g_Src

! Solid continuity MMS source term
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MMS_ROP_s_Src
! Solid momentum MMS source terms
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MMS_U_s_Src
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MMS_V_s_Src
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MMS_W_s_Src

! Solid energy equation MMS source term
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MMS_T_s_Src

! Granular energy MMS source term
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MMS_theta_m_Src

! Gemporary variable for pressure shifting while plotting and
! discretization error norm calculation.
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  P_G_Sh

! Index for pressure shifting !
      INTEGER :: IJK_Sh

!  x, y, z locations of top-right corner of a cell
        DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  xtr, ytr, ztr

      contains



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name:CALCULATE_MMS                                           !
!  Purpose: Calculate manufactured solutions                           !
!                                                                      !
!  Author: Aniruddha Choudhary                        Date: 17-Oct-11  !
!  email: anirudd@vt.edu                                               !
!                                                                      !
!  Reviewer: J.Musser                                 Date: 04-Dec-13  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALCULATE_MMS

      USE param
      USE param1
      USE parallel
      USE constant
      USE run
      USE toleranc
      USE geometry
      USE indices
      USE compar
      USE sendrecv
      USE fldvar

      IMPLICIT NONE

! Indices
        INTEGER  :: I, J, K, IJK, II, JJ, KK
        INTEGER  :: AllocateStatus

! MMS_Function(xt, yt, zt, index)
! index: 1=Pg, 2=Ug, 3=Vg, 4=Wg 5=Us 6=Vs 7=Ws 8=Tg 9=Ts &
!        10=Ep_g 11=ROP_s 12=Theta_m !
        DOUBLE PRECISION :: xt, yt, zt

! temporary variables holding x, y z locations
       INCLUDE 'function.inc'

       if(myPE == pe_io) write(*,"(3x,'Calculating MMS')")

        !! Allocate MMS and plotting variables !!
!
       If(.not.allocated(MMS_Ep_G)) ALLOCATE (MMS_Ep_G(DIMENSION_3), STAT = AllocateStatus)
!       If (AllocateStatus /= 0) STOP "*** Not enough memory ***"
!
       If(.not.allocated(MMS_P_G)) ALLOCATE (MMS_P_G(DIMENSION_3), STAT = AllocateStatus)
!       If (AllocateStatus /= 0) STOP "*** Not enough memory ***"
!
       If(.not.allocated(MMS_U_G)) ALLOCATE (MMS_U_G(DIMENSION_3), STAT = AllocateStatus)
!       If (AllocateStatus /= 0) STOP "*** Not enough memory ***"
!
       If(.not.allocated(MMS_V_G)) ALLOCATE (MMS_V_G(DIMENSION_3), STAT = AllocateStatus)
!       If (AllocateStatus /= 0) STOP "*** Not enough memory ***"
!
       If(.not.allocated(MMS_W_G)) ALLOCATE (MMS_W_G(DIMENSION_3), STAT = AllocateStatus)
!       If (AllocateStatus /= 0) STOP "*** Not enough memory ***"
!
       If(.not.allocated(MMS_T_G)) ALLOCATE (MMS_T_G(DIMENSION_3), STAT = AllocateStatus)
!       If (AllocateStatus /= 0) STOP "*** Not enough memory ***"
!
       If(.not.allocated(MMS_ROP_s)) ALLOCATE (MMS_ROP_s(DIMENSION_3), STAT = AllocateStatus)
!       If (AllocateStatus /= 0) STOP "*** Not enough memory ***"
!
       If(.not.allocated(MMS_U_S)) ALLOCATE (MMS_U_S(DIMENSION_3), STAT = AllocateStatus)
!       If (AllocateStatus /= 0) STOP "*** Not enough memory ***"
!
       If(.not.allocated(MMS_V_S)) ALLOCATE (MMS_V_S(DIMENSION_3), STAT = AllocateStatus)
!       If (AllocateStatus /= 0) STOP "*** Not enough memory ***"
!
       If(.not.allocated(MMS_W_S)) ALLOCATE (MMS_W_S(DIMENSION_3), STAT = AllocateStatus)
!       If (AllocateStatus /= 0) STOP "*** Not enough memory ***"
!
       If(.not.allocated(MMS_T_S)) ALLOCATE (MMS_T_S(DIMENSION_3), STAT = AllocateStatus)
!       If (AllocateStatus /= 0) STOP "*** Not enough memory ***"
!
       If(.not.allocated(MMS_THETA_M)) ALLOCATE (MMS_THETA_M(DIMENSION_3), STAT = AllocateStatus)
!       If (AllocateStatus /= 0) STOP "*** Not enough memory ***"
!
       If(.not.allocated(P_G_Sh)) ALLOCATE (P_G_Sh(DIMENSION_3), STAT = AllocateStatus)
!       If (AllocateStatus /= 0) STOP "*** Not enough memory ***"
!
       If(.not.allocated(xtr)) ALLOCATE (xtr(DIMENSION_3), STAT = AllocateStatus)
!       If (AllocateStatus /= 0) STOP "*** Not enough memory ***"
!
       If(.not.allocated(ytr)) ALLOCATE (ytr(DIMENSION_3), STAT = AllocateStatus)
!       If (AllocateStatus /= 0) STOP "*** Not enough memory ***"
!
       If(.not.allocated(ztr)) ALLOCATE (ztr(DIMENSION_3), STAT = AllocateStatus)
!       If (AllocateStatus /= 0) STOP "*** Not enough memory ***"
!

      !! Set reference point for shifting pressure 
      IJK_Sh = FUNIJK_GL(IMAX1/2+1,JMAX1/2+1,KMAX/2+1)

      !! Generate grid locations for plotting and MMS calculations
      DO IJK = ijkstart3, ijkend3

          I = I_OF(IJK)
          J = J_OF(IJK)
          K = K_OF(IJK)

          xt = ZERO - DX(1)
          yt = ZERO - DY(1)
          zt = ZERO - DZ(1)

          DO II = 1, I
            xt = xt + DX(II)
          END DO

          DO JJ = 1, J
            yt = yt + DY(JJ)
          END DO

          DO KK = 1, K
            zt = zt + DZ(KK)
          END DO

          xtr(IJK) = xt
          ytr(IJK) = yt
          ztr(IJK) = zt

      ENDDO

      call send_recv(xtr,2)
      call send_recv(ytr,2)
      call send_recv(ztr,2)


      !! Set MMS variables
      DO IJK = ijkstart3, ijkend3

         I = I_OF(IJK)
         J = J_OF(IJK)
         K = K_OF(IJK)

! Scalar variables !
         xt = xtr(IJK) - DX(I)*HALF
         yt = ytr(IJK) - DY(J)*HALF
         zt = ztr(IJK) - DZ(K)*HALF
         MMS_P_G(IJK) = MMS_Function(xt, yt, zt, 1)
         MMS_T_G(IJK) = MMS_Function(xt, yt, zt, 8)
         MMS_T_S(IJK) = MMS_Function(xt, yt, zt, 9)
         MMS_EP_G(IJK) = MMS_Function(xt, yt, zt, 10)
         MMS_ROP_S(IJK) = MMS_Function(xt, yt, zt, 11)
         MMS_THETA_M(IJK) = MMS_Function(xt, yt, zt, 12)

! U-velocity MMS !
         xt = xtr(IJK)
         yt = ytr(IJK) - DY(J)*HALF
         zt = ztr(IJK) - DZ(K)*HALF
         MMS_U_G(IJK) = MMS_Function(xt, yt, zt, 2)
         MMS_U_S(IJK) = MMS_Function(xt, yt, zt, 5)

! V-velocity MMS !
         xt = xtr(IJK) - DX(I)*HALF
         yt = ytr(IJK)
         zt = ztr(IJK) - DZ(K)*HALF
         MMS_V_G(IJK) = MMS_Function(xt, yt, zt, 3)
         MMS_V_S(IJK) = MMS_Function(xt, yt, zt, 6)

! W-velocity MMS !
         xt = xtr(IJK) - DX(I)*HALF
         yt = ytr(IJK) - DY(J)*HALF
         zt = ztr(IJK)
         MMS_W_G(IJK) = MMS_Function(xt, yt, zt, 4)
         MMS_W_S(IJK) = MMS_Function(xt, yt, zt, 7)

      ENDDO


      call send_recv(MMS_P_G,2)
      call send_recv(MMS_U_G,2)
      call send_recv(MMS_V_G,2)
      call send_recv(MMS_W_G,2)
      call send_recv(MMS_U_S,2)
      call send_recv(MMS_V_S,2)
      call send_recv(MMS_W_S,2)
      call send_recv(MMS_T_G,2)
      call send_recv(MMS_T_S,2)
      call send_recv(MMS_EP_G,2)
      call send_recv(MMS_ROP_S,2)
      call send_recv(MMS_THETA_M,2)

      END SUBROUTINE CALCULATE_MMS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Function: MMS_Function                                              !
!  Purpose:                                                            !
!                                                                      !
!  Author: Aniruddha Choudhary                        Date: 17-Oct-11  !
!  email: anirudd@vt.edu                                               !
!                                                                      !
!  Reviewer: J.Musser                                 Date: 04-Dec-13  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      DOUBLE PRECISION Function MMS_Function(xt, yt, zt, idx)

      USE constant
      USE physprop

      IMPLICIT NONE

! temporary coordinates
      DOUBLE PRECISION :: xt, yt, zt

! index for the variable
! index: 1=Pg, 2=Ug, 3=Vg, 4=Wg 5=Us 6=Vs 7=Ws 8=Tg 9=Ts &
!        10=Ep_g 11=ROP_s 12=Theta_m !
      INTEGER :: idx

! Set coefficients in MMS variables
! Note that for solid volume fraction, all outer coefficients except es0 are set 
! as 0.0d0 which makes this a constant volume fraction case.
!      DOUBLE PRECISION :: pg0=1.0d+5, pgx=0.2d+5, pgy=-0.5d+5, pgz=0.2d+5, &
!        pgxy=-0.25d+5, pgyz=-0.1d+5, pgzx=0.1d+5, &
!        apgx=0.4d0, apgy=0.45d0, apgz=0.85d0, apgxy=0.75d0, apgyz=0.7d0,&
!        apgzx=0.8d0   
      DOUBLE PRECISION :: pg0=100.0d0, pgx=20.0d0, pgy=-50.0d0, pgz=20.0d0, &
        pgxy=-25.0d0, pgyz=-10.0d0, pgzx=10.0d0, &
        apgx=0.4d0, apgy=0.45d0, apgz=0.85d0, apgxy=0.75d0, apgyz=0.7d0,&
        apgzx=0.8d0
      DOUBLE PRECISION :: ug0=7.0d0, ugx=3.0d0, ugy=-4.0d0, ugz=-3.0d0, &
        ugxy=2.0d0, ugyz=1.5d0, ugzx=-2.0d0, &
        augx=0.5d0, augy=0.85d0, augz=0.4d0, augxy=0.6d0, augyz=0.8d0, augzx=0.9d0  
      DOUBLE PRECISION :: vg0=9.0d0, vgx=-5.0d0, vgy=4.0d0, vgz=5.0d0, &
        vgxy=-3.0d0, vgyz=2.5d0, vgzx=3.5d0, &
        avgx=0.8d0, avgy=0.8d0, avgz=0.5d0, avgxy=0.9d0, avgyz=0.4d0, avgzx=0.6d0
      DOUBLE PRECISION :: wg0=8.0d0, wgx=-4.0d0, wgy=3.5d0, wgz=4.2d0, &
        wgxy=-2.2d0, wgyz=2.1d0, wgzx=2.5d0, &
        awgx=0.85d0, awgy=0.9d0, awgz=0.5d0, awgxy=0.4d0, awgyz=0.8d0, awgzx=0.75d0
      DOUBLE PRECISION :: us0=5.0d0
      DOUBLE PRECISION :: vs0=5.0d0
      DOUBLE PRECISION :: ws0=5.0d0
      DOUBLE PRECISION :: Tg0=350.0d0, Tgx=10.0d0, Tgy=-30.0d0, Tgz=20.0d0, &
        Tgxy=-12.0d0, Tgyz=10.0d0, Tgzx=8.0d0, &
        aTgx=0.75d0, aTgy=1.25d0, aTgz=0.8d0, aTgxy=0.65d0, aTgyz=0.5d0, aTgzx=0.6d0
      DOUBLE PRECISION :: Ts0=300.0d0, Tsx=15.0d0, Tsy=-20.0d0, Tsz=15.0d0, &
        Tsxy=-10.0d0, Tsyz=12.0d0, Tszx=10.0d0, &
        aTsx=0.5d0, aTsy=0.9d0, aTsz=0.8d0, aTsxy=0.5d0, aTsyz=0.65d0, aTszx=0.4d0
      DOUBLE PRECISION :: es0=0.0d0, esx=0.0d0, esy=0.0d0, esz=0.0d0, &
        esxy=0.0d0, esyz=0.0d0, eszx=0.0d0, &
        aesx=0.5d0, aesy=0.5d0, aesz=0.5d0, aesxy=0.4d0, aesyz=0.4d0, aeszx=0.4d0
!      DOUBLE PRECISION :: es0=0.1d0, esx=0.0d0, esy=0.0d0, esz=0.0d0, &
!        esxy=0.0d0, esyz=0.0d0, eszx=0.0d0, &
!        aesx=0.5d0, aesy=0.5d0, aesz=0.5d0, aesxy=0.4d0, aesyz=0.4d0, aeszx=0.4d0
      DOUBLE PRECISION :: Ths0=100.0d0, Thsx=5.0d0, Thsy=-10.0d0, Thsz=12.0d0, &
        Thsxy=-8.0d0, Thsyz=10.0d0, Thszx=7.0d0, &
        aThsx=0.8d0, aThsy=1.25d0, aThsz=0.7d0, aThsxy=0.5d0, aThsyz=0.6d0, aThszx=0.7d0

      DOUBLE PRECISION :: ros
  
      ros= RO_s0(1)
 
      ! NSW
      ug0=10.0d0
      ugx=1.0d0
      vg0=10.0d0
      pgx=100.0d0

      ! FSW
!      ug0=10.0d0
!      vg0=5.0d0

      !!FLAGMMS
      !write(*,*) "MMS_Function: ros= ", ros
  
      SELECT CASE(idx)
      CASE(1)
      !Pg!
        MMS_Function = pg0 + pgx*xt**2*Cos(Pi*(xt + yt))
      CASE(2)
      !Ug!
        MMS_Function = ug0*xt**2*(1 + Sin(Pi*(xt + yt))**2)
      CASE(3)
      !Vg!
        MMS_Function = xt*(5*ug0 - 3*ug0*yt + (ug0*xt*Cos(2*Pi*xt)*Cos(2*Pi*yt))/2. + &
         (ug0*Cos(2*Pi*yt)*Sin(2*Pi*xt))/(2.*Pi) + (ug0*Cos(2*Pi*xt)*Sin(2*Pi*yt))/(2.*Pi) - &
         (ug0*xt*Sin(2*Pi*xt)*Sin(2*Pi*yt))/2.)
      CASE(4)
      !Wg!
        MMS_Function = ZERO
      Case(5)
      !Us!
        MMS_Function = ZERO
      Case(6)
      !Vs!
        MMS_Function = ZERO
      Case(7)
      !Ws!
        MMS_Function = ZERO
      Case(8)
      !Tg!
        MMS_Function = ZERO
      Case(9)
      !Ts!
        MMS_Function = ZERO
      Case(10)
      !Ep_g!
        MMS_Function = 1.0d0
      Case(11)
      !ROP_s!
        MMS_Function = ZERO
      Case(12)
      !Theta_m!
        MMS_Function = ZERO
      END SELECT

      End Function MMS_Function





!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Function: CALCULATE_MMS_SOURCE                                      !
!  Purpose:                                                            !
!                                                                      !
!  Author: Aniruddha Choudhary                        Date: 17-Oct-11  !
!  email: anirudd@vt.edu                                               !
!                                                                      !
!  Reviewer: J.Musser                                 Date: 04-Dec-13  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALCULATE_MMS_SOURCE

      USE param
      USE param1
      USE parallel
      USE constant
      USE run
      USE toleranc
      USE geometry
      USE indices
      USE compar
      USE sendrecv
      USE fldvar

      IMPLICIT NONE
!                      indices
        INTEGER  :: I, J, K, IJK
        INTEGER  :: AllocateStatus

        ! X, Y cell-center locations (not including ghost cells) !
        DOUBLE PRECISION :: xt, yt, zt
          ! temporary variables holding x, y, z locations

       INCLUDE 'function.inc'

        !! Allocate MMS source term variables !!
!
       If(.not.allocated(MMS_ROP_G_SRC)) ALLOCATE (MMS_ROP_G_SRC(DIMENSION_3), STAT = AllocateStatus)
!       If (AllocateStatus /= 0) STOP "*** Not enough memory ***"
!
       If(.not.allocated(MMS_U_G_SRC)) ALLOCATE (MMS_U_G_SRC(DIMENSION_3), STAT = AllocateStatus)
!       If (AllocateStatus /= 0) STOP "*** Not enough memory ***"
!
       If(.not.allocated(MMS_V_G_SRC)) ALLOCATE (MMS_V_G_SRC(DIMENSION_3), STAT = AllocateStatus)
!       If (AllocateStatus /= 0) STOP "*** Not enough memory ***"
!
       If(.not.allocated(MMS_W_G_SRC)) ALLOCATE (MMS_W_G_SRC(DIMENSION_3), STAT = AllocateStatus)
!       If (AllocateStatus /= 0) STOP "*** Not enough memory ***"
!
       If(.not.allocated(MMS_T_G_SRC)) ALLOCATE (MMS_T_G_SRC(DIMENSION_3), STAT = AllocateStatus)
!       If (AllocateStatus /= 0) STOP "*** Not enough memory ***"
!
       If(.not.allocated(MMS_ROP_S_SRC)) ALLOCATE (MMS_ROP_S_SRC(DIMENSION_3), STAT = AllocateStatus)
!       If (AllocateStatus /= 0) STOP "*** Not enough memory ***"
!
       If(.not.allocated(MMS_U_S_SRC)) ALLOCATE (MMS_U_S_SRC(DIMENSION_3), STAT = AllocateStatus)
!       If (AllocateStatus /= 0) STOP "*** Not enough memory ***"
!
       If(.not.allocated(MMS_V_S_SRC)) ALLOCATE (MMS_V_S_SRC(DIMENSION_3), STAT = AllocateStatus)
!       If (AllocateStatus /= 0) STOP "*** Not enough memory ***"
!
       If(.not.allocated(MMS_W_S_SRC)) ALLOCATE (MMS_W_S_SRC(DIMENSION_3), STAT = AllocateStatus)
!       If (AllocateStatus /= 0) STOP "*** Not enough memory ***"
!
       If(.not.allocated(MMS_T_S_SRC)) ALLOCATE (MMS_T_S_SRC(DIMENSION_3), STAT = AllocateStatus)
!       If (AllocateStatus /= 0) STOP "*** Not enough memory ***"
!
       If(.not.allocated(MMS_THETA_M_SRC)) ALLOCATE (MMS_THETA_M_SRC(DIMENSION_3), STAT = AllocateStatus)
!       If (AllocateStatus /= 0) STOP "*** Not enough memory ***"

      DO IJK = ijkstart3, ijkend3

         I = I_OF(IJK)
         J = J_OF(IJK)
         K = K_OF(IJK)

! Scalar variables
         xt = xtr(IJK) - DX(I)*HALF
         yt = ytr(IJK) - DY(J)*HALF
         zt = ztr(IJK) - DZ(K)*HALF
         MMS_T_G_SRC(IJK) = MMS_Source(xt, yt, zt, 8)
         MMS_T_S_SRC(IJK) = MMS_Source(xt, yt, zt, 9)
         MMS_ROP_G_SRC(IJK) = MMS_Source(xt, yt, zt, 10)
         MMS_ROP_S_SRC(IJK) = MMS_Source(xt, yt, zt, 11)
         MMS_THETA_M_SRC(IJK) = MMS_Source(xt, yt, zt, 12)

! U-velocity MMS !
         xt = xtr(IJK)
         yt = ytr(IJK) - DY(J)*HALF
         zt = ztr(IJK) - DZ(K)*HALF
        MMS_U_g_Src(IJK) = MMS_Source(xt, yt, zt, 2)
        MMS_U_s_Src(IJK) = MMS_Source(xt, yt, zt, 5)


! V-velocity MMS !
         xt = xtr(IJK) - DX(I)*HALF
         yt = ytr(IJK)
         zt = ztr(IJK) - DZ(K)*HALF
        MMS_V_g_Src(IJK) = MMS_Source(xt, yt, zt, 3)
        MMS_V_s_Src(IJK) = MMS_Source(xt, yt, zt, 6)

! W-velocity MMS !
         xt = xtr(IJK) - DX(I)*HALF
         yt = ytr(IJK) - DY(J)*HALF
         zt = ztr(IJK)
        MMS_W_g_Src(IJK) = MMS_Source(xt, yt, zt, 4)
        MMS_W_s_Src(IJK) = MMS_Source(xt, yt, zt, 7)

      ENDDO


      call send_recv(MMS_U_g_Src,2)
      call send_recv(MMS_V_g_Src,2)
      call send_recv(MMS_W_g_Src,2)
      call send_recv(MMS_U_s_Src,2)
      call send_recv(MMS_V_s_Src,2)
      call send_recv(MMS_W_s_Src,2)
      call send_recv(MMS_T_G_SRC,2)
      call send_recv(MMS_T_S_SRC,2)
      call send_recv(MMS_ROP_G_SRC,2)
      call send_recv(MMS_ROP_S_SRC,2)
      call send_recv(MMS_THETA_M_SRC,2)

      End SUBROUTINE CALCULATE_MMS_SOURCE


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Function: MMS_Source                                                !
!  Purpose:                                                            !
!                                                                      !
!  Author: Aniruddha Choudhary                        Date: 17-Oct-11  !
!  email: anirudd@vt.edu                                               !
!                                                                      !
!  Reviewer: J.Musser                                 Date: 04-Dec-13  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      DOUBLE PRECISION Function MMS_Source(xt, yt, zt, idx)

      USE constant
      USE physprop

      IMPLICIT NONE

      DOUBLE PRECISION :: xt, yt, zt
      ! temporary coordinates

      INTEGER :: idx
      ! index for the variable
      ! index: 1=--, 2=Ug_src, 3=Vg_src, 4=Wg_Src 5=Us_src 6=Vs_src 7=Ws_src
      ! 8=Tg_src 9=Ts_src 10=rop_g_src 11=rop_s_src 12=Theta_m_src !

      !! Set coefficients in MMS variables
      !Note that for solid volume fraction, all outer coefficients except es0 are set 
      !as 0.0d0 which makes this a constant volume fraction case.
!      DOUBLE PRECISION :: pg0=1.0d+5, pgx=0.2d+5, pgy=-0.5d+5, pgz=0.2d+5, &
!        pgxy=-0.25d+5, pgyz=-0.1d+5, pgzx=0.1d+5, &
!        apgx=0.4d0, apgy=0.45d0, apgz=0.85d0, apgxy=0.75d0, apgyz=0.7d0,&
!        apgzx=0.8d0   
      DOUBLE PRECISION :: pg0=100.0d0, pgx=20.0d0, pgy=-50.0d0, pgz=20.0d0, &
        pgxy=-25.0d0, pgyz=-10.0d0, pgzx=10.0d0, &
        apgx=0.4d0, apgy=0.45d0, apgz=0.85d0, apgxy=0.75d0, apgyz=0.7d0,&
        apgzx=0.8d0  
      DOUBLE PRECISION :: ug0=7.0d0, ugx=3.0d0, ugy=-4.0d0, ugz=-3.0d0, &
        ugxy=2.0d0, ugyz=1.5d0, ugzx=-2.0d0, &
        augx=0.5d0, augy=0.85d0, augz=0.4d0, augxy=0.6d0, augyz=0.8d0, augzx=0.9d0  
      DOUBLE PRECISION :: vg0=9.0d0, vgx=-5.0d0, vgy=4.0d0, vgz=5.0d0, &
        vgxy=-3.0d0, vgyz=2.5d0, vgzx=3.5d0, &
        avgx=0.8d0, avgy=0.8d0, avgz=0.5d0, avgxy=0.9d0, avgyz=0.4d0, avgzx=0.6d0
      DOUBLE PRECISION :: wg0=8.0d0, wgx=-4.0d0, wgy=3.5d0, wgz=4.2d0, &
        wgxy=-2.2d0, wgyz=2.1d0, wgzx=2.5d0, &
        awgx=0.85d0, awgy=0.9d0, awgz=0.5d0, awgxy=0.4d0, awgyz=0.8d0, awgzx=0.75d0
      DOUBLE PRECISION :: us0=5.0d0
      DOUBLE PRECISION :: vs0=5.0d0
      DOUBLE PRECISION :: ws0=5.0d0
      DOUBLE PRECISION :: Tg0=350.0d0, Tgx=10.0d0, Tgy=-30.0d0, Tgz=20.0d0, &
        Tgxy=-12.0d0, Tgyz=10.0d0, Tgzx=8.0d0, &
        aTgx=0.75d0, aTgy=1.25d0, aTgz=0.8d0, aTgxy=0.65d0, aTgyz=0.5d0, aTgzx=0.6d0
      DOUBLE PRECISION :: Ts0=300.0d0, Tsx=15.0d0, Tsy=-20.0d0, Tsz=15.0d0, &
        Tsxy=-10.0d0, Tsyz=12.0d0, Tszx=10.0d0, &
        aTsx=0.5d0, aTsy=0.9d0, aTsz=0.8d0, aTsxy=0.5d0, aTsyz=0.65d0, aTszx=0.4d0
      DOUBLE PRECISION :: es0=0.0d0, esx=0.0d0, esy=0.0d0, esz=0.0d0, &
        esxy=0.0d0, esyz=0.0d0, eszx=0.0d0, &
        aesx=0.5d0, aesy=0.5d0, aesz=0.5d0, aesxy=0.4d0, aesyz=0.4d0, aeszx=0.4d0
!      DOUBLE PRECISION :: es0=0.1d0, esx=0.0d0, esy=0.0d0, esz=0.0d0, &
!        esxy=0.0d0, esyz=0.0d0, eszx=0.0d0, &
!        aesx=0.5d0, aesy=0.5d0, aesz=0.5d0, aesxy=0.4d0, aesyz=0.4d0, aeszx=0.4d0
      DOUBLE PRECISION :: Ths0=100.0d0, Thsx=5.0d0, Thsy=-10.0d0, Thsz=12.0d0, &
        Thsxy=-8.0d0, Thsyz=10.0d0, Thszx=7.0d0, &
        aThsx=0.8d0, aThsy=1.25d0, aThsz=0.7d0, aThsxy=0.5d0, aThsyz=0.6d0, aThszx=0.7d0

      ! temporary variables for physical constants
      DOUBLE PRECISION  :: mug, MW, Rg, ros, mus, Cpg, kg, Cps, ks, rog
  
      ! physical constant values
      mug= MU_g0
      rog= RO_g0
      MW= MW_AVG
      Rg= Gas_Const
      ros= RO_s0(1)
      mus= MU_s0
      Cpg= C_pg0
      kg= K_g0
      Cps= C_ps0
      ks= K_s0

      ! NSW
      ug0=10.0d0
      ugx=1.0d0
      vg0=10.0d0
      pgx=100.0d0

      ! FSW
      !ug0=10.0d0
      !vg0=5.0d0

      !!FLAGMMS
      !write(*,*) "MMS_Source: ros= ", ros
  
      SELECT CASE(idx)
      CASE(1)
      !pgscr -> no source terms for pressure correction equation !
        Write(*,*) "Wrong index in MMS_Source. Stop."
        Stop
      CASE(2)
      !ugsrc!
        MMS_Source = -(mug*ug0) + ug0**2*xt**3 - mug*ug0*Cos(2*Pi*xt)*Cos(2*Pi*yt) + &
       2*mug*Pi**2*ug0*xt**2*Cos(2*Pi*xt)*Cos(2*Pi*yt) + &
       ug0**2*xt**3*Cos(2*Pi*xt)*Cos(2*Pi*yt) + 2*pgx*xt*Cos(Pi*(xt + yt)) - &
       6*mug*Pi**2*ug0*xt**2*Cos(Pi*(xt + yt))**2 + &
       4*mug*Pi*ug0*xt*Cos(2*Pi*yt)*Sin(2*Pi*xt) - &
       Pi*ug0**2*xt**4*Cos(2*Pi*yt)*Sin(2*Pi*xt) + &
       4*mug*Pi*ug0*xt*Cos(2*Pi*xt)*Sin(2*Pi*yt) - &
       Pi*ug0**2*xt**4*Cos(2*Pi*xt)*Sin(2*Pi*yt) + mug*ug0*Sin(2*Pi*xt)*Sin(2*Pi*yt) - &
       2*mug*Pi**2*ug0*xt**2*Sin(2*Pi*xt)*Sin(2*Pi*yt) - &
       ug0**2*xt**3*Sin(2*Pi*xt)*Sin(2*Pi*yt) - pgx*Pi*xt**2*Sin(Pi*(xt + yt)) - &
       16*mug*Pi*ug0*xt*Cos(Pi*(xt + yt))*Sin(Pi*(xt + yt)) + &
       10*Pi*ug0**2*xt**3*Cos(Pi*(xt + yt))*Sin(Pi*(xt + yt)) + &
       4*Pi*ug0**2*xt**4*Cos(Pi*(xt + yt))*Sin(Pi*(xt + yt)) - &
       6*Pi*ug0**2*xt**3*yt*Cos(Pi*(xt + yt))*Sin(Pi*(xt + yt)) + &
       Pi*ug0**2*xt**4*Cos(2*Pi*xt)*Cos(2*Pi*yt)*Cos(Pi*(xt + yt))*Sin(Pi*(xt + yt)) + &
       ug0**2*xt**3*Cos(2*Pi*yt)*Cos(Pi*(xt + yt))*Sin(2*Pi*xt)*Sin(Pi*(xt + yt)) + &
       ug0**2*xt**3*Cos(2*Pi*xt)*Cos(Pi*(xt + yt))*Sin(2*Pi*yt)*Sin(Pi*(xt + yt)) - &
       Pi*ug0**2*xt**4*Cos(Pi*(xt + yt))*Sin(2*Pi*xt)*Sin(2*Pi*yt)*Sin(Pi*(xt + yt)) - &
       4*mug*ug0*Sin(Pi*(xt + yt))**2 + 6*mug*Pi**2*ug0*xt**2*Sin(Pi*(xt + yt))**2 + &
       5*ug0**2*xt**3*Sin(Pi*(xt + yt))**2 + &
       ug0**2*xt**3*Cos(2*Pi*xt)*Cos(2*Pi*yt)*Sin(Pi*(xt + yt))**2 - &
       Pi*ug0**2*xt**4*Cos(2*Pi*yt)*Sin(2*Pi*xt)*Sin(Pi*(xt + yt))**2 - &
       Pi*ug0**2*xt**4*Cos(2*Pi*xt)*Sin(2*Pi*yt)*Sin(Pi*(xt + yt))**2 - &
       ug0**2*xt**3*Sin(2*Pi*xt)*Sin(2*Pi*yt)*Sin(Pi*(xt + yt))**2 + &
       4*Pi*ug0**2*xt**4*Cos(Pi*(xt + yt))*Sin(Pi*(xt + yt))**3 + &
       4*ug0**2*xt**3*Sin(Pi*(xt + yt))**4
      CASE(3)
      !vgsrc!
        MMS_Source =  -15*ug0**2*xt**2 + 9*ug0**2*xt**2*yt - 3*mug*ug0*Cos(2*Pi*xt)*Cos(2*Pi*yt) + &
       6*mug*Pi**2*ug0*xt**2*Cos(2*Pi*xt)*Cos(2*Pi*yt) + &
       10*ug0**2*xt**2*Cos(2*Pi*xt)*Cos(2*Pi*yt) - &
       6*ug0**2*xt**2*yt*Cos(2*Pi*xt)*Cos(2*Pi*yt) + &
       ug0**2*xt**3*Cos(2*Pi*xt)**2*Cos(2*Pi*yt)**2 - &
       2*mug*Pi**2*ug0*xt**2*Cos(Pi*(xt + yt))**2 + &
       10*mug*Pi*ug0*xt*Cos(2*Pi*yt)*Sin(2*Pi*xt) - &
       (3*ug0**2*xt**2*Cos(2*Pi*yt)*Sin(2*Pi*xt))/(2.*Pi) - &
       10*Pi*ug0**2*xt**3*Cos(2*Pi*yt)*Sin(2*Pi*xt) - &
       Pi*ug0**2*xt**4*Cos(2*Pi*yt)*Sin(2*Pi*xt) + &
       6*Pi*ug0**2*xt**3*yt*Cos(2*Pi*yt)*Sin(2*Pi*xt) + &
       (ug0**2*xt**2*Cos(2*Pi*xt)*Cos(2*Pi*yt)**2*Sin(2*Pi*xt))/Pi - &
       Pi*ug0**2*xt**4*Cos(2*Pi*xt)*Cos(2*Pi*yt)**2*Sin(2*Pi*xt) - &
       ug0**2*xt**3*Cos(2*Pi*yt)**2*Sin(2*Pi*xt)**2 + &
       10*mug*Pi*ug0*xt*Cos(2*Pi*xt)*Sin(2*Pi*yt) - &
       (3*ug0**2*xt**2*Cos(2*Pi*xt)*Sin(2*Pi*yt))/(2.*Pi) - &
       10*Pi*ug0**2*xt**3*Cos(2*Pi*xt)*Sin(2*Pi*yt) - &
       Pi*ug0**2*xt**4*Cos(2*Pi*xt)*Sin(2*Pi*yt) + &
       6*Pi*ug0**2*xt**3*yt*Cos(2*Pi*xt)*Sin(2*Pi*yt) + &
       (ug0**2*xt**2*Cos(2*Pi*xt)**2*Cos(2*Pi*yt)*Sin(2*Pi*yt))/Pi - &
       Pi*ug0**2*xt**4*Cos(2*Pi*xt)**2*Cos(2*Pi*yt)*Sin(2*Pi*yt) + &
       3*mug*ug0*Sin(2*Pi*xt)*Sin(2*Pi*yt) - &
       6*mug*Pi**2*ug0*xt**2*Sin(2*Pi*xt)*Sin(2*Pi*yt) - &
       10*ug0**2*xt**2*Sin(2*Pi*xt)*Sin(2*Pi*yt) + &
       6*ug0**2*xt**2*yt*Sin(2*Pi*xt)*Sin(2*Pi*yt) - &
       4*ug0**2*xt**3*Cos(2*Pi*xt)*Cos(2*Pi*yt)*Sin(2*Pi*xt)*Sin(2*Pi*yt) - &
       (ug0**2*xt**2*Cos(2*Pi*yt)*Sin(2*Pi*xt)**2*Sin(2*Pi*yt))/Pi + &
       Pi*ug0**2*xt**4*Cos(2*Pi*yt)*Sin(2*Pi*xt)**2*Sin(2*Pi*yt) - &
       ug0**2*xt**3*Cos(2*Pi*xt)**2*Sin(2*Pi*yt)**2 - &
       (ug0**2*xt**2*Cos(2*Pi*xt)*Sin(2*Pi*xt)*Sin(2*Pi*yt)**2)/Pi + &
       Pi*ug0**2*xt**4*Cos(2*Pi*xt)*Sin(2*Pi*xt)*Sin(2*Pi*yt)**2 + &
       ug0**2*xt**3*Sin(2*Pi*xt)**2*Sin(2*Pi*yt)**2 - pgx*Pi*xt**2*Sin(Pi*(xt + yt)) - &
       4*mug*Pi*ug0*xt*Cos(Pi*(xt + yt))*Sin(Pi*(xt + yt)) + &
       10*Pi*ug0**2*xt**3*Cos(Pi*(xt + yt))*Sin(Pi*(xt + yt)) - &
       6*Pi*ug0**2*xt**3*yt*Cos(Pi*(xt + yt))*Sin(Pi*(xt + yt)) + &
       Pi*ug0**2*xt**4*Cos(2*Pi*xt)*Cos(2*Pi*yt)*Cos(Pi*(xt + yt))*Sin(Pi*(xt + yt)) + &
       ug0**2*xt**3*Cos(2*Pi*yt)*Cos(Pi*(xt + yt))*Sin(2*Pi*xt)*Sin(Pi*(xt + yt)) + &
       ug0**2*xt**3*Cos(2*Pi*xt)*Cos(Pi*(xt + yt))*Sin(2*Pi*yt)*Sin(Pi*(xt + yt)) - &
       Pi*ug0**2*xt**4*Cos(Pi*(xt + yt))*Sin(2*Pi*xt)*Sin(2*Pi*yt)*Sin(Pi*(xt + yt)) + &
       2*mug*Pi**2*ug0*xt**2*Sin(Pi*(xt + yt))**2 + &
       15*ug0**2*xt**2*Sin(Pi*(xt + yt))**2 - 9*ug0**2*xt**2*yt*Sin(Pi*(xt + yt))**2 + &
       3*ug0**2*xt**3*Cos(2*Pi*xt)*Cos(2*Pi*yt)*Sin(Pi*(xt + yt))**2 + &
       (3*ug0**2*xt**2*Cos(2*Pi*yt)*Sin(2*Pi*xt)*Sin(Pi*(xt + yt))**2)/(2.*Pi) - &
       Pi*ug0**2*xt**4*Cos(2*Pi*yt)*Sin(2*Pi*xt)*Sin(Pi*(xt + yt))**2 + &
       (3*ug0**2*xt**2*Cos(2*Pi*xt)*Sin(2*Pi*yt)*Sin(Pi*(xt + yt))**2)/(2.*Pi) - &
       Pi*ug0**2*xt**4*Cos(2*Pi*xt)*Sin(2*Pi*yt)*Sin(Pi*(xt + yt))**2 - &
       3*ug0**2*xt**3*Sin(2*Pi*xt)*Sin(2*Pi*yt)*Sin(Pi*(xt + yt))**2
      CASE(4)
      !wgsrc!
        MMS_Source = ZERO 
      CASE(5)
      !ussrc!
        MMS_Source = ZERO
      CASE(6)
      !vssrc!
        MMS_Source = ZERO
      CASE(7)
      !wssrc!
        MMS_Source = ZERO
      CASE(8)
      !Tgsrc!
        MMS_Source = ZERO
      CASE(9)
      !Tssrc!
        MMS_Source = ZERO
      CASE(10)
      !ropgsrc!
        ! continuity not solved for incompressible, const. volume fraction case
        MMS_Source = ZERO
      CASE(11)
      !ropssrc!
        ! continuity not solved for incompressible, const. volume fraction case
        MMS_Source = ZERO
      CASE(12)
      !Thssrc!
        MMS_Source = ZERO
      END SELECT

      End Function MMS_Source

      END MODULE MMS
