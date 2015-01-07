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
!      DOUBLE PRECISION :: pg0=0.0d0, pgx=0.2d0, pgy=-0.5d0, pgz=0.2d0, &
!        pgxy=-0.25d0, pgyz=-0.1d0, pgzx=0.1d0, &
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
  
      !!FLAGMMS
      !write(*,*) "MMS_Function: ros= ", ros
 
      !ug0=ZERO
      !vg0=ZERO
      !wg0=ZERO

      SELECT CASE(idx)
      CASE(1)
      !Pg!
        MMS_Function = pg0 + pgx*Cos(apgx*Pi*xt) + pgy*Cos(apgy*Pi*yt) + pgxy*Cos(apgxy*Pi*xt*yt) + &
         pgzx*Cos(apgzx*Pi*xt*zt) + pgz*Sin(apgz*Pi*zt) + pgyz*Sin(apgyz*Pi*yt*zt)
      CASE(2)
      !Ug!
        MMS_Function = awgy*Pi*wgy*Cos(awgy*Pi*yt) + awgxy*Pi*wgxy*xt*Cos(awgxy*Pi*xt*yt) - &
         avgyz*Pi*vgyz*yt*Cos(avgyz*Pi*yt*zt) + awgyz*Pi*wgyz*zt*Cos(awgyz*Pi*yt*zt) + &
         avgz*Pi*vgz*Sin(avgz*Pi*zt) + avgzx*Pi*vgzx*xt*Sin(avgzx*Pi*xt*zt)
      CASE(3)
      !Vg!
        MMS_Function = -(awgxy*Pi*wgxy*yt*Cos(awgxy*Pi*xt*yt)) + &
         augyz*Pi*ugyz*yt*Cos(augyz*Pi*yt*zt) + awgx*Pi*wgx*Sin(awgx*Pi*xt) - &
         augz*Pi*ugz*Sin(augz*Pi*zt) - augzx*Pi*ugzx*xt*Sin(augzx*Pi*xt*zt) + &
         awgzx*Pi*wgzx*zt*Sin(awgzx*Pi*xt*zt)
      CASE(4)
      !Wg!
        MMS_Function = avgx*Pi*vgx*Cos(avgx*Pi*xt) - augyz*Pi*ugyz*zt*Cos(augyz*Pi*yt*zt) + &
         augy*Pi*ugy*Sin(augy*Pi*yt) + augxy*Pi*ugxy*xt*Sin(augxy*Pi*xt*yt) - &
         avgxy*Pi*vgxy*yt*Sin(avgxy*Pi*xt*yt) - avgzx*Pi*vgzx*zt*Sin(avgzx*Pi*xt*zt)
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
!      DOUBLE PRECISION :: pg0=0.0d0, pgx=0.2d0, pgy=-0.5d0, pgz=0.2d0, &
!        pgxy=-0.25d0, pgyz=-0.1d0, pgzx=0.1d0, &
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
  
      !!FLAGMMS
      !write(*,*) "MMS_Source: ros= ", ros
  
      SELECT CASE(idx)
      CASE(1)
      !pgscr -> no source terms for pressure correction equation !
        Write(*,*) "Wrong index in MMS_Source. Stop."
        Stop
      CASE(2)
      !ugsrc!
        MMS_Source = -(mug*(awgxy**3*Pi**3*wgxy*xt*yt**2*Cos(awgxy*Pi*xt*yt) + &
              2*awgxy**2*Pi**2*wgxy*yt*Sin(awgxy*Pi*xt*yt))) - &
         mug*(-(awgxy**3*Pi**3*wgxy*xt*yt**2*Cos(awgxy*Pi*xt*yt)) + &
            2*avgzx**2*Pi**2*vgzx*zt*Cos(avgzx*Pi*xt*zt) - &
            2*awgxy**2*Pi**2*wgxy*yt*Sin(awgxy*Pi*xt*yt) - &
            avgzx**3*Pi**3*vgzx*xt*zt**2*Sin(avgzx*Pi*xt*zt)) - &
         mug*(-2*avgzx**2*Pi**2*vgzx*zt*Cos(avgzx*Pi*xt*zt) + &
            avgzx**3*Pi**3*vgzx*xt*zt**2*Sin(avgzx*Pi*xt*zt)) + &
         (-(apgx*pgx*Pi*Sin(apgx*Pi*xt)) - apgxy*pgxy*Pi*yt*Sin(apgxy*Pi*xt*yt) - &
            apgzx*pgzx*Pi*zt*Sin(apgzx*Pi*xt*zt))*&
          (1. - es0 - esx*Cos(aesx*Pi*xt) - esy*Cos(aesy*Pi*yt) - &
            esxy*Cos(aesxy*Pi*xt*yt) - eszx*Cos(aeszx*Pi*xt*zt) - &
            esz*Sin(aesz*Pi*zt) - esyz*Sin(aesyz*Pi*yt*zt)) - &
         mug*(-(awgy**3*Pi**3*wgy*Cos(awgy*Pi*yt)) - &
            awgxy**3*Pi**3*wgxy*xt**3*Cos(awgxy*Pi*xt*yt) + &
            avgyz**3*Pi**3*vgyz*yt**3*Cos(avgyz*Pi*yt*zt) - &
            awgyz**3*Pi**3*wgyz*zt**3*Cos(awgyz*Pi*yt*zt) + &
            awgxy*Pi*wgxy*(-(awgxy**2*Pi**2*xt*yt**2*Cos(awgxy*Pi*xt*yt)) - &
               2*awgxy*Pi*yt*Sin(awgxy*Pi*xt*yt)) - &
            avgz**3*Pi**3*vgz*Sin(avgz*Pi*zt) - &
            avgzx**3*Pi**3*vgzx*xt**3*Sin(avgzx*Pi*xt*zt) + &
            avgzx*Pi*vgzx*(2*avgzx*Pi*zt*Cos(avgzx*Pi*xt*zt) - &
               avgzx**2*Pi**2*xt*zt**2*Sin(avgzx*Pi*xt*zt)) - &
            avgyz*Pi*(-(avgyz**2*Pi**2*vgyz*yt*zt**2*Cos(avgyz*Pi*yt*zt)) - &
               2*avgyz*Pi*vgyz*zt*Sin(avgyz*Pi*yt*zt)) + &
            awgyz*Pi*wgyz*(-(awgyz**2*Pi**2*yt**2*zt*Cos(awgyz*Pi*yt*zt)) - &
               2*awgyz*Pi*yt*Sin(awgyz*Pi*yt*zt))) + &
         rog*(1. - es0 - esx*Cos(aesx*Pi*xt) - esy*Cos(aesy*Pi*yt) - &
            esxy*Cos(aesxy*Pi*xt*yt) - eszx*Cos(aeszx*Pi*xt*zt) - &
            esz*Sin(aesz*Pi*zt) - esyz*Sin(aesyz*Pi*yt*zt))*&
          (2*(awgxy*Pi*wgxy*Cos(awgxy*Pi*xt*yt) + &
               avgzx**2*Pi**2*vgzx*xt*zt*Cos(avgzx*Pi*xt*zt) - &
               awgxy**2*Pi**2*wgxy*xt*yt*Sin(awgxy*Pi*xt*yt) + &
               avgzx*Pi*vgzx*Sin(avgzx*Pi*xt*zt))*&
             (awgy*Pi*wgy*Cos(awgy*Pi*yt) + awgxy*Pi*wgxy*xt*Cos(awgxy*Pi*xt*yt) - &
               avgyz*Pi*vgyz*yt*Cos(avgyz*Pi*yt*zt) + &
               awgyz*Pi*wgyz*zt*Cos(awgyz*Pi*yt*zt) + avgz*Pi*vgz*Sin(avgz*Pi*zt) + &
               avgzx*Pi*vgzx*xt*Sin(avgzx*Pi*xt*zt)) + &
            (awgy*Pi*wgy*Cos(awgy*Pi*yt) + awgxy*Pi*wgxy*xt*Cos(awgxy*Pi*xt*yt) - &
               avgyz*Pi*vgyz*yt*Cos(avgyz*Pi*yt*zt) + &
               awgyz*Pi*wgyz*zt*Cos(awgyz*Pi*yt*zt) + avgz*Pi*vgz*Sin(avgz*Pi*zt) + &
               avgzx*Pi*vgzx*xt*Sin(avgzx*Pi*xt*zt))*&
             (-(awgxy*Pi*wgxy*Cos(awgxy*Pi*xt*yt)) + &
               augyz*Pi*ugyz*Cos(augyz*Pi*yt*zt) + &
               awgxy**2*Pi**2*wgxy*xt*yt*Sin(awgxy*Pi*xt*yt) - &
               augyz**2*Pi**2*ugyz*yt*zt*Sin(augyz*Pi*yt*zt)) + &
            (awgy*Pi*wgy*Cos(awgy*Pi*yt) + awgxy*Pi*wgxy*xt*Cos(awgxy*Pi*xt*yt) - &
               avgyz*Pi*vgyz*yt*Cos(avgyz*Pi*yt*zt) + &
               awgyz*Pi*wgyz*zt*Cos(awgyz*Pi*yt*zt) + avgz*Pi*vgz*Sin(avgz*Pi*zt) + &
               avgzx*Pi*vgzx*xt*Sin(avgzx*Pi*xt*zt))*&
             (-(avgzx**2*Pi**2*vgzx*xt*zt*Cos(avgzx*Pi*xt*zt)) - &
               augyz*Pi*ugyz*Cos(augyz*Pi*yt*zt) - &
               avgzx*Pi*vgzx*Sin(avgzx*Pi*xt*zt) + &
               augyz**2*Pi**2*ugyz*yt*zt*Sin(augyz*Pi*yt*zt)) + &
            (avgx*Pi*vgx*Cos(avgx*Pi*xt) - augyz*Pi*ugyz*zt*Cos(augyz*Pi*yt*zt) + &
               augy*Pi*ugy*Sin(augy*Pi*yt) + augxy*Pi*ugxy*xt*Sin(augxy*Pi*xt*yt) - &
               avgxy*Pi*vgxy*yt*Sin(avgxy*Pi*xt*yt) - &
               avgzx*Pi*vgzx*zt*Sin(avgzx*Pi*xt*zt))*&
             (avgz**2*Pi**2*vgz*Cos(avgz*Pi*zt) + &
               avgzx**2*Pi**2*vgzx*xt**2*Cos(avgzx*Pi*xt*zt) + &
               awgyz*Pi*wgyz*Cos(awgyz*Pi*yt*zt) + &
               avgyz**2*Pi**2*vgyz*yt**2*Sin(avgyz*Pi*yt*zt) - &
               awgyz**2*Pi**2*wgyz*yt*zt*Sin(awgyz*Pi*yt*zt)) + &
            (-(awgxy*Pi*wgxy*yt*Cos(awgxy*Pi*xt*yt)) + &
               augyz*Pi*ugyz*yt*Cos(augyz*Pi*yt*zt) + awgx*Pi*wgx*Sin(awgx*Pi*xt) - &
               augz*Pi*ugz*Sin(augz*Pi*zt) - augzx*Pi*ugzx*xt*Sin(augzx*Pi*xt*zt) + &
               awgzx*Pi*wgzx*zt*Sin(awgzx*Pi*xt*zt))*&
             (-(avgyz*Pi*vgyz*Cos(avgyz*Pi*yt*zt)) - &
               awgy**2*Pi**2*wgy*Sin(awgy*Pi*yt) - &
               awgxy**2*Pi**2*wgxy*xt**2*Sin(awgxy*Pi*xt*yt) + &
               avgyz**2*Pi**2*vgyz*yt*zt*Sin(avgyz*Pi*yt*zt) - &
               awgyz**2*Pi**2*wgyz*zt**2*Sin(awgyz*Pi*yt*zt)))
      CASE(3)
      !vgsrc!
        MMS_Source = -(mug*(-(awgxy**3*Pi**3*wgxy*xt**2*yt*Cos(awgxy*Pi*xt*yt)) - &
              2*awgxy**2*Pi**2*wgxy*xt*Sin(awgxy*Pi*xt*yt))) + &
         (apgyz*pgyz*Pi*zt*Cos(apgyz*Pi*yt*zt) - apgy*pgy*Pi*Sin(apgy*Pi*yt) - &
            apgxy*pgxy*Pi*xt*Sin(apgxy*Pi*xt*yt))*&
          (1. - es0 - esx*Cos(aesx*Pi*xt) - esy*Cos(aesy*Pi*yt) - &
            esxy*Cos(aesxy*Pi*xt*yt) - eszx*Cos(aeszx*Pi*xt*zt) - &
            esz*Sin(aesz*Pi*zt) - esyz*Sin(aesyz*Pi*yt*zt)) - &
         mug*(awgxy**3*Pi**3*wgxy*xt**2*yt*Cos(awgxy*Pi*xt*yt) - &
            augyz**3*Pi**3*ugyz*yt*zt**2*Cos(augyz*Pi*yt*zt) + &
            2*awgxy**2*Pi**2*wgxy*xt*Sin(awgxy*Pi*xt*yt) - &
            2*augyz**2*Pi**2*ugyz*zt*Sin(augyz*Pi*yt*zt)) - &
         mug*(augyz**3*Pi**3*ugyz*yt*zt**2*Cos(augyz*Pi*yt*zt) + &
            2*augyz**2*Pi**2*ugyz*zt*Sin(augyz*Pi*yt*zt)) - &
         mug*(awgxy**3*Pi**3*wgxy*yt**3*Cos(awgxy*Pi*xt*yt) - &
            augyz**3*Pi**3*ugyz*yt**3*Cos(augyz*Pi*yt*zt) - &
            awgx**3*Pi**3*wgx*Sin(awgx*Pi*xt) - &
            awgxy*Pi*(-(awgxy**2*Pi**2*wgxy*xt**2*yt*Cos(awgxy*Pi*xt*yt)) - &
               2*awgxy*Pi*wgxy*xt*Sin(awgxy*Pi*xt*yt)) + &
            augz**3*Pi**3*ugz*Sin(augz*Pi*zt) + &
            augzx**3*Pi**3*ugzx*xt**3*Sin(augzx*Pi*xt*zt) - &
            augzx*Pi*(2*augzx*Pi*ugzx*zt*Cos(augzx*Pi*xt*zt) - &
               augzx**2*Pi**2*ugzx*xt*zt**2*Sin(augzx*Pi*xt*zt)) - &
            awgzx**3*Pi**3*wgzx*zt**3*Sin(awgzx*Pi*xt*zt) + &
            awgzx*Pi*wgzx*(2*awgzx*Pi*xt*Cos(awgzx*Pi*xt*zt) - &
               awgzx**2*Pi**2*xt**2*zt*Sin(awgzx*Pi*xt*zt)) + &
            augyz*Pi*ugyz*(-(augyz**2*Pi**2*yt*zt**2*Cos(augyz*Pi*yt*zt)) - &
               2*augyz*Pi*zt*Sin(augyz*Pi*yt*zt))) + &
         rog*(1. - es0 - esx*Cos(aesx*Pi*xt) - esy*Cos(aesy*Pi*yt) - &
            esxy*Cos(aesxy*Pi*xt*yt) - eszx*Cos(aeszx*Pi*xt*zt) - &
            esz*Sin(aesz*Pi*zt) - esyz*Sin(aesyz*Pi*yt*zt))*&
          ((awgx**2*Pi**2*wgx*Cos(awgx*Pi*xt) - &
               augzx**2*Pi**2*ugzx*xt*zt*Cos(augzx*Pi*xt*zt) + &
               awgzx**2*Pi**2*wgzx*zt**2*Cos(awgzx*Pi*xt*zt) + &
               awgxy**2*Pi**2*wgxy*yt**2*Sin(awgxy*Pi*xt*yt) - &
               augzx*Pi*ugzx*Sin(augzx*Pi*xt*zt))*&
             (awgy*Pi*wgy*Cos(awgy*Pi*yt) + awgxy*Pi*wgxy*xt*Cos(awgxy*Pi*xt*yt) - &
               avgyz*Pi*vgyz*yt*Cos(avgyz*Pi*yt*zt) + &
               awgyz*Pi*wgyz*zt*Cos(awgyz*Pi*yt*zt) + avgz*Pi*vgz*Sin(avgz*Pi*zt) + &
               avgzx*Pi*vgzx*xt*Sin(avgzx*Pi*xt*zt)) + &
            (awgxy*Pi*wgxy*Cos(awgxy*Pi*xt*yt) + &
               avgzx**2*Pi**2*vgzx*xt*zt*Cos(avgzx*Pi*xt*zt) - &
               awgxy**2*Pi**2*wgxy*xt*yt*Sin(awgxy*Pi*xt*yt) + &
               avgzx*Pi*vgzx*Sin(avgzx*Pi*xt*zt))*&
             (-(awgxy*Pi*wgxy*yt*Cos(awgxy*Pi*xt*yt)) + &
               augyz*Pi*ugyz*yt*Cos(augyz*Pi*yt*zt) + awgx*Pi*wgx*Sin(awgx*Pi*xt) - &
               augz*Pi*ugz*Sin(augz*Pi*zt) - augzx*Pi*ugzx*xt*Sin(augzx*Pi*xt*zt) + &
               awgzx*Pi*wgzx*zt*Sin(awgzx*Pi*xt*zt)) + &
            (avgx*Pi*vgx*Cos(avgx*Pi*xt) - augyz*Pi*ugyz*zt*Cos(augyz*Pi*yt*zt) + &
               augy*Pi*ugy*Sin(augy*Pi*yt) + augxy*Pi*ugxy*xt*Sin(augxy*Pi*xt*yt) - &
               avgxy*Pi*vgxy*yt*Sin(avgxy*Pi*xt*yt) - &
               avgzx*Pi*vgzx*zt*Sin(avgzx*Pi*xt*zt))*&
             (-(augz**2*Pi**2*ugz*Cos(augz*Pi*zt)) - &
               augzx**2*Pi**2*ugzx*xt**2*Cos(augzx*Pi*xt*zt) + &
               awgzx**2*Pi**2*wgzx*xt*zt*Cos(awgzx*Pi*xt*zt) + &
               awgzx*Pi*wgzx*Sin(awgzx*Pi*xt*zt) - &
               augyz**2*Pi**2*ugyz*yt**2*Sin(augyz*Pi*yt*zt)) + &
            2*(-(awgxy*Pi*wgxy*yt*Cos(awgxy*Pi*xt*yt)) + &
               augyz*Pi*ugyz*yt*Cos(augyz*Pi*yt*zt) + awgx*Pi*wgx*Sin(awgx*Pi*xt) - &
               augz*Pi*ugz*Sin(augz*Pi*zt) - augzx*Pi*ugzx*xt*Sin(augzx*Pi*xt*zt) + &
               awgzx*Pi*wgzx*zt*Sin(awgzx*Pi*xt*zt))*&
             (-(awgxy*Pi*wgxy*Cos(awgxy*Pi*xt*yt)) + &
               augyz*Pi*ugyz*Cos(augyz*Pi*yt*zt) + &
               awgxy**2*Pi**2*wgxy*xt*yt*Sin(awgxy*Pi*xt*yt) - &
               augyz**2*Pi**2*ugyz*yt*zt*Sin(augyz*Pi*yt*zt)) + &
            (-(awgxy*Pi*wgxy*yt*Cos(awgxy*Pi*xt*yt)) + &
               augyz*Pi*ugyz*yt*Cos(augyz*Pi*yt*zt) + awgx*Pi*wgx*Sin(awgx*Pi*xt) - &
               augz*Pi*ugz*Sin(augz*Pi*zt) - augzx*Pi*ugzx*xt*Sin(augzx*Pi*xt*zt) + &
               awgzx*Pi*wgzx*zt*Sin(awgzx*Pi*xt*zt))*&
             (-(avgzx**2*Pi**2*vgzx*xt*zt*Cos(avgzx*Pi*xt*zt)) - &
               augyz*Pi*ugyz*Cos(augyz*Pi*yt*zt) - &
               avgzx*Pi*vgzx*Sin(avgzx*Pi*xt*zt) + &
               augyz**2*Pi**2*ugyz*yt*zt*Sin(augyz*Pi*yt*zt)))
      CASE(4)
      !wgsrc!
        MMS_Source = -(mug*(2*avgzx**2*Pi**2*vgzx*xt*Cos(avgzx*Pi*xt*zt) - &
              avgzx**3*Pi**3*vgzx*xt**2*zt*Sin(avgzx*Pi*xt*zt))) + &
         (apgz*pgz*Pi*Cos(apgz*Pi*zt) + apgyz*pgyz*Pi*yt*Cos(apgyz*Pi*yt*zt) - &
            apgzx*pgzx*Pi*xt*Sin(apgzx*Pi*xt*zt))*&
          (1. - es0 - esx*Cos(aesx*Pi*xt) - esy*Cos(aesy*Pi*yt) - &
            esxy*Cos(aesxy*Pi*xt*yt) - eszx*Cos(aeszx*Pi*xt*zt) - &
            esz*Sin(aesz*Pi*zt) - esyz*Sin(aesyz*Pi*yt*zt)) - &
         mug*(-(augyz**3*Pi**3*ugyz*yt**2*zt*Cos(augyz*Pi*yt*zt)) - &
            2*augyz**2*Pi**2*ugyz*yt*Sin(augyz*Pi*yt*zt)) - &
         mug*(-2*avgzx**2*Pi**2*vgzx*xt*Cos(avgzx*Pi*xt*zt) + &
            augyz**3*Pi**3*ugyz*yt**2*zt*Cos(augyz*Pi*yt*zt) + &
            avgzx**3*Pi**3*vgzx*xt**2*zt*Sin(avgzx*Pi*xt*zt) + &
            2*augyz**2*Pi**2*ugyz*yt*Sin(augyz*Pi*yt*zt)) - &
         mug*(-(avgx**3*Pi**3*vgx*Cos(avgx*Pi*xt)) + &
            augyz**3*Pi**3*ugyz*zt**3*Cos(augyz*Pi*yt*zt) - &
            augy**3*Pi**3*ugy*Sin(augy*Pi*yt) - &
            augxy**3*Pi**3*ugxy*xt**3*Sin(augxy*Pi*xt*yt) + &
            augxy*Pi*ugxy*(2*augxy*Pi*yt*Cos(augxy*Pi*xt*yt) - &
               augxy**2*Pi**2*xt*yt**2*Sin(augxy*Pi*xt*yt)) + &
            avgxy**3*Pi**3*vgxy*yt**3*Sin(avgxy*Pi*xt*yt) - &
            avgxy*Pi*(2*avgxy*Pi*vgxy*xt*Cos(avgxy*Pi*xt*yt) - &
               avgxy**2*Pi**2*vgxy*xt**2*yt*Sin(avgxy*Pi*xt*yt)) + &
            avgzx**3*Pi**3*vgzx*zt**3*Sin(avgzx*Pi*xt*zt) - &
            avgzx*Pi*(2*avgzx*Pi*vgzx*xt*Cos(avgzx*Pi*xt*zt) - &
               avgzx**2*Pi**2*vgzx*xt**2*zt*Sin(avgzx*Pi*xt*zt)) - &
            augyz*Pi*(-(augyz**2*Pi**2*ugyz*yt**2*zt*Cos(augyz*Pi*yt*zt)) - &
               2*augyz*Pi*ugyz*yt*Sin(augyz*Pi*yt*zt))) + &
         rog*(1. - es0 - esx*Cos(aesx*Pi*xt) - esy*Cos(aesy*Pi*yt) - &
            esxy*Cos(aesxy*Pi*xt*yt) - eszx*Cos(aeszx*Pi*xt*zt) - &
            esz*Sin(aesz*Pi*zt) - esyz*Sin(aesyz*Pi*yt*zt))*&
          ((augxy**2*Pi**2*ugxy*xt*yt*Cos(augxy*Pi*xt*yt) - &
               avgxy**2*Pi**2*vgxy*yt**2*Cos(avgxy*Pi*xt*yt) - &
               avgzx**2*Pi**2*vgzx*zt**2*Cos(avgzx*Pi*xt*zt) - &
               avgx**2*Pi**2*vgx*Sin(avgx*Pi*xt) + augxy*Pi*ugxy*Sin(augxy*Pi*xt*yt))*&
             (awgy*Pi*wgy*Cos(awgy*Pi*yt) + awgxy*Pi*wgxy*xt*Cos(awgxy*Pi*xt*yt) - &
               avgyz*Pi*vgyz*yt*Cos(avgyz*Pi*yt*zt) + &
               awgyz*Pi*wgyz*zt*Cos(awgyz*Pi*yt*zt) + avgz*Pi*vgz*Sin(avgz*Pi*zt) + &
               avgzx*Pi*vgzx*xt*Sin(avgzx*Pi*xt*zt)) + &
            (awgxy*Pi*wgxy*Cos(awgxy*Pi*xt*yt) + &
               avgzx**2*Pi**2*vgzx*xt*zt*Cos(avgzx*Pi*xt*zt) - &
               awgxy**2*Pi**2*wgxy*xt*yt*Sin(awgxy*Pi*xt*yt) + &
               avgzx*Pi*vgzx*Sin(avgzx*Pi*xt*zt))*&
             (avgx*Pi*vgx*Cos(avgx*Pi*xt) - augyz*Pi*ugyz*zt*Cos(augyz*Pi*yt*zt) + &
               augy*Pi*ugy*Sin(augy*Pi*yt) + augxy*Pi*ugxy*xt*Sin(augxy*Pi*xt*yt) - &
               avgxy*Pi*vgxy*yt*Sin(avgxy*Pi*xt*yt) - &
               avgzx*Pi*vgzx*zt*Sin(avgzx*Pi*xt*zt)) + &
            (avgx*Pi*vgx*Cos(avgx*Pi*xt) - augyz*Pi*ugyz*zt*Cos(augyz*Pi*yt*zt) + &
               augy*Pi*ugy*Sin(augy*Pi*yt) + augxy*Pi*ugxy*xt*Sin(augxy*Pi*xt*yt) - &
               avgxy*Pi*vgxy*yt*Sin(avgxy*Pi*xt*yt) - &
               avgzx*Pi*vgzx*zt*Sin(avgzx*Pi*xt*zt))*&
             (-(awgxy*Pi*wgxy*Cos(awgxy*Pi*xt*yt)) + &
               augyz*Pi*ugyz*Cos(augyz*Pi*yt*zt) + &
               awgxy**2*Pi**2*wgxy*xt*yt*Sin(awgxy*Pi*xt*yt) - &
               augyz**2*Pi**2*ugyz*yt*zt*Sin(augyz*Pi*yt*zt)) + &
            2*(avgx*Pi*vgx*Cos(avgx*Pi*xt) - augyz*Pi*ugyz*zt*Cos(augyz*Pi*yt*zt) + &
               augy*Pi*ugy*Sin(augy*Pi*yt) + augxy*Pi*ugxy*xt*Sin(augxy*Pi*xt*yt) - &
               avgxy*Pi*vgxy*yt*Sin(avgxy*Pi*xt*yt) - &
               avgzx*Pi*vgzx*zt*Sin(avgzx*Pi*xt*zt))*&
             (-(avgzx**2*Pi**2*vgzx*xt*zt*Cos(avgzx*Pi*xt*zt)) - &
               augyz*Pi*ugyz*Cos(augyz*Pi*yt*zt) - &
               avgzx*Pi*vgzx*Sin(avgzx*Pi*xt*zt) + &
               augyz**2*Pi**2*ugyz*yt*zt*Sin(augyz*Pi*yt*zt)) + &
            (-(awgxy*Pi*wgxy*yt*Cos(awgxy*Pi*xt*yt)) + &
               augyz*Pi*ugyz*yt*Cos(augyz*Pi*yt*zt) + awgx*Pi*wgx*Sin(awgx*Pi*xt) - &
               augz*Pi*ugz*Sin(augz*Pi*zt) - augzx*Pi*ugzx*xt*Sin(augzx*Pi*xt*zt) + &
               awgzx*Pi*wgzx*zt*Sin(awgzx*Pi*xt*zt))*&
             (augy**2*Pi**2*ugy*Cos(augy*Pi*yt) + &
               augxy**2*Pi**2*ugxy*xt**2*Cos(augxy*Pi*xt*yt) - &
               avgxy**2*Pi**2*vgxy*xt*yt*Cos(avgxy*Pi*xt*yt) - &
               avgxy*Pi*vgxy*Sin(avgxy*Pi*xt*yt) + &
               augyz**2*Pi**2*ugyz*zt**2*Sin(augyz*Pi*yt*zt)))
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
