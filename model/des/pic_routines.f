!  Author: R. Garg                                                     !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!


! contains the following subroutines:
!   mppic_compute_ps_grad,
!   mppic_bc_u_s, mppic_bc_v_s, mppic_bc_w_s
!   mppic_apply_ps_grad,
!   mppic_avg_eps
!   mppic_time_march


      SUBROUTINE PIC_TIME_MARCH
      USE param
      USE param1
      USE run
      USE output
      USE funits
      USE discretelement
      USE constant
      USE sendrecv
      USE des_bc
      USE cutcell
      USE mfix_pic
      USE pic_bc
      USE error_manager
      USE fldvar, only: P_g
      IMPLICIT NONE
!------------------------------------------------
! Local variables
!------------------------------------------------
      INTEGER NN, FACTOR, NP, IJK, I, J, K, BCV_I, LL, PC, TIME_LOOP_COUNT, IJK_BOT, IJK_TOP

!     Local variables to keep track of time when dem restart and des
!     write data need to be written when des_continuum_coupled is F
      DOUBLE PRECISION DES_RES_TIME, DES_SPX_TIME

!     Temporary variables when des_continuum_coupled is T to track
!     changes in solid time step
      DOUBLE PRECISION TMP_DTS, DTSOLID_TMP
      CHARACTER*5 FILENAME


!     Logical to see whether this is the first entry to this routine
      LOGICAL,SAVE:: FIRST_PASS = .TRUE.

!     temp
      DOUBLE PRECISION xpos,ypos, zpos, NORM_CF(3), DIST


      DOUBLE PRECISION :: DES_KE_VEC(DIMN)

      !time till which the PIC loop will be run
      double precision :: TEND_PIC_LOOP
      !number of PIC time steps
      Integer :: PIC_ITERS
! Identifies that the indicated particle is of interest for debugging
      LOGICAL FOCUS

      INCLUDE '../function.inc'
      INCLUDE '../fun_avg1.inc'
      INCLUDE '../fun_avg2.inc'

      CALL INIT_ERR_MSG("PIC_TIME_MARCH")

      S_TIME = TIME
      TIME_LOOP_COUNT = 0
      !compute the gas-phase pressure gradient at the beginning of the
      !des loop as the gas-phase pressure field will not change during
      !des calls
      IF(DES_CONTINUUM_COUPLED)   CALL COMPUTE_PG_GRAD

      TEND_PIC_LOOP = MERGE(TIME+DT, TSTOP, DES_CONTINUUM_COUPLED)
      PIC_ITERS = 0
      DO WHILE(S_TIME.LT.TEND_PIC_LOOP)
         ! If the current time in the discrete loop exceeds the current time in
         ! the continuum simulation, exit the lagrangian loop

         !DTPIC_MAX = MIN( 1e-04, DTPIC_MAX)
         DTSOLID = MERGE(MIN(DTPIC_MAX, DT), DTPIC_MAX, DES_CONTINUUM_COUPLED)
         DTSOLID_ORIG  = DTSOLID
         IF(MOD(PIC_ITERS, 10).eq.0) then
            IF(DES_CONTINUUM_COUPLED) then
               WRITE(ERR_MSG, 2000) DTSOLID, DTPIC_CFL, DTPIC_TAUP, DT

2000           FORMAT(/5x, &
                    & 'DTSOLID CURRENT  = ', g17.8, /5x,  &
                    & 'DTPIC_CFL        = ', g17.8, /5x,  &
                    & 'DTPIC TAUP       = ', g17.8, /5x, &
                    & 'DT FLOW          = ', g17.8)
            ELSE

               WRITE(ERR_MSG, 2001) S_TIME, DTSOLID, DTPIC_CFL, DTPIC_TAUP, DT

2001           FORMAT(/5x, &
                    & 'TIME             = ', g17.8, /5x,  &
                    & 'DTSOLID CURRENT  = ', g17.8, /5x,  &
                    & 'DTPIC_CFL        = ', g17.8, /5x,  &
                    & 'DTPIC TAUP       = ', g17.8, /5x, &
                    & 'DT FLOW          = ', g17.8)
            ENDIF
            CALL flush_err_msg(header = .false., footer = .false.)
         ENDIF

         PIC_ITERS  = PIC_ITERS + 1

         IF(S_TIME + DTSOLID.GT.TEND_PIC_LOOP) then
            ! If next time step in the discrete loop will exceed the current time
            ! in the continuum simulation, modify the discrete time step so final
            ! time will match
            WRITE(ERR_MSG, 2010)  DTSOLID, TIME + DT - S_TIME
            CALL flush_err_msg(header = .false., footer = .false.)

2010        FORMAT(/5X, &
                 & 'REDUCING DTSOLID TO ENSURE STIME + DTSOLID LE TIME + DT', /5x, &
                 & 'DTSOLID ORIG         = ', g17.8, /5x, &
                 & 'DTSOLID ACTUAL       = ', g17.8)

            DTSOLID = TEND_PIC_LOOP - S_TIME
         ENDIF

         CALL PARTICLES_IN_CELL

         CALL MPPIC_COMPUTE_PS_GRAD
         IF(DES_CONTINUUM_COUPLED)   CALL CALC_DES_DRAG_GS
         CALL CFUPDATEOLD

         CALL CFNEWVALUES

         ! Impose the wall-particle boundary condition for pic case
         ! The same routine also applies the mass inflow/outflow BCs as well
         CALL PIC_APPLY_WALLBC_STL


         ! Update time to reflect changes
         S_TIME = S_TIME + DTSOLID


         DTPIC_MAX = MIN(DTPIC_CFL, DTPIC_TAUP)

         ! When coupled, all write calls are made in time_march (the continuum
         ! portion) according to user settings for spx_time and res_time.
         ! The following section targets data writes for DEM only cases:
         IF(.NOT.DES_CONTINUUM_COUPLED) THEN
! Keep track of TIME for DEM simulations
            TIME = S_TIME

! Write data using des_spx_time and des_res_time; note the time will
! reflect current position of particles
            IF(PRINT_DES_DATA) THEN
               IF ( (S_TIME+0.1d0*DTSOLID >= DES_SPX_TIME) .OR. &
                    (S_TIME+0.1d0*DTSOLID >= TSTOP)) then
                  DES_SPX_TIME = &
                     ( INT((S_TIME+0.1d0*DTSOLID)/DES_SPX_DT) &
                     + 1 )*DES_SPX_DT
                  CALL WRITE_DES_DATA
                  IF(DMP_LOG) WRITE(UNIT_LOG,'(3X,A,X,ES15.5)') &
                     'DES data file written at time =', S_TIME
               ENDIF
            ENDIF

            IF ( (S_TIME+0.1d0*DTSOLID >= DES_RES_TIME) .OR. &
                 (S_TIME+0.1d0*DTSOLID >= TSTOP) .OR. &
                 (NN == FACTOR) ) THEN
               DES_RES_TIME = &
                  ( INT((S_TIME+0.1d0*DTSOLID)/DES_RES_DT) &
                  + 1 )*DES_RES_DT
                  CALL WRITE_RES0_DES
! Write RES1 here since it won't be called in time_march.  This will
! also keep track of TIME
               CALL WRITE_RES1
               IF(DMP_LOG) WRITE(UNIT_LOG,'(3X,A,X,ES15.5)') &
               'DES.RES and .RES files written at time =', S_TIME
            ENDIF
         ENDIF  ! end if (.not.des_continuum_coupled)




      ENDDO

      IF(DMP_LOG)       WRITE(UNIT_LOG,'(/5x, A, 2(2x, i10))') 'NUMBER OF TIMES MPPIC LOOP WAS CALLED AND PARTICLE COUNT = ', PIC_ITERS, PIP
      IF(mype.eq.pe_IO) WRITE(*,'(/5x, A, 2(2x, i10))') 'NUMBER OF TIMES MPPIC LOOP WAS CALLED AND PARTICLE COUNT = ', PIC_ITERS, PIP

      IJK_BOT = funijk(imin1, 2,kmin1)
      IJK_TOP = funijk(imin1, jmax1, kmin1)
      WRITE(*,'(/5X, A, 3(2x,g17.8))') 'MPPIC: PRES BOTTOM, TOP, AND DIFF KPA', P_G(IJK_BOT)/10000.d0, P_G(IJK_TOP)/10000.d0, (P_G(IJK_BOT) -  P_G(IJK_TOP))/10000.d0
      !WRITE(*,'(/5X, A, 3(2x,g17.8))') 'PRES BOTTOM, TOP, AND DIFF ', P_G(IJK_BOT), P_G(IJK_TOP), P_G(IJK_BOT) -  P_G(IJK_TOP)
      IF(.NOT.DES_CONTINUUM_COUPLED)then
         if(dmp_log)write(unit_log,'(1X,A)')&
         '<---------- END MPPIC_TIME_MARCH ----------'
!Pradeep call send recv for variables
      else
         !call send_recv(ep_g,2)
         !call send_recv(rop_g,2)
         !call send_recv(des_u_s,2)
         !call send_recv(des_v_s,2)
         !if(dimn.eq.3) call send_recv(des_w_s,2)
         !call send_recv(rop_s,2)
! now the above are communicated in comp_mean_fields_interp itself.
! so no need to communicate them here.
      END IF


      CALL FINL_ERR_MSG
    end SUBROUTINE PIC_TIME_MARCH



      SUBROUTINE MPPIC_COMP_EULERIAN_VELS_NON_CG
      USE param
      USE param1
      USE parallel
      USE constant
      USE geometry
      USE indices
      USE discretelement
      use desmpi
      !USE cutcell
      USE mfix_pic
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! local variable used for debugging
      LOGICAL :: FOCUS
! general i, j, k indices
      INTEGER :: I, J, K, IJK, IPJK, IJPK, IJKP, IMJK, IJMK, IJKM,&
                 IPJPK, IPJKP, IJPKP, IPJPKP, II, JJ, KK

! index of solid phase that particle NP belongs to
      INTEGER :: M
!-----------------------------------------------
! Include statement functions
!-----------------------------------------------
      INCLUDE '../function.inc'


      DO IJK = ijkstart3, ijkend3
            I = I_OF(IJK)
            J = J_OF(IJK)
            K = K_OF(IJK)
            !U_so(IJK, :) = U_s(IJK, :)
            !V_so(IJK, :) = V_s(IJK, :)
            !W_so(IJK, :) = W_s(IJK, :)
            PIC_U_S(IJK, :) = ZERO
            PIC_V_S(IJK, :) = ZERO
            PIC_W_S(IJK, :) = ZERO
            IF (.NOT.IS_ON_myPE_plus1layer(I,J,K)) CYCLE
            IF(.NOT.FLUID_AT(IJK)) CYCLE

            IF(I.GE.IMIN1.AND.I.LT.IMAX1) then !because we don't want to
               !calculate solids velocity at the wall cells.
               IPJK = IP_OF(IJK)
               DO M = 1, DES_MMAX
                  PIC_U_S(IJK,M) = 0.5d0*(DES_U_S(IJK, M) + DES_U_S(IPJK,M))
               ENDDO
            ENDIF

            if(J.GE.JMIN1.AND.J.LT.JMAX1) then !because we don't want to
               !calculate solids velocity at the wall cells.
               IJPK = JP_OF(IJK)
               DO M = 1, DES_MMAX
                  PIC_V_S(IJK,M) = 0.5d0*(DES_V_S(IJK, M) + DES_V_S(IJPK,M))
               ENDDO
            ENDIF


            if(K.GE.KMIN1.AND.K.LT.KMAX1.AND.DO_K) then !because we don't want to
               !calculate solids velocity at the wall cells.
               IJKP = KP_OF(IJK)
               DO M = 1, DES_MMAX
                  PIC_W_S(IJK,M) = 0.5d0*(DES_W_S(IJK, M) + DES_W_S(IJKP,M))
               ENDDO
            ENDIF
         ENDDO
         !CALL SET_WALL_BC(IER)
         !the above routine will apply noslip or free slip BC as per the mfix  convention.
         !currently, this implies NSW or FSW wall BC's will be re-applied to gas-phase
         !field as well. This can be changed later on to be more specific to MPPIC case
         !CALL WRITE_MPPIC_VEL_S
         CALL MPPIC_BC_U_S
         CALL MPPIC_BC_V_S
         IF(DO_K) CALL MPPIC_BC_W_S

      END SUBROUTINE MPPIC_COMP_EULERIAN_VELS_NON_CG

      SUBROUTINE MPPIC_COMP_EULERIAN_VELS_CG
      USE param
      USE param1
      USE parallel
      USE constant
      USE geometry
      USE indices
      USE discretelement
      use desmpi
      USE cutcell
      USE mfix_pic
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! local variable used for debugging
      LOGICAL :: FOCUS
! general i, j, k indices
      INTEGER :: I, J, K, IJK, IPJK, IJPK, IJKP, IMJK, IJMK, IJKM,&
                 IPJPK, IPJKP, IJPKP, IPJPKP, II, JJ, KK

! index of solid phase that particle NP belongs to
      INTEGER :: M
!-----------------------------------------------
! Include statement functions
!-----------------------------------------------
      INCLUDE '../function.inc'

      DO IJK = ijkstart3, ijkend3
         I = I_OF(IJK)
         J = J_OF(IJK)
         K = K_OF(IJK)
         !U_so(IJK, :) = U_s(IJK, :)
         !V_so(IJK, :) = V_s(IJK, :)
         !W_so(IJK, :) = W_s(IJK, :)
         PIC_U_S(IJK, :) = ZERO
         PIC_V_S(IJK, :) = ZERO
         PIC_W_S(IJK, :) = ZERO

         IF (.NOT.IS_ON_myPE_plus1layer(I,J,K)) CYCLE

         IF(WALL_U_AT(IJK)) THEN
            PIC_U_S(IJK, :) = ZERO
            !currently only No slip BC is being set on this mean
            !solid's velocity field. Later this wall part can be
            !treated separately and U_S set only for scalar cells
            !where FLUID_AT(IJK) is true.
         ELSE
            if(.not.FLUID_AT(IJK)) cycle

            IPJK = IP_OF(IJK)
            IF(FLUID_AT(IPJK)) THEN
               DO M = 1, DES_MMAX
                  PIC_U_S(IJK,M) = 0.5d0*(DES_U_S(IJK, M) + DES_U_S(IPJK,M))
               ENDDO
            ELSE
               PIC_U_S(IJK,:) = DES_U_S(IJK, :)
            ENDIF
         ENDIF

         IF(WALL_V_AT(IJK)) THEN
            PIC_V_S(IJK, :) = ZERO
         ELSE
            if(.not.FLUID_AT(IJK)) cycle
            IJPK = JP_OF(IJK)
            IF(FLUID_AT(IJPK)) THEN
               DO M = 1, DES_MMAX
                  PIC_V_S(IJK,M) = 0.5d0*(DES_V_S(IJK, M) + DES_V_S(IJPK,M))
               ENDDO
            ELSE
               PIC_V_S(IJK,:) = DES_V_S(IJK, :)
            ENDIF
         ENDIF

         IF(DO_K) THEN
            IF(WALL_W_AT(IJK)) THEN
               PIC_W_S(IJK, :) = ZERO
            ELSE
               if(.not.FLUID_AT(IJK)) cycle
               IJKP = KP_OF(IJK)
               IF(FLUID_AT(IJKP)) THEN
                  DO M = 1, DES_MMAX
                     PIC_W_S(IJK,M) = 0.5d0*(DES_W_S(IJK, M) + DES_W_S(IJKP,M))
                  ENDDO
               ELSE
                  PIC_W_S(IJK,:) = DES_W_S(IJK, :)
               ENDIF
            ENDIF
         ENDIF
      ENDDO

      END SUBROUTINE MPPIC_COMP_EULERIAN_VELS_CG

      SUBROUTINE MPPIC_APPLY_PS_GRAD_PART(L)

      USE param
      USE param1
      USE parallel
      USE constant
      USE discretelement
      USE mpi_utility
      USE mfix_pic
      USE cutcell
      USE fldvar, only: ep_g
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: L
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
      INTEGER M, IDIM
      INTEGER I, J, K, IJK, IJK_C

      DOUBLE PRECISION D(DIMN), DIST,  DP_BAR, COEFF_EN, COEFF_EN2

      DOUBLE PRECISION DELUP(DIMN), UPRIMETAU(DIMN), UPRIMETAU_INT(DIMN), PS_FORCE(DIMN), VEL_ORIG(DIMN)
! index to track accounted for particles
      INTEGER PC , epg_min_loc(1)

! Logical for local debug warnings
      LOGICAL DES_LOC_DEBUG

! maximum distance particles can move in MPPIC
      DOUBLE PRECISION MAXDIST_PIC, UPRIMEMOD, UPRIMEMODNEW, signvel

! dt's in each direction  based on cfl_pic for the mppic case

      DOUBLE PRECISION  MEANUS(DIMN, DES_MMAX), RELVEL(DIMN), MEANVEL(DIMN), VEL_NEW(DIMN)
!      INTEGER :: TOT_CASE, case1_count, case2_count, case3_count, case4_count

      LOGICAL :: INSIDE_DOMAIN
!-----------------------------------------------
! Functions
!-----------------------------------------------
      DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT

!-----------------------------------------------

      INCLUDE '../function.inc'
      INCLUDE '../fun_avg1.inc'
      INCLUDE '../fun_avg2.inc'

      M = PIJK(L,5)
      IJK = PIJK(L,4)
      COEFF_EN  = MPPIC_COEFF_EN1
      COEFF_EN2 = MPPIC_COEFF_EN2

      VEL_ORIG(1:DIMN) = DES_VEL_NEW(:,L)
      VEL_NEW (1:DIMN) = DES_VEL_NEW(:,L)
      !IF(L.eq.1) WRITE(*,*) 'MPPIC COEFFS = ', COEFF_EN, COEFF_EN2
      IF(L.EQ.FOCUS_PARTICLE) THEN

         WRITE(*,'(A20,2x,3(2x,i4))') 'CELL ID = ', PIJK(L,1:3)
         WRITE(*,'(A20,2x,3(2x,g17.8))') 'EPS = ', 1.d0 - EP_g(PIJK(L,4))
         WRITE(*,'(A20,2x,3(2x,g17.8))') 'DES_VEL ORIG = ', DES_VEL_NEW(:,L)

         WRITE(*,'(A20,2x,3(2x,g17.8))') 'FC = ', FC(:,L)
      ENDIF

      MEANVEL(1) = DES_U_S(IJK,M)
      MEANVEL(2) = DES_V_S(IJK,M)
      IF(DO_K) MEANVEL(3) = DES_W_S(IJK,M)

      PS_FORCE(:) = PS_GRAD(L, :)
      !IF(ABS(PS_FORCE(2)).GT.ZERO)  WRITE(*,*) 'PS_FORCE = ', PS_FORCE
      DELUP(:) = -PS_FORCE(:)

      MEANUS(:,M) =  AVGSOLVEL_P (:,L)
      !MEANUS(:,M) = MEANVEL(:)
      RELVEL(:) = DES_VEL_NEW(:,L) - MEANUS(:,M)

      !IF(EPg_P(L).gt.1.2d0*ep_star) RETURN

      DO IDIM = 1, DIMN

         IF(ABS(PS_FORCE(IDIM)).eq.zero) cycle

         IF(VEL_ORIG(IDIM)*MEANUS(IDIM,M).GT.ZERO) THEN

            IF(VEL_ORIG(IDIM)*DELUP(IDIM).GT.ZERO) THEN

               IF(ABS(MEANUS(IDIM,M)) .GT. ABS(VEL_ORIG(IDIM))) THEN
                  IJK_C = IJK
                  IF(IDIM.eq.1) then
                     if(VEL_ORIG(IDIM).LT.ZERO) IJK_C = IM_OF(IJK)
                     if(VEL_ORIG(IDIM).GT.ZERO) IJK_C = IP_OF(IJK)
                  ELSEIF(IDIM.eq.2) then
                     if(VEL_ORIG(IDIM).LT.ZERO) IJK_C = JM_OF(IJK)
                     if(VEL_ORIG(IDIM).GT.ZERO) IJK_C = JP_OF(IJK)
                  ELSEIF(IDIM.eq.3) then
                     if(VEL_ORIG(IDIM).LT.ZERO) IJK_C = KM_OF(IJK)
                     if(VEL_ORIG(IDIM).GT.ZERO) IJK_C = KP_OF(IJK)
                  ENDIF
                  INSIDE_DOMAIN = .false.
                  INSIDE_DOMAIN = fluid_at(IJK_C)!and.(.not.cut_cell_at(IJK_C))

                  if(INSIDE_DOMAIN) then
                     VEL_NEW(IDIM) = MEANUS(IDIM,M) - COEFF_EN*RELVEL(IDIM)
                  endif
!                  case4_count = case4_count + 1
               ELSE
                  !do nothing
               ENDIF
            ELSE
               IF(ABS(VEL_ORIG(IDIM)).GT.ABS(MEANUS(IDIM,M))) then
                  !VEL_NEW(IDIM) = MEANUS(IDIM,M) - COEFF_EN*RELVEL(IDIM)
                  IJK_C = IJK
                  IF(IDIM.eq.1) then
                     if(VEL_ORIG(IDIM).LT.ZERO) IJK_C = IM_OF(IJK)
                     if(VEL_ORIG(IDIM).GT.ZERO) IJK_C = IP_OF(IJK)
                  ELSEIF(IDIM.eq.2) then
                     if(VEL_ORIG(IDIM).LT.ZERO) IJK_C = JM_OF(IJK)
                     if(VEL_ORIG(IDIM).GT.ZERO) IJK_C = JP_OF(IJK)
                  ELSEIF(IDIM.eq.3) then
                     if(VEL_ORIG(IDIM).LT.ZERO) IJK_C = KM_OF(IJK)
                     if(VEL_ORIG(IDIM).GT.ZERO) IJK_C = KP_OF(IJK)
                  ENDIF
                  !VEL_NEW(IDIM) = (MEANUS(IDIM,M) - COEFF_EN*RELVEL(IDIM))

                  INSIDE_DOMAIN = .false.
                  INSIDE_DOMAIN = fluid_at(IJK_C)!.and.(.not.cut_cell_at(IJK_C))

                  if(INSIDE_DOMAIN) then
                     VEL_NEW(IDIM) = (MEANUS(IDIM,M) - COEFF_EN*RELVEL(IDIM))
                  else
                     VEL_NEW(IDIM) = -COEFF_EN*VEL_ORIG(IDIM)
                  endif

                  IJK_C = IJK
                  IF(IDIM.eq.1) then
                     if(VEL_NEW(IDIM).LT.ZERO) IJK_C = IM_OF(IJK)
                     if(VEL_NEW(IDIM).GT.ZERO) IJK_C = IP_OF(IJK)
                  ELSEIF(IDIM.eq.2) then
                     if(VEL_NEW(IDIM).LT.ZERO) IJK_C = JM_OF(IJK)
                     if(VEL_NEW(IDIM).GT.ZERO) IJK_C = JP_OF(IJK)
                  ELSEIF(IDIM.eq.3) then
                     if(VEL_NEW(IDIM).LT.ZERO) IJK_C = KM_OF(IJK)
                     if(VEL_NEW(IDIM).GT.ZERO) IJK_C = KP_OF(IJK)
                  ENDIF
                  INSIDE_DOMAIN = .false.
                  INSIDE_DOMAIN = fluid_at(IJK_C) !.and.(.not.cut_cell_at(IJK_C))

                  !if(.not.INSIDE_DOMAIN) then
                  !   VEL_NEW(IDIM) = -COEFF_EN*VEL_ORIG(IDIM)
                  !ENDIF

!                  case1_count = case1_count + 1
               ELSE
                  !do nothing
                  VEL_NEW(IDIM) = COEFF_EN2 * VEL_ORIG(IDIM)
                  !turning on the above would make the model uncondtionally stable
!                  case1_count = case1_count + 1

               ENDIF
            ENDIF
         ELSE
            IF(MEANUS(IDIM,M)*DELUP(IDIM).GT.ZERO) THEN
               VEL_NEW(IDIM) = MEANUS(IDIM,M) - COEFF_EN*RELVEL(IDIM)
!               case2_count = case2_count + 1
            ELSE
!               case3_count = case3_count + 1
               !DO NOTHING
            ENDIF
         ENDIF

         IF(MPPIC_GRAV_TREATMENT) THEN
            IF(DELUP(IDIM)*GRAV(IDIM).LT.ZERO.AND.VEL_ORIG(IDIM)*GRAV(IDIM).GT.ZERO) THEN
               VEL_NEW(IDIM) = -COEFF_EN*VEL_ORIG(IDIM)
               !MPPIC_VPTAU(L, IDIM) = VEL_NEW(IDIM) -  DES_VEL_NEW(IDIM, L)
            ENDIF
         ENDIF
         !MPPIC_VPTAU(L, IDIM) = VEL_NEW(IDIM) -  DES_VEL_NEW(IDIM, L)
         DES_VEL_NEW(IDIM, L) = VEL_NEW(IDIM)
      ENDDO

         !
      if(L.eq.FOCUS_PARTICLE) THEN
         !iF((IJK.eq.epg_min_loc(1).or.IJK_OLD.eq.epg_min_loc(1)).and.epg_min2.lt.0.38) then
         !if(j.ne.2) cycle

         WRITE(*,'(A20,2x,3(2x,i5))') 'PIJK I, J, K =', I_OF(IJK),J_OF(IJK),K_OF(IJK)
         WRITE(*,'(A20,2x,3(2x,g17.8))') 'DES_VEL ORIG = ', VEL_ORIG(:)
         WRITE(*,'(A20,2x,3(2x,g17.8))') 'DES_VEL NEW = ', DES_VEL_NEW(:,L)
         WRITE(*,'(A20,2x,3(2x,g17.8))') 'MEANUS = ', MEANUS(:,1)
         WRITE(*,'(A20,2x,3(2x,g17.8))') 'DES_POS_NEW = ', DES_POS_NEW(:,L)
         WRITE(*,'(A20,2x,3(2x,g17.8))') 'GRAD PS = ', PS_FORCE(:)
         WRITE(*,'(A20,2x,3(2x,g17.8))') 'DELUP =  ', DELUP(:)

         !WRITE(*,'(A20,2x,3(2x,g17.8))') 'UPRIMETAU_INT = ', UPRIMETAU_INT(:)
         !WRITE(*,'(A20,2x,3(2x,g17.8))') 'UPRIMETAU = ', UPRIMETAU(:)
         !IF(UPRIMEMOD*DTSOLID.GT.MEAN_FREE_PATH) WRITE(*,'(A20,2x,3(2x,g17.8))') 'U*DT, MFP =', UPRIMEMOD*DTSOLID, MEAN_FREE_PATH
         read(*,*)
      ENDIF

      !TOT_CASE = case1_count + case2_count + case3_count + case4_count
      !IF(TOT_CASE.GT.0) THEN
      !WRITE(*,'(A, 4(2x,i10))') 'CASE COUNT NUMBERS  = ', case1_count ,case2_count ,case3_count ,case4_count
      !WRITE(*,'(A, 4(2x,g12.7))') 'CASE COUNT %AGE = ', real(case1_count)*100./real(tot_case),real(case2_count)*100./real(tot_case), real(case3_count)*100./real(tot_case), real(case4_count)*100./real(tot_case)
      !ENDIF
      RETURN

      !MEANUS(:,M) = MEANVEL(:)
      !RELVEL(:) = DES_VEL_NEW(:,L) - MEANUS(:,M)
      !DO IDIM = 1, DIMN
      !    IF(ABS(PS_FORCE(IDIM)).eq.zero) cycle

      !   IF(RELVEL(IDIM)*DELUP(IDIM).GT.ZERO) THEN
      !do nothing
      !    ELSE
      !       DES_VEL_NEW(IDIM,L) = MEANUS(IDIM,M) - 0.4d0*RELVEL(IDIM)

      !IF(DES_VEL_NEW(IDIM,L)*DELUP(IDIM).LT.ZERO) DES_VEL_NEW(IDIM,L) = -0.5d0*DES_VEL_NEW(IDIM,L)
      !    ENDIF
      ! ENDDO
      ! CYCLE


      END SUBROUTINE MPPIC_APPLY_PS_GRAD_PART


      SUBROUTINE MPPIC_COMPUTE_PS_GRAD
      USE param
      USE param1
      USE parallel
      USE physprop
      USE fldvar, only: ep_g, u_g, v_g, w_g
      USE run
      USE geometry
      USE indices
      USE bc
      USE compar
      USE sendrecv
      USE discretelement
      USE constant
      USE cutcell
      USE interpolation
      USE mfix_pic
      implicit none

      ! general i, j, k indices
      INTEGER I, J, K, IJK, IPJK, IJPK, IJKP, IMJK, IJMK, IJKM, IDIM, M

      ! temporary variables used to calculate pressure at scalar cell edge
      DOUBLE PRECISION TEMP1, TEMP2, avg_factor, VOL_TOT_VEC(3), VOL_TOT_SCAL

      integer :: korder, iw,ie,js,jn,kb,ktp, onew, pcell(3), cur_ijk, NP, nindx

      integer :: ii,jj,kk, ipjpk, ijpkp, ipjkp, ipjpkp, I1, J1, K1

      double precision :: vol_ijk, vol_ipjk, vol_ijpk, vol_ipjpk
      double precision :: vol_ijkp, vol_ipjkp, vol_ijpkp, vol_ipjpkp

      INCLUDE '../function.inc'
      INCLUDE '../fun_avg1.inc'
      INCLUDE '../fun_avg2.inc'

      if(MPPIC_SOLID_STRESS_SNIDER) then

         DO IJK = IJKSTART3, IJKEND3
            IF(FLUID_AT(IJK)) THEN
               PIC_P_S(IJK,1) = PSFAC_FRIC_PIC * ((1.d0 - EP_G(IJK))**FRIC_EXP_PIC)/ MAX(EP_G(IJK) - EP_STAR, FRIC_NON_SING_FAC*EP_G(IJK))
               !write(102,'(2(2x,i4),2(2x,g17.8))') I_OF(IJK), J_OF(IJK), EP_S(IJK,1), P_STAR(IJK)
            ELSE
               !So that ghost cells have higher pressure
               PIC_P_S(IJK,1) = PSFAC_FRIC_PIC * ((1.d0 - 0.1d0)**FRIC_EXP_PIC)/ MAX(0.1d0 - EP_STAR, FRIC_NON_SING_FAC*0.1d0)

            ENDIF
         ENDDO


      ELSE
         DO IJK = IJKSTART3, IJKEND3
            PS_FORCE_PIC(IJK,:) = ZERO

            IF(FLUID_AT(IJK)) THEN
               if(EP_G(IJK).lt.ep_star) then
                  PIC_P_S(IJK,1) = one*(one-ep_g(ijk))
               else

                  PIC_P_S(IJK,1) = ZERO
               endif

            ELSE
               PIC_P_S(IJK,1) = 1.!\*0.d0
            ENDIF
         ENDDO

      ENDIF

      CALL SEND_RECV(PIC_P_S,1)

      !Since EP_G is already shared across the processors, the pressure gradient calculation
      !can be made a function call so that the extra communication of P_S can be avoided.

      !DO k = kstart2, kend1
      !do j = jstart2, jend1
      !do i = istart2, iend1
      DO IJK = IJKSTART3, IJKEND3
         PS_FORCE_PIC(IJK,:) = ZERO
         IF(.NOT.FLUID_AT(IJK)) CYCLE

         I = I_OF(IJK)
         J = J_OF(IJK)
         K = K_OF(IJK)
         IPJK = IP_OF(IJK)
         IJPK = JP_OF(IJK)
         IJKP = KP_OF(IJK)

         if(fluid_at(ipjk)) then
            PS_FORCE_PIC(IJK,1) = 2.d0*(PIC_P_S(IPJK,1) - PIC_P_S(IJK,1))/(DX(I)+DX(I_of(ipjk)))
         else
            if(PIC_P_S(IJK,1).GT.ZERO) then
               !this will always be true for Snider's case
               PS_FORCE_PIC(IJK,1) = 2.d0*(PIC_P_S(IPJK,1) - PIC_P_S(IJK,1))/(DX(I)+DX(I_of(ipjk)))
            ELSE
               PS_FORCE_PIC(IJK,1) = zero
            endif
         ENDIF

         IF(FLUID_AT(IJPK)) THEN
            PS_FORCE_PIC(IJK,2) = 2.d0*(PIC_P_S(IJPK,1) - PIC_P_S(IJK,1))/(DY(j)+Dy(j_of(ijpk)))
         ELSE
            IF(PIC_P_S(IJK,1).GT.ZERO) then
               PS_FORCE_PIC(IJK,2) = 2.d0*(PIC_P_S(IJPK,1) - PIC_P_S(IJK,1))/(DY(j)+Dy(j_of(ijpk)))
            ELSE
               PS_FORCE_PIC(IJK,2) = zero
            ENDIF
         ENDIF

         IF(DO_K) THEN
            IF(FLUID_AT(IJKP)) then
               PS_FORCE_PIC(IJK,3) = 2.d0*(PIC_P_S(IJKP,1) - PIC_P_S(IJK,1))/(Dz(k)+Dz(k_of(ijkp)))
            ELSE
               IF(PIC_P_S(IJK,1).GT.ZERO) then
                  PS_FORCE_PIC(IJK,3) = 2.d0*(PIC_P_S(IJKP,1) - PIC_P_S(IJK,1))/(Dz(k)+Dz(k_of(ijkp)))
               ELSE
                  PS_FORCE_PIC(IJK,3) = ZERO
               ENDIF
            ENDIF
         ENDIF
      ENDDO
      !Rahul:
      !the above will not compute pressure gradients normal to the east, south and bottom faces
      !which are very important
      I1 = 1
      DO K1 = KSTART3, KEND3
         DO J1 = JSTART3, JEND3
!//
            IF(I1.NE.ISTART2)   EXIT
            IF (.NOT.IS_ON_myPE_owns(I1,J1,K1)) CYCLE
            IJK = FUNIJK(I1,J1,K1)
            I = I_OF(IJK)
            J = J_OF(IJK)
            K = K_OF(IJK)

            IPJK = IP_OF(IJK)
            IF(PIC_P_S(IPJK,1).GT.ZERO) then
               PS_FORCE_PIC(IJK,1) = 2.d0*(PIC_P_S(IPJK,1) - PIC_P_S(IJK,1))/(DX(I)+DX(I_of(ipjk)))
            else
               PS_FORCE_PIC(IJK,1) = zero
            endif

         ENDDO
      ENDDO
      J1 = 1
      DO K1 = KSTART3, KEND3
         DO I1 = ISTART3, IEND3
!//
            IF(J1.NE.JSTART2)   EXIT
            IF (.NOT.IS_ON_myPE_owns(I1,J1,K1)) CYCLE
            IJK = FUNIJK(I1,J1,K1)
            I = I_OF(IJK)
            J = J_OF(IJK)
            K = K_OF(IJK)


            IJPK = JP_OF(IJK)

            IF(PIC_P_S(IJPK,1).GT.ZERO) then

               PS_FORCE_PIC(IJK,2) = 2.d0*(PIC_P_S(IJPK,1) - PIC_P_S(IJK,1))/(DY(j)+Dy(j_of(ijpk)))
            ELSE
               PS_FORCE_PIC(IJK,2) = ZERO
            ENDIF

         END DO
      END DO

      IF(DO_K) then
         K1 = 1
         DO J1 = JSTART3, JEND3
            DO I1 = ISTART3, IEND3
               IF(K1.NE.KSTART2)   EXIT
               IF (.NOT.IS_ON_myPE_owns(I1,J1,K1)) CYCLE
               IJK = FUNIJK(I1,J1,K1)
               I = I_OF(IJK)
               J = J_OF(IJK)
               K = K_OF(IJK)


               IJKP = KP_OF(IJK)

               IF(PIC_P_S(IJKP,1).GT.ZERO) then

                  PS_FORCE_PIC(IJK,3) = 2.d0*(PIC_P_S(IJKP,1) - PIC_P_S(IJK,1))/(Dz(k)+Dz(k_of(ijkp)))
               ELSE
                  PS_FORCE_PIC(IJK, 3) = ZERO
               ENDIF

            END DO
         END DO
      ENDIF

      DO IDIM = 1, merge(2,3,NO_K)
         CALL SEND_RECV(PS_FORCE_PIC(:,IDIM),1)
      ENDDO

      CALL SET_INTERPOLATION_SCHEME(2)

      KORDER = merge ( 1, 2, no_k) !1+(DIMN-2)

! avg_factor=0.25 (in 3D) or =0.5 (in 2D)
      !avg_factor = 0.25d0*(dimn-2) + 0.5d0*(3-dimn)
      AVG_FACTOR = merge(0.5d0, 0.25D0, NO_K)

      do ijk = ijkstart3,ijkend3

         if(.not.fluid_at(ijk) .or. pinc(ijk).eq.0) cycle
         i = i_of(ijk)
         j = j_of(ijk)
         k = k_of(ijk)
         pcell(1) = i-1
         pcell(2) = j-1
         pcell(3) = merge(1, k-1, no_k) ! =k-1 (in 3d) or =1 (in 2d)
         call set_interpolation_stencil(pcell,iw,ie,js,jn,kb,&
         ktp,interp_scheme,merge(2,3,no_k),ordernew = onew)

!Compute velocity at grid nodes and set the geometric stencil
         do k = 1, merge(1, onew, no_K)!  (3-dimn)*1+(dimn-2)*onew
            do j = 1,onew
               do i = 1,onew
                  ii = iw + i-1
                  jj = js + j-1
                  kk = kb + k-1
                  cur_ijk = funijk(imap_c(ii),jmap_c(jj),kmap_c(kk))
                  ipjk    = funijk(imap_c(ii+1),jmap_c(jj),kmap_c(kk))
                  ijpk    = funijk(imap_c(ii),jmap_c(jj+1),kmap_c(kk))
                  ipjpk   = funijk(imap_c(ii+1),jmap_c(jj+1),kmap_c(kk))

                  vol_ijk = zero
                  vol_ipjk = zero
                  vol_ijpk = zero
                  vol_ipjpk = zero

                  vol_ijkp = zero
                  vol_ipjkp = zero
                  vol_ijpkp = zero
                  vol_ipjpkp = zero

                  if(fluid_at(cur_ijk)) vol_ijk   = vol(cur_ijk)
                  if(fluid_at(ipjk))    vol_ipjk  = vol(ipjk)
                  if(fluid_at(ijpk))    vol_ijpk  = vol(ijpk)
                  if(fluid_at(ipjpk))   vol_ipjpk = vol(ipjpk)


                  if(DO_K) then
                     ijkp    = funijk(imap_c(ii),jmap_c(jj),kmap_c(kk+1))
                     ijpkp   = funijk(imap_c(ii),jmap_c(jj+1),kmap_c(kk+1))
                     ipjkp   = funijk(imap_c(ii+1),jmap_c(jj),kmap_c(kk+1))
                     ipjpkp  = funijk(imap_c(ii+1),jmap_c(jj+1),kmap_c(kk+1))


                     if(fluid_at(ijkp))     vol_ijkp   = vol(ijkp)
                     if(fluid_at(ipjkp))    vol_ipjkp  = vol(ipjkp)
                     if(fluid_at(ijpkp))    vol_ijpkp  = vol(ijpkp)
                     if(fluid_at(ipjpkp))   vol_ipjpkp = vol(ipjpkp)

                  endif
                  gstencil(i,j,k,1) = xe(ii)
                  gstencil(i,j,k,2) = yn(jj)
                  gstencil(i,j,k,3) = merge(DZ(1), ZT(KK), NO_K)

                  VOL_TOT_SCAL = ZERO

                  VOL_TOT_SCAL = vol_ijk + vol_ipjk + vol_ijpk + vol_ipjpk + &
                  & vol_ijkp + vol_ipjkp + vol_ijpkp + vol_ipjpkp

                  VOL_TOT_VEC = ZERO

                  VOL_TOT_VEC(1) = VOL(CUR_IJK) + VOL(IJPK)
                  VOL_TOT_VEC(2) = VOL(CUR_IJK) + VOL(IPJK)

                  DO M = 1, DES_MMAX
                     VEL_SOL_STENCIL(i,j,k,1,M) = pic_u_s(cur_ijk,M)*vol(cur_ijk) + pic_u_s(ijpk,M)*vol(ijpk)

                     VEL_SOL_STENCIL(i,j,k,2,M) = pic_v_s(cur_ijk,M)*vol(cur_ijk) + pic_v_s(ipjk,M)*vol(ipjk)
                  ENDDO
                  !ep_g(cur_ijk)*vol_ijk+ ep_g(ipjk)*vol_ipjk+ &
                  !   &  ep_g(ijpk)*vol_ijpk + ep_g(ipjpk)*vol_ipjpk,
                  sstencil(i,j,k) = ep_g(cur_ijk)*vol_ijk+ ep_g(ipjk)*vol_ipjk+ &
                  &  ep_g(ijpk)*vol_ijpk + ep_g(ipjpk)*vol_ipjpk

                  psgradstencil(i,j,k,1) = PS_FORCE_PIC(cur_ijk,1)*VOL(CUR_IJK) + PS_FORCE_PIC(ijpk,1)*VOL(IJPK)

                  psgradstencil(i,j,k,2) = PS_FORCE_PIC(cur_ijk,2)*VOL(CUR_IJK) + PS_FORCE_PIC(ipjk,2)*VOL(IPJK)

                  vstencil(i,j,k,1) = u_g(cur_ijk)*vol(cur_ijk) + u_g(ijpk)*vol(ijpk)
                  vstencil(i,j,k,2) = v_g(cur_ijk)*vol(cur_ijk) + v_g(ipjk)*vol(ipjk)
                  if(DO_K) then
                     VOL_TOT_VEC(1) = VOL_TOT_VEC(1) + VOL(IJKP) + VOL(IJPKP)
                     VOL_TOT_VEC(2) = VOL_TOT_VEC(2) + VOL(IJKP) + VOL(IPJKP)
                     VOL_TOT_VEC(3) = VOL(CUR_IJK) + VOL(IPJK) + VOL(IJPK) + VOL(IPJPK)
                     sstencil(i,j,k) = sstencil(i,j,k) + ep_g(ijkp)*vol_ijkp + ep_g(ipjkp)*vol_ipjkp &
                     & + ep_g(ijpkp)*vol_ijpkp + ep_g(ipjpkp)*vol_ipjpkp

                     psgradstencil(i,j,k,1) = psgradstencil(i,j,k,1)+PS_FORCE_PIC(ijkp,1)*VOL(ijkp) + PS_FORCE_PIC(ijpkp,1)*vol(ijpkp)
                     psgradstencil(i,j,k,2) = psgradstencil(i,j,k,2)+PS_FORCE_PIC(ijkp,2)*vol(ijkp) + PS_FORCE_PIC(ipjkp,2)*vol(ipjkp)
                     psgradstencil(i,j,k,3) = PS_FORCE_PIC(cur_ijk,3)*vol(cur_ijk)+&
                     & PS_FORCE_PIC(ijpk,3)*vol(ijpk)+PS_FORCE_PIC(ipjk,3)*vol(ipjk)+&
                     & PS_FORCE_PIC(ipjpk,3)*vol(ipjpk)

                     vstencil(i,j,k,1) = vstencil(i,j,k,1) + u_g(ijkp)*vol(ijkp) + u_g(ijpkp)*vol(ijpkp)

                     vstencil(i,j,k,2) = vstencil(i,j,k,2) + v_g(ijkp)*vol(ijkp) + v_g(ipjkp)*vol(ipjkp)
                     vstencil(i,j,k,3) = w_g(cur_ijk)*vol(cur_ijk)+&
                     & w_g(ijpk)*vol(ijpk)+w_g(ipjk)*vol(ipjk)+w_g(ipjpk)*vol(ipjpk)

                     DO M = 1, DES_MMAX
                        VEL_SOL_STENCIL(i,j,k,1, M) = VEL_SOL_STENCIL(i,j,k,1,M) &
                        & + PIC_u_s(ijkp,M)*vol(ijkp) + PIC_u_s(ijpkp,M)*vol(ijpkp)
                        VEL_SOL_STENCIL(i,j,k,2, M) = VEL_SOL_STENCIL(i,j,k,2,M) &
                        & + PIC_v_s(ijkp,M)*vol(ijkp) + PIC_v_s(ipjkp,M)*vol(ipjkp)
                        VEL_SOL_STENCIL(i,j,k,3, M) = PIC_w_s(cur_ijk,M)*vol(cur_ijk) +&
                        PIC_w_s(ijpk,M)*vol(ijpk)+PIC_w_s(ipjk,M)*vol(ipjk)+PIC_w_s(ipjpk,M)*vol(ipjpk)
                     ENDDO
                  else
                     psgradstencil(i,j,k,3) = 0.d0
                     VEL_SOL_STENCIL(i,j,k,3, 1:DES_MMAX) = 0.d0
                     vstencil(i,j,k,3) = 0.d0

                  endif
                  DO IDIM = 1, merge(2,3,NO_K)
                     IF(VOL_TOT_VEC(IDIM).GT.ZERO)  THEN
                        psgradstencil(i,j,k,idim) = psgradstencil(i,j,k,idim)/VOL_TOT_VEC(idim)

                        VEL_SOL_STENCIL(i,j,k,idim, 1:DES_MMAX) = VEL_SOL_STENCIL(i,j,k,idim, 1:DES_MMAX)/VOL_TOT_VEC(idim)

                        vstencil(i,j,k,idim) = vstencil(i,j,k,idim)/VOL_TOT_VEC(idim)

                        !no need for if as sum of positive numbers can only be zero
                        !if and only if each one of them are zero
                     ENDIF

                  ENDDO


                  !write(*,*) 'epg*vol     = ', ep_g(cur_ijk)*vol_ijk, ep_g(ipjk)*vol_ipjk, &
                  !&  ep_g(ijpk)*vol_ijpk , ep_g(ipjpk)*vol_ipjpk,  ep_g(cur_ijk)*vol_ijk+ ep_g(ipjk)*vol_ipjk+ &
                  !&  ep_g(ijpk)*vol_ijpk + ep_g(ipjpk)*vol_ipjpk,sstencil(i,j,k)


                  if(VOL_TOT_SCAL.gt.zero) sstencil(i,j,k) = sstencil(i,j,k)/VOL_TOT_SCAL


               enddo
            enddo
         enddo

                  !loop through particles in the cell
         do nindx = 1,pinc(ijk)
            np = pic(ijk)%p(nindx)
            m = pijk(np,5)

            if(NO_K) then !2-D
               call interpolator(gstencil(1:onew,1:onew,1,1:dimn), &
               psgradstencil(1:onew,1:onew,1,1:dimn), &
               des_pos_new(1:dimn,np),PS_GRAD(np,1:dimn),  &
               onew,interp_scheme,weightp)
            else !3-D, diff in psgradstencil size
               call interpolator(gstencil(1:onew,1:onew,1:onew,1:dimn), &
               psgradstencil(1:onew,1:onew,1:onew,1:dimn), &
               des_pos_new(1:dimn,np),PS_GRAD(np,1:dimn),  &
               onew,interp_scheme,weightp)
            endif

            do idim = 1,  merge(2,3,NO_K)
               AVGSOLVEL_P(IDIM,NP) = ARRAY_DOT_PRODUCT(VEL_SOL_STENCIL(:,:,:,IDIM,M),WEIGHTP(:,:,:))
               VEL_FP(IDIM,NP) = ARRAY_DOT_PRODUCT(VSTENCIL(:,:,:,IDIM),WEIGHTP(:,:,:))
            ENDDO

            EPG_P(NP) = ARRAY_DOT_PRODUCT(SSTENCIL(:,:,:),WEIGHTP(:,:,:))

         END DO
      END DO


      END SUBROUTINE MPPIC_COMPUTE_PS_GRAD


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: MPPIC_BC_U_S                                            !
!  Purpose:                                                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE MPPIC_BC_U_S

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE parallel
      USE constant
      USE geometry
      USE indices
      USE compar
      USE mfix_pic
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Error index
      INTEGER :: IER
! Boundary condition
      INTEGER ::  L
! Indices
      INTEGER :: I, J, K, IM, I1, I2, J1, J2, K1, K2, IJK,&
                 JM, KM, IJKW, IMJK, IPJK, IP, IJK_WALL
!-----------------------------------------------
! Include statement functions
!-----------------------------------------------
      INCLUDE '../function.inc'
!-----------------------------------------------


! Set the default boundary conditions
      IF (DO_K) THEN
         K1 = 1
         DO J1 = jmin3,jmax3
            DO I1 = imin3, imax3
               IF (.NOT.IS_ON_myPE_owns(I1,J1,K1)) CYCLE
               IJK_WALL = FUNIJK(I1,J1,K1)
               IJK = FUNIJK(I1,J1,K1+1)
               PIC_U_S(IJK_WALL, :) = PIC_U_S(IJK,:)
            END DO
         END DO

         K1 = KMAX2
         DO J1 = jmin3,jmax3
            DO I1 = imin3, imax3
               IF (.NOT.IS_ON_myPE_owns(I1,J1,K1)) CYCLE
               IJK_WALL = FUNIJK(I1,J1,K1)
               IJK = FUNIJK(I1,J1,K1-1)
               PIC_U_S(IJK_WALL, :) = PIC_U_S(IJK,:)
            END DO
         END DO
      ENDIF

      J1 = 1
      DO K1 = kmin3, kmax3
         DO I1 = imin3, imax3
            IF (.NOT.IS_ON_myPE_owns(I1,J1,K1)) CYCLE
            IJK_WALL = FUNIJK(I1,J1,K1)
            IJK = FUNIJK(I1,J1+1,K1)
            PIC_U_S(IJK_WALL, :) = PIC_U_S(IJK,:)
         END DO
      END DO
      J1 = JMAX2
      DO K1 = kmin3, kmax3
         DO I1 = imin3, imax3
            IF (.NOT.IS_ON_myPE_owns(I1,J1,K1)) CYCLE
            IJK_WALL = FUNIJK(I1,J1,K1)
            IJK = FUNIJK(I1,J1-1,K1)
            PIC_U_S(IJK_WALL, :) = PIC_U_S(IJK,:)
         END DO
      END DO

      RETURN
      END SUBROUTINE MPPIC_BC_U_S



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: MPPIC_BC_V_S                                            !
!  Purpose:                                                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE MPPIC_BC_V_S


!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE parallel
      USE constant
      USE geometry
      USE indices
      USE compar
      USE mfix_pic
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Error index
      INTEGER          IER
! Boundary condition
      INTEGER          L
! Indices
      INTEGER          I,  J, K, JM, I1, I2, J1, J2, K1, K2, IJK,&
                       IM, KM, IJKS, IJMK, IJPK, IJK_WALL
!-----------------------------------------------
! Include statement functions
!-----------------------------------------------
      INCLUDE '../function.inc'
!-----------------------------------------------

! Set the default boundary conditions
      IF (DO_K) THEN
         K1 = 1
         DO J1 = jmin3,jmax3
            DO I1 = imin3, imax3
               IF (.NOT.IS_ON_myPE_owns(I1,J1,K1)) CYCLE
               IJK_WALL = FUNIJK(I1,J1,K1)
               IJK = FUNIJK(I1,J1,K1+1)
               PIC_V_S(IJK_WALL, :) = PIC_V_S(IJK,:)
            END DO
         END DO
         K1 = KMAX2
         DO J1 = jmin3,jmax3
            DO I1 = imin3, imax3
               IF (.NOT.IS_ON_myPE_owns(I1,J1,K1)) CYCLE
               IJK_WALL = FUNIJK(I1,J1,K1)
               IJK = FUNIJK(I1,J1,K1-1)
               PIC_V_S(IJK_WALL, :) = PIC_V_S(IJK,:)
            END DO
         END DO
      ENDIF

      I1 = 1
      DO K1 = kmin3, kmax3
         DO J1 = jmin3, jmax3
            IF (.NOT.IS_ON_myPE_owns(I1,J1,K1)) CYCLE
            IJK_WALL = FUNIJK(I1,J1,K1)
            IJK = FUNIJK(I1+1,J1,K1)
            PIC_V_S(IJK_WALL, :) = PIC_V_S(IJK,:)

         END DO
      END DO
      I1 = IMAX2
      DO K1 = kmin3, kmax3
         DO J1 = jmin3, jmax3
            IF (.NOT.IS_ON_myPE_owns(I1,J1,K1)) CYCLE
            IJK_WALL = FUNIJK(I1,J1,K1)
            IJK = FUNIJK(I1-1,J1,K1)
            PIC_V_S(IJK_WALL, :) = PIC_V_S(IJK,:)

         END DO
      END DO
      END SUBROUTINE MPPIC_BC_V_S



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: MPPIC_BC_W_S                                            !
!  Purpose:                                                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE MPPIC_BC_W_S

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE parallel
      USE constant
      USE geometry
      USE indices
      USE compar
      USE mfix_pic
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Error index
      INTEGER :: IER
! Boundary condition
      INTEGER :: L
! Indices
      INTEGER :: I, J, K, KM, I1, I2, J1, J2, K1, K2, IJK,&
                 IM, JM, IJKB, IJKM, IJKP, IJK_WALL
!-----------------------------------------------
! Include statement functions
!-----------------------------------------------
      INCLUDE '../function.inc'
!-----------------------------------------------

! Set the default boundary conditions
      J1 = 1
      DO K1 = kmin3,kmax3
         DO I1 = imin3,imax3
            IF (.NOT.IS_ON_myPE_owns(I1,J1,K1)) CYCLE
            IJK = FUNIJK(I1,J1+1,K1)
            IJK_WALL = FUNIJK(I1,J1,K1)
            PIC_W_S(IJK_WALL,:) = PIC_W_S(IJK,:)
         END DO
      END DO
      J1 = JMAX2
      DO K1 = kmin3, kmax3
         DO I1 = imin3, imax3
           IF (.NOT.IS_ON_myPE_owns(I1,J1,K1)) CYCLE
            IJK = FUNIJK(I1,J1-1,K1)
            IJK_WALL = FUNIJK(I1,J1,K1)
            PIC_W_S(IJK_WALL,:) = PIC_W_S(IJK,:)
         END DO
      END DO
      I1 = 1
      DO K1 = kmin3, kmax3
         DO J1 = jmin3, jmax3
            IF (.NOT.IS_ON_myPE_owns(I1,J1,K1)) CYCLE
            IJK = FUNIJK(I1+1,J1,K1)
            IJK_WALL = FUNIJK(I1,J1,K1)
            PIC_W_S(IJK_WALL,:) = PIC_W_S(IJK,:)
         END DO
      END DO
      I1 = IMAX2
      DO K1 = kmin3,kmax3
         DO J1 = jmin3,jmax3
            IF (.NOT.IS_ON_myPE_owns(I1,J1,K1)) CYCLE
            IJK = FUNIJK(I1-1,J1,K1)
            IJK_WALL = FUNIJK(I1,J1,K1)
            PIC_W_S(IJK_WALL,:) = PIC_W_S(IJK,:)
         END DO
      END DO


      END SUBROUTINE MPPIC_BC_W_S


      SUBROUTINE WRITE_NODEDATA(bufin, funit)
      USE param
      USE param1
      USE parallel
      USE constant
      USE run
      USE geometry
      USE indices
      USE compar
      USE discretelement
      use desmpi

      IMPLICIT NONE
      integer, intent(in) :: funit
      double precision, dimension(:), intent(in)  :: bufin

      integer :: ijk, i, j,k
      INCLUDE '../function.inc'
      INCLUDE '../fun_avg1.inc'
      INCLUDE '../fun_avg2.inc'

      write(funit,*)'VARIABLES= ',' "I" ',' "J" ',' "K" ',' "DES_ROPS_NODE" '

      write(funit, *)'ZONE F=POINT, I=', (IEND1-ISTART2)+1,  ', J=', JEND1-JSTART2+1, ', K=', KEND1-KSTART2 + 1

      DO K = KSTART2, KEND1
         DO J = JSTART2, JEND1
            DO I = ISTART2, IEND1
               IJK = funijk(I, J, K)
               !WRITE(*,*) 'IJK = ', IJK, I, J, K , SIZE(BUFIN,1)

               WRITE(funit, '(3(2x, i10), 3x, g17.8)') I, J, K , DES_ROPS_NODE(IJK,1)
            ENDDO
         ENDDO
      ENDDO
      CLOSE(funit, status = 'keep')
      end SUBROUTINE WRITE_NODEDATA

      SUBROUTINE WRITE_MPPIC_VEL_S
      USE param
      USE param1
      USE parallel
      USE constant
      USE run
      USE geometry
      USE indices
      USE compar
      USE fldvar, only : ep_g
      USE discretelement
      USE mfix_pic
      implicit none
      integer :: i, j, k, ijk, fluid_ind, LL, PC, IDIM
      double precision :: zcor
      character*100 :: filename
      logical finish
      INCLUDE '../function.inc'

!      INCLUDE '../ep_s1.inc'
!      INCLUDE '../ep_s2.inc'

      WRITE(filename,'(A,"_",I5.5,".dat")') TRIM(RUN_NAME)//'_U_S_',myPE
      OPEN(1000, file = TRIM(filename), form ='formatted', status='unknown')
      IF(DIMN.eq.2) then
         write(1000,*)'VARIABLES= ',' "X" ',' "Y" ',' "Z" ',' "EP_s " ', ' "U_S" ', ' "V_S" ',' "DES_U_s" ', ' "DES_V_s" '!, ' "P_S_FOR1" ', ' "P_S_FOR2" '
         write(1000,*)'ZONE F=POINT, I=', (IEND3-ISTART3)+1,  ', J=', JEND3-JSTART3+1, ', K=', KEND3-KSTART3 + 1
      else
         write(1000,*)'VARIABLES= ',' "X" ',' "Y" ',' "Z" ',' "EP_s " ', ' "U_S" ', ' "V_S" ', ' "W_S" ',' "DES_U_s" ', ' "DES_V_s" ', ' "DES_W_s" '!, &
        ! & ' "P_S_FOR1" ', ' "P_S_FOR2" ', ' "P_S_FOR3" '
         write(1000,*)'ZONE F=POINT, I=', (IEND3-ISTART3)+1,  ', J=', JEND3-JSTART3+1, ', K=', KEND3-KSTART3 + 1
      ENDIF
      DO K=KSTART3, KEND3
         DO J=JSTART3, JEND3
            DO I=ISTART3, IEND3
               IJK  = FUNIJK(I,J,K)
               IF(FLUID_AT(IJK)) THEN
                  FLUID_IND = 1
               ELSE
                  FLUID_IND = 0
               END IF
               IF(DIMN.EQ.2) THEN
                  ZCOR = ZT(K)
                  WRITE(1000,'(3(2X,G17.8),4( 2X, G17.8))')  XE(I-1)+DX(I), YN(J-1)+DY(J),ZCOR, &
                     1.D0 - EP_G(IJK), PIC_U_S(IJK,1), PIC_V_S(IJK,1), DES_U_S(IJK,1), &
                     DES_V_S(IJK,1) !, PS_FORCE_PIC(IJK,1), PS_FORCE_PIC(IJK,2)
               ELSE
                  ZCOR = ZT(K-1) + DZ(K)
                  WRITE(1000,'(3(2X,G17.8),4( 2X, G17.8))')  XE(I-1)+DX(I), YN(J-1)+DY(J),ZCOR, &
                     1.D0 - EP_G(IJK), PIC_U_S(IJK,1), PIC_V_S(IJK,1), PIC_W_S(IJK,1), DES_U_S(IJK,1), &
                    DES_V_S(IJK,1), DES_W_S(IJK,1)!, PS_FORCE_PIC(IJK,1), PS_FORCE_PIC(IJK,2),  PS_FORCE_PIC(IJK,3)
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      CLOSE(1000, STATUS='KEEP')

      return

      WRITE(FILENAME,'(A,"_",I5.5,".DAT")') TRIM(RUN_NAME)//'_PS_FORCE_',myPE
      OPEN(1000, file = TRIM(filename), form ='formatted', status='unknown')

      IF(DIMN.eq.3) write(1000,*)'VARIABLES= ',' "X" ',' "Y" ',' "Z" ',' "DELPX" ', '"DELPY"', '"DELPZ" ',' "US_part" ', '"VS_part"' , '"WS_part"', '"EP_s_part"'
      IF(DIMN.eq.2) write(1000,*)'VARIABLES= ',' "X" ',' "Y" ',' "DELPX" ', '"DELPY"', ' "US_part" ', '"VS_part"' , '"EP_S_part"'

      PC = 1
      DO LL = 1, MAX_PIP

         IF(PC .GT. PIP) EXIT
         IF(.NOT.PEA(LL,1)) CYCLE
         pc = pc+1
         IF(PEA(LL,4)) CYCLE

         WRITE(1000,'(10( 2x, g17.8))') (DES_POS_NEW(IDIM, LL), IDIM = 1, DIMN), (PS_GRAD(LL, IDIM) , IDIM = 1, DIMN), (AVGSOLVEL_P (IDIM, LL) , IDIM = 1, DIMN), 1-EPg_P(LL)
      ENDDO
               close(1000, status='keep')

      !write(*,*) 'want to quit ?', LL, mAX_PIP, PIP
      !read(*,*) finish
      !if(finish) STOp
      END SUBROUTINE WRITE_MPPIC_VEL_S





