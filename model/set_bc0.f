!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SET_BC0                                                 C
!  Purpose: This module does the initial setting of all boundary       C
!           conditions. The user specifications of the boundary        C
!           conditions are checked for veracity in check_data_07       C
!                                                                      C
!  Author: M. Syamlal                                 Date: 29-JAN-92  C
!  Reviewer: P. Nicoletti, W. Rogers, S. Venkatesan   Date: 29-JAN-92  C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: Check whether the velocity components have the correct     C
!           sign at MASS_INFLOW cells                                  C
!  Author: M. Syamlal                                 Date: 08-JUL-92  C
!  Reviewer: M. Syamlal                               Date: 11-DEC-92  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: BC_DEFINED, BC_TYPE, BC_DT_0, TIME, BC_Jet_g0,C
!                        BC_K_b, BC_K_t, BC_J_s, BC_J_n, BC_I_w,       C
!                        BC_I_e, BC_PLANE, BC_EP_g, BC_P_g, BC_T_g,    C
!                        BC_T_s,  BC_U_g, BC_V_g, BC_W_g,              C
!                        MMAX, BC_ROP_s, BC_U_s, BC_V_s, BC_W_s        C
!  Variables modified: BC_TIME, BC_V_g, I, J, K, IJK, EP_g, P_g, T_g,  C
!                      T_s, U_g, V_g, W_g, ROP_s, U_s, V_s, W_s,       C
!                      M                                               C
!                                                                      C
!  Local variables: L, IJK1, IJK2, IJK3                                C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE SET_BC0

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE geometry
      USE compar
      USE mpi_utility
      USE physprop
      USE constant
      USE bc
      USE fldvar
      USE indices
      USE run
      USE funits
      USE scales
      USE scalars
      USE boundfunijk
      USE toleranc
      USE sendrecv
      USE mms
      USE functions

      IMPLICIT NONE

!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Local index for boundary condition
      INTEGER ::  L
! indices
      INTEGER :: I, J, K, IJK, M, N
! Local index for setting U velocity b.c.
      INTEGER :: IJK1
! Local index for setting V velocity b.c.
      INTEGER :: IJK2
! Local index for setting W velocity b.c.
      INTEGER :: IJK3
! number densities for use in GHD theory only
      DOUBLE PRECISION :: nM, nTOT
!-----------------------------------------------
! External functions
!-----------------------------------------------
      DOUBLE PRECISION, EXTERNAL :: EOSG, CALC_MW
!----------------------------------------------


! Incompressible cases require that Ppg specified for one cell.
! The following attempts to pick an appropriate cell.
      CALL SET_IJK_P_G

      IF(USE_MMS) THEN
! IJK_P_G is set as UNDEFINED for MMS since current IJK_P_G setting
! is not a second order operation.
         IJK_P_G = UNDEFINED_I
! Calculate MMS variables. Better place might be inside interate before
! every time-step. (?)
         CALL CALCULATE_MMS
         CALL CALCULATE_MMS_SOURCE
      END IF

! Set the boundary conditions - Defining the field variables at the
! boundaries according to the user specifications. These are not the
! real values in the wall cells, only initial guesses.
      DO L = 1, DIMENSION_BC
         IF (BC_DEFINED(L)) THEN

! setting wall (nsw, fsw, psw) boundary conditions
!----------------------------------------------------------------->>>
            IF (BC_TYPE(L)=='FREE_SLIP_WALL' .OR.                      &
                BC_TYPE(L)=='NO_SLIP_WALL' .OR.                        &
                BC_TYPE(L)=='PAR_SLIP_WALL') THEN

               DO K = BC_K_B(L), BC_K_T(L)
               DO J = BC_J_S(L), BC_J_N(L)
               DO I = BC_I_W(L), BC_I_E(L)

                  IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                        IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
                  IJK = BOUND_FUNIJK(I,J,K)

                  IF (WALL_AT(IJK)) THEN
! this conditional is probably unnecessary since this section already
! falls under the 'wall' if/else branch for the bc_type check...
                     IF(BC_Tw_g(L) /= UNDEFINED) T_g(IJK) = BC_Tw_g(L)

                     IF(NMAX(0) > 0)                                   &
                        WHERE (BC_Xw_G(L,:NMAX(0)) /= UNDEFINED)       &
                        X_G(IJK,:NMAX(0)) = BC_Xw_G(L,:NMAX(0))

                     IF(SMAX > 0)                                      &
                        WHERE(BC_Tw_s(L,:SMAX) /= UNDEFINED)           &
                        T_s(IJK,:SMAX) = BC_Tw_s(L,:SMAX)

                     IF(MMAX > 0)                                      &
                        WHERE(BC_Thetaw_m(L,:MMAX) /= UNDEFINED)       &
                        Theta_m(IJK,:MMAX) = BC_Thetaw_m(L,:MMAX)

                     DO M = 1, SMAX
                        IF(NMAX(M) > 0)                                &
                           WHERE (BC_Xw_s(L,M,:NMAX(M)) /= UNDEFINED)  &
                           X_s(IJK,M,:NMAX(M)) = BC_Xw_s(L,M,:NMAX(M))
                     ENDDO

                     IF(NScalar > 0)                                   &
                        WHERE (BC_ScalarW(L,:NScalar) /= UNDEFINED)    &
                          Scalar(IJK,:NScalar) = BC_ScalarW(L,:NScalar)
                  ENDIF   ! end if (wall_at(ijk))

               ENDDO   ! end do (i=bc_i_w(l),bc_i_e(l))
               ENDDO   ! end do (j=bc_j_s(l),bc_j_n(l))
               ENDDO   ! end do (k=bc_k_b(l),bc_k_t(l))
! end setting (FSW, NSW or PSW)
!-----------------------------------------------------------------<<<

            ELSE   ! else branch of if(FSW, NSW or PSW)


! setting all other 'non-wall' boundary conditions
!----------------------------------------------------------------->>>


! initializing for time dependent mass inflow
               BC_JET_G(L) = UNDEFINED
               IF (BC_DT_0(L) /= UNDEFINED) THEN
                  BC_TIME(L) = TIME + BC_DT_0(L)
                  BC_OUT_N(L) = 0
                  BC_MOUT_G(L) = ZERO
                  BC_VOUT_G(L) = ZERO
                  IF (SMAX > 0) THEN
                     BC_MOUT_S(L,:SMAX) = ZERO
                     BC_VOUT_S(L,:SMAX) = ZERO
                  ENDIF
                  BC_JET_G(L) = BC_JET_G0(L)
                  IF (BC_JET_G(L) /= UNDEFINED) THEN
                     SELECT CASE (TRIM(BC_PLANE(L)))
                     CASE ('W')
                        BC_U_G(L) = BC_JET_G(L)
                     CASE ('E')
                        BC_U_G(L) = BC_JET_G(L)
                     CASE ('S')
                        BC_V_G(L) = BC_JET_G(L)
                     CASE ('N')
                        BC_V_G(L) = BC_JET_G(L)
                     CASE ('B')
                        BC_W_G(L) = BC_JET_G(L)
                     CASE ('T')
                        BC_W_G(L) = BC_JET_G(L)
                     END SELECT
                  ENDIF
               ELSE
                  BC_TIME(L) = UNDEFINED
               ENDIF


               DO K = BC_K_B(L), BC_K_T(L)
               DO J = BC_J_S(L), BC_J_N(L)
               DO I = BC_I_W(L), BC_I_E(L)

                  IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                        IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
                  IJK = BOUND_FUNIJK(I,J,K)

! this conditional may still be necessary to avoid cyclic boundaries.
! cyclic boundaries would fall under this branch but are considered
! 'wall_at'
                  IF (.NOT.WALL_AT(IJK)) THEN

! setting p_outflow_at or outflow_at boundary conditions
! ---------------------------------------------------------------->>>
                     IF(P_OUTFLOW_AT(IJK) .OR. OUTFLOW_AT(IJK)) THEN

! Why check if the bc quantity for the field variable is defined before
! performing the assignment?:
! For a new run the field variable will be undefined in these locations
! so any check is irrelevant. However, for a restart run the field
! variable in the boundary cell may have an existing value based on the
! flow field. In this case, the check acts to maintain the existing value
! unless given a user defined BC value. This is unlike the PI, MI and MO
! treatment where the values of field variables in boundary cells are
! strictly overwritten regardless if the BC is defined. In addition,
! check_data_07 basically ensures such quantities will be defined for
! MI or PI if they are required.

                        P_STAR(IJK) = ZERO
                        P_G(IJK) = SCALE(BC_P_G(L))
                        IF (BC_EP_G(L) /= UNDEFINED)                   &
                           EP_G(IJK) = BC_EP_G(L)
! unless gas flow is entering back into the domain and ep_g in the
! boundary cell is defined, t_g will be set to its adjacent fluid cell
! value in set_bc1. so why do this here? if during the initial start-up,
! the gas has backflow and ep_g is defined then t_g may be undefined in
! the boundary.. otherwise this seems unnecessary?

                        T_G(IJK)= merge(BC_T_G(L), TMIN,               &
                           BC_T_G(L) /= UNDEFINED)

                        IF (NMAX(0) > 0)                               &
                           WHERE (BC_X_G(L,:NMAX(0)) /= UNDEFINED)     &
                           X_G(IJK,:NMAX(0)) = BC_X_G(L,:NMAX(0))

                        IF (NScalar > 0)                               &
                           WHERE (BC_Scalar(L,:NScalar) /= UNDEFINED)  &
                           Scalar(IJK,:NScalar) = BC_Scalar(L,:NScalar)

                        IF (K_Epsilon) THEN
                           IF (BC_K_Turb_G(L) /= UNDEFINED)            &
                              K_Turb_G(IJK) = BC_K_Turb_G(L)
                           IF (BC_E_Turb_G(L) /= UNDEFINED)            &
                              E_Turb_G(IJK) = BC_E_Turb_G(L)
                        ENDIF

                        DO M = 1, SMAX
                           IF (BC_ROP_S(L,M) /= UNDEFINED)             &
                              ROP_S(IJK,M) = BC_ROP_S(L,M)
                           IF(BC_T_S(L,M) /= UNDEFINED)                &
                              T_S(IJK,M)=BC_T_S(L,M)
                           IF (BC_THETA_M(L,M) /= UNDEFINED)           &
                              THETA_M(IJK,M) = BC_THETA_M(L,M)

                           IF (NMAX(M) > 0)                            &
                              WHERE (BC_X_S(L,M,:NMAX(M)) /= UNDEFINED)&
                              X_S(IJK,M,:NMAX(M)) = BC_X_S(L,M,:NMAX(M))
                        ENDDO

! for GHD theory to compute mixture BC of velocity and density
                        IF(KT_TYPE_ENUM == GHD_2007) THEN
                           ROP_S(IJK,MMAX) = ZERO
                           nTOT = ZERO
                           THETA_M(IJK,MMAX) = ZERO
                           DO M = 1, SMAX
                              IF (BC_ROP_S(L,M) /= UNDEFINED) THEN
                                 ROP_S(IJK,MMAX) = ROP_S(IJK,MMAX) +   &
                                    BC_ROP_S(L,M)
                                 nM = BC_ROP_S(L,M)*6d0 /              &
                                    (PI*D_p(IJK,M)**3*RO_S(IJK,M))
                                 nTOT = nTOT + nM
                              IF (BC_THETA_M(L,M) /= UNDEFINED)        &
                                 THETA_M(IJK,MMAX) = THETA_M(IJK,MMAX)+&
                                    nM*BC_THETA_M(L,M)
                            ENDIF
                          ENDDO
                          IF(ROP_S(IJK,MMAX) > ZERO) THETA_M(IJK,MMAX)=&
                             THETA_M(IJK,MMAX) / nTOT

! Set MMS BCs when PO boundary condition is used.
                        IF (USE_MMS) THEN
                           P_G(IJK) = SCALE(MMS_P_G(IJK))
                           EP_G(IJK) = MMS_EP_G(IJK)
                           T_G(IJK) = MMS_T_G(IJK)

                           M = 1 ! Single solid phase for MMS cases
                           ROP_S(IJK,M) = MMS_ROP_S(IJK)
                           T_S(IJK,M) = MMS_T_S(IJK)
                           THETA_M(IJK,M) = MMS_THETA_M(IJK)
                        ENDIF ! end if(USE_MMS)

                        ENDIF
! end setting P_outflow_at or outflow_at boundary conditions
! ----------------------------------------------------------------<<<

                     ELSE   ! else branch of if(P_outflow_at .or. outflow_at)

! setting for p_inflow, mass_inflow or mass_outflow boundary
! conditions
! ---------------------------------------------------------------->>>

! Unlike PO or O boundary condition branch above, in this branch no
! checks are made to determine if the BC is defined before it is
! assigned to the field variable. However, check_data_07 generally
! ensures such BC quantities will be defined if they are needed for
! PI and MI boundaries.
                        P_STAR(IJK) = ZERO
                        P_G(IJK) = SCALE(BC_P_G(L))
                        EP_G(IJK) = BC_EP_G(L)
                        T_G(IJK) = BC_T_G(L)

                        IF (NMAX(0) > 0)                               &
                          X_G(IJK,:NMAX(0)) = BC_X_G(L,:NMAX(0))

                        IF (NScalar > 0)                               &
                          Scalar(IJK,:NScalar) = BC_Scalar(L,:NScalar)

                        IF (K_Epsilon) THEN
                          K_Turb_G(IJK) = BC_K_Turb_G(L)
                          E_Turb_G(IJK) = BC_E_Turb_G(L)
                        ENDIF

                        DO M = 1, SMAX
                          ROP_S(IJK,M) = BC_ROP_S(L,M)
                          T_S(IJK,M) = BC_T_S(L,M)
                          THETA_M(IJK,M) = BC_THETA_M(L,M)
                          IF (NMAX(M) > 0)                             &
                             X_S(IJK,M,:NMAX(M)) = BC_X_S(L,M,:NMAX(M))
                        ENDDO

! Changed to make consistent approach
                        IJK1 = IJK
                        IJK2 = IJK
                        IJK3 = IJK
                        SELECT CASE (TRIM(BC_PLANE(L)))
                          CASE ('W')
                            IJK1 = BOUND_FUNIJK(IM1(I),J,K)
                          CASE ('S')
                            IJK2 = BOUND_FUNIJK(I,JM1(J),K)
                          CASE ('B')
                            IJK3 = BOUND_FUNIJK(I,J,KM1(K))
                        END SELECT

! When the boundary plane is W, S, or B the velocity components
! need to be set for the boundary cell and adjacent fluid cell.
                        U_G(IJK) = BC_U_G(L)
                        V_G(IJK) = BC_V_G(L)
                        W_G(IJK) = BC_W_G(L)
                        U_G(IJK1) = BC_U_G(L)
                        V_G(IJK2) = BC_V_G(L)
                        W_G(IJK3) = BC_W_G(L)

                        IF (MMAX > 0) THEN
                           U_S(IJK,:SMAX) = BC_U_S(L,:SMAX)
                           V_S(IJK,:SMAX) = BC_V_S(L,:SMAX)
                           W_S(IJK,:SMAX) = BC_W_S(L,:SMAX)
                           U_S(IJK1,:SMAX) = BC_U_S(L,:SMAX)
                           V_S(IJK2,:SMAX) = BC_V_S(L,:SMAX)
                           W_S(IJK3,:SMAX) = BC_W_S(L,:SMAX)
                        ENDIF

! for GHD theory to compute mixture BC of velocity and density
                        IF(KT_TYPE_ENUM == GHD_2007) THEN
                           ROP_S(IJK,MMAX) = ZERO
                           nTOT = ZERO
                           THETA_M(IJK,MMAX) = ZERO
                           U_S(IJK,MMAX) =  ZERO
                           V_S(IJK,MMAX) =  ZERO
                           W_S(IJK,MMAX) =  ZERO
                           U_S(IJK1,MMAX) =  ZERO
                           V_S(IJK2,MMAX) =  ZERO
                           W_S(IJK3,MMAX) =  ZERO
                           DO M = 1, SMAX
                              ROP_S(IJK,MMAX) = ROP_S(IJK,MMAX) + BC_ROP_S(L,M)
                              nM = BC_ROP_S(L,M)*6d0/(PI*D_p(IJK,M)**3*RO_S(IJK,M))
                              nTOT = nTOT + nM
                              THETA_M(IJK,MMAX) = THETA_M(IJK,MMAX) + nM*BC_THETA_M(L,M)
                              U_S(IJK,MMAX) = U_S(IJK,MMAX) + BC_ROP_S(L,M)*BC_U_S(L,M)
                              V_S(IJK,MMAX) = V_S(IJK,MMAX) + BC_ROP_S(L,M)*BC_V_S(L,M)
                              W_S(IJK,MMAX) = W_S(IJK,MMAX) + BC_ROP_S(L,M)*BC_W_S(L,M)
                              SELECT CASE (TRIM(BC_PLANE(L)))
                                CASE ('W')
                                  U_S(IJK1,MMAX) = U_S(IJK1,MMAX) + BC_ROP_S(L,M)*BC_U_S(L,M)
                                CASE ('S')
                                  V_S(IJK2,MMAX) = V_S(IJK2,MMAX) + BC_ROP_S(L,M)*BC_V_S(L,M)
                                CASE ('B')
                                  W_S(IJK3,MMAX) = W_S(IJK3,MMAX) + BC_ROP_S(L,M)*BC_W_S(L,M)
                              END SELECT
                           ENDDO
                           IF(ROP_S(IJK,MMAX) > ZERO) THEN
                              THETA_M(IJK,MMAX) = THETA_M(IJK,MMAX) / nTOT
                              U_S(IJK,MMAX) = U_S(IJK,MMAX) / ROP_S(IJK,MMAX)
                              V_S(IJK,MMAX) = V_S(IJK,MMAX) / ROP_S(IJK,MMAX)
                              W_S(IJK,MMAX) = W_S(IJK,MMAX) / ROP_S(IJK,MMAX)
                              SELECT CASE (TRIM(BC_PLANE(L)))
                                CASE ('W')
                                  U_S(IJK1,MMAX) = U_S(IJK1,MMAX) / ROP_S(IJK,MMAX)
                                CASE ('S')
                                  V_S(IJK2,MMAX) = V_S(IJK2,MMAX) / ROP_S(IJK,MMAX)
                                CASE ('B')
                                  W_S(IJK3,MMAX) = W_S(IJK3,MMAX) / ROP_S(IJK,MMAX)
                              END SELECT
                           ENDIF
                        ENDIF  ! end if (trim(kt_type) ='ghd')

! Set MMS BCs when MI boundary condition is used.
                        IF (USE_MMS) THEN
                           P_G(IJK) = SCALE(MMS_P_G(IJK))
                           EP_G(IJK) = MMS_EP_G(IJK)
                           T_G(IJK) = MMS_T_G(IJK)

                           DO M = 1, SMAX
                             ROP_S(IJK,M) = MMS_ROP_S(IJK)
                             T_S(IJK,M) = MMS_T_S(IJK)
                             THETA_M(IJK,M) = MMS_THETA_M(IJK)
                           ENDDO

                           IJK1 = IJK
                           IJK2 = IJK
                           IJK3 = IJK
                           SELECT CASE (TRIM(BC_PLANE(L)))
                             CASE ('W')
                               IJK1 = BOUND_FUNIJK(IM1(I),J,K)
                             CASE ('S')
                               IJK2 = BOUND_FUNIJK(I,JM1(J),K)
                             CASE ('B')
                               IJK3 = BOUND_FUNIJK(I,J,KM1(K))
                           END SELECT
! When the boundary plane is W, S, or B the velocity components
! need to be set for both sides of the boundary cell.
                          U_G(IJK) = MMS_U_G(IJK)
                          V_G(IJK) = MMS_V_G(IJK)
                          W_G(IJK) = MMS_W_G(IJK)
                          U_G(IJK1) = MMS_U_G(IJK1)
                          V_G(IJK2) = MMS_V_G(IJK2)
                          W_G(IJK3) = MMS_W_G(IJK3)

                          IF (MMAX > 0) THEN
                             U_S(IJK,:SMAX) = MMS_U_S(IJK)
                             V_S(IJK,:SMAX) = MMS_V_S(IJK)
                             W_S(IJK,:SMAX) = MMS_W_S(IJK)
                             U_S(IJK1,:SMAX) = MMS_U_S(IJK1)
                             V_S(IJK2,:SMAX) = MMS_V_S(IJK2)
                             W_S(IJK3,:SMAX) = MMS_W_S(IJK3)
                          ENDIF
                        ENDIF ! end if(USE_MMS)

                     ENDIF     ! end if/else branch of if(P_outflow_at or outflow_at)
! end setting for p_inflow, mass_inflow or mass_outflow boundary
! conditions
! ----------------------------------------------------------------<<<
                  ENDIF       ! end if (not wall_at)

               ENDDO   ! do i
               ENDDO   ! do j
               ENDDO   ! do k

            ENDIF      ! end if/else branch if(FSW, NSW or PSW)
! end setting all boundary conditions
!-----------------------------------------------------------------<<<

         ENDIF         ! if (bc_defined)
      ENDDO            ! do dimension_bc

! FIX AEOLUS Sofiane's bug fix to make T_g nonzero in k=0,1 ghost layers
! when k-decomposition employed
      call send_recv(T_G,2)
      call send_recv(P_G,2)
      call send_recv(X_G,2)

      RETURN
      END SUBROUTINE SET_BC0




!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SET_IJK_P_G                                             !
!  Purpose: Pick an appropriate control volume to specify Ppg.         !
!                                                                      !
!  Author: J. Musser                                  Date: 07-Nov-13  !
!  Reviewer:                                          Date:            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE SET_IJK_P_G

! IJK location where Ppg is fixed.
      use bc, only: IJK_P_G
! Specified constant gas density.
      use physprop, only: RO_G0

      use geometry, only: CYCLIC
      use geometry, only: CYCLIC_X, CYCLIC_X_PD, CYCLIC_X_MF
      use geometry, only: CYCLIC_Y, CYCLIC_Y_PD, CYCLIC_Y_MF
      use geometry, only: CYCLIC_Z, CYCLIC_Z_PD, CYCLIC_Z_MF

      use geometry, only: iMAX1, iMin1
      use geometry, only: jMAX1, jMin1
      use geometry, only: kMAX1, kMin1

      use geometry, only: do_K

      use funits, only: DMP_LOG

      use bc, only: BC_I_w, BC_I_e
      use bc, only: BC_J_s, BC_J_n
      use bc, only: BC_K_b, BC_K_t

      use bc, only: BC_DEFINED
      use bc, only: BC_TYPE
      use bc, only: BC_PLANE

! MFIX Runtime parameters:
      use param, only: DIMENSION_BC
      use param1, only: UNDEFINED
      use param1, only: UNDEFINED_I

      use mpi_utility

      implicit none

      INTEGER :: BCV

      CHARACTER(len=7) :: Map

      CHARACTER(len=128) :: lMsg

      INTEGER :: l3
      INTEGER :: l2, u2
      INTEGER :: l1, u1

      LOGICAL, parameter :: setDBG = .FALSE.
      LOGICAL :: dFlag

      INTEGER :: iErr

      EXTERNAL JKI_MAP, IKJ_MAP, KIJ_MAP

      dFlag = (DMP_LOG .AND. setDBG)

! Initialize.
      iErr = 0
      IJK_P_G = UNDEFINED_I

! This is not needed for compressible cases.
      IF(RO_G0 == UNDEFINED) THEN
         IF(dFlag) write(*,"(3x,A)")                                   &
            'Compressible: IJK_P_g remaining undefined.'
         return
      ENDIF

! If there are no cyclic boundaries, look for a pressure outflow.
      lpBCV: DO BCV = 1, DIMENSION_BC
         IF (.NOT.BC_DEFINED(BCV)) cycle lpBCV
         IF (BC_TYPE(BCV) == 'P_OUTFLOW') THEN
            IF(dFlag) write(*,"(3x,A)")                                &
               'Outflow PC defiend: IJK_P_g remaining undefined.'
            RETURN
         ENDIF
      ENDDO lpBCV

! Initialize.
         l3 = UNDEFINED_I
         l2 = UNDEFINED_I; u2=l2
         l1 = UNDEFINED_I; u1=l1

! If there are cyclic boundaries, flag a cell along the positive
! domain extreme in the cyclic direction (e.g., JMAX1).
      IF(CYCLIC_Y .OR. CYCLIC_Y_PD .OR. CYCLIC_Y_MF) THEN

         Map = 'JKI_MAP'
         l3 = JMAX1
         l2 = KMIN1;  u2 = KMAX1
         l1 = IMIN1;  u1 = IMAX1
         lMsg='Cyclic in Y'

      ELSEIF(CYCLIC_X .OR. CYCLIC_X_PD .OR. CYCLIC_X_MF) THEN

         Map = 'IKJ_MAP'
         l3 = IMAX1
         l2 = KMIN1;  u2 = KMAX1
         l1 = JMIN1;  u1 = JMAX1
         lMsg='Cyclic in X'

      ELSEIF(CYCLIC_Z .OR. CYCLIC_Z_PD .OR. CYCLIC_Z_MF) THEN

         Map = 'KIJ_MAP'
         l3 = KMAX1
         l2 = IMIN1;  u2 = IMAX1
         l1 = JMIN1;  u1 = JMAX1
         lMsg='Cyclic in Z'

      ENDIF

! No cyclic boundaries or pressure outflows. The IJ plane is used in
! this case to maximize search region for 2D problems.
      IF(l3 == UNDEFINED_I) THEN
         Map = 'KIJ_MAP'
         l3 = merge((KMAX1 - KMIN1)/2 + 1, KMIN1, do_K)
         l2 = IMIN1;  u2 = IMAX1
         l1 = JMIN1;  u1 = JMAX1
         lMsg='Center of domain'
      ENDIF

! Debugging messages.
      IF(dFlag) THEN
         write(*,"(3/,3x,'Map: ',A)") Map
         write(*,"(/5x,'l3:',2x,I4)") l3
         write(*,"( 5x,'l2:',2(2x,I4))") l2, u2
         write(*,"( 5x,'l1:',2(2x,I4))") l1, u1
         write(*,"( 5x,'Msg: ',A)") trim(lMsg)
      ENDIF

! Invoke the search routine.
      SELECT CASE (Map)
      CASE ('JKI_MAP')
         CALL IJK_Pg_SEARCH(l3, l2, u2, l1, u1, JKI_MAP, dFlag, iErr)
      CASE ('IKJ_MAP')
         CALL IJK_Pg_SEARCH(l3, l2, u2, l1, u1, IKJ_MAP, dFlag, iErr)
      CASE ('KIJ_MAP')
         CALL IJK_Pg_SEARCH(l3, l2, u2, l1, u1, KIJ_MAP, dFlag, iErr)
      CASE DEFAULT
         iErr = 1001
      END SELECT

      IF(iErr == 0) RETURN

! Error management.
      IF(DMP_LOG) THEN
         SELECT CASE (iErr)
         CASE ( 1001);  WRITE(UNIT_LOG, 1001); WRITE(*,1001)
         CASE ( 2000);  WRITE(UNIT_LOG, 2000); WRITE(*,2000)
         CASE ( 2001);  WRITE(UNIT_LOG, 2001); WRITE(*,2001)
         CASE ( 2002);  WRITE(UNIT_LOG, 2002); WRITE(*,2002)
         CASE DEFAULT
            WRITE(UNIT_LOG, 1000) iErr
            WRITE(*,1000) iErr
         END SELECT

         WRITE(UNIT_LOG, 9000) MAP(1:1), l3, MAP(2:2),                 &
            l2, u2, MAP(3:3), l1, u1
         WRITE(*, 9000) MAP(1:1), l3, MAP(2:2),                        &
            l2, u2, MAP(3:3), l1, u1

         WRITE(*, 9999)
         WRITE(UNIT_LOG, 9999)

      ENDIF


      CALL MFIX_EXIT(myPE)


 1000 FORMAT(//1X,70('*')/' From: SET_IJK_Pg',/,                       &
         ' Error 1000: Unknown error reported. x', I4.4)

 1001 FORMAT(//1X,70('*')/' From: SET_IJK_Pg',/,                       &
         ' Error 1001: Invalid mapping function.')

 2000 FORMAT(//1X,70('*')/' From: SET_IJK_Pg > IJK_Pg_SEARCH',/,       &
         ' Error 2000: Unknown error reported from IJK_Pg_SEARCH.')

 2001 FORMAT(//1X,70('*')/' From: SET_IJK_Pg > IJK_Pg_SEARCH',/,       &
         ' Error 2001: Unable to locate fluid cell in search region.')

 2002 FORMAT(//1X,70('*')/' From: SET_IJK_Pg > IJK_Pg_SEARCH',/,       &
         ' Error 2002: Unable to locate fluid cell owner.')

 9000 FORMAT(/' Search plane information:',/,3x,A1,': ',I8,            &
          2(/3x,A1,': ',I8,' x ',I8))

 9999 FORMAT(/' Fatal Error --> Invoking MFIX_EXIT',/1x,70('*'),2/)


      END SUBROUTINE SET_IJK_P_G


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SET_IJK_P_G                                             !
!  Purpose: Pick an appropriate control volume to specify Ppg.         !
!                                                                      !
!  Author: J. Musser                                  Date: 07-Nov-13  !
!  Reviewer:                                          Date:            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE IJK_Pg_SEARCH(ll3, ll2, lu2, ll1, lu1,  lMAP,         &
         ldFlag, iErr)

! IJK location where Ppg is fixed.
      use bc, only: IJK_P_g

      use indices
      use mpi_utility
      use functions

      implicit none

      INTEGER, INTENT(IN)  :: ll3
      INTEGER, INTENT(IN)  :: ll2, lu2
      INTEGER, INTENT(IN)  :: ll1, lu1

      LOGICAL, intent(in) :: ldFlag

      INTEGER, INTENT(OUT)  :: iErr

      EXTERNAL lMAP

      INTEGER :: lc2, lS2, lE2
      INTEGER :: lc1, lS1, lE1

      INTEGER :: I, J, K, IJK
      LOGICAL :: recheck

      INTEGER :: IJK_Pg_Owner, proc
      INTEGER :: gIJK(0:numPEs-1)

      INTEGER :: I_J_K_Pg(3)

      INTEGER :: lpCnt

      CHARACTER(len=32) :: cInt

! Initialize Error Flag
      iErr = 2000

! Initialize the Owner ID
      IJK_Pg_Owner = UNDEFINED_I

! Set the initial search region, a single cell.
      lS1 = ll1 + (lu1-ll1)/2 + 1; lE1 = lS1
      lS2 = ll2 + (lu2-ll2)/2 + 1; lE2 = lS2

      lpCnt = 1
      recheck = .TRUE.
      do while(recheck)

! Initialize the global IJK array to zero. Resetting this array inside
! this do-loop is most likely overkill. This loop should only cycle
! if gIJK is zero.
         gIJK = 0

! Report debugging information for the search region.
         if(ldFlag) then
            write(*,"(/5x,'Pass: ',I4)") lpCnt
            write(*,"( 5x,'lp2 bounds:',2(2x,I4))")lS2, lE2
            write(*,"( 5x,'lp1 bounds:',2(2x,I4))")lS1, lE1
         endif

         lp2: do lc2 = lS2, lE2
         lp1: do lc1 = lS1, lE1
! Map the loop counters to I/J/K indices.
            CALL lMAP(lc1, lc2, ll3, I, J, K)

! Only the rank that owns this I/J/K proceeds.
            if(.NOT.IS_ON_myPE_owns(I,J,K)) cycle
! Calculate the triple loop index.
            IJK = funijk(I,J,K)
! If there is fluid at this location, store the IJK and exit loops.
            if(fluid_at(IJK)) then
               gIJK(myPE) = IJK
               exit lp2
            endif
         enddo lp1
         enddo lp2

! Sync gIJK across all processes. Select the lowest ranked process that
! has gIJK defined. This choice is arbitray and doesn't really matter.
! It just needs to be consistent.
         CALL global_all_sum(gIJK)
         proc_lp: do proc=0, numPEs-1
            if(gIJK(proc) /= 0) then
               IJK_P_g = gIJK(proc)
               IJK_Pg_Owner = proc
               recheck = .FALSE.
               exit proc_lp
            endif
         enddo proc_lp

! If the proceeding section did not find a fluid cell, expand the search
! area and try again.
         if(recheck) then
            if(lS1 > ll1 .OR. lE1 < lu1 .OR.                           &
               lS2 > ll2 .OR. lE2 < lu2) then
! Expand the 1-axis
               lS1 = max((lS1-1), ll1)
               lE1 = min((lE1+1), lu1)
! Expand the 2-axis
               lS2 = max((lS2-1), ll2)
               lE2 = min((lE2+1), lu2)
! The entire seach plane was checked with no fluid cell identified.
! Force IJK_P_g to undefined for later error checking.
            else
               recheck = .FALSE.
               IJK_P_g = UNDEFINED_I
            endif
         endif
      enddo

! Verify that one fluid cell was detected. Otherwise flag the possible
! errors and return.
      if(IJK_P_G == UNDEFINED_I) then
         iErr = 2001
         return
      elseif(IJK_Pg_Owner == UNDEFINED_I) then
         iErr = 2002
         return
      endif


! The owner if the IJK_Pg gets the global I/J/K values and sends
! them to all ranks.
      if(myPE == IJK_Pg_Owner) then
         I_J_K_Pg(1) = I_OF(IJK_P_G)
         I_J_K_Pg(2) = J_OF(IJK_P_G)
         I_J_K_Pg(3) = K_OF(IJK_P_G)
      endif
      CALL BCAST(I_J_K_Pg, IJK_Pg_Owner)

      I = I_J_K_Pg(1)
      J = I_J_K_Pg(2)
      K = I_J_K_Pg(3)

! If debugging, have PE_IO report some information before the
! data is overwritten.
      if(ldFlag) then
         write(*,"(/3x,'IJK_P_g successfully identified!')")
         cInt=''; write(cInt,*) IJK_Pg_Owner
         write(*,"( 5x,'Owner Rank: ',A)")trim(adjustl(cInt))
         cInt=''; write(cInt,*) IJK_P_G
         write(*,"(5x, 'IJK: ',A)", advance='no') trim(adjustl(cInt))
         write(*,"(' :: ')", advance='no')
         cInt=''; write(cInt,*) I
         write(*,"('(',A)",advance='no') trim(adjustl(cInt))
         cInt=''; write(cInt,*) J
         write(*,"(',',A)",advance='no') trim(adjustl(cInt))
         cInt=''; write(cInt,*) K
         write(*,"(',',A,')',2/)") trim(adjustl(cInt))
      endif

! Ranks that 'see' IJK_P_g store their local IJK value. Everyone else
! resets IJK_P_g to UNDEFINED_I. This removes the need for getting
! I/J/K values later on in source_PPg.
!      IJK_P_g = merge(funijk(I,J,K), UNDEFINED_I,                      &
!         IS_ON_myPE_plus2layers(I,J,K))

      IF(IS_ON_myPE_plus2layers(I,J,K)) THEN
         IJK_P_g = funijk(I,J,K)
      ELSE
         IJK_P_g = UNDEFINED_I
      ENDIF

      IERR = 0
      RETURN
      END SUBROUTINE IJK_Pg_SEARCH


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: NegROs_LOG                                              !
!  Purpose: Record information about the location and conditions that  !
!           resulted in a negative solids phase density.               !
!                                                                      !
!  Author: J. Musser                                  Date: 09-Oct-13  !
!  Reviewer:                                          Date:            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE JKI_MAP(in1, in2, in3, lI, lJ, lK)

      implicit none

      INTEGER, intent(in) :: in1, in2, in3
      INTEGER, intent(out) :: lI, lJ, lK

      lI=in1; lJ=in3; lK=in2; return
      return
      END SUBROUTINE JKI_MAP



      SUBROUTINE IKJ_MAP(in1, in2, in3, lI, lJ, lK)

      implicit none

      INTEGER, intent(in) :: in1, in2, in3
      INTEGER, intent(out) :: lI, lJ, lK

      lI=in3; lJ=in1; lK=in2; return
      return
      END SUBROUTINE IKJ_MAP



      SUBROUTINE KIJ_MAP(in1, in2, in3, lI, lJ, lK)

      implicit none

      INTEGER, intent(in) :: in1, in2, in3
      INTEGER, intent(out) :: lI, lJ, lK

      lI=in2; lJ=in1; lK=in3; return
      return
      END SUBROUTINE KIJ_MAP
