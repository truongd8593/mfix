!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module: kintheory                                                   C
!  Purpose: Common block containing constants, variables, functions    C
!  used by various kinetic theory models                               C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Garzo, V., Tenneti, S., Subramaniam, S., and Hrenya, C. M.,         C
!      "Enskog kinetic theory for monodisperse gas-solid flows", JFM,  C
!      Vol. 712, 2012, pp. 129-168                                     C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      MODULE kintheory


! coefficient terms needed for stress (in addition to
! mu_s and lambda_s which are defined in visc_s_mod, allocated
! in allocate_arrays and initialized in set_constprop; and
! P_s which is defined in fldvar_mod, allocated in
! allocate_arrays, and initialized in init_fvars)
!     stress term with gradient in particle M velocity
      DOUBLE PRECISION, DIMENSION(:, :, :), ALLOCATABLE :: MU_sM_ip
!     stress term with gradient in particle L velocity
      DOUBLE PRECISION, DIMENSION(:, :, :), ALLOCATABLE :: MU_sL_ip
!     stress term with trace in particle M velocity
      DOUBLE PRECISION, DIMENSION(:, :, :), ALLOCATABLE :: XI_sM_ip
!     stress term with trace in particle L velocity
      DOUBLE PRECISION, DIMENSION(:, :, :), ALLOCATABLE :: XI_sL_ip


! coefficient terms needed for momentum source (in addition to
! F_SS which is defined in drag_mod, allocated in allocate_arrays
! and initialized in set_constprop)
!     momentum source term with gradient in number density
      DOUBLE PRECISION, DIMENSION(:, :, :), ALLOCATABLE :: Fnu_s_ip
!     momentum source term with gradient in mixture temperature or
!     with the gradient in temperature of species M
      DOUBLE PRECISION, DIMENSION(:, :, :), ALLOCATABLE :: FT_sM_ip
!     momentum source term with gradient in temperature of species L
      DOUBLE PRECISION, DIMENSION(:, :, :), ALLOCATABLE :: FT_sL_ip


! coefficient terms needed for heat flux (in addition to
! kth_s which is defined in set_constprop, allocated
! allocate_arrays and initialized in init_fvars)
!     heat flux term with gradient in granular temperature of species L
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: Kth_sL_ip
!     heat flux term with gradient in number density of species M
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: Knu_sM_ip
!     heat flux term with gradient in number density of species L
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: Knu_sL_ip
!     heat flux term with velocity difference
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: Kvel_s_ip


! coefficient terms needed for energy dissipation
!     energy dissipation with difference in species granular
!     temperature: transfer between solid solid phases
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: ED_ss_ip
!     energy dissipation term
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: EDT_s_ip
!     energy dissipation with divergence of velocity of species M
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: EDvel_sM_ip
!     energy dissipation with divergence of velocity of species L
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: EDvel_sL_ip
!     coefficient A2, xsi used in multiple places in GTSH theory
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: A2_gtsh
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: xsi_gtsh

! Solids source terms needed for Iddir & Arastoopour (2005)
! kinetic theory model
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  KTMOM_U_s
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  KTMOM_V_s
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  KTMOM_W_s

! parameter in the theory of GTSH that is related to length scale
! of lubrication effects. For details see GTSH, 2012.
      DOUBLE PRECISION, PARAMETER :: EpM = 0.01d0



      CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Purpose: Calculate the magnitude of the gas-solids relative         C
!  velocity at i, j, k                                                 C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      DOUBLE PRECISION FUNCTION KT_RVEL(IJK, m)

! Modules
!---------------------------------------------------------------------//
      use fldvar, only: u_s, v_s, w_s
      use fldvar, only: u_g, v_g, w_g
      use run, only: shear
      use vshear, only: vsh
      use indices, only: i_of
      use functions, only: im_of, jm_of, km_of
      use functions, only: fluid_at
      use fun_avg, only: avg_x_e, avg_y_n, avg_z_t
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! ijk index
      INTEGER, INTENT(IN) :: ijk
! solids phase index
      INTEGER, INTENT(IN) :: M

! Local variables
!---------------------------------------------------------------------//
! cell indices
      INTEGER :: I, IMJK, IJMK, IJKM
! Cell center value of solids and gas velocities
      DOUBLE PRECISION :: USCM, VSCM, WSCM, &
                          UGC, VGC, WGC
! y-component of velocity with 'shear' applied
      DOUBLE PRECISION :: vs_j, vs_jm
!---------------------------------------------------------------------//
      I = I_OF(IJK)
      IMJK  = IM_OF(IJK)
      IJMK  = JM_OF(IJK)
      IJKM  = KM_OF(IJK)

! Awkward here but captures what was done in early version of
! calc_mu_s. Otherwise any call to this routine must be done with
! 'shear' already applied to v_s.
      vs_j = V_S(IJK,M)
      vs_jm = V_S(IJMK,M)
      IF (SHEAR) THEN
         vs_j = V_S(IJK,M)+VSH(IJK)
         IF(FLUID_AT(IJMK)) THEN
            vs_jm = V_S(IJMK,M)+VSH(IJMK)
         ELSE
            vs_jm = v_s(IJMK,M)
         ENDIF
      ENDIF

! Start calculation for relative velocity
! Calculate velocity components at i, j, k
      UGC = AVG_X_E(U_G(IMJK),U_G(IJK),I)
      VGC = AVG_Y_N(V_G(IJMK),V_G(IJK))
      WGC = AVG_Z_T(W_G(IJKM),W_G(IJK))

      USCM = AVG_X_E(U_S(IMJK,M),U_S(IJK,M),I)
      VSCM = AVG_Y_N(vs_jm, vs_j)
      WSCM = AVG_Z_T(W_S(IJKM,M),W_S(IJK,M))

! Magnitude of gas-solids relative velocity
      KT_RVEL = SQRT((UGC - USCM)**2 + (VGC - VSCM)**2 + &
                     (WGC - WSCM)**2)

      RETURN
      END FUNCTION KT_RVEL

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      DOUBLE PRECISION FUNCTION KT_COS_THETA(ijk, M)

! Modules
!---------------------------------------------------------------------//
      USE param1, only: zero, small_number
      USE fldvar, only: u_g, v_g, w_g
      USE fldvar, only: u_s, v_s, w_s
      USE fldvar, only: ep_s
      USE run, only: shear
      USE vshear, only: vsh
      USE toleranc, only: zero_ep_s
      USE indices, only: i_of
      USE functions, only: im_of, jm_of, km_of
      USE functions, only: fluid_at
      USE fun_avg, only: avg_x_e, avg_y_n, avg_z_t
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! ijk index
      INTEGER, INTENT(IN) :: ijk
! solids phase index
      INTEGER, INTENT(IN) :: M

! Local variables
!---------------------------------------------------------------------//
! cell indices
      INTEGER :: I, IMJK, IJMK, IJKM
! Cell center value of solids and gas velocities
      DOUBLE PRECISION :: USCM, VSCM, WSCM, &
                          UGC, VGC, WGC
      DOUBLE PRECISION :: vs_j, vs_jm, speed, rvel
!---------------------------------------------------------------------//

      I = I_OF(IJK)
      IMJK  = IM_OF(IJK)
      IJMK  = JM_OF(IJK)
      IJKM  = KM_OF(IJK)

! Awkward here but captures what was done in early version of
! calc_mu_s. Otherwise any call to this routine must be done with
! 'shear' already applied to v_s.
      vs_j = V_S(IJK,M)
      vs_jm = V_S(IJMK,M)
      IF (SHEAR) THEN
         vs_j = V_S(IJK,M)+VSH(IJK)
         IF(FLUID_AT(IJMK)) THEN
            vs_jm = V_S(IJMK,M)+VSH(IJMK)
         ELSE
            vs_jm = v_s(IJMK,M)
         ENDIF
      ENDIF

      UGC = AVG_X_E(U_G(IMJK),U_G(IJK),I)
      VGC = AVG_Y_N(V_G(IJMK),V_G(IJK))
      WGC = AVG_Z_T(W_G(IJKM),W_G(IJK))

      USCM = AVG_X_E(U_S(IMJK,M),U_S(IJK,M),I)
      VSCM = AVG_Y_N(vs_jm, vs_j)
      WSCM = AVG_Z_T(W_S(IJKM,M),W_S(IJK,M))

      RVEL =  SQRT((UGC - USCM)**2 + (VGC - VSCM)**2 + &
                   (WGC - WSCM)**2)

      SPEED = SQRT(USCM**2+VSCM**2+WSCM**2)

! motion viewed by the particles (crossing trajectory effect)
      IF(SPEED > Small_Number .AND. RVEL > Small_Number .AND. &
         EP_S(IJK,M) > ZERO_EP_S) THEN
         KT_Cos_Theta = ( (UGC-USCM)*USCM + (VGC-VSCM)*VSCM + &
                       (WGC-WSCM)*WSCM )/(RVEL - SPEED)
      ELSE
         KT_Cos_Theta = ZERO
      ENDIF

      RETURN
      END FUNCTION  KT_COS_THETA


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Purpose: Calculate a single particle drag coefficient/ep_s, which   C
!  is used in evaluting limiting values of certain granular kinetic    C
!  theory terms                                                        C
!                                                                      C
!  Comments: This is currently based on wen-yu but in the the future   C
!  we may want to make this consistent with drag_gs (that is, replace  C
!  this function call and assign appropriate variable within drag_gs   C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      DOUBLE PRECISION FUNCTION KT_DGA(IJK, m)

! Modules
!---------------------------------------------------------------------//
      use param1, only: zero, one, small_number, large_number
      use fldvar, only: rop_g
      use fldvar, only: d_p
      use physprop, only: mu_g
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! ijk index
      INTEGER, INTENT(IN) :: IJK
! solids phase index
      INTEGER, INTENT(IN) :: M

! Local variables
!---------------------------------------------------------------------//
! single particle drag coefficient, reynolds number
      DOUBLE PRECISION :: C_d, Re
! local value for relative velocity
      DOUBLE PRECISION :: rvel
!---------------------------------------------------------------------//

! initialization
      kt_dga = zero
      rvel = zero
      RVEL = KT_RVEL(IJK, M)

! Defining single particle drag coefficient dgA based on Wen-Yu
! correlation
      RE = D_p(IJK,M)*RVEL*ROP_G(IJK)/(MU_G(IJK) + SMALL_NUMBER)
      IF(RE .LE. 1000.d0)THEN
         C_d = (24.d0/(Re+SMALL_NUMBER)) * &
            (ONE + 0.15d0 * Re**0.687D0)
      ELSE
         C_d = 0.44d0
      ENDIF

! dga_s is local to this routine as is lrvel
      kt_dga = 0.75d0 * C_d * RVEL * ROP_g(IJK) / D_p(IJK,M)

! set value for 1st iteration and 1st time step
      IF(RVEL == ZERO) kt_dga = LARGE_NUMBER

      RETURN
      END FUNCTION KT_DGA

      END MODULE kintheory
