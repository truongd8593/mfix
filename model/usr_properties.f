!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: USR_PROPERTIES                                          C
!  Purpose: Hook for user defined physical properties and transport    C
!  coefficients including interphase exchange coefficients for all     C
!  phases. The exact quantities are listed here for clarity.           C
!                                                                      C
!         Quantity      GAS (M=0)     SOLIDS (M=1toMMAX)               C
!         density       RO_G(:)        RO_S(:,M)                       C
!      specific heat    C_PG(:)        C_PS(:,M)                       C
!       conductivity    K_G(:)         K_S(:,M)                        C
!       diffusivity     DIF_G(:,N)     DIF_S(:,M,N)                    C
!        viscosity      MU_G(:)        MU_S(:,M)                       C
!     gas-solids drag                  F_GS(:,M)                       C
!   solids-solids drag                 F_SS(:,LM)                      C
!   gas-solids heat tr.                GAMA(:,M)                       C
!                                                                      C
!  Comments:                                                           C
!  - gas-solids drag is momentum transfer due to relative velocity     C
!    differences (skin friction and form drag) between the gas phase   C
!    (M=0) and each solids phase (M=1 to MMAX).                        C
!  - solids-solids drag is momentum transfer due to relative velocity  C
!    differences between solids phases M and L, where M and L range    C
!    from 1 to MMAX and M!=L.                                          C
!  - gas-solids heat transfer is heat transfer due to relative         C
!    temperature differences between the gas phase phase (M=0) and     C
!    each solids phase (M=1 to MMAX).                                  C
!  - No solids-solids heat transfer is allowed. To account for such    C
!    would require appropriate closures and additional code            C
!    development.                                                      C
!  - the specific heat assigned in this routine only applies to the    C
!    mixture average specific heat invoked in the general energy       C
!    equations. reacting flow simulations require values for species   C
!    specific heat and are taken from the Burcat database or read      C
!    in that format from the mfix.dat file. therefore inconsitencies   C
!    may arise in calculations involving species specific heat. this   C
!    requires further development to fully address.                    C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE USR_PROPERTIES(lprop, IJK, M, N)

! Modules 
!-----------------------------------------------
      use constant, only: pi, gas_const, gravity
      use error_manager

      use fldvar, only: u_g, v_g, w_g, ep_g
      use fldvar, only: u_s, v_s, w_s, ep_s
      use fldvar, only: p_g, rop_g, ro_g, T_g, X_g
      use fldvar, only: p_s, rop_s, ro_s, T_s, X_s
      use fldvar, only: d_p, theta_m

      use fldvar, only: scalar
      use functions
      use geometry

      use indices, only: i_of, j_of, k_of
      use indices, only: im1, ip1, jm1, jp1, km1, kp1

      use param1, only: zero, one, half, undefined, undefined_i
      use physprop, only: k_g, c_pg, dif_g, mu_g
      use physprop, only: k_s, c_ps, dif_s
      use scalars, only: phase4scalar
      use visc_s, only: mu_s, lambda_s
      use usr_prop
      IMPLICIT NONE

! Dummy arguments
!-----------------------------------------------
! reference equation 
      INTEGER, INTENT(IN) :: lprop
! index
      INTEGER, INTENT(IN) :: IJK
! Phase index  
      INTEGER, INTENT(IN) :: M
! Species index or second solids phase index for solids-solids drag
! (if applicable otherwise undefined_i)
      INTEGER, INTENT(IN) :: N

! Local variables
!-----------------------------------------------
! Error flag
      INTEGER :: IER
      CHARACTER(len=40):: err_prop

!-----------------------------------------------
! initialize
      ier = undefined_i

! in each case the ier flag is set to ensure that if a user defined
! quantity is invoked that the associated user defined quantity is
! actually specified. this flag must be removed by the user or the
! code will fail.

      SELECT CASE(lprop)

! 
      CASE (GAS_DENSITY)
! if using this quantity then remove definition of ier
         ier = 1
         write(err_prop, '("gas density")')
!        RO_G(IJK) =


      CASE (SOLIDS_DENSITY)
! if using this quantity then remove definition of ier
         ier = 1
         write(err_prop, '("solids phase ",I2," density")') M
!        RO_S(IJK,M) =

      CASE (GAS_SPECIFICHEAT)
! if using this quantity then remove definition of ier
         ier = 1
         write(err_prop, '("gas specific heat")')
!        C_PG(IJK,M) = 


      CASE (SOLIDS_SPECIFICHEAT)
! if using this quantity then remove definition of ier
         ier = 1
         write(err_prop, '("solids phase ",I2," specific heat")') M
!        C_PS(IJK,M) = 


      CASE (GAS_CONDUCTIVITY)
! if using this quantity then remove definition of ier
         ier = 1
         write(err_prop, '("gas conductivity")')
!        K_G(IJK) =  


      CASE (SOLIDS_CONDUCTIVITY)
! if using this quantity then remove definition of ier
         ier = 1
         write(err_prop, '("solids phase ",I2," conductivity")') M
!        K_S(IJK,M) =  


      CASE (GAS_DIFFUSIVITY)
! if using this quantity then remove definition of ier
         ier = 1
         write(err_prop, '("gas diffusivity")')
!        DIF_G(IJK,N) = 

      CASE (SOLIDS_DIFFUSIVITY)
! if using this quantity then remove definition of ier
         ier = 1
         write(err_prop, '("solids phase",I2," diffusivity")') M
!        DIF_S(IJK,M,N) = 


! assign gas viscosity value mu_g. bulk viscosity is taken as zero.
! second viscosity (lambda_g) is automatically defined as
! lambda_g = -2/3mu_g
      CASE (GAS_VISCOSITY)
! if using this quantity then remove definition of ier
         ier = 1
         write(err_prop, '("gas viscosity")')
!        MU_G(IJK) = 


! assign solids phase M viscosity mu_s. 
! assign second viscosity (lambda_s): lambda_s = mu_sbulk - 2/3mu_s
! assign solids pressure p_s. 
      CASE (SOLIDS_VISCOSITY)
! if using this quantity then remove definition of ier
         ier = 1
         write(err_prop, '("solids phase ", I2," viscosity")') M
!        MU_S(IJK,M) =
!        LAMBDA_S(IJK,M) = 
!        P_S(IJK,M) = 


      CASE (GASSOLIDS_HEATTRANSFER)
! if using this quantity then remove definition of ier
         ier = 1
         write(err_prop, &
               '("gas-solids phase ",I2," heat transfer coefficient")') M
!       GAMA_GS(IJK,M) = 


      CASE (GASSOLIDS_DRAG)
! if using this quantity then remove definition of ier
         ier = 1
         write(err_prop, &
               '("gas-solids phase ",I2," drag coefficient")') M
!        F_GS(IJK,M) =


      CASE (SOLIDSSOLIDS_DRAG)
! if using this quantity then remove definition of ier
         ier = 1
         write(err_prop, &
            '("solids-solids phases ",I2," & ",I2," drag coefficient")') M, N
!        F_SS(IJK,LM) = 


      END SELECT

      IF (IER /= UNDEFINED_I) THEN
         CALL INIT_ERR_MSG('USR_PROPERTIES') 
         WRITE(ERR_MSG,9999) trim(err_prop)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

 9999 FORMAT('ERROR 9999: The user-defined properties routine was ',&
         'invoked for',/,A,' but this generic error',/,'message exits.',&
         'Either choose a different model or correct',/,'mfix/model/',&
         'usr_properties.f')
      ENDIF

      RETURN
      END SUBROUTINE USR_PROPERTIES



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: USR_PROP_ROg                                            !
!  Purpose: User hook for calculating the gas phase density.           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE USR_PROP_ROg(IJK)

! Modules
!---------------------------------------------------------------------//
      use error_manager, only: init_err_msg, err_msg, flush_err_msg
      use fldvar, only: t_g, x_g, p_g, ro_g
      use param1, only: undefined_i, zero, one, half
      use physprop, only: mw_g, mw_avg, nmax
      implicit none

! Dummy arguments
!---------------------------------------------------------------------//
! index
      INTEGER, INTENT(IN) :: IJK

! Local Variables:
!---------------------------------------------------------------------//
! error flag
      INTEGER :: IER = undefined_i
      CHARACTER(LEN=40) :: err_prop
!......................................................................!

! if using this quantity then remove definition of ier
      ier = 1

! Assign the fluid density 
      RO_G(IJK) = ZERO


      IF (IER /= UNDEFINED_I) THEN
         write(err_prop, '("gas density")')
         CALL INIT_ERR_MSG('USR_PROP_ROg') 
         WRITE(ERR_MSG,9999) trim(err_prop)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
 9999 FORMAT('ERROR 9999: The user-defined properties routine was ',&
         'invoked for',/,A,' but this generic error',/,'message exi',&
         'sts. Either choose a different model or correct',/,'mfix/,'&
         'model/usr_properties.f')
      ENDIF

      RETURN
      END SUBROUTINE USR_PROP_ROg


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: USR_PROP_ROs                                            !
!  Purpose: User hook for calculating solids phase density.            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE USR_PROP_ROs(IJK,M)

! Modules
!---------------------------------------------------------------------//
      use error_manager, only: init_err_msg, err_msg, flush_err_msg
      use fldvar, only: ro_s, T_s, X_s
      use param1, only: undefined_i, zero, one, half
      use physprop, only: mw_s, nmax
      implicit none

! Dummy arguments
!---------------------------------------------------------------------//
! index
      INTEGER, INTENT(IN) :: IJK
! solids phase index
      INTEGER, INTENT(IN) :: M

! Local Variables:
!---------------------------------------------------------------------//
! error flag
      INTEGER :: IER = undefined_i
      CHARACTER(LEN=40) :: err_prop
!......................................................................!
! if using this quantity then remove definition of ier
      ier = 1


! Assign the fluid density 
      RO_S(IJK,M) = ZERO


      IF (IER /= UNDEFINED_I) THEN
         write(err_prop, '("solids phase ",I2," density")') M
         CALL INIT_ERR_MSG('USR_PROP_ROS') 
         WRITE(ERR_MSG,9999) trim(err_prop)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
 9999 FORMAT('ERROR 9999: The user-defined properties routine was ',&
         'invoked for',/,A,' but this generic error',/,'message exi',&
         'sts. Either choose a different model or correct',/,'mfix/,'&
         'model/usr_properties.f')
      ENDIF

      RETURN
      END SUBROUTINE USR_PROP_ROs


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: USR_PROP_CPg                                            !
!  Purpose: User hook for calculating the gas phase constant pressure  !
!  specific heat.                                                      !
!                                                                      !
!  Comments:                                                           !
!  - The specific heat assigned in this routine only applies to the    !
!    mixture average specific heat invoked in the general energy       !
!    equations. MFIX has no global variable representing species       !
!    specific heats. The reason for this limitation is explained.      !
!                                                                      !
!    Species specific heat values are locally evaluated based on a     !
!    specific polynominal format set by the Burcat database. Values    !
!    for the polynominal coefficients are read from either the         !
!    database or the mfix.dat. Species specific heats are needed by    !
!    reacting flow simulations.                                        !
!                                                                      !
!    Inconsistencies may arise in reacting flow systems if a user      !
!    specifies a phase average specific heat that is not consistent    !
!    with the values of the species specific heats that comprise that  !
!    phase.                                                            !
!    - IT may be possible to circumvent the matter long as the formula !
!      for the specific heat follows the burcat database form and      !
!      the quantites Thigh(M,N), Tlow(M,N), Tcom, Alow(M,N) and        !
!      Ahigh(M,N) are all appropriately assigned. However, the same    !
!      can be achieved by simply entering the data into the mfix.dat   !
!      as indicated by the user guide.                                 !
!    - To permit different forms of the polynominal for species        !
!      specific heats would require more code development to ensure    !
!      reacting flow simulations evaluate/reference the species        !
!      specific heats accordingly.                                     !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE USR_PROP_CPg(IJK)

! Modules
!---------------------------------------------------------------------//
      use error_manager, only: init_err_msg, err_msg, flush_err_msg
      use constant, only: RGAS => GAS_CONST_cal
      use fldvar, only: t_g, x_g
      use param1, only: undefined_i, zero, one, half
      use physprop, only: mw_g, c_pg, nmax
      use read_thermochemical, only: calc_CpoR
      use run, only: UNITS
      implicit none

! Dummy arguments
!---------------------------------------------------------------------//
! index
      INTEGER, INTENT(IN) :: IJK

! Local Variables:
!---------------------------------------------------------------------//
! error flag
      INTEGER :: IER = undefined_i
      CHARACTER(LEN=40) :: err_prop
!......................................................................!

! if using this quantity then remove definition of ier
      ier = 1

! Assign the fluid density 
      C_PG(IJK) = ZERO


      IF (IER /= UNDEFINED_I) THEN
         write(err_prop, '("gas specific heat")')
         CALL INIT_ERR_MSG('USR_PROP_CPg') 
         WRITE(ERR_MSG,9999) trim(err_prop)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
 9999 FORMAT('ERROR 9999: The user-defined properties routine was ',&
         'invoked for',/,A,' but this generic error',/,'message exi',&
         'sts. Either choose a different model or correct',/,'mfix/,'&
         'model/usr_properties.f')
      ENDIF

      RETURN
      END SUBROUTINE USR_PROP_CPg


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: USR_PROP_CPs                                            !
!  Purpose: User hook for calculating solids phase constant pressure   !
!  specific heat.                                                      !
!                                                                      !
!  Comments:                                                           !
!  - See comments under USER_PROP_CPg                                  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE USR_PROP_CPs(IJK, M)
! Modules
!---------------------------------------------------------------------//
      use error_manager, only: init_err_msg, err_msg, flush_err_msg
      use constant, only: RGAS => GAS_CONST_cal
      use fldvar, only: t_s, x_s
      use param1, only: undefined_i, zero, one, half
      use physprop, only: mw_s, c_ps, nmax
      use read_thermochemical, only: calc_CpoR
      use run, only: UNITS
      implicit none

! Dummy arguments
!---------------------------------------------------------------------//
! index
      INTEGER, INTENT(IN) :: IJK
! solids phase index
      INTEGER, INTENT(IN) :: M

! Local Variables:
!---------------------------------------------------------------------//
! error flag
      INTEGER :: IER = undefined_i
      CHARACTER(LEN=40) :: err_prop
!......................................................................!

! if using this quantity then remove definition of ier
      ier = 1

! Assign the fluid density 
      C_PS(IJK,M) = ZERO


      IF (IER /= UNDEFINED_I) THEN
         write(err_prop, '("solids phase ",I2," specific heat")') M
         CALL INIT_ERR_MSG('USR_PROP_CPs') 
         WRITE(ERR_MSG,9999) trim(err_prop)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
 9999 FORMAT('ERROR 9999: The user-defined properties routine was ',&
         'invoked for',/,A,' but this generic error',/,'message exi',&
         'sts. Either choose a different model or correct',/,'mfix/,'&
         'model/usr_properties.f')
      ENDIF

      RETURN
      END SUBROUTINE USR_PROP_CPs
