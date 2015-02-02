MODULE ChiScheme
!  Purpose: Determine factors for Chi-scheme of Darwish and Moukalled
!           to ensure consistency of equation sets (e.g., species mass
!           fractions add up to 1)
!
!  Author: M. Syamlal                                 Date: 1-AUG-03
!  Reviewer:                                          Date:
!
!
!  Literature/Document References: Darwish, M. and Moukalled, F.,
!   "The Chi-shemes: a new consistent high-resolution formulation based on
!    the normalized variable methodology," Comput. Methods Appl. Mech. Engrg.,
!    192, 1711-1730 (2003)
!

! To initiate Chi-Scheme
!   call set_chi( ...)
!
! and to terminate Chi-Scheme
!   call unset_chi(ier)
!

      USE param
      USE param1
      USE run
      USE geometry
      USE indices
      USE compar
      USE sendrecv
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Chi_e
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Chi_n
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Chi_t


   LOGICAL :: ChiScheme_allocated = .false.
   LOGICAL :: Chi_flag = .false.

   CONTAINS
      SUBROUTINE Set_Chi(DISCR, PHI, Nmax, U, V, W)
!
!                      discretization method
      INTEGER          DISCR
!
!                      Second dimension of Phi array
      INTEGER          NMax
!
!                      convected quantity
      DOUBLE PRECISION PHI(DIMENSION_3, Nmax)
!
!                      Velocity components
      DOUBLE PRECISION U(DIMENSION_3), V(DIMENSION_3), W(DIMENSION_3)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Chi_e_temp
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Chi_n_temp
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Chi_t_temp
!
!                      index
      INTEGER          IJK, N

      if(.not.ChiScheme_allocated)then
        Allocate( Chi_e(DIMENSION_3) , &
                  Chi_n(DIMENSION_3) , &
                  Chi_t(DIMENSION_3)     )
        ChiScheme_allocated = .true.
      endif

      if(Chi_flag)then
        ! Error: Chi-Scheme is already active.  This routine cannot be called
        ! again before unsetting the flag
        Print *, 'Module ChiScheme: Cannot call Set_Chi again, before Unset_chi'
        call Mfix_Exit(0)

      else
      ! Set Chi_flag to indicate that future calls to calc_Xsi will use
      ! the Chi-Scheme for discretization
        Chi_e = large_number
        Chi_n = large_number
        Chi_t = large_number
        Chi_flag = .true.
      Endif

      Allocate( Chi_e_temp(DIMENSION_3) , &
                Chi_n_temp(DIMENSION_3) , &
                Chi_t_temp(DIMENSION_3)  )

!  Start Chi calculations
      DO N = 1, Nmax

        CALL CALC_CHI(DISCR, PHI(1,N), U, V, W, Chi_e_temp, Chi_n_temp, Chi_t_temp)

!!!$omp    parallel do private(IJK)
        DO IJK = ijkstart3, ijkend3
          Chi_e(IJK) = MIN(Chi_e(IJK), Chi_e_temp(IJK))
          Chi_n(IJK) = MIN(Chi_n(IJK), Chi_n_temp(IJK))
          Chi_t(IJK) = MIN(Chi_t(IJK), Chi_t_temp(IJK))
        ENDDO

      ENDDO

!  End Chi calculations

      call send_recv(CHI_E,2)
      call send_recv(CHI_N,2)
      call send_recv(CHI_T,2)

      Deallocate( Chi_e_temp , &
                  Chi_n_temp , &
                  Chi_t_temp  )


      RETURN
      END SUBROUTINE Set_Chi


      SUBROUTINE Unset_Chi()
        Chi_flag = .false.
      RETURN
      END SUBROUTINE Unset_Chi

!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine name: CALC_CHI(DISCR, PHI, U, V, W, CHI_e, CHI_n, CHI_t)     C
!  Purpose: Determine CHI factors for higher order discretization.
!  Author: M. Syamlal                                 Date: 4-AUG-03   C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References: Darwish and Moukalled (2003)
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CALC_CHI(DISCR, PHI, U, V, W, CHI_E, CHI_N, CHI_T)
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE run
      USE geometry
      USE indices
      USE vshear
      USE compar
      USE sendrecv
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

!                      discretization method
      INTEGER          DISCR
!
!                      convected quantity
      DOUBLE PRECISION PHI(DIMENSION_3)
!
!                      Velocity components
      DOUBLE PRECISION U(DIMENSION_3), V(DIMENSION_3), W(DIMENSION_3)
!
!                      Convection weighting factors
      DOUBLE PRECISION CHI_e(DIMENSION_3), CHI_n(DIMENSION_3),&
                       CHI_t(DIMENSION_3)
!
!                      Indices
      INTEGER          IJK, IJKC, IJKD, IJKU

!
!                      Error message
      CHARACTER(LEN=80)     LINE(1)
!
!
      DOUBLE PRECISION PHI_C
!
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      DOUBLE PRECISION , EXTERNAL :: PHI_C_OF, CHI4SMART, CHI4MUSCL

        IF (SHEAR) THEN
! calculate CHI_E,CHI_N,CHI_T when periodic shear BCs are used

!       call CXS(incr,DISCR,U,V,W,PHI,CHI_E,CHI_N,CHI_T)  !need implementation
         print *,'From CALC_CHI:  "Shear" option not implemented'
         Call MFIX_EXIT(0)


        ELSE
!
!
      SELECT CASE (DISCR)                        !first order upwinding
      CASE (:1)
!
!!!$omp    parallel do private(IJK)
         DO IJK = ijkstart3, ijkend3
            CHI_E(IJK) = ZERO
            CHI_N(IJK) = ZERO
            IF (DO_K) CHI_T(IJK) = ZERO
         END DO
!      CASE (2)                                   !Superbee
      CASE (3)                                   !SMART
!
!!!$omp    parallel do private(IJK, IJKC,IJKD,IJKU, PHI_C)
         DO IJK = ijkstart3, ijkend3
          IF (.NOT.WALL_AT(IJK)) THEN ! no need to do these calculations for walls
            IF (U(IJK) >= ZERO) THEN
               IJKC = IJK
               IJKD = EAST_OF(IJK)
               IJKU = WEST_OF(IJKC)
            ELSE
               IJKC = EAST_OF(IJK)
               IJKD = IJK
               IJKU = EAST_OF(IJKC)
            ENDIF
!
            PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD))
!
            CHI_E(IJK) = CHI4SMART(PHI_C,PHI(IJKU),PHI(IJKC),PHI(IJKD))
!
            IF (V(IJK) >= ZERO) THEN
               IJKC = IJK
               IJKD = NORTH_OF(IJK)
               IJKU = SOUTH_OF(IJKC)
            ELSE
               IJKC = NORTH_OF(IJK)
               IJKD = IJK
               IJKU = NORTH_OF(IJKC)
            ENDIF
!
            PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD))
!
            CHI_N(IJK) = CHI4SMART(PHI_C,PHI(IJKU),PHI(IJKC),PHI(IJKD))
!
            IF (DO_K) THEN
               IF (W(IJK) >= ZERO) THEN
                  IJKC = IJK
                  IJKD = TOP_OF(IJK)
                  IJKU = BOTTOM_OF(IJKC)
               ELSE
                  IJKC = TOP_OF(IJK)
                  IJKD = IJK
                  IJKU = TOP_OF(IJKC)
               ENDIF
!
               PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD))
!
               CHI_T(IJK) = CHI4SMART(PHI_C,PHI(IJKU),PHI(IJKC),PHI(IJKD))
            ENDIF
          ELSE
            CHI_E(IJK) = ZERO
            CHI_N(IJK) = ZERO
            CHI_T(IJK) = ZERO
          ENDIF
         END DO
!      CASE (4)                                   !ULTRA-QUICK
!      CASE (5)                                   !QUICKEST
      CASE (6)                                   !MUSCL

!!!$omp    parallel do private(IJK, IJKC,IJKD,IJKU, PHI_C )
         DO IJK = ijkstart3, ijkend3
          IF (.NOT.WALL_AT(IJK)) THEN ! no need to do these calculations for walls
            IF (U(IJK) >= ZERO) THEN
               IJKC = IJK
               IJKD = EAST_OF(IJK)
               IJKU = WEST_OF(IJKC)
            ELSE
               IJKC = EAST_OF(IJK)
               IJKD = IJK
               IJKU = EAST_OF(IJKC)
            ENDIF
!
            PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD))
!
            CHI_E(IJK) = CHI4MUSCL(PHI_C,PHI(IJKU),PHI(IJKC),PHI(IJKD))
!
            IF (V(IJK) >= ZERO) THEN
               IJKC = IJK
               IJKD = NORTH_OF(IJK)
               IJKU = SOUTH_OF(IJKC)
            ELSE
               IJKC = NORTH_OF(IJK)
               IJKD = IJK
               IJKU = NORTH_OF(IJKC)
            ENDIF
!
            PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD))
!
            CHI_N(IJK) = CHI4MUSCL(PHI_C,PHI(IJKU),PHI(IJKC),PHI(IJKD))
!
            IF (DO_K) THEN
               IF (W(IJK) >= ZERO) THEN
                  IJKC = IJK
                  IJKD = TOP_OF(IJK)
                  IJKU = BOTTOM_OF(IJKC)
               ELSE
                  IJKC = TOP_OF(IJK)
                  IJKD = IJK
                  IJKU = TOP_OF(IJKC)
               ENDIF
!
               PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD))
!
               CHI_T(IJK) = CHI4MUSCL(PHI_C,PHI(IJKU),PHI(IJKC),PHI(IJKD))
            ENDIF
          ELSE
            CHI_E(IJK) = ZERO
            CHI_N(IJK) = ZERO
            CHI_T(IJK) = ZERO
          ENDIF ! for walls
         END DO
!      CASE (7)                                   !Van Leer
!      CASE (8)                                   !Minmod
      CASE DEFAULT                               !Error
         WRITE (LINE, '(A,I2,A)') 'Chi-Scheme for DISCRETIZE = ', DISCR, ' not supported.'
         CALL WRITE_ERROR ('CALC_CHI', LINE, 1)
         CALL MFIX_EXIT(myPE)
      END SELECT

      ENDIF

      call send_recv(CHI_E,2)
      call send_recv(CHI_N,2)
      call send_recv(CHI_T,2)

      RETURN
      END SUBROUTINE CALC_CHI

END MODULE CHIScheme

