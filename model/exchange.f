!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  SUBROUTINE: EXCHANGE                                                C
!  Purpose: Calls routines to drive calculations of the interphase     C
!           mass, momentum, and energy exchange coefficients/terms     C
!           if directed to do so by the corresponding flags            C
!                                                                      C
!  Author: M. Syamlal                                 Date: 25-APR-96  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE EXCHANGE(IER)

! Global Variables
!-----------------------------------------------------------------------
! Flags for calculating drag coefficient.
      use coeff, only: DRAGCOEF
! Flags for calculating heat transfer coefficient.
      use coeff, only: HEAT_TR
      use physprop, only: smax
      use run, only: granular_energy
      use run, only: kt_type_enum, ia_2005
      implicit none

! Dummy arguments
!-----------------------------------------------------------------------
      INTEGER, INTENT(INOUT) :: IER ! Error index

! Local variables 
!-----------------------------------------------------------------------
! loop counter
      INTEGER :: M, L
!-----------------------------------------------------------------------

! Calculate drag coefficients
      CALL CALC_DRAG (DRAGCOEF, IER)

! Calculate additional interphase interaction coefficient
      IF (GRANULAR_ENERGY) THEN
         SELECT CASE(KT_TYPE_ENUM)
            CASE(IA_2005)
               DO M=1,SMAX
                  DO L=1,SMAX
                     CALL COLL_MOMENTUM_COEFF_IA(L, M)
                  ENDDO
               ENDDO
            CASE DEFAULT
         END SELECT
      ENDIF

! Calculate interphase heat transfer coefficients
      CALL CALC_GAMA (HEAT_TR)


      return
      END SUBROUTINE EXCHANGE
