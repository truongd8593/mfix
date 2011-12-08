!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: DES_RADIATION                                          !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!  Author: J.Musser                                   Date: 25-Jun-10  !
!                                                                      !
!  CommenT_ENV:                                                           !
!                                                                      !
!  REF: Batchelor and O'Brien, "Thermal or electrical conduction       !
!       through a granular material," Proceedings of the Royal Society !
!       of London. Series A, Mathematical and Physical Sciences,       !
!       Vol. 355, no 1682, pp. 313-333, July 1977.                     !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DES_RADIATION(I, T_ENV, FOCUS)

      Use constant
      Use des_thermo
      Use discretelement
      Use param1

      IMPLICIT NONE

! Index value of particle I
      INTEGER, INTENT(IN) :: I
! Average solids temperature in radiation domain
      DOUBLE PRECISION, INTENT(IN) :: T_ENV
! Logical used for debugging
      LOGICAL, INTENT(IN) :: FOCUS

!-----------------------------------------------
! Local variables
!-----------------------------------------------
! file name
      CHARACTER*64 :: FNAME
! logical used for testing is the data file already exists
      LOGICAL :: F_EXISTS
! Surface area of particle
      DOUBLE PRECISION A_S

! Calculate the surface area of the particle
      A_S = 4.0d0 * Pi * DES_RADIUS(I)**2

      Qrd(I) = SB_CONST * A_s * DES_EM(PIJK(I,5)) * &
         (T_ENV**4 - (DES_T_s_NEW(I))**4)

      IF(DEBUG_DES .AND. FOCUS)THEN
         WRITE(*,"(//5X,A)")'From: DES_RADIATION -'
         WRITE(*,"(8X,A,D12.6)")'Tp: ',DES_T_s_NEW(I)
         WRITE(*,"(8X,A,D12.6)")'T_ENV: ',T_ENV
         WRITE(*,"(8X,A,D12.6)")'Qrd: ',Qrd(I)
         WRITE(*,"(5X,25('-')/)")
      ENDIF

      RETURN

      END SUBROUTINE DES_RADIATION
