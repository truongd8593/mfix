!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: DISPLAY_RESID(NIT, IER)                             C
!  Purpose: Display residuals                                          C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 8-JUL-96   C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE DISPLAY_RESID(NIT, IER)
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE physprop
      USE residual
      USE fldvar
      USE compar
      USE geometry
      USE scalars, only : NScalar
      USE run
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!
!                      iteration number
      INTEGER          NIT
!
!                      Error index
      INTEGER          IER
!
!                      residual number
      INTEGER          L
!
!                      Print Location of Max_Resid
      LOGICAL,PARAMETER:: Print_ijk=.FALSE.
!

!
      if (myPE.ne.PE_IO) return
!
      IF(GROUP_RESID) THEN

         IF (NIT == 1) THEN
            WRITE (*, '(A)',ADVANCE="NO") '  Nit'
            WRITE (*, '(2X, A7)',ADVANCE="NO") RESID_GRP_STRING(HYDRO_GRP)
            IF(GRANULAR_ENERGY) WRITE (*, '(2X, A7)',ADVANCE="NO") RESID_GRP_STRING(THETA_GRP)
            IF(ENERGY_EQ) WRITE (*, '(2X, A7)',ADVANCE="NO") RESID_GRP_STRING(ENERGY_GRP)
            IF(ANY_SPECIES_EQ) WRITE (*, '(2X, A7)',ADVANCE="NO") RESID_GRP_STRING(SPECIES_GRP)
            IF(NScalar > 0) WRITE (*, '(2X, A7)',ADVANCE="NO") RESID_GRP_STRING(SCALAR_GRP)
            IF(K_EPSILON) WRITE (*, '(2X, A7)',ADVANCE="NO") RESID_GRP_STRING(KE_GRP)
            WRITE (*, '(2X, A)') 'Max res'
         ENDIF

         WRITE (*, '(I5)',ADVANCE="NO") NIT
         WRITE (*, '(2X, 1PG8.1)',ADVANCE="NO") RESID_GRP(HYDRO_GRP)
         IF(GRANULAR_ENERGY) WRITE (*, '(2X, 1PG8.1)',ADVANCE="NO") RESID_GRP(THETA_GRP)
         IF(ENERGY_EQ) WRITE (*, '(2X, 1PG8.1)',ADVANCE="NO") RESID_GRP(ENERGY_GRP)
         IF(ANY_SPECIES_EQ) WRITE (*, '(2X, 1PG8.1)',ADVANCE="NO") RESID_GRP(SPECIES_GRP)
         IF(NScalar > 0) WRITE (*, '(2X, 1PG8.1)',ADVANCE="NO") RESID_GRP(SCALAR_GRP)
         IF(K_EPSILON) WRITE (*, '(2X, 1PG8.1)',ADVANCE="NO") RESID_GRP(KE_GRP)
         WRITE (*, '(2X, A)') RESID_STRING(8)


      ELSE

         IF (NIT == 1) THEN
            WRITE (*, '(A)',ADVANCE="NO") '  Nit'
            DO L = 1, 8
               IF (RESID_INDEX(L,1) /= UNDEFINED_I) WRITE (*, '(5X, A4)',ADVANCE="NO") &
                  RESID_STRING(L)
            END DO
            IF (RESID_INDEX(8,1) == UNDEFINED_I) THEN
               WRITE (*, '(2X, A)') 'Max res'
            ELSE
               WRITE (*, *)
            ENDIF
         ENDIF
!
         WRITE (*, '(I5)',ADVANCE="NO") NIT
         DO L = 1, 8
            IF (RESID_INDEX(L,1) /= UNDEFINED_I) WRITE (*, '(2X, 1PG8.1)',ADVANCE="NO") &
               RESID(RESID_INDEX(L,1),RESID_INDEX(L,2))
         END DO
         IF (RESID_INDEX(8,1) == UNDEFINED_I) THEN
            WRITE (*, '(2X, A)') RESID_STRING(8)
         ELSE
            WRITE (*, *)
         ENDIF

      ENDIF

!
!
!
!     Display maximum values of residuals
      IF(PRINT_IJK) WRITE(*,'(A, G12.3, 3I6, A, G12.3, 3I6, A, G12.3)') &
      & " Max Res/IJK: P_g: ", MAX_RESID(RESID_P, 0), &
      & I_OF_G(IJK_RESID(RESID_P, 0)), &
      & J_OF_G(IJK_RESID(RESID_P, 0)), &
      & K_OF_G(IJK_RESID(RESID_P, 0)), &
      & " P_s: ", MAX_RESID(RESID_p, 1), &
      & I_OF_G(IJK_RESID(RESID_p, 1)), &
      & J_OF_G(IJK_RESID(RESID_p, 1)), &
      & K_OF_G(IJK_RESID(RESID_p, 1)), &
      & " P_star=",  P_star(IJK_RESID(RESID_p, 1))
!
      RETURN

      contains

!       ----------------------------------------------------------
!       functions for determining i, j, k for global ijk_resid
!       -----------------------------------------------------------
!        integer i_of_g,j_of_g,k_of_g, ijk

        integer function k_of_g(ijk)
          integer ijk
          k_of_g = int( (ijk-1)/( (imax3-imin3+1)*(jmax3-jmin3+1) ) ) + kmin3
        end function

        integer function i_of_g(ijk)
          integer ijk
          i_of_g = int( ( (ijk-  (k_of_g(ijk)-kmin3)*((imax3-imin3+1)*(jmax3-jmin3+1))) &
                      - 1)/(jmax3-jmin3+1)) + imin3
        end function

        integer function j_of_g(ijk)
          integer ijk
          j_of_g = ijk - (i_of_g(ijk)-imin3)*(jmax3-jmin3+1) - &
               (k_of_g(ijk)-kmin3)*((imax3-imin3+1)*(jmax3-jmin3+1)) - 1 + jmin3
        end function


      END SUBROUTINE DISPLAY_RESID

!// Comments on the modifications for DMP version implementation
!// 001 Include header file and common declarations for parallelization
