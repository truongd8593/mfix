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
!                      Phase index 
      INTEGER          M 
! 
  
!
!
      IF (NIT == 1) THEN 
         WRITE (*, '(A, $)') '  Nit' 
         DO L = 1, 8 
            IF (RESID_INDEX(L,1) /= UNDEFINED_I) WRITE (*, '(5X, A4, $)') &
               RESID_STRING(L) 
         END DO 
         IF (RESID_INDEX(8,1) == UNDEFINED_I) THEN 
            WRITE (*, '(2X, A)') 'Max res' 
         ELSE 
            WRITE (*, *) 
         ENDIF 
      ENDIF 
!
      WRITE (*, '(I5, $)') NIT 
      DO L = 1, 8 
         IF (RESID_INDEX(L,1) /= UNDEFINED_I) WRITE (*, '(2X, 1PG7.1, $)') &
            RESID(RESID_INDEX(L,1),RESID_INDEX(L,2)) 
      END DO 
      IF (RESID_INDEX(8,1) == UNDEFINED_I) THEN 
         WRITE (*, '(2X, A)') RESID_STRING(8) 
      ELSE 
         WRITE (*, *) 
      ENDIF 
!
!
!
!     Display maximum values of residuals
!      WRITE(*,'(A, G12.3, I6, A, G12.3, I6, A, G12.3)')
!     & " Max Res/IJK: P_g: ", MAX_RESID(RESID_P, 0),
!     & IJK_RESID(RESID_P, 0),
!     & " P_s: ", MAX_RESID(RESID_p, 1),
!     & IJK_RESID(RESID_p, 1),
!     & " P_star=",  P_star(IJK_RESID(RESID_p, 1))
!
      RETURN  
      END SUBROUTINE DISPLAY_RESID 
