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
!                      Print Location of Max_Resid
      LOGICAL,PARAMETER:: Print_ijk=.FALSE.
! 
!       ----------------------------------------------------------
!       inline functions for determining i, j, k for global ijk_resid
!       -----------------------------------------------------------
        integer i_of_g,j_of_g,k_of_g, ijk

        k_of_g(ijk) = int( (ijk-1)/( (imax3-imin3+1)*(jmax3-jmin3+1) ) ) + kmin3
        i_of_g(ijk) = int( ( (ijk-  (k_of_g(ijk)-kmin3)*((imax3-imin3+1)*(jmax3-jmin3+1))) &
                      - 1)/(jmax3-jmin3+1)) + imin3
        j_of_g(ijk) = ijk - (i_of_g(ijk)-imin3)*(jmax3-jmin3+1) - &
                      (k_of_g(ijk)-kmin3)*((imax3-imin3+1)*(jmax3-jmin3+1)) - 1 + jmin3

  
!
      if (myPE.ne.PE_IO) return   
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
      END SUBROUTINE DISPLAY_RESID 

!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization
