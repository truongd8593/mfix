!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                         C  
!     Module name: CALC_JACOBIAN                                          C
!     Purpose: Provide the Jacobian matrix of source terms of ODEs        C
!                                                                         C
!     Author: Nan Xie                                   Date: 02-Aug-04   C
!     Reviewer:                                         Date:             C
!                                                                         C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CALC_JACOBIAN(Y, NX, Ja, IJK)
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param1
      USE usr
      IMPLICIT NONE
!-----------------------------------------------
!     G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!     L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!     L o c a l   V a r i a b l e s
!----------------------------------------------- 
!
!                      Number of odes
      INTEGER          NX
!
      INTEGER          IJK
!
!                      Values of odes
      DOUBLE PRECISION Y(NX)
!
!                      Jacobian matrix 
      DOUBLE PRECISION Ja(NX,NX)
!  
!                      ADIFOR
      DOUBLE PRECISION G_X(NX,NX), G_Y(NX,NX),Y_source(NX)
!                     
!                      Loop indices
      INTEGER          NL, NM

!
!     Using ADIFOR to calulate the Jacobian 
!              
      DO NL = 1, NX
         DO NM = 1, NX
            G_X(NL, NM) = ZERO
         END DO
      END DO

      DO NL = 1, NX
         DO NM = 1, NX
            IF (NL .EQ. NM) THEN
               G_X(NL, NM) = ONE
            END IF
         END DO
      END DO   
      
      CALL G_DERIVS(NX, Y, G_X, NX, Y_source, G_Y, NX, N_sh(IJK,1))
!
!     transport g_y
!
      DO NL = 1, NX
         DO NM = 1, NX
            Ja(NM, NL) = G_Y(NL, NM)
         END DO
      END DO

      RETURN
      END 
