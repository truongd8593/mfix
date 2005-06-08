!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                         C  
!     Module name: EXPONENTIAL                                            C
!     Purpose: partial G/ partial t = Ja*G  G(0) = I (unit matrix)        C
!              G(dt) = G(0)exp(Ja*dt)    where G and Ja are matrix        C 
!                                                                         C
!     Author: Nan Xie                                   Date: 02-Aug-04   C
!     Reviewer:                                         Date:             C
!                                                                         C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE EXPONENTIAL(Ga, Ja, NX, ODE_dt)
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------   
      IMPLICIT NONE
!-----------------------------------------------
!     G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!     L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: ideg = 6
!-----------------------------------------------
!     L o c a l   V a r i a b l e s
!----------------------------------------------- 
!
!                      Number of odes
      INTEGER          NX
!
!                      Variables for dgpadm
      DOUBLE PRECISION ODE_dt
!
!                      Ga will return values for isatab
      DOUBLE PRECISION Ga(NX, NX)
!
!                      Ja is the mid-point reaction source Jacobian
      DOUBLE PRECISION Ja(NX, NX)




      RETURN
      END

