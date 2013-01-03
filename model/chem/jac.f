!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                         C  
!     Module name: JAC                                                    C
!     Purpose: Provide the Jacobian matrix of source terms of ODEs        C
!              in ODEPACK solver PD(I, J) = df(i)/dy(j)                   C 
!                                                                         C
!     Author: Nan Xie                                   Date: 02-Aug-04   C
!     Reviewer:                                         Date:             C
!                                                                         C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE JAC(NEQ, T, Y, ML, MU, PD, NROWPD)
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      Use run
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
!                      NEQ(1)Number of equations, NEQ(2) IJK of the cell      
      INTEGER          NEQ(2)
      INTEGER          IJK
      INTEGER          NX
!
!                      half-bandwidth parameters
      INTEGER          ML, MU 
!
! 
      DOUBLE PRECISION T
!
!                      variable
      DOUBLE PRECISION Y(NEQ(1))
!
!                      Jacobian matrix
      INTEGER          NROWPD
      DOUBLE PRECISION PD(NROWPD, NEQ(1))
!
!                      Loop indice
      INTEGER          NL, NM

!begin by Q.X
      IF(SOLID_RO_V) RETURN
!end by Q.X
!
      NX  = NEQ(1)
      IJK = NEQ(2)
!
!     Get jacobian matrix and store in PD
!
      CALL CALC_JACOBIAN(Y, NX, PD, IJK)


      RETURN
      END

