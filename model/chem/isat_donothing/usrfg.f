!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                         C  
!     Module name: USRFG                                                  C
!     Purpose:  This subroutine is provided to evaluate f(x), g(x)        C 
!                                                                         C
!     Author: Nan Xie                                   Date: 02-Aug-04   C
!     Reviewer:                                         Date:             C
!                                                                         C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE USRFG(NEED, NX, X, NF, NH, IUSR, RUSR, Fa, Ga, Ha) 
!
!-----------------------------------------------
!   M o d u l e s 
!----------------------------------------------- 
      USE param1
      USE run
      USE mchem
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
!     need       integer(3),  input: the values of need indicate whether or not
!                f(x) and g(x) need to be evaluated.
!                     need(1) = 1 f(x) needs to be evaluated
!                     need(2) = 1 g(x) needs to be evaluated
!                     need(3) = 1 h(x) needs to be evaluated
!                     need() = 0  not evaluated
      INTEGER          NEED(3)
!
!                      Dimension of odes
      INTEGER          NX, NF, NH
!
!                      User-defined integer array passed through isatab
      INTEGER          IUSR(1)
!
!                      User-defined real array passed through isatab
      INTEGER          RUSR
!
!                      Components of f(x) which must be returned if need(1)=1
      DOUBLE PRECISION Fa(NF)
!
!                      Components of g(x) which must be returned if need(2)=1 
!                      Ga(i,j) = partial fi/ partial xj
      DOUBLE PRECISION Ga(NF, NX)   
!
!                      Components of h(x) which must be returned if need(3)=1
      DOUBLE PRECISION Ha(NH)
!
!                      ODEs values at t
      DOUBLE PRECISION X(NX)
!
!
!


      RETURN
      END

