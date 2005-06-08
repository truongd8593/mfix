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
      INTEGER          IJK
!
      DOUBLE PRECISION Y(NX)
!
!                      Jacobian matrix at t and t+dt
      DOUBLE PRECISION Ja(NX,NX), Ja1(NX,NX)
!                     
!                      Loop indices
      INTEGER          NL, NM
!
!                      Time for integration
      DOUBLE PRECISION ODE_Dt
!
!                      Controlling parameters for odepack
      DOUBLE PRECISION RWORK(22+9*NX+NX**2+10), ODSpec(NX), T1
      INTEGER          IWORK(20+NX+10), ITASK, ISTATE
      INTEGER          LRW, LIW, NEQ(2)

       external FEX, JAC
!
!
       NEED(1) = 1
       NEED(2) = 1
       NEED(3) = 0
!
       IJK = IUSR(1)
!
!      Allocate values
!
       DO NL = 1, NX
         Y(NL) = X(NL)
       END DO 
!
!     Jacobian matrix at t = 0
!
      CALL CALC_JACOBIAN(Y, NX, Ja, IJK)
!
!     Call ode solver to get the new values for reaction progress 
!
      DO NL = 1, NX
         ODSpec(NL) = X(NL)
      END DO
!
!     Contants for odepack
!
      NEQ(1) = NX
      NEQ(2) = IJK
!
!     Integration time
!
      ODE_dt = ISATdt
!
      T1 = ZERO
!
!     Controlling parameters for ODEPACK
      ISTATE = 1
      ITASK  = 1
      LRW    = 22+9*NX+NX**2+10
      LIW    = 20+NX+10
!
      CALL DLSODA(FEX, NEQ, ODSPEC, T1, ODE_dt, ITOL, RTOL, ATOL, &
         ITASK, ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, JT)
!         
!     Get fa
! 
      DO NL = 1, NX
         Fa(NL) = ODSpec(NL)
      END DO
!      
!     Get ga  
!     Jacobian at t = dt 
      DO NL = 1, NX
         Y(NL) = Fa(NL)
      END DO

      CALL CALC_JACOBIAN(Y, NX, Ja1, IJK)
!
!       Ja is equal to mid-point value
!
      DO NL = 1, NX
         DO NM = 1, NX
            Ja(NL, NM) = (Ja(NL, NM) + Ja1(NL, NM))*HALF
         END DO
      END DO
!
!     get G(dt) = G(0)exp(Ja*dt) G(0) is a unit matrix.
!     use the matrix exponential method to calculate
      CALL EXPONENTIAL(Ga, Ja, NX, ODE_dt)

      RETURN
      END

