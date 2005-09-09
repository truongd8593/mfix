!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                         C  
!     Module name: REACT                                                  C
!     Purpose: user CFD code (MFIX) to call ISATAB or DI subroutines      C 
!                                                                         C
!     Author: Nan Xie                                   Date: 02-Aug-04   C
!     Reviewer:                                         Date:             C
!                                                                         C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE REACT(Spec0, R_tmp1, IJK, NSpec1, ODE_dt)
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param1
      USE toleranc
      USE fldvar
      USE physprop
      USE rxns
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
!
!                      Dimensions of ODEs solved in ISAT or DI
      INTEGER          NSpec1
!
!                      ODEs solved in ISAT or DI   
      DOUBLE PRECISION Spec0(NSpec1)
! 
!                      Mass transfer 
      Double precision R_tmp1(0:MMAX,0:MMAX)
!
!
      INTEGER          IJK
!
!                      for ISATAB table
      INTEGER          NX, NF
      INTEGER          Nh, Nhd
      DOUBLE PRECISION Ha(1)
!
!
      DOUBLE PRECISION Stats(100)
!
!                      User-defined integer passed to usrfg
      INTEGER          IUSR(1)
!
!                      Usr-defined real array passed to usrfg
      DOUBLE PRECISION RUSR
!
!                      ODEs solved in ISAT or DI(after Integration)
      DOUBLE PRECISION Spect(NSpec1)
!
!                      ISATAB marching matrix
      DOUBLE PRECISION Ga(NSpec1, NSpec1)
!
!                      For calculating mass transfer 
!                      (1 indicate at t=0, 2 at t=dt)
      INTEGER          N1(2), N2(2)
      DOUBLE PRECISION X1, X2, YDOT1(NSpec1), YDOT2(NSpec1)
      DOUBLE PRECISION R1_g(DIM_N_g), R2_g(DIM_N_g)
      DOUBLE PRECISION R1_s(DIM_M, DIM_N_s), R2_s(DIM_M, DIM_N_s)    
!
!                      Molecular weight for mixture
      DOUBLE PRECISION MW_MIX_g_ISAT 
!
!                      Loop indices
      INTEGER          NL, NM  
!
!                      Time for integration
      DOUBLE PRECISION ODE_Dt
!
!                      Controlling parameters for odepack
      DOUBLE PRECISION RWORK(22+9*NSpec1+NSpec1**2+10), ODSpec(NSpec1), T1
      INTEGER          IWORK(20+NSpec1+10), ITASK, ISTATE
      INTEGER          LRW, LIW, NEQ(2)    
!
!
      external USRFG, FEX, JAC
!
!     Initial value for arrays
!
      R1_g = ZERO
      R2_g = ZERO
      R1_s = ZERO
      R2_s = ZERO
      SUM_R_G(IJK) = ZERO
      DO NL = 1, MMAX
         SUM_R_S(IJK, NL) = ZERO
      END DO
      DO NL = 1, MMAX
         R_tmp1(0, MMAX) = ZERO
      END DO
!
!     Dimensions of ISATAB
!
      NX = NSpec1
      NF = NSpec1
!
!     Initial values for ISATAB table 
!
      Nh = 0
      Nhd = 1
!
!     The integer variables transferred to usrfg
!
      IUSR(1) = IJK
!
!     Mass transfer at t = 0.0
!
      N1(1) = NSpec1
      N1(2) = IJK
      X1 = ZERO
      CALL FEX(N1, X1, SPEC0, YDOT1) 
      DO NL = 1, NMAX(0)
         R1_g(NL) = RXN_source_g(NL)
      END DO
      DO NL = 1, MMAX
         DO NM = 1, NMAX(NL)
            R1_s(NL, NM) = RXN_source_s(NL, NM)
         END DO
      END DO

      IF (CALL_DI) THEN
!
!     Call ode solver to get the new values for reaction progress 
!
         DO NL = 1, NSpec1
            ODSpec(NL) = SPEC0(NL)
         END DO
!
!     Contants for odepack
!
         NEQ(1) = NSpec1
         NEQ(2) = IJK
!
!
         T1 = ZERO
!
!     Controlling parameters for ODEPACK
         ISTATE = 1
         ITASK  = 1
         LRW    = 22+9*NSpec1+NSpec1**2+10
         LIW    = 20+NSpec1+10
!
         CALL DLSODA(FEX, NEQ, ODSPEC, T1, ODE_dt, ITOL, RTOL, ATOL, &
         ITASK, ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, JT)
!         
!     Get Spec0
! 
         DO NL = 1, NSpec1
            Spec0(NL) = ODSpec(NL)
         END DO         

      END IF
!
!
!      
!      IF (CALL_ISAT) THEN
!        
!     Call isatab subroutine
!
!         CALL isatab(idtab,mode,nx,spec0,nf,nh,nhd,usrfg,iusr,rusr,info, &
!     &        rinfo,spect,ga,ha,stats)       
!
!     Transfer data back to spec0
!
!         DO NL = 1, NX
!            Spec0(NL) = Spect(NL)
!         END DO
!      END IF
!
!     Mass transfer at t = dt
!
      N2(1) = NSpec1
      N2(2) = IJK
      X2 = ZERO
      CALL FEX(N2, X2, SPEC0, YDOT2) 
      DO NL = 1, NMAX(0)
         R2_g(NL) = RXN_source_g(NL)
      END DO
      DO NL = 1, MMAX
         DO NM = 1, NMAX(NL)
            R2_s(NL, NM) = RXN_source_s(NL, NM)
         END DO
      END DO
!
!     Calculate P_G(IJK)
!
      MW_MIX_g_ISAT = ZERO
      DO NL = 1, NMAX(0)
         MW_MIX_g_ISAT = MW_MIX_g_ISAT + Spec0(NL+2)/MW_g(NL)
      END DO
      MW_MIX_g_ISAT = ONE/MAX(MW_MIX_g_ISAT,OMW_MAX)
      P_G(IJK) = Spec0(1)*(8314.56D4*Spec0(2))/MW_MIX_g_ISAT
!
!     Total change of gas phase in continuity eqn
! 
      SUM_R_G(IJK) = HALF * (SUM(R1_g(1:NMAX(0))) + SUM(R2_g(1:NMAX(0))))
!
!     Total change of solid phase in continuity eqn
!
      DO NL = 1, MMAX
         SUM_R_S(IJK, NL) = HALF * (SUM(R1_s(NL, 1:NMAX(NL))) + SUM(R2_s(NL, 1:NMAX(NL))))
      END DO
!
!      Store R_tmp values in an array
!
      DO NL = 1, MMAX
         R_tmp1(0, NL) = -HALF * (SUM(R1_s(NL, 1:NMAX(NL))) + SUM(R2_s(NL, 1:NMAX(NL))))
      END DO

      
      RETURN
      END

