!DECK DGMRES
      SUBROUTINE DGMRES(N, B, X, NELT, IA, JA, A, ISYM, MATVEC, MSOLVE, ITOL, &
         TOL, ITMAX, ITER, ERR, IERR, IUNIT, SB, SX, RGWK, LRGW, IGWK, LIGW, &
         RWORK, IWORK) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  15:55:30  12/05/98  
!...Switches: -xf
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER N, NELT, ISYM, ITOL, ITMAX, ITER, IERR, IUNIT, LRGW, LIGW 
      DOUBLE PRECISION TOL, ERR 
      INTEGER, DIMENSION(NELT) :: IA, JA 
      INTEGER, DIMENSION(LIGW) :: IGWK 
      INTEGER, DIMENSION(*) :: IWORK 
      DOUBLE PRECISION, DIMENSION(N) :: B, X 
      DOUBLE PRECISION, DIMENSION(NELT) :: A 
      DOUBLE PRECISION, DIMENSION(N) :: SB, SX 
      DOUBLE PRECISION, DIMENSION(LRGW) :: RGWK 
      DOUBLE PRECISION, DIMENSION(*) :: RWORK 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: JPRE, KMP, MAXL, NMS, MAXLP1, NMSL, NRSTS, NRMAX, I, IFLAG, LR&
         , LDL, LHES, LGMR, LQ, LV, LW, JSCAL, LXL, LZ, LZM1 
      DOUBLE PRECISION :: BNRM, RHOL, SUM 
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      DOUBLE PRECISION , EXTERNAL :: D1MACH, DNRM2 
      EXTERNAL MATVEC, MSOLVE 
!-----------------------------------------------
!
!***FIRST EXECUTABLE STATEMENT  DGMRES
      IERR = 0 
!   ------------------------------------------------------------------
!         Load method parameters with user values or defaults.
!   ------------------------------------------------------------------
      MAXL = IGWK(1) 
      IF (MAXL == 0) MAXL = 10 
      MAXL = MIN0(N,MAXL) 
      KMP = IGWK(2) 
      IF (KMP == 0) KMP = MAXL 
      KMP = MIN0(MAXL,KMP) 
      JSCAL = IGWK(3) 
      JPRE = IGWK(4) 
!         Check for consistent values of ITOL and JPRE.
      IF (ITOL/=1 .OR. JPRE>=0) THEN 
         IF (ITOL/=2 .OR. JPRE<0) THEN 
            NRMAX = IGWK(5) 
            IF (NRMAX == 0) NRMAX = 10 
!         If NRMAX .eq. -1, then set NRMAX = 0 to turn off restarting.
            IF (NRMAX == (-1)) NRMAX = 0 
!         If input value of TOL is zero, set it to its default value.
            IF (TOL == 0.0D0) TOL = 500.0*D1MACH(3) 
!
!         Initialize counters.
            ITER = 0 
            NMS = 0 
            NRSTS = 0 
!   ------------------------------------------------------------------
!         Form work array segment pointers.
!   ------------------------------------------------------------------
            MAXLP1 = MAXL + 1 
            LV = 1 
            LR = LV + N*MAXLP1 
            LHES = LR + N + 1 
            LQ = LHES + MAXL*MAXLP1 
            LDL = LQ + 2*MAXL 
            LW = LDL + N 
            LXL = LW + N 
            LZ = LXL + N 
!
!         Load igwk(6) with required minimum length of the rgwk array.
            IGWK(6) = LZ + N - 1 
            IF (LZ + N - 1 <= LRGW) THEN 
!   ------------------------------------------------------------------
!         Calculate scaled-preconditioned norm of RHS vector b.
!   ------------------------------------------------------------------
               IF (JPRE < 0) THEN 
                  CALL MSOLVE(N,B,RGWK(LR),NELT,IA,JA,A,ISYM,RWORK,IWORK) 
                  NMS = NMS + 1 
               ELSE 
                  CALL DCOPY (N, B, 1, RGWK(LR), 1) 
               ENDIF 
               IF (JSCAL==2 .OR. JSCAL==3) THEN 
                  SUM = 0.D0 
                  DO I = 1, N 
                     SUM = SUM + (RGWK(LR-1+I)*SB(I))**2 
                  END DO 
                  BNRM = DSQRT(SUM) 
               ELSE 
                  BNRM = DNRM2(N,RGWK(LR),1) 
               ENDIF 
!   ------------------------------------------------------------------
!         Calculate initial residual.
!   ------------------------------------------------------------------
               CALL MATVEC (N, X, RGWK(LR), NELT, IA, JA, A, ISYM) 
               I = 1 
               IF (N > 0) THEN 
                  RGWK(LR:N-1+LR) = B(:N) - RGWK(LR:N-1+LR) 
                  I = N + 1 
               ENDIF 
  100          CONTINUE 
               IF (NRSTS <= NRMAX) THEN 
!         Copy the curr residual to different loc in the Rgwk array.
                  IF (NRSTS > 0) CALL DCOPY (N, RGWK(LDL), 1, RGWK(LR), 1) 
!   ------------------------------------------------------------------
!         Use the DPIGMR algorithm to solve the linear system A*Z = R.
!   ------------------------------------------------------------------
                  CALL DPIGMR (N, RGWK(LR), SB, SX, JSCAL, MAXL, MAXLP1, KMP, &
                     NRSTS, JPRE, MATVEC, MSOLVE, NMSL, RGWK(LZ), RGWK(LV), &
                     RGWK(LHES), RGWK(LQ), LGMR, RWORK, IWORK, RGWK(LW), RGWK(&
                     LDL), RHOL, NRMAX, B, BNRM, X, RGWK(LXL), ITOL, TOL, NELT&
                     , IA, JA, A, ISYM, IUNIT, IFLAG, ERR) 
                  ITER = ITER + LGMR 
                  NMS = NMS + NMSL 
!
!         Increment X by the current approximate solution Z of A*Z = R.
!
                  LZM1 = LZ - 1 
                  I = 1 
                  IF (N > 0) THEN 
                     X(:N) = X(:N) + RGWK(LZM1+1:N+LZM1) 
                     I = N + 1 
                  ENDIF 
                  IF (IFLAG == 0) GO TO 600 
                  IF (IFLAG == 1) THEN 
                     NRSTS = NRSTS + 1 
                     GO TO 100 
                  ENDIF 
                  IF (IFLAG == 2) GO TO 620 
!   ------------------------------------------------------------------
!         All returns are made through this section.
!   ------------------------------------------------------------------
!         The iteration has converged.
!
  600             CONTINUE 
                  IGWK(7) = NMS 
                  RGWK(1) = RHOL 
                  IERR = 0 
                  RETURN  
!
!         Max number((NRMAX+1)*MAXL) of linear iterations performed.
               ENDIF 
               IGWK(7) = NMS 
               RGWK(1) = RHOL 
               IERR = 1 
               RETURN  
!
!         GMRES failed to reduce last residual in MAXL iterations.
!         The iteration has stalled.
  620          CONTINUE 
               IGWK(7) = NMS 
               RGWK(1) = RHOL 
               IERR = 2 
               RETURN  
!         Error return.  Insufficient length for Rgwk array.
            ENDIF 
            ERR = TOL 
            IERR = -1 
            RETURN  
!         Error return.  Inconsistent ITOL and JPRE values.
         ENDIF 
      ENDIF 
      ERR = TOL 
      IERR = -2 
      RETURN  
!------------- LAST LINE OF DGMRES FOLLOWS ----------------------------
      END SUBROUTINE DGMRES 
!DECK DSDGMR
      SUBROUTINE DSDGMR(N, B, X, NELT, IA, JA, A, ISYM, NSAVE, ITOL, TOL, ITMAX&
         , ITER, ERR, IERR, IUNIT, RWORK, LENW, IWORK, LENIW) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  15:55:30  12/05/98  
!...Switches: -xf
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER N, NELT, ISYM, NSAVE, ITOL, ITMAX, ITER, IERR, IUNIT, LENW, LENIW 
      DOUBLE PRECISION TOL, ERR 
      INTEGER, DIMENSION(NELT) :: IA, JA 
      INTEGER, DIMENSION(LENIW) :: IWORK 
      DOUBLE PRECISION, DIMENSION(N) :: B, X 
      DOUBLE PRECISION, DIMENSION(NELT) :: A 
      DOUBLE PRECISION, DIMENSION(LENW) :: RWORK 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: LOCRB = 1 
      INTEGER, PARAMETER :: LOCIB = 11 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: LOCIGW, LOCIW, LOCDIN, LOCRGW, LOCW, MYITOL 
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      DOUBLE PRECISION , EXTERNAL :: DSMV, DSDI 
!-----------------------------------------------
      IERR = 0 
      ERR = 0.0 
      IF (NSAVE <= 1) THEN 
         IERR = 3 
         RETURN  
      ENDIF 
      CALL DS2Y (N, NELT, IA, JA, A, ISYM) 
!
!         Set up the workspace.  We assume MAXL=KMP=NSAVE.
!         Compute the inverse of the diagonal of the matrix.
      LOCIGW = LOCIB 
      LOCIW = LOCIGW + 20 
!
      LOCDIN = LOCRB 
      LOCRGW = LOCDIN + N 
      LOCW = LOCRGW + 1 + N*(NSAVE + 6) + NSAVE*(NSAVE + 3) 
!
      IWORK(4) = LOCDIN 
      IWORK(9) = LOCIW 
      IWORK(10) = LOCW 
!
!         Check the workspace allocations.
      CALL DCHKW ('DSDGMR', LOCIW, LENIW, LOCW, LENW, IERR, ITER, ERR) 
      IF (IERR /= 0) RETURN  
!
      CALL DSDS (N, NELT, IA, JA, A, ISYM, RWORK(LOCDIN)) 
!
!         Perform the Diagonaly Scaled Generalized Minimum
!         Residual iteration algorithm.  The following DGMRES
!         defaults are used MAXL = KMP = NSAVE, JSCAL = 0,
!         JPRE = -1, NRMAX = ITMAX/NSAVE
      IWORK(LOCIGW) = NSAVE 
      IWORK(LOCIGW+1) = NSAVE 
      IWORK(LOCIGW+2) = 0 
      IWORK(LOCIGW+3) = -1 
      IWORK(LOCIGW+4) = ITMAX/NSAVE 
      MYITOL = 0 
!
      CALL DGMRES (N, B, X, NELT, IA, JA, A, ISYM, DSMV, DSDI, MYITOL, TOL, &
         ITMAX, ITER, ERR, IERR, IUNIT, RWORK, RWORK, RWORK(LOCRGW), LENW - &
         LOCRGW, IWORK(LOCIGW), 20, RWORK, IWORK) 
!
      IF (ITER > ITMAX) IERR = 2 
      RETURN  
!------------- LAST LINE OF DSDGMR FOLLOWS ----------------------------
      END SUBROUTINE DSDGMR 
!DECK DSLUGM
      SUBROUTINE DSLUGM(N, B, X, NELT, IA, JA, A, ISYM, NSAVE, ITOL, TOL, ITMAX&
         , ITER, ERR, IERR, IUNIT, RWORK, LENW, IWORK, LENIW) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  15:55:30  12/05/98  
!...Switches: -xf
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER N, NELT, ISYM, NSAVE, ITOL, ITMAX, ITER, IERR, IUNIT, LENW, LENIW 
      DOUBLE PRECISION TOL, ERR 
      INTEGER, DIMENSION(NELT) :: IA, JA 
      INTEGER, DIMENSION(LENIW) :: IWORK 
      DOUBLE PRECISION, DIMENSION(N) :: B, X 
      DOUBLE PRECISION, DIMENSION(NELT) :: A 
      DOUBLE PRECISION, DIMENSION(LENW) :: RWORK 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: LOCRB = 1 
      INTEGER, PARAMETER :: LOCIB = 11 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NL, NU, ICOL, JBGN, JEND, J, LOCIGW, LOCIL, LOCJL, LOCIU, &
         LOCJU, LOCNR, LOCNC, LOCIW, LOCL, LOCDIN, LOCU, LOCRGW, LOCW, MYITOL 
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      DOUBLE PRECISION , EXTERNAL :: DSMV, DSLUI 
!-----------------------------------------------
      IERR = 0 
      ERR = 0.0 
      IF (NSAVE <= 1) THEN 
         IERR = 3 
         RETURN  
      ENDIF 
      CALL DS2Y (N, NELT, IA, JA, A, ISYM) 
!
!         Count number of Non-Zero elements preconditioner ILU matrix.
!         Then set up the work arrays.  We assume MAXL=KMP=NSAVE.
      NL = 0 
      NU = 0 
      DO ICOL = 1, N 
!         Don't count diagonal.
         JBGN = JA(ICOL) + 1 
         JEND = JA(ICOL+1) - 1 
         IF (JBGN <= JEND) THEN 
!VD$ NOVECTOR
            DO J = JBGN, JEND 
               IF (IA(J) > ICOL) THEN 
                  NL = NL + 1 
                  IF (ISYM /= 0) NU = NU + 1 
               ELSE 
                  NU = NU + 1 
               ENDIF 
            END DO 
         ENDIF 
      END DO 
      LOCIGW = LOCIB 
      LOCIL = LOCIGW + 20 
      LOCJL = LOCIL + N + 1 
      LOCIU = LOCJL + NL 
      LOCJU = LOCIU + NU 
      LOCNR = LOCJU + N + 1 
      LOCNC = LOCNR + N 
      LOCIW = LOCNC + N 
!
      LOCL = LOCRB 
      LOCDIN = LOCL + NL 
      LOCU = LOCDIN + N 
      LOCRGW = LOCU + NU 
      LOCW = LOCRGW + 1 + N*(NSAVE + 6) + NSAVE*(NSAVE + 3) 
!
!         Check the workspace allocations.
      CALL DCHKW ('DSLUGM', LOCIW, LENIW, LOCW, LENW, IERR, ITER, ERR) 
      IF (IERR /= 0) RETURN  
!
      IWORK(1) = LOCIL 
      IWORK(2) = LOCJL 
      IWORK(3) = LOCIU 
      IWORK(4) = LOCJU 
      IWORK(5) = LOCL 
      IWORK(6) = LOCDIN 
      IWORK(7) = LOCU 
      IWORK(9) = LOCIW 
      IWORK(10) = LOCW 
!
!         Compute the Incomplete LU decomposition.
      CALL DSILUS (N, NELT, IA, JA, A, ISYM, NL, IWORK(LOCIL), IWORK(LOCJL), &
         RWORK(LOCL), RWORK(LOCDIN), NU, IWORK(LOCIU), IWORK(LOCJU), RWORK(LOCU&
         ), IWORK(LOCNR), IWORK(LOCNC)) 
!
!         Perform the Incomplet LU Preconditioned Generalized Minimum
!         Residual iteration algorithm.  The following DGMRES
!         defaults are used MAXL = KMP = NSAVE, JSCAL = 0,
!         JPRE = -1, NRMAX = ITMAX/NSAVE
      IWORK(LOCIGW) = NSAVE 
      IWORK(LOCIGW+1) = NSAVE 
      IWORK(LOCIGW+2) = 0 
      IWORK(LOCIGW+3) = -1 
      IWORK(LOCIGW+4) = ITMAX/NSAVE 
      MYITOL = 0 
!
      CALL DGMRES (N, B, X, NELT, IA, JA, A, ISYM, DSMV, DSLUI, MYITOL, TOL, &
         ITMAX, ITER, ERR, IERR, IUNIT, RWORK, RWORK, RWORK(LOCRGW), LENW - &
         LOCRGW, IWORK(LOCIGW), 20, RWORK, IWORK) 
!
      IF (ITER > ITMAX) IERR = 2 
      RETURN  
!------------- LAST LINE OF DSLUGM FOLLOWS ----------------------------
      END SUBROUTINE DSLUGM 
!DECK DHELS
      SUBROUTINE DHELS(A, LDA, N, Q, B) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  15:55:30  12/05/98  
!...Switches: -xf
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER LDA, N 
      DOUBLE PRECISION, DIMENSION(LDA,*) :: A 
      DOUBLE PRECISION, DIMENSION(*) :: Q, B 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IQ, K, KB, KP1 
      DOUBLE PRECISION :: C, S, T, T1, T2 
!-----------------------------------------------
!
!         Local Variables.
!
!
!         minimize(B-A*X,B-A*X).  First form Q*B.
!
      DO K = 1, N 
         KP1 = K + 1 
         IQ = 2*(K - 1) + 1 
         C = Q(IQ) 
         S = Q(IQ+1) 
         T1 = B(K) 
         T2 = B(KP1) 
         B(K) = C*T1 - S*T2 
         B(KP1) = S*T1 + C*T2 
      END DO 
      DO KB = 1, N 
         K = N + 1 - KB 
         B(K) = B(K)/A(K,K) 
         T = -B(K) 
         CALL DAXPY (K - 1, T, A(1,K), 1, B(1), 1) 
      END DO 
      RETURN  
!------------- LAST LINE OF DHELS FOLLOWS ----------------------------
      END SUBROUTINE DHELS 
!DECK DHEQR
      SUBROUTINE DHEQR(A, LDA, N, Q, INFO, IJOB) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  15:55:30  12/05/98  
!...Switches: -xf
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER LDA, N, INFO, IJOB 
      DOUBLE PRECISION, DIMENSION(LDA,*) :: A 
      DOUBLE PRECISION, DIMENSION(*) :: Q 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, IQ, J, K, KM1, KP1, NM1 
      DOUBLE PRECISION :: C, S, T, T1, T2 
!-----------------------------------------------
!
!         Local Variables.
!
!
!***FIRST EXECUTABLE STATEMENT  DHEQR
      IF (IJOB <= 1) THEN 
!   -------------------------------------------------------------------
!         A new facorization is desired.
!   -------------------------------------------------------------------
!         QR decomposition without pivoting.
!
         INFO = 0 
         DO K = 1, N 
            KM1 = K - 1 
            KP1 = K + 1 
!
!           Compute K-th column of R.
!           First, multiply the K-th column of a by the previous
!           K-1 Givens rotations.
!
            IF (KM1 >= 1) THEN 
               DO J = 1, KM1 
                  I = 2*(J - 1) + 1 
                  T1 = A(J,K) 
                  T2 = A(J+1,K) 
                  C = Q(I) 
                  S = Q(I+1) 
                  A(J,K) = C*T1 - S*T2 
                  A(J+1,K) = S*T1 + C*T2 
               END DO 
            ENDIF 
            IQ = 2*KM1 + 1 
            T1 = A(K,K) 
            T2 = A(KP1,K) 
            IF (T2 == 0.0D0) THEN 
               C = 1.0D0 
               S = 0.0D0 
            ELSE IF (ABS(T2) >= ABS(T1)) THEN 
               T = T1/T2 
               S = -1.0D0/DSQRT(1.0D0 + T*T) 
               C = -S*T 
            ELSE 
               T = T2/T1 
               C = 1.0D0/DSQRT(1.0D0 + T*T) 
               S = -C*T 
            ENDIF 
            Q(IQ) = C 
            Q(IQ+1) = S 
            A(K,K) = C*T1 - S*T2 
            IF (A(K,K) == 0.0D0) INFO = K 
         END DO 
         RETURN  
!   -------------------------------------------------------------------
!         The old factorization of a will be updated.  A row and a
!         column has been added to the matrix A.  N by N-1 is now
!         the old size of the matrix.
!   -------------------------------------------------------------------
      ENDIF 
      NM1 = N - 1 
!   -------------------------------------------------------------------
!         Multiply the new column by the N previous Givens rotations.
!   -------------------------------------------------------------------
      DO K = 1, NM1 
         I = 2*(K - 1) + 1 
         T1 = A(K,N) 
         T2 = A(K+1,N) 
         C = Q(I) 
         S = Q(I+1) 
         A(K,N) = C*T1 - S*T2 
         A(K+1,N) = S*T1 + C*T2 
      END DO 
      INFO = 0 
      T1 = A(N,N) 
      T2 = A(N+1,N) 
      IF (T2 == 0.0D0) THEN 
         C = 1.0D0 
         S = 0.0D0 
      ELSE IF (ABS(T2) >= ABS(T1)) THEN 
         T = T1/T2 
         S = -1.0D0/DSQRT(1.0D0 + T*T) 
         C = -S*T 
      ELSE 
         T = T2/T1 
         C = 1.0D0/DSQRT(1.0D0 + T*T) 
         S = -C*T 
      ENDIF 
      IQ = 2*N - 1 
      Q(IQ) = C 
      Q(IQ+1) = S 
      A(N,N) = C*T1 - S*T2 
      IF (A(N,N) == 0.0D0) INFO = N 
      RETURN  
!------------- LAST LINE OF DHEQR FOLLOWS ----------------------------
      END SUBROUTINE DHEQR 
!DECK DORTH
      SUBROUTINE DORTH(VNEW, V, HES, N, LL, LDHES, KMP, SNORMW) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  15:55:30  12/05/98  
!...Switches: -xf
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER N, LL, LDHES, KMP 
      DOUBLE PRECISION SNORMW 
      DOUBLE PRECISION, DIMENSION(*) :: VNEW 
      DOUBLE PRECISION, DIMENSION(N,*) :: V 
      DOUBLE PRECISION, DIMENSION(LDHES,*) :: HES 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, I0 
      DOUBLE PRECISION :: ARG, SUMDSQ, TEM, VNRM 
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      DOUBLE PRECISION , EXTERNAL :: DNRM2, DDOT 
!-----------------------------------------------
!
!         Internal variables.
!
!
!         Get norm of unaltered VNEW for later use.
!***FIRST EXECUTABLE STATEMENT  DORTH
      VNRM = DNRM2(N,VNEW,1) 
!   -------------------------------------------------------------------
!         Perform the modified gram-schmidt procedure on VNEW =A*V(LL).
!         Scaled inner products give new column of HES.
!         Projections of earlier vectors are subtracted from VNEW.
!   -------------------------------------------------------------------
      I0 = MAX0(1,LL - KMP + 1) 
      DO I = I0, LL 
         HES(I,LL) = DDOT(N,V(1,I),1,VNEW,1) 
         TEM = -HES(I,LL) 
         CALL DAXPY (N, TEM, V(1,I), 1, VNEW, 1) 
      END DO 
      SNORMW = DNRM2(N,VNEW,1) 
      IF (VNRM + 0.001D0*SNORMW /= VNRM) RETURN  
      SUMDSQ = 0.0D0 
      DO I = I0, LL 
         TEM = -DDOT(N,V(1,I),1,VNEW,1) 
         IF (HES(I,LL) + 0.001D0*TEM /= HES(I,LL)) THEN 
            HES(I,LL) = HES(I,LL) - TEM 
            CALL DAXPY (N, TEM, V(1,I), 1, VNEW, 1) 
            SUMDSQ = SUMDSQ + TEM**2 
         ENDIF 
      END DO 
      IF (SUMDSQ == 0.0D0) RETURN  
      ARG = MAX(0.0D0,SNORMW**2 - SUMDSQ) 
      SNORMW = DSQRT(ARG) 
!
      RETURN  
!------------- LAST LINE OF DORTH FOLLOWS ----------------------------
      END SUBROUTINE DORTH 
!DECK DPIGMR
      SUBROUTINE DPIGMR(N, R0, SR, SZ, JSCAL, MAXL, MAXLP1, KMP, NRSTS, JPRE, &
         MATVEC, MSOLVE, NMSL, Z, V, HES, Q, LGMR, RPAR, IPAR, WK, DL, RHOL, &
         NRMAX, B, BNRM, X, XL, ITOL, TOL, NELT, IA, JA, A, ISYM, IUNIT, IFLAG&
         , ERR) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  15:55:30  12/05/98  
!...Switches: -xf
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER N, JSCAL, MAXL, MAXLP1, KMP, NRSTS, JPRE, NMSL, LGMR, NRMAX, ITOL&
         , NELT, ISYM, IUNIT, IFLAG 
      DOUBLE PRECISION RHOL, BNRM, TOL, ERR 
      INTEGER, DIMENSION(*) :: IPAR 
      INTEGER, DIMENSION(NELT) :: IA, JA 
      DOUBLE PRECISION, DIMENSION(*) :: R0, SR, SZ, Z 
      DOUBLE PRECISION, DIMENSION(N,*) :: V 
      DOUBLE PRECISION, DIMENSION(MAXLP1,*) :: HES 
      DOUBLE PRECISION, DIMENSION(*) :: Q, RPAR, WK, DL, B, X, XL 
      DOUBLE PRECISION, DIMENSION(NELT) :: A 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, INFO, IP1, I2, J, K, LL, LLP1, ITMAX, ITER 
      DOUBLE PRECISION :: R0NRM, C, DLNRM, PROD, RHO, S, SNORMW, TEM 
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      INTEGER , EXTERNAL :: ISDGMR 
      DOUBLE PRECISION , EXTERNAL :: DNRM2 
      EXTERNAL MATVEC, MSOLVE 
!-----------------------------------------------
!
!         Local variables.
!
!
!         Zero out the z array.
!***FIRST EXECUTABLE STATEMENT  DPIGMR
      I = 1 
      IF (N > 0) THEN 
         Z(:N) = 0.0D0 
         I = N + 1 
      ENDIF 
      IFLAG = 0 
      LGMR = 0 
      NMSL = 0 
!         Load ITMAX, the maximum number of iterations.
      ITMAX = (NRMAX + 1)*MAXL 
!   -------------------------------------------------------------------
!         The initial residual is the vector R0.
!         Apply left precon. if JPRE < 0 and this is not a restart.
!         Apply scaling to R0 if JSCAL = 2 or 3.
!   -------------------------------------------------------------------
      IF (JPRE<0 .AND. NRSTS==0) THEN 
         CALL DCOPY (N, R0, 1, WK, 1) 
         CALL MSOLVE (N, WK, R0, NELT, IA, JA, A, ISYM, RPAR, IPAR) 
         NMSL = NMSL + 1 
      ENDIF 
      IF ((JSCAL==2 .OR. JSCAL==3) .AND. NRSTS==0) THEN 
         I = 1 
         IF (N > 0) THEN 
            V(:N,1) = R0(:N)*SR(:N) 
            I = N + 1 
         ENDIF 
      ELSE 
         I = 1 
         IF (N > 0) THEN 
            V(:N,1) = R0(:N) 
            I = N + 1 
         ENDIF 
      ENDIF 
      R0NRM = DNRM2(N,V,1) 
      ITER = NRSTS*MAXL 
!
!         Call stopping routine ISDGMR.
!
      IF (ISDGMR(N,B,X,XL,NELT,IA,JA,A,ISYM,MSOLVE,NMSL,ITOL,TOL,ITMAX,ITER,ERR&
         ,IUNIT,V(1,1),Z,WK,RPAR,IPAR,R0NRM,BNRM,SR,SZ,JSCAL,KMP,LGMR,MAXL,&
         MAXLP1,V,Q,SNORMW,PROD,R0NRM,HES,JPRE) /= 0) RETURN  
      TEM = 1.0D0/R0NRM 
      CALL DSCAL (N, TEM, V(1,1), 1) 
!
!         Zero out the HES array.
!
      J = 1 
      IF (MAXL > 0) THEN 
         I = 1 
         IF (MAXLP1 > 0) THEN 
            HES(:MAXLP1,:MAXL) = 0.0D0 
            I = MAXLP1 + 1 
         ENDIF 
         J = MAXL + 1 
      ENDIF 
      PROD = 1.0D0 
      DO LL = 1, MAXL 
         LGMR = LL 
!   -------------------------------------------------------------------
!        Unscale  the  current V(LL)  and store  in WK.  Call routine
!        msolve    to   compute(M-inverse)*WK,   where    M   is  the
!        preconditioner matrix.  Save the answer in Z.   Call routine
!        MATVEC to compute  VNEW  = A*Z,  where  A is  the the system
!        matrix.  save the answer in  V(LL+1).  Scale V(LL+1).   Call
!        routine DORTH  to  orthogonalize the    new vector VNEW   =
!        V(*,LL+1).  Call routine DHEQR to update the factors of HES.
!   -------------------------------------------------------------------
         IF (JSCAL==1 .OR. JSCAL==3) THEN 
            I = 1 
            IF (N > 0) THEN 
               WK(:N) = V(:N,LL)/SZ(:N) 
               I = N + 1 
            ENDIF 
         ELSE 
            CALL DCOPY (N, V(1,LL), 1, WK, 1) 
         ENDIF 
         IF (JPRE > 0) THEN 
            CALL MSOLVE (N, WK, Z, NELT, IA, JA, A, ISYM, RPAR, IPAR) 
            NMSL = NMSL + 1 
            CALL MATVEC (N, Z, V(1,LL+1), NELT, IA, JA, A, ISYM) 
         ELSE 
            CALL MATVEC (N, WK, V(1,LL+1), NELT, IA, JA, A, ISYM) 
         ENDIF 
         IF (JPRE < 0) THEN 
            CALL DCOPY (N, V(1,LL+1), 1, WK, 1) 
            CALL MSOLVE (N, WK, V(1,LL+1), NELT, IA, JA, A, ISYM, RPAR, IPAR) 
            NMSL = NMSL + 1 
         ENDIF 
         IF (JSCAL==2 .OR. JSCAL==3) THEN 
            I = 1 
            IF (N > 0) THEN 
               V(:N,LL+1) = V(:N,LL+1)*SR(:N) 
               I = N + 1 
            ENDIF 
         ENDIF 
         CALL DORTH (V(1,LL+1), V, HES, N, LL, MAXLP1, KMP, SNORMW) 
         HES(LL+1,LL) = SNORMW 
         CALL DHEQR (HES, MAXLP1, LL, Q, INFO, LL) 
         IF (INFO == LL) GO TO 120 
!   -------------------------------------------------------------------
!         Update RHO, the estimate of the norm of the residual R0-A*ZL.
!         If KMP <  MAXL, then the vectors V(*,1),...,V(*,LL+1) are not
!         necessarily orthogonal for LL > KMP.  The vector DL must then
!         be computed, and its norm used in the calculation of RHO.
!   -------------------------------------------------------------------
         PROD = PROD*Q(2*LL) 
         RHO = ABS(PROD*R0NRM) 
         IF (LL>KMP .AND. KMP<MAXL) THEN 
            IF (LL == KMP + 1) THEN 
               CALL DCOPY (N, V(1,1), 1, DL, 1) 
               DO I = 1, KMP 
                  IP1 = I + 1 
                  I2 = I*2 
                  S = Q(I2) 
                  C = Q(I2-1) 
                  K = 1 
                  IF (N > 0) THEN 
                     DL(:N) = S*DL(:N) + C*V(:N,IP1) 
                     K = N + 1 
                  ENDIF 
               END DO 
            ENDIF 
            S = Q(2*LL) 
            C = Q(2*LL-1)/SNORMW 
            LLP1 = LL + 1 
            K = 1 
            IF (N > 0) THEN 
               DL(:N) = S*DL(:N) + C*V(:N,LLP1) 
               K = N + 1 
            ENDIF 
            DLNRM = DNRM2(N,DL,1) 
            RHO = RHO*DLNRM 
         ENDIF 
         RHOL = RHO 
!   -------------------------------------------------------------------
!         Test for convergence.  If passed, compute approximation ZL.
!         If failed and LL < MAXL, then continue iterating.
!   -------------------------------------------------------------------
         ITER = NRSTS*MAXL + LGMR 
         IF (ISDGMR(N,B,X,XL,NELT,IA,JA,A,ISYM,MSOLVE,NMSL,ITOL,TOL,ITMAX,ITER,&
            ERR,IUNIT,DL,Z,WK,RPAR,IPAR,RHOL,BNRM,SR,SZ,JSCAL,KMP,LGMR,MAXL,&
            MAXLP1,V,Q,SNORMW,PROD,R0NRM,HES,JPRE) /= 0) GO TO 200 
         IF (LL == MAXL) EXIT  
!   -------------------------------------------------------------------
!         Rescale so that the norm of V(1,LL+1) is one.
!   -------------------------------------------------------------------
         TEM = 1.0D0/SNORMW 
         CALL DSCAL (N, TEM, V(1,LL+1), 1) 
      END DO 
      IF (RHO < R0NRM) GO TO 150 
  120 CONTINUE 
      IFLAG = 2 
!
!         Load approximate solution with zero.
!
      I = 1 
      IF (N > 0) THEN 
         Z(:N) = 0.D0 
         I = N + 1 
      ENDIF 
      RETURN  
  150 CONTINUE 
      IFLAG = 1 
!
!        If performing restarting (NRMAX > 0)  calculate the residual
!        vector RL and  store it in the DL  array.  If the incomplete
!        version is being used (KMP < MAXL) then DL has  already been
!        calculated up to a scaling factor.   Use DRLCAL to calculate
!        the scaled residual vector.
!
      IF(NRMAX>0)CALL DRLCAL(N,KMP,MAXL,MAXL,V,Q,DL,SNORMW,PROD,R0NRM) 
!   -------------------------------------------------------------------
!         Compute the approximation ZL to the solution.  Since the
!         vector Z was used as work space, and the initial guess
!         of the linear iteration is zero, Z must be reset to zero.
!   -------------------------------------------------------------------
  200 CONTINUE 
      LL = LGMR 
      LLP1 = LL + 1 
      K = 1 
      IF (LLP1 > 0) THEN 
         R0(:LLP1) = 0.0D0 
         K = LLP1 + 1 
      ENDIF 
      R0(1) = R0NRM 
      CALL DHELS (HES, MAXLP1, LL, Q, R0) 
      K = 1 
      IF (N > 0) THEN 
         Z(:N) = 0.0D0 
         K = N + 1 
      ENDIF 
      DO I = 1, LL 
         CALL DAXPY (N, R0(I), V(1,I), 1, Z, 1) 
      END DO 
      IF (JSCAL==1 .OR. JSCAL==3) THEN 
         I = 1 
         IF (N > 0) THEN 
            Z(:N) = Z(:N)/SZ(:N) 
            I = N + 1 
         ENDIF 
      ENDIF 
      IF (JPRE > 0) THEN 
         CALL DCOPY (N, Z, 1, WK, 1) 
         CALL MSOLVE (N, WK, Z, NELT, IA, JA, A, ISYM, RPAR, IPAR) 
         NMSL = NMSL + 1 
      ENDIF 
      RETURN  
!------------- LAST LINE OF DPIGMR FOLLOWS ----------------------------
      END SUBROUTINE DPIGMR 
!DECK DRLCAL
      SUBROUTINE DRLCAL(N, KMP, LL, MAXL, V, Q, RL, SNORMW, PROD, R0NRM) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  15:55:30  12/05/98  
!...Switches: -xf
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER N, KMP, LL, MAXL 
      DOUBLE PRECISION SNORMW, PROD, R0NRM 
      DOUBLE PRECISION, DIMENSION(N,*) :: V 
      DOUBLE PRECISION, DIMENSION(*) :: Q 
      DOUBLE PRECISION, DIMENSION(N) :: RL 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, IP1, I2, K, LLM1, LLP1 
      DOUBLE PRECISION :: S, C, TEM 
!-----------------------------------------------
!
!         Internal Variables.
!
!
!***FIRST EXECUTABLE STATEMENT  DRLCAL
      IF (KMP == MAXL) THEN 
!
!         calculate RL.  Start by copying V(*,1) into RL.
!
         CALL DCOPY (N, V(1,1), 1, RL, 1) 
         LLM1 = LL - 1 
         DO I = 1, LLM1 
            IP1 = I + 1 
            I2 = I*2 
            S = Q(I2) 
            C = Q(I2-1) 
            K = 1 
            IF (N > 0) THEN 
               RL(:N) = S*RL(:N) + C*V(:N,IP1) 
               K = N + 1 
            ENDIF 
         END DO 
         S = Q(2*LL) 
         C = Q(2*LL-1)/SNORMW 
         LLP1 = LL + 1 
         K = 1 
         IF (N > 0) THEN 
            RL(:N) = S*RL(:N) + C*V(:N,LLP1) 
            K = N + 1 
         ENDIF 
      ENDIF 
!
!         When KMP < MAXL, RL vector already partially calculated.
!         Scale RL by R0NRM*PROD to obtain the residual RL.
!
      TEM = R0NRM*PROD 
      CALL DSCAL (N, TEM, RL, 1) 
      RETURN  
!------------- LAST LINE OF DRLCAL FOLLOWS ----------------------------
      END SUBROUTINE DRLCAL 
!DECK DXLCAL
      SUBROUTINE DXLCAL(N, LGMR, X, XL, ZL, HES, MAXLP1, Q, V, R0NRM, WK, SZ, &
         JSCAL, JPRE, MSOLVE, NMSL, RPAR, IPAR, NELT, IA, JA, A, ISYM) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  15:55:30  12/05/98  
!...Switches: -xf
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER N, LGMR, MAXLP1, JSCAL, JPRE, NMSL, NELT, ISYM 
      DOUBLE PRECISION R0NRM 
      INTEGER, DIMENSION(*) :: IPAR 
      INTEGER, DIMENSION(NELT) :: IA, JA 
      DOUBLE PRECISION, DIMENSION(N) :: X, XL, ZL 
      DOUBLE PRECISION, DIMENSION(MAXLP1,*) :: HES 
      DOUBLE PRECISION, DIMENSION(*) :: Q 
      DOUBLE PRECISION, DIMENSION(N,*) :: V 
      DOUBLE PRECISION, DIMENSION(N) :: WK 
      DOUBLE PRECISION, DIMENSION(*) :: SZ, RPAR 
      DOUBLE PRECISION, DIMENSION(NELT) :: A 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, K, LL, LLP1 
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      EXTERNAL MSOLVE 
!-----------------------------------------------
!
!         Internal variables.
!
!
!***FIRST EXECUTABLE STATEMENT  DXLCAL
      LL = LGMR 
      LLP1 = LL + 1 
      K = 1 
      IF (LLP1 > 0) THEN 
         WK(:LLP1) = 0.0D0 
         K = LLP1 + 1 
      ENDIF 
      WK(1) = R0NRM 
      CALL DHELS (HES, MAXLP1, LL, Q, WK) 
      K = 1 
      IF (N > 0) THEN 
         ZL(:N) = 0.0D0 
         K = N + 1 
      ENDIF 
      DO I = 1, LL 
         CALL DAXPY (N, WK(I), V(1,I), 1, ZL, 1) 
      END DO 
      IF (JSCAL==1 .OR. JSCAL==3) THEN 
         K = 1 
         IF (N > 0) THEN 
            ZL(:N) = ZL(:N)/SZ(:N) 
            K = N + 1 
         ENDIF 
      ENDIF 
      IF (JPRE > 0) THEN 
         CALL DCOPY (N, ZL, 1, WK, 1) 
         CALL MSOLVE (N, WK, ZL, NELT, IA, JA, A, ISYM, RPAR, IPAR) 
         NMSL = NMSL + 1 
      ENDIF 
!         calculate XL from X and ZL.
      K = 1 
      IF (N > 0) THEN 
         XL(:N) = X(:N) + ZL(:N) 
         K = N + 1 
      ENDIF 
      RETURN  
!------------- LAST LINE OF DXLCAL FOLLOWS ----------------------------
      END SUBROUTINE DXLCAL 
!DECK ISDGMR
      INTEGER FUNCTION ISDGMR (N, B, X, XL, NELT, IA, JA, A, ISYM, MSOLVE, NMSL&
         , ITOL, TOL, ITMAX, ITER, ERR, IUNIT, R, Z, DZ, RWORK, IWORK, RNRM, &
         BNRM, SB, SX, JSCAL, KMP, LGMR, MAXL, MAXLP1, V, Q, SNORMW, PROD, &
         R0NRM, HES, JPRE) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  15:55:30  12/05/98  
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE SOLBLK 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER N, NELT, ISYM, NMSL, ITOL, ITMAX, ITER, IUNIT, JSCAL, KMP, LGMR, &
         MAXL, MAXLP1, JPRE 
      DOUBLE PRECISION TOL, ERR, RNRM, BNRM, SNORMW, PROD, R0NRM 
      DOUBLE PRECISION,DIMENSION(*)::B,X,XL,IA,JA,A,R,Z,DZ,RWORK,IWORK,SB,SX 
      DOUBLE PRECISION, DIMENSION(N,*) :: V 
      DOUBLE PRECISION, DIMENSION(*) :: Q 
      DOUBLE PRECISION, DIMENSION(MAXLP1,MAXL) :: HES 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IELMAX, I 
      DOUBLE PRECISION :: DXNRM, SOLNRM, TEM, FUZZ, RATMAX, RAT 

      SAVE SOLNRM 
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      DOUBLE PRECISION , EXTERNAL :: DNRM2, D1MACH
      EXTERNAL MSOLVE 
!-----------------------------------------------
!
!***FIRST EXECUTABLE STATEMENT ISDGMR
      ISDGMR = 0 
!
!       Use input from DPIGMR to determine if stop conditions are met.
!
      IF (ITOL == 0) ERR = RNRM/BNRM 
      IF (ITOL>0 .AND. ITOL<=3) THEN 
!
!       Use DRLCAL to calculate the scaled residual vector.
!       Store answer in R.
!
         IF(LGMR/=0)CALL DRLCAL(N,KMP,LGMR,MAXL,V,Q,R,SNORMW,PROD,R0NRM) 
         IF (ITOL <= 2) THEN 
!         err = ||Residual||/||RightHandSide||(2-Norms).
            ERR = DNRM2(N,R,1)/BNRM 
!
!         Unscale R by R0NRM*PROD when KMP < MAXL.
!
            IF (KMP<MAXL .AND. LGMR/=0) THEN 
               TEM = 1.0D0/(R0NRM*PROD) 
               CALL DSCAL (N, TEM, R, 1) 
            ENDIF 
         ELSE IF (ITOL == 3) THEN 
!         err = Max |(Minv*Residual)(i)/x(i)|
!         When jpre .lt. 0, r already contains Minv*Residual.
            IF (JPRE > 0) THEN 
               CALL MSOLVE (N, R, DZ, NELT, IA, JA, A, ISYM, RWORK, IWORK) 
               NMSL = NMSL + 1 
            ENDIF 
!
!         Unscale R by R0NRM*PROD when KMP < MAXL.
!
            IF (KMP<MAXL .AND. LGMR/=0) THEN 
               TEM = 1.0D0/(R0NRM*PROD) 
               CALL DSCAL (N, TEM, R, 1) 
            ENDIF 
!
            FUZZ = D1MACH(1) 
            IELMAX = 1 
            RATMAX = ABS(DZ(1))/MAX(ABS(X(1)),FUZZ) 
            DO I = 2, N 
               RAT = ABS(DZ(I))/MAX(ABS(X(I)),FUZZ) 
               IF (RAT > RATMAX) THEN 
                  IELMAX = I 
                  RATMAX = RAT 
               ENDIF 
            END DO 
            ERR = RATMAX 
            IF (RATMAX <= TOL) ISDGMR = 1 
            IF (IUNIT > 0) WRITE (IUNIT, 1020) ITER, IELMAX, RATMAX 
            RETURN  
         ENDIF 
      ENDIF 
      IF (ITOL == 11) THEN 
!
!       Use DXLCAL to calculate the approximate solution XL.
!
         IF (LGMR/=0 .AND. ITER>0) THEN 
            CALL DXLCAL (N, LGMR, X, XL, XL, HES, MAXLP1, Q, V, R0NRM, DZ, SX, &
               JSCAL, JPRE, MSOLVE, NMSL, RWORK, IWORK, NELT, IA, JA, A, ISYM) 
         ELSE IF (ITER == 0) THEN 
!         Copy X to XL to check if initial guess is good enough.
            CALL DCOPY (N, X, 1, XL, 1) 
         ELSE 
!         Return since this is the first call to DPIGMR on a restart.
            RETURN  
         ENDIF 
!
         IF (JSCAL==0 .OR. JSCAL==2) THEN 
!         err = ||x-TrueSolution||/||TrueSolution||(2-Norms).
            IF (ITER == 0) SOLNRM = DNRM2(N,SOLN,1) 
            I = 1 
            IF (N > 0) THEN 
               DZ(:N) = XL(:N) - SOLN(:N) 
               I = N + 1 
            ENDIF 
            ERR = DNRM2(N,DZ,1)/SOLNRM 
         ELSE 
            IF (ITER == 0) THEN 
               I = 1 
               IF (N > 0) THEN 
                  SOLNRM = SUM((SX(:N)*SOLN(:N))**2) 
                  I = N + 1 
               ENDIF 
               SOLNRM = DSQRT(SOLNRM) 
            ENDIF 
            I = 1 
            IF (N > 0) THEN 
               DXNRM = SUM((SX(:N)*(XL(:N)-SOLN(:N)))**2) 
               I = N + 1 
            ENDIF 
            DXNRM = DSQRT(DXNRM) 
!         err = ||SX*(x-TrueSolution)||/||SX*TrueSolution|| (2-Norms).
            ERR = DXNRM/SOLNRM 
         ENDIF 
      ENDIF 
!
      IF (IUNIT /= 0) THEN 
         IF (ITER == 0) WRITE (IUNIT, 1000) N, ITOL, MAXL, KMP 
         WRITE (IUNIT, 1010) ITER, RNRM/BNRM, ERR 
      ENDIF 
      IF (ERR <= TOL) ISDGMR = 1 
!
      RETURN  
 1000 FORMAT(' Generalized Minimum Residual(',I3,I3,') for ','N, ITOL = ',I5,I5&
         ,/' ITER','   Natral Err Est','   Error Estimate') 
 1010 FORMAT(1X,I4,1X,E16.7,1X,E16.7) 
 1020 FORMAT(1X,' ITER = ',I5,' IELMAX = ',I5,' |R(IELMAX)/X(IELMAX)| = ',E12.5&
         ) 
!------------- LAST LINE OF ISDGMR FOLLOWS ----------------------------
      END FUNCTION ISDGMR 
