!DECK DCGS
      SUBROUTINE DCGS(N, B, X, NELT, IA, JA, A, ISYM, MATVEC, MSOLVE, ITOL, TOL&
         , ITMAX, ITER, ERR, IERR, IUNIT, R, R0, P, Q, U, V1, V2, RWORK, IWORK) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  15:55:30  12/05/98  
!...Switches: -xf
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER N, NELT, ISYM, ITOL, ITMAX, ITER, IERR, IUNIT 
      DOUBLE PRECISION TOL, ERR 
      INTEGER, DIMENSION(NELT) :: IA, JA 
      INTEGER, DIMENSION(*) :: IWORK 
      DOUBLE PRECISION, DIMENSION(N) :: B, X 
      DOUBLE PRECISION, DIMENSION(NELT) :: A 
      DOUBLE PRECISION, DIMENSION(N) :: R, R0, P, Q, U, V1, V2 
      DOUBLE PRECISION, DIMENSION(*) :: RWORK 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, K 
      DOUBLE PRECISION::TOLMIN,AK,BK,BNRM,SOLNRM,FUZZ,RHONM1,RHON,SIGMA,AKM 
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      INTEGER , EXTERNAL :: ISDCGS 
      DOUBLE PRECISION , EXTERNAL :: D1MACH, DDOT 
      EXTERNAL MATVEC, MSOLVE 
!-----------------------------------------------
!
!         Check some of the input data.
!***FIRST EXECUTABLE STATEMENT  DCGS
      ITER = 0 
      IERR = 0 
      IF (N < 1) THEN 
         IERR = 3 
         RETURN  
      ENDIF 
      TOLMIN = 500.0*D1MACH(3) 
      IF (TOL < TOLMIN) THEN 
         TOL = TOLMIN 
         IERR = 4 
      ENDIF 
!
!         Calculate initial residual and pseudo-residual, and check
!         stopping criterion.
      CALL MATVEC (N, X, R, NELT, IA, JA, A, ISYM) 
      I = 1 
      IF (N > 0) THEN 
         V1(:N) = R(:N) - B(:N) 
         I = N + 1 
      ENDIF 
      CALL MSOLVE (N, V1, R, NELT, IA, JA, A, ISYM, RWORK, IWORK) 
!
      IF (ISDCGS(N,B,X,NELT,IA,JA,A,ISYM,MATVEC,MSOLVE,ITOL,TOL,ITMAX,ITER,ERR,&
         IERR,IUNIT,R,R0,P,Q,U,V1,V2,RWORK,IWORK,AK,BK,BNRM,SOLNRM) == 0) THEN 
         IF (IERR /= 0) RETURN  
!
!         Set initial values.
!
         FUZZ = D1MACH(3)**2 
         I = 1 
         IF (N > 0) THEN 
            R0(:N) = R(:N) 
            I = N + 1 
         ENDIF 
         RHONM1 = 1.0 
!
!         ***** ITERATION LOOP *****
!
         DO K = 1, ITMAX 
            ITER = K 
!
!         Calculate coefficient BK and direction vectors U, V and P.
            RHON = DDOT(N,R0,1,R,1) 
            IF (ABS(RHONM1) < FUZZ) GO TO 998 
            BK = RHON/RHONM1 
            IF (ITER == 1) THEN 
               I = 1 
               IF (N > 0) THEN 
                  U(:N) = R(:N) 
                  P(:N) = R(:N) 
                  I = N + 1 
               ENDIF 
            ELSE 
               I = 1 
               IF (N > 0) THEN 
                  U(:N) = R(:N) + BK*Q(:N) 
                  V1(:N) = Q(:N) + BK*P(:N) 
                  I = N + 1 
               ENDIF 
               I = 1 
               IF (N > 0) THEN 
                  P(:N) = U(:N) + BK*V1(:N) 
                  I = N + 1 
               ENDIF 
            ENDIF 
!
!         Calculate coefficient AK, new iterate X, Q
            CALL MATVEC (N, P, V2, NELT, IA, JA, A, ISYM) 
            CALL MSOLVE (N, V2, V1, NELT, IA, JA, A, ISYM, RWORK, IWORK) 
            SIGMA = DDOT(N,R0,1,V1,1) 
            IF (ABS(SIGMA) < FUZZ) GO TO 999 
            AK = RHON/SIGMA 
            AKM = -AK 
            I = 1 
            IF (N > 0) THEN 
               Q(:N) = U(:N) + AKM*V1(:N) 
               I = N + 1 
            ENDIF 
            I = 1 
            IF (N > 0) THEN 
               V1(:N) = U(:N) + Q(:N) 
               I = N + 1 
            ENDIF 
            CALL DAXPY (N, AKM, V1, 1, X, 1) 
!                     -1
!         R = R - ak*M  *A*V1
            CALL MATVEC (N, V1, V2, NELT, IA, JA, A, ISYM) 
            CALL MSOLVE (N, V2, V1, NELT, IA, JA, A, ISYM, RWORK, IWORK) 
            CALL DAXPY (N, AKM, V1, 1, R, 1) 
!
!         check stopping criterion.
            IF (ISDCGS(N,B,X,NELT,IA,JA,A,ISYM,MATVEC,MSOLVE,ITOL,TOL,ITMAX,&
               ITER,ERR,IERR,IUNIT,R,R0,P,Q,U,V1,V2,RWORK,IWORK,AK,BK,BNRM,&
               SOLNRM) /= 0) GO TO 200 
!
!         Update RHO.
            RHONM1 = RHON 
         END DO 
         ITER = ITMAX + 1 
         IERR = 2 
      ENDIF 
  200 CONTINUE 
      RETURN  
!
!         Breakdown of method detected.
  998 CONTINUE 
      IERR = 5 
      RETURN  
!
!         Stagnation of method detected.
  999 CONTINUE 
      IERR = 6 
      RETURN  
!------------- LAST LINE OF DCGS FOLLOWS ----------------------------
      END SUBROUTINE DCGS 
!DECK DSDCGS
      SUBROUTINE DSDCGS(N, B, X, NELT, IA, JA, A, ISYM, ITOL, TOL, ITMAX, ITER&
         , ERR, IERR, IUNIT, RWORK, LENW, IWORK, LENIW) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  15:55:30  12/05/98  
!...Switches: -xf
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER N, NELT, ISYM, ITOL, ITMAX, ITER, IERR, IUNIT, LENW, LENIW 
      DOUBLE PRECISION TOL, ERR 
      INTEGER, DIMENSION(NELT) :: IA, JA 
      INTEGER, DIMENSION(LENIW) :: IWORK 
      DOUBLE PRECISION, DIMENSION(N) :: B, X, A 
      DOUBLE PRECISION, DIMENSION(LENW) :: RWORK 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: LOCRB = 1 
      INTEGER, PARAMETER :: LOCIB = 11 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER::LOCIW,LOCDIN,LOCR,LOCR0,LOCP,LOCQ,LOCU,LOCV1,LOCV2,LOCW 
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      DOUBLE PRECISION , EXTERNAL :: DSMV, DSDI 
!-----------------------------------------------
      IERR = 0 
      IF (N<1 .OR. NELT<1) THEN 
         IERR = 3 
         RETURN  
      ENDIF 
      CALL DS2Y (N, NELT, IA, JA, A, ISYM) 
!
!         Set up the workspace.  Compute the inverse of the
!         diagonal of the matrix.
      LOCIW = LOCIB 
!
      LOCDIN = LOCRB 
      LOCR = LOCDIN + N 
      LOCR0 = LOCR + N 
      LOCP = LOCR0 + N 
      LOCQ = LOCP + N 
      LOCU = LOCQ + N 
      LOCV1 = LOCU + N 
      LOCV2 = LOCV1 + N 
      LOCW = LOCV2 + N 
!
!         Check the workspace allocations.
      CALL DCHKW ('DSDCGS', LOCIW, LENIW, LOCW, LENW, IERR, ITER, ERR) 
      IF (IERR /= 0) RETURN  
!
      IWORK(4) = LOCDIN 
      IWORK(9) = LOCIW 
      IWORK(10) = LOCW 
!
      CALL DSDS (N, NELT, IA, JA, A, ISYM, RWORK(LOCDIN)) 
!
!         Perform the Diagonally Scaled
!         BiConjugate Gradient Squared algorithm.
      CALL DCGS (N, B, X, NELT, IA, JA, A, ISYM, DSMV, DSDI, ITOL, TOL, ITMAX, &
         ITER, ERR, IERR, IUNIT, RWORK(LOCR), RWORK(LOCR0), RWORK(LOCP), RWORK(&
         LOCQ), RWORK(LOCU), RWORK(LOCV1), RWORK(LOCV2), RWORK(1), IWORK(1)) 
      RETURN  
!------------- LAST LINE OF DSDCGS FOLLOWS ----------------------------
      END SUBROUTINE DSDCGS 
!DECK DSLUCS
      SUBROUTINE DSLUCS(N, B, X, NELT, IA, JA, A, ISYM, ITOL, TOL, ITMAX, ITER&
         , ERR, IERR, IUNIT, RWORK, LENW, IWORK, LENIW) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  15:55:30  12/05/98  
!...Switches: -xf
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER N, NELT, ISYM, ITOL, ITMAX, ITER, IERR, IUNIT, LENW, LENIW 
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
      INTEGER :: NL, NU, ICOL, JBGN, JEND, J, LOCIL, LOCJL, LOCIU, LOCJU, LOCNR&
         , LOCNC, LOCIW, LOCL, LOCDIN, LOCUU, LOCR, LOCR0, LOCP, LOCQ, LOCU, &
         LOCV1, LOCV2, LOCW 
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      DOUBLE PRECISION , EXTERNAL :: DSMV, DSLUI 
!-----------------------------------------------
      IERR = 0 
      IF (N<1 .OR. NELT<1) THEN 
         IERR = 3 
         RETURN  
      ENDIF 
      CALL DS2Y (N, NELT, IA, JA, A, ISYM) 
!
!         Count number of Non-Zero elements preconditioner ILU matrix.
!         Then set up the work arrays.
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
      LOCIL = LOCIB 
      LOCJL = LOCIL + N + 1 
      LOCIU = LOCJL + NL 
      LOCJU = LOCIU + NU 
      LOCNR = LOCJU + N + 1 
      LOCNC = LOCNR + N 
      LOCIW = LOCNC + N 
!
      LOCL = LOCRB 
      LOCDIN = LOCL + NL 
      LOCUU = LOCDIN + N 
      LOCR = LOCUU + NU 
      LOCR0 = LOCR + N 
      LOCP = LOCR0 + N 
      LOCQ = LOCP + N 
      LOCU = LOCQ + N 
      LOCV1 = LOCU + N 
      LOCV2 = LOCV1 + N 
      LOCW = LOCV2 + N 
!
!         Check the workspace allocations.
      CALL DCHKW ('DSLUCS', LOCIW, LENIW, LOCW, LENW, IERR, ITER, ERR) 
      IF (IERR /= 0) RETURN  
!
      IWORK(1) = LOCIL 
      IWORK(2) = LOCJL 
      IWORK(3) = LOCIU 
      IWORK(4) = LOCJU 
      IWORK(5) = LOCL 
      IWORK(6) = LOCDIN 
      IWORK(7) = LOCUU 
      IWORK(9) = LOCIW 
      IWORK(10) = LOCW 
!
!         Compute the Incomplete LU decomposition.
      CALL DSILUS (N, NELT, IA, JA, A, ISYM, NL, IWORK(LOCIL), IWORK(LOCJL), &
         RWORK(LOCL), RWORK(LOCDIN), NU, IWORK(LOCIU), IWORK(LOCJU), RWORK(&
         LOCUU), IWORK(LOCNR), IWORK(LOCNC)) 
!
!         Perform the incomplete LU preconditioned
!         BiConjugate Gradient Squared algorithm.
      CALL DCGS (N, B, X, NELT, IA, JA, A, ISYM, DSMV, DSLUI, ITOL, TOL, ITMAX&
         , ITER, ERR, IERR, IUNIT, RWORK(LOCR), RWORK(LOCR0), RWORK(LOCP), &
         RWORK(LOCQ), RWORK(LOCU), RWORK(LOCV1), RWORK(LOCV2), RWORK, IWORK) 
      RETURN  
!------------- LAST LINE OF DSLUCS FOLLOWS ----------------------------
      END SUBROUTINE DSLUCS 
!DECK ISDCGS
      INTEGER FUNCTION ISDCGS (N, B, X, NELT, IA, JA, A, ISYM, MATVEC, MSOLVE, &
         ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, R, R0, P, Q, U, V1, V2, &
         RWORK, IWORK, AK, BK, BNRM, SOLNRM) 
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
      INTEGER N, NELT, ISYM, ITOL, ITMAX, ITER, IERR, IUNIT 
      DOUBLE PRECISION TOL, ERR, AK, BK, BNRM, SOLNRM 
      INTEGER, DIMENSION(NELT) :: IA, JA 
      INTEGER, DIMENSION(*) :: IWORK 
      DOUBLE PRECISION, DIMENSION(N) :: B, X 
      DOUBLE PRECISION, DIMENSION(NELT) :: A 
      DOUBLE PRECISION, DIMENSION(N) :: R, R0, P, Q, U, V1, V2 
      DOUBLE PRECISION, DIMENSION(*) :: RWORK 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I 
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      DOUBLE PRECISION , EXTERNAL :: DNRM2 
      EXTERNAL MATVEC, MSOLVE 
!-----------------------------------------------
!
!***FIRST EXECUTABLE STATEMENT  ISDCGS
      ISDCGS = 0 
!
      SELECT CASE (ITOL)  
      CASE (1)  
!         err = ||Residual||/||RightHandSide|| (2-Norms).
         IF (ITER == 0) BNRM = DNRM2(N,B,1) 
         CALL MATVEC (N, X, V2, NELT, IA, JA, A, ISYM) 
         I = 1 
         IF (N > 0) THEN 
            V2(:N) = V2(:N) - B(:N) 
            I = N + 1 
         ENDIF 
         ERR = DNRM2(N,V2,1)/BNRM 
      CASE (2)  
!                  -1              -1
!         err = ||M  Residual||/||M  RightHandSide|| (2-Norms).
         IF (ITER == 0) THEN 
            CALL MSOLVE (N, B, V2, NELT, IA, JA, A, ISYM, RWORK, IWORK) 
            BNRM = DNRM2(N,V2,1) 
         ENDIF 
         ERR = DNRM2(N,R,1)/BNRM 
      CASE (11)  
!         err = ||x-TrueSolution||/||TrueSolution|| (2-Norms).
         IF (ITER == 0) SOLNRM = DNRM2(N,SOLN,1) 
         I = 1 
         IF (N > 0) THEN 
            V2(:N) = X(:N) - SOLN(:N) 
            I = N + 1 
         ENDIF 
         ERR = DNRM2(N,V2,1)/SOLNRM 
      CASE DEFAULT 
!
!         If we get here ITOL is not one of the acceptable values.
         ERR = 1.0E10 
         IERR = 3 
      END SELECT 
!
!         Print the error and Coeficients AK, BK on each step,
!         if desired.
      IF (IUNIT /= 0) THEN 
         IF (ITER == 0) WRITE (IUNIT, 1000) N, ITOL 
         WRITE (IUNIT, 1010) ITER, ERR, AK, BK 
      ENDIF 
      IF (ERR <= TOL) ISDCGS = 1 
!
      RETURN  
 1000 FORMAT(' Preconditioned BiConjugate Gradient Squared for ','N, ITOL = ',&
         I5,I5,/' ITER','   Error Estimate','            Alpha',&
         '             Beta') 
 1010 FORMAT(1X,I4,1X,E16.7,1X,E16.7,1X,E16.7) 
!------------- LAST LINE OF ISDCGS FOLLOWS ----------------------------
      END FUNCTION ISDCGS 
