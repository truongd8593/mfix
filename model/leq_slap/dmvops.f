!DECK DSMV
      SUBROUTINE DSMV(N, X, Y, NELT, IA, JA, A, ISYM) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  15:55:30  12/05/98  
!...Switches: -xf
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER N, NELT, ISYM 
      INTEGER, DIMENSION(NELT) :: IA, JA 
      DOUBLE PRECISION, DIMENSION(N) :: X, Y 
      DOUBLE PRECISION, DIMENSION(NELT) :: A 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, ICOL, IBGN, IEND, IROW, JBGN, JEND, J 
!-----------------------------------------------
!
!         Zero out the result vector.
!***FIRST EXECUTABLE STATEMENT  DSMV
      I = 1 
      IF (N > 0) THEN 
         Y = 0.0D0 
         I = N + 1 
      ENDIF 
      DO ICOL = 1, N 
         IBGN = JA(ICOL) 
         IEND = JA(ICOL+1) - 1 
!LLL. OPTION ASSERT (NOHAZARD)
!DIR$ IVDEP
!VD$ NODEPCHK
         I = IBGN 
         IF (IEND - IBGN + 1 > 0) THEN 
            Y(IA(IBGN:IEND)) = Y(IA(IBGN:IEND)) + A(IBGN:IEND)*X(ICOL) 
            I = IEND + 1 
         ENDIF 
      END DO 
      IF (ISYM == 1) THEN 
!
!         The matrix is non-symmetric.  Need to get the other half in...
!         This loops assumes that the diagonal is the first entry in
!         each column.
!
!$omp    parallel do private(irow,jbgn,jend,j)
         DO IROW = 1, N 
            JBGN = JA(IROW) + 1 
            JEND = JA(IROW+1) - 1 
            IF (JBGN <= JEND) THEN 
               J = JBGN 
               IF (JEND - JBGN + 1 > 0) THEN 
                  Y(IROW) = Y(IROW) + SUM(A(JBGN:JEND)*X(IA(JBGN:JEND))) 
                  J = JEND + 1 
               ENDIF 
            ENDIF 
         END DO 
      ENDIF 
      RETURN  
!------------- LAST LINE OF DSMV FOLLOWS ----------------------------
      END SUBROUTINE DSMV 
!DECK DSMTV
      SUBROUTINE DSMTV(N, X, Y, NELT, IA, JA, A, ISYM) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  15:55:30  12/05/98  
!...Switches: -xf
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER N, NELT, ISYM 
      INTEGER, DIMENSION(NELT) :: IA, JA 
      DOUBLE PRECISION, DIMENSION(N) :: X, Y 
      DOUBLE PRECISION, DIMENSION(NELT) :: A 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, IROW, IBGN, IEND, ICOL, JBGN, JEND, J 
!-----------------------------------------------
!
!         Zero out the result vector.
!***FIRST EXECUTABLE STATEMENT  DSMTV
      I = 1 
      IF (N > 0) THEN 
         Y = 0.0D0 
         I = N + 1 
      ENDIF 
      DO IROW = 1, N 
         IBGN = JA(IROW) 
         IEND = JA(IROW+1) - 1 
!VD$ ASSOC
         I = IBGN 
         IF (IEND - IBGN + 1 > 0) THEN 
            Y(IROW) = Y(IROW) + SUM(A(IBGN:IEND)*X(IA(IBGN:IEND))) 
            I = IEND + 1 
         ENDIF 
      END DO 
      IF (ISYM == 1) THEN 
!
!         The matrix is non-symmetric.  Need to get the other half in...
!         This loops assumes that the diagonal is the first entry in
!         each column.
!
         DO ICOL = 1, N 
            JBGN = JA(ICOL) + 1 
            JEND = JA(ICOL+1) - 1 
            IF (JBGN <= JEND) THEN 
!LLL. OPTION ASSERT (NOHAZARD)
!DIR$ IVDEP
!VD$ NODEPCHK
               J = JBGN 
               IF (JEND - JBGN + 1 > 0) THEN 
                  Y(IA(JBGN:JEND)) = Y(IA(JBGN:JEND)) + A(JBGN:JEND)*X(ICOL) 
                  J = JEND + 1 
               ENDIF 
            ENDIF 
         END DO 
      ENDIF 
      RETURN  
!------------- LAST LINE OF DSMTV FOLLOWS ----------------------------
      END SUBROUTINE DSMTV 
!DECK DSDI
      SUBROUTINE DSDI(N, B, X, NELT, IA, JA, A, ISYM, RWORK, IWORK) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  15:55:30  12/05/98  
!...Switches: -xf
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER N, NELT, ISYM 
      INTEGER, DIMENSION(NELT) :: IA, JA 
      INTEGER, DIMENSION(10) :: IWORK 
      DOUBLE PRECISION, DIMENSION(N) :: B, X 
      DOUBLE PRECISION, DIMENSION(NELT) :: A 
      DOUBLE PRECISION, DIMENSION(*) :: RWORK 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: LOCD, I 
!-----------------------------------------------
!
!         Determine where the inverse of the diagonal
!         is in the work array and then scale by it.
!***FIRST EXECUTABLE STATEMENT  DSDI
      LOCD = IWORK(4) - 1 
      I = 1 
      IF (N > 0) THEN 
         X = RWORK(LOCD+1:N+LOCD)*B 
         I = N + 1 
      ENDIF 
      RETURN  
!------------- LAST LINE OF DSDI FOLLOWS ----------------------------
      END SUBROUTINE DSDI 
!DECK DSLI
      SUBROUTINE DSLI(N, B, X, NELT, IA, JA, A, ISYM, RWORK, IWORK) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  15:55:30  12/05/98  
!...Switches: -xf
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER N, NELT, ISYM 
      INTEGER, DIMENSION(NELT) :: IA, JA 
      INTEGER, DIMENSION(10) :: IWORK 
      DOUBLE PRECISION, DIMENSION(N) :: B, X 
      DOUBLE PRECISION, DIMENSION(NELT) :: A 
      DOUBLE PRECISION, DIMENSION(*) :: RWORK 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NEL, LOCIEL, LOCJEL, LOCEL 
!-----------------------------------------------
!***FIRST EXECUTABLE STATEMENT  DSLI
!
      NEL = IWORK(1) 
      LOCIEL = IWORK(2) 
      LOCJEL = IWORK(3) 
      LOCEL = IWORK(4) 
      CALL DSLI2 (N, B, X, NEL, IWORK(LOCIEL), IWORK(LOCJEL), RWORK(LOCEL)) 
!
      RETURN  
!------------- LAST LINE OF DSLI FOLLOWS ----------------------------
      END SUBROUTINE DSLI 
!DECK DSLI2
      SUBROUTINE DSLI2(N, B, X, NEL, IEL, JEL, EL) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  15:55:30  12/05/98  
!...Switches: -xf
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER N, NEL 
      INTEGER, DIMENSION(NEL) :: IEL, JEL 
      DOUBLE PRECISION, DIMENSION(N) :: B, X 
      DOUBLE PRECISION, DIMENSION(NEL) :: EL 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, ICOL, JBGN, JEND, J 
!-----------------------------------------------
!
!         Initialize the solution by copying the right hands side
!         into it.
!***FIRST EXECUTABLE STATEMENT  DSLI2
      I = 1 
      IF (N > 0) THEN 
         X = B 
         I = N + 1 
      ENDIF 
      DO ICOL = 1, N 
         X(ICOL) = X(ICOL)/EL(JEL(ICOL)) 
         JBGN = JEL(ICOL) + 1 
         JEND = JEL(ICOL+1) - 1 
         IF (JBGN <= JEND) THEN 
!LLL. OPTION ASSERT (NOHAZARD)
!DIR$ IVDEP
!VD$ NOCONCUR
!VD$ NODEPCHK
            J = JBGN 
            IF (JEND - JBGN + 1 > 0) THEN 
               X(IEL(JBGN:JEND)) = X(IEL(JBGN:JEND)) - EL(JBGN:JEND)*X(ICOL) 
               J = JEND + 1 
            ENDIF 
         ENDIF 
      END DO 
      RETURN  
!------------- LAST LINE OF DSLI2 FOLLOWS ----------------------------
      END SUBROUTINE DSLI2 
!DECK DSLLTI
      SUBROUTINE DSLLTI(N, B, X, NELT, IA, JA, A, ISYM, RWORK, IWORK) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  15:55:30  12/05/98  
!...Switches: -xf
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER N, NELT, ISYM 
      INTEGER, DIMENSION(NELT) :: IA, JA 
      INTEGER, DIMENSION(*) :: IWORK 
      DOUBLE PRECISION, DIMENSION(*) :: B, X 
      DOUBLE PRECISION, DIMENSION(NELT) :: A 
      DOUBLE PRECISION, DIMENSION(*) :: RWORK 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NEL, LOCIEL, LOCJEL, LOCEL, LOCDIN 
!-----------------------------------------------
!
!***FIRST EXECUTABLE STATEMENT  DSLLTI
      NEL = IWORK(1) 
      LOCIEL = IWORK(3) 
      LOCJEL = IWORK(2) 
      LOCEL = IWORK(4) 
      LOCDIN = IWORK(5) 
      CALL SLLTI2 (N, B, X, NEL, IWORK(LOCIEL), IWORK(LOCJEL), RWORK(LOCEL), &
         RWORK(LOCDIN)) 
!
      RETURN  
!------------- LAST LINE OF DSLLTI FOLLOWS ----------------------------
      END SUBROUTINE DSLLTI 
!DECK SLLTI2
      SUBROUTINE SLLTI2(N, B, X, NEL, IEL, JEL, EL, DINV) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  15:55:30  12/05/98  
!...Switches: -xf
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER N, NEL 
      INTEGER, DIMENSION(NEL) :: IEL 
      INTEGER, DIMENSION(*) :: JEL 
      DOUBLE PRECISION, DIMENSION(N) :: B, X 
      DOUBLE PRECISION, DIMENSION(NEL) :: EL 
      DOUBLE PRECISION, DIMENSION(N) :: DINV 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, IROW, IBGN, IEND 
!-----------------------------------------------
!
!         solve  l*y = b,  storing result in x.
!***FIRST EXECUTABLE STATEMENT  SLLTI2
      I = 1 
      IF (N > 0) THEN 
         X = B 
         I = N + 1 
      ENDIF 
      DO IROW = 1, N 
         IBGN = IEL(IROW) + 1 
         IEND = IEL(IROW+1) - 1 
         IF (IBGN <= IEND) THEN 
!LLL. OPTION ASSERT (NOHAZARD)
!DIR$ IVDEP
!VD$ NOCONCUR
!VD$ NODEPCHK
            I = IBGN 
            IF (IEND - IBGN + 1 > 0) THEN 
               X(IROW) = X(IROW) - SUM(EL(IBGN:IEND)*X(JEL(IBGN:IEND))) 
               I = IEND + 1 
            ENDIF 
         ENDIF 
      END DO 
      I = 1 
      IF (N > 0) THEN 
         X = X*DINV 
         I = N + 1 
      ENDIF 
      DO IROW = N, 2, -1 
         IBGN = IEL(IROW) + 1 
         IEND = IEL(IROW+1) - 1 
         IF (IBGN <= IEND) THEN 
!LLL. OPTION ASSERT (NOHAZARD)
!DIR$ IVDEP
!VD$ NOCONCUR
!VD$ NODEPCHK
            I = IBGN 
            IF (IEND - IBGN + 1 > 0) THEN 
               X(JEL(IBGN:IEND)) = X(JEL(IBGN:IEND)) - EL(IBGN:IEND)*X(IROW) 
               I = IEND + 1 
            ENDIF 
         ENDIF 
      END DO 
      RETURN  
!------------- LAST LINE OF SLTI2 FOLLOWS ----------------------------
      END SUBROUTINE SLLTI2 
!DECK DSLUI
      SUBROUTINE DSLUI(N, B, X, NELT, IA, JA, A, ISYM, RWORK, IWORK) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  15:55:30  12/05/98  
!...Switches: -xf
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER N, NELT, ISYM 
      INTEGER, DIMENSION(NELT) :: IA, JA 
      INTEGER, DIMENSION(10) :: IWORK 
      DOUBLE PRECISION, DIMENSION(N) :: B, X 
      DOUBLE PRECISION, DIMENSION(NELT) :: A 
      DOUBLE PRECISION, DIMENSION(*) :: RWORK 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: LOCIL, LOCJL, LOCIU, LOCJU, LOCL, LOCDIN, LOCU 
!-----------------------------------------------
!
!         Pull out the locations of the arrays holding the ILU
!         factorization.
!***FIRST EXECUTABLE STATEMENT  DSLUI
      LOCIL = IWORK(1) 
      LOCJL = IWORK(2) 
      LOCIU = IWORK(3) 
      LOCJU = IWORK(4) 
      LOCL = IWORK(5) 
      LOCDIN = IWORK(6) 
      LOCU = IWORK(7) 
!
!         Solve the system LUx = b
      CALL DSLUI2 (N, B, X, IWORK(LOCIL), IWORK(LOCJL), RWORK(LOCL), RWORK(&
         LOCDIN), IWORK(LOCIU), IWORK(LOCJU), RWORK(LOCU)) 
!
      RETURN  
!------------- LAST LINE OF DSLUI FOLLOWS ----------------------------
      END SUBROUTINE DSLUI 
!DECK DSLUI2
      SUBROUTINE DSLUI2(N, B, X, IL, JL, L, DINV, IU, JU, U) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  15:55:30  12/05/98  
!...Switches: -xf
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER N 
      INTEGER, DIMENSION(*) :: IL, JL, IU, JU 
      DOUBLE PRECISION, DIMENSION(N) :: B, X 
      DOUBLE PRECISION, DIMENSION(*) :: L 
      DOUBLE PRECISION, DIMENSION(N) :: DINV 
      DOUBLE PRECISION, DIMENSION(*) :: U 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, IROW, JBGN, JEND, J, ICOL 
!-----------------------------------------------
!
!         Solve  L*Y = B,  storing result in X, L stored by rows.
!***FIRST EXECUTABLE STATEMENT  DSLUI2
      I = 1 
      IF (N > 0) THEN 
         X = B 
         I = N + 1 
      ENDIF 
      DO IROW = 2, N 
         JBGN = IL(IROW) 
         JEND = IL(IROW+1) - 1 
         IF (JBGN <= JEND) THEN 
!LLL. OPTION ASSERT (NOHAZARD)
!DIR$ IVDEP
!VD$ ASSOC
!VD$ NODEPCHK
            J = JBGN 
            IF (JEND - JBGN + 1 > 0) THEN 
               X(IROW) = X(IROW) - SUM(L(JBGN:JEND)*X(JL(JBGN:JEND))) 
               J = JEND + 1 
            ENDIF 
         ENDIF 
      END DO 
      I = 1 
      IF (N > 0) THEN 
         X = X*DINV 
         I = N + 1 
      ENDIF 
      DO ICOL = N, 2, -1 
         JBGN = JU(ICOL) 
         JEND = JU(ICOL+1) - 1 
         IF (JBGN <= JEND) THEN 
!LLL. OPTION ASSERT (NOHAZARD)
!DIR$ IVDEP
!VD$ NODEPCHK
            J = JBGN 
            IF (JEND - JBGN + 1 > 0) THEN 
               X(IU(JBGN:JEND)) = X(IU(JBGN:JEND)) - U(JBGN:JEND)*X(ICOL) 
               J = JEND + 1 
            ENDIF 
         ENDIF 
      END DO 
      RETURN  
!------------- LAST LINE OF DSLUI2 FOLLOWS ----------------------------
      END SUBROUTINE DSLUI2 
!DECK DSLUTI
      SUBROUTINE DSLUTI(N, B, X, NELT, IA, JA, A, ISYM, RWORK, IWORK) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  15:55:30  12/05/98  
!...Switches: -xf
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER N, NELT, ISYM 
      INTEGER, DIMENSION(NELT) :: IA, JA 
      INTEGER, DIMENSION(10) :: IWORK 
      DOUBLE PRECISION, DIMENSION(N) :: B, X, A 
      DOUBLE PRECISION, DIMENSION(*) :: RWORK 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: LOCIL, LOCJL, LOCIU, LOCJU, LOCL, LOCDIN, LOCU 
!-----------------------------------------------
!
!         Pull out the pointers to the L, D and U matricies and call
!         the workhorse routine.
!***FIRST EXECUTABLE STATEMENT  DSLUTI
      LOCIL = IWORK(1) 
      LOCJL = IWORK(2) 
      LOCIU = IWORK(3) 
      LOCJU = IWORK(4) 
      LOCL = IWORK(5) 
      LOCDIN = IWORK(6) 
      LOCU = IWORK(7) 
!
      CALL DSLUI4 (N, B, X, IWORK(LOCIL), IWORK(LOCJL), RWORK(LOCL), RWORK(&
         LOCDIN), IWORK(LOCIU), IWORK(LOCJU), RWORK(LOCU)) 
!
      RETURN  
!------------- LAST LINE OF DSLUTI FOLLOWS ----------------------------
      END SUBROUTINE DSLUTI 
!DECK DSLUI4
      SUBROUTINE DSLUI4(N, B, X, IL, JL, L, DINV, IU, JU, U) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  15:55:30  12/05/98  
!...Switches: -xf
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER N 
      INTEGER, DIMENSION(*) :: IL, JL, IU, JU 
      DOUBLE PRECISION, DIMENSION(N) :: B, X 
      DOUBLE PRECISION, DIMENSION(*) :: L 
      DOUBLE PRECISION, DIMENSION(N) :: DINV 
      DOUBLE PRECISION, DIMENSION(*) :: U 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, IROW, JBGN, JEND, J, ICOL 
!-----------------------------------------------
!
!***FIRST EXECUTABLE STATEMENT  DSLUI4
      I = 1 
      IF (N > 0) THEN 
         X = B 
         I = N + 1 
      ENDIF 
      DO IROW = 2, N 
         JBGN = JU(IROW) 
         JEND = JU(IROW+1) - 1 
         IF (JBGN <= JEND) THEN 
!LLL. OPTION ASSERT (NOHAZARD)
!DIR$ IVDEP
!VD$ ASSOC
!VD$ NODEPCHK
            J = JBGN 
            IF (JEND - JBGN + 1 > 0) THEN 
               X(IROW) = X(IROW) - SUM(U(JBGN:JEND)*X(IU(JBGN:JEND))) 
               J = JEND + 1 
            ENDIF 
         ENDIF 
      END DO 
      I = 1 
      IF (N > 0) THEN 
         X = X*DINV 
         I = N + 1 
      ENDIF 
      DO ICOL = N, 2, -1 
         JBGN = IL(ICOL) 
         JEND = IL(ICOL+1) - 1 
         IF (JBGN <= JEND) THEN 
!LLL. OPTION ASSERT (NOHAZARD)
!DIR$ IVDEP
!VD$ NODEPCHK
            J = JBGN 
            IF (JEND - JBGN + 1 > 0) THEN 
               X(JL(JBGN:JEND)) = X(JL(JBGN:JEND)) - L(JBGN:JEND)*X(ICOL) 
               J = JEND + 1 
            ENDIF 
         ENDIF 
      END DO 
      RETURN  
!------------- LAST LINE OF DSLUI4 FOLLOWS ----------------------------
      END SUBROUTINE DSLUI4 
!DECK DSMMTI
      SUBROUTINE DSMMTI(N, B, X, NELT, IA, JA, A, ISYM, RWORK, IWORK) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  15:55:30  12/05/98  
!...Switches: -xf
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER N, NELT, ISYM 
      INTEGER, DIMENSION(NELT) :: IA, JA 
      INTEGER, DIMENSION(10) :: IWORK 
      DOUBLE PRECISION, DIMENSION(N) :: B, X 
      DOUBLE PRECISION, DIMENSION(NELT) :: A 
      DOUBLE PRECISION, DIMENSION(*) :: RWORK 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: LOCIL, LOCJL, LOCIU, LOCJU, LOCL, LOCDIN, LOCU 
!-----------------------------------------------
!
!         Pull out the locations of the arrays holding the ILU
!         factorization.
!***FIRST EXECUTABLE STATEMENT  DSMMTI
      LOCIL = IWORK(1) 
      LOCJL = IWORK(2) 
      LOCIU = IWORK(3) 
      LOCJU = IWORK(4) 
      LOCL = IWORK(5) 
      LOCDIN = IWORK(6) 
      LOCU = IWORK(7) 
!
      CALL DSMMI2 (N, B, X, IWORK(LOCIL), IWORK(LOCJL), RWORK(LOCL), RWORK(&
         LOCDIN), IWORK(LOCIU), IWORK(LOCJU), RWORK(LOCU)) 
!
      RETURN  
!------------- LAST LINE OF DSMMTI FOLLOWS ----------------------------
      END SUBROUTINE DSMMTI 
!DECK DSMMI2
      SUBROUTINE DSMMI2(N, B, X, IL, JL, L, DINV, IU, JU, U) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  15:55:30  12/05/98  
!...Switches: -xf
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER N 
      INTEGER, DIMENSION(1) :: IL, JL, IU, JU 
      DOUBLE PRECISION, DIMENSION(N) :: B, X 
      DOUBLE PRECISION, DIMENSION(*) :: L 
      DOUBLE PRECISION, DIMENSION(N) :: DINV, U 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, IROW, JBGN, JEND, J, ICOL 
!-----------------------------------------------
!
!         Solve  L*Y = B,  storing result in X, L stored by rows.
!***FIRST EXECUTABLE STATEMENT  DSMMI2
      I = 1 
      IF (N > 0) THEN 
         X = B 
         I = N + 1 
      ENDIF 
      DO IROW = 2, N 
         JBGN = IL(IROW) 
         JEND = IL(IROW+1) - 1 
         IF (JBGN <= JEND) THEN 
!LLL. OPTION ASSERT (NOHAZARD)
!DIR$ IVDEP
!VD$ ASSOC
!VD$ NODEPCHK
            J = JBGN 
            IF (JEND - JBGN + 1 > 0) THEN 
               X(IROW) = X(IROW) - SUM(L(JBGN:JEND)*X(JL(JBGN:JEND))) 
               J = JEND + 1 
            ENDIF 
         ENDIF 
      END DO 
      I = 1 
      IF (N > 0) THEN 
         X = X*DINV 
         I = N + 1 
      ENDIF 
      DO ICOL = N, 2, -1 
         JBGN = JU(ICOL) 
         JEND = JU(ICOL+1) - 1 
         IF (JBGN <= JEND) THEN 
!LLL. OPTION ASSERT (NOHAZARD)
!DIR$ IVDEP
!VD$ NODEPCHK
            J = JBGN 
            IF (JEND - JBGN + 1 > 0) THEN 
               X(IU(JBGN:JEND)) = X(IU(JBGN:JEND)) - U(JBGN:JEND)*X(ICOL) 
               J = JEND + 1 
            ENDIF 
         ENDIF 
      END DO 
      DO IROW = 2, N 
         JBGN = JU(IROW) 
         JEND = JU(IROW+1) - 1 
         IF (JBGN <= JEND) THEN 
!LLL. OPTION ASSERT (NOHAZARD)
!DIR$ IVDEP
!VD$ ASSOC
!VD$ NODEPCHK
            J = JBGN 
            IF (JEND - JBGN + 1 > 0) THEN 
               X(IROW) = X(IROW) - SUM(U(JBGN:JEND)*X(IU(JBGN:JEND))) 
               J = JEND + 1 
            ENDIF 
         ENDIF 
      END DO 
      I = 1 
      IF (N > 0) THEN 
         X = X*DINV 
         I = N + 1 
      ENDIF 
      DO ICOL = N, 2, -1 
         JBGN = IL(ICOL) 
         JEND = IL(ICOL+1) - 1 
         IF (JBGN <= JEND) THEN 
!LLL. OPTION ASSERT (NOHAZARD)
!DIR$ IVDEP
!VD$ NODEPCHK
            J = JBGN 
            IF (JEND - JBGN + 1 > 0) THEN 
               X(JL(JBGN:JEND)) = X(JL(JBGN:JEND)) - L(JBGN:JEND)*X(ICOL) 
               J = JEND + 1 
            ENDIF 
         ENDIF 
      END DO 
      RETURN  
!------------- LAST LINE OF DSMMI2 FOLLOWS ----------------------------
      END SUBROUTINE DSMMI2 
