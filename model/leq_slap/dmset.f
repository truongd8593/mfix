!DECK DSDS
      SUBROUTINE DSDS(N, NELT, IA, JA, A, ISYM, DINV) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  15:55:30  12/05/98  
!...Switches: -xf
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER N, NELT, ISYM 
      INTEGER, DIMENSION(NELT) :: IA, JA 
      DOUBLE PRECISION, DIMENSION(NELT) :: A 
      DOUBLE PRECISION, DIMENSION(N) :: DINV 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: ICOL 
!-----------------------------------------------
!
!         Assume the Diagonal elements are the first in each column.
!         This loop should *VECTORIZE*.  If it does not you may have
!         to add a compiler directive.  We do not check for a zero
!         (or near zero) diagonal element since this would interfere
!         with vectorization.  If this makes you nervous put a check
!         in!  It will run much slower.
!***FIRST EXECUTABLE STATEMENT  DSDS
      ICOL = 1 
      IF (N > 0) THEN 
         DINV = 1.0D0/A(JA(:N)) 
         ICOL = N + 1 
      ENDIF 
      RETURN  
!------------- LAST LINE OF DSDS FOLLOWS ----------------------------
      END SUBROUTINE DSDS 
!DECK DSDSCL
      SUBROUTINE DSDSCL(N, NELT, IA, JA, A, ISYM, X, B, DINV, JOB, ITOL) 
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
      INTEGER N, NELT, ISYM, JOB, ITOL 
      INTEGER, DIMENSION(NELT) :: IA, JA 
      DOUBLE PRECISION, DIMENSION(NELT) :: A 
      DOUBLE PRECISION, DIMENSION(N) :: X, B, DINV 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: ICOL, JBGN, JEND, J 
      DOUBLE PRECISION :: DI 
!-----------------------------------------------
!
!         SCALING...
!
      IF (JOB /= 0) THEN 
         ICOL = 1 
         IF (N > 0) THEN 
            DINV = 1.0D0/SQRT(A(JA(:N))) 
            ICOL = N + 1 
         ENDIF 
      ELSE 
!
!         UNSCALING...
!
         ICOL = 1 
         IF (N > 0) THEN 
            DINV = 1.0D0/DINV 
            ICOL = N + 1 
         ENDIF 
      ENDIF 
!
      DO ICOL = 1, N 
         JBGN = JA(ICOL) 
         JEND = JA(ICOL+1) - 1 
         DI = DINV(ICOL) 
         J = JBGN 
         IF (JEND - JBGN + 1 > 0) THEN 
            A(JBGN:JEND) = DINV(IA(JBGN:JEND))*A(JBGN:JEND)*DI 
            J = JEND + 1 
         ENDIF 
      END DO 
      ICOL = 1 
      IF (N > 0) THEN 
         B = B*DINV 
         X = X/DINV 
         ICOL = N + 1 
      ENDIF 
      IF (ITOL == 11) THEN 
         ICOL = 1 
         IF (N > 0) THEN 
            SOLN(:N) = SOLN(:N)/DINV 
            ICOL = N + 1 
         ENDIF 
      ENDIF 
!
      RETURN  
      END SUBROUTINE DSDSCL 
!DECK DSD2S
      SUBROUTINE DSD2S(N, NELT, IA, JA, A, ISYM, DINV) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  15:55:30  12/05/98  
!...Switches: -xf
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER N, NELT, ISYM 
      INTEGER, DIMENSION(NELT) :: IA, JA 
      DOUBLE PRECISION, DIMENSION(NELT) :: A 
      DOUBLE PRECISION, DIMENSION(N) :: DINV 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, KBGN, KEND, K 
!-----------------------------------------------
!
!***FIRST EXECUTABLE STATEMENT  DSD2S
      I = 1 
      IF (N > 0) THEN 
         DINV = 0. 
         I = N + 1 
      ENDIF 
      DO I = 1, N 
         KBGN = JA(I) 
         KEND = JA(I+1) - 1 
!
!         Add in the contributions for each row that has a non-zero
!         in this column.
!LLL. OPTION ASSERT (NOHAZARD)
!DIR$ IVDEP
!VD$ NODEPCHK
         K = KBGN 
         IF (KEND - KBGN + 1 > 0) THEN 
            DINV(IA(KBGN:KEND)) = DINV(IA(KBGN:KEND)) + A(KBGN:KEND)**2 
            K = KEND + 1 
         ENDIF 
         IF (ISYM == 1) THEN 
!
!         Lower triangle stored by columns => upper triangle stored by
!         rows with Diagonal being the first entry.  Loop across the
!         rest of the row.
            KBGN = KBGN + 1 
            IF (KBGN <= KEND) THEN 
               K = KBGN 
               IF (KEND - KBGN + 1 > 0) THEN 
                  DINV(I) = DINV(I) + SUM(A(KBGN:KEND)**2) 
                  K = KEND + 1 
               ENDIF 
            ENDIF 
         ENDIF 
      END DO 
      I = 1 
      IF (N > 0) THEN 
         DINV = 1./DINV 
         I = N + 1 
      ENDIF 
      RETURN  
!------------- LAST LINE OF DSD2S FOLLOWS ----------------------------
      END SUBROUTINE DSD2S 
!DECK DS2LT
      SUBROUTINE DS2LT(N, NELT, IA, JA, A, ISYM, NEL, IEL, JEL, EL) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  15:55:30  12/05/98  
!...Switches: -xf
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER N, NELT, ISYM, NEL 
      INTEGER, DIMENSION(NELT) :: IA, JA 
      INTEGER, DIMENSION(NEL) :: IEL, JEL 
      DOUBLE PRECISION, DIMENSION(NELT) :: A, EL 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: ICOL, JBGN, JEND, J, I 
!-----------------------------------------------
!***FIRST EXECUTABLE STATEMENT  DS2LT
      IF (ISYM == 0) THEN 
!
!         The matrix is stored non-symmetricly.  Pick out the lower
!         triangle.
!
         NEL = 0 
         DO ICOL = 1, N 
            JEL(ICOL) = NEL + 1 
            JBGN = JA(ICOL) 
            JEND = JA(ICOL+1) - 1 
!VD$ NOVECTOR
            DO J = JBGN, JEND 
               IF (IA(J) >= ICOL) THEN 
                  NEL = NEL + 1 
                  IEL(NEL) = IA(J) 
                  EL(NEL) = A(J) 
               ENDIF 
            END DO 
         END DO 
         JEL(N+1) = NEL + 1 
      ELSE 
!
!         The matrix is symmetric and only the lower triangle is
!         stored.  Copy it to IEL, JEL, EL.
!
         NEL = NELT 
         I = 1 
         IF (NELT > 0) THEN 
            IEL(:NELT) = IA 
            EL = A 
            I = NELT + 1 
         ENDIF 
         I = 1 
         IF (N + 1 > 0) THEN 
            JEL(:N+1) = JA(:N+1) 
            I = N + 2 
         ENDIF 
      ENDIF 
      RETURN  
!------------- LAST LINE OF DS2LT FOLLOWS ----------------------------
      END SUBROUTINE DS2LT 
!DECK DSICS
      SUBROUTINE DSICS(N,NELT,IA,JA,A,ISYM,NEL,IEL,JEL,EL,D,R,IWARN) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  15:55:30  12/05/98  
!...Switches: -xf
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER N, NELT, ISYM, NEL, IWARN 
      INTEGER, DIMENSION(NELT) :: IA, JA 
      INTEGER, DIMENSION(NEL) :: IEL, JEL 
      DOUBLE PRECISION, DIMENSION(NELT) :: A 
      DOUBLE PRECISION, DIMENSION(NEL) :: EL 
      DOUBLE PRECISION, DIMENSION(N) :: D, R 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IROW, ICBGN, ICEND, IC, ICOL, JBGN, JEND, J, IBGN, IEND, I, &
         JELTMP, IRBGN, IREND, IRR, IR 
      DOUBLE PRECISION :: ELTMP 
!-----------------------------------------------
!
!         Set the lower triangle in IEL, JEL, EL
!***FIRST EXECUTABLE STATEMENT  DSICS
      IWARN = 0 
!
!         All matrix elements stored in IA, JA, A.  Pick out the lower
!         triangle (making sure that the Diagonal of EL is one) and
!         store by rows.
!
      NEL = 1 
      IEL(1) = 1 
      JEL(1) = 1 
      EL(1) = 1.0D0 
      D(1) = A(1) 
!VD$R NOCONCUR
      DO IROW = 2, N 
!         Put in the Diagonal.
         NEL = NEL + 1 
         IEL(IROW) = NEL 
         JEL(NEL) = IROW 
         EL(NEL) = 1.0D0 
         D(IROW) = A(JA(IROW)) 
!
!         Look in all the lower triangle columns for a matching row.
!         Since the matrix is symmetric, we can look across the
!         irow-th row by looking down the irow-th column (if it is
!         stored ISYM=0)...
         IF (ISYM == 0) THEN 
            ICBGN = JA(IROW) 
            ICEND = JA(IROW+1) - 1 
         ELSE 
            ICBGN = 1 
            ICEND = IROW - 1 
         ENDIF 
         L20: DO IC = ICBGN, ICEND 
            IF (ISYM == 0) THEN 
               ICOL = IA(IC) 
               IF (ICOL >= IROW) CYCLE  L20 
            ELSE 
               ICOL = IC 
            ENDIF 
            JBGN = JA(ICOL) + 1 
            JEND = JA(ICOL+1) - 1 
            IF (JBGN<=JEND .AND. IA(JEND)>=IROW) THEN 
!VD$ NOVECTOR
               DO J = JBGN, JEND 
                  IF (IA(J) == IROW) THEN 
                     NEL = NEL + 1 
                     JEL(NEL) = ICOL 
                     EL(NEL) = A(J) 
                     CYCLE  L20 
                  ENDIF 
               END DO 
            ENDIF 
         END DO L20 
      END DO 
      IEL(N+1) = NEL + 1 
!
!         Sort ROWS of lower triangle into descending order (count out
!         along rows out from Diagonal).
!
      DO IROW = 2, N 
         IBGN = IEL(IROW) + 1 
         IEND = IEL(IROW+1) - 1 
         IF (IBGN < IEND) THEN 
            DO I = IBGN, IEND - 1 
!VD$ NOVECTOR
               DO J = I + 1, IEND 
                  IF (JEL(I) > JEL(J)) THEN 
                     JELTMP = JEL(J) 
                     JEL(J) = JEL(I) 
                     JEL(I) = JELTMP 
                     ELTMP = EL(J) 
                     EL(J) = EL(I) 
                     EL(I) = ELTMP 
                  ENDIF 
               END DO 
            END DO 
         ENDIF 
      END DO 
      IRBGN = JA(1) + 1 
      IREND = JA(2) - 1 
      DO IRR = IRBGN, IREND 
         IR = IA(IRR) 
!         Find the index into EL for EL(1,IR).
!         Hint: it's the second entry.
         I = IEL(IR) + 1 
         EL(I) = EL(I)/D(1) 
      END DO 
      DO IROW = 2, N 
!
!         Update the IROW-th diagonal.
!
         I = 1 
         IF (IROW - 1 > 0) THEN 
            R(:IROW-1) = 0.0D0 
            I = IROW 
         ENDIF 
         IBGN = IEL(IROW) + 1 
         IEND = IEL(IROW+1) - 1 
         IF (IBGN <= IEND) THEN 
!LLL. OPTION ASSERT (NOHAZARD)
!DIR$ IVDEP
!VD$ NODEPCHK
            I = IBGN 
            IF (IEND - IBGN + 1 > 0) THEN 
               R(JEL(IBGN:IEND)) = EL(IBGN:IEND)*D(JEL(IBGN:IEND)) 
               D(IROW) = D(IROW) - SUM(EL(IBGN:IEND)*R(JEL(IBGN:IEND))) 
               I = IEND + 1 
            ENDIF 
            IF (D(IROW) <= 0.0D0) THEN 
               IF (IWARN == 0) IWARN = IROW 
               D(IROW) = 1.0D0 
            ENDIF 
         ENDIF 
!
!         Update each EL(IROW+1:N,IROW), if there are any.
!         Use the structure of A to determine the Non-zero elements
!         of the IROW-th column of EL.
!
         IRBGN = JA(IROW) 
         IREND = JA(IROW+1) - 1 
         L100: DO IRR = IRBGN, IREND 
            IR = IA(IRR) 
            IF (IR > IROW) THEN 
!         Find the index into EL for EL(IR,IROW)
               IBGN = IEL(IR) + 1 
               IEND = IEL(IR+1) - 1 
               IF (JEL(IBGN) <= IROW) THEN 
                  DO I = IBGN, IEND 
                     IF (JEL(I) == IROW) THEN 
                        ICEND = IEND 
   91                   CONTINUE 
                        IF (JEL(ICEND) >= IROW) THEN 
                           ICEND = ICEND - 1 
                           GO TO 91 
                        ENDIF 
!         Sum up the EL(IR,1:IROW-1)*R(1:IROW-1) contributions.
!LLL. OPTION ASSERT (NOHAZARD)
!DIR$ IVDEP
!VD$ NODEPCHK
                        IC = IBGN 
                        IF (ICEND - IBGN + 1 > 0) THEN 
                           EL(I)=EL(I)-SUM(EL(IBGN:ICEND)*R(JEL(IBGN:ICEND))) 
                           IC = ICEND + 1 
                        ENDIF 
                        EL(I) = EL(I)/D(IROW) 
                        CYCLE  L100 
                     ENDIF 
                  END DO 
                  CALL XERRWV (&
                     'DSICS -- A and EL data structure mismatch in row (i1)', &
                     53, 1, 2, 1, IROW, 0, 0, 0.0, 0.0) 
               ENDIF 
            ENDIF 
         END DO L100 
      END DO 
      I = 1 
      IF (N > 0) THEN 
         D = 1.0D0/D 
         I = N + 1 
      ENDIF 
      RETURN  
!------------- LAST LINE OF DSICS FOLLOWS ----------------------------
      END SUBROUTINE DSICS 
!DECK DSILUS
      SUBROUTINE DSILUS(N, NELT, IA, JA, A, ISYM, NL, IL, JL, L, DINV, NU, IU, &
         JU, U, NROW, NCOL) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  15:55:30  12/05/98  
!...Switches: -xf
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER N, NELT, ISYM, NL, NU 
      INTEGER, DIMENSION(NELT) :: IA, JA 
      INTEGER, DIMENSION(NL) :: IL, JL 
      INTEGER, DIMENSION(NU) :: IU, JU 
      INTEGER, DIMENSION(N) :: NROW, NCOL 
      DOUBLE PRECISION, DIMENSION(NELT) :: A 
      DOUBLE PRECISION, DIMENSION(NL) :: L 
      DOUBLE PRECISION, DIMENSION(N) :: DINV 
      DOUBLE PRECISION, DIMENSION(NU) :: U 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, ICOL, JBGN, JEND, J, IROW, K, ITEMP, IBGN, IEND, JTEMP, &
         INDX1, INDX2, INDX, INDXR1, INDXR2, INDXC1, INDXC2, KR, KC 
      DOUBLE PRECISION :: TEMP 
!-----------------------------------------------
!
!         Count number of elements in each row of the lower triangle.
!***FIRST EXECUTABLE STATEMENT  DSILUS
      I = 1 
      IF (N > 0) THEN 
         NROW = 0 
         NCOL = 0 
         I = N + 1 
      ENDIF 
!VD$R NOVECTOR
      DO ICOL = 1, N 
         JBGN = JA(ICOL) + 1 
         JEND = JA(ICOL+1) - 1 
         IF (JBGN <= JEND) THEN 
            DO J = JBGN, JEND 
               IF (IA(J) < ICOL) THEN 
                  NCOL(ICOL) = NCOL(ICOL) + 1 
               ELSE 
                  NROW(IA(J)) = NROW(IA(J)) + 1 
                  IF (ISYM /= 0) NCOL(IA(J)) = NCOL(IA(J)) + 1 
               ENDIF 
            END DO 
         ENDIF 
      END DO 
      JU(1) = 1 
      IL(1) = 1 
      DO ICOL = 1, N 
         IL(ICOL+1) = IL(ICOL) + NROW(ICOL) 
         JU(ICOL+1) = JU(ICOL) + NCOL(ICOL) 
         NROW(ICOL) = IL(ICOL) 
         NCOL(ICOL) = JU(ICOL) 
      END DO 
      DO ICOL = 1, N 
         DINV(ICOL) = A(JA(ICOL)) 
         JBGN = JA(ICOL) + 1 
         JEND = JA(ICOL+1) - 1 
         IF (JBGN <= JEND) THEN 
            DO J = JBGN, JEND 
               IROW = IA(J) 
               IF (IROW < ICOL) THEN 
!         Part of the upper triangle.
                  IU(NCOL(ICOL)) = IROW 
                  U(NCOL(ICOL)) = A(J) 
                  NCOL(ICOL) = NCOL(ICOL) + 1 
               ELSE 
!         Part of the lower triangle (stored by row).
                  JL(NROW(IROW)) = ICOL 
                  L(NROW(IROW)) = A(J) 
                  NROW(IROW) = NROW(IROW) + 1 
                  IF (ISYM /= 0) THEN 
!         Symmetric...Copy lower triangle into upper triangle as well.
                     IU(NCOL(IROW)) = ICOL 
                     U(NCOL(IROW)) = A(J) 
                     NCOL(IROW) = NCOL(IROW) + 1 
                  ENDIF 
               ENDIF 
            END DO 
         ENDIF 
      END DO 
      DO K = 2, N 
         JBGN = JU(K) 
         JEND = JU(K+1) - 1 
         IF (JBGN < JEND) THEN 
            DO J = JBGN, JEND - 1 
               DO I = J + 1, JEND 
                  IF (IU(J) > IU(I)) THEN 
                     ITEMP = IU(J) 
                     IU(J) = IU(I) 
                     IU(I) = ITEMP 
                     TEMP = U(J) 
                     U(J) = U(I) 
                     U(I) = TEMP 
                  ENDIF 
               END DO 
            END DO 
         ENDIF 
         IBGN = IL(K) 
         IEND = IL(K+1) - 1 
         IF (IBGN < IEND) THEN 
            DO I = IBGN, IEND - 1 
               DO J = I + 1, IEND 
                  IF (JL(I) > JL(J)) THEN 
                     JTEMP = JU(I) 
                     JU(I) = JU(J) 
                     JU(J) = JTEMP 
                     TEMP = L(I) 
                     L(I) = L(J) 
                     L(J) = TEMP 
                  ENDIF 
               END DO 
            END DO 
         ENDIF 
      END DO 
      DO I = 2, N 
!
!           I-th row of L
         INDX1 = IL(I) 
         INDX2 = IL(I+1) - 1 
         IF (INDX1 <= INDX2) THEN 
            DO INDX = INDX1, INDX2 
               IF (INDX /= INDX1) THEN 
                  INDXR1 = INDX1 
                  INDXR2 = INDX - 1 
                  INDXC1 = JU(JL(INDX)) 
                  INDXC2 = JU(JL(INDX)+1) - 1 
                  IF (INDXC1 <= INDXC2) THEN 
  160                CONTINUE 
                     KR = JL(INDXR1) 
  170                CONTINUE 
                     KC = IU(INDXC1) 
                     IF(KR > KC) THEN
                        INDXC1 = INDXC1 + 1 
                        IF (INDXC1 <= INDXC2) GO TO 170 
                     ELSEIF(KR < KC) THEN
                        INDXR1 = INDXR1 + 1 
                        IF (INDXR1 <= INDXR2) GO TO 160 
                     ELSEIF(KR == KC) THEN
                        L(INDX) = L(INDX) - L(INDXR1)*DINV(KC)*U(INDXC1) 
                        INDXR1 = INDXR1 + 1 
                        INDXC1 = INDXC1 + 1 
                        IF (INDXR1<=INDXR2 .AND. INDXC1<=INDXC2) GO TO 160 
                     ENDIF 
                  ENDIF 
               ENDIF 
               L(INDX) = L(INDX)/DINV(JL(INDX)) 
            END DO 
         ENDIF 
         INDX1 = JU(I) 
         INDX2 = JU(I+1) - 1 
         IF (INDX1 <= INDX2) THEN 
            DO INDX = INDX1, INDX2 
               IF (INDX /= INDX1) THEN 
                  INDXC1 = INDX1 
                  INDXC2 = INDX - 1 
                  INDXR1 = IL(IU(INDX)) 
                  INDXR2 = IL(IU(INDX)+1) - 1 
                  IF (INDXR1 <= INDXR2) THEN 
  210                CONTINUE 
                     KR = JL(INDXR1) 
  220                CONTINUE 
                     KC = IU(INDXC1) 
                     IF(KR > KC) THEN
                        INDXC1 = INDXC1 + 1 
                        IF (INDXC1 <= INDXC2) GO TO 220 
                     ELSEIF(KR < KC) THEN
                        INDXR1 = INDXR1 + 1 
                        IF (INDXR1 <= INDXR2) GO TO 210 
                     ELSEIF(KR == KC) THEN
                        U(INDX) = U(INDX) - L(INDXR1)*DINV(KC)*U(INDXC1) 
                        INDXR1 = INDXR1 + 1 
                        INDXC1 = INDXC1 + 1 
                        IF (INDXR1<=INDXR2 .AND. INDXC1<=INDXC2) GO TO 210 
                     ENDIF 
                  ENDIF 
               ENDIF 
               U(INDX) = U(INDX)/DINV(IU(INDX)) 
            END DO 
         ENDIF 
         INDXR1 = IL(I) 
         INDXR2 = IL(I+1) - 1 
         IF (INDXR1 <= INDXR2) THEN 
            INDXC1 = JU(I) 
            INDXC2 = JU(I+1) - 1 
            IF (INDXC1 <= INDXC2) THEN 
  270          CONTINUE 
               KR = JL(INDXR1) 
  280          CONTINUE 
               KC = IU(INDXC1) 
               IF(KR > KC) THEN
                  INDXC1 = INDXC1 + 1 
                  IF (INDXC1 <= INDXC2) GO TO 280 
               ELSEIF(KR < KC) THEN
                  INDXR1 = INDXR1 + 1 
                  IF (INDXR1 <= INDXR2) GO TO 270 
               ELSEIF(KR == KC) THEN
                  DINV(I) = DINV(I) - L(INDXR1)*DINV(KC)*U(INDXC1) 
                  INDXR1 = INDXR1 + 1 
                  INDXC1 = INDXC1 + 1 
                  IF (INDXR1<=INDXR2 .AND. INDXC1<=INDXC2) GO TO 270 
               ENDIF
!
            ENDIF 
         ENDIF 
      END DO 
      I = 1 
      IF (N > 0) THEN 
         DINV = 1./DINV 
         I = N + 1 
      ENDIF 
      RETURN  
!------------- LAST LINE OF DSILUS FOLLOWS ----------------------------
      END SUBROUTINE DSILUS 
