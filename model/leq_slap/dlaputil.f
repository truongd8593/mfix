!DECK DBHIN
      SUBROUTINE DBHIN(N, NELT, IA, JA, A, ISYM, SOLN, RHS, IUNIT, JOB) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  15:55:30  12/05/98  
!...Switches: -xf
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER N, NELT, ISYM, IUNIT, JOB 
      INTEGER, DIMENSION(NELT) :: IA, JA 
      DOUBLE PRECISION, DIMENSION(NELT) :: A 
      DOUBLE PRECISION, DIMENSION(N) :: SOLN, RHS 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NLINE, NPLS, NRILS, NNVLS, NRHSLS, NROW, NCOL, NIND, NELE, &
         JOBRET, I, ICOL, IBGN, IEND, ITEMP, J 
      DOUBLE PRECISION :: TEMP 
      CHARACTER :: TITLE*80, CODE*3, PNTFMT*16, RINFMT*16, NVLFMT*20, RHSFMT*20 
!-----------------------------------------------
!
!         Local Variables
!
!
!
!         Read Matrices In BOEING-HARWELL format.
!
! NLINE  Number of Data (after the header) lines in the file.
! NPLS   Number of lines for the Column Pointer data in the file.
! NRILS  Number of lines for the Row indicies in the data file.
! NNVLS  Number of lines for the Matrix elements in the data file.
! NRHSLS Number of lines for the RHS in the data file.
!
!***FIRST EXECUTABLE STATEMENT  DBHIN
      READ (IUNIT, 9000) TITLE 
      READ (IUNIT, 9010) NLINE, NPLS, NRILS, NNVLS, NRHSLS 
      READ (IUNIT, 9020) CODE, NROW, NCOL, NIND, NELE 
      READ (IUNIT, 9030) PNTFMT, RINFMT, NVLFMT, RHSFMT 
!
      IF (NROW > N) THEN 
         N = NROW 
         JOBRET = -1 
         GO TO 999 
      ENDIF 
      IF (NIND > NELT) THEN 
         NELT = NIND 
         JOBRET = -2 
         GO TO 999 
      ENDIF 
!
!         Set the parameters.
!
      N = NROW 
      NELT = NIND 
      IF (CODE == 'RUA') THEN 
         ISYM = 0 
      ELSE IF (CODE == 'RSA') THEN 
         ISYM = 1 
      ELSE 
         JOBRET = -3 
         GO TO 999 
      ENDIF 
      READ (IUNIT, PNTFMT) (JA(I),I=1,N + 1) 
      READ (IUNIT, RINFMT) (IA(I),I=1,NELT) 
      JOBRET = 10 
      IF (NNVLS > 0) THEN 
         READ (IUNIT, NVLFMT) (A(I),I=1,NELT) 
         JOBRET = 0 
      ENDIF 
      IF (NRHSLS>0 .AND. MOD(JOB,2)==1) THEN 
         READ (5, RHSFMT) (RHS(I),I=1,N) 
         JOBRET = JOBRET + 1 
      ENDIF 
!
!         Now loop thru the IA(i) array making sure that the Diagonal
!         matrix element appears first in the column.  Then sort the
!         rest of the column in ascending order.
!
!VD$R NOCONCUR
!VD$R NOVECTOR
      DO ICOL = 1, N 
         IBGN = JA(ICOL) 
         IEND = JA(ICOL+1) - 1 
         DO I = IBGN, IEND 
            IF (IA(I) == ICOL) THEN 
!         Swap the diag element with the first element in the column.
               ITEMP = IA(I) 
               IA(I) = IA(IBGN) 
               IA(IBGN) = ITEMP 
               TEMP = A(I) 
               A(I) = A(IBGN) 
               A(IBGN) = TEMP 
               EXIT  
            ENDIF 
         END DO 
         IBGN = IBGN + 1 
         IF (IBGN < IEND) THEN 
            DO I = IBGN, IEND 
               DO J = I + 1, IEND 
                  IF (IA(I) > IA(J)) THEN 
                     ITEMP = IA(I) 
                     IA(I) = IA(J) 
                     IA(J) = ITEMP 
                     TEMP = A(I) 
                     A(I) = A(J) 
                     A(J) = TEMP 
                  ENDIF 
               END DO 
            END DO 
         ENDIF 
      END DO 
  999 CONTINUE 
      JOB = JOBRET 
      RETURN  
 9000 FORMAT(A80) 
 9010 FORMAT(5I14) 
 9020 FORMAT(A3,11X,4I14) 
 9030 FORMAT(2A16,2A20) 
!------------- LAST LINE OF DBHIN FOLLOWS ------------------------------
      END SUBROUTINE DBHIN 
!DECK DCHKW
      SUBROUTINE DCHKW(NAME, LOCIW, LENIW, LOCW, LENW, IERR, ITER, ERR) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  15:55:30  12/05/98  
!...Switches: -xf
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER LOCIW, LENIW, LOCW, LENW, IERR, ITER 
      DOUBLE PRECISION ERR 
      CHARACTER NAME*(*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      CHARACTER :: MESG*72 
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      DOUBLE PRECISION , EXTERNAL :: D1MACH 
      EXTERNAL XERRWV 
!-----------------------------------------------
!
!         Check the Integer workspace situation.
!***FIRST EXECUTABLE STATEMENT  DCHKW
      IERR = 0 
      IF (LOCIW > LENIW) THEN 
         IERR = 1 
         ITER = 0 
         ERR = D1MACH(2) 
         MESG = NAME//': INTEGER work array too short. '//&
            ' IWORK needs i1: have allocated i2.' 
         CALL XERRWV (MESG, LEN(MESG), 1, 1, 2, LOCIW, LENIW, 0, 0.0, 0.0) 
      ENDIF 
!
!         Check the Double Precision workspace situation.
      IF (LOCW > LENW) THEN 
         IERR = 1 
         ITER = 0 
         ERR = D1MACH(2) 
         MESG = NAME//': DOUBLE PRECISION work array too short. '//&
            ' RWORK needs i1: have allocated i2.' 
         CALL XERRWV (MESG, LEN(MESG), 1, 1, 2, LOCW, LENW, 0, 0.0, 0.0) 
      ENDIF 
      RETURN  
!------------- LAST LINE OF DCHKW FOLLOWS ----------------------------
      END SUBROUTINE DCHKW 
!DECK QS2I1D
      SUBROUTINE QS2I1D(IA, JA, A, N, KFLAG) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  15:55:30  12/05/98  
!...Switches: -xf
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER N, KFLAG 
      INTEGER, DIMENSION(N) :: IA, JA 
      DOUBLE PRECISION, DIMENSION(N) :: A 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, DIMENSION(21) :: IL, IU 
      INTEGER :: IT, IIT, JT, JJT, NN, KK, I, M, J, K, IJ, L 
      DOUBLE PRECISION :: TA, TTA, R 
!-----------------------------------------------
!VD$R NOVECTOR
!VD$R NOCONCUR
!
!***FIRST EXECUTABLE STATEMENT  QS2I1D
      NN = N 
      IF (NN < 1) THEN 
         CALL XERROR (&
            'QS2I1D- the number of values to be sorted was noT POSITIVE.', 59, &
            1, 1) 
         RETURN  
      ENDIF 
      IF (N == 1) RETURN  
      KK = IABS(KFLAG) 
      IF (KK /= 1) THEN 
         CALL XERROR ('QS2I1D- the sort control parameter, k, was not 1 OR -1.'&
            , 55, 2, 1) 
         RETURN  
      ENDIF 
!
!     Alter array IA to get decreasing order if needed.
!
      IF (KFLAG < 1) THEN 
         DO I = 1, NN 
            IA(I) = -IA(I) 
         END DO 
      ENDIF 
!
!     Sort IA and carry JA and A along.
!     And now...Just a little black magic...
      M = 1 
      I = 1 
      J = NN 
      R = .375 
  210 CONTINUE 
      IF (R <= 0.5898437) THEN 
         R = R + 3.90625E-2 
      ELSE 
         R = R - .21875 
      ENDIF 
  225 CONTINUE 
      K = I 
!
!     Select a central element of the array and save it in location
!     it, jt, at.
!
      IJ = I + IDINT(DBLE(J - I)*R) 
      IT = IA(IJ) 
      JT = JA(IJ) 
      TA = A(IJ) 
!
!     If first element of array is greater than it, interchange with it.
!
      IF (IA(I) > IT) THEN 
         IA(IJ) = IA(I) 
         IA(I) = IT 
         IT = IA(IJ) 
         JA(IJ) = JA(I) 
         JA(I) = JT 
         JT = JA(IJ) 
         A(IJ) = A(I) 
         A(I) = TA 
         TA = A(IJ) 
      ENDIF 
      L = J 
!
!     If last element of array is less than it, swap with it.
!
      IF (IA(J) < IT) THEN 
         IA(IJ) = IA(J) 
         IA(J) = IT 
         IT = IA(IJ) 
         JA(IJ) = JA(J) 
         JA(J) = JT 
         JT = JA(IJ) 
         A(IJ) = A(J) 
         A(J) = TA 
         TA = A(IJ) 
!
!     If first element of array is greater than it, swap with it.
!
         IF (IA(I) > IT) THEN 
            IA(IJ) = IA(I) 
            IA(I) = IT 
            IT = IA(IJ) 
            JA(IJ) = JA(I) 
            JA(I) = JT 
            JT = JA(IJ) 
            A(IJ) = A(I) 
            A(I) = TA 
            TA = A(IJ) 
         ENDIF 
      ENDIF 
!
!     Find an element in the second half of the array which is
!     smaller than it.
!
  240 CONTINUE 
      L = L - 1 
 1003 CONTINUE 
      IF (IA(L) <= IT) GO TO 1002 
      L = L - 1 
      IF (IA(L) <= IT) GO TO 1002 
      L = L - 1 
      IF (IA(L) <= IT) GO TO 1002 
      L = L - 1 
      IF (IA(L) <= IT) GO TO 1002 
      L = L - 1 
      IF (IA(L) <= IT) GO TO 1002 
      L = L - 1 
      IF (IA(L) <= IT) GO TO 1002 
      L = L - 1 
      IF (IA(L) <= IT) GO TO 1002 
      L = L - 1 
      IF (IA(L) <= IT) GO TO 1002 
      L = L - 1 
      GO TO 1003 
 1002 CONTINUE 
!
!     Find an element in the first half of the array which is
!     greater than it.
!
      K = K + 1 
 1005 CONTINUE 
      IF (IA(K) >= IT) GO TO 1004 
      K = K + 1 
      IF (IA(K) >= IT) GO TO 1004 
      K = K + 1 
      IF (IA(K) >= IT) GO TO 1004 
      K = K + 1 
      IF (IA(K) >= IT) GO TO 1004 
      K = K + 1 
      IF (IA(K) >= IT) GO TO 1004 
      K = K + 1 
      IF (IA(K) >= IT) GO TO 1004 
      K = K + 1 
      IF (IA(K) >= IT) GO TO 1004 
      K = K + 1 
      IF (IA(K) >= IT) GO TO 1004 
      K = K + 1 
      GO TO 1005 
 1004 CONTINUE 
!
!     Interchange these elements.
!
      IF (K <= L) THEN 
         IIT = IA(L) 
         IA(L) = IA(K) 
         IA(K) = IIT 
         JJT = JA(L) 
         JA(L) = JA(K) 
         JA(K) = JJT 
         TTA = A(L) 
         A(L) = A(K) 
         A(K) = TTA 
         GO TO 240 
      ENDIF 
!
!     Save upper and lower subscripts of the array yet to be sorted.
!
      IF (L - I > J - K) THEN 
         IL(M) = I 
         IU(M) = L 
         I = K 
         M = M + 1 
      ELSE 
         IL(M) = K 
         IU(M) = J 
         J = L 
         M = M + 1 
      ENDIF 
      GO TO 260 
!
!     Begin again on another portion of the unsorted array.
!
  255 CONTINUE 
      M = M - 1 
      IF (M == 0) GO TO 300 
      I = IL(M) 
      J = IU(M) 
  260 CONTINUE 
      IF (J - I >= 1) GO TO 225 
      IF (I == J) GO TO 255 
      IF (I == 1) GO TO 210 
      I = I - 1 
  265 CONTINUE 
      I = I + 1 
      IF (I == J) GO TO 255 
      IT = IA(I+1) 
      JT = JA(I+1) 
      TA = A(I+1) 
      DO WHILE(IA(I) <= IT) 
         I = I + 1 
         IF (I == J) GO TO 255 
         IT = IA(I+1) 
         JT = JA(I+1) 
         TA = A(I+1) 
      END DO 
      K = I 
      IA(K+1) = IA(K) 
      JA(K+1) = JA(K) 
      A(K+1) = A(K) 
      K = K - 1 
      DO WHILE(IT < IA(K)) 
         IA(K+1) = IA(K) 
         JA(K+1) = JA(K) 
         A(K+1) = A(K) 
         K = K - 1 
      END DO 
      IA(K+1) = IT 
      JA(K+1) = JT 
      A(K+1) = TA 
      GO TO 265 
!
!     Clean up, if necessary.
!
  300 CONTINUE 
      IF (KFLAG < 1) THEN 
         DO I = 1, NN 
            IA(I) = -IA(I) 
         END DO 
      ENDIF 
      RETURN  
!------------- LAST LINE OF QS2I1D FOLLOWS ----------------------------
      END SUBROUTINE QS2I1D 
!DECK DS2Y
      SUBROUTINE DS2Y(N, NELT, IA, JA, A, ISYM) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  15:55:30  12/05/98  
!...Switches: -xf
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER N, NELT, ISYM 
      INTEGER, DIMENSION(NELT) :: IA, JA 
      DOUBLE PRECISION, DIMENSION(NELT) :: A 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: ICOL, J, IBGN, IEND, I, ITEMP 
      DOUBLE PRECISION :: TEMP 
!-----------------------------------------------
!
!         Check to see if the (IA,JA,A) arrays are in SLAP Column
!         format.  If it's not then transform from SLAP Triad.
!***FIRST EXECUTABLE STATEMENT  DS2LT
      IF (JA(N+1) == NELT + 1) RETURN  
!
!         Sort into ascending order by COLUMN (on the ja array).
!         This will line up the columns.
!
      CALL QS2I1D (JA, IA, A, NELT, 1) 
!
!         Loop over each column to see where the column indicies change
!         in the column index array ja.  This marks the beginning of the
!         next column.
!
!VD$R NOVECTOR
      JA(1) = 1 
      L20: DO ICOL = 1, N - 1 
         DO J = JA(ICOL) + 1, NELT 
            IF (JA(J) /= ICOL) THEN 
               JA(ICOL+1) = J 
               CYCLE  L20 
            ENDIF 
         END DO 
      END DO L20 
      JA(N+1) = NELT + 1 
!
!         Mark the n+2 element so that future calls to a SLAP routine
!         utilizing the YSMP-Column storage format will be able to tell.
!
      JA(N+2) = 0 
!
!         Now loop thru the ia(i) array making sure that the Diagonal
!         matrix element appears first in the column.  Then sort the
!         rest of the column in ascending order.
!
      DO ICOL = 1, N 
         IBGN = JA(ICOL) 
         IEND = JA(ICOL+1) - 1 
         DO I = IBGN, IEND 
            IF (IA(I) == ICOL) THEN 
!         Swap the diag element with the first element in the column.
               ITEMP = IA(I) 
               IA(I) = IA(IBGN) 
               IA(IBGN) = ITEMP 
               TEMP = A(I) 
               A(I) = A(IBGN) 
               A(IBGN) = TEMP 
               EXIT  
            ENDIF 
         END DO 
         IBGN = IBGN + 1 
         IF (IBGN < IEND) THEN 
            DO I = IBGN, IEND 
               DO J = I + 1, IEND 
                  IF (IA(I) > IA(J)) THEN 
                     ITEMP = IA(I) 
                     IA(I) = IA(J) 
                     IA(J) = ITEMP 
                     TEMP = A(I) 
                     A(I) = A(J) 
                     A(J) = TEMP 
                  ENDIF 
               END DO 
            END DO 
         ENDIF 
      END DO 
      RETURN  
!------------- LAST LINE OF DS2Y FOLLOWS ----------------------------
      END SUBROUTINE DS2Y 
!DECK DCPPLT
      SUBROUTINE DCPPLT(N, NELT, IA, JA, A, ISYM, IUNIT) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  15:55:30  12/05/98  
!...Switches: -xf
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER N, NELT, ISYM, IUNIT 
      INTEGER, DIMENSION(NELT) :: IA, JA 
      DOUBLE PRECISION, DIMENSION(NELT) :: A 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NMAX, I, ICOL, JBGN, JEND, J, IROW 
      CHARACTER, DIMENSION(225) :: CHMAT*225 
!-----------------------------------------------
!
!         Set up the character matrix...
!***FIRST EXECUTABLE STATEMENT  DCPPLT
      NMAX = MIN(225,N) 
      DO I = 1, NMAX 
         CHMAT(I)(1:NMAX) = ' ' 
      END DO 
      DO ICOL = 1, NMAX 
         JBGN = JA(ICOL) 
         JEND = JA(ICOL+1) - 1 
         DO J = JBGN, JEND 
            IROW = IA(J) 
            IF (IROW <= NMAX) THEN 
               IF (ISYM /= 0) THEN 
!         Put in non-dym part as well...
                  IF (A(J) == 0.0D0) THEN 
                     CHMAT(IROW)(ICOL:ICOL) = '0' 
                  ELSE IF (A(J) > 0.0D0) THEN 
                     CHMAT(IROW)(ICOL:ICOL) = '#' 
                  ELSE 
                     CHMAT(IROW)(ICOL:ICOL) = '*' 
                  ENDIF 
               ENDIF 
               IF (IROW == ICOL) THEN 
!         Diagonal entry.
                  IF (A(J) == 0.0D0) THEN 
                     CHMAT(IROW)(ICOL:ICOL) = '0' 
                  ELSE IF (A(J) > 0.0D0) THEN 
                     CHMAT(IROW)(ICOL:ICOL) = 'D' 
                  ELSE 
                     CHMAT(IROW)(ICOL:ICOL) = 'N' 
                  ENDIF 
               ELSE 
!         Off-Diagonal entry
                  IF (A(J) == 0.0D0) THEN 
                     CHMAT(IROW)(ICOL:ICOL) = '0' 
                  ELSE IF (A(J) > 0.0D0) THEN 
                     CHMAT(IROW)(ICOL:ICOL) = '#' 
                  ELSE 
                     CHMAT(IROW)(ICOL:ICOL) = '*' 
                  ENDIF 
               ENDIF 
            ENDIF 
         END DO 
      END DO 
      WRITE (IUNIT, 1000) N, NELT, FLOAT(NELT)/FLOAT(N*N) 
      WRITE (IUNIT, 1010) (MOD(I,10),I=1,NMAX) 
!
!         Write out the character representations matrix elements.
      DO IROW = 1, NMAX 
         WRITE (IUNIT, 1020) IROW, CHMAT(IROW)(1:NMAX) 
      END DO 
      RETURN  
 1000 FORMAT(/'**** Picture of Column SLAP matrix follows ****'/&
         ' N, NELT and Density = ',2I10,E16.7) 
 1010 FORMAT(4X,255(I1)) 
 1020 FORMAT(1X,I3,A) 
!------------- LAST LINE OF DCPPLT FOLLOWS ----------------------------
      END SUBROUTINE DCPPLT 
!DECK DTOUT
      SUBROUTINE DTOUT(N, NELT, IA, JA, A, ISYM, SOLN, RHS, IUNIT, JOB) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  15:55:30  12/05/98  
!...Switches: -xf
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER N, NELT, ISYM, IUNIT, JOB 
      INTEGER, DIMENSION(NELT) :: IA, JA 
      DOUBLE PRECISION, DIMENSION(NELT) :: A 
      DOUBLE PRECISION, DIMENSION(N) :: SOLN, RHS 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IRHS, ISOLN, I 
!-----------------------------------------------
!
!         Local variables.
!
!
!         If RHS and SOLN are to be printed also.
!         Write out the information heading.
!***FIRST EXECUTABLE STATEMENT  DTOUT
      IRHS = 0 
      ISOLN = 0 
      IF (JOB==1 .OR. JOB==3) IRHS = 1 
      IF (JOB > 1) ISOLN = 1 
      WRITE (IUNIT, 1000) N, NELT, ISYM, IRHS, ISOLN 
!
!         Write out the matrix non-zeros in Triad format.
      DO I = 1, NELT 
         WRITE (IUNIT, 1010) IA(I), JA(I), A(I) 
      END DO 
      IF (IRHS == 1) WRITE (IUNIT, 1020) (RHS(I),I=1,N) 
      IF (ISOLN == 1) WRITE (IUNIT, 1020) (SOLN(I),I=1,N) 
      RETURN  
 1000 FORMAT(5I10) 
 1010 FORMAT(1X,I5,1X,I5,1X,E16.7) 
 1020 FORMAT(1X,E16.7) 
!------------- LAST LINE OF DTOUT FOLLOWS ----------------------------
      END SUBROUTINE DTOUT 
!DECK DTIN
      SUBROUTINE DTIN(N, NELT, IA, JA, A, ISYM, SOLN, RHS, IUNIT, JOB) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  15:55:30  12/05/98  
!...Switches: -xf
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER N, NELT, ISYM, IUNIT, JOB 
      INTEGER, DIMENSION(NELT) :: IA, JA 
      DOUBLE PRECISION, DIMENSION(NELT) :: A 
      DOUBLE PRECISION, DIMENSION(N) :: SOLN, RHS 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IRHS, ISOLN, I, NELTMAX, JOBRET 
!-----------------------------------------------
!
!         Local variables.
!
!
!         Read in the information heading.
!***FIRST EXECUTABLE STATEMENT  DTIN
      NELTMAX = NELT 
      READ (IUNIT, 1000) N, NELT, ISYM, IRHS, ISOLN 
      NELT = MIN(NELT,NELTMAX) 
!
!         Read in the matrix non-zeros in Triad format.
      DO I = 1, NELT 
         READ (IUNIT, 1010) IA(I), JA(I), A(I) 
      END DO 
      JOBRET = 0 
      IF (JOB==1 .OR. JOB==3) THEN 
!
!         Check to see if rhs is in the file.
         IF (IRHS == 1) THEN 
            JOBRET = 1 
            READ (IUNIT, 1020) (RHS(I),I=1,N) 
         ELSE 
            I = 1 
            IF (N > 0) THEN 
               RHS(:N) = 0.0D0 
               I = N + 1 
            ENDIF 
         ENDIF 
      ENDIF 
!
!         If requested, read in the soln.
      IF (JOB > 1) THEN 
!
!         Check to see if soln is in the file.
         IF (ISOLN == 1) THEN 
            JOBRET = JOBRET + 2 
            READ (IUNIT, 1020) (SOLN(I),I=1,N) 
         ELSE 
            I = 1 
            IF (N > 0) THEN 
               SOLN(:N) = 0.0D0 
               I = N + 1 
            ENDIF 
         ENDIF 
      ENDIF 
!
      JOB = JOBRET 
      RETURN  
 1000 FORMAT(5I10) 
 1010 FORMAT(1X,I5,1X,I5,1X,E16.7) 
 1020 FORMAT(1X,E16.7) 
!------------- LAST LINE OF DTIN FOLLOWS ----------------------------
      END SUBROUTINE DTIN 
