      REAL FUNCTION R1MACH (I) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  15:55:30  12/05/98  
!...Switches: -xf
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER I 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, DIMENSION(2) :: SMALL, LARGE, RIGHT, DIVER, LOG10 
      REAL, DIMENSION(5) :: RMACH 
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      INTEGER , EXTERNAL :: I1MACH 
!-----------------------------------------------
!
!  SINGLE-PRECISION MACHINE CONSTANTS
!
!  R1MACH(1) = B**(EMIN-1), THE SMALLEST POSITIVE MAGNITUDE.
!
!  R1MACH(2) = B**EMAX*(1 - B**(-T)), THE LARGEST MAGNITUDE.
!
!  R1MACH(3) = B**(-T), THE SMALLEST RELATIVE SPACING.
!
!  R1MACH(4) = B**(1-T), THE LARGEST RELATIVE SPACING.
!
!  R1MACH(5) = LOG10(B)
!
!  TO ALTER THIS FUNCTION FOR A PARTICULAR ENVIRONMENT,
!  THE DESIRED SET OF DATA STATEMENTS SHOULD BE ACTIVATED BY
!  REMOVING THE C FROM COLUMN 1.
!  ON RARE MACHINES A STATIC STATEMENT MAY NEED TO BE ADDED.
!  (BUT PROBABLY MORE SYSTEMS PROHIBIT IT THAN REQUIRE IT.)
!
!  WHERE POSSIBLE, OCTAL OR HEXADECIMAL CONSTANTS HAVE BEEN USED
!  TO SPECIFY THE CONSTANTS EXACTLY WHICH HAS IN SOME CASES
!  REQUIRED THE USE OF EQUIVALENT INTEGER ARRAYS.
!
!
!
      EQUIVALENCE (RMACH(1), SMALL(1)) 
      EQUIVALENCE (RMACH(2), LARGE(1)) 
      EQUIVALENCE (RMACH(3), RIGHT(1)) 
      EQUIVALENCE (RMACH(4), DIVER(1)) 
      EQUIVALENCE (RMACH(5), LOG10(1)) 
      DATA RMACH(1)/1.17549435E-38/ 
!      DATA RMACH(2) / 3.40282347E+38 /
      DATA RMACH(2)/3.402823E+38/ 
      DATA RMACH(3)/5.96016605E-08/ 
      DATA RMACH(4)/1.19203321E-07/ 
      DATA RMACH(5)/3.01030010E-01/ 
!
!     MACHINE CONSTANTS FOR THE Alliant FX/8 UNIX Fortran COMPILER
!     WITH THE -r8 COMMAND LINE OPTION.  This option causes all variables
!     declared with 'real' to be of type 'real*8' or double precision.
!     This option does not override the 'real*4' declarations.  These
!     R1MACH numbers below and the coresponding I1MACH are simply the double
!     precision or 'real*8' numbers.  If you use the -r8 your whole code
!     (and the user libraries you link with, the system libraries are taken
!     care of automagicly) must be compiled with this option.
!
!$$$      DATA RMACH(1) / 2.22507385850721D-308 /
!$$$      DATA RMACH(2) / 1.79769313486231D+308 /
!$$$      DATA RMACH(3) / 1.1101827117665D-16 /
!$$$      DATA RMACH(4) / 2.2203654423533D-16 /
!$$$      DATA RMACH(5) / 3.01029995663981E-1 /
!
!     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES.
!
!      DATA RMACH(1) / O000400000000 /
!      DATA RMACH(2) / O377777777777 /
!      DATA RMACH(3) / O146400000000 /
!      DATA RMACH(4) / O147400000000 /
!      DATA RMACH(5) / O177464202324 /
!
!     MACHINE CONSTANTS FOR THE VAX UNIX F77 COMPILER
!
!      DATA SMALL(1) /       128 /
!      DATA LARGE(1) /    -32769 /
!      DATA RIGHT(1) /     13440 /
!      DATA DIVER(1) /     13568 /
!      DATA LOG10(1) / 547045274 /
!
!     MACHINE CONSTANTS FOR THE VAX-11 WITH
!     FORTRAN IV-PLUS COMPILER
!
!      DATA RMACH(1) / Z00000080 /
!      DATA RMACH(2) / ZFFFF7FFF /
!      DATA RMACH(3) / Z00003480 /
!      DATA RMACH(4) / Z00003500 /
!      DATA RMACH(5) / Z209B3F9A /
!
!     MACHINE CONSTANTS FOR VAX/VMS VERSION 2.2
!
!      DATA RMACH(1) /       '80'X /
!      DATA RMACH(2) / 'FFFF7FFF'X /
!      DATA RMACH(3) /     '3480'X /
!      DATA RMACH(4) /     '3500'X /
!      DATA RMACH(5) / '209B3F9A'X /
!
!     MACHINE CONSTANTS FOR THE SEQUENT BALANCE 8000 AND SVS FORTRAN ON
!     THE AT&T 7300 (UNIX PC)
!
!      DATA SMALL(1) / $00800000 /
!      DATA LARGE(1) / $7F7FFFFF /
!      DATA RIGHT(1) / $33800000 /
!      DATA DIVER(1) / $34000000 /
!      DATA LOG10(1) / $3E9A209B /
!
!     MACHINE CONSTANTS FOR RM FORTRAN (ON THE AT&T 7300)
!
!      DATA RMACH(1) / Z'00800000' /
!      DATA RMACH(2) / Z'7F7FFFFF' /
!      DATA RMACH(3) / Z'33800000' /
!      DATA RMACH(4) / Z'34000000' /
!      DATA RMACH(5) / Z'3E9A209B' /
!
!
      IF (I>=1 .AND. I<=5) THEN 
         R1MACH = RMACH(I) 
         RETURN  
      ENDIF 
      WRITE (I1MACH(2), 1999) I 
 1999 FORMAT(' R1MACH - I OUT OF BOUNDS',I10) 
      STOP  
      END FUNCTION R1MACH 
      DOUBLE PRECISION FUNCTION D1MACH (I) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  15:55:30  12/05/98  
!...Switches: -xf
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER I 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, DIMENSION(4) :: SMALL, LARGE, RIGHT, DIVER, LOG10 
      DOUBLE PRECISION, DIMENSION(5) :: DMACH 
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      INTEGER , EXTERNAL :: I1MACH 
!-----------------------------------------------
!
!  DOUBLE-PRECISION MACHINE CONSTANTS
!
!  D1MACH( 1) = B**(EMIN-1), THE SMALLEST POSITIVE MAGNITUDE.
!
!  D1MACH( 2) = B**EMAX*(1 - B**(-T)), THE LARGEST MAGNITUDE.
!
!  D1MACH( 3) = B**(-T), THE SMALLEST RELATIVE SPACING.
!
!  D1MACH( 4) = B**(1-T), THE LARGEST RELATIVE SPACING.
!
!  D1MACH( 5) = LOG10(B)
!
!  TO ALTER THIS FUNCTION FOR A PARTICULAR ENVIRONMENT,
!  THE DESIRED SET OF DATA STATEMENTS SHOULD BE ACTIVATED BY
!  REMOVING THE C FROM COLUMN 1.
!  ON RARE MACHINES A STATIC STATEMENT MAY NEED TO BE ADDED.
!  (BUT PROBABLY MORE SYSTEMS PROHIBIT IT THAN REQUIRE IT.)
!
!  WHERE POSSIBLE, OCTAL OR HEXADECIMAL CONSTANTS HAVE BEEN USED
!  TO SPECIFY THE CONSTANTS EXACTLY WHICH HAS IN SOME CASES
!  REQUIRED THE USE OF EQUIVALENT INTEGER ARRAYS.
!
!
!
      EQUIVALENCE (DMACH(1), SMALL(1)) 
      EQUIVALENCE (DMACH(2), LARGE(1)) 
      EQUIVALENCE (DMACH(3), RIGHT(1)) 
      EQUIVALENCE (DMACH(4), DIVER(1)) 
      EQUIVALENCE (DMACH(5), LOG10(1)) 
      DATA DMACH(1)/2.22507385850720D-308/ 
      DATA DMACH(2)/1.79769313486231D+308/ 
      DATA DMACH(3)/1.1101827117665D-16/ 
      DATA DMACH(4)/2.2203654423533D-16/ 
      DATA DMACH(5)/3.01029995663981E-1/ 
!
!     MACHINE CONSTANTS FOR THE ALLIANT FX/8 UNIX FORTRAN COMPILER.
!
!$$$      DATA DMACH(1) / 2.22507385850721D-308 /
!$$$      DATA DMACH(2) / 1.79769313486231D+308 /
!$$$      DATA DMACH(3) / 1.1101827117665D-16 /
!$$$      DATA DMACH(4) / 2.2203654423533D-16 /
!$$$      DATA DMACH(5) / 3.01029995663981E-1 /
!
!     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES.
!
!      DATA SMALL(1),SMALL(2) / O000040000000, O000000000000 /
!      DATA LARGE(1),LARGE(2) / O377777777777, O777777777777 /
!      DATA RIGHT(1),RIGHT(2) / O170540000000, O000000000000 /
!      DATA DIVER(1),DIVER(2) / O170640000000, O000000000000 /
!      DATA LOG10(1),LOG10(2) / O177746420232, O411757177572 /
!
!     MACHINE CONSTANTS FOR THE VAX UNIX F77 COMPILER
!
!      DATA SMALL(1),SMALL(2) /        128,           0 /
!      DATA LARGE(1),LARGE(2) /     -32769,          -1 /
!      DATA RIGHT(1),RIGHT(2) /       9344,           0 /
!      DATA DIVER(1),DIVER(2) /       9472,           0 /
!      DATA LOG10(1),LOG10(2) /  546979738,  -805796613 /
!
!     MACHINE CONSTANTS FOR THE VAX-11 WITH
!     FORTRAN IV-PLUS COMPILER
!
!      DATA SMALL(1),SMALL(2) / Z00000080, Z00000000 /
!      DATA LARGE(1),LARGE(2) / ZFFFF7FFF, ZFFFFFFFF /
!      DATA RIGHT(1),RIGHT(2) / Z00002480, Z00000000 /
!      DATA DIVER(1),DIVER(2) / Z00002500, Z00000000 /
!      DATA LOG10(1),LOG10(2) / Z209A3F9A, ZCFF884FB /
!
!     MACHINE CONSTANTS FOR VAX/VMS VERSION 2.2
!
!      DATA SMALL(1),SMALL(2) /       '80'X,        '0'X /
!      DATA LARGE(1),LARGE(2) / 'FFFF7FFF'X, 'FFFFFFFF'X /
!      DATA RIGHT(1),RIGHT(2) /     '2480'X,        '0'X /
!      DATA DIVER(1),DIVER(2) /     '2500'X,        '0'X /
!      DATA LOG10(1),LOG10(2) / '209A3F9A'X, 'CFF884FB'X /
!
!     MACHINE CONSTANTS FOR THE SEQUENT BALANCE 8000
!
!      DATA SMALL(1),SMALL(2) / $00000000,  $00100000 /
!      DATA LARGE(1),LARGE(2) / $FFFFFFFF,  $7FEFFFFF /
!      DATA RIGHT(1),RIGHT(2) / $00000000,  $3CA00000 /
!      DATA DIVER(1),DIVER(2) / $00000000,  $3CB00000 /
!      DATA LOG10(1),LOG10(2) / $509F79FF,  $3FD34413 /
!
!     MACHINE CONSTANTS FOR SVS FORTRAN ON THE AT&T 7300 (UNIX PC)
!
!      DATA SMALL(1),SMALL(2) / $00100000, $00000000 /
!      DATA LARGE(1),LARGE(2) / $7FEFFFFF, $FFFFFFFF /
!      DATA RIGHT(1),RIGHT(2) / $3CA00000, $00000000 /
!      DATA DIVER(1),DIVER(2) / $3CB00000, $00000000 /
!      DATA LOG10(1),LOG10(2) / $3FD34413, $509F79FF /
!
!     MACHINE CONSTANTS FOR THE RM FORTRAN ON THE AT&T 7300 (UNIX PC)
!
!      DATA SMALL(1),SMALL(2) / Z'00100000', Z'00000000' /
!      DATA LARGE(1),LARGE(2) / Z'7FEFFFFF', Z'FFFFFFFF' /
!      DATA RIGHT(1),RIGHT(2) / Z'3CA00000', Z'00000000' /
!      DATA DIVER(1),DIVER(2) / Z'3CB00000', Z'00000000' /
!      DATA LOG10(1),LOG10(2) / Z'3FD34413', Z'509F79FF' /
!
      IF (I>=1 .AND. I<=5) THEN 
         D1MACH = DMACH(I) 
         RETURN  
      ENDIF 
      WRITE (I1MACH(2), 1999) I 
 1999 FORMAT(' D1MACH - I OUT OF BOUNDS',I10) 
      STOP  
      END FUNCTION D1MACH 
      INTEGER FUNCTION I1MACH (I) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  15:55:30  12/05/98  
!...Switches: -xf
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER I 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, DIMENSION(16) :: IMACH 
      INTEGER :: OUTPUT 
!-----------------------------------------------
!
!  I/O UNIT NUMBERS.
!
!    I1MACH( 1) = THE STANDARD INPUT UNIT.
!
!    I1MACH( 2) = THE STANDARD OUTPUT UNIT.
!
!    I1MACH( 3) = THE STANDARD PUNCH UNIT.
!
!    I1MACH( 4) = THE STANDARD ERROR MESSAGE UNIT.
!
!  WORDS.
!
!    I1MACH( 5) = THE NUMBER OF BITS PER INTEGER STORAGE UNIT.
!
!    I1MACH( 6) = THE NUMBER OF CHARACTERS PER INTEGER STORAGE UNIT.
!
!  INTEGERS.
!
!    ASSUME INTEGERS ARE REPRESENTED IN THE S-DIGIT, BASE-A FORM
!
!               SIGN ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )
!
!               WHERE 0 .LE. X(I) .LT. A FOR I=0,...,S-1.
!
!    I1MACH( 7) = A, THE BASE.
!
!    I1MACH( 8) = S, THE NUMBER OF BASE-A DIGITS.
!
!    I1MACH( 9) = A**S - 1, THE LARGEST MAGNITUDE.
!
!  FLOATING-POINT NUMBERS.
!
!    ASSUME FLOATING-POINT NUMBERS ARE REPRESENTED IN THE T-DIGIT,
!    BASE-B FORM
!
!               SIGN (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
!
!               WHERE 0 .LE. X(I) .LT. B FOR I=1,...,T,
!               0 .LT. X(1), AND EMIN .LE. E .LE. EMAX.
!
!    I1MACH(10) = B, THE BASE.
!
!  SINGLE-PRECISION
!
!    I1MACH(11) = T, THE NUMBER OF BASE-B DIGITS.
!
!    I1MACH(12) = EMIN, THE SMALLEST EXPONENT E.
!
!    I1MACH(13) = EMAX, THE LARGEST EXPONENT E.
!
!  DOUBLE-PRECISION
!
!    I1MACH(14) = T, THE NUMBER OF BASE-B DIGITS.
!
!    I1MACH(15) = EMIN, THE SMALLEST EXPONENT E.
!
!    I1MACH(16) = EMAX, THE LARGEST EXPONENT E.
!
!  TO ALTER THIS FUNCTION FOR A PARTICULAR ENVIRONMENT,
!  THE DESIRED SET OF DATA STATEMENTS SHOULD BE ACTIVATED BY
!  REMOVING THE C FROM COLUMN 1.  ALSO, THE VALUES OF
!  I1MACH(1) - I1MACH(4) SHOULD BE CHECKED FOR CONSISTENCY
!  WITH THE LOCAL OPERATING SYSTEM.
!  ON RARE MACHINES A STATIC STATEMENT MAY NEED TO BE ADDED.
!  (BUT PROBABLY MORE SYSTEMS PROHIBIT IT THAN REQUIRE IT.)
!
!
      EQUIVALENCE (IMACH(4), OUTPUT) 
      DATA IMACH(1)/5/ 
      DATA IMACH(2)/6/ 
      DATA IMACH(3)/6/ 
      DATA IMACH(4)/0/ 
      DATA IMACH(5)/32/ 
      DATA IMACH(6)/4/ 
      DATA IMACH(7)/2/ 
      DATA IMACH(8)/32/ 
      DATA IMACH(9)/2147483647/ 
      DATA IMACH(10)/2/ 
      DATA IMACH(11)/24/ 
      DATA IMACH(12)/-126/ 
      DATA IMACH(13)/128/ 
      DATA IMACH(14)/53/ 
      DATA IMACH(15)/-1022/ 
      DATA IMACH(16)/1024/ 
!
!     MACHINE CONSTANTS FOR THE ALLIANT FX/8 UNIX FORTRAN COMPILER.
!
!$$$      DATA IMACH( 1) /     5 /
!$$$      DATA IMACH( 2) /     6 /
!$$$      DATA IMACH( 3) /     6 /
!$$$      DATA IMACH( 4) /     0 /
!$$$      DATA IMACH( 5) /    32 /
!$$$      DATA IMACH( 6) /     4 /
!$$$      DATA IMACH( 7) /     2 /
!$$$      DATA IMACH( 8) /    32 /
!$$$      DATA IMACH( 9) /2147483647/
!$$$      DATA IMACH(10) /     2 /
!$$$      DATA IMACH(11) /    24 /
!$$$      DATA IMACH(12) /  -126 /
!$$$      DATA IMACH(13) /   128 /
!$$$      DATA IMACH(14) /    53 /
!$$$      DATA IMACH(15) / -1022 /
!$$$      DATA IMACH(16) /  1024 /
!
!     MACHINE CONSTANTS FOR THE ALLIANT FX/8 UNIX FORTRAN COMPILER.
!     WITH THE -r8 COMMAND LINE OPTION.
!
!$$$      DATA IMACH( 1) /     5 /
!$$$      DATA IMACH( 2) /     6 /
!$$$      DATA IMACH( 3) /     6 /
!$$$      DATA IMACH( 4) /     0 /
!$$$      DATA IMACH( 5) /    32 /
!$$$      DATA IMACH( 6) /     4 /
!$$$      DATA IMACH( 7) /     2 /
!$$$      DATA IMACH( 8) /    32 /
!$$$      DATA IMACH( 9) /2147483647/
!$$$      DATA IMACH(10) /     2 /
!$$$      DATA IMACH(11) /    53 /
!$$$      DATA IMACH(12) / -1022 /
!$$$      DATA IMACH(13) /  1024 /
!$$$      DATA IMACH(14) /    53 /
!$$$      DATA IMACH(15) / -1022 /
!$$$      DATA IMACH(16) /  1024 /
!
!     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES.
!
!     NOTE THAT THE PUNCH UNIT, I1MACH(3), HAS BEEN SET TO 7
!     WHICH IS APPROPRIATE FOR THE UNIVAC-FOR SYSTEM.
!     IF YOU HAVE THE UNIVAC-FTN SYSTEM, SET IT TO 1.
!
!      DATA IMACH( 1) /    5 /
!      DATA IMACH( 2) /    6 /
!      DATA IMACH( 3) /    7 /
!      DATA IMACH( 4) /    6 /
!      DATA IMACH( 5) /   36 /
!      DATA IMACH( 6) /    6 /
!      DATA IMACH( 7) /    2 /
!      DATA IMACH( 8) /   35 /
!      DATA IMACH( 9) / O377777777777 /
!      DATA IMACH(10) /    2 /
!      DATA IMACH(11) /   27 /
!      DATA IMACH(12) / -128 /
!      DATA IMACH(13) /  127 /
!      DATA IMACH(14) /   60 /
!      DATA IMACH(15) /-1024 /
!      DATA IMACH(16) / 1023 /
!
!     MACHINE CONSTANTS FOR VAX
!
!      DATA IMACH( 1) /    5 /
!      DATA IMACH( 2) /    6 /
!      DATA IMACH( 3) /    7 /
!      DATA IMACH( 4) /    6 /
!      DATA IMACH( 5) /   32 /
!      DATA IMACH( 6) /    4 /
!      DATA IMACH( 7) /    2 /
!      DATA IMACH( 8) /   31 /
!      DATA IMACH( 9) / 2147483647 /
!      DATA IMACH(10) /    2 /
!      DATA IMACH(11) /   24 /
!      DATA IMACH(12) / -127 /
!      DATA IMACH(13) /  127 /
!      DATA IMACH(14) /   56 /
!      DATA IMACH(15) / -127 /
!      DATA IMACH(16) /  127 /
!
!     MACHINE CONSTANTS FOR THE SEQUENT BALANCE 8000 AND SVS FORTRAN ON
!     THE AT&T 7300 (UNIX PC)
!
!      DATA IMACH( 1) /     0 /
!      DATA IMACH( 2) /     0 /
!      DATA IMACH( 3) /     7 /
!      DATA IMACH( 4) /     0 /
!      DATA IMACH( 5) /    32 /
!      DATA IMACH( 6) /     1 /
!      DATA IMACH( 7) /     2 /
!      DATA IMACH( 8) /    31 /
!      DATA IMACH( 9) /  2147483647 /
!      DATA IMACH(10) /     2 /
!      DATA IMACH(11) /    24 /
!      DATA IMACH(12) /  -125 /
!      DATA IMACH(13) /   128 /
!      DATA IMACH(14) /    53 /
!      DATA IMACH(15) / -1021 /
!      DATA IMACH(16) /  1024 /
!
!     MACHINE CONSTANTS FOR THE RM FORTRAN ON THE AT&T 7300 (UNIX PC)
!
!      DATA IMACH( 1) /     5 /
!      DATA IMACH( 2) /     6 /
!      DATA IMACH( 3) /     7 /
!      DATA IMACH( 4) /     6 /
!      DATA IMACH( 5) /    32 /
!      DATA IMACH( 6) /     1 /
!      DATA IMACH( 7) /     2 /
!      DATA IMACH( 8) /    31 /
!      DATA IMACH( 9) /  2147483647 /
!      DATA IMACH(10) /     2 /
!      DATA IMACH(11) /    24 /
!      DATA IMACH(12) /  -125 /
!      DATA IMACH(13) /   128 /
!      DATA IMACH(14) /    53 /
!      DATA IMACH(15) / -1021 /
!      DATA IMACH(16) /  1024 /
!
      IF (I>=1 .AND. I<=16) THEN 
         I1MACH = IMACH(I) 
         RETURN  
      ENDIF 
      WRITE (OUTPUT, 1999) I 
 1999 FORMAT(' I1MACH - I OUT OF BOUNDS',I10) 
      STOP  
      END FUNCTION I1MACH 
