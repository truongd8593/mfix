!********************* MODIFIED IGCG ROUTINES ************************
!*                                                                   *
!*   THIS SUBROUTINE COMPUTES THE SPLITTING MATRIX.  THE BANDS       *
!*   OF MATRIX 'B' ARE GIVEN NAMES ACCORDING TO THE FOLLOWING RULE:  *
!*                                                                   *
!*   <BD00>: DIAGONAL OF THE ORIGINAL UNNORMALIZED MATRIX 'B'        *
!*   <BUX> : SIDE-BAND NO. 'X' (COUNTED FROM DIAGONAL) OF UPPER      *
!*           TRIANGULAR PART OF 'B'.                                 *
!*   <BLX> : SIDE-BAND NO. 'X' (COUNTED FROM DIAGONAL) OF LOWER      *
!*           TRIANGULAR PART OF 'B'.                                 *
!*                                                                   *
!*   THE MATRIX ELEMENTS OF THE FACTORIZATION HAVE THE SAME NAMES    *
!*   WITH A POSTFIX 'N'                                              *
!*                                                                   *
!*   ENTRY-POINTS: <LUZERO> --- PRESETS SPLITTING MATRIX WITH ZEROS  *
!*                    <ILU> --- COMPUTES THE INCOMPLETE FACTORIZATION*
!*                                                                   *
!*********************************************************************
!
!
      SUBROUTINE LUZERO(BL01N, BL03N, BL09N, BU01N, BU03N, BU09N) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:24:36  12/09/98  
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE IGCG_I 
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      DOUBLE PRECISION, DIMENSION(DIMENSION_3 + DIM_I1) :: BL01N 
      DOUBLE PRECISION, DIMENSION(DIMENSION_3 + DIM_I3) :: BL03N 
      DOUBLE PRECISION, DIMENSION(DIMENSION_3 + DIM_I9) :: BL09N 
      DOUBLE PRECISION, DIMENSION(1 - DIM_I1:DIMENSION_3) :: BU01N 
      DOUBLE PRECISION, DIMENSION(1 - DIM_I3:DIMENSION_3) :: BU03N 
      DOUBLE PRECISION, DIMENSION(1 - DIM_I9:DIMENSION_3) :: BU09N 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I 
!-----------------------------------------------
!
!    ARGUMENT LIST
!
!
!
!     LOCAL VARIABLES
!
!
!     PRESET ALL AUXILIARY ELEMENTS OF SPLITTING MATRIX WITH ZERO
!
      I = 1 
      IF (I1 > 0) THEN 
         BL01N(NVOL+1:I1+NVOL) = 0. 
         BU01N(:1-I1:(-1)) = 0. 
         I = I1 + 1 
      ENDIF 
      I = 1 
      IF (I3 > 0) THEN 
         BL03N(NVOL+1:I3+NVOL) = 0. 
         BU03N(0:1-I3:(-1)) = 0. 
         I = I3 + 1 
      ENDIF 
      I = 1 
      IF (I9 > 0) THEN 
         BL09N(NVOL+1:I9+NVOL) = 0. 
         BU09N(0:1-I9:(-1)) = 0. 
         I = I9 + 1 
      ENDIF 
      RETURN  
      END SUBROUTINE LUZERO 
!
!
!
      SUBROUTINE ILU(BL01, BL03, BL09, BU01, BU03, BU09, BL01N, BL03N, BL09N, &
         BU01N, BU03N, BU09N, BD00N) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:24:36  12/09/98  
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE IGCG_I 
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      DOUBLE PRECISION, DIMENSION(NVOL) :: BL01, BL03, BL09, BU01, BU03, BU09 
      DOUBLE PRECISION, DIMENSION(DIMENSION_3 + DIM_I1) :: BL01N 
      DOUBLE PRECISION, DIMENSION(DIMENSION_3 + DIM_I3) :: BL03N 
      DOUBLE PRECISION, DIMENSION(DIMENSION_3 + DIM_I9) :: BL09N 
      DOUBLE PRECISION, DIMENSION(1 - DIM_I1:DIMENSION_3) :: BU01N 
      DOUBLE PRECISION, DIMENSION(1 - DIM_I3:DIMENSION_3) :: BU03N 
      DOUBLE PRECISION, DIMENSION(1 - DIM_I9:DIMENSION_3) :: BU09N 
      DOUBLE PRECISION, DIMENSION(NVOL) :: BD00N 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, IND 
!-----------------------------------------------
!
!    ARGUMENT LIST
!
!
!
!
!     LOCAL VARIABLES
!
!
!  THE MATRIX 'B' IS DECOMPOSED INTO TWO TRIANGULAR MATRICES 'L' AND
!  'U'.  ONLY ELEMENTS ARE SET WHICH ARE NON-ZERO IN 'B'.  THIS MEANS
!  B = LU + E WITH SOME ERROR MATRIX E, WHICH HAS VANISHING ELEMENTS
!  IN THE NON-ZERO BANDS OF 'B'.  THE DECOMPOSITION IS DONE COLUMN BY
!  COLUMN.  7 TYPES OF COLUMNS HAVE TO BE CONSIDERED.
!
      IND = 1 
      IF (NVOL > 0) THEN 
         BU09N(1:NVOL) = BU09 
         BU03N(1:NVOL) = BU03 
         BU01N(1:NVOL) = BU01 
         BL09N(:NVOL) = BL09 
         BL03N(:NVOL) = BL03 
         BL01N(:NVOL) = BL01 
         BD00N = 1. 
         IND = NVOL + 1 
      ENDIF 
      DO I = 1, NVOL 
         BD00N(I) = BD00N(I) - BL01N(I)*BU01N(I-I1) - BL03N(I)*BU03N(I-I3) - &
            BL09N(I)*BU09N(I-I9) 
         BD00N(I) = 1./BD00N(I) 
         BL01N(I+I1) = BL01N(I+I1)*BD00N(I) 
         BL03N(I+I3) = BL03N(I+I3)*BD00N(I) 
         BL09N(I+I9) = BL09N(I+I9)*BD00N(I) 
      END DO 
      IND = 1 
      IF (NVOL > 0) THEN 
         BU09N(1:NVOL) = BU09N(1:NVOL)*BD00N 
         BU03N(1:NVOL) = BU03N(1:NVOL)*BD00N 
         BU01N(1:NVOL) = BU01N(1:NVOL)*BD00N 
         IND = NVOL + 1 
      ENDIF 
      RETURN  
      END SUBROUTINE ILU 
!
!
!
      SUBROUTINE IGCG(BL01, BL03, BL09, BU01, BU03, BU09, BL01N, BL03N, BL09N, &
         BU01N, BU03N, BU09N, BD00N, FPRES, XPRES, ITMAX, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:24:36  12/09/98  
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE IGCG_I 
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER ITMAX, IER 
      DOUBLE PRECISION,DIMENSION(DIMENSION_3)::BL01,BL03,BL09,BU01,BU03,BU09 
      DOUBLE PRECISION, DIMENSION(DIMENSION_3 + DIM_I1) :: BL01N 
      DOUBLE PRECISION, DIMENSION(DIMENSION_3 + DIM_I3) :: BL03N 
      DOUBLE PRECISION, DIMENSION(DIMENSION_3 + DIM_I9) :: BL09N 
      DOUBLE PRECISION, DIMENSION(1 - DIM_I1:DIMENSION_3) :: BU01N 
      DOUBLE PRECISION, DIMENSION(1 - DIM_I3:DIMENSION_3) :: BU03N 
      DOUBLE PRECISION, DIMENSION(1 - DIM_I9:DIMENSION_3) :: BU09N 
      DOUBLE PRECISION, DIMENSION(DIMENSION_3) :: BD00N, FPRES, XPRES 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      DOUBLE PRECISION, PARAMETER :: RESMAX = 1.0E-20 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: ITER, IND 
      DOUBLE PRECISION,DIMENSION(DIMENSION_3)::DE,DLBPP,LBDE,LBPP,PP,UDE 
      DOUBLE PRECISION :: ALF, FACLAM, FNORM, PROD, RES, RESVGL, VAR 
!-----------------------------------------------
!
!  ARGUMENT LIST
!
!
!
!
!
!
!
!
!  LOCAL VARIABLES
!
!
!
      IER = 0 
      ITER = 0 
      FNORM = 0. 
      CALL LIMUL (DE, FPRES, BL01N, BL03N, BL09N) 
      CALL UIMUL (UDE, DE, BU01N, BU03N, BU09N, BD00N) 
      IND = 1 
      IF (NVOL > 0) THEN 
         FNORM = DOT_PRODUCT(UDE(:NVOL),UDE(:NVOL)) 
         IND = NVOL + 1 
      ENDIF 
      IF (FNORM > 1.0) THEN 
         FNORM = 1./FNORM 
      ELSE 
         FNORM = 1.0 
      ENDIF 
      CALL BMUL (DE, XPRES, BL01, BL03, BL09, BU01, BU03, BU09) 
      IND = 1 
      IF (NVOL > 0) THEN 
         DE(:NVOL) = FPRES(:NVOL) - DE(:NVOL) 
         IND = NVOL + 1 
      ENDIF 
      CALL LIMUL (UDE, DE, BL01N, BL03N, BL09N) 
      CALL UIMUL (DE, UDE, BU01N, BU03N, BU09N, BD00N) 
      IND = 1 
      IF (NVOL > 0) THEN 
         RES = DOT_PRODUCT(DE(:NVOL),DE(:NVOL)) 
         IND = NVOL + 1 
      ENDIF 
      RES = RES*FNORM 
      IF (RES > RESMAX) THEN 
         RES = 0. 
         IND = 1 
         IF (NVOL > 0) THEN 
            RES = DOT_PRODUCT(XPRES(:NVOL),XPRES(:NVOL)) 
            IND = NVOL + 1 
         ENDIF 
         RESVGL = RES 
         CALL BMUL (FPRES, DE, BL01, BL03, BL09, BU01, BU03, BU09) 
         CALL LIMUL (LBDE, FPRES, BL01N, BL03N, BL09N) 
         IND = 1 
         IF (NVOL > 0) THEN 
            PP(:NVOL) = DE(:NVOL) 
            LBPP(:NVOL) = LBDE(:NVOL) 
            DLBPP(:NVOL) = LBDE(:NVOL)*BD00N(:NVOL) 
            IND = NVOL + 1 
         ENDIF 
         VAR = 0. 
         IND = 1 
         IF (NVOL > 0) THEN 
            VAR = DOT_PRODUCT(DLBPP(:NVOL),LBPP(:NVOL)) 
            IND = NVOL + 1 
         ENDIF 
         PROD = VAR 
         DO ITER = 1, ITMAX 
            IND = 1 
            IF (NVOL > 0) THEN 
               FACLAM = SUM(UDE(:NVOL)*DLBPP(:NVOL)) 
               IND = NVOL + 1 
            ENDIF 
            IF (ABS(PROD) < 1.0E-32) GO TO 999 
            FACLAM = FACLAM/PROD 
            IND = 1 
            IF (NVOL > 0) THEN 
               XPRES(:NVOL) = XPRES(:NVOL) + FACLAM*PP(:NVOL) 
               UDE(:NVOL) = UDE(:NVOL) - FACLAM*LBPP(:NVOL) 
               IND = NVOL + 1 
            ENDIF 
            CALL UIMUL (DE, UDE, BU01N, BU03N, BU09N, BD00N) 
            IND = 1 
            IF (NVOL > 0) THEN 
               RES = SUM(XPRES(:NVOL)*XPRES(:NVOL)) 
               IND = NVOL + 1 
            ENDIF 
            RES = RES*FNORM 
            IF (ABS(RES - RESVGL) <= RESMAX) GO TO 999 
            CALL BMUL (FPRES, DE, BL01, BL03, BL09, BU01, BU03, BU09) 
            CALL LIMUL (LBDE, FPRES, BL01N, BL03N, BL09N) 
            VAR = 0. 
            IND = 1 
            IF (NVOL > 0) THEN 
               VAR = SUM(LBDE(:NVOL)*DLBPP(:NVOL)) 
               IND = NVOL + 1 
            ENDIF 
            IF (ABS(PROD) < 1.0E-32) GO TO 999 
            ALF = -VAR/PROD 
            IND = 1 
            IF (NVOL > 0) THEN 
               PP(:NVOL) = DE(:NVOL) + ALF*PP(:NVOL) 
               LBPP(:NVOL) = LBDE(:NVOL) + ALF*LBPP(:NVOL) 
               IND = NVOL + 1 
            ENDIF 
            IND = 1 
            IF (NVOL > 0) THEN 
               DLBPP(:NVOL) = LBPP(:NVOL)*BD00N(:NVOL) 
               IND = NVOL + 1 
            ENDIF 
            VAR = 0. 
            IND = 1 
            IF (NVOL > 0) THEN 
               VAR = SUM(LBPP(:NVOL)*DLBPP(:NVOL)) 
               IND = NVOL + 1 
            ENDIF 
            PROD = VAR 
            RESVGL = RES 
         END DO 
         IER = 1 
      ENDIF 
  999 CONTINUE 
      RETURN  
      END SUBROUTINE IGCG 
!
!
!***********************************************************************
!*                                                                     *
!*  THIS SECTION CONTAINS THREE SUBROUTINES FOR MULTIPLICATION         *
!*  OF MATRICES WITH VECTORS.                                          *
!*                                                                     *
!*  ENTRIES:  <LIMUL> MULTIPLIES THE INVERSE OF 'L' WITH SOME          *
!*                    VECTOR 'WK' TO OBTAIN 'WKS'                      *
!*            <UIMUL> MULTIPLIES THE INVERSE OF 'U' WITH SOME          *
!*                    VECTOR 'WK' TO OBTAIN 'WKS'                      *
!*            <BMUL>  MULTIPLIES 'B' WITH SOME VECTOR 'WK' TO          *
!*                    OBTAIN 'WKS'                                     *
!*                                                                     *
!***********************************************************************
      SUBROUTINE LIMUL(WKS, WK, BL01N, BL03N, BL09N) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:24:36  12/09/98  
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE IGCG_I 
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      DOUBLE PRECISION, DIMENSION(DIMENSION_3) :: WKS, WK 
      DOUBLE PRECISION, DIMENSION(DIMENSION_3 + DIM_I1) :: BL01N 
      DOUBLE PRECISION, DIMENSION(DIMENSION_3 + DIM_I3) :: BL03N 
      DOUBLE PRECISION, DIMENSION(DIMENSION_3 + DIM_I9) :: BL09N 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I 
!-----------------------------------------------
!
!     ARGUMENT LIST
!
!
!     LOCAL VARIABLES
!
!
!  'L' IS THE LOWER TRIANGULAR MATRIX DERIVED FROM 27-POINT-STAR.
!  IT HAS UNIT DIAGONAL
!
      WKS(1) = WK(1) 
      DO I = I1 + 1, I3 
         WKS(I) = WK(I) - BL01N(I)*WKS(I-I1) 
      END DO 
      DO I = I3 + 1, I9 
         WKS(I) = WK(I) - BL01N(I)*WKS(I-I1) - BL03N(I)*WKS(I-I3) 
      END DO 
      DO I = I9 + 1, NVOL 
         WKS(I) = WK(I) - BL01N(I)*WKS(I-I1) - BL03N(I)*WKS(I-I3) - BL09N(I)*&
            WKS(I-I9) 
      END DO 
      RETURN  
      END SUBROUTINE LIMUL 
!
!
      SUBROUTINE UIMUL(WKS, WK, BU01N, BU03N, BU09N, BD00N) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:24:36  12/09/98  
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE IGCG_I 
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      DOUBLE PRECISION, DIMENSION(DIMENSION_3) :: WKS, WK 
      DOUBLE PRECISION, DIMENSION(1 - DIM_I1:DIMENSION_3) :: BU01N 
      DOUBLE PRECISION, DIMENSION(1 - DIM_I3:DIMENSION_3) :: BU03N 
      DOUBLE PRECISION, DIMENSION(1 - DIM_I9:DIMENSION_3) :: BU09N 
      DOUBLE PRECISION, DIMENSION(DIMENSION_3) :: BD00N 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I 
!-----------------------------------------------
!
!     ARGUMENT LIST
!
!
!
!     LOCAL VARIABLES
!
!
!  'U' IS THE UPPER TRIANGULAR MATRIX DERIVED FROM 27-POINT-STAR
!
      I = 1 
      IF (NVOL > 0) THEN 
         WKS(:NVOL) = WK(:NVOL)*BD00N(:NVOL) 
         I = NVOL + 1 
      ENDIF 
      DO I = NVOL - I1, NVOL - I3 + 1, -1 
         WKS(I) = WKS(I) - BU01N(I)*WKS(I+I1) 
      END DO 
      DO I = NVOL - I3, NVOL - I9 + 1, -1 
         WKS(I) = WKS(I) - BU01N(I)*WKS(I+I1) - BU03N(I)*WKS(I+I3) 
      END DO 
      DO I = NVOL - I9, 1, -1 
         WKS(I) = WKS(I) - BU01N(I)*WKS(I+I1) - BU03N(I)*WKS(I+I3) - BU09N(I)*&
            WKS(I+I9) 
      END DO 
      RETURN  
      END SUBROUTINE UIMUL 
!
!
      SUBROUTINE BMUL(WKS, WK, BL01, BL03, BL09, BU01, BU03, BU09) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:24:36  12/09/98  
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE IGCG_I 
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      DOUBLE PRECISION, DIMENSION(DIMENSION_3) :: WKS, WK, BL01, BL03, BL09, &
         BU01, BU03, BU09 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I 
!-----------------------------------------------
!
!     ARGUMENT LIST
!
!
!     LOCAL VARIABLES
!
!
      I = 1 
      IF (NVOL > 0) THEN 
         WKS(:NVOL) = WK(:NVOL) 
         I = NVOL + 1 
      ENDIF 
      I = 1 
      IF (NVOL - I1 > 0) THEN 
         WKS(:NVOL-I1) = WKS(:NVOL-I1) + BU01(:NVOL-I1)*WK(1+I1:NVOL) 
         I = NVOL - I1 + 1 
      ENDIF 
      I = 1 
      IF (NVOL - I3 > 0) THEN 
         WKS(:NVOL-I3) = WKS(:NVOL-I3) + BU03(:NVOL-I3)*WK(1+I3:NVOL) 
         I = NVOL - I3 + 1 
      ENDIF 
      I = 1 
      IF (NVOL - I9 > 0) THEN 
         WKS(:NVOL-I9) = WKS(:NVOL-I9) + BU09(:NVOL-I9)*WK(1+I9:NVOL) 
         I = NVOL - I9 + 1 
      ENDIF 
      I = I1 + 1 
      IF (NVOL - I1 > 0) THEN 
         WKS(I1+1:NVOL) = WKS(I1+1:NVOL) + BL01(I1+1:NVOL)*WK(:NVOL-I1) 
         I = NVOL + 1 
      ENDIF 
      I = I3 + 1 
      IF (NVOL - I3 > 0) THEN 
         WKS(I3+1:NVOL) = WKS(I3+1:NVOL) + BL03(I3+1:NVOL)*WK(:NVOL-I3) 
         I = NVOL + 1 
      ENDIF 
      I = I9 + 1 
      IF (NVOL - I9 > 0) THEN 
         WKS(I9+1:NVOL) = WKS(I9+1:NVOL) + BL09(I9+1:NVOL)*WK(:NVOL-I9) 
         I = NVOL + 1 
      ENDIF 
      RETURN  
      END SUBROUTINE BMUL 
!
!------------------------------------------------------------------------------
!
!  subroutine for initializing I's
!  Added for interfacing with MFIX
!
      SUBROUTINE IGCG_INIT 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:24:36  12/09/98  
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE geometry
      USE IGCG_I 
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NCOL, NAREA 
!-----------------------------------------------
!
!
!
!
      NCOL = IMAX2 
      NAREA = IJMAX2 
      NVOL = IJKMAX2 
!
      I1 = 1 
      I2 = NCOL - 1 
      I3 = NCOL 
      I4 = NCOL + 1 
      I5 = NAREA - NCOL - 1 
      I6 = NAREA - NCOL 
      I7 = NAREA - NCOL + 1 
      I8 = NAREA - 1 
      I9 = NAREA 
      I10 = NAREA + 1 
      I11 = NAREA + NCOL - 1 
      I12 = NAREA + NCOL 
      I13 = NAREA + NCOL + 1 
!
      RETURN  
      END SUBROUTINE IGCG_INIT 
 
 
