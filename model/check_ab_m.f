!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: Check_Ab_m(A_m, b_m, M, src, IER)                      C                     C
!  Purpose: Check the entries of the sparse matrix coefficients and theC
!           source vector, if src is set true.                         C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 16-MAY-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CHECK_AB_M(A_M, B_M, M, SRC, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE matrix 
      USE geometry
      USE indices
      USE compar   !//d
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Error index
      INTEGER          IER
!
!                      Phase index
      INTEGER          M
!
!                      Check source.  Check source term only for
!                      variables, such as temperatures and
!                      mass fractions, that are always positive.
      LOGICAL          SRC
!
!                      Error message
      CHARACTER*80     LINE(1)
!
!                      cell index
      INTEGER          IJK
!
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
!
!                      Source vector
      DOUBLE PRECISION b_m(DIMENSION_3, 0:DIMENSION_M)
!-----------------------------------------------
      INCLUDE 'function.inc'
!
      IER = 0 
      DO IJK = 1, IJKMAX2 
         IF (.NOT.WALL_AT(IJK)) THEN 
            IF (A_M(IJK,B,M) < ZERO) THEN 
               IF (ABS(A_M(IJK,B,M)) > SMALL_NUMBER) THEN 
                  WRITE (LINE(1), *) 'Error: Diagonal-b < 0. Phase = ', M, &
                     ' IJK = ', IJK 
                  CALL WRITE_ERROR ('CHECK_Ab_m', LINE, 1) 
                  GO TO 500 
               ELSE 
                  A_M(IJK,B,M) = ZERO 
               ENDIF 
            ENDIF 
            IF (A_M(IJK,S,M) < ZERO) THEN 
               IF (ABS(A_M(IJK,S,M)) > SMALL_NUMBER) THEN 
                  WRITE (LINE(1), *) 'Error: Diagonal-s < 0. Phase = ', M, &
                     ' IJK = ', IJK 
                  CALL WRITE_ERROR ('CHECK_Ab_m', LINE, 1) 
                  GO TO 500 
               ELSE 
                  A_M(IJK,S,M) = ZERO 
               ENDIF 
            ENDIF 
            IF (A_M(IJK,W,M) < ZERO) THEN 
               IF (ABS(A_M(IJK,W,M)) > SMALL_NUMBER) THEN 
                  WRITE (LINE(1), *) 'Error: Diagonal-w < 0. Phase = ', M, &
                     ' IJK = ', IJK 
                  CALL WRITE_ERROR ('CHECK_Ab_m', LINE, 1) 
                  GO TO 500 
               ELSE 
                  A_M(IJK,W,M) = ZERO 
               ENDIF 
            ENDIF 
            IF (A_M(IJK,E,M) < ZERO) THEN 
               IF (ABS(A_M(IJK,E,M)) > SMALL_NUMBER) THEN 
                  WRITE (LINE(1), *) 'Error: Diagonal-e < 0. Phase = ', M, &
                     ' IJK = ', IJK 
                  CALL WRITE_ERROR ('CHECK_Ab_m', LINE, 1) 
                  GO TO 500 
               ELSE 
                  A_M(IJK,E,M) = ZERO 
               ENDIF 
            ENDIF 
            IF (A_M(IJK,N,M) < ZERO) THEN 
               IF (ABS(A_M(IJK,N,M)) > SMALL_NUMBER) THEN 
                  WRITE (LINE(1), *) 'Error: Diagonal-n < 0. Phase = ', M, &
                     ' IJK = ', IJK 
                  CALL WRITE_ERROR ('CHECK_Ab_m', LINE, 1) 
                  GO TO 500 
               ELSE 
                  A_M(IJK,N,M) = ZERO 
               ENDIF 
            ENDIF 
            IF (A_M(IJK,T,M) < ZERO) THEN 
               IF (ABS(A_M(IJK,T,M)) > SMALL_NUMBER) THEN 
                  WRITE (LINE(1), *) 'Error: Diagonal-t < 0. Phase = ', M, &
                     ' IJK = ', IJK 
                  CALL WRITE_ERROR ('CHECK_Ab_m', LINE, 1) 
                  GO TO 500 
               ELSE 
                  A_M(IJK,T,M) = ZERO 
               ENDIF 
            ENDIF 
            IF (A_M(IJK,0,M) >= ZERO) THEN 
               WRITE (LINE(1), *) 'Error: Main Diagonal >= 0. Phase = ', M, &
                  ' IJK = ', IJK 
               CALL WRITE_ERROR ('CHECK_Ab_m', LINE, 1) 
               GO TO 500 
            ENDIF 
            IF (SRC) THEN 
               IF (B_M(IJK,M) > ZERO) THEN 
                  IF (B_M(IJK,M) > SMALL_NUMBER) THEN 
                     WRITE (LINE(1), *) 'Error: Source term >0. Phase = ', M, &
                        ' IJK = ', IJK 
                     CALL WRITE_ERROR ('CHECK_Ab_m', LINE, 1) 
                     GO TO 500 
                  ELSE 
                     B_M(IJK,M) = ZERO 
                  ENDIF 
               ENDIF 
            ENDIF 
         ENDIF 
      END DO 
      RETURN  
!
!     error condition
  500 CONTINUE 
      IER = 1 
      CALL WRITE_AB_M (A_M, B_M, IJKMAX2, M, IER) 
      STOP  
      END SUBROUTINE CHECK_AB_M 
