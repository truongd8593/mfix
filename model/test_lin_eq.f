!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: TEST_LIN_EQ(IJKMAX2, IJMAX2, IMAX2, A_m, TEST, DO_K,   C
!                  IER)                                                C
!  Purpose: Routine for testing the accuracy of linear equation solver C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 4-JUN-96   C
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
!
!
      SUBROUTINE TEST_LIN_EQ(IJKMAX2, IJMAX2, IMAX2, A_M, TEST, DO_K, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER IJKMAX2, IJMAX2, IMAX2, TEST, IER 
      LOGICAL DO_K 
      DOUBLE PRECISION, DIMENSION(IJKMAX2,-3:3) :: A_M 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      DOUBLE PRECISION, PARAMETER :: TOL = 1.0D-12 
      INTEGER, PARAMETER :: ITMAX = 1000 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: ISEED, IJK, IJKERR 
      DOUBLE PRECISION, DIMENSION(DIMENSION_3,-3:3) :: A 
      DOUBLE PRECISION, DIMENSION(DIMENSION_3) :: B, X, X_SOL 
      DOUBLE PRECISION :: ERR, ERRMAX, ERRSUM, XSUM 
      CHARACTER, DIMENSION(7) :: LINE*80 
      REAL  :: Harvest 
!-----------------------------------------------
!
!
!  Initialize the random number generator
!
      CALL RANDOM_SEED 
!
!  Fill the A and x arrays with random numbers, but ensuring that
!  the matrix is diagonally dominant
!
      DO IJK = 1, IJKMAX2 
         CALL RANDOM_NUMBER(HARVEST)
         X(IJK) = DBLE(HARVEST) 
         X_SOL(IJK) = 0.0 
         IF (TEST == 0) THEN 
            A(IJK,-3) = A_M(IJK,-3) 
            A(IJK,-2) = A_M(IJK,-2) 
            A(IJK,-1) = A_M(IJK,-1) 
            A(IJK,0) = A_M(IJK,0) 
            A(IJK,1) = A_M(IJK,1) 
            A(IJK,2) = A_M(IJK,2) 
            A(IJK,3) = A_M(IJK,3) 
         ELSE 
            CALL RANDOM_NUMBER(HARVEST)
            A(IJK,-3) = DBLE(HARVEST) 
            CALL RANDOM_NUMBER(HARVEST)
            A(IJK,-2) = DBLE(HARVEST) 
            CALL RANDOM_NUMBER(HARVEST)
            A(IJK,-1) = DBLE(HARVEST) 
            CALL RANDOM_NUMBER(HARVEST)
            A(IJK,0) = -DBLE(MAX(HARVEST,0.1))*70. 
            CALL RANDOM_NUMBER(HARVEST)
            A(IJK,1) = DBLE(HARVEST) 
            CALL RANDOM_NUMBER(HARVEST)
            A(IJK,2) = DBLE(HARVEST) 
            CALL RANDOM_NUMBER(HARVEST)
            A(IJK,3) = DBLE(HARVEST) 
         ENDIF 
      END DO 
      IJK = 1 
      IF (IJKMAX2 > 0) THEN 
         B(:IJKMAX2) = A(:IJKMAX2,0)*X(:IJKMAX2) 
         IJK = IJKMAX2 + 1 
      ENDIF 
      IJK = 2 
      IF (IJKMAX2 - 1 > 0) THEN 
         B(2:IJKMAX2) = B(2:IJKMAX2) + A(2:IJKMAX2,-1)*X(:IJKMAX2-1) 
         IJK = IJKMAX2 + 1 
      ENDIF 
      IJK = 1 
      IF (IJKMAX2 - 1 > 0) THEN 
         B(:IJKMAX2-1) = B(:IJKMAX2-1) + A(:IJKMAX2-1,1)*X(2:IJKMAX2) 
         IJK = IJKMAX2 
      ENDIF 
      IJK = IMAX2 + 1 
      IF (IJKMAX2 - IMAX2 > 0) THEN 
         B(IMAX2+1:IJKMAX2) = B(IMAX2+1:IJKMAX2) + A(IMAX2+1:IJKMAX2,-2)*X(:&
            IJKMAX2-IMAX2) 
         IJK = IJKMAX2 + 1 
      ENDIF 
      IJK = 1 
      IF (IJKMAX2 - IMAX2 > 0) THEN 
         B(:IJKMAX2-IMAX2) = B(:IJKMAX2-IMAX2) + A(:IJKMAX2-IMAX2,2)*X(1+IMAX2:&
            IJKMAX2) 
         IJK = IJKMAX2 - IMAX2 + 1 
      ENDIF 
      IF (DO_K) THEN 
         IJK = IJMAX2 + 1 
         IF (IJKMAX2 - IJMAX2 > 0) THEN 
            B(IJMAX2+1:IJKMAX2) = B(IJMAX2+1:IJKMAX2) + A(IJMAX2+1:IJKMAX2,-3)*&
               X(:IJKMAX2-IJMAX2) 
            IJK = IJKMAX2 + 1 
         ENDIF 
         IJK = 1 
         IF (IJKMAX2 - IJMAX2 > 0) THEN 
            B(:IJKMAX2-IJMAX2) = B(:IJKMAX2-IJMAX2) + A(:IJKMAX2-IJMAX2,3)*X(1+&
               IJMAX2:IJKMAX2) 
            IJK = IJKMAX2 - IJMAX2 + 1 
         ENDIF 
      ENDIF 
!
!  Solve the linear equation
!
      CALL SOLVE_LIN_EQ ('Test', X_SOL, A, B, 0, ITMAX, 3, 'II' , 1.0E-4, IER) 
!
!  Check the solution
!
      ERRSUM = 0.0 
      XSUM = 0.0 
      ERRMAX = 0.0 
      IJKERR = 0 
      DO IJK = 1, IJKMAX2 
         IF (X(IJK) /= 0.0) THEN 
            ERR = ABS(X_SOL(IJK)-X(IJK))/X(IJK) 
         ELSE IF (X_SOL(IJK) == 0.0) THEN 
            ERR = 0.0 
         ELSE 
            ERR = 1.0D32 
         ENDIF 
         IF (ERR > ERRMAX) THEN 
            ERRMAX = ERR 
            IJKERR = IJK 
         ENDIF 
         ERRSUM = ERRSUM + ABS(X_SOL(IJK)-X(IJK)) 
         XSUM = XSUM + ABS(X(IJK)) 
      END DO 
      IF (XSUM /= 0.0) THEN 
         ERR = ERRSUM/XSUM 
      ELSE IF (ERRSUM == 0.0) THEN 
         ERR = 0.0 
      ELSE 
         ERR = 1.0D32 
      ENDIF 
!
      IF (ERR < TOL) THEN 
         IER = 0 
         LINE(1) = 'Message: Lin equation solution satisfies tolerance.' 
      ELSE 
         IER = 1 
         LINE(1) = 'Error: Lin equation solution does not satisfy tolerance!' 
      ENDIF 
!
      WRITE (LINE(2), *) 'Average normalized error = ', ERR 
      WRITE (LINE(3), *) 'Max normalized error = ', ERRMAX 
      WRITE (LINE(4), *) 'Location of max error = ', IJKERR 
      WRITE (LINE(5), *) 'Sample values of actual (Xa) and solution (Xs):' 
      WRITE (LINE(6), '(A,G12.5, A, I6, A, G12.5, A, I6, A, G12.5)') 'Xa(1)=', &
         X(1), '  Xa(', IJKMAX2/2, ')=', X(IJKMAX2/2), '  Xa(', IJKMAX2, ')=', &
         X(IJKMAX2) 
      WRITE (LINE(7), '(A,G12.5, A, I6, A, G12.5, A, I6, A, G12.5)') 'Xs(1)=', &
         X_SOL(1), '  Xs(', IJKMAX2/2, ')=', X_SOL(IJKMAX2/2), '  Xs(', IJKMAX2&
         , ')=', X_SOL(IJKMAX2) 
!
      CALL WRITE_ERROR ('TEST_LIN_EQ', LINE, 7) 
!
      RETURN  
      END SUBROUTINE TEST_LIN_EQ 
 
    
      
      
