!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: LEQ_DSLUCS(Vname, Var, A_m, B_m, M, ITMAX, IER)        C
!  Purpose: Incomplete factorization preconditioner with BCGS: DSLUCS  C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 31-MAY-96  C
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
      SUBROUTINE LEQ_DSLUCS(VNAME, VAR, A_M, B_M, M, ITMAX, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE geometry
      USE indices
      Use dslucs_a
      USE compar        !//d
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER M, ITMAX, IER 
      CHARACTER VNAME*(*) 
      DOUBLE PRECISION, DIMENSION(DIMENSION_3) :: VAR 
      DOUBLE PRECISION, DIMENSION(DIMENSION_3,-3:3,0:DIMENSION_M) :: A_M 
      DOUBLE PRECISION, DIMENSION(DIMENSION_3,0:DIMENSION_M) :: B_M 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      DOUBLE PRECISION, PARAMETER :: TOL = 1.D-12 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: L, LNOW, IJK, N, NELT 
      INTEGER :: ISYM, ITOL, ITER, IERR, IUNIT 
      DOUBLE PRECISION :: ERR 
      CHARACTER :: LINE*80 
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      DOUBLE PRECISION , EXTERNAL :: DNRM2 
!-----------------------------------------------
      INCLUDE 'function.inc'
!
!
!     If the norm of rhs vector is zero do not call linear equation
!     solver.  Set the solution vector to 0 and return.
!
      IF (DNRM2(IJKMAX2,B_M(1,M),1) == ZERO) THEN 
         L = 1 
         IF (IJKMAX2 > 0) THEN 
            VAR(:IJKMAX2) = ZERO 
            L = IJKMAX2 + 1 
         ENDIF 
         RETURN  
      ENDIF 
!
!
!   Initialize
!
!
! *Arguments:
! N      :IN       Integer.
!         Order of the Matrix.
      N = IJKMAX2 
!
! B      :IN       Double Precision B(N). = B_m(1, M)
!         Right-hand side vector.
!
! X      :INOUT    Double Precision X(N). = Var
!         On input X is your initial guess for solution vector.
!         On output X is the final approximate solution.
!
! NELT   :IN       Integer.
!         Number of Non-Zeros stored in A.  -> set after filling A, IA,
! IA     :INOUT    Integer IA(NELT).
! JA     :INOUT    Integer JA(NELT).
! A      :INOUT    Double Precision A(NELT).
!         These arrays should hold the matrix A in either the SLAP
!         Triad format or the SLAP Column format.  See "Description",
!         below.  If the SLAP Triad format is chosen it is changed
!         internally to the SLAP Column format.
      DO L = 1, IJKMAX2 
         A(L) = A_M(L,0,M) 
         IA(L) = L 
         JA(L) = L 
      END DO 
      LNOW = IJKMAX2 
!
      DO L = 1, IJKMAX2 - 1                      !ijk = L + 1 
         IF (A_M(L+1,-1,M) /= ZERO) THEN 
            LNOW = LNOW + 1 
            A(LNOW) = A_M(L+1,-1,M) 
            IA(LNOW) = L + 1 
            JA(LNOW) = IM_OF(L + 1) 
         ENDIF 
      END DO 
      DO L = 1, IJKMAX2 - 1                      !ijk = L 
         IF (A_M(L,1,M) /= ZERO) THEN 
            LNOW = LNOW + 1 
            A(LNOW) = A_M(L,1,M) 
            IA(LNOW) = L 
            JA(LNOW) = IP_OF(L) 
         ENDIF 
      END DO 
      DO L = 1, IJKMAX2 - IMAX2                  !ijk = L+IMAX2 
         IF (A_M(L+IMAX2,-2,M) /= ZERO) THEN 
            LNOW = LNOW + 1 
            A(LNOW) = A_M(L+IMAX2,-2,M) 
            IA(LNOW) = L + IMAX2 
            JA(LNOW) = JM_OF(L + IMAX2) 
         ENDIF 
      END DO 
      DO L = 1, IJKMAX2 - IMAX2                  !ijk = L 
         IF (A_M(L,2,M) /= ZERO) THEN 
            LNOW = LNOW + 1 
            A(LNOW) = A_M(L,2,M) 
            IA(LNOW) = L 
            JA(LNOW) = JP_OF(L) 
         ENDIF 
      END DO 
      IF (DO_K) THEN 
!
         DO L = 1, IJKMAX2 - IJMAX2              !ijk = L+IJMAX2 
            IF (A_M(L+IJMAX2,-3,M) /= ZERO) THEN 
               LNOW = LNOW + 1 
               A(LNOW) = A_M(L+IJMAX2,-3,M) 
               IA(LNOW) = L + IJMAX2 
               JA(LNOW) = KM_OF(L + IJMAX2) 
            ENDIF 
         END DO 
         DO L = 1, IJKMAX2 - IJMAX2              !ijk = L 
            IF (A_M(L,3,M) /= ZERO) THEN 
               LNOW = LNOW + 1 
               A(LNOW) = A_M(L,3,M) 
               IA(LNOW) = L 
               JA(LNOW) = KP_OF(L) 
            ENDIF 
         END DO 
      ENDIF 
!
!
      NELT = LNOW 
! ISYM   :IN       Integer.
!         Flag to indicate symmetric storage format.
!         If ISYM=0, all nonzero entries of the matrix are stored.
!         If ISYM=1, the matrix is symmetric, and only the upper
!         or lower triangle of the matrix is stored.
      ISYM = 0 
!
! ITOL   :IN       Integer.
!         Flag to indicate type of convergence criterion.
!         If ITOL=1, iteration stops when the 2-norm of the residual
!         divided by the 2-norm of the right-hand side is less than TOL.
!         This routine must calculate the residual from R = A*X - B.
!         This is un-natural and hence expensive for this type of iter-
!         ative method.  ITOL=2 is *STRONGLY* recommended.
!         If ITOL=2, iteration stops when the 2-norm of M-inv times the
!         residual divided by the 2-norm of M-inv times the right hand
!         side is less than tol, where M-inv time a vector is the pre-
!         conditioning step.  This is the *NATURAL* stopping for this
!         iterative method and is *STRONGLY* recommended.
      ITOL = 2 
!
! TOL    :IN       Double Precision. -> set as a parameter
!         Convergence criterion, as described above.
!
! ITMAX  :IN       Integer. -> set as a parameter
!         Maximum number of iterations.
!
! ITER   :OUT      Integer.
!         Number of iterations required to reach convergence, or
!         ITMAX+1 if convergence criterion could not be achieved in
!         ITMAX iterations.
! ERR    :OUT      Double Precision.
!         Error estimate of error in final approximate solution, as
!         defined by ITOL.
! IERR   :OUT      Integer.
!         Return error flag.
!           IERR = 0 => All went well.
!           IERR = 1 => Insufficient storage allocated
!                       for WORK or IWORK.
!           IERR = 2 => Method failed to converge in
!                       ITMAX steps.
!           IERR = 3 => Error in user input.  Check input
!                       value of N, ITOL.
!           IERR = 4 => User error tolerance set too tight.
!                       Reset to 500.0*D1MACH(3).  Iteration proceeded.
!           IERR = 5 => Breakdown of the method detected.
!                       $(r0,r) approximately 0.0$.
!           IERR = 6 => Stagnation of the method detected.
!                        $(r0,v) approximately 0.0$.
!           IERR = 7 => Incomplete factorization broke down
!                       and was fudged.  Resulting preconditioning may
!                       be less than the best.
! IUNIT  :IN       Integer.
!         Unit number on which to write the error at each iteration,
!         if this is desired for monitoring convergence.  If unit
!         number is 0, no writing will occur.
      IUNIT = 0 
!
! RWORK  :WORK     Double Precision RWORK(LENW).
!         Double Precision array used for workspace.  NEL is the
!         number of non-
!         zeros in the lower triangle of the matrix (including the
!         diagonal).  NU is the number of nonzeros in the upper
!         triangle of the matrix (including the diagonal).
! LENW   :IN       Integer.   -> set as a parameter
!         Length of the double precision workspace, RWORK.
!         LENW >= NEL+NU+8*N.
! IWORK  :WORK     Integer IWORK(LENIW).
!         Integer array used for workspace.  NEL is the number of non-
!         zeros in the lower triangle of the matrix (including the
!         diagonal).  NU is the number of nonzeros in the upper
!         triangle of the matrix (including the diagonal).
!         Upon return the following locations of IWORK hold information
!         which may be of use to the user:
!         IWORK(9)  Amount of Integer workspace actually used.
!         IWORK(10) Amount of Double Precision workspace actually used.
! LENIW  :IN       Integer.   -> set as a parameter
!         Length of the integer workspace, IWORK.
!         LENIW >= NEL+NU+4*N+12.
!
      CALL DSLUCS (N, B_M(1,M), VAR, NELT, IA, JA, A, ISYM, ITOL, TOL, ITMAX, &
         ITER, ERR, IERR, IUNIT, RWORK, LENW, IWORK, LENIW) 
!
!      IF(IERR .EQ. 1) THEN
!        LINE = "Error: Insufficient work space allocated for LEq Solver"
!        CALL WRITE_ERROR('LEQ_DSLUCS', LINE, 1)
!        STOP
!      ELSEIF(IERR .EQ. 2) THEN
!        LINE = "Warning: Max # of iterations exceeded in LEq Solver"
!        CALL WRITE_ERROR('LEQ_DSLUCS', LINE, 1)
!      ELSEIF(IERR .EQ. 3)THEN
!        WRITE(LINE, '(A, I6)')
!     &    "Error: Check input value of N and ITOL in LEq Solver"
!        CALL WRITE_ERROR('LEQ_DSLUCS', LINE, 1)
!        STOP
!      ELSEIF(IERR .EQ. 4) THEN
!        LINE = "Warning: Tolerance too tight. Tolerance reset."
!        CALL WRITE_ERROR('LEQ_DSLUCS', LINE, 1)
!      ELSEIF(IERR .EQ. 5) THEN
!        LINE = "Error: LEq solver method broke down."
!        CALL WRITE_ERROR('LEQ_DSLUCS', LINE, 1)
!        STOP
!      ELSEIF(IERR .EQ. 6) THEN
!        LINE = "Warning: LEq iterations stagnate."
!        CALL WRITE_ERROR('LEQ_DSLUCS', LINE, 1)
!      ELSEIF(IERR .EQ. 7) THEN
!        LINE = "Warning: Incomplete factorization broke down."
!        CALL WRITE_ERROR('LEQ_DSLUCS', LINE, 1)
!      ELSEIF(IERR .NE. 0)THEN
!        LINE = "Warning: Unknown error from LEq Solver"
!        CALL WRITE_ERROR('LEQ_DSLUCS', LINE, 1)
!      ENDIF
!     Incomplete factorization preconditioner with BCGS: DSLUCS
!------------------------------------------------------------------------------
      RETURN  
      END SUBROUTINE LEQ_DSLUCS 
