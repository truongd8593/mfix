!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: LEQ_DSLUGM(Vname, Var, A_m, B_m, M, ITMAX, IER)        C
!  Purpose: Incomplete factorization preconditioner with GMRES: DSLUGM C
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
      SUBROUTINE LEQ_DSLUGM(VNAME, VAR, A_M, B_M, M, ITMAX, IER) 
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
      Use dslugm_a
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
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: L, LNOW, IJK, N, NELT 
      INTEGER :: ISYM, ITOL, ITER, IERR, IUNIT 
      DOUBLE PRECISION :: ERR , TOL
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
! NELT   :IN       Integer.  -> set after filling A, IA, and JA
!         Number of Non-Zeros stored in A.
! IA     :IN       Integer IA(NELT).
! JA     :IN       Integer JA(NELT).
! A      :IN       Double Precision A(NELT).
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
!
! ISYM   :IN       Integer.
!         Flag to indicate symmetric storage format.
!         If ISYM=0, all nonzero entries of the matrix are stored.
!         If ISYM=1, the matrix is symmetric, and only the upper
!         or lower triangle of the matrix is stored.
      ISYM = 0 
!
! NSAVE  :IN       Integer.  -> set as a parameter
!         Number of direction vectors to save and orthogonalize against.
!         Must be greater than 1.
!
! ITOL   :IN       Integer.
!         Flag to indicate the type of convergence criterion used.
!         ITOL=0  Means the  iteration stops when the test described
!                 below on  the  residual RL  is satisfied.  This is
!                 the  "Natural Stopping Criteria" for this routine.
!                 Other values  of   ITOL  cause  extra,   otherwise
!                 unnecessary, computation per iteration and     are
!                 therefore  much less  efficient.  See  ISDGMR (the
!                 stop test routine) for more information.
!         ITOL=1  Means   the  iteration stops   when the first test
!                 described below on  the residual RL  is satisfied,
!                 and there  is either right  or  no preconditioning
!                 being used.
!         ITOL=2  Implies     that   the  user    is   using    left
!                 preconditioning, and the second stopping criterion
!                 below is used.
!         ITOL=3  Means the  iteration stops   when  the  third test
!                 described below on Minv*Residual is satisfied, and
!                 there is either left  or no  preconditioning begin
!                 used.
!         ITOL=11 is    often  useful  for   checking  and comparing
!                 different routines.  For this case, the  user must
!                 supply  the  "exact" solution or  a  very accurate
!                 approximation (one with  an  error much less  than
!                 TOL) through a common block,
!                     COMMON /SOLBLK/ SOLN(1)
!                 if ITOL=11, iteration stops when the 2-norm of the
!                 difference between the iterative approximation and
!                 the user-supplied solution  divided by the  2-norm
!                 of the  user-supplied solution  is  less than TOL.
!                 Note that this requires  the  user to  set up  the
!                 "COMMON     /SOLBLK/ SOLN(LENGTH)"  in the calling
!                 routine.  The routine with this declaration should
!                 be loaded before the stop test so that the correct
!                 length is used by  the loader.  This procedure  is
!                 not standard Fortran and may not work correctly on
!                 your   system (although  it  has  worked  on every
!                 system the authors have tried).  If ITOL is not 11
!                 then this common block is indeed standard Fortran.
      ITOL = 0 
!
! TOL    :INOUT    Double Precision.
!         Convergence criterion, as described below.  If TOL is set
!         to zero on input, then a default value of 500*(the smallest
!         positive magnitude, machine epsilon) is used.
      TOL = 0. 
!
! ITMAX  :IN       Integer.
!         Maximum number of iterations.  This routine uses the default
!         of NRMAX = ITMAX/NSAVE to determine the when each restart
!         should occur.  See the description of NRMAX and MAXL in
!         DGMRES for a full and frightfully interesting discussion of
!         this topic.
!       ITMAX = 1000
!
! ITER   :OUT      Integer.
!         Number of iterations required to reach convergence, or
!         ITMAX+1 if convergence criterion could not be achieved in
!         ITMAX iterations.
! ERR    :OUT      Double Precision.
!         Error estimate of error in final approximate solution, as
!         defined by ITOL.  Letting norm() denote the Euclidean
!         norm, ERR is defined as follows...
!         If ITOL=0, then ERR = norm(SB*(B-A*X(L)))/norm(SB*B),
!                               for right or no preconditioning, and
!                         ERR = norm(SB*(M-inverse)*(B-A*X(L)))/
!                                norm(SB*(M-inverse)*B),
!                               for left preconditioning.
!         If ITOL=1, then ERR = norm(SB*(B-A*X(L)))/norm(SB*B),
!                               since right or no preconditioning
!                               being used.
!         If ITOL=2, then ERR = norm(SB*(M-inverse)*(B-A*X(L)))/
!                                norm(SB*(M-inverse)*B),
!                               since left preconditioning is being
!                               used.
!         If ITOL=3, then ERR =  Max  |(Minv*(B-A*X(L)))(i)/x(i)|
!                               i=1,n
!         If ITOL=11, then ERR = norm(SB*(X(L)-SOLN))/norm(SB*SOLN).
! IERR   :OUT      Integer.
!         Return error flag.
!               IERR = 0 => All went well.
!               IERR = 1 => Insufficient storage allocated for
!                           RGWK or IGWK.
!               IERR = 2 => Routine DPIGMR failed to reduce the norm
!                           of the current residual on its last call,
!                           and so the iteration has stalled.  In
!                           this case, X equals the last computed
!                           approximation.  The user must either
!                           increase MAXL, or choose a different
!                           initial guess.
!               IERR =-1 => Insufficient length for RGWK array.
!                           IGWK(6) contains the required minimum
!                           length of the RGWK array.
!               IERR =-2 => Inconsistent ITOL and JPRE values.
!         For IERR <= 2, RGWK(1) = RHOL, which is the norm on the
!         left-hand-side of the relevant stopping test defined
!         below associated with the residual for the current
!         approximation X(L).
! IUNIT  :IN       Integer.
!         Unit number on which to write the error at each iteration,
!         if this is desired for monitoring convergence.  If unit
!         number is 0, no writing will occur.
      IUNIT = 0 
!
! RWORK  :WORK    Double Precision RWORK(LENW).
!         Double Precision array of size LENW.
!
! LENW   :IN       Integer.   -> set as a parameter
!         Length of the double precision workspace, RWORK.
!         LENW >= 1 + N*(NSAVE+7) +  NSAVE*(NSAVE+3)+NEL+NU.
!         For the recommended values,  RWORK
!         has size at least 131 + 17*N + NEL + NU.  Where  NEL is  the
!         number of non- zeros  in  the  lower triangle of  the matrix
!         (including the diagonal).  NU is the  number  of nonzeros in
!         the upper triangle of the matrix (including the diagonal).
!
! IWORK  :INOUT    Integer IWORK(LENIW).
!         Used to hold pointers into the RWORK array.
!         Upon return the following locations of IWORK hold information
!         which may be of use to the user:
!         IWORK(9)  Amount of Integer workspace actually used.
!         IWORK(10) Amount of Double Precision workspace actually used.
! LENIW  :IN       Integer.
!         Length of the integer workspace, IWORK.
!         LENIW >= NEL+NU+4*N+32.
!
      CALL DSLUGM (N, B_M(1,M), VAR, NELT, IA, JA, A, ISYM, NSAVE, ITOL, TOL, &
         ITMAX, ITER, ERR, IERR, IUNIT, RWORK, LENW, IWORK, LENIW) 
!
!      IF(IERR .EQ. 1) THEN
!        LINE = "Error: Insufficient work space allocated for LEq Solver"
!        CALL WRITE_ERROR('LEQ_DSLUGM: '//VName, LINE, 1)
!        STOP
!      ELSEIF(IERR .EQ. 2) THEN
!c        LINE = "Warning: LEq Solver iteration has stalled"
!c        CALL WRITE_ERROR('LEQ_DSLUGM: '//VName, LINE, 1)
!      ELSEIF(IERR .EQ. -1)THEN
!        WRITE(LINE, '(A, I6)')
!     &    "Error: Need to increase the length RWORK array to ", IWORK(6)
!        CALL WRITE_ERROR('LEQ_DSLUGM: '//VName, LINE, 1)
!        STOP
!      ELSEIF(IERR .EQ. -2) THEN
!        LINE = "Warning: Inconsistent ITOL and JPRE values"
!        CALL WRITE_ERROR('LEQ_DSLUGM: '//VName, LINE, 1)
!      ELSEIF(IERR .NE. 0)THEN
!        LINE = "Warning: Unknown error from LEq Solver"
!        CALL WRITE_ERROR('LEQ_DSLUGM: '//VName, LINE, 1)
!      ENDIF
      RETURN  
      END SUBROUTINE LEQ_DSLUGM 
