!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: LEQ_IGCG(Vname, Var, A_m, B_m, M, ITMAX, IER)          C
!  Purpose: Conjugate gradient type preconditioned by an incomplete    C
!           factorization.  by Kapitza & Eppel, 1987                   C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 9-AUG-96   C
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
      SUBROUTINE LEQ_IGCG(VNAME, VAR, A_M, B_M, M, ITMAX, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE geometry
      USE matrix 
      USE igcg_a
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!                      Error index
      INTEGER          IER
!
!                      phase index
      INTEGER          M
!
!                      maximum number of iterations
      INTEGER          ITMAX
!
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
!
!                      Vector b_m
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M)
!
!                      Variable name
      CHARACTER*(*)    Vname
!
!                      Variable
      DOUBLE PRECISION Var(DIMENSION_3)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IJK 
      DOUBLE PRECISION :: BD 
      CHARACTER :: LINE*80 
!-----------------------------------------------
!
!
!
      CALL LUZERO (BL01N, BL03N, BL09N, BU01N, BU03N, BU09N) 
!
!
!  NORMALIZE THE EQUATION SET
!
!$omp parallel do private(IJK,BD)
      DO IJK = 1, IJKMAX2 
         BD = A_M(IJK,0,M) 
         A_M(IJK,E,M) = A_M(IJK,E,M)/BD 
         A_M(IJK,W,M) = A_M(IJK,W,M)/BD 
         A_M(IJK,N,M) = A_M(IJK,N,M)/BD 
         A_M(IJK,S,M) = A_M(IJK,S,M)/BD 
         A_M(IJK,T,M) = A_M(IJK,T,M)/BD 
         A_M(IJK,B,M) = A_M(IJK,B,M)/BD 
         B_M(IJK,M) = B_M(IJK,M)/BD 
      END DO 
      CALL ILU (A_M(1,W,M), A_M(1,S,M), A_M(1,B,M), A_M(1,E,M), A_M(1,N,M), A_M&
         (1,T,M), BL01N, BL03N, BL09N, BU01N, BU03N, BU09N, BD00N) 
!
!
!
!  OBTAIN Solution
!
      CALL IGCG (A_M(1,W,M), A_M(1,S,M), A_M(1,B,M), A_M(1,E,M), A_M(1,N,M), &
         A_M(1,T,M), BL01N, BL03N, BL09N, BU01N, BU03N, BU09N, BD00N, B_M(1,M)&
         , VAR, ITMAX, IER) 
!
!      IF(IER .NE. 0) THEN
!        LINE = "Warning: No convergence in LEq Solver"
!        CALL WRITE_ERROR('LEQ_IGCG: '//VName, LINE, 1)
!      ENDIF
!
      RETURN  
      END SUBROUTINE LEQ_IGCG 
