!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SOLVE_LIN_EQ(Vname, Vno, Var, A_m, B_m, M, ITMAX, METHOD,   C
!                            SWEEP, TOL1, PC,            IER)           C
!  Purpose: Interface for linear equation solver                       C
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
      SUBROUTINE SOLVE_LIN_EQ(VNAME, Vno, VAR, A_M, B_M, M, ITMAX, METHOD,&
                              SWEEP, TOL1, PC, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE geometry
      USE compar  
      USE residual  
      USE toleranc  
      USE leqsol  
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
!                      Maximum number of iterations
      INTEGER          ITMAX
!
!                      Linear equation solver method
      INTEGER          METHOD
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
!                      Variable number
      INTEGER          Vno
!
!                      Adjust LEQ Tolerance flag
      LOGICAL, PARAMETER :: adjust_leq_tol = .FALSE.
      LOGICAL, PARAMETER :: leq_tol_scheme1 = .FALSE.
      DOUBLE PRECISION max_resid_local, tol_resid_max
!
!                      Variable
      DOUBLE PRECISION Var(DIMENSION_3)

!                      convergence tolerance for leq solver
      DOUBLE PRECISION  :: TOL1, TOL  
!                      Preconditioner
      CHARACTER*4       :: PC  
      
!                      sweep direction of leq solver
      CHARACTER*4	:: SWEEP
      
      CHARACTER*80     LINE0, LINE1
!-----------------------------------------------
!      DOUBLE PRECISION, PARAMETER :: TOL = 1.0D-4
      INTEGER, PARAMETER :: MAX_IT = 1
!
      double precision, allocatable, dimension(:,:) :: A_mt
      integer :: ii, ijk 

!     Adjusting the tolerances
!     
      if(adjust_leq_tol) then 
         max_resid_local = maxval(resid(:,M),1)
         tol_resid_max   = max(tol_resid, tol_resid_t, tol_resid_th, tol_resid_x)
         if(leq_tol_scheme1.and.resid(Vno,M).lt.1.0e-1) then
            if(Vno.le.5) then
               TOL = MAX(TOL1,TOL1*RESID(Vno,M)/TOL_RESID)
            elseif (Vno.eq.6) then 
               TOL = MAX(TOL1,TOL1*RESID(Vno,M)/TOL_RESID_T)
            elseif (Vno.eq.7) then 
               TOL = MAX(TOL1,TOL1*RESID(Vno,M)/TOL_RESID_Th)
            elseif (Vno.eq.10) then 
               TOL = MAX(TOL1,TOL1*RESID(Vno,M)/TOL_RESID_X)
            endif 
            Write(*,*) 'Adjusting LEQ_Tolerance', Vname, tol, resid(Vno,M)
         else if(max_resid_local.lt.1.0e-1) then
            TOL = MAX(TOL1,TOL1*max_resid_local/TOL_RESID_max)
            Write(*,*) 'Adjusting LEQ_Tolerance', Vname, tol, max_resid_local
         endif 
      else
        TOL = TOL1
      endif

!
      SELECT CASE (METHOD)  
      CASE (1)  
!       Successive over relaxation method from Templates
        CALL LEQ_SOR (VNAME, VNO, VAR, A_M, B_M, M, ITMAX, IER) 
!
      CASE (2)  

      if(do_transpose) then	
        allocate( A_mt(-3:3, ijkstart3:ijkend3 ))
!!!$omp parallel do private(ijk,ii)
        do ijk=ijkstart3,ijkend3
          do ii=-3,3
               A_mt(ii,ijk) = A_m(ijk,ii,M)
          enddo
        enddo
      endif
!       bicgstab

      if(do_transpose) then
        call leq_bicgst( VNAME,VNO, VAR, A_Mt(:,:), B_M(:,M), SWEEP, TOL, PC, ITMAX,IER)
      else 
        call leq_bicgs( VNAME,VNO, VAR, A_M(:,:,M), B_M(:,M), SWEEP, TOL, PC, ITMAX,IER)
      endif

      if(do_transpose) then	
        deallocate( A_mt )
      endif

      case(3)
        call leq_gmres( VNAME, VNO, VAR, A_M(:,:,M), B_M(:,M), &
                        SWEEP, TOL, ITMAX, MAX_IT, IER )

      CASE (4)  
        IER = 0
        call leq_bicgs( VNAME,VNO, VAR, A_M(:,:,M), B_M(:,M), SWEEP, TOL, PC, ITMAX,IER)
        if (IER .eq. -2) then
          IER = 0
          print*,'calling leq_gmres', Vname
          call leq_gmres( VNAME, VNO, VAR, A_M(:,:,M), B_M(:,M), &
                        SWEEP, TOL, ITMAX, MAX_IT, IER )
        endif

      CASE (5)  
        call leq_cg( VNAME,VNO, VAR, A_M(:,:,M), B_M(:,M), SWEEP, TOL, PC, ITMAX,IER)

!     CASE (6) - Disabled
!       CALL LEQ_LSOR( VNAME, VAR, A_M(:,:,M), B_M(:,M), ITMAX, IER)

      CASE DEFAULT 
        LINE0(1:14) = 'SOLVE_LIN_EQ: '
	LINE0(15:80)= VName
        WRITE(LINE1,'(A, I2, A)') &
             'Error: LEQ_METHOD = ', METHOD, ' is invalid'
        CALL WRITE_ERROR(LINE0, LINE1, 1)
        call mfix_exit(myPE)
      END SELECT 
      RETURN  
      END SUBROUTINE SOLVE_LIN_EQ 
