!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SOLVE_LIN_EQ(Vname, Var, A_m, B_m, M, ITMAX, METHOD,   C
!                            SWEEP, TOL,                IER)           C
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
      SUBROUTINE SOLVE_LIN_EQ(VNAME, VAR, A_M, B_M, M, ITMAX, METHOD,&
                              SWEEP, TOL, IER) 
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
!                      Variable
      DOUBLE PRECISION Var(DIMENSION_3)

!                           convergence tolerance for leq solver
      DOUBLE PRECISION  :: TOL  
      
!			sweep direction of leq solver
      CHARACTER*4	:: SWEEP
      
      CHARACTER*80     LINE0, LINE1
!-----------------------------------------------
!      DOUBLE PRECISION, PARAMETER :: TOL = 1.0D-4
      INTEGER, PARAMETER :: MAX_IT = 1
!
      double precision, allocatable, dimension(:,:) :: A_mt
      logical, parameter :: do_transpose = .true. 
      integer :: ii, ijk 

!
!
      SELECT CASE (METHOD)  
      CASE (1)  
!       Successive over relaxation method from Templates
        CALL LEQ_SOR (VNAME, VAR, A_M, B_M, M, ITMAX, IER) 
!
      CASE (2)  

      if(do_transpose) then	
        allocate( A_mt(-3:3, ijkstart3:ijkend3 ))
!$omp parallel do private(ijk,ii)
        do ijk=ijkstart3,ijkend3
          do ii=-3,3
               A_mt(ii,ijk) = A_m(ijk,ii,M)
          enddo
        enddo
      endif
!       bicgstab
      if(do_transpose) then
        call leq_bicgst( VNAME,VAR, A_Mt(:,:), B_M(:,M), SWEEP, TOL, ITMAX,IER)
      else 
        call leq_bicgs( VNAME,VAR, A_M(:,:,M), B_M(:,M), SWEEP, TOL, ITMAX,IER)
      endif

      if(do_transpose) then	
        deallocate( A_mt )
      endif

      case(3)
        call leq_gmres( VNAME, VAR, A_M(:,:,M), B_M(:,M), &
                        SWEEP, TOL, ITMAX, MAX_IT, IER )

      CASE (4)  
        IER = 0
        call leq_bicgs( VNAME,VAR, A_M(:,:,M), B_M(:,M), SWEEP, TOL, ITMAX,IER)
        if (IER .eq. -2) then
          IER = 0
          print*,'calling leq_gmres', Vname
          call leq_gmres( VNAME, VAR, A_M(:,:,M), B_M(:,M), &
                        SWEEP, TOL, ITMAX, MAX_IT, IER )
        endif
	
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
