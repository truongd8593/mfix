!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: LEQ_CG(Vname, Var, A_m, B_m,                           C
!                         cmethod, TOL, ITMAX, IER )                   C
!  Purpose: Compute residual of linear system                          C
!                                                                      C
!  Author: S. Pannala                                 Date: 18-JUL-07  C
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
      SUBROUTINE LEQ_CG(VNAME, VNO, VAR, A_M, B_m,  cmethod, TOL, PC, ITMAX,IER)
      
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE matrix
      USE geometry
      USE compar
      USE indices
      USE leqsol
      USE funits
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Error indicator
      INTEGER ::          IER
!                      maximum number of iterations
      INTEGER ::          ITMAX
!                      variable number
      INTEGER ::          VNO
!                      convergence tolerance
      DOUBLE PRECISION ::  TOL
!                      Preconditioner
      CHARACTER*4   ::  PC
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION, DIMENSION(ijkstart3:ijkend3,-3:3) :: A_m
!                      Vector b_m
      DOUBLE PRECISION, DIMENSION(ijkstart3:ijkend3) :: B_m
!                      Variable name
      CHARACTER*(*) ::    Vname
!                      Variable
      DOUBLE PRECISION, DIMENSION(ijkstart3:ijkend3) :: Var
!                    sweep direction
      CHARACTER*(*) :: CMETHOD
!
!-------------------------------------------------
      DOUBLE PRECISION DNRM2
      EXTERNAL LEQ_MATVEC, LEQ_MSOLVE, LEQ_MSOLVE0, LEQ_MSOLVE1


!--------------------------------------------------

      if(PC.eq.'LINE') then
         call LEQ_CG0( Vname, Vno, Var, A_m, B_m,                        &
         cmethod, TOL, ITMAX, LEQ_MATVEC, LEQ_MSOLVE, IER )
      elseif(PC.eq.'DIAG') then
         call LEQ_CG0( Vname, Vno, Var, A_m, B_m,                        &
         cmethod, TOL, ITMAX, LEQ_MATVEC, LEQ_MSOLVE1, IER )
      elseif(PC.eq.'NONE') then
         call LEQ_CG0( Vname, Vno, Var, A_m, B_m,                        &
         cmethod, TOL, ITMAX, LEQ_MATVEC, LEQ_MSOLVE0, IER )
      else
         IF(DMP_LOG)WRITE (UNIT_LOG,*) 'preconditioner option not found - check mfix.dat and readme'
         call mfix_exit(myPE)
      endif

      return
      END SUBROUTINE LEQ_CG




!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: LEQ_CG0(Vname, Var, A_m, B_m,                          C
!                         cmethod, TOL, ITMAX, MATVEC, MSOLVE, IER )   C
!  Purpose: Compute residual of linear system                          C
!                                                                      C
!  Author: S. Pannala                                 Date: 18-JUL-07  C
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
      SUBROUTINE LEQ_CG0(VNAME, VNO, VAR, A_M, B_m,  cmethod, TOL, ITMAX,  &
                            MATVEC, MSOLVE, IER ) 
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE matrix
      USE geometry
      USE compar
      USE mpi_utility
      USE sendrecv
      USE indices
      USE leqsol
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!
!                      Error indicator
      INTEGER ::          IER
!                      maximum number of iterations
      INTEGER ::          ITMAX
!                      variable number
      INTEGER ::          VNO
!                      convergence tolerance
      DOUBLE PRECISION ::  TOL
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION, DIMENSION(ijkstart3:ijkend3,-3:3) :: A_m
!                      Vector b_m
      DOUBLE PRECISION, DIMENSION(ijkstart3:ijkend3) :: B_m
!                      Variable name
      CHARACTER*(*) ::    Vname
!                      Variable
      DOUBLE PRECISION, DIMENSION(ijkstart3:ijkend3) :: Var
!                    sweep direction
      CHARACTER*(*) :: CMETHOD
!
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, PARAMETER :: idebugl = 1
      DOUBLE PRECISION :: ratiotol = 0.2

      DOUBLE PRECISION, DIMENSION(ijkstart3:ijkend3) :: Xinit, R,P,Zvec,Q
      DOUBLE PRECISION, DIMENSION(0:ITMAX+1) :: alpha,beta,rho
      DOUBLE PRECISION :: RxZ, PxQ, aijmax, oam, Rnorm=0, Rnorm0, TOLMIN, pnorm
      LOGICAL :: isconverged
      INTEGER :: i, ii, j, k, ijk, itemp, iter

!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
!     DOUBLE PRECISION , EXTERNAL :: DOT_PRODUCT_PAR
      EXTERNAL  MATVEC, MSOLVE

      INTERFACE
         DOUBLE PRECISION FUNCTION DOT_PRODUCT_PAR( R1, R2 )
         use compar
         DOUBLE PRECISION, INTENT(IN), DIMENSION(ijkstart3:ijkend3) :: R1,R2
         END FUNCTION DOT_PRODUCT_PAR
      END INTERFACE

      logical, parameter :: do_unit_scaling = .true.

!-----------------------------------------------
      INCLUDE 'function.inc'
!     
      is_serial = numPEs.eq.1.and.is_serial

      alpha(:)  = zero
      beta(:)   = zero
      rho(:)    = zero

!     
!     ---------------------------------------------
!     zero out R,Zvec, P and Q
!     ---------------------------------------------
      if (use_doloop) then

!$omp  parallel do private(ijk)
         do ijk=ijkstart3,ijkend3
            R(ijk) = zero
            Zvec(ijk) = zero
            P(ijk) = zero
            Q(ijk) = zero
         enddo

      else

         R(:) = zero
         Zvec(:) = zero
         P(:) = zero
         Q(:) = zero

      endif

      TOLMIN = EPSILON( one )

      if (do_unit_scaling) then
!     
!     Scale matrix to have unit diagonal
!     
!$omp parallel do private(ijk,i,j,k,oam,aijmax)
         do k = kstart2,kend2
            do i = istart2,iend2
               do j = jstart2,jend2

                  IJK = funijk(i,j,k)


!                 aijmax = maxval(abs(A_M(ijk,:)) )
              
!                 if(aijmax.ne.abs(A_M(ijk,0))) & 
!                 write(*,*) 'Not positive definite', k,i,j,(A_M(ijk,:))

                  OAM = one/A_M(ijk,0)
                  
                  A_M(IJK,:) = A_M(IJK,:)*OAM

                  B_M(IJK) = B_M(IJK)*OAM

               enddo
            enddo
         enddo
      endif

!     
!     Compute initial residual, assume initial guess in Var + some small random number
!     r = b - A*x : Line 1

!     call random_number(Xinit(:))

!     if (use_doloop) then

!$omp   parallel do private(ijk)
!        do ijk=ijkstart3,ijkend3
!           Xinit(ijk) = Var(ijk)*(ONE + (2.0d0*Xinit(ijk)-1.0d0)*1.0d-6)
!        enddo
!     else
!        Xinit(:) = Var(:)* (ONE + (2.0d0*Xinit(:)-1.0d0)*1.0d-6)
!     endif

!     Xinit(:) = Zero

      if (idebugl >= 1) then
         if(myPE.eq.0) print*,'leq_bicgs, initial: ', Vname,' resid ', Rnorm0
      endif

      call MATVEC( Vname, Var, A_M, R )


      if (use_doloop) then

!$omp   parallel do private(ijk)
         do ijk=ijkstart3,ijkend3
            R(ijk) = B_m(ijk) - R(ijk)
         enddo
      else
         R(:) = B_m(:) - R(:)
      endif

      if(is_serial) then
         Rnorm0 = zero
         if (use_doloop) then

!$omp          parallel do private(ijk) reduction(+:Rnorm0)
            do ijk=ijkstart3,ijkend3
               Rnorm0 = Rnorm0 + R(ijk)*R(ijk)
            enddo
         else
            Rnorm0 = dot_product(R,R)
         endif
         Rnorm0 = sqrt( Rnorm0 )
      else
         Rnorm0 = sqrt( dot_product_par( R, R ) )
      endif

      if (idebugl >= 1) then
         if(myPE.eq.0) print*,'leq_bicgs, initial: ', Vname,' resid ', Rnorm0
      endif
!     
!     Main loop : Line 2
!     
      iter = 1
      do i=1,itmax
!     
!     Solve M Zvec(:) = R(:) : Line 3
!     
         call MSOLVE( Vname, R, A_m, Zvec, CMETHOD)

!     
!     Solve Rho = RxZ : Line 4
!     
         if(is_serial) then
            if (use_doloop) then
               RxZ = zero
!$omp        parallel do private(ijk) reduction(+:RxZ)
               do ijk=ijkstart3,ijkend3
                  RxZ = RxZ + R(ijk) * Zvec(ijk)
               enddo
               rho(i-1) = RxZ
            else
               rho(i-1) = dot_product( R, Zvec )
            endif
         else
            rho(i-1) = dot_product_par( R, Zvec )
         endif                  ! is_serial


         if (rho(i-1) .eq. zero) then
            if(i /= 1)then
!     ------------
!     Method fails
!     ------------
               ier = -2
            else
!     ------------
!     converged.  residual is already zero
!     ------------
               ier = 0
            endif
            call send_recv(var,2)
            return
         endif                  ! rho(i-1).eq.0

         if (i .eq. 1) then

!
!     P_1 = Z_0 | Line 6
!

            if (use_doloop) then
!$omp        parallel do private(ijk)
               do ijk=ijkstart3,ijkend3
                  P(ijk) = Zvec(ijk)
               enddo
            else
               P(:) = Zvec(:)
            endif
         else
!
!     beta = rho(i-1)/rho(i-2) | Line 8
!     P = Z + beta*P | Line 9
!
            beta(i-1) = ( rho(i-1)/rho(i-2) )
            if (use_doloop) then
!$omp        parallel do private(ijk)
               do ijk=ijkstart3,ijkend3
                  P(ijk) = Zvec(ijk) + beta(i-1)* P(ijk)
               enddo
            else
               P(:) = Zvec(:) + beta(i-1)*P(:)
            endif
         endif                  ! i.eq.1

!     Q(:) = A*P(:) : Line 10
!     
         call MATVEC( Vname, P, A_m, Q )
         
         if(is_serial) then
            if (use_doloop) then
               PxQ = zero
!$omp         parallel do private(ijk) reduction(+:PxQ)
               do ijk=ijkstart3,ijkend3
                  PxQ = PxQ + P(ijk) * Q(ijk)
               enddo
            else
               PxQ = dot_product( P, Q )
            endif
         else
            PxQ = dot_product_par( P, Q )
         endif                  ! is_serial

!     alpha = rho/PxQ : Line 11

         alpha(i) = rho(i-1)/PxQ

!     x = x + alpha*p : Line 12
!     r = r - alpha*q : Line 13

         if (use_doloop) then
!$omp     parallel do private(ijk)
            do ijk=ijkstart3,ijkend3
               R(ijk) = R(ijk) - alpha(i) * Q(ijk)
               Var(ijk) = Var(ijk) + alpha(i) * P(ijk)
            enddo
         else
            R(:) = R(:) - alpha(i) * Q(:)
            Var(:) = Var(:) + alpha(i) * P(:)
         endif                  ! use_doloop

!     
!     Check norm of R(:); if small enough, Exit
!     
         if(is_serial) then
            if (use_doloop) then
               Rnorm = zero
!$omp       parallel do private(ijk) reduction(+:Rnorm)
               do ijk=ijkstart3,ijkend3
                  Rnorm = Rnorm + R(ijk) * R(ijk)
               enddo
            else
               Rnorm = dot_product( R, R )
            endif
            Rnorm = sqrt( Rnorm )
         else
            Rnorm = sqrt( dot_product_par( R, R ) )
         endif                  ! is_serial


         if (idebugl >= 1) then
            print*,'leq_bicgs, initial: ', Vname,' Vnorm ', Rnorm
         endif 
         

         if (idebugl.ge.1) then
            if (myPE.eq.PE_IO) then
               print*,'iter, Rnorm ', iter, Rnorm
               print*,'PxQ, rho(i-1) ', PxQ, rho(i-1)
               print*,'alpha(i), beta(i-1) ', alpha(i), beta(i-1)
            endif
         endif
!     
         isconverged = (Rnorm <= TOL*Rnorm0)

         if (isconverged) then
            iter_tot(vno) = iter_tot(vno) + iter + 1
            EXIT
         endif


!     Advance the iteration count
         iter = iter + 1

      enddo

      if (idebugl >= 1) then
         call MATVEC( Vname, Var, A_m, R )
         if (use_doloop) then
!$omp  parallel do private(ijk)
            do ijk=ijkstart3,ijkend3
               R(ijk) = R(ijk) - B_m(ijk)
            enddo
         else
            R(:) = R(:) - B_m(:)
         endif

         if(is_serial) then
            if (use_doloop) then
               Rnorm = zero
!$omp         parallel do private(ijk) reduction(+:Rnorm)
               do ijk=ijkstart3,ijkend3
                  Rnorm = Rnorm + R(ijk) * R(ijk)
               enddo
            else
               Rnorm = dot_product( R,R)
            endif
            Rnorm = sqrt( Rnorm )
         else
            Rnorm = sqrt( dot_product_par( R,R) )
         endif

         if(myPE.eq.0) print*,'leq_bicgs: final Rnorm ', Rnorm

         if(myPE.eq.0)  print*,'leq_bicgs ratio : ', Vname,' ',iter,     &
         ' L-2', Rnorm/Rnorm0
      endif 

      IER = 0
      if (.not.isconverged) then
         IER = -1
         if (real(Rnorm) >= ratiotol*real(Rnorm0)) then
            IER = -2
         endif
      endif

      call send_recv(var,2)
      
      return
      end subroutine LEQ_CG0
