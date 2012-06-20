!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine LEQ_BICGS                                                C
!  Purpose: Solve system of linear system using BICGS method           C
!           Biconjugate gradients stabilized                           C
!                                                                      C
!  Author: Ed D'Azevedo                               Date: 21-JAN-99  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE LEQ_BICGS(VNAME, VNO, VAR, A_M, B_m, cmethod, &
                           TOL, PC, ITMAX, IER)
      
!-----------------------------------------------
! Modules
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
! Dummy arguments
!-----------------------------------------------
! variable name
      CHARACTER*(*), INTENT(IN) :: Vname
! variable number (not really used here; see calling subroutine)
      INTEGER, INTENT(IN) :: VNO
! variable 
!     e.g., pp_g, epp, rop_g, rop_s, u_g, u_s, v_g, v_s, w_g,
!     w_s, T_g, T_s, x_g, x_s, Theta_m, scalar, K_Turb_G, 
!     e_Turb_G
      DOUBLE PRECISION, DIMENSION(ijkstart3:ijkend3), INTENT(INOUT) :: Var
! Septadiagonal matrix A_m
      DOUBLE PRECISION, DIMENSION(ijkstart3:ijkend3,-3:3), INTENT(INOUT) :: A_m
! Vector b_m
      DOUBLE PRECISION, DIMENSION(ijkstart3:ijkend3), INTENT(INOUT) :: B_m
! Sweep direction of leq solver (leq_sweep)
!     e.g., options = 'isis', 'rsrs' (default), 'asas'
! Note: this setting only seems to matter when leq_pc='line'      
      CHARACTER*(*), INTENT(IN) :: CMETHOD
! convergence tolerance (generally leq_tol)
      DOUBLE PRECISION, INTENT(IN) :: TOL
! preconditioner (leq_pc)
!     options = 'line' (default), 'diag', 'none'      
      CHARACTER*4, INTENT(IN) ::  PC
! maximum number of iterations (generally leq_it)
      INTEGER, INTENT(IN) :: ITMAX      
! error indicator
      INTEGER, INTENT(INOUT) :: IER
!-------------------------------------------------
! Local Variables      
!-------------------------------------------------
!------------------------------------------------- 
! External subroutines 
!------------------------------------------------- 
! These procedures are effectively dummy arguments (procedures as
! arguments within the subroutine leq_bicgs0)
      EXTERNAL LEQ_MATVEC, LEQ_MSOLVE, LEQ_MSOLVE0, LEQ_MSOLVE1
!--------------------------------------------------


      if(PC.eq.'LINE') then   ! default
         call LEQ_BICGS0( Vname, Vno, Var, A_m, B_m,  &
            cmethod, TOL, ITMAX, LEQ_MATVEC, LEQ_MSOLVE, IER )
      elseif(PC.eq.'DIAG') then
         call LEQ_BICGS0( Vname, Vno, Var, A_m, B_m,   &
            cmethod, TOL, ITMAX, LEQ_MATVEC, LEQ_MSOLVE1, IER )
      elseif(PC.eq.'NONE') then
         call LEQ_BICGS0( Vname, Vno, Var, A_m, B_m,   &
            cmethod, TOL, ITMAX, LEQ_MATVEC, LEQ_MSOLVE0, IER )
      else
         IF(DMP_LOG)WRITE (UNIT_LOG,*) &
           'preconditioner option not found - check mfix.dat and readme'
         call mfix_exit(myPE)
      endif

      return
      END SUBROUTINE LEQ_BICGS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: LEQ_BICGS0                                              C
!  Purpose: Compute residual of linear system                          C
!                                                                      C
!  Author: Ed D'Azevedo                               Date: 21-JAN-99  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE LEQ_BICGS0(VNAME, VNO, VAR, A_M, B_m, cmethod, &
                            TOL, ITMAX, MATVEC, MSOLVE, IER )

!-----------------------------------------------
! Modules
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
! Dummy arguments/procedure
!-----------------------------------------------
! variable name
      CHARACTER*(*), INTENT(IN) :: Vname
! variable number (not really used here-see calling subroutine)
      INTEGER, INTENT(IN) :: VNO
! variable 
!     e.g., pp_g, epp, rop_g, rop_s, u_g, u_s, v_g, v_s, w_g,
!     w_s, T_g, T_s, x_g, x_s, Theta_m, scalar, K_Turb_G, 
!     e_Turb_G
      DOUBLE PRECISION, INTENT(INOUT) :: Var(ijkstart3:ijkend3)
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(ijkstart3:ijkend3,-3:3)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(ijkstart3:ijkend3)
! Sweep direction of leq solver (leq_sweep)
!     e.g., options = 'isis', 'rsrs' (default), 'asas'
      CHARACTER*(*), INTENT(IN) :: CMETHOD
! convergence tolerance (generally leq_tol)
      DOUBLE PRECISION, INTENT(IN) :: TOL
! maximum number of iterations (generally leq_it)
      INTEGER, INTENT(IN) :: ITMAX      
! error indicator
      INTEGER, INTENT(INOUT) :: IER
! dummy arguments/procedures set as indicated
!     matvec->leq_matvec      
! for preconditioner (leq_pc)
!    'line' msolve->leq_msolve  (default)
!    'diag' msolve->leq_msolve1
!    'none' msolve->leq_msolve0
      EXTERNAL  MATVEC, MSOLVE
!-----------------------------------------------
! Local parameters
!-----------------------------------------------
      INTEGER, PARAMETER :: idebugl = 0
      DOUBLE PRECISION, PARAMETER :: ratiotol = 0.2
      LOGICAL, PARAMETER :: do_unit_scaling = .true.      
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      DOUBLE PRECISION, DIMENSION(ijkstart3:ijkend3) ::   &
                        R, Rtilde, P, Phat, Svec, Shat, Tvec, V
      DOUBLE PRECISION, DIMENSION(0:ITMAX+1) :: &
                        alpha, beta, omega, rho
      DOUBLE PRECISION :: TxS, TxT, RtildexV, RtildexR, &
                          aijmax, oam
      DOUBLE PRECISION :: Rnorm, Rnorm0, Snorm, TOLMIN, pnorm
      LOGICAL :: isconverged
      INTEGER :: i, j, k, ijk
      INTEGER :: iter
      DOUBLE PRECISION, DIMENSION(2) :: TxS_TxT
!-----------------------------------------------
! External subroutines/functions
!-----------------------------------------------

      INTERFACE
         DOUBLE PRECISION FUNCTION DOT_PRODUCT_PAR( R1, R2 )
         use compar
         DOUBLE PRECISION, INTENT(IN), DIMENSION(ijkstart3:ijkend3) :: R1, R2
         END FUNCTION DOT_PRODUCT_PAR
      END INTERFACE

      INTERFACE
         FUNCTION DOT_PRODUCT_PAR2( R1, R2, R3, R4 )
         use compar
         DOUBLE PRECISION, INTENT(IN), DIMENSION(ijkstart3:ijkend3) :: &
                                       R1, R2, R3, R4
         DOUBLE PRECISION, DIMENSION(2) :: DOT_PRODUCT_PAR2
         END FUNCTION DOT_PRODUCT_PAR2
      END INTERFACE

!-----------------------------------------------
! Include statement functions      
!-----------------------------------------------
      INCLUDE 'function.inc'
!-----------------------------------------------

      is_serial = numPEs.eq.1.and.is_serial

! these scalars should not be necessary to initialize but done as failsafe      
      rnorm = ZERO
      rnorm0 = ZERO
      snorm = ZERO
      pnorm = ZERO

! initializing      
      alpha(:)  = zero
      beta(:)   = zero
      omega(:)  = zero
      rho(:)    = zero



! zero out R, Rtilde, P, Phat, Svec, Shat, Tvec, V
! --------------------------------
      if (use_doloop) then   ! mfix.dat keyword default=false
!$omp  parallel do private(ijk)
         do ijk=ijkstart3,ijkend3
            R(ijk) = zero
            Rtilde(ijk) = zero
            P(ijk) = zero
            Phat(ijk) = zero
            Svec(ijk) = zero
            Shat(ijk) = zero
            Tvec(ijk) = zero
            V(ijk) = zero
         enddo
      else
         R(:) = zero
         Rtilde(:) = zero
         P(:) = zero
         Phat(:) = zero
         Svec(:) = zero
         Shat(:) = zero
         Tvec(:) = zero
         V(:) = zero
      endif

      TOLMIN = EPSILON( one )

! Scale matrix to have unit diagonal
! ---------------------------------------------------------------->>>
      if (do_unit_scaling) then
!$omp parallel do private(ijk,i,j,k,oam,aijmax)
         do k = kstart2,kend2
            do i = istart2,iend2
               do j = jstart2,jend2
                  IJK = funijk(i,j,k)
                  aijmax = maxval(abs(A_M(ijk,:)) )
                  OAM = one/aijmax
                  A_M(IJK,:) = A_M(IJK,:)*OAM
                  B_M(IJK) = B_M(IJK)*OAM
               enddo
            enddo
         enddo
      endif
! ----------------------------------------------------------------<<<


! Compute initial residual (R = b-A*x) for Ax=b
!    assume initial guess in Var
!    rtilde = r
! ---------------------------------------------------------------->>>
      call MATVEC(Vname, Var, A_M, R)   ! returns R=A*Var

      if (use_doloop) then   ! mfix.dat keyword default=false
!!$omp   parallel do private(ijk)
         do ijk=ijkstart3,ijkend3
            R(ijk) = B_m(ijk) - R(ijk)
         enddo
      else
         R(:) = B_m(:) - R(:)
      endif

      if(is_serial) then
         Rnorm0 = zero
         if (use_doloop) then   ! mfix.dat keyword default=false
!!$omp          parallel do private(ijk) reduction(+:Rnorm0)
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

! determine an initial guess for the residual = residual + small random
! number (so it could be set to anything). note that since random_number 
! is used to supply the guess, this line could potentially be the source
! of small differences between runs.  the random number is shifted below
! between -1 and 1 and then scaled by factor 1.0D-6*Rnorm0 
      call random_number(Rtilde(:))

      if (use_doloop) then   ! mfix.dat keyword default=false
!!$omp   parallel do private(ijk)
         do ijk=ijkstart3,ijkend3
            Rtilde(ijk) = R(ijk) + (2.0d0*Rtilde(ijk)-1.0d0)*1.0d-6*Rnorm0
         enddo
      else
         Rtilde(:) = R(:) + (2.0d0*Rtilde(:)-1.0d0)*1.0d-6*Rnorm0
      endif

      if (idebugl >= 1) then
         if(myPE.eq.0) print*,'leq_bicgs, initial: ', Vname,' resid ', Rnorm0
      endif
! ----------------------------------------------------------------<<<


! Main loop
! ---------------------------------------------------------------->>>
      iter = 1
      do i=1,itmax

         if(is_serial) then
            if (use_doloop) then   ! mfix.dat keyword default=false
               RtildexR = zero
!!$omp        parallel do private(ijk) reduction(+:RtildexR)
               do ijk=ijkstart3,ijkend3
                  RtildexR = RtildexR + Rtilde(ijk) * R(ijk)
               enddo
               rho(i-1) = RtildexR
            else
               rho(i-1) = dot_product( Rtilde, R )
            endif
         else
            rho(i-1) = dot_product_par( Rtilde, R )
         endif ! is_serial

!         print*,'leq_bicgs, initial: ', Vname,' rho(i-1) ', rho(i-1)

         if (rho(i-1) .eq. zero) then
            if(i /= 1)then
! Method fails
! --------------------------------
!               print*, 'leq_bicgs,',Vname,': rho(i-1) == 0 '
               ier = -2
            else
! Method converged.  residual is already zero
! --------------------------------
               ier = 0
            endif
            call send_recv(var,2)
            return
         endif ! rho(i-1).eq.0

         if (i .eq. 1) then
            if (use_doloop) then
!!$omp        parallel do private(ijk)
               do ijk=ijkstart3,ijkend3
                  P(ijk) = R(ijk)
               enddo
            else
               P(:) = R(:)
            endif
         else
            beta(i-1) = ( rho(i-1)/rho(i-2) )*( alpha(i-1) / omega(i-1) )
            if (use_doloop) then
!!$omp        parallel do private(ijk)
               do ijk=ijkstart3,ijkend3
                  P(ijk) = R(ijk) + beta(i-1)*( P(ijk) - omega(i-1)*V(ijk) )
               enddo
            else
               P(:) = R(:) + beta(i-1)*( P(:) - omega(i-1)*V(:) )
            endif
         endif ! i.eq.1

     
! Solve A*Phat(:) = P(:)
! V(:) = A*Phat(:)  
! --------------------------------
         call MSOLVE(Vname, P, A_m, Phat, CMETHOD) ! returns Phat  
         call MATVEC(Vname, Phat, A_m, V)   ! returns V=A*Phat
         
         if(is_serial) then
            if (use_doloop) then
               RtildexV = zero
!!$omp         parallel do private(ijk) reduction(+:RtildexV)
               do ijk=ijkstart3,ijkend3
                  RtildexV = RtildexV + Rtilde(ijk) * V(ijk)
               enddo
            else
               RtildexV = dot_product( Rtilde, V )
            endif
         else
            RtildexV = dot_product_par( Rtilde, V )
         endif ! is_serial

!         print*,'leq_bicgs, initial: ', Vname,' RtildexV ', RtildexV

! compute alpha
! --------------------------------
         alpha(i) = rho(i-1) / RtildexV

! compute Svec
! --------------------------------
         if (use_doloop) then
!!$omp     parallel do private(ijk)
            do ijk=ijkstart3,ijkend3
               Svec(ijk) = R(ijk) - alpha(i) * V(ijk)
            enddo
         else
            Svec(:) = R(:) - alpha(i) * V(:)
         endif ! use_doloop



! Check norm of Svec(:); if small enough:
! set X(:) = X(:) + alpha(i)*Phat(:) and stop
! --------------------------------
         if(.not.minimize_dotproducts) then
            if(is_serial) then
               if (use_doloop) then
                  Snorm = zero
!!$omp       parallel do private(ijk) reduction(+:Snorm)
                  do ijk=ijkstart3,ijkend3
                     Snorm = Snorm + Svec(ijk) * Svec(ijk)
                  enddo
               else
                  Snorm = dot_product( Svec, Svec )
               endif
               Snorm = sqrt( Snorm )
            else
               Snorm = sqrt( dot_product_par( Svec, Svec ) )
            endif               ! is_serial
!            print*,'leq_bicgs, initial: ', Vname,' Snorm ', real(Snorm)


            if (Snorm <= TOLMIN) then
               if (use_doloop) then
!!$omp          parallel do private(ijk)
                  do ijk=ijkstart3,ijkend3
                     Var(ijk) = Var(ijk) + alpha(i)*Phat(ijk)
                  enddo
               else
                  Var(:) = Var(:) + alpha(i)*Phat(:)
               endif            ! use_doloop

! Recompute residual norm 
! --------------------------------
               if (idebugl >= 1) then
                  call MATVEC(Vname, Var, A_m, R)   ! returns R=A*Var
!                  Rnorm = sqrt( dot_product_par( Var, Var ) )
!                  print*,'leq_bicgs, initial: ', Vname,' Vnorm ', Rnorm

                  if (use_doloop) then
!!$omp          parallel do private(ijk)
                     do ijk=ijkstart3,ijkend3
                        R(ijk) = B_m(ijk) - R(ijk)
                     enddo
                  else
                     R(:) = B_m(:) - R(:)
                  endif

                  if(is_serial) then
                     if (use_doloop) then
                        Rnorm = zero
!!$omp            parallel do private(ijk) reduction(+:Rnorm)
                        do ijk=ijkstart3,ijkend3
                           Rnorm = Rnorm + R(ijk)*R(ijk)
                        enddo
                     else
                        Rnorm =  dot_product( R, R ) 
                     endif
                     Rnorm = sqrt( Rnorm )
                  else
                     Rnorm = sqrt( dot_product_par( R, R ) )
                  endif
!                  print*,'leq_bicgs, initial: ', Vname,' Rnorm ', Rnorm
               endif            ! idebugl >= 1

               EXIT
            endif               ! end if (Snorm <= TOLMIN)
         endif                  ! end if (.not.minimize_dotproducts)



! Solve A*Shat(:) = Svec(:)
! Tvec(:) = A*Shat(:)
! --------------------------------
         call MSOLVE( Vname, Svec, A_m, Shat, CMETHOD)  ! returns Shat
         call MATVEC( Vname, Shat, A_m, Tvec )   ! returns Tvec=A*Shat

         if(is_serial) then
            if (use_doloop) then
               TxS = zero
               TxT = zero
!!$omp  parallel do private(ijk) reduction(+:TxS,TxT)
               do ijk=ijkstart3,ijkend3
                  TxS = TxS + Tvec(ijk)  * Svec(ijk)
                  TxT = TxT + Tvec(ijk)  * Tvec(ijk)
               enddo
            else
               TxS = dot_product( Tvec, Svec )
               TxT = dot_product( Tvec, Tvec )
            endif
         else
            if(.not.minimize_dotproducts) then
               TxS = dot_product_par( Tvec, Svec )
               TxT = dot_product_par( Tvec, Tvec )
            else
               TxS_TxT = dot_product_par2(Tvec, Svec, Tvec, Tvec )
               TxS = TxS_TxT(1)
               TxT = TxS_TxT(2)
            endif
         endif

         IF(TxT.eq.Zero) TxT = SMALL_NUMBER

! compute omega
! --------------------------------
         omega(i) = TxS / TxT

! compute new guess for Var
! --------------------------------
         if (use_doloop) then
!!$omp    parallel do private(ijk)
            do ijk=ijkstart3,ijkend3
               Var(ijk) = Var(ijk) +  &
                  alpha(i)*Phat(ijk) + omega(i)*Shat(ijk)
                  R(ijk) = Svec(ijk) - omega(i)*Tvec(ijk)
            enddo
         else
            Var(:) = Var(:) +  &
               alpha(i)*Phat(:) + omega(i)*Shat(:)
               R(:) = Svec(:) - omega(i)*Tvec(:)
         endif


! --------------------------------
         if(.not.minimize_dotproducts.or.(mod(iter,icheck_bicgs).eq.0)) then
            if(is_serial) then
               if (use_doloop) then
                  Rnorm = zero
!!$omp       parallel do private(ijk) reduction(+:Rnorm)
                  do ijk=ijkstart3,ijkend3
                     Rnorm = Rnorm + R(ijk) * R(ijk)
                  enddo
               else
                  Rnorm =  dot_product(R, R )
               endif
               Rnorm = sqrt( Rnorm )
            else
               Rnorm = sqrt( dot_product_par(R, R) )
            endif               ! is_serial

            if (idebugl.ge.1) then
               if (myPE.eq.PE_IO) then
                  print*,'iter, Rnorm ', iter, Rnorm, Snorm
                  print*,'alpha(i), omega(i) ', alpha(i), omega(i)
                  print*,'TxS, TxT ', TxS, TxT
                  print*,'RtildexV, rho(i-1) ', RtildexV, rho(i-1)
               endif
            endif

!           call mfix_exit(myPE)

! Check convergence; continue if necessary
! for continuation, it is necessary that omega(i) .ne. 0
            isconverged = (Rnorm <= TOL*Rnorm0)

            if (isconverged) then
               iter_tot(vno) = iter_tot(vno) + iter + 1
               EXIT
            endif
         endif                  ! end if(.not.minimize_dotproducts)

! Advance the iteration count
         iter = iter + 1

      enddo   ! end do i=1,itmax
! end of linear solver loop 
! ----------------------------------------------------------------<<<


      if (idebugl >= 1) then
         call MATVEC(Vname, Var, A_m, R)   ! returns R=A*Var
         if (use_doloop) then
!!$omp  parallel do private(ijk)
            do ijk=ijkstart3,ijkend3
               R(ijk) = R(ijk) - B_m(ijk)
            enddo
         else
            R(:) = R(:) - B_m(:)
         endif

         if(is_serial) then
            if (use_doloop) then
               Rnorm = zero
!!$omp         parallel do private(ijk) reduction(+:Rnorm)
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

         if(myPE.eq.0)  print*,'leq_bicgs ratio : ', Vname,' ',iter,  &
         ' L-2', Rnorm/Rnorm0
      endif   ! end if(idebugl >=1)

      isconverged = (real(Rnorm) <= TOL*Rnorm0);
!     write(*,*) '***',iter, isconverged, Rnorm, TOL, Rnorm0, myPE
      IER = 0
      if (.not.isconverged) then
         IER = -1
         iter_tot(vno) = iter_tot(vno) + iter
         if (real(Rnorm) >= ratiotol*real(Rnorm0)) then
            IER = -2
         endif
      endif

      call send_recv(var,2)
      
      return
      end subroutine LEQ_BICGS0


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: LEQ_ISWEEP                                              C
!  Purpose: Perform line sweep at coordinate I                         C
!                                                                      C
!  Author: Ed D'Azevedo                               Date: 21-JAN-99  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE LEQ_ISWEEP(I, Vname, VAR, A_M, B_M)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE matrix
      USE geometry
      USE compar
      USE indices
      USE funits
      USE sendrecv
      USE mpi_utility
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
!  Line position
      INTEGER, INTENT(IN) :: I
! Variable name
      CHARACTER*(*), INTENT(IN) :: Vname
! Variable
      DOUBLE PRECISION, INTENT(INOUT) :: Var(ijkstart3:ijkend3)
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(IN) :: A_m(ijkstart3:ijkend3, -3:3)
! Vector b_m
      DOUBLE PRECISION, INTENT(IN) :: B_m(ijkstart3:ijkend3)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      DOUBLE PRECISION, DIMENSION (JSTART:JEND) :: CC, DD, EE, BB
      INTEGER :: NSTART, NEND, INFO       
      INTEGER :: IJK, J, K, IM1JK, IP1JK
!-----------------------------------------------
! Include statement functions
!-----------------------------------------------
      INCLUDE 'function.inc'
!-----------------------------------------------

      NEND = JEND
      NSTART = JSTART
      K = 1

      DO J=NSTART, NEND
!         IJK = FUNIJK(IMAP_C(I),JMAP_C(J),KMAP_C(K))
         IJK = FUNIJK(I,J,K)
         IM1JK = IM_OF(IJK)
         IP1JK = IP_OF(IJK)
         DD(J) = A_M(IJK,  0)
         CC(J) = A_M(IJK, -2)
         EE(J) = A_M(IJK,  2)
         BB(J) = B_M(IJK) -  A_M(IJK,-1) * Var( IM1JK )  &
                          -  A_M(IJK, 1) * Var( IP1JK )
      ENDDO

      CC(NSTART) = ZERO
      EE(NEND) = ZERO
      INFO = 0
!     CALL DGTSL( JEND-JSTART+1, CC, DD, EE, BB, INFO )
      CALL DGTSV( JEND-JSTART+1, 1, CC(JSTART+1), DD, EE, BB, JEND-JSTART+1, INFO )
      
      IF (INFO.NE.0) THEN
         RETURN
      ENDIF
      
      DO J=NSTART, NEND
         IJK = FUNIJK(I,J,K)
         Var(IJK) =  BB(J) 
      ENDDO

      RETURN
      END SUBROUTINE LEQ_ISWEEP


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: LEQ_IKSWEEP                                             C
!  Purpose: Perform line sweep at coordinate I, K                      C
!                                                                      C
!  Author: Ed D'Azevedo                               Date: 21-JAN-99  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE LEQ_IKSWEEP(I, K, Vname, VAR, A_M, B_M)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE matrix
      USE geometry
      USE compar
      USE funits
      USE indices
      USE sendrecv
      USE mpi_utility
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Line position
      INTEGER, INTENT(IN) :: I, K
! Variable name
      CHARACTER*(*), INTENT(IN) :: Vname
! Variable
      DOUBLE PRECISION, INTENT(INOUT) :: Var(ijkstart3:ijkend3)
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(IN) :: A_m(ijkstart3:ijkend3, -3:3)
! Vector b_m
      DOUBLE PRECISION, INTENT(IN) :: B_m(ijkstart3:ijkend3)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      DOUBLE PRECISION, DIMENSION(JSTART:JEND) :: CC, DD, EE, BB
      INTEGER :: NSTART, NEND, INFO 
      INTEGER :: IJK, J, IM1JK, IP1JK, IJKM1, IJKP1
!-----------------------------------------------
! Include statement functions
!-----------------------------------------------
      INCLUDE 'function.inc'
!-----------------------------------------------

      NEND = JEND
      NSTART = JSTART

!!$omp parallel do private(j,ijk,im1jk,ip1jk,ijkm1,ijkp1)
      DO J=NSTART, NEND
!         IJK = FUNIJK(IMAP_C(I),JMAP_C(J),KMAP_C(K))
         IJK = FUNIJK(I,J,K)
         IM1JK = IM_OF(IJK)
         IP1JK = IP_OF(IJK)
         IJKM1 = KM_OF(IJK)
         IJKP1 = KP_OF(IJK)
         DD(J) = A_M(IJK,  0)
         CC(J) = A_M(IJK, -2)
         EE(J) = A_M(IJK,  2)
         BB(J) = B_M(IJK) -  A_M(IJK,-1) * Var( IM1JK ) &
                          -  A_M(IJK, 1) * Var( IP1JK ) &
                          -  A_M(IJK,-3) * Var( IJKM1 ) &
                          -  A_M(IJK, 3) * Var( IJKP1 )
      ENDDO

      CC(NSTART) = ZERO
      EE(NEND) = ZERO
      INFO = 0
!     CALL DGTSL( JEND-JSTART+1, CC, DD, EE, BB, INFO )
      CALL DGTSV(NEND-NSTART+1, 1, CC(NSTART+1), DD, EE, BB, NEND-NSTART+1, INFO)
      
      IF (INFO.NE.0) THEN
         write(*,*) 'leq_iksweep',INFO, myPE
         IF(DMP_LOG)WRITE (UNIT_LOG,*) 'ROUTINE = ', ' IKSWEEP'
         RETURN
      ENDIF
      
      DO J=NSTART, NEND
         IJK = FUNIJK(I,J,K)
         Var(IJK) = BB(J)
      ENDDO
      
      RETURN
      END SUBROUTINE LEQ_IKSWEEP


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: LEQ_JKSWEEP                                             C
!  Purpose: Perform line sweep at coordinate I, K                      C
!                                                                      C
!  Author: Ed D'Azevedo                               Date: 21-JAN-99  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE LEQ_JKSWEEP(J, K, Vname, VAR, A_M, B_M)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE matrix
      USE geometry
      USE funits
      USE compar
      USE indices
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Line position
      INTEGER, INTENT(IN) :: J, K
! Variable name
      CHARACTER*(*), INTENT(IN) :: Vname
! Variable
      DOUBLE PRECISION, INTENT(INOUT) :: Var(ijkstart3:ijkend3)
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(IN) :: A_m(ijkstart3:ijkend3, -3:3)
! Vector b_m
      DOUBLE PRECISION, INTENT(IN) :: B_m(ijkstart3:ijkend3)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      DOUBLE PRECISION, DIMENSION (ISTART:IEND) :: CC, DD, EE, BB
      INTEGER :: NSTART, NEND, INFO, IJK, I 
!-----------------------------------------------
! Include statement functions
!-----------------------------------------------
      INCLUDE 'function.inc'
!-----------------------------------------------
      
      NEND = IEND
      NSTART = ISTART

      DO I=NSTART,NEND
         IJK = FUNIJK(I,J,K)
         DD(I) = A_M(IJK,  0)
         CC(I) = A_M(IJK, -1)
         EE(I) = A_M(IJK,  1)
         BB(I) = B_M(IJK)    -  A_M(IJK,-2) * Var( JM_OF(IJK) ) &
                             -  A_M(IJK, 2) * Var( JP_OF(IJK) ) &
                             -  A_M(IJK,-3) * Var( KM_OF(IJK) ) &
                             -  A_M(IJK, 3) * Var( KP_OF(IJK) )
      ENDDO

      CC(NSTART) = ZERO
      EE(NEND) = ZERO
      INFO = 0
      CALL DGTSV(NEND-NSTART+1, 1, CC(NSTART+1), DD, EE, BB, NEND-NSTART+1, INFO)

      IF (INFO.NE.0) THEN
         IF(DMP_LOG)WRITE (UNIT_LOG,*) 'VNAME = ', VNAME
         IF(DMP_LOG)WRITE (UNIT_LOG,*) 'ROUTINE = ', ' JKSWEEP'
         IF(DMP_LOG)WRITE (UNIT_LOG,*) 'DGTSV RETURNS INFO = ', INFO
         call mfix_exit(myPE)
      ENDIF

      DO I=NSTART, NEND
         IJK = FUNIJK(I,J,K)
         Var(IJK) = BB(I)
      ENDDO

      RETURN
      END SUBROUTINE LEQ_JKSWEEP


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: LEQ_IJSWEEP                                             C
!  Purpose: Perform line sweep at coordinate I, K                      C
!                                                                      C
!  Author: Ed D'Azevedo                               Date: 21-JAN-99  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE LEQ_IJSWEEP(I, J, Vname, VAR, A_M, B_M)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE matrix
      USE geometry
      USE funits
      USE compar
      USE indices
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Line position
      INTEGER, INTENT(IN) :: I, J
! Variable name
      CHARACTER*(*), INTENT(IN) :: Vname
! Variable
      DOUBLE PRECISION, INTENT(INOUT) :: Var(ijkstart3:ijkend3)
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(IN) :: A_m(ijkstart3:ijkend3, -3:3)
! Vector b_m
      DOUBLE PRECISION, INTENT(IN) :: B_m(ijkstart3:ijkend3)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      DOUBLE PRECISION, DIMENSION (KSTART:KEND) :: CC, DD, EE, BB
      INTEGER :: NEND, NSTART, INFO, IJK, K
!-----------------------------------------------
! Include statement functions
!-----------------------------------------------
      INCLUDE 'function.inc'
!-----------------------------------------------

      NEND = KEND
      NSTART = KSTART

      DO K=NSTART, NEND
         IJK = FUNIJK(I,J,K)
         DD(K) = A_M(IJK,  0)
         CC(K) = A_M(IJK, -3)
         EE(K) = A_M(IJK,  3)
         BB(K) = B_M(IJK)    -  A_M(IJK,-2) * Var( JM_OF(IJK) ) &
                             -  A_M(IJK, 2) * Var( JP_OF(IJK) ) &
                             -  A_M(IJK,-1) * Var( IM_OF(IJK) ) &
                             -  A_M(IJK, 1) * Var( IP_OF(IJK) )
      ENDDO

      CC(NSTART) = ZERO
      EE(NEND) = ZERO
      INFO = 0
      CALL DGTSV(NEND-NSTART+1, 1, CC(NSTART+1), DD, EE, BB, NEND-NSTART+1, INFO)

      IF (INFO.NE.0) THEN
         IF(DMP_LOG)WRITE (UNIT_LOG,*) 'VNAME = ', VNAME
         IF(DMP_LOG)WRITE (UNIT_LOG,*) 'ROUTINE = ', ' IJSWEEP'
         IF(DMP_LOG)WRITE (UNIT_LOG,*) 'DGTSV RETURNS INFO = ', INFO
         call mfix_exit(myPE)
      ENDIF

      DO K=NSTART, NEND
         IJK = FUNIJK(I,J,K)
         Var(IJK) = BB(K)
      ENDDO

      RETURN
      END SUBROUTINE LEQ_IJSWEEP


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: LEQ_MATVEC                                              C
!  Purpose: Compute matrix vector multiplication                       C
!           (for linear equation Ax=b compute Ax                       C
!                                                                      C
!  Author: Ed D'Azevedo                               Date: 21-JAN-99  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE LEQ_MATVEC(VNAME, VAR, A_M, Avar)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE matrix
      USE geometry
      USE compar
      USE indices
      USE sendrecv
      USE mpi_utility
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Variable name
      CHARACTER*(*), INTENT(IN) :: Vname
! Variable
      DOUBLE PRECISION, INTENT(IN) :: Var(ijkstart3:ijkend3)
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(IN) :: A_m(ijkstart3:ijkend3, -3:3)
! Vector AVar
      DOUBLE PRECISION, INTENT(OUT) :: AVar(ijkstart3:ijkend3)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Variable
      INTEGER :: I, J, K, IJK
      integer :: im1jk, ip1jk, ijm1k, ijp1k, ijkm1, ijkp1
!-----------------------------------------------
! Include statement functions
!-----------------------------------------------
      INCLUDE 'function.inc'
!-----------------------------------------------

      if (do_k) then
!$omp    parallel  do &
!$omp&   private(     &
!$omp&           ijk,i,j,k, &
!$omp&           im1jk,ip1jk,ijm1k,ijp1k,ijkm1,ijkp1) collapse (3)
         do k = kstart,kend
            do i = istart,iend
               do j = jstart,jend
                  IJK = funijk(i,j,k)
                  im1jk = im_of(ijk)
                  ip1jk = ip_of(ijk)
                  ijm1k = jm_of(ijk)
                  ijp1k = jp_of(ijk)
                  ijkm1 = km_of(ijk)
                  ijkp1 = kp_of(ijk)
                  AVar(ijk) =  A_m(ijk,-3) * Var(ijkm1)   &
                             + A_m(ijk,-2) * Var(ijm1k)   &
                             + A_m(ijk,-1) * Var(im1jk)   &
                             + A_m(ijk, 0) * Var(ijk)     &
                             + A_m(ijk, 1) * Var(ip1jk)   &
                             + A_m(ijk, 2) * Var(ijp1k)   &
                             + A_m(ijk, 3) * Var(ijkp1)
               enddo
            enddo
         enddo
      else
         k = 1
!$omp parallel do private(i,j,ijk,im1jk,ip1jk,ijm1k,ijp1k) collapse (2)
         do i = istart,iend
            do j = jstart,jend
               IJK = funijk(i,j,k)
               im1jk = im_of(ijk)
               ip1jk = ip_of(ijk)
               ijm1k = jm_of(ijk)
               ijp1k = jp_of(ijk)
               AVar(ijk) =  A_m(ijk,-2) * Var(ijm1k)   &
                          + A_m(ijk,-1) * Var(im1jk)   &
                          + A_m(ijk, 0) * Var(ijk)     &
                          + A_m(ijk, 1) * Var(ip1jk)   &
                          + A_m(ijk, 2) * Var(ijp1k)
            enddo
         enddo
      endif

      call send_recv(Avar,nlayers_bicgs)
      RETURN
      END SUBROUTINE LEQ_MATVEC


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: LEQ_MSOLVE                                              C
!  Purpose:                                                            C
!  Notes: if leq_method is biggs or cg then this subroutine is         C
!         invoked when leq_pc='line'. if leq_method is gmres then      C
!         this subroutine is invoked (leq_pc setting does not matter)  C
!                                                                      C
!  Author: Ed D'Azevedo                               Date: 21-JAN-99  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE LEQ_MSOLVE(VNAME, B_m, A_M, Var, CMETHOD)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE matrix
      USE geometry
      USE compar
      USE indices
      USE sendrecv
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Variable name
      CHARACTER*(*), INTENT(IN) :: Vname
! Vector b_m
      DOUBLE PRECISION, INTENT(IN) :: B_m(ijkstart3:ijkend3)
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(IN) :: A_m(ijkstart3:ijkend3, -3:3)
! Variable
      DOUBLE PRECISION, INTENT(INOUT) :: Var(ijkstart3:ijkend3)
! Sweep direction of leq solver (leq_sweep)
!     e.g., options = 'isis', 'rsrs' (default), 'asas'
      CHARACTER*4, INTENT(IN) :: CMETHOD
!-----------------------------------------------
! Local parameters
!-----------------------------------------------
      LOGICAL, PARAMETER :: USE_IKLOOP = .FALSE.
      LOGICAL, PARAMETER :: SETGUESS = .TRUE.
!-----------------------------------------------
! Local variables
!-----------------------------------------------
!     
      INTEGER :: ITER, NITER
      INTEGER :: IJK, I , J, K
      INTEGER :: I1, J1, K1, I2, J2, K2, IK, JK, IJ
      INTEGER :: ISIZE, JSIZE, KSIZE 
      INTEGER :: ICASE
      
!     CHARACTER*4, PARAMETER :: CMETHOD = 'II'
      CHARACTER :: CH
      LOGICAL :: DO_ISWEEP, DO_JSWEEP, DO_KSWEEP
      LOGICAL :: DO_SENDRECV, DO_REDBLACK, DO_ALL
!-----------------------------------------------
! Include statement functions
!-----------------------------------------------
      INCLUDE 'function.inc'
!-----------------------------------------------
!!$      double precision omp_start, omp_end
!!$      double precision omp_get_wtime	      
!       by Tingwen      
!!$      omp_start=omp_get_wtime()

      IF (SETGUESS) THEN
!$omp   parallel do private(i,j,k,ijk) collapse (3)
         do k = kstart3,kend3
            do i = istart3,iend3
               do j = jstart3,jend3
                  IJK = funijk(i,j,k)
                  VAR(IJK) = B_M(IJK)
               enddo
            enddo
         enddo
         call send_recv(var,nlayers_bicgs)
      ENDIF

      NITER = LEN( CMETHOD )

      DO ITER=1,NITER
     
! Perform sweeps
         CH = CMETHOD( ITER:ITER )
         DO_ISWEEP = (CH .EQ. 'I') .OR. (CH .EQ. 'i')
         DO_JSWEEP = (CH .EQ. 'J') .OR. (CH .EQ. 'j')  
         DO_KSWEEP = (CH .EQ. 'K') .OR. (CH .EQ. 'k') 
         DO_ALL = (CH .EQ. 'A') .OR. (CH .EQ. 'a')         
         DO_REDBLACK = (CH .EQ. 'R') .OR. (CH .EQ. 'r')
         DO_SENDRECV = (CH .EQ. 'S') .OR. (CH .EQ. 's')

         IF (NO_K) THEN   ! two dimensional
! 2D run no need to enable openmp parallel
            IF ( DO_ISWEEP ) THEN
!!$omp   parallel do private(I)   
               DO I=istart,iend,1
                  CALL LEQ_ISWEEP( I, Vname, Var, A_m, B_m )
               ENDDO        
            ENDIF

         ELSE   ! three dimensional


! do_all true only for leq_pc='asas'
! ---------------------------------------------------------------->>>
            IF(DO_ALL) THEN        ! redblack for all sweeps, not used by default            
! JK Loop
! --------------------------------
               j1 = jstart
               k1 = kstart
               j2 = jend
               k2 = kend
               jsize = j2-j1+1
               ksize = k2-k1+1
               DO icase = 1, 2
!$omp   parallel do private(K,J,JK)
                  DO JK=icase, ksize*jsize, 2
                     if (mod(jk,jsize).ne.0) then
                        k = int( jk/jsize ) + k1
                     else
                        k = int( jk/jsize ) + k1 -1
                     endif
                     j = (jk-1-(k-k1)*jsize) + j1 + mod(k,2)
                     if(j.gt.j2) j=j-j2 + j1 -1                  
                     CALL LEQ_JKSWEEP(J, K, Vname, Var, A_m, B_m)
                  ENDDO
               ENDDO
               call send_recv(var,nlayers_bicgs)

! IJ Loop
! --------------------------------
               i1 = istart
               j1 = jstart
               i2 = iend
               j2 = jend
               isize = i2-i1+1
               jsize = j2-j1+1
               DO icase = 1, 2
!$omp   parallel do private(J,I,IJ)
                  DO IJ=icase, jsize*isize, 2
                     if (mod(ij,isize).ne.0) then
                        j = int( ij/isize ) + j1
                     else
                        j = int( ij/isize ) + j1 -1
                     endif
                     i = (ij-1-(j-j1)*isize) + i1 + mod(j,2)
                     if(i.gt.i2) i=i-i2 + i1 -1  
                     CALL LEQ_IJSWEEP(I, J, Vname, Var, A_m, B_m)
                  ENDDO
               ENDDO
               call send_recv(var,nlayers_bicgs)

! IK Loop
! --------------------------------
               i1 = istart
               k1 = kstart
               i2 = iend
               k2 = kend
               isize = i2-i1+1
               ksize = k2-k1+1

               DO icase = 1, 2
!$omp   parallel do private(K,I,IK)
                  DO IK=icase, ksize*isize, 2
                     if (mod(ik,isize).ne.0) then
                        k = int( ik/isize ) + k1
                     else
                        k = int( ik/isize ) + k1 -1
                     endif
                     i = (ik-1-(k-k1)*isize) + i1 + mod(k,2)
                     if(i.gt.i2) i=i-i2 + i1 -1      
                     CALL LEQ_IKSWEEP(I, K, Vname, Var, A_m, B_m)
                  ENDDO
               ENDDO
            ENDIF ! end DO_ALL
! ----------------------------------------------------------------<<<


! do_redblack only true leq_pc='rsrs'
! ---------------------------------------------------------------->>>
            IF(DO_REDBLACK) THEN
               i1 = istart
               k1 = kstart
               i2 = iend
               k2 = kend
               isize = i2-i1+1
               ksize = k2-k1+1
               DO icase = 1, 2
!$omp   parallel do private(K,I,IK)
                  DO IK=icase, ksize*isize, 2
                     if (mod(ik,isize).ne.0) then
                        k = int( ik/isize ) + k1
                     else
                        k = int( ik/isize ) + k1 -1
                     endif
                     i = (ik-1-(k-k1)*isize) + i1 + mod(k,2)
                     if(i.gt.i2) i=i-i2 + i1 -1
                     CALL LEQ_IKSWEEP(I, K, Vname, Var, A_m, B_m)
                  ENDDO
               ENDDO
            ENDIF   ! end if(do_redblack)
! ----------------------------------------------------------------<<<

!  Not sure the purpose of us_ikloop
!  The SMP directives below need review                        !Tingwen Jan 2012
            IF(USE_IKLOOP) THEN  
! use_ikloop is currently hard-wired to false (so goto else branch)
! ---------------------------------------------------------------->>>
               i1 = istart
               k1 = kstart
               i2 = iend
               k2 = kend
               isize = i2-i1+1
               ksize = k2-k1+1
               IF (DO_ISWEEP) THEN
!!$omp   parallel do private(K,I,IK)
                  DO IK=1, ksize*isize
                     if (mod(ik,isize).ne.0) then
                        k = int( ik/isize ) + k1
                     else
                        k = int( ik/isize ) + k1 -1
                     endif
                     i = (ik-1-(k-k1)*isize) + i1
                     CALL LEQ_IKSWEEP(I, K, Vname, Var, A_m, B_m)
                  ENDDO
               ENDIF
               IF (DO_KSWEEP) THEN
!!$omp   parallel do private(K,I,IK)
                  DO IK=1, ksize*isize
                     if (mod(ik,ksize).ne.0) then
                        i = int( ik/ksize ) + i1
                     else
                        i = int( ik/ksize ) + i1 -1
                     endif
                     k = (ik-1-(i-i1)*ksize) + k1
                     CALL LEQ_IKSWEEP(I, K, Vname, Var, A_m, B_m)
                  ENDDO
               ENDIF
! ----------------------------------------------------------------<<<               
            ELSE   ! else branch of if(use_ikloop)
!  Not sure the purpose of us_ikloop
!  The SMP directives below need review                        !Tingwen Jan 2012
! ---------------------------------------------------------------->>>
               IF (DO_ISWEEP) THEN
!!$omp   parallel do private(K,I)
                  DO K=kstart,kend
                     DO I=istart,iend
                        CALL LEQ_IKSWEEP(I, K, Vname, Var, A_m, B_m)
                     ENDDO
                  ENDDO
               ENDIF
               IF (DO_KSWEEP) THEN
!!$omp   parallel do private(K,I)
                  DO I=istart,iend
                     DO K=kstart,kend
                        CALL LEQ_IKSWEEP(I, K, Vname, Var, A_m, B_m)
                     ENDDO
                  ENDDO
               ENDIF
            ENDIF   ! end if/else (use(ikloop)
! ----------------------------------------------------------------<<<

         ENDIF   ! end if/else (do_k)


! this is called for all settings of leq_pc
         IF (DO_SENDRECV) call send_recv(var,nlayers_bicgs)


      ENDDO   ! end do iter=1,niter
!!$      omp_end=omp_get_wtime()
!!$      write(*,*)'leq_msolve:',omp_end - omp_start	

      RETURN
      END SUBROUTINE LEQ_MSOLVE


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: LEQ_MSOLVE0                                             C
!  Notes: do nothing or no preconditioning (leq_pc='none')             C
!                                                                      C
!  Author: Ed D'Azevedo                               Date: 21-JAN-99  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C 

      SUBROUTINE LEQ_MSOLVE0(VNAME, B_m, A_M, Var, CMETHOD)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE matrix
      USE geometry
      USE compar
      USE indices
      USE sendrecv
      use parallel
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Variable name
      CHARACTER*(*), INTENT(IN) :: Vname
! Vector b_m
      DOUBLE PRECISION, INTENT(IN) :: B_m(ijkstart3:ijkend3)
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(IN) :: A_m(ijkstart3:ijkend3, -3:3)
! Variable
      DOUBLE PRECISION, INTENT(OUT) :: Var(ijkstart3:ijkend3)
! sweep direction
      CHARACTER*4, INTENT(IN) :: CMETHOD
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      integer :: ijk
!-----------------------------------------------

! do nothing or no preconditioning
      if (use_doloop) then   ! mfix.dat keyword default=false
!!$omp  parallel do private(ijk)
         do ijk=ijkstart3,ijkend3
            var(ijk) = b_m(ijk)
         enddo
      else
         var(:) = b_m(:)
      endif
      call send_recv(var,nlayers_bicgs)

      return
      end subroutine leq_msolve0


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: LEQ_MSOLVE1                                             C
!  Notes: diagonal scaling (leq_pc='diag')                             C
!                                                                      C
!  Author: Ed D'Azevedo                               Date: 21-JAN-99  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE LEQ_MSOLVE1(VNAME, B_m, A_M, Var, CMETHOD)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE matrix
      USE geometry
      USE compar
      USE indices
      USE sendrecv
      use parallel
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Variable name
      CHARACTER*(*), INTENT(IN) :: Vname
! Vector b_m
      DOUBLE PRECISION, INTENT(IN) :: B_m(ijkstart3:ijkend3)
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(IN) :: A_m(ijkstart3:ijkend3, -3:3)
! Variable
      DOUBLE PRECISION, INTENT(OUT) :: Var(ijkstart3:ijkend3)
! sweep direction
      CHARACTER*4, INTENT(IN) :: CMETHOD
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      integer :: i,j,k, ijk 
!-----------------------------------------------
! Include statement functions
!-----------------------------------------------
      include 'function.inc'
!-----------------------------------------------

      if (use_doloop) then   ! mfix.dat keyword default=false
!!$omp    parallel do private(ijk)
         do ijk=ijkstart3,ijkend3
            var(ijk) = zero
         enddo
      else
         var(:) = ZERO
      endif

! diagonal scaling
!$omp   parallel do private(i,j,k,ijk) collapse (3)
      do k=kstart2,kend2
         do i=istart2,iend2
            do j=jstart2,jend2
               ijk = funijk( i,j,k )
               var(ijk) = b_m(ijk)/A_m(ijk,0)
            enddo
         enddo
      enddo

      call send_recv(var,nlayers_bicgs)

      return
      end subroutine leq_msolve1


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      double precision function dot_product_par(r1,r2)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      use mpi_utility
      use geometry
      use compar
      use indices
      implicit none
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
      double precision, intent(in), dimension(ijkstart3:ijkend3) :: r1,r2
!-----------------------------------------------
! Local parameters
!-----------------------------------------------
      logical, parameter :: do_global_sum = .true.
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      DOUBLE PRECISION, allocatable, Dimension(:) :: r1_g, r2_g
      double precision :: prod
      integer :: i, j, k, ijk
!-----------------------------------------------
! Include statement functions
!-----------------------------------------------
      include 'function.inc'
!-----------------------------------------------

      if(do_global_sum) then
         prod = 0.0d0
!$omp parallel do private(i,j,k,ijk) reduction(+:prod) collapse (3)
         do k = kstart1, kend1
            do i = istart1, iend1
               do j = jstart1, jend1
                  ijk = funijk (imap_c(i),jmap_c(j),kmap_c(k))
!                  ijk = funijk (i,j,k)
                  prod = prod + r1(ijk)*r2(ijk)
               enddo
            enddo
         enddo
         call global_all_sum(prod, dot_product_par)

      else
         if(myPE.eq.root) then
            allocate (r1_g(1:ijkmax3))
            allocate (r2_g(1:ijkmax3))
         else
            allocate (r1_g(10))
            allocate (r2_g(10))
         endif
         call gather(r1,r1_g)
         call gather(r2,r2_g)
        
         if(myPE.eq.root) then
            prod = 0.0d0
            
!$omp parallel do private(i,j,k,ijk) reduction(+:prod) collapse (3)
            do k = kmin1, kmax1
               do i = imin1, imax1
                  do j = jmin1, jmax1
                     ijk = funijk_gl (imap_c(i),jmap_c(j),kmap_c(k))
!                     ijk = funijk_gl (i,j,k)
                     prod = prod + r1_g(ijk)*r2_g(ijk)
                  enddo
               enddo
            enddo
         endif
         call bcast( prod)
         
         dot_product_par = prod

         deallocate (r1_g)
         deallocate (r2_g)

      endif
      
      end function dot_product_par


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      function dot_product_par2(r1,r2,r3,r4)

!-----------------------------------------------
! Modules
!-----------------------------------------------      
      use mpi_utility
      use geometry
      use compar
      use indices
      implicit none
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
      double precision, intent(in), dimension(ijkstart3:ijkend3) :: r1,r2,r3,r4
!-----------------------------------------------
! Local parameters
!-----------------------------------------------
      logical, parameter :: do_global_sum = .true.
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      DOUBLE PRECISION, allocatable, Dimension(:,:) :: r_temp, rg_temp
      double precision, Dimension(2) :: prod, dot_product_par2
      integer :: i, j, k, ijk
!-----------------------------------------------
! Include statement functions
!-----------------------------------------------
      include 'function.inc'
!-----------------------------------------------

      if(do_global_sum) then
         
         prod(:) = 0.0d0
         
!$omp parallel do private(i,j,k,ijk) reduction(+:prod) collapse (3)
         do k = kstart1, kend1
            do i = istart1, iend1
               do j = jstart1, jend1
                  ijk = funijk (imap_c(i),jmap_c(j),kmap_c(k))
!                  ijk = funijk (i,j,k)
                  prod(1) = prod(1) + r1(ijk)*r2(ijk)
                  prod(2) = prod(2) + r3(ijk)*r4(ijk)
               enddo
            enddo
         enddo
         call global_all_sum(prod, dot_product_par2)

      else
         allocate (r_temp(ijkstart3:ijkend3,4))
         r_temp(:,1) = r1
         r_temp(:,2) = r2
         r_temp(:,3) = r3
         r_temp(:,4) = r4

         if(myPE.eq.root) then
            allocate (rg_temp(1:ijkmax3,4))
         else
            allocate (rg_temp(10,4))
         endif
         call gather(r_temp,rg_temp)
         
         if(myPE.eq.root) then
            prod = 0.0d0
!$omp parallel do private(i,j,k,ijk) reduction(+:prod) collapse (3)
            do k = kmin1, kmax1
               do i = imin1, imax1
                  do j = jmin1, jmax1
                     ijk = funijk_gl (imap_c(i),jmap_c(j),kmap_c(k))
!                     ijk = funijk_gl (i,j,k)
                     prod(1) = prod(1) + rg_temp(ijk,1)*rg_temp(ijk,2)
                     prod(2) = prod(2) + rg_temp(ijk,3)*rg_temp(ijk,4)
                  enddo
               enddo
            enddo
         endif
         call bcast( prod)
         
         dot_product_par2 = prod

         deallocate (r_temp)
         deallocate (rg_temp)

      endif
      
      end function dot_product_par2

    

