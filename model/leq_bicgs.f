!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: LEQ_BICGS(Vname, Var, A_m, B_m,                        C
!                         cmethod, TOL, ITMAX, IER )
!  Purpose: Compute residual of linear system                          C
!                                                                      C
!  Author: Ed D'Azevedo                               Date: 21-JAN-99  C
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
      SUBROUTINE LEQ_BICGS(VNAME, VAR, A_M, B_m,  cmethod, TOL, ITMAX,IER)
      
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE matrix
      USE geometry
      USE compar
      USE indices
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Error indicator
      INTEGER, INTENT(INOUT) ::          IER
!                      maximum number of iterations
      INTEGER, INTENT(IN) ::          ITMAX
!                      convergence tolerance
      DOUBLE PRECISION, INTENT(IN) ::  TOL
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(IN), DIMENSION(:,:) :: A_m
!                      Vector b_m
      DOUBLE PRECISION, INTENT(IN), DIMENSION(:) :: B_m
!                      Variable name
      CHARACTER*(*), INTENT(IN) ::    Vname
!                      Variable
      DOUBLE PRECISION, INTENT(INOUT), DIMENSION(:) :: Var
!                    sweep direction
      CHARACTER*(*), INTENT(INOUT) :: CMETHOD
!
      INTEGER          IER
      INTEGER          ITMAX
!
!-------------------------------------------------
      DOUBLE PRECISION DNRM2
      EXTERNAL LEQ_MATVEC, LEQ_MSOLVE


!--------------------------------------------------


      cmethod = 'RSRS'
      call LEQ_BICGS0( Vname, Var, A_m, B_m,                        &
                       cmethod, TOL, ITMAX, LEQ_MATVEC, LEQ_MSOLVE, IER )

      return
      END SUBROUTINE LEQ_BICGS




!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: LEQ_BICGS0(Vname, Var, A_m, B_m,                       C
!                         cmethod, TOL, ITMAX, MATVEC, MSOLVE, IER )   C
!  Purpose: Compute residual of linear system                          C
!                                                                      C
!  Author: Ed D'Azevedo                               Date: 21-JAN-99  C
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
      SUBROUTINE LEQ_BICGS0(VNAME, VAR, A_M, B_m,  cmethod, TOL, ITMAX,  &
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
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!
!                      Error indicator
      INTEGER          IER
!
!                      maximum number of iterations
      INTEGER          ITMAX
!
!
!                      convergence tolerance
      DOUBLE PRECISION TOL
!
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(ijkstart3:ijkend3, -3:3)
!
!                      Vector b_m
      DOUBLE PRECISION B_m(ijkstart3:ijkend3)
!
!                      Variable name
      CHARACTER*(*)    Vname
!
!                      Variable
      DOUBLE PRECISION Var(ijkstart3:ijkend3)

!                    sweep direction
      CHARACTER*4 :: CMETHOD

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, PARAMETER :: idebugl = 0
      DOUBLE PRECISION :: ratiotol = 0.2

      DOUBLE PRECISION, DIMENSION(ijkstart3:ijkend3) ::                       &
                                R,Rtilde, P,Phat, Svec, Shat, Tvec,V
      DOUBLE PRECISION, DIMENSION(0:ITMAX+1) :: alpha,beta,omega,rho
      DOUBLE PRECISION :: TxS, TxT, oam,RtildexV,                   &
		      RtildexR, aijmax, Rnorm, Rnorm0, Snorm, TOLMIN, pnorm
      LOGICAL :: isconverged
      INTEGER :: i, ii, j, k, ijk, itemp, iter

!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
!     DOUBLE PRECISION , EXTERNAL :: DOT_PRODUCT_PAR
      EXTERNAL  MATVEC, MSOLVE

        INTERFACE
        DOUBLE PRECISION FUNCTION DOT_PRODUCT_PAR( R1, R2 )
        DOUBLE PRECISION, INTENT(IN), DIMENSION(:) :: R1,R2
        END FUNCTION DOT_PRODUCT_PAR
        END INTERFACE

!//AIKE 0201
        logical, parameter :: do_unit_scaling = .true.
!       logical, parameter :: do_unit_scaling = .false.

!-----------------------------------------------
      INCLUDE 'function.inc'


      alpha(:)  = zero
      beta(:)   = zero
      omega(:)  = zero
      rho(:)    = zero

      TOLMIN = EPSILON( one )

      if (do_unit_scaling) then
!
!     Scale matrix to have unit diagonal
!
!$omp parallel do private(ijk,oam,aijmax)
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

!
!    Compute initial residual, assume initial guess in Var
!    r = b - A*x
!    rtilde = r


    call MATVEC( Vname, Var, A_M, R )

!      call out_array(var,'var')
!      call out_array(r,'r')

	R(ijkstart3:ijkend3) = B_m(ijkstart3:ijkend3) - R(ijkstart3:ijkend3)
!	call send_recv(R,2)

	Rtilde(:) = R(:)
	
	Rnorm0 = sqrt( dot_product_par( R, R ) )

!      print*,'leq_bicgs, initial: ', Vname,' resid ', real(Rnorm0)
    if (idebugl >= 1) then
       print*,'leq_bicgs, initial: ', Vname,' resid ', real(Rnorm0)
    endif
!
!   Main loop
!
    iter = 1
    do i=1,itmax

        rho(i-1) = dot_product_par( Rtilde, R )
!      print*,'leq_bicgs, initial: ', Vname,' rho(i-1) ', real(rho(i-1))

        if (rho(i-1) .eq. zero) then
	  if(i /= 1)then
!           ------------
!           Method fails
!           ------------
            !print*, 'leq_bicgs,',Vname,': rho(i-1) == 0 '
	    ier = -2
	  else
!           ------------
!           converged.  residual is already zero
!           ------------
            ier = 0
	  endif

          call send_recv(var,2)
	  return
        endif

        if (i .eq. 1) then

	   P(:) = R(:)

        else
           beta(i-1) = ( rho(i-1)/rho(i-2) )*( alpha(i-1) / omega(i-1) )

           P(:) = R(:) + beta(i-1)*( P(:) - omega(i-1)*V(:) )


        endif

!
!       Solve M Phat(:) = P(:)
!       V(:) = A*Phat(:)
!       

        call MSOLVE( Vname, P, A_m, Phat, CMETHOD)

        call MATVEC( Vname, Phat, A_m, V )

!      call out_array(p,'p')
!      call out_array(phat,'ph')
!      call out_array(v,'v')
!     call write_ab_m(a_m, b_m, ijkmax2, 0, ier)

            
        RtildexV = dot_product_par( Rtilde, V )
!      print*,'leq_bicgs, initial: ', Vname,' RtildexV ', real(RtildexV)

        alpha(i) = rho(i-1) / RtildexV

        Svec(:) = R(:) - alpha(i) * V(:)


!
!       Check norm of Svec(:); if small enough:
!       set X(:) = X(:) + alpha(i)*Phat(:) and stop
!
        Snorm = sqrt( dot_product_par( Svec, Svec ) )
!      print*,'leq_bicgs, initial: ', Vname,' Snorm ', real(Snorm)

        if (Snorm <= TOLMIN) then

!	write(91,*) alpha(i)
!	write(91,*) Var(ijkstart3:ijkend3)
!	write(91,*) Phat(ijkstart3:ijkend3)
!	write(*,*) '**************************'

             Var(ijkstart3:ijkend3) = Var(ijkstart3:ijkend3) + alpha(i)*Phat(ijkstart3:ijkend3)
!            if (idebugl >= 1) then
!
!               Recompute residual norm
!          
             call MATVEC( Vname, Var, A_m, R )
             Rnorm = sqrt( dot_product_par( Var, Var ) )
!      print*,'leq_bicgs, initial: ', Vname,' Vnorm ', real(Rnorm)
             R(ijkstart3:ijkend3) = B_m(ijkstart3:ijkend3) - R(ijkstart3:ijkend3)
!            endif
             Rnorm = sqrt( dot_product_par( R, R ) )
!      print*,'leq_bicgs, initial: ', Vname,' Rnorm ', real(Rnorm)

            EXIT
        endif

!
!       Solve M Shat(:) = Svec(:)
!       Tvec(:) = A * Shat(:)
!
        call MSOLVE( Vname, Svec, A_m, Shat, CMETHOD)

        
        call MATVEC( Vname, Shat, A_m, Tvec )

!       do k = kmin2, kmax2
!       do j = jmin2, jmax2
!       do ii = imin2, imax2

!       IJK = funijk(ii,j,k)
!       write(90,*) 'initial', ii,j,k,Svec(ijk),Tvec(ijk)

!       enddo
!       enddo
!       enddo


        TxS = dot_product_par( Tvec(:), Svec(:) )
        TxT = dot_product_par( Tvec(:), Tvec(:) )
        omega(i) = TxS / TxT


        Var(ijkstart3:ijkend3) = Var(ijkstart3:ijkend3) +                           &
              alpha(i)*Phat(ijkstart3:ijkend3) + omega(i)*Shat(ijkstart3:ijkend3)


        R(:) = Svec(:) - omega(i)*Tvec(:)

        Rnorm = sqrt( dot_product_par(R, R) )

        if (idebugl.ge.1) then
           if (myPE.eq.PE_IO) then
             print*,'iter, Rnorm ', iter, Rnorm, Snorm
             print*,'alpha(i), omega(i) ', alpha(i), omega(i)
             print*,'TxS, TxT ', TxS, TxT
             print*,'RtildexV, rho(i-1) ', RtildexV, rho(i-1)
           endif
        endif

!       Check convergence; continue if necessary
!       for continuation, it is necessary that omega(i) .ne. 0
!

        isconverged = (Rnorm <= TOL*Rnorm0)

        if (isconverged) then
            EXIT
        endif

        iter = iter + 1
        enddo


        if (idebugl >= 1) then

          call MATVEC( Vname, Var, A_m, R )
          R(:) = R(:) - B_m(:)
          Rnorm = sqrt( dot_product_par( R,R) )

          print*,'leq_bicgs: final Rnorm ', Rnorm

           print*,'leq_bicgs ratio : ', Vname,' ',iter,     &
                   ' L-2', real(Rnorm/Rnorm0)
        endif 

        isconverged = (real(Rnorm) <= TOL*Rnorm0);
!	write(*,*) '***',iter, isconverged, Rnorm, TOL, Rnorm0, myPE
        IER = 0
        if (.not.isconverged) then
            IER = -1
            if (real(Rnorm) >= ratiotol*real(Rnorm0)) then
                  IER = -2
            endif
        endif

        call send_recv(var,2)
        
        return
        end subroutine LEQ_BICGS0


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: LEQ_ISWEEP(I, Vname, Var, A_m, B_m )                   C
!  Purpose: Perform line sweep at coordiante I                         C
!                                                                      C
!                                                                      C
!  Author: Ed D'Azevedo                               Date: 21-JAN-99  C
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
      SUBROUTINE LEQ_ISWEEP(I,Vname, VAR, A_M, B_M)

!-----------------------------------------------
!   M o d u l e s
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
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!                      Line position
      INTEGER          I
!
!
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(ijkstart3:ijkend3, -3:3)
!
!                      Vector b_m
      DOUBLE PRECISION B_m(ijkstart3:ijkend3)
!
!                      Variable name
      CHARACTER*(*)    Vname

!
!                      Variable
      DOUBLE PRECISION Var(ijkstart3:ijkend3)

!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!


!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

      INTEGER :: NSTART, NEND 
      DOUBLE PRECISION, DIMENSION (JSTART:JEND) :: CC,DD,EE,BB
      INTEGER :: INFO, IJK, J, K, IM1JK, IP1JK

      INCLUDE 'function.inc'

      NEND = JEND
      NSTART = JSTART
      K = 1

      DO J=NSTART, NEND
         IJK = FUNIJK(IMAP_C(I),JMAP_C(J),KMAP_C(K))
         IM1JK = IM_OF(IJK)
         IP1JK = IP_OF(IJK)

         DD(J) = A_M(IJK,  0)
         CC(J) = A_M(IJK, -2)
         EE(J) = A_M(IJK,  2)
         BB(J) = B_M(IJK) -  A_M(IJK,-1) * Var( IM1JK )         &
                          -  A_M(IJK, 1) * Var( IP1JK )

     ENDDO

     CC(NSTART) = ZERO
     EE(NEND) = ZERO
     INFO = 0

!    CALL DGTSL( JEND-JSTART+1, CC, DD, EE, BB, INFO )
     CALL DGTSV( JEND-JSTART+1, 1, CC(JSTART+1), DD, EE, BB,  JEND-JSTART+1, INFO )

     IF (INFO.NE.0) THEN
	RETURN
     ENDIF

     DO J=NSTART, NEND
        IJK = FUNIJK(I,J,K)
        Var(IJK) =  BB(J) 
     ENDDO

     RETURN
     END SUBROUTINE  LEQ_ISWEEP
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: LEQ_IKSWEEP(I, K, Vname, Var, A_m, B_m )               C
!  Purpose: Perform line sweep at coordiante I,K                       C
!                                                                      C
!                                                                      C
!  Author: Ed D'Azevedo                               Date: 21-JAN-99  C
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
      SUBROUTINE LEQ_IKSWEEP(I,K,Vname, VAR, A_M, B_M )

!-----------------------------------------------
!   M o d u l e s
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
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!                      Line position
      INTEGER          I,K
!
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(ijkstart3:ijkend3, -3:3)
!
!                      Vector b_m
      DOUBLE PRECISION B_m(ijkstart3:ijkend3)
!
!                      Variable name
      CHARACTER*(*)    Vname

!
!                      Variable
      DOUBLE PRECISION Var(ijkstart3:ijkend3)
!

!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!


!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

      DOUBLE PRECISION, DIMENSION (JSTART:JEND) :: CC,DD,EE, BB
      INTEGER :: NSTART, NEND, INFO, IJK, J,  IM1JK, IP1JK, IJKM1, IJKP1

      INCLUDE 'function.inc'

      NEND = JEND
      NSTART = JSTART

!!$omp parallel do private(j,ijk,im1jk,ip1jk,ijkm1,ijkp1)
      DO J=NSTART, NEND

         IJK = FUNIJK(IMAP_C(I),JMAP_C(J),KMAP_C(K))
         IM1JK = IM_OF(IJK)
         IP1JK = IP_OF(IJK)
         IJKM1 = KM_OF(IJK)
         IJKP1 = KP_OF(IJK)

         DD(J) = A_M(IJK,  0)
         CC(J) = A_M(IJK, -2)
         EE(J) = A_M(IJK,  2)
         BB(J) = B_M(IJK) -  A_M(IJK,-1) * Var( IM1JK )         &
                             -  A_M(IJK, 1) * Var( IP1JK )         &
                             -  A_M(IJK,-3) * Var( IJKM1 )         &
                             -  A_M(IJK, 3) * Var( IJKP1 )

     ENDDO

     CC(NSTART) = ZERO
     EE(NEND) = ZERO
     INFO = 0
!    CALL DGTSL( JEND-JSTART+1, CC, DD, EE, BB, INFO )
     CALL DGTSV( JEND-JSTART+1, 1, CC(JSTART+1), DD, EE, BB,  JEND-JSTART+1, INFO )

     IF (INFO.NE.0) THEN
        write(*,*) 'leq_iksweep',INFO, myPE
        RETURN
     ENDIF

      DO J=NSTART, NEND

        IJK = FUNIJK(I,J,K)
        Var(IJK) = BB(J)

     ENDDO

     RETURN
     END SUBROUTINE  LEQ_IKSWEEP


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: LEQ_MATVEC(Vname, Var, A_m, B_m )                      C
!  Purpose: Compute residual of linear system                          C
!                                                                      C
!  Author: Ed D'Azevedo                               Date: 21-JAN-99  C
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
      SUBROUTINE LEQ_MATVEC(VNAME, VAR, A_M, Avar )
!
!-----------------------------------------------
!   M o d u l e s
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
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(ijkstart3:ijkend3, -3:3)
!
!                      Vector AVar
      DOUBLE PRECISION AVar(ijkstart3:ijkend3)
!
!                      Variable name
      CHARACTER*(*)    Vname
!
!                      Variable
      DOUBLE PRECISION Var(ijkstart3:ijkend3)
!                      Variable
      DOUBLE PRECISION, allocatable, Dimension(:) :: Var_g
!
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

!
      INTEGER          I,  J, K, IJK, ITER 
      INTEGER          II,  JJ, KK
      DOUBLE PRECISION oAm

      integer :: im1jk,ip1jk, ijm1k,ijp1k, ijkm1,ijkp1
      logical, parameter :: use_send_recv = .true.
      logical, parameter :: need_distribute_Avar = .true.
      logical, parameter :: use_funijk = .true.

      integer :: i1,i2, j1,j2, k1,k2, isize,jsize
      integer :: cc, ci, cj, ck

      integer :: ijk_c, im1jk_c,ip1jk_c, ijm1k_c,ijp1k_c, ijkm1_c,ijkp1_c


!-----------------------------------------------
      INCLUDE 'function.inc'



!     Initialize the residual

      Avar(:) = 0.0
    
!
!     Calculate residual
!

!     if (use_send_recv) then
!       call send_recv(var,2)
!     else
!       if(myPE.eq.root) then
! 	     allocate (var_g(1:ijkmax3))
!     else
!     allocate (var_g(1:10))
!       endif
!       call gather(var,var_g)
!       call scatter(var,var_g)

!       call MPI_Barrier(MPI_COMM_WORLD,mpierr)
!     endif

        if (use_funijk) then
        i = istart3
        j = jstart3
        k = kstart3
        cj = funijk( i,j+1,k) - funijk(i,j,k)
        if (no_k) then
            ck = 0
        else
            ck = funijk(i,j,k+1)-funijk(i,j,k)
        endif
        ci = funijk(i+1,j,k) - funijk(i,j,k)
        cc = funijk(i,j,k) - (ci*i + cj*j  + ck*k)
	endif

      if (do_k) then

!$omp    parallel  do &
!$omp&   private(ijk_c,im1jk_c,ip1jk_c,ijm1k_c,ijp1k_c,ijkm1_c,&
!$omp&           ijkp1_c,ijk,k,i, &
!$omp&           im1jk,ip1jk,ijm1k,ijp1k,ijkm1,ijkp1)
        do k = kstart,kend
        do i = istart,iend

            ijk_c   = (cc + ck*k + ci*i)
            im1jk_c = (cc + ck*k + ci*im1(i))
            ip1jk_c = (cc + ck*k + ci*ip1(i))

            ijm1k_c = (cc + ck*k + ci*i)
            ijp1k_c = (cc + ck*k + ci*i)

            ijkm1_c = (cc + ck*km1(k) + ci*i)
            ijkp1_c = (cc + ck*kp1(k) + ci*i)

        do j = jstart,jend

        if (use_funijk) then


            ijk   = ijk_c   + cj*j

            im1jk = im1jk_c + cj*j
            ip1jk = ip1jk_c + cj*j

            ijm1k = ijm1k_c + cj*jm1(j)
            ijp1k = ijp1k_c + cj*jp1(j)

            ijkm1 = ijkm1_c + cj*j
            ijkp1 = ijkp1_c + cj*j

        else
           IJK = funijk(i,j,k)

           im1jk = im_of(ijk)
           ip1jk = ip_of(ijk)
           ijm1k = jm_of(ijk)
           ijp1k = jp_of(ijk)
!
           ijkm1 = km_of(ijk)
           ijkp1 = kp_of(ijk)
        endif


           AVar(ijk) =      A_m(ijk,-3) * Var(ijkm1)   &
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
!$omp parallel do private(i,j,k,ijk,   im1jk,ip1jk,ijm1k,ijp1k,ijkm1,ijkp1)
        do i = istart,iend
        do j = jstart,jend

        k = 1

        IJK = funijk(i,j,k)


           im1jk = im_of(ijk)
           ip1jk = ip_of(ijk)
           ijm1k = jm_of(ijk)
           ijp1k = jp_of(ijk)
           AVar(ijk) =      A_m(ijk,-2) * Var(ijm1k)   &
                          + A_m(ijk,-1) * Var(im1jk)   &
                          + A_m(ijk, 0) * Var(ijk)     &
                          + A_m(ijk, 1) * Var(ip1jk)   &
                          + A_m(ijk, 2) * Var(ijp1k)

        enddo
        enddo


      endif

      if (need_distribute_Avar) then

      if (use_send_recv) then
 	     call send_recv(Avar,1)
      else
	      call gather(Avar,var_g)
	      call scatter(Avar,var_g)
              deallocate( var_g )
      endif
      endif

      return
      END SUBROUTINE LEQ_MATVEC
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: LEQ_MSOLVE(Vname, B_m, A_m, Var, CMETHOD)
!  Purpose: Successive line over-relaxation method -- Cyclic bc        C
!                                                                      C
!  Author: Ed D'Azevedo                               Date: 21-JAN-99  C
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
      SUBROUTINE LEQ_MSOLVE(VNAME, B_m, A_M, Var, CMETHOD)
!
!-----------------------------------------------
!   M o d u l e s
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
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(ijkstart3:ijkend3, -3:3)
!
!                      Vector b_m
      DOUBLE PRECISION B_m(ijkstart3:ijkend3)
!
!                      Variable name
      CHARACTER*(*)    Vname
!
!                      Variable
      DOUBLE PRECISION Var(ijkstart3:ijkend3)


!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

!
      INTEGER ::   IJK, I , J, K , ITER, NITER
      INTEGER ::   I1 , K1 , I2, K2, IK, ISIZE, KSIZE
      INTEGER ::   ICASE
      
!      CHARACTER*4, PARAMETER :: CMETHOD = 'II'
!                    sweep direction
      CHARACTER*4 :: CMETHOD
      CHARACTER :: CH
      LOGICAL :: DO_ISWEEP, DO_JSWEEP, DO_KSWEEP
      LOGICAL :: DO_SENDRECV, DO_REDBLACK
      LOGICAL, PARAMETER :: USE_IKLOOP = .FALSE.

      LOGICAL, PARAMETER :: SETGUESS = .TRUE.


!-----------------------------------------------
      INCLUDE 'function.inc'

     IF (SETGUESS) THEN

!$omp   parallel do private(ijk)
        do k = kstart3,kend3
        do i = istart3,iend3
        do j = jstart3,jend3

        IJK = funijk(i,j,k)

           VAR(IJK) = B_M(IJK)

        enddo
        enddo
        enddo

!    call send_recv(var,1)

     ENDIF

     NITER = LEN( CMETHOD )

     DO ITER=1,NITER
!
!     Perform sweeps
!
      CH = CMETHOD( ITER:ITER )
      DO_ISWEEP = (CH .EQ. 'I') .OR. (CH .EQ. 'i')
      DO_JSWEEP = (CH .EQ. 'J') .OR. (CH .EQ. 'j')
      DO_KSWEEP = (CH .EQ. 'K') .OR. (CH .EQ. 'k')
      DO_SENDRECV = (CH .EQ. 'S') .OR. (CH .EQ. 's')
      DO_REDBLACK = (CH .EQ. 'R') .OR. (CH .EQ. 'r')

     IF (NO_K) THEN

       IF ( DO_ISWEEP ) THEN

!$omp   parallel do private(I)
        DO I=istart,iend
           CALL LEQ_ISWEEP( I, Vname, Var, A_m, B_m )
        ENDDO

      ENDIF

     ELSE

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
        i = (ik-1-(k-k1)*isize) + i1

             CALL LEQ_IKSWEEP( I,K, Vname, Var, A_m, B_m )
        ENDDO
        ENDDO

      ENDIF

      IF(USE_IKLOOP) THEN

      i1 = istart
      k1 = kstart
      i2 = iend
      k2 = kend
      isize = i2-i1+1
      ksize = k2-k1+1

      IF (DO_ISWEEP) THEN

!$omp   parallel do private(K,I,IK)
        DO IK=1, ksize*isize
        if (mod(ik,isize).ne.0) then
                k = int( ik/isize ) + k1
        else
                k = int( ik/isize ) + k1 -1
        endif
        i = (ik-1-(k-k1)*isize) + i1

             CALL LEQ_IKSWEEP( I,K, Vname, Var, A_m, B_m )
        ENDDO
      ENDIF

      IF (DO_KSWEEP) THEN

!$omp   parallel do private(K,I,IK)
        DO IK=1, ksize*isize
        if (mod(ik,ksize).ne.0) then
                i = int( ik/ksize ) + i1
        else
                i = int( ik/ksize ) + i1 -1
        endif
        k = (ik-1-(i-i1)*ksize) + k1

             CALL LEQ_IKSWEEP( I,K, Vname, Var, A_m, B_m )
        ENDDO
      ENDIF

      ELSE

      IF (DO_ISWEEP) THEN
!$omp   parallel do private(K,I)
        DO K=kstart,kend
          DO I=istart,iend
             CALL LEQ_IKSWEEP( I,K, Vname, Var, A_m, B_m )
          ENDDO
        ENDDO
      ENDIF

      IF (DO_KSWEEP) THEN

!$omp   parallel do private(K,I)
        DO I=istart,iend
          DO K=kstart,kend
             CALL LEQ_IKSWEEP( I,K, Vname, Var, A_m, B_m )
          ENDDO
        ENDDO
      ENDIF

      ENDIF

      IF (DO_SENDRECV) call send_recv(var,1)

     ENDIF

     ENDDO

    RETURN
    END SUBROUTINE LEQ_MSOLVE


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: LEQ_JKSWEEP(J, K, Vname, Var, A_m, B_m )               C
!  Purpose: Perform line sweep at coordiante I,K                       C
!                                                                      C
!                                                                      C
!  Author: Ed D'Azevedo                               Date: 21-JAN-99  C
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
      SUBROUTINE LEQ_JKSWEEP(J,K,Vname, VAR, A_M, B_M )

!-----------------------------------------------
!   M o d u l e s
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
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!                      Line position
      INTEGER          J,K
!
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3)
!
!                      Vector b_m
      DOUBLE PRECISION B_m(DIMENSION_3)
!
!                      Variable name
      CHARACTER*(*)    Vname
!
!                      Variable
      DOUBLE PRECISION Var(DIMENSION_3)

!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!


!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

      DOUBLE PRECISION, DIMENSION (IMAX2) :: CC,DD,EE,BB
      INTEGER :: NN, INFO, IJK, I

      INCLUDE 'function.inc'

      NN = IMAX2

      DO I=1,NN
         IJK = FUNIJK(I,J,K)

         DD(I) = A_M(IJK,  0)
         CC(I) = A_M(IJK, -1)
         EE(I) = A_M(IJK,  1)
         BB(I) = B_M(IJK)    -  A_M(IJK,-2) * Var( JM_OF(IJK) )         &
                             -  A_M(IJK, 2) * Var( JP_OF(IJK) )         &
                             -  A_M(IJK,-3) * Var( KM_OF(IJK) )         &
                             -  A_M(IJK, 3) * Var( KP_OF(IJK) )

     ENDDO

     CC(1) = ZERO
     EE(NN) = ZERO
!    DL(1:NEND-1) = CC(2:NEND)
     INFO = 0
     CALL DGTSL( NN, CC, DD, EE, BB, INFO )
!    CALL DGTSV( JEND-JSTART+1, 1, DL, DD, EE, BB,  JEND-JSTART+1, INFO )

     IF (INFO.NE.0) THEN
        write(unit_log,*) 'VNAME = ', VNAME
        write(unit_log,*) 'DGTSV RETURNS INFO = ', INFO
        call mfix_exit(myPE)
     ENDIF

     DO I=1,NN
        IJK = FUNIJK(I,J,K)
        Var(IJK) = BB(I)
     ENDDO

     RETURN
     END SUBROUTINE  LEQ_JKSWEEP
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: LEQ_IJSWEEP(I,J, Vname, Var, A_m, B_m )               C
!  Purpose: Perform line sweep at coordiante I,K                       C
!                                                                      C
!                                                                      C
!  Author: Ed D'Azevedo                               Date: 21-JAN-99  C
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
      SUBROUTINE LEQ_IJSWEEP(I,J,Vname, VAR, A_M, B_M )

!-----------------------------------------------
!   M o d u l e s
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
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!                      Line position
      INTEGER          I,J
!
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3)
!
!                      Vector b_m
      DOUBLE PRECISION B_m(DIMENSION_3)
!
!                      Variable name
      CHARACTER*(*)    Vname

!
!                      Variable
      DOUBLE PRECISION Var(DIMENSION_3)

!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!


!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

      DOUBLE PRECISION, DIMENSION (KMAX2) :: CC,DD,EE,BB
      INTEGER :: NN, INFO, IJK, K

      INCLUDE 'function.inc'

      NN = KMAX2

      DO K=1,NN
         IJK = FUNIJK(I,J,K)

         DD(K) = A_M(IJK,  0)
         CC(K) = A_M(IJK, -3)
         EE(K) = A_M(IJK,  3)
         BB(K) = B_M(IJK)    -  A_M(IJK,-2) * Var( JM_OF(IJK) )         &
                             -  A_M(IJK, 2) * Var( JP_OF(IJK) )         &
                             -  A_M(IJK,-1) * Var( IM_OF(IJK) )         &
                             -  A_M(IJK, 1) * Var( IP_OF(IJK) )

     ENDDO

     CC(1) = ZERO
     EE(NN) = ZERO
!    DL(1:NEND-1) = CC(2:NEND)
     INFO = 0
     CALL DGTSL( NN, CC, DD, EE, BB, INFO )
!    CALL DGTSV( JEND-JSTART+1, 1, DL, DD, EE, BB,  JEND-JSTART+1, INFO )

     IF (INFO.NE.0) THEN
        write(unit_log,*) 'VNAME = ', VNAME
        write(unit_log,*) 'DGTSV RETURNS INFO = ', INFO
        call mfix_exit(myPE)
     ENDIF

     DO K=1,NN
        IJK = FUNIJK(I,J,K)
        Var(IJK) = BB(K)
     ENDDO

     RETURN
     END SUBROUTINE  LEQ_IJSWEEP


    Double Precision function dot_product_par(r1,r2)

    use mpi_utility
    use geometry
    use compar
    use indices

    implicit none
    double precision, intent(in), dimension(:) :: r1,r2
!                         Variable
    DOUBLE PRECISION, allocatable, Dimension(:) :: r1_g, r2_g
    double precision :: prod, prod_gl
    integer :: i, j, k, ijk
    logical, parameter :: do_global_sum = .true.

    include 'function.inc'

    if(do_global_sum) then
  
      prod = 0.0d0
   
      do k = kstart, kend
        do i = istart, iend
          do j = jstart, jend
   
            ijk = funijk (imap_c(i),jmap_c(j),kmap_c(k))

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
  
!$omp parallel do private(i,j,k,ijk) reduction(+:prod)
        do k = kmin2, kmax2
          do i = imin2, imax2
            do j = jmin2, jmax2
  
              ijk = funijk_gl (imap_c(i),jmap_c(j),kmap_c(k))
  
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
    
      SUBROUTINE LEQ_MSOLVE0(VNAME, B_m, A_M, Var, CMETHOD )
!
!-----------------------------------------------
!   M o d u l e s
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
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(ijkstart3:ijkend3, -3:3)
!
!                      Vector b_m
      DOUBLE PRECISION B_m(ijkstart3:ijkend3)
!
!                      Variable name
      CHARACTER*(*)    Vname
!
!                      Variable
      DOUBLE PRECISION Var(ijkstart3:ijkend3)

      CHARACTER*4 :: CMETHOD



! do nothing or no preconditioning

     var(:) = b_m(:)
     call send_recv(var,1)

     return
     end subroutine leq_msolve0

      SUBROUTINE LEQ_msolve1(VNAME, B_m, A_M, Var, CMETHOD )
!
!-----------------------------------------------
!   M o d u l e s
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
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(ijkstart3:ijkend3, -3:3)
!
!                      Vector b_m
      DOUBLE PRECISION B_m(ijkstart3:ijkend3)
!
!                      Variable name
      CHARACTER*(*)    Vname
!
!                      Variable
      DOUBLE PRECISION Var(ijkstart3:ijkend3)

      CHARACTER*4 :: CMETHOD

	integer :: i,j,k, ijk 

    include 'function.inc'

	var(:) = ZERO

! diagonal scaling

	do k=kstart2,kend2
	do i=istart2,iend2
	do j=jstart2,jend2

	 ijk = funijk( i,j,k )
	 var(ijk) = b_m(ijk)/A_m(ijk,0)

	enddo
	enddo
	enddo

     call send_recv(var,1)

     return
     end subroutine leq_msolve1

      SUBROUTINE DGTSV( N, NRHS, DL, D, DU, B, LDB, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1999
!
!     .. Scalar Arguments ..
      INTEGER            INFO, LDB, N, NRHS
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   B( LDB, * ), D( * ), DL( * ), DU( * )
!     ..
!
!  Purpose
!  =======
!
!  DGTSV  solves the equation
!
!     A*X = B,
!
!  where A is an n by n tridiagonal matrix, by Gaussian elimination with
!  partial pivoting.
!
!  Note that the equation  A'*X = B  may be solved by interchanging the
!  order of the arguments DU and DL.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.
!
!  DL      (input/output) DOUBLE PRECISION array, dimension (N-1)
!          On entry, DL must contain the (n-1) sub-diagonal elements of
!          A.
!
!          On exit, DL is overwritten by the (n-2) elements of the
!          second super-diagonal of the upper triangular matrix U from
!          the LU factorization of A, in DL(1), ..., DL(n-2).
!
!  D       (input/output) DOUBLE PRECISION array, dimension (N)
!          On entry, D must contain the diagonal elements of A.
!
!          On exit, D is overwritten by the n diagonal elements of U.
!
!  DU      (input/output) DOUBLE PRECISION array, dimension (N-1)
!          On entry, DU must contain the (n-1) super-diagonal elements
!          of A.
!
!          On exit, DU is overwritten by the (n-1) elements of the first
!          super-diagonal of U.
!
!  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
!          On entry, the N by NRHS matrix of right hand side matrix B.
!          On exit, if INFO = 0, the N by NRHS solution matrix X.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!          > 0: if INFO = i, U(i,i) is exactly zero, and the solution
!               has not been computed.  The factorization has not been
!               completed unless i = N.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, J
      DOUBLE PRECISION   FACT, TEMP
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     ..
!     .. Executable Statements ..
!
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGTSV ', -INFO )
         RETURN
      END IF
!
      IF( N.EQ.0 ) RETURN
!
      IF( NRHS.EQ.1 ) THEN
         DO 10 I = 1, N - 2
            IF( ABS( D( I ) ).GE.ABS( DL( I ) ) ) THEN
!
!              No row interchange required
!
               IF( D( I ).NE.ZERO ) THEN
                  FACT = DL( I ) / D( I )
                  D( I+1 ) = D( I+1 ) - FACT*DU( I )
                  B( I+1, 1 ) = B( I+1, 1 ) - FACT*B( I, 1 )
               ELSE
                  INFO = I
                  RETURN
               END IF
               DL( I ) = ZERO
            ELSE
!
!              Interchange rows I and I+1
!
               FACT = D( I ) / DL( I )
               D( I ) = DL( I )
               TEMP = D( I+1 )
               D( I+1 ) = DU( I ) - FACT*TEMP
               DL( I ) = DU( I+1 )
               DU( I+1 ) = -FACT*DL( I )
               DU( I ) = TEMP
               TEMP = B( I, 1 )
               B( I, 1 ) = B( I+1, 1 )
               B( I+1, 1 ) = TEMP - FACT*B( I+1, 1 )
            END IF
   10    CONTINUE
         IF( N.GT.1 ) THEN
            I = N - 1
            IF( ABS( D( I ) ).GE.ABS( DL( I ) ) ) THEN
               IF( D( I ).NE.ZERO ) THEN
                  FACT = DL( I ) / D( I )
                  D( I+1 ) = D( I+1 ) - FACT*DU( I )
                  B( I+1, 1 ) = B( I+1, 1 ) - FACT*B( I, 1 )
               ELSE
                  INFO = I
                  RETURN
               END IF
            ELSE
               FACT = D( I ) / DL( I )
               D( I ) = DL( I )
               TEMP = D( I+1 )
               D( I+1 ) = DU( I ) - FACT*TEMP
               DU( I ) = TEMP
               TEMP = B( I, 1 )
               B( I, 1 ) = B( I+1, 1 )
               B( I+1, 1 ) = TEMP - FACT*B( I+1, 1 )
            END IF
         END IF
         IF( D( N ).EQ.ZERO ) THEN
            INFO = N
            RETURN
         END IF
      ELSE
         DO 40 I = 1, N - 2
            IF( ABS( D( I ) ).GE.ABS( DL( I ) ) ) THEN
!
!              No row interchange required
!
               IF( D( I ).NE.ZERO ) THEN
                  FACT = DL( I ) / D( I )
                  D( I+1 ) = D( I+1 ) - FACT*DU( I )
                  DO 20 J = 1, NRHS
                     B( I+1, J ) = B( I+1, J ) - FACT*B( I, J )
   20             CONTINUE
               ELSE
                  INFO = I
                  RETURN
               END IF
               DL( I ) = ZERO
            ELSE
!
!              Interchange rows I and I+1
!
               FACT = D( I ) / DL( I )
               D( I ) = DL( I )
               TEMP = D( I+1 )
               D( I+1 ) = DU( I ) - FACT*TEMP
               DL( I ) = DU( I+1 )
               DU( I+1 ) = -FACT*DL( I )
               DU( I ) = TEMP
               DO 30 J = 1, NRHS
                  TEMP = B( I, J )
                  B( I, J ) = B( I+1, J )
                  B( I+1, J ) = TEMP - FACT*B( I+1, J )
   30          CONTINUE
            END IF
   40    CONTINUE
         IF( N.GT.1 ) THEN
            I = N - 1
            IF( ABS( D( I ) ).GE.ABS( DL( I ) ) ) THEN
               IF( D( I ).NE.ZERO ) THEN
                  FACT = DL( I ) / D( I )
                  D( I+1 ) = D( I+1 ) - FACT*DU( I )
                  DO 50 J = 1, NRHS
                     B( I+1, J ) = B( I+1, J ) - FACT*B( I, J )
   50             CONTINUE
               ELSE
                  INFO = I
                  RETURN
               END IF
            ELSE
               FACT = D( I ) / DL( I )
               D( I ) = DL( I )
               TEMP = D( I+1 )
               D( I+1 ) = DU( I ) - FACT*TEMP
               DU( I ) = TEMP
               DO 60 J = 1, NRHS
                  TEMP = B( I, J )
                  B( I, J ) = B( I+1, J )
                  B( I+1, J ) = TEMP - FACT*B( I+1, J )
   60          CONTINUE
            END IF
         END IF
         IF( D( N ).EQ.ZERO ) THEN
            INFO = N
            RETURN
         END IF
      END IF
!
!     Back solve with the matrix U from the factorization.
!
      IF( NRHS.LE.2 ) THEN
         J = 1
   70    CONTINUE
         B( N, J ) = B( N, J ) / D( N )
         IF( N.GT.1 ) &
            B( N-1, J ) = ( B( N-1, J )-DU( N-1 )*B( N, J ) ) / D( N-1 )
         DO 80 I = N - 2, 1, -1
            B( I, J ) = ( B( I, J )-DU( I )*B( I+1, J )-DL( I )* &
                        B( I+2, J ) ) / D( I )
   80    CONTINUE
         IF( J.LT.NRHS ) THEN
            J = J + 1
            GO TO 70
         END IF
      ELSE
         DO 100 J = 1, NRHS
            B( N, J ) = B( N, J ) / D( N )
            IF( N.GT.1 ) &
               B( N-1, J ) = ( B( N-1, J )-DU( N-1 )*B( N, J ) ) / &
                             D( N-1 )
            DO 90 I = N - 2, 1, -1
               B( I, J ) = ( B( I, J )-DU( I )*B( I+1, J )-DL( I )*&
                           B( I+2, J ) ) / D( I )
   90       CONTINUE
  100    CONTINUE
      END IF
!
      RETURN
!
!     End of DGTSV
!
      END SUBROUTINE DGTSV

!// Comments on the modifications for DMP version implementation      
!// 400 Send Receive two updated ghostlayers of Var to every processor
