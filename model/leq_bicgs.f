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
      USE indices
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!
!                      Error indicator
      INTEGER          IER
!                      maximum number of iterations
      INTEGER          ITMAX
!
!
!                      convergence tolerance
      DOUBLE PRECISION TOL
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

!                    sweep direction
      CHARACTER*4 :: CMETHOD
!-------------------------------------------------
      DOUBLE PRECISION DNRM2
      EXTERNAL LEQ_MATVEC, LEQ_MSOLVE


!--------------------------------------------------
      
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
      
!                    sweep direction
      CHARACTER*4 :: CMETHOD

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, PARAMETER :: idebug = 0
      DOUBLE PRECISION :: ratiotol = 0.2

      DOUBLE PRECISION, DIMENSION(IJKMAX2) ::                       &
                                R,Rtilde, P,Phat, Svec, Shat, Tvec,V
      DOUBLE PRECISION, DIMENSION(0:ITMAX+1) :: alpha,beta,omega,rho
      DOUBLE PRECISION :: TxS, TxT, oam,RtildexV,                   &
		      RtildexR, aijmax, Rnorm, Rnorm0, Snorm, TOLMIN
      LOGICAL :: isconverged
      INTEGER :: i,ijk, itemp, iter

!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      EXTERNAL  MATVEC, MSOLVE
!-----------------------------------------------
      INCLUDE 'function.inc'


      TOLMIN = EPSILON( one )

!
!     Scale matrix to have unit diagonal
!
!$omp parallel do private(ijk,oam,aijmax)
      DO IJK = 1, IJKMAX2

         aijmax =                                                 &
             max(                                                 &
                 max(                                             &
		      max( abs(A_M(ijk,-1)), abs(A_M(ijk,1)) ),   &
		      max( abs(A_M(ijk,-2)), abs(A_M(ijk,2)) )    &
                     ),                                           &
                 max(                                             &
		      max( abs(A_M(ijk,-3)), abs(A_M(ijk,3)) ),   &
		      abs(A_M(ijk,0))                             &
                    )                                             &
               )
         OAM = one/aijmax
         
         A_M(IJK,:) = A_M(IJK,:)*OAM

         B_M(IJK) = B_M(IJK)*OAM
      END DO

!
!    Compute initial residual, assume initial guess in Var
!    r = b - A*x
!    rtilde = r

    call MATVEC( Vname, Var, A_M, R )

    if (USE_DOLOOP) then

       Rnorm0 = ZERO
!$omp  parallel do private(ijk) reduction(+:Rnorm0)
       do ijk=1,ijkmax2
            R(ijk) = B_m(ijk) - R(ijk)
            Rtilde(ijk) = R(ijk)
            Rnorm0 = Rnorm0 + R(ijk)*R(ijk)
       enddo
       Rnorm0 = sqrt( Rnorm0 )
    else

	R(1:ijkmax2) = B_m(1:ijkmax2) - R(1:ijkmax2)
	Rtilde(:) = R(:)
	Rnorm0 = sqrt( dot_product( R(:), R(:) ) )
    endif

    if (idebug >= 1) then
       print*,'leq_bicgs, initial: ', Vname,' resid ', real(Rnorm0)
    endif
!
!   Main loop
!
    iter = 1
    do i=1,itmax
      if (USE_DOLOOP) then
         RtildexR = ZERO

!$omp    parallel do private(ijk) reduction(+:RtildexR)
         do ijk = 1, ijkmax2
           RtildexR = RtildexR + Rtilde(ijk)*R(ijk)
         enddo
         rho(i-1) = RtildexR

      else
        rho(i-1) = dot_product( Rtilde(:), R(:) )
      endif

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
	  return
        endif

        if (i .eq. 1) then
           if (USE_DOLOOP) then
!$omp          parallel do private(ijk)
               do ijk=1,ijkmax2
                   P(ijk) = R(ijk)
               enddo
           else
	       P(:) = R(:)
           endif
        else
           beta(i-1) = ( rho(i-1)/rho(i-2) )*( alpha(i-1) / omega(i-1) )
           if (USE_DOLOOP) then
!$omp         parallel do private(ijk)
              do ijk=1,ijkmax2
                 P(ijk) = R(ijk) + beta(i-1)*( P(ijk) - omega(i-1)*V(ijk) )
              enddo
           else
              P(:) = R(:) + beta(i-1)*( P(:) - omega(i-1)*V(:) )
           endif
        endif

!
!       Solve M Phat(:) = P(:)
!       V(:) = A*Phat(:)
!       
        call MSOLVE( Vname, P, A_m, Phat, CMETHOD )

        call MATVEC( Vname, Phat, A_m, V )

            
        if (USE_DOLOOP) then
           RtildexV = ZERO

!$omp      parallel do private(ijk) reduction(+:RtildexV)
           do ijk=1,ijkmax2
               RtildexV = RtildexV + Rtilde(ijk)*V(ijk) 
           enddo

        else
           RtildexV = dot_product( Rtilde(:), V(:) )
        endif

        alpha(i) = rho(i-1) / RtildexV

        if (USE_DOLOOP) then

!$omp       parallel do private(ijk)
            do ijk=1,ijkmax2
              Svec(ijk) = R(ijk) - alpha(i) * V(ijk)
            enddo

        else
            Svec(:) = R(:) - alpha(i) * V(:)
        endif

!
!       Check norm of Svec(:); if small enough:
!       set X(:) = X(:) + alpha(i)*Phat(:) and stop
!
        if (USE_DOLOOP) then
           Snorm = ZERO

!$omp      parallel do private(ijk) reduction(+:Snorm)
           do ijk=1,ijkmax2
               Snorm = Snorm + Svec(ijk)*Svec(ijk)
           enddo

        else
           Snorm = sqrt( dot_product( Svec(:), Svec(:) ) )
        endif

        if (Snorm <= TOLMIN) then
             if (USE_DOLOOP) then
!$omp            parallel do private(ijk)
                 do ijk=1,ijkmax2
                   Var(ijk) = Var(ijk) + alpha(i)*Phat(ijk)
                 enddo
             else
                Var(1:ijkmax2) = Var(1:ijkmax2) + alpha(i)*Phat(1:ijkmax2)
             endif
             if (idebug >= 1) then
!
!               Recompute residual norm
!          
             call MATVEC( Vname, Var, A_m, R )
             if (USE_DOLOOP) then
                Rnorm = ZERO

!$omp           parallel do private(ijk) reduction(+:Rnorm)
                do ijk=1,ijkmax2
                   R(ijk) = B_m(ijk) - R(ijk)
                   Rnorm = Rnorm + R(ijk)*R(ijk)
                enddo

                Rnorm = sqrt( Rnorm )
             else
                R(1:ijkmax2) = B_m(1:ijkmax2) - R(1:ijkmax2)
                Rnorm = sqrt( dot_product( R(:), R(:) ) )
             endif
             endif

            EXIT
        endif

!
!       Solve M Shat(:) = Svec(:)
!       Tvec(:) = A * Shat(:)
!
        call MSOLVE( Vname, Svec, A_m, Shat, CMETHOD )
        
        call MATVEC( Vname, Shat, A_m, Tvec )

        
        if (USE_DOLOOP) then
           TxS = ZERO
           TxT = ZERO

!$omp      parallel do private(ijk) reduction(+:TxS,TxT)
           do ijk=1,ijkmax2
             TxS = TxS + Tvec(ijk)*Svec(ijk)
             TxT = TxT + Tvec(ijk)*Tvec(ijk)
           enddo

        else
           TxS = dot_product( Tvec(:), Svec(:) )
           TxT = dot_product( Tvec(:), Tvec(:) )
        endif
        omega(i) = TxS / TxT


        if (USE_DOLOOP) then

           Rnorm = ZERO
!$omp      parallel do private(ijk) reduction(+:Rnorm)
           do ijk=1,ijkmax2
              Var(ijk) = Var(ijk) + alpha(i)*Phat(ijk) + omega(i)*Shat(ijk)
              R(ijk) = Svec(ijk) - omega(i)*Tvec(ijk)
              Rnorm = Rnorm + R(ijk)*R(ijk)
           enddo
           Rnorm = sqrt( Rnorm )

        else
           Var(1:ijkmax2) = Var(1:ijkmax2) +                           &
                 alpha(i)*Phat(1:ijkmax2) + omega(i)*Shat(1:ijkmax2)
           R(:) = Svec(:) - omega(i)*Tvec(:)
           Rnorm = sqrt( dot_product(R(:), R(:)) )
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


        if (idebug >= 1) then
           print*,'leq_bicgs ratio : ', Vname,' ',iter,     &
                   ' L-2', real(Rnorm/Rnorm0)
           if (real(Rnorm/Rnorm0) .gt. ratiotol) then
              call leq_dump( Vname, Var, A_m, B_m )
	      print*,'Rnorm ', real(Rnorm), ' Rnorm0 ', real(Rnorm0)
           endif
        endif 

        isconverged = (real(Rnorm) <= TOL*Rnorm0);
        IER = 0
        if (.not.isconverged) then
            IER = -1
            if (real(Rnorm) >= ratiotol*real(Rnorm0)) then
                  IER = -2
            endif
        endif
        
        return
        end subroutine LEQ_BICGS0
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: LEQ_LSOR(Vname, Var, A_m, B_m,  ITMAX, IER)          C
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
      SUBROUTINE LEQ_LSOR(VNAME, VAR, A_M, B_M,  ITMAX, IER)
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE matrix
      USE geometry
      USE indices
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
!                      maximum number of iterations
      INTEGER          ITMAX
!
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
!                      OVERRELAXATION FACTOR
      DOUBLE PRECISION, PARAMETER :: OMEGA = 1.2
      logical, parameter :: need_scale = .false.
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

      integer          idebug
      parameter(idebug=0)
!
      INTEGER          I,  J, K, IJK, ITER ,itemp
      DOUBLE PRECISION oAm
      double precision :: resid1,resid2, rmax1,rmax2,aijmax
      double precision, parameter :: ratiotol = 0.9


!-----------------------------------------------
      INCLUDE 'function.inc'

!
!     Scale matrix to have unit diagonal
!
    if (need_scale) then
      DO IJK = 1, IJKMAX2

         OAM = one/A_M(ijk,0)

         A_M(ijk,:) = A_M(ijk,:) * OAM
         B_M(IJK) = B_M(IJK)*OAM
      END DO
     endif

!
!     Calculate residual
!
     if (idebug >= 1) then
       call leq_resid( Vname, Var, A_m, B_m,  resid1, rmax1 )
       print*,'leq_lsor, initial: ',Vname, ' resid1, rmax1 ',      &
              real(resid1), real(rmax1)
     endif





     DO ITER=1,ITMAX
!
!     Perform symmetric sweeps
!

      IF (NO_K) THEN
        DO I=1,IMAX2
           CALL LEQ_ISWEEP( I, Vname, Var, A_m, B_m )
        ENDDO

        DO I=IMAX2,1,-1
           CALL LEQ_ISWEEP( I, Vname, Var, A_m, B_m )
        ENDDO
      ELSE
        DO K=1,KMAX2
          DO I=1,IMAX2
             CALL LEQ_IKSWEEP( I,K, Vname, Var, A_m, B_m )
          ENDDO
        ENDDO

        DO K=KMAX2,1,-1
           DO I=IMAX2,1,-1
              CALL LEQ_IKSWEEP( I,K, Vname, Var, A_m, B_m )
           ENDDO
        ENDDO
      ENDIF
    ENDDO


!
!     Calculate residual
!
     if (idebug >= 1) then
       call leq_resid( Vname, Var, A_m, B_m,  resid2, rmax2 )
       print*,'leq_lsor, final: ',Vname, ' resid2, rmax2 ',       &
               real(resid2), real(rmax2)

       print*,'leq_lsor ratio : ',Vname,' L-2 ',                  &
               real(resid2/resid1), ' L-inf ',real(rmax2/rmax1)

       if (min( real(resid2/resid1), real(rmax2/rmax1) ) .gt. ratiotol ) then
           call leq_dump(  Vname, Var, A_m, B_m )
           print*,'leq_lsor ', Vname
           print*,' resid1, rmax1 ', resid1,rmax1
           print*,' resid2, rmax2 ', resid2,rmax2
           stop '** error ** '
       endif

      endif

    RETURN
    END SUBROUTINE LEQ_LSOR


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
      SUBROUTINE LEQ_ISWEEP(I,Vname, VAR, A_M, B_M )

!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE matrix
      USE geometry
      USE indices
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

      INTEGER :: NN 
      DOUBLE PRECISION, DIMENSION (JMAX2) :: CC,DD,EE,BB
      INTEGER :: INFO, IJK, J, K, IM1JK, IP1JK

      INCLUDE 'function.inc'


      NN = JMAX2
      K = 1

      DO J=1,NN
         IJK = FUNIJK(I,J,K)
         IM1JK = IM_OF(IJK)
         IP1JK = IP_OF(IJK)

         DD(J) = A_M(IJK,  0)
         CC(J) = A_M(IJK, -2)
         EE(J) = A_M(IJK,  2)
         BB(J) = B_M(IJK) -  A_M(IJK,-1) * Var( IM1JK )         &
                          -  A_M(IJK, 1) * Var( IP1JK )

     ENDDO

     CC(1) = ZERO
     EE(NN) = ZERO
     INFO = 0
     CALL DGTSL( NN, CC, DD, EE, BB, INFO )

     IF (INFO.NE.0) THEN
        PRINT *,'VNAME = ', VNAME
        PRINT*,'DGTSL RETURNS INFO = ', INFO
        STOP 'ERROR'
     ENDIF

     DO J=1,NN
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
      USE indices
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

      INTEGER :: NN 
      DOUBLE PRECISION, DIMENSION (JMAX2) :: CC,DD,EE,BB
      INTEGER :: NN, INFO, IJK, J,  IM1JK, IP1JK, IJKM1, IJKP1

      INCLUDE 'function.inc'

      NN = JMAX2

!$omp parallel do private(j,ijk,im1jk,ip1jk,ijkm1,ijkp1)
      DO J=1,NN
         IJK = FUNIJK(I,J,K)
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

     CC(1) = ZERO
     EE(NN) = ZERO
     INFO = 0
     CALL DGTSL( NN, CC, DD, EE, BB, INFO )

     IF (INFO.NE.0) THEN
        PRINT *,'VNAME = ', VNAME
        PRINT*,'DGTSL RETURNS INFO = ', INFO
        STOP 'ERROR'
     ENDIF

     DO J=1,NN
        IJK = FUNIJK(I,J,K)
        Var(IJK) = BB(J)
     ENDDO

     RETURN
     END SUBROUTINE  LEQ_IKSWEEP


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: LEQ_RESID(Vname, Var, A_m, B_m,  resid2, rmax)       C
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
      SUBROUTINE LEQ_resid(VNAME, VAR, A_M, B_M,  resid2,rmax)
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE matrix
      USE geometry
      USE indices
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
!                      maximum number of iterations
      INTEGER          ITMAX
!
!                      phase index
      INTEGER          M
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
!
!     L-2 norm of residual vector
      DOUBLE PRECISION RESID2
!
!     max norm of residual vector
! 
      DOUBLE PRECISION RMAX

!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!
!                      OVERRELAXATION FACTOR
      DOUBLE PRECISION, PARAMETER :: OMEGA = 1.2
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

!
      INTEGER          I,  K, IJK, ITER 

      double precision resid_ijk
      integer :: im1jk,ip1jk, ijm1k,ijp1k, ijkm1,ijkp1


!-----------------------------------------------
      INCLUDE 'function.inc'


    
!
!     Calculate residual
!

      rmax = ZERO
      resid2 = ZERO
     

      if (do_k) then

!$omp   parallel  do &
!$omp&  private(im1jk,ip1jk,ijm1k,ijp1k,ijkm1,ijkp1,resid_ijk) &
!$omp&  reduction(MAX:rmax) reduction(+:resid2)
        do ijk=1,ijkmax2

           im1jk = im_of(ijk)
           ip1jk = ip_of(ijk)
           ijm1k = jm_of(ijk)
           ijp1k = jp_of(ijk)
 
           ijkm1 = km_of(ijk)
           ijkp1 = kp_of(ijk)

           resid_ijk  = abs( B_m(ijk)                  &
                          - A_m(ijk, 0) * Var(ijk)     &
                          - A_m(ijk,-1) * Var(im1jk)   &
                          - A_m(ijk, 1) * Var(ip1jk)   &
                          - A_m(ijk,-2) * Var(ijm1k)   &
                          - A_m(ijk, 2) * Var(ijp1k)   &
                          - A_m(ijk,-3) * Var(ijkm1)   &
                          - A_m(ijk, 3) * Var(ijkp1)   &
                          )
          rmax = MAX( rmax, resid_ijk )
          resid2 = resid2 + resid_ijk*resid_ijk
        enddo
        resid2 = sqrt( resid2 )

      else
!$omp   parallel do &
!$omp&  private(im1jk,ip1jk,ijm1k,ijp1k,ijkm1,ijkp1,resid_ijk) &
!$omp&  reduction(MAX:rmax) reduction(+:resid2)
        do ijk=1,ijkmax2

           im1jk = im_of(ijk)
           ip1jk = ip_of(ijk)
           ijm1k = jm_of(ijk)
           ijp1k = jp_of(ijk)
 

           resid_ijk  = abs( B_m(ijk)                  &
                          - A_m(ijk, 0) * Var(ijk)     &
                          - A_m(ijk,-1) * Var(im1jk)   &
                          - A_m(ijk, 1) * Var(ip1jk)   &
                          - A_m(ijk,-2) * Var(ijm1k)   &
                          - A_m(ijk, 2) * Var(ijp1k)   &
                          )
          rmax = MAX( rmax, resid_ijk )
          resid2 = resid2 + resid_ijk*resid_ijk
        enddo

        resid2 = sqrt( resid2 )


      endif


      return
      END SUBROUTINE LEQ_RESID


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: LEQ_dump( Vname, Var, A_m, B_m)                        C
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
      SUBROUTINE LEQ_dump(VNAME, VAR, A_M, B_M )


!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE matrix
      USE geometry
      USE indices
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3)
!
!                      Vector b_m
      DOUBLE PRECISION B_m(DIMENSION_3)
!
!
!
!                      Variable name
      CHARACTER*(*)    Vname
!
!                      Variable
      DOUBLE PRECISION Var(DIMENSION_3)
!
!     L-2 norm of residual vector
      DOUBLE PRECISION RESID2
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!
!                      I/O device
      integer, parameter :: iodev = 2
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      character*80 :: filename
      integer :: i,j,k,ijk

      include 'function.inc'


      filename = 'B_m.m'
      open(unit=iodev,file=filename,form='formatted',   &
           access='sequential')
      rewind(iodev)

      write(iodev,*) '% ', Vname
      do ijk=1,ijkmax2
        write(iodev,9001) ijk,B_m(ijk)
 9001   format('B(',i6,') = ', e14.4, ';')
      enddo

      write(iodev,*) 'B = B(:);'
      close(iodev)

      
      filename = 'Var.m'
      open(unit=iodev,file=filename,form='formatted',   &
           access='sequential')
      rewind(iodev)


      write(iodev,*) '% ', Vname
      do ijk=1,ijkmax2
        write(iodev,9002) ijk, Var(ijk)
 9002   format('X(',i6,') = ', e14.4, ';' )
      enddo

      write(iodev,*) 'X = X(:); '
      close(iodev)



      filename = 'A_m.m'
      open(unit=iodev,file=filename,form='formatted',   &
           access='sequential')
      rewind(iodev)
      
      write(iodev,*) '% ',Vname

      write(iodev,*) 'A = sparse( ',ijkmax2,',',ijkmax2,');';
      do k=1,kmax2
      do j=1,jmax2
      do i=1,imax2
        ijk = funijk(i,j,k)

        if (i-1.ge.1) then 
           write(iodev,9003) ijk, im_of(ijk), A_m(ijk,-1)
        endif
        if (i+1.le.imax2) then
          write(iodev,9003) ijk, ip_of(ijk), A_m(ijk, 1)
        endif
        if (j-1.ge.1) then
          write(iodev,9003) ijk, jm_of(ijk), A_m(ijk,-2)
        endif
        if (j+1.le.jmax2) then
          write(iodev,9003) ijk, jp_of(ijk), A_m(ijk, 2)
        endif


        if (do_k) then
          if (k-1.ge.1) then
             write(iodev,9003) ijk, km_of(ijk), A_m(ijk,-3)
          endif
          if (k+1.le.kmax2) then
             write(iodev,9003) ijk, kp_of(ijk), A_m(ijk, 3)
          endif
        endif

        write(iodev,9003) ijk, ijk, A_m(ijk,0)

 9003   format('A(',i6,',',i6,') = ',e14.4,';')
      enddo
      enddo
      enddo

      close(iodev)

   
      return
      end subroutine leq_dump
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
      USE indices
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
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3)
!
!                      Vector AVar
      DOUBLE PRECISION AVar(DIMENSION_3)
!
!                      Variable name
      CHARACTER*(*)    Vname
!
!                      Variable
      DOUBLE PRECISION Var(DIMENSION_3)
!
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

!
      INTEGER          I,  K, IJK, ITER 
      DOUBLE PRECISION oAm

      double precision, dimension (ijkmax2) :: resid
      integer :: im1jk,ip1jk, ijm1k,ijp1k, ijkm1,ijkp1


!-----------------------------------------------
      INCLUDE 'function.inc'


    
!
!     Calculate residual
!

     

      if (do_k) then

!$omp parallel  do private(im1jk,ip1jk,ijm1k,ijp1k,ijkm1,ijkp1)
        do ijk=1,ijkmax2

           im1jk = im_of(ijk)
           ip1jk = ip_of(ijk)
           ijm1k = jm_of(ijk)
           ijp1k = jp_of(ijk)
 
           ijkm1 = km_of(ijk)
           ijkp1 = kp_of(ijk)

           AVar(ijk) =      A_m(ijk,-3) * Var(ijkm1)   &
                          + A_m(ijk,-2) * Var(ijm1k)   &
                          + A_m(ijk,-1) * Var(im1jk)   &
                          + A_m(ijk, 0) * Var(ijk)     &
                          + A_m(ijk, 1) * Var(ip1jk)   &
                          + A_m(ijk, 2) * Var(ijp1k)   &
                          + A_m(ijk, 3) * Var(ijkp1)    
        enddo

      else
!$omp parallel do private(im1jk,ip1jk,ijm1k,ijp1k,ijkm1,ijkp1)
        do ijk=1,ijkmax2

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


      endif

      return
      END SUBROUTINE LEQ_MATVEC
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: LEQ_MSOLVE(Vname, B_m, A_m, Var, CMETHOD )
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
      SUBROUTINE LEQ_MSOLVE(VNAME, B_m, A_M, Var, CMETHOD )
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE matrix
      USE geometry
      USE indices
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
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

!
      INTEGER ::   IJK, I , J, K , ITER, NITER
      
!      CHARACTER*4, PARAMETER :: CMETHOD = 'II'
!                    sweep direction
      CHARACTER*4 :: CMETHOD
      CHARACTER :: CH
      LOGICAL :: DO_ISWEEP, DO_JSWEEP, DO_KSWEEP

      LOGICAL, PARAMETER :: SETGUESS = .TRUE.


!-----------------------------------------------
      INCLUDE 'function.inc'


     IF (SETGUESS) THEN

!$omp   parallel do private(ijk)
        DO IJK=1,IJKMAX2
           VAR(IJK) = B_M(IJK)
        ENDDO
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

      IF (NO_K) THEN

       IF ( DO_ISWEEP ) THEN

!$omp   parallel do private(I)
        DO I=1,IMAX2
           CALL LEQ_ISWEEP( I, Vname, Var, A_m, B_m  )
        ENDDO

      ENDIF
      IF (DO_JSWEEP) THEN

!$omp   parallel do private(J)
        DO J=1,JMAX2
           CALL LEQ_JSWEEP( J, Vname, Var, A_m, B_m )
        ENDDO

      ENDIF

     ELSE

      IF (DO_ISWEEP) THEN

!$omp   parallel do private(K,I)
        DO K=1,KMAX2
          DO I=1,IMAX2
             CALL LEQ_IKSWEEP( I,K, Vname, Var, A_m, B_m  )
          ENDDO
        ENDDO
      ENDIF
      IF (DO_JSWEEP) THEN

!$omp   parallel do private(J,K)
        DO K=1,KMAX2
          DO J=1,JMAX2
             CALL LEQ_JKSWEEP( J,K, Vname, Var, A_m, B_m )
          ENDDO
        ENDDO

      ENDIF
      IF (DO_KSWEEP) THEN

!$omp    parallel do private(I,J)
         DO J=1,JMAX2
           DO I=1,IMAX2
             CALL LEQ_IJSWEEP( I,J, Vname, Var, A_m, B_m )
           ENDDO
         ENDDO

      ENDIF

      ENDIF

      ENDDO


    RETURN
    END SUBROUTINE LEQ_MSOLVE
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: LEQ_MSOLVE22(Vname, B_m, A_m, Var )                    C
!  Purpose: Successive relaxation method over hyperplanes              C
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
      SUBROUTINE LEQ_MSOLVE2(VNAME, B_m, A_M, Var  )
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE matrix
      USE geometry
      USE indices
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
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

!
!                      OVERRELAXATION FACTOR
      DOUBLE PRECISION, PARAMETER :: OMEGA = 1.2

!
      INTEGER ::   IJK, I , J, K , ITER
      INTEGER, PARAMETER :: NITER = 2
      
      logical, parameter :: setguess = .true.


      integer :: Lstart,Lend,Linc,L,iplane
      integer :: ncount,ilevel,nlevel, istart,iend,inode, iposition
      integer, allocatable, dimension(:) :: xlist,list

      integer, parameter :: idebug = 0

      logical :: isfirst = .true.
      logical :: isodd, isvalid
      save
!-----------------------------------------------
      INCLUDE 'function.inc'


     if (setguess) then

!$omp   parallel do private(ijk)
        do ijk=1,ijkmax2
           Var(ijk) = B_m(ijk)
        enddo
     endif

     if (isfirst) then
          isfirst = .false.

          nlevel = (imax2+jmax2+kmax2)-3+1
          allocate( xlist(nlevel+1) )
          allocate( list(ijkmax2) )
         

         ilevel = 1
         xlist(ilevel) = 1
         do iplane=3,imax2+jmax2+kmax2

            ncount = 0

            do k=1,min(kmax2, iplane-2)
             do i=1,min(imax2, iplane-k-1)

                  j = iplane - i - k
                  isvalid = (1 <= j) .and. (j <= jmax2)
                  if (isvalid) then

                     ijk = funijk(i,j,k)
                     ncount = ncount + 1
                     iposition = xlist(ilevel)-1+ncount
                     list(iposition) = ijk
                  endif

              enddo
            enddo

            xlist(ilevel+1) = xlist(ilevel) + ncount
            ilevel = ilevel + 1
          enddo

          if (idebug >= 1) then
	    print*,'nlevel ', nlevel
	    do ilevel=1,nlevel
		print*,'ilevel: ', ilevel
		istart = xlist(ilevel)
                iend = xlist(ilevel+1)-1
                print*, 'istart,iend ', istart,iend

                do inode=istart,iend
                    ijk = list(inode)
		    i = i_of(ijk)
		    j = j_of(ijk)
		    k = k_of(ijk)
           
                    isvalid = (ijk .eq. funijk(i,j,k) )
		    if ((idebug >= 2) .or. (.not. isvalid)) then
                       print*,'i,j,k, ',i,j,k, '  ijk ',ijk
                    endif
                enddo
             enddo
           endif


      endif
        


     DO ITER=1,NITER
!
!     Perform sweeps
!
      isodd = (mod(iter,2) .ne. 0)
      if (isodd) then

!        forward sweep

         Lstart = 1
         Lend = nlevel
         Linc = 1
      else
!        backward sweep 
         Lstart = nlevel
         Lend = 1
         Linc = -1
      endif

      IF (NO_K) THEN

         do L=Lstart,Lend,Linc
            istart = xlist(L)
            iend = xlist(L+1)-1

!$omp       parallel do private(inode,ijk)
            do inode=istart,iend
               ijk = list(inode)
               Var(ijk) = Var(ijk) + &
                  omega*(                                        &
			  B_m(ijk) - A_m(ijk,-1)*Var( im_of(ijk) ) &
                                   - A_m(ijk, 1)*Var( ip_of(ijk) ) &
                                   - A_m(ijk,-2)*Var( jm_of(ijk) ) &
                                   - A_m(ijk, 2)*Var( jp_of(ijk) ) &
                                   - Var(ijk)                    &
                          )
            enddo
          enddo
               


      ELSE


         do L=Lstart,Lend,Linc
            istart = xlist(L)
            iend = xlist(L+1)-1

!$omp       parallel do private(inode,ijk)
            do inode=istart,iend
               ijk = list(inode)
               Var(ijk) = Var(ijk) + &
                  omega*(                                        &
			  B_m(ijk) - A_m(ijk,-1)*Var( im_of(ijk) ) &
                                   - A_m(ijk, 1)*Var( ip_of(ijk) ) &
                                   - A_m(ijk,-2)*Var( jm_of(ijk) ) &
                                   - A_m(ijk, 2)*Var( jp_of(ijk) ) &
                                   - A_m(ijk,-3)*Var( km_of(ijk) ) &
                                   - A_m(ijk, 3)*Var( kp_of(ijk) ) &
                                   - Var(ijk)                    &
                        )
            enddo
          enddo


      ENDIF

      ENDDO


    RETURN
    END SUBROUTINE LEQ_MSOLVE2
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: LEQ_JSWEEP(I, Vname, Var, A_m, B_m )                   C
!  Purpose: Perform line sweep at coordiante J                         C
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
      SUBROUTINE LEQ_JSWEEP(J,Vname, VAR, A_M, B_M )

!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE matrix
      USE geometry
      USE indices
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!                      Line position
      INTEGER          J
!
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

      INTEGER :: NN 
      DOUBLE PRECISION, DIMENSION (IMAX2) :: CC,DD,EE,BB
      INTEGER :: INFO, IJK, I, K

      INCLUDE 'function.inc'


      NN = IMAX2
      K = 1

      DO I=1,NN
         IJK = FUNIJK(I,J,K)


         DD(I) = A_M(IJK,  0)
         CC(I) = A_M(IJK, -1)
         EE(I) = A_M(IJK,  1)
         BB(I) = B_M(IJK) -  A_M(IJK,-2) * Var( JM_OF(IJK) )         &
                          -  A_M(IJK, 2) * Var( JP_OF(IJK) )

     ENDDO

     CC(1) = ZERO
     EE(NN) = ZERO
     INFO = 0
     CALL DGTSL( NN, CC, DD, EE, BB, INFO )

     IF (INFO.NE.0) THEN
        PRINT *,'VNAME = ', VNAME
        PRINT*,'DGTSL RETURNS INFO = ', INFO
        STOP 'ERROR'
     ENDIF

     DO I=1,NN
        IJK = FUNIJK(I,J,K)
        Var(IJK) =  BB(I) 
     ENDDO

     RETURN
     END SUBROUTINE  LEQ_JSWEEP
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

      INTEGER :: NN 
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
     INFO = 0
     CALL DGTSL( NN, CC, DD, EE, BB, INFO )

     IF (INFO.NE.0) THEN
        PRINT *,'VNAME = ', VNAME
        PRINT*,'DGTSL RETURNS INFO = ', INFO
        STOP 'ERROR'
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

      INTEGER :: NN 
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
     INFO = 0
     CALL DGTSL( NN, CC, DD, EE, BB, INFO )

     IF (INFO.NE.0) THEN
        PRINT *,'VNAME = ', VNAME
        PRINT*,'DGTSL RETURNS INFO = ', INFO
        STOP 'ERROR'
     ENDIF

     DO K=1,NN
        IJK = FUNIJK(I,J,K)
        Var(IJK) = BB(K)
     ENDDO

     RETURN
     END SUBROUTINE  LEQ_IJSWEEP
