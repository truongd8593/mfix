!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: LEQ_SOR(Vname, Var, A_m, B_m, M, ITMAX, IER)          C
!  Purpose: Successive over-relaxation method -- Cyclic bc             C
!                                                                      C
!  Author: M. Syamlal                                 Date: 19-AUG-96  C
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
      SUBROUTINE LEQ_SOR(VNAME, VAR, A_M, B_M, M, ITMAX, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
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
!
!                      OVERRELAXATION FACTOR
      DOUBLE PRECISION, PARAMETER :: OMEGA = 1.2 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer          idebug
      parameter( idebug = 0 )
 
! 
      INTEGER          I, J, K, IJK, ITER, IJK01, IJK02, IJK11, IJK12 
      DOUBLE PRECISION oAm 
  
      double precision :: resid1,resid2,rmax1,rmax2

!-----------------------------------------------
      INCLUDE 'function.inc'
!
!
!     Successive Over relaxation method
!
!$omp parallel do private(IJK,OAM)
      DO IJK = 1, IJKMAX2 
!
         OAM = ONE/A_M(IJK,0,M) 

         A_M(IJK,0,M) = ONE

         A_M(IJK,-2,M) = A_M(IJK,-2,M)*OAM 
         A_M(IJK,-1,M) = A_M(IJK,-1,M)*OAM 
         A_M(IJK,1,M) = A_M(IJK,1,M)*OAM 
         A_M(IJK,2,M) = A_M(IJK,2,M)*OAM 
         A_M(IJK,-3,M) = A_M(IJK,-3,M)*OAM 
         A_M(IJK,3,M) = A_M(IJK,3,M)*OAM 
         B_M(IJK,M) = B_M(IJK,M)*OAM 
      END DO 

!
!     Calculate residual
!
     if (idebug >= 1) then
       call leq_resid( Vname, Var, A_m(:,:,M), B_m(:,M),  resid1, rmax1 )
       print*,'LEQ_SOR, initial: ',Vname, ' resid1, rmax1 ',      &
              real(resid1), real(rmax1)
     endif


      DO ITER = 1, ITMAX 
!
!  SOR procedure
!
         IF (DO_K) THEN 

!$omp parallel do private(IJK)
            DO IJK = 1, IJKMAX2 
!
               VAR(IJK) = VAR(IJK) + OMEGA*(B_M(IJK,M)-A_M(IJK,-2,M)*VAR(JM_OF(&
                  IJK))-A_M(IJK,-1,M)*VAR(IM_OF(IJK))-A_M(IJK,1,M)*VAR(IP_OF(&
                  IJK))-A_M(IJK,2,M)*VAR(JP_OF(IJK))-A_M(IJK,-3,M)*VAR(KM_OF(&
                  IJK))-A_M(IJK,3,M)*VAR(KP_OF(IJK))-VAR(IJK)) 
            END DO 
         ELSE 

!$omp parallel do private(IJK)
            DO IJK = 1, IJKMAX2 
!
               VAR(IJK) = VAR(IJK) + OMEGA*(B_M(IJK,M)-A_M(IJK,-2,M)*VAR(JM_OF(&
                  IJK))-A_M(IJK,-1,M)*VAR(IM_OF(IJK))-A_M(IJK,1,M)*VAR(IP_OF(&
                  IJK))-A_M(IJK,2,M)*VAR(JP_OF(IJK))-VAR(IJK)) 
            END DO 
         ENDIF 
      END DO 
!
!     Calculate residual
!
      if (idebug >= 1) then
       call leq_resid( Vname, Var, A_m(:,:,M), B_m(:,M), resid2, rmax2 )
       print*,'LEQ_SOR, final: ',Vname, ' resid2, rmax2 ',       &
               real(resid2), real(rmax2)
       print*,'LEQ_SOR ratio: ',Vname,' L-2 ', &
              real(resid2/resid1), ' L-inf ',real(rmax2/rmax1)
      endif

      RETURN  
      END SUBROUTINE LEQ_SOR 
!
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: LEQ_SORX(Vname, Var, A_m, B_m, M, ITMAX, IER)
!  Purpose: Successive over-relaxation method                          C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 19-AUG-96  C
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
      SUBROUTINE LEQ_SORX(VNAME, VAR, A_M, B_M, M, ITMAX, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
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
!
!                      OVERRELAXATION FACTOR
      DOUBLE PRECISION, PARAMETER :: OMEGA = 1.2 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IJK, ITER 
      DOUBLE PRECISION :: OAM 
      double precision :: resid2,rmax

!-----------------------------------------------
      include 'function.inc'
!
!
!
!     Successive Over relaxation method
!
      DO IJK = 1, IJKMAX2 
!
         OAM = ONE/A_M(IJK,0,M) 
         A_M(IJK,-2,M) = A_M(IJK,-2,M)*OAM 
         A_M(IJK,-1,M) = A_M(IJK,-1,M)*OAM 
         A_M(IJK,1,M) = A_M(IJK,1,M)*OAM 
         A_M(IJK,2,M) = A_M(IJK,2,M)*OAM 
         A_M(IJK,-3,M) = A_M(IJK,-3,M)*OAM 
         A_M(IJK,3,M) = A_M(IJK,3,M)*OAM 
         B_M(IJK,M) = B_M(IJK,M)*OAM 
      END DO 
      DO ITER = 1, ITMAX 
!
!
!  SOR procedure
!
         IF (DO_K) THEN 
            IJK = 1 
            IF (IJKMAX2 > 0) THEN 
!
               DO IJK = 1, IJKMAX2 
                  VAR(IJK) = VAR(IJK) + OMEGA*(B_M(IJK,M)-A_M(IJK,-2,M)*VAR(1-&
                     IMAX2+IJK-1)-A_M(IJK,-1,M)*VAR(-1+IJK)-A_M(IJK,1,M)*VAR(1+&
                     IJK)-A_M(IJK,2,M)*VAR(1+IMAX2+IJK-1)-A_M(IJK,-3,M)*VAR(1-&
                     IJMAX2+IJK-1)-A_M(IJK,3,M)*VAR(1+IJMAX2+IJK-1)-VAR(IJK)) 
               END DO 
               IJK = IJKMAX2 + 1 
            ENDIF 
         ELSE 
            IJK = 1 
            IF (IJKMAX2 > 0) THEN 
!
               DO IJK = 1, IJKMAX2 
                  VAR(IJK) = VAR(IJK) + OMEGA*(B_M(IJK,M)-A_M(IJK,-2,M)*VAR(IJK&
                     -IMAX2)-A_M(IJK,-1,M)*VAR(IJK-1)-A_M(IJK,1,M)*VAR(1+IJK)-&
                     A_M(IJK,2,M)*VAR(IMAX2+IJK)-VAR(IJK)) 
               END DO 
               IJK = IJKMAX2 + 1 
            ENDIF 
         ENDIF 
      END DO 

      RETURN  
      END SUBROUTINE LEQ_SORX 
