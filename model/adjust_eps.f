!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: ADJUST_EPS                                             C
!  Purpose: Eliminate the solids phases that occupy only very small    C
!           fractions of the computational cell volume                 C
!                                                                      C
!  Author: M. Syamlal                                 Date: 28-JAN-92  C
!  Reviewer: P. Nicoletti, W. Rogers, S. Venkatesan   Date: 29-JAN-92  C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: Include CONSTANTS.INC                                      C
!  Author: M. Syamlal                                 Date: 7-FEB-92   C
!  Reviewer: P. Nicoletti                             Date: 11-DEC-92  C
!                                                                      C
!  Literature/Document References: None                                C
!                                                                      C
!  Variables referenced: FLAG, RO_s, ROP_s, U_s, V_s, W_s, EP_g        C
!                        IMAX1, JMAX1, KMAX1, MMAX, IMIN1, JMIN1, KMIN1C
!  Variables modified: ROP_s, U_s, V_s, W_s, EP_g, ROP_g, I, J, K, IJK,C
!                      RO_s                                            C
!                                                                      C
!  Local variables: EPSUM                                              C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE ADJUST_EPS 

!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE toleranc 
      USE constant
      USE fldvar
      USE geometry
      USE indices
      USE physprop
      USE run
      USE compar     
      USE sendrecv
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!
!                      Indices
      INTEGER          I, J, K, IJK
!
!                      Solids phase
      INTEGER          M
!                      Sum of (very small) solids volume fractions that
!                      are set to zero.
      DOUBLE PRECISION EPSUM, epsMix, epSolid
!-----------------------------------------------
      INCLUDE 'function.inc'
!
!

      DO K = Kstart1, Kend1 
         DO J = Jstart1, Jend1 
            DO I = Istart1, Iend1 

               IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells

               IJK = FUNIJK(I,J,K) 
               IF (FLUID_AT(IJK)) THEN 
                  EPSUM = ZERO 
		  epsMix = ZERO
                  DO M = 1, SMAX 
!QX
                     epSolid = ROP_S(IJK,M)/RO_S(IJK,M)
		     epsMix = epsMix +  epSolid
		     IF (epSolid < ZERO_EP_S) THEN 

!  Remove solids in very small quantities and set solids velocity to zero
!  if there is outflow from the present cell.

                        EPSUM = EPSUM + epSolid
                        IF(TRIM(KT_TYPE) == 'GHD') ROP_S(IJK,MMAX) = ROP_S(IJK,MMAX) - ROP_S(IJK,M)
                        ROP_S(IJK,M) = ZERO 
                        U_S(IJK,M) = MIN(U_S(IJK,M),ZERO) 
                        V_S(IJK,M) = MIN(V_S(IJK,M),ZERO) 
                        W_S(IJK,M) = MIN(W_S(IJK,M),ZERO) 
                        U_S(IM_OF(IJK),M) = MAX(U_S(IM_OF(IJK),M),ZERO) 
                        V_S(JM_OF(IJK),M) = MAX(V_S(JM_OF(IJK),M),ZERO) 
                        W_S(KM_OF(IJK),M) = MAX(W_S(KM_OF(IJK),M),ZERO) 
                     ENDIF 
                  END DO 
                  epsMix = epsMix - EPSUM
		  IF(TRIM(KT_TYPE) == 'GHD' .AND. epsMix < ZERO_EP_S) THEN
                    U_S(IJK,MMAX) = MIN(U_S(IJK,MMAX),ZERO) 
                    V_S(IJK,MMAX) = MIN(V_S(IJK,MMAX),ZERO) 
                    W_S(IJK,MMAX) = MIN(W_S(IJK,MMAX),ZERO) 
                    U_S(IM_OF(IJK),MMAX) = MAX(U_S(IM_OF(IJK),MMAX),ZERO) 
                    V_S(JM_OF(IJK),MMAX) = MAX(V_S(JM_OF(IJK),MMAX),ZERO) 
                    W_S(KM_OF(IJK),MMAX) = MAX(W_S(KM_OF(IJK),MMAX),ZERO)  
                  ENDIF
                  EP_G(IJK) = EP_G(IJK) + EPSUM 
                  ROP_G(IJK) = RO_G(IJK)*EP_G(IJK) 
               ENDIF 
            END DO 
         END DO 
      END DO            
     
!// Communicate field variables calculated in the do i,j,k loop
      call send_recv(ROP_S,2)
      call send_recv(U_S,2) 
      call send_recv(V_S,2) 
      call send_recv(W_S,2) 
      call send_recv(EP_G,2) 
      call send_recv(ROP_G,2)

      RETURN  
      END SUBROUTINE ADJUST_EPS
