!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_K_s(M, IER)                                       C
!  Purpose: Calculate the effective conductivity of solids phases      C
!                                                                      C
!  Author:M. Syamlal                                  Date: 24-APR-96  C
!  Reviewer:                                          Date: dd-mmm-yy  C
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
      SUBROUTINE CALC_K_S(M, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE physprop
      USE fldvar
      USE geometry
      USE indices
      USE constant
      USE toleranc 
      USE compar     !//d
      USE sendrecv   !// 400
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!
!                      conductivity of solids in cal/s.cm.K
      DOUBLE PRECISION K_p
      PARAMETER (K_p = 0.2)
!
!                      Constant in conductivity equation
      DOUBLE PRECISION PHI_k
      PARAMETER (PHI_k = 7.26E-3)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

!                      Error index
      INTEGER          IER
!
!                      Indices
      INTEGER          IJK
!
!                      Solids phase
      INTEGER          M
!
!                      Quantities in solids conductivity formula
      DOUBLE PRECISION B, R_km, BoR, L_rm
!-----------------------------------------------
      INCLUDE 'function.inc'
!
!
      IF (K_S0 /= UNDEFINED) RETURN  

!!$omp parallel do private(IJK,B,R_km,BoR,L_rm) &
!!$omp& schedule(dynamic,chunk_size)


!// 350 1112 MTP changed do loop limits 1,ijkmax2 ==> ijkstart3, ijkend3
      DO IJK = IJKSTART3, IJKEND3            
!
!     Solids conductivity in cal/s.cm.K
!     Bauer & Schlunder's (1978) theory
!          B = 1.25 * ((ONE - EP_g(IJK))/EP_g(IJK))**(10./9.)
!          IF( (ONE - EP_g(IJK)) .GT. DIL_EP_s) THEN
!            R_km = K_p/K_g(IJK)
!            BoR  = B/R_km
!            L_rm = -(2./(ONE - BoR))
!     &           * ( (R_km - ONE) * BoR / (ONE - BoR)**2 * LOG(BoR)
!     &               + (B - ONE)/(ONE - BoR) + (B + ONE)/2.         )
!            K_s(IJK, M) = ( Phi_k * R_km + (ONE - Phi_k) * L_rm )
!     &                      * K_g(IJK) / SQRT(ONE - EP_g(IJK))
!          ELSE
!            K_s(IJK, M) = ZERO
!          ENDIF
!     An approximate average value for the solids conductivity is 2.5*K_g
         IF (.NOT.WALL_AT(IJK)) K_S(IJK,M) = 2.5*K_G(IJK) 
      END DO 

!//S 1113 try to move this COMM to the end of transport_prop to do all COMMs
!//       at certain locations, provided that no data dependency in between.

!// 400 1113 MTP communicate boundaries
      CALL SEND_RECV(K_S, 2)     
    
      RETURN  
      END SUBROUTINE CALC_K_S 
