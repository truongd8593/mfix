!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_trD_g(trD_g, IER)                                 C
!  Purpose: Calculate the trace of gas phase rate of strain tensor     C
!                                                                      C
!  Author: M. Syamlal                                 Date: 19-DEC-96  C
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
      SUBROUTINE CALC_TRD_G(TRD_G, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!     Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE geometry
      USE fldvar
      USE indices
      USE compar   !//d
      USE sendrecv   !//SP
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Error index
      INTEGER          IER
!                      Indices
      INTEGER          I, J, K, IJK, IMJK, IJMK, IJKM, IM
!
!                      Strain rate tensor components for mth solids phase
      DOUBLE PRECISION trD_g(DIMENSION_3)
!-----------------------------------------------
      INCLUDE 'function.inc'
!
!//SP
      TRD_G(:) = ZERO
!
!!$omp  parallel do private( IJK, I,J,K, IM,IMJK,IJMK,IJKM ) &
!!$omp& schedule(dynamic,chunk_size)
!//SP
      DO IJK = IJKSTART3, IJKEND3 
         IF (FLUID_AT(IJK)) THEN 
            I = I_OF(IJK) 
            J = J_OF(IJK) 
            K = K_OF(IJK) 
            IM = IM1(I) 
            IMJK = IM_OF(IJK) 
            IJMK = JM_OF(IJK) 
            IJKM = KM_OF(IJK) 
            TRD_G(IJK) = (X_E(I)*U_G(IJK)-X_E(IM)*U_G(IMJK))*OX(I)*ODX(I) + (&
               V_G(IJK)-V_G(IJMK))*ODY(J) + (W_G(IJK)-W_G(IJKM))*(OX(I)*ODZ(K)) 
         ENDIF 
      END DO 
      RETURN  
!//SP
      call send_recv(TRD_G,2)
      END SUBROUTINE CALC_TRD_G 
