!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_trD_s(trD_s, IER)                                 C
!  Purpose: Calculate the trace of solids phase rate of strain tensor  C
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
      SUBROUTINE CALC_TRD_S(TRD_S, IER) 
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
      USE physprop
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
      INTEGER          I, J, K, IJK, IMJK, IJMK, IJKM, IM, M
!
!                      Strain rate tensor components for mth solids phase
      DOUBLE PRECISION trD_s(DIMENSION_3, DIMENSION_M)
!-----------------------------------------------
      INCLUDE 'function.inc'
!
!
      DO M = 1, MMAX 
!!$omp    parallel do private(ijk,i,j,k,im,imjk,ijmk,ijkm)
!//SP
      DO IJK = IJKSTART3, IJKEND3
            IF (.NOT.WALL_AT(IJK)) THEN 
               I = I_OF(IJK) 
               J = J_OF(IJK) 
               K = K_OF(IJK) 
               IM = IM1(I) 
               IMJK = IM_OF(IJK) 
               IJMK = JM_OF(IJK) 
               IJKM = KM_OF(IJK) 
               TRD_S(IJK,M) = (X_E(I)*U_S(IJK,M)-X_E(IM)*U_S(IMJK,M))*OX(I)*ODX&
                  (I) + (V_S(IJK,M)-V_S(IJMK,M))*ODY(J) + (W_S(IJK,M)-W_S(IJKM,&
                  M))*(OX(I)*ODZ(K)) 
            ENDIF 
         END DO 
      END DO 
      RETURN  
!//SP
      call send_recv(TRD_S,2)
      END SUBROUTINE CALC_TRD_S 
