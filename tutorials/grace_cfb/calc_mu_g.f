!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_MU_g(IER)                                         C
!  Purpose: Calculate the effective viscosity for a turbulent flow,    C
!           which is the sum of molecular and eddy viscosities         C
!                                                                      C
!  Author: W. Sams/M. Syamlal                         Date: 18-JUL-94  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: MFIX 2.0 mods (previous name CALC_MU_gt)                   C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
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
      SUBROUTINE CALC_MU_G(IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE physprop
      USE geometry
      USE fldvar
      USE visc_g
      USE visc_s
      USE indices
      USE constant
      USE compar    
      USE sendrecv 
      USE scalars
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
      DOUBLE PRECISION, PARAMETER :: F2O3 = 2./3. 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!
!                      Error index
      INTEGER          IER
!                      Indices
      INTEGER          I, J, K, IJK, IMJK, IPJK, IJMK, IJPK, IJKM, IJKP, &
                      IM, JM, KM
      INTEGER          IMJPK, IMJMK, IMJKP, IMJKM, IPJKM, IPJMK, IJMKP, &
                      IJMKM, IJPKM

!
!                      Strain rate tensor components for mth solids phase
      DOUBLE PRECISION D_g(3,3)
!
!                      U_g at the north face of the THETA cell-(i, j+1/2, k)
      DOUBLE PRECISION U_g_N
!
!                      U_g at the south face of the THETA cell-(i, j-1/2, k)
      DOUBLE PRECISION U_g_S
!
!                      U_g at the top face of the THETA cell-(i, j, k+1/2)
      DOUBLE PRECISION U_g_T
!
!                      U_g at the bottom face of the THETA cell-(i, j, k-1/2)
      DOUBLE PRECISION U_g_B
!
!                      U_g at the center of the THETA cell-(i, j, k)
!                      Calculated for Cylindrical coordinates only.
      DOUBLE PRECISION U_g_C
!
!                      V_g at the east face of the THETA cell-(i+1/2, j, k)
      DOUBLE PRECISION V_g_E
!
!                      V_g at the west face of the THETA cell-(i-1/2, j, k)
      DOUBLE PRECISION V_g_W
!
!                      V_g at the top face of the THETA cell-(i, j, k+1/2)
      DOUBLE PRECISION V_g_T
!
!                      V_g at the bottom face of the THETA cell-(i, j, k-1/2)
      DOUBLE PRECISION V_g_B
!
!                      W_g at the east face of the THETA cell-(i+1/2, j, k)
      DOUBLE PRECISION W_g_E
!
!                      W_g at the west face of the THETA cell-(1-1/2, j, k)
      DOUBLE PRECISION W_g_W
!
!                      W_g at the north face of the THETA cell-(i, j+1/2, k)
      DOUBLE PRECISION W_g_N
!
!                      W_g at the south face of the THETA cell-(i, j-1/2, k)
      DOUBLE PRECISION W_g_S
!
!                      W_g at the center of the THETA cell-(i, j, k).
!                      Calculated for Cylindrical coordinates only.
      DOUBLE PRECISION W_g_C
!
!                      Second invariant of the deviator of D_g
      DOUBLE PRECISION I2_devD_g
      DOUBLE PRECISION C_MU
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'ep_s2.inc'
      INCLUDE 'fun_avg2.inc'
!
!
!!$omp parallel do private(ijk) schedule(dynamic,chunk_size)
      DO IJK = ijkstart3, ijkend3 
         IF (FLUID_AT(IJK)) THEN 
!
!  Molecular viscosity
!
            IF (MU_G0 == UNDEFINED) MU_G(IJK) = 1.7D-4*(T_G(IJK)/273.0)**1.5*(&
               383./(T_G(IJK)+110.)) 
            MU_GT(IJK) = MU_G(IJK) 
            LAMBDA_GT(IJK) = -F2O3*MU_GT(IJK)

	 ELSE
            MU_G(IJK)  = ZERO 
            MU_GT(IJK) = ZERO 
            LAMBDA_GT(IJK) = ZERO
         ENDIF 
      END DO 

!!$omp parallel do &
!!$omp$ schedule(dynamic,chunk_size) &
!!$omp$ private(IJK, I,J,K,IM,JM,KM, &
!!$omp& IMJK,IPJK,IJMK,IJPK,IJKM,IJKP,IMJPK,IMJMK,IMJKP, &
!!$omp& IMJKM,IPJKM,IPJMK,IJMKP,IJMKM,IJPKM, &
!!$omp& U_G_N,U_G_S,U_G_T,U_G_B,V_G_E,V_G_W,V_G_T,V_G_B, &
!!$omp$ W_G_N,W_G_S,W_G_E,W_G_W,  U_G_C,W_G_C, D_G,I2_DEVD_G )
      DO IJK = ijkstart3, ijkend3
!//SP
!! When K and Epsilon are computed, no need to compute strain again.
!SOF
!
         IF ( FLUID_AT(IJK) .AND. L_SCALE(IJK)/=ZERO) THEN 
	 IF (SCALAR(IJK,1).GT. SMALL_NUMBER .AND. SCALAR(IJK,2) .GT. SMALL_NUMBER) THEN	 
!
	C_MU = 9D-02
            MU_GT(IJK) = MU_G(IJK) +  RO_G(IJK)*C_MU*Scalar(IJK, 1)**2	&
	                  /SCALAR(IJK, 2)
	ELSE
	    MU_GT(IJK) = MU_G(IJK)
	ENDIF  !for scalars with very small values
!	    		  
	    MU_GT(IJK) = MIN(MU_GMAX, MU_GT(IJK))
            LAMBDA_GT(IJK) = -F2O3*MU_GT(IJK)
!
	    
         ENDIF 
      END DO 

      RETURN  
      END SUBROUTINE CALC_MU_G 

!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization 
!// 350 Changed do loop limits: 1,ijkmax2-> ijkstart3, ijkend3
