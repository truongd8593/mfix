!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: TRANSPORT_COEFF_GHD (M, IER)                           C
!  Purpose: Calculate all transport coefficients for use in GHD theory C
!                                                                      C
!  Author: S. Benyahia                              Date: 03-MAR-09    C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!     C. Hrenya handnotes and Garzo, Hrenya, Dufty papers (PRE, 2007)  C
!                                                                      C
!  Variables referenced:                                               C
!                                                                      C
!  Variables modified:  All Transport Coefficients and species temp.   C
!                                                                      C
!     Local variables:                                                 C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE TRANSPORT_COEFF_GHD (M, IER)
!
!-----------------------------------------------
!     Modules 
!-----------------------------------------------  
      USE param
      USE param1
      USE geometry
      USE compar
      USE fldvar
      USE indices
      USE visc_s
      USE ghdtheory
      USE physprop
      USE run
      USE constant
      USE toleranc
      IMPLICIT NONE
!-----------------------------------------------
!     Local variables
!-----------------------------------------------  
!                      Index
      INTEGER          IJK, I, J, K
!     
!                      Solids phase
      INTEGER          M, L 
!     
!     Error index
      INTEGER          IER  
!     
!     particles properties and GHD transport coefficients
      DOUBLE PRECISION SIGMAI(smax), Mi(smax), phii(smax), TMix, tmpZeta0, tmpZetaU, &
                       Ti(smax), tmpP, tmpKappa, tmpEta, tmpDT(smax), tmpLambda, &
		       tmpLij(smax,smax), tmpDij(smax,smax), tmpDF(smax,smax), &
		       tmpDijQ(smax,smax)
!
!----------------------------------------------- 
!     Function subroutines
!----------------------------------------------- 

!-----------------------------------------------
!     Include statement functions
!-----------------------------------------------
      INCLUDE 'function.inc'
!-----------------------------------------------   
     
      DO 200 IJK = ijkstart3, ijkend3
          I = I_OF(IJK)
          J = J_OF(IJK)
          K = K_OF(IJK)

     
          IF ( FLUID_AT(IJK) ) THEN

               DO M = 1, SMAX
                 SIGMAI(M) = D_P(IJK,M)
                 Mi(M)     = (PI/6.d0)*SIGMAI(M)**3 * RO_S(M)
                 phii(M)   = ROP_S(IJK,M) / RO_S(M)
               ENDDO
               TMix = THETA_M(IJK,MMAX)
               
               CALL GHD(SMAX, SIGMAI, r_p(:smax,:smax), Mi, phii, TMix, tmpZeta0, &
                        tmpZetaU, Ti, tmpP, tmpKappa, tmpEta, tmpDT, tmpDF, &
			tmpLambda, tmpLij, tmpDij, tmpDijQ)

!              zeroth order cooling rate
               Zeta0(IJK) =  tmpZeta0

!              first order cooling rate transport coefficient
               ZetaU(IJK) =  tmpZetaU

!              species temperature Ti and GHD-specific transport properties
               DO M = 1, SMAX
                 THETA_M(IJK,M) = Ti(M)
!                Thermal diffusivity
                 DiT(IJK,M) = tmpDT(M)
                 DO L = 1, SMAX
!                mass mobility coefficient
                   DijF(IJK,M,L) = tmpDF(M,L)
!                thermal mobility coefficient
                   Lij(IJK,M,L) = tmpLij(M,L)
!                ordinary diffucsion coefficient
                   Dij(IJK,M,L) = tmpDij(M,L)
!                Dufour coefficient:	   
                   DijQ(IJK,M,L) = tmpDijQ(M,L)
		 ENDDO
               ENDDO

!              mixture solids pressure
               P_s_v(IJK) = tmpP

!              mixture bulk viscosity
               Mu_b_v(IJK) = tmpKappa

!              mixture kinetic-collisional viscosity
               Mu_s_v(IJK) = tmpEta

!              mixture total bulk viscosity
               LAMBDA_S_V(IJK) = Mu_b_v(IJK) - (2.d0/3.d0)*Mu_s_v(IJK)

!              granular conductivity in Mth solids phase
               Kth_s(IJK,M) = tmpLambda

!              Boyle-Massoudi stress coefficient
               ALPHA_s_v(IJK) = ZERO

          ENDIF     ! Fluid_at
 200  CONTINUE     ! outer IJK loop

      RETURN
      END SUBROUTINE TRANSPORT_COEFF_GHD
