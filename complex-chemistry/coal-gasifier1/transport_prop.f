!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: TRANSPORT_PROP(VISC, COND, DIFF, IER)                  C
!  Purpose: Calculate transport properties that vary with time         C
!                                                                      C
!  Author: M. Syamlal                                 Date: 17-JUL-92  C
!  Reviewer: P. Nicoletti                             Date: 11-DEC-92  C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: Mods for MFIX 2.0 (old name CALC_PHYSPROP)                 C
!  Author: M. Syamlal                                 Date: 23-APR-96  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: None                                          C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables: None                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE TRANSPORT_PROP(VISC, COND, DIFF, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE fldvar
      USE physprop
      USE geometry
      USE indices
      USE run
      USE toleranc 
      USE visc_g
      USE compar
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                       Error index
      INTEGER          IER
!
!                      Indices
      INTEGER          IJK, M, N
!
!                      Flags to tell whether to calculate or not
      LOGICAL          VISC(0:DIMENSION_M), COND(0:DIMENSION_M), &
                       DIFF(0:DIMENSION_M)
      INCLUDE 'function.inc'
!                 
!-----------------------------------------------
!
! 1.2.1 Fluid viscosity  (effective viscosity is modified later in this routine)
      IF (VISC(0)) CALL CALC_MU_g (IER) 
      
!
! 1.2.2 Solids viscosity
      DO M = 1, MMAX 
         IF (VISC(M)) CALL CALC_MU_S (M, IER) 
      END DO 
      
      M = 1
!
!$omp  parallel do private(IJK)  
         DO IJK = IJKSTART3, IJKEND3 
!----------------------------------------------------------------------	  
!
! 1.2.1    Fluid viscosity
!         Mu_gt(IJK) = 0.24 !effective turbulent viscosity based on a Fluent run
	 
!         IF (.NOT.WALL_AT(IJK)) THEN 
           IF (FLUID_AT(IJK)) THEN
        Mu_gt(IJK) = 0.24 !effective turbulent viscosity based on a Fluentrun
! 1.2.3    Fluid conductivity
!
!          assuming turbulent Prandtl number of 0.85
           K_g(IJK) = 0.85 * C_pg(IJK) * MU_gt(IJK)
	   
	   
!
! 1.2.4    Solids conductivity
!          Solids conductivity in cal/s.cm.K
!          modified Gort's model to account for in-bed radiation
!
           K_s(IJK, M) = 2. * 1.355E-12 * d_p(M) * T_s(IJK,M)**3
	   
	   
!
! 1.2.5    Fluid diffusivity
!          assuming turbulent Schmidt number of 0.7
           DO N = 1, NMAX (0)
               DIF_g(IJK,N) = 0.7 * MU_gt(IJK)
           END DO 
           
	   
!
! 1.2.6    Solids diffusivity
           DO N = 1, NMAX (M)
               DIF_S(IJK,M,N) = ZERO 
           END DO 
	   
	ELSE
           Mu_gt(IJK) = ZERO
           K_g(IJK) = ZERO
           K_s(IJK, M) = ZERO
           DO N = 1, NMAX (0)
               DIF_g(IJK,N) = ZERO
           END DO 
           DO N = 1, NMAX (M)
               DIF_S(IJK,M,N) = ZERO 
           END DO 
	
	ENDIF
      END DO
      RETURN  
      END SUBROUTINE TRANSPORT_PROP 
