!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: calc_external_forces(IER)                              C
!  Purpose: Calculate the 3-components of external forces at cell      C
!           faces                                                      C
!                                                                      C
!  Author: S. Benyahia                              Date: 13-MAY-09    C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!     C. Hrenya handnotes and Garzo, Hrenya, Dufty papers (PRE, 2007)  C
!                                                                      C
!  Variables referenced:                                               C
!                                                                      C
!  Variables modified:   JoiX  JoiY   JoiZ                             C
!                                                                      C
!     Local variables: all terms in mass flux                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CALC_EXTERNAL_FORCES(IER)
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
      USE ghdtheory
      USE physprop
      USE run
      USE constant
      USE drag
      USE bc
      use scales
      IMPLICIT NONE
!-----------------------------------------------
!     Local variables
!-----------------------------------------------  
!                      Index
      INTEGER          IJK, I, J, K   
! 
!                      Index
      INTEGER          IJKE, IJKN, IJKT
!     
!                      Solids phase
      INTEGER          M 
!     
!     Error index
      INTEGER          IER  
!     
!     number densities to compute del(Nj)
      DOUBLE PRECISION NjC, NjE, NjN, NjT  
!     
!     mass, volume of species
      DOUBLE PRECISION Mj, Vj
!     
!     mixture density at cell faces
      DOUBLE PRECISION ropsE, ropsN, ropsT
!     
!     drag force on a particle
      DOUBLE PRECISION dragFc, dragFe, dragFn, dragFt, dragFx, dragFy, dragFz
!     
!     pressure terms in mass mobility
      DOUBLE PRECISION PGE, PGN, PGT, SDPx, SDPy, SDPz
!
!----------------------------------------------- 
!     Function subroutines
!----------------------------------------------- 

!-----------------------------------------------
!     Include statement functions
!-----------------------------------------------
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'b_force1.inc'
      INCLUDE 'b_force2.inc'
!-----------------------------------------------   
     
      DO 200 IJK = ijkstart3, ijkend3
          I = I_OF(IJK)
          J = J_OF(IJK)
          K = K_OF(IJK)

     
          IF ( FLUID_AT(IJK) ) THEN      
               
	       IJKE = EAST_OF(IJK)
	       IJKN = NORTH_OF(IJK)
               IJKT = TOP_OF(IJK)

! pressure term (no need to have it inside smax loop)
               PGE = P_G(IJKE) 
               PGN = P_G(IJKN)
               PGT = P_G(IJKT) 
               IF (CYCLIC_X_PD) THEN 
                 IF (CYCLIC_AT_E(IJK)) PGE = P_G(IJKE) - DELP_X 
               ENDIF 
               IF (CYCLIC_Y_PD) THEN 
                 IF (CYCLIC_AT_N(IJK)) PGN = P_G(IJKN) - DELP_Y 
               ENDIF 
               IF (CYCLIC_Z_PD) THEN 
                 IF (CYCLIC_AT_T(IJK)) PGT = P_G(IJKT) - DELP_Z 
               ENDIF
               SDPx = -P_SCALE*(PGE - P_G(IJK))  * oDX_E(I)
               SDPy = -P_SCALE*(PGN - P_G(IJK))  * oDY_N(J)
               SDPz = -P_SCALE*(PGT - P_G(IJK))  * (oX_E(I)*oDZ_T(K))

	       DO M = 1, SMAX
                 Vj = (PI/6.d0)*D_P(IJK,M)**3 ! particle volume
                 Mj = Vj * RO_S(M)            ! particle mass

                 NjC = ROP_s(IJK,M) / Mj
                 NjE = ROP_S(IJKE,M) / Mj
                 NjN = ROP_S(IJKN,M) / Mj
                 NjT = ROP_S(IJKT,M) / Mj

! drag force on a particle in -x -y -z directions
                 dragFc = zero
		 dragFe = zero 
		 dragFn = zero
		 dragFt = zero
                 if(NjC > zero) dragFc = F_GS(IJK ,M)/NjC
		 if(NjE > zero) dragFe = F_GS(IJKE,M)/NjE 
		 if(NjN > zero) dragFn = F_GS(IJKN,M)/NjN
		 if(NjT > zero) dragFt = F_GS(IJKT,M)/NjT
		 
		 dragFx = AVG_X(dragFc,dragFe,I) * (U_g(IJK) - U_s(IJK,M))
		 dragFy = AVG_Y(dragFc,dragFn,J) * (V_g(IJK) - V_s(IJK,M))
		 dragFz = AVG_Z(dragFc,dragFt,K) * (W_g(IJK) - W_s(IJK,M))
		 
		 FiX(IJK,M) =  (Mj * BFX_S(IJK,M) + dragFx + Vj*SDPx)
		 FiY(IJK,M) =  (Mj * BFY_S(IJK,M) + dragFy + Vj*SDPy)
		 FiZ(IJK,M) =  (Mj * BFZ_S(IJK,M) + dragFz + Vj*SDPz)
               ENDDO
          ENDIF     ! Fluid_at
 200  CONTINUE     ! outer IJK loop

      RETURN
      END SUBROUTINE CALC_EXTERNAL_FORCES
