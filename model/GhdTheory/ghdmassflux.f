!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: GHDMASSFLUX (IER)                                      C
!  Purpose: Calculate the species mass flux 3-components of Joi at cellC
!           faces to compute species velocities and source terms in T  C
!           equation.                                                  C
!                                                                      C
!  Author: S. Benyahia                              Date: 13-MAR-09    C
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
      SUBROUTINE GHDMASSFLUX (IER)
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
      USE drag
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
      INTEGER          M, L 
!     
!     Error index
      INTEGER          IER  
!     
!     number densities to compute del(Nj)
      DOUBLE PRECISION NjC, NjE, NjN, NjT  
!     
!     mass of species
      DOUBLE PRECISION Mi, Mj
!     
!     mixture density and temperature at cell faces
      DOUBLE PRECISION ropsE, ropsN, ropsT, ThetaE, ThetaN, ThetaT
      DOUBLE PRECISION EPSA
!     
!     transport coefficient at cell faces
      DOUBLE PRECISION DiTE, DiTN, DiTT
      DOUBLE PRECISION DijE, DijN, DijT
      DOUBLE PRECISION DijFE, DijFN, DijFT
!     
!     Terms in the calculation of Joi-X,Y,Z
      DOUBLE PRECISION ordinDiffTermX, ordinDiffTermY, ordinDiffTermZ
      DOUBLE PRECISION massMobilityTermX, massMobilityTermY, massMobilityTermZ
      DOUBLE PRECISION thermalDiffTermX, thermalDiffTermY, thermalDiffTermZ
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
!-----------------------------------------------   
    DO M = 1, SMAX 
      DO 200 IJK = ijkstart3, ijkend3
          I = I_OF(IJK)
          J = J_OF(IJK)
          K = K_OF(IJK)

     
          IF ( FLUID_AT(IJK) ) THEN      
               
	       IJKE = EAST_OF(IJK)
	       IJKN = NORTH_OF(IJK)
               IJKT = TOP_OF(IJK)

               Mi = (PI/6.d0)*D_P(IJK,M)**3 * RO_S(M)

! mixture density and temperature evaluated at cell faces
	       ropsE = AVG_X(ROP_S(IJK,MMAX),ROP_S(IJKE,MMAX),I)
	       ropsN = AVG_Y(ROP_S(IJK,MMAX),ROP_S(IJKN,MMAX),J)
               ropsT = AVG_Z(ROP_S(IJK,MMAX),ROP_S(IJKT,MMAX),K) 

	       ThetaE = AVG_X(THETA_M(IJK,MMAX),THETA_M(IJKE,MMAX),I)
	       ThetaN = AVG_Y(THETA_M(IJK,MMAX),THETA_M(IJKN,MMAX),J)
               ThetaT = AVG_Z(THETA_M(IJK,MMAX),THETA_M(IJKT,MMAX),K) 

! Thermal diffusion evaluated at cell faces (all used transport coef. will be evaluated this way)
	       DiTE = AVG_X(DiT(IJK,M),DiT(IJKE,M),I)
	       DiTN = AVG_Y(DiT(IJK,M),DiT(IJKN,M),J)
	       DiTT = AVG_Z(DiT(IJK,M),DiT(IJKT,M),K)

! initializing variables for summation over L
	       ordinDiffTermX = ZERO
               ordinDiffTermY = ZERO
	       ordinDiffTermZ = ZERO 
	       massMobilityTermX = ZERO
	       massMobilityTermY = ZERO
	       massMobilityTermZ = ZERO

	       DO L = 1, SMAX
                 Mj  = (PI/6.d0)*D_P(IJK,L)**3 * RO_S(L)

                 NjC = ROP_s(IJK,L) / Mj
                 NjE = ROP_S(IJKE,L) / Mj
                 NjN = ROP_S(IJKN,L) / Mj
                 NjT = ROP_S(IJKT,L) / Mj

		 DijE = AVG_X(Dij(IJK,M,L),Dij(IJKE,M,L),I)
		 DijN = AVG_Y(Dij(IJK,M,L),Dij(IJKN,M,L),J)
		 DijT = AVG_Z(Dij(IJK,M,L),Dij(IJKT,M,L),K)

		 DijFE = AVG_X(DijF(IJK,M,L),DijF(IJKE,M,L),I)
		 DijFN = AVG_Y(DijF(IJK,M,L),DijF(IJKN,M,L),J)
		 DijFT = AVG_Z(DijF(IJK,M,L),DijF(IJKT,M,L),K)
		 
		 ordinDiffTermX = ordinDiffTermX + Mj * DijE * (NjE - NjC) * oDX_E(I)
		 ordinDiffTermY = ordinDiffTermY + Mj * DijN * (NjN - NjC) * oDY_N(J)
		 ordinDiffTermZ = ordinDiffTermZ + Mj * DijT * (NjT - NjC) * (oX_E(I)*oDZ_T(K))
		 
		 massMobilityTermX = massMobilityTermX + DijFE * FiX(IJK,L)
		 massMobilityTermY = massMobilityTermY + DijFN * FiY(IJK,L)
		 massMobilityTermZ = massMobilityTermZ + DijFT * FiZ(IJK,L)
               ENDDO  

	       
	       IF(ropsE > ZERO) THEN  ! just to make sure we're not / zero
	         ordinDiffTermX = ordinDiffTermX * Mi/ropsE
	       ELSE
	         ordinDiffTermX = zero
	       ENDIF
	       IF(ropsN > ZERO) THEN
	         ordinDiffTermY = ordinDiffTermY * Mi/ropsN
	       ELSE
	         ordinDiffTermY = zero
	       ENDIF
	       IF(ropsT > ZERO) THEN
	         ordinDiffTermZ = ordinDiffTermZ * Mi/ropsT
	       ELSE
	         ordinDiffTermZ = zero
	       ENDIF
	       
	       thermalDiffTermX = ropsE*DiTE/ThetaE * ( THETA_M(IJKE,MMAX) - THETA_M(IJK,MMAX) )  * oDX_E(I)
	       thermalDiffTermY = ropsN*DiTN/ThetaN * ( THETA_M(IJKN,MMAX) - THETA_M(IJK,MMAX) )  * oDY_N(J)
	       thermalDiffTermZ = ropsT*DiTT/ThetaT * ( THETA_M(IJKT,MMAX) - THETA_M(IJK,MMAX) )  * (oX_E(I)*oDZ_T(K))
	       
	       JoiX(IJK,M) = -(ordinDiffTermX + thermalDiffTermX + massMobilityTermX)
	       JoiY(IJK,M) = -(ordinDiffTermY + thermalDiffTermY + massMobilityTermY)
	       JoiZ(IJK,M) = -(ordinDiffTermZ + thermalDiffTermZ + massMobilityTermZ)

! set fluxes to zero in case very dilute conditions	       
	       EPSA = AVG_X(ROP_S(IJK,M),ROP_S(IJKE,M),I) / RO_S(M)
	       IF(EPSA <= ZERO_EP_S) JoiX(IJK,M) = ZERO
	       EPSA = AVG_Y(ROP_S(IJK,M),ROP_S(IJKN,M),J) / RO_S(M)
	       IF(EPSA <= ZERO_EP_S) JoiY(IJK,M) = ZERO
	       EPSA = AVG_Z(ROP_S(IJK,M),ROP_S(IJKT,M),K) / RO_S(M)
	       IF(EPSA <= ZERO_EP_S) JoiZ(IJK,M) = ZERO

! set flux components to zero in case of walls        
               IF (IP_AT_E(IJK))  JoiX(IJK,M) = ZERO 
               IF (IP_AT_N(IJK))  JoiY(IJK,M) = ZERO 
               IF (IP_AT_T(IJK))  JoiZ(IJK,M) = ZERO

          ELSE
	       JoiX(IJK,M) = ZERO
	       JoiY(IJK,M) = ZERO
	       JoiZ(IJK,M) = ZERO

          ENDIF     ! Fluid_at
 200  CONTINUE     ! outer IJK loop  
    ENDDO ! for m = 1,smax

      RETURN
      END SUBROUTINE GHDMASSFLUX

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name:UpdateSpeciesVelocities(IER)                            C
!  Purpose: Update solids velocities at celll faces based on the       C
!           formula Ui = U + Joi/(mi ni); also calculate averaged      C
!           velocities for dilute conditions.                          C
!                                                                      C
!  Author: S. Benyahia                              Date: 6-MAR-09     C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!     C. Hrenya handnotes and Garzo, Hrenya, Dufty papers (PRE, 2007)  C
!                                                                      C
!  Variables referenced:                                               C
!                                                                      C
!  Variables modified: solid species velocity components: Us, Vs, Ws   C
!                                                                      C
!     Local variables:                                                 C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE UpdateSpeciesVelocities (IER)
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
      USE is
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
!                      Index
      INTEGER          IJKE, IJKN, IJKT, IJKW, IJKS, IJKB, IMJK, IJMK, IJKM, ISV
!     
!                      Solids phase
      INTEGER          M, L 
!     
!     Error index
      INTEGER          IER  
!     
!     species density at cell faces
      DOUBLE PRECISION ropsE, ropsN, ropsT
      DOUBLE PRECISION EPSA
      DOUBLE PRECISION EPSw, EPSe, EPSn, EPSs, EPSt, EPSb
      DOUBLE PRECISION tmpVel, counter
!
!----------------------------------------------- 
!     Function subroutines
!----------------------------------------------- 

!-----------------------------------------------
!     Include statement functions
!-----------------------------------------------
      INCLUDE 'function.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'
!-----------------------------------------------   
    DO M = 1, SMAX     
      DO 200 IJK = ijkstart3, ijkend3
          I = I_OF(IJK)
          J = J_OF(IJK)
          K = K_OF(IJK)

     
          IF ( FLUID_AT(IJK) ) THEN      
               
	    IJKE = EAST_OF(IJK)
	    IJKW = WEST_OF(IJK)
	    IJKN = NORTH_OF(IJK)
	    IJKS = SOUTH_OF(IJK)
            IJKT = TOP_OF(IJK)    
            IJKB = BOTTOM_OF(IJK)    
            IMJK = IM_OF(IJK)      
            IJMK = JM_OF(IJK)  
            IJKM = KM_OF(IJK)     

            EPSw = EP_S(IJKW,M)
            EPSe = EP_S(IJKE,M)
            EPSn = EP_S(IJKN,M)
            EPSs = EP_S(IJKS,M)
            IF(.NOT. NO_K) THEN
               EPSt = EP_S(IJKT,M)
               EPSb = EP_S(IJKB,M)
            ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!      First Compute Us        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	    EPSA = AVG_X(EP_S(IJK,M),EP_S(IJKE,M),I)                 
            IF (IP_AT_E(IJK)) THEN 
	       U_s(IJK, M) = ZERO
            ELSEIF (SIP_AT_E(IJK)) THEN 
               ISV = IS_ID_AT_E(IJK) 
	       U_s(IJK, M) = IS_VEL_S(ISV,M)
! Dilute flow
            ELSEIF (EPSA <= DIL_EP_S) THEN 
! using the average boundary cell values to compute U_s
               tmpVel = ZERO
	       counter = ZERO
	       IF (EPSw > DIL_EP_S .AND. .NOT.IS_AT_E(IMJK)) THEN
	          tmpVel = tmpVel + U_s(IJKW,M)
		  counter = counter + ONE
	       ENDIF
               IF (EPSe > DIL_EP_S .AND. .NOT.IS_AT_E(IJK)) THEN
	          tmpVel = tmpVel + U_s(IJKE,M)
		  counter = counter + ONE
	       ENDIF
               IF (EPSs > DIL_EP_S .AND. .NOT.IS_AT_N(IJMK)) THEN
	          tmpVel = tmpVel + U_s(IJKS,M)
		  counter = counter + ONE
	       ENDIF
               IF (EPSn > DIL_EP_S .AND. .NOT.IS_AT_N(IJK)) THEN
	          tmpVel = tmpVel + U_s(IJKN,M)
		  counter = counter + ONE
	       ENDIF
               IF(.NOT. NO_K) THEN
                  IF (EPSb > DIL_EP_S .AND. .NOT.IS_AT_T(IJKM)) THEN
	            tmpVel = tmpVel + U_s(IJKB,M)
		    counter = counter + ONE
	          ENDIF
                  IF (EPSt > DIL_EP_S .AND. .NOT.IS_AT_T(IJK)) THEN
	            tmpVel = tmpVel + U_s(IJKT,M)
		    counter = counter + ONE
	          ENDIF
               ENDIF
! compute Us as averaged of neighbours non-dilute velocities.
               IF (counter >= ONE) THEN
	          U_S(IJK,M) = tmpVel / counter
               ELSE
	          U_S(IJK,M) = U_S(IJK,M) ! same as in source_u_s, no change in Us values.
               ENDIF
! Normal case
            ELSE

! species density evaluated at cell faces
	       ropsE = AVG_X(ROP_S(IJK,M),ROP_S(IJKE,M),I)
! Calculate species velocity Us at cell faces
	       U_s(IJK, M) = U_s(IJK, MMAX) + JoiX(IJK,M) / ropsE
            ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!      Now Compute Vs        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	    EPSA =  AVG_Y(EP_S(IJK,M),EP_S(IJKN,M),J)         
            IF (IP_AT_N(IJK)) THEN 
	       V_s(IJK, M) = ZERO
            ELSEIF (SIP_AT_N(IJK)) THEN 
               ISV = IS_ID_AT_N(IJK) 
	       V_s(IJK, M) = IS_VEL_S(ISV,M) 
! Dilute flow
            ELSEIF (EPSA <= DIL_EP_S) THEN 
! using the average boundary cell values to compute V_s
               tmpVel = ZERO
	       counter = ZERO
	       IF (EPSw > DIL_EP_S .AND. .NOT.IS_AT_E(IMJK)) THEN
	          tmpVel = tmpVel + V_s(IJKW,M)
		  counter = counter + ONE
	       ENDIF
               IF (EPSe > DIL_EP_S .AND. .NOT.IS_AT_E(IJK)) THEN
	          tmpVel = tmpVel + V_s(IJKE,M)
		  counter = counter + ONE
	       ENDIF
               IF (EPSs > DIL_EP_S .AND. .NOT.IS_AT_N(IJMK)) THEN
	          tmpVel = tmpVel + V_s(IJKS,M)
		  counter = counter + ONE
	       ENDIF
               IF (EPSn > DIL_EP_S .AND. .NOT.IS_AT_N(IJK)) THEN
	          tmpVel = tmpVel + V_s(IJKN,M)
		  counter = counter + ONE
	       ENDIF
               IF(.NOT. NO_K) THEN
                  IF (EPSb > DIL_EP_S .AND. .NOT.IS_AT_T(IJKM)) THEN
	            tmpVel = tmpVel + V_s(IJKB,M)
		    counter = counter + ONE
	          ENDIF
                  IF (EPSt > DIL_EP_S .AND. .NOT.IS_AT_T(IJK)) THEN
	            tmpVel = tmpVel + V_s(IJKT,M)
		    counter = counter + ONE
	          ENDIF
               ENDIF
! compute Us as averaged of neighbours non-dilute velocities.
               IF (counter >= ONE) THEN
	          V_S(IJK,M) = tmpVel / counter
               ELSE
	          V_S(IJK,M) = V_S(IJK,M) ! same as in source_v_s, no change in Vs values.
               ENDIF
! Normal case
            ELSE

! species density evaluated at cell faces
	       ropsN = AVG_Y(ROP_S(IJK,M),ROP_S(IJKN,M),J)
! Calculate species velocity Us at cell faces
	       V_s(IJK, M) = V_s(IJK, MMAX) + JoiY(IJK,M) / ropsN
            ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!      Finaly Compute Ws        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            IF(.NOT. NO_K) THEN
	      EPSA =  AVG_Z(EP_S(IJK,M),EP_S(IJKT,M),K)       
              IF (IP_AT_T(IJK)) THEN 
	         W_s(IJK, M) = ZERO
              ELSEIF (SIP_AT_T(IJK)) THEN 
                 ISV = IS_ID_AT_T(IJK)
	         W_s(IJK, M) = IS_VEL_S(ISV,M) 
! Dilute flow
              ELSEIF (EPSA <= DIL_EP_S) THEN 
! using the average boundary cell values to compute W_s
                 tmpVel = ZERO
	         counter = ZERO
	         IF (EPSw > DIL_EP_S .AND. .NOT.IS_AT_E(IMJK)) THEN
	            tmpVel = tmpVel + W_s(IJKW,M)
		    counter = counter + ONE
	         ENDIF
                 IF (EPSe > DIL_EP_S .AND. .NOT.IS_AT_E(IJK)) THEN
	            tmpVel = tmpVel + W_s(IJKE,M)
		    counter = counter + ONE
	         ENDIF
                 IF (EPSs > DIL_EP_S .AND. .NOT.IS_AT_N(IJMK)) THEN
	            tmpVel = tmpVel + W_s(IJKS,M)
		    counter = counter + ONE
	         ENDIF
                 IF (EPSn > DIL_EP_S .AND. .NOT.IS_AT_N(IJK)) THEN
	            tmpVel = tmpVel + W_s(IJKN,M)
		    counter = counter + ONE
	         ENDIF
                 IF (EPSb > DIL_EP_S .AND. .NOT.IS_AT_T(IJKM)) THEN
	           tmpVel = tmpVel + W_s(IJKB,M)
		   counter = counter + ONE
	         ENDIF
                 IF (EPSt > DIL_EP_S .AND. .NOT.IS_AT_T(IJK)) THEN
	           tmpVel = tmpVel + W_s(IJKT,M)
		   counter = counter + ONE
	         ENDIF
! compute Ws as averaged of neighbours non-dilute velocities.
                 IF (counter >= ONE) THEN
	            W_S(IJK,M) = tmpVel / counter
                 ELSE
	            W_S(IJK,M) = W_S(IJK,M) ! same as in source_W_s, no change in Ws values.
                 ENDIF
! Normal case
              ELSE

! species density evaluated at cell faces
	         ropsT = AVG_Z(ROP_S(IJK,M),ROP_S(IJKT,M),K) 
! Calculate species velocity Us at cell faces
	         W_s(IJK, M) = W_s(IJK, MMAX) + JoiZ(IJK,M) / ropsT
              ENDIF
            ENDIF

! if .not. fluid_at(ijk), do nothing (keep old values of velocities). 

          ENDIF     ! Fluid_at
 200  CONTINUE     ! outer IJK loop   
    ENDDO ! for m = 1, smax

      RETURN
      END SUBROUTINE UpdateSpeciesVelocities
