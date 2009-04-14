!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SOURCE_GHD_GRANULAR_ENERGY(sourcelhs,sourcerhs,IJK,IER)
!  Purpose: Calculate the source terms in the granular energy equation C
!           for GHD theory                                             C
!                                                                      C
!  Author: S. Benyahia, J. Galvin, C. Hrenya        Date: 22-JAN-09    C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!     C. Hrenya handnotes and Garzo, Hrenya, Dufty papers (PRE, 2007)  C
!                                                                      C
!  Variables referenced:                                               C
!                                                                      C
!  Variables modified:                                                 C
!                                                                      C
!     Local variables: sourcelhs, sourcerhs                            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE SOURCE_GHD_GRANULAR_ENERGY(SOURCELHS, SOURCERHS, IJK, IER) 
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
      USE physprop
      USE run
      USE drag
      USE geometry
      USE fldvar
      USE visc_s
      USE ghdtheory
      USE trace
      USE indices
      USE constant
      USE compar        !//d
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
!     
!                      Indices
      INTEGER          IJK, I, J, K, IM, JM, KM, IMJK, IJMK, IJKM,&
                       IJKE, IJKW, IJKN, IJKS, IJKT, IJKB
!
!                      number densities
      DOUBLE PRECISION NiE, NiW, NiN, NiS, NiT,&
                       NiB, Nip, Ntotal
!
!                      mass of particles & non-zero gran. temp.
      DOUBLE PRECISION Mi, NonZeroTheta
!
!                      Thermal mobility-related source terms
      DOUBLE PRECISION ThermMobilityX, ThermMobilityY, ThermMobilityZ
!
!                      del.Joi and Fi.Joi terms
      DOUBLE PRECISION DelDotJoi, FiDotJoi, JoiXC, JoiYC, JoiZC
      DOUBLE PRECISION UGC, VGC, WGC, USCM, VSCM, WSCM, dragFx, dragFy, dragFz
!
!                      phase index 
      INTEGER          M, L
!
!                      Dufour-related source terms
      DOUBLE PRECISION DufourX, DufourY, DufourZ, DijQTerm, &
                       DijQTermE, DijQTermW, DijQTermN, DijQTermS, DijQTermT, DijQTermB
!                       
!                      Source terms to be kept on RHS
      DOUBLE PRECISION SOURCERHS, PressureRhs, ShearProduction, BulkViscRhs, DissDivURhs
!
!                      Source terms to be kept on LHS
      DOUBLE PRECISION SOURCELHS, PressureLhs, CollDissipation, BulkViscLhs, DissDivULhs
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------

!----------------------------------------------- 
!     Include statement functions
!-----------------------------------------------
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'ep_s2.inc'
      INCLUDE 'b_force1.inc'
      INCLUDE 'b_force2.inc'
!-----------------------------------------------

      I = I_OF(IJK) 
      J = J_OF(IJK) 
      K = K_OF(IJK) 
      IM = Im1(I)
      JM = Jm1(J)
      KM = Km1(K)
      IMJK = IM_OF(IJK)
      IJMK = JM_OF(IJK)
      IJKM = KM_OF(IJK) 
      IJKE = EAST_OF(IJK) 
      IJKW = WEST_OF(IJK)
      IJKN = NORTH_OF(IJK) 
      IJKS = SOUTH_OF(IJK)
      IJKT = TOP_OF(IJK) 
      IJKB = BOTTOM_OF(IJK) 

      NonZeroTheta = MAX(THETA_M(IJK,MMAX), SMALL_NUMBER)
      
      Ntotal = ZERO       
      DO M = 1,SMAX
          Ntotal = Ntotal + ROP_S(IJK,M)*6.d0/(Pi*D_P(IJK,M)**3 * RO_S(M))
      ENDDO          

! Production by shear: (S:grad(v))
!     p_s*tr(D)     
      PressureRhs = P_S_C(IJK,MMAX)*ZMAX(( -TRD_S_C(IJK,MMAX) ))
      PressureLhs = P_S_C(IJK,MMAX)*ZMAX((  TRD_S_C(IJK,MMAX) ))

!     mu_s*tr(D^2)   
      ShearProduction = 2.d0*MU_S_C(IJK,MMAX)*TRD_S2(IJK,MMAX)

!     lambda_s*tr(D)^2  
      BulkViscRhs = ZMAX(  LAMBDA_S_C(IJK,MMAX) ) * TRD_S_C(IJK,MMAX)**2
      BulkViscLhs = ZMAX( -LAMBDA_S_C(IJK,MMAX) ) * TRD_S_C(IJK,MMAX)**2


! Energy dissipation by collisions: (3/2)*n*T*zeta
!     (3/2)*n*T*zeta0; zeroth order cooling rate term    
      CollDissipation = 1.5d0*Ntotal*Zeta0(IJK)

!     (3/2)*n*T*zetaU*div(U) : 
      DissDivURhs = 1.5d0*Ntotal*Theta_m(IJK,MMAX)* ZMAX( -ZetaU(IJK)*TRD_S_C(IJK,MMAX) )
      DissDivULhs = 1.5d0*Ntotal* ZMAX(  ZetaU(IJK)*TRD_S_C(IJK,MMAX) )


      DufourX = ZERO
      DufourY = ZERO
      DufourZ = ZERO
      ThermMobilityX = ZERO
      ThermMobilityY = ZERO
      ThermMobilityZ = ZERO
      DelDotJoi = ZERO
      FiDotJoi  = ZERO
      DO M = 1,SMAX

! Part of heat flux: div (q)
!     Sum_ij [ div( T^2*DijQ/nj*grad(nj)) ] -> Dufour term
 
          DO L = 1,SMAX      
              Mi = (Pi/6.d0)*D_P(IJK,L)**3 * RO_S(L)
              Nip = ROP_S(IJK,L) /Mi
              NiE = ROP_S(IJKE,L)/Mi
              NiW = ROP_S(IJKW,L)/Mi
              NiN = ROP_S(IJKN,L)/Mi
              NiS = ROP_S(IJKS,L)/Mi

              DijQTerm = Theta_m(IJK,MMAX)**2 * DijQ(IJK,M,L) / Nip
              DijQTermE = AVG_X_S(DijQTerm, &
	                          Theta_m(IJKE,MMAX)**2*DijQ(IJKE,M,L)/NiE, I)
              DijQTermW = AVG_X_S(Theta_m(IJKW,MMAX)**2*DijQ(IJKW,M,L)/NiW, &
	                          DijQTerm, IM)
              DijQTermN = AVG_Y_S(DijQTerm, &
	                          Theta_m(IJKN,MMAX)**2*DijQ(IJKN,M,L)/NiN, J)
              DijQTermS = AVG_Y_S(Theta_m(IJKS,MMAX)**2*DijQ(IJKS,M,L)/NiS, &
	                          DijQTerm, JM) 
	      DufourX = DufourX + ( DijQTermE*(NiE-Nip)*ODX_E(I)*AYZ(IJK) - &
	                            DijQTermW*(Nip-NiW)*ODX_E(IM)*AYZ(IMJK) )
	      DufourY = DufourY + ( DijQTermN*(NiN-Nip)*ODY_N(J)*AXZ(IJK) - &
	                            DijQTermS*(Nip-NiS)*ODY_N(JM)*AXZ(IJMK) )

	      IF(.NOT. NO_K) THEN 
                 NiT = ROP_S(IJKT,L)/Mi
                 NiB = ROP_S(IJKB,L)/Mi
                
		 DijQTermT = AVG_Z_S(DijQTerm, &
	                          Theta_m(IJKT,MMAX)**2*DijQ(IJKT,M,L)/NiT, K)
                 DijQTermB = AVG_Z_S(Theta_m(IJKB,MMAX)**2*DijQ(IJKB,M,L)/NiB, &
	                          DijQTerm, KM)
	         
		 DufourZ = DufourZ + ( DijQTermT*(NiT-Nip)*ODZ_T(K)*OX(I)*AXY(IJK) - &
		                       DijQTermB*(Nip-NiB)*ODZ_T(KM)*OX(I)*AXY(IJKM) )
	      ENDIF

!     Sum_ij [ div( Lij*Fj) ]; thermal mobility term
!     Where Fj = Body Force
              
	       ThermMobilityX = ThermMobilityX + Mi * BFX_S(IJK,L) * ( &
	                        AVG_X_S(Lij(IJK,M,L), Lij(IJKE,M,L), I)*AYZ(IJK) - &
				AVG_X_S(Lij(IJKW,M,L), Lij(IJK,M,L), IM)*AYZ(IMJK) )
              
	       ThermMobilityY = ThermMobilityY + Mi * BFY_S(IJK,L) * ( &
	                        AVG_Y_S(Lij(IJK,M,L), Lij(IJKN,M,L), J)*AXZ(IJK) - &
				AVG_Y_S(Lij(IJKS,M,L), Lij(IJK,M,L), JM)*AXZ(IJMK) )

	      IF(.NOT. NO_K) ThermMobilityZ = ThermMobilityZ + Mi * BFZ_S(IJK,L) * ( &
	                        AVG_Z_S(Lij(IJK,M,L), Lij(IJKT,M,L), K) *AXY(IJK) -&
				AVG_Z_S(Lij(IJKB,M,L), Lij(IJK,M,L), KM) *AXY(IJKM) )

          ENDDO ! for L = 1, smax

! Additional term arising from subtraction of 3/2*T*continuity
!     + (3/2)*T* Sum_i [ div (joi/mi) ]

          Mi = (Pi/6.d0)*D_P(IJK,M)**3 * RO_S(M)
      
          DelDotJoi = DelDotJoi + 1.5d0*THETA_M(IJK,MMAX)/Mi * ( &
	              JoiX(IJK,M)*AYZ(IJK) - JoiX(IMJK,M)*AYZ(IMJK) + JoiY(IJK,M)*&
                      AXZ(IJK) - JoiY(IJMK,M)*AXZ(IJMK) + JoiZ(IJK,M)*AXY(IJK) - JoiZ(&
                      IJKM,M)*AXY(IJKM) )



! Species force dot species mass flux 
!     Sum_i [ Fi dot joi/mi ]  
            
!     Calculate species mass flux components at cell center
          JoiXC = AVG_X_E(JoiX(IMJK,M),JoiX(IJK,M),I) 
          JoiYC = AVG_Y_N(JoiY(IJMK,M),JoiY(IJK,M)) 
          JoiZC = AVG_Z_T(JoiZ(IJKM,M),JoiZ(IJK,M))

! drag force on a particle
          UGC = AVG_X_E(U_G(IMJK),U_G(IJK),I)
          VGC = AVG_Y_N(V_G(IJMK),V_G(IJK)) 
	  WGC = AVG_Z_T(W_G(IJKM),W_G(IJK))
	  USCM = AVG_X_E(U_S(IMJK,M),U_S(IJK,M),I)
	  VSCM = AVG_Y_N(V_S(IJMK,M),V_S(IJK,M))
	  WSCM = AVG_Z_T(W_S(IJKM,M),W_S(IJK,M))
	  
	  dragFx = F_GS(IJK ,M)/ROP_S(IJK,M) * (UGC - USCM)
	  dragFy = F_GS(IJK ,M)/ROP_S(IJK,M) * (VGC - VSCM)
	  dragFz = F_GS(IJK ,M)/ROP_S(IJK,M) * (WGC - WSCM)
	  
	  FiDotJoi  = FiDotJoi  + JoiXC*BFX_S(IJK,M) + JoiYC*BFY_S(IJK,M) + JoiZC*BFZ_S(IJK,M)
          


      ENDDO ! for M = 1, smax


      SOURCERHS = (PressureRhs + ShearProduction + BulkViscRhs + DissDivURhs)*VOL(IJK) &
                 + ZMAX(DufourX)+ZMAX(DufourY)+ZMAX(DufourZ) &
                 + ZMAX(ThermMobilityX)+ZMAX(ThermMobilityY)+ZMAX(ThermMobilityZ)      &
		 + ZMAX(DelDotJoi) + ZMAX(FiDotJoi)*VOL(IJK)

      SOURCELHS = ( (PressureLhs + BulkViscLhs)/NonZeroTheta   + &
                  (CollDissipation + DissDivULhs) ) * VOL(IJK) + &
		  ( ZMAX(-DufourX)+ZMAX(-DufourY)+ZMAX(-DufourZ)    + &
                 + ZMAX(-ThermMobilityX)+ZMAX(-ThermMobilityY)+ZMAX(-ThermMobilityZ)   &
		 + ZMAX(-DelDotJoi) + ZMAX(-FiDotJoi)*VOL(IJK) )/ NonZeroTheta
      
      RETURN  
      END SUBROUTINE SOURCE_GHD_GRANULAR_ENERGY 
