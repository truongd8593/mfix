!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: DRAG_gs (M, IER)                                       C
!  Purpose: Calculate the gas-solids drag coefficient                  C
!                                                                      C
!  Author: M. Syamlal                                 Date: 20-MAY-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: To use Tsuji drag correlation also if needed in DES        C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number: 2                                                  C
!  Purpose: Added switch function to Gidaspow - GIDASPOW_BLEND         C
!  Author: Charles E.A. Finney                        Date: 23-Mar-06  C
!  Reviewer: Sreekanth Pannala                        Date: 24-Mar-06  C
!                                                                      C
!  Literature/Document References:                                     C
!    Syamlal-O'Brien:                                                  C
!      Syamlal M, O'Brien TJ (1988). International Journal of          C
!        Multiphase Flow 14: 473-481.                                  C
!    Gidaspow:                                                         C
!      Ding J, Gidaspow D (1990). AIChE Journal 36: 523-538.           C
!    Gidaspow,blended (original source unknown):                       C
!      Lathouwers D, Bellan J (2000). Proceedings of the 2000 U.S. DOE C
!        Hydrogen Program Review NREL/CP-570-28890. Available from     C
!      http://www.eere.energy.gov/hydrogenandfuelcells/pdfs/28890k.pdf.C
!    Wen-Yu:                                                           C
!      Wen CY, Yu YH (1966). Chemical Engineering Progress Symposium   C
!        Series 62: 100-111.                                           C
!    Koch-Hill (modified):                                             C
!      Benyahia S, Syamlal M, O'Brien TJ (2006). Powder Technology     C
!        162: 166-174.                                                 C
!      Hill RJ, Koch DL, Ladd JC (2001). Journal of Fluid Mechanics    C
!        448: 213-241.                                                 C
!      Hill RJ, Koch DL, Ladd JC (2001). Journal of Fluid Mechanics    C
!        448: 243-278.                                                 C
!                                                                      C
!  Variables referenced: EP_g, RO_g, MU_g, D_p                         C
!  Variables modified: DRAG_gs                                         C
!                                                                      C
!  Local variables: A, B, V_rm, Re                                     C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!     
      SUBROUTINE DRAG_GS(M, IER) 
!...  Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...  Switches: -xf
!-----------------------------------------------
!     M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE fldvar
      USE geometry
      USE indices
      USE physprop
      USE run
      USE constant
      USE compar  
      USE drag  
      USE sendrecv 
      USE discretelement
      USE ur_facs 
      
      IMPLICIT NONE
!-----------------------------------------------
!     G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!     D u m m y   A r g u m e n t s
!-----------------------------------------------
!     
!     Solids phase 
      INTEGER          M 
!     
!     Error index 
      INTEGER          IER 
!     
!-----------------------------------------------
!     L o c a l   P a r a m e t e r s
!-----------------------------------------------
!     Parameters in the Cluster-effect model
!     PARAMETER (a1 = 250.)   !for G_s = 98 kg/m^2.s
!     PARAMETER (a1 = 1500.)  !for G_s = 147 kg/m^2.s
!     a1 depends upon solids flux.  It has been represented by C(1) 
!     defined in the data file.
      DOUBLE PRECISION, PARAMETER :: A2 = 0.005D0 
      DOUBLE PRECISION, PARAMETER :: A3 = 90.0D0 
      DOUBLE PRECISION, PARAMETER :: RE_C = 5.D0 
      DOUBLE PRECISION, PARAMETER :: EP_C = 0.92D0 
!-----------------------------------------------
!     L o c a l   V a r i a b l e s
!-----------------------------------------------
!     
!     Indices 
      INTEGER          I,  IJK, IMJK, IJMK, IJKM 
!     
!     Cell center value of U_sm 
      DOUBLE PRECISION USCM 
!     
!     Cell center value of U_g 
      DOUBLE PRECISION UGC 
!     
!     Cell center value of V_sm 
      DOUBLE PRECISION VSCM 
!     
!     Cell center value of V_g 
      DOUBLE PRECISION VGC 
!     
!     Cell center value of W_sm 
      DOUBLE PRECISION WSCM 
!     
!     Cell center value of W_g 
      DOUBLE PRECISION WGC 
!     
!     Magnitude of gas-solids relative velocity 
      DOUBLE PRECISION VREL 
!     
!     Reynolds number 
      DOUBLE PRECISION Re 
!     
!     Ratio of settling velocity of a multiparticle 
!     system to that of a single particle 
      DOUBLE PRECISION V_rm 
!     
!     Function of EP_g 
      DOUBLE PRECISION A 
!     
!     Function of EP_g 
      DOUBLE PRECISION B 
!     
!     Single sphere drag coefficient x Re 
      DOUBLE PRECISION C_DsxRe, C_DsxReT 
!     
!     single sphere drag coefficient 
      DOUBLE PRECISION C_d 
!     
!     drag coefficient 
      DOUBLE PRECISION DgA  

! --- Gidaspow switch function variables [ceaf 2006-03-23]
      DOUBLE PRECISION Ergun
      DOUBLE PRECISION WenYu
      DOUBLE PRECISION PHI_gs
! --- end Gidaspow switch function variables

!     
!     Gas Laminar viscosity redefined here to set
!     viscosity at pressure boundaries
      DOUBLE PRECISION Mu
!     
!     Gidaspow Reynolds number
      DOUBLE PRECISION Re_g
!
!***********************************************************
!     Declaration of variables relevant to the Koch and Hill
!     drag correlation, sof
!***********************************************************     
!     Stokes Drag Force
      DOUBLE PRECISION F_STOKES
!
!      zero Re function for low Reynolds number
       DOUBLE PRECISION F_0
!
!      inertial function for low Reynolds number
       DOUBLE PRECISION F_1
!
!      zero Re function for high Reynolds number
       DOUBLE PRECISION F_2
!
!      inertial function for high Reynolds number
       DOUBLE PRECISION F_3
!      
!      dimensionless drag force F
       DOUBLE PRECISION F
!      
!      transition Reynolds numbers
       DOUBLE PRECISION Re_Trans_1, Re_Trans_2
!      
!      solids volume fraction
       DOUBLE PRECISION phis
!      
!      weighting factor to compute F_0 and F_2
       DOUBLE PRECISION w
!     
!     Hill and Koch Reynolds number
      DOUBLE PRECISION Re_kh
!
!     End of Koch and Hill variables declaration, sof
!***********************************************************
!     
!     
!     Current value of F_gs (i.e., without underrelaxation)
      DOUBLE PRECISION F_gstmp
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'
!     
      C_DSXRET(RE) = 24D0*(1 + 0.15D0*RE**0.687D0)/RE ! Tsuji drag
      C_DSXRE(RE) = (0.63D0*SQRT(RE) + 4.8D0)**2 ! Dalla Valle (1948) 
!     C_DsxRe (Re) = 24.D0 * (1.D0 + 0.173D0 * Re**0.657D0)      ! Turton and
!     &          + 0.413D0 * Re**2.09D0 / (Re**1.09D0 + 16300.D0) ! Levenspiel (1986)
!     
!     
!     
!     !$omp  parallel do private( I,  IJK, IMJK, IJMK, IJKM, &
!     !$omp&  USCM, VSCM, WSCM, &
!     !$omp&  VREL, UGC, VGC, WGC, Re, V_rm, A, B) &
!     !$omp&  schedule(static)
      DO IJK = ijkstart3, ijkend3
         IF (FLUIDorP_FLOW_AT(IJK)) THEN 
!     
            I = I_OF(IJK) 
            IMJK = IM_OF(IJK) 
            IJMK = JM_OF(IJK) 
            IJKM = KM_OF(IJK) 
            
!     Calculate velocity components at i, j, k
            UGC = AVG_X_E(U_G(IMJK),U_G(IJK),I) 
            VGC = AVG_Y_N(V_G(IJMK),V_G(IJK)) 
            WGC = AVG_Z_T(W_G(IJKM),W_G(IJK)) 
            USCM = AVG_X_E(U_S(IMJK,M),U_S(IJK,M),I) 
            VSCM = AVG_Y_N(V_S(IJMK,M),V_S(IJK,M)) 
            WSCM = AVG_Z_T(W_S(IJKM,M),W_S(IJK,M)) 
!     
!     magnitude of gas-solids relative velocity
!     
            VREL = SQRT((UGC - USCM)**2 + (VGC - VSCM)**2 + (WGC - WSCM)**2) 
!     
!     ! Laminar viscosity at a pressure boundary is given the value of the fluid cell next to it.
!     ! This applies just to the calculation of the drag, in other routines the value of viscosity
!     ! at a pressure boundary has always a zero value.
!     
            IF (P_OUTFLOW_AT(IJK)) THEN
               IF( FLUID_AT(EAST_OF(IJK) )) THEN
                  Mu = MU_G(EAST_OF(IJK))
               ELSE IF ( FLUID_AT(WEST_OF(IJK)) ) THEN
                  Mu = MU_G(WEST_OF(IJK))
               ELSE IF ( FLUID_AT(NORTH_OF(IJK)) ) THEN
                  Mu = MU_G(NORTH_OF(IJK))
               ELSE IF ( FLUID_AT(SOUTH_OF(IJK)) ) THEN
                  Mu = MU_G(SOUTH_OF(IJK))
               ELSE IF ( FLUID_AT(TOP_OF(IJK)) ) THEN
                  Mu = MU_G(TOP_OF(IJK))
               ELSE IF ( FLUID_AT(BOTTOM_OF(IJK)) ) THEN
                  Mu = MU_G(BOTTOM_OF(IJK))
               ENDIF
            ELSE
               Mu = MU_G(IJK)
            ENDIF
            
            

!     Reynolds number
            if(Mu > ZERO)then
               RE = D_P(IJK,M)*VREL*RO_G(IJK)/Mu

!     Note the presence of gas volume fraction in ROP_G
               RE_G = D_P(IJK,M)*VREL*ROP_G(IJK)/Mu

!     Note Reynolds' number for Hill and Koch has an additional factor of 1/2 & ep_g
               RE_kh = 0.5D0*D_P(IJK,M)*VREL*ROP_G(IJK)/Mu
            else 
               RE = LARGE_NUMBER 
               RE_G = LARGE_NUMBER
               RE_kh = LARGE_NUMBER
            endif
!     
!     To select one of the following models uncomment (delete) lower
!     case c's.
!     
!---------------Begin Syamlal and O'Brien ---------------------------
!     
!     Calculate V_rm
!     
            IF(TRIM(DRAG_TYPE).EQ.'SYAM_OBRIEN') then

               IF (EP_s(IJK,M) <= ZERO) THEN 
                  F_gstmp = ZERO 
               ELSE IF (EP_G(IJK) == ZERO) THEN 
                  F_gstmp = ZERO 
               ELSE 
                  A = EP_G(IJK)**4.14D0 
                  IF (EP_G(IJK) <= 0.85D0) THEN 
                     B = drag_c1*EP_G(IJK)**1.28D0 
                  ELSE 
                     B = EP_G(IJK)**drag_d1
                  ENDIF 
                  V_RM=HALF*(A-0.06D0*RE+SQRT(3.6D-3*RE*RE+0.12D0*RE*(2.D0*B-A)+A*A)) 
!------------------Begin cluster correction --------------------------
!     uncomment the following four lines and comment the above line V_RM=... 
!     V_RM=HALF*(A-0.06D0*RE+SQRT(3.6D-3*RE*RE+0.12D0*RE*(2.D0*B-A)+A*A)) & 
!     * ( ONE + C(1) * exp( -a2*(Re - Re_c)**2 &
!     - a3*(EP_g(IJK)-ep_c)**2 &
!     )       * Re * (1. - EP_g(IJK))                )
!------------------End cluster correction ----------------------------
!     
!     Calculate the drag coefficient (Model B coeff = Model A coeff/EP_g)
!     
                  IF(TSUJI_DRAG) THEN
                     IF(EP_G(IJK).LE.0.8D0) THEN
                        F_gstmp = (Mu*EP_S(IJK,M)/(D_P(IJK,M)**2))*&
                        (150D0*(EP_S(IJK,M)/EP_G(IJK)) + 1.75D0*RE)
                     ELSE IF(EP_G(IJK).GT.0.8D0) THEN
                        IF(RE*EP_G(IJK).GT.1000D0) THEN
                           F_gstmp = 0.75D0*0.43D0*Mu*EP_S(IJK,M)*RE/(D_P(IJK,M)**2 *&
                           EP_G(IJK)**1.7D0)
                        ELSE IF(RE*EP_G(IJK).LE.1000D0) THEN
                           F_gstmp = 0.75D0*C_DSXRET(RE*EP_G(IJK))*Mu*EP_S(IJK,M)*&
                           RE/(D_P(IJK,M)**2 *EP_G(IJK)**1.7D0)
                        END IF
                     END IF 
                  ELSE IF(MODEL_B) THEN 
                     F_gstmp = 0.75D0*Mu*EP_S(IJK,M)*C_DSXRE(RE/V_RM)/(&
                     V_RM*D_P(IJK,M)*D_P(IJK,M)) 
                  ELSE
                     F_gstmp = 0.75D0*Mu*EP_S(IJK,M)*EP_G(IJK)*C_DSXRE(RE&
                     /V_RM)/(V_RM*D_P(IJK,M)*D_P(IJK,M)) 
                  ENDIF 
               ENDIF 
!---------------End Syamlal and O'Brien ---------------------------
!     
!--------------------------Begin Gidaspow --------------------------
            ELSE IF(TRIM(DRAG_TYPE).EQ.'GIDASPOW') then
               IF(EP_g(IJK) .LE. 0.8D0) THEN
                  DgA = 150D0 * (ONE - EP_g(IJK)) * Mu &
                  / ( EP_g(IJK) * D_p(IJK,M)**2 ) &
                  + 1.75D0 * RO_g(IJK) * VREL / D_p(IJK,M)
               ELSE
                  IF(Re_G .LE. 1000D0)THEN
                     C_d = (24.D0/(Re_G+SMALL_NUMBER)) * (ONE + 0.15D0 * Re_G**0.687D0)
                  ELSE
                     C_d = 0.44D0
                  ENDIF
                  DgA = 0.75D0 * C_d * VREL * ROP_g(IJK) * EP_g(IJK)**(-2.65D0) &
                  /D_p(IJK,M)
               ENDIF
               
!              Calculate the drag coefficient (Model B coeff = Model A coeff/EP_g)
               IF(Model_B)THEN
                  F_gstmp = DgA * EP_s(IJK, M)/EP_g(IJK)
               ELSE
                  F_gstmp = DgA * EP_s(IJK, M)
               ENDIF
               
!--------------------------End Gidaspow --------------------------
!     
!-----------------------Begin Gidaspow_blend ---------------------
            ELSE IF(TRIM(DRAG_TYPE).EQ.'GIDASPOW_BLEND') then
!              Dense phase - EP_g < 0.8
               Ergun = 150D0 * (ONE - EP_g(IJK)) * Mu &
               / ( EP_g(IJK) * D_p(IJK,M)**2 ) &
               + 1.75D0 * RO_g(IJK) * VREL / D_p(IJK,M)
!
!              Dilute phase - EP_g >= 0.8
               IF(Re_G .LE. 1000D0)THEN
                  C_d = (24.D0/(Re_G+SMALL_NUMBER)) * (ONE + 0.15D0 * Re_G**0.687D0)
               ELSE
                  C_d = 0.44D0
               ENDIF
               WenYu = 0.75D0 * C_d * VREL * ROP_g(IJK) * EP_g(IJK)**(-2.65D0) &
               /D_p(IJK,M)
!
!              Switch function
               PHI_gs = ATAN(150D0 * 1.75D0 * (EP_g(IJK) - 0.8D0)) / PI + 0.5D0
!              Blend the models
               DgA = (1D0 - PHI_gs) * Ergun + PHI_gs * WenYu
               
!              Calculate the drag coefficient (Model B coeff = Model A coeff/EP_g)
               IF(Model_B)THEN
                  F_gstmp = DgA * EP_s(IJK, M)/EP_g(IJK)
               ELSE
                  F_gstmp = DgA * EP_s(IJK, M)
               ENDIF
               
!-----------------------End Gidaspow_blend -----------------------
!     
!--------------------------Begin WEN_YU --------------------------
            ELSE IF(TRIM(DRAG_TYPE).EQ.'WEN_YU') then
                IF(Re_G .LE. 1000D0)THEN
                   C_d = (24.D0/(Re_G+SMALL_NUMBER)) * (ONE + 0.15D0 * Re_G**0.687D0)
                ELSE
                   C_d = 0.44D0
                ENDIF
                DgA = 0.75D0 * C_d * VREL * ROP_g(IJK) * EP_g(IJK)**(-2.65D0) &
                  /D_p(IJK,M)
               
!              Calculate the drag coefficient (Model B coeff = Model A coeff/EP_g)
               IF(Model_B)THEN
                  F_gstmp = DgA * EP_s(IJK, M)/EP_g(IJK)
               ELSE
                  F_gstmp = DgA * EP_s(IJK, M)
               ENDIF
               
!--------------------------End WEN_YU ----------------------------
!     
!--------------------Begin Koch & Hill (2001) --------------------
!     
!!!   Added by Clay Sutton (Lehigh University) 7-14-04
!!!   
!!!   MODIFICATIONS:
!!!   
!!!   1) Declared new variables F_STOKES, F_0, F_1, F_3
!!!   
!!!   2) Added new drag closure lines
!!!
!!!  Clay's implementation was modified by Sof (01-21-2005)
!!!  for a report explaining these changes contact sof@fluent.com
!
            ELSE IF(TRIM(DRAG_TYPE).EQ.'KOCH_HILL') then
!     
              F_STOKES = 18D0*MU_g(IJK)*EP_g(IJK)*EP_g(IJK)/D_p(IJK,M)**2
	       
	      phis = EP_s(IJK,M)
	      w = EXP(-10.0D0*(0.4D0-phis)/phis)
	   
	      IF(phis > 0.01D0 .AND. phis < 0.4D0) THEN
	        F_0 = (1.0D0-w) *                                           &
	              (1.0D0 + 3.0D0*dsqrt(phis/2.0D0) + 135.0D0/64.0D0*phis    &
	                *LOG(phis) + 17.14D0*phis) / (1.0D0 + 0.681D0*      &
	  	        phis - 8.48D0*phis*phis + 8.16D0*phis**3) + w *   &
			10.0D0*phis/(1.0D0-phis)**3
	               
	      ELSE IF(phis >= 0.4D0) THEN
	        F_0 = 10.0D0*phis/(1.0D0-phis)**3
	      ENDIF
	   
	      IF(phis > 0.01D0 .AND. phis <= 0.1D0) THEN
	        F_1 = dsqrt(2.0D0/phis) / 40.0D0
	      ELSE IF(phis > 0.1D0) THEN
	        F_1 = 0.11D0 + 5.1D-04 * exp(11.6D0*phis)
	      ENDIF
	   
	      IF(phis < 0.4D0) THEN
	        F_2 = (1.0D0-w) *                                           &
	              (1.0D0 + 3.0D0*dsqrt(phis/2.0D0) + 135.0D0/64.0D0*phis    &
	                *LOG(phis) + 17.89D0*phis) / (1.0D0 + 0.681D0*      &
	  	        phis - 11.03D0*phis*phis + 15.41D0*phis**3)+ w *  &
			10.0D0*phis/(1.0D0-phis)**3
	   
	      ELSE
	        F_2 = 10.0D0*phis/(1.0D0-phis)**3
	      ENDIF
	   
	      IF(phis < 0.0953D0) THEN
	        F_3 = 0.9351D0*phis + 0.03667D0
	      ELSE
	        F_3 = 0.0673D0 + 0.212D0*phis +0.0232D0/(1.0-phis)**5
	      ENDIF
	   
	      Re_Trans_1 = (F_2 - 1.0D0)/(3.0D0/8.0D0 - F_3)
	      Re_Trans_2 = (F_3 + dsqrt(F_3*F_3 - 4.0D0*F_1 &
	                 *(F_0-F_2))) / (2.0D0*F_1)
	   
	      IF(phis <= 0.01D0 .AND. Re_kh <= Re_Trans_1) THEN
	        F = 1.0D0 + 3.0D0/8.0D0*Re_kh
	   
	      ELSE IF(phis > 0.01D0 .AND. Re_kh <= Re_Trans_2) THEN
	        F = F_0 + F_1*Re_kh*Re_kh
	   
	  
	      ELSE IF(phis <= 0.01D0 .AND. Re_kh > Re_Trans_1 .OR.         &
	           phis >  0.01D0 .AND. Re_kh > Re_Trans_2) THEN
	        F = F_2 + F_3*Re_kh
	  
	      ELSE
	        F = zero
	      ENDIF
	   
!  This is a check for phis (or eps_(ijk,m)) to be within physical range
	      IF(phis <= ZERO .OR. phis > ONE) F = zero
	   
	      DgA = F * F_STOKES
!!!   
!!!   Calculate the drag coefficient (Model B coeff = Model A coeff/EP_g)
              IF(Model_B)THEN
                  F_gstmp = DgA * EP_s(IJK, M)/EP_g(IJK)
              ELSE
                  F_gstmp = DgA * EP_s(IJK, M)
              ENDIF
!     
!---------------------End Koch & Hill (2001) ----------------------
! 
            ELSE
              CALL START_LOG 
              IF(.not.DMP_LOG)call open_pe_log(ier)
	      if(mype == pe_io) WRITE (*, '(A,A)') 'Unknown DRAG_TYPE: ', DRAG_TYPE
              WRITE (UNIT_LOG, '(A,A)') 'Unknown DRAG_TYPE: ', DRAG_TYPE
              CALL END_LOG 
              call mfix_exit(myPE)  
            ENDIF
      
            F_gs(IJK, M) = (ONE - UR_F_gs) * F_gs(IJK, M) + UR_F_gs * F_gstmp
         
	 ELSE 
            F_gs(IJK, M) = ZERO 
         ENDIF 
      END DO
      
      RETURN  
      END SUBROUTINE DRAG_GS 

!//   Comments on the modifications for DMP version implementation      
!//   001 Include header file and common declarations for parallelization
!//   350 Changed do loop limits: 1,ijkmax2-> ijkstart3, ijkend3
!     
