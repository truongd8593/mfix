!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: DRAG_gs                                                 C
!  Purpose: Calculate the gas-solids drag coefficient                  C
!                                                                      C
!  Author: M. Syamlal                                 Date: 20-MAY-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
     
      SUBROUTINE DRAG_GS(M, IER) !DISCRETE_FLAG, IER) 

!-----------------------------------------------
! Modules
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
      USE funits
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Solids phase index
      INTEGER, INTENT(IN) :: M 
! Error index 
      INTEGER, INTENT(INOUT) :: IER 
! Flag used only when the hybrid model is invoked and notifies the
! routine that the solid phase index M refers to the indice of a
! discrete 'phase' not a continuous phase so that the appropriate
! variables are referenced.
!      LOGICAL, INTENT(IN) :: DISCRETE_FLAG
       LOGICAL :: DISCRETE_FLAG      
!-----------------------------------------------
! Local parameters
!-----------------------------------------------

!-----------------------------------------------
! Local variables
!-----------------------------------------------
! indices 
      INTEGER :: I, IJK, IMJK, IJMK, IJKM
! cell center value of U_sm, V_sm, W_sm
      DOUBLE PRECISION :: USCM, VSCM, WSCM
! cell center value of x, y and z-particle velocity 
      DOUBLE PRECISION :: USCM_HYS, VSCM_HYS, WSCM_HYS
! cell center value of U_g, V_g, W_g
      DOUBLE PRECISION :: UGC, VGC, WGC
! magnitude of gas-solids relative velocity 
      DOUBLE PRECISION :: VREL 
! gas laminar viscosity redefined here to set viscosity at pressure
! boundaries
      DOUBLE PRECISION :: Mu
! drag coefficient 
      DOUBLE PRECISION :: DgA  
! current value of F_gs (i.e., without underrelaxation)
      DOUBLE PRECISION F_gstmp 
! indices of solids phases (continuous, discrete)
      INTEGER :: CM, DM, L
! temporary shift of total number of solids phases to account for both
! discrete and continuous solids phases used for the hybrid mdoel
      INTEGER :: MAXM
! tmp local variable for the particle diameter of solids 
! phase M (continuous or discrete)
      DOUBLE PRECISION :: DP_loc(2*DIM_M)
! tmp local variable for the solids volume fraction of solids
! phase M (continuous or discrete)
      DOUBLE PRECISION :: EPs_loc(2*DIM_M)
! correction factors for implementing polydisperse drag model 
! proposed by van der Hoef et al. (2005)
      DOUBLE PRECISION :: F_cor, tmp_sum, tmp_fac
! average particle diameter in polydisperse systems
      DOUBLE PRECISION :: DPA
! diameter ratio in polydisperse systems
      DOUBLE PRECISION :: Y_i
! total solids volume fraction
      DOUBLE PRECISION :: phis
! temporary local variables to use for dummy arguments in subroutines
! void fraction, gas density, gas bulk density, solids volume fraction
! particle diameter      
      DOUBLE PRECISION :: EPG, ROg, ROPg, EP_SM, DPM
!-----------------------------------------------
! Include statement functions
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'
!-----------------------------------------------

!!$omp  parallel do private( I,  IJK, IMJK, IJMK, IJKM, &
!!$omp&  USCM, VSCM, WSCM, &
!!$omp&  VREL, UGC, VGC, WGC, Re, V_rm, A, B) &
!!$omp&  schedule(static)

      discrete_flag =.false.

      DO IJK = ijkstart3, ijkend3

         IF (FLUIDorP_FLOW_AT(IJK)) THEN 

            I = I_OF(IJK) 
            IMJK = IM_OF(IJK) 
            IJMK = JM_OF(IJK) 
            IJKM = KM_OF(IJK) 

! assign local variables: DP_loc, EPs_loc, and MAXM.
! These represent arrays for the particle diameter, solids volume
! fraction, and number of particle types                 
            IF (.NOT.DES_CONTINUUM_HYBRID) THEN
               MAXM = SMAX
               DO CM = 1,SMAX
                  DP_loc(CM) = D_p(IJK,CM)
                  EPs_loc(CM) = EP_S(IJK,CM)
               ENDDO
            ELSE   ! des_continuum_hybrid branch
! for the hybrid model store the diameters of both discrete and
! continuous phases in a single new local variable.  likewise with the
! solids volume fraction.  also change the loop limit to include all
! solids phases (discrete+continuum)
               MAXM = SMAX + DES_MMAX
               IF (DISCRETE_FLAG) THEN
! the subroutine is being called to calculate the drag coefficient on
! a discrete particle 'phase' so populate DP, EPS starting with discrete
! phases
                  DO DM = 1,DES_MMAX
                     DP_loc(DM) = DES_D_p0(DM)         
                     EPs_loc(DM) = DES_ROP_S(IJK,DM)/DES_RO_S(DM)
                  ENDDO
                  L=DES_MMAX 
                  DO CM = 1,SMAX
                     L = L + 1
                     DP_loc(L) = D_P(IJK,CM)
                     EPs_loc(L) = EP_S(IJK,CM)
                  ENDDO
               ELSE
! the subroutine is being called to calculate the drag coefficient on
! a continous solids phase so populate DP, EPS starting with continuous
! phases
                  DO CM = 1,SMAX
                     DP_loc(CM) = D_p(IJK,CM)         
                     EPs_loc(CM) = EP_S(IJK,CM)
                  ENDDO
                  L=SMAX 
                  DO DM = 1,DES_MMAX
                     L = L + 1
                     DP_loc(L) = DES_D_p0(DM)
                     EPs_loc(L) = DES_ROP_S(IJK,DM)/DES_RO_S(DM)
                  ENDDO          
               ENDIF
            ENDIF   ! end if/else (.not.des_continuum_hybrid)

! Calculate velocity components at i, j, k
            UGC = AVG_X_E(U_G(IMJK),U_G(IJK),I) 
            VGC = AVG_Y_N(V_G(IJMK),V_G(IJK)) 
            WGC = AVG_Z_T(W_G(IJKM),W_G(IJK)) 
            IF(DES_CONTINUUM_COUPLED) THEN
               USCM = DES_U_S(IJK,M)
               VSCM = DES_V_S(IJK,M)
               WSCM = DES_W_S(IJK,M)
            ELSE               
               USCM = AVG_X_E(U_S(IMJK,M),U_S(IJK,M),I) 
               VSCM = AVG_Y_N(V_S(IJMK,M),V_S(IJK,M)) 
               WSCM = AVG_Z_T(W_S(IJKM,M),W_S(IJK,M))
            ENDIF
! magnitude of gas-solids relative velocity
            VREL = SQRT((UGC - USCM)**2 + (VGC - VSCM)**2 + (WGC - WSCM)**2) 

! Laminar viscosity at a pressure boundary is given the value of the
! fluid cell next to it. This applies just to the calculation of the
! drag, in other routines the value of viscosity at a pressure boundary
! always has a zero value.
            IF (P_FLOW_AT(IJK)) THEN
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

! calculate the total solids volume fraction
            phis = ZERO
            DO L = 1, MAXM
! this is slightly /= one-ep_g due to round-off
               phis = phis + EPs_loc(L)
            ENDDO 

! calculate the average paricle diameter and particle ratio            
            DPA = ZERO
            tmp_sum = ZERO
            tmp_fac = ZERO
            DO L = 1, MAXM
               IF (phis .GT. ZERO) THEN
                  tmp_fac = EPs_loc(L)/phis
                  tmp_sum = tmp_sum + tmp_fac/DP_loc(L)
                ELSE
                  tmp_sum = tmp_sum + ONE/DP_loc(L) ! not important, but will avoid NaN's in empty cells
                ENDIF
            ENDDO 
            DPA = ONE / tmp_sum
            Y_i = DP_loc(M)/DPA

! assign variables for short dummy arguments
            EPg = EP_G(IJK)
            ROg = RO_G(IJK)
            ROPg = ROP_G(IJK)
            EP_SM = EPs_loc(M)
            DPM = DP_loc(M)

! determine the drag coefficient
            SELECT CASE(TRIM(DRAG_TYPE))
            CASE ('SYAM_OBRIEN')
               CALL DRAG_SYAM_OBRIEN(DgA,EPG,Mu,ROg,VREL,&
                   DPM)
            CASE ('GIDASPOW') 
               CALL DRAG_GIDASPOW(DgA,EPg,Mu,ROg,ROPg,VREL,&
                    DPM)
            CASE ('GIDASPOW_PCF')
               CALL DRAG_GIDASPOW(DgA,EPg,Mu,ROg,ROPg,VREL,&
                    DPA)
            CASE ('GIDASPOW_BLEND')
               CALL DRAG_GIDASPOW_BLEND(DgA,EPg,Mu,ROg,ROPg,VREL,&
                    DPM)
            CASE ('GIDASPOW_BLEND_PCF')
               CALL DRAG_GIDASPOW_BLEND(DgA,EPg,Mu,ROg,ROPg,VREL,&
                    DPA)
            CASE ('WEN_YU')
               CALL DRAG_WEN_YU(DgA,EPg,Mu,ROPg,VREL,&
                    DPM)
            CASE ('WEN_YU_PCF')
               CALL DRAG_WEN_YU(DgA,EPg,Mu,ROPg,VREL,&
                    DPA)
            CASE ('KOCH_HILL')
               CALL DRAG_KOCH_HILL(DgA,EPg,Mu,ROPg,VREL,&
                    DPM,DPM,phis)
            CASE ('KOCH_HILL_PCF')
               CALL DRAG_KOCH_HILL(DgA,EPg,Mu,ROPg,VREL,&
                    DPM,DPA,phis)
            CASE ('BVK')
               CALL DRAG_BVK(DgA,EPg,Mu,ROPg,VREL,&
                    DPM,DPA,phis)
            CASE ('HYS')
! calculate velocity components of each solids phase
               USCM_HYS = ZERO
               VSCM_HYS = ZERO
               WSCM_HYS = ZERO
               IF(phis > ZERO) THEN
                  DO L = 1, MAXM
                     USCM_HYS = USCM_HYS + EPs_loc(L)*(UGC - &
                        AVG_X_E(U_S(IMJK,L),U_S(IJK,L),I))
                     VSCM_HYS = VSCM_HYS + EPs_loc(L)*(VGC - &
                        AVG_Y_N(V_S(IJMK,L),V_S(IJK,L)))
                     WSCM_HYS = WSCM_HYS + EPs_loc(L)*(WGC - &
                        AVG_Z_T(W_S(IJKM,L),W_S(IJK,L)))
                  ENDDO 
                  USCM_HYS = USCM_HYS/phis
                  VSCM_HYS = VSCM_HYS/phis
                  WSCM_HYS = WSCM_HYS/phis
               ENDIF
! magnitude of gas-solids relative velocity
               VREL = SQRT(USCM_HYS**2 +VSCM_HYS**2 +WSCM_HYS**2)
               CALL DRAG_HYS(DgA,EPg,Mu,ROPg,VREL,&
                    DP_loc(:),DPA,Y_i,EPs_loc(:),phis,M,IJK,MAXM)
            CASE DEFAULT
               CALL START_LOG 
               IF(.NOT.DMP_LOG) call open_pe_log(ier)
               IF(DMP_LOG) WRITE (*, '(A,A)') &
                  'Unknown DRAG_TYPE: ', DRAG_TYPE
               WRITE (UNIT_LOG, '(A,A)') 'Unknown DRAG_TYPE: ', DRAG_TYPE
               CALL END_LOG 
               CALL mfix_exit(myPE)  
            END SELECT   ! end selection of drag_type


! Modify drag coefficient to account for possible corrections and
! for differences between Model B and Model A
            IF(TRIM(DRAG_TYPE) == 'HYS') THEN
! this drag model is handled differently than the others
               IF(Model_B)THEN
                  F_gstmp = DgA/EPg
               ELSE
                  F_gstmp = DgA
               ENDIF
            ELSE
               IF (EP_SM <= ZERO) THEN 
                  DgA = ZERO 
               ELSEIF (EPg == ZERO) THEN  
! this case will already be caught in most drag subroutines whenever
! RE==0 (for correlations in which RE includes EPg)
                  DgA = ZERO   
               ENDIF 
               IF(TRIM(DRAG_TYPE) == 'GIDASPOW_PCF' .OR. &
                  TRIM(DRAG_TYPE) == 'GIDASPOW_BLEND_PCF' .OR. &
                  TRIM(DRAG_TYPE) == 'WEN_YU_PCF' .OR. &
                  TRIM(DRAG_TYPE) == 'KOCH_HILL_PCF' .OR. &
                  TRIM(DRAG_TYPE) == 'BVK') THEN
! see erratum by Beetstra et al. (2007) : the correction factor differs
! for model A versus model B. 
! application of the correction factor for model A is found from
! the correction factor for model B and neglects the Y_i**3 term
                  IF(Model_B) THEN
                     IF (M == 1) THEN
                        F_cor = (EPg*Y_i + phis*Y_i**2) 
                     ELSE
                        F_cor = (EPg*Y_i + phis*Y_i**2 + &
                           0.064d0*EPg*Y_i**3) 
                     ENDIF
                  ELSE
                     F_cor = Y_i
                  ENDIF
                  DgA = ONE/(Y_i*Y_i) * DgA * F_cor 
               ENDIF
! Calculate the drag coefficient (Model B coeff = Model A coeff/EP_g)
               IF(Model_B)THEN
                  F_gstmp = DgA * EP_SM/EPg
               ELSE
                  F_gstmp = DgA * EP_SM
               ENDIF
            ENDIF   !end if/else trim(drag_type=='hys')

            F_gs(IJK, M) = (ONE - UR_F_gs)*F_gs(IJK,M) + &
               UR_F_gs*F_gstmp

            IF(TRIM(KT_TYPE) == 'GHD') THEN
              IF(M==1) THEN
                F_gs(IJK,MMAX) = F_gs(IJK,M)
              ELSE
                F_gs(IJK,MMAX) = F_gs(IJK,MMAX) + F_gs(IJK,M)
              ENDIF
            ENDIF
         
         ELSE   ! .not.(fluidorp_flow_at(ijk)) branch

            F_gs(IJK,M) = ZERO 
            IF(TRIM(KT_TYPE) == 'GHD') F_gs(IJK, MMAX) = ZERO

         ENDIF   ! end if (fluidorp_flow_at(ijk))

      ENDDO   ! end do (ijk=ijkstart3,ijkend3)

      RETURN  
      END SUBROUTINE DRAG_GS 



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Function(s): C_DsxRe                                                C
!  Purpose:                                                            C
!     Calculate single sphere drag correlation multiplied by           C
!     the Reynolds number or                                           C
!     Calculate the single sphere drag correlation                     C      
!                                                                      C
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC

! Dalla Valle (1948)
!----------------------------------------------------------------->>>
      DOUBLE PRECISION FUNCTION C_DSXRE_DV(RE)
      USE param
      USE param1
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: RE ! Reynolds number 

      C_DSXRE_DV = (0.63D0*SQRT(RE) + 4.8D0)**2
      RETURN
      END FUNCTION C_DSXRE_DV
!-----------------------------------------------------------------<<<

! Turton and Levenspiel (1986)
!----------------------------------------------------------------->>>
      DOUBLE PRECISION FUNCTION C_DSXRE_TL(RE)
      USE param
      USE param1
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: RE ! Reynolds number 

      C_DSXRE_TL = 24.D0*(1.D0 + 0.173D0*RE**0.657D0) + &
         0.413D0*RE**2.09D0/(RE**1.09D0 + 16300.D0)
      RETURN
      END FUNCTION C_DSXRE_TL
!-----------------------------------------------------------------<<<


! Schiller and Naumann (1933)
!----------------------------------------------------------------->>>
      DOUBLE PRECISION FUNCTION C_DS_SN(RE)
      USE param
      USE param1
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: RE ! Reynolds number 

      C_DS_SN = 24.D0*(1.D0 + 0.15D0*RE**0.687D0)/(RE+SMALL_NUMBER)
      RETURN
      END FUNCTION C_DS_SN
!-----------------------------------------------------------------<<<



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: DRAG_SYAM_OBRIEN                                        C
!  Purpose: Calculate the gas-solids drag coefficient                  C
!                                                                      C
!  Literature/Document References:                                     C
!     Syamlal M, O'Brien TJ (1988). International Journal of           C
!        Multiphase Flow 14: 473-481.                                  C
!                                                                      C
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC

      SUBROUTINE DRAG_SYAM_OBRIEN(lDgA,EPg,Mug,ROg,VREL,&
                 DPM)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param 
      USE param1
      USE constant
      use run
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! drag coefficient
      DOUBLE PRECISION, INTENT(OUT) :: ldGA
! gas volume fraction 
      DOUBLE PRECISION, INTENT(IN) :: EPg
! gas laminar viscosity 
      DOUBLE PRECISION, INTENT(IN) :: Mug
! gas density
      DOUBLE PRECISION, INTENT(IN) :: ROg
! Magnitude of gas-solids relative velocity 
      DOUBLE PRECISION, INTENT(IN) :: VREL
! particle diameter of solids phase M
      DOUBLE PRECISION, INTENT(IN) :: DPM
!-----------------------------------------------
! Local parameters
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
! Local variables
!-----------------------------------------------
! Variables which are function of EP_g
      DOUBLE PRECISION :: A, B
! Ratio of settling velocity of a multiparticle system to 
! that of a single particle
      DOUBLE PRECISION :: V_rm 
! Reynolds number 
      DOUBLE PRECISION :: RE
!-----------------------------------------------
! External functions
!-----------------------------------------------
! Single sphere drag coefficient x Re 
      DOUBLE PRECISION, EXTERNAL :: C_DSXRE_DV
!-----------------------------------------------
     
      IF(Mug > ZERO) THEN
         RE = DPM*VREL*ROg/Mug
      ELSE
         RE = LARGE_NUMBER 
      ENDIF

! Calculate V_rm
      A = EPg**4.14D0 
      IF (EPg <= 0.85D0) THEN 
         B = drag_c1*EPg**1.28D0 
      ELSE 
        B = EPg**drag_d1
      ENDIF

      V_RM=HALF*(A-0.06D0*RE+&
           SQRT((3.6D-3)*RE*RE+0.12D0*RE*(2.D0*B-A)+A*A) )

!------------------Begin cluster correction --------------------------
! uncomment the following four lines ... 
!       V_RM = V_RM * (ONE + C(1)*&
!                      EXP(-A2*(RE-RE_C)**2 - A3*(EPg-EP_C)**2)* &
!                      RE*(1. - EPg)) 
!------------------End cluster correction ----------------------------

      lDgA = 0.75D0*Mug*EPg*C_DSXRE_DV(RE/V_RM)/(V_RM*DPM*DPM)
    
      IF (RE == ZERO) lDgA = ZERO 

      RETURN
      END SUBROUTINE DRAG_SYAM_OBRIEN

      
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: DRAG_GIDASPOW                                           C
!  Purpose: Calculate the gas-solids drag coefficient                  C
!                                                                      C
!  Literature/Document References:                                     C
!     Ding J, Gidaspow D (1990). AIChE Journal 36: 523-538.            C
!                                                                      C
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC

      SUBROUTINE DRAG_GIDASPOW(lDgA,EPg,Mug,ROg,ROPg,VREL,&
                 DPM)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param 
      USE param1
      USE constant
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! drag coefficient
      DOUBLE PRECISION, INTENT(OUT) :: lDgA
! gas volume fraction 
      DOUBLE PRECISION, INTENT(IN) :: EPg
! gas laminar viscosity 
      DOUBLE PRECISION, INTENT(IN) :: Mug
! gas density
      DOUBLE PRECISION, INTENT(IN) :: ROg
! gas density*EP_g
      DOUBLE PRECISION, INTENT(IN) :: ROPg      
! Magnitude of gas-solids relative velocity 
      DOUBLE PRECISION, INTENT(IN) :: VREL
! particle diameter of solids phase M or
! average particle diameter if PCF
      DOUBLE PRECISION, INTENT(IN) :: DPM
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Reynolds number 
      DOUBLE PRECISION :: RE
! Single sphere drag coefficient 
      DOUBLE PRECISION :: C_d 
!-----------------------------------------------

      IF(Mug > ZERO) THEN
! Note the presence of gas volume fraction in ROPG
         RE = DPM*VREL*ROPg/Mug
      ELSE
         RE = LARGE_NUMBER 
      ENDIF

! Dense phase 
      IF(EPg <= 0.8D0) THEN
         lDgA = 150D0*(ONE-EPg)*Mug / (EPg*DPM**2) + &
                1.75D0*ROg*VREL/DPM
      ELSE
! Dilute phase - EP_g >= 0.8
         IF(RE <= 1000D0)THEN
! this could be replaced with the function C_DS_SN 
            C_d = (24.D0/(RE+SMALL_NUMBER)) * &
                  (ONE + 0.15D0 * RE**0.687D0)
         ELSE
            C_d = 0.44D0
         ENDIF
         lDgA = 0.75D0*C_d*VREL*ROPg*EPg**(-2.65D0) / DPM
      ENDIF

      IF (RE == ZERO) lDgA = ZERO      

      RETURN
      END SUBROUTINE DRAG_GIDASPOW



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: DRAG_GIDASPOW_BLEND                                     C
!  Purpose: Calculate the gas-solids drag coefficient                  C
!                                                                      C
!  Author: Charles E.A. Finney                        Date: 23-Mar-06  C
!  Reviewer: Sreekanth Pannala                        Date: 24-Mar-06  C
!                                                                      C
!  Literature/Document References:                                     C
!     original source unknown:                                         C
!     Lathouwers D, Bellan J (2000). Proceedings of the 2000 U.S. DOE  C
!        Hydrogen Program Review NREL/CP-570-28890. Available from     C
!     http://www.eere.energy.gov/hydrogenandfuelcells/pdfs/28890k.pdf. C
!                                                                      C
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC

      SUBROUTINE DRAG_GIDASPOW_BLEND(lDgA,EPg,Mug,ROg,ROPg,VREL,&
                 DPM)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param 
      USE param1
      USE constant
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! drag coefficient
      DOUBLE PRECISION, INTENT(OUT) :: lDgA
! gas volume fraction 
      DOUBLE PRECISION, INTENT(IN) :: EPg
! gas laminar viscosity 
      DOUBLE PRECISION, INTENT(IN) :: Mug
! gas density
      DOUBLE PRECISION, INTENT(IN) :: ROg
! gas density*EP_g
      DOUBLE PRECISION, INTENT(IN) :: ROPg      
! Magnitude of gas-solids relative velocity 
      DOUBLE PRECISION, INTENT(IN) :: VREL
! particle diameter of solids phase M or
! average particle diameter if PCF
      DOUBLE PRECISION, INTENT(IN) :: DPM
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Reynolds number 
      DOUBLE PRECISION :: RE
! Single sphere drag coefficient 
      DOUBLE PRECISION :: C_d 
! Gidaspow switch function variables
      DOUBLE PRECISION :: Ergun
      DOUBLE PRECISION :: WenYu
      DOUBLE PRECISION :: PHI_gs      
!-----------------------------------------------

      IF(Mug > ZERO) THEN
! Note the presence of gas volume fraction in ROPG
         RE = DPM*VREL*ROPg/Mug
      ELSE
         RE = LARGE_NUMBER 
      ENDIF

! Dense phase - EP_g <= 0.8
      Ergun = 150D0*(ONE-EPg)*Mug / (EPg*DPM**2) + &
               1.75D0*ROg*VREL/DPM
! Dilute phase - EP_g >= 0.8
      IF(RE <= 1000D0)THEN
         C_d = (24.D0/(RE+SMALL_NUMBER)) * &
           (ONE + 0.15D0 * RE**0.687D0)
      ELSE
         C_d = 0.44D0
      ENDIF
      WenYu = 0.75D0*C_d*VREL*ROPg*EPg**(-2.65D0) / DPM

! Switch function
      PHI_gs = ATAN(150.D0*1.75D0*(EPg - 0.8D0))/PI + 0.5D0

! Blend the models
      lDgA = (1.D0-PHI_gs)*Ergun + PHI_gs*WenYu
      IF (RE == ZERO) lDgA = ZERO

      RETURN
      END SUBROUTINE DRAG_GIDASPOW_BLEND



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: DRAG_WEN_YU                                             C
!  Purpose: Calculate the gas-solids drag coefficient                  C
!                                                                      C
!  Literature/Document References:                                     C
!     Wen CY, Yu YH (1966). Chemical Engineering Progress Symposium    C
!        Series 62: 100-111.                                           C
!                                                                      C
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC

      SUBROUTINE DRAG_WEN_YU(lDgA,EPg,Mug,ROPg,VREL,&
                 DPM)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param 
      USE param1
      USE constant
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! drag coefficient
      DOUBLE PRECISION, INTENT(OUT) :: lDgA
! gas volume fraction 
      DOUBLE PRECISION, INTENT(IN) :: EPg
! gas laminar viscosity 
      DOUBLE PRECISION, INTENT(IN) :: Mug
! gas density*EP_g
      DOUBLE PRECISION, INTENT(IN) :: ROPg      
! Magnitude of gas-solids relative velocity 
      DOUBLE PRECISION, INTENT(IN) :: VREL
! particle diameter of solids phase M or
! average particle diameter if PCF
      DOUBLE PRECISION, INTENT(IN) :: DPM
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Reynolds number 
      DOUBLE PRECISION :: RE
! Single sphere drag coefficient 
      DOUBLE PRECISION :: C_d 
!-----------------------------------------------

      IF(Mug > ZERO) THEN
! Note the presence of gas volume fraction in ROPG
         RE = DPM*VREL*ROPg/Mug
      ELSE
         RE = LARGE_NUMBER 
      ENDIF      

      IF(RE <= 1000D0)THEN
         C_d = (24.D0/(RE+SMALL_NUMBER)) * &
               (ONE + 0.15D0*RE**0.687D0)
      ELSE
         C_d = 0.44D0
      ENDIF

      lDgA = 0.75D0 * C_d*VREL*ROPg*EPg**(-2.65D0) / DPM
      IF (RE == ZERO) lDgA = ZERO
               
      RETURN
      END SUBROUTINE DRAG_WEN_YU



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: DRAG_KOCH_HILL                                          C
!  Purpose: Calculate the gas-solids drag coefficient                  C
!                                                                      C
!  Author: Clay Sutton (Lehigh University)            Date: 14-Jul-04  C
!                                                                      C
!  Revision: 1                                                         C
!  Author: Sofiane Benyahia                           Date: 21-Jan-05  C
!                                                                      C
!  Literature/Document References:                                     C
!     Benyahia S, Syamlal M, O'Brien TJ (2006). Powder Technology      C
!        162: 166-174.                                                 C
!     Hill RJ, Koch DL, Ladd JC (2001). Journal of Fluid Mechanics     C
!        448: 213-241.                                                 C
!     Hill RJ, Koch DL, Ladd JC (2001). Journal of Fluid Mechanics     C
!        448: 243-278.                                                 C
!                                                                      C
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC

      SUBROUTINE DRAG_KOCH_HILL(lDgA,EPg,Mug,ROPg,VREL,&
                 DPM,DPA,PHIS)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param 
      USE param1
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! drag coefficient
      DOUBLE PRECISION, INTENT(OUT) :: lDgA
! gas volume fraction 
      DOUBLE PRECISION, INTENT(IN) :: EPg
! gas laminar viscosity 
      DOUBLE PRECISION, INTENT(IN) :: Mug
! gas density*EP_g
      DOUBLE PRECISION, INTENT(IN) :: ROPg      
! Magnitude of gas-solids relative velocity 
      DOUBLE PRECISION, INTENT(IN) :: VREL
! particle diameter of solids phase M 
      DOUBLE PRECISION, INTENT(IN) :: DPM
! average particle diameter if pcf otherwise DPM again
      DOUBLE PRECISION, INTENT(IN) :: DPA
! total solids volume fraction of solids phases 
      DOUBLE PRECISION, INTENT(IN) :: PHIS      
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Reynolds number 
      DOUBLE PRECISION :: RE
! transition Reynolds numbers
      DOUBLE PRECISION :: Re_Trans_1, Re_Trans_2
! Stokes Drag Force
      DOUBLE PRECISION :: F_STOKES
! zero Re function for low Reynolds number
      DOUBLE PRECISION :: F_0
! inertial function for low Reynolds number
      DOUBLE PRECISION :: F_1
! zero Re function for high Reynolds number
      DOUBLE PRECISION :: F_2
! inertial function for high Reynolds number
      DOUBLE PRECISION :: F_3
! dimensionless drag force F
      DOUBLE PRECISION :: F
! weighting factor to compute F_0 and F_2
      DOUBLE PRECISION :: w
!-----------------------------------------------


      IF(Mug > ZERO) THEN
! Note the presence of gas volume fraction in ROPG and factor of 1/2
         RE = 0.5D0*DPA*VREL*ROPG/Mug        ! if pcf DPA otherwise DPM
      ELSE
         RE = LARGE_NUMBER 
      ENDIF      

      F_STOKES = 18.D0*Mug*EPg*EPg/DPM**2    ! use DPM
      w = EXP(-10.0D0*(0.4D0-phis)/phis)

      IF(phis > 0.01D0 .AND. phis < 0.4D0) THEN
         F_0 = (1.0D0-w) * (1.0D0 + 3.0D0*dsqrt(phis/2.0D0) + &
            135.0D0/64.0D0*phis*LOG(phis) + 17.14D0*phis) / &
            (1.0D0 + 0.681D0*phis - 8.48D0*phis*phis + &
            8.16D0*phis**3) + w*10.0D0*phis/(1.0D0-phis)**3
      ELSEIF(phis >= 0.4D0) THEN
         F_0 = 10.0D0*phis/(1.0D0-phis)**3
      ENDIF

      IF(phis > 0.01D0 .AND. phis <= 0.1D0) THEN
        F_1 = dsqrt(2.0D0/phis) / 40.0D0
      ELSE IF(phis > 0.1D0) THEN
        F_1 = 0.11D0 + 5.1D-04 * exp(11.6D0*phis)
      ENDIF

      IF(phis < 0.4D0) THEN
        F_2 = (1.0D0-w) * (1.0D0 + 3.0D0*dsqrt(phis/2.0D0) + &
           135.0D0/64.0D0*phis*LOG(phis) + 17.89D0*phis) / &
           (1.0D0 + 0.681D0*phis - 11.03D0*phis*phis + &
           15.41D0*phis**3)+ w*10.0D0*phis/(1.0D0-phis)**3
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

      IF(phis <= 0.01D0 .AND. RE <= Re_Trans_1) THEN
         F = 1.0D0 + 3.0D0/8.0D0*RE
      ELSEIF(phis > 0.01D0 .AND. RE <= Re_Trans_2) THEN
         F = F_0 + F_1*RE*RE
      ELSEIF(phis <= 0.01D0 .AND. RE > Re_Trans_1 .OR.   &
         phis >  0.01D0 .AND. RE > Re_Trans_2) THEN
         F = F_2 + F_3*RE
      ELSE
         F = zero
      ENDIF

      lDgA = F * F_STOKES
      IF (RE == ZERO) lDgA = ZERO
  
      RETURN
      END SUBROUTINE DRAG_KOCH_HILL
       


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: DRAG_BVK                                                C
!  Purpose: Calculate the gas-solids drag coefficient                  C
!                                                                      C
!  Literature/Document References:                                     C
!     Beetstra, van der Hoef, Kuipers, Chem. Eng. Science 62           C
!     (Jan 2007)                                                       C
!                                                                      C
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC

      SUBROUTINE DRAG_BVK(lDgA,EPg,Mug,ROPg,VREL,&
                 DPM,DPA,PHIS)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param 
      USE param1
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! drag coefficient
      DOUBLE PRECISION, INTENT(OUT) :: lDgA
! gas volume fraction 
      DOUBLE PRECISION, INTENT(IN) :: EPg
! gas laminar viscosity 
      DOUBLE PRECISION, INTENT(IN) :: Mug
! gas density*EP_g
      DOUBLE PRECISION, INTENT(IN) :: ROPg      
! magnitude of gas-solids relative velocity 
      DOUBLE PRECISION, INTENT(IN) :: VREL
! particle diameter of solids phase M or
      DOUBLE PRECISION, INTENT(IN) :: DPM
! average particle diameter 
      DOUBLE PRECISION, INTENT(IN) :: DPA
! total solids volume fraction of solids phases 
      DOUBLE PRECISION, INTENT(IN) :: PHIS      
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Reynolds number 
      DOUBLE PRECISION :: RE
! Stokes Drag Force
      DOUBLE PRECISION :: F_STOKES
! dimensionless drag force F
      DOUBLE PRECISION :: F      
!-----------------------------------------------

      IF(Mug > ZERO) THEN
! Note the presence of gas volume fraction in ROPG
         RE = DPA*VREL*ROPg/Mug        ! use DPA
      ELSE
         RE = LARGE_NUMBER 
      ENDIF      

! eq(9) BVK J. fluid. Mech. 528, 2005 
! (this F_Stokes is /= of Koch_Hill by a factor of ep_g)
      F_STOKES = 18D0*Mug*EPg/DPM**2   ! use DPM

      F = 10d0*phis/EPg**2 + EPg**2*(ONE+1.5d0*DSQRT(phis))
      F = F + 0.413d0*RE/(24.d0*EPg**2) * &
             (ONE/EPg + 3d0*EPg*phis + 8.4d0/RE**0.343) / &
             (ONE+10.d0**(3d0*phis)/RE**(0.5+2.d0*phis))

      lDgA = F*F_STOKES
      IF (RE == ZERO) lDgA = ZERO      
   
      RETURN
      END SUBROUTINE DRAG_BVK
       


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: DRAG_HYS                                                C
!  Purpose: Calculate the gas-solids drag coefficient                  C
!                                                                      C
!  Literature/Document References:                                     C
!     Yin, X, Sundaresan, S. (2009). AIChE Journal 55: no 6, 1352-     C
!     1368                                                             C
!                                                                      C
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC

      SUBROUTINE DRAG_HYS(lDgA,EPg,Mug,ROPg,VREL,&
                 DPM,DPA,Y_I,EP_sM,PHIS,M,MAXM,IJK)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param 
      USE param1
      USE drag
      USE run
      USE physprop
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! drag coefficient
      DOUBLE PRECISION, INTENT(OUT) :: lDgA
! gas volume fraction 
      DOUBLE PRECISION, INTENT(IN) :: EPg
! gas laminar viscosity 
      DOUBLE PRECISION, INTENT(IN) :: Mug
! gas density*EP_g
      DOUBLE PRECISION, INTENT(IN) :: ROPg
! magnitude of gas-solids relative velocity 
      DOUBLE PRECISION, INTENT(IN) :: VREL
! local variable for the particle diameter 
      DOUBLE PRECISION :: DPM(2*DIM_M)      
! average particle diameter 
      DOUBLE PRECISION, INTENT(IN) :: DPA
! diameter ratio in polydisperse systems
      DOUBLE PRECISION, INTENT(IN) :: Y_i
! local variable for the solids volume fraction 
      DOUBLE PRECISION :: EP_sM(2*DIM_M)
! total solids volume fraction of solids phases 
      DOUBLE PRECISION, INTENT(IN) :: PHIS 
! current solids phase index and fluid cell index
      INTEGER, INTENT(IN) :: M, IJK
! maximum number of solids phases
      INTEGER, INTENT(IN) :: MAXM      
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! minimum particle diameter in mixture
      DOUBLE PRECISION :: DPmin
! Index for particles of other species
      INTEGER :: L
! Reynolds number 
      DOUBLE PRECISION :: RE
! Stokes Drag Force
      DOUBLE PRECISION :: F_STOKES
! dimensionless drag force F
      DOUBLE PRECISION :: F
! Polydisperse correction factor for YS drag relation
      DOUBLE PRECISION :: a_YS
! Lubrication interaction prefactor in YS drag relation	
      DOUBLE PRECISION :: alpha_YS
! Friction coefficient for a particle of type i (HYS drag relation)
      DOUBLE PRECISION :: beta_i_HYS
! Friction coefficient for a particle of type j (HYS drag relation)
      DOUBLE PRECISION :: beta_j_HYS
! Stokes drag of a particle of type j
      DOUBLE PRECISION :: FSTOKES_j
! Diameter ratio for particle of type j 
      DOUBLE PRECISION :: Y_i_J
! Variable for Beetstra et. al. drag relation
      DOUBLE PRECISION :: F_D_BVK
! Variable for YS drag relation
      DOUBLE PRECISION :: F_YS
!-----------------------------------------------

      IF (Mug > ZERO) THEN
! Note the presence of gas volume fraction in ROPG              
         RE = DPA*VREL*ROPg/Mug   ! use DPA
      ELSE
         RE = LARGE_NUMBER
      ENDIF

! (this F_Stokes is /= of Koch_Hill by a factor of ep_g/ep_sm)
! (this F_Stokes is /= of BVK by a factor of 1/ep_sm)      
      F_STOKES = 18D0*Mug*EPg*EP_SM(M)/DPM(M)**2   ! use DPM

! Find smallest diameter if number of particle types is greater than 1
      Dpmin= DPM(1)               
      IF (MAXM > 1) THEN
         DO L=2,MAXM
            Dpmin = MIN(Dpmin,DPM(L))
         ENDDO
      ENDIF

      a_YS = 1d0 - 2.66d0*phis + 9.096d0*phis**2 - 11.338d0*phis**3 

! Calculate the prefactor of the off-diagonal friction coefficient
! Use default value of lamdba if there are no particle asparities
      alpha_YS = 1.313d0*LOG10(DPmin/lam_HYS) - 1.249d0


! Beetstra correction for monodisperse drag
      F_D_BVK = ZERO
      F = 10.d0*phis/EPg**2 + EPg**2*(ONE+1.5d0*DSQRT(phis))

      IF(RE > ZERO) F_D_BVK = F + 0.413d0*RE/(24.d0*EPg**2)*&
         (ONE/EPg + 3.d0*EPg*phis + 8.4d0/RE**0.343d0) / &
         (ONE+10**(3.d0*phis)/RE**(0.5d0+2.d0*phis))

! YS correction for polydisperse drag
      F_YS = 1d0/EPg + (F_D_BVK - 1d0/EPg)*&
                       (a_YS*Y_i+(1d0-a_YS)*Y_i**2)
      F_YS = F_YS*F_STOKES
      beta_i_HYS = F_YS
              
      DO L= 1,MAXM
         IF (L /= M) THEN
            Y_i_J = DPM(L)/DPA
            beta_j_HYS = 1.d0/EPg + (F_D_BVK - 1.d0/EPg) * &
               (a_YS*Y_i_J + (1d0-a_YS)*Y_i_J**2)
            FSTOKES_j = 18.D0*Mug*EP_sM(L)*EPg/&
               DPM(L)**2 

            beta_j_HYS = beta_j_HYS*FSTOKES_j
                      
! Calculate off-diagonal friction coefficient
            beta_ij(IJK,M,L) = ZERO
                      
! This if statement prevents NaN values from appearing for beta_ij
            IF (EP_sM(M) > ZERO .AND. EP_SM(L) > ZERO) &
               beta_ij(IJK,M,L) = (2.d0*alpha_YS*EP_sM(M)*EP_sM(L))/ &
                  (EP_sM(M)/beta_i_HYS + EP_sM(L)/beta_j_HYS)
            F_YS = F_YS + beta_ij(IJK,M,L)

         ENDIF   ! end if (J/=M)
      ENDDO   ! end do (J=1,MAXM)


      lDgA = F_YS

      RETURN
      END SUBROUTINE DRAG_HYS


 
     
