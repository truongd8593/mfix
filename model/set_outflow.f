!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SET_OUTFLOW                                             C
!  Purpose: Set specified pressure outflow bc for a specified range    C
!           of cells. This routine is also called for mass_outlow or   C
!           outflow bcs                                                C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-JAN-92  C
!  Reviewer:M. Syamlal, S. Venkatesan, P. Nicoletti,  Date: 29-JAN-92  C
!           W. Rogers                                                  C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables modified: RO_g, ROP_g, EP_g                               C
!                      ROP_s, U_g, U_s, V_g, V_s,                      C
!                      W_g, W_s, T_g, T_s                              C
!                                                                      C
!  Local variables: IJK, I, J, K, M, LFLUID                            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE SET_OUTFLOW(BCV, I1, I2, J1, J2, K1, K2) 

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param 
      USE param1 
      USE bc
      USE fldvar
      USE geometry
      USE indices
      USE physprop
      USE constant
      USE scalars
      USE run
      USE compar
      USE mflux
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Boundary condition number 
      INTEGER, INTENT(IN) :: BCV 
! Starting and ending I index 
      INTEGER, INTENT(IN) :: I1, I2
! Starting and ending J index 
      INTEGER, INTENT(IN) :: J1, J2 
! Starting and ending K index 
      INTEGER, INTENT(IN) :: K1, K2
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! indices 
      INTEGER :: I, J, K, M, N 
! Local index for boundary cell 
      INTEGER :: IJK 
! Local index for a fluid cell near the boundary cell 
      INTEGER :: LFLUID 
!-----------------------------------------------
! External functions
!-----------------------------------------------
      DOUBLE PRECISION, EXTERNAL :: EOSG 
!-----------------------------------------------
! Include statement functions
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'ep_s2.inc'
!-----------------------------------------------

      DO K = K1, K2 
         DO J = J1, J2 
            DO I = I1, I2 
! Check if current i,j,k resides on this PE
               IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
               IJK = FUNIJK(I,J,K) 


! Fluid cell at West
! ---------------------------------------------------------------->>>
               IF (FLUID_AT(IM_OF(IJK))) THEN 
                  LFLUID = IM_OF(IJK) 
                  IF (U_G(LFLUID)>=ZERO .OR. EP_G(IJK)==UNDEFINED) THEN 
                     IF (BC_TYPE(BCV) /= 'P_OUTFLOW') P_G(IJK) = P_G(LFLUID) 
                     T_G(IJK) = T_G(LFLUID) 
                     N = 1 
                     IF (NMAX(0) > 0) THEN 
                        X_G(IJK,:NMAX(0)) = X_G(LFLUID,:NMAX(0)) 
                        N = NMAX(0) + 1 
                     ENDIF 
                     MW_MIX_G(IJK) = MW_MIX_G(LFLUID) 
                     IF (RO_G0 == UNDEFINED) RO_G(IJK) = EOSG(MW_MIX_G(IJK),P_G&
                        (IJK),T_G(IJK)) 
                  ENDIF 
                  P_STAR(IJK) = P_STAR(LFLUID) 
! if bc_ep_g undefined, set ep_g to 1 
                  IF (BC_EP_G(BCV) == UNDEFINED) EP_G(IJK) = ONE 
  
                  DO N = 1, NScalar
                     M = Phase4Scalar(N)
                     IF(M == 0)Then
                        IF (U_G(LFLUID)>=ZERO) THEN 
                           Scalar(IJK, N) = Scalar(LFLUID, N)
                        ENDIF
                     Else
                        IF (U_s(LFLUID, M)>=ZERO) THEN 
                           Scalar(IJK, N) = Scalar(LFLUID, N)
                        ENDIF
                     Endif
                  ENDDO

                  IF(K_Epsilon) THEN
                     IF (U_G(LFLUID) >= ZERO) THEN 
                        K_Turb_G(IJK) = K_Turb_G(LFLUID)
                        E_Turb_G(IJK) = E_Turb_G(LFLUID)
                     ENDIF
                  ENDIF

                  IF (TRIM(KT_TYPE) == 'GHD') THEN 
                     P_S(IJK,MMAX) = P_S(LFLUID,MMAX) 
                     ROP_S(IJK,MMAX) = ZERO
                  ENDIF

                  DO M = 1, SMAX 
                     P_S(IJK,M) = P_S(LFLUID,M)
! check if solids are entering boundary cell from adjacent fluid cell,
! if so set boundary cell values of rop_s, t_s and theta_m according to
! adjacent fluid cell values
                     IF (U_S(LFLUID,M) >= ZERO) THEN 
                        ROP_S(IJK,M) = ROP_S(LFLUID,M) 
                        T_S(IJK,M) = T_S(LFLUID,M)
                        THETA_M(IJK,M) =  THETA_M(LFLUID,M)
                     ELSE 
                        ROP_S(IJK,M) = ZERO 
                        THETA_M(IJK,M) =  THETA_M(LFLUID,M)  ! cannot have theta_m = zero anywhere
                     ENDIF 

! if bc_rop_s is defined at the boundary, set rop_s accordingly
                     IF(BC_ROP_S(BCV,M)/=UNDEFINED)ROP_S(IJK,M)=BC_ROP_S(BCV,M) 

                     IF (TRIM(KT_TYPE) == 'GHD') THEN
                        ROP_S(IJK,MMAX) = ROP_S(IJK,MMAX) + ROP_S(IJK,M)
                        IF(M==SMAX) THETA_M(IJK,MMAX) =  THETA_M(LFLUID,MMAX)
                     ENDIF

! if bc_ep_g is undefined, set ep_g at boundary to be consistent with 1=ep_s+ep_g
                     IF(BC_EP_G(BCV)==UNDEFINED)EP_G(IJK)=EP_G(IJK)-EP_S(IJK,M) 

                     N = 1 
                     IF (NMAX(M) > 0) THEN 
                        X_S(IJK,M,:NMAX(M)) = X_S(LFLUID,M,:NMAX(M)) 
                        N = NMAX(M) + 1 
                     ENDIF 
                  ENDDO 
                  ROP_G(IJK) = RO_G(IJK)*EP_G(IJK) 
                  IF (ROP_G(IJK) > ZERO) THEN 
                     U_G(IJK) = ROP_G(LFLUID)*U_G(LFLUID)/ROP_G(IJK) 
                  ELSE 
                     U_G(IJK) = ZERO 
                  ENDIF 
! set boundary cell values of velocity according to adjacent fluid cell
! values
                  V_G(IJK) = V_G(LFLUID) 
                  W_G(IJK) = W_G(LFLUID) 
                  Flux_gE(IJK) = Flux_gE(LFLUID)
                  Flux_gN(IJK) = Flux_gN(LFLUID)
                  Flux_gT(IJK) = Flux_gT(LFLUID) 
                  IF(ADDED_MASS) THEN
                    Flux_gSE(IJK) = Flux_gSE(LFLUID)
                    Flux_gSN(IJK) = Flux_gSN(LFLUID)
                    Flux_gST(IJK) = Flux_gST(LFLUID) 
                  ENDIF

! for GHD theory
                  IF (TRIM(KT_TYPE) == 'GHD') THEN 
                     Flux_nE(IJK) = Flux_nE(LFLUID)
                     Flux_nN(IJK) = Flux_nN(LFLUID)
                     Flux_nT(IJK) = Flux_nT(LFLUID) 
                  ENDIF
! end of GHD theory

                  IF (MMAX > 0) THEN 
                     WHERE (ROP_S(IJK,:MMAX) > ZERO)  
                        U_S(IJK,:MMAX) = ROP_S(LFLUID,:MMAX)*U_S(LFLUID,:MMAX)/&
                           ROP_S(IJK,:MMAX) 
                     ELSEWHERE 
                        U_S(IJK,:MMAX) = ZERO 
                     END WHERE 
                     V_S(IJK,:MMAX) = V_S(LFLUID,:MMAX) 
                     W_S(IJK,:MMAX) = W_S(LFLUID,:MMAX) 
                     Flux_sE(IJK,:MMAX) = Flux_sE(LFLUID,:MMAX)
                     Flux_sN(IJK,:MMAX) = Flux_sN(LFLUID,:MMAX)
                     Flux_sT(IJK,:MMAX) = Flux_sT(LFLUID,:MMAX)  
                     IF(ADDED_MASS) THEN 
                       Flux_sSE(IJK) = Flux_sSE(LFLUID)
                       Flux_sSN(IJK) = Flux_sSN(LFLUID)
                       Flux_sST(IJK) = Flux_sST(LFLUID)  
                     ENDIF
                  ENDIF 
               ENDIF          ! if (FLUID_AT(IM_OF(IJK)))
! ----------------------------------------------------------------<<<


! Fluid cell at East
! ---------------------------------------------------------------->>>
               IF (FLUID_AT(IP_OF(IJK))) THEN 
                  LFLUID = IP_OF(IJK) 
                  IF (U_G(IJK)<=ZERO .OR. EP_G(IJK)==UNDEFINED) THEN 
                     IF (BC_TYPE(BCV) /= 'P_OUTFLOW') P_G(IJK) = P_G(LFLUID) 
                     T_G(IJK) = T_G(LFLUID) 
                     N = 1 
                     IF (NMAX(0) > 0) THEN 
                        X_G(IJK,:NMAX(0)) = X_G(LFLUID,:NMAX(0)) 
                        N = NMAX(0) + 1 
                     ENDIF 
                     MW_MIX_G(IJK) = MW_MIX_G(LFLUID) 
                     IF (RO_G0 == UNDEFINED) RO_G(IJK) = EOSG(MW_MIX_G(IJK),P_G&
                        (IJK),T_G(IJK)) 
                  ENDIF 
                  P_STAR(IJK) = P_STAR(LFLUID) 
                  IF (BC_EP_G(BCV) == UNDEFINED) EP_G(IJK) = ONE 

                  DO N = 1, NScalar
                     M = Phase4Scalar(N)
                     IF(M == 0)Then
                        IF (U_G(LFLUID) <= ZERO) THEN 
                           Scalar(IJK, N) = Scalar(LFLUID, N)
                        ENDIF
                     Else
                        IF (U_s(LFLUID, M) <= ZERO) THEN 
                           Scalar(IJK, N) = Scalar(LFLUID, N)
                        ENDIF
                     Endif
                  ENDDO

                  IF(K_Epsilon) THEN
                     IF (U_G(LFLUID) <= ZERO) THEN 
                        K_Turb_G(IJK) = K_Turb_G(LFLUID)
                        E_Turb_G(IJK) = E_Turb_G(LFLUID)
                     ENDIF
                  ENDIF

                  IF (TRIM(KT_TYPE) == 'GHD') THEN 
                     P_S(IJK,MMAX) = P_S(LFLUID,MMAX) 
                     ROP_S(IJK,MMAX) = ZERO
                  ENDIF

                  DO M = 1, SMAX 
                     P_S(IJK,M) = P_S(LFLUID,M) 
                     IF (U_S(IJK,M) <= ZERO) THEN 
                        ROP_S(IJK,M) = ROP_S(LFLUID,M) 
                        T_S(IJK,M) = T_S(LFLUID,M) 
                        THETA_M(IJK,M) =  THETA_M(LFLUID,M)
                     ELSE 
                        ROP_S(IJK,M) = ZERO 
                        THETA_M(IJK,M) =  THETA_M(LFLUID,M)
                     ENDIF 

                     IF(BC_ROP_S(BCV,M)/=UNDEFINED)ROP_S(IJK,M)=BC_ROP_S(BCV,M)   

                     IF (TRIM(KT_TYPE) == 'GHD') THEN
                        ROP_S(IJK,MMAX) = ROP_S(IJK,MMAX) + ROP_S(IJK,M)
                        IF(M==SMAX) THETA_M(IJK,MMAX) =  THETA_M(LFLUID,MMAX)
                     ENDIF

                     IF(BC_EP_G(BCV)==UNDEFINED)EP_G(IJK)=EP_G(IJK)-EP_S(IJK,M) 

                     N = 1 
                     IF (NMAX(M) > 0) THEN 
                        X_S(IJK,M,:NMAX(M)) = X_S(LFLUID,M,:NMAX(M)) 
                        N = NMAX(M) + 1 
                     ENDIF 
                  ENDDO 
                  ROP_G(IJK) = RO_G(IJK)*EP_G(IJK) 
                  IF (U_G(IJK) == UNDEFINED) THEN 
                     IF (ROP_G(IJK) > ZERO) THEN 
                        U_G(IJK) = ROP_G(LFLUID)*U_G(LFLUID)/ROP_G(IJK) 
                     ELSE 
                        U_G(IJK) = ZERO 
                     ENDIF 
                  ENDIF 
                  V_G(IJK) = V_G(LFLUID) 
                  W_G(IJK) = W_G(LFLUID) 
                  Flux_gE(IJK) = Flux_gE(LFLUID)
                  Flux_gN(IJK) = Flux_gN(LFLUID)
                  Flux_gT(IJK) = Flux_gT(LFLUID)  
                  IF(ADDED_MASS) THEN
                    Flux_gSE(IJK) = Flux_gSE(LFLUID)
                    Flux_gSN(IJK) = Flux_gSN(LFLUID)
                    Flux_gST(IJK) = Flux_gST(LFLUID) 
                  ENDIF

! for GHD theory
                  IF (TRIM(KT_TYPE) == 'GHD') THEN 
                     Flux_nE(IJK) = Flux_nE(LFLUID)
                     Flux_nN(IJK) = Flux_nN(LFLUID)
                     Flux_nT(IJK) = Flux_nT(LFLUID) 
                  ENDIF
! end of GHD theory

                  DO M = 1, MMAX 
                     IF (U_S(IJK,M) == UNDEFINED) THEN 
                        IF (ROP_S(IJK,M) > ZERO) THEN 
                           U_S(IJK,M) = ROP_S(LFLUID,M)*U_S(LFLUID,M)/ROP_S(IJK&
                              ,M) 
                        ELSE 
                           U_S(IJK,M) = ZERO 
                        ENDIF 
                     ENDIF 
                     V_S(IJK,M) = V_S(LFLUID,M) 
                     W_S(IJK,M) = W_S(LFLUID,M) 
                     Flux_sE(IJK,M) = Flux_sE(LFLUID,M)
                     Flux_sN(IJK,M) = Flux_sN(LFLUID,M)
                     Flux_sT(IJK,M) = Flux_sT(LFLUID,M) 
                     IF(ADDED_MASS) THEN 
                       Flux_sSE(IJK) = Flux_sSE(LFLUID)
                       Flux_sSN(IJK) = Flux_sSN(LFLUID)
                       Flux_sST(IJK) = Flux_sST(LFLUID)  
                     ENDIF
                  END DO 
               ENDIF          ! if (FLUID_AT(IP_OF(IJK)))
! ----------------------------------------------------------------<<<


! Fluid cell at South
! ---------------------------------------------------------------->>>
               IF (FLUID_AT(JM_OF(IJK))) THEN 
                  LFLUID = JM_OF(IJK) 
                  IF (V_G(LFLUID)>=ZERO .OR. EP_G(IJK)==UNDEFINED) THEN 
                     IF (BC_TYPE(BCV) /= 'P_OUTFLOW') P_G(IJK) = P_G(LFLUID) 
                     T_G(IJK) = T_G(LFLUID) 
                     N = 1 
                     IF (NMAX(0) > 0) THEN 
                        X_G(IJK,:NMAX(0)) = X_G(LFLUID,:NMAX(0)) 
                        N = NMAX(0) + 1 
                     ENDIF 
                     MW_MIX_G(IJK) = MW_MIX_G(LFLUID) 
                     IF (RO_G0 == UNDEFINED) RO_G(IJK) = EOSG(MW_MIX_G(IJK),P_G&
                        (IJK),T_G(IJK)) 
                  ENDIF 
                  P_STAR(IJK) = P_STAR(LFLUID) 
                  IF (BC_EP_G(BCV) == UNDEFINED) EP_G(IJK) = ONE 

                  DO N = 1, NScalar
                     M = Phase4Scalar(N)
                     IF(M == 0)Then
                        IF (V_G(LFLUID) >= ZERO) THEN 
                           Scalar(IJK, N) = Scalar(LFLUID, N)
                        ENDIF
                     Else
                        IF (V_s(LFLUID, M) >= ZERO) THEN 
                           Scalar(IJK, N) = Scalar(LFLUID, N)
                        ENDIF
                     Endif
                  ENDDO

                  IF(K_Epsilon) THEN
                    IF (V_G(LFLUID) >= ZERO) THEN 
                       K_Turb_G(IJK) = K_Turb_G(LFLUID)
                       E_Turb_G(IJK) = E_Turb_G(LFLUID)
                    ENDIF
                  ENDIF

                  IF (TRIM(KT_TYPE) == 'GHD') THEN 
                     P_S(IJK,MMAX) = P_S(LFLUID,MMAX) 
                     ROP_S(IJK,MMAX) = ZERO
                  ENDIF

                  DO M = 1, SMAX  
                     P_S(IJK,M) = P_S(LFLUID,M) 
                     IF (V_S(LFLUID,M) >= 0.) THEN 
                        ROP_S(IJK,M) = ROP_S(LFLUID,M) 
                        T_S(IJK,M) = T_S(LFLUID,M) 
                        THETA_M(IJK,M) =  THETA_M(LFLUID,M)
                     ELSE 
                        ROP_S(IJK,M) = ZERO 
                        THETA_M(IJK,M) =  THETA_M(LFLUID,M)
                     ENDIF 

! overwrite calculation of rop_s with defined value (bc_rop_s)
                     IF(BC_ROP_S(BCV,M)/=UNDEFINED)ROP_S(IJK,M)=BC_ROP_S(BCV,M)   
                     IF (TRIM(KT_TYPE) == 'GHD') THEN
                        ROP_S(IJK,MMAX) = ROP_S(IJK,MMAX) + ROP_S(IJK,M)
                        IF(M==SMAX) THETA_M(IJK,MMAX) =  THETA_M(LFLUID,MMAX)
                     ENDIF

! make ep_g consistent with calculated rop_s
                     IF(BC_EP_G(BCV)==UNDEFINED)EP_G(IJK)=EP_G(IJK)-EP_S(IJK,M) 

                     N = 1 
                     IF (NMAX(M) > 0) THEN 
                        X_S(IJK,M,:NMAX(M)) = X_S(LFLUID,M,:NMAX(M)) 
                        N = NMAX(M) + 1 
                     ENDIF 
                  ENDDO 

                  ROP_G(IJK) = RO_G(IJK)*EP_G(IJK) 
                  U_G(IJK) = U_G(LFLUID) 
                  IF (ROP_G(IJK) > ZERO) THEN 
                     V_G(IJK) = ROP_G(LFLUID)*V_G(LFLUID)/ROP_G(IJK) 
                  ELSE 
                     V_G(IJK) = ZERO 
                  ENDIF 
                  W_G(IJK) = W_G(LFLUID) 
                  Flux_gE(IJK) = Flux_gE(LFLUID)
                  Flux_gN(IJK) = Flux_gN(LFLUID)
                  Flux_gT(IJK) = Flux_gT(LFLUID)  
                  IF(ADDED_MASS) THEN
                    Flux_gSE(IJK) = Flux_gSE(LFLUID)
                    Flux_gSN(IJK) = Flux_gSN(LFLUID)
                    Flux_gST(IJK) = Flux_gST(LFLUID) 
                  ENDIF
! for GHD theory
                  IF (TRIM(KT_TYPE) == 'GHD') THEN 
                     Flux_nE(IJK) = Flux_nE(LFLUID)
                     Flux_nN(IJK) = Flux_nN(LFLUID)
                     Flux_nT(IJK) = Flux_nT(LFLUID) 
                  ENDIF
! end of GHD theory

                  IF (MMAX > 0) THEN 
                     U_S(IJK,:MMAX) = U_S(LFLUID,:MMAX) 
                     V_S(IJK,:MMAX) = V_S(LFLUID,:MMAX) 
                     W_S(IJK,:MMAX) = W_S(LFLUID,:MMAX) 
                     Flux_sE(IJK,:MMAX) = Flux_sE(LFLUID,:MMAX)
                     Flux_sN(IJK,:MMAX) = Flux_sN(LFLUID,:MMAX)
                     Flux_sT(IJK,:MMAX) = Flux_sT(LFLUID,:MMAX) 
                     IF(ADDED_MASS) THEN 
                       Flux_sSE(IJK) = Flux_sSE(LFLUID)
                       Flux_sSN(IJK) = Flux_sSN(LFLUID)
                       Flux_sST(IJK) = Flux_sST(LFLUID)  
                     ENDIF
                  ENDIF 
               ENDIF          ! if FLUID_AT(JM_OF(IJK)))
! ----------------------------------------------------------------<<<


! Fluid cell at North
! ---------------------------------------------------------------->>>
               IF (FLUID_AT(JP_OF(IJK))) THEN 
                  LFLUID = JP_OF(IJK) 
                  IF (V_G(IJK)<=ZERO .OR. EP_G(IJK)==UNDEFINED) THEN 
                     IF (BC_TYPE(BCV) /= 'P_OUTFLOW') P_G(IJK) = P_G(LFLUID) 
                     T_G(IJK) = T_G(LFLUID) 
                     N = 1 
                     IF (NMAX(0) > 0) THEN 
                        X_G(IJK,:NMAX(0)) = X_G(LFLUID,:NMAX(0)) 
                        N = NMAX(0) + 1 
                     ENDIF 
                     MW_MIX_G(IJK) = MW_MIX_G(LFLUID) 
                     IF (RO_G0 == UNDEFINED) RO_G(IJK) = EOSG(MW_MIX_G(IJK),P_G&
                        (IJK),T_G(IJK)) 
                  ENDIF 
                  P_STAR(IJK) = P_STAR(LFLUID) 
                  IF (BC_EP_G(BCV) == UNDEFINED) EP_G(IJK) = ONE 

                  DO N = 1, NScalar
                     M = Phase4Scalar(N)
                     IF(M == 0)Then
                        IF (V_G(LFLUID) <= ZERO) THEN 
                           Scalar(IJK, N) = Scalar(LFLUID, N)
                        ENDIF
                     Else
                        IF (V_s(LFLUID, M) <= ZERO) THEN 
                           Scalar(IJK, N) = Scalar(LFLUID, N)
                        ENDIF
                     Endif
                  ENDDO

                  IF(K_Epsilon) THEN
                    IF (V_G(LFLUID) <= ZERO) THEN 
                       K_Turb_G(IJK) = K_Turb_G(LFLUID)
                       E_Turb_G(IJK) = E_Turb_G(LFLUID)
                    ENDIF
                  ENDIF

                  IF (TRIM(KT_TYPE) == 'GHD') THEN 
                     P_S(IJK,MMAX) = P_S(LFLUID,MMAX) 
                     ROP_S(IJK,MMAX) = ZERO
                  ENDIF

                  DO M = 1, SMAX 
                     P_S(IJK,M) = P_S(LFLUID,M) 
                     IF (V_S(IJK,M) <= ZERO) THEN 
                        ROP_S(IJK,M) = ROP_S(LFLUID,M) 
                        T_S(IJK,M) = T_S(LFLUID,M) 
                        THETA_M(IJK,M) =  THETA_M(LFLUID,M)
                     ELSE 
                        ROP_S(IJK,M) = ZERO 
                        THETA_M(IJK,M) =  THETA_M(LFLUID,M)
                     ENDIF 

                     IF(BC_ROP_S(BCV,M)/=UNDEFINED)ROP_S(IJK,M)=BC_ROP_S(BCV,M)   
                     IF (TRIM(KT_TYPE) == 'GHD') THEN
                        ROP_S(IJK,MMAX) = ROP_S(IJK,MMAX) + ROP_S(IJK,M)
                        IF(M==SMAX) THETA_M(IJK,MMAX) =  THETA_M(LFLUID,MMAX)
                     ENDIF

                     IF(BC_EP_G(BCV)==UNDEFINED)EP_G(IJK)=EP_G(IJK)-EP_S(IJK,M) 

                     N = 1 
                     IF (NMAX(M) > 0) THEN 
                        X_S(IJK,M,:NMAX(M)) = X_S(LFLUID,M,:NMAX(M)) 
                        N = NMAX(M) + 1 
                     ENDIF 
                  ENDDO 
                  ROP_G(IJK) = RO_G(IJK)*EP_G(IJK) 
                  U_G(IJK) = U_G(LFLUID) 
                  IF (V_G(IJK) == UNDEFINED) THEN 
                     IF (ROP_G(IJK) > ZERO) THEN 
                        V_G(IJK) = ROP_G(LFLUID)*V_G(LFLUID)/ROP_G(IJK) 
                     ELSE 
                        V_G(IJK) = ZERO 
                     ENDIF 
                  ENDIF 
                  W_G(IJK) = W_G(LFLUID) 
                  Flux_gE(IJK) = Flux_gE(LFLUID)
                  Flux_gN(IJK) = Flux_gN(LFLUID)
                  Flux_gT(IJK) = Flux_gT(LFLUID) 
                  IF(ADDED_MASS) THEN
                    Flux_gSE(IJK) = Flux_gSE(LFLUID)
                    Flux_gSN(IJK) = Flux_gSN(LFLUID)
                    Flux_gST(IJK) = Flux_gST(LFLUID) 
                  ENDIF

! for GHD theory
                  IF (TRIM(KT_TYPE) == 'GHD') THEN 
                     Flux_nE(IJK) = Flux_nE(LFLUID)
                     Flux_nN(IJK) = Flux_nN(LFLUID)
                     Flux_nT(IJK) = Flux_nT(LFLUID) 
                  ENDIF
! end of GHD theory

                  IF (MMAX > 0) THEN 
                     U_S(IJK,:MMAX) = U_S(LFLUID,:MMAX) 
                     WHERE (V_S(IJK,:MMAX) == UNDEFINED) V_S(IJK,:MMAX) = V_S(&
                        LFLUID,:MMAX) 
                     W_S(IJK,:MMAX) = W_S(LFLUID,:MMAX) 
                     Flux_sE(IJK,:MMAX) = Flux_sE(LFLUID,:MMAX)
                     Flux_sN(IJK,:MMAX) = Flux_sN(LFLUID,:MMAX)
                     Flux_sT(IJK,:MMAX) = Flux_sT(LFLUID,:MMAX) 
                     IF(ADDED_MASS) THEN 
                       Flux_sSE(IJK) = Flux_sSE(LFLUID)
                       Flux_sSN(IJK) = Flux_sSN(LFLUID)
                       Flux_sST(IJK) = Flux_sST(LFLUID)  
                     ENDIF
                  ENDIF 
               ENDIF          ! if (FLUID_AT(JP_OF(IJK)))
! ----------------------------------------------------------------<<<


! Fluid cell at Bottom
! ---------------------------------------------------------------->>>
               IF (FLUID_AT(KM_OF(IJK))) THEN 
                  LFLUID = KM_OF(IJK) 
                  IF (W_G(LFLUID)>=ZERO .OR. EP_G(IJK)==UNDEFINED) THEN 
                     IF (BC_TYPE(BCV) /= 'P_OUTFLOW') P_G(IJK) = P_G(LFLUID) 
                     T_G(IJK) = T_G(LFLUID) 
                     N = 1 
                     IF (NMAX(0) > 0) THEN 
                        X_G(IJK,:NMAX(0)) = X_G(LFLUID,:NMAX(0)) 
                        N = NMAX(0) + 1 
                     ENDIF 
                     MW_MIX_G(IJK) = MW_MIX_G(LFLUID) 
                     IF (RO_G0 == UNDEFINED) RO_G(IJK) = EOSG(MW_MIX_G(IJK),P_G&
                        (IJK),T_G(IJK)) 
                  ENDIF 
                  P_STAR(IJK) = P_STAR(LFLUID) 
                  IF (BC_EP_G(BCV) == UNDEFINED) EP_G(IJK) = ONE 

                  DO N = 1, NScalar
                     M = Phase4Scalar(N)
                     IF(M == 0)Then
                        IF (W_G(LFLUID) >= ZERO) THEN 
                           Scalar(IJK, N) = Scalar(LFLUID, N)
                        ENDIF
                     Else
                        IF (W_s(LFLUID, M) >= ZERO) THEN 
                           Scalar(IJK, N) = Scalar(LFLUID, N)
                        ENDIF
                     Endif
                  END DO

                  IF(K_Epsilon) THEN
                     IF (W_G(LFLUID) >= ZERO) THEN 
                        K_Turb_G(IJK) = K_Turb_G(LFLUID)
                        E_Turb_G(IJK) = E_Turb_G(LFLUID)
                     ENDIF
                  ENDIF

                  IF (TRIM(KT_TYPE) == 'GHD') THEN 
                     P_S(IJK,MMAX) = P_S(LFLUID,MMAX) 
                     ROP_S(IJK,MMAX) = ZERO
                  ENDIF
                  DO M = 1, SMAX 
                     P_S(IJK,M) = P_S(LFLUID,M) 
                     IF (W_S(LFLUID,M) >= 0.) THEN 
                        ROP_S(IJK,M) = ROP_S(LFLUID,M) 
                        T_S(IJK,M) = T_S(LFLUID,M) 
                        THETA_M(IJK,M) =  THETA_M(LFLUID,M)
                     ELSE 
                        ROP_S(IJK,M) = ZERO 
                        THETA_M(IJK,M) =  THETA_M(LFLUID,M)
                     ENDIF 

                     IF(BC_ROP_S(BCV,M)/=UNDEFINED)ROP_S(IJK,M)=BC_ROP_S(BCV,M)  

                     IF (TRIM(KT_TYPE) == 'GHD') THEN
                        ROP_S(IJK,MMAX) = ROP_S(IJK,MMAX) + ROP_S(IJK,M)
                        IF(M==SMAX) THETA_M(IJK,MMAX) =  THETA_M(LFLUID,MMAX)
                     ENDIF

                     IF(BC_EP_G(BCV)==UNDEFINED)EP_G(IJK)=EP_G(IJK)-EP_S(IJK,M) 
                     N = 1 
                     IF (NMAX(M) > 0) THEN 
                        X_S(IJK,M,:NMAX(M)) = X_S(LFLUID,M,:NMAX(M)) 
                        N = NMAX(M) + 1 
                     ENDIF 
                  END DO 
                  ROP_G(IJK) = RO_G(IJK)*EP_G(IJK) 
                  U_G(IJK) = U_G(LFLUID) 
                  V_G(IJK) = V_G(LFLUID) 
                  IF (ROP_G(IJK) > ZERO) THEN 
                     W_G(IJK) = ROP_G(LFLUID)*W_G(LFLUID)/ROP_G(IJK) 
                  ELSE 
                     W_G(IJK) = ZERO 
                  ENDIF 
                  Flux_gE(IJK) = Flux_gE(LFLUID)
                  Flux_gN(IJK) = Flux_gN(LFLUID)
                  Flux_gT(IJK) = Flux_gT(LFLUID) 
                  IF(ADDED_MASS) THEN
                    Flux_gSE(IJK) = Flux_gSE(LFLUID)
                    Flux_gSN(IJK) = Flux_gSN(LFLUID)
                    Flux_gST(IJK) = Flux_gST(LFLUID) 
                  ENDIF

! for GHD theory
                  IF (TRIM(KT_TYPE) == 'GHD') THEN 
                     Flux_nE(IJK) = Flux_nE(LFLUID)
                     Flux_nN(IJK) = Flux_nN(LFLUID)
                     Flux_nT(IJK) = Flux_nT(LFLUID) 
                  ENDIF
! end of GHD theory

                  IF (MMAX > 0) THEN 
                     U_S(IJK,:MMAX) = U_S(LFLUID,:MMAX) 
                     V_S(IJK,:MMAX) = V_S(LFLUID,:MMAX) 
                     W_S(IJK,:MMAX) = W_S(LFLUID,:MMAX) 
                     Flux_sE(IJK,:MMAX) = Flux_sE(LFLUID,:MMAX)
                     Flux_sN(IJK,:MMAX) = Flux_sN(LFLUID,:MMAX)
                     Flux_sT(IJK,:MMAX) = Flux_sT(LFLUID,:MMAX) 
                     IF(ADDED_MASS) THEN 
                       Flux_sSE(IJK) = Flux_sSE(LFLUID)
                       Flux_sSN(IJK) = Flux_sSN(LFLUID)
                       Flux_sST(IJK) = Flux_sST(LFLUID)  
                     ENDIF
                  ENDIF 
               ENDIF          ! if (FLUID_AT(KM_OF(IJK)))
! ----------------------------------------------------------------<<<


! Fluid cell at Top
! ---------------------------------------------------------------->>>
               IF (FLUID_AT(KP_OF(IJK))) THEN 
                  LFLUID = KP_OF(IJK) 
                  IF (W_G(IJK)<=ZERO .OR. EP_G(IJK)==UNDEFINED) THEN 
                     IF (BC_TYPE(BCV) /= 'P_OUTFLOW') P_G(IJK) = P_G(LFLUID) 
                     T_G(IJK) = T_G(LFLUID) 
                     N = 1 
                     IF (NMAX(0) > 0) THEN 
                        X_G(IJK,:NMAX(0)) = X_G(LFLUID,:NMAX(0)) 
                        N = NMAX(0) + 1 
                     ENDIF 
                     MW_MIX_G(IJK) = MW_MIX_G(LFLUID) 
                     IF (RO_G0 == UNDEFINED) RO_G(IJK) = EOSG(MW_MIX_G(IJK),P_G&
                        (IJK),T_G(IJK)) 
                  ENDIF 
                  P_STAR(IJK) = P_STAR(LFLUID) 
                  IF (BC_EP_G(BCV) == UNDEFINED) EP_G(IJK) = ONE 

                  DO N = 1, NScalar
                     M = Phase4Scalar(N)
                     IF(M == 0)Then
                       IF (W_G(LFLUID) <= ZERO) THEN 
                          Scalar(IJK, N) = Scalar(LFLUID, N)
                        ENDIF
                     Else
                        IF (W_s(LFLUID, M) <= ZERO) THEN 
                           Scalar(IJK, N) = Scalar(LFLUID, N)
                        ENDIF
                     Endif
                  END DO

                  IF(K_Epsilon) THEN
                     IF (W_G(LFLUID) <= ZERO) THEN 
                        K_Turb_G(IJK) = K_Turb_G(LFLUID)
                        E_Turb_G(IJK) = E_Turb_G(LFLUID)
                     ENDIF
                  ENDIF

                  IF (TRIM(KT_TYPE) == 'GHD') THEN 
                     P_S(IJK,MMAX) = P_S(LFLUID,MMAX) 
                     ROP_S(IJK,MMAX) = ZERO
                  ENDIF
                  DO M = 1, SMAX 
                     P_S(IJK,M) = P_S(LFLUID,M) 
                     IF (W_S(IJK,M) <= ZERO) THEN 
                        ROP_S(IJK,M) = ROP_S(LFLUID,M) 
                        T_S(IJK,M) = T_S(LFLUID,M) 
                        THETA_M(IJK,M) =  THETA_M(LFLUID,M)
                     ELSE 
                        ROP_S(IJK,M) = ZERO 
                        THETA_M(IJK,M) =  THETA_M(LFLUID,M)
                     ENDIF 

                     IF(BC_ROP_S(BCV,M)/=UNDEFINED)ROP_S(IJK,M)=BC_ROP_S(BCV,M)   

                     IF (TRIM(KT_TYPE) == 'GHD') THEN
                        ROP_S(IJK,MMAX) = ROP_S(IJK,MMAX) + ROP_S(IJK,M)
                        IF(M==SMAX) THETA_M(IJK,MMAX) =  THETA_M(LFLUID,MMAX)
                     ENDIF

                     IF(BC_EP_G(BCV)==UNDEFINED)EP_G(IJK)=EP_G(IJK)-EP_S(IJK,M) 

                     N = 1 
                     IF (NMAX(M) > 0) THEN 
                        X_S(IJK,M,:NMAX(M)) = X_S(LFLUID,M,:NMAX(M)) 
                        N = NMAX(M) + 1 
                     ENDIF 
                  END DO 
                  ROP_G(IJK) = RO_G(IJK)*EP_G(IJK) 
                  U_G(IJK) = U_G(LFLUID) 
                  V_G(IJK) = V_G(LFLUID) 
                  IF (W_G(IJK) == UNDEFINED) THEN 
                     IF (ROP_G(IJK) > ZERO) THEN 
                        W_G(IJK) = ROP_G(LFLUID)*W_G(LFLUID)/ROP_G(IJK) 
                     ELSE 
                        W_G(IJK) = ZERO 
                     ENDIF 
                  ENDIF 
                  Flux_gE(IJK) = Flux_gE(LFLUID)
                  Flux_gN(IJK) = Flux_gN(LFLUID)
                  Flux_gT(IJK) = Flux_gT(LFLUID) 
                  IF(ADDED_MASS) THEN
                    Flux_gSE(IJK) = Flux_gSE(LFLUID)
                    Flux_gSN(IJK) = Flux_gSN(LFLUID)
                    Flux_gST(IJK) = Flux_gST(LFLUID) 
                  ENDIF

! for GHD theory
                  IF (TRIM(KT_TYPE) == 'GHD') THEN 
                     Flux_nE(IJK) = Flux_nE(LFLUID)
                     Flux_nN(IJK) = Flux_nN(LFLUID)
                     Flux_nT(IJK) = Flux_nT(LFLUID) 
                  ENDIF
! end of GHD theory

                  IF (MMAX > 0) THEN 
                     U_S(IJK,:MMAX) = U_S(LFLUID,:MMAX) 
                     V_S(IJK,:MMAX) = V_S(LFLUID,:MMAX) 
                     WHERE (W_S(IJK,:MMAX) == UNDEFINED) W_S(IJK,:MMAX) = W_S(&
                        LFLUID,:MMAX) 
                     Flux_sE(IJK,:MMAX) = Flux_sE(LFLUID,:MMAX)
                     Flux_sN(IJK,:MMAX) = Flux_sN(LFLUID,:MMAX)
                     Flux_sT(IJK,:MMAX) = Flux_sT(LFLUID,:MMAX) 
                     IF(ADDED_MASS) THEN 
                       Flux_sSE(IJK) = Flux_sSE(LFLUID)
                       Flux_sSN(IJK) = Flux_sSN(LFLUID)
                       Flux_sST(IJK) = Flux_sST(LFLUID)  
                     ENDIF
                  ENDIF 
               ENDIF          ! if (FLUID_AT(KP_OF(IJK)))
! ----------------------------------------------------------------<<<

            END DO 
         END DO 
      END DO 

      RETURN  
      END SUBROUTINE SET_OUTFLOW 
