!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SET_IC                                                 C
!  Purpose: This module sets all the initial conditions.               C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-JAN-92  C
!  Reviewer:M. Syamlal, S. Venkatesan, P. Nicoletti,  Date: 29-JAN-92  C
!           W. Rogers                                                  C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: Modifications for initializing solids pressure             C
!  Author: M. Syamlal                                 Date: 11-MAY-93  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: IC_DEFINED, IC_P_g, FLAG, IC_I_e, IC_I_w,     C
!                        IC_EP_g, IC_T_g, IC_T_s,  IC_U_g,    C
!                        IC_V_g, IC_W_g, IC_ROP_s, MMAX, IC_U_s,       C
!                        IC_V_s, IC_W_s, IC_K_b, IC_K_t, IC_J_s,       C
!                        IC_J_n                                        C
!                                                                      C
!  Variables modified: M, I, J, K, IJK, EP_g, P_g, T_g, T_s,     C
!                      U_g, V_g, W_g, ROP_s, U_s, V_s, W_s, P_star     C
!                                                                      C
!  Local variables: EPGX, PGX, TGX, TS1X, TS2X, UGX, VGX, WGX, ROPSX,  C
!                   USX, VSX, WSX, L                                  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE SET_IC 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE geometry
      USE constant
      USE physprop
      USE ic
      USE fldvar
      USE visc_g
      USE indices
      USE scales 
      USE energy
      USE scalars
      USE compar 
      USE run
      USE sendrecv 
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
! 
!                      indices 
      INTEGER          I, J, K, IJK, M, N 
! 
!                      Local index for initial condition 
      INTEGER          L 
! 
!                      Temporary variable for storing IC_EP_g 
      DOUBLE PRECISION EPGX 
! 
!                      Temporary variable for storing IC_P_g 
      DOUBLE PRECISION PGX 
! 
!                      Temporary variable for storing P_star 
      DOUBLE PRECISION PSX 
! 
!                      Temporary variable for storing IC_T_g 
      DOUBLE PRECISION TGX 
! 
!                      Temporary variable for storing IC_U_g 
      DOUBLE PRECISION UGX 
! 
!                      Temporary variable for storing IC_V_g 
      DOUBLE PRECISION VGX 
! 
!                      Temporary variable for storing IC_W_g 
      DOUBLE PRECISION WGX 
! 
!                      number densities for use in GHD theory only
      DOUBLE PRECISION ni(DIMENSION_M) 
! 
!                      Temporary variable for storing IC_ROP_s 
      DOUBLE PRECISION ROPSX (DIMENSION_M) 
! 
!                      Temporary variable for storing IC_T_s 
      DOUBLE PRECISION TSX (DIMENSION_M) 
! 
!                      Temporary variable for storing IC_U_s 
      DOUBLE PRECISION USX (DIMENSION_M) 
! 
!                      Temporary variable for storing IC_V_s 
      DOUBLE PRECISION VSX (DIMENSION_M) 
! 
!                      Temporary variable for storing IC_W_s 
      DOUBLE PRECISION WSX (DIMENSION_M) 
! 
!                      number density for GHD theory 
      DOUBLE PRECISION nM, nTOT 
! 
!-----------------------------------------------
      INCLUDE 'sc_p_g1.inc'
      INCLUDE 's_pr1.inc'
      INCLUDE 'function.inc'
      INCLUDE 's_pr2.inc'
      INCLUDE 'sc_p_g2.inc'

!  Set the initial conditions.

      DO L = 1, DIMENSION_IC 
         IF (IC_DEFINED(L)) THEN 
            EPGX = IC_EP_G(L) 
            PGX = IC_P_G(L) 
            PSX = IC_P_STAR(L) 
            IF (PSX==UNDEFINED .AND. IC_TYPE(L)/='PATCH') PSX = ZERO 
            TGX = IC_T_G(L) 
            UGX = IC_U_G(L) 
            VGX = IC_V_G(L) 
            WGX = IC_W_G(L) 
            IF (IC_L_SCALE(L) /= UNDEFINED) RECALC_VISC_G = .TRUE. 
            IF(K_Epsilon) RECALC_VISC_G = .TRUE.

            M = 1 
            IF (MMAX > 0) THEN 
              ROPSX(:MMAX) = IC_ROP_S(L,:MMAX) 
              TSX(:MMAX) = IC_T_S(L,:MMAX) 
              USX(:MMAX) = IC_U_S(L,:MMAX) 
              VSX(:MMAX) = IC_V_S(L,:MMAX) 
              WSX(:MMAX) = IC_W_S(L,:MMAX) 
              M = MMAX + 1 
            ENDIF 

            DO K = IC_K_B(L), IC_K_T(L) 
            DO J = IC_J_S(L), IC_J_N(L) 
            DO I = IC_I_W(L), IC_I_E(L) 
              IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
              IJK = FUNIJK(I,J,K) 

              IF (FLAG(IJK) == 1) THEN 
                  IF (EPGX /= UNDEFINED) EP_G(IJK) = EPGX 
                  IF (IC_TYPE(L) == 'PATCH') THEN 
                      IF (PGX /= UNDEFINED) P_G(IJK) = SCALE(PGX) 
                  ELSE 
                      IF (PGX /= UNDEFINED) THEN 
! If P_g is UNDEFINED then it is
                          P_G(IJK) = SCALE(PGX) 
                      ELSE                  !  reset in SET_FLUIDBED_P 
                          P_G(IJK) = UNDEFINED 
                      ENDIF 
                  ENDIF 
                  IF (PSX /= UNDEFINED) P_STAR(IJK) = PSX 
                  IF (TGX /= UNDEFINED) T_G(IJK) = TGX 
                  IF (IC_L_SCALE(L) /= UNDEFINED) L_SCALE(IJK) = &
                      IC_L_SCALE(L)

                  N = 1 
                  IF (NMAX(0) > 0) THEN 
                      WHERE (IC_X_G(L,:NMAX(0)) /= UNDEFINED) X_G(IJK,:&
                          NMAX(0)) = IC_X_G(L,:NMAX(0)) 
                      N = NMAX(0) + 1 
                  ENDIF 

                  IF (NScalar > 0) THEN 
                      WHERE (IC_Scalar(L,:NScalar) /= UNDEFINED)&
                          Scalar(IJK,:NScalar) = IC_Scalar(L,:NScalar) 
                  ENDIF  

                  IF (K_Epsilon) THEN 
                      IF (IC_K_Turb_G(L) /= UNDEFINED)&
                          K_Turb_G(IJK) = IC_K_Turb_G(L) 
                      IF (IC_E_Turb_G(L) /= UNDEFINED)&
                          E_Turb_G(IJK) = IC_E_Turb_G(L)
                  ENDIF

                  IF (UGX /= UNDEFINED) U_G(IJK) = UGX 
                  IF (VGX /= UNDEFINED) V_G(IJK) = VGX 
                  IF (WGX /= UNDEFINED) W_G(IJK) = WGX 

                  GAMA_RG(IJK) = IC_GAMA_RG(L) 
                  IF (IC_T_RG(L) /= UNDEFINED) THEN 
                      T_RG(IJK) = IC_T_RG(L) 
                  ELSE 
                      T_RG(IJK) = ZERO 
                  ENDIF 

                  DO M = 1, MMAX 
                    IF (ROPSX(M) /= UNDEFINED) ROP_S(IJK,M) = ROPSX(M)
                    IF (TSX(M) /= UNDEFINED) T_S(IJK,M) = TSX(M) 
                    IF (IC_THETA_M(L,M) /= UNDEFINED) THETA_M(IJK,M) = &
                        IC_THETA_M(L,M) 
                    IF (USX(M) /= UNDEFINED) U_S(IJK,M) = USX(M) 
                    IF (VSX(M) /= UNDEFINED) V_S(IJK,M) = VSX(M) 
                    IF (WSX(M) /= UNDEFINED) W_S(IJK,M) = WSX(M) 

                    GAMA_RS(IJK,M) = IC_GAMA_RS(L,M) 
                    IF (IC_T_RS(L,M) /= UNDEFINED) THEN 
                        T_RS(IJK,M) = IC_T_RS(L,M) 
                    ELSE 
                      T_RS(IJK,M) = ZERO 
                    ENDIF 

                    N = 1 
                    IF (NMAX(M) > 0) THEN 
                        WHERE (IC_X_S(L,M,:NMAX(M)) /= UNDEFINED) X_S(IJK&
                            ,M,:NMAX(M)) = IC_X_S(L,M,:NMAX(M)) 
                        N = NMAX(M) + 1 
                    ENDIF 
                  ENDDO 

! for GHD theory to compute mixture IC of velocity and density
                  IF(TRIM(KT_TYPE) == 'GHD') THEN
                    ROP_S(IJK,MMAX) = ZERO
                    U_S(IJK,MMAX) = ZERO
                    V_S(IJK,MMAX) = ZERO
                    W_S(IJK,MMAX) = ZERO
                    THETA_M(IJK,MMAX) = ZERO
                    nTOT = ZERO
                    nM = ZERO
                    DO M = 1, SMAX 
                      IF (ROPSX(M) /= UNDEFINED) THEN
                        ROP_S(IJK,MMAX) = ROP_S(IJK,MMAX) + ROPSX(M) 
                        nM = ROPSX(M)*6d0/(PI*D_p(IJK,M)**3*RO_S(M))  
                        nTOT = nTOT + nM 
                      ENDIF
                      IF (IC_THETA_M(L,M) /= UNDEFINED) THETA_M(IJK,MMAX) = &
                        THETA_M(IJK,MMAX) + nM*IC_THETA_M(L,M) 
                      IF (USX(M) /= UNDEFINED .AND. ROPSX(M) /= UNDEFINED) &
                        U_S(IJK,MMAX) = U_S(IJK,MMAX) + ROPSX(M)*USX(M) 
                      IF (VSX(M) /= UNDEFINED .AND. ROPSX(M) /= UNDEFINED) &
                        V_S(IJK,MMAX) = V_S(IJK,MMAX) +  ROPSX(M)*VSX(M) 
                      IF (WSX(M) /= UNDEFINED .AND. ROPSX(M) /= UNDEFINED) &
                        W_S(IJK,MMAX) = W_S(IJK,MMAX) +  ROPSX(M)*WSX(M)
                    ENDDO 

                    IF(ROP_S(IJK,MMAX) > ZERO) THEN  ! if ropsTotal > 0 then RoN_T > 0
                      U_S(IJK,MMAX) = U_S(IJK,MMAX) / ROP_S(IJK,MMAX)
                      V_S(IJK,MMAX) = V_S(IJK,MMAX) / ROP_S(IJK,MMAX)
                      W_S(IJK,MMAX) = W_S(IJK,MMAX) / ROP_S(IJK,MMAX) 
                      THETA_M(IJK,MMAX) = THETA_M(IJK,MMAX) / nTOT 
                    ELSE ! for initially empty bed
                      U_S(IJK,MMAX) = U_S(IJK,MMAX)
                      V_S(IJK,MMAX) = V_S(IJK,MMAX)
                      W_S(IJK,MMAX) = W_S(IJK,MMAX)
                      DO M = 1, SMAX 
		        THETA_M(IJK,MMAX) = THETA_M(IJK,M) ! set T > 0 in case Ti > 0
                      ENDDO
                      IF(THETA_M(IJK,MMAX)==ZERO) THETA_M(IJK,MMAX) = small_number
                    ENDIF
                  ENDIF
! end of modifications for GHD theory

              ENDIF     ! if (flag = 1)
            ENDDO       ! over i
            ENDDO       ! over j
            ENDDO       ! over k
         ENDIF          ! if (ic_defined)
      ENDDO             ! over dimension_ic
      
      CALL SEND_RECV(L_SCALE,2)
      RETURN  
      END SUBROUTINE SET_IC 

!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization
!// 350 Changed do loop limits: 1,kmax2->kmin3,kmax3      
!// 400 Added sendrecv module and send_recv calls for COMMunication
