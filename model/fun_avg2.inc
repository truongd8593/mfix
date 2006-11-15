! Arithmetic averages
!     i+1/2, j, k. Use IJKE for scalar
!          F(IJK) F(IPJK)  I
      AVG_X(XXXm, XXXp,    xL) = FX_E(xL) * XXXm + FX_E_bar(xL) * XXXp
!
!     i+1, j, k
!          F(IJK) F(IPJK)  IP
      AVG_X_E(XXXm, XXXp,    xL) = FX(xL) * XXXm + FX_bar(xL) * XXXp
!
!     i, j+1/2, k. Use IJKN for scalar
!           F(IJK) F(IJPK) J  
      AVG_Y(XXXm,  XXXp,   xL) = FY_N(xL) * XXXm + FY_N_bar(xL) * XXXp
!
!     i, j+1, k
!           F(IJK) F(IJPK)  
      AVG_Y_N(XXXm,  XXXp) = HALF *( XXXm + XXXp )
!
!     i, j, k+1/2. Use IJKT for scalar
!           F(IJK) F(IJKP) K  
      AVG_Z(XXXm,  XXXp,   xL) = FZ_T(xL) * XXXm + FZ_T_bar(xL) * XXXp
!
!     i, j, k+1
!           F(IJK) F(IJKP)
      AVG_Z_T(XXXm,  XXXp) = HALF *( XXXm + XXXp )
!
! Harmonic averages
!
!     i+1/2, j, k. Use IJKE for scalar
!             F(IJK) F(IPJK)  I      
      AVG_X_h(XXXm,  XXXp,    xL) = XXXm * XXXp / &
       MAX(SMALL_NUMBER, ( FX_E(xL) * XXXm + FX_E_bar(xL) * XXXp) )
!
!
!     i, j+1/2, k. Use IJKN for scalar
!             F(IJK) F(IJPK) J  
      AVG_Y_h(XXXm,  XXXp,   xL) = XXXm * XXXp / &
       MAX(SMALL_NUMBER, ( FY_N(xL) * XXXm + FY_N_bar(xL) * XXXp) )
!
!
!     i, j, k+1/2. Use IJKT for scalar
!             F(IJK) F(IJKP) K  
      AVG_Z_h(XXXm,  XXXp,   xL) = XXXm * XXXp / &
       MAX(SMALL_NUMBER, ( FZ_T(xL) * XXXm + FZ_T_bar(xL) * XXXp) )
!
! Harmonic averages for possibly negative scalars (sof, Aug 22 2006)
!
!     i+1/2, j, k. Use IJKE for scalar
!             F(IJK) F(IPJK)  I      
      AVG_X_S(XXXm,  XXXp,    xL) = XXXm * XXXp / &
       (1D-30 + ( FX_E(xL) * XXXm + FX_E_bar(xL) * XXXp) )
!
!
!     i, j+1/2, k. Use IJKN for scalar
!             F(IJK) F(IJPK) J  
      AVG_Y_S(XXXm,  XXXp,   xL) = XXXm * XXXp / &
       (1D-30 + ( FY_N(xL) * XXXm + FY_N_bar(xL) * XXXp) )
!
!
!     i, j, k+1/2. Use IJKT for scalar
!             F(IJK) F(IJKP) K  
      AVG_Z_S(XXXm,  XXXp,   xL) = XXXm * XXXp / &
       (1D-30 + ( FZ_T(xL) * XXXm + FZ_T_bar(xL) * XXXp) )

