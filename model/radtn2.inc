!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!  Source terms for the energy equations.  The source term is linearized as
!  S = S_c - S_p * T,  where S_c and S_p must be positive.
!
!  By default the source terms have been coded for radiation sources.

!
!  S_p for gas phase at i, j, k
      S_Rpg(IJK) = 4.d0 * GAMA_Rg(IJK) *  T_g(IJK)**3
!
!  S_c for gas phase at i, j, k
      S_Rcg(IJK) = GAMA_Rg(IJK) * ( T_Rg(IJK)**4 + 3.d0 * T_g(IJK)**4 )

!
!  S_p for solids phase at i, j, k
      S_Rps(IJK, M) = 4.d0 * GAMA_Rs(IJK, M) *  T_s(IJK, M)**3
!
!  S_c for solids phase at i, j, k
      S_Rcs(IJK, M) = GAMA_Rs(IJK, M) * ( T_Rs(IJK, M)**4 &
                                       + 3.d0 * T_s(IJK, M)**4 )
