C
C                      Gas-phase heat of reaction
      DOUBLE PRECISION HOR_g (DIMENSION_3)
C
C                      Solids-phase heat of reaction
      DOUBLE PRECISION HOR_s (DIMENSION_3, DIMENSION_M)
C
C                      Gas-solids heat transfer coefficient
      DOUBLE PRECISION GAMA_gs (DIMENSION_3, DIMENSION_M)
C
C                      Gas-phase radiation coefficient
      DOUBLE PRECISION GAMA_Rg (DIMENSION_3)
C
C                      Solids-phase radiation coefficient
      DOUBLE PRECISION GAMA_Rs (DIMENSION_3, DIMENSION_M)
C
C                      Gas-phase radiation temperature
      DOUBLE PRECISION T_Rg (DIMENSION_3)
C
C                      Solids-phase radiation temperature
      DOUBLE PRECISION T_Rs (DIMENSION_3, DIMENSION_M)
C
      COMMON /energy_DP/
     &  HOR_g, HOR_s,
     &  GAMA_gs, GAMA_Rg, GAMA_Rs, T_Rg, T_Rs
