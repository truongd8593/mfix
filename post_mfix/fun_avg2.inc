CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: Averaging functions                                    C
C  Purpose:                                                            C
C                                                                      C
C                                                                      C
C  Author: M. Syamlal                                 Date: 1-MAY-96   C
C  Reviewer:                                          Date:            C
C                                                                      C
C                                                                      C
C  Literature/Document References:                                     C
C                                                                      C
C  Variables referenced:                                               C
C  Variables modified:                                                 C
C                                                                      C
C  Local variables:                                                    C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
C Arithmetic averages
C     i+1/2, j, k. Use IJKE for scalar
C          F(IJK) F(IPJK)  I
      AVG_X(XXXm, XXXp,    xL) = FX_E(xL) * XXXm + FX_E_bar(xL) * XXXp
C
C     i+1, j, k
C          F(IJK) F(IPJK)  IP
      AVG_X_E(XXXm, XXXp,    xL) = FX(xL) * XXXm + FX_bar(xL) * XXXp
C
C     i, j+1/2, k. Use IJKN for scalar
C           F(IJK) F(IJPK) J  
      AVG_Y(XXXm,  XXXp,   xL) = FY_N(xL) * XXXm + FY_N_bar(xL) * XXXp
C
C     i, j+1, k
C           F(IJK) F(IJPK)  
      AVG_Y_N(XXXm,  XXXp) = HALF *( XXXm + XXXp )
C
C     i, j, k+1/2. Use IJKT for scalar
C           F(IJK) F(IJKP) K  
      AVG_Z(XXXm,  XXXp,   xL) = FZ_T(xL) * XXXm + FZ_T_bar(xL) * XXXp
C
C     i, j, k+1
C           F(IJK) F(IJKP)
      AVG_Z_T(XXXm,  XXXp) = HALF *( XXXm + XXXp )
C
C Harmonic averages
C
C     i+1/2, j, k. Use IJKE for scalar
C             F(IJK) F(IPJK)  I      
      AVG_X_h(XXXm,  XXXp,    xL) = XXXm * XXXp /
     &  MAX(SMALL_NUMBER, ( FX_E(xL) * XXXm + FX_E_bar(xL) * XXXp) )
C
C
C     i, j+1/2, k. Use IJKN for scalar
C             F(IJK) F(IJPK) J  
      AVG_Y_h(XXXm,  XXXp,   xL) = XXXm * XXXp /
     &  MAX(SMALL_NUMBER, ( FY_N(xL) * XXXm + FY_N_bar(xL) * XXXp) )
C
C
C     i, j, k+1/2. Use IJKT for scalar
C             F(IJK) F(IJKP) K  
      AVG_Z_h(XXXm,  XXXp,   xL) = XXXm * XXXp /
     &  MAX(SMALL_NUMBER, ( FZ_T(xL) * XXXm + FZ_T_bar(xL) * XXXp) )
