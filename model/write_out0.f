!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_OUT0                                             C
!  Purpose: echo user input                                            C
!                                                                      C
!  Author: P. Nicoletti, M. Syamlal                   Date: 04-DEC-91  C
!  Reviewer: W. Rogers, M. Syamlal, S. Venkatesan     Date: 31-JAN-92  C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: add node, version                                          C
!  Author: P.Nicoletti                                Date: 07-FEB-92  C
!  Reviewer: S. Venkatesan                            Date: 11-DEC-92  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: ID_MONTH, ID_DAY, ID_YEAR, RUN_NAME, UNITS    C
!                        DESCRIPTION, RUN_TYPE, DX, IMAX1, IMAX, DY    C
!                        JMAX1, YLENGTH, XLENGTH, JMAX, DZ, KMAX1      C
!                        KMAX, ZLENGTH, MMAX, D_p, RO_s, EP_star, MU_g0C
!                        MW_AVG, IC_DEFINED, IC_X_w, IC_X_e, IC_Y_s    C
!                        IC_Z_b, IC_Z_t, IC_I_w, IC_I_e, IC_J_s        C
!                        IC_J_n, IC_K_b, IC_K_t, IC_EP_g, IC_P_g       C
!                        IC_U_g, IC_V_g, IC_W_g, IC_ROP_s, IC_T_s      C
!                        IC_U_s, IC_V_s, IC_W_s, BC_DEFINED            C
!                        BC_TYPE, BC_X_w, BC_X_e, BC_Y_s, BC_Y_n       C
!                        BC_Z_b, BC_Z_t, BC_I_w, BC_I_e, BC_J_s        C
!                        BC_J_n, BC_Z_b, BC_Z_t, BC_EP_g, BC_P_g       C
!                        BC_T_g, BC_U_g, BC_V_g, BC_W_g, BC_ROP_s      C
!                        BC_T_s, BC_U_s, BC_V_s, BC_W_s                C
!                        ICBC_FLAG, IMIN1, JMIN1, KMIN1, ID_NODE       C
!                        ID_VERSION, RO_g0                             C
!  Variables modified: M                                               C
!                                                                      C
!  Local variables: L, LOC                                             C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE WRITE_OUT0 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE run
      USE output
      USE physprop
      USE geometry
      USE ic
      USE bc
      USE is
      USE fldvar
      USE constant
      USE indices
      USE funits 
      USE toleranc 
      USE scales 
      USE scalars
      USE ur_facs 
      USE leqsol 
      USE compar         !//d
      USE mpi_utility    !//d
      USE sendrecv    !//d
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
      INTEGER :: L, M, N 
      DOUBLE PRECISION, DIMENSION(6) :: LOC 
! 
!                      Coefficient of restitution (old symbol) 
      DOUBLE PRECISION :: E 
      CHARACTER, DIMENSION(3) :: LEGEND*3 
      CHARACTER, DIMENSION(0:8) :: DISCR_NAME*12 
      CHARACTER, DIMENSION(0:8) :: DISCR_NAME1*12 
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      DOUBLE PRECISION , EXTERNAL :: LOCATION 
!-----------------------------------------------

!
!
!
!                      Coefficient of restitution (old symbol)
      DATA DISCR_NAME/'FOUP', 'FOUP', 'Superbee', 'Smart', 'Ultra-Quick', &
         'QUICKEST', 'Muscl', 'VanLeer', 'Minmod'/ 
      DATA DISCR_NAME1/'FOUP', 'FOUP', 'Fourth Order', 'Smart', 'Ultra-Quick', &
         'QUICKEST', 'Muscl', 'VanLeer', 'Minmod'/ 

      if (myPE.ne.PE_IO) return

!
!  Write Headers for .OUT file
!
      WRITE(UNIT_OUT,1000)ID_VERSION,ID_HOUR,ID_MINUTE,ID_MONTH,ID_DAY,ID_YEAR 
      WRITE (UNIT_OUT, 1010) ID_NODE(1:50) 
!
!  Echo input data
!
!  Run control section
!
      WRITE (UNIT_OUT, 1100) 
      WRITE (UNIT_OUT, 1110) RUN_NAME 
      WRITE (UNIT_OUT, 1120) DESCRIPTION 
      WRITE (UNIT_OUT, 1130) UNITS 
      IF (DT /= UNDEFINED) THEN 
         WRITE (UNIT_OUT, 1135) TIME, TSTOP, DT, DT_MAX, DT_MIN, DT_FAC 
      ELSE 
         WRITE (UNIT_OUT, 1136) 
      ENDIF 
      WRITE (UNIT_OUT, 1137) RUN_TYPE 
      IF (RUN_TYPE == 'NEW') THEN 
         WRITE (UNIT_OUT, 1138) 
      ELSE IF (RUN_TYPE == 'RESTART_1') THEN 
         WRITE (UNIT_OUT, 1139) 
      ENDIF 
      IF (MOMENTUM_X_EQ(0)) THEN 
         WRITE (UNIT_OUT, 1140) 'X', ' ' 
      ELSE 
         WRITE (UNIT_OUT, 1140) 'X', ' NOT ' 
      ENDIF 
      IF (MOMENTUM_Y_EQ(0)) THEN 
         WRITE (UNIT_OUT, 1140) 'Y', ' ' 
      ELSE 
         WRITE (UNIT_OUT, 1140) 'Y', ' NOT ' 
      ENDIF 
      IF (MOMENTUM_Z_EQ(0)) THEN 
         WRITE (UNIT_OUT, 1140) 'Z', ' ' 
      ELSE 
         WRITE (UNIT_OUT, 1140) 'Z', ' NOT ' 
      ENDIF 
      DO M = 1, MMAX 
         IF (MOMENTUM_X_EQ(M)) THEN 
            WRITE (UNIT_OUT, 1141) M, 'X', ' ' 
         ELSE 
            WRITE (UNIT_OUT, 1141) M, 'X', ' NOT ' 
         ENDIF 
         IF (MOMENTUM_Y_EQ(M)) THEN 
            WRITE (UNIT_OUT, 1141) M, 'Y', ' ' 
         ELSE 
            WRITE (UNIT_OUT, 1141) M, 'Y', ' NOT ' 
         ENDIF 
         IF (MOMENTUM_Z_EQ(M)) THEN 
            WRITE (UNIT_OUT, 1141) M, 'Z', ' ' 
         ELSE 
            WRITE (UNIT_OUT, 1141) M, 'Z', ' NOT ' 
         ENDIF 
      END DO 
      IF (ENERGY_EQ) THEN 
         WRITE (UNIT_OUT, 1143) 
      ELSE 
         WRITE (UNIT_OUT, 1144) 
      ENDIF 
      IF (SPECIES_EQ(0)) THEN 
         WRITE (UNIT_OUT, 1145) 
      ELSE 
         WRITE (UNIT_OUT, 1146) 
      ENDIF 
      DO M = 1, MMAX 
         IF (SPECIES_EQ(M)) THEN 
            WRITE (UNIT_OUT, 1147) M 
         ELSE 
            WRITE (UNIT_OUT, 1148) M 
         ENDIF 
      END DO 
      IF (CALL_USR) THEN 
         WRITE (UNIT_OUT, 1149) ' ' 
      ELSE 
         WRITE (UNIT_OUT, 1149) ' NOT ' 
      ENDIF 
      IF (MODEL_B) WRITE (UNIT_OUT, 1101) 
      IF (Nscalar /= 0)THEN
        WRITE (UNIT_OUT, 1102)NScalar
        DO L = 1, NScalar
          WRITE (UNIT_OUT, 1103)L, Phase4Scalar(L)
        END DO 
      ENDIF
      IF (K_Epsilon) WRITE (UNIT_OUT, 1104)    
!
!  Physical and numerical parameters
!
      WRITE (UNIT_OUT, 1150) 
      IF (C_E /= UNDEFINED) WRITE (UNIT_OUT, 1151) C_E 
      IF (C_F /= UNDEFINED) WRITE (UNIT_OUT, 1152) C_F 
      IF (PHI /= UNDEFINED) WRITE (UNIT_OUT, 1153) PHI 
      IF (PHI_W /= UNDEFINED) WRITE (UNIT_OUT, 1154) PHI_W 
      WRITE (UNIT_OUT, 1155) L_SCALE0, MU_GMAX 
      IF (V_EX /= ZERO) WRITE (UNIT_OUT, 1156) V_EX 
      WRITE (UNIT_OUT, 1157) P_REF, P_SCALE, GRAVITY 
      WRITE (UNIT_OUT, 1158) 
      IF(FPFOI) THEN
         WRITE (UNIT_OUT, 1159) (UR_FAC(L),LEQ_IT(L),LEQ_METHOD(L),&
	                     LEQ_SWEEP(L), LEQ_TOL(L),&
			     DISCR_NAME1(DISCRETIZE(L)),L=1,9) 
      ELSE
         WRITE (UNIT_OUT, 1159) (UR_FAC(L),LEQ_IT(L),LEQ_METHOD(L),&
	                     LEQ_SWEEP(L), LEQ_TOL(L),&
			     DISCR_NAME(DISCRETIZE(L)),L=1,9) 
      ENDIF

      DO L = 1, DIMENSION_C 
         IF (C(L) /= UNDEFINED) WRITE (UNIT_OUT, 1190) C_NAME(L), L, C(L) 
      END DO 
      WRITE (UNIT_OUT, 1200) 
      WRITE (UNIT_OUT, 1201) COORDINATES 
      IF (CYCLIC_X_PD) THEN 
         WRITE (UNIT_OUT, 1202) 'X', ' with pressure drop' 
         WRITE (UNIT_OUT, 1203) 'X', DELP_X 
      ELSE IF (CYCLIC_X) THEN 
         WRITE (UNIT_OUT, 1202) 'X' 
      ENDIF 
      IF (CYCLIC_Y_PD) THEN 
         WRITE (UNIT_OUT, 1202) 'Y', ' with pressure drop' 
         WRITE (UNIT_OUT, 1203) 'Y', DELP_Y 
      ELSE IF (CYCLIC_Y) THEN 
         WRITE (UNIT_OUT, 1202) 'Y' 
      ENDIF 
      IF (CYCLIC_Z_PD) THEN 
         WRITE (UNIT_OUT, 1202) 'Z', ' with pressure drop' 
         WRITE (UNIT_OUT, 1203) 'Z', DELP_Z 
      ELSE IF (CYCLIC_Z) THEN 
         WRITE (UNIT_OUT, 1202) 'Z' 
      ENDIF 
      WRITE (UNIT_OUT, 1210) 
      LEGEND(1) = '  I' 
      LEGEND(2) = ' DX' 
      LEGEND(3) = 'X_E' 
      CALL WRITE_TABLE (LEGEND, DX, XMIN, 1, IMAX2) 
      IF (XMIN /= ZERO) WRITE (UNIT_OUT, 1211) XMIN 
      WRITE (UNIT_OUT, 1212) IMAX 
      WRITE (UNIT_OUT, 1213) XLENGTH 
      WRITE (UNIT_OUT, 1220) 
      LEGEND(1) = '  J' 
      LEGEND(2) = ' DY' 
      LEGEND(3) = 'Y_N' 
      CALL WRITE_TABLE (LEGEND, DY, ZERO, 1, JMAX2) 
      WRITE (UNIT_OUT, 1221) JMAX 
      WRITE (UNIT_OUT, 1222) YLENGTH 
      WRITE (UNIT_OUT, 1230) 
      LEGEND(1) = '  K' 
      LEGEND(2) = ' DZ' 
      LEGEND(3) = 'Z_T' 
      CALL WRITE_TABLE (LEGEND, DZ, ZERO, 1, KMAX2) 
      WRITE (UNIT_OUT, 1231) KMAX 
      WRITE (UNIT_OUT, 1232) ZLENGTH 
!
!  Gas Section
!
      WRITE (UNIT_OUT, 1300) 
      IF (RO_G0 /= UNDEFINED) WRITE (UNIT_OUT, 1305) RO_G0 
      IF (MU_G0 /= UNDEFINED) WRITE (UNIT_OUT, 1310) MU_G0 
      IF (SPECIES_EQ(0)) THEN 
         WRITE (UNIT_OUT, 1315) NMAX(0) 
         WRITE (UNIT_OUT, 1316) 
         DO N = 1, NMAX(0) 
            WRITE (UNIT_OUT, 1317) N, MW_G(N) 
         END DO 
      ENDIF 
      IF (MW_AVG /= UNDEFINED) WRITE (UNIT_OUT, 1320) MW_AVG 
!
!  Particle Section
!
      WRITE (UNIT_OUT, 1400) 
      IF (MU_S0 /= UNDEFINED) WRITE (UNIT_OUT, 1405) MU_S0 
      WRITE (UNIT_OUT, 1410) MMAX 
      WRITE (UNIT_OUT, 1420) 
      DO M = 1, MMAX 
         WRITE (UNIT_OUT, 1421) M, D_P(M), RO_S(M), CLOSE_PACKED(M) 
      END DO 
      DO M = 1, MMAX 
         IF (SPECIES_EQ(M)) THEN 
            WRITE (UNIT_OUT, 1422) M, M, NMAX(M) 
            WRITE (UNIT_OUT, 1423) 
            DO N = 1, NMAX(M) 
               WRITE (UNIT_OUT, 1424) N, MW_S(M,N) 
            END DO 
         ENDIF 
      END DO 
      WRITE (UNIT_OUT, 1430) EP_STAR 
!
!  Initial Conditions Section
!
      WRITE (UNIT_OUT, 1500) 
      DO L = 1, DIMENSION_IC 
         IF (IC_DEFINED(L)) THEN 
            WRITE (UNIT_OUT, 1510) L 
            LOC(1) = LOCATION(IC_I_W(L),XMIN,DX) - HALF*DX(IC_I_W(L)) 
            LOC(2) = LOCATION(IC_I_E(L),XMIN,DX) + HALF*DX(IC_I_E(L)) 
            LOC(3) = LOCATION(IC_J_S(L),ZERO,DY) - HALF*DY(IC_J_S(L)) 
            LOC(4) = LOCATION(IC_J_N(L),ZERO,DY) + HALF*DY(IC_J_N(L)) 
            LOC(5) = LOCATION(IC_K_B(L),ZERO,DZ) - HALF*DZ(IC_K_B(L)) 
            LOC(6) = LOCATION(IC_K_T(L),ZERO,DZ) + HALF*DZ(IC_K_T(L)) 
            WRITE (UNIT_OUT, 1520) IC_X_W(L), LOC(1), IC_X_E(L), LOC(2), IC_Y_S&
               (L), LOC(3), IC_Y_N(L), LOC(4), IC_Z_B(L), LOC(5), IC_Z_T(L), &
               LOC(6) 
            WRITE (UNIT_OUT, 1530) IC_I_W(L), IC_I_E(L), IC_J_S(L), IC_J_N(L), &
               IC_K_B(L), IC_K_T(L) 
            WRITE (UNIT_OUT, 1540) IC_EP_G(L) 
            IF (IC_P_G(L) /= UNDEFINED) WRITE (UNIT_OUT, 1541) IC_P_G(L) 
            WRITE (UNIT_OUT, 1542) IC_T_G(L) 
            IF (SPECIES_EQ(0)) THEN 
               WRITE (UNIT_OUT, 1543) 
               DO N = 1, NMAX(0) 
                  WRITE (UNIT_OUT, 1544) N, IC_X_G(L,N) 
               END DO 
            ENDIF 
            IF (IC_GAMA_RG(L) /= ZERO) WRITE (UNIT_OUT, 1545) IC_GAMA_RG(L), &
               IC_T_RG(L) 
!
            WRITE (UNIT_OUT, 1550) IC_U_G(L), IC_V_G(L), IC_W_G(L) 
            DO M = 1, MMAX 
               WRITE (UNIT_OUT, 1560) M, IC_ROP_S(L,M) 
               WRITE (UNIT_OUT, 1561) M, IC_T_S(L,M) 
            END DO 
            DO M = 1, MMAX 
               IF (SPECIES_EQ(M)) THEN 
                  WRITE (UNIT_OUT, 1563) M 
                  DO N = 1, NMAX(M) 
                     WRITE (UNIT_OUT, 1564) N, IC_X_S(L,M,N) 
                  END DO 
               ENDIF 
            END DO 
            DO M = 1, MMAX 
               IF (IC_GAMA_RS(L,M) /= ZERO) WRITE (UNIT_OUT, 1565) M, &
                  IC_GAMA_RS(L,M), IC_T_RS(L,M) 
!
               WRITE(UNIT_OUT,1570)M,IC_U_S(L,M),M,IC_V_S(L,M),M,IC_W_S(L,M) 
            END DO 
            IF (IC_P_STAR(L) /= UNDEFINED) WRITE (UNIT_OUT, 1574) IC_P_STAR(L) 
            IF(IC_L_SCALE(L)/=UNDEFINED)WRITE(UNIT_OUT,1575)IC_L_SCALE(L) 
         ENDIF 
      END DO 
      WRITE (UNIT_OUT, 1600) 
      IF (U_G0 /= UNDEFINED) WRITE (UNIT_OUT, 1601) 'U_g (U_g0) = ', U_G0 
      IF (V_G0 /= UNDEFINED) WRITE (UNIT_OUT, 1601) 'V_g (V_g0) = ', V_G0 
      IF (W_G0 /= UNDEFINED) WRITE (UNIT_OUT, 1601) 'W_g (W_g0) = ', W_G0 
      DO M = 1, MMAX 
         IF (U_S0(M) /= UNDEFINED) WRITE (UNIT_OUT, 1602) 'U_s (U_s0[', M, &
            ']) = ', U_S0(M) 
         IF (V_S0(M) /= UNDEFINED) WRITE (UNIT_OUT, 1602) 'V_s (V_s0[', M, &
            ']) = ', V_S0(M) 
         IF (W_S0(M) /= UNDEFINED) WRITE (UNIT_OUT, 1602) 'W_s (W_s0[', M, &
            ']) = ', W_S0(M) 
      END DO 
      DO L = 1, DIMENSION_BC 
         IF (BC_DEFINED(L)) THEN 
            WRITE (UNIT_OUT, 1610) L 
            WRITE (UNIT_OUT, 1611) BC_TYPE(L) 
            SELECT CASE (TRIM(BC_TYPE(L)))
            CASE ('MASS_INFLOW')  
               WRITE (UNIT_OUT, 1612) 
            CASE ('MASS_OUTFLOW')  
               WRITE (UNIT_OUT, 1613) 
            CASE ('P_INFLOW')  
               WRITE (UNIT_OUT, 1614) 
            CASE ('P_OUTFLOW')  
               WRITE (UNIT_OUT, 1615) 
            CASE ('FREE_SLIP_WALL')  
               WRITE (UNIT_OUT, 1616) 
            CASE ('NO_SLIP_WALL')  
               WRITE (UNIT_OUT, 1617) 
            CASE ('PAR_SLIP_WALL')  
               WRITE (UNIT_OUT, 1618) 
            CASE ('OUTFLOW')  
               WRITE (UNIT_OUT, 1619) 
            END SELECT 
            LOC(1) = LOCATION(BC_I_W(L),XMIN,DX) - HALF*DX(BC_I_W(L)) 
            LOC(2) = LOCATION(BC_I_E(L),XMIN,DX) + HALF*DX(BC_I_E(L)) 
            LOC(3) = LOCATION(BC_J_S(L),ZERO,DY) - HALF*DY(BC_J_S(L)) 
            LOC(4) = LOCATION(BC_J_N(L),ZERO,DY) + HALF*DY(BC_J_N(L)) 
            LOC(5) = LOCATION(BC_K_B(L),ZERO,DZ) - HALF*DZ(BC_K_B(L)) 
            LOC(6) = LOCATION(BC_K_T(L),ZERO,DZ) + HALF*DZ(BC_K_T(L)) 
            WRITE (UNIT_OUT, 1620) BC_X_W(L), LOC(1), BC_X_E(L), LOC(2), BC_Y_S&
               (L), LOC(3), BC_Y_N(L), LOC(4), BC_Z_B(L), LOC(5), BC_Z_T(L), &
               LOC(6) 
            WRITE (UNIT_OUT, 1630) BC_I_W(L), BC_I_E(L), BC_J_S(L), BC_J_N(L), &
               BC_K_B(L), BC_K_T(L) 
            IF (BC_EP_G(L) /= UNDEFINED) WRITE (UNIT_OUT, 1640) BC_EP_G(L) 
            IF (BC_P_G(L) /= UNDEFINED) WRITE (UNIT_OUT, 1641) BC_P_G(L) 
            IF (BC_T_G(L) /= UNDEFINED) WRITE (UNIT_OUT, 1642) BC_T_G(L) 
            IF (SPECIES_EQ(0) .AND. BC_X_G(L,1)/=UNDEFINED) THEN 
               WRITE (UNIT_OUT, 1643) 
               DO N = 1, NMAX(0) 
                  WRITE (UNIT_OUT, 1644) N, BC_X_G(L,N) 
               END DO 
            ENDIF 
            IF (BC_MASSFLOW_G(L) /= UNDEFINED) WRITE (UNIT_OUT, 1648) &
               BC_MASSFLOW_G(L) 
            IF (BC_VOLFLOW_G(L) /= UNDEFINED) WRITE (UNIT_OUT, 1649) &
               BC_VOLFLOW_G(L) 
            IF (BC_U_G(L) /= UNDEFINED) WRITE (UNIT_OUT, 1650) BC_U_G(L) 
            IF (BC_V_G(L) /= UNDEFINED) WRITE (UNIT_OUT, 1651) BC_V_G(L) 
            IF (BC_W_G(L) /= UNDEFINED) WRITE (UNIT_OUT, 1652) BC_W_G(L) 
            IF (BC_DT_0(L) /= UNDEFINED) THEN 
               IF (BC_JET_G0(L) /= UNDEFINED) THEN 
                  WRITE (UNIT_OUT, 1655) BC_DT_0(L), BC_JET_G0(L), BC_DT_L(L), &
                     BC_JET_GL(L), BC_DT_H(L), BC_JET_GH(L) 
               ELSE 
                  WRITE (UNIT_OUT, 1656) BC_DT_0(L) 
               ENDIF 
            ENDIF 
            DO M = 1, MMAX 
               IF (BC_ROP_S(L,M) /= UNDEFINED) THEN 
                  WRITE (UNIT_OUT, 1660) M, BC_ROP_S(L,M) 
                  WRITE (UNIT_OUT, 1661) M, BC_T_S(L,M) 
               ENDIF 
            END DO 
            DO M = 1, MMAX 
               IF (SPECIES_EQ(M) .AND. BC_X_S(L,M,1)/=UNDEFINED) THEN 
                  WRITE (UNIT_OUT, 1663) M 
                  DO N = 1, NMAX(M) 
                     WRITE (UNIT_OUT, 1664) N, BC_X_S(L,M,N) 
                  END DO 
               ENDIF 
            END DO 
            DO M = 1, MMAX 
               IF (BC_MASSFLOW_S(L,M) /= UNDEFINED) WRITE (UNIT_OUT, 1668) M, &
                  BC_MASSFLOW_S(L,M) 
               IF (BC_VOLFLOW_S(L,M) /= UNDEFINED) WRITE (UNIT_OUT, 1669) M, &
                  BC_VOLFLOW_S(L,M) 
               IF(BC_U_S(L,M)/=UNDEFINED)WRITE(UNIT_OUT,1670)M,BC_U_S(L,M) 
               IF(BC_V_S(L,M)/=UNDEFINED)WRITE(UNIT_OUT,1671)M,BC_V_S(L,M) 
               IF(BC_W_S(L,M)/=UNDEFINED)WRITE(UNIT_OUT,1672)M,BC_W_S(L,M) 
            END DO 
            IF (BC_TYPE(L) == 'PAR_SLIP_WALL') THEN 
               WRITE (UNIT_OUT, 1675) BC_HW_G(L), BC_UW_G(L), BC_VW_G(L), &
                  BC_WW_G(L) 
               DO M = 1, MMAX 
                  WRITE (UNIT_OUT, 1676) M, BC_HW_S(L,M), BC_UW_S(L,M), BC_VW_S&
                     (L,M), BC_WW_S(L,M) 
               END DO 
            ENDIF 
         ENDIF 
      END DO 
      WRITE (UNIT_OUT, 1700) 
      DO L = 1, DIMENSION_IS 
         IF (IS_DEFINED(L)) THEN 
            WRITE (UNIT_OUT, 1710) L 
            WRITE (UNIT_OUT, 1711) IS_TYPE(L) 
            IF (IS_TYPE(L)=='IMPERMEABLE' .OR. IS_TYPE(L)(3:13)=='IMPERMEABLE'&
               ) THEN 
               WRITE (UNIT_OUT, 1712) 
            ELSE IF (IS_TYPE(L)=='SEMIPERMEABLE' .OR. IS_TYPE(L)(3:15)==&
                  'SEMIPERMEABLE') THEN 
               WRITE (UNIT_OUT, 1713) 
            ENDIF 
            LOC(1) = LOCATION(IS_I_W(L),XMIN,DX) - HALF*DX(IS_I_W(L)) 
            LOC(2) = LOCATION(IS_I_E(L),XMIN,DX) + HALF*DX(IS_I_E(L)) 
            LOC(3) = LOCATION(IS_J_S(L),ZERO,DY) - HALF*DY(IS_J_S(L)) 
            LOC(4) = LOCATION(IS_J_N(L),ZERO,DY) + HALF*DY(IS_J_N(L)) 
            LOC(5) = LOCATION(IS_K_B(L),ZERO,DZ) - HALF*DZ(IS_K_B(L)) 
            LOC(6) = LOCATION(IS_K_T(L),ZERO,DZ) + HALF*DZ(IS_K_T(L)) 
            WRITE (UNIT_OUT, 1720) IS_X_W(L), LOC(1), IS_X_E(L), LOC(2), IS_Y_S&
               (L), LOC(3), IS_Y_N(L), LOC(4), IS_Z_B(L), LOC(5), IS_Z_T(L), &
               LOC(6) 
            WRITE (UNIT_OUT, 1730) IS_I_W(L), IS_I_E(L), IS_J_S(L), IS_J_N(L), &
               IS_K_B(L), IS_K_T(L) 
            IF (IS_PC(L,1) /= UNDEFINED) WRITE (UNIT_OUT, 1740) IS_PC(L,1) 
            IF (IS_PC(L,2) /= UNDEFINED) WRITE (UNIT_OUT, 1741) IS_PC(L,2) 
            DO M = 1, MMAX 
               WRITE (UNIT_OUT, 1742) M, IS_VEL_S(L,M) 
            END DO 
         ENDIF 
      END DO 
      WRITE (UNIT_OUT, 1800) OUT_DT, RES_DT, SPX_DT 
!
!  Print out tolerance values from TOLERANCE.INC
!
      WRITE (UNIT_OUT, 1900) 
      WRITE (UNIT_OUT, 1901) ZERO_EP_S 
      WRITE (UNIT_OUT, 1904) TOL_RESID, TOL_RESID_T, TOL_RESID_X, TOL_DIVERGE 
      WRITE (UNIT_OUT, 1905) TOL_COM
      IF(NScalar /= 0)WRITE (UNIT_OUT, 1906) TOL_RESID_Scalar
      IF(K_Epsilon)WRITE (UNIT_OUT, 1907) TOL_RESID_K_Epsilon
!
!  Echo user defined input data
!
      WRITE (UNIT_OUT, '(/,1X,1A1)') CHAR(12) 
      IF (CALL_USR) CALL USR_WRITE_OUT0 
!
      RETURN  
 1000 FORMAT(17X,'MM      MM  FFFFFFFFFF    IIIIII    XX      XX',/17X,&
         'MM      MM  FFFFFFFFFF    IIIIII    XX      XX',/17X,&
         'MMMM  MMMM  FF              II      XX      XX',/17X,&
         'MMMM  MMMM  FF              II      XX      XX',/17X,&
         'MM  MM  MM  FF              II        XX  XX  ',/17X,&
         'MM  MM  MM  FF              II        XX  XX  ',/17X,&
         'MM      MM  FFFFFFFF        II          XX    ',/17X,&
         'MM      MM  FFFFFFFF        II          XX    ',/17X,&
         'MM      MM  FF              II        XX  XX  ',/17X,&
         'MM      MM  FF              II        XX  XX  ',/17X,&
         'MM      MM  FF              II      XX      XX',/17X,&
         'MM      MM  FF              II      XX      XX',/17X,&
         'MM      MM  FF            IIIIII    XX      XX',/17X,&
         'MM      MM  FF            IIIIII    XX      XX',2/20X,&
         'Multiphase Flow with Interphase eXchanges'/34X,'Version: ',A,/20X,&
         'Time: ',I2,':',I2,20X,'Date: ',I2,'-',I2,'-',I4) 
 1010 FORMAT(/7X,'Computer : ',A50,/,1X,79('_')) 
 1100 FORMAT(//,3X,'1. RUN CONTROL',/) 
 1101 FORMAT(/7X,'* Model B momentum equations are solved')
 1102 FORMAT(/7X,'Number of scalars = ', I4,&
             /7X,'Scalar No.        Carrier Phase (Phase4Scalar)')
 1103 FORMAT(/7X, I4,'               ',I4)
 1104 FORMAT(/7X,'* K and Epsilon equations are solved')
 1110 FORMAT(7X,'Run name(RUN_NAME): ',A60) 
 1120 FORMAT(7X,'Brief description of the run (DESCRIPTION) :',/9X,A60) 
 1130 FORMAT(7X,'Units (UNITS) : ',A16) 
 1135 FORMAT(7X,'Start-time (TIME) = ',G12.5,/7X,'Stop_time (TSTOP) = ',G12.5,/&
         7X,'Time step (DT) = ',G12.5,/7X,'Max time step (DT_MAX) = ',G12.5,/7X&
         ,'Min time step (DT_MIN) = ',G12.5,/7X,&
         'Time step adjustment factor (DT_FAC) = ',G12.5) 
 1136 FORMAT(7X,'* Steady state simulation.') 
 1137 FORMAT(7X,'Type of run (RUN_TYPE) : ',A16) 
 1138 FORMAT(30X,'(Initial conditions from the input (.DAT) file)') 
 1139 FORMAT(30X,'(Initial conditions from the restart (.RES) file)') 
 1140 FORMAT(/7X,'* Gas momentum equation-',A,' is',A,'solved.') 
 1141 FORMAT(/7X,'* Solids-',I1,' momentum equation-',A,' is',A,'solved.') 
 1143 FORMAT(/7X,'* Energy equations are solved.') 
 1144 FORMAT(/7X,'* Energy equations are NOT solved.') 
 1145 FORMAT(/7X,'* Gas Species equations are solved.') 
 1146 FORMAT(/7X,'* Gas Species equations are NOT solved.') 
 1147 FORMAT(/7X,'* Solids-',I1,' Species equations are solved.') 
 1148 FORMAT(/7X,'* Solids-',I1,' Species equations are NOT solved.') 
 1149 FORMAT(/7X,'* User-defined subroutines are',A,'called.') 
!
 1150 FORMAT(//,3X,'2. PHYSICAL AND NUMERICAL PARAMETERS',/) 
 1151 FORMAT(7X,'Coefficient of restitution (C_e) = ',G12.5) 
 1152 FORMAT(7X,'Coefficient of friction (C_f) = ',G12.5) 
 1153 FORMAT(7X,'Angle of internal friction (Phi) = ',G12.5) 
 1154 FORMAT(7X,'Angle of wall_particle friction (Phi_w) = ',G12.5) 
 1155 FORMAT(7X,'Default turbulence length scale (L_scale0) = ',G12.5,/7X,&
         'Maximum turbulent viscosity (MU_gmax) = ',G12.5) 
 1156 FORMAT(7X,'Excluded volume for B-M stress term (V_ex) = ',G12.5) 
 1157 FORMAT(7X,'Reference pressure (P_ref) = ',G12.5,/7X,&
         'Pressure scale-factor (P_scale) = ',G12.5,/7X,&
         'Gravitational acceleration (GRAVITY) = ',G12.5) 
 1158 FORMAT(7X,'Under relaxation (UR_FAC) and',&
         ' Iterations in Leq solver (LEQ_IT):'/,9X,&
         '                        UR_FAC',2X,'LEQ_IT','  LEQ_METHOD',&
         '  LEQ_SWEEP', '  LEQ_TOL', '  DISCRETIZE') 
 1159 FORMAT(9X,&
         'Fluid cont. and P_g   = ',F5.3,2X,I4,5X,I4,9x,A4,4X,G11.4,2X,A12/9X,&
         'Solids cont. and P_s  = ',F5.3,2X,I4,5X,I4,9x,A4,4X,G11.4,2X,A12/9X,&
         'U velocity            = ',F5.3,2X,I4,5X,I4,9x,A4,4X,G11.4,2X,A12/9X,&
         'V velocity            = ',F5.3,2X,I4,5X,I4,9x,A4,4X,G11.4,2X,A12/9X,&
         'W velocity            = ',F5.3,2X,I4,5X,I4,9x,A4,4X,G11.4,2X,A12/9X,&
         'Energy                = ',F5.3,2X,I4,5X,I4,9x,A4,4X,G11.4,2X,A12/9X,&
         'Species               = ',F5.3,2X,I4,5X,I4,9x,A4,4X,G11.4,2X,A12/9X,&
         'Granular Energy       = ',F5.3,2X,I4,5X,I4,9x,A4,4X,G11.4,2X,A12/9X,&
         'User scalar           = ',F5.3,2X,I4,5X,I4,9x,A4,4X,G11.4,2X,A12/) 
 1190 FORMAT(7X,1A20,'- C(',I2,') = ',G12.5) 
!
 1200 FORMAT(//,3X,'3. GEOMETRY AND DISCRETIZATION',/) 
 1201 FORMAT(7X,'Coordinates: ',1A16/) 
 1202 FORMAT(7X,'Cyclic boundary conditions in ',A,' direction',A) 
 1203 FORMAT(7X,'Pressure drop (DELP_',A,') = ',G12.5) 
 1210 FORMAT(7X,'X-direction cell sizes (DX) and East face locations:') 
 1211 FORMAT(7X,'Minimum value of X, or R (XMIN) =',G12.5) 
 1212 FORMAT(7X,'Number of cells in X, or R, direction (IMAX) = ',I4) 
 1213 FORMAT(7X,'Reactor length in X, or R, direction (XLENGTH) =',G12.5//) 
 1220 FORMAT(7X,'Y-direction cell sizes (DY) and North face locations:') 
 1221 FORMAT(7X,'Number of cells in Y direction (JMAX) = ',I4) 
 1222 FORMAT(7X,'Reactor length in Y direction (YLENGTH) =',G12.5//) 
 1230 FORMAT(7X,'Z-direction cell sizes (DZ) and Top face locations:') 
 1231 FORMAT(7X,'Number of cells in Z, or theta, direction (KMAX) = ',I4) 
 1232 FORMAT(7X,'Reactor length in Z, or theta, direction (ZLENGTH) =',G12.5) 
!
 1300 FORMAT(//,3X,'4. GAS PHASE',/) 
 1305 FORMAT(7X,'Gas density (RO_g0) = ',G12.5,&
         '  (A constant value is used everywhere)') 
 1310 FORMAT(7X,'Viscosity (MU_g0) = ',G12.5,&
         '  (A constant value is used everywhere)') 
 1315 FORMAT(7X,'Number of gas species (NMAX(0)) = ',I3) 
 1316 FORMAT(7X,'Gas species',5X,'Molecular weight (MW_g)') 
 1317 FORMAT(7X,3X,I3,15X,G12.5) 
 1320 FORMAT(7X,'Average molecular weight (MW_avg) = ',G12.5,&
         '  (A constant value is used everywhere)') 
!
 1400 FORMAT(//,3X,'5. SOLIDS PHASE',/) 
 1405 FORMAT(7X,'Viscosity (MU_s0) = ',G12.5,&
         '  (A constant value is used everywhere)') 
 1410 FORMAT(7X,'Number of particulate phases (MMAX) = ',I2) 
 1420 FORMAT(/7X,'M',5X,'Diameter (D_p)',T35,'Density (RO_s)',T50,&
         'Close_Packed') 
 1421 FORMAT(7X,I1,5X,G12.5,T35,G12.5,T55,L1) 
 1422 FORMAT(7X,'Number of solids-',I1,' species (NMAX(',I1,')) = ',I3) 
 1423 FORMAT(7X,'Solids species',5X,'Molecular weight (MW_s)') 
 1424 FORMAT(7X,3X,I3,18X,G12.5) 
 1430 FORMAT(/7X,'Void fraction at maximum packing (EP_star) = ',G12.5) 
!
 1500 FORMAT(//,3X,'6. INITIAL CONDITIONS') 
 1510 FORMAT(/7X,'Initial condition no : ',I4) 
 1520 FORMAT(9X,39X,' Specified  ',5X,' Simulated  ',/9X,&
         'X coordinate of west face   (IC_X_w) = ',G12.5,5X,G12.5/,9X,&
         'X coordinate of east face   (IC_X_e) = ',G12.5,5X,G12.5/,9X,&
         'Y coordinate of south face  (IC_Y_s) = ',G12.5,5X,G12.5/,9X,&
         'Y coordinate of north face  (IC_Y_n) = ',G12.5,5X,G12.5/,9X,&
         'Z coordinate of bottom face (IC_Z_b) = ',G12.5,5X,G12.5/,9X,&
         'Z coordinate of top face    (IC_Z_t) = ',G12.5,5X,G12.5) 
 1530 FORMAT(9X,'I index of cell at west   (IC_I_w) = ',I4,/,9X,&
         'I index of cell at east   (IC_I_e) = ',I4,/,9X,&
         'J index of cell at south  (IC_J_s) = ',I4,/,9X,&
         'J index of cell at north  (IC_J_n) = ',I4,/,9X,&
         'K index of cell at bottom (IC_K_b) = ',I4,/,9X,&
         'K index of cell at top    (IC_K_t) = ',I4) 
 1540 FORMAT(9X,'Void fraction (IC_EP_g) = ',G12.5) 
 1541 FORMAT(9X,'Gas pressure (IC_P_g) = ',G12.5) 
 1542 FORMAT(9X,'Gas temperature (IC_T_g) = ',G12.5) 
 1543 FORMAT(9X,'Gas species',5X,'Mass fraction (IC_X_g)') 
 1544 FORMAT(9X,3X,I3,15X,G12.5) 
 1545 FORMAT(9X,'Gas radiation coefficient   (IC_GAMA_Rg) = ',G12.5,/,9X,&
         'Gas radiation temperature   (IC_T_Rg) = ',G12.5) 
 1550 FORMAT(9X,'X-component of gas velocity (IC_U_g) = ',G12.5,/9X,&
         'Y-component of gas velocity (IC_V_g) = ',G12.5,/9X,&
         'Z-component of gas velocity (IC_W_g) = ',G12.5) 
 1560 FORMAT(9X,'Solids phase-',I1,' Density x Volume fr. (IC_ROP_s) = ',G12.5) 
 1561 FORMAT(9X,'Solids phase-',I1,' temperature (IC_T_s) = ',G12.5) 
 1563 FORMAT(9X,'Solids-',I1,' species',5X,'Mass fraction (IC_X_s)') 
 1564 FORMAT(9X,3X,I3,20X,G12.5) 
 1565 FORMAT(9X,'Solids phase-',I1,' radiation coefficient (IC_GAMA_Rs)',' =',&
         G12.5,/9X,'Solids phase-',I1,' radiation temperature (IC_T_Rs) =',&
         G12.5) 
 1570 FORMAT(9X,'X-component of solids phase-',I1,' velocity (IC_U_s) =',G12.5,&
         /9X,'Y-component of solids phase-',I1,' velocity (IC_V_s) =',G12.5,/9X&
         ,'Z-component of solids phase-',I1,' velocity (IC_W_s) =',G12.5) 
 1574 FORMAT(9X,'Solids pressure (IC_P_star) = ',G12.5) 
 1575 FORMAT(9X,'Turbulence length scale (IC_L_scale) = ',G12.5) 
!
 1600 FORMAT(//,3X,'7. BOUNDARY CONDITIONS') 
 1601 FORMAT(/7X,'Average value of ',A,G12.5) 
 1602 FORMAT(/7X,'Average value of ',A,I2,A,G12.5) 
 1610 FORMAT(/7X,'Boundary condition no : ',I4) 
 1611 FORMAT(9X,'Type of boundary condition : ',A16) 
 1612 FORMAT(11X,'(Inlet with specified gas and solids mass flux)') 
 1613 FORMAT(11X,'(Outlet with specified gas and solids mass flux)') 
 1614 FORMAT(11X,'(Inlet with specified gas pressure)') 
 1615 FORMAT(11X,'(Outlet with specified gas pressure)') 
 1616 FORMAT(11X,'(Gradients of parallel velocity components are zero)') 
 1617 FORMAT(11X,'(Velocity is zero at wall)') 
 1618 FORMAT(11X,'(Partial slip condition at wall)') 
 1619 FORMAT(11X,'(Outflow condition)') 
 1620 FORMAT(9X,39X,' Specified  ',5X,' Simulated  ',/9X,&
         'X coordinate of west face   (BC_X_w) = ',G12.5,5X,G12.5/,9X,&
         'X coordinate of east face   (BC_X_e) = ',G12.5,5X,G12.5/,9X,&
         'Y coordinate of south face  (BC_Y_s) = ',G12.5,5X,G12.5/,9X,&
         'Y coordinate of north face  (BC_Y_n) = ',G12.5,5X,G12.5/,9X,&
         'Z coordinate of bottom face (BC_Z_b) = ',G12.5,5X,G12.5/,9X,&
         'Z coordinate of top face    (BC_Z_t) = ',G12.5,5X,G12.5) 
 1630 FORMAT(9X,'I index of cell at west   (BC_I_w) = ',I4,/,9X,&
         'I index of cell at east   (BC_I_e) = ',I4,/,9X,&
         'J index of cell at south  (BC_J_s) = ',I4,/,9X,&
         'J index of cell at north  (BC_J_n) = ',I4,/,9X,&
         'K index of cell at bottom (BC_K_b) = ',I4,/,9X,&
         'K index of cell at top    (BC_K_t) = ',I4) 
 1640 FORMAT(9X,'Void fraction (BC_EP_g) = ',G12.5) 
 1641 FORMAT(9X,'Gas pressure (BC_P_g) = ',G12.5) 
 1642 FORMAT(9X,'Gas temperature (BC_T_g) = ',G12.5) 
 1643 FORMAT(9X,'Gas species',5X,'Mass fraction (BC_X_g)') 
 1644 FORMAT(9X,3X,I3,15X,G12.5) 
 1648 FORMAT(9X,'Gas mass flow rate (BC_MASSFLOW_g) = ',G12.5) 
 1649 FORMAT(9X,'Gas volumetric flow rate (BC_VOLFLOW_g) = ',G12.5) 
 1650 FORMAT(9X,'X-component of gas velocity (BC_U_g) = ',G12.5) 
 1651 FORMAT(9X,'Y-component of gas velocity (BC_V_g) = ',G12.5) 
 1652 FORMAT(9X,'Z-component of gas velocity (BC_W_g) = ',G12.5) 
 1655 FORMAT(9X,'Initial interval when jet vel= BC_Jet_g0 (BC_DT_0) = ',G12.5,/&
         9X,'Initial jet velocity (BC_Jet_g0) = ',G12.5,/9X,&
         'Interval when jet vel= BC_Jet_gl (BC_DT_l) = ',G12.5,/9X,&
         'Low value of jet velocity (BC_Jet_gl) = ',G12.5,/9X,&
         'Interval when jet vel = BC_Jet_gh (BC_DT_h) = ',G12.5,/9X,&
         'High value of jet velocity (BC_Jet_gh) = ',G12.5) 
 1656 FORMAT(9X,'Interval for averaging outflow rates= (BC_DT_0) = ',G12.5) 
 1660 FORMAT(9X,'Solids phase-',I1,' Density x Volume fr. (BC_ROP_s) = ',G12.5) 
 1661 FORMAT(9X,'Solids phase-',I1,' temperature (BC_T_s) = ',G12.5) 
 1663 FORMAT(9X,'Solids-',I1,' species',5X,'Mass fraction (BC_X_s)') 
 1664 FORMAT(9X,3X,I3,20X,G12.5) 
 1668 FORMAT(9X,'Solids phase-',I1,' mass flow rate (BC_MASSFLOW_s) =',G12.5) 
 1669 FORMAT(9X,'Solids phase-',I1,' volumetric flow rate (BC_VOLFLOW_s) =',&
         G12.5) 
 1670 FORMAT(9X,'X-component of solids phase-',I1,' velocity (BC_U_s) =',G12.5) 
 1671 FORMAT(9X,'Y-component of solids phase-',I1,' velocity (BC_V_s) =',G12.5) 
 1672 FORMAT(9X,'Z-component of solids phase-',I1,' velocity (BC_W_s) =',G12.5) 
 1675 FORMAT(9X,'Partial slip coefficient   (BC_hw_g) = ',G12.5,/,9X,&
         'Slip velociity U at wall   (BC_Uw_g) = ',G12.5,/,9X,&
         'Slip velociity V at wall   (BC_Vw_g) = ',G12.5,/,9X,&
         'Slip velociity W at wall   (BC_Ww_g) = ',G12.5) 
 1676 FORMAT(9X,'Solids phase: ',I1,/,11X,&
         'Partial slip coefficient   (BC_hw_s) = ',G12.5,/,11X,&
         'Slip velociity U at wall   (BC_Uw_s) = ',G12.5,/,11X,&
         'Slip velociity V at wall   (BC_Vw_s) = ',G12.5,/,11X,&
         'Slip velociity W at wall   (BC_Ww_s) = ',G12.5) 
!
 1700 FORMAT(//,3X,'8. INTERNAL SURFACES') 
 1710 FORMAT(/7X,'Internal surface no : ',I4) 
 1711 FORMAT(9X,'Type of internal surface : ',A16) 
 1712 FORMAT(11X,'(No gas or solids flow through the surface)') 
 1713 FORMAT(11X,'(Only gas flows through the surface)') 
 1720 FORMAT(9X,39X,' Specified  ',5X,' Simulated  ',/9X,&
         'X coordinate of west face   (IS_X_w) = ',G12.5,5X,G12.5/,9X,&
         'X coordinate of east face   (IS_X_e) = ',G12.5,5X,G12.5/,9X,&
         'Y coordinate of south face  (IS_Y_s) = ',G12.5,5X,G12.5/,9X,&
         'Y coordinate of north face  (IS_Y_n) = ',G12.5,5X,G12.5/,9X,&
         'Z coordinate of bottom face (IS_Z_b) = ',G12.5,5X,G12.5/,9X,&
         'Z coordinate of top face    (IS_Z_t) = ',G12.5,5X,G12.5) 
 1730 FORMAT(9X,'I index of cell at west   (IS_I_w) = ',I4,/,9X,&
         'I index of cell at east   (IS_I_e) = ',I4,/,9X,&
         'J index of cell at south  (IS_J_s) = ',I4,/,9X,&
         'J index of cell at north  (IS_J_n) = ',I4,/,9X,&
         'K index of cell at bottom (IS_K_b) = ',I4,/,9X,&
         'K index of cell at top    (IS_K_t) = ',I4) 
 1740 FORMAT(9X,'Permeability (IS_PC1) = ',G12.5) 
 1741 FORMAT(9X,'Inertial resistance factor (IS_PC2) = ',G12.5) 
 1742 FORMAT(9X,'Solids phase-',I2,' Velocity (IS_VEL_s) = ',G12.5) 
!
 1800 FORMAT(//,3X,'9. OUTPUT DATA FILES:',/7X,'Extension',T18,'Description',&
         T59,'Interval for writing',/7X,'.OUT',T18,'This file (ASCII)',T61,&
         G12.5,/7X,'.LOG',T18,'Log file containing messages (ASCII)',&
         /7X,'.RES',&
         T18,'Restart file (Binary)',T61,G12.5,/7X,'.SP1',T18,&
         'EP_g (Binary, single precision)',T61,G12.5,/7X,'.SP2',T18,&
         'P_g, P_star (Binary, single precision)',T61,G12.5,/7X,'.SP3',T18,&
         'U_g, V_g, W_g (Binary, single precision)',T61,G12.5,/7X,'.SP4',T18,&
         'U_s, V_s, W_s (Binary, single precision)',T61,G12.5,/7X,'.SP5',T18,&
         'ROP_s (Binary, single precision)',T61,G12.5,/7X,'.SP6',T18,&
         'T_g, T_s (Binary, single precision)',T61,G12.5,/7X,'.SP7',T18,&
         'X_g, X_s (Binary, single precision)',T61,G12.5,/7X,'.SP8',T18,&
         'Theta_m (Binary, single precision)',T61,G12.5,/7X,'.SP9',T18,&
         'User Scalar (Binary, single precision)',T61,G12.5) 
!
 1900 FORMAT(//,3X,'10. TOLERANCES',/7X,&
         'The following values are specified in the file TOLERANCE.INC.') 
 1901 FORMAT(/7X,'Minimum value of EP_s tracked (ZERO_EP_s) = ',G12.5) 
 1904 FORMAT(7X,'Maximum average residual (TOL_RESID) = ',G12.5,/7X,&
         'Maximum average residual (TOL_RESID_T) = ',G12.5,/7X,&
         'Maximum average residual (TOL_RESID_X) = ',G12.5,/7X,&
         'Minimum residual at divergence (TOL_DIVERGE) = ',G12.5) 
 1905 FORMAT(7X,'Tolerance for species and energy balances (TOL_COM) = ',G12.5) 
 1906 FORMAT(7X,'Tolerance for scalar mass balances (TOL_RESID_Scalar) = ',G12.5)  
 1907 FORMAT(7X,'Tolerance for K-Epsilon balances (TOL_RESID_K_Epsilon) = ',G12.5)
!
      END SUBROUTINE WRITE_OUT0 
      
      SUBROUTINE WRITE_FLAGS
      USE param
      USE param1
      USE funits
      USE geometry
      USE indices
      USE compar         !//d
      USE mpi_utility    !//d
      USE sendrecv    !//d
      IMPLICIT NONE
      integer ijk
!
      character*3, allocatable :: array1(:)   !//d
      character*4, dimension(:), allocatable :: array2, array3
      include 'function.inc'
      
  
      if (myPE .eq. PE_IO) then
         allocate (array1(ijkmax3))
         allocate (array2(dimension_3))
         allocate (array3(ijkmax3))
      else
         allocate (array1(1))
         allocate (array2(dimension_3))
         allocate (array3(1))
      end if

!write(*,*) 'ijkmax3', ijkmax3, dimension_3
      
!//SP Filling the processor ghost layer with the correct values

      call gather (icbc_flag,array1,PE_IO)
      call scatter (icbc_flag,array1,PE_IO)
      
!
!  Superimpose internal surface flags on Initial and boundary condition flags
!
      DO ijk = IJKSTART3, IJKEND3
        array2(ijk) = '    '
        array2(ijk)(1:3) = icbc_flag(ijk)(1:3)
        IF (IP_AT_E(IJK)) THEN 
           array2(IJK)(4:4) = 'E' 
        ELSE IF (SIP_AT_E(IJK)) THEN 
           array2(IJK)(4:4) = 'e' 
        ENDIF 
!
        IF (IP_AT_N(IJK)) THEN 
           array2(IJK)(4:4) = 'N' 
        ELSE IF (SIP_AT_N(IJK)) THEN 
           array2(IJK)(4:4) = 'n' 
        ENDIF 
!
        IF (IP_AT_T(IJK)) THEN 
           array2(IJK)(4:4) = 'T' 
        ELSE IF (SIP_AT_T(IJK)) THEN 
           array2(IJK)(4:4) = 't' 
        ENDIF 
      ENDDO
      call gather (array2,array3,PE_IO)
      
      if(myPE.eq.PE_IO) then
        WRITE (UNIT_OUT, 2000) CHAR(12) 
        CALL OUT_ARRAY_C (array3, 'BC/IC condition flags') 
        WRITE (UNIT_OUT, *)
      ENDIF
      

      deallocate (array1) 
      deallocate (array2) 
      deallocate (array3) 
!
 2000 FORMAT(//,3X,'11. INITIAL AND BOUNDARY CONDITION FLAGS',/7X,&
         'The initial and boundary conditions specified are shown in',/7X,&
         'the following map. Each computational cell is represented',/7X,&
         'by a string of three characters.  The first character',/7X,&
         'represents the type of cell, and the last two characters',/7X,&
         'give a number that identifies a boundary or initial condi-',/7X,&
         'tion.  For example, .02 indicates a cell where Initial',/7X,&
         'Condition No. 2 will be specified. Only the last two digits'/7X,&
         'are written.  Hence, for example, Condition No. 12, 112, 212'/7X,&
         'etc. will be represented only as 12.',/7X,&
         '  First Character       Description'/7X,&
         '       .                Initial condition'/7X,&
         '       W                No slip wall'/7X,&
         '       S                Free-slip wall'/7X,&
         '       s                Partial-slip wall'/7X,&
         '       c                Cyclic boundary'/7X,&
         '       C                Cyclic boundary with pressure drop'/7X,&
         '       I                Specified mass-flux inflow cell'/7X,&
         '       O                Outflow cell'/7X,&
         '       p                Specified pressure inflow cell'/7X,&
         '       P                Specified pressure outflow cell'/7X,&
         '                                                       '/7X,&
         'Internal surfaces at East, North or Top of each cell is',/7X,&
         'is represented by the following letters to the right of the',/7X,&
         'three-character string:',/7X,&
         '  Side          Impermeable           Semipermeable',/7X,&
         '  East             E                       e       ',/7X,&
         '  North            N                       n       ',/7X,&
         '  Top              T                       t       ',/7X,&
         'For cells with internal surfaces on more than one side',/7X,&
         'the characters will be over-written in the above order',/1X,A1) 
	 RETURN
	 END SUBROUTINE WRITE_FLAGS

      
