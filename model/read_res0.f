!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: READ_RES0                                              C
!  Purpose: read the initial restart records (namelist data)           C
!                                                                      C
!  Author: P. Nicoletti                               Date: 13-DEC-91  C
!  Reviewer: P. Nicoletti, W. Rogers, M. Syamlal      Date:            C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables modified : RUN_NAME  ,  ID_MONTH  ,  ID_DAY , ID_YEAR     C
!                        ID_HOUR, ID_MINUTE, ID_SECOND, IMAX, JMAX     C
!                        KMAX, IMAX1, JMAX1, KMAX1, IMAX2, JMAX2,KMAX2 C
!                        IJMAX2, IJKMAX2, MMAX, DT, XLENGTH, YLENGTH   C
!                        ZLENGTH, DX, DY, DZ, RUN_NAME, DESCRIPTION    C
!                        UNITS, RUN_TYPE, CORDINATES, D_p, RO_s,       C
!                        EP_star, MU_g, MW_AVG, IC_X_w, IC_X_e, IC_Y_s C
!                        IC_Y_n, IC_Z_b, IC_Z_t, IC_I_w, IC_I_e        C
!                        IC_J_s, IC_J_n, IC_K_b, IC_K_t, IC_EP_g       C
!                        IC_P_g, IC_T_g, IC_T_s,  IC_U_g      C
!                        IC_V_g, IC_W_g, IC_ROP_s, IC_U_s, IC_V_s      C
!                        IC_W_s, BC_X_w, BC_X_e, BC_Y_s, BC_Y_n        C
!                        BC_Z_b, BC_Z_t, BC_I_w, BC_I_e, BC_J_s        C
!                        BC_K_b, BC_K_t, BC_EP_g, BC_P_g, BC_T_g       C
!                        BC_T_s,  BC_U_g, BC_V_g, BC_W_g      C
!                        BC_RO_g, BC_ROP_g, BC_VOLFLOW_g,BC_MASSFLOW_g C
!                        BC_ROP_s, BC_U_s, BC_V_s, BC_VOLFLOW_s        C
!                        BC_MASSFLOW_s, BC_TYPE, FLAG                  C
!  Variables referenced: None                                          C
!                                                                      C
!  Local variables: LC, L, NEXT_RECA, VERSION                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE READ_RES0 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE geometry
      USE physprop
      USE run
      USE ic
      USE bc
      USE is
      USE constant
      USE funits 
      USE output
      USE scales 
      USE ur_facs 
      USE toleranc 
      USE leqsol 
      USE tmp_array
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
!                loop counters
      INTEGER    LC, L , N, M
!
!                Pointer to the next record
      INTEGER    NEXT_RECA
!
!                file version id
      CHARACTER  VERSION*512
!
!                version number
      REAL       VERSION_NUMBER
!
!                      Temporary arrays
      DOUBLE PRECISION IC_Tmp(DIMENSION_IC), BC_Tmp(DIMENSION_BC)
!
      INTEGER    DIM_IC , DIM_BC , DIM_C , DIM_IS
!-----------------------------------------------

      call lock_tmp_array

!
!
!  1) Check to ensure that this subroutine was updated.
!  2) Initialize missing constants from earlier versions.
!  3) Add new read statements at the end of the file.
!
      READ (UNIT_RES, REC=1) VERSION 
      READ (VERSION(6:512), *) VERSION_NUMBER 
      IF (VERSION_NUMBER > 1.2) THEN 
         WRITE (*, *) ' Update Subroutine read_res0' 
         CALL SLUMBER 
         STOP  
      ENDIF 
!
!  Initialize required constants missing from earlier versions
      P_REF = ZERO 
      P_SCALE = ONE 
      DIM_IC = 5 
      DIM_BC = 5 
      DIM_C = 5 
      DIM_IS = 5 
      C_E = 1.0 
      C_F = 0.0 
      PHI = 0.0 
      PHI_W = 0.0 
!
!
      READ (UNIT_RES, REC=2) RUN_NAME, ID_MONTH, ID_DAY, ID_YEAR, ID_HOUR, &
         ID_MINUTE, ID_SECOND 
      READ (UNIT_RES, REC=3) NEXT_RECA 
      IF (VERSION == 'RES = 01.00') THEN 
         READ (UNIT_RES, REC=4) IMIN1, JMIN1, KMIN1, IMAX, JMAX, KMAX, IMAX1, &
            JMAX1, KMAX1, IMAX2, JMAX2, KMAX2, IJMAX2, IJKMAX2, MMAX, DT, &
            XLENGTH, YLENGTH, ZLENGTH 
      ELSE IF (VERSION=='RES = 01.01' .OR. VERSION=='RES = 01.02') THEN 
         READ (UNIT_RES, REC=4) IMIN1, JMIN1, KMIN1, IMAX, JMAX, KMAX, IMAX1, &
            JMAX1, KMAX1, IMAX2, JMAX2, KMAX2, IJMAX2, IJKMAX2, MMAX, DIM_IC, &
            DIM_BC, DT, XLENGTH, YLENGTH, ZLENGTH 
      ELSE IF (VERSION == 'RES = 01.03') THEN 
         READ (UNIT_RES, REC=4) IMIN1, JMIN1, KMIN1, IMAX, JMAX, KMAX, IMAX1, &
            JMAX1, KMAX1, IMAX2, JMAX2, KMAX2, IJMAX2, IJKMAX2, MMAX, DIM_IC, &
            DIM_BC, DT, XMIN, XLENGTH, YLENGTH, ZLENGTH 
      ELSE IF (VERSION == 'RES = 01.04') THEN 
         READ (UNIT_RES, REC=4) IMIN1, JMIN1, KMIN1, IMAX, JMAX, KMAX, IMAX1, &
            JMAX1, KMAX1, IMAX2, JMAX2, KMAX2, IJMAX2, IJKMAX2, MMAX, DIM_IC, &
            DIM_BC, DIM_C, DT, XMIN, XLENGTH, YLENGTH, ZLENGTH 
      ELSE IF (VERSION == 'RES = 01.05') THEN 
         READ (UNIT_RES, REC=4) IMIN1, JMIN1, KMIN1, IMAX, JMAX, KMAX, IMAX1, &
            JMAX1, KMAX1, IMAX2, JMAX2, KMAX2, IJMAX2, IJKMAX2, MMAX, DIM_IC, &
            DIM_BC, DIM_C, DIM_IS, DT, XMIN, XLENGTH, YLENGTH, ZLENGTH 
      ELSE 
         READ (UNIT_RES, REC=4) IMIN1, JMIN1, KMIN1, IMAX, JMAX, KMAX, IMAX1, &
            JMAX1, KMAX1, IMAX2, JMAX2, KMAX2, IJMAX2, IJKMAX2, MMAX, DIM_IC, &
            DIM_BC, DIM_C, DIM_IS, DT, XMIN, XLENGTH, YLENGTH, ZLENGTH, C_E, &
            C_F, PHI, PHI_W 
      ENDIF 
!
!
! CHECK DIMENSIONS
!
      IF (IMAX2 <= DIMENSION_I) THEN 
         IF (JMAX2 <= DIMENSION_J) THEN 
            IF (KMAX2 <= DIMENSION_K) THEN 
               IF (IJKMAX2 <= DIMENSION_3) THEN 
                  IF (MMAX <= DIMENSION_M) THEN 
                     IF (DIM_IC <= DIMENSION_IC) THEN 
                        IF (DIM_BC <= DIMENSION_BC) THEN 
                           IF (DIM_C <= DIMENSION_C) THEN 
                              IF (DIM_IS <= DIMENSION_IS) THEN 
                                 M = 0 
                                 IF (MMAX + 1 > 0) THEN 
                                    NMAX(:MMAX) = 1 
                                    M = MMAX + 1 
                                 ENDIF 
                                 NEXT_RECA = 5 
!
                                 IF (VERSION_NUMBER >= 1.04) THEN 
                                    CALL IN_BIN_512 (UNIT_RES, C, DIM_C, &
                                       NEXT_RECA) 
!                                                ! work around for -O3 compiler bug
                                    NEXT_RECA = 1 + NEXT_RECA 
                                    NEXT_RECA = NEXT_RECA - 1 
                                    DO LC = 1, DIM_C 
                                    READ (UNIT_RES, REC=NEXT_RECA) C_NAME(LC) 
                                    NEXT_RECA = NEXT_RECA + 1 
                                    END DO 
                                    IF (VERSION_NUMBER < 1.12) THEN 
                                    CALL IN_BIN_512I (UNIT_RES, NMAX, MMAX + 1&
                                       , NEXT_RECA) 
                                    ELSE 
                                    READ (UNIT_RES, REC=NEXT_RECA) (NMAX(L),L=0&
                                       ,MMAX) 
                                    NEXT_RECA = NEXT_RECA + 1 
                                    ENDIF 
                                 ENDIF 
!
                                 IF (NMAX(0) <= DIMENSION_N_G) THEN 
                                    DO M = 1, MMAX 
                                    IF (NMAX(M) > DIMENSION_N_S) GO TO 900 
                                    END DO 
                                    CALL IN_BIN_512 (UNIT_RES, DX, IMAX2, &
                                       NEXT_RECA) 
                                    CALL IN_BIN_512 (UNIT_RES, DY, JMAX2, &
                                       NEXT_RECA) 
                                    CALL IN_BIN_512 (UNIT_RES, DZ, KMAX2, &
                                       NEXT_RECA) 
!
                                    READ (UNIT_RES, REC=NEXT_RECA) RUN_NAME, &
                                       DESCRIPTION, UNITS, RUN_TYPE, &
                                       COORDINATES 
                                    NEXT_RECA = NEXT_RECA + 1 
                                    IF (VERSION=='RES = 01.00' .OR. VERSION==&
                                       'RES = 01.01') THEN 
                                    READ (UNIT_RES, REC=NEXT_RECA) (D_P(L),L=1,&
                                       MMAX), (RO_S(L),L=1,MMAX), EP_STAR, &
                                       MU_G0, MW_AVG 
                                    ELSE IF (VERSION == 'RES = 01.02') THEN 
                                    READ (UNIT_RES, REC=NEXT_RECA) (D_P(L),L=1,&
                                       MMAX), (RO_S(L),L=1,MMAX), EP_STAR, &
                                       RO_G0, MU_G0, MW_AVG 
                                    ELSE IF (VERSION == 'RES = 01.03') THEN 
                                    READ (UNIT_RES, REC=NEXT_RECA) (D_P(L),L=1,&
                                       MMAX), (RO_S(L),L=1,MMAX), EP_STAR, &
                                       RO_G0, MU_G0, MW_AVG 
                                    ELSE IF (VERSION_NUMBER >= 1.04) THEN 
                                    READ (UNIT_RES, REC=NEXT_RECA) (D_P(L),L=1,&
                                       MMAX), (RO_S(L),L=1,MMAX), EP_STAR, &
                                       RO_G0, MU_G0, MW_AVG 
                                    ENDIF 
                                    NEXT_RECA = NEXT_RECA + 1 
!
                                    IF (VERSION_NUMBER >= 1.04) THEN 
                                    CALL IN_BIN_512 (UNIT_RES, MW_G, NMAX(0), &
                                       NEXT_RECA) 
                                    DO LC = 1, MMAX 
                                    READ (UNIT_RES, REC=NEXT_RECA) (MW_S(LC,N),&
                                       N=1,NMAX(LC)) 
                                    NEXT_RECA = NEXT_RECA + 1 
                                    END DO 
                                    ENDIF 
!
                                    CALL IN_BIN_512 (UNIT_RES, IC_X_W, DIM_IC, &
                                       NEXT_RECA) 
                                    CALL IN_BIN_512 (UNIT_RES, IC_X_E, DIM_IC, &
                                       NEXT_RECA) 
                                    CALL IN_BIN_512 (UNIT_RES, IC_Y_S, DIM_IC, &
                                       NEXT_RECA) 
                                    CALL IN_BIN_512 (UNIT_RES, IC_Y_N, DIM_IC, &
                                       NEXT_RECA) 
                                    CALL IN_BIN_512 (UNIT_RES, IC_Z_B, DIM_IC, &
                                       NEXT_RECA) 
                                    CALL IN_BIN_512 (UNIT_RES, IC_Z_T, DIM_IC, &
                                       NEXT_RECA) 
                                    CALL IN_BIN_512I (UNIT_RES, IC_I_W, DIM_IC&
                                       , NEXT_RECA) 
                                    CALL IN_BIN_512I (UNIT_RES, IC_I_E, DIM_IC&
                                       , NEXT_RECA) 
                                    CALL IN_BIN_512I (UNIT_RES, IC_J_S, DIM_IC&
                                       , NEXT_RECA) 
                                    CALL IN_BIN_512I (UNIT_RES, IC_J_N, DIM_IC&
                                       , NEXT_RECA) 
                                    CALL IN_BIN_512I (UNIT_RES, IC_K_B, DIM_IC&
                                       , NEXT_RECA) 
                                    CALL IN_BIN_512I (UNIT_RES, IC_K_T, DIM_IC&
                                       , NEXT_RECA) 
                                    CALL IN_BIN_512 (UNIT_RES, IC_EP_G, DIM_IC&
                                       , NEXT_RECA) 
                                    CALL IN_BIN_512 (UNIT_RES, IC_P_G, DIM_IC, &
                                       NEXT_RECA) 
                                    CALL IN_BIN_512 (UNIT_RES, IC_T_G, DIM_IC, &
                                       NEXT_RECA) 
                                    IF (VERSION_NUMBER < 1.15) THEN 
                                    CALL IN_BIN_512 (UNIT_RES, IC_T_S(1,1), &
                                       DIM_IC, NEXT_RECA) 
                                    IF (MMAX >= 2) THEN 
                                    CALL IN_BIN_512 (UNIT_RES, IC_T_S(1,2), &
                                       DIM_IC, NEXT_RECA) 
                                    ELSE 
                                    CALL IN_BIN_512 (UNIT_RES, IC_TMP, DIM_IC, &
                                       NEXT_RECA) 
                                    ENDIF 
                                    ENDIF 
!
                                    IF (VERSION_NUMBER >= 1.04) THEN 
                                    DO N = 1, NMAX(0) 
                                    CALL IN_BIN_512 (UNIT_RES, IC_X_G(1,N), &
                                       DIM_IC, NEXT_RECA) 
                                    END DO 
                                    ENDIF 
!
                                    CALL IN_BIN_512 (UNIT_RES, IC_U_G, DIM_IC, &
                                       NEXT_RECA) 
                                    CALL IN_BIN_512 (UNIT_RES, IC_V_G, DIM_IC, &
                                       NEXT_RECA) 
                                    CALL IN_BIN_512 (UNIT_RES, IC_W_G, DIM_IC, &
                                       NEXT_RECA) 
                                    DO LC = 1, MMAX 
                                    CALL IN_BIN_512 (UNIT_RES, IC_ROP_S(1,LC), &
                                       DIM_IC, NEXT_RECA) 
                                    CALL IN_BIN_512 (UNIT_RES, IC_U_S(1,LC), &
                                       DIM_IC, NEXT_RECA) 
                                    CALL IN_BIN_512 (UNIT_RES, IC_V_S(1,LC), &
                                       DIM_IC, NEXT_RECA) 
                                    CALL IN_BIN_512 (UNIT_RES, IC_W_S(1,LC), &
                                       DIM_IC, NEXT_RECA) 
                                    IF (VERSION_NUMBER >= 1.15) CALL IN_BIN_512&
                                        (UNIT_RES, IC_T_S(1,LC), DIM_IC, &
                                       NEXT_RECA) 
!
                                    IF (VERSION_NUMBER >= 1.04) THEN 
                                    DO N = 1, NMAX(LC) 
                                    CALL IN_BIN_512 (UNIT_RES, IC_X_S(1,LC,N), &
                                       DIM_IC, NEXT_RECA) 
                                    END DO 
                                    ENDIF 
                                    END DO 
                                    CALL IN_BIN_512 (UNIT_RES, BC_X_W, DIM_BC, &
                                       NEXT_RECA) 
                                    CALL IN_BIN_512 (UNIT_RES, BC_X_E, DIM_BC, &
                                       NEXT_RECA) 
                                    CALL IN_BIN_512 (UNIT_RES, BC_Y_S, DIM_BC, &
                                       NEXT_RECA) 
                                    CALL IN_BIN_512 (UNIT_RES, BC_Y_N, DIM_BC, &
                                       NEXT_RECA) 
                                    CALL IN_BIN_512 (UNIT_RES, BC_Z_B, DIM_BC, &
                                       NEXT_RECA) 
                                    CALL IN_BIN_512 (UNIT_RES, BC_Z_T, DIM_BC, &
                                       NEXT_RECA) 
                                    CALL IN_BIN_512I (UNIT_RES, BC_I_W, DIM_BC&
                                       , NEXT_RECA) 
                                    CALL IN_BIN_512I (UNIT_RES, BC_I_E, DIM_BC&
                                       , NEXT_RECA) 
                                    CALL IN_BIN_512I (UNIT_RES, BC_J_S, DIM_BC&
                                       , NEXT_RECA) 
                                    CALL IN_BIN_512I (UNIT_RES, BC_J_N, DIM_BC&
                                       , NEXT_RECA) 
                                    CALL IN_BIN_512I (UNIT_RES, BC_K_B, DIM_BC&
                                       , NEXT_RECA) 
                                    CALL IN_BIN_512I (UNIT_RES, BC_K_T, DIM_BC&
                                       , NEXT_RECA) 
                                    CALL IN_BIN_512 (UNIT_RES, BC_EP_G, DIM_BC&
                                       , NEXT_RECA) 
                                    CALL IN_BIN_512 (UNIT_RES, BC_P_G, DIM_BC, &
                                       NEXT_RECA) 
                                    CALL IN_BIN_512 (UNIT_RES, BC_T_G, DIM_BC, &
                                       NEXT_RECA) 
                                    IF (VERSION_NUMBER < 1.15) THEN 
                                    CALL IN_BIN_512 (UNIT_RES, BC_T_S(1,1), &
                                       DIM_BC, NEXT_RECA) 
                                    IF (MMAX >= 2) THEN 
                                    CALL IN_BIN_512 (UNIT_RES, BC_T_S(1,2), &
                                       DIM_BC, NEXT_RECA) 
                                    ELSE 
                                    CALL IN_BIN_512 (UNIT_RES, BC_TMP, DIM_BC, &
                                       NEXT_RECA) 
                                    ENDIF 
                                    ENDIF 
!
                                    IF (VERSION_NUMBER >= 1.04) THEN 
                                    DO N = 1, NMAX(0) 
                                    CALL IN_BIN_512 (UNIT_RES, BC_X_G(1,N), &
                                       DIM_BC, NEXT_RECA) 
                                    END DO 
                                    ENDIF 
!
                                    CALL IN_BIN_512 (UNIT_RES, BC_U_G, DIM_BC, &
                                       NEXT_RECA) 
                                    CALL IN_BIN_512 (UNIT_RES, BC_V_G, DIM_BC, &
                                       NEXT_RECA) 
                                    CALL IN_BIN_512 (UNIT_RES, BC_W_G, DIM_BC, &
                                       NEXT_RECA) 
                                    CALL IN_BIN_512 (UNIT_RES, BC_RO_G, DIM_BC&
                                       , NEXT_RECA) 
                                    CALL IN_BIN_512 (UNIT_RES, BC_ROP_G, DIM_BC&
                                       , NEXT_RECA) 
                                    CALL IN_BIN_512 (UNIT_RES, BC_VOLFLOW_G, &
                                       DIM_BC, NEXT_RECA) 
                                    CALL IN_BIN_512 (UNIT_RES, BC_MASSFLOW_G, &
                                       DIM_BC, NEXT_RECA) 
                                    DO LC = 1, MMAX 
                                    CALL IN_BIN_512 (UNIT_RES, BC_ROP_S(1,LC), &
                                       DIM_BC, NEXT_RECA) 
                                    CALL IN_BIN_512 (UNIT_RES, BC_U_S(1,LC), &
                                       DIM_BC, NEXT_RECA) 
                                    CALL IN_BIN_512 (UNIT_RES, BC_V_S(1,LC), &
                                       DIM_BC, NEXT_RECA) 
!
! Note : previous versions did not write out BC_W_s
!
                                    IF (VERSION_NUMBER >= 1.04) THEN 
                                    CALL IN_BIN_512 (UNIT_RES, BC_W_S(1,LC), &
                                       DIM_BC, NEXT_RECA) 
                                    IF (VERSION_NUMBER >= 1.15) CALL IN_BIN_512&
                                        (UNIT_RES, BC_T_S(1,LC), DIM_BC, &
                                       NEXT_RECA) 
!
                                    DO N = 1, NMAX(LC) 
                                    CALL IN_BIN_512 (UNIT_RES, BC_X_S(1,LC,N), &
                                       DIM_BC, NEXT_RECA) 
                                    END DO 
                                    ENDIF 
!
                                    CALL IN_BIN_512 (UNIT_RES, BC_VOLFLOW_S(1,&
                                       LC), DIM_BC, NEXT_RECA) 
                                    CALL IN_BIN_512 (UNIT_RES, BC_MASSFLOW_S(1,&
                                       LC), DIM_BC, NEXT_RECA) 
                                    END DO 
                                    IF (VERSION == 'RES = 01.00') THEN 
                                    L = 10 
                                    ELSE 
                                    L = DIM_BC 
                                    ENDIF 
                                    DO LC = 1, L 
                                    READ (UNIT_RES, REC=NEXT_RECA) BC_TYPE(LC) 
                                    NEXT_RECA = NEXT_RECA + 1 
                                    END DO 
                                    CALL IN_BIN_512I (UNIT_RES, array1i, IJKMAX2, &
                                       NEXT_RECA) 
                                    call convert_from_io_i(array1i,flag,ijkmax2)
!
                                    IF (VERSION_NUMBER >= 1.04) THEN 
                                    CALL IN_BIN_512 (UNIT_RES, IS_X_W, DIM_IS, &
                                       NEXT_RECA) 
                                    CALL IN_BIN_512 (UNIT_RES, IS_X_E, DIM_IS, &
                                       NEXT_RECA) 
                                    CALL IN_BIN_512 (UNIT_RES, IS_Y_S, DIM_IS, &
                                       NEXT_RECA) 
                                    CALL IN_BIN_512 (UNIT_RES, IS_Y_N, DIM_IS, &
                                       NEXT_RECA) 
                                    CALL IN_BIN_512 (UNIT_RES, IS_Z_B, DIM_IS, &
                                       NEXT_RECA) 
                                    CALL IN_BIN_512 (UNIT_RES, IS_Z_T, DIM_IS, &
                                       NEXT_RECA) 
                                    CALL IN_BIN_512I (UNIT_RES, IS_I_W, DIM_IS&
                                       , NEXT_RECA) 
                                    CALL IN_BIN_512I (UNIT_RES, IS_I_E, DIM_IS&
                                       , NEXT_RECA) 
                                    CALL IN_BIN_512I (UNIT_RES, IS_J_S, DIM_IS&
                                       , NEXT_RECA) 
                                    CALL IN_BIN_512I (UNIT_RES, IS_J_N, DIM_IS&
                                       , NEXT_RECA) 
                                    CALL IN_BIN_512I (UNIT_RES, IS_K_B, DIM_IS&
                                       , NEXT_RECA) 
                                    CALL IN_BIN_512I (UNIT_RES, IS_K_T, DIM_IS&
                                       , NEXT_RECA) 
                                    CALL IN_BIN_512 (UNIT_RES, IS_PC(1,1), &
                                       DIM_IS, NEXT_RECA) 
                                    CALL IN_BIN_512 (UNIT_RES, IS_PC(1,2), &
                                       DIM_IS, NEXT_RECA) 
                                    IF (VERSION_NUMBER >= 1.07) THEN 
                                    DO LC = 1, MMAX 
                                    CALL IN_BIN_512 (UNIT_RES, IS_VEL_S(1,LC), &
                                       DIM_IS, NEXT_RECA) 
                                    END DO 
                                    ENDIF 
                                    DO LC = 1, DIM_IS 
                                    READ (UNIT_RES, REC=NEXT_RECA) IS_TYPE(LC) 
                                    NEXT_RECA = NEXT_RECA + 1 
                                    END DO 
                                    ENDIF 
!
!  Additions from new versions of .RES file
!
                                    IF (VERSION_NUMBER >= 1.08) THEN 
                                    READ (UNIT_RES, REC=NEXT_RECA) CYCLIC_X, &
                                       CYCLIC_Y, CYCLIC_Z, CYCLIC_X_PD, &
                                       CYCLIC_Y_PD, CYCLIC_Z_PD, DELP_X, DELP_Y&
                                       , DELP_Z, U_G0, U_S0, V_G0, V_S0, W_G0, &
                                       W_S0 
                                    NEXT_RECA = NEXT_RECA + 1 
                                    ENDIF 
!
                                    IF (VERSION_NUMBER >= 1.09) THEN 
                                    READ (UNIT_RES, REC=NEXT_RECA) TIME, TSTOP&
                                       , ENERGY_EQ, RES_DT, OUT_DT, NLOG, &
                                       L_SCALE0, NO_I, NO_J, NO_K, CALL_USR 
                                    NEXT_RECA = NEXT_RECA + 1 
                                    DO LC = 1, N_SPX 
                                    READ (UNIT_RES, REC=NEXT_RECA) SPX_DT(LC) 
                                    NEXT_RECA = NEXT_RECA + 1 
                                    END DO 
                                    DO LC = 0, MMAX 
                                    READ (UNIT_RES, REC=NEXT_RECA) SPECIES_EQ(&
                                       LC) 
                                    NEXT_RECA = NEXT_RECA + 1 
                                    END DO 
                                    CALL IN_BIN_512 (UNIT_RES, USR_DT, &
                                       DIMENSION_USR, NEXT_RECA) 
                                    CALL IN_BIN_512 (UNIT_RES, USR_X_W, &
                                       DIMENSION_USR, NEXT_RECA) 
                                    CALL IN_BIN_512 (UNIT_RES, USR_X_E, &
                                       DIMENSION_USR, NEXT_RECA) 
                                    CALL IN_BIN_512 (UNIT_RES, USR_Y_S, &
                                       DIMENSION_USR, NEXT_RECA) 
                                    CALL IN_BIN_512 (UNIT_RES, USR_Y_N, &
                                       DIMENSION_USR, NEXT_RECA) 
                                    CALL IN_BIN_512 (UNIT_RES, USR_Z_B, &
                                       DIMENSION_USR, NEXT_RECA) 
                                    CALL IN_BIN_512 (UNIT_RES, USR_Z_T, &
                                       DIMENSION_USR, NEXT_RECA) 
                                    DO LC = 1, DIMENSION_USR 
                                    READ (UNIT_RES, REC=NEXT_RECA) USR_FORMAT(&
                                       LC), USR_EXT(LC), USR_TYPE(LC), USR_VAR(&
                                       LC) 
                                    NEXT_RECA = NEXT_RECA + 1 
                                    END DO 
                                    CALL IN_BIN_512 (UNIT_RES, IC_P_STAR, &
                                       DIM_IC, NEXT_RECA) 
                                    CALL IN_BIN_512 (UNIT_RES, IC_L_SCALE, &
                                       DIM_IC, NEXT_RECA) 
                                    DO LC = 1, DIM_IC 
                                    READ (UNIT_RES, REC=NEXT_RECA) IC_TYPE(LC) 
                                    NEXT_RECA = NEXT_RECA + 1 
                                    END DO 
                                    CALL IN_BIN_512 (UNIT_RES, BC_DT_0, DIM_BC&
                                       , NEXT_RECA) 
                                    CALL IN_BIN_512 (UNIT_RES, BC_JET_G0, &
                                       DIM_BC, NEXT_RECA) 
                                    CALL IN_BIN_512 (UNIT_RES, BC_DT_H, DIM_BC&
                                       , NEXT_RECA) 
                                    CALL IN_BIN_512 (UNIT_RES, BC_JET_GH, &
                                       DIM_BC, NEXT_RECA) 
                                    CALL IN_BIN_512 (UNIT_RES, BC_DT_L, DIM_BC&
                                       , NEXT_RECA) 
                                    CALL IN_BIN_512 (UNIT_RES, BC_JET_GL, &
                                       DIM_BC, NEXT_RECA) 
                                    ENDIF 
!
                                    IF (VERSION_NUMBER >= 1.10) THEN 
                                    READ (UNIT_RES, REC=NEXT_RECA) MU_GMAX 
                                    NEXT_RECA = NEXT_RECA + 1 
                                    ENDIF 
!
                                    IF (VERSION_NUMBER >= 1.11) THEN 
                                    READ (UNIT_RES, REC=NEXT_RECA) V_EX, &
                                       MODEL_B 
                                    NEXT_RECA = NEXT_RECA + 1 
                                    ENDIF 
!
                                    IF (VERSION_NUMBER >= 1.12) THEN 
                                    READ (UNIT_RES, REC=NEXT_RECA) P_REF, &
                                       P_SCALE, UR_FAC, TOL_RESID, DT_MAX, &
                                       DT_MIN, DT_FAC, CLOSE_PACKED, GRAVITY, &
                                       MU_S0 
                                    NEXT_RECA = NEXT_RECA + 1 
                                    READ (UNIT_RES, REC=NEXT_RECA) LEQ_IT, &
                                       LEQ_METHOD 
                                    NEXT_RECA = NEXT_RECA + 1 
                                    CALL IN_BIN_512 (UNIT_RES, BC_HW_G, DIM_BC&
                                       , NEXT_RECA) 
                                    CALL IN_BIN_512 (UNIT_RES, BC_UW_G, DIM_BC&
                                       , NEXT_RECA) 
                                    CALL IN_BIN_512 (UNIT_RES, BC_VW_G, DIM_BC&
                                       , NEXT_RECA) 
                                    CALL IN_BIN_512 (UNIT_RES, BC_WW_G, DIM_BC&
                                       , NEXT_RECA) 
                                    DO LC = 1, MMAX 
                                    CALL IN_BIN_512 (UNIT_RES, BC_HW_S(1,LC), &
                                       DIM_BC, NEXT_RECA) 
                                    CALL IN_BIN_512 (UNIT_RES, BC_UW_S(1,LC), &
                                       DIM_BC, NEXT_RECA) 
                                    CALL IN_BIN_512 (UNIT_RES, BC_VW_S(1,LC), &
                                       DIM_BC, NEXT_RECA) 
                                    CALL IN_BIN_512 (UNIT_RES, BC_WW_S(1,LC), &
                                       DIM_BC, NEXT_RECA) 
                                    END DO 
                                    ENDIF 
!
                                    LC = 0 
                                    IF (MMAX + 1 > 0) THEN 
                                    MOMENTUM_X_EQ(:MMAX) = .TRUE. 
                                    MOMENTUM_Y_EQ(:MMAX) = .TRUE. 
                                    MOMENTUM_Z_EQ(:MMAX) = .TRUE. 
                                    LC = MMAX + 1 
                                    ENDIF 
                                    TOL_DIVERGE = 1.E+4 
                                    IF (VERSION_NUMBER >= 1.13) THEN 
                                    READ (UNIT_RES, REC=NEXT_RECA) &
                                       MOMENTUM_X_EQ, MOMENTUM_Y_EQ, &
                                       MOMENTUM_Z_EQ, TOL_DIVERGE, DISCRETIZE, &
                                       FULL_LOG 
                                    NEXT_RECA = NEXT_RECA + 1 
                                    ENDIF 
!
                                    IF (VERSION_NUMBER >= 1.14) THEN 
                                    READ (UNIT_RES, REC=NEXT_RECA) DETECT_STALL 
                                    NEXT_RECA = NEXT_RECA + 1 
                                    ENDIF 
!
                                    IF (VERSION_NUMBER >= 1.15) THEN 
                                    READ (UNIT_RES, REC=NEXT_RECA) K_G0, K_S0, &
                                       C_PG0, C_PS0, TOL_RESID_T, TOL_RESID_X 
                                    NEXT_RECA = NEXT_RECA + 1 
                                    CALL IN_BIN_512 (UNIT_RES, IC_GAMA_RG, &
                                       DIM_IC, NEXT_RECA) 
                                    CALL IN_BIN_512 (UNIT_RES, IC_T_RG, DIM_IC&
                                       , NEXT_RECA) 
                                    DO LC = 1, MMAX 
                                    CALL IN_BIN_512 (UNIT_RES, IC_GAMA_RS(1,LC)&
                                       , DIM_IC, NEXT_RECA) 
                                    CALL IN_BIN_512 (UNIT_RES, IC_T_RS(1,LC), &
                                       DIM_IC, NEXT_RECA) 
                                    END DO 
                                    ENDIF 
!
                                    IF (VERSION_NUMBER >= 1.2) THEN 
                                    READ (UNIT_RES, REC=NEXT_RECA) NORM_G, &
                                       NORM_S 
                                    NEXT_RECA = NEXT_RECA + 1 
                                    ENDIF 
!
!  Add new read statements above this line.  Remember to update NEXT_RECA.
!  Remember to update the version number check near begining of this subroutine.
!------------------------------------------------------------------------------
!
                                    READ (UNIT_RES, REC=3) NEXT_RECA 
!
!  Since the value of UNDEFINED was changed ...
!
                                    IF (RO_G0 >= 1E30) RO_G0 = UNDEFINED 
                                    IF (MU_G0 >= 1E30) MU_G0 = UNDEFINED 
                                    IF (MW_AVG >= 1E30) MW_AVG = UNDEFINED 
                                    IF (C_E >= 1E30) C_E = UNDEFINED 
!
                                    call unlock_tmp_array
                                    RETURN  
!
! HERE IF DIMENSION ERROR
!
                                 ENDIF 
                              ENDIF 
                           ENDIF 
                        ENDIF 
                     ENDIF 
                  ENDIF 
               ENDIF 
            ENDIF 
         ENDIF 
      ENDIF 
  900 CONTINUE 
      WRITE (*, *) ' ' 
      WRITE (*, *) ' **************************************' 
      WRITE (*, *) ' From: READ_RES0' 
      WRITE (*, *) ' DIMENSION ERROR ---' 
      WRITE (*, *) ' ' 
      WRITE (*, *) ' DIMENSION_I  = ', DIMENSION_I, ' IMAX2        = ', IMAX2 
      WRITE (*, *) ' DIMENSION_J  = ', DIMENSION_J, ' JMAX2        = ', JMAX2 
      WRITE (*, *) ' DIMENSION_K  = ', DIMENSION_K, ' KMAX2        = ', KMAX2 
      WRITE (*, *) ' DIMENSION_3  = ', DIMENSION_3, ' IJKMAX2      = ', IJKMAX2 
      WRITE (*, *) ' DIMENSION_M  = ', DIMENSION_M, ' MMAX         = ', MMAX 
      WRITE (*, *) ' DIMENSION_IC = ', DIMENSION_IC, ' DIM_IC       = ', DIM_IC 
      WRITE (*, *) ' DIMENSION_BC = ', DIMENSION_BC, ' DIM_BC       = ', DIM_BC 
      WRITE (*, *) ' DIMENSION_IS = ', DIMENSION_IS, ' DIM_IS       = ', DIM_IS 
      WRITE (*, *) ' DIMENSION_C  = ', DIMENSION_C, ' DIM_C        = ', DIM_C 
      WRITE (*, *) ' DIMENSION_N_g= ', DIMENSION_N_G, ' NMAX(0)      = ', NMAX(&
         0) 
      WRITE (*, *) ' DIMENSION_N_s= ', DIMENSION_N_S 
      DO M = 1, MMAX 
         WRITE (*, '(A, I2, A, I4)') ' NMAX(', M, ') = ', NMAX(M) 
      END DO 
      WRITE (*, *) ' ' 
!
      call unlock_tmp_array
      STOP  
      END SUBROUTINE READ_RES0 
