!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_RES0                                             C
!  Purpose: write out the initial restart records (namelist data)      C
!                                                                      C
!  Author: P. Nicoletti                               Date: 13-DEC-91  C
!  Reviewer: P. Nicoletti, W. Rogers, M. Syamlal      Date: 24-JAN-92  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: RUN_NAME  ,  ID_MONTH  ,  ID_DAY , ID_YEAR    C
!                        ID_HOUR, ID_MINUTE, ID_SECOND, IMAX, JMAX     C
!                        KMAX, IMAX1, JMAX1, KMAX1, IMAX2, JMAX2,KMAX2 C
!                        IJMAX2, IJKMAX2, MMAX, DT, XLENGTH, YLENGTH   C
!                        ZLENGTH, DX, DY, DZ, RUN_NAME, DESCRIPTION    C
!                        UNITS, RUN_TYPE, CORDINATES, D_p0, RO_s,      C
!                        EP_star, MU_g0, MW_AVG, IC_X_w, IC_X_e, IC_Y_sC
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
!                        BC_MASSFLOW_s, BC_TYPE, FLAG, RO_g0           C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables: LC, L, N, NEXT_RECA, VERSION                       C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE WRITE_RES0 
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
      USE is
      USE bc
      USE constant
      USE funits 
      USE output
      USE scales 
      USE scalars
      USE rxns
      USE ur_facs 
      USE leqsol 
      USE toleranc 
      USE cdist
      USE compar           !//
      USE mpi_utility      !// for gather
      USE sendrecv         !// for filling the boundary information
      USE stiff_chem
!//       USE tmp_array    !// no longer using these arrays
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
!//  ... temporary arrays
!
      integer, allocatable :: arr1(:)
      integer, allocatable :: arr2(:)

      integer :: work_around(100)
!
!                loop counters
      INTEGER :: LC, L, N, M
!
!                Pointer to the next record
      INTEGER :: NEXT_RECA 
!
!                file version id
      CHARACTER :: VERSION*512 
!-----------------------------------------------
!
!
      NEXT_RECA = 5 
      
      ! note: the value below for version must be the same as
      !       the value of version in mfix.f
      !       If you change it in this subroutine, you must
      !       change it in mfix.f also.
      !   
      !       The value should also be consistent with the check
      !       in read_res0
      
      VERSION = 'RES = 01.8' 
!
!     Add new data entries at the end of the file and identify version no.

      if (myPE.eq.PE_IO) then
         allocate (arr1(ijkmax3))        !// 
         allocate (arr2(ijkmax2))        !// 
      else
         allocate (arr1(1))              !// 
         allocate (arr2(1))              !// 
         goto 1100
      end if

!------------------------------------------------------------------------
!//SP - Temporarily disabled so that the binary files can be diffed.....
!
      WRITE (UNIT_RES, REC=1) VERSION 
      WRITE (UNIT_RES, REC=2) RUN_NAME, ID_MONTH, ID_DAY, ID_YEAR, ID_HOUR, &
         ID_MINUTE, ID_SECOND 
      WRITE (UNIT_RES, REC=3) NEXT_RECA 
      WRITE (UNIT_RES, REC=4) IMIN1, JMIN1, KMIN1, IMAX, JMAX, KMAX, IMAX1, &
         JMAX1, KMAX1, IMAX2, JMAX2, KMAX2, IJMAX2, IJKMAX2, MMAX, DIMENSION_IC&
         , DIMENSION_BC, DIMENSION_C, DIMENSION_IS, DT, XMIN, XLENGTH, YLENGTH&
         , ZLENGTH, C_E, C_F, PHI, PHI_W
      CALL OUT_BIN_512 (UNIT_RES, C, DIMENSION_C, NEXT_RECA)
      NEXT_RECA = 1 + NEXT_RECA                  ! work around for -O3 compiler bug
      NEXT_RECA = NEXT_RECA - 1 
      DO LC = 1, DIMENSION_C 
         WRITE (UNIT_RES, REC=NEXT_RECA) C_NAME(LC) 
         NEXT_RECA = NEXT_RECA + 1 
      END DO
      write (*,*) 'next_reca = ' , next_reca
      do l = 0,mmax
         write (*,*) L,nmax(l)
	 work_around(L+1) = nmax(L)
      end do
 !     next_reca = next_reca + 1
 !     CALL OUT_BIN_512I (UNIT_RES, work_around, MMAX+1, NEXT_RECA)
!      WRITE (UNIT_RES, REC=NEXT_RECA) (work_around(L),L=1,MMAX+1)
      WRITE (UNIT_RES, REC=NEXT_RECA) (NMAX(L),L=0,MMAX)
 !      do l =0,mmax
 !         write (unit_res,rec=next_reca) nmax(l)
!	  next_reca = next_reca + 1
 !      end do
!      write (unit_res,rec=next_reca) nmax(0) , nmax(1)
      NEXT_RECA = NEXT_RECA + 1
      write (*,*) ' next_reca = ' , next_reca
!
      CALL OUT_BIN_512 (UNIT_RES, DX(1), IMAX2, NEXT_RECA) 
      CALL OUT_BIN_512 (UNIT_RES, DY(1), JMAX2, NEXT_RECA) 
      CALL OUT_BIN_512 (UNIT_RES, DZ(1), KMAX2, NEXT_RECA) 
      WRITE (UNIT_RES, REC=NEXT_RECA) RUN_NAME, DESCRIPTION, UNITS, RUN_TYPE, &
         COORDINATES 
      NEXT_RECA = NEXT_RECA + 1 
      WRITE (UNIT_RES, REC=NEXT_RECA) (D_P0(L),L=1,MMAX), (RO_S0(L),L=1,MMAX), &
         EP_STAR, RO_G0, MU_G0, MW_AVG 
      NEXT_RECA = NEXT_RECA + 1 
      CALL OUT_BIN_512 (UNIT_RES, MW_G, NMAX(0), NEXT_RECA) 
      DO LC = 1, MMAX 
         WRITE (UNIT_RES, REC=NEXT_RECA) (MW_S(LC,N),N=1,NMAX(LC)) 
         NEXT_RECA = NEXT_RECA + 1 
      END DO 
      CALL OUT_BIN_512 (UNIT_RES, IC_X_W, DIMENSION_IC, NEXT_RECA) 
      CALL OUT_BIN_512 (UNIT_RES, IC_X_E, DIMENSION_IC, NEXT_RECA) 
      CALL OUT_BIN_512 (UNIT_RES, IC_Y_S, DIMENSION_IC, NEXT_RECA) 
      CALL OUT_BIN_512 (UNIT_RES, IC_Y_N, DIMENSION_IC, NEXT_RECA) 
      CALL OUT_BIN_512 (UNIT_RES, IC_Z_B, DIMENSION_IC, NEXT_RECA) 
      CALL OUT_BIN_512 (UNIT_RES, IC_Z_T, DIMENSION_IC, NEXT_RECA) 
      CALL OUT_BIN_512I (UNIT_RES, IC_I_W, DIMENSION_IC, NEXT_RECA) 
      CALL OUT_BIN_512I (UNIT_RES, IC_I_E, DIMENSION_IC, NEXT_RECA) 
      CALL OUT_BIN_512I (UNIT_RES, IC_J_S, DIMENSION_IC, NEXT_RECA) 
      CALL OUT_BIN_512I (UNIT_RES, IC_J_N, DIMENSION_IC, NEXT_RECA) 
      CALL OUT_BIN_512I (UNIT_RES, IC_K_B, DIMENSION_IC, NEXT_RECA)
      CALL OUT_BIN_512I (UNIT_RES, IC_K_T, DIMENSION_IC, NEXT_RECA)
      CALL OUT_BIN_512 (UNIT_RES, IC_EP_G, DIMENSION_IC, NEXT_RECA) 
      CALL OUT_BIN_512 (UNIT_RES, IC_P_G, DIMENSION_IC, NEXT_RECA) 
      CALL OUT_BIN_512 (UNIT_RES, IC_T_G, DIMENSION_IC, NEXT_RECA) 
      DO N = 1, NMAX(0) 
         CALL OUT_BIN_512 (UNIT_RES, IC_X_G(1,N), DIMENSION_IC, NEXT_RECA) 
      END DO 
      CALL OUT_BIN_512 (UNIT_RES, IC_U_G, DIMENSION_IC, NEXT_RECA) 
      CALL OUT_BIN_512 (UNIT_RES, IC_V_G, DIMENSION_IC, NEXT_RECA) 
      CALL OUT_BIN_512 (UNIT_RES, IC_W_G, DIMENSION_IC, NEXT_RECA) 
      DO LC = 1, MMAX 
         CALL OUT_BIN_512 (UNIT_RES, IC_ROP_S(1,LC), DIMENSION_IC, NEXT_RECA) 
         CALL OUT_BIN_512 (UNIT_RES, IC_U_S(1,LC), DIMENSION_IC, NEXT_RECA) 
         CALL OUT_BIN_512 (UNIT_RES, IC_V_S(1,LC), DIMENSION_IC, NEXT_RECA) 
         CALL OUT_BIN_512 (UNIT_RES, IC_W_S(1,LC), DIMENSION_IC, NEXT_RECA) 
         CALL OUT_BIN_512 (UNIT_RES, IC_T_S(1,LC), DIMENSION_IC, NEXT_RECA) 
         DO N = 1, NMAX(LC) 
            CALL OUT_BIN_512(UNIT_RES,IC_X_S(1,LC,N),DIMENSION_IC,NEXT_RECA) 
         END DO 
      END DO 
      CALL OUT_BIN_512 (UNIT_RES, BC_X_W, DIMENSION_BC, NEXT_RECA) 
      CALL OUT_BIN_512 (UNIT_RES, BC_X_E, DIMENSION_BC, NEXT_RECA) 
      CALL OUT_BIN_512 (UNIT_RES, BC_Y_S, DIMENSION_BC, NEXT_RECA) 
      CALL OUT_BIN_512 (UNIT_RES, BC_Y_N, DIMENSION_BC, NEXT_RECA) 
      CALL OUT_BIN_512 (UNIT_RES, BC_Z_B, DIMENSION_BC, NEXT_RECA) 
      CALL OUT_BIN_512 (UNIT_RES, BC_Z_T, DIMENSION_BC, NEXT_RECA) 
      
      CALL OUT_BIN_512I (UNIT_RES, BC_I_W, DIMENSION_BC, NEXT_RECA) 
      CALL OUT_BIN_512I (UNIT_RES, BC_I_E, DIMENSION_BC, NEXT_RECA) 
      CALL OUT_BIN_512I (UNIT_RES, BC_J_S, DIMENSION_BC, NEXT_RECA) 
      CALL OUT_BIN_512I (UNIT_RES, BC_J_N, DIMENSION_BC, NEXT_RECA) 
      CALL OUT_BIN_512I (UNIT_RES, BC_K_B, DIMENSION_BC, NEXT_RECA) 
      CALL OUT_BIN_512I (UNIT_RES, BC_K_T, DIMENSION_BC, NEXT_RECA) 
      CALL OUT_BIN_512 (UNIT_RES, BC_EP_G, DIMENSION_BC, NEXT_RECA) 
      CALL OUT_BIN_512 (UNIT_RES, BC_P_G, DIMENSION_BC, NEXT_RECA) 
      CALL OUT_BIN_512 (UNIT_RES, BC_T_G, DIMENSION_BC, NEXT_RECA) 
      DO N = 1, NMAX(0) 
         CALL OUT_BIN_512 (UNIT_RES, BC_X_G(1,N), DIMENSION_BC, NEXT_RECA) 
      END DO 
      CALL OUT_BIN_512 (UNIT_RES, BC_U_G, DIMENSION_BC, NEXT_RECA) 
      CALL OUT_BIN_512 (UNIT_RES, BC_V_G, DIMENSION_BC, NEXT_RECA) 
      CALL OUT_BIN_512 (UNIT_RES, BC_W_G, DIMENSION_BC, NEXT_RECA) 
      CALL OUT_BIN_512 (UNIT_RES, BC_RO_G, DIMENSION_BC, NEXT_RECA) 
      CALL OUT_BIN_512 (UNIT_RES, BC_ROP_G, DIMENSION_BC, NEXT_RECA) 
      CALL OUT_BIN_512 (UNIT_RES, BC_VOLFLOW_G, DIMENSION_BC, NEXT_RECA) 
      CALL OUT_BIN_512 (UNIT_RES, BC_MASSFLOW_G, DIMENSION_BC, NEXT_RECA) 
      DO LC = 1, MMAX 
         CALL OUT_BIN_512 (UNIT_RES, BC_ROP_S(1,LC), DIMENSION_BC, NEXT_RECA) 
         CALL OUT_BIN_512 (UNIT_RES, BC_U_S(1,LC), DIMENSION_BC, NEXT_RECA) 
         CALL OUT_BIN_512 (UNIT_RES, BC_V_S(1,LC), DIMENSION_BC, NEXT_RECA) 
         CALL OUT_BIN_512 (UNIT_RES, BC_W_S(1,LC), DIMENSION_BC, NEXT_RECA) 
         CALL OUT_BIN_512 (UNIT_RES, BC_T_S(1,LC), DIMENSION_BC, NEXT_RECA) 
         DO N = 1, NMAX(LC) 
            CALL OUT_BIN_512(UNIT_RES,BC_X_S(1,LC,N),DIMENSION_BC,NEXT_RECA) 
         END DO 
         CALL OUT_BIN_512 (UNIT_RES, BC_VOLFLOW_S(1,LC), DIMENSION_BC, &
            NEXT_RECA) 
         CALL OUT_BIN_512 (UNIT_RES, BC_MASSFLOW_S(1,LC), DIMENSION_BC, &
            NEXT_RECA) 
      END DO 
      DO LC = 1, DIMENSION_BC 
         WRITE (UNIT_RES, REC=NEXT_RECA) BC_TYPE(LC) 
         NEXT_RECA = NEXT_RECA + 1 
      END DO 
!
 1100 continue
!     call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
!     call send_recv (flag,2)
      call gather (flag,arr1,root)                               !// 
!//SP To take care of filling the processor ghost layers with the correct values
      call scatter (flag,arr1,root)                               !// 
!     call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
      if (myPE .ne. PE_IO) goto 1200
      call convert_to_io_i(arr1,arr2,ijkmax2)                    !//
      CALL OUT_BIN_512I (UNIT_RES, arr2, IJKMAX2, NEXT_RECA)     !//
!
      CALL OUT_BIN_512 (UNIT_RES, IS_X_W, DIMENSION_IS, NEXT_RECA) 
      CALL OUT_BIN_512 (UNIT_RES, IS_X_E, DIMENSION_IS, NEXT_RECA) 
      CALL OUT_BIN_512 (UNIT_RES, IS_Y_S, DIMENSION_IS, NEXT_RECA) 
      CALL OUT_BIN_512 (UNIT_RES, IS_Y_N, DIMENSION_IS, NEXT_RECA) 
      CALL OUT_BIN_512 (UNIT_RES, IS_Z_B, DIMENSION_IS, NEXT_RECA) 
      CALL OUT_BIN_512 (UNIT_RES, IS_Z_T, DIMENSION_IS, NEXT_RECA) 
      CALL OUT_BIN_512I (UNIT_RES, IS_I_W, DIMENSION_IS, NEXT_RECA) 
      CALL OUT_BIN_512I (UNIT_RES, IS_I_E, DIMENSION_IS, NEXT_RECA) 
      CALL OUT_BIN_512I (UNIT_RES, IS_J_S, DIMENSION_IS, NEXT_RECA) 
      CALL OUT_BIN_512I (UNIT_RES, IS_J_N, DIMENSION_IS, NEXT_RECA) 
      CALL OUT_BIN_512I (UNIT_RES, IS_K_B, DIMENSION_IS, NEXT_RECA) 
      CALL OUT_BIN_512I (UNIT_RES, IS_K_T, DIMENSION_IS, NEXT_RECA) 
      CALL OUT_BIN_512 (UNIT_RES, IS_PC(1,1), DIMENSION_IS, NEXT_RECA) 
      CALL OUT_BIN_512 (UNIT_RES, IS_PC(1,2), DIMENSION_IS, NEXT_RECA) 
      DO LC = 1, MMAX 
         CALL OUT_BIN_512 (UNIT_RES, IS_VEL_S(1,LC), DIMENSION_IS, NEXT_RECA) 
      END DO 
      DO LC = 1, DIMENSION_IS 
         WRITE (UNIT_RES, REC=NEXT_RECA) IS_TYPE(LC) 
         NEXT_RECA = NEXT_RECA + 1 
      END DO 
      WRITE (UNIT_RES, REC=NEXT_RECA) CYCLIC_X, CYCLIC_Y, CYCLIC_Z, CYCLIC_X_PD&
         , CYCLIC_Y_PD, CYCLIC_Z_PD, DELP_X, DELP_Y, DELP_Z, U_G0, U_S0, V_G0, &
         V_S0, W_G0, W_S0 
      NEXT_RECA = NEXT_RECA + 1 
!
!     Version 01.09
      WRITE (UNIT_RES, REC=NEXT_RECA) TIME, TSTOP, ENERGY_EQ, RES_DT, OUT_DT, &
         NLOG, L_SCALE0, NO_I, NO_J, NO_K, CALL_USR 
      NEXT_RECA = NEXT_RECA + 1
      write (unit_res,rec=next_reca) n_spx 
      NEXT_RECA = NEXT_RECA + 1
      DO LC = 1, N_SPX 
         WRITE (UNIT_RES, REC=NEXT_RECA) SPX_DT(LC) 
         NEXT_RECA = NEXT_RECA + 1 
      END DO 
      DO LC = 0, MMAX 
         WRITE (UNIT_RES, REC=NEXT_RECA) SPECIES_EQ(LC) 
         NEXT_RECA = NEXT_RECA + 1 
      END DO 
      CALL OUT_BIN_512 (UNIT_RES, USR_DT, DIMENSION_USR, NEXT_RECA) 
      CALL OUT_BIN_512 (UNIT_RES, USR_X_W, DIMENSION_USR, NEXT_RECA) 
      CALL OUT_BIN_512 (UNIT_RES, USR_X_E, DIMENSION_USR, NEXT_RECA) 
      CALL OUT_BIN_512 (UNIT_RES, USR_Y_S, DIMENSION_USR, NEXT_RECA) 
      CALL OUT_BIN_512 (UNIT_RES, USR_Y_N, DIMENSION_USR, NEXT_RECA) 
      CALL OUT_BIN_512 (UNIT_RES, USR_Z_B, DIMENSION_USR, NEXT_RECA) 
      CALL OUT_BIN_512 (UNIT_RES, USR_Z_T, DIMENSION_USR, NEXT_RECA) 
      DO LC = 1, DIMENSION_USR 
         WRITE (UNIT_RES, REC=NEXT_RECA) USR_FORMAT(LC), USR_EXT(LC), USR_TYPE(&
            LC), USR_VAR(LC) 
         NEXT_RECA = NEXT_RECA + 1 
      END DO 
      CALL OUT_BIN_512 (UNIT_RES, IC_P_STAR, DIMENSION_IC, NEXT_RECA) 
      CALL OUT_BIN_512 (UNIT_RES, IC_L_SCALE, DIMENSION_IC, NEXT_RECA) 
      DO LC = 1, DIMENSION_IC 
         WRITE (UNIT_RES, REC=NEXT_RECA) IC_TYPE(LC) 
         NEXT_RECA = NEXT_RECA + 1 
      END DO 
      CALL OUT_BIN_512 (UNIT_RES, BC_DT_0, DIMENSION_BC, NEXT_RECA) 
      CALL OUT_BIN_512 (UNIT_RES, BC_JET_G0, DIMENSION_BC, NEXT_RECA) 
      CALL OUT_BIN_512 (UNIT_RES, BC_DT_H, DIMENSION_BC, NEXT_RECA) 
      CALL OUT_BIN_512 (UNIT_RES, BC_JET_GH, DIMENSION_BC, NEXT_RECA) 
      CALL OUT_BIN_512 (UNIT_RES, BC_DT_L, DIMENSION_BC, NEXT_RECA) 
      CALL OUT_BIN_512 (UNIT_RES, BC_JET_GL, DIMENSION_BC, NEXT_RECA) 
!
!     Version 01.10
      WRITE (UNIT_RES, REC=NEXT_RECA) MU_GMAX 
      NEXT_RECA = NEXT_RECA + 1 
!
!     Version 01.11
      WRITE (UNIT_RES, REC=NEXT_RECA) V_EX, MODEL_B 
      NEXT_RECA = NEXT_RECA + 1 
!
!     Version 01.12
      WRITE (UNIT_RES, REC=NEXT_RECA) P_REF, P_SCALE, UR_FAC, TOL_RESID, DT_MAX&
         , DT_MIN, DT_FAC, CLOSE_PACKED, GRAVITY, MU_S0 
      NEXT_RECA = NEXT_RECA + 1 
      WRITE (UNIT_RES, REC=NEXT_RECA) LEQ_IT, LEQ_METHOD 
      NEXT_RECA = NEXT_RECA + 1 
      CALL OUT_BIN_512 (UNIT_RES, BC_HW_G, DIMENSION_BC, NEXT_RECA) 
      CALL OUT_BIN_512 (UNIT_RES, BC_UW_G, DIMENSION_BC, NEXT_RECA) 
      CALL OUT_BIN_512 (UNIT_RES, BC_VW_G, DIMENSION_BC, NEXT_RECA) 
      CALL OUT_BIN_512 (UNIT_RES, BC_WW_G, DIMENSION_BC, NEXT_RECA) 
      DO LC = 1, MMAX 
         CALL OUT_BIN_512 (UNIT_RES, BC_HW_S(1,LC), DIMENSION_BC, NEXT_RECA) 
         CALL OUT_BIN_512 (UNIT_RES, BC_UW_S(1,LC), DIMENSION_BC, NEXT_RECA) 
         CALL OUT_BIN_512 (UNIT_RES, BC_VW_S(1,LC), DIMENSION_BC, NEXT_RECA) 
         CALL OUT_BIN_512 (UNIT_RES, BC_WW_S(1,LC), DIMENSION_BC, NEXT_RECA) 
      END DO 
      WRITE (UNIT_RES, REC=NEXT_RECA) MOMENTUM_X_EQ, MOMENTUM_Y_EQ, &
         MOMENTUM_Z_EQ, TOL_DIVERGE, DISCRETIZE, FULL_LOG 
      NEXT_RECA = NEXT_RECA + 1 
!
!     Version 01.14
      WRITE (UNIT_RES, REC=NEXT_RECA) DETECT_STALL 
      NEXT_RECA = NEXT_RECA + 1 
!
!     Version 01.15
      WRITE (UNIT_RES, REC=NEXT_RECA) K_G0, K_S0, C_PG0, C_PS0, TOL_RESID_T, &
         TOL_RESID_X 
      NEXT_RECA = NEXT_RECA + 1 
      CALL OUT_BIN_512 (UNIT_RES, IC_GAMA_RG, DIMENSION_IC, NEXT_RECA) 
      CALL OUT_BIN_512 (UNIT_RES, IC_T_RG, DIMENSION_IC, NEXT_RECA) 
      DO LC = 1, MMAX 
         CALL OUT_BIN_512 (UNIT_RES, IC_GAMA_RS(1,LC), DIMENSION_IC, NEXT_RECA) 
         CALL OUT_BIN_512 (UNIT_RES, IC_T_RS(1,LC), DIMENSION_IC, NEXT_RECA) 
      END DO 
!
!     Version 01.2
      WRITE (UNIT_RES, REC=NEXT_RECA) NORM_G, NORM_S 
      NEXT_RECA = NEXT_RECA + 1 
!
!     Version 01.3
      WRITE (UNIT_RES, REC=NEXT_RECA) NScalar, TOL_RESID_Scalar, DIM_SCALAR  
      NEXT_RECA = NEXT_RECA + 1 
      CALL OUT_BIN_512I (UNIT_RES, Phase4Scalar, DIM_SCALAR, NEXT_RECA)
!
!     Version 1.4 -- write radiation variables in write_res1 

!
!     Version 1.5 -- write nRR
      WRITE (UNIT_RES, REC=NEXT_RECA) nRR  
      NEXT_RECA = NEXT_RECA + 1 

!
!     Version 1.6 -- write k and epsilon in write_res1 and spx1
      WRITE (UNIT_RES, REC=NEXT_RECA) K_epsilon  
      NEXT_RECA = NEXT_RECA + 1 

!
!     Version 1.7 -- write STIFF_CHEMISTRY and CALL_ISAT
      WRITE (UNIT_RES, REC=NEXT_RECA) STIFF_CHEMISTRY, CALL_ISAT
      NEXT_RECA = NEXT_RECA + 1 

!
!     Version 1.8 -- write SOLID_RO_V and densities of each solids species
      WRITE (UNIT_RES, REC=NEXT_RECA) SOLID_RO_V
      NEXT_RECA = NEXT_RECA + 1 
      DO LC = 1, MMAX 
         WRITE (UNIT_RES, REC=NEXT_RECA) (RO_SS(LC,N),N=1,NMAX(LC))
         NEXT_RECA = NEXT_RECA + 1 
      END DO 


!
!  Add new write statements above this line.  Remember to update NEXT_RECA.
!  Remember to change the version number near begining of this subroutine.
!  Also modify READ_RES0.  The routines such as OUT_BIN_512 etc. writes
!  arrays dimensioned ARRAY(DIM).  So arrays dimensioned ARRAY(DIM1:DIM2)
!  should be passed as ARRAY(DIM1) and array length as DIM2-DIM1+1.
!---------------------------------------------------------------------------
      WRITE (UNIT_RES, REC=3) NEXT_RECA 
      CALL FLUSH (UNIT_RES) 

 1200 continue

      if (bDist_IO .and. myPE.ne.PE_IO) then
         WRITE (UNIT_RES, REC=1) VERSION 
         WRITE (UNIT_RES, REC=2) RUN_NAME, ID_MONTH, ID_DAY, ID_YEAR, ID_HOUR, &
             ID_MINUTE, ID_SECOND 
         WRITE (UNIT_RES, REC=3) 4 
         CALL FLUSH (UNIT_RES) 
      end if

      deallocate (arr1)              !//
      deallocate (arr2)              !//

      RETURN  
      END SUBROUTINE WRITE_RES0 
