!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_RES1                                             C
!  Purpose: write out the time-dependent restart records               C
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
!  Variables referenced: TIME, NSTEP, EP_g, P_g, P_star, RO_g, ROP_g   C
!                        T_g, T_s, U_g, V_g, W_g, ROP_s, U_s    C
!                        V_s, W_s, IJKMAX2                             C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables: LC, N, NEXT_REC                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE WRITE_RES1 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE fldvar
      USE geometry
      USE physprop
      USE run
      USE funits 
      USE output
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
!
!
!
!             loop counter
      INTEGER :: LC, N
!
!             pointer to first time-dependent record in restart file
      INTEGER :: NEXT_REC 
!-----------------------------------------------
!
      call lock_tmp_array

      READ (UNIT_RES, REC=3) NEXT_REC 
      WRITE (UNIT_RES, REC=NEXT_REC) TIME, DT, NSTEP 
      NEXT_REC = NEXT_REC + 1 
!
      call convert_to_io_dp(EP_g,array1,ijkmax2)
      CALL OUT_BIN_512 (UNIT_RES, array1, IJKMAX2, NEXT_REC) 
      call convert_to_io_dp(P_g,array1,ijkmax2)
      CALL OUT_BIN_512 (UNIT_RES, array1, IJKMAX2, NEXT_REC) 
      call convert_to_io_dp(P_STAR,array1,ijkmax2)
      CALL OUT_BIN_512 (UNIT_RES, array1, IJKMAX2, NEXT_REC) 
      call convert_to_io_dp(RO_g,array1,ijkmax2)
      CALL OUT_BIN_512 (UNIT_RES, array1, IJKMAX2, NEXT_REC) 
      call convert_to_io_dp(ROP_g,array1,ijkmax2)
      CALL OUT_BIN_512 (UNIT_RES, array1, IJKMAX2, NEXT_REC) 
      call convert_to_io_dp(T_g,array1,ijkmax2)
      CALL OUT_BIN_512 (UNIT_RES, array1, IJKMAX2, NEXT_REC) 
      DO N = 1, NMAX(0) 
         call convert_to_io_dp(X_G(1,N),array1,ijkmax2)
         CALL OUT_BIN_512 (UNIT_RES, array1, IJKMAX2, NEXT_REC) 
      END DO 
      call convert_to_io_dp(U_g,array1,ijkmax2)
      CALL OUT_BIN_512 (UNIT_RES, array1, IJKMAX2, NEXT_REC) 
      call convert_to_io_dp(V_g,array1,ijkmax2)
      CALL OUT_BIN_512 (UNIT_RES, array1, IJKMAX2, NEXT_REC) 
      call convert_to_io_dp(W_g,array1,ijkmax2)
      CALL OUT_BIN_512 (UNIT_RES, array1, IJKMAX2, NEXT_REC) 
!
      DO LC = 1, MMAX 
         call convert_to_io_dp(ROP_S(1,LC),array1,ijkmax2)
         CALL OUT_BIN_512 (UNIT_RES,array1 , IJKMAX2, NEXT_REC) 
         call convert_to_io_dp(T_S(1,LC),array1,ijkmax2)
         CALL OUT_BIN_512 (UNIT_RES,array1 , IJKMAX2, NEXT_REC) 
         call convert_to_io_dp(U_S(1,LC),array1,ijkmax2)
         CALL OUT_BIN_512 (UNIT_RES,array1, IJKMAX2, NEXT_REC) 
         call convert_to_io_dp(V_S(1,LC),array1,ijkmax2)
         CALL OUT_BIN_512 (UNIT_RES,array1 , IJKMAX2, NEXT_REC) 
         call convert_to_io_dp(W_S(1,LC),array1,ijkmax2)
         CALL OUT_BIN_512 (UNIT_RES,array1, IJKMAX2, NEXT_REC) 
         call convert_to_io_dp(THETA_M(1,LC),array1,ijkmax2)
         CALL OUT_BIN_512 (UNIT_RES,array1 , IJKMAX2, NEXT_REC) 
         DO N = 1, NMAX(LC) 
            call convert_to_io_dp(X_S(1,LC,N),array1,ijkmax2)
            CALL OUT_BIN_512 (UNIT_RES,array1, IJKMAX2, NEXT_REC) 
         END DO 
      END DO 
      CALL FLUSH (UNIT_RES) 
!
      call unlock_tmp_array
!
      RETURN  
      END SUBROUTINE WRITE_RES1 
