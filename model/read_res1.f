!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: READ_RES1                                              C
!  Purpose: read in the time-dependent restart records                 C
!                                                                      C
!  Author: P. Nicoletti                               Date: 03-JAN-92  C
!  Reviewer: P. Nicoletti, W. Rogers, M. Syamlal      Date: 24-JAN-92  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: IJKMAX2, MMAX, DT                             C
!  Variables modified: TIME, NSTEP, EP_g, P_g, P_star, RO_g            C
!                      ROP_g, T_g, T_s,  U_g, V_g, W_g, ROP_s    C
!                      U_s, V_s, W_s                                   C
!                                                                      C
!  Local variables: TIME_READ, LC, NEXT_REC                            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE READ_RES1 
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
!             loop counter
      INTEGER LC 
! 
!                      Local species index
      INTEGER          N
!
!             pointer to the next record
      INTEGER NEXT_REC
!
!                file version id
      CHARACTER  VERSION*512
!
!                version number
      REAL       VERSION_NUMBER
!
!                      Dummy array
      DOUBLE PRECISION Tmp(DIMENSION_3)
      DOUBLE PRECISION DT_SAVE
!-----------------------------------------------
!
      call lock_tmp_array
!
!     Use DT from data file if DT_FAC is set to 1.0
      IF (DT_FAC == ONE) DT_SAVE = DT 
!
!
      READ (UNIT_RES, REC=1) VERSION 
      READ (VERSION(6:512), *) VERSION_NUMBER 
!
      READ (UNIT_RES, REC=3) NEXT_REC 
      IF (VERSION_NUMBER >= 1.12) THEN 
         READ (UNIT_RES, REC=NEXT_REC) TIME, DT, NSTEP 
      ELSE 
         READ (UNIT_RES, REC=NEXT_REC) TIME, NSTEP 
      ENDIF 
!
      NEXT_REC = NEXT_REC + 1 
!
      CALL IN_BIN_512 (UNIT_RES, array1, IJKMAX2, NEXT_REC) 
      CALL convert_from_io_dp(array1, EP_G, IJKMAX2) 

      CALL IN_BIN_512 (UNIT_RES, array1, IJKMAX2, NEXT_REC) 
      CALL convert_from_io_dp(array1, P_G, IJKMAX2) 

      CALL IN_BIN_512 (UNIT_RES, array1, IJKMAX2, NEXT_REC) 
      CALL convert_from_io_dp(array1, P_STAR, IJKMAX2) 

      CALL IN_BIN_512 (UNIT_RES, array1, IJKMAX2, NEXT_REC) 
      CALL convert_from_io_dp(array1, RO_G, IJKMAX2) 

      CALL IN_BIN_512 (UNIT_RES, array1, IJKMAX2, NEXT_REC) 
      CALL convert_from_io_dp(array1, ROP_G, IJKMAX2) 

      CALL IN_BIN_512 (UNIT_RES, array1, IJKMAX2, NEXT_REC) 
      CALL convert_from_io_dp(array1, T_G, IJKMAX2) 

      IF (VERSION_NUMBER < 1.15) THEN 
         CALL IN_BIN_512 (UNIT_RES, array1, IJKMAX2, NEXT_REC) 
         CALL convert_from_io_dp(array1, T_S(1,1), IJKMAX2) 
         IF (MMAX >= 2) THEN 
            CALL IN_BIN_512 (UNIT_RES, array1, IJKMAX2, NEXT_REC) 
            CALL convert_from_io_dp(array1, T_S(1,2), IJKMAX2) 
         ELSE 
            CALL IN_BIN_512 (UNIT_RES, array1, IJKMAX2, NEXT_REC) 
            CALL convert_from_io_dp(array1, tmp, IJKMAX2) 
         ENDIF 
      ENDIF 
      IF (VERSION_NUMBER >= 1.05) THEN 
         DO N = 1, NMAX(0) 
            CALL IN_BIN_512 (UNIT_RES, array1, IJKMAX2, NEXT_REC) 
            CALL convert_from_io_dp(array1, X_G(1,N), IJKMAX2) 
         END DO 
      ENDIF 
      CALL IN_BIN_512 (UNIT_RES, array1, IJKMAX2, NEXT_REC) 
      CALL convert_from_io_dp(array1, U_G, IJKMAX2) 
      CALL IN_BIN_512 (UNIT_RES, array1, IJKMAX2, NEXT_REC) 
      CALL convert_from_io_dp(array1, V_G, IJKMAX2) 
      CALL IN_BIN_512 (UNIT_RES, array1, IJKMAX2, NEXT_REC) 
      CALL convert_from_io_dp(array1, W_G, IJKMAX2) 
!
      DO LC = 1, MMAX 
         CALL IN_BIN_512 (UNIT_RES, array1, IJKMAX2, NEXT_REC) 
         CALL convert_from_io_dp(array1, ROP_S(1,LC), IJKMAX2) 
         IF (VERSION_NUMBER >= 1.15) then
            CALL IN_BIN_512 (UNIT_RES, array1, &
                         IJKMAX2, NEXT_REC) 
            CALL convert_from_io_dp(array1, T_S(1,LC), IJKMAX2) 
         end if
         CALL IN_BIN_512 (UNIT_RES, array1, IJKMAX2, NEXT_REC) 
         CALL convert_from_io_dp(array1, U_S(1,LC), IJKMAX2) 
         CALL IN_BIN_512 (UNIT_RES, array1, IJKMAX2, NEXT_REC) 
         CALL convert_from_io_dp(array1, V_S(1,LC), IJKMAX2) 
         CALL IN_BIN_512 (UNIT_RES, array1, IJKMAX2, NEXT_REC) 
         CALL convert_from_io_dp(array1, W_S(1,LC), IJKMAX2) 
         IF (VERSION_NUMBER >= 1.2) then
              CALL IN_BIN_512 (UNIT_RES, array1, &
                       IJKMAX2, NEXT_REC) 
              CALL convert_from_io_dp(array1, THETA_M(1,LC), IJKMAX2) 
         end if
         IF (VERSION_NUMBER >= 1.05) THEN 
            DO N = 1, NMAX(LC) 
               CALL IN_BIN_512 (UNIT_RES, array1, IJKMAX2, NEXT_REC) 
               CALL convert_from_io_dp(array1, X_S(1,LC,N), IJKMAX2) 
            END DO 
         ENDIF 
      END DO 
      IF (DT_FAC == ONE) DT = DT_SAVE 
!
      call unlock_tmp_array

      RETURN  
      END SUBROUTINE READ_RES1 
