!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_SPX1                                             C
!  Purpose: write out the time-dependent restart records (REAL)        C
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
!  Variables referenced: TIME, NSTEP, EP_g, P_g, P_star, U_g           C
!                        V_g, W_g, U_s, V_s, W_s, ROP_s, T_g, T_s      C
!                        IJKMAX2, MMAX                                 C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables:  LC, N, NEXT_REC, NUM_REC                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE WRITE_SPX1(L) 
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
      USE compar           !//
      USE mpi_utility      !//
!//       USE tmp_array
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!             flag whether to write a particular SPx file 
      INTEGER L
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!
! local variables
!
!//
      double precision, allocatable :: array1(:)     !//
      double precision, allocatable :: array2(:)     !//

!             loop counters
      INTEGER LC, N
!
!             Pointer to the next record
      INTEGER NEXT_REC
!
!              Number of records written each time step
      INTEGER  NUM_REC
!-----------------------------------------------

      if (myPE.ne.PE_IO) return    !//
!//     call lock_tmp_array
!
      allocate (array1(ijkmax2))   !//
      allocate (array2(ijkmax2))   !//
!
! ".SP1" FILE         EP_g    [ ROP_g, RO_g  must be calculated ...
!                                        not written out ]
!
      SELECT CASE (L)  
      CASE (1)  
         READ (UNIT_SPX + L, REC=3) NEXT_REC, NUM_REC 
         NUM_REC = NEXT_REC 
         WRITE (UNIT_SPX + L, REC=NEXT_REC) REAL(TIME), NSTEP 
         NEXT_REC = NEXT_REC + 1 
         call gather (EP_g,array2,root)   !//
         call convert_to_io_dp(array2,array1,ijkmax2)
         CALL OUT_BIN_R (UNIT_SPX + L, array1, IJKMAX2, NEXT_REC) 
         NUM_REC = NEXT_REC - NUM_REC 
         WRITE (UNIT_SPX + L, REC=3) NEXT_REC, NUM_REC 
         CALL FLUSH (UNIT_SPX + L) 
!
! ".SP2" FILE         P_g , P_star
!
      CASE (2)  
         READ (UNIT_SPX + L, REC=3) NEXT_REC, NUM_REC 
         NUM_REC = NEXT_REC 
         WRITE (UNIT_SPX + L, REC=NEXT_REC) REAL(TIME), NSTEP 
         NEXT_REC = NEXT_REC + 1 
         call gather (P_g,array2,root)   !//
         call convert_to_io_dp(array2,array1,ijkmax2)
         CALL OUT_BIN_R (UNIT_SPX + L, array1, IJKMAX2, NEXT_REC) 
         call gather (P_star,array2,root)   !//
         call convert_to_io_dp(array2,array1,ijkmax2)
         CALL OUT_BIN_R (UNIT_SPX + L, array1, IJKMAX2, NEXT_REC) 
         NUM_REC = NEXT_REC - NUM_REC 
         WRITE (UNIT_SPX + L, REC=3) NEXT_REC, NUM_REC 
         CALL FLUSH (UNIT_SPX + L) 
!
! ".SP3" FILE         U_g , V_g , W_g
!
      CASE (3)  
         READ (UNIT_SPX + L, REC=3) NEXT_REC, NUM_REC 
         NUM_REC = NEXT_REC 
         WRITE (UNIT_SPX + L, REC=NEXT_REC) REAL(TIME), NSTEP 
         NEXT_REC = NEXT_REC + 1 
         call gather (U_g,array2,root)   !//
         call convert_to_io_dp(array2,array1,ijkmax2)
         CALL OUT_BIN_R (UNIT_SPX + L, array1, IJKMAX2, NEXT_REC) 
         call gather (V_g,array2,root)   !//
         call convert_to_io_dp(array2,array1,ijkmax2)
         CALL OUT_BIN_R (UNIT_SPX + L, array1, IJKMAX2, NEXT_REC) 
         call gather (W_g,array2,root)   !//
         call convert_to_io_dp(array2,array1,ijkmax2)
         CALL OUT_BIN_R (UNIT_SPX + L,array1, IJKMAX2, NEXT_REC) 
         NUM_REC = NEXT_REC - NUM_REC 
         WRITE (UNIT_SPX + L, REC=3) NEXT_REC, NUM_REC 
         CALL FLUSH (UNIT_SPX + L) 
!
! ".SP4" FILE         U_s , V_s , W_s
!
      CASE (4)  
         READ (UNIT_SPX + L, REC=3) NEXT_REC, NUM_REC 
         NUM_REC = NEXT_REC 
         WRITE (UNIT_SPX + L, REC=NEXT_REC) REAL(TIME), NSTEP 
         NEXT_REC = NEXT_REC + 1 
         DO LC = 1, MMAX 
            call gather (U_s(:,LC),array2,root)   !//
            call convert_to_io_dp(array2,array1,ijkmax2)
            CALL OUT_BIN_R (UNIT_SPX + L, array1, IJKMAX2, NEXT_REC) 
            call gather (V_s(:,LC),array2,root)   !//
            call convert_to_io_dp(array2,array1,ijkmax2)
            CALL OUT_BIN_R (UNIT_SPX + L,array1 , IJKMAX2, NEXT_REC) 
            call gather (W_s(:,LC),array2,root)   !//
            call convert_to_io_dp(array2,array1,ijkmax2)
            CALL OUT_BIN_R (UNIT_SPX + L,array1, IJKMAX2, NEXT_REC) 
         END DO 
         NUM_REC = NEXT_REC - NUM_REC 
         WRITE (UNIT_SPX + L, REC=3) NEXT_REC, NUM_REC 
         CALL FLUSH (UNIT_SPX + L) 
!
! ".SP5" FILE         ROP_s
!
      CASE (5)  
         READ (UNIT_SPX + L, REC=3) NEXT_REC, NUM_REC 
         NUM_REC = NEXT_REC 
         WRITE (UNIT_SPX + L, REC=NEXT_REC) REAL(TIME), NSTEP 
         NEXT_REC = NEXT_REC + 1 
         DO LC = 1, MMAX 
            call gather (ROP_s(:,LC),array2,root)   !//
            call convert_to_io_dp(array2,array1,ijkmax2)
            CALL OUT_BIN_R (UNIT_SPX + L, array1, IJKMAX2, NEXT_REC) 
         END DO 
         NUM_REC = NEXT_REC - NUM_REC 
         WRITE (UNIT_SPX + L, REC=3) NEXT_REC, NUM_REC 
         CALL FLUSH (UNIT_SPX + L) 
!
! ".SP6" FILE         T_g  , T_s
!
      CASE (6)  
         READ (UNIT_SPX + L, REC=3) NEXT_REC, NUM_REC 
         NUM_REC = NEXT_REC 
         WRITE (UNIT_SPX + L, REC=NEXT_REC) REAL(TIME), NSTEP 
         NEXT_REC = NEXT_REC + 1 
         call gather (T_g,array2,root)   !//
         call convert_to_io_dp(array2,array1,ijkmax2)
         CALL OUT_BIN_R (UNIT_SPX + L, array1, IJKMAX2, NEXT_REC) 
         DO LC = 1, MMAX 
            call gather (T_s(:,LC),array2,root)   !//
            call convert_to_io_dp(array2,array1,ijkmax2)
            CALL OUT_BIN_R (UNIT_SPX + L, array1, IJKMAX2, NEXT_REC) 
         END DO 
         NUM_REC = NEXT_REC - NUM_REC 
         WRITE (UNIT_SPX + L, REC=3) NEXT_REC, NUM_REC 
         CALL FLUSH (UNIT_SPX + L) 
!
! ".SP7" FILE         X_g, X_s
!
      CASE (7)  
         READ (UNIT_SPX + L, REC=3) NEXT_REC, NUM_REC 
         NUM_REC = NEXT_REC 
         WRITE (UNIT_SPX + L, REC=NEXT_REC) REAL(TIME), NSTEP 
         NEXT_REC = NEXT_REC + 1 
         DO N = 1, NMAX(0) 
            call gather (X_G(:,N),array2,root)   !//
            call convert_to_io_dp(array2,array1,ijkmax2)
            CALL OUT_BIN_R (UNIT_SPX + L,array1 , IJKMAX2, NEXT_REC) 
         END DO 
         DO LC = 1, MMAX 
            DO N = 1, NMAX(LC) 
               call gather (X_s(:,LC,N),array2,root)  !//
               call convert_to_io_dp(array2,array1,ijkmax2)
               CALL OUT_BIN_R (UNIT_SPX + L,array1, IJKMAX2, NEXT_REC) 
            END DO 
         END DO 
         NUM_REC = NEXT_REC - NUM_REC 
         WRITE (UNIT_SPX + L, REC=3) NEXT_REC, NUM_REC 
         CALL FLUSH (UNIT_SPX + L) 
!
! ".SP8" FILE         THETA_m
!
      CASE (8)  
         READ (UNIT_SPX + L, REC=3) NEXT_REC, NUM_REC 
         NUM_REC = NEXT_REC 
         WRITE (UNIT_SPX + L, REC=NEXT_REC) REAL(TIME), NSTEP 
         NEXT_REC = NEXT_REC + 1 
         DO LC = 1, MMAX 
            call gather (THETA_m(:,LC),array2,root)   !//
            call convert_to_io_dp(array2,array1,ijkmax2)
            CALL OUT_BIN_R (UNIT_SPX + L,array1 , IJKMAX2, NEXT_REC) 
         END DO 
         NUM_REC = NEXT_REC - NUM_REC 
         WRITE (UNIT_SPX + L, REC=3) NEXT_REC, NUM_REC 
         CALL FLUSH (UNIT_SPX + L) 
!
!
      END SELECT 

!//      call unlock_tmp_array
!
      deallocate (array1)    !//
      deallocate (array2)    !//
!
      RETURN  
      END SUBROUTINE WRITE_SPX1 
