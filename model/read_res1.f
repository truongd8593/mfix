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
      USE compar      !// 001 Include MPI header file
      USE mpi_utility !//
!//TD      USE dbg_utility !//PARDBG
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

!//PAR_I/O 0814 declare global scratch arrays
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: dGTEMP, dGTEMP2  !//PAR_I/O declare real*8 Global SCRatch array
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: dGTEMP3 !//PAR_I/O declare real*8 Global SCRatch array
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: dGTEMP4 !//PAR_I/O declare real*8 Global SCRatch array
      INTEGER allocstatus, i
!-----------------------------------------------
!
      call lock_tmp_array
!
!     Use DT from data file if DT_FAC is set to 1.0
      IF (DT_FAC == ONE) DT_SAVE = DT 
!
!

!//PAR_I/O 0813 only PE_IO reads the restart file
    if (myPE == PE_IO ) then
      READ (UNIT_RES, REC=1) VERSION 
      READ (VERSION(6:512), *) VERSION_NUMBER 
!
!//TD need bcast_0c   call bcast(VERSION, PE_IO)        !//PAR_I/O BCAST0c
!//TD need bcast_0c	  call bcast(VERSION_NUMBER, PE_IO) !//PAR_I/O BCAST0r

      READ (UNIT_RES, REC=3) NEXT_REC 
      IF (VERSION_NUMBER >= 1.12) THEN 
         READ (UNIT_RES, REC=NEXT_REC) TIME, DT, NSTEP 
      ELSE 
         READ (UNIT_RES, REC=NEXT_REC) TIME, NSTEP 
      ENDIF 
      call bcast(TIME, PE_IO)        !//PAR_I/O BCAST0d
      call bcast(NSTEP, PE_IO)       !//PAR_I/O BCAST0i
      if (VERSION_NUMBER >= 1.12) &
         call bcast(DT, PE_IO)       !//PAR_I/O BCAST0d	     

!
      NEXT_REC = NEXT_REC + 1 
!
!//? are we preserving the index ijkmax2 as global ijkmax2?
!//S if array1 is not used somewhere else instead of defining a new array, dGTEMP 
!//S (see read_res0 for array1 usage)
!//S could we use array1 as the global scratch array for reading in the record
      Allocate( dGTEMP(IJKMAX2), STAT=allocstatus) !//PAR_I/O ALLOCate R*8 Global scratch
      Allocate( dGTEMP2(IJKMAX2), STAT=allocstatus) !//PAR_I/O ALLOCate R*8 Global scratch
      Allocate( dGTEMP3(IJKMAX2,DIMENSION_M), STAT=allocstatus) !//PAR_I/O ALLOCate R*8 Global scratch
      Allocate( dGTEMP4(IJKMAX2,DIMENSION_M,DIMENSION_N_s), STAT=allocstatus) !//PAR_I/O ALLOCate R*8 Global scratch
!      CALL IN_BIN_512 (UNIT_RES, array1, IJKMAX2, NEXT_REC) 
      CALL IN_BIN_512 (UNIT_RES, dGTEMP, IJKMAX2, NEXT_REC) 
!      CALL convert_from_io_dp(array1, EP_G, IJKMAX2) 
      CALL convert_from_io_dp(dGTEMP, dGTEMP2, IJKMAX2) 
      call scatter(EP_G,dGTEMP2,PE_IO)        !//PAR_I/O SCATTER1d

!//? should we initialize dGTEMP, dGTEMP2 to zero before reading next rec for safety purposes?

      CALL IN_BIN_512 (UNIT_RES, dGTEMP, IJKMAX2, NEXT_REC) 
      CALL convert_from_io_dp(dGTEMP, dGTEMP2, IJKMAX2) 
      call scatter(P_G,dGTEMP2,PE_IO)        !//PAR_I/O SCATTER1d

      CALL IN_BIN_512 (UNIT_RES, dGTEMP, IJKMAX2, NEXT_REC) 
      CALL convert_from_io_dp(dGTEMP, dGTEMP2, IJKMAX2) 
      call scatter(P_STAR,dGTEMP2,PE_IO)     !//PAR_I/O SCATTER1d

      CALL IN_BIN_512 (UNIT_RES, dGTEMP, IJKMAX2, NEXT_REC) 
      CALL convert_from_io_dp(dGTEMP, dGTEMP2, IJKMAX2) 
      call scatter(RO_G,dGTEMP2,PE_IO)       !//PAR_I/O SCATTER1d

      CALL IN_BIN_512 (UNIT_RES, dGTEMP, IJKMAX2, NEXT_REC) 
      CALL convert_from_io_dp(dGTEMP, dGTEMP2, IJKMAX2) 
      call scatter(ROP_G,dGTEMP2,PE_IO)      !//PAR_I/O SCATTER1d

      CALL IN_BIN_512 (UNIT_RES, dGTEMP, IJKMAX2, NEXT_REC) 
      CALL convert_from_io_dp(dGTEMP, dGTEMP2, IJKMAX2) 
      call scatter(T_G,dGTEMP2,PE_IO)        !//PAR_I/O SCATTER1d

!//? make sure 2D array mapping to 1D scratch is done properly
      IF (VERSION_NUMBER < 1.15) THEN 
         CALL IN_BIN_512 (UNIT_RES, dGTEMP, IJKMAX2, NEXT_REC) 
         CALL convert_from_io_dp(dGTEMP, dGTEMP3(1,1), IJKMAX2) 

         IF (MMAX >= 2) THEN 
            CALL IN_BIN_512 (UNIT_RES, dGTEMP, IJKMAX2, NEXT_REC) 
            CALL convert_from_io_dp(dGTEMP, dGTEMP3(1,2), IJKMAX2) 
         ELSE 
            CALL IN_BIN_512 (UNIT_RES, dGTEMP, IJKMAX2, NEXT_REC) 
            CALL convert_from_io_dp(dGTEMP, dGTEMP2, IJKMAX2) 
            call scatter(tmp,dGTEMP2,PE_IO)    !//PAR_I/O SCATTER2d
         ENDIF 

         call scatter(T_S,dGTEMP3,PE_IO)       !//PAR_I/O SCATTER2d
      ENDIF 

      IF (VERSION_NUMBER >= 1.05) THEN 
         DO N = 1, NMAX(0) 
            CALL IN_BIN_512 (UNIT_RES, dGTEMP, IJKMAX2, NEXT_REC) 
            CALL convert_from_io_dp(dGTEMP, dGTEMP3(1,N), IJKMAX2) 
         END DO 
         call scatter(X_G, dGTEMP3, PE_IO)     !//PAR_I/O SCATTER2d
      ENDIF 

      CALL IN_BIN_512 (UNIT_RES, dGTEMP, IJKMAX2, NEXT_REC) 
      CALL convert_from_io_dp(dGTEMP, dGTEMP2, IJKMAX2) 
      call scatter(U_G, dGTEMP2, PE_IO)        !//PAR_I/O SCATTER2d

      CALL IN_BIN_512 (UNIT_RES, dGTEMP, IJKMAX2, NEXT_REC) 
      CALL convert_from_io_dp(dGTEMP, dGTEMP2, IJKMAX2) 
      call scatter(V_G, dGTEMP2, PE_IO)        !//PAR_I/O SCATTER2d

      CALL IN_BIN_512 (UNIT_RES, dGTEMP, IJKMAX2, NEXT_REC) 
      CALL convert_from_io_dp(dGTEMP, dGTEMP2, IJKMAX2) 
      call scatter(W_G, dGTEMP2, PE_IO)        !//PAR_I/O SCATTER2d

!
      DO LC = 1, MMAX 
         CALL IN_BIN_512 (UNIT_RES, dGTEMP, IJKMAX2, NEXT_REC) 
!         CALL convert_from_io_dp(array1, ROP_S(1,LC), IJKMAX2) 
         CALL convert_from_io_dp(dGTEMP, dGTEMP3(1,LC), IJKMAX2) 
!//? due to the structure of this DO loop, can not do the scatter out of the loop
!//? on the other hand (1,LC) needs to be scattered but all full array will be 
!//? done if implemented like this
         call scatter(ROP_S, dGTEMP3, PE_IO)        !//PAR_I/O SCATTER2d

         IF (VERSION_NUMBER >= 1.15) THEN
            CALL IN_BIN_512 (UNIT_RES, dGTEMP, &
                         IJKMAX2, NEXT_REC) 
            CALL convert_from_io_dp(dGTEMP, dGTEMP3(1,LC), IJKMAX2) 
            call scatter(T_S, dGTEMP3, PE_IO)        !//PAR_I/O SCATTER2d
         END IF

         CALL IN_BIN_512 (UNIT_RES, dGTEMP, IJKMAX2, NEXT_REC) 
         CALL convert_from_io_dp(dGTEMP, dGTEMP3(1,LC), IJKMAX2) 
         call scatter(U_S, dGTEMP3, PE_IO)        !//PAR_I/O SCATTER2d

         CALL IN_BIN_512 (UNIT_RES, dGTEMP, IJKMAX2, NEXT_REC) 
         CALL convert_from_io_dp(dGTEMP, dGTEMP3(1,LC), IJKMAX2) 
         call scatter(V_S, dGTEMP3, PE_IO)        !//PAR_I/O SCATTER2d

         CALL IN_BIN_512 (UNIT_RES, dGTEMP, IJKMAX2, NEXT_REC) 
         CALL convert_from_io_dp(dGTEMP, dGTEMP3(1,LC), IJKMAX2) 
         call scatter(W_S, dGTEMP3, PE_IO)        !//PAR_I/O SCATTER2d

         IF (VERSION_NUMBER >= 1.2) then
              CALL IN_BIN_512 (UNIT_RES, dGTEMP, &
                       IJKMAX2, NEXT_REC) 
              CALL convert_from_io_dp(dGTEMP, dGTEMP3(1,LC), IJKMAX2) 
              call scatter(THETA_M, dGTEMP3, PE_IO)        !//PAR_I/O SCATTER2d
         end if
         IF (VERSION_NUMBER >= 1.05) THEN 
            DO N = 1, NMAX(LC) 
               CALL IN_BIN_512 (UNIT_RES, dGTEMP, IJKMAX2, NEXT_REC) 
!               CALL convert_from_io_dp(array1, X_S(1,LC,N), IJKMAX2) 
               CALL convert_from_io_dp(dGTEMP, dGTEMP4(1,LC,N), IJKMAX2) 
               call scatter(X_S, dGTEMP4, PE_IO)        !//PAR_I/O SCATTER2d
            END DO 
         ENDIF 
      END DO 


      DeAllocate( dGTEMP, dGTEMP2, dGTEMP3, dGTEMP4)

    else

!//TD need bcast_0c   call bcast(VERSION, PE_IO)        !//PAR_I/O BCAST0c (recv)
!//TD need bcast_0c	  call bcast(VERSION_NUMBER, PE_IO) !//PAR_I/O BCAST0r (recv)
      call bcast(TIME, PE_IO)        !//PAR_I/O BCAST0d (recv)
      call bcast(NSTEP, PE_IO)       !//PAR_I/O BCAST0i (recv)
      if (VERSION_NUMBER >= 1.12) &
         call bcast(DT, PE_IO)       !//PAR_I/O BCAST0d	(recv)

      call gather(EP_G,dGTEMP2,PE_IO)   !//PAR_I/O GATHER1d
      call gather(P_G,dGTEMP2,PE_IO)    !//PAR_I/O GATHER1d
      call gather(P_STAR,dGTEMP2,PE_IO) !//PAR_I/O GATHER1d
      call gather(RO_G,dGTEMP2,PE_IO)   !//PAR_I/O GATHER1d
      call gather(ROP_G,dGTEMP2,PE_IO)  !//PAR_I/O GATHER1d
      call gather(T_G,dGTEMP2,PE_IO)    !//PAR_I/O GATHER1d
      if (MMAX < 2) then 
        call gather(tmp,dGTEMP2,PE_IO)  !//PAR_I/O GATHER1d
      endif
      call gather(T_S,dGTEMP3,PE_IO)    !//PAR_I/O GATHER2d

      if (VERSION_NUMBER >= 1.05) then
         call gather(X_G, dGTEMP3, PE_IO)     !//PAR_I/O GATHER2d
      endif

      call gather(U_G,dGTEMP2,PE_IO)   !//PAR_I/O GATHER1d
      call gather(V_G,dGTEMP2,PE_IO)   !//PAR_I/O GATHER1d
      call gather(W_G,dGTEMP2,PE_IO)   !//PAR_I/O GATHER1d

      DO LC = 1, MMAX 
         call gather(ROP_S, dGTEMP3, PE_IO)        !//PAR_I/O GATHER2d

         IF (VERSION_NUMBER >= 1.15) THEN
            call gather(T_S, dGTEMP3, PE_IO)        !//PAR_I/O GATHER2d
         END IF
         call gather(U_S, dGTEMP3, PE_IO)        !//PAR_I/O GATHER2d
         call gather(V_S, dGTEMP3, PE_IO)        !//PAR_I/O GATHER2d
         call gather(W_S, dGTEMP3, PE_IO)        !//PAR_I/O GATHER2d

         IF (VERSION_NUMBER >= 1.2) then
              call gather(THETA_M, dGTEMP3, PE_IO)     !//PAR_I/O GATHER2d
         end if
         IF (VERSION_NUMBER >= 1.05) THEN 
            DO N = 1, NMAX(LC) 
               call gather(X_S, dGTEMP4, PE_IO)        !//PAR_I/O GATHER2d
            END DO 
         ENDIF 
      END DO

    endif





      IF (DT_FAC == ONE) DT = DT_SAVE 
!
      call unlock_tmp_array

      RETURN  
      END SUBROUTINE READ_RES1 
