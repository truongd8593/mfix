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
      USE compar      !// 001 Include MPI header file
      USE mpi_utility !//
      USE sendrecv    !//
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
      DOUBLE PRECISION DT_SAVE

!//PAR_I/O 0814 declare global scratch arrays
      double precision, allocatable :: array1(:)
      double precision, allocatable :: array2(:)
      INTEGER allocstatus, i
!-----------------------------------------------
!
      if (myPE .eq. PE_IO) then
         allocate (array1(ijkmax2))
         allocate (array2(ijkmax3))
      else
         allocate (array1(1))
         allocate (array2(1))
      end if

!//SP
      array1(:) = Undefined
      array2(:) = Undefined

      call MPI_barrier(MPI_COMM_WORLD,mpierr)
!
!     Use DT from data file if DT_FAC is set to 1.0
      IF (DT_FAC == ONE) DT_SAVE = DT 
!
!

!//PAR_I/O 0813 only PE_IO reads the restart file
      if (myPE == PE_IO ) then
         READ (UNIT_RES, REC=1) VERSION 
         READ (VERSION(6:512), *) VERSION_NUMBER 

         READ (UNIT_RES, REC=3) NEXT_REC 
         IF (VERSION_NUMBER >= 1.12) THEN 
            READ (UNIT_RES, REC=NEXT_REC) TIME, DT, NSTEP 
         ELSE 
            READ (UNIT_RES, REC=NEXT_REC) TIME, NSTEP 
         ENDIF 
      end if
!
      call MPI_barrier(MPI_COMM_WORLD,mpierr)
      call bcast(VERSION, PE_IO)        !//PAR_I/O BCAST0c
      call bcast(VERSION_NUMBER, PE_IO) !//PAR_I/O BCAST0r
      call bcast(TIME, PE_IO)           !//PAR_I/O BCAST0d
      call bcast(NSTEP, PE_IO)          !//PAR_I/O BCAST0i
      if (VERSION_NUMBER >= 1.12) call bcast(DT, PE_IO)   !//PAR_I/O BCAST0d	     
      call MPI_barrier(MPI_COMM_WORLD,mpierr)

!
      if (myPE == PE_IO) then
        NEXT_REC = NEXT_REC + 1 
        CALL IN_BIN_512 (UNIT_RES, array1, IJKMAX2, NEXT_REC) 
        CALL convert_from_io_dp(array1, array2, IJKMAX2) 
      end if
      call MPI_barrier(MPI_COMM_WORLD,mpierr)
      call scatter(EP_G,array2,PE_IO)
      call MPI_barrier(MPI_COMM_WORLD,mpierr)

      if (myPE == PE_IO) then
        CALL IN_BIN_512 (UNIT_RES, array1, IJKMAX2, NEXT_REC) 
        CALL convert_from_io_dp(array1, array2, IJKMAX2) 
      end if
      call MPI_barrier(MPI_COMM_WORLD,mpierr)
      call scatter(P_G,array2,PE_IO)
      call MPI_barrier(MPI_COMM_WORLD,mpierr)

      if (myPE == PE_IO) then
        CALL IN_BIN_512 (UNIT_RES, array1, IJKMAX2, NEXT_REC) 
        CALL convert_from_io_dp(array1, array2, IJKMAX2) 
      end if
      call MPI_barrier(MPI_COMM_WORLD,mpierr)
      call scatter(P_STAR,array2,PE_IO)
      call MPI_barrier(MPI_COMM_WORLD,mpierr)

      if (myPE == PE_IO) then
        CALL IN_BIN_512 (UNIT_RES, array1, IJKMAX2, NEXT_REC) 
        CALL convert_from_io_dp(array1, array2, IJKMAX2) 
      end if
      call MPI_barrier(MPI_COMM_WORLD,mpierr)
      call scatter(RO_G,array2,PE_IO)
      call MPI_barrier(MPI_COMM_WORLD,mpierr)

      if (myPE == PE_IO) then
         CALL IN_BIN_512 (UNIT_RES, array1, IJKMAX2, NEXT_REC) 
         CALL convert_from_io_dp(array1, array2, IJKMAX2) 
      end if
      call MPI_barrier(MPI_COMM_WORLD,mpierr)
      call scatter(ROP_G,array2,PE_IO)
      call MPI_barrier(MPI_COMM_WORLD,mpierr)

      if (myPE == PE_IO) then
         CALL IN_BIN_512 (UNIT_RES, array1, IJKMAX2, NEXT_REC) 
         CALL convert_from_io_dp(array1, array2, IJKMAX2) 
      end if
      call MPI_barrier(MPI_COMM_WORLD,mpierr)
      call scatter(T_G,array2,PE_IO)
      call MPI_barrier(MPI_COMM_WORLD,mpierr)

!

      IF (VERSION_NUMBER < 1.15) THEN 
         if (myPE == PE_IO) then
            CALL IN_BIN_512 (UNIT_RES, array1, IJKMAX2, NEXT_REC) 
            CALL convert_from_io_dp(array1, array2, IJKMAX2)
         end if 
         call MPI_barrier(MPI_COMM_WORLD,mpierr)
         call scatter (T_s(:,1),array2,PE_IO)
         call MPI_barrier(MPI_COMM_WORLD,mpierr)

         IF (MMAX >= 2) THEN 
            if (myPE == PE_IO) then
               CALL IN_BIN_512 (UNIT_RES, array1, IJKMAX2, NEXT_REC) 
               CALL convert_from_io_dp(array1, array2, IJKMAX2) 
            end if
            call MPI_barrier(MPI_COMM_WORLD,mpierr)
            call scatter (T_s(:,2),array2,PE_IO)
            call MPI_barrier(MPI_COMM_WORLD,mpierr)
         ELSE 
            if (myPE == PE_IO) &
                   CALL IN_BIN_512 (UNIT_RES, array1, IJKMAX2, NEXT_REC) 
         ENDIF 
      ENDIF 

      IF (VERSION_NUMBER >= 1.05) THEN 
         DO N = 1, NMAX(0) 
            if (myPE == PE_IO) then
               CALL IN_BIN_512 (UNIT_RES, array1, IJKMAX2, NEXT_REC) 
               CALL convert_from_io_dp(array1, array2, IJKMAX2)
            end if 
            call MPI_barrier(MPI_COMM_WORLD,mpierr)
            call scatter (X_g(:,n),array2,PE_IO)
            call MPI_barrier(MPI_COMM_WORLD,mpierr)
         END DO 
      ENDIF 

      if (myPE == PE_IO) then
         CALL IN_BIN_512 (UNIT_RES, array1, IJKMAX2, NEXT_REC) 
         CALL convert_from_io_dp(array1, array2, IJKMAX2) 
      end if
      call MPI_barrier(MPI_COMM_WORLD,mpierr)
      call scatter(U_G, array2, PE_IO)
      call MPI_barrier(MPI_COMM_WORLD,mpierr)

      if (myPE == PE_IO) then
         CALL IN_BIN_512 (UNIT_RES, array1, IJKMAX2, NEXT_REC) 
         CALL convert_from_io_dp(array1, array2, IJKMAX2) 
      end if
      call MPI_barrier(MPI_COMM_WORLD,mpierr)
      call scatter(V_G, array2, PE_IO)
      call MPI_barrier(MPI_COMM_WORLD,mpierr)

      if (myPE == PE_IO) then
         CALL IN_BIN_512 (UNIT_RES, array1, IJKMAX2, NEXT_REC) 
         CALL convert_from_io_dp(array1, array2, IJKMAX2) 
      end if
      call MPI_barrier(MPI_COMM_WORLD,mpierr)
      call scatter(W_G, array2, PE_IO)
      call MPI_barrier(MPI_COMM_WORLD,mpierr)

!
      DO LC = 1, MMAX 

         if (myPE == PE_IO) then
            CALL IN_BIN_512 (UNIT_RES, array1, IJKMAX2, NEXT_REC) 
            CALL convert_from_io_dp(array1, array2, IJKMAX2) 
         end if
         call MPI_barrier(MPI_COMM_WORLD,mpierr)
         call scatter(ROP_S(:,LC), array2, PE_IO)
         call MPI_barrier(MPI_COMM_WORLD,mpierr)

         IF (VERSION_NUMBER >= 1.15) THEN
            if (myPE == PE_IO) then
               CALL IN_BIN_512 (UNIT_RES, array1, IJKMAX2, NEXT_REC) 
               CALL convert_from_io_dp(array1, array2, IJKMAX2) 
            end if
            call MPI_barrier(MPI_COMM_WORLD,mpierr)
            call scatter(T_S(:,LC), array2, PE_IO)
            call MPI_barrier(MPI_COMM_WORLD,mpierr)
         END IF

         if (myPE == PE_IO) then
            CALL IN_BIN_512 (UNIT_RES, array1, IJKMAX2, NEXT_REC) 
            CALL convert_from_io_dp(array1, array2, IJKMAX2) 
         end if
         call MPI_barrier(MPI_COMM_WORLD,mpierr)
         call scatter(U_S(:,LC), array2, PE_IO)
         call MPI_barrier(MPI_COMM_WORLD,mpierr)

         if (myPE == PE_IO) then
            CALL IN_BIN_512 (UNIT_RES, array1, IJKMAX2, NEXT_REC) 
            CALL convert_from_io_dp(array1, array2, IJKMAX2)
         end if
         call MPI_barrier(MPI_COMM_WORLD,mpierr)
         call scatter(V_S(:,LC), array2, PE_IO)
         call MPI_barrier(MPI_COMM_WORLD,mpierr)

         if (myPE == PE_IO) then
            CALL IN_BIN_512 (UNIT_RES, array1, IJKMAX2, NEXT_REC) 
            CALL convert_from_io_dp(array1, array2, IJKMAX2) 
         end if
         call MPI_barrier(MPI_COMM_WORLD,mpierr)
         call scatter(W_S(:,LC), array2, PE_IO)
         call MPI_barrier(MPI_COMM_WORLD,mpierr)

         IF (VERSION_NUMBER >= 1.2) then
            if (myPE == PE_IO) then
               CALL IN_BIN_512 (UNIT_RES, array1, IJKMAX2, NEXT_REC) 
               CALL convert_from_io_dp(array1, array2, IJKMAX2) 
            end if
            call MPI_barrier(MPI_COMM_WORLD,mpierr)
            call scatter(THETA_M(:,LC), array2, PE_IO)
            call MPI_barrier(MPI_COMM_WORLD,mpierr)
         end if
         IF (VERSION_NUMBER >= 1.05) THEN 
            DO N = 1, NMAX(LC) 
               if (myPE == PE_IO) then
                  CALL IN_BIN_512 (UNIT_RES, array1, IJKMAX2, NEXT_REC) 
                  CALL convert_from_io_dp(array1, array2, IJKMAX2) 
               end if
               call MPI_barrier(MPI_COMM_WORLD,mpierr)
               call scatter(X_S(:,LC,N), array2, PE_IO)
               call MPI_barrier(MPI_COMM_WORLD,mpierr)
            END DO 
         ENDIF 
      END DO 


      call MPI_barrier(MPI_COMM_WORLD,mpierr)
      deallocate( array1 )
      deallocate( array2 )
      call MPI_barrier(MPI_COMM_WORLD,mpierr)

!    else

!//TD need bcast_0c   call bcast(VERSION, PE_IO)        !//PAR_I/O BCAST0c (recv)
!//TD need bcast_0c	  call bcast(VERSION_NUMBER, PE_IO) !//PAR_I/O BCAST0r (recv)
!      call bcast(TIME, PE_IO)        !//PAR_I/O BCAST0d (recv)
!      call bcast(NSTEP, PE_IO)       !//PAR_I/O BCAST0i (recv)
!      if (VERSION_NUMBER >= 1.12) &
!         call bcast(DT, PE_IO)       !//PAR_I/O BCAST0d	(recv)

!      call gather(EP_G,dGTEMP2,PE_IO)   !//PAR_I/O GATHER1d
!      call gather(P_G,dGTEMP2,PE_IO)    !//PAR_I/O GATHER1d
!      call gather(P_STAR,dGTEMP2,PE_IO) !//PAR_I/O GATHER1d
!      call gather(RO_G,dGTEMP2,PE_IO)   !//PAR_I/O GATHER1d
!      call gather(ROP_G,dGTEMP2,PE_IO)  !//PAR_I/O GATHER1d
!      call gather(T_G,dGTEMP2,PE_IO)    !//PAR_I/O GATHER1d
!      if (MMAX < 2) then 
!        call gather(tmp,dGTEMP2,PE_IO)  !//PAR_I/O GATHER1d
!      endif
!      call gather(T_S,dGTEMP3,PE_IO)    !//PAR_I/O GATHER2d

!      if (VERSION_NUMBER >= 1.05) then
!         call gather(X_G, dGTEMP3, PE_IO)     !//PAR_I/O GATHER2d
!      endif

!      call gather(U_G,dGTEMP2,PE_IO)   !//PAR_I/O GATHER1d
!      call gather(V_G,dGTEMP2,PE_IO)   !//PAR_I/O GATHER1d
!      call gather(W_G,dGTEMP2,PE_IO)   !//PAR_I/O GATHER1d

!      DO LC = 1, MMAX 
 !        call gather(ROP_S, dGTEMP3, PE_IO)        !//PAR_I/O GATHER2d

!         IF (VERSION_NUMBER >= 1.15) THEN
!            call gather(T_S, dGTEMP3, PE_IO)        !//PAR_I/O GATHER2d
 !        END IF
!         call gather(U_S, dGTEMP3, PE_IO)        !//PAR_I/O GATHER2d
!         call gather(V_S, dGTEMP3, PE_IO)        !//PAR_I/O GATHER2d
!         call gather(W_S, dGTEMP3, PE_IO)        !//PAR_I/O GATHER2d

!         IF (VERSION_NUMBER >= 1.2) then
!              call gather(THETA_M, dGTEMP3, PE_IO)     !//PAR_I/O GATHER2d
!         end if
!         IF (VERSION_NUMBER >= 1.05) THEN 
!            DO N = 1, NMAX(LC) 
!               call gather(X_S, dGTEMP4, PE_IO)        !//PAR_I/O GATHER2d
!            END DO 
!         ENDIF 
!      END DO

!    endif


!//SP
      call send_recv(rop_g)
      call send_recv(ro_g)
      call send_recv(rop_s)


      IF (DT_FAC == ONE) DT = DT_SAVE 
!

      RETURN  
      END SUBROUTINE READ_RES1 
