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
      SUBROUTINE WRITE_SPX1(L, unit_add) 
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
      USE scalars
      USE output
      USE rxns
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

!              offset for use in post_mfix
      INTEGER  unit_add
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
      
      INTEGER  uspx   ! UNIT_SPX + offset from post_mfix
      CHARACTER, DIMENSION(1) :: LINE*50   !error message
!-----------------------------------------------
      uspx = UNIT_SPX + unit_add

!
      if (myPE .eq.PE_IO) then
         allocate (array1(ijkmax2))   !//
         allocate (array2(ijkmax3))   !//
      else
         allocate (array1(1))   !//
         allocate (array2(1))   !//
      end if
!
! ".SP1" FILE         EP_g    [ ROP_g, RO_g  must be calculated ...
!                                        not written out ]
!
      SELECT CASE (L)  
      CASE (1)  
         if (myPE.eq.PE_IO) then 
            READ (uspx + L, REC=3) NEXT_REC, NUM_REC 
            NUM_REC = NEXT_REC 
            WRITE (uspx + L, REC=NEXT_REC) REAL(TIME), NSTEP 
            NEXT_REC = NEXT_REC + 1 
         end if
         call gatherWriteSpx (EP_g,array2, array1, uspx+L, NEXT_REC)   !//
         if (myPE .eq. PE_IO) then
            NUM_REC = NEXT_REC - NUM_REC 
            WRITE (uspx + L, REC=3) NEXT_REC, NUM_REC 
            if(unit_add == 0) CALL FLUSH_bin (uspx + L) 
         end if
!        call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
!
! ".SP2" FILE         P_g , P_star
!
      CASE (2)  
         if (myPE.eq.PE_IO) then 
            READ (uspx + L, REC=3) NEXT_REC, NUM_REC 
            NUM_REC = NEXT_REC 
            WRITE (uspx + L, REC=NEXT_REC) REAL(TIME), NSTEP 
            NEXT_REC = NEXT_REC + 1 
         end if
         call gatherWriteSpx (P_g,array2, array1, uspx+L, NEXT_REC)   !//
         call gatherWriteSpx (P_star,array2, array1, uspx+L, NEXT_REC)   !//
         if (myPE.eq.PE_IO) then 
            NUM_REC = NEXT_REC - NUM_REC 
            WRITE (uspx + L, REC=3) NEXT_REC, NUM_REC 
            if(unit_add == 0) CALL FLUSH_bin (uspx + L) 
         end if
!        call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
!
! ".SP3" FILE         U_g , V_g , W_g
!
      CASE (3)  
         if (myPE.eq.PE_IO) then 
            READ (uspx + L, REC=3) NEXT_REC, NUM_REC 
            NUM_REC = NEXT_REC 
            WRITE (uspx + L, REC=NEXT_REC) REAL(TIME), NSTEP 
            NEXT_REC = NEXT_REC + 1
         end if 
         call gatherWriteSpx (U_g,array2, array1, uspx+L, NEXT_REC)   !//
         call gatherWriteSpx (V_g,array2, array1, uspx+L, NEXT_REC)   !//
         call gatherWriteSpx (W_g,array2, array1, uspx+L, NEXT_REC)   !//
         if (myPE.eq.PE_IO) then
            NUM_REC = NEXT_REC - NUM_REC 
            WRITE (uspx + L, REC=3) NEXT_REC, NUM_REC 
            if(unit_add == 0) CALL FLUSH_bin (uspx + L) 
         end if
!        call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
!
! ".SP4" FILE         U_s , V_s , W_s
!
      CASE (4)  
         if (myPE.eq.PE_IO) then
            READ (uspx + L, REC=3) NEXT_REC, NUM_REC 
            NUM_REC = NEXT_REC 
            WRITE (uspx + L, REC=NEXT_REC) REAL(TIME), NSTEP 
            NEXT_REC = NEXT_REC + 1 
         end if
         DO LC = 1, MMAX 
            call gatherWriteSpx (U_s(:,LC),array2, array1, uspx+L, NEXT_REC)
            call gatherWriteSpx (V_s(:,LC),array2, array1, uspx+L, NEXT_REC)
            call gatherWriteSpx (W_s(:,LC),array2, array1, uspx+L, NEXT_REC)
         END DO 
         if (myPE.eq.PE_IO) then
            NUM_REC = NEXT_REC - NUM_REC 
            WRITE (uspx + L, REC=3) NEXT_REC, NUM_REC 
            if(unit_add == 0) CALL FLUSH_bin (uspx + L) 
         end if
!        call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
!
! ".SP5" FILE         ROP_s
!
      CASE (5)  
         if (myPE.eq.PE_IO) then
            READ (uspx + L, REC=3) NEXT_REC, NUM_REC 
            NUM_REC = NEXT_REC 
            WRITE (uspx + L, REC=NEXT_REC) REAL(TIME), NSTEP 
            NEXT_REC = NEXT_REC + 1 
         end if
         DO LC = 1, MMAX 
            call gatherWriteSpx (ROP_s(:,LC),array2, array1, uspx+L, NEXT_REC)
         END DO 
         if (myPE.eq.PE_IO) then
            NUM_REC = NEXT_REC - NUM_REC 
            WRITE (uspx + L, REC=3) NEXT_REC, NUM_REC 
            if(unit_add == 0) CALL FLUSH_bin (uspx + L) 
         end if
!        call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
!
! ".SP6" FILE         T_g  , T_s
!
      CASE (6)  
         if (myPE.eq.PE_IO) then
            READ (uspx + L, REC=3) NEXT_REC, NUM_REC 
            NUM_REC = NEXT_REC 
            WRITE (uspx + L, REC=NEXT_REC) REAL(TIME), NSTEP 
            NEXT_REC = NEXT_REC + 1 
         end if
         call gatherWriteSpx (T_g,array2, array1, uspx+L, NEXT_REC)   !//
         DO LC = 1, MMAX 
            call gatherWriteSpx (T_s(:,LC),array2, array1, uspx+L, NEXT_REC)
         END DO 
         if (myPE.eq.PE_IO) then
            NUM_REC = NEXT_REC - NUM_REC 
            WRITE (uspx + L, REC=3) NEXT_REC, NUM_REC 
            if(unit_add == 0) CALL FLUSH_bin (uspx + L) 
         end if
!
! ".SP7" FILE         X_g, X_s
!
      CASE (7)  
         if (myPE.eq.PE_IO) then
            READ (uspx + L, REC=3) NEXT_REC, NUM_REC 
            NUM_REC = NEXT_REC 
            WRITE (uspx + L, REC=NEXT_REC) REAL(TIME), NSTEP 
            NEXT_REC = NEXT_REC + 1 
         end if
         DO N = 1, NMAX(0) 
            call gatherWriteSpx (X_G(:,N),array2, array1, uspx+L, NEXT_REC)
         END DO 
         DO LC = 1, MMAX 
            DO N = 1, NMAX(LC) 
               call gatherWriteSpx (X_s(:,LC,N),array2, array1, uspx+L, NEXT_REC)
            END DO 
         END DO 
         if (myPE.eq.PE_IO) then
            NUM_REC = NEXT_REC - NUM_REC 
            WRITE (uspx + L, REC=3) NEXT_REC, NUM_REC 
            if(unit_add == 0) CALL FLUSH_bin (uspx + L) 
         end if
!        call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
!
! ".SP8" FILE         THETA_m
!
      CASE (8)  
         if (myPE.eq.PE_IO) then
            READ (uspx + L, REC=3) NEXT_REC, NUM_REC 
            NUM_REC = NEXT_REC 
            WRITE (uspx + L, REC=NEXT_REC) REAL(TIME), NSTEP 
            NEXT_REC = NEXT_REC + 1
         end if 
         DO LC = 1, MMAX 
            call gatherWriteSpx (THETA_m(:,LC),array2, array1, uspx+L, NEXT_REC)
         END DO 
         if (myPE.eq.PE_IO) then
            NUM_REC = NEXT_REC - NUM_REC 
            WRITE (uspx + L, REC=3) NEXT_REC, NUM_REC 
            if(unit_add == 0) CALL FLUSH_bin (uspx + L) 
         end if
!        call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
!
! ".SP9" FILE         Scalar
!
      CASE (9)  
         if (myPE.eq.PE_IO) then
            READ (uspx + L, REC=3) NEXT_REC, NUM_REC 
            NUM_REC = NEXT_REC 
            WRITE (uspx + L, REC=NEXT_REC) REAL(TIME), NSTEP 
            NEXT_REC = NEXT_REC + 1
         end if 
         DO LC = 1, Nscalar 
            call gatherWriteSpx (Scalar(:,LC),array2, array1, uspx+L, NEXT_REC) 
         END DO 
         if (myPE.eq.PE_IO) then
            NUM_REC = NEXT_REC - NUM_REC 
            WRITE (uspx + L, REC=3) NEXT_REC, NUM_REC 
            if(unit_add == 0) CALL FLUSH_bin (uspx + L) 
         end if
!        call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
!
      CASE (10)

         if (myPE.eq.PE_IO) then
            READ (uspx + L, REC=3) NEXT_REC, NUM_REC 
            NUM_REC = NEXT_REC 
            WRITE (uspx + L, REC=NEXT_REC) REAL(TIME), NSTEP 
            NEXT_REC = NEXT_REC + 1
         end if 
         DO LC = 1, nRR
            call gatherWriteSpx (ReactionRates(:,LC),array2, array1, uspx+L, NEXT_REC) 
         END DO 
         if (myPE.eq.PE_IO) then
            NUM_REC = NEXT_REC - NUM_REC 
            WRITE (uspx + L, REC=3) NEXT_REC, NUM_REC 
            if(unit_add == 0) CALL FLUSH_bin (uspx + L) 
         end if
!
! ".SP11" FILE         turbulence
!
      CASE (11)  
         if (myPE.eq.PE_IO) then
            READ (uspx + L, REC=3) NEXT_REC, NUM_REC 
            NUM_REC = NEXT_REC 
            WRITE (uspx + L, REC=NEXT_REC) REAL(TIME), NSTEP 
            NEXT_REC = NEXT_REC + 1
         end if 
	 if (K_Epsilon) then
            call gatherWriteSpx (K_Turb_G,array2, array1, uspx+L, NEXT_REC) 
            call gatherWriteSpx (E_Turb_G,array2, array1, uspx+L, NEXT_REC)
	 end if
	    
         if (myPE.eq.PE_IO) then
            NUM_REC = NEXT_REC - NUM_REC 
            WRITE (uspx + L, REC=3) NEXT_REC, NUM_REC 
            if(unit_add == 0) CALL FLUSH_bin (uspx + L) 
         end if
!        call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
!
!
      CASE DEFAULT
            LINE(1) = 'Unknown SPx file index' 
            CALL WRITE_ERROR ('WRITE_SPX1', LINE, 1) 
            CALL MFIX_EXIT(myPE)
      END SELECT 

!//      call unlock_tmp_array
!
      deallocate (array1)    !//
      deallocate (array2)    !//
!
      RETURN  
      END SUBROUTINE WRITE_SPX1 
      
      subroutine gatherWriteSpx(VAR, array2, array1, uspxL, NEXT_REC)
        USE geometry
        USE compar           !//
        USE mpi_utility      !//d pnicol : for gatherWriteSpx
        USE sendrecv         !//d pnicol : for gatherWriteSpx
        IMPLICIT NONE
	integer uspxL, NEXT_REC
        double precision, dimension(ijkmax2) :: array1       
        double precision, dimension(ijkmax3) :: array2     
        double precision, dimension(DIMENSION_3) :: VAR    

!       call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
        call gather (VAR,array2,root)  !//d pnicol
!       call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
        if (myPE.eq.PE_IO) then
           call convert_to_io_dp(array2,array1,ijkmax2)  
           CALL OUT_BIN_R (uspxL, array1, IJKMAX2, NEXT_REC) 
        end if
      
      End subroutine gatherWriteSpx
      
