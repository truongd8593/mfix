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
      USE scalars
      USE rxns
      USE funits 
      USE output
      USE energy
      USE compar           !//
      USE mpi_utility      !//d pnicol : for gather
      USE sendrecv         !//d pnicol : for gather
!//d pnicol  ... not needed    USE tmp_array
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
!//d pnicol : allocate arrays for gather/convert_to_io_dp
      double precision, allocatable :: array1(:)
      double precision, allocatable :: array2(:)


!             loop counter
      INTEGER :: LC, N
!
!             pointer to first time-dependent record in restart file
      INTEGER :: NEXT_REC 
!-----------------------------------------------
!

!//d pnicol      
      if (myPE.eq.PE_IO) then
         allocate (array1(ijkmax2)) 
         allocate (array2(ijkmax3))  
      else
         allocate (array1(1)) 
         allocate (array2(1))  
      end if


      if (myPE.eq.PE_IO) then
         READ (UNIT_RES, REC=3) NEXT_REC 
         WRITE (UNIT_RES, REC=NEXT_REC) TIME, DT, NSTEP 
         NEXT_REC = NEXT_REC + 1 
      end if
!
!\\SP Local Send Receive - need to be moved to source later!!
      call send_recv(EP_g,2)
      call send_recv(P_g,2)
      call send_recv(P_star,2)
      call send_recv(RO_g,2)
      call send_recv(ROP_g,2)
      call send_recv(X_g,2)
      call send_recv(T_g,2)
      call send_recv(U_g,2)
      call send_recv(V_g,2)
      call send_recv(W_g,2)
      call send_recv(ROP_S,2)
      call send_recv(T_S,2)
      call send_recv(U_S,2)
      call send_recv(V_S,2)
      call send_recv(W_S,2)
      call send_recv(THETA_M,2)
      call send_recv(X_S,2)
      if(NScalar > 0)call send_recv(Scalar,2)
      if(K_Epsilon) THEN
        call send_recv(K_Turb_G,2)
        call send_recv(E_Turb_G,2)
      endif
      call send_recv(GAMA_RG,2)
      call send_recv(T_RG,2)
      call send_recv(GAMA_RS,2)
      call send_recv(T_RS,2)
      if(nRR > 0)call send_recv(ReactionRates,2)


      call gatherWriteRes (EP_g,array2, array1, NEXT_REC)  !//d pnicol
!
      call gatherWriteRes (P_g,array2, array1, NEXT_REC)  !//d pnicol
!
      call gatherWriteRes (P_star,array2, array1, NEXT_REC)  !//d pnicol
!
      call gatherWriteRes (RO_g,array2, array1, NEXT_REC)  !//d pnicol
!
      call gatherWriteRes (ROP_g,array2, array1, NEXT_REC)  !//d pnicol
!
      call gatherWriteRes (T_g,array2, array1, NEXT_REC)  !//d pnicol
!
      DO N = 1, NMAX(0) 
            call gatherWriteRes (X_g(:,n),array2, array1, NEXT_REC)  !//d pnicol
      END DO 
!
      call gatherWriteRes (U_g,array2, array1, NEXT_REC)  !//d pnicol
!
      call gatherWriteRes (V_g,array2, array1, NEXT_REC)  !//d pnicol
!
      call gatherWriteRes (W_g,array2, array1, NEXT_REC)  !//d pnicol
!
      DO LC = 1, MMAX 
!
        call gatherWriteRes (ROP_s(:,LC),array2, array1, NEXT_REC)  !//d pnicol
        
!
        call gatherWriteRes (T_s(:,LC),array2, array1, NEXT_REC)  !//d pnicol
!
        call gatherWriteRes (U_s(:,LC),array2, array1, NEXT_REC)  !//d pnicol
!
        call gatherWriteRes (V_s(:,LC),array2, array1, NEXT_REC)  !//d pnicol
!
        call gatherWriteRes (W_s(:,LC),array2, array1, NEXT_REC)  !//d pnicol
!
        call gatherWriteRes (THETA_M(:,LC),array2, array1, NEXT_REC)  !//d pnicol
!
         DO N = 1, NMAX(LC) 
            call gatherWriteRes (X_s(:,LC,N),array2, array1, NEXT_REC)  !//d pnicol
         END DO 
      END DO 
!
!     Version 1.3

      DO LC = 1, NScalar 
            call gatherWriteRes (Scalar(:,LC),array2, array1, NEXT_REC)  !//d pnicol
      END DO
!
!     Version 1.4 -- write radiation variables in write_res1 
      call gatherWriteRes (GAMA_RG,array2, array1, NEXT_REC)  !//d pnicol
 
      call gatherWriteRes (T_RG,array2, array1, NEXT_REC)  !//d pnicol

      DO LC = 1, MMAX 
        call gatherWriteRes (GAMA_RS(1,LC),array2, array1, NEXT_REC)  !//d pnicol

        call gatherWriteRes (T_RS(1,LC),array2, array1, NEXT_REC)  !//d pnicol
      ENDDO 

!
!     Version 1.5
      DO LC = 1, nRR 
            call gatherWriteRes (ReactionRates(:,LC),array2, array1, NEXT_REC)  !//d pnicol
      END DO 
!
!     Version 1.6

      if (K_epsilon) then
            call gatherWriteRes (K_turb_G,array2, array1, NEXT_REC)  !//d pnicol
	    call gatherWriteRes (E_turb_G,array2, array1, NEXT_REC)
      endif
!--------------------------------------------------------------------- 
 
      if (myPE.eq.PE_IO) CALL FLUSH_res (UNIT_RES) 
!     call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
!
      deallocate (array1)  !//d pnicol
      deallocate (array2)  !//d pnicol
!     call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
!
      RETURN  
      END SUBROUTINE WRITE_RES1 
      
      subroutine gatherWriteRes(VAR, array2, array1, NEXT_REC)
        USE geometry
        USE funits 
        USE compar           !//
        USE mpi_utility      !//d pnicol : for gather
        USE sendrecv         !//d pnicol : for gather
        IMPLICIT NONE
        double precision, dimension(ijkmax2) :: array1       
        double precision, dimension(ijkmax3) :: array2     
        double precision, dimension(DIMENSION_3) :: VAR    
        INTEGER :: NEXT_REC 

!       call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
        call gather (VAR,array2,root)  !//d pnicol
!       call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
        if (myPE.eq.PE_IO) then
           call convert_to_io_dp(array2,array1,ijkmax2)  
           CALL OUT_BIN_512 (UNIT_RES, array1, IJKMAX2, NEXT_REC) 
        end if
      
      End subroutine gatherWriteRes
      
