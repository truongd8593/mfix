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
      USE rxns
      USE scalars
      USE funits 
      USE energy
      USE compar     
      USE cdist 
      USE mpi_utility 
      USE sendrecv    

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

!//PAR_I/O declare global scratch arrays
      double precision, allocatable :: array1(:)
      double precision, allocatable :: array2(:)
      INTEGER allocstatus, i
!-----------------------------------------------
!
      if (myPE .eq. PE_IO .or. .not.bStart_with_one_res) then
         allocate (array1(ijkmax2))
         allocate (array2(ijkmax3))
      else
         allocate (array1(1))
         allocate (array2(1))
      end if

!// Reset global scratch arrays
      array1(:) = Undefined
      array2(:) = Undefined

!      call MPI_barrier(MPI_COMM_WORLD,mpierr)
!
!     Use DT from data file if DT_FAC is set to 1.0
      IF (DT_FAC == ONE) DT_SAVE = DT 
!

!//PAR_I/O only PE_IO reads the restart file
      if (myPE == PE_IO .or. (bDist_IO .and. .not.bStart_with_one_RES)) then
         READ (UNIT_RES, REC=1) VERSION 
         READ (VERSION(6:512), *) VERSION_NUMBER 

         READ (UNIT_RES, REC=3) NEXT_REC 
         IF (VERSION_NUMBER >= 1.12) THEN 
            READ (UNIT_RES, REC=NEXT_REC) TIME, DT, NSTEP 
         ELSE 
            READ (UNIT_RES, REC=NEXT_REC) TIME, NSTEP 
         ENDIF 
         NEXT_REC = NEXT_REC + 1 
      end if
!
      if (.not.bDist_IO  .or. bStart_with_one_RES) then
!        call MPI_barrier(MPI_COMM_WORLD,mpierr)
      
         call bcast(VERSION, PE_IO)        !//PAR_I/O BCAST0c
         call bcast(VERSION_NUMBER, PE_IO) !//PAR_I/O BCAST0r
         call bcast(TIME, PE_IO)           !//PAR_I/O BCAST0d
         call bcast(NSTEP, PE_IO)          !//PAR_I/O BCAST0i
         if (VERSION_NUMBER >= 1.12) call bcast(DT, PE_IO)   !//PAR_I/O BCAST0d	  
	end if   
!      call MPI_barrier(MPI_COMM_WORLD,mpierr)

!AE TIME 091501 Store the timestep counter level at the begin of RESTART run
      NSTEPRST = NSTEP

!
      call readScatterRes(EP_G,array2, array1, NEXT_REC)

      call readScatterRes(P_G,array2, array1, NEXT_REC)

      call readScatterRes(P_STAR,array2, array1, NEXT_REC)

      call readScatterRes(RO_G,array2, array1, NEXT_REC)

      call readScatterRes(ROP_G,array2, array1, NEXT_REC)

      call readScatterRes(T_G,array2, array1, NEXT_REC)
!

      IF (VERSION_NUMBER < 1.15) THEN 
         call readScatterRes (T_s(:,1),array2, array1, NEXT_REC)

         IF (MMAX >= 2) THEN 
            call readScatterRes (T_s(:,2),array2, array1, NEXT_REC)
         ELSE 
            if (myPE == PE_IO) &
                   CALL IN_BIN_512 (UNIT_RES, array1, IJKMAX2, NEXT_REC) 
         ENDIF 
      ENDIF 

      IF (VERSION_NUMBER >= 1.05) THEN 
         DO N = 1, NMAX(0) 
            call readScatterRes (X_g(:,n),array2, array1, NEXT_REC)
         END DO 
      ENDIF 

      call readScatterRes(U_G, array2, array1, NEXT_REC)
      call readScatterRes(V_G, array2, array1, NEXT_REC)
      call readScatterRes(W_G, array2, array1, NEXT_REC)
!
      DO LC = 1, MMAX 
         call readScatterRes(ROP_S(:,LC), array2, array1, NEXT_REC)

         IF (VERSION_NUMBER >= 1.15) THEN
            call readScatterRes(T_S(:,LC), array2, array1, NEXT_REC)
         END IF
         call readScatterRes(U_S(:,LC), array2, array1, NEXT_REC)
         call readScatterRes(V_S(:,LC), array2, array1, NEXT_REC)
         call readScatterRes(W_S(:,LC), array2, array1, NEXT_REC)

         IF (VERSION_NUMBER >= 1.2) then
            call readScatterRes(THETA_M(:,LC), array2, array1, NEXT_REC)
         end if
         IF (VERSION_NUMBER >= 1.05) THEN 
            DO N = 1, NMAX(LC) 
               call readScatterRes(X_S(:,LC,N), array2, array1, NEXT_REC)
            END DO 
         ENDIF 
      END DO 

      IF (VERSION_NUMBER >= 1.3) THEN
        DO N = 1, NScalar 
          call readScatterRes(Scalar(:,N), array2, array1, NEXT_REC)
        END DO 
      ENDIF

      IF (VERSION_NUMBER >= 1.4) THEN
        call readScatterRes(GAMA_RG, array2, array1, NEXT_REC)
 
        call readScatterRes(T_RG, array2, array1, NEXT_REC)

        DO LC = 1, MMAX 
          call readScatterRes(GAMA_RS(1,LC), array2, array1, NEXT_REC)

          call readScatterRes(T_RS(1,LC), array2, array1, NEXT_REC)

        ENDDO 
      ELSE
        GAMA_RG(:)   = ZERO
        T_RG (:)     = ZERO
        GAMA_RS(:,:) = ZERO
        T_RS(:,:)    = ZERO
      ENDIF


      IF (VERSION_NUMBER >= 1.5) THEN
        DO N = 1, nRR 
          call readScatterRes(ReactionRates(:,N), array2, array1, NEXT_REC)
        END DO 
      ENDIF 

      IF (VERSION_NUMBER >= 1.6 .AND. K_Epsilon) THEN 
          call readScatterRes(K_Turb_G, array2, array1, NEXT_REC)
          call readScatterRes(E_Turb_G, array2, array1, NEXT_REC)
      ENDIF
!------------------------------------------------------------------------
!


!      call MPI_barrier(MPI_COMM_WORLD,mpierr)
      deallocate( array1 )
      deallocate( array2 )
!      call MPI_barrier(MPI_COMM_WORLD,mpierr)

      if (.not.bDist_IO .or. bStart_with_one_RES) then
         call send_recv(rop_g)
         call send_recv(ro_g)
         call send_recv(rop_s)
      end if


      IF (DT_FAC == ONE) DT = DT_SAVE 
!

      RETURN  
      END SUBROUTINE READ_RES1 
      
      subroutine readScatterRes(VAR, array2, array1, NEXT_REC)
        USE geometry
        USE funits 
        USE compar      
	  USE cdist     
        USE mpi_utility      
        USE sendrecv         
        IMPLICIT NONE
        double precision, dimension(ijkmax2) :: array1       
        double precision, dimension(ijkmax3) :: array2     
        double precision, dimension(DIMENSION_3) :: VAR
        INTEGER :: NEXT_REC 
	   
	if (.not.bDist_IO .or. bStart_with_one_RES) then 
         if (myPE == PE_IO) then
            CALL IN_BIN_512 (UNIT_RES, array1, IJKMAX2, NEXT_REC) 
            CALL convert_from_io_dp(array1, array2, IJKMAX2) 
         end if
!        call MPI_barrier(MPI_COMM_WORLD,mpierr)
         call scatter(VAR, array2, PE_IO)
!        call MPI_barrier(MPI_COMM_WORLD,mpierr)
      else
         CALL IN_BIN_512 (UNIT_RES, var, size(var) , NEXT_REC) 
	end if
      
      End subroutine readScatterRes

!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization
!// 020 New local variables for parallelization: array1, array2
!// 400 Added sendrecv module and send_recv calls for COMMunication
!// 400 Added mpi_utility module and other global reduction (bcast) calls
