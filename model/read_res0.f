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
      USE compar      !// 001 Include MPI header file
      USE mpi_utility !//
!     USE dbg_util !//PARDBG


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

!//PAR_I/O 0814 declare global scratch arrays
      INTEGER, ALLOCATABLE, DIMENSION(:) :: iGTEMP,iGTEMP2    !//PAR_I/O declare integer Global SCRatch array
      REAL, ALLOCATABLE, DIMENSION(:) :: rGTEMP    !//PAR_I/O declare real*4 Global SCRatch array
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: dGTEMP  !//PAR_I/O declare real*8 Global SCRatch array
!//PAR_I/O 0815 declaring following arrays to pack scalar variables when BCASTing to reduce the number of BCAST calls
      INTEGER, ALLOCATABLE, DIMENSION(:) :: INTPACK           !//PAR_I/O packing array for integers
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: DBLPACK  !//PAR_I/O packing array for doubles    

!//? 0906 temporarily setting KMAX3 here, need to find a better place to initialize it
!      KMAX3 = KMAX2 + 1
      KMAX3 = KMAX2
!-----------------------------------------------

!      call lock_tmp_array
!      write(*,"('(PE ',I2,'): begin READ_RES0')") myPE  !//AIKEPARDBG

!
!
!  1) Check to ensure that this subroutine was updated.
!  2) Initialize missing constants from earlier versions.
!  3) Add new read statements at the end of the file.
!
!//PAR_I/O 0813 only PE_IO reads the restart file
    if (myPE == PE_IO ) then
      READ (UNIT_RES, REC=1) VERSION 
      READ (VERSION(6:512), *) VERSION_NUMBER 
      IF (VERSION_NUMBER > 1.2) THEN 
         WRITE (*, *) ' Update Subroutine read_res0' 
         CALL SLUMBER 
!         STOP  
         call exitMPI(myPE)  !// 990 0807 Abort all PEs, not only the current one
      ENDIF 
    endif

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
!//PAR_I/O 0813 !//PAR_I/O only PE_IO reads the restart
    if (myPE == PE_IO ) then
      READ (UNIT_RES, REC=2) RUN_NAME, ID_MONTH, ID_DAY, ID_YEAR, ID_HOUR, &
         ID_MINUTE, ID_SECOND 
      READ (UNIT_RES, REC=3) NEXT_RECA 

!      write(*,"('(PE ',I2,'): M2B31 ',A20)") myPE,version   !//AIKEPARDBG   
!      write(*,"('(PE ',I2,'): M2B33 ',I7)") myPE, NEXT_RECA !//AIKEPARDBG   
      
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
    endif
!d ver 01.01 & 01.02 adds DIM_IC , DIM_BC in to .RES header
!d ver 01.03         adds XMIN in addition to above
!d ver 01.04         adds DIM_C in addition to above
!d ver 01.05         adds DIM_IS in addition to above
!d ver>01.05         adds C_E, C_F, PHI, PHI_W

    Allocate( INTPACK(30))  !//PAR_I/O ALLOCate packing array
    Allocate( DBLPACK(20))  !//PAR_I/O ALLOCate packing array
    INTPACK = 0
    DBLPACK = 0.d0
!//PAR_I/O 0814 Root PE (PE_IO) : Pack and broadcast ; others unpack
!//WARNING: This implementation assumes VERSION > 1.05, filtering needed for earlier versions
    if (myPE == PE_IO ) then
       call bcast(VERSION_NUMBER, PE_IO)  !//PAR_I/O BCAST0r
       INTPACK(1) = ID_MONTH
       INTPACK(2) = ID_DAY
       INTPACK(3) = ID_YEAR
       INTPACK(4) = ID_HOUR
       INTPACK(5) = ID_MINUTE
       INTPACK(6) = ID_SECOND
       INTPACK(7) = IMIN1       
       INTPACK(8) = JMIN1
       INTPACK(9) = KMIN1
       INTPACK(10) = IMAX
       INTPACK(11) = JMAX
       INTPACK(12) = KMAX
       INTPACK(13) = IMAX1
       INTPACK(14) = JMAX1
       INTPACK(15) = KMAX1
       INTPACK(16) = IMAX2
       INTPACK(17) = JMAX2
       INTPACK(18) = KMAX2
       INTPACK(19) = IJMAX2
       INTPACK(20) = IJKMAX2
       INTPACK(21) = MMAX
       INTPACK(22) = DIM_IC
       INTPACK(23) = DIM_BC
       INTPACK(24) = DIM_C
       INTPACK(25) = DIM_IS
!//S In spite of the overhead, using MPI_Pack/MPI_Unpack may be worth using here
       call bcast(INTPACK(1:25),PE_IO)  !//PAR_I/O BCAST1i
       DBLPACK(1) = DT
       DBLPACK(2) = XMIN
       DBLPACK(3) = XLENGTH
       DBLPACK(4) = YLENGTH
       DBLPACK(5) = ZLENGTH
       DBLPACK(6) = C_E
       DBLPACK(7) = C_F
       DBLPACK(8) = PHI
       DBLPACK(9) = PHI_W 
       call bcast(DBLPACK,PE_IO) !//PAR_I/O BCAST1d
       call bcast(VERSION,PE_IO) !//PAR_I/O BCAST0c
       call bcast(RUN_NAME,PE_IO) !//PAR_I/O BCAST0c
!//S need a count variable as arg in bcast so that we can arrange the packing based
!    on version number and send with a fixed number of arrays
    else
       call bcast(VERSION_NUMBER, PE_IO)  !//PAR_I/O BCAST0r    
       call bcast(INTPACK(1:25),PE_IO)    !//PAR_I/O BCAST1i (receive)
       ID_MONTH = INTPACK(1)
       ID_DAY = INTPACK(2)
       ID_YEAR = INTPACK(3) 
       ID_HOUR = INTPACK(4)
       ID_MINUTE = INTPACK(5)
       ID_SECOND = INTPACK(6)
       IMIN1 = INTPACK(7)
       JMIN1 = INTPACK(8)
       KMIN1 = INTPACK(9)
       IMAX = INTPACK(10)
       JMAX = INTPACK(11)
       KMAX = INTPACK(12)
       IMAX1 = INTPACK(13)
       JMAX1 = INTPACK(14)
       KMAX1 = INTPACK(15)
       IMAX2 = INTPACK(16)
       JMAX2 = INTPACK(17)
       KMAX2 = INTPACK(18)
       IJMAX2 = INTPACK(19)
       IJKMAX2 = INTPACK(20)
       MMAX = INTPACK(21)
       DIM_IC = INTPACK(22)
       DIM_BC = INTPACK(23)
       DIM_C = INTPACK(24)
       DIM_IS = INTPACK(25)
       call bcast(DBLPACK,PE_IO)  !//PAR_I/O BCAST1d (recv)
       DT = DBLPACK(1) 
       XMIN = DBLPACK(2)
       XLENGTH = DBLPACK(3)
       YLENGTH = DBLPACK(4)
       ZLENGTH = DBLPACK(5)
       C_E = DBLPACK(6)
       C_F = DBLPACK(7)
       PHI = DBLPACK(8)
       PHI_W = DBLPACK(9)
       call bcast(VERSION,PE_IO)  !//PAR_I/O BCAST0c (recv)
       call bcast(RUN_NAME,PE_IO) !//PAR_I/O BCAST0c (recv)
    endif

!//AIKEPARDBG 0816 DEBUG STOP
!    call dbgprn(INTPACK,25)   !//AIKEPARDBG
!    call dbgprn(DBLPACK,9)    !//AIKEPARDBG
!     write(UNIT_LOG,*) 'INTPACK',INTPACK   !//AIKEPARDBG
!     write(UNIT_LOG,*) 'DBLPACK',DBLPACK   !//AIKEPARDBG  
!    write(*,"('(PE ',I2,'): VERSION=',A80)") myPE,VERSION !//AIKEPARDBG
!    write(*,"('(PE ',I2,'): aft Bcast of INT in read_res0')") myPE !//AIKEPARDBG
!    call exitMPI(myPE) !//AIKEPARDBG
!//PAR_I/O deaallocate the packing arrays used to pack scalars for reducing Bcast calls
    DeAllocate (INTPACK)

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

!//AIKEPARDBG Block1
                                if (myPE == PE_IO) then

                                  IF (VERSION_NUMBER >= 1.04) THEN 

                                    CALL IN_BIN_512 (UNIT_RES, C, DIM_C, &
                                       NEXT_RECA) 
                                    call bcast(C,PE_IO)   !//PAR_I/O BCAST1d user defined constants,C
!                                                ! work around for -O3 compiler bug
                                    NEXT_RECA = 1 + NEXT_RECA 
                                    NEXT_RECA = NEXT_RECA - 1 
                                    DO LC = 1, DIM_C 
                                    READ (UNIT_RES, REC=NEXT_RECA) C_NAME(LC) 
                                    NEXT_RECA = NEXT_RECA + 1 
                                    END DO 
                                    call bcast(C_NAME,PE_IO)  !//PAR_I/O BCAST1c user defined constant names                                    
                                    call bcast(MMAX,PE_IO)    !//PAR_I/O BCAST0i # of solid phases
                                    IF (VERSION_NUMBER < 1.12) THEN 
                                      CALL IN_BIN_512I (UNIT_RES, NMAX, MMAX + 1&
                                       , NEXT_RECA) 
                                    ELSE 
                                      READ (UNIT_RES, REC=NEXT_RECA) (NMAX(L),L=0&
                                       ,MMAX) 
                                      NEXT_RECA = NEXT_RECA + 1 
                                    ENDIF 
                                    call bcast(NMAX,PE_IO)   !//PAR_I/O BCAST1i total # of gas OR solid species
                                  ENDIF

                                else    ! else of if (myPE == PE_IO)

                                  IF (VERSION_NUMBER >= 1.04) THEN 

                                    call bcast(C,PE_IO)      !//PAR_I/O BCAST1d (recv)   
                                    call bcast(C_NAME,PE_IO) !//PAR_I/O BCAST1c (recv)
                                    call bcast(MMAX,PE_IO)   !//PAR_I/O BCAST0i (recv)
                                    call bcast(NMAX,PE_IO)   !//PAR_I/O BCAST1i (recv)
                                  ENDIF
                                endif   ! end of if (myPE == PE_IO)

!//AIKEPARDBG 0906 DEBUGSTOP
!    write(*,"('(PE ',I2,'): MMAX:',I4)") myPE,MMAX !//AIKEPARDBG
!    write(*,"('(PE ',I2,'): aft Block1 in read_res0')") myPE !//AIKEPARDBG
!    call exitMPI(myPE) !//AIKEPARDBG


                                 IF (NMAX(0) <= DIMENSION_N_G) THEN 
                                    DO M = 1, MMAX 
                                      IF (NMAX(M) > DIMENSION_N_S) GO TO 900 
                                    END DO 

                                   if (myPE == PE_IO) then    !//PAR_I/O only PE_IO reads RES file
                                      CALL IN_BIN_512 (UNIT_RES, DX, IMAX2, &
                                        NEXT_RECA) 
                                      call bcast(dx, PE_IO)   !//PAR_I/O BCAST1d

                                      CALL IN_BIN_512 (UNIT_RES, DY, JMAX2, &
                                        NEXT_RECA) 
                                      call bcast(dy, PE_IO)   !//PAR_I/O BCAST1d
!//AIKEPARDBGSTOP 0907
    write(*,"('(PE ',I2,'): aft dy in Block2 in read_res0')") myPE !//AIKEPARDBG
    write(UNIT_LOG,*) dy    !//AIKEPARDBG
    call exitMPI(myPE)      !//AIKEPARDBG

!//? 0906 kmax2 or kmax3 is the index for the global dummy array? 
                                      Allocate( dGTEMP(KMAX3)) !//PAR_I/O ALLOCate R*8 Global scratch
                                      CALL IN_BIN_512 (UNIT_RES, dGTEMP, KMAX3, &
                                        NEXT_RECA) 
                                      call scatter(dz,dGTEMP,PE_IO)   !//PAR_I/O SCATTER1d
                                      DeAllocate( dGTEMP)
				      
!    call MPI_Barrier(MPI_COMM_WORLD,mpierr)	
    write(*,"('(PE ',I2,'): KMAX3',I2)") myPE, KMAX3 !//AIKEPARDBG
    write(UNIT_LOG,"('(PE ',I2,'): aft dz in Block2 in read_res0')") myPE !//AIKEPARDBG
    write(UNIT_LOG,*) dz    !//AIKEPARDBG
    call exitMPI(myPE)      !//AIKEPARDBG

                                      READ (UNIT_RES, REC=NEXT_RECA) RUN_NAME, &
                                         DESCRIPTION, UNITS, RUN_TYPE, &
                                         COORDINATES 
                                      NEXT_RECA = NEXT_RECA + 1
                                      call bcast(RUN_NAME,PE_IO)    !//PAR_I/O BCAST0c
                                      call bcast(DESCRIPTION,PE_IO) !//PAR_I/O BCAST0c
                                      call bcast(UNITS,PE_IO)       !//PAR_I/O BCAST0c
                                      call bcast(RUN_TYPE,PE_IO)    !//PAR_I/O BCAST0c
                                      call bcast(COORDINATES,PE_IO) !//PAR_I/O BCAST0c

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

                                      call bcast(D_P, PE_IO)       !//PAR_I/O BCAST1d
                                      call bcast(RO_S, PE_IO)      !//PAR_I/O BCAST1d
                                      call bcast(EP_STAR, PE_IO)   !//PAR_I/O BCAST0d
                                      if (VERSION /= 'RES = 01.00' .OR. VERSION /= &
                                         'RES = 01.01') &
                                      call bcast(RO_G0, PE_IO)     !//PAR_I/O BCAST0d
                                      call bcast(MU_G0, PE_IO)     !//PAR_I/O BCAST0d
                                      call bcast(MW_AVG, PE_IO)    !//PAR_I/O BCAST0d


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
!//? adding ncount var'ble to bcast may be better instead of bcasts with full array size
                                        call bcast(MW_G, PE_IO)     !//PAR_I/O BCAST1d
                                        call bcast(MW_S, PE_IO)     !//PAR_I/O BCAST2d
                                      ENDIF 

                                      CALL IN_BIN_512 (UNIT_RES, IC_X_W, DIM_IC, &
                                         NEXT_RECA)
                                      CALL IN_BIN_512 (UNIT_RES, IC_X_E, DIM_IC, &
                                         NEXT_RECA) 
                                      CALL IN_BIN_512 (UNIT_RES, IC_Y_S, DIM_IC, &
                                         NEXT_RECA) 
                                      CALL IN_BIN_512 (UNIT_RES, IC_Y_N, DIM_IC, &
                                         NEXT_RECA) 
                                      call bcast(IC_X_W, PE_IO)     !//PAR_I/O BCAST1d
                                      call bcast(IC_X_E, PE_IO)     !//PAR_I/O BCAST1d
                                      call bcast(IC_Y_S, PE_IO)     !//PAR_I/O BCAST1d
                                      call bcast(IC_Y_N, PE_IO)     !//PAR_I/O BCAST1d
                                      CALL IN_BIN_512 (UNIT_RES, IC_Z_B, DIM_IC, &
                                         NEXT_RECA) 
                                      CALL IN_BIN_512 (UNIT_RES, IC_Z_T, DIM_IC, &
                                         NEXT_RECA) 
!//? do we need to bcast to all or only IC_Z_T to last PE, does the interior PE need this info
                                      call bcast(IC_Z_B, PE_IO)     !//PAR_I/O BCAST1d
                                      call bcast(IC_Z_T, PE_IO)     !//PAR_I/O BCAST1d

                                      CALL IN_BIN_512I (UNIT_RES, IC_I_W, DIM_IC&
                                         , NEXT_RECA) 
                                      CALL IN_BIN_512I (UNIT_RES, IC_I_E, DIM_IC&
                                         , NEXT_RECA) 
                                      CALL IN_BIN_512I (UNIT_RES, IC_J_S, DIM_IC&
                                         , NEXT_RECA) 
                                      CALL IN_BIN_512I (UNIT_RES, IC_J_N, DIM_IC&
                                         , NEXT_RECA)
                                      call bcast(IC_I_W, PE_IO)     !//PAR_I/O BCAST1i
                                      call bcast(IC_I_E, PE_IO)     !//PAR_I/O BCAST1i
                                      call bcast(IC_J_S, PE_IO)     !//PAR_I/O BCAST1i
                                      call bcast(IC_J_N, PE_IO)     !//PAR_I/O BCAST1i

                                      CALL IN_BIN_512I (UNIT_RES, IC_K_B, DIM_IC&
                                         , NEXT_RECA) 
                                      CALL IN_BIN_512I (UNIT_RES, IC_K_T, DIM_IC&
                                         , NEXT_RECA) 
!//? do we need to bcast to all or only IC_Z_T to last PE, does the interior PE need this info
                                      call bcast(IC_K_B, PE_IO)     !//PAR_I/O BCAST1i
                                      call bcast(IC_K_T, PE_IO)     !//PAR_I/O BCAST1i

                                      CALL IN_BIN_512 (UNIT_RES, IC_EP_G, DIM_IC&
                                         , NEXT_RECA) 
                                      CALL IN_BIN_512 (UNIT_RES, IC_P_G, DIM_IC, &
                                         NEXT_RECA) 
                                      CALL IN_BIN_512 (UNIT_RES, IC_T_G, DIM_IC, &
                                         NEXT_RECA) 
                                      call bcast(IC_EP_G, PE_IO)    !//PAR_I/O BCAST1d
                                      call bcast(IC_P_G, PE_IO)     !//PAR_I/O BCAST1d
                                      call bcast(IC_T_G, PE_IO)     !//PAR_I/O BCAST1d

                                      IF (VERSION_NUMBER < 1.15) THEN 
                                        CALL IN_BIN_512 (UNIT_RES, IC_T_S(1,1), &
                                        DIM_IC, NEXT_RECA) 
                                        IF (MMAX >= 2) THEN 
                                          CALL IN_BIN_512 (UNIT_RES, IC_T_S(1,2), &
                                           DIM_IC, NEXT_RECA) 
                                        ELSE 
                                          CALL IN_BIN_512 (UNIT_RES, IC_TMP, DIM_IC, &
                                          NEXT_RECA) 
                                          call bcast(IC_TMP, PE_IO)    !//PAR_I/O BCAST1d
                                        ENDIF 
                                        call bcast(IC_T_S, PE_IO)    !//PAR_I/O BCAST2d

                                      ENDIF 

                                      IF (VERSION_NUMBER >= 1.04) THEN 
                                        DO N = 1, NMAX(0) 
                                          CALL IN_BIN_512 (UNIT_RES, IC_X_G(1,N), &
                                          DIM_IC, NEXT_RECA) 
                                        END DO 
                                        call bcast(IC_X_G, PE_IO)    !//PAR_I/O BCAST2d
                                      ENDIF 
                                      CALL IN_BIN_512 (UNIT_RES, IC_U_G, DIM_IC, &
                                         NEXT_RECA) 
                                      CALL IN_BIN_512 (UNIT_RES, IC_V_G, DIM_IC, &
                                         NEXT_RECA) 
                                      CALL IN_BIN_512 (UNIT_RES, IC_W_G, DIM_IC, &
                                         NEXT_RECA) 
                                      call bcast(IC_U_G, PE_IO)    !//PAR_I/O BCAST1d
                                      call bcast(IC_V_G, PE_IO)    !//PAR_I/O BCAST1d
                                      call bcast(IC_W_G, PE_IO)    !//PAR_I/O BCAST1d

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
                                      call bcast(IC_ROP_S, PE_IO)  !//PAR_I/O BCAST2d
                                      call bcast(IC_U_S, PE_IO)    !//PAR_I/O BCAST2d
                                      call bcast(IC_V_S, PE_IO)    !//PAR_I/O BCAST2d
                                      call bcast(IC_W_S, PE_IO)    !//PAR_I/O BCAST2d
                                      if (VERSION_NUMBER >= 1.15) &
                                        call bcast(IC_T_S, PE_IO)  !//PAR_I/O BCAST2d
                                      if (VERSION_NUMBER >= 1.04) &
                                        call bcast(IC_X_S, PE_IO)  !//PAR_I/O BCAST3d

                                    CALL IN_BIN_512 (UNIT_RES, BC_X_W, DIM_BC, &
                                       NEXT_RECA) 
                                    CALL IN_BIN_512 (UNIT_RES, BC_X_E, DIM_BC, &
                                       NEXT_RECA) 
                                    CALL IN_BIN_512 (UNIT_RES, BC_Y_S, DIM_BC, &
                                       NEXT_RECA) 
                                    CALL IN_BIN_512 (UNIT_RES, BC_Y_N, DIM_BC, &
                                       NEXT_RECA) 

                                    call bcast(BC_X_W, PE_IO)    !//PAR_I/O BCAST1d
                                    call bcast(BC_X_E, PE_IO)    !//PAR_I/O BCAST1d
                                    call bcast(BC_Y_S, PE_IO)    !//PAR_I/O BCAST1d
                                    call bcast(BC_Y_N, PE_IO)    !//PAR_I/O BCAST1d

                                    CALL IN_BIN_512 (UNIT_RES, BC_Z_B, DIM_BC, &
                                       NEXT_RECA) 
                                    CALL IN_BIN_512 (UNIT_RES, BC_Z_T, DIM_BC, &
                                       NEXT_RECA) 
!//? do we need to bcast to all or only BC_Z_T to last PE, does the interior PE need this info
                                    call bcast(BC_Z_B, PE_IO)    !//PAR_I/O BCAST1d
                                    call bcast(BC_Z_T, PE_IO)    !//PAR_I/O BCAST1d

                                    CALL IN_BIN_512I (UNIT_RES, BC_I_W, DIM_BC&
                                       , NEXT_RECA) 
                                    CALL IN_BIN_512I (UNIT_RES, BC_I_E, DIM_BC&
                                       , NEXT_RECA) 
                                    CALL IN_BIN_512I (UNIT_RES, BC_J_S, DIM_BC&
                                       , NEXT_RECA) 
                                    CALL IN_BIN_512I (UNIT_RES, BC_J_N, DIM_BC&
                                       , NEXT_RECA) 

                                    call bcast(BC_I_W, PE_IO)    !//PAR_I/O BCAST1i
                                    call bcast(BC_I_E, PE_IO)    !//PAR_I/O BCAST1i
                                    call bcast(BC_J_S, PE_IO)    !//PAR_I/O BCAST1i
                                    call bcast(BC_J_N, PE_IO)    !//PAR_I/O BCAST1i

                                    CALL IN_BIN_512I (UNIT_RES, BC_K_B, DIM_BC&
                                       , NEXT_RECA) 
                                    CALL IN_BIN_512I (UNIT_RES, BC_K_T, DIM_BC&
                                       , NEXT_RECA) 
!//? do we need to bcast to all or only BC_K_T to last PE, does the interior PE need this info
                                    call bcast(BC_K_B, PE_IO)    !//PAR_I/O BCAST1i
                                    call bcast(BC_K_T, PE_IO)    !//PAR_I/O BCAST1i

                                    CALL IN_BIN_512 (UNIT_RES, BC_EP_G, DIM_BC&
                                       , NEXT_RECA) 
                                    CALL IN_BIN_512 (UNIT_RES, BC_P_G, DIM_BC, &
                                       NEXT_RECA) 
                                    CALL IN_BIN_512 (UNIT_RES, BC_T_G, DIM_BC, &
                                       NEXT_RECA) 

                                    call bcast(BC_EP_G, PE_IO)   !//PAR_I/O BCAST1d
                                    call bcast(BC_P_G, PE_IO)    !//PAR_I/O BCAST1d
                                    call bcast(BC_T_G, PE_IO)    !//PAR_I/O BCAST1d

                                    IF (VERSION_NUMBER < 1.15) THEN 
                                      CALL IN_BIN_512 (UNIT_RES, BC_T_S(1,1), &
                                       DIM_BC, NEXT_RECA) 
                                      IF (MMAX >= 2) THEN 
                                        CALL IN_BIN_512 (UNIT_RES, BC_T_S(1,2), &
                                             DIM_BC, NEXT_RECA) 
                                      ELSE 
                                        CALL IN_BIN_512 (UNIT_RES, BC_TMP, DIM_BC, &
                                             NEXT_RECA) 
                                        call bcast(BC_TMP, PE_IO)    !//PAR_I/O BCAST1d
                                      ENDIF 
                                      call bcast(BC_T_S, PE_IO)   !//PAR_I/O BCAST2d
                                    ENDIF 
!


                                    IF (VERSION_NUMBER >= 1.04) THEN 
                                      DO N = 1, NMAX(0) 
                                        CALL IN_BIN_512 (UNIT_RES, BC_X_G(1,N), &
                                             DIM_BC, NEXT_RECA) 
                                      END DO 
                                      call bcast(BC_X_G, PE_IO)   !//PAR_I/O BCAST2d
                                    ENDIF 

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

                                    call bcast(BC_U_G, PE_IO)       !//PAR_I/O BCAST1d
                                    call bcast(BC_V_G, PE_IO)       !//PAR_I/O BCAST1d
                                    call bcast(BC_W_G, PE_IO)       !//PAR_I/O BCAST1d
                                    call bcast(BC_RO_G, PE_IO)      !//PAR_I/O BCAST1d
                                    call bcast(BC_ROP_G, PE_IO)     !//PAR_I/O BCAST1d
                                    call bcast(BC_VOLFLOW_G, PE_IO) !//PAR_I/O BCAST1d
                                    call bcast(BC_MASSFLOW_G, PE_IO)!//PAR_I/O BCAST1d

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

                                    call bcast(BC_ROP_S, PE_IO) !//PAR_I/O BCAST2d
                                    call bcast(BC_U_S, PE_IO)   !//PAR_I/O BCAST2d
                                    call bcast(BC_V_S, PE_IO)   !//PAR_I/O BCAST2d

                                    if (VERSION_NUMBER >= 1.04) then
                                       call bcast(BC_W_S, PE_IO)   !//PAR_I/O BCAST2d
                                       if (VERSION_NUMBER >= 1.15)  &
                                          call bcast(BC_T_S, PE_IO)   !//PAR_I/O BCAST2d
!
                                       call bcast(BC_X_S, PE_IO)   !//PAR_I/O BCAST2d
                                    endif 
                                    call bcast(BC_VOLFLOW_S, PE_IO)   !//PAR_I/O BCAST2d
                                    call bcast(BC_MASSFLOW_S, PE_IO)  !//PAR_I/O BCAST2d

                                    IF (VERSION == 'RES = 01.00') THEN 
                                      L = 10 
                                    ELSE 
                                      L = DIM_BC 
                                    ENDIF 
                                    DO LC = 1, L 
                                      READ (UNIT_RES, REC=NEXT_RECA) BC_TYPE(LC) 
                                      NEXT_RECA = NEXT_RECA + 1 
                                    END DO 

                                    call bcast(BC_TYPE, PE_IO)   !//PAR_I/O BCAST1c

!//? is ijkmax2 or ijkmax3 the global index? if not change followings
                                    Allocate(iGTEMP(IJKMAX2))   !//PAR_I/O ALLOCate INT Global scratch
                                    Allocate(iGTEMP2(IJKMAX2))   !//PAR_I/O ALLOCate INT Global scratch

                                    CALL IN_BIN_512I (UNIT_RES, iGTEMP, IJKMAX2, &
                                       NEXT_RECA) 
                                    call convert_from_io_i(iGTEMP,iGTEMP2,ijkmax2)

                                    call scatter(array1i,iGTEMP2,PE_IO)        !//PAR_I/O SCATTER1d
                                    DeAllocate (iGTEMP, iGTEMP2)


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

                                      call bcast(IS_X_W, PE_IO)   !//PAR_I/O BCAST1d
                                      call bcast(IS_X_E, PE_IO)   !//PAR_I/O BCAST1d
                                      call bcast(IS_Y_S, PE_IO)   !//PAR_I/O BCAST1d
                                      call bcast(IS_Y_N, PE_IO)   !//PAR_I/O BCAST1d
!//? do we need to bcast to all or only IS_Z_T to last PE, does the interior PE need this info
                                      call bcast(IS_Z_B, PE_IO)   !//PAR_I/O BCAST1d
                                      call bcast(IS_Z_T, PE_IO)   !//PAR_I/O BCAST1d

                                      call bcast(IS_I_W, PE_IO)   !//PAR_I/O BCAST1i
                                      call bcast(IS_I_E, PE_IO)   !//PAR_I/O BCAST1i
                                      call bcast(IS_J_S, PE_IO)   !//PAR_I/O BCAST1i
                                      call bcast(IS_J_N, PE_IO)   !//PAR_I/O BCAST1i


                                      IF (VERSION_NUMBER >= 1.07) THEN 
                                        DO LC = 1, MMAX 
                                          CALL IN_BIN_512 (UNIT_RES, IS_VEL_S(1,LC), &
                                                DIM_IS, NEXT_RECA) 
                                        END DO 
                                        call bcast(IS_VEL_S, PE_IO)   !//PAR_I/O BCAST2d
                                      ENDIF 
                                      DO LC = 1, DIM_IS 
                                        READ (UNIT_RES, REC=NEXT_RECA) IS_TYPE(LC) 
                                        NEXT_RECA = NEXT_RECA + 1 
                                      END DO 
                                      call bcast(IS_TYPE, PE_IO)   !//PAR_I/O BCAST1c
                                    ENDIF 
!
!
!  Additions from new versions of .RES file
!
                                    IF (VERSION_NUMBER >= 1.08) THEN 
                                      READ (UNIT_RES, REC=NEXT_RECA) CYCLIC_X, &
                                        CYCLIC_Y, CYCLIC_Z, CYCLIC_X_PD, &
                                        CYCLIC_Y_PD, CYCLIC_Z_PD, DELP_X, DELP_Y&
                                        , DELP_Z, U_G0, U_S0, V_G0, V_S0, W_G0, &
                                        W_S0 
                                      call bcast(CYCLIC_X,PE_IO)     !//PAR_I/O BCAST0l
                                      call bcast(CYCLIC_Y,PE_IO)     !//PAR_I/O BCAST0l
                                      call bcast(CYCLIC_Z,PE_IO)     !//PAR_I/O BCAST0l
                                      call bcast(CYCLIC_X_PD,PE_IO)  !//PAR_I/O BCAST0l
                                      call bcast(CYCLIC_Y_PD,PE_IO)  !//PAR_I/O BCAST0l
                                      call bcast(CYCLIC_Z_PD,PE_IO)  !//PAR_I/O BCAST0l
                                      
                                      DBLPACK(1) = DELP_X
                                      DBLPACK(2) = DELP_Y
                                      DBLPACK(3) = DELP_Z
                                      DBLPACK(4) = U_G0
                                      DBLPACK(5) = V_G0
                                      DBLPACK(6) = W_G0
                                      call bcast(DBLPACK,PE_IO) !//PAR_I/O BCAST1d
                                      call bcast(U_S0,PE_IO)    !//PAR_I/O BCAST1d
                                      call bcast(V_S0,PE_IO)    !//PAR_I/O BCAST1d
                                      call bcast(W_S0,PE_IO)    !//PAR_I/O BCAST1d

                                      NEXT_RECA = NEXT_RECA + 1 
                                    ENDIF 
                                    
!
                                    IF (VERSION_NUMBER >= 1.09) THEN 
                                      READ (UNIT_RES, REC=NEXT_RECA) TIME, TSTOP&
                                       , ENERGY_EQ, RES_DT, OUT_DT, NLOG, &
                                       L_SCALE0, NO_I, NO_J, NO_K, CALL_USR 

                                      call bcast(TIME,PE_IO)    !//PAR_I/O BCAST0d
                                      call bcast(TSTOP,PE_IO)   !//PAR_I/O BCAST0d
                                      call bcast(RES_DT,PE_IO)  !//PAR_I/O BCAST0d
                                      call bcast(OUT_DT,PE_IO)  !//PAR_I/O BCAST0d
                                      call bcast(L_SCALE0,PE_IO)!//PAR_I/O BCAST0d
                                      call bcast(NLOG,PE_IO)    !//PAR_I/O BCAST0i
                                      call bcast(ENERGY_EQ,PE_IO) !//PAR_I/O BCAST0l 
                                      call bcast(NO_I,PE_IO)      !//PAR_I/O BCAST0l
                                      call bcast(NO_J,PE_IO)      !//PAR_I/O BCAST0l
                                      call bcast(NO_K,PE_IO)       !//PAR_I/O BCAST0l
                                      call bcast(CALL_USR,PE_IO)  !//PAR_I/O BCAST0l

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

                                      call bcast(SPX_DT,PE_IO) !//PAR_I/O BCAST1d
                                      call bcast(SPECIES_EQ,PE_IO) !//PAR_I/O BCAST1l (recv)

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

                                      call bcast(USR_DT,PE_IO)  !//PAR_I/O BCAST1d
                                      call bcast(USR_X_W,PE_IO) !//PAR_I/O BCAST1d
                                      call bcast(USR_X_E,PE_IO) !//PAR_I/O BCAST1d
                                      call bcast(USR_Y_S,PE_IO) !//PAR_I/O BCAST1d
                                      call bcast(USR_Y_N,PE_IO) !//PAR_I/O BCAST1d
!//? do we need to bcast to all or only USR_Z_T to last PE, does the interior PE need this info
                                      call bcast(USR_Z_B,PE_IO) !//PAR_I/O BCAST1d
                                      call bcast(USR_Z_T,PE_IO) !//PAR_I/O BCAST1d

                                      DO LC = 1, DIMENSION_USR 
                                        READ (UNIT_RES, REC=NEXT_RECA) USR_FORMAT(&
                                             LC), USR_EXT(LC), USR_TYPE(LC), USR_VAR(&
                                             LC) 
                                        NEXT_RECA = NEXT_RECA + 1 
                                      END DO 

                                      call bcast(USR_FORMAT,PE_IO) !//PAR_I/O BCAST1c 
                                      call bcast(USR_EXT,PE_IO) !//PAR_I/O BCAST1c 
                                      call bcast(USR_TYPE,PE_IO) !//PAR_I/O BCAST1c
                                      call bcast(USR_VAR,PE_IO) !//PAR_I/O BCAST1c

                                      CALL IN_BIN_512 (UNIT_RES, IC_P_STAR, &
                                           DIM_IC, NEXT_RECA) 
                                      CALL IN_BIN_512 (UNIT_RES, IC_L_SCALE, &
                                           DIM_IC, NEXT_RECA) 

                                      call bcast(IC_P_STAR,PE_IO) !//PAR_I/O BCAST1d
                                      call bcast(IC_L_SCALE,PE_IO) !//PAR_I/O BCAST1d

                                      DO LC = 1, DIM_IC 
                                        READ (UNIT_RES, REC=NEXT_RECA) IC_TYPE(LC) 
                                        NEXT_RECA = NEXT_RECA + 1 
                                      END DO 

                                      call bcast(IC_TYPE,PE_IO) !//PAR_I/O BCAST1c

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

                                      call bcast(BC_DT_0,PE_IO) !//PAR_I/O BCAST1d
                                      call bcast(BC_DT_H,PE_IO) !//PAR_I/O BCAST1d
                                      call bcast(BC_DT_L,PE_IO) !//PAR_I/O BCAST1d
                                      call bcast(BC_JET_G0,PE_IO) !//PAR_I/O BCAST1d
                                      call bcast(BC_JET_GH,PE_IO) !//PAR_I/O BCAST1d
                                      call bcast(BC_JET_GL,PE_IO) !//PAR_I/O BCAST1d
                                    ENDIF 

!
                                    IF (VERSION_NUMBER >= 1.10) THEN 
                                    READ (UNIT_RES, REC=NEXT_RECA) MU_GMAX 
                                    call bcast(MU_GMAX,PE_IO) !//PAR_I/O BCAST0d
                                    NEXT_RECA = NEXT_RECA + 1 
                                    ENDIF 

!
                                    IF (VERSION_NUMBER >= 1.11) THEN 
                                    READ (UNIT_RES, REC=NEXT_RECA) V_EX, &
                                       MODEL_B 
                                    call bcast(V_EX,PE_IO) !//PAR_I/O BCAST0d 
                                    call bcast(MODEL_B,PE_IO) !//PAR_I/O BCAST0l 
                                    NEXT_RECA = NEXT_RECA + 1 
                                    ENDIF 
!
                                    IF (VERSION_NUMBER >= 1.12) THEN 
                                    READ (UNIT_RES, REC=NEXT_RECA) P_REF, &
                                       P_SCALE, UR_FAC, TOL_RESID, DT_MAX, &
                                       DT_MIN, DT_FAC, CLOSE_PACKED, GRAVITY, &
                                       MU_S0 
!//S may be worth to pack them and then bcast
                                    call bcast(P_REF,PE_IO) !//PAR_I/O BCAST0d
                                    call bcast(P_SCALE,PE_IO) !//PAR_I/O BCAST0d
                                    call bcast(UR_FAC,PE_IO) !//PAR_I/O BCAST0d
                                    call bcast(TOL_RESID,PE_IO) !//PAR_I/O BCAST0d
                                    call bcast(DT_MAX,PE_IO) !//PAR_I/O BCAST0d
                                    call bcast(DT_MIN,PE_IO) !//PAR_I/O BCAST0d
                                    call bcast(DT_FAC,PE_IO) !//PAR_I/O BCAST0d
                                    call bcast(GRAVITY,PE_IO) !//PAR_I/O BCAST0d
                                    call bcast(MU_S0,PE_IO) !//PAR_I/O BCAST0d
                                    call bcast(CLOSE_PACKED,PE_IO) !//PAR_I/O BCAST0l

                                    NEXT_RECA = NEXT_RECA + 1 

                                    READ (UNIT_RES, REC=NEXT_RECA) LEQ_IT, &
                                       LEQ_METHOD 
                                    NEXT_RECA = NEXT_RECA + 1 
                                    call bcast(LEQ_IT,PE_IO)     !//PAR_I/O BCAST1i 
                                    call bcast(LEQ_METHOD,PE_IO) !//PAR_I/O BCAST1i 

                                    CALL IN_BIN_512 (UNIT_RES, BC_HW_G, DIM_BC&
                                       , NEXT_RECA) 
                                    CALL IN_BIN_512 (UNIT_RES, BC_UW_G, DIM_BC&
                                       , NEXT_RECA) 
                                    CALL IN_BIN_512 (UNIT_RES, BC_VW_G, DIM_BC&
                                       , NEXT_RECA) 
                                    CALL IN_BIN_512 (UNIT_RES, BC_WW_G, DIM_BC&
                                       , NEXT_RECA) 

                                    call bcast(BC_HW_G,PE_IO) !//PAR_I/O BCAST1d 
                                    call bcast(BC_UW_G,PE_IO) !//PAR_I/O BCAST1d 
                                    call bcast(BC_VW_G,PE_IO) !//PAR_I/O BCAST1d 
                                    call bcast(BC_WW_G,PE_IO) !//PAR_I/O BCAST1d 


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

                                    call bcast(BC_HW_S,PE_IO) !//PAR_I/O BCAST2d 
                                    call bcast(BC_UW_S,PE_IO) !//PAR_I/O BCAST2d 
                                    call bcast(BC_VW_S,PE_IO) !//PAR_I/O BCAST2d 
                                    call bcast(BC_WW_S,PE_IO) !//PAR_I/O BCAST2d 
                                    ENDIF 

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

                                    call bcast(TOL_DIVERGE,PE_IO) !//PAR_I/O BCAST0d 
                                    call bcast(DISCRETIZE,PE_IO) !//PAR_I/O BCAST1i 
                                    call bcast(MOMENTUM_X_EQ,PE_IO) !//PAR_I/O BCAST1l 
                                    call bcast(MOMENTUM_Y_EQ,PE_IO) !//PAR_I/O BCAST1l 
                                    call bcast(MOMENTUM_Z_EQ,PE_IO) !//PAR_I/O BCAST1l 
                                    call bcast(FULL_LOG,PE_IO) !//PAR_I/O BCAST0l 

                                    NEXT_RECA = NEXT_RECA + 1 
                                    ENDIF 
!
                                    IF (VERSION_NUMBER >= 1.14) THEN 
                                      READ (UNIT_RES, REC=NEXT_RECA) DETECT_STALL 
                                      call bcast(DETECT_STALL,PE_IO) !//PAR_I/O BCAST0l 
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
!//S may be worth to pack them and then bcast
                                      call bcast(K_G0,PE_IO) !//PAR_I/O BCAST0d 
                                      call bcast(K_S0,PE_IO) !//PAR_I/O BCAST0d 
                                      call bcast(C_PG0,PE_IO) !//PAR_I/O BCAST0d 
                                      call bcast(C_PS0,PE_IO) !//PAR_I/O BCAST0d 
                                      call bcast(TOL_RESID_T,PE_IO) !//PAR_I/O BCAST0d
                                      call bcast(TOL_RESID_X,PE_IO) !//PAR_I/O BCAST0d 

                                      call bcast(IC_GAMA_RG,PE_IO) !//PAR_I/O BCAST1d 
                                      call bcast(IC_T_RG,PE_IO) !//PAR_I/O BCAST1d 

                                      DO LC = 1, MMAX 
                                       CALL IN_BIN_512 (UNIT_RES, IC_GAMA_RS(1,LC)&
                                       , DIM_IC, NEXT_RECA) 
                                       CALL IN_BIN_512 (UNIT_RES, IC_T_RS(1,LC), &
                                       DIM_IC, NEXT_RECA) 
                                      END DO 
                                      call bcast(IC_GAMA_RS,PE_IO) !//PAR_I/O BCAST2d 
                                      call bcast(IC_T_RS,PE_IO) !//PAR_I/O BCAST2d 

                                    ENDIF 
!
                                    IF (VERSION_NUMBER >= 1.2) THEN 
                                      READ (UNIT_RES, REC=NEXT_RECA) NORM_G, &
                                        NORM_S 
                                      call bcast(NORM_G,PE_IO) !//PAR_I/O BCAST0d 
                                      call bcast(NORM_S,PE_IO) !//PAR_I/O BCAST0d 

                                      NEXT_RECA = NEXT_RECA + 1 
                                    ENDIF 
!
!  Add new read statements above this line.  Remember to update NEXT_RECA.
!  Remember to update the version number check near begining of this subroutine.
!------------------------------------------------------------------------------
!
                                    READ (UNIT_RES, REC=3) NEXT_RECA 


                               else    ! else of if (myPE == PE_IO) from Block 2
			       			       
                                    call bcast(dx, PE_IO)         !//PAR_I/O BCAST1d (recv)

                                    call bcast(dy, PE_IO)         !//PAR_I/O BCAST1d (recv)
!//AIKEPARDBGSTOP 0907
    write(*,"('(PE ',I2,'): aft dy in Block2 in read_res0')") myPE !//AIKEPARDBG
    write(UNIT_LOG,*) dy  !//AIKEPARDBG  
    call exitMPI(myPE)    !//AIKEPARDBG

                                    call gather(dz,dGTEMP,PE_IO)  !//PAR_I/O GATHER1d  (recv)
!    call MPI_Barrier(MPI_COMM_WORLD,mpierr)					    
    write(*,"('(PE ',I2,'): KMAX3',I2)") myPE, KMAX3 !//AIKEPARDBG
    write(UNIT_LOG,"('(PE ',I2,'): aft dz in Block2 in read_res0')") myPE !//AIKEPARDBG
    write(UNIT_LOG,*) dz    !//AIKEPARDBG
    call exitMPI(myPE)      !//AIKEPARDBG
				    
                                    call bcast(RUN_NAME,PE_IO)    !//PAR_I/O BCAST0c (recv)
                                    call bcast(DESCRIPTION,PE_IO) !//PAR_I/O BCAST0c (recv)
                                    call bcast(UNITS,PE_IO)       !//PAR_I/O BCAST0c (recv)
                                    call bcast(RUN_TYPE,PE_IO)    !//PAR_I/O BCAST0c (recv)
                                    call bcast(COORDINATES,PE_IO) !//PAR_I/O BCAST0c (recv)
                                    call bcast(D_P, PE_IO)       !//PAR_I/O BCAST1d (recv)
                                    call bcast(RO_S, PE_IO)      !//PAR_I/O BCAST1d (recv)
                                    call bcast(EP_STAR, PE_IO)   !//PAR_I/O BCAST0d (recv)
                                    if (VERSION /= 'RES = 01.00' .OR. VERSION /= &
                                         'RES = 01.01') &
                                    call bcast(RO_G0, PE_IO)     !//PAR_I/O BCAST0d (recv)
                                    call bcast(MU_G0, PE_IO)     !//PAR_I/O BCAST0d (recv)
                                    call bcast(MW_AVG, PE_IO)    !//PAR_I/O BCAST0d (recv)
                                    if (VERSION_NUMBER >= 1.04) then
                                       call bcast(MW_G, PE_IO)     !//PAR_I/O BCAST1d (recv)
                                       call bcast(MW_S, PE_IO)     !//PAR_I/O BCAST2d (recv)
                                    endif
                                    call bcast(IC_X_W, PE_IO)     !//PAR_I/O BCAST1d (recv)
                                    call bcast(IC_X_E, PE_IO)     !//PAR_I/O BCAST1d (recv)
                                    call bcast(IC_Y_S, PE_IO)     !//PAR_I/O BCAST1d (recv)
                                    call bcast(IC_Y_N, PE_IO)     !//PAR_I/O BCAST1d (recv)
!//? do we need to bcast to all or only Z_T to last PE
                                    call bcast(IC_Z_B, PE_IO)     !//PAR_I/O BCAST1d (recv)
                                    call bcast(IC_Z_T, PE_IO)     !//PAR_I/O BCAST1d (recv)

                                    call bcast(IC_I_W, PE_IO)     !//PAR_I/O BCAST1i (recv)
                                    call bcast(IC_I_E, PE_IO)     !//PAR_I/O BCAST1i (recv)
                                    call bcast(IC_J_S, PE_IO)     !//PAR_I/O BCAST1i (recv)
                                    call bcast(IC_J_N, PE_IO)     !//PAR_I/O BCAST1i (recv)

!//? do we need to bcast to all or only IC_Z_T to last PE, does the interior PE need this info
                                    call bcast(IC_K_B, PE_IO)     !//PAR_I/O BCAST1i (recv)
                                    call bcast(IC_K_T, PE_IO)     !//PAR_I/O BCAST1i (recv)

                                    call bcast(IC_EP_G, PE_IO)    !//PAR_I/O BCAST1d (recv)
                                    call bcast(IC_P_G, PE_IO)     !//PAR_I/O BCAST1d (recv)
                                    call bcast(IC_T_G, PE_IO)     !//PAR_I/O BCAST1d (recv)

                                    if (VERSION_NUMBER < 1.15) then
                                       IF (MMAX < 2) THEN 
                                         call bcast(IC_TMP, PE_IO) !//PAR_I/O BCAST1d (recv)
                                       ENDIF 
                                       call bcast(IC_T_S, PE_IO)   !//PAR_I/O BCAST2d (recv)
                                    endif

                                    if (VERSION_NUMBER >= 1.04) then
                                       call bcast(IC_X_G, PE_IO)   !//PAR_I/O BCAST2d (recv)
                                    endif

                                    call bcast(IC_U_G, PE_IO)    !//PAR_I/O BCAST1d (recv)
                                    call bcast(IC_V_G, PE_IO)    !//PAR_I/O BCAST1d (recv)
                                    call bcast(IC_W_G, PE_IO)    !//PAR_I/O BCAST1d (recv)

                                    call bcast(IC_ROP_S, PE_IO)  !//PAR_I/O BCAST2d (recv)
                                    call bcast(IC_U_S, PE_IO)    !//PAR_I/O BCAST2d (recv)
                                    call bcast(IC_V_S, PE_IO)    !//PAR_I/O BCAST2d (recv)
                                    call bcast(IC_W_S, PE_IO)    !//PAR_I/O BCAST2d (recv)
                                    if (VERSION_NUMBER >= 1.15) &
                                       call bcast(IC_T_S, PE_IO)  !//PAR_I/O BCAST2d (recv)
                                    if (VERSION_NUMBER >= 1.04) &
                                       call bcast(IC_X_S, PE_IO)  !//PAR_I/O BCAST3d (recv)

                                    call bcast(BC_X_W, PE_IO)    !//PAR_I/O BCAST1d (recv)
                                    call bcast(BC_X_E, PE_IO)    !//PAR_I/O BCAST1d (recv)
                                    call bcast(BC_Y_S, PE_IO)    !//PAR_I/O BCAST1d (recv)
                                    call bcast(BC_Y_N, PE_IO)    !//PAR_I/O BCAST1d (recv)

!//? do we need to bcast to all or only BC_Z_T to last PE, does the interior PE need this info
                                    call bcast(BC_Z_B, PE_IO)    !//PAR_I/O BCAST1d (recv)
                                    call bcast(BC_Z_T, PE_IO)    !//PAR_I/O BCAST1d (recv)

                                    call bcast(BC_I_W, PE_IO)    !//PAR_I/O BCAST1i (recv)
                                    call bcast(BC_I_E, PE_IO)    !//PAR_I/O BCAST1i (recv)
                                    call bcast(BC_J_S, PE_IO)    !//PAR_I/O BCAST1i (recv)
                                    call bcast(BC_J_N, PE_IO)    !//PAR_I/O BCAST1i (recv)

!//? do we need to bcast to all or only BC_K_T to last PE, does the interior PE need this info
                                    call bcast(BC_K_B, PE_IO)    !//PAR_I/O BCAST1i (recv)
                                    call bcast(BC_K_T, PE_IO)    !//PAR_I/O BCAST1i (recv)

                                    call bcast(BC_EP_G, PE_IO)   !//PAR_I/O BCAST1d (recv)
                                    call bcast(BC_P_G, PE_IO)    !//PAR_I/O BCAST1d (recv)
                                    call bcast(BC_T_G, PE_IO)    !//PAR_I/O BCAST1d (recv)


                                    IF (VERSION_NUMBER < 1.15) THEN 
                                      IF (MMAX < 2) THEN 
                                        call bcast(BC_TMP, PE_IO)    !//PAR_I/O BCAST1d (recv)
                                      ENDIF 
                                      call bcast(BC_T_S, PE_IO)   !//PAR_I/O BCAST2d (recv)
                                    ENDIF 
!
                                    IF (VERSION_NUMBER >= 1.04) THEN 
                                      call bcast(BC_X_G, PE_IO)   !//PAR_I/O BCAST2d (recv)
                                    ENDIF 

                                    call bcast(BC_U_G, PE_IO)   !//PAR_I/O BCAST1d (recv)
                                    call bcast(BC_V_G, PE_IO)   !//PAR_I/O BCAST1d (recv)
                                    call bcast(BC_W_G, PE_IO)   !//PAR_I/O BCAST1d (recv)
                                    call bcast(BC_RO_G, PE_IO)   !//PAR_I/O BCAST1d (recv)
                                    call bcast(BC_ROP_G, PE_IO)   !//PAR_I/O BCAST1d (recv)
                                    call bcast(BC_VOLFLOW_G, PE_IO)   !//PAR_I/O BCAST1d (recv)
                                    call bcast(BC_MASSFLOW_G, PE_IO)   !//PAR_I/O BCAST1d (recv)

                                    call bcast(BC_ROP_S, PE_IO) !//PAR_I/O BCAST2d (recv)
                                    call bcast(BC_U_S, PE_IO)   !//PAR_I/O BCAST2d (recv)
                                    call bcast(BC_V_S, PE_IO)   !//PAR_I/O BCAST2d (recv)

                                    if (VERSION_NUMBER >= 1.04) then
                                       call bcast(BC_W_S, PE_IO)   !//PAR_I/O BCAST2d (recv)
                                       if (VERSION_NUMBER >= 1.15)  &
                                          call bcast(BC_T_S, PE_IO)   !//PAR_I/O BCAST2d (recv)
!
                                       call bcast(BC_X_S, PE_IO)   !//PAR_I/O BCAST2d (recv)
                                    endif 
                                    call bcast(BC_VOLFLOW_S, PE_IO)   !//PAR_I/O BCAST2d (recv)
                                    call bcast(BC_MASSFLOW_S, PE_IO)   !//PAR_I/O BCAST2d (recv)

                                    call bcast(BC_TYPE, PE_IO)   !//PAR_I/O BCAST1c (recv)

                                    call gather(array1i,iGTEMP2,PE_IO)  !//PAR_I/O GATHER1d

                                    IF (VERSION_NUMBER >= 1.04) THEN 

                                      call bcast(IS_X_W, PE_IO)   !//PAR_I/O BCAST1d (recv)
                                      call bcast(IS_X_E, PE_IO)   !//PAR_I/O BCAST1d (recv)
                                      call bcast(IS_Y_S, PE_IO)   !//PAR_I/O BCAST1d (recv)
                                      call bcast(IS_Y_N, PE_IO)   !//PAR_I/O BCAST1d (recv)
!//? do we need to bcast to all or only IS_Z_T to last PE, does the interior PE need this info
                                      call bcast(IS_Z_B, PE_IO)   !//PAR_I/O BCAST1d (recv)
                                      call bcast(IS_Z_T, PE_IO)   !//PAR_I/O BCAST1d (recv)

                                      call bcast(IS_I_W, PE_IO)   !//PAR_I/O BCAST1i (recv)
                                      call bcast(IS_I_E, PE_IO)   !//PAR_I/O BCAST1i (recv)
                                      call bcast(IS_J_S, PE_IO)   !//PAR_I/O BCAST1i (recv)
                                      call bcast(IS_J_N, PE_IO)   !//PAR_I/O BCAST1i (recv)

                                      IF (VERSION_NUMBER >= 1.07) & 
                                        call bcast(IS_VEL_S, PE_IO)!//PAR_I/O BCAST2d  (recv)
                                      
                                      call bcast(IS_TYPE, PE_IO)   !//PAR_I/O BCAST1c (recv)
                                    ENDIF 

                                    if (VERSION_NUMBER >= 1.08) then
                                      call bcast(CYCLIC_X,PE_IO)     !//PAR_I/O BCAST0l
                                      call bcast(CYCLIC_Y,PE_IO)     !//PAR_I/O BCAST0l
                                      call bcast(CYCLIC_Z,PE_IO)     !//PAR_I/O BCAST0l
                                      call bcast(CYCLIC_X_PD,PE_IO)  !//PAR_I/O BCAST0l
                                      call bcast(CYCLIC_Y_PD,PE_IO)  !//PAR_I/O BCAST0l
                                      call bcast(CYCLIC_Z_PD,PE_IO)  !//PAR_I/O BCAST0l                                    
                                      call bcast(DBLPACK,PE_IO) !//PAR_I/O BCAST1d (recv)
                                      DELP_X = DBLPACK(1)
                                      DELP_Y = DBLPACK(2) 
                                      DELP_Z = DBLPACK(3)
                                      U_G0 = DBLPACK(4)
                                      V_G0 = DBLPACK(5)
                                      W_G0 = DBLPACK(6)

                                      call bcast(U_S0,PE_IO) !//PAR_I/O BCAST1d (recv)
                                      call bcast(V_S0,PE_IO) !//PAR_I/O BCAST1d (recv)
                                      call bcast(W_S0,PE_IO) !//PAR_I/O BCAST1d (recv)
                                    endif


                                    IF (VERSION_NUMBER >= 1.09) THEN 
                                      call bcast(TIME,PE_IO)    !//PAR_I/O BCAST0d (recv)
                                      call bcast(TSTOP,PE_IO)   !//PAR_I/O BCAST0d (recv)
                                      call bcast(RES_DT,PE_IO)  !//PAR_I/O BCAST0d (recv)
                                      call bcast(OUT_DT,PE_IO)  !//PAR_I/O BCAST0d (recv)
                                      call bcast(L_SCALE0,PE_IO)!//PAR_I/O BCAST0d (recv)
                                      call bcast(NLOG,PE_IO)    !//PAR_I/O BCAST0i (recv)
                                      call bcast(ENERGY_EQ,PE_IO)!//PAR_I/O BCAST0l (recv)
                                      call bcast(NO_I,PE_IO)     !//PAR_I/O BCAST0l (recv)
                                      call bcast(NO_J,PE_IO)     !//PAR_I/O BCAST0l (recv)
                                      call bcast(NO_K,PE_IO)      !//PAR_I/O BCAST0l (recv)
                                      call bcast(CALL_USR,PE_IO) !//PAR_I/O BCAST0l (recv)

                                      call bcast(SPX_DT,PE_IO) !//PAR_I/O BCAST1d (recv)
                                      call bcast(SPECIES_EQ,PE_IO) !//PAR_I/O BCAST1l (recv)

                                      call bcast(USR_DT,PE_IO)  !//PAR_I/O BCAST1d (recv)
                                      call bcast(USR_X_W,PE_IO) !//PAR_I/O BCAST1d (recv)
                                      call bcast(USR_X_E,PE_IO) !//PAR_I/O BCAST1d (recv)
                                      call bcast(USR_Y_S,PE_IO) !//PAR_I/O BCAST1d (recv)
                                      call bcast(USR_Y_N,PE_IO) !//PAR_I/O BCAST1d (recv)
!//? do we need to bcast to all or only USR_Z_T to last PE, does the interior PE need this info
                                      call bcast(USR_Z_B,PE_IO) !//PAR_I/O BCAST1d (recv)
                                      call bcast(USR_Z_T,PE_IO) !//PAR_I/O BCAST1d (recv)


                                      call bcast(USR_FORMAT,PE_IO) !//PAR_I/O BCAST1c (recv)
                                      call bcast(USR_EXT,PE_IO) !//PAR_I/O BCAST1c (recv)
                                      call bcast(USR_TYPE,PE_IO) !//PAR_I/O BCAST1c (recv)
                                      call bcast(USR_VAR,PE_IO) !//PAR_I/O BCAST1c (recv)


                                      call bcast(IC_P_STAR,PE_IO) !//PAR_I/O BCAST1d (recv)
                                      call bcast(IC_L_SCALE,PE_IO) !//PAR_I/O BCAST1d (recv)

                                      call bcast(IC_TYPE,PE_IO) !//PAR_I/O BCAST1c

                                      call bcast(BC_DT_0,PE_IO) !//PAR_I/O BCAST1d (recv)
                                      call bcast(BC_DT_H,PE_IO) !//PAR_I/O BCAST1d (recv)
                                      call bcast(BC_DT_L,PE_IO) !//PAR_I/O BCAST1d (recv)
                                      call bcast(BC_JET_G0,PE_IO) !//PAR_I/O BCAST1d (recv)
                                      call bcast(BC_JET_GH,PE_IO) !//PAR_I/O BCAST1d (recv)
                                      call bcast(BC_JET_GL,PE_IO) !//PAR_I/O BCAST1d (recv)
                                    ENDIF 

                                    IF (VERSION_NUMBER >= 1.10) &
                                       call bcast(MU_GMAX,PE_IO) !//PAR_I/O BCAST0d (recv)
!
                                    IF (VERSION_NUMBER >= 1.11) THEN 
                                    call bcast(V_EX,PE_IO) !//PAR_I/O BCAST0d (recv)
                                    call bcast(MODEL_B,PE_IO) !//PAR_I/O BCAST0l (recv)
                                    ENDIF 
!
                                    IF (VERSION_NUMBER >= 1.12) THEN 
!//S may be worth to pack them and then bcast
                                     call bcast(P_REF,PE_IO) !//PAR_I/O BCAST0d (recv)
                                     call bcast(P_SCALE,PE_IO) !//PAR_I/O BCAST0d (recv)
                                     call bcast(UR_FAC,PE_IO) !//PAR_I/O BCAST0d (recv)
                                     call bcast(TOL_RESID,PE_IO) !//PAR_I/O BCAST0d (recv)
                                     call bcast(DT_MAX,PE_IO) !//PAR_I/O BCAST0d (recv)
                                     call bcast(DT_MIN,PE_IO) !//PAR_I/O BCAST0d (recv)
                                     call bcast(DT_FAC,PE_IO) !//PAR_I/O BCAST0d (recv)
                                     call bcast(GRAVITY,PE_IO) !//PAR_I/O BCAST0d (recv)
                                     call bcast(MU_S0,PE_IO) !//PAR_I/O BCAST0d (recv)
                                     call bcast(CLOSE_PACKED,PE_IO) !//PAR_I/O BCAST0l (recv)

                                     call bcast(LEQ_IT,PE_IO)     !//PAR_I/O BCAST1i (recv) 
                                     call bcast(LEQ_METHOD,PE_IO) !//PAR_I/O BCAST1i (recv) 

                                     call bcast(BC_HW_G,PE_IO) !//PAR_I/O BCAST1d (recv) 
                                     call bcast(BC_UW_G,PE_IO) !//PAR_I/O BCAST1d (recv) 
                                     call bcast(BC_VW_G,PE_IO) !//PAR_I/O BCAST1d (recv) 
                                     call bcast(BC_WW_G,PE_IO) !//PAR_I/O BCAST1d (recv)

                                     call bcast(BC_HW_S,PE_IO) !//PAR_I/O BCAST2d (recv) 
                                     call bcast(BC_UW_S,PE_IO) !//PAR_I/O BCAST2d (recv) 
                                     call bcast(BC_VW_S,PE_IO) !//PAR_I/O BCAST2d (recv) 
                                     call bcast(BC_WW_S,PE_IO) !//PAR_I/O BCAST2d (recv) 
                                    ENDIF 

                                    LC = 0 
                                    IF (MMAX + 1 > 0) THEN 
                                    MOMENTUM_X_EQ(:MMAX) = .TRUE. 
                                    MOMENTUM_Y_EQ(:MMAX) = .TRUE. 
                                    MOMENTUM_Z_EQ(:MMAX) = .TRUE. 
                                    LC = MMAX + 1 
                                    ENDIF 
                                    TOL_DIVERGE = 1.E+4 

                                    IF (VERSION_NUMBER >= 1.13) THEN 
                                    call bcast(TOL_DIVERGE,PE_IO) !//PAR_I/O BCAST0d (recv) 
                                    call bcast(DISCRETIZE,PE_IO) !//PAR_I/O BCAST1i (recv) 
                                    call bcast(MOMENTUM_X_EQ,PE_IO) !//PAR_I/O BCAST1l (recv) 
                                    call bcast(MOMENTUM_Y_EQ,PE_IO) !//PAR_I/O BCAST1l (recv) 
                                    call bcast(MOMENTUM_Z_EQ,PE_IO) !//PAR_I/O BCAST1l (recv) 
                                    call bcast(FULL_LOG,PE_IO) !//PAR_I/O BCAST0l (recv) 
                                    ENDIF 
!
                                    IF (VERSION_NUMBER >= 1.14) THEN 
                                      call bcast(DETECT_STALL,PE_IO) !//PAR_I/O BCAST0l (recv) 
                                    ENDIF 
!
                                    IF (VERSION_NUMBER >= 1.15) THEN 
!//S may be worth to pack them and then bcast
                                      call bcast(K_G0,PE_IO) !//PAR_I/O BCAST0d (recv) 
                                      call bcast(K_S0,PE_IO) !//PAR_I/O BCAST0d (recv) 
                                      call bcast(C_PG0,PE_IO) !//PAR_I/O BCAST0d (recv) 
                                      call bcast(C_PS0,PE_IO) !//PAR_I/O BCAST0d (recv) 
                                      call bcast(TOL_RESID_T,PE_IO) !//PAR_I/O BCAST0d (recv)
                                      call bcast(TOL_RESID_X,PE_IO) !//PAR_I/O BCAST0d (recv) 

                                      call bcast(IC_GAMA_RG,PE_IO) !//PAR_I/O BCAST1d (recv) 
                                      call bcast(IC_T_RG,PE_IO) !//PAR_I/O BCAST1d (recv)

                                      call bcast(IC_GAMA_RS,PE_IO) !//PAR_I/O BCAST2d (recv) 
                                      call bcast(IC_T_RS,PE_IO) !//PAR_I/O BCAST2d (recv) 
                                    ENDIF 
!
                                    IF (VERSION_NUMBER >= 1.2) THEN 
                                      call bcast(NORM_G,PE_IO) !//PAR_I/O BCAST0d (recv) 
                                      call bcast(NORM_S,PE_IO) !//PAR_I/O BCAST0d (recv) 
                                    ENDIF 


                                    endif  ! end of if (myPE == PE_IO)

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
      WRITE (*, "('(PE ',I3,'): From: READ_RES0')") myPE 
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
!//PAR_I/O 0815 write the error message into the XXX.log files
      WRITE (UNIT_LOG, *) ' ' 
      WRITE (UNIT_LOG, *) ' **************************************' 
      WRITE (UNIT_LOG, *) ' From: READ_RES0' 
      WRITE (UNIT_LOG, *) ' DIMENSION ERROR ---' 
      WRITE (UNIT_LOG, *) ' ' 
      WRITE (UNIT_LOG, *) ' DIMENSION_I  = ', DIMENSION_I, ' IMAX2        = ', IMAX2 
      WRITE (UNIT_LOG, *) ' DIMENSION_J  = ', DIMENSION_J, ' JMAX2        = ', JMAX2 
      WRITE (UNIT_LOG, *) ' DIMENSION_K  = ', DIMENSION_K, ' KMAX2        = ', KMAX2 
      WRITE (UNIT_LOG, *) ' DIMENSION_3  = ', DIMENSION_3, ' IJKMAX2      = ', IJKMAX2 
      WRITE (UNIT_LOG, *) ' DIMENSION_M  = ', DIMENSION_M, ' MMAX         = ', MMAX 
      WRITE (UNIT_LOG, *) ' DIMENSION_IC = ', DIMENSION_IC, ' DIM_IC       = ', DIM_IC 
      WRITE (UNIT_LOG, *) ' DIMENSION_BC = ', DIMENSION_BC, ' DIM_BC       = ', DIM_BC 
      WRITE (UNIT_LOG, *) ' DIMENSION_IS = ', DIMENSION_IS, ' DIM_IS       = ', DIM_IS 
      WRITE (UNIT_LOG, *) ' DIMENSION_C  = ', DIMENSION_C, ' DIM_C        = ', DIM_C 
      WRITE (UNIT_LOG, *) ' DIMENSION_N_g= ', DIMENSION_N_G, ' NMAX(0)      = ', NMAX(&
         0) 
      WRITE (UNIT_LOG, *) ' DIMENSION_N_s= ', DIMENSION_N_S 
      DO M = 1, MMAX 
         WRITE (UNIT_LOG, '(A, I2, A, I4)') ' NMAX(', M, ') = ', NMAX(M) 
      END DO 
      WRITE (UNIT_LOG, *) ' ' 
!
      call unlock_tmp_array

!// 990 0807 Abort all PEs, not only the current one
!      STOP  
      call exitMPI(myPE)

      END SUBROUTINE READ_RES0 
