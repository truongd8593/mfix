!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_VTU_FILE                                         C
!  Purpose: Writes the cut cell grid in VTK format (Unstructured VTU)  C
!           Binary format                                              C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE WRITE_VTU_FILE  
    
      USE param
      USE param1
      USE parallel
      USE constant
      USE run
      USE toleranc
      USE geometry
      USE indices  
      USE compar
      USE sendrecv
      USE quadric
      USE cutcell
      USE fldvar
      USE visc_s
      USE physprop
      USE pgcor
      USE vtk
      USE rxns      
      USE output
      USE scalars
      USE stl 

      USE mpi_utility     
      USE parallel_mpi


      USE pgcor
      USE pscor
      USE discretelement, Only : DES_CELLWISE_BCDATA, DISCRETE_ELEMENT
      USE mfix_pic

      IMPLICIT NONE
      DOUBLE PRECISION:: Xw,Xe,Yn,Ys
      INTEGER :: I,J,K,L,M,N,R,IM,JM,KM,IP,JP,KP,IJK
      INTEGER :: IMJK,IJMK,IJKM,IMJMK,IMJKM,IJMKM,IMJMKM

      INTEGER sw,se,ne,nw
      INTEGER, DIMENSION(10) :: additional_node
      DOUBLE PRECISION, DIMENSION(2*DIMENSION_3) ::  X_OF
      DOUBLE PRECISION, DIMENSION(2*DIMENSION_3) ::  Y_OF
      DOUBLE PRECISION, DIMENSION(2*DIMENSION_3) ::  Z_OF
      INTEGER, DIMENSION(DIMENSION_3) ::  INDEX_OF_E_ADD_NODE
      INTEGER, DIMENSION(DIMENSION_3) ::  INDEX_OF_N_ADD_NODE
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  FACET_COUNT_DES, NEIGHBORING_FACET

      INTEGER :: SPECIES_COUNTER,LT

      CHARACTER (LEN=32) :: SUBM,SUBN,SUBR
      CHARACTER (LEN=64) :: VAR_NAME

      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  DP_BC_ID, COUNT_DES_BC,IJK_ARRAY

      INTEGER :: PASS
      INTEGER :: WRITE_HEADER = 1
      INTEGER :: WRITE_DATA   = 2


      include "function.inc"

      IF(.NOT.CARTESIAN_GRID) RETURN


      DX(IEND3+1) = DX(IEND3)
      DY(JEND3+1) = DY(JEND3)
      DZ(KEND3+1) = DZ(KEND3)

!     Location of U-momentum cells for original (uncut grid)
      IF (DO_I) THEN 
        XG_E(1) = ZERO
        DO I = IMIN1, IMAX2 
           XG_E(I) = XG_E(I-1) + DX(I) 
        END DO 
      ENDIF

!     Location of V-momentum cells for original (uncut grid)
      IF (DO_J) THEN 
        YG_N(1) = ZERO
        DO J = JMIN1, JMAX2 
           YG_N(J) = YG_N(J-1) + DY(J) 
        END DO 
      ENDIF

!     Location of W-momentum cells for original (uncut grid)
      IF (DO_K) THEN 
        ZG_T(1) = ZERO
        DO K = KMIN1, KMAX2 
           ZG_T(K) = ZG_T(K-1) + DZ(K) 
        END DO 
      ELSE
         ZG_T = ZERO
      ENDIF

      CALL OPEN_VTU_FILE_BIN
      CALL OPEN_PVD_FILE

      DO PASS=WRITE_HEADER,WRITE_DATA


         CALL WRITE_GEOMETRY_IN_VTU_BIN(PASS)

         DO L = 1, DIM_VTK_VAR

            SELECT CASE (VTK_VAR(L))
              
               CASE (1)
                  CALL WRITE_SCALAR_IN_VTU_BIN('EP_G',EP_G,PASS)

                
               CASE (2)
                  CALL WRITE_SCALAR_IN_VTU_BIN('P_G',P_G,PASS)
                  CALL WRITE_SCALAR_IN_VTU_BIN('P_S',P_S,PASS)


               CASE (3)
                  CALL WRITE_VECTOR_IN_VTU_BIN('Gas_Velocity',U_G,V_G,W_G,PASS)

               
               CASE (4)
                  DO M = 1,MMAX
                     WRITE(SUBM,*)M
                     CALL WRITE_VECTOR_IN_VTU_BIN('Solids_Velocity_'//ADJUSTL(SUBM),U_S(:,M),V_S(:,M),W_S(:,M),PASS)
                  END DO

               
               CASE (5)
                  DO M = 1,MMAX
                     WRITE(SUBM,*)M
                     CALL WRITE_SCALAR_IN_VTU_BIN('Solids_density_'//ADJUSTL(SUBM),ROP_S(:,M),PASS)
                  END DO

               
               CASE (6)
                  CALL WRITE_SCALAR_IN_VTU_BIN('Gas_temperature',T_g,PASS)
                  DO M = 1,MMAX
                     WRITE(SUBM,*)M
                     CALL WRITE_SCALAR_IN_VTU_BIN('Solids_temperature_'//ADJUSTL(SUBM),T_S(:,M),PASS)
                  END DO

               
               CASE (7)
                  SPECIES_COUNTER = 0
                  DO N = 1,NMAX(0)
                     WRITE(SUBN,*)N
                     IF(USE_RRATES) THEN
                        SPECIES_COUNTER = SPECIES_COUNTER + 1
                        VAR_NAME = ADJUSTL(SPECIES_NAME(SPECIES_COUNTER))
                        LT = LEN_TRIM(ADJUSTL(SPECIES_NAME(SPECIES_COUNTER)))
                     ELSE
                        VAR_NAME = ADJUSTL(SPECIES_ALIAS_g(N))
                        LT = LEN_TRIM(ADJUSTL(SPECIES_ALIAS_g(N)))
                     ENDIF
                     VAR_NAME = VAR_NAME(1:LT)//'_Gas_mass_fractions_'//ADJUSTL(SUBN)
                     CALL WRITE_SCALAR_IN_VTU_BIN(VAR_NAME,X_g(:,N),PASS)
                  END DO

                  DO M = 1, MMAX 
                     WRITE(SUBM,*)M
                     DO N = 1,NMAX(M)
                        WRITE(SUBN,*)N
                        IF(USE_RRATES) THEN
                           SPECIES_COUNTER = SPECIES_COUNTER + 1
                           VAR_NAME = ADJUSTL(SPECIES_NAME(SPECIES_COUNTER))
                           LT = LEN_TRIM(ADJUSTL(SPECIES_NAME(SPECIES_COUNTER)))
                        ELSE
                           VAR_NAME = ADJUSTL(SPECIES_ALIAS_s(M,N))
                           LT = LEN_TRIM(ADJUSTL(SPECIES_ALIAS_s(M,N)))
                        ENDIF
                        VAR_NAME = VAR_NAME(1:LT)//'_Solids_mass_fractions_'//TRIM(ADJUSTL(SUBM))//'_'//ADJUSTL(SUBN)
                     CALL WRITE_SCALAR_IN_VTU_BIN(VAR_NAME,X_s(:,M,N),PASS)
                     END DO
                  END DO  

               
               CASE (8)
                  DO M = 1,MMAX
                     WRITE(SUBM,*)M
                     CALL WRITE_SCALAR_IN_VTU_BIN('Granular_temperature_'//ADJUSTL(SUBM),Theta_m(:,M),PASS)
                  END DO

               
               CASE (9)
                  DO N = 1,NSCALAR
                     WRITE(SUBN,*)N
                     VAR_NAME = 'Scalar_'//ADJUSTL(SUBN)
                     CALL WRITE_SCALAR_IN_VTU_BIN(VAR_NAME,Scalar(:,N),PASS)
                  END DO


               CASE (10)
                  DO R = 1,nRR
                     WRITE(SUBR,*)R
                     VAR_NAME = 'RRates_'//ADJUSTL(SUBR)
                     CALL WRITE_SCALAR_IN_VTU_BIN(VAR_NAME,ReactionRates(:, R),PASS)
                  END DO


               CASE (11)
                  IF(K_EPSILON) THEN
                     CALL WRITE_SCALAR_IN_VTU_BIN('K_Turb_G',K_Turb_G,PASS)                
                     CALL WRITE_SCALAR_IN_VTU_BIN('E_Turb_G',E_Turb_G,PASS)                                

                  ENDIF

               CASE (12)
                  CALL CALC_VORTICITY

                  CALL WRITE_SCALAR_IN_VTU_BIN('VORTICITY_MAG',VORTICITY,PASS)                
                  CALL WRITE_SCALAR_IN_VTU_BIN('LAMBDA_2',LAMBDA2,PASS)                
                

               
               CASE (100)
                  IF(DISCRETE_ELEMENT.AND.MPPIC) THEN
                     ALLOCATE( COUNT_DES_BC(DIMENSION_3))
                     DO IJK = IJKSTART3, IJKEND3
                        COUNT_DES_BC (IJK) = 0.d0
                        COUNT_DES_BC (IJK) = DES_CELLWISE_BCDATA(IJK)%COUNT_DES_BC
                     ENDDO
                     CALL WRITE_SCALAR_IN_VTU_BIN('COUNT_BC',COUNT_DES_BC,PASS)
                     DEALLOCATE(COUNT_DES_BC)
                  ELSE
                     CALL WRITE_SCALAR_IN_VTU_BIN('PARTITION',PARTITION,PASS)
                  ENDIF


               CASE (101)
                  Allocate(DP_BC_ID(DIMENSION_3))
                  DP_BC_ID = DFLOAT(BC_ID)
                  CALL WRITE_SCALAR_IN_VTU_BIN('BC_ID',DP_BC_ID,PASS)

                  DeAllocate(DP_BC_ID)

               CASE (102)
                  CALL WRITE_SCALAR_IN_VTU_BIN('DISTANCE_TO_WALL',DWALL,PASS)
    
               CASE (103)
                  IF(DISCRETE_ELEMENT.AND.USE_STL) THEN
                     ALLOCATE(FACET_COUNT_DES(DIMENSION_3))
                     
                     DO IJK = IJKSTART3, IJKEND3 
                        FACET_COUNT_DES(IJK) = LIST_FACET_AT_DES(IJK)%COUNT_FACETS 
                     ENDDO

                     CALL WRITE_SCALAR_IN_VTU_BIN('FACET_COUNT', FACET_COUNT_DES,PASS)
                     DEALLOCATE(FACET_COUNT_DES)
                  ENDIF

               CASE (104)
                  IF(DISCRETE_ELEMENT.AND.USE_STL) THEN
                     ALLOCATE(NEIGHBORING_FACET(DIMENSION_3))
                     
                     DO IJK = IJKSTART3, IJKEND3 
                        NEIGHBORING_FACET(IJK) = 1.0
                        IF(NO_NEIGHBORING_FACET_DES(IJK)) NEIGHBORING_FACET(IJK) = 0.0
                     ENDDO

                     CALL WRITE_SCALAR_IN_VTU_BIN('NEIGH_FACET', NEIGHBORING_FACET,PASS)
                     DEALLOCATE(NEIGHBORING_FACET)
                  ENDIF


               CASE(999)
                  Allocate(IJK_ARRAY(DIMENSION_3))
                  DO IJK = IJKSTART3, IJKEND3
                     IJK_ARRAY(IJK) = DFLOAT(IJK)
                  ENDDO
                  CALL WRITE_SCALAR_IN_VTU_BIN('IJK',IJK_ARRAY,PASS)
                  DeAllocate(IJK_ARRAY)

               CASE(1000)
                  CALL WRITE_VECTOR_IN_VTU_BIN('Scalar normal',NORMAL_S(:,1),NORMAL_S(:,2),NORMAL_S(:,3),PASS)

               CASE (1001)
                  CALL WRITE_SCALAR_IN_VTU_BIN('DEBUG_1',DEBUG_CG(:,1),PASS)

               CASE (1002)
                  CALL WRITE_SCALAR_IN_VTU_BIN('DEBUG_2',DEBUG_CG(:,2),PASS)

               CASE (1003)
                  CALL WRITE_SCALAR_IN_VTU_BIN('DEBUG_3',DEBUG_CG(:,3),PASS)

               CASE (1004)
                  CALL WRITE_SCALAR_IN_VTU_BIN('DEBUG_4',DEBUG_CG(:,4),PASS)

               CASE (1005)
                  CALL WRITE_SCALAR_IN_VTU_BIN('DEBUG_5',DEBUG_CG(:,5),PASS)

               CASE (1006)
                  CALL WRITE_SCALAR_IN_VTU_BIN('DEBUG_6',DEBUG_CG(:,6),PASS)

               CASE (1007)
                  CALL WRITE_SCALAR_IN_VTU_BIN('DEBUG_7',DEBUG_CG(:,7),PASS)

               CASE (1008)
                  CALL WRITE_SCALAR_IN_VTU_BIN('DEBUG_8',DEBUG_CG(:,8),PASS)

               CASE (1009)
                  CALL WRITE_SCALAR_IN_VTU_BIN('DEBUG_9',DEBUG_CG(:,9),PASS)

               CASE (1010)
                  CALL WRITE_SCALAR_IN_VTU_BIN('DEBUG_10',DEBUG_CG(:,10),PASS)

               CASE (1011)
                  CALL WRITE_SCALAR_IN_VTU_BIN('DEBUG_11',DEBUG_CG(:,11),PASS)

               CASE (1012)
                  CALL WRITE_SCALAR_IN_VTU_BIN('DEBUG_12',DEBUG_CG(:,12),PASS)

               CASE (1013)
                  CALL WRITE_SCALAR_IN_VTU_BIN('DEBUG_13',DEBUG_CG(:,13),PASS)

               CASE (1014)
                  CALL WRITE_SCALAR_IN_VTU_BIN('DEBUG_14',DEBUG_CG(:,14),PASS)

               CASE (1015)
                  CALL WRITE_SCALAR_IN_VTU_BIN('DEBUG_15',DEBUG_CG(:,15),PASS)



               CASE (0) ! do nothing

               CASE (UNDEFINED_I) ! do nothing

               CASE DEFAULT

                  WRITE(*,30) ' Unknown VTK variable flag ',L,':',VTK_VAR(L)
                  WRITE(*,30) ' Available flags are : '
                  WRITE(*,30) ' 1 : Void fraction (EP_g)'
                  WRITE(*,30) ' 2 : Gas pressure, solids pressure (P_g, P_star)'
                  WRITE(*,30) ' 3 : Gas velocity (U_g, V_g, W_g)'
                  WRITE(*,30) ' 4 : Solids velocity (U_s, V_s, W_s)'
                  WRITE(*,30) ' 5 : Solids density (ROP_s)'
                  WRITE(*,30) ' 6 : Gas and solids temperature (T_g, T_s1, T_s2)'
                  WRITE(*,30) ' 7 : Gas and solids mass fractions (X_g, X-s)'
                  WRITE(*,30) ' 8 : Granular temperature (G)'
                  write(*,30) ' 9 : User defined scalars'
                  write(*,30) '10 : Reaction Rates'
                  write(*,30) '11 : Turbulence quantities (k and ε)'
                  write(*,30) '12 : Gas Vorticity magnitude and Lambda_2 (VORTICITY, LAMBDA_2)'
                  write(*,30) '100: Processor assigned to scalar cell (Partition)'
                  write(*,30) '101: Boundary condition flag for scalar cell (BC_ID)'
                  write(*,30) 'MFiX will exit now.'
                  CALL MFIX_EXIT(myPE) 

               END SELECT

         END DO


      ENDDO ! PASS LOOP, EITHER HEADER OR DATA


      CALL CLOSE_VTU_FILE_BIN
      CALL UPDATE_AND_CLOSE_PVD_FILE

      call MPI_barrier(MPI_COMM_WORLD,mpierr)

      IF (myPE == PE_IO.AND.TIME_DEPENDENT_FILENAME) THEN
         OPEN(UNIT = VTU_FRAME_UNIT, FILE = TRIM(VTU_FRAME_FILENAME))
         WRITE(VTU_FRAME_UNIT,*)FRAME
         CLOSE(VTU_FRAME_UNIT)
      ENDIF


     IF (FULL_LOG.AND.myPE == PE_IO) WRITE(*,20)' DONE.'

10    FORMAT(A,$)
20    FORMAT(A,1X/)
30    FORMAT(1X,A)    
      RETURN
      
      END SUBROUTINE WRITE_VTU_FILE 



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: OPEN_VTU_FILE                                          C
!  Purpose: Open a vtu file and writes the header                      C
!           Binary format                                              C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE OPEN_VTU_FILE_BIN
    
      USE param
      USE param1
      USE parallel
      USE constant
      USE run
      USE toleranc
      USE geometry
      USE indices  
      USE compar
      USE sendrecv
      USE output
      USE quadric
      USE cutcell
      USE fldvar
      USE vtk
      use cdist

     
      IMPLICIT NONE
      DOUBLE PRECISION:: Xw,Xe,Yn,Ys
      INTEGER :: I,J,K,L,IM,JM,KM,IP,JP,KP,IJK
      INTEGER :: IMJK,IJMK,IJKM,IMJMK,IMJKM,IJMKM,IMJMKM
      LOGICAL :: VTU_FRAME_FILE_EXISTS
      INTEGER :: ISTAT 

      include "function.inc"

! Only open the file form head node when not using distributed I/O
      IF (myPE /= PE_IO.AND.(.NOT.BDIST_IO)) RETURN 


      IF(TIME_DEPENDENT_FILENAME) THEN
         INQUIRE(FILE=VTU_FRAME_FILENAME,EXIST=VTU_FRAME_FILE_EXISTS)
         IF(VTU_FRAME_FILE_EXISTS) THEN
            OPEN(UNIT = VTU_FRAME_UNIT, FILE = TRIM(VTU_FRAME_FILENAME))
            READ(VTU_FRAME_UNIT,*)FRAME
            CLOSE(VTU_FRAME_UNIT)
         ENDIF
         IF(RESET_FRAME_AT_TIME_ZERO.AND.TIME==ZERO) FRAME=-1
         FRAME = FRAME + 1 
      ENDIF




      IF (BDIST_IO) THEN 
! For distributed I/O, define the file name for each processor
         IF(TIME_DEPENDENT_FILENAME) THEN
            WRITE(VTU_FILENAME,20) TRIM(RUN_NAME),FRAME,MYPE
         ELSE
            WRITE(VTU_FILENAME,25) TRIM(RUN_NAME),MYPE
         ENDIF
      ELSE 
         IF(MYPE.EQ.PE_IO) THEN 
            IF(TIME_DEPENDENT_FILENAME) THEN
               WRITE(VTU_FILENAME,30) TRIM(RUN_NAME),FRAME
            ELSE
               WRITE(VTU_FILENAME,35) TRIM(RUN_NAME)
            ENDIF
         END IF  
      END IF 


! Add the VTU directory path if necessary

      IF(TRIM(VTU_DIR)/='.') VTU_FILENAME='./'//TRIM(VTU_DIR)//'/'//VTU_FILENAME

! Echo 
      IF (FULL_LOG) THEN
         IF (.NOT.BDIST_IO) THEN
            WRITE(*,10)' WRITING VTU FILE : ', TRIM(VTU_FILENAME),' .'
         ELSE
            IF(myPE==PE_IO) WRITE(*,15)' EACH PROCESOR IS WRITING ITS OWN VTU FILE.'
         ENDIF
      ENDIF

! Open File
!      OPEN(UNIT = VTU_UNIT, FILE = TRIM(VTU_FILENAME),FORM='BINARY',IOSTAT=ISTAT)


      OPEN(UNIT     = VTU_UNIT,           &
           FILE     = TRIM(VTU_FILENAME), &
           FORM     = 'UNFORMATTED',      &  ! works with gfortran 4.3.4 and ifort 10.1 but may not be supported by all compilers
                                             ! use 'BINARY' if 'UNFORMATTED' is not supported 
           ACCESS   = 'STREAM',           &  ! works with gfortran 4.3.4 and ifort 10.1 but may not be supported by all compilers
                                             ! use 'SEQUENTIAL' if 'STREAM' is not supported 
           ACTION   = 'WRITE',            &
           CONVERT  = 'BIG_ENDIAN',       &
           IOSTAT=ISTAT)


      IF (ISTAT /= 0) THEN
         IF(DMP_LOG) WRITE(UNIT_LOG, 1001) VTU_FILENAME, VTU_UNIT,VTU_DIR
         IF(FULL_LOG.AND.myPE == PE_IO) WRITE(*, 1001) VTU_FILENAME, VTU_UNIT,VTU_DIR
         CALL MFIX_EXIT(myPE)
      ENDIF


 1001 FORMAT(/1X,70('*')//, ' From: OPEN_VTU_FILE',/,' Message: ',          &
         'Error opening vtu file. Terminating run.',/10X,          &
         'File name:  ',A,/10X,                                         &
         'DES_UNIT :  ',i4, /10X,                                       &
         'PLEASE VERIFY THAT VTU_DIR EXISTS: ', A, &
         /1X,70('*')/)


! Write file Header
      BUFFER='<?xml version="1.0"?>'
      WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC


      WRITE(BUFFER,110)'<!-- Time =',TIME,' sec. -->'
      WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

      BUFFER='<VTKFile type="UnstructuredGrid" version="0.1" byte_order="BigEndian">'
      WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

      BUFFER='  <UnstructuredGrid>'
      WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC






! For distributed I/O, open .pvtu file that combines all *.vtu files for a given FRAME
! this is a simple ASCII file

      IF (myPE == PE_IO.AND.BDIST_IO) THEN    

         IF(TIME_DEPENDENT_FILENAME) THEN
            WRITE(PVTU_FILENAME,40) TRIM(RUN_NAME),FRAME
         ELSE
            WRITE(PVTU_FILENAME,45) TRIM(RUN_NAME)
         ENDIF

         IF(TRIM(VTU_DIR)/='.') PVTU_FILENAME='./'//TRIM(VTU_DIR)//'/'//PVTU_FILENAME

         OPEN(UNIT = PVTU_UNIT, FILE = TRIM(PVTU_FILENAME))

         WRITE(PVTU_UNIT,100) '<?xml version="1.0"?>'
         WRITE(PVTU_UNIT,110) '<!-- Time =',TIME,' sec. -->'
         WRITE(PVTU_UNIT,120) '<VTKFile type="PUnstructuredGrid"',&
                  ' version="0.1" byte_order="BigEndian">'

         WRITE(PVTU_UNIT,100) '  <PUnstructuredGrid GhostLevel="0">'
         WRITE(PVTU_UNIT,100) '      <PPoints>'
         WRITE(PVTU_UNIT,100) '        <PDataArray type="Float32" Name="coordinates" NumberOfComponents="3" format="appended" offset=" 0" />'
         WRITE(PVTU_UNIT,100) '      </PPoints>'
         WRITE(PVTU_UNIT,100) ''
         WRITE(PVTU_UNIT,100) '      <PCellData Scalars="scalars">'


      ENDIF



100   FORMAT(A)
110   FORMAT(A,E14.8,A)
120   FORMAT(A,A)
10    FORMAT(/1X,3A,$)
15    FORMAT(/1X,A,$)
20    FORMAT(A,"_",I4.4,"_",I5.5,".vtu")
25    FORMAT(A,"_",I5.5,".vtu")
30    FORMAT(A,"_",I4.4,".vtu")
35    FORMAT(A,".vtu")
40    FORMAT(A,"_",I4.4,".pvtu")
45    FORMAT(A,".pvtu")

      RETURN

      END SUBROUTINE OPEN_VTU_FILE_BIN





!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_GEOMETRY_IN_VTU_BIN                              C
!  Purpose: Write Geometry and connectivity in a vtu file              C
!           Binary format                                              C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE WRITE_GEOMETRY_IN_VTU_BIN(PASS)
    
      USE param
      USE param1
      USE parallel
      USE constant
      USE run
      USE toleranc
      USE geometry
      USE indices  
      USE compar
      USE mpi_utility 
      USE sendrecv
      USE quadric
      USE cutcell
      USE fldvar
      USE vtk
      USE cdist

     
      IMPLICIT NONE

      INTEGER :: IJK,I,J,K,L
      INTEGER :: IJK_OFFSET,OFFSET

      INTEGER :: iproc,IERR
      INTEGER, DIMENSION(0:numPEs-1) :: disp,rcount
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: SHIFTED_CONNECTIVITY
      INTEGER :: CELL_TYPE

      REAL*4 :: float,SP_X,SP_Y,SP_Z
      INTEGER :: int

      INTEGER ::     nbytes_xyz,nbytes_connectivity,nbytes_offset,nbytes_type 
      INTEGER ::     offset_xyz,offset_connectivity,offset_offset,offset_type 

      INTEGER :: PASS
      INTEGER :: WRITE_HEADER = 1
      INTEGER :: WRITE_DATA   = 2

      include "function.inc"


! First a series of tags is written for the geometry (PASS=WRITE_HEADER)
!  - Coordinates
!  - Connectivity
!  - Connectivity offset
!  - cell types
!
! Since the data is appended (i.e., written after all tags), the offset, in number of bytes must be specified.
! The offset includes the size of the data for each field, plus the size of the integer that stores the number of bytes.
! this is why the offset of a field equals the offset of the previous field plus sizeof(int) plus the number of bytes of the field.

! Next, the actual data is written for the geometry (PASS=WRITE_DATA)
! The DATA is converted to single precision to save memory.

      IF (myPE == PE_IO.AND.(.NOT.BDIST_IO)) THEN

         NUMBER_OF_VTK_CELLS = NUMBER_OF_CELLS - NUMBER_OF_BLOCKED_CELLS

! Number of bytes of each field
         nbytes_xyz          = NUMBER_OF_POINTS * 3 * sizeof(float)

         nbytes_connectivity = 0
         DO IJK = 1,IJKMAX3
            IF (GLOBAL_INTERIOR_CELL_AT(IJK))      THEN
               IF (.NOT.GLOBAL_BLOCKED_CELL_AT(IJK)) THEN
                  nbytes_connectivity = nbytes_connectivity + GLOBAL_NUMBER_OF_NODES(IJK)
               ENDIF
            ENDIF
         END DO
         nbytes_connectivity = nbytes_connectivity * sizeof(int)

         
         nbytes_offset       = NUMBER_OF_VTK_CELLS * sizeof(int)

         nbytes_type         = NUMBER_OF_VTK_CELLS * sizeof(int)


! Offset of each field
         offset_xyz = 0
         offset_connectivity = offset_xyz          + sizeof(int) + nbytes_xyz
         offset_offset       = offset_connectivity + sizeof(int) + nbytes_connectivity
         offset_type         = offset_offset       + sizeof(int) + nbytes_offset


         IF(PASS==WRITE_HEADER) THEN

            WRITE(BUFFER,100)'    <Piece NumberOfPoints="',NUMBER_OF_POINTS,'" NumberOfCells="',NUMBER_OF_VTK_CELLS,'" >'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

            WRITE(BUFFER,110)'      <Points>'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

            WRITE(BUFFER,100)'        <DataArray type="Float32" Name="coordinates" NumberOfComponents="3" format="appended" offset="',offset_xyz,'" />'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

            WRITE(BUFFER,110)'      </Points>'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

            WRITE(BUFFER,110)'      <Cells>'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

            WRITE(BUFFER,100)'        <DataArray type="Int32" Name="connectivity" format="appended" offset="',offset_connectivity,'" />'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

            WRITE(BUFFER,100)'        <DataArray type="Int32" Name="offsets" format="appended" offset="',offset_offset,'" />'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

            WRITE(BUFFER,100)'        <DataArray type="Int32" Name="types" format="appended" offset="',offset_type,'" />'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

            WRITE(BUFFER,110)'      </Cells>'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

            VTU_offset =  offset_type       + sizeof(int) + nbytes_type  ! Store offset for first variable to be written

            WRITE(BUFFER,110)'      <CellData>'                          ! Preparing CellData
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC




         ELSEIF(PASS==WRITE_DATA) THEN

            WRITE(BUFFER,110)'      </CellData>'                          
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

            WRITE(BUFFER,110)'    </Piece>'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC


            WRITE(BUFFER,110)'  </UnstructuredGrid>'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC


            WRITE(BUFFER,110)'  <AppendedData encoding="raw">'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC


! Starting raw binary data with an underscore

            WRITE(BUFFER,110)'_'
            WRITE(VTU_UNIT)TRIM(BUFFER)

! X,Y,Z coordinates
            WRITE(VTU_UNIT) nbytes_xyz,(REAL(XG_E(GLOBAL_I_OF(IJK))),REAL(YG_N(GLOBAL_J_OF(IJK))),REAL(ZG_T(GLOBAL_K_OF(IJK))), IJK = 1,IJKMAX3), &
                                       (REAL(GLOBAL_X_NEW_POINT(IJK)),REAL(GLOBAL_Y_NEW_POINT(IJK)),REAL(GLOBAL_Z_NEW_POINT(IJK)),IJK = 1,&
                                        GLOBAL_NUMBER_OF_NEW_POINTS)

! Conectivity
            WRITE(VTU_UNIT) nbytes_connectivity

            DO IJK = 1,IJKMAX3
               IF (GLOBAL_INTERIOR_CELL_AT(IJK))      THEN
                  IF (.NOT.GLOBAL_BLOCKED_CELL_AT(IJK)) WRITE(VTU_UNIT) (GLOBAL_CONNECTIVITY(IJK,L)-1,L=1,GLOBAL_NUMBER_OF_NODES(IJK))
               ENDIF
            END DO

! Offsets
            WRITE(VTU_UNIT) nbytes_offset

            OFFSET = 0
            DO IJK = 1,IJKMAX3
               IF (GLOBAL_INTERIOR_CELL_AT(IJK))      THEN
                  IF (.NOT.GLOBAL_BLOCKED_CELL_AT(IJK)) THEN
                     OFFSET = OFFSET + GLOBAL_NUMBER_OF_NODES(IJK)
                     WRITE(VTU_UNIT) OFFSET
                  ENDIF
               ENDIF
            END DO

! Types
            WRITE(VTU_UNIT)nbytes_type

            IF(NO_K) THEN
               CELL_TYPE = 7
            ELSE
               CELL_TYPE = 41
            ENDIF

            DO IJK = 1,IJKMAX3
               IF (GLOBAL_INTERIOR_CELL_AT(IJK))      THEN
                  IF (.NOT.GLOBAL_BLOCKED_CELL_AT(IJK)) WRITE(VTU_UNIT) CELL_TYPE
               ENDIF
            END DO


         ENDIF


      ELSEIF(BDIST_IO) THEN

! For distributed I/O, it works the same as above, except, the data is local to each processor
! First compute local number of cells and points

         NUMBER_OF_CELLS = 0
         NUMBER_OF_BLOCKED_CELLS  = 0

         DO IJK = 1, IJKEND3   
            IF (INTERIOR_CELL_AT(IJK)) THEN
               NUMBER_OF_CELLS = NUMBER_OF_CELLS + 1
               IF (BLOCKED_CELL_AT(IJK))  NUMBER_OF_BLOCKED_CELLS  = NUMBER_OF_BLOCKED_CELLS + 1
            ENDIF
         END DO


         NUMBER_OF_POINTS = IJKEND3 + NUMBER_OF_NEW_POINTS
         NUMBER_OF_VTK_CELLS = NUMBER_OF_CELLS - NUMBER_OF_BLOCKED_CELLS

! Number of bytes of each field
         nbytes_xyz          = NUMBER_OF_POINTS * 3 * sizeof(float)

         nbytes_connectivity = 0
         DO IJK = 1,IJKEND3
            IF (INTERIOR_CELL_AT(IJK))      THEN
               IF (.NOT.BLOCKED_CELL_AT(IJK)) THEN
                  nbytes_connectivity = nbytes_connectivity + NUMBER_OF_NODES(IJK)
               ENDIF
            ENDIF
         END DO
         nbytes_connectivity = nbytes_connectivity * sizeof(int)
         
         nbytes_offset       = NUMBER_OF_VTK_CELLS * sizeof(int)

         nbytes_type         = NUMBER_OF_VTK_CELLS * sizeof(int)


! Offset of each field
         offset_xyz = 0
         offset_connectivity = offset_xyz          + sizeof(int) + nbytes_xyz
         offset_offset       = offset_connectivity + sizeof(int) + nbytes_connectivity
         offset_type         = offset_offset       + sizeof(int) + nbytes_offset



         IF(PASS==WRITE_HEADER) THEN

            WRITE(BUFFER,100)'    <Piece NumberOfPoints="',NUMBER_OF_POINTS,'" NumberOfCells="',NUMBER_OF_VTK_CELLS,'" >'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

            WRITE(BUFFER,110)'      <Points>'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

            WRITE(BUFFER,100)'        <DataArray type="Float32" Name="coordinates" NumberOfComponents="3" format="appended" offset="',offset_xyz,'" />'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

            WRITE(BUFFER,110)'      </Points>'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC


            WRITE(BUFFER,110)'      <Cells>'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

            WRITE(BUFFER,100)'        <DataArray type="Int32" Name="connectivity" format="appended" offset="',offset_connectivity,'" />'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

            WRITE(BUFFER,100)'        <DataArray type="Int32" Name="offsets" format="appended" offset="',offset_offset,'" />'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

            WRITE(BUFFER,100)'        <DataArray type="Int32" Name="types" format="appended" offset="',offset_type,'" />'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

            WRITE(BUFFER,110)'      </Cells>'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

            VTU_offset =  offset_type       + sizeof(int) + nbytes_type  ! Store offset for first variable to be written

            WRITE(BUFFER,110)'      <CellData>'                          ! Preparing CellData
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC




         ELSEIF(PASS==WRITE_DATA) THEN

            WRITE(BUFFER,110)'      </CellData>'                          
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

            WRITE(BUFFER,110)'    </Piece>'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC


            WRITE(BUFFER,110)'  </UnstructuredGrid>'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC


            WRITE(BUFFER,110)'  <AppendedData encoding="raw">'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC


! Starting raw binary data with an underscore

            WRITE(BUFFER,110)'_'
            WRITE(VTU_UNIT)TRIM(BUFFER)

! X,Y,Z coordinates
            WRITE(VTU_UNIT) nbytes_xyz,(REAL(XG_E(I_OF(IJK))),REAL(YG_N(J_OF(IJK))),REAL(ZG_T(K_OF(IJK))), IJK = 1,IJKEND3), &
                                       (REAL(X_NEW_POINT(IJK)),REAL(Y_NEW_POINT(IJK)),REAL(Z_NEW_POINT(IJK)),IJK = 1,&
                                        NUMBER_OF_NEW_POINTS)

! Conectivity
            WRITE(VTU_UNIT) nbytes_connectivity

            DO IJK = 1,IJKEND3
               IF (INTERIOR_CELL_AT(IJK))      THEN
                  IF (.NOT.BLOCKED_CELL_AT(IJK)) WRITE(VTU_UNIT) (CONNECTIVITY(IJK,L)-1,L=1,NUMBER_OF_NODES(IJK))
               ENDIF
            END DO

! Offsets
            WRITE(VTU_UNIT) nbytes_offset

            OFFSET = 0
            DO IJK = 1,IJKEND3
               IF (INTERIOR_CELL_AT(IJK))      THEN
                  IF (.NOT.BLOCKED_CELL_AT(IJK)) THEN
                     OFFSET = OFFSET + NUMBER_OF_NODES(IJK)
                     WRITE(VTU_UNIT) OFFSET
                  ENDIF
               ENDIF
            END DO

! Types
            WRITE(VTU_UNIT)nbytes_type

            IF(NO_K) THEN
               CELL_TYPE = 7
            ELSE
               CELL_TYPE = 41
            ENDIF

            DO IJK = 1,IJKEND3
               IF (INTERIOR_CELL_AT(IJK))      THEN
                  IF (.NOT.BLOCKED_CELL_AT(IJK)) WRITE(VTU_UNIT) CELL_TYPE
               ENDIF
            END DO


         ENDIF


      ENDIF


100   FORMAT(A,I12,A,I12,A)
110   FORMAT(A)
120   FORMAT(10X,3(E16.8E3,2X))
130   FORMAT(10X,15(I12,2X))

      RETURN
      
      END SUBROUTINE WRITE_GEOMETRY_IN_VTU_BIN


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_SCALAR_IN_VTU_BIN                                C
!  Purpose: Write Scalar variable in a vtu file                        C
!           Binary format                                              C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE WRITE_SCALAR_IN_VTU_BIN(VAR_NAME,VAR,PASS)
    
      USE param
      USE param1
      USE parallel
      USE constant
      USE run
      USE toleranc
      USE geometry
      USE indices  
      USE compar
      USE mpi_utility 
      USE sendrecv
      USE quadric
      USE cutcell
      USE fldvar
      USE vtk
      USE cdist
      USE output

      IMPLICIT NONE
      INTEGER :: I,IJK,L

      CHARACTER (*) :: VAR_NAME
      DOUBLE PRECISION, DIMENSION(DIMENSION_3) ::  VAR
      DOUBLE PRECISION, ALLOCATABLE :: GLOBAL_VAR(:)


      INTEGER :: int
      REAL*4 :: float

      INTEGER :: nbytes_scalar

      INTEGER :: PASS
      INTEGER :: WRITE_HEADER = 1
      INTEGER :: WRITE_DATA   = 2


      include "function.inc"

      IF (.NOT.BDIST_IO) THEN

! For each scalar, write a tag, with corresponding offset

         nbytes_scalar = NUMBER_OF_VTK_CELLS * sizeof(float)
 
         IF(PASS==WRITE_HEADER) THEN
!           For each scalar, write a tag, with corresponding offset

            DO I = 1,LEN_TRIM(VAR_NAME)
               IF(VAR_NAME(I:I) == ' ') VAR_NAME(I:I) = '_'
            ENDDO

            WRITE(BUFFER,90)'        <DataArray type="Float32" Name="',TRIM(VAR_NAME),'" format="appended" offset="',VTU_offset,'" />'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

            VTU_offset = VTU_offset + sizeof(float) + nbytes_scalar


         ELSEIF(PASS==WRITE_DATA) THEN
!           and write the data, always preceded by its size in number of bytes

            IF (myPE == PE_IO) THEN
               allocate (GLOBAL_VAR(ijkmax3))     
            ELSE
               allocate (GLOBAL_VAR(1))     
            ENDIF

            call gather (VAR,GLOBAL_VAR,root) 


            IF (myPE /= PE_IO) RETURN


            WRITE(VTU_UNIT) nbytes_scalar

            DO IJK = 1,IJKMAX3
               IF (GLOBAL_INTERIOR_CELL_AT(IJK))      THEN
                  IF (.NOT.GLOBAL_BLOCKED_CELL_AT(IJK))   WRITE(VTU_UNIT) REAL(GLOBAL_VAR(IJK))
               ENDIF
            ENDDO


            Deallocate (GLOBAL_VAR) 

         ENDIF


      ELSE ! BDIST_IO=.TRUE.


         nbytes_scalar = NUMBER_OF_VTK_CELLS * sizeof(float)
 
         IF(PASS==WRITE_HEADER) THEN
!           For each scalar, write a tag, with corresponding offset

            DO I = 1,LEN_TRIM(VAR_NAME)
               IF(VAR_NAME(I:I) == ' ') VAR_NAME(I:I) = '_'
            ENDDO

            WRITE(BUFFER,90)'        <DataArray type="Float32" Name="',TRIM(VAR_NAME),'" format="appended" offset="',VTU_offset,'" />'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC


            VTU_offset = VTU_offset + sizeof(float) + nbytes_scalar


         ELSEIF(PASS==WRITE_DATA) THEN
!           and write the data, always preceded by its size in number of bytes

            WRITE(VTU_UNIT) nbytes_scalar

            DO IJK = 1,IJKEND3
               IF (INTERIOR_CELL_AT(IJK))      THEN
                  IF (.NOT.BLOCKED_CELL_AT(IJK))   WRITE(VTU_UNIT) REAL(VAR(IJK))
               ENDIF
            ENDDO

         ENDIF


         IF (myPE == PE_IO) THEN       ! Update pvtu file with variable name
            WRITE(PVTU_UNIT,90) '        <DataArray type="Float32" Name="',TRIM(VAR_NAME),'" format="appended" offset="',VTU_offset,'" />'
         ENDIF


      ENDIF


      IF (PASS==WRITE_DATA.AND.FULL_LOG.AND.myPE == PE_IO) WRITE(*,10)'.'

10    FORMAT(A,$)
90    FORMAT(A,A,A,I12,A)
100   FORMAT(A)
110   FORMAT(A,A,A)
120   FORMAT(10X,E16.8E3)



      RETURN
      
      END SUBROUTINE WRITE_SCALAR_IN_VTU_BIN




!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_VECTOR_IN_VTU                                    C
!  Purpose: Write Vector variable in a vtu file                        C
!           Binary format                                              C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE WRITE_VECTOR_IN_VTU_BIN(VAR_NAME,VARX,VARY,VARZ,PASS)
    
      USE param
      USE param1
      USE parallel
      USE constant
      USE run
      USE toleranc
      USE geometry
      USE indices  
      USE compar
      USE mpi_utility 
      USE sendrecv
      USE quadric
      USE cutcell
      USE fldvar
      USE vtk
      USE cdist
      USE output
      
      IMPLICIT NONE
      INTEGER :: IJK,L

      CHARACTER (*) :: VAR_NAME
      DOUBLE PRECISION, DIMENSION(DIMENSION_3) ::  VARX,VARY,VARZ
      DOUBLE PRECISION, ALLOCATABLE :: GLOBAL_VARX(:),GLOBAL_VARY(:),GLOBAL_VARZ(:)


      INTEGER :: int
      REAL*4 :: float

      INTEGER :: nbytes_vector

      INTEGER :: PASS
      INTEGER :: WRITE_HEADER = 1
      INTEGER :: WRITE_DATA   = 2


      include "function.inc"

      IF (.NOT.BDIST_IO) THEN

         nbytes_vector = NUMBER_OF_VTK_CELLS * 3 * sizeof(float)
 
         IF(PASS==WRITE_HEADER) THEN
!           For each vector, write a tag, with corresponding offset

            WRITE(BUFFER,90)'        <DataArray type="Float32" Name="',TRIM(VAR_NAME),'"  NumberOfComponents="3" format="appended" offset="',VTU_offset,'" />'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC


            VTU_offset = VTU_offset + sizeof(float) + nbytes_vector


         ELSEIF(PASS==WRITE_DATA) THEN
!           and write the data, always preceded by its size in number of bytes

            IF (myPE == PE_IO) THEN
               allocate (GLOBAL_VARX(ijkmax3))
               allocate (GLOBAL_VARY(ijkmax3))     
               allocate (GLOBAL_VARZ(ijkmax3))          
            ELSE
               allocate (GLOBAL_VARX(1))
               allocate (GLOBAL_VARY(1))     
               allocate (GLOBAL_VARZ(1))          
            ENDIF

            call gather (VARX,GLOBAL_VARX,root)
            call gather (VARY,GLOBAL_VARY,root) 
            call gather (VARZ,GLOBAL_VARZ,root)  

            IF (myPE /= PE_IO) RETURN


            WRITE(VTU_UNIT) nbytes_vector

            DO IJK = 1,IJKMAX3
               IF (GLOBAL_INTERIOR_CELL_AT(IJK))      THEN
                  IF (.NOT.GLOBAL_BLOCKED_CELL_AT(IJK))   WRITE(VTU_UNIT) REAL(GLOBAL_VARX(IJK)),REAL(GLOBAL_VARY(IJK)),REAL(GLOBAL_VARZ(IJK))
               ENDIF
            ENDDO


            Deallocate (GLOBAL_VARX)
            Deallocate (GLOBAL_VARY)   
            Deallocate (GLOBAL_VARZ) 

         ENDIF


      ELSE ! BDIST_IO=.TRUE.


         nbytes_vector = NUMBER_OF_VTK_CELLS * 3 * sizeof(float)
 
         IF(PASS==WRITE_HEADER) THEN
!           For each vector, write a tag, with corresponding offset


            WRITE(BUFFER,90)'        <DataArray type="Float32" Name="',TRIM(VAR_NAME),'"  NumberOfComponents="3" format="appended" offset="',VTU_offset,'" />'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC


            VTU_offset = VTU_offset + sizeof(float) + nbytes_vector


         ELSEIF(PASS==WRITE_DATA) THEN
!           and write the data, always preceded by its size in number of bytes

            WRITE(VTU_UNIT) nbytes_vector

            DO IJK = 1,IJKEND3
               IF (INTERIOR_CELL_AT(IJK))      THEN
                  IF (.NOT.BLOCKED_CELL_AT(IJK))   WRITE(VTU_UNIT) REAL(VARX(IJK)),REAL(VARY(IJK)),REAL(VARZ(IJK))
               ENDIF
            ENDDO

         ENDIF


         IF (myPE == PE_IO) THEN       ! Update pvtu file with variable name
            WRITE(PVTU_UNIT,90)'        <DataArray type="Float32" Name="',TRIM(VAR_NAME),'"  NumberOfComponents="3" format="appended" offset="',VTU_offset,'" />'
         ENDIF

      ENDIF


      IF (PASS==WRITE_DATA.AND.FULL_LOG.AND.myPE == PE_IO) WRITE(*,10)'.'

10    FORMAT(A,$)
90    FORMAT(A,A,A,I12,A)
100   FORMAT(A)
110   FORMAT(A,A,A)
120   FORMAT(10X,3(E16.8E3,2X))


      RETURN
      
      END SUBROUTINE WRITE_VECTOR_IN_VTU_BIN




!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CLOSE_VTU_FILE_BIN                                     C
!  Purpose: Close a vtu file                                           C
!           Binary format                                              C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE CLOSE_VTU_FILE_BIN

      USE compar
      Use run
      USE vtk
      use cdist

      IMPLICIT NONE
    
      INTEGER:: N
      CHARACTER (LEN=32)  :: VTU_NAME

      IF (myPE /= PE_IO.AND.(.NOT.BDIST_IO)) RETURN 

! Write last tags and close the vtu file
       WRITE(BUFFER,110)'  </AppendedData>'
       WRITE(VTU_UNIT)END_REC//TRIM(BUFFER)//END_REC

       WRITE(BUFFER,110)'</VTKFile>'
       WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

      CLOSE(VTU_UNIT)

! Update pvtu file and close
      IF (myPE == PE_IO.AND.BDIST_IO) THEN       
         WRITE(PVTU_UNIT,100) '      </PCellData>'

         DO N = 0,NumPEs-1
            IF(TIME_DEPENDENT_FILENAME) THEN
               WRITE(VTU_NAME,20) TRIM(RUN_NAME),FRAME,N
            ELSE
               WRITE(VTU_NAME,25) TRIM(RUN_NAME),N
            ENDIF

            WRITE(PVTU_UNIT,110) '      <Piece Source="',TRIM(VTU_NAME),'"/>'
         ENDDO


         WRITE(PVTU_UNIT,100) '  </PUnstructuredGrid>'
         WRITE(PVTU_UNIT,100) '</VTKFile>'
         CLOSE(PVTU_UNIT)
      ENDIF


20    FORMAT(A,"_",I4.4,"_",I5.5,".vtu")
25    FORMAT(A,"_",I5.5,".vtu")

100   FORMAT(A)
110   FORMAT(A,A,A)

      RETURN

      END SUBROUTINE CLOSE_VTU_FILE_BIN




!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_VTU_FILE                                         C
!  Purpose: Writes the cut cell grid in VTK format (Unstructured VTU)  C
!           ASCII format                                               C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE WRITE_VTU_FILE_ASCII
    
      USE param
      USE param1
      USE parallel
      USE constant
      USE run
      USE toleranc
      USE geometry
      USE indices  
      USE compar
      USE sendrecv
      USE quadric
      USE cutcell
      USE fldvar
      USE visc_s
      USE physprop
      USE pgcor
      USE vtk
      USE rxns      
      USE output
      USE scalars

      USE mpi_utility     
      USE parallel_mpi


      USE pgcor
      USE pscor
      USE discretelement, Only : DES_CELLWISE_BCDATA, DISCRETE_ELEMENT
      USE mfix_pic

      IMPLICIT NONE
      DOUBLE PRECISION:: Xw,Xe,Yn,Ys
      INTEGER :: I,J,K,L,M,N,IM,JM,KM,IP,JP,KP,IJK
      INTEGER :: IMJK,IJMK,IJKM,IMJMK,IMJKM,IJMKM,IMJMKM

      INTEGER sw,se,ne,nw
      INTEGER, DIMENSION(10) :: additional_node
      DOUBLE PRECISION, DIMENSION(2*DIMENSION_3) ::  X_OF
      DOUBLE PRECISION, DIMENSION(2*DIMENSION_3) ::  Y_OF
      DOUBLE PRECISION, DIMENSION(2*DIMENSION_3) ::  Z_OF
      INTEGER, DIMENSION(DIMENSION_3) ::  INDEX_OF_E_ADD_NODE
      INTEGER, DIMENSION(DIMENSION_3) ::  INDEX_OF_N_ADD_NODE
      INTEGER :: SPECIES_COUNTER,LT

      CHARACTER (LEN=32) :: SUBM,SUBN
      CHARACTER (LEN=64) :: VAR_NAME

      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  DP_BC_ID, COUNT_DES_BC,IJK_ARRAY

      include "function.inc"

      IF(.NOT.CARTESIAN_GRID) RETURN

      DX(IEND3+1) = DX(IEND3)
      DY(JEND3+1) = DY(JEND3)
      DZ(KEND3+1) = DZ(KEND3)

!     Location of U-momentum cells for original (uncut grid)
      IF (DO_I) THEN 
        XG_E(1) = ZERO
        DO I = IMIN1, IMAX2 
           XG_E(I) = XG_E(I-1) + DX(I) 
        END DO 
      ENDIF

!     Location of V-momentum cells for original (uncut grid)
      IF (DO_J) THEN 
        YG_N(1) = ZERO
        DO J = JMIN1, JMAX2 
           YG_N(J) = YG_N(J-1) + DY(J) 
        END DO 
      ENDIF

!     Location of W-momentum cells for original (uncut grid)
      IF (DO_K) THEN 
        ZG_T(1) = ZERO
        DO K = KMIN1, KMAX2 
           ZG_T(K) = ZG_T(K-1) + DZ(K) 
        END DO 
      ELSE
         ZG_T = ZERO
      ENDIF

      CALL OPEN_VTU_FILE_ASCII
      CALL OPEN_PVD_FILE

      CALL WRITE_GEOMETRY_IN_VTU_ASCII

      DO L = 1, DIM_VTK_VAR

         SELECT CASE (VTK_VAR(L))
           
            CASE (1)
               CALL WRITE_SCALAR_IN_VTU_ASCII('EP_G',EP_G)
               IF (FULL_LOG.AND.myPE == PE_IO) WRITE(*,10)'.'
             
            CASE (2)
               CALL WRITE_SCALAR_IN_VTU_ASCII('P_G',P_G)
               CALL WRITE_SCALAR_IN_VTU_ASCII('P_S',P_S)
               IF (FULL_LOG.AND.myPE == PE_IO) WRITE(*,10)'.'

            CASE (3)
               CALL WRITE_VECTOR_IN_VTU_ASCII('Gas_Velocity',U_G,V_G,W_G)
               IF (FULL_LOG.AND.myPE == PE_IO) WRITE(*,10)'.'
            
            CASE (4)
               DO M = 1,MMAX
                  WRITE(SUBM,*)M
                  CALL WRITE_VECTOR_IN_VTU_ASCII('Solids_Velocity_'//ADJUSTL(SUBM),U_S(:,M),V_S(:,M),W_S(:,M))
               END DO
               IF (FULL_LOG.AND.myPE == PE_IO) WRITE(*,10)'.'
            
            CASE (5)
               DO M = 1,MMAX
                  WRITE(SUBM,*)M
                  CALL WRITE_SCALAR_IN_VTU_ASCII('Solids_density_'//ADJUSTL(SUBM),ROP_S(:,M))
               END DO
               IF (FULL_LOG.AND.myPE == PE_IO) WRITE(*,10)'.'
            
            CASE (6)
               CALL WRITE_SCALAR_IN_VTU_ASCII('Gas_temperature',T_g)
               DO M = 1,MMAX
                  WRITE(SUBM,*)M
                  CALL WRITE_SCALAR_IN_VTU_ASCII('Solids_temperature_'//ADJUSTL(SUBM),T_S(:,M))
               END DO
               IF (FULL_LOG.AND.myPE == PE_IO) WRITE(*,10)'.'
            
            CASE (7)
               SPECIES_COUNTER = 0
               DO N = 1,NMAX(0)
                  WRITE(SUBN,*)N
                  SPECIES_COUNTER = SPECIES_COUNTER + 1
                  VAR_NAME = ADJUSTL(SPECIES_NAME(SPECIES_COUNTER))
                  LT = LEN_TRIM(ADJUSTL(SPECIES_NAME(SPECIES_COUNTER)))
                  VAR_NAME = VAR_NAME(1:LT)//'_Gas_mass_fractions_'//ADJUSTL(SUBN)
                  CALL WRITE_SCALAR_IN_VTU_ASCII(VAR_NAME,X_g(:,N))
               END DO

               DO M = 1, MMAX 
                  WRITE(SUBM,*)M
                  DO N = 1,NMAX(M)
                     WRITE(SUBN,*)N
                     SPECIES_COUNTER = SPECIES_COUNTER + 1
                     VAR_NAME = ADJUSTL(SPECIES_NAME(SPECIES_COUNTER))
                     LT = LEN_TRIM(ADJUSTL(SPECIES_NAME(SPECIES_COUNTER)))
                     VAR_NAME = VAR_NAME(1:LT)//'_Solids_mass_fractions_'//TRIM(ADJUSTL(SUBM))//'_'//ADJUSTL(SUBN)
                  CALL WRITE_SCALAR_IN_VTU_ASCII(VAR_NAME,X_g(:,N))
                  END DO
               END DO  
               IF (FULL_LOG.AND.myPE == PE_IO) WRITE(*,10)'.'
            
            CASE (8)
               DO M = 1,MMAX
                  WRITE(SUBM,*)M
                  CALL WRITE_SCALAR_IN_VTU_ASCII('Granular_temperature_'//ADJUSTL(SUBM),Theta_m(:,M))
               END DO
               IF (FULL_LOG.AND.myPE == PE_IO) WRITE(*,10)'.'
            
            CASE (9)
               SPECIES_COUNTER = 0
               DO N = 1,NSCALAR
                  WRITE(SUBN,*)N
                  SPECIES_COUNTER = SPECIES_COUNTER + 1
                  VAR_NAME = 'Scalar_'//ADJUSTL(SUBN)
                  CALL WRITE_SCALAR_IN_VTU_ASCII(VAR_NAME,Scalar(:,N))
               END DO
               IF (FULL_LOG.AND.myPE == PE_IO) WRITE(*,10)'.'

            CASE (11)
               IF(K_EPSILON) THEN
                  CALL WRITE_SCALAR_IN_VTU_ASCII('K_Turb_G',K_Turb_G)                
                  CALL WRITE_SCALAR_IN_VTU_ASCII('E_Turb_G',E_Turb_G)                                
                  IF (FULL_LOG.AND.myPE == PE_IO) WRITE(*,10)'.'
               ENDIF

            CASE (12)
               CALL CALC_VORTICITY

               CALL WRITE_SCALAR_IN_VTU_ASCII('VORTICITY_MAG',VORTICITY)                
               CALL WRITE_SCALAR_IN_VTU_ASCII('LAMBDA_2',LAMBDA2)                
             
               IF (FULL_LOG.AND.myPE == PE_IO) WRITE(*,10)'.'
            
            CASE (100)
               IF(DISCRETE_ELEMENT.AND.MPPIC) THEN
                  ALLOCATE( COUNT_DES_BC(DIMENSION_3))
                  DO IJK = IJKSTART3, IJKEND3
                     COUNT_DES_BC (IJK) = 0.d0
                     COUNT_DES_BC (IJK) = DES_CELLWISE_BCDATA(IJK)%COUNT_DES_BC
                  ENDDO
                  CALL WRITE_SCALAR_IN_VTU_ASCII('COUNT_BC',COUNT_DES_BC)
                  DEALLOCATE(COUNT_DES_BC)
               ELSE
                  CALL WRITE_SCALAR_IN_VTU_ASCII('PARTITION',PARTITION)
               ENDIF
               IF (FULL_LOG.AND.myPE == PE_IO) WRITE(*,10)'.'

            CASE (101)
               Allocate(DP_BC_ID(DIMENSION_3))
               DP_BC_ID = DFLOAT(BC_ID)
               CALL WRITE_SCALAR_IN_VTU_ASCII('BC_ID',DP_BC_ID)
               IF (FULL_LOG.AND.myPE == PE_IO) WRITE(*,10)'.'
               DeAllocate(DP_BC_ID)
 
            CASE(999)
               Allocate(IJK_ARRAY(DIMENSION_3))
               DO IJK = IJKSTART3, IJKEND3
                  IJK_ARRAY(IJK) = DFLOAT(IJK)
               ENDDO
               CALL WRITE_SCALAR_IN_VTU_ASCII('IJK',IJK_ARRAY)
               CALL WRITE_SCALAR_IN_VTU_ASCII('F_AT',F_AT)
               DeAllocate(IJK_ARRAY)

            CASE(1000)
               CALL WRITE_VECTOR_IN_VTU_ASCII('Scalar normal',NORMAL_S(:,1),NORMAL_S(:,2),NORMAL_S(:,3))

            CASE (1001)
               CALL WRITE_SCALAR_IN_VTU_ASCII('DEBUG_1',DEBUG_CG(:,1))

            CASE (1002)
               CALL WRITE_SCALAR_IN_VTU_ASCII('DEBUG_2',DEBUG_CG(:,2))

            CASE (1003)
               CALL WRITE_SCALAR_IN_VTU_ASCII('DEBUG_3',DEBUG_CG(:,3))

            CASE (1004)
               CALL WRITE_SCALAR_IN_VTU_ASCII('DEBUG_4',DEBUG_CG(:,4))

            CASE (1005)
               CALL WRITE_SCALAR_IN_VTU_ASCII('DEBUG_5',DEBUG_CG(:,5))

            CASE (1006)
               CALL WRITE_SCALAR_IN_VTU_ASCII('DEBUG_6',DEBUG_CG(:,6))

            CASE (1007)
               CALL WRITE_SCALAR_IN_VTU_ASCII('DEBUG_7',DEBUG_CG(:,7))

            CASE (1008)
               CALL WRITE_SCALAR_IN_VTU_ASCII('DEBUG_8',DEBUG_CG(:,8))

            CASE (1009)
               CALL WRITE_SCALAR_IN_VTU_ASCII('DEBUG_9',DEBUG_CG(:,9))

            CASE (1010)
               CALL WRITE_SCALAR_IN_VTU_ASCII('DEBUG_10',DEBUG_CG(:,10))

            CASE (1011)
               CALL WRITE_SCALAR_IN_VTU_ASCII('DEBUG_11',DEBUG_CG(:,11))

            CASE (1012)
               CALL WRITE_SCALAR_IN_VTU_ASCII('DEBUG_12',DEBUG_CG(:,12))

            CASE (1013)
               CALL WRITE_SCALAR_IN_VTU_ASCII('DEBUG_13',DEBUG_CG(:,13))

            CASE (1014)
               CALL WRITE_SCALAR_IN_VTU_ASCII('DEBUG_14',DEBUG_CG(:,14))

            CASE (1015)
               CALL WRITE_SCALAR_IN_VTU_ASCII('DEBUG_15',DEBUG_CG(:,15))



            CASE (0) ! do nothing

            CASE (UNDEFINED_I) ! do nothing

            CASE DEFAULT

               WRITE(*,30) ' Unknown VTK variable flag ',L,':',VTK_VAR(L)
               WRITE(*,30) ' Available flags are : '
               WRITE(*,30) ' 1 : Void fraction (EP_g)'
               WRITE(*,30) ' 2 : Gas pressure, solids pressure (P_g, P_star)'
               WRITE(*,30) ' 3 : Gas velocity (U_g, V_g, W_g)'
               WRITE(*,30) ' 4 : Solids velocity (U_s, V_s, W_s)'
               WRITE(*,30) ' 5 : Solids density (ROP_s)'
               WRITE(*,30) ' 6 : Gas and solids temperature (T_g, T_s1, T_s2)'
               WRITE(*,30) ' 7 : Gas and solids mass fractions (X_g, X-s)'
               WRITE(*,30) ' 8 : Granular temperature (G)'
               write(*,30) ' 9 : User defined scalars'
!               write(*,30) '10 : Reaction Rates'
               write(*,30) '11 : Turbulence quantities (k and ε)'
               write(*,30) '12 : Gas Vorticity magnitude and Lambda_2 (VORTICITY, LAMBDA_2)'
               write(*,30) '100: Processor assigned to scalar cell (Partition)'
               write(*,30) '101: Boundary condition flag for scalar cell (BC_ID)'
               write(*,30) 'MFiX will exit now.'
               CALL MFIX_EXIT(myPE) 

            END SELECT

      END DO

!      print*,'calling CLOSE_VTU_FILE from myPE=',myPE
      CALL CLOSE_VTU_FILE_ASCII
      CALL UPDATE_AND_CLOSE_PVD_FILE

      call MPI_barrier(MPI_COMM_WORLD,mpierr)

      IF (myPE == PE_IO.AND.TIME_DEPENDENT_FILENAME) THEN
         OPEN(UNIT = VTU_FRAME_UNIT, FILE = TRIM(VTU_FRAME_FILENAME))
         WRITE(VTU_FRAME_UNIT,*)FRAME
         CLOSE(VTU_FRAME_UNIT)
      ENDIF


     IF (FULL_LOG.AND.myPE == PE_IO) WRITE(*,20)' DONE.'

10    FORMAT(A,$)
20    FORMAT(A,1X/)
30    FORMAT(1X,A)    
      RETURN
      
      END SUBROUTINE WRITE_VTU_FILE_ASCII



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: OPEN_VTU_FILE                                          C
!  Purpose: Open a vtu file and writes the header                      C
!           ASCII format                                               C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE OPEN_VTU_FILE_ASCII
    
      USE param
      USE param1
      USE parallel
      USE constant
      USE run
      USE toleranc
      USE geometry
      USE indices  
      USE compar
      USE sendrecv
      USE output
      USE quadric
      USE cutcell
      USE fldvar
      USE vtk
      use cdist

     
      IMPLICIT NONE
      DOUBLE PRECISION:: Xw,Xe,Yn,Ys
      INTEGER :: I,J,K,L,IM,JM,KM,IP,JP,KP,IJK
      INTEGER :: IMJK,IJMK,IJKM,IMJMK,IMJKM,IJMKM,IMJMKM
      LOGICAL :: VTU_FRAME_FILE_EXISTS
      INTEGER :: ISTAT 

      include "function.inc"


!      print*,'Entering Open_vtu_file from myPE=',MyPE, (myPE /= PE_IO.AND.(.NOT.BDIST_IO))

      IF (myPE /= PE_IO.AND.(.NOT.BDIST_IO)) RETURN 


      IF(TIME_DEPENDENT_FILENAME) THEN
         INQUIRE(FILE=VTU_FRAME_FILENAME,EXIST=VTU_FRAME_FILE_EXISTS)
         IF(VTU_FRAME_FILE_EXISTS) THEN
            OPEN(UNIT = VTU_FRAME_UNIT, FILE = TRIM(VTU_FRAME_FILENAME))
            READ(VTU_FRAME_UNIT,*)FRAME
            CLOSE(VTU_FRAME_UNIT)
         ENDIF
         IF(RESET_FRAME_AT_TIME_ZERO.AND.TIME==ZERO) FRAME=-1
         FRAME = FRAME + 1 
      ENDIF



! Define file name
      IF (BDIST_IO) THEN 
         IF(TIME_DEPENDENT_FILENAME) THEN
            WRITE(VTU_FILENAME,20) TRIM(RUN_NAME),FRAME,MYPE
         ELSE
            WRITE(VTU_FILENAME,25) TRIM(RUN_NAME),MYPE
         ENDIF
      ELSE 
         IF(MYPE.EQ.PE_IO) THEN 
            IF(TIME_DEPENDENT_FILENAME) THEN
               WRITE(VTU_FILENAME,30) TRIM(RUN_NAME),FRAME
            ELSE
               WRITE(VTU_FILENAME,35) TRIM(RUN_NAME)
            ENDIF
         END IF  
      END IF 


      IF(TRIM(VTU_DIR)/='.') VTU_FILENAME='./'//TRIM(VTU_DIR)//'/'//VTU_FILENAME

!      print*,'From Open_vtu_file:',MyPE,BDIST_IO, TRIM(VTU_FILENAME)

      IF (FULL_LOG) THEN
         WRITE(*,10)' WRITING VTU FILE : ', TRIM(VTU_FILENAME),' .'
      ENDIF


      OPEN(UNIT = VTU_UNIT, FILE = TRIM(VTU_FILENAME),IOSTAT=ISTAT)

      IF (ISTAT /= 0) THEN
         IF(DMP_LOG) WRITE(UNIT_LOG, 1001) VTU_FILENAME, VTU_UNIT,VTU_DIR
         IF(FULL_LOG.AND.myPE == PE_IO) WRITE(*, 1001) VTU_FILENAME, VTU_UNIT,VTU_DIR
         CALL MFIX_EXIT(myPE)
      ENDIF


 1001 FORMAT(/1X,70('*')//, ' From: OPEN_VTU_FILE',/,' Message: ',          &
         'Error opening vtu file. Terminating run.',/10X,          &
         'File name:  ',A,/10X,                                         &
         'DES_UNIT :  ',i4, /10X,                                       &
         'PLEASE VERIFY THAT VTU_DIR EXISTS: ', A, &
         /1X,70('*')/)


      WRITE(VTU_UNIT,100) '<?xml version="1.0"?>'
      WRITE(VTU_UNIT,110) '<!-- Time =',TIME,' sec. -->'
      WRITE(VTU_UNIT,120) '<VTKFile type="UnstructuredGrid"',&
               ' version="0.1" byte_order="LittleEndian">'

      WRITE(VTU_UNIT,100) '  <UnstructuredGrid>'




      IF (myPE == PE_IO.AND.BDIST_IO) THEN    ! Open .pvtu file that combines all *.vtu files for a given FRAME

         IF(TIME_DEPENDENT_FILENAME) THEN
            WRITE(PVTU_FILENAME,40) TRIM(RUN_NAME),FRAME
         ELSE
            WRITE(PVTU_FILENAME,45) TRIM(RUN_NAME)
         ENDIF

         IF(TRIM(VTU_DIR)/='.') PVTU_FILENAME='./'//TRIM(VTU_DIR)//'/'//PVTU_FILENAME

         OPEN(UNIT = PVTU_UNIT, FILE = TRIM(PVTU_FILENAME))

         WRITE(PVTU_UNIT,100) '<?xml version="1.0"?>'
         WRITE(PVTU_UNIT,110) '<!-- Time =',TIME,' sec. -->'
         WRITE(PVTU_UNIT,120) '<VTKFile type="PUnstructuredGrid"',&
                  ' version="0.1" byte_order="LittleEndian">'

         WRITE(PVTU_UNIT,100) '  <PUnstructuredGrid GhostLevel="0">'
         WRITE(PVTU_UNIT,100) '      <PPoints>'
         WRITE(PVTU_UNIT,100) '        <PDataArray type="Float32" NumberOfComponents="3" format="ascii"/>'
         WRITE(PVTU_UNIT,100) '      </PPoints>'
         WRITE(PVTU_UNIT,100) ''
         WRITE(PVTU_UNIT,100) '      <PCellData Scalars="scalars">'


      ENDIF




100   FORMAT(A)
110   FORMAT(A,E14.8,A)
120   FORMAT(A,A)
10    FORMAT(/1X,3A,$)
20    FORMAT(A,"_",I4.4,"_",I5.5,".vtu")
25    FORMAT(A,"_",I5.5,".vtu")
30    FORMAT(A,"_",I4.4,".vtu")
35    FORMAT(A,".vtu")
40    FORMAT(A,"_",I4.4,".pvtu")
45    FORMAT(A,".pvtu")

      RETURN

      END SUBROUTINE OPEN_VTU_FILE_ASCII

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: OPEN_PVD_FILE                                          C
!  Purpose: Open a PVD file and writes the header                      C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE OPEN_PVD_FILE
    
      USE param
      USE param1
      USE parallel
      USE constant
      USE run
      USE toleranc
      USE geometry
      USE indices  
      USE compar
      USE sendrecv
      USE output
      USE quadric
      USE cutcell
      USE fldvar
      USE vtk

     
      IMPLICIT NONE
      DOUBLE PRECISION:: Xw,Xe,Yn,Ys
      INTEGER :: I,J,K,L,IM,JM,KM,IP,JP,KP,IJK
      INTEGER :: IMJK,IJMK,IJKM,IMJMK,IMJKM,IJMKM,IMJMKM
      LOGICAL :: PVD_EXISTS,VTU_FRAME_FILE_EXISTS

      include "function.inc"

      IF (myPE /= PE_IO) RETURN 


      PVD_FILENAME = TRIM(RUN_NAME) // '.pvd'


! First, check if the file already exists.

      INQUIRE(FILE=PVD_FILENAME,EXIST=PVD_EXISTS)

! The first time this subroutine is executed, properly initialize the pvd file

      IF(.NOT.PVD_FILE_INITIALIZED) THEN

         IF(RUN_TYPE == 'NEW')THEN
            ! For a new run, the pvd file should not exist, and is created with appropriate header
            IF (.NOT.PVD_EXISTS) THEN
               OPEN(UNIT = PVD_UNIT, FILE = TRIM(PVD_FILENAME))
               WRITE(PVD_UNIT,100) '<?xml version="1.0"?>'
               WRITE(PVD_UNIT,100) '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">'
               WRITE(PVD_UNIT,100) '<Collection>'
!               CALL UPDATE_AND_CLOSE_PVD_FILE
               PVD_FILE_INITIALIZED=.TRUE.
            ELSE ! If the pvd file exists, print error message and exits
               WRITE(*,1002) TRIM(PVD_FILENAME)
               WRITE(UNIT_LOG, 1002) TRIM(PVD_FILENAME)
               CALL MFIX_EXIT(myPE)
            ENDIF
         ELSE
            ! For a restart run, the pvd file must exist
            IF (.NOT.PVD_EXISTS) THEN
               ! If the pvd file does not exist, print error message and exits
               WRITE(*,1003) TRIM(PVD_FILENAME)
               WRITE(UNIT_LOG, 1003) TRIM(PVD_FILENAME)
               CALL MFIX_EXIT(myPE)
            ELSE 
           ! If it already exists, go to the bottom of the file and prepare to append data (remove last two lines)
               OPEN(UNIT=PVD_UNIT,FILE = TRIM(PVD_FILENAME),POSITION="APPEND",STATUS='OLD')
               BACKSPACE(PVD_UNIT)
               BACKSPACE(PVD_UNIT)
               PVD_FILE_INITIALIZED=.TRUE.
            ENDIF
         ENDIF
      ELSE ! When properly initialized, open the file and go to the bottom of the file and prepare to append data (remove last two lines)
         OPEN(UNIT=PVD_UNIT,FILE = TRIM(PVD_FILENAME),POSITION="APPEND",STATUS='OLD')
         BACKSPACE(PVD_UNIT)
         BACKSPACE(PVD_UNIT)
      ENDIF


100   FORMAT(A)

1002  FORMAT(/1X,70('*')/,' From: OPEN_PVD_FILE',/,' Message: ',       &
         A,' already exists in the run directory.',/10X,               &
           'This is not allowed for a new run.',/10X,                   &
         'Terminating run.',/1X,70('*')/)

1003  FORMAT(/1X,70('*')/,' From: OPEN_PVD_FILE',/,' Message: ',       &
         A,' is missing from the  the run directory,',/10X,            &
           ' and must be present for a restart run.',/10X,              &
         'Terminating run.',/1X,70('*')/)

1004  FORMAT(/1X,70('*')/,' From: OPEN_PVD_FILE',/,' Message: ',       &
         ' Current VTU frame is ',I10,/10X,                              &
         ' (from ',A,').',/1X,70('*')/)

1005  FORMAT(/1X,70('*')/,' From: OPEN_PVD_FILE',/,' Message: ',       &
         ' Current VTU frame is ',I10,/10X,                              &
           ' (from mfix.dat).',/1X,70('*')/)



      RETURN

      END SUBROUTINE OPEN_PVD_FILE



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_GEOMETRY_IN_VTU                                  C
!  Purpose: Write Geometry and connectivity in a vtu file              C
!           ASCII format                                               C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE WRITE_GEOMETRY_IN_VTU_ASCII
    
      USE param
      USE param1
      USE parallel
      USE constant
      USE run
      USE toleranc
      USE geometry
      USE indices  
      USE compar
      USE mpi_utility 
      USE sendrecv
      USE quadric
      USE cutcell
      USE fldvar
      USE vtk
      USE cdist

     
      IMPLICIT NONE

      INTEGER :: IJK,I,J,K,L
      INTEGER :: IJK_OFFSET,OFFSET

      INTEGER :: iproc,IERR
      INTEGER, DIMENSION(0:numPEs-1) :: disp,rcount
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: SHIFTED_CONNECTIVITY
      INTEGER :: CELL_TYPE

      include "function.inc"

      IF (myPE == PE_IO.AND.(.NOT.BDIST_IO)) THEN

         NUMBER_OF_VTK_CELLS = NUMBER_OF_CELLS - NUMBER_OF_BLOCKED_CELLS

         WRITE(VTU_UNIT,100) '    <Piece NumberOfPoints="',NUMBER_OF_POINTS,'" NumberOfCells="',NUMBER_OF_VTK_CELLS,'" >'
         WRITE(VTU_UNIT,110) '      <Points>'
         WRITE(VTU_UNIT,110) '        <DataArray type="Float32" NumberOfComponents="3" format="ascii">'

         DO IJK = 1,IJKMAX3
            WRITE(VTU_UNIT,120)XG_E(GLOBAL_I_OF(IJK)),YG_N(GLOBAL_J_OF(IJK)),ZG_T(GLOBAL_K_OF(IJK))
         ENDDO
         DO IJK = 1,GLOBAL_NUMBER_OF_NEW_POINTS
            WRITE(VTU_UNIT,120)GLOBAL_X_NEW_POINT(IJK),GLOBAL_Y_NEW_POINT(IJK),GLOBAL_Z_NEW_POINT(IJK)
         ENDDO

         WRITE(VTU_UNIT,110) '        </DataArray>'
         WRITE(VTU_UNIT,110) ''
         WRITE(VTU_UNIT,110) '      </Points>'
         WRITE(VTU_UNIT,110) ''
         WRITE(VTU_UNIT,110) '      <Cells>'
         WRITE(VTU_UNIT,110) '        <DataArray type="Int32" Name="connectivity" format="ascii">'

         DO IJK = 1,IJKMAX3
            IF (GLOBAL_INTERIOR_CELL_AT(IJK))      THEN
               IF (.NOT.GLOBAL_BLOCKED_CELL_AT(IJK)) WRITE(VTU_UNIT,130) (GLOBAL_CONNECTIVITY(IJK,L)-1,L=1,GLOBAL_NUMBER_OF_NODES(IJK))
            ENDIF
         END DO

         WRITE(VTU_UNIT,110) '        </DataArray>'
         WRITE(VTU_UNIT,110) ''
         WRITE(VTU_UNIT,110) '        <DataArray type="Int32" Name="offsets" format="ascii">'

         OFFSET = 0
         DO IJK = 1,IJKMAX3
            IF (GLOBAL_INTERIOR_CELL_AT(IJK))      THEN
               IF (.NOT.GLOBAL_BLOCKED_CELL_AT(IJK)) THEN
                  OFFSET = OFFSET + GLOBAL_NUMBER_OF_NODES(IJK)
                  WRITE(VTU_UNIT,140) OFFSET
               ENDIF
            ENDIF
         END DO

         WRITE(VTU_UNIT,110) '        </DataArray>'
         WRITE(VTU_UNIT,110) ''
         WRITE(VTU_UNIT,110) '        <DataArray type="UInt8" Name="types" format="ascii">'

         IF(NO_K) THEN
            CELL_TYPE = 7
         ELSE
            CELL_TYPE = 41
         ENDIF

         DO IJK = 1,IJKMAX3
            IF (GLOBAL_INTERIOR_CELL_AT(IJK))  WRITE(VTU_UNIT,140) CELL_TYPE
         END DO
     

      ELSEIF(BDIST_IO) THEN

         NUMBER_OF_VTK_CELLS = NUMBER_OF_CELLS - NUMBER_OF_BLOCKED_CELLS

         WRITE(VTU_UNIT,100) '    <Piece NumberOfPoints="',NUMBER_OF_POINTS,'" NumberOfCells="',NUMBER_OF_VTK_CELLS,'" >'
         WRITE(VTU_UNIT,110) '      <Points>'
         WRITE(VTU_UNIT,110) '        <DataArray type="Float32" NumberOfComponents="3" format="ascii">'

         DO IJK = 1,IJKEND3
            WRITE(VTU_UNIT,120)XG_E(I_OF(IJK)),YG_N(J_OF(IJK)),ZG_T(K_OF(IJK))
         ENDDO
         DO IJK = 1,NUMBER_OF_NEW_POINTS
            WRITE(VTU_UNIT,120)X_NEW_POINT(IJK),Y_NEW_POINT(IJK),Z_NEW_POINT(IJK)
         ENDDO

         WRITE(VTU_UNIT,110) '        </DataArray>'
         WRITE(VTU_UNIT,110) ''
         WRITE(VTU_UNIT,110) '      </Points>'
         WRITE(VTU_UNIT,110) ''
         WRITE(VTU_UNIT,110) '      <Cells>'
         WRITE(VTU_UNIT,110) '        <DataArray type="Int32" Name="connectivity" format="ascii">'

         DO IJK = 1,IJKEND3
            IF (INTERIOR_CELL_AT(IJK))      THEN
               IF (.NOT.BLOCKED_CELL_AT(IJK)) WRITE(VTU_UNIT,130) (CONNECTIVITY(IJK,L)-1,L=1,NUMBER_OF_NODES(IJK))
            ENDIF
         END DO

         WRITE(VTU_UNIT,110) '        </DataArray>'
         WRITE(VTU_UNIT,110) ''
         WRITE(VTU_UNIT,110) '        <DataArray type="Int32" Name="offsets" format="ascii">'

         OFFSET = 0
         DO IJK = 1,IJKEND3
            IF (INTERIOR_CELL_AT(IJK))      THEN
               IF (.NOT.BLOCKED_CELL_AT(IJK)) THEN
                  OFFSET = OFFSET + NUMBER_OF_NODES(IJK)
                  WRITE(VTU_UNIT,140) OFFSET
               ENDIF
            ENDIF
         END DO

         WRITE(VTU_UNIT,110) '        </DataArray>'
         WRITE(VTU_UNIT,110) ''
         WRITE(VTU_UNIT,110) '        <DataArray type="UInt8" Name="types" format="ascii">'

         IF(NO_K) THEN
            CELL_TYPE = 7
         ELSE
            CELL_TYPE = 41
         ENDIF

         DO IJK = 1,IJKMAX3
            IF (INTERIOR_CELL_AT(IJK))  WRITE(VTU_UNIT,140) CELL_TYPE
         END DO
     


      ENDIF



      WRITE(VTU_UNIT,110) '        </DataArray>'
      WRITE(VTU_UNIT,110) ''
      WRITE(VTU_UNIT,110) '      </Cells>'
      WRITE(VTU_UNIT,100) '      <CellData Scalars="scalars">'



100   FORMAT(A,I8,A,I8,A)
110   FORMAT(A)
120   FORMAT(10X,3(E16.8E3,2X))
130   FORMAT(10X,15(I8,2X))
140   FORMAT(10X,I8)





      RETURN
      
      END SUBROUTINE WRITE_GEOMETRY_IN_VTU_ASCII




!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_SCALAR_IN_VTU                                    C
!  Purpose: Write Scalar variable in a vtu file                        C
!           ASCII format                                               C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE WRITE_SCALAR_IN_VTU_ASCII(VAR_NAME,VAR)
    
      USE param
      USE param1
      USE parallel
      USE constant
      USE run
      USE toleranc
      USE geometry
      USE indices  
      USE compar
      USE mpi_utility 
      USE sendrecv
      USE quadric
      USE cutcell
      USE fldvar
      USE vtk
      USE cdist

      IMPLICIT NONE
      INTEGER :: I,IJK,L

      CHARACTER (*) :: VAR_NAME
      DOUBLE PRECISION, DIMENSION(DIMENSION_3) ::  VAR
      DOUBLE PRECISION, ALLOCATABLE :: GLOBAL_VAR(:)

      include "function.inc"

      IF (.NOT.BDIST_IO) THEN

         IF (myPE == PE_IO) THEN
            allocate (GLOBAL_VAR(ijkmax3))     
         ELSE
            allocate (GLOBAL_VAR(1))     
         ENDIF

         call gather (VAR,GLOBAL_VAR,root) 


         IF (myPE /= PE_IO) RETURN


         DO I = 1,LEN_TRIM(VAR_NAME)
            IF(VAR_NAME(I:I) == ' ') VAR_NAME(I:I) = '_'
         ENDDO


         WRITE(VTU_UNIT,110) '        <DataArray type="Float32" Name="',TRIM(VAR_NAME),'" format="ascii">'
         DO IJK = 1,IJKMAX3
            IF (GLOBAL_INTERIOR_CELL_AT(IJK))      THEN
               IF (.NOT.GLOBAL_BLOCKED_CELL_AT(IJK))   WRITE(VTU_UNIT,120) GLOBAL_VAR(IJK)
            ENDIF
         ENDDO
         WRITE(VTU_UNIT,100) '        </DataArray>'

         Deallocate (GLOBAL_VAR)   


      ELSE ! BDIST_IO=.TRUE.


         DO I = 1,LEN_TRIM(VAR_NAME)
            IF(VAR_NAME(I:I) == ' ') VAR_NAME(I:I) = '_'
         ENDDO

         WRITE(VTU_UNIT,110) '        <DataArray type="Float32" Name="',TRIM(VAR_NAME),'" format="ascii">'
         DO IJK = 1,IJKEND3
            IF (INTERIOR_CELL_AT(IJK))      THEN
               IF (.NOT.BLOCKED_CELL_AT(IJK))   WRITE(VTU_UNIT,120) VAR(IJK)
            ENDIF
         ENDDO
         WRITE(VTU_UNIT,100) '        </DataArray>'

         IF (myPE == PE_IO) THEN       ! Update pvtu file with variable name
            WRITE(PVTU_UNIT,110) '        <PDataArray type="Float32" Name="',TRIM(VAR_NAME),'" format="ascii"/>'
         ENDIF


      ENDIF


100   FORMAT(A)
110   FORMAT(A,A,A)
120   FORMAT(10X,E16.8E3)



      RETURN
      
      END SUBROUTINE WRITE_SCALAR_IN_VTU_ASCII



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_VECTOR_IN_VTU                                    C
!  Purpose: Write Vector variable in a vtu file                        C
!           ASCII format                                               C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE WRITE_VECTOR_IN_VTU_ASCII(VAR_NAME,VARX,VARY,VARZ)
    
      USE param
      USE param1
      USE parallel
      USE constant
      USE run
      USE toleranc
      USE geometry
      USE indices  
      USE compar
      USE mpi_utility 
      USE sendrecv
      USE quadric
      USE cutcell
      USE fldvar
      USE vtk
      USE cdist
      
      IMPLICIT NONE
      INTEGER :: IJK,L

      CHARACTER (*) :: VAR_NAME
      DOUBLE PRECISION, DIMENSION(DIMENSION_3) ::  VARX,VARY,VARZ
      DOUBLE PRECISION, ALLOCATABLE :: GLOBAL_VARX(:),GLOBAL_VARY(:),GLOBAL_VARZ(:)

      include "function.inc"

      IF (.NOT.BDIST_IO) THEN

         IF (myPE == PE_IO) THEN
            allocate (GLOBAL_VARX(ijkmax3))
            allocate (GLOBAL_VARY(ijkmax3))     
            allocate (GLOBAL_VARZ(ijkmax3))          
         ELSE
            allocate (GLOBAL_VARX(1))
            allocate (GLOBAL_VARY(1))     
            allocate (GLOBAL_VARZ(1))          
         ENDIF

         call gather (VARX,GLOBAL_VARX,root)
         call gather (VARY,GLOBAL_VARY,root) 
         call gather (VARZ,GLOBAL_VARZ,root)  

         IF (myPE /= PE_IO) RETURN


         WRITE(VTU_UNIT,110) '        <DataArray type="Float32" Name="',TRIM(VAR_NAME),'" NumberOfComponents="3"  format="ascii">'
         DO IJK = 1,IJKMAX3
            IF (GLOBAL_INTERIOR_CELL_AT(IJK))      THEN
               IF (.NOT.GLOBAL_BLOCKED_CELL_AT(IJK))   WRITE(VTU_UNIT,120) GLOBAL_VARX(IJK),GLOBAL_VARY(IJK),GLOBAL_VARZ(IJK)
            ENDIF
         ENDDO
         WRITE(VTU_UNIT,100) '        </DataArray>'

         Deallocate (GLOBAL_VARX)
         Deallocate (GLOBAL_VARY)   
         Deallocate (GLOBAL_VARZ) 


      ELSE ! BDIST_IO=.TRUE.


         WRITE(VTU_UNIT,110) '        <DataArray type="Float32" Name="',TRIM(VAR_NAME),'" NumberOfComponents="3"  format="ascii">'
         DO IJK = 1,IJKEND3
            IF (INTERIOR_CELL_AT(IJK))      THEN
               IF (.NOT.BLOCKED_CELL_AT(IJK))   WRITE(VTU_UNIT,120) VARX(IJK),VARY(IJK),VARZ(IJK)
            ENDIF
         ENDDO
         WRITE(VTU_UNIT,100) '        </DataArray>'

         IF (myPE == PE_IO) THEN       ! Update pvtu file with variable name
            WRITE(PVTU_UNIT,110) '        <PDataArray type="Float32" Name="',TRIM(VAR_NAME),'" NumberOfComponents="3"  format="ascii"/>'
         ENDIF

      ENDIF


100   FORMAT(A)
110   FORMAT(A,A,A)
120   FORMAT(10X,3(E16.8E3,2X))


      RETURN
      
      END SUBROUTINE WRITE_VECTOR_IN_VTU_ASCII



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CLOSE_VTU_FILE                                         C
!  Purpose: Close a vtu file                                           C
!           ASCII format                                               C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE CLOSE_VTU_FILE_ASCII

      USE compar
      Use run
      USE vtk
      use cdist

      IMPLICIT NONE
    
      INTEGER:: N
      CHARACTER (LEN=32)  :: VTU_NAME

      IF (myPE /= PE_IO.AND.(.NOT.BDIST_IO)) RETURN 

!      print*,'VTU_UNIT,MyPE=',VTU_UNIT,MyPE


      WRITE(VTU_UNIT,100) '      </CellData>'
      WRITE(VTU_UNIT,100) '    </Piece>'
      WRITE(VTU_UNIT,100) '  </UnstructuredGrid>'
      WRITE(VTU_UNIT,100) '</VTKFile>'
      CLOSE(VTU_UNIT)

      IF (myPE == PE_IO.AND.BDIST_IO) THEN       ! Update pvtu file and close
         WRITE(PVTU_UNIT,100) '      </PCellData>'


         DO N = 0,NumPEs-1
            IF(TIME_DEPENDENT_FILENAME) THEN
               WRITE(VTU_NAME,20) TRIM(RUN_NAME),FRAME,N
            ELSE
               WRITE(VTU_NAME,25) TRIM(RUN_NAME),N
            ENDIF

            WRITE(PVTU_UNIT,110) '      <Piece Source="',TRIM(VTU_NAME),'"/>'
         ENDDO


         WRITE(PVTU_UNIT,100) '  </PUnstructuredGrid>'
         WRITE(PVTU_UNIT,100) '</VTKFile>'
      ENDIF




20    FORMAT(A,"_",I4.4,"_",I5.5,".vtu")
25    FORMAT(A,"_",I5.5,".vtu")

100   FORMAT(A)
110   FORMAT(A,A,A)



      RETURN

      END SUBROUTINE CLOSE_VTU_FILE_ASCII

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: UPDATE_AND_CLOSE_PVD_FILE                              C
!  Purpose: Updates and close a pvd file                               C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE UPDATE_AND_CLOSE_PVD_FILE
    
      USE compar
      use cdist
      USE run
      USE vtk

      IMPLICIT NONE

      CHARACTER (LEN=32)  :: FILENAME

      IF (myPE /= PE_IO) RETURN 


      IF(.NOT.BDIST_IO) THEN
         FILENAME=VTU_FILENAME
      ELSE
         IF(TIME_DEPENDENT_FILENAME) THEN
            WRITE(FILENAME,40) TRIM(RUN_NAME),FRAME
         ELSE
            WRITE(FILENAME,45) TRIM(RUN_NAME)
         ENDIF
         IF(TRIM(VTU_DIR)/='.') FILENAME='./'//TRIM(VTU_DIR)//'/'//FILENAME
      ENDIF


! Write the data to the file
         WRITE(PVD_UNIT,100)&
         '<DataSet timestep="',TIME,'" ',& ! simulation time
         'group="" part="0" ',& ! necessary file data
         'file="',TRIM(FILENAME),'"/>' ! file name of vtp

! Write the closing tags
         WRITE(PVD_UNIT,110)'</Collection>'
         WRITE(PVD_UNIT,110)'</VTKFile>'
         
         CLOSE(PVD_UNIT)


40    FORMAT(A,"_",I4.4,".pvtu")
45    FORMAT(A,".pvtu")
100   FORMAT(6X,A,E14.8,5A)
110   FORMAT(A)


      RETURN

      END SUBROUTINE UPDATE_AND_CLOSE_PVD_FILE


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_VTK_FILE                                         C
!  Purpose: Writes the cut cell grid in VTK format                     C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE WRITE_VTK_FILE
    
      USE param
      USE param1
      USE parallel
      USE constant
      USE run
      USE toleranc
      USE geometry
      USE indices  
      USE compar
      USE sendrecv
      USE quadric
      USE cutcell
      USE fldvar
      USE visc_s
      USE physprop
      USE pgcor
      USE vtk
      USE rxns      
      USE output
      USE scalars

      USE pgcor
      USE pscor
      USE discretelement, Only : DES_CELLWISE_BCDATA, DISCRETE_ELEMENT
      USE mfix_pic

      IMPLICIT NONE
      DOUBLE PRECISION:: Xw,Xe,Yn,Ys
      INTEGER :: I,J,K,L,M,N,IM,JM,KM,IP,JP,KP,IJK
      INTEGER :: IMJK,IJMK,IJKM,IMJMK,IMJKM,IJMKM,IMJMKM

      INTEGER sw,se,ne,nw
      INTEGER, DIMENSION(10) :: additional_node
      DOUBLE PRECISION, DIMENSION(2*DIMENSION_3) ::  X_OF
      DOUBLE PRECISION, DIMENSION(2*DIMENSION_3) ::  Y_OF
      DOUBLE PRECISION, DIMENSION(2*DIMENSION_3) ::  Z_OF
      INTEGER, DIMENSION(DIMENSION_3) ::  INDEX_OF_E_ADD_NODE
      INTEGER, DIMENSION(DIMENSION_3) ::  INDEX_OF_N_ADD_NODE
      INTEGER :: SPECIES_COUNTER,LT

      CHARACTER (LEN=32) :: SUBM,SUBN
      CHARACTER (LEN=64) :: VAR_NAME

      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  DP_BC_ID, COUNT_DES_BC

      include "function.inc"

      IF(.NOT.CARTESIAN_GRID) RETURN

      DX(IMAX3+1) = DX(IMAX3)
      DY(JMAX3+1) = DY(JMAX3)
      DZ(KMAX3+1) = DZ(KMAX3)

!     Location of U-momentum cells for original (uncut grid)
      IF (DO_I) THEN 
        XG_E(1) = ZERO
        DO I = IMIN1, IMAX2 
           XG_E(I) = XG_E(I-1) + DX(I) 
        END DO 
      ENDIF

!     Location of V-momentum cells for original (uncut grid)
      IF (DO_J) THEN 
        YG_N(1) = ZERO
        DO J = JMIN1, JMAX2 
           YG_N(J) = YG_N(J-1) + DY(J) 
        END DO 
      ENDIF

!     Location of W-momentum cells for original (uncut grid)
      IF (DO_K) THEN 
        ZG_T(1) = ZERO
        DO K = KMIN1, KMAX2 
           ZG_T(K) = ZG_T(K-1) + DZ(K) 
        END DO 
      ELSE
         ZG_T = ZERO
      ENDIF

      IF(WRITE_ANI_CUTCELL) THEN
         CALL OPEN_VTK_FILE
         CALL WRITE_GEOMETRY_IN_VTK
         CALL CLOSE_VTK_FILE
         IF (FULL_LOG) THEN
            WRITE(*,30)'WROTE VTK FILE : ani_cutcell.vtk'
         ENDIF
         WRITE_ANI_CUTCELL = .FALSE.
         RETURN
      ENDIF

      CALL OPEN_VTK_FILE

      CALL WRITE_GEOMETRY_IN_VTK

      DO L = 1, DIM_VTK_VAR

         SELECT CASE (VTK_VAR(L))
           
            CASE (1)
               CALL WRITE_SCALAR_IN_VTK('EP_G',EP_G)
               IF (FULL_LOG.AND.myPE == PE_IO) WRITE(*,10)'.'
             
            CASE (2)
               CALL WRITE_SCALAR_IN_VTK('P_G',P_G)
               CALL WRITE_SCALAR_IN_VTK('P_S',P_S)
               IF (FULL_LOG.AND.myPE == PE_IO) WRITE(*,10)'.'

            CASE (3)
               CALL WRITE_VECTOR_IN_VTK('Gas_Velocity',U_G,V_G,W_G)
               IF (FULL_LOG.AND.myPE == PE_IO) WRITE(*,10)'.'
            
            CASE (4)
               DO M = 1,MMAX
                  WRITE(SUBM,*)M
                  CALL WRITE_VECTOR_IN_VTK('Solids_Velocity_'//ADJUSTL(SUBM),U_S(:,M),V_S(:,M),W_S(:,M))
               END DO
               IF (FULL_LOG.AND.myPE == PE_IO) WRITE(*,10)'.'
            
            CASE (5)
               DO M = 1,MMAX
                  WRITE(SUBM,*)M
                  CALL WRITE_SCALAR_IN_VTK('Solids_density_'//ADJUSTL(SUBM),ROP_S(:,M))
               END DO
               IF (FULL_LOG.AND.myPE == PE_IO) WRITE(*,10)'.'
            
            CASE (6)
               CALL WRITE_SCALAR_IN_VTK('Gas_temperature',T_g)
               DO M = 1,MMAX
                  WRITE(SUBM,*)M
                  CALL WRITE_SCALAR_IN_VTK('Solids_temperature_'//ADJUSTL(SUBM),T_S(:,M))
               END DO
               IF (FULL_LOG.AND.myPE == PE_IO) WRITE(*,10)'.'
            
            CASE (7)
               SPECIES_COUNTER = 0
               DO N = 1,NMAX(0)
                  WRITE(SUBN,*)N
                  SPECIES_COUNTER = SPECIES_COUNTER + 1
                  VAR_NAME = ADJUSTL(SPECIES_NAME(SPECIES_COUNTER))
                  LT = LEN_TRIM(ADJUSTL(SPECIES_NAME(SPECIES_COUNTER)))
                  VAR_NAME = VAR_NAME(1:LT)//'_Gas_mass_fractions_'//ADJUSTL(SUBN)
                  CALL WRITE_SCALAR_IN_VTK(VAR_NAME,X_g(:,N))
               END DO

               DO M = 1, MMAX 
                  WRITE(SUBM,*)M
                  DO N = 1,NMAX(M)
                     WRITE(SUBN,*)N
                     SPECIES_COUNTER = SPECIES_COUNTER + 1
                     VAR_NAME = ADJUSTL(SPECIES_NAME(SPECIES_COUNTER))
                     LT = LEN_TRIM(ADJUSTL(SPECIES_NAME(SPECIES_COUNTER)))
                     VAR_NAME = VAR_NAME(1:LT)//'_Solids_mass_fractions_'//TRIM(ADJUSTL(SUBM))//'_'//ADJUSTL(SUBN)
                     CALL WRITE_SCALAR_IN_VTK(VAR_NAME,X_s(:,M,N))
                  END DO
               END DO  
               IF (FULL_LOG.AND.myPE == PE_IO) WRITE(*,10)'.'
            
            CASE (8)
               DO M = 1,MMAX
                  WRITE(SUBM,*)M
                  CALL WRITE_SCALAR_IN_VTK('Granular_temperature_'//ADJUSTL(SUBM),Theta_m(:,M))
               END DO
               IF (FULL_LOG.AND.myPE == PE_IO) WRITE(*,10)'.'
            
            CASE (9)
               SPECIES_COUNTER = 0
               DO N = 1,NSCALAR
                  WRITE(SUBN,*)N
                  SPECIES_COUNTER = SPECIES_COUNTER + 1
                  VAR_NAME = 'Scalar_'//ADJUSTL(SUBN)
                  CALL WRITE_SCALAR_IN_VTK(VAR_NAME,Scalar(:,N))
               END DO
               IF (FULL_LOG.AND.myPE == PE_IO) WRITE(*,10)'.'

            CASE (11)
               IF(K_EPSILON) THEN
                  CALL WRITE_SCALAR_IN_VTK('K_Turb_G',K_Turb_G)                
                  CALL WRITE_SCALAR_IN_VTK('E_Turb_G',E_Turb_G)                
                  IF (FULL_LOG.AND.myPE == PE_IO) WRITE(*,10)'.'
               ENDIF

            CASE (12)
               CALL CALC_VORTICITY

               CALL WRITE_SCALAR_IN_VTK('VORTICITY_MAG',VORTICITY)                
               CALL WRITE_SCALAR_IN_VTK('LAMBDA_2',LAMBDA2)                
               IF (FULL_LOG.AND.myPE == PE_IO) WRITE(*,10)'.'
            
            CASE (100)
               IF(DISCRETE_ELEMENT.AND.MPPIC) THEN
                  ALLOCATE( COUNT_DES_BC(DIMENSION_3))
                  DO IJK = IJKSTART3, IJKEND3
                     COUNT_DES_BC (IJK) = 0.d0
                     COUNT_DES_BC (IJK) = DES_CELLWISE_BCDATA(IJK)%COUNT_DES_BC
                  ENDDO
                  CALL WRITE_SCALAR_IN_VTK('COUNT_BC',COUNT_DES_BC)
                  DEALLOCATE(COUNT_DES_BC)
               ELSE
                  CALL WRITE_SCALAR_IN_VTK('PARTITION',PARTITION)
               ENDIF
               IF (FULL_LOG.AND.myPE == PE_IO) WRITE(*,10)'.'
            CASE (101)
       
               Allocate(DP_BC_ID(DIMENSION_3))
               DP_BC_ID = DFLOAT(BC_ID)
!               CALL WRITE_SCALAR_IN_VTK('BC_ID',DFLOAT(BC_ID))
               CALL WRITE_SCALAR_IN_VTK('BC_ID',DP_BC_ID)
               IF (FULL_LOG.AND.myPE == PE_IO) WRITE(*,10)'.'
               DeAllocate(DP_BC_ID)

            CASE (0) ! do nothing

            CASE (UNDEFINED_I) ! do nothing

            CASE DEFAULT

               WRITE(*,30) ' Unknown VTK variable flag ',L,':',VTK_VAR(L)
               WRITE(*,30) ' Available flags are : '
               WRITE(*,30) ' 1 : Void fraction (EP_g)'
               WRITE(*,30) ' 2 : Gas pressure, solids pressure (P_g, P_star)'
               WRITE(*,30) ' 3 : Gas velocity (U_g, V_g, W_g)'
               WRITE(*,30) ' 4 : Solids velocity (U_s, V_s, W_s)'
               WRITE(*,30) ' 5 : Solids density (ROP_s)'
               WRITE(*,30) ' 6 : Gas and solids temperature (T_g, T_s1, T_s2)'
               WRITE(*,30) ' 7 : Gas and solids mass fractions (X_g, X-s)'
               WRITE(*,30) ' 8 : Granular temperature (G)'
!               write(*,30) ' 9 : User defined scalars'
!               write(*,30) '10 : Reaction Rates'
               write(*,30) '11 : Turbulence quantities (k and ε)'
               write(*,30) '12 : Gas Vorticity magnitude and Lambda_2 (VORTICITY, LAMBDA_2)'
               write(*,30) '100: Processor assigned to scalar cell (Partition)'
               write(*,30) '101: Boundary condition flag for scalar cell (BC_ID)'
               write(*,30) 'MFiX will exit now.'
               CALL MFIX_EXIT(myPE) 

            END SELECT

      END DO

      CALL CLOSE_VTK_FILE


     IF (FULL_LOG.AND.myPE == PE_IO) WRITE(*,20)' DONE.'

10    FORMAT(A,$)
20    FORMAT(A,1X/)
30    FORMAT(1X,A)    
      RETURN
      
      END SUBROUTINE WRITE_VTK_FILE

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: OPEN_VTK_FILE                                          C
!  Purpose: Open a vtk file and writes the header                      C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE OPEN_VTK_FILE
    
      USE param
      USE param1
      USE parallel
      USE constant
      USE run
      USE toleranc
      USE geometry
      USE indices  
      USE compar
      USE sendrecv
      USE output
      USE quadric
      USE cutcell
      USE fldvar
      USE vtk
      
      IMPLICIT NONE
      DOUBLE PRECISION:: Xw,Xe,Yn,Ys
      INTEGER :: I,J,K,L,IM,JM,KM,IP,JP,KP,IJK
      INTEGER :: IMJK,IJMK,IJKM,IMJMK,IMJKM,IJMKM,IMJMKM

      include "function.inc"

      IF (myPE /= PE_IO) RETURN 

      IF(.NOT.WRITE_ANI_CUTCELL) THEN

         VTK_FILENAME = TRIM(RUN_NAME)

         IF(TIME_DEPENDENT_FILENAME) THEN
            FRAME = FRAME + 1 
            WRITE(FRAME_CHAR,*) FRAME
            FRAME_CHAR = ADJUSTL(FRAME_CHAR)
            VTK_FILENAME = TRIM(VTK_FILENAME) // '_' // TRIM(FRAME_CHAR) // '.vtk'
         ELSE
            VTK_FILENAME = TRIM(VTK_FILENAME) // '.vtk'
         ENDIF

         IF (FULL_LOG) THEN
            WRITE(*,10)' WRITING VTK FILE : ', TRIM(VTK_FILENAME),' .'
         ENDIF

      ELSE
  
         VTK_FILENAME = 'ani_cutcell.vtk'

      ENDIF

      VTK_UNIT   = 123

      OPEN(UNIT     = VTK_UNIT,           &
           FILE     = TRIM(VTK_FILENAME), &
           FORM     = 'UNFORMATTED',    &  ! works with gfortran 4.3.4 and ifort 10.1 but may not be supported by all compilers
                                           ! use 'BINARY' if 'UNFORMATTED' is not supported 
           ACCESS   = 'STREAM',   &        ! works with gfortran 4.3.4 and ifort 10.1 but may not be supported by all compilers
                                           ! use 'SEQUENTIAL' if 'STREAM' is not supported 
           ACTION   = 'WRITE', &
           CONVERT  = 'BIG_ENDIAN')

      WRITE(UNIT=VTK_UNIT)'# vtk DataFile Version 2.0'//END_REC
      WRITE(BUFFER,FMT='(A,A,E14.8)')TRIM(RUN_NAME),', Time = ',TIME
      WRITE(UNIT=VTK_UNIT)TRIM(BUFFER)//END_REC
      WRITE(UNIT=VTK_UNIT)TRIM('BINARY')//END_REC
      IF(NO_K) THEN
         WRITE(UNIT=VTK_UNIT)'DATASET POLYDATA'//END_REC
      ELSE
         WRITE(UNIT=VTK_UNIT)'DATASET UNSTRUCTURED_GRID'//END_REC
      ENDIF

10    FORMAT(/1X,3A,$)

      RETURN

      END SUBROUTINE OPEN_VTK_FILE


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_GEOMETRY_IN_VTK                                  C
!  Purpose: Write Geometry and connectivity in a vtk file              C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE WRITE_GEOMETRY_IN_VTK
    
      USE param
      USE param1
      USE parallel
      USE constant
      USE run
      USE toleranc
      USE geometry
      USE indices  
      USE compar
      USE mpi_utility 
      USE sendrecv
      USE quadric
      USE cutcell
      USE fldvar
      USE vtk

     
      IMPLICIT NONE

      INTEGER :: IJK,I,J,K,L
      INTEGER :: IJK_OFFSET

      INTEGER :: iproc,IERR
      INTEGER, DIMENSION(0:numPEs-1) :: disp,rcount
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: SHIFTED_CONNECTIVITY

      include "function.inc"

      IF (myPE /= PE_IO) RETURN

         NUMBER_OF_VTK_CELLS = NUMBER_OF_CELLS - NUMBER_OF_BLOCKED_CELLS

         WRITE(BUFFER,FMT='(A,I8,A)')'POINTS ',NUMBER_OF_POINTS,' double'
         WRITE(UNIT=VTK_UNIT)TRIM(BUFFER)//END_REC

         WRITE(UNIT=VTK_UNIT) (XG_E(GLOBAL_I_OF(IJK)),YG_N(GLOBAL_J_OF(IJK)),ZG_T(GLOBAL_K_OF(IJK)), IJK = 1,IJKMAX3), &
                              (GLOBAL_X_NEW_POINT(IJK),GLOBAL_Y_NEW_POINT(IJK),GLOBAL_Z_NEW_POINT(IJK),IJK = 1,&
                              GLOBAL_NUMBER_OF_NEW_POINTS)
         WRITE(UNIT=VTK_UNIT) END_REC

         IF(NO_K) THEN
            WRITE(BUFFER,FMT='(A,2(I8,2X))')'POLYGONS ',NUMBER_OF_VTK_CELLS,POLY_COUNTER
         ELSE
            WRITE(BUFFER,FMT='(A,2(I8,2X))')'CELLS ',NUMBER_OF_VTK_CELLS,POLY_COUNTER
         ENDIF

         WRITE(UNIT=VTK_UNIT)TRIM(BUFFER)//END_REC

         DO IJK = 1,IJKMAX3
            IF (GLOBAL_INTERIOR_CELL_AT(IJK))      THEN
              IF (.NOT.GLOBAL_BLOCKED_CELL_AT(IJK)) WRITE(UNIT=VTK_UNIT) GLOBAL_NUMBER_OF_NODES(IJK),&
              (GLOBAL_CONNECTIVITY(IJK,L)-1,L=1,GLOBAL_NUMBER_OF_NODES(IJK))
            ENDIF
         END DO
         WRITE(UNIT=VTK_UNIT) END_REC

         IF(DO_K) THEN
            WRITE(BUFFER,FMT='(A,I8)')'CELL_TYPES ',NUMBER_OF_VTK_CELLS
            WRITE(UNIT=VTK_UNIT)TRIM(BUFFER)//END_REC

            DO IJK = 1,IJKMAX3
               IF (GLOBAL_INTERIOR_CELL_AT(IJK))      THEN
                  IF (.NOT.GLOBAL_BLOCKED_CELL_AT(IJK)) THEN
                     IF(GLOBAL_CUT_CELL_AT(IJK)) THEN
                        WRITE(UNIT=VTK_UNIT) 41
                     ELSE
                        WRITE(UNIT=VTK_UNIT) 11
                     ENDIF
                  ENDIF
               ENDIF
            END DO
            WRITE(UNIT=VTK_UNIT) END_REC
         ENDIF

         WRITE(BUFFER,FMT='(A,I8)') 'CELL_DATA ',NUMBER_OF_VTK_CELLS
         WRITE(UNIT=VTK_UNIT)TRIM(BUFFER)//END_REC

      RETURN
      
      END SUBROUTINE WRITE_GEOMETRY_IN_VTK



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_SCALAR_IN_VTK                                    C
!  Purpose: Write Scalar variable in a vtk file                        C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE WRITE_SCALAR_IN_VTK(VAR_NAME,VAR)
    
      USE param
      USE param1
      USE parallel
      USE constant
      USE run
      USE toleranc
      USE geometry
      USE indices  
      USE compar
      USE mpi_utility 
      USE sendrecv
      USE quadric
      USE cutcell
      USE fldvar
      USE vtk
      
      IMPLICIT NONE
      INTEGER :: I,IJK,L

      CHARACTER (*) :: VAR_NAME
      DOUBLE PRECISION, DIMENSION(DIMENSION_3) ::  VAR
      DOUBLE PRECISION, ALLOCATABLE :: GLOBAL_VAR(:)

      include "function.inc"


      IF (myPE == PE_IO) THEN
         allocate (GLOBAL_VAR(ijkmax3))     
      ELSE
         allocate (GLOBAL_VAR(1))     
      ENDIF

      call gather (VAR,GLOBAL_VAR,root) 


      IF (myPE /= PE_IO) RETURN


      DO I = 1,LEN_TRIM(VAR_NAME)
         IF(VAR_NAME(I:I) == ' ') VAR_NAME(I:I) = '_'
      ENDDO


      WRITE(BUFFER,FMT='(A)')'SCALARS '//TRIM(VAR_NAME)//' double 1'
      WRITE(UNIT=VTK_UNIT)TRIM(BUFFER)//END_REC
      WRITE(BUFFER,FMT='(A)')'LOOKUP_TABLE default'
      WRITE(UNIT=VTK_UNIT)TRIM(BUFFER)//END_REC

      DO IJK = 1,IJKMAX3
         IF (GLOBAL_INTERIOR_CELL_AT(IJK))      THEN
            IF (.NOT.GLOBAL_BLOCKED_CELL_AT(IJK))   WRITE(UNIT=VTK_UNIT) GLOBAL_VAR(IJK)
         ENDIF
      ENDDO
      WRITE(UNIT=VTK_UNIT)END_REC

      Deallocate (GLOBAL_VAR)   

      RETURN
      
      END SUBROUTINE WRITE_SCALAR_IN_VTK


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_VECTOR_IN_VTK                                    C
!  Purpose: Write Vector variable in a vtk file                        C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE WRITE_VECTOR_IN_VTK(VAR_NAME,VARX,VARY,VARZ)
    
      USE param
      USE param1
      USE parallel
      USE constant
      USE run
      USE toleranc
      USE geometry
      USE indices  
      USE compar
      USE mpi_utility 
      USE sendrecv
      USE quadric
      USE cutcell
      USE fldvar
      USE vtk
      
      IMPLICIT NONE
      INTEGER :: IJK,L

      CHARACTER (*) :: VAR_NAME
      DOUBLE PRECISION, DIMENSION(DIMENSION_3) ::  VARX,VARY,VARZ
      DOUBLE PRECISION, ALLOCATABLE :: GLOBAL_VARX(:),GLOBAL_VARY(:),GLOBAL_VARZ(:)

      include "function.inc"


      IF (myPE == PE_IO) THEN
         allocate (GLOBAL_VARX(ijkmax3))
         allocate (GLOBAL_VARY(ijkmax3))     
         allocate (GLOBAL_VARZ(ijkmax3))          
      ELSE
         allocate (GLOBAL_VARX(1))
         allocate (GLOBAL_VARY(1))     
         allocate (GLOBAL_VARZ(1))          
      ENDIF

      call gather (VARX,GLOBAL_VARX,root)
      call gather (VARY,GLOBAL_VARY,root) 
      call gather (VARZ,GLOBAL_VARZ,root)  

      IF (myPE /= PE_IO) RETURN


      WRITE(BUFFER,FMT='(A)')'VECTORS '//TRIM(VAR_NAME)//' double'
      WRITE(UNIT=VTK_UNIT)TRIM(BUFFER)//END_REC

      DO IJK = 1,IJKMAX3
         IF (GLOBAL_INTERIOR_CELL_AT(IJK))      THEN
            IF (.NOT.GLOBAL_BLOCKED_CELL_AT(IJK))   WRITE(UNIT=VTK_UNIT) GLOBAL_VARX(IJK),GLOBAL_VARY(IJK),GLOBAL_VARZ(IJK)
         ENDIF
      ENDDO
      WRITE(UNIT=VTK_UNIT)END_REC


      Deallocate (GLOBAL_VARX)
      Deallocate (GLOBAL_VARY)   
      Deallocate (GLOBAL_VARZ)      

      RETURN
      
      END SUBROUTINE WRITE_VECTOR_IN_VTK


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CLOSE_VTK_FILE                                         C
!  Purpose: Close a vtk file                                           C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE CLOSE_VTK_FILE
    
      USE vtk

      IF (myPE /= PE_IO) RETURN 

      CLOSE(VTK_UNIT)

      RETURN

      END SUBROUTINE CLOSE_VTK_FILE





!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_CUT_SURFACE_VTK                                  C
!  Purpose: Writes the cut cell surface in VTK format                  C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE WRITE_CUT_SURFACE_VTK
    
      USE param
      USE param1
      USE parallel
      USE constant
      USE run
      USE toleranc
      USE geometry
      USE indices  
      USE compar
      USE sendrecv
      USE quadric
      USE cutcell
      USE fldvar
      USE vtk
      USE polygon
      USE stl
      
      IMPLICIT NONE

      INTEGER :: I,J,K,L,IM,JM,KM,IP,JP,KP,IJK,NODE
      INTEGER :: IMJK,IJMK,IJKM,IMJMK,IMJKM,IJMKM,IMJMKM
      INTEGER :: POINT_ID,POLY_COUNT,FACE_ID,Q_ID,Q_ID2
      INTEGER :: N_CUT_FACE_NODES,BCID2

      INTEGER NUMBER_OF_FACES
      INTEGER NUMBER_OF_SURFACE_POINTS

      DOUBLE PRECISION, DIMENSION(15,3) :: COORD_CUT_FACE_NODES
      DOUBLE PRECISION, DIMENSION(3)    :: NORMAL

      INTEGER, DIMENSION(DIMENSION_MAX_CUT_CELL,6) ::  FACE_CONNECTIVITY
      INTEGER, DIMENSION(DIMENSION_MAX_CUT_CELL)   ::  NUMBER_OF_CUT_FACE_POINTS

      DOUBLE PRECISION, DIMENSION(DIMENSION_MAX_CUT_CELL) ::  X_FACE_POINT
      DOUBLE PRECISION, DIMENSION(DIMENSION_MAX_CUT_CELL) ::  Y_FACE_POINT
      DOUBLE PRECISION, DIMENSION(DIMENSION_MAX_CUT_CELL) ::  Z_FACE_POINT

      DOUBLE PRECISION :: X_COPY,Y_COPY,Z_COPY,F_COPY,F2

      LOGICAL :: CLIP_FLAG,INTERSECT_FLAG,PRINT_FLAG

      CHARACTER (LEN=32) :: FILENAME

      LOGICAL :: CORNER_POINT
      INTEGER :: NODE_OF_CORNER

      include "function.inc"

      IF(myPE/=0) RETURN

!======================================================================
!  Set-up connectivity for each cell, i.e., regular cells and cut cells
!======================================================================

      POLY_COUNT = 0

      NUMBER_OF_SURFACE_POINTS = 0

      NUMBER_OF_FACES = 0

      DO IJK = 1,IJKMAX3

         IF(GLOBAL_CUT_CELL_AT(IJK)) THEN

!======================================================================
!  Filter the connectivity to identify nodes belonging to cut face
!======================================================================


            NUMBER_OF_FACES = NUMBER_OF_FACES + 1
     
            N_CUT_FACE_NODES = 0

            CALL GET_GLOBAL_CELL_NODE_COORDINATES(IJK,'SCALAR')

            DO L = 1, GLOBAL_NUMBER_OF_NODES(IJK)
               IF(GLOBAL_CONNECTIVITY(IJK,L)>IJKMAX3) THEN   ! One of the new point          
                  X_COPY = GLOBAL_X_NEW_POINT(GLOBAL_CONNECTIVITY(IJK,L)-IJKMAX3)
                  Y_COPY = GLOBAL_Y_NEW_POINT(GLOBAL_CONNECTIVITY(IJK,L)-IJKMAX3)
                  Z_COPY = GLOBAL_Z_NEW_POINT(GLOBAL_CONNECTIVITY(IJK,L)-IJKMAX3)
                  CORNER_POINT = .FALSE.
               ELSE                                   ! An existing point          
                  DO NODE = 1,8
                  CORNER_POINT = .TRUE.
                     IF(GLOBAL_CONNECTIVITY(IJK,L) == IJK_OF_NODE(NODE)) THEN
                        NODE_OF_CORNER = NODE
                        X_COPY = X_NODE(NODE)
                        Y_COPY = Y_NODE(NODE)
                        Z_COPY = Z_NODE(NODE)

                        IF (GLOBAL_SNAP(IJK_OF_NODE(NODE))) THEN ! One of the snapped corner point which now belongs to the cut face
                           N_CUT_FACE_NODES = N_CUT_FACE_NODES + 1
                           COORD_CUT_FACE_NODES(N_CUT_FACE_NODES,1) = X_COPY
                           COORD_CUT_FACE_NODES(N_CUT_FACE_NODES,2) = Y_COPY
                           COORD_CUT_FACE_NODES(N_CUT_FACE_NODES,3) = Z_COPY
                        ENDIF
                     ENDIF
                  END DO

               ENDIF




               IF(CORNER_POINT) THEN
                  Q_ID = 1

                  CALL EVAL_F('QUADRIC',X_COPY,Y_COPY,Z_COPY,Q_ID,F_COPY,CLIP_FLAG)

                  CALL EVAL_F('POLYGON',X_COPY,Y_COPY,Z_COPY,N_POLYGON,F_COPY,CLIP_FLAG)

                  CALL EVAL_F('USR_DEF',X_COPY,Y_COPY,Z_COPY,N_USR_DEF,F_COPY,CLIP_FLAG)

                  IF(USE_STL.OR.USE_MSH) F_COPY = GLOBAL_F_AT(IJK_OF_NODE(NODE_OF_CORNER))

!                  CALL EVAL_STL_FCT_AT('SCALAR',IJK,NODE_OF_CORNER,F_COPY,CLIP_FLAG,BCID2)
               ELSE
                  F_COPY = ZERO
               ENDIF

               IF (ABS(F_COPY) < TOL_F ) THEN ! belongs to cut face
                  N_CUT_FACE_NODES = N_CUT_FACE_NODES + 1
                  COORD_CUT_FACE_NODES(N_CUT_FACE_NODES,1) = X_COPY
                  COORD_CUT_FACE_NODES(N_CUT_FACE_NODES,2) = Y_COPY
                  COORD_CUT_FACE_NODES(N_CUT_FACE_NODES,3) = Z_COPY
               ENDIF

            END DO

            CALL REORDER_POLYGON(N_CUT_FACE_NODES,COORD_CUT_FACE_NODES,NORMAL)

            NUMBER_OF_CUT_FACE_POINTS(NUMBER_OF_FACES) = N_CUT_FACE_NODES 
            POLY_COUNT = POLY_COUNT + N_CUT_FACE_NODES + 1
            DO NODE = 1,N_CUT_FACE_NODES
               NUMBER_OF_SURFACE_POINTS = NUMBER_OF_SURFACE_POINTS + 1

               IF(NUMBER_OF_SURFACE_POINTS>=DIMENSION_MAX_CUT_CELL) THEN
                  WRITE(*,3000) 'ERROR IN SUBROUTINE WRITE_3DCUT_SURFACE_VTK:'
                  WRITE(*,3000) 'NUMBER_OF_SURFACE_POINTS>=DIMENSION_MAX_CUT_CELL:'
                  WRITE(*,3000) 'INCREASE VALUE OF FAC_DIM_MAX_CUT_CELL.'
                  WRITE(*,3010) 'CURRENT VALUE OF FAC_DIM_MAX_CUT_CELL =',FAC_DIM_MAX_CUT_CELL
                  WRITE(*,3020) 'CURRENT VALUE OF DIMENSION_MAX_CUT_CELL =',DIMENSION_MAX_CUT_CELL
                  WRITE(*,3000) 'MFiX will exit now.'
                  CALL MFIX_EXIT(myPE)
               ENDIF

               X_FACE_POINT(NUMBER_OF_SURFACE_POINTS) = COORD_CUT_FACE_NODES(NODE,1)
               Y_FACE_POINT(NUMBER_OF_SURFACE_POINTS) = COORD_CUT_FACE_NODES(NODE,2)
               Z_FACE_POINT(NUMBER_OF_SURFACE_POINTS) = COORD_CUT_FACE_NODES(NODE,3)
               FACE_CONNECTIVITY(NUMBER_OF_FACES,NODE) = NUMBER_OF_SURFACE_POINTS
            ENDDO

         ENDIF

      END DO



      FILENAME= TRIM(RUN_NAME) // '_boundary.vtk'
      FILENAME = TRIM(FILENAME)
      OPEN(UNIT = 123, FILE = FILENAME)
      WRITE(123,1001)'# vtk DataFile Version 2.0'
      WRITE(123,1001)'3D CUT-CELL SURFACE'
      WRITE(123,1001)'ASCII'

      IF(NO_K) THEN   ! 2D GEOMETRY
         WRITE(123,1001)'DATASET UNSTRUCTURED_GRID'      
      ELSE            ! 3D GEOMETRY
         WRITE(123,1001)'DATASET POLYDATA'      
      ENDIF

      WRITE(123,1010)'POINTS ',NUMBER_OF_SURFACE_POINTS,' float'

      DO POINT_ID = 1,NUMBER_OF_SURFACE_POINTS
         WRITE(123,1020) X_FACE_POINT(POINT_ID),Y_FACE_POINT(POINT_ID),Z_FACE_POINT(POINT_ID)
      ENDDO
     
      IF(NO_K) THEN   ! 2D GEOMETRY

         WRITE(123,1030)'CELLS ',NUMBER_OF_FACES,POLY_COUNT
         DO FACE_ID = 1 , NUMBER_OF_FACES
            WRITE(123,1040) NUMBER_OF_CUT_FACE_POINTS(FACE_ID),(FACE_CONNECTIVITY(FACE_ID,L)-1,&
            L=1,NUMBER_OF_CUT_FACE_POINTS(FACE_ID))         
         ENDDO
         WRITE(123,1030)'CELL_TYPES ',NUMBER_OF_FACES
         DO FACE_ID = 1 , NUMBER_OF_FACES
            WRITE(123,1040) 3
         ENDDO

      ELSE            ! 3D GEOMETRY
      
         WRITE(123,1030)'POLYGONS ',NUMBER_OF_FACES,POLY_COUNT
         DO FACE_ID = 1 , NUMBER_OF_FACES
            WRITE(123,1040) NUMBER_OF_CUT_FACE_POINTS(FACE_ID),(FACE_CONNECTIVITY(FACE_ID,L)-1,&
            L=1,NUMBER_OF_CUT_FACE_POINTS(FACE_ID))         
         ENDDO

      ENDIF

1001  FORMAT(A)
1010  FORMAT(A,I8,A)
1020  FORMAT(3(E16.8,2X))
1030  FORMAT(A,2(I8,2X))
1040  FORMAT(20(I8,2X))
1050  FORMAT(A,I8)
1060  FORMAT(E16.8)
1070  FORMAT(3(E16.8,2X))
1080  FORMAT(I5)
3000  FORMAT(1X,A) 
3010  FORMAT(1X,A,F8.4) 
3020  FORMAT(1X,A,I8) 
3030  FORMAT(1X,A,A) 
      CLOSE (123)


      WRITE(*,3030)'WROTE BOUNDARY IN VTK FILE : ',FILENAME
      RETURN

      
      END SUBROUTINE WRITE_CUT_SURFACE_VTK


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: GATHER_DATA                                            C
!  Purpose: Gather data from all processes in preparation of           C
!           Writing VTK files                                          C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE GATHER_DATA
    
      USE param
      USE param1
      USE parallel
      USE constant
      USE run
      USE toleranc
      USE geometry
      USE indices  
      USE compar
      USE mpi_utility 
      USE sendrecv
      USE quadric
      USE cutcell
      USE fldvar
      USE vtk

     
      IMPLICIT NONE

      INTEGER :: IJK,I,J,K,L
      INTEGER :: IJK_OFFSET

      INTEGER :: iproc,IERR
      INTEGER, DIMENSION(0:numPEs-1) :: disp,rcount
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: SHIFTED_CONNECTIVITY

      include "function.inc"


!======================================================================
!  parallel processing
!======================================================================

      CALL allgather_1i (NUMBER_OF_NEW_POINTS,rcount,IERR)

      IF (myPE == 0) THEN
         IJK_OFFSET = 0
      ELSE
         IJK_OFFSET = 0
         DO iproc=0,myPE-1
            IJK_OFFSET = IJK_OFFSET + rcount(iproc)
         ENDDO
      ENDIF

      CALL allgather_1i (IJK_OFFSET,disp,IERR)

      IF(.NOT.GLOBAL_VAR_ALLOCATED) THEN

         IF (myPE == PE_IO) THEN
            allocate (GLOBAL_I_OF(ijkmax3))     
            allocate (GLOBAL_J_OF(ijkmax3))     
            allocate (GLOBAL_K_OF(ijkmax3))     
            allocate (GLOBAL_CONNECTIVITY(ijkmax3,15)) 
            allocate (GLOBAL_NUMBER_OF_NODES(ijkmax3)) 
            allocate (GLOBAL_INTERIOR_CELL_AT(ijkmax3))
            allocate (GLOBAL_BLOCKED_CELL_AT(ijkmax3)) 
            allocate (GLOBAL_STANDARD_CELL_AT(ijkmax3))
            allocate (GLOBAL_CUT_CELL_AT(ijkmax3))    
            allocate (GLOBAL_SNAP(ijkmax3))    
            allocate (GLOBAL_X_NEW_POINT(ijkmax3))     
            allocate (GLOBAL_Y_NEW_POINT(ijkmax3))     
            allocate (GLOBAL_Z_NEW_POINT(ijkmax3))
            allocate (GLOBAL_F_AT(ijkmax3))         

         ELSE
            allocate (GLOBAL_I_OF(1))
            allocate (GLOBAL_J_OF(1))
            allocate (GLOBAL_K_OF(1))
            allocate (GLOBAL_CONNECTIVITY(1,15))     
            allocate (GLOBAL_NUMBER_OF_NODES(1))    
            allocate (GLOBAL_INTERIOR_CELL_AT(1))   
            allocate (GLOBAL_BLOCKED_CELL_AT(1))    
            allocate (GLOBAL_STANDARD_CELL_AT(1))   
            allocate (GLOBAL_CUT_CELL_AT(1))    
            allocate (GLOBAL_SNAP(1))     
            allocate (GLOBAL_X_NEW_POINT(1))     
            allocate (GLOBAL_Y_NEW_POINT(1))     
            allocate (GLOBAL_Z_NEW_POINT(1))     
            allocate (GLOBAL_F_AT(1))     
         ENDIF

         GLOBAL_VAR_ALLOCATED = .TRUE.
 
      ENDIF
      
      IF(numPEs==1) THEN  ! Serial run
         GLOBAL_X_NEW_POINT(1:NUMBER_OF_NEW_POINTS) =  X_NEW_POINT(1:NUMBER_OF_NEW_POINTS)
         GLOBAL_Y_NEW_POINT(1:NUMBER_OF_NEW_POINTS) =  Y_NEW_POINT(1:NUMBER_OF_NEW_POINTS)
         IF(DO_K) GLOBAL_Z_NEW_POINT(1:NUMBER_OF_NEW_POINTS) =  Z_NEW_POINT(1:NUMBER_OF_NEW_POINTS)
      ELSE !Parallel run
         call gatherv_1d( X_NEW_POINT, NUMBER_OF_NEW_POINTS, GLOBAL_X_NEW_POINT, rcount, disp, PE_IO, ierr )
         call gatherv_1d( Y_NEW_POINT, NUMBER_OF_NEW_POINTS, GLOBAL_Y_NEW_POINT, rcount, disp, PE_IO, ierr )
         call gatherv_1d( Z_NEW_POINT, NUMBER_OF_NEW_POINTS, GLOBAL_Z_NEW_POINT, rcount, disp, PE_IO, ierr )
      ENDIF

      call global_sum(NUMBER_OF_NEW_POINTS, GLOBAL_NUMBER_OF_NEW_POINTS,  PE_IO, ierr )

      Allocate(  SHIFTED_CONNECTIVITY  (DIMENSION_3,15) )

      SHIFTED_CONNECTIVITY = CONNECTIVITY

      WHERE (SHIFTED_CONNECTIVITY > IJKEND3)
         SHIFTED_CONNECTIVITY = SHIFTED_CONNECTIVITY - IJKEND3 + IJKMAX3 + disp(myPE)
      END WHERE

      DO IJK = IJKSTART3,IJKEND3
         DO L=1,NUMBER_OF_NODES(IJK)
            IF(CONNECTIVITY(IJK,L) <= IJKEND3) THEN
               I = I_OF(CONNECTIVITY(IJK,L))
               J = J_OF(CONNECTIVITY(IJK,L))
               K = K_OF(CONNECTIVITY(IJK,L))
               SHIFTED_CONNECTIVITY(IJK,L) = funijk_gl(I,J,K) 
            ENDIF
         ENDDO
      ENDDO

      call gather (I_OF,GLOBAL_I_OF,root)    
      call gather (J_OF,GLOBAL_J_OF,root)    
      call gather (K_OF,GLOBAL_K_OF,root)    
      call gather (SHIFTED_CONNECTIVITY,GLOBAL_CONNECTIVITY,root)    
      call gather (NUMBER_OF_NODES,GLOBAL_NUMBER_OF_NODES,root)    
      call gather (INTERIOR_CELL_AT,GLOBAL_INTERIOR_CELL_AT,root)  
      call gather (BLOCKED_CELL_AT,GLOBAL_BLOCKED_CELL_AT,root)    
      call gather (STANDARD_CELL_AT,GLOBAL_STANDARD_CELL_AT,root)  
      call gather (CUT_CELL_AT,GLOBAL_CUT_CELL_AT,root)   
      call gather (SNAP,GLOBAL_SNAP,root)   
      call gather (F_AT,GLOBAL_F_AT,root)   

      Deallocate(  SHIFTED_CONNECTIVITY )

      IF (myPE == PE_IO) THEN

         POLY_COUNTER = 0
  
         NUMBER_OF_CELLS = 0

         NUMBER_OF_CUT_CELLS = 0

         NUMBER_OF_BLOCKED_CELLS = 0

         NUMBER_OF_STANDARD_CELLS = 0

         DO IJK = 1, IJKMAX3   

            IF (GLOBAL_INTERIOR_CELL_AT(IJK)) THEN

               NUMBER_OF_CELLS = NUMBER_OF_CELLS + 1

               IF (GLOBAL_BLOCKED_CELL_AT(IJK))  NUMBER_OF_BLOCKED_CELLS  = NUMBER_OF_BLOCKED_CELLS + 1
               IF (GLOBAL_STANDARD_CELL_AT(IJK)) NUMBER_OF_STANDARD_CELLS = NUMBER_OF_STANDARD_CELLS + 1
               IF (GLOBAL_CUT_CELL_AT(IJK))      NUMBER_OF_CUT_CELLS = NUMBER_OF_CUT_CELLS + 1
    
               IF (.NOT.GLOBAL_BLOCKED_CELL_AT(IJK))   POLY_COUNTER = POLY_COUNTER + GLOBAL_NUMBER_OF_NODES(IJK) + 1

            ENDIF

         END DO


         NUMBER_OF_POINTS = IJKMAX3 + GLOBAL_NUMBER_OF_NEW_POINTS

      ENDIF

      RETURN

      
      END SUBROUTINE GATHER_DATA


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: PRINT_GRID_STATISTICS                                  C
!  Purpose: PRINT_GRID_STATISTICS ON SCREEN                            C
!                                                                      C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE PRINT_GRID_STATISTICS
    
      USE param
      USE param1
      USE parallel
      USE constant
      USE run
      USE toleranc
      USE geometry
      USE indices  
      USE compar
      USE mpi_utility 
      USE sendrecv
      USE quadric
      USE cutcell
      USE fldvar
      USE vtk

     
      IMPLICIT NONE

      INTEGER :: IJK,I,J,K,L
      INTEGER :: IJK_OFFSET

      INTEGER :: iproc,IERR

      DOUBLE PRECISION :: MIN_VOL, MAX_VOL, GLOBAL_MIN_VOL,GLOBAL_MAX_VOL
      DOUBLE PRECISION :: MIN_AYZ, MAX_AYZ, GLOBAL_MIN_AYZ,GLOBAL_MAX_AYZ
      DOUBLE PRECISION :: MIN_AXZ, MAX_AXZ, GLOBAL_MIN_AXZ,GLOBAL_MAX_AXZ
      DOUBLE PRECISION :: MIN_AXY, MAX_AXY, GLOBAL_MIN_AXY,GLOBAL_MAX_AXY
      DOUBLE PRECISION :: MIN_CUT, MAX_CUT, GLOBAL_MIN_CUT,GLOBAL_MAX_CUT
      DOUBLE PRECISION :: LOCAL_MIN_Q,LOCAL_MAX_Q, GLOBAL_MIN_Q,GLOBAL_MAX_Q


      include "function.inc"

      IF (myPE == PE_IO) THEN

         IF(.NOT.GRID_INFO_PRINTED_ON_SCREEN) THEN
            WRITE(*,5) 'GRID STATISTICS:'
            WRITE(*,5) 'NUMBER OF CELLS          = ', NUMBER_OF_CELLS 
            WRITE(*,10)'NUMBER OF STANDARD CELLS = ', &
                        NUMBER_OF_STANDARD_CELLS,DFLOAT(NUMBER_OF_STANDARD_CELLS) / DFLOAT(NUMBER_OF_CELLS) * 100.0D0
            WRITE(*,10)'NUMBER OF CUT CELLS      = ', &
                        NUMBER_OF_CUT_CELLS,DFLOAT(NUMBER_OF_CUT_CELLS) / DFLOAT(NUMBER_OF_CELLS) * 100.0D0
            WRITE(*,10)'NUMBER OF BLOCKED CELLS  = ', &
                        NUMBER_OF_BLOCKED_CELLS,DFLOAT(NUMBER_OF_BLOCKED_CELLS) / DFLOAT(NUMBER_OF_CELLS) * 100.0D0

5           FORMAT(1X,A,I8)
10          FORMAT(1X,A,I8,' (',F6.2,' % of Total)')

         ENDIF

         GRID_INFO_PRINTED_ON_SCREEN = .TRUE.

      ENDIF


!======================================================================
!  Scalar Cell volumes and areas
!======================================================================

      MIN_VOL =   LARGE_NUMBER
      MAX_VOL = - LARGE_NUMBER
      MIN_AYZ =   LARGE_NUMBER 
      MAX_AYZ = - LARGE_NUMBER 
      MIN_AXZ =   LARGE_NUMBER 
      MAX_AXZ = - LARGE_NUMBER
      MIN_AXY =   LARGE_NUMBER
      MAX_AXY = - LARGE_NUMBER

      DO IJK = IJKSTART3, IJKEND3
         IF(STANDARD_CELL_AT(IJK)) THEN              ! STANDARD CELLS

            MIN_VOL =   DMIN1(MIN_VOL,VOL(IJK))
            MAX_VOL =   DMAX1(MAX_VOL,VOL(IJK))
            MIN_AYZ =   DMIN1(MIN_AYZ,AYZ(IJK))
            MAX_AYZ =   DMAX1(MAX_AYZ,AYZ(IJK))
            MIN_AXZ =   DMIN1(MIN_AXZ,AXZ(IJK))
            MAX_AXZ =   DMAX1(MAX_AXZ,AXZ(IJK))
            MIN_AXY =   DMIN1(MIN_AXY,AXY(IJK))
            MAX_AXY =   DMAX1(MAX_AXY,AXY(IJK))

         ENDIF
      END DO

      call global_min(MIN_VOL, GLOBAL_MIN_VOL,  PE_IO, ierr )
      call global_max(MAX_VOL, GLOBAL_MAX_VOL,  PE_IO, ierr )
      call global_min(MIN_AYZ, GLOBAL_MIN_AYZ,  PE_IO, ierr )
      call global_max(MAX_AYZ, GLOBAL_MAX_AYZ,  PE_IO, ierr )
      call global_min(MIN_AXZ, GLOBAL_MIN_AXZ,  PE_IO, ierr )
      call global_max(MAX_AXZ, GLOBAL_MAX_AXZ,  PE_IO, ierr )
      call global_min(MIN_AXY, GLOBAL_MIN_AXY,  PE_IO, ierr )
      call global_max(MAX_AXY, GLOBAL_MAX_AXY,  PE_IO, ierr )

      IF (myPE == PE_IO) THEN
         WRITE(UNIT_CUT_CELL_LOG,1000)  '################################################################'
         WRITE(UNIT_CUT_CELL_LOG,1000)  '                       CELLS STATISTICS                         '
         WRITE(UNIT_CUT_CELL_LOG,1000)  '################################################################'
         WRITE(UNIT_CUT_CELL_LOG,1000)  'SCALAR STANDARD CELLS:'
         WRITE(UNIT_CUT_CELL_LOG,1000)  'RANGE OF AXY                   = ', GLOBAL_MIN_AXY,GLOBAL_MAX_AXY
         WRITE(UNIT_CUT_CELL_LOG,1000)  'RANGE OF AXZ                   = ', GLOBAL_MIN_AXZ,GLOBAL_MAX_AXZ
         WRITE(UNIT_CUT_CELL_LOG,1000)  'RANGE OF AYZ                   = ', GLOBAL_MIN_AYZ,GLOBAL_MAX_AYZ
         WRITE(UNIT_CUT_CELL_LOG,1000)  'RANGE OF VOLUME                = ', GLOBAL_MIN_VOL,GLOBAL_MAX_VOL
      ENDIF


      MIN_VOL =   LARGE_NUMBER
      MAX_VOL = - LARGE_NUMBER
      MIN_AYZ =   LARGE_NUMBER 
      MAX_AYZ = - LARGE_NUMBER 
      MIN_AXZ =   LARGE_NUMBER 
      MAX_AXZ = - LARGE_NUMBER
      MIN_AXY =   LARGE_NUMBER
      MAX_AXY = - LARGE_NUMBER

      DO IJK = IJKSTART3, IJKEND3
         IF(CUT_CELL_AT(IJK)) THEN                   ! CUT CELLS

            MIN_VOL =   DMIN1(MIN_VOL,VOL(IJK))
            MAX_VOL =   DMAX1(MAX_VOL,VOL(IJK))
            MIN_AYZ =   DMIN1(MIN_AYZ,AYZ(IJK))
            MAX_AYZ =   DMAX1(MAX_AYZ,AYZ(IJK))
            MIN_AXZ =   DMIN1(MIN_AXZ,AXZ(IJK))
            MAX_AXZ =   DMAX1(MAX_AXZ,AXZ(IJK))
            MIN_AXY =   DMIN1(MIN_AXY,AXY(IJK))
            MAX_AXY =   DMAX1(MAX_AXY,AXY(IJK))

         ENDIF
      END DO

      call global_min(MIN_VOL, GLOBAL_MIN_VOL,  PE_IO, ierr )
      call global_max(MAX_VOL, GLOBAL_MAX_VOL,  PE_IO, ierr )
      call global_min(MIN_AYZ, GLOBAL_MIN_AYZ,  PE_IO, ierr )
      call global_max(MAX_AYZ, GLOBAL_MAX_AYZ,  PE_IO, ierr )
      call global_min(MIN_AXZ, GLOBAL_MIN_AXZ,  PE_IO, ierr )
      call global_max(MAX_AXZ, GLOBAL_MAX_AXZ,  PE_IO, ierr )
      call global_min(MIN_AXY, GLOBAL_MIN_AXY,  PE_IO, ierr )
      call global_max(MAX_AXY, GLOBAL_MAX_AXY,  PE_IO, ierr )

      IF (myPE == PE_IO) THEN
         WRITE(UNIT_CUT_CELL_LOG,1000)  'SCALAR CUT CELLS:'
         WRITE(UNIT_CUT_CELL_LOG,1000)  'RANGE OF AXY                   = ', GLOBAL_MIN_AXY,GLOBAL_MAX_AXY
         WRITE(UNIT_CUT_CELL_LOG,1000)  'RANGE OF AXZ                   = ', GLOBAL_MIN_AXZ,GLOBAL_MAX_AXZ
         WRITE(UNIT_CUT_CELL_LOG,1000)  'RANGE OF AYZ                   = ', GLOBAL_MIN_AYZ,GLOBAL_MAX_AYZ
         WRITE(UNIT_CUT_CELL_LOG,1000)  'RANGE OF VOLUME                = ', GLOBAL_MIN_VOL,GLOBAL_MAX_VOL
         WRITE(UNIT_CUT_CELL_LOG,1010)  'NUMBER OF SMALL SCALAR CELLS   = ', NUMBER_OF_SMALL_CELLS 
         WRITE(UNIT_CUT_CELL_LOG,1000)  '################################################################'
      ENDIF


1000 FORMAT(A,E14.8,2X,E14.8)
1010 FORMAT(A,I8)

!======================================================================
!  U-Momentum Cell volumes and areas
!======================================================================



      MIN_VOL =   LARGE_NUMBER
      MAX_VOL = - LARGE_NUMBER
      MIN_AYZ =   LARGE_NUMBER 
      MAX_AYZ = - LARGE_NUMBER 
      MIN_AXZ =   LARGE_NUMBER 
      MAX_AXZ = - LARGE_NUMBER
      MIN_AXY =   LARGE_NUMBER
      MAX_AXY = - LARGE_NUMBER

      DO IJK = IJKSTART3, IJKEND3
         IF(STANDARD_U_CELL_AT(IJK)) THEN              ! STANDARD CELLS

            MIN_VOL =   DMIN1(MIN_VOL,VOL_U(IJK))
            MAX_VOL =   DMAX1(MAX_VOL,VOL_U(IJK))
            MIN_AYZ =   DMIN1(MIN_AYZ,AYZ_U(IJK))
            MAX_AYZ =   DMAX1(MAX_AYZ,AYZ_U(IJK))
            MIN_AXZ =   DMIN1(MIN_AXZ,AXZ_U(IJK))
            MAX_AXZ =   DMAX1(MAX_AXZ,AXZ_U(IJK))
            MIN_AXY =   DMIN1(MIN_AXY,AXY_U(IJK))
            MAX_AXY =   DMAX1(MAX_AXY,AXY_U(IJK))

         ENDIF
      END DO

      call global_min(MIN_VOL, GLOBAL_MIN_VOL,  PE_IO, ierr )
      call global_max(MAX_VOL, GLOBAL_MAX_VOL,  PE_IO, ierr )
      call global_min(MIN_AYZ, GLOBAL_MIN_AYZ,  PE_IO, ierr )
      call global_max(MAX_AYZ, GLOBAL_MAX_AYZ,  PE_IO, ierr )
      call global_min(MIN_AXZ, GLOBAL_MIN_AXZ,  PE_IO, ierr )
      call global_max(MAX_AXZ, GLOBAL_MAX_AXZ,  PE_IO, ierr )
      call global_min(MIN_AXY, GLOBAL_MIN_AXY,  PE_IO, ierr )
      call global_max(MAX_AXY, GLOBAL_MAX_AXY,  PE_IO, ierr )

      IF (myPE == PE_IO) THEN
         WRITE(UNIT_CUT_CELL_LOG,1000)  'U-MOMENTUM STANDARD CELLS:'
         WRITE(UNIT_CUT_CELL_LOG,1000)  'RANGE OF AXY                   = ', GLOBAL_MIN_AXY,GLOBAL_MAX_AXY
         WRITE(UNIT_CUT_CELL_LOG,1000)  'RANGE OF AXZ                   = ', GLOBAL_MIN_AXZ,GLOBAL_MAX_AXZ
         WRITE(UNIT_CUT_CELL_LOG,1000)  'RANGE OF AYZ                   = ', GLOBAL_MIN_AYZ,GLOBAL_MAX_AYZ
         WRITE(UNIT_CUT_CELL_LOG,1000)  'RANGE OF VOLUME                = ', GLOBAL_MIN_VOL,GLOBAL_MAX_VOL
      ENDIF

      MIN_VOL =   LARGE_NUMBER
      MAX_VOL = - LARGE_NUMBER
      MIN_AYZ =   LARGE_NUMBER 
      MAX_AYZ = - LARGE_NUMBER 
      MIN_AXZ =   LARGE_NUMBER 
      MAX_AXZ = - LARGE_NUMBER
      MIN_AXY =   LARGE_NUMBER
      MAX_AXY = - LARGE_NUMBER
      MIN_CUT =   LARGE_NUMBER
      MAX_CUT = - LARGE_NUMBER

      DO IJK = IJKSTART3, IJKEND3
         IF(CUT_U_CELL_AT(IJK).AND.(.NOT.WALL_U_AT(IJK))) THEN      ! CUT CELLS

            MIN_VOL =   DMIN1(MIN_VOL,VOL_U(IJK))
            MAX_VOL =   DMAX1(MAX_VOL,VOL_U(IJK))
            MIN_AYZ =   DMIN1(MIN_AYZ,AYZ_U(IJK))
            MAX_AYZ =   DMAX1(MAX_AYZ,AYZ_U(IJK))
            MIN_AXZ =   DMIN1(MIN_AXZ,AXZ_U(IJK))
            MAX_AXZ =   DMAX1(MAX_AXZ,AXZ_U(IJK))
            MIN_AXY =   DMIN1(MIN_AXY,AXY_U(IJK))
            MAX_AXY =   DMAX1(MAX_AXY,AXY_U(IJK))
            MIN_CUT =   DMIN1(MIN_CUT,AREA_U_CUT(IJK))
            MAX_CUT =   DMAX1(MAX_CUT,AREA_U_CUT(IJK))

         ENDIF
      END DO

      call global_min(MIN_VOL, GLOBAL_MIN_VOL,  PE_IO, ierr )
      call global_max(MAX_VOL, GLOBAL_MAX_VOL,  PE_IO, ierr )
      call global_min(MIN_AYZ, GLOBAL_MIN_AYZ,  PE_IO, ierr )
      call global_max(MAX_AYZ, GLOBAL_MAX_AYZ,  PE_IO, ierr )
      call global_min(MIN_AXZ, GLOBAL_MIN_AXZ,  PE_IO, ierr )
      call global_max(MAX_AXZ, GLOBAL_MAX_AXZ,  PE_IO, ierr )
      call global_min(MIN_AXY, GLOBAL_MIN_AXY,  PE_IO, ierr )
      call global_max(MAX_AXY, GLOBAL_MAX_AXY,  PE_IO, ierr )
      call global_min(MIN_CUT, GLOBAL_MIN_CUT,  PE_IO, ierr )
      call global_max(MAX_CUT, GLOBAL_MAX_CUT,  PE_IO, ierr )

      IF (myPE == PE_IO) THEN
         WRITE(UNIT_CUT_CELL_LOG,1000)  'U-MOMENTUM CUT CELLS:'
         WRITE(UNIT_CUT_CELL_LOG,1000)  'RANGE OF AXY                   = ', GLOBAL_MIN_AXY,GLOBAL_MAX_AXY
         WRITE(UNIT_CUT_CELL_LOG,1000)  'RANGE OF AXZ                   = ', GLOBAL_MIN_AXZ,GLOBAL_MAX_AXZ
         WRITE(UNIT_CUT_CELL_LOG,1000)  'RANGE OF AYZ                   = ', GLOBAL_MIN_AYZ,GLOBAL_MAX_AYZ
         WRITE(UNIT_CUT_CELL_LOG,1000)  'RANGE OF CUT AREA              = ', GLOBAL_MIN_CUT,GLOBAL_MAX_CUT
         WRITE(UNIT_CUT_CELL_LOG,1000)  'RANGE OF VOLUME                = ', GLOBAL_MIN_VOL,GLOBAL_MAX_VOL
         WRITE(UNIT_CUT_CELL_LOG,1010)  'NUMBER OF U WALL CELLS         = ', NUMBER_OF_U_WALL_CELLS 
         WRITE(UNIT_CUT_CELL_LOG,1000)  '################################################################'
      ENDIF
!======================================================================
!  V-Momentum Cell volumes and areas
!======================================================================


      MIN_VOL =   LARGE_NUMBER
      MAX_VOL = - LARGE_NUMBER
      MIN_AYZ =   LARGE_NUMBER 
      MAX_AYZ = - LARGE_NUMBER 
      MIN_AXZ =   LARGE_NUMBER 
      MAX_AXZ = - LARGE_NUMBER
      MIN_AXY =   LARGE_NUMBER
      MAX_AXY = - LARGE_NUMBER

      DO IJK = IJKSTART3, IJKEND3
         IF(STANDARD_V_CELL_AT(IJK)) THEN              ! STANDARD CELLS

            MIN_VOL =   DMIN1(MIN_VOL,VOL_V(IJK))
            MAX_VOL =   DMAX1(MAX_VOL,VOL_V(IJK))
            MIN_AYZ =   DMIN1(MIN_AYZ,AYZ_V(IJK))
            MAX_AYZ =   DMAX1(MAX_AYZ,AYZ_V(IJK))
            MIN_AXZ =   DMIN1(MIN_AXZ,AXZ_V(IJK))
            MAX_AXZ =   DMAX1(MAX_AXZ,AXZ_V(IJK))
            MIN_AXY =   DMIN1(MIN_AXY,AXY_V(IJK))
            MAX_AXY =   DMAX1(MAX_AXY,AXY_V(IJK))

         ENDIF
      END DO

      call global_min(MIN_VOL, GLOBAL_MIN_VOL,  PE_IO, ierr )
      call global_max(MAX_VOL, GLOBAL_MAX_VOL,  PE_IO, ierr )
      call global_min(MIN_AYZ, GLOBAL_MIN_AYZ,  PE_IO, ierr )
      call global_max(MAX_AYZ, GLOBAL_MAX_AYZ,  PE_IO, ierr )
      call global_min(MIN_AXZ, GLOBAL_MIN_AXZ,  PE_IO, ierr )
      call global_max(MAX_AXZ, GLOBAL_MAX_AXZ,  PE_IO, ierr )
      call global_min(MIN_AXY, GLOBAL_MIN_AXY,  PE_IO, ierr )
      call global_max(MAX_AXY, GLOBAL_MAX_AXY,  PE_IO, ierr )

      IF (myPE == PE_IO) THEN
         WRITE(UNIT_CUT_CELL_LOG,1000)  'V-MOMENTUM STANDARD CELLS:'
         WRITE(UNIT_CUT_CELL_LOG,1000)  'RANGE OF AXY                   = ', GLOBAL_MIN_AXY,GLOBAL_MAX_AXY
         WRITE(UNIT_CUT_CELL_LOG,1000)  'RANGE OF AXZ                   = ', GLOBAL_MIN_AXZ,GLOBAL_MAX_AXZ
         WRITE(UNIT_CUT_CELL_LOG,1000)  'RANGE OF AYZ                   = ', GLOBAL_MIN_AYZ,GLOBAL_MAX_AYZ
         WRITE(UNIT_CUT_CELL_LOG,1000)  'RANGE OF VOLUME                = ', GLOBAL_MIN_VOL,GLOBAL_MAX_VOL
      ENDIF

      MIN_VOL =   LARGE_NUMBER
      MAX_VOL = - LARGE_NUMBER
      MIN_AYZ =   LARGE_NUMBER 
      MAX_AYZ = - LARGE_NUMBER 
      MIN_AXZ =   LARGE_NUMBER 
      MAX_AXZ = - LARGE_NUMBER
      MIN_AXY =   LARGE_NUMBER
      MAX_AXY = - LARGE_NUMBER
      MIN_CUT =   LARGE_NUMBER
      MAX_CUT = - LARGE_NUMBER

      DO IJK = IJKSTART3, IJKEND3
         IF(CUT_V_CELL_AT(IJK).AND.(.NOT.WALL_V_AT(IJK))) THEN      ! CUT CELLS

            MIN_VOL =   DMIN1(MIN_VOL,VOL_V(IJK))
            MAX_VOL =   DMAX1(MAX_VOL,VOL_V(IJK))
            MIN_AYZ =   DMIN1(MIN_AYZ,AYZ_V(IJK))
            MAX_AYZ =   DMAX1(MAX_AYZ,AYZ_V(IJK))
            MIN_AXZ =   DMIN1(MIN_AXZ,AXZ_V(IJK))
            MAX_AXZ =   DMAX1(MAX_AXZ,AXZ_V(IJK))
            MIN_AXY =   DMIN1(MIN_AXY,AXY_V(IJK))
            MAX_AXY =   DMAX1(MAX_AXY,AXY_V(IJK))
            MIN_CUT =   DMIN1(MIN_CUT,AREA_V_CUT(IJK))
            MAX_CUT =   DMAX1(MAX_CUT,AREA_V_CUT(IJK))

         ENDIF
      END DO

      call global_min(MIN_VOL, GLOBAL_MIN_VOL,  PE_IO, ierr )
      call global_max(MAX_VOL, GLOBAL_MAX_VOL,  PE_IO, ierr )
      call global_min(MIN_AYZ, GLOBAL_MIN_AYZ,  PE_IO, ierr )
      call global_max(MAX_AYZ, GLOBAL_MAX_AYZ,  PE_IO, ierr )
      call global_min(MIN_AXZ, GLOBAL_MIN_AXZ,  PE_IO, ierr )
      call global_max(MAX_AXZ, GLOBAL_MAX_AXZ,  PE_IO, ierr )
      call global_min(MIN_AXY, GLOBAL_MIN_AXY,  PE_IO, ierr )
      call global_max(MAX_AXY, GLOBAL_MAX_AXY,  PE_IO, ierr )
      call global_min(MIN_CUT, GLOBAL_MIN_CUT,  PE_IO, ierr )
      call global_max(MAX_CUT, GLOBAL_MAX_CUT,  PE_IO, ierr )

      IF (myPE == PE_IO) THEN
         WRITE(UNIT_CUT_CELL_LOG,1000)  'V-MOMENTUM CUT CELLS:'
         WRITE(UNIT_CUT_CELL_LOG,1000)  'RANGE OF AXY                   = ', GLOBAL_MIN_AXY,GLOBAL_MAX_AXY
         WRITE(UNIT_CUT_CELL_LOG,1000)  'RANGE OF AXZ                   = ', GLOBAL_MIN_AXZ,GLOBAL_MAX_AXZ
         WRITE(UNIT_CUT_CELL_LOG,1000)  'RANGE OF AYZ                   = ', GLOBAL_MIN_AYZ,GLOBAL_MAX_AYZ
         WRITE(UNIT_CUT_CELL_LOG,1000)  'RANGE OF CUT AREA              = ', GLOBAL_MIN_CUT,GLOBAL_MAX_CUT
         WRITE(UNIT_CUT_CELL_LOG,1000)  'RANGE OF VOLUME                = ', GLOBAL_MIN_VOL,GLOBAL_MAX_VOL
         WRITE(UNIT_CUT_CELL_LOG,1010)  'NUMBER OF V WALL CELLS         = ', NUMBER_OF_V_WALL_CELLS 
         WRITE(UNIT_CUT_CELL_LOG,1000)  '################################################################'
      ENDIF

!======================================================================
!  W-Momentum Cell volumes and areas
!======================================================================


      IF(DO_K) THEN

         MIN_VOL =   LARGE_NUMBER
         MAX_VOL = - LARGE_NUMBER
         MIN_AYZ =   LARGE_NUMBER 
         MAX_AYZ = - LARGE_NUMBER 
         MIN_AXZ =   LARGE_NUMBER 
         MAX_AXZ = - LARGE_NUMBER
         MIN_AXY =   LARGE_NUMBER
         MAX_AXY = - LARGE_NUMBER

         DO IJK = IJKSTART3, IJKEND3
            IF(STANDARD_W_CELL_AT(IJK)) THEN              ! STANDARD CELLS

               MIN_VOL =   DMIN1(MIN_VOL,VOL_W(IJK))
               MAX_VOL =   DMAX1(MAX_VOL,VOL_W(IJK))
               MIN_AYZ =   DMIN1(MIN_AYZ,AYZ_W(IJK))
               MAX_AYZ =   DMAX1(MAX_AYZ,AYZ_W(IJK))
               MIN_AXZ =   DMIN1(MIN_AXZ,AXZ_W(IJK))
               MAX_AXZ =   DMAX1(MAX_AXZ,AXZ_W(IJK))
               MIN_AXY =   DMIN1(MIN_AXY,AXY_W(IJK))
               MAX_AXY =   DMAX1(MAX_AXY,AXY_W(IJK))

            ENDIF
         END DO

         call global_min(MIN_VOL, GLOBAL_MIN_VOL,  PE_IO, ierr )
         call global_max(MAX_VOL, GLOBAL_MAX_VOL,  PE_IO, ierr )
         call global_min(MIN_AYZ, GLOBAL_MIN_AYZ,  PE_IO, ierr )
         call global_max(MAX_AYZ, GLOBAL_MAX_AYZ,  PE_IO, ierr )
         call global_min(MIN_AXZ, GLOBAL_MIN_AXZ,  PE_IO, ierr )
         call global_max(MAX_AXZ, GLOBAL_MAX_AXZ,  PE_IO, ierr )
         call global_min(MIN_AXY, GLOBAL_MIN_AXY,  PE_IO, ierr )
         call global_max(MAX_AXY, GLOBAL_MAX_AXY,  PE_IO, ierr )

         IF (myPE == PE_IO) THEN
            WRITE(UNIT_CUT_CELL_LOG,1000)  'W-MOMENTUM STANDARD CELLS:'
            WRITE(UNIT_CUT_CELL_LOG,1000)  'RANGE OF AXY                   = ', GLOBAL_MIN_AXY,GLOBAL_MAX_AXY
            WRITE(UNIT_CUT_CELL_LOG,1000)  'RANGE OF AXZ                   = ', GLOBAL_MIN_AXZ,GLOBAL_MAX_AXZ
            WRITE(UNIT_CUT_CELL_LOG,1000)  'RANGE OF AYZ                   = ', GLOBAL_MIN_AYZ,GLOBAL_MAX_AYZ
            WRITE(UNIT_CUT_CELL_LOG,1000)  'RANGE OF VOLUME                = ', GLOBAL_MIN_VOL,GLOBAL_MAX_VOL
         ENDIF

         MIN_VOL =   LARGE_NUMBER
         MAX_VOL = - LARGE_NUMBER
         MIN_AYZ =   LARGE_NUMBER 
         MAX_AYZ = - LARGE_NUMBER 
         MIN_AXZ =   LARGE_NUMBER 
         MAX_AXZ = - LARGE_NUMBER
         MIN_AXY =   LARGE_NUMBER
         MAX_AXY = - LARGE_NUMBER
         MIN_CUT =   LARGE_NUMBER
         MAX_CUT = - LARGE_NUMBER

         DO IJK = IJKSTART3, IJKEND3
            IF(CUT_W_CELL_AT(IJK).AND.(.NOT.WALL_W_AT(IJK))) THEN      ! CUT CELLS

               MIN_VOL =   DMIN1(MIN_VOL,VOL_W(IJK))
               MAX_VOL =   DMAX1(MAX_VOL,VOL_W(IJK))
               MIN_AYZ =   DMIN1(MIN_AYZ,AYZ_W(IJK))
               MAX_AYZ =   DMAX1(MAX_AYZ,AYZ_W(IJK))
               MIN_AXZ =   DMIN1(MIN_AXZ,AXZ_W(IJK))
               MAX_AXZ =   DMAX1(MAX_AXZ,AXZ_W(IJK))
               MIN_AXY =   DMIN1(MIN_AXY,AXY_W(IJK))
               MAX_AXY =   DMAX1(MAX_AXY,AXY_W(IJK))
               MIN_CUT =   DMIN1(MIN_CUT,AREA_W_CUT(IJK))
               MAX_CUT =   DMAX1(MAX_CUT,AREA_W_CUT(IJK))

            ENDIF
         END DO

         call global_min(MIN_VOL, GLOBAL_MIN_VOL,  PE_IO, ierr )
         call global_max(MAX_VOL, GLOBAL_MAX_VOL,  PE_IO, ierr )
         call global_min(MIN_AYZ, GLOBAL_MIN_AYZ,  PE_IO, ierr )
         call global_max(MAX_AYZ, GLOBAL_MAX_AYZ,  PE_IO, ierr )
         call global_min(MIN_AXZ, GLOBAL_MIN_AXZ,  PE_IO, ierr )
         call global_max(MAX_AXZ, GLOBAL_MAX_AXZ,  PE_IO, ierr )
         call global_min(MIN_AXY, GLOBAL_MIN_AXY,  PE_IO, ierr )
         call global_max(MAX_AXY, GLOBAL_MAX_AXY,  PE_IO, ierr )
         call global_min(MIN_CUT, GLOBAL_MIN_CUT,  PE_IO, ierr )
         call global_max(MAX_CUT, GLOBAL_MAX_CUT,  PE_IO, ierr )

         IF (myPE == PE_IO) THEN
            WRITE(UNIT_CUT_CELL_LOG,1000)  'W-MOMENTUM CUT CELLS:'
            WRITE(UNIT_CUT_CELL_LOG,1000)  'RANGE OF AXY                   = ', GLOBAL_MIN_AXY,GLOBAL_MAX_AXY
            WRITE(UNIT_CUT_CELL_LOG,1000)  'RANGE OF AXZ                   = ', GLOBAL_MIN_AXZ,GLOBAL_MAX_AXZ
            WRITE(UNIT_CUT_CELL_LOG,1000)  'RANGE OF AYZ                   = ', GLOBAL_MIN_AYZ,GLOBAL_MAX_AYZ
            WRITE(UNIT_CUT_CELL_LOG,1000)  'RANGE OF CUT AREA              = ', GLOBAL_MIN_CUT,GLOBAL_MAX_CUT
            WRITE(UNIT_CUT_CELL_LOG,1000)  'RANGE OF VOLUME                = ', GLOBAL_MIN_VOL,GLOBAL_MAX_VOL
            WRITE(UNIT_CUT_CELL_LOG,1010)  'NUMBER OF W WALL CELLS         = ', NUMBER_OF_W_WALL_CELLS 
            WRITE(UNIT_CUT_CELL_LOG,1000)  '################################################################'
         ENDIF

      ENDIF



      LOCAL_MIN_Q = MINVAL(Alpha_Ue_c)
      LOCAL_MAX_Q = MAXVAL(Alpha_Ue_c)
      call global_min(LOCAL_MIN_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
      call global_max(LOCAL_MAX_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
      IF (myPE == PE_IO)  WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF U-MOMENTUM Alpha_Ue_c = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q

      LOCAL_MIN_Q = MINVAL(Alpha_Un_c)
      LOCAL_MAX_Q = MAXVAL(Alpha_Un_c)
      call global_min(LOCAL_MIN_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
      call global_max(LOCAL_MAX_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
      IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF U-MOMENTUM Alpha_Un_c = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q

      LOCAL_MIN_Q = MINVAL(Alpha_Ut_c)
      LOCAL_MAX_Q = MAXVAL(Alpha_Ut_c)
      call global_min(LOCAL_MIN_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
      call global_max(LOCAL_MAX_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
      IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF U-MOMENTUM Alpha_Ut_c = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q
      IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  

      LOCAL_MIN_Q = MINVAL(Theta_Ue)
      LOCAL_MAX_Q = MAXVAL(Theta_Ue)
      call global_min(LOCAL_MIN_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
      call global_max(LOCAL_MAX_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
      IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF U-MOMENTUM Theta_Ue   = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q

      LOCAL_MIN_Q = MINVAL(Theta_Un)
      LOCAL_MAX_Q = MAXVAL(Theta_Un)
      call global_min(LOCAL_MIN_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
      call global_max(LOCAL_MAX_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
      IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF U-MOMENTUM Theta_Un   = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q

      LOCAL_MIN_Q = MINVAL(Theta_Ut)
      LOCAL_MAX_Q = MAXVAL(Theta_Ut)
      call global_min(LOCAL_MIN_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
      call global_max(LOCAL_MAX_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
      IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF U-MOMENTUM Theta_Ut   = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q
      IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  

      LOCAL_MIN_Q = MINVAL(Theta_U_ne)
      LOCAL_MAX_Q = MAXVAL(Theta_U_ne)
      call global_min(LOCAL_MIN_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
      call global_max(LOCAL_MAX_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
      IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF U-MOMENTUM Theta_U_ne = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q

      LOCAL_MIN_Q = MINVAL(Theta_U_te)
      LOCAL_MAX_Q = MAXVAL(Theta_U_te)
      call global_min(LOCAL_MIN_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
      call global_max(LOCAL_MAX_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
      IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF U-MOMENTUM Theta_U_te = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q
      IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  

      LOCAL_MIN_Q = MINVAL(NOC_U_E)
      LOCAL_MAX_Q = MAXVAL(NOC_U_E)
      call global_min(LOCAL_MIN_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
      call global_max(LOCAL_MAX_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
      IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF U-MOMENTUM NOC_U_E    = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q

      LOCAL_MIN_Q = MINVAL(NOC_U_N)
      LOCAL_MAX_Q = MAXVAL(NOC_U_N)
      call global_min(LOCAL_MIN_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
      call global_max(LOCAL_MAX_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
      IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF U-MOMENTUM NOC_U_N    = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q

      LOCAL_MIN_Q = MINVAL(NOC_U_T)
      LOCAL_MAX_Q = MAXVAL(NOC_U_T)
      call global_min(LOCAL_MIN_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
      call global_max(LOCAL_MAX_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
      IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF U-MOMENTUM NOC_U_T    = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q
      IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  

      LOCAL_MIN_Q = MINVAL(DELH_U)
      LOCAL_MAX_Q = MAXVAL(DELH_U)
      call global_min(LOCAL_MIN_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
      call global_max(LOCAL_MAX_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
      IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF U-MOMENTUM DELH_U     = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q
      IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  '################################################################'



      LOCAL_MIN_Q = MINVAL(Alpha_Ve_c)
      LOCAL_MAX_Q = MAXVAL(Alpha_Ve_c)
      call global_min(LOCAL_MIN_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
      call global_max(LOCAL_MAX_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
      IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF V-MOMENTUM Alpha_Ve_c = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q
      LOCAL_MIN_Q = MINVAL(Alpha_Vn_c)
      LOCAL_MAX_Q = MAXVAL(Alpha_Vn_c)
      call global_min(LOCAL_MIN_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
      call global_max(LOCAL_MAX_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
      IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF V-MOMENTUM Alpha_Vn_c = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q
      LOCAL_MIN_Q = MINVAL(Alpha_Vt_c)
      LOCAL_MAX_Q = MAXVAL(Alpha_Vt_c)
      call global_min(LOCAL_MIN_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
      call global_max(LOCAL_MAX_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
      IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF V-MOMENTUM Alpha_Vt_c = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q
      IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  
      LOCAL_MIN_Q = MINVAL(Theta_Ve)
      LOCAL_MAX_Q = MAXVAL(Theta_Ve)
      call global_min(LOCAL_MIN_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
      call global_max(LOCAL_MAX_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
      IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF V-MOMENTUM Theta_Ve   = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q
      LOCAL_MIN_Q = MINVAL(Theta_Vn)
      LOCAL_MAX_Q = MAXVAL(Theta_Vn)
      call global_min(LOCAL_MIN_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
      call global_max(LOCAL_MAX_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
      IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF V-MOMENTUM Theta_Vn   = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q
      LOCAL_MIN_Q = MINVAL(Theta_Vt)
      LOCAL_MAX_Q = MAXVAL(Theta_Vt)
      call global_min(LOCAL_MIN_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
      call global_max(LOCAL_MAX_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
      IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF V-MOMENTUM Theta_Vt   = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q
      IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  
      LOCAL_MIN_Q = MINVAL(Theta_V_ne)
      LOCAL_MAX_Q = MAXVAL(Theta_V_ne)
      call global_min(LOCAL_MIN_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
      call global_max(LOCAL_MAX_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
      IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF V-MOMENTUM Theta_V_ne = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q
      LOCAL_MIN_Q = MINVAL(Theta_V_nt)
      LOCAL_MAX_Q = MAXVAL(Theta_V_nt)
      call global_min(LOCAL_MIN_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
      call global_max(LOCAL_MAX_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
      IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF V-MOMENTUM Theta_V_nt = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q
      IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  
      LOCAL_MIN_Q = MINVAL(NOC_V_E)
      LOCAL_MAX_Q = MAXVAL(NOC_V_E)
      call global_min(LOCAL_MIN_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
      call global_max(LOCAL_MAX_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
      IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF V-MOMENTUM NOC_V_E    = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q
      LOCAL_MIN_Q = MINVAL(NOC_V_N)
      LOCAL_MAX_Q = MAXVAL(NOC_V_N)
      call global_min(LOCAL_MIN_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
      call global_max(LOCAL_MAX_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
      IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF V-MOMENTUM NOC_V_N    = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q
      LOCAL_MIN_Q = MINVAL(NOC_V_T)
      LOCAL_MAX_Q = MAXVAL(NOC_V_T)
      call global_min(LOCAL_MIN_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
      call global_max(LOCAL_MAX_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
      IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF V-MOMENTUM NOC_V_T    = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q
      IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  
      LOCAL_MIN_Q = MINVAL(DELH_V)
      LOCAL_MAX_Q = MAXVAL(DELH_V)
      call global_min(LOCAL_MIN_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
      call global_max(LOCAL_MAX_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
      IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF V-MOMENTUM DELH_V     = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q
      IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  '################################################################'


      IF(DO_K) THEN

         LOCAL_MIN_Q = MINVAL(Alpha_We_c)
         LOCAL_MAX_Q = MAXVAL(Alpha_We_c)
         call global_min(LOCAL_MIN_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
         call global_max(LOCAL_MAX_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
         IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF W-MOMENTUM Alpha_We_c = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q
         LOCAL_MIN_Q = MINVAL(Alpha_Wn_c)
         LOCAL_MAX_Q = MAXVAL(Alpha_Wn_c)
         call global_min(LOCAL_MIN_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
         call global_max(LOCAL_MAX_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
         IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF W-MOMENTUM Alpha_Wn_c = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q
         LOCAL_MIN_Q = MINVAL(Alpha_Wt_c)
         LOCAL_MAX_Q = MAXVAL(Alpha_Wt_c)
         call global_min(LOCAL_MIN_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
         call global_max(LOCAL_MAX_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
         IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF W-MOMENTUM Alpha_Wt_c = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q
         IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  
         LOCAL_MIN_Q = MINVAL(Theta_We)
         LOCAL_MAX_Q = MAXVAL(Theta_We)
         call global_min(LOCAL_MIN_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
         call global_max(LOCAL_MAX_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
         IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF W-MOMENTUM Theta_We   = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q
         LOCAL_MIN_Q = MINVAL(Theta_Wn)
         LOCAL_MAX_Q = MAXVAL(Theta_Wn)
         call global_min(LOCAL_MIN_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
         call global_max(LOCAL_MAX_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
         IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF W-MOMENTUM Theta_Wn   = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q
         LOCAL_MIN_Q = MINVAL(Theta_Wt)
         LOCAL_MAX_Q = MAXVAL(Theta_Wt)
         call global_min(LOCAL_MIN_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
         call global_max(LOCAL_MAX_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
         IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF W-MOMENTUM Theta_Wt   = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q
         IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  
         LOCAL_MIN_Q = MINVAL(Theta_W_te)
         LOCAL_MAX_Q = MAXVAL(Theta_W_te)
         call global_min(LOCAL_MAX_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
         call global_max(LOCAL_MIN_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
         IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF W-MOMENTUM Theta_W_te = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q
         LOCAL_MIN_Q = MINVAL(Theta_W_tn)
         LOCAL_MAX_Q = MAXVAL(Theta_W_tn)
         call global_min(LOCAL_MIN_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
         call global_max(LOCAL_MAX_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
         IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF W-MOMENTUM Theta_W_tn = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q
         IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  
         LOCAL_MIN_Q = MINVAL(NOC_W_E)
         LOCAL_MAX_Q = MAXVAL(NOC_W_E)
         call global_min(LOCAL_MIN_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
         call global_max(LOCAL_MAX_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
         IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF W-MOMENTUM NOC_W_E    = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q
         LOCAL_MIN_Q = MINVAL(NOC_W_N)
         LOCAL_MAX_Q = MAXVAL(NOC_W_N)
         call global_min(LOCAL_MIN_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
         call global_max(LOCAL_MAX_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
         IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF W-MOMENTUM NOC_W_N    = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q
         LOCAL_MIN_Q = MINVAL(NOC_W_T)
         LOCAL_MAX_Q = MAXVAL(NOC_W_T)
         call global_min(LOCAL_MIN_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
         call global_max(LOCAL_MAX_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
         IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF W-MOMENTUM NOC_W_T    = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q
         IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  
         LOCAL_MIN_Q = MINVAL(DELH_W)
         LOCAL_MAX_Q = MAXVAL(DELH_W)
         call global_min(LOCAL_MIN_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
         call global_max(LOCAL_MAX_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
         IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF W-MOMENTUM DELH_W     = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q
         IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  '################################################################'
   
      ENDIF

      RETURN

      
      END SUBROUTINE PRINT_GRID_STATISTICS

