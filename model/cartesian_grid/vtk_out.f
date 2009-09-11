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

      USE pgcor
      USE pscor


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
               CALL WRITE_SCALAR_IN_VTK('PARTITION',PARTITION)
               IF (FULL_LOG.AND.myPE == PE_IO) WRITE(*,10)'.'
            
            CASE (101)
               CALL WRITE_SCALAR_IN_VTK('BC_ID',DFLOAT(BC_ID))
               IF (FULL_LOG.AND.myPE == PE_IO) WRITE(*,10)'.'
            

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
           FORM     = 'BINARY',    &
           ACCESS   = 'SEQUENTIAL',   &
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
      
      IMPLICIT NONE

      INTEGER :: I,J,K,L,IM,JM,KM,IP,JP,KP,IJK,NODE
      INTEGER :: IMJK,IJMK,IJKM,IMJMK,IMJKM,IJMKM,IMJMKM
      INTEGER :: POINT_ID,POLY_COUNT,FACE_ID,Q_ID,Q_ID2
      INTEGER :: N_CUT_FACE_NODES

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
               ELSE                                   ! An existing point          
                  DO NODE = 1,8
                     IF(GLOBAL_CONNECTIVITY(IJK,L) == IJK_OF_NODE(NODE)) THEN
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

               Q_ID = 1

               CALL EVAL_F('QUADRIC',X_COPY,Y_COPY,Z_COPY,Q_ID,F_COPY,CLIP_FLAG)

               CALL EVAL_F('POLYGON',X_COPY,Y_COPY,Z_COPY,N_POLYGON,F_COPY,CLIP_FLAG)

               CALL EVAL_F('USR_DEF',X_COPY,Y_COPY,Z_COPY,N_USR_DEF,F_COPY,CLIP_FLAG)


         
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

         ELSE
            allocate (GLOBAL_I_OF(1))
            allocate (GLOBAL_J_OF(1))
            allocate (GLOBAL_K_OF(1))
            allocate (GLOBAL_CONNECTIVITY(1,1))     
            allocate (GLOBAL_NUMBER_OF_NODES(1))    
            allocate (GLOBAL_INTERIOR_CELL_AT(1))   
            allocate (GLOBAL_BLOCKED_CELL_AT(1))    
            allocate (GLOBAL_STANDARD_CELL_AT(1))   
            allocate (GLOBAL_CUT_CELL_AT(1))    
            allocate (GLOBAL_SNAP(1))     
            allocate (GLOBAL_X_NEW_POINT(1))     
            allocate (GLOBAL_Y_NEW_POINT(1))     
            allocate (GLOBAL_Z_NEW_POINT(1))     
         ENDIF

         GLOBAL_VAR_ALLOCATED = .TRUE.
 
      ENDIF

      call gatherv_1d( X_NEW_POINT, NUMBER_OF_NEW_POINTS, GLOBAL_X_NEW_POINT, rcount, disp, PE_IO, ierr )
      call gatherv_1d( Y_NEW_POINT, NUMBER_OF_NEW_POINTS, GLOBAL_Y_NEW_POINT, rcount, disp, PE_IO, ierr )
      call gatherv_1d( Z_NEW_POINT, NUMBER_OF_NEW_POINTS, GLOBAL_Z_NEW_POINT, rcount, disp, PE_IO, ierr )

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

