!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_DES_DATA                                         C
!  Purpose: Writing DES output in Paraview format                      
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 26-Jul-06  C
!  Reviewer: Sreekanth Pannala                        Date: 31-Oct-06  C
!                                                                      !
!  Reviewer: J. Musser                                Date: 20-Apr-10  !
!  Comments: Split original subroutine into one for ParaView *.vtp     !
!  files, and a second for TECPLOT files *.dat.                        !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE WRITE_DES_DATA

      USE run
      USE discretelement
      USE des_bc
      IMPLICIT NONE

!-----------------------------------------------
! Local Variables
!-----------------------------------------------

!-----------------------------------------------
! Functions 
!-----------------------------------------------

!-----------------------------------------------

      IF (TRIM(DES_OUTPUT_TYPE) .EQ. 'TECPLOT') THEN
         CALL WRITE_DES_TECPLOT
      ELSE
         CALL WRITE_DES_VTP
      ENDIF

! Invoke at own risk      
      IF (.FALSE.) CALL WRITE_DES_THETA
      IF (.FALSE.) CALL WRITE_DES_BEDHEIGHT

      RETURN
      END SUBROUTINE WRITE_DES_DATA
!----------------------------------------------- 



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!
!  Module name: WRITE_DES_VTP
!  Purpose: Writing DES output in Paraview format
!
!
!  Reviewer: Rahul Garg                               Dare: 01-Aug-07
!  Comments: Added one more output file containing averaged bed height
!     
!  NOTE: If the system starts with zero particles, ParaView may have
!  trouble showing the results. To view the results in the current
!  version of ParaView, Version 3.6.1:
!    i - load the *.vtp files
!   ii - filter with glyph (Filters > Common > Glyph)
!        a - change glyph to sphere
!        b - change scale mode to scalar
!        c - check the "Edit" Set Scale Factor Box
!        d - change the value to 1.0
!     
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      SUBROUTINE WRITE_DES_VTP
    
      USE run
      USE compar
      USE funits
      USE discretelement
      USE des_bc
      IMPLICIT NONE

!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! logical that identifies that the data file has been created
! and is already opened (initial checks/writes)
      LOGICAL, SAVE :: FIRST_PASS = .TRUE.

! logical used for testing is the data file already exists
      LOGICAL :: F_EXISTS      

! file unit for ParaView *.vtp data      
      INTEGER, PARAMETER :: DES_UNIT = 2000

! file unit for ParaView *.vbd data      
      INTEGER, PARAMETER :: PVD_UNIT = 2050  

! index used when writing DES_*.vtp file
      CHARACTER*5 F_INDEX

! formatted file name
      CHARACTER*64 :: FNAME_VTP = ''
      CHARACTER*64 :: FNAME_PVD = ''      

! formatted solids time
      CHARACTER*12 :: S_TIME_CHAR = ''      

! index to track accounted for particles
      INTEGER PC

! dummy index values
      INTEGER L, I, J, K, M, IJK

! dummy values to maintain format for dimn=2
      REAL POS_Z, VEL_W 
 
! check whether an error occurs in opening a file
      INTEGER ISTAT      
!-----------------------------------------------
! Functions 
!-----------------------------------------------
!-----------------------------------------------


! Convert the index VTP_FINDEX from an integer to a string and force
! leading zeros


! Force file name format
      WRITE (F_INDEX,"(I5.5)") VTP_FINDEX
      FNAME_VTP=TRIM(RUN_NAME)//'_DES_'//TRIM(F_INDEX)//'.vtp'

      OPEN(UNIT=DES_UNIT,FILE=FNAME_VTP,STATUS='NEW',IOSTAT=ISTAT)
      IF (ISTAT /= 0) THEN
         WRITE(*,999) 
         WRITE(UNIT_LOG, 999)
         CALL MFIX_EXIT(myPE)
      ENDIF

! Dummy values to maintain format for 2D runs
      POS_Z = 0.0
      VEL_W = 0.0

! Bengin writing the *.vtp file
      WRITE(DES_UNIT,"(A)") '<?xml version="1.0"?>'
! Write the S_TIME as a comment for reference
      WRITE(DES_UNIT,"(A,ES24.16,A)") '<!-- Time =',S_TIME,'s -->'

! Write the necessary header information for a PolyData file type
      WRITE(DES_UNIT,"(A,A)") '<VTKFile type="PolyData"',&
         ' version="0.1" byte_order="LittleEndian">'
      WRITE(DES_UNIT,"(3X,A)") '<PolyData>'

! Write piece tag and identify the number of particles in the system.
      WRITE(DES_UNIT,"(6X,A,I6.6,A,A)")&
         '<Piece NumberOfPoints="',PIS,'" NumberOfVerts="0" ',&
         'NumberOfLines="0" NumberOfStrips="0" NumberOfPolys="0">'
      WRITE(DES_UNIT,"(9X,A)")&
         '<PointData Scalars="Diameter" Vectors="Velocity">'

!-----------------------
! Write the diameter data.         
      WRITE(DES_UNIT,"(12X,A)")&
         '<DataArray type="Float32" Name="Diameter" format="ascii">'
      PC = 1
      DO L = 1, MAX_PIS
         IF(PC .GT. PIS) EXIT
         IF(.NOT.PEA(L,1)) CYCLE
         WRITE (DES_UNIT,"(15X,ES12.6)") (real(2.d0*DES_RADIUS(L)))
         PC = PC + 1
      ENDDO
! Write end tag
      WRITE(DES_UNIT,"(12X,A)") '</DataArray>'


!-----------------------
! Write velocity data. Force to three dimensions. So for 2D runs, a 
! dummy value of zero is supplied as the 3rd point.
      WRITE(DES_UNIT,"(12X,A,A)") '<DataArray type="Float32" ',&
         'Name="Velocity" NumberOfComponents="3" format="ascii">'
      IF(DIMN.EQ.2) THEN
         PC = 1
         DO L = 1, MAX_PIS
            IF(PC .GT. PIS) EXIT
            IF(.NOT.PEA(L,1)) CYCLE
            WRITE (DES_UNIT,"(15X,3(ES13.6,3X))")&
               (real(DES_VEL_NEW(L,K)),K=1,DIMN),VEL_W
            PC = PC + 1
         ENDDO
      ELSE ! 3D 
         PC = 1
         DO L = 1, MAX_PIS
            IF(PC .GT. PIS) EXIT
            IF(.NOT.PEA(L,1)) CYCLE
            WRITE (DES_UNIT,"(15X,3(ES13.6,3X))")&
               (real(DES_VEL_NEW(L,K)),K=1,DIMN)
            PC = PC + 1
         ENDDO
      ENDIF
! Write end tags
      WRITE(DES_UNIT,"(12X,A,/9X,A)") '</DataArray>','</PointData>'


!-----------------------      
! Skip CellData tag, no data. (vtp format style)
      WRITE(DES_UNIT,"(9X,A)") '<CellData></CellData>'


!-----------------------      
! Write position data. Point data must be supplied in 3 dimensions. So
! for 2D runs, a dummy value of zero is supplied as the 3rd point.
      WRITE(DES_UNIT,"(9X,A)") '<Points>'
      WRITE(DES_UNIT,"(12X,A,A)") '<DataArray type="Float32" ',&
         'NAME="Position" NumberOfComponents="3" format="ascii">'
      IF(DIMN.EQ.2) THEN
         PC = 1
         DO L = 1, MAX_PIS
            IF(PC .GT. PIS) EXIT
            IF(.NOT.PEA(L,1)) CYCLE
            WRITE (DES_UNIT,"(15X,3(ES13.6,3X))")&
               (real(DES_POS_NEW(L,K)),K=1,DIMN),POS_Z 
            PC = PC + 1
         ENDDO
      ELSE
         PC = 1
         DO L = 1, MAX_PIS
            IF(PC .GT. PIS) EXIT
            IF(.NOT.PEA(L,1)) CYCLE
            WRITE (DES_UNIT,"(15X,3(ES13.6,3X))")&
               (real(DES_POS_NEW(L,K)),K=1,DIMN) 
            PC = PC + 1
         ENDDO
      ENDIF
! Write end tags
      WRITE(DES_UNIT,"(12X,A,/9X,A)")'</DataArray>','</Points>'


!-----------------------      
! Write tags for data not included (vtp format style)
      WRITE(DES_UNIT,"(9X,A,/9X,A,/9X,A,/9X,A)")'<Verts></Verts>',&
          '<Lines></Lines>','<Strips></Strips>','<Polys></Polys>'
! Write end tags
      WRITE(DES_UNIT,"(6X,A,/3X,A,/A)")&
         '</Piece>','</PolyData>','</VTKFile>'
    
      CLOSE(DES_UNIT)

! Index output file count
      VTP_FINDEX = VTP_FINDEX+1



!-----------------------      
! Construct the file that contains all the file names and solids time
! in *.vbd format. This file can be read into ParaView in place of the
! *.vtp files while providing the S_TIME data

! Determine if the "RUN_NAME"_DES_DATA.dat file exists
      F_EXISTS = .FALSE.
      FNAME_PVD = TRIM(RUN_NAME)//'_DES.pvd'

      IF(FIRST_PASS) THEN
         INQUIRE(FILE=FNAME_PVD,EXIST=F_EXISTS)

         IF (.NOT.F_EXISTS) THEN
! If the file does not exist, then create it with the necessary
! header information.
            OPEN(UNIT=PVD_UNIT,FILE=FNAME_PVD,STATUS='NEW')
            WRITE(PVD_UNIT,"(A)")'<?xml version="1.0"?>'
            WRITE(PVD_UNIT,"(A,A)")'<VTKFile type="Collection" ',&
               'version="0.1" byte_order="LittleEndian">'
            WRITE(PVD_UNIT,"(3X,A)")'<Collection>'
! write two generic lines that will be removed later
            WRITE(PVD_UNIT,*)"SCRAP LINE 1"
            WRITE(PVD_UNIT,*)"SCRAP LINE 2"
         ELSE   
! The file exists but first_pass is also true so most likely an existing
! file from an earlier/other run is present in the directory which it
! should not be if the run is NEW
            IF(RUN_TYPE .EQ. 'NEW') THEN
! Exit to prevent overwriting existing file accidently
               WRITE(*,997) FNAME_PVD
               WRITE(UNIT_LOG, 997) FNAME_PVD
               CALL MFIX_EXIT(myPE)
            ELSE
! Open the file for appending of new data (RESTART_1 Case)
               OPEN(UNIT=PVD_UNIT,FILE=FNAME_PVD,&
                  POSITION="APPEND",STATUS='OLD',IOSTAT=ISTAT)
               IF (ISTAT /= 0) THEN
                  WRITE(*,998) 
                  WRITE(UNIT_LOG, 998)
                  CALL MFIX_EXIT(myPE)
               ENDIF
            ENDIF
         ENDIF
! Identify that the files has been created and opened for next pass
         FIRST_PASS = .FALSE.         
      ELSE
         OPEN(UNIT=PVD_UNIT,FILE=FNAME_PVD,&
            POSITION="APPEND",STATUS='OLD',IOSTAT=ISTAT)
         IF (ISTAT /= 0) THEN
            WRITE(*,998) 
            WRITE(UNIT_LOG, 998)
            CALL MFIX_EXIT(myPE)
         ENDIF                 
      ENDIF

! Remove the last two lines written so that additional data can be added
      BACKSPACE(PVD_UNIT)
      BACKSPACE(PVD_UNIT)

! Force time formatting #####.######  (Forcing leading zeros)
      IF(S_TIME .LT. 1.0d0)THEN
         WRITE (S_TIME_CHAR,"(A,F7.6)")"00000",S_TIME
      ELSEIF(S_TIME .LT. 10.0d0) THEN
         WRITE (S_TIME_CHAR,"(A,F8.6)")"0000",S_TIME
      ELSEIF(S_TIME .LT. 100.0d0) THEN
         WRITE (S_TIME_CHAR,"(A,F9.6)")"000",S_TIME
      ELSEIF(S_TIME .LT. 1000.0d0) THEN
         WRITE (S_TIME_CHAR,"(A,F10.6)")"00",S_TIME
      ELSEIF(S_TIME .LT. 10000.0d0)THEN
         WRITE (S_TIME_CHAR,"(A,F11.6)")"0",S_TIME
      ELSE
         WRITE (S_TIME_CHAR,"(F12.6)"),S_TIME
      ENDIF

! Write the data to the file
      WRITE(PVD_UNIT,"(6X,A,A,A,A,A,A,A)")&
         '<DataSet timestep="',TRIM(S_TIME_CHAR),'" ',& ! simulation time
         'group="" part="0" ',&   ! necessary file data
         'file="',TRIM(FNAME_VTP),'"/>'    ! file name of vtp

! Write the closing tags
      WRITE(PVD_UNIT,"(3X,A)")'</Collection>'
      WRITE(PVD_UNIT,"(A)")'</VTKFile>'

      CLOSE(PVD_UNIT)


      RETURN


  997 FORMAT(/1X,70('*')//, ' From: WRITE_DES_VTP',/,' Message: ',&
         A, ' already exists in the run directory.',/10X,&
         'Terminating run.',/1X,70('*')/)

  998 FORMAT(/1X,70('*')//, ' From: WRITE_DES_VTP',/,' Message: ',&
         'Error opening DES vbd file. Terminating run.',/1X,70('*')/)

  999 FORMAT(/1X,70('*')//, ' From: WRITE_DES_VTP',/,' Message: ',&
         'Error opening DES vtp file. Terminating run.',/1X,70('*')/)


      END SUBROUTINE WRITE_DES_VTP 
!----------------------------------------------- 




!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!
!  Module name: WRITE_DES_TECPLOT
!  Purpose: Writing DES output in TECPLOT format
!
!
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      SUBROUTINE WRITE_DES_TECPLOT

      USE compar
      USE fldvar
      USE funits
      USE geometry
      USE indices      
      USE physprop
      USE run
      USE discretelement
      USE des_bc
      IMPLICIT NONE

!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! logical that identifies that the data file has been created
! and is already opened (initial checks/writes)
      LOGICAL, SAVE :: FIRST_PASS = .TRUE.

! logical used for testing is the data file already exists
      LOGICAL :: F_EXISTS      

! file units for 
      INTEGER, PARAMETER :: DES_DATA = 2001  ! Tecplot particle data
      INTEGER, PARAMETER :: DES_EX   = 2002  ! Tecplot extra data
      INTEGER, PARAMETER :: DES_EPS  = 2003  ! Tecplot solids fraction

! output file for basic DES variables including: position, velocity,
! radius, density, mark (flag)
      CHARACTER*50     :: FNAME_DATA

! output file for extra DES variables including:
! solids time (S_TIME), maximum neighbor count, maximum overlap
! granular energy and granular temperature
      CHARACTER*50     :: FNAME_EXTRA

! output file for axial solids volume fraction and granular temp
      CHARACTER*50     :: FNAME_EPS

! tmp character value
      CHARACTER*150    :: TMP_CHAR

! dummy indices
      INTEGER L, I, J, K, M, IJK 

! index to track accounted for particles
      INTEGER PC

! tmp variables for calculations of axial solids volume fraction, and
! granular temperature
      DOUBLE PRECISION :: AVG_EPS(JMAX2, DES_MMAX), &
         AVG_THETA(JMAX2,DES_MMAX)

!-----------------------------------------------
! Functions 
!-----------------------------------------------
!-----------------------------------------------

      INCLUDE 'function.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'ep_s2.inc'

! Set the file names for the output TECPLOT data
      FNAME_DATA = TRIM(RUN_NAME)//'_DES_DATA.dat'
      FNAME_EXTRA = TRIM(RUN_NAME)//'_DES_EXTRA.dat'
      FNAME_EPS = TRIM(RUN_NAME)//'_AVG_EPS.dat'

      IF(TECPLOT_FINDEX .GT. 0)THEN

         IF(FIRST_PASS) THEN
! If the files have not been opened since the start of the simulation
! (either because it is the first pass of a new run or a RESTART_1) 
! create the files or flag with errors and exit.
                 

! Check "RUN_NAME"_DES_DATA.dat
! the file associated with particle position, velocity, radius,
! denisty, and tagged particles.                 
!-----------------------------------------------
! Determine if the "RUN_NAME"_DES_DATA.dat file exists
            F_EXISTS = .FALSE.
            INQUIRE(FILE=FNAME_DATA,EXIST=F_EXISTS)

            IF (.NOT.F_EXISTS) THEN
! If the file does not exist, then create it with the necessary
! header information.

               OPEN(UNIT=DES_DATA,FILE=FNAME_DATA,status='new')
               IF(DIMN.EQ.3) THEN 
                  WRITE (DES_DATA, '(9(A,3X),A)')&
                  'VARIABLES = ', '"x"', '"y"', '"z"', '"vx"', '"vy"',&
                  '"vz"', '"rad"', '"den"', '"mark"'
               ELSE 
                  WRITE (DES_DATA, '(8(A,3X),A)') &
                  'VARIABLES = ', '"x"', '"y"', '"vx"', '"vy"',&
                  '"omega"','"rad"', '"den"', '"mark"'
               ENDIF

            ELSE   
! The file exists but first_pass is also true.  Thus, for the file to 
! exist it is likely an existing file from a earlier/other run (should
! not be in directory if is is a NEW run)
               IF(RUN_TYPE .EQ. 'NEW') THEN
! To prevent overwriting existing file accidently, exit if the file
! exists and this is a NEW run.                      
                  WRITE(*,3100) FNAME_DATA
                  WRITE(UNIT_LOG, 3100) FNAME_DATA
                  CALL MFIX_EXIT(myPE)
               ELSE
! Open the file for appending of new data (RESTART_1 Case)
                  OPEN(UNIT=DES_DATA,FILE=FNAME_DATA,POSITION="append")
               ENDIF
            ENDIF
!-----------------------------------------------


! Check "RUN_NAME"_DES_EXTRA.dat
! the file associated with extra simulation data: (at current time step)
!   MAX_NEIGH: Maximum number of neighbors on any particle
!   MAX_OVERLAP: Maximum overlap between any two particles
!   GRAN_ENERGY/TEMP: Global granular energy/temperature obtained 
!     by averaging over all the particles                        
!-----------------------------------------------
! Determine if the "RUN_NAME"_DES_EXTRA.dat file exists
            F_EXISTS = .FALSE.
            INQUIRE(FILE=FNAME_EXTRA,EXIST=F_EXISTS)

            IF (.NOT.F_EXISTS) THEN
! If the file does not exist, then create it with the necessary
! header information.

               OPEN(unit=DES_EX,FILE=FNAME_EXTRA, status='new')
               WRITE(DES_EX,"(5(A,3X),A)")&
               'VARIABLES = ', '"t"', '"MAX_NEIGH"', '"MAX_OVERLAP"',&
               '"GRAN_ENERGY"','"GRAN_TEMP"'

            ELSE
! To prevent overwriting existing files accidently, exit if the file
! exists and this is a NEW run.
               IF(RUN_TYPE .EQ. 'NEW') THEN
                  WRITE(*,3100) FNAME_EXTRA
                  WRITE(UNIT_LOG, 3100) FNAME_EXTRA
                  CALL MFIX_EXIT(myPE)
               ELSE
! Open the file for appending of new data (RESTART_1 Case)
                  OPEN(UNIT=DES_EX, FILE=FNAME_EXTRA, POSITION="append")
               ENDIF
            ENDIF
!-----------------------------------------------


! Check "RUN_NAME"_AVG_EPS.dat
! the file associated with axial profile in solids volume fraction and
! granular temperature
!----------------------------------------------------------------------
! Determine if the "RUN_NAME"_AVG_EPS.dat file exists
            F_EXISTS = .FALSE. 
            INQUIRE(FILE=FNAME_EPS,EXIST=F_EXISTS)

            IF (.NOT.F_EXISTS) THEN
! If the file does not exist, then create it with the necessary
! header information.

               OPEN(UNIT=DES_EPS,FILE=FNAME_EPS,status='NEW')
               WRITE(TMP_CHAR, *) ""
               DO M=1, DES_MMAX
                  WRITE(TMP_CHAR,"(A,A,I2.2,A)")&
                  TRIM(TMP_CHAR),'"`e_S_,_', M, '",'
               ENDDO
               DO M=1, DES_MMAX 
                  WRITE(TMP_CHAR,"(A,A,I2.2,A)")&
                  TRIM(TMP_CHAR), '"`Q_S_,_',M, '",'
               ENDDO
               WRITE(DES_EPS,"(A,3X,A,3X,A)")&
               "VARIABLES =", TRIM(TMP_CHAR), '"y"'

            ELSE 
! To prevent overwriting existing files accidently, exit if the file
! exists and this is a NEW run.
               IF(RUN_TYPE .EQ. 'NEW') THEN
                  WRITE(*,3100) FNAME_EPS
                  WRITE(UNIT_LOG, 3100) FNAME_EPS
                  CALL MFIX_EXIT(myPE)
               ELSE
! Open the file for appending of new data (RESTART_1 Case)
                  OPEN(UNIT=DES_EPS,FILE=FNAME_EPS,POSITION="append")
               ENDIF
            ENDIF
!-----------------------------------------------
            
! Identify that the files has been created and opened for next pass
            FIRST_PASS = .FALSE.
         ELSE 
! Open each file and mark for appending
            OPEN(UNIT=DES_DATA,  FILE=FNAME_DATA,  POSITION="append")
            OPEN(UNIT=DES_EX,    FILE=FNAME_EXTRA, POSITION="append")
            OPEN(UNIT=DES_EPS,   FILE=FNAME_EPS,   POSITION="append")
         ENDIF ! end if(FIRST_PASS)/else


! Write to the "RUN_NAME"_DES_DATA.dat file: particle position, velocity, 
! radius, denisty, and tagged particles.        
         WRITE (DES_DATA, "(A,ES24.16,A)") 'ZONE T = "', S_TIME, '"'
         PC = 1
         DO L = 1, MAX_PIS 
            IF(PC .GT. PIS) EXIT
            IF(.NOT.PEA(L,1)) CYCLE 
            IF(DIMN.EQ.3) THEN
               WRITE (DES_DATA, '(10(2X,G12.5))')&
                  (DES_POS_NEW(L, K), K = 1,DIMN), &  ! particle L position
                  (DES_VEL_NEW(L, K), K = 1,DIMN), &  ! particle L velocity
                  DES_RADIUS(L), Ro_Sol(L), &  ! radius and density
                  MARK_PART(L)  ! flagged by user
            ELSE
               WRITE (DES_DATA, '(10(2X,G12.5))')&
                  (DES_POS_NEW(L, K), K = 1,DIMN), &  ! particle L position
                  (DES_VEL_NEW(L, K), K = 1,DIMN), &  ! particle L position
                  OMEGA_NEW(L,1), & ! particle L rotation
                  DES_RADIUS(L), Ro_Sol(L), &  ! radius and density
                  MARK_PART(L)  ! flagged by user
            ENDIF
            PC = PC + 1
         ENDDO
! Close the file
         CLOSE(DES_DATA)


! Write to the "RUN_NAME"_DES_EXTRA.dat file: extra simulation data at 
! current time ste including MAX_NEIGH, MAX_OVERLAP, GRAN_ENERGY,
! GRAN_TEMP
         WRITE(DES_EX,"(G12.5,2X,I4,2X,3(G12.5,2X))")&
            S_TIME, NEIGH_MAX, OVERLAP_MAX, &
            SUM(GLOBAL_GRAN_ENERGY(1:DIMN)), &
            SUM(GLOBAL_GRAN_TEMP(1:DIMN))/(DIMN)
! Close the file
         CLOSE(DES_EX)


! Write to the "RUN_NAME"_AVG_EPS.dat file: axial profiles in solids 
! volume fraction and granular temperature of each phase
         WRITE(DES_EPS, "(A,I5,A,A,A,I5,A,ES24.16)")&
            'ZONE T = "', TECPLOT_FINDEX, '",',& ! Number of passes of routine
            'DATAPACKING = POINT,',&
            'J=',JMAX, &
            'SOLUTIONTIME=', S_TIME
          
! loop over fluid cells          
         DO J = JMIN1, JMAX1
            AVG_EPS(J,:) = ZERO 
            AVG_THETA(J,:) = ZERO
            DO K = KMIN1, KMAX1
               DO I = IMIN1, IMAX1
                  IJK = FUNIJK(I,J,K)
                  DO M = 1, DES_MMAX
                     AVG_EPS(J,M) =  AVG_EPS(J,M) + EP_S(IJK,M)
                     AVG_THETA(J,M) =  AVG_THETA(J,M) + DES_THETA(IJK,M)
                  ENDDO
               ENDDO
            ENDDO
            AVG_EPS(J,:) = AVG_EPS(J,:)/(IMAX*KMAX)
            AVG_THETA(J,:) = AVG_THETA(J,:)/(IMAX*KMAX)

            WRITE(DES_EPS,"(20(G12.5,2X))")&
               (AVG_EPS(J,M), M = 1, DES_MMAX) , &  ! average solids fraction
               (AVG_THETA(J,M), M = 1, DES_MMAX), & !
               0.5d0*(YN(J)+YN(J-1))
         ENDDO
! Close the file
         CLOSE(DES_EPS)


      ENDIF ! if ifi > 0

! Index file iteration count
      TECPLOT_FINDEX = TECPLOT_FINDEX+1
      
      RETURN

 3100 FORMAT(/1X,70('*')//, ' From: WRITE_DES_TECPLOT',/,' Message: ',&
         A, ' already exists in the run',/10X,&
         'directory. Run terminated to prevent accidental overwriting',&
         /10X,'of files.',/1X,70('*')/)

      END SUBROUTINE WRITE_DES_TECPLOT 
!-----------------------------------------------



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!
!  Module name: WRITE_DES_BEDHEIGHT
!  Purpose: Writing DES output on bed height. 

!  WARNING: This code is out-of-date and should be modified for consistency
!  with current DEM version.  Also this routine will be fairly specific
!  to a user needs and should probably be tailored as such

!     
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      SUBROUTINE WRITE_DES_BEDHEIGHT

      USE compar
      USE funits
      USE physprop
      USE run
      USE discretelement
      USE des_bc
      IMPLICIT NONE

!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! logical that identifies that the data file has been created
! and is already opened (initial checks/writes)
      LOGICAL, SAVE :: FIRST_PASS = .TRUE.
! logical used for testing is the data files already exists 
      LOGICAL :: F_EXISTS
! output file for the bed height data
      CHARACTER*50     :: FNAME_BH
! file unit for the bed height data      
      INTEGER, PARAMETER :: BH_UNIT = 2010  
! dummy index values
      INTEGER I, M
! variables for bed height calculation
      INTEGER, SAVE :: tcount = 1
      DOUBLE PRECISION :: height_avg, height_rms
      DOUBLE PRECISION, PARAMETER :: tmin = 5.d0
      DOUBLE PRECISION, DIMENSION(5000), SAVE :: bed_height_time, dt_time

!-----------------------------------------------
! Functions 
!-----------------------------------------------
!-----------------------------------------------


! after tmin start storing bed height. after enough measurements
! have been taken (i.e. tcount > 20) start to calculate a running
! average bed height and running rms bed height for solids phase 1 only
      height_avg = zero
      height_rms = zero
      
      if(time.gt.tmin) then 
         if(tcount.le.5000)  then 
            bed_height_time(tcount) = bed_height(1)
            !dt_time(tcount) = DT
            tcount = tcount + 1
            
            if(tcount.gt.20)  then
               do i = 1, tcount-1,1
                  height_avg = height_avg + bed_height_time(i)!*dt_time(i)
               enddo
               height_avg = height_avg/(tcount-1)
               do i = 1, tcount-1,1
                  height_rms = height_rms + ((bed_height_time(i)&
                       &-height_avg)**2)!*dt_time(i)
               enddo
               
               height_rms = sqrt(height_rms/(tcount-1))
            endif
         endif
      endif

      FNAME_BH = TRIM(RUN_NAME)//'_DES_BEDHEIGHT.dat'      
      IF(FIRST_PASS) THEN 
         F_EXISTS = .FALSE.
         INQUIRE(FILE=FNAME_BH,EXIST=F_EXISTS)
! If the file does not exist, then create it with the necessary
! header information.
         IF (.NOT.F_EXISTS) THEN
            OPEN(UNIT=BH_UNIT,FILE=FNAME_BH,&
               FORM="formatted",STATUS="new")
         ELSE
! To prevent overwriting existing files accidently, exit if the file
! exists and this is a NEW run.
            IF(RUN_TYPE .EQ. 'NEW') THEN
               WRITE(*,1000)
               WRITE(UNIT_LOG, 1000)
               CALL MFIX_EXIT(myPE)
            ELSE
! Open the file for appending of new data (RESTART_1 Case)
               OPEN(UNIT=BH_UNIT,FILE=FNAME_BH,POSITION="append")
            ENDIF
         ENDIF
         FIRST_PASS = .FALSE.
      ELSE
! Open the file and mark for appending              
         OPEN(UNIT=BH_UNIT,FILE=FNAME_BH,POSITION="append")
      ENDIF

      WRITE(BH_UNIT, '(10(2X,E20.12))') s_time, &
         (bed_height(M), M=1,DES_MMAX), height_avg, height_rms
! Close the file
      CLOSE(BH_UNIT)

      RETURN

 1000 FORMAT(/1X,70('*')//, ' From: WRITE_DES_BEDHEIGHT',/,&
         ' Message: bed_height.dat already exists in the run',&
         ' directory.',/10X, 'Run terminated to prevent',&
         ' accidental overwriting of files.',/1X,70('*')/)

      END SUBROUTINE WRITE_DES_BEDHEIGHT
!-----------------------------------------------



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!
!  Module name: WRITE_DES_THETA
!  Purpose: The following code writes out des_theta to a file for each
!  ijk cell in the system each time des_granular_temperature is called.
!
!
!     
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      SUBROUTINE WRITE_DES_THETA

      USE compar
      USE funits
      USE geometry
      USE indices
      USE physprop
      USE run
      USE discretelement
      USE des_bc
      IMPLICIT NONE

!-----------------------------------------------
! Local Variables
!-----------------------------------------------      
! indices
      INTEGER I, J, K, IJK
! 
      INTEGER M, LL, NP
! logical that identifies that the data file has been created
! and is already opened (initial checks/writes)
      LOGICAL, SAVE :: FIRST_PASS = .TRUE.
! logical used for testing is the data file already exists
      LOGICAL :: F_EXISTS      
! file unit for the granular temperature data
      INTEGER, PARAMETER :: GT_UNIT = 2020
! output file for the granular temperature data
      CHARACTER*50  :: FNAME_GT      
!-----------------------------------------------      

      INCLUDE 'function.inc'


      FNAME_GT = TRIM(RUN_NAME)//'_DES_THETA.dat'
      IF (FIRST_PASS) THEN
         F_EXISTS = .FALSE.               
         INQUIRE(FILE=FNAME_GT,EXIST=F_EXISTS)

         IF (.NOT.F_EXISTS) THEN
! If the file does not exist, then create it with the necessary
! header information.
            OPEN(UNIT=GT_UNIT,FILE=FNAME_GT,STATUS='NEW')
         ELSE
            IF(RUN_TYPE .EQ. 'NEW') THEN
! If the run is new and the GT file already exists replace it with a
! new file.                   
!               OPEN(UNIT=GT_UNIT,FILE=FNAME_GT,STATUS='REPLACE')
! Prevent overwriting an existing file by exiting if the file exists
! and this is a NEW run.
               WRITE(*,1001) FNAME_GT
               WRITE(UNIT_LOG,1001) FNAME_GT
               CALL MFIX_EXIT(myPE)                       
            ELSE
! Open the file for appending of new data (RESTART_1 Case)
               OPEN(UNIT=GT_UNIT, FILE=FNAME_GT, POSITION='APPEND')
            ENDIF
         ENDIF
         FIRST_PASS =  .FALSE.
      ELSE 
! Open file and mark for appending              
         OPEN(UNIT=GT_UNIT,FILE=FNAME_GT,POSITION='APPEND') 
      ENDIF   ! endif (first_pass)

      WRITE(GT_UNIT,*) ''
      WRITE(GT_UNIT,'(A6,ES24.16)') 'Time=', S_TIME
      WRITE(GT_UNIT,'(A6,2X,3(A6,2X),A8,$)') 'IJK', &
         'I', 'J', 'K', 'NP'
      DO M = 1,DES_MMAX
         WRITE(GT_UNIT,'(7X,A6,I1,$)') 'THETA_',M
      ENDDO
      WRITE(GT_UNIT,*) ''
      DO IJK = IJKSTART3, IJKEND3
         IF(FLUID_AT(IJK)) THEN
            I = I_OF(IJK)
            J = J_OF(IJK)
            K = K_OF(IJK)
            NP = PINC(IJK)
            WRITE(GT_UNIT,'(I6,2X,3(I6,2X),I8,(2X,ES15.5))') &
               IJK, I, J, K, NP, (DES_THETA(IJK,M), M = 1,DES_MMAX)
         ENDIF
      ENDDO

! Close the file
      CLOSE(GT_UNIT)

      RETURN

 1001 FORMAT(/1X,70('*')//, ' From: WRITE_DES_THETA',/,&
         ' Message: ', A, ' already exists in the run',/10X,&
         'directory. Run terminated to prevent accidental overwriting',&
         /10X,'of files.',/1X,70('*')/)

      END SUBROUTINE WRITE_DES_THETA
!-----------------------------------------------
