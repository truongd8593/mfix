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

      USE param
      USE param1
      USE parallel
      USE fldvar
      USE discretelement
      USE run
      USE geometry
      USE physprop
      USE sendrecv
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
! NOTE: If the system starts with zero particles, ParaView may have
! trouble showing the results. To view the results in the current
! version of ParaView, Version 3.6.1:
!   i - load the *.vtp files
!  ii - filter with glyph (Filters > Common > Glyph)
!       a - change glyph to sphere
!       b - change scale mode to scalar
!       c - check the "Edit" Set Scale Factor Box
!       d - change the value to 1.0
!     
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      SUBROUTINE WRITE_DES_VTP

      USE param
      USE param1
      USE parallel
      USE fldvar
      USE discretelement
      USE run
      USE geometry
      USE physprop
      USE sendrecv
      USE des_bc
      IMPLICIT NONE

!-----------------------------------------------
! Local Variables
!-----------------------------------------------
      INTEGER, PARAMETER :: DES_UNIT = 2000  ! ParaView *.vtp data

! Logicals used for testing is the data files already exist 
      LOGICAL :: F_EXISTS

! Output file for bed height calculations
      CHARACTER*50     :: FNAME_BH
      INTEGER, PARAMETER :: BH_UNIT = 2010  ! bed height data

! This logical identifies that a data file has been created and
! already opened.
      LOGICAL, SAVE :: FIRST_PASS = .TRUE.

! Index used when writing DES_*.vtp file
      CHARACTER*5 F_INDEX

! index to track accounted for particles
      INTEGER PC

! Dummy index values
      INTEGER L, I, J, K, M, IJK

! dummy values to maintain format for dimn=2
      REAL POS_Z, VEL_W 
 
! variables for bed height calculation
      INTEGER, SAVE :: tcount = 1
      DOUBLE PRECISION :: height_avg, height_rms
      DOUBLE PRECISION, PARAMETER :: tmin = 5.d0
      DOUBLE PRECISION, DIMENSION(5000), SAVE :: bed_height_time, dt_time

!-----------------------------------------------
! Functions 
!-----------------------------------------------
!-----------------------------------------------

      INCLUDE 'function.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'ep_s2.inc'

! Convert the index VTP_FINDEX from an integer to a string and force
! leading zeros
      WRITE (F_INDEX,"(I5.5)") VTP_FINDEX
      OPEN(UNIT=DES_UNIT,FILE=TRIM(RUN_NAME)//'_DES_'//F_INDEX//'.vtp',&
         STATUS='NEW',ERR=999)

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
      END DO
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

      RETURN


! Following code is not currently accessible to this subroutine (after
! return statement). This code is also out-of-date and should be
! modified for consistency with current DEM version
!######################################################################

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

      FNAME_BH = 'bed_height.dat'      
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
         (bed_height(j), j=1,MMAX), height_avg, height_rms
! Close the file and keep
      CLOSE(BH_UNIT, STATUS="KEEP")

! End code not currently accessible to subroutine      
!######################################################################




  999 WRITE(*,"(/1X,70('*'),//,A,/,A,/1X,70('*'))")&
         ' From: WRITE_DES_VTP ',&
         ' Message: Error opening DES vtp file. Terminating run.'
      CALL MFIX_EXIT(myPE)

 1000 FORMAT(/1X,70('*')//, ' From: WRITE_DES_VTP',/,&
         ' Message: bed_height.dat already exists in the run',&
         ' directory.',/10X, 'Run terminated to prevent',&
         ' accidental overwriting of files.',/1X,70('*')/)



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

      USE param
      USE param1
      USE parallel
      USE fldvar
      USE discretelement
      USE run
      USE geometry
      USE physprop
      USE sendrecv
      USE des_bc
      IMPLICIT NONE

!-----------------------------------------------
! Local Variables
!-----------------------------------------------
      INTEGER, PARAMETER :: DES_DATA = 2001  ! Tecplot particle data
      INTEGER, PARAMETER :: DES_EX   = 2002  ! Tecplot extra data
      INTEGER, PARAMETER :: DES_EPS  = 2003  ! Tecplot solids fraction

! This logical identifies that the data files have been created and
! already opened.
      LOGICAL, SAVE :: FIRST_PASS = .TRUE.

! Logicals used for testing is the data files already exist and/or
! if the file is currently open.
      LOGICAL :: F_EXISTS

! Output file for basic DES variables including: position, velocity,
! radius, density, mark (flag)
      CHARACTER*50     :: FNAME_DATA

! Output file for extra DES variables including:
! solids time (S_TIME), maximum neighbor count, maximum overlap
! granular energy and granular temperature
      CHARACTER*50     :: FNAME_EXTRA

! Output file for axial solids volume fraction and granular temp
      CHARACTER*50     :: FNAME_EPS

! Tmp character value
      CHARACTER*150    :: TMP_CHAR

! Dummy indices
      INTEGER L, I, J, K, M, IJK 

! Index to track accounted for particles
      INTEGER PC

! Tmp variables for calculationsof axial solids volume fraction, and
! granular temperature
      DOUBLE PRECISION :: AVG_EPS(JMAX2, MMAX), AVG_THETA(JMAX2, MMAX)

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
                 

! Check"RUN_NAME"_DES_DATA.dat
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


! Check"RUN_NAME"_DES_EXTRA.dat
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


! Check"RUN_NAME"_AVG_EPS.dat
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
               DO M=1, MMAX
                  WRITE(TMP_CHAR,"(A,A,I2.2,A)")&
                  TRIM(TMP_CHAR),'"`e_S_,_', M, '",'
               ENDDO
               DO M=1, MMAX 
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
               END IF
            ENDIF
!-----------------------------------------------
            
! Identify that the files has been created and opened for next pass
            FIRST_PASS = .FALSE.
         ELSE 
! Open each file and mark for appending
            OPEN(UNIT= DES_DATA, FILE=FNAME_DATA,  POSITION="append")
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
! Close the file and keep
         CLOSE(DES_DATA, STATUS = "keep")


! Write to the "RUN_NAME"_DES_EXTRA.dat file: extra simulation data at 
! current time ste including MAX_NEIGH, MAX_OVERLAP, GRAN_ENERGY,
! GRAN_TEMP
         WRITE(DES_EX,"(G12.5,2X,I4,2X,3(G12.5,2X))")&
            S_TIME, NEIGH_MAX, OVERLAP_MAX, &
            SUM(GLOBAL_GRAN_ENERGY(1:DIMN)), &
            SUM(GLOBAL_GRAN_TEMP(1:DIMN))/(DIMN)
! Close the file and keep
         CLOSE(DES_EX, STATUS = "keep")


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
                  DO M = 1, MMAX
                     AVG_EPS(J,M) =  AVG_EPS(J,M) + EP_S(IJK,M)
                     AVG_THETA(J,M) =  AVG_THETA(J,M) + DES_THETA(IJK,M)
                  ENDDO
               ENDDO
            ENDDO
            AVG_EPS(J,:) = AVG_EPS(J,:)/(IMAX*KMAX)
            AVG_THETA(J,:) = AVG_THETA(J,:)/(IMAX*KMAX)

            WRITE(DES_EPS,"(20(G12.5,2X))")&
               (AVG_EPS(J,M), M = 1, MMAX) , &  ! average solids fraction
               (AVG_THETA(J,M), M = 1, MMAX), & !
               0.5d0*(YN(J)+YN(J-1))
         ENDDO
! Close the file and keep
         CLOSE(DES_EPS, STATUS = "keep")


      ENDIF ! if ifi > 0

! Index file iteration count
! Note: This value is also indexed in the subroutine WRITE_DES_DATA. 
! Currently these routines can not be called at the same time so this
! should not create any problems.
      TECPLOT_FINDEX = TECPLOT_FINDEX+1
      
      RETURN

 3100 FORMAT(/1X,70('*')//, ' From: WRITE_DES_TECPLOT',/,' Message: ',&
         A, ' already exists in the run',/10X,&
         'directory. Run terminated to prevent accidental overwriting',&
         /10X,'of files.',/1X,70('*')/)

      END SUBROUTINE WRITE_DES_TECPLOT 
!-----------------------------------------------

