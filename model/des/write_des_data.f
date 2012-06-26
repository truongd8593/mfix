!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: WRITE_DES_DATA                                          !
!  Purpose: Writing DES output in Paraview format                      !
!                                                                      !
!                                                                      !
!  Author: Jay Boyalakuntla                           Date: 26-Jul-06  !
!  Reviewer: Sreekanth Pannala                        Date: 31-Oct-06  !
!                                                                      !
!  Reviewer: J. Musser                                Date: 20-Apr-10  !
!  Comments: Split original subroutine into one for ParaView *.vtp     !
!  files, and a second for TECPLOT files *.dat.                        !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE WRITE_DES_DATA

!-----------------------------------------------
! Modules
!-----------------------------------------------
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

! Invoke at own risk      
      IF (.FALSE.) CALL WRITE_DES_THETA
      IF (.FALSE.) CALL WRITE_DES_BEDHEIGHT

      !if(vtp_findex.eq.2)  stop
      RETURN
      END SUBROUTINE WRITE_DES_DATA
!----------------------------------------------- 



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: WRITE_DES_VTP                                           !
!  Purpose: Writing DES output in Paraview format                      !
!                                                                      !
!                                                                      !
!  Reviewer: Rahul Garg                               Dare: 01-Aug-07  !
!  Comments: Added one more output file containing averaged bed height !
!                                                                      !
!  Revision : For parallel runs added cumulative and parallel IO       !
!  Author   : Pradeep G.                                               !
!                                                                      !
!  NOTE: If the system starts with zero particles, ParaView may have   !
!  trouble showing the results. To view the results in the current     !
!  version of ParaView, Version 3.6.1:                                 !
!    i - load the *.vtp files                                          !
!   ii - filter with glyph (Filters > Common > Glyph)                  !
!        a - change glyph to sphere                                    !
!        b - change scale mode to scalar                               !
!        c - check the "Edit" Set Scale Factor Box                     !
!        d - change the value to 1.0                                   !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE WRITE_DES_VTP

!-----------------------------------------------
! Modules
!-----------------------------------------------
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
      use mpi_utility
      use compar 
      use desmpi
      use cdist
      Use des_thermo
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! logical that identifies that the data file has been created
! and is already opened (initial checks/writes)
      LOGICAL, SAVE :: FIRST_PASS = .TRUE.

! logical used for testing is the data file already exists
      LOGICAL :: F_EXISTS      

! file unit for ParaView *.vtp data      
!      INTEGER, PARAMETER :: DES_UNIT = 2000
      INTEGER :: DES_UNIT = 2000

! file unit for ParaView *.vbd data      
      INTEGER, PARAMETER :: PVD_UNIT = 2050
      
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


! Variables related to gather 
      integer llocalcnt,lglocnt,lgathercnts(0:numpes-1),lproc
      real,dimension(:,:), allocatable :: ltemp_array

! check whether an error occurs in opening a file
      INTEGER ISTAT      
!-----------------------------------------------
! Include statement functions 
!-----------------------------------------------
      INCLUDE 'function.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'ep_s2.inc'
!-----------------------------------------------

! set the file name and unit number and open file 
      if (bdist_io) then 
         write(fname_vtp,'(A,"_DES",I4.4,"_",I5.5,".vtp")') trim(run_name),vtp_findex,mype
         CALL OPEN_VTP
      else 
         if(mype.eq.pe_io) then 
            write(fname_vtp,'(A,"_DES_",I5.5,".vtp")') trim(run_name),vtp_findex 
            CALL OPEN_VTP
         end if  
      end if 


! Dummy values to maintain format for 2D runs
      POS_Z = 0.0
      VEL_W = 0.0

!Write VTP file either in distributed mode or single file  
      if (bdist_io) then  
         write(des_unit,"(a)") '<?xml version="1.0"?>'
! Write the S_TIME as a comment for reference
         write(des_unit,"(a,es24.16,a)") '<!-- Time =',s_time,'s -->'

! Write the necessary header information for a PolyData file type
         write(des_unit,"(a,a)") '<VTKFile type="PolyData"',&
            ' version="0.1" byte_order="LittleEndian">'
         write(des_unit,"(3x,a)") '<PolyData>'

! Write Piece tag and identify the number of particles in the system.
         write(des_unit,"(6x,a,i10.10,a,a)")&
            '<Piece NumberOfPoints="',pip-ighost_cnt,&
            '" NumberOfVerts="0" ',&
            'NumberOfLines="0" NumberOfStrips="0" NumberOfPolys="0">'
         write(des_unit,"(9x,a)")&
            '<PointData Scalars="Diameter" Vectors="Velocity">'

! Write the diameter data.         
         write(des_unit,"(12x,a)")&
            '<DataArray type="Float32" Name="Diameter" format="ascii">'
         pc = 1
         do l = 1,max_pip
            if(pc.gt.pip) exit
            if(.not.pea(l,1)) cycle 
            pc = pc+1
            if(pea(l,4)) cycle 
            write (des_unit,"(15x,es12.6)") (real(2.d0*des_radius(l)))
         end do
         write(des_unit,"(12x,a)") '</DataArray>'

! Write the mag. of cohesive force scaled with the weight of the particle.
         if(use_cohesion) then
           write(des_unit,"(12x,a)")&
        '<DataArray type="Float32" Name="cohesiveForce" format="ascii">'
           pc = 1
           do l = 1,max_pip
              if(pc.gt.pip) exit
              if(.not.pea(l,1)) cycle 
              pc = pc+1
              if(pea(l,4)) cycle 
              write (des_unit,"(15x,es12.6)") (PostCohesive(l))
           end do
           write(des_unit,"(12x,a)") '</DataArray>'        
         endif

! Write velocity data. Force to three dimensions. So for 2D runs, a 
! dummy value of zero is supplied as the 3rd point.
         write(des_unit,"(12x,a,a)") '<DataArray type="Float32" ',&
            'Name="Velocity" NumberOfComponents="3" format="ascii">'
         if(dimn.eq.2) then
            pc = 1
            do l = 1,max_pip
               if(pc.gt.pip) exit
               if(.not.pea(l,1)) cycle 
               pc = pc+1
               if(pea(l,4)) cycle 
               write (des_unit,"(15x,3(es13.6,3x))")&
                  (real(des_vel_new(l,k)),k=1,dimn),vel_w
            enddo
         else ! 3d 
            pc = 1
            do l = 1,max_pip
               if(pc.gt.pip) exit
               if(.not.pea(l,1)) cycle 
               pc = pc+1
               if(pea(l,4)) cycle 
               write (des_unit,"(15x,3(es13.6,3x))")&
                  (real(des_vel_new(l,k)),k=1,dimn)
            enddo
         endif
         write(des_unit,"(12x,a,/9x,a)") '</DataArray>','</PointData>'

! Skip CellData tag, no data. (vtp format style)
         write(des_unit,"(9x,a)") '<CellData></CellData>'

! Write position data. Point data must be supplied in 3 dimensions. So
! for 2D runs, a dummy value of zero is supplied as the 3rd point.
         write(des_unit,"(9x,a)") '<Points>'
         write(des_unit,"(12x,a,a)") '<DataArray type="Float32" ',&
            'Name="Position" NumberOfComponents="3" format="ascii">'
         if(dimn.eq.2) then
            pc = 1
            do l = 1,max_pip
               if(pc.gt.pip) exit
               if(.not.pea(l,1)) cycle 
               pc = pc+1
               if(pea(l,4)) cycle 
               write (des_unit,"(15x,3(es13.6,3x))")&
                  (real(des_pos_new(l,k)),k=1,dimn),pos_z 
            enddo
         else
            pc = 1
            do l = 1,max_pip
               if(pc.gt.pip) exit
               if(.not.pea(l,1)) cycle 
               pc = pc+1
               if(pea(l,4)) cycle 
               write (des_unit,"(15x,3(es13.6,3x))")&
                  (real(des_pos_new(l,k)),k=1,dimn) 
            enddo
         endif
         write(des_unit,"(12x,a,/9x,a)")'</DataArray>','</Points>'

! Write tags for data not included (vtp format style)
         write(des_unit,"(9x,a,/9x,a,/9x,a,/9x,a)")'<Verts></Verts>',&
             '<Lines></Lines>','<Strips></Strips>','<Polys></Polys>'
         write(des_unit,"(6x,a,/3x,a,/a)")&
            '</Piece>','</PolyData>','</VTKFile>'
      else 
! write a single file at PE_IO and set parameters for gather 
         lglocnt = 10
         llocalcnt = pip - ighost_cnt
         call global_sum(llocalcnt,lglocnt) 
         allocate (dprocbuf(llocalcnt),drootbuf(lglocnt),iprocbuf(llocalcnt),irootbuf(lglocnt))
         allocate (ltemp_array(lglocnt,3)) 
         igath_sendcnt = llocalcnt 
         lgathercnts = 0
         lgathercnts(mype) = llocalcnt
         call global_sum(lgathercnts,igathercnts)
         idispls(0) = 0 
         do lproc = 1,numpes-1 
            idispls(lproc) = idispls(lproc-1) + igathercnts(lproc-1)  
         end do 

! write diameter 
         call des_gather(des_radius)
         if (mype.eq.pe_io) then 
            write(des_unit,"(a)") '<?xml version="1.0"?>'
            write(des_unit,"(a,es24.16,a)") '<!-- time =',s_time,'s -->'
            write(des_unit,"(a,a)") '<VTKFile type="PolyData"',&
               ' version="0.1" byte_order="LittleEndian">'
            write(des_unit,"(3x,a)") '<PolyData>'
            write(des_unit,"(6x,a,i10.10,a,a)")&
               '<Piece NumberOfPoints="',lglocnt,'" NumberOfVerts="0" ',&
               'NumberOfLines="0" NumberOfStrips="0" NumberOfPolys="0">'
            write(des_unit,"(9x,a)")&
               '<PointData Scalars="Diameter" Vectors="Velocity">'
            write(des_unit,"(12x,a)")&
               '<DataArray type="Float32" Name="Diameter" format="ascii">'
            write (des_unit,"(15x,es12.6)") (real(2.d0*drootbuf(l)),l=1,lglocnt)
            write(des_unit,"(12x,a)") '</DataArray>'
         endif         

! write scaled mag. of cohesive force 
         if(use_cohesion) then
           call des_gather(PostCohesive)
           if (mype.eq.pe_io) then 
              write(des_unit,"(12x,a)")&
                 '<DataArray type="Float32" Name="Fcohesive" format="ascii">'
              write (des_unit,"(15x,es12.6)") ((drootbuf(l)),l=1,lglocnt)
              write(des_unit,"(12x,a)") '</DataArray>'
           endif    
         endif    ! for cohesion

!-----------------------
! Write the temperature data.
         IF((MYPE .EQ. PE_IO) .AND. DES_ENERGY_EQ) THEN
            WRITE(DES_UNIT,"(12X,A)")&
               '<DataArray type="Float32" Name="Temperature" format="ascii">'
            PC = 1
            DO L = 1, MAX_PIP
               IF(PC .GT. PIP) EXIT
               IF(.NOT.PEA(L,1)) CYCLE
               WRITE (DES_UNIT,"(15X,ES12.6)") (real(DES_T_s_NEW(L)))
               PC = PC + 1
            END DO
! Write end tag
            WRITE(DES_UNIT,"(12X,A)") '</DataArray>'
        ENDIF

! Write velocity data.
         ltemp_array = 0.0 
         do k = 1,dimn
            call des_gather(des_vel_new(:,k))
            ltemp_array(:,k) = drootbuf(:)
         end do
         if (mype.eq.pe_io) then 
            write(des_unit,"(12x,a,a)") '<DataArray type="Float32" ',&
               'Name="Velocity" NumberOfComponents="3" format="ascii">'
            write (des_unit,"(15x,3(es13.6,3x))")&
               ((ltemp_array(l,k),k=1,3),l=1,lglocnt)
            write(des_unit,"(12x,a,/9x,a)") '</DataArray>','</PointData>'
! skip cell data
            write(des_unit,"(9x,a)") '<CellData></CellData>'
         end if 


!write position data 
         ltemp_array = 0.0 
         do k = 1,dimn
            call des_gather(des_pos_new(:,k))
            ltemp_array(:,k) = drootbuf(:)
         end do
         if (mype.eq.pe_io) then 
            write(des_unit,"(9x,a)") '<Points>'
            write(des_unit,"(12x,a,a)") '<DataArray type="Float32" ',&
               'Name="Position" NumberOfComponents="3" format="ascii">'
            write (des_unit,"(15x,3(es13.6,3x))")&
               ((ltemp_array(l,k),k=1,3),l=1,lglocnt)
            write(des_unit,"(12x,a,/9x,a)")'</DataArray>','</Points>'
! Write tags for data not included (vtp format style)
            write(des_unit,"(9x,a,/9x,a,/9x,a,/9x,a)")'<Verts></Verts>',&
                '<Lines></Lines>','<Strips></Strips>','<Polys></Polys>'
            write(des_unit,"(6x,a,/3x,a,/a)")&
               '</Piece>','</PolyData>','</VTKFile>'
         endif
         deallocate (dprocbuf,drootbuf,iprocbuf,irootbuf,ltemp_array)
      end if  

      VTP_FINDEX=VTP_FINDEX+1
 
      if (bdist_io .or. (myPE .eq.pe_IO)) close(des_unit)


      IF(.not.bdist_io.and.myPE.eq.pe_IO) THEN 
!-----------------------      
! Construct the file that contains all the file names and solids time
! in *.vbd format. This file can be read into ParaView in place of the
! *.vtp files while providing the S_TIME data

! Obtain the file name and open the pvd file
         FNAME_PVD = TRIM(RUN_NAME)//'_DES.pvd'
         CALL OPEN_PVD

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
         'group="" part="0" ',& ! necessary file data
         'file="',TRIM(FNAME_VTP),'"/>' ! file name of vtp

! Write the closing tags
         WRITE(PVD_UNIT,"(3X,A)")'</Collection>'
         WRITE(PVD_UNIT,"(A)")'</VTKFile>'
         
         CLOSE(PVD_UNIT)
      ENDIF
      
      
      RETURN

      CONTAINS
!......................................................................!
! SUBROUTINE: OPEN_VTP                                                 !
!                                                                      !
! Purpose: This routine opens a vtp file.                              !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE OPEN_VTP

!-----------------------------------------------
! Local Variables
!----------------------------------------------- 
! logical used for testing is the data file already exists
      LOGICAL :: EXISTS_VTP
! status of the vtp file to be written
      CHARACTER*8 :: STATUS_VTP
!-----------------------------------------------

! Check to see if the file already exists.
      INQUIRE(FILE=FNAME_VTP,EXIST=EXISTS_VTP)
! The given file should not exist if the run type is NEW.
      IF(RUN_TYPE == 'NEW' .AND. EXISTS_VTP)THEN
         !WRITE(*,1000) FNAME_VTP
         IF(DMP_LOG) WRITE(UNIT_LOG, 1000) FNAME_VTP
! A RESTART_1 case may need to over-write vtp files created during the
! previous run.
      ELSEIF(RUN_TYPE == 'RESTART_1' .AND. EXISTS_VTP)THEN
         STATUS_VTP = 'REPLACE'
! Set the status to NEW so that any problems in opening the file
! will be accurately reported.
      ELSEIF(RUN_TYPE == 'NEW' .AND. .NOT.EXISTS_VTP)THEN
         STATUS_VTP = 'NEW'
      ELSEIF(RUN_TYPE == 'RESTART_1' .AND. .NOT.EXISTS_VTP)THEN
         STATUS_VTP = 'NEW'
      ENDIF

      OPEN(UNIT=DES_UNIT,FILE=FNAME_VTP,STATUS=STATUS_VTP,IOSTAT=ISTAT)
      IF (ISTAT /= 0) THEN
         IF(DMP_LOG) WRITE(UNIT_LOG, 1001) FNAME_VTP, DES_UNIT, STATUS_VTP
         IF(PRINT_DES_SCREEN ) WRITE(*, 1001) FNAME_VTP, DES_UNIT, STATUS_VTP
         CALL MFIX_EXIT(myPE)
      ENDIF


 1000 FORMAT(/1X,70('*')/,' From: OPEN_VTP',/,' Message: ',            &
         'The following vtp file was found in the run directory ',     &
         'for a',/10X,'run classified as NEW. Please correct.',/10X,   &
         'File name: ',A,/1X,70('*')/)

 1001 FORMAT(/1X,70('*')//, ' From: OPEN_VTP',/,' Message: ',          &
         'Error opening DES vtp file. Terminating run.',/10X,          &
         'File name:  ',A,/10X,                                         &
         'DES_UNIT :  ',i4, /10X,                                       &
         'STATUS_VTP: ', A, &
         /1X,70('*')/)

      END SUBROUTINE OPEN_VTP


!......................................................................!
! SUBROUTINE: OPEN_PVD                                                 !
!                                                                      !
! Purpose: This routine opens the pvd file.                            !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE OPEN_PVD

!-----------------------------------------------
! Local Variables
!-----------------------------------------------       
! Index position of desired character
      INTEGER IDX_f, IDX_b
! logical used for testing is the data file already exists
      LOGICAL :: EXISTS_PVD
! status of the vtp file to be written
      CHARACTER*8 :: STATUS_PVD
! Generic input limited to 256 characters
      CHARACTER*256 INPUT
!----------------------------------------------- 

! Check to see if the file already exists.
      INQUIRE(FILE=FNAME_PVD,EXIST=EXISTS_PVD)

      IF(FIRST_PASS) THEN
         IF(RUN_TYPE == 'NEW')THEN
! If the file does not exist, then create it with the necessary
! header information.
            IF (.NOT.EXISTS_PVD) THEN
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
! Exit to prevent overwriting existing file accidently
               WRITE(*,1002) TRIM(FNAME_PVD)
               WRITE(UNIT_LOG, 1002) TRIM(FNAME_PVD)
               STOP
            ENDIF
! This is the first pass of a restart run. Extra care is needed to make
! sure that the pvd file is ready to accept new data.
         ELSE ! a restart run
            IF (.NOT.EXISTS_PVD) THEN
! For a restart run, there should be a pvd file already in the run
! directory. If there is not a pvd file, notifiy the user and exit.
               WRITE(*,1003) TRIM(FNAME_PVD)
               WRITE(UNIT_LOG, 1003) TRIM(FNAME_PVD)
               STOP
! Modify the pvd file.
            ELSE ! a pvd file does exist
! Open the file at the beginning.
               OPEN(UNIT=PVD_UNIT,FILE=FNAME_PVD,&
                  POSITION="REWIND",STATUS='OLD',IOSTAT=ISTAT)
               IF (ISTAT /= 0) THEN
                  WRITE(*,1004) 
                  WRITE(UNIT_LOG, 1004)
                  STOP
               ENDIF

! Loop over the entries in the PVD file, looking for a match to the
! file that is being written. If no match is found, the data will be
! appended to the end of the pvd file, otherwise, the old data will
! be over-written.
               DO
! Read in the entires of the PVD file.
                  READ(PVD_UNIT,"(A)",IOSTAT=ISTAT)INPUT
                  IF(ISTAT > 0)THEN
                     WRITE(*,1005)
                     STOP
                  ELSEIF(ISTAT<0)THEN
! The end of the pvd file has been reached without finding an entry 
! matching the current record. Exit the loop.
                     BACKSPACE(PVD_UNIT)
                     EXIT
                  ENDIF
! Find the first instances of file=" and "/> in the read data.
                  IDX_f = INDEX(INPUT,'file="')
                  IDX_b = INDEX(INPUT,'"/>')
! Skip rows that do not contain file data
                  IF(IDX_f == 0 .AND. IDX_b == 0) CYCLE
! Truncate the file name from the read data
                  WRITE (INPUT,"(A)") INPUT(IDX_f+6:IDX_b-1)
! If the file name matches the current VTP record, return to the calling
! routine to over-write this record.
                  IF(TRIM(FNAME_VTP) == TRIM(INPUT)) THEN
                     BACKSPACE(PVD_UNIT)
                     RETURN
                  ENDIF
               ENDDO

            ENDIF ! a pvd file does not/does exist
         ENDIF ! run_type new or restart

! Identify that the files has been created and opened for next pass
         FIRST_PASS = .FALSE.

      ELSE ! not FIRST_PASS
         OPEN(UNIT=PVD_UNIT,FILE=FNAME_PVD,&
            POSITION="APPEND",STATUS='OLD',IOSTAT=ISTAT)
         IF (ISTAT /= 0) THEN
            WRITE(*,1004) 
            WRITE(UNIT_LOG, 1004)
            STOP
         ENDIF                 
      ENDIF

! Remove the last two lines written so that additional data can be added
      BACKSPACE(PVD_UNIT)
      BACKSPACE(PVD_UNIT)

! Return to the calling routine
      RETURN

 1002 FORMAT(/1X,70('*')/,' From: OPEN_PVD',/,' Message: ',            &
         A,' already exists in the run directory.',/10X,               &
         'Terminating run.',/1X,70('*')/)

 1003 FORMAT(/1X,70('*')/,' From: OPEN_PVD',/,' Message: ',            &
         A,' is missing from the  the run directory.',/10X,            &
         'Terminating run.',/1X,70('*')/)

 1004 FORMAT(/1X,70('*')/,' From: OPEN_PVD',/,' Message: ',            &
         'Error opening DES pvd file. Terminating run.',/1X,70('*')/)

 1005 FORMAT(/1X,70('*')/,' From: OPEN_PVD',/,' Message: ',            &
         'Error reading DES pvd file. Terminating run.',/1X,70('*')/)

      END SUBROUTINE OPEN_PVD


      END SUBROUTINE WRITE_DES_VTP 
!----------------------------------------------- 



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!
!  Module name: WRITE_DES_TECPLOT
!  Purpose: Writing DES output in TECPLOT format
!
!  Revision: For parallel runs added distributed and single IO
!  Comment: In earlier version the time instances are keep appended to 
!           the tecplot file. This will make tecplot file so large for 
!           large simulations. Hence seperate files are written for 
!           each instances
!  Author : Pradeep G.
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
      use compar 
      use cdist 
      use desmpi
      use mpi_utility
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
! Pradeep remove parameter from following variables
      INTEGER:: DES_DATA,DES_EX,DES_EPS

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
      DOUBLE PRECISION :: AVG_EPS(JMAX2, DES_MMAX),&
                          AVG_THETA(JMAX2, DES_MMAX)

! Variables related to gathering info at PE_IO 
      integer llocalcnt,lglocnt,lgathercnts(0:numpes-1),lproc,ltotvar,lcount
      real,dimension(:,:), allocatable :: ltemp_array
!-----------------------------------------------
! Functions 
!-----------------------------------------------
!-----------------------------------------------

      INCLUDE 'function.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'ep_s2.inc'

! set the total variable based on dimension 
      if (dimn.eq.3) then  
        ltotvar = 9 
      else
        ltotvar = 8 
      end if 

! set the file name and unit number and open file 
      des_data = 2000
      des_ex = 2100
      des_eps = 2200
      if (bdist_io) then 
         write(fname_data,'(A,"_DES_DATA",I4.4,"_",I4.4,".dat")') trim(run_name),tecplot_findex,mype 
         write(fname_extra,'(A,"_DES_EXTRA",I4.4,"_",I4.4,".dat")') trim(run_name),tecplot_findex,mype 
         write(fname_eps,'(A,"_DES_EPS",I4.4,"_",I4.4,".dat")') trim(run_name),tecplot_findex,mype
         open(unit=des_data,file=fname_data,status='new',err=999)
         open(unit=des_ex,file=fname_extra,status='new',err=999)
         open(unit=des_eps,file=fname_eps,status='new',err=999)
      else 
         if(mype.eq.pe_io) then 
            write(fname_data,'(A,"_DES_DATA_",I4.4,".dat")') trim(run_name),tecplot_findex
            write(fname_extra,'(A,"_DES_EXTRA_",I4.4,".dat")') trim(run_name),tecplot_findex
            write(fname_eps,'(A,"_DES_EPS_",I4.4,".dat")') trim(run_name),tecplot_findex
            open(unit=des_data,file=fname_data,status='new',err=999)
            open(unit=des_ex,file=fname_extra,status='new',err=999)
            open(unit=des_eps,file=fname_eps,status='new',err=999)
         end if  
      end if 
      tecplot_findex = tecplot_findex + 1
       
! write header 
      if (bdist_io .or. mype .eq. pe_io) then 
         if(dimn.eq.3) then 
            write (des_data,'(A)') & 
            'variables = "x" "y" "z" "vx" "vy" "vz" "rad" "den" "mark"'
         else 
            write (des_data, '(A)') &
            'variables = "x" "y" "vx" "vy" "omega" "rad" "den" "mark"'
         endif
         write (des_data, "(A,F15.7,A)")'zone T ="', s_time, '"'
      end if 

! write data in ditributed mode or as a single file  
      if(bdist_io) then 
         pc = 1
         do l = 1,max_pip
            if(pc.gt.pip) exit
            if(.not.pea(l,1)) cycle 
            pc = pc+1
            if(pea(l,4)) cycle 
            if(dimn.eq.3) then
               write (des_data, '(8(2x,es12.5),I5)')&
                  (des_pos_new(l,k),k=1,dimn),(des_vel_new(l,k),k=1,dimn), &
                   des_radius(l),ro_sol(l),mark_part(l) 
            else
               write (des_data, '(7(2x,es12.5),I5)')&
                  (des_pos_new(l,k),k=1,dimn), (des_vel_new(l,k),k=1,dimn), & 
                  omega_new(l,1), des_radius(l), ro_sol(l),mark_part(l) 
            endif
        end do 

      else ! if bdist_io 
! set parameters required for gathering info at PEIO and write as single file   
         lglocnt = 10
         llocalcnt = pip - ighost_cnt
         call global_sum(llocalcnt,lglocnt) 
         allocate (dprocbuf(llocalcnt),drootbuf(lglocnt),iprocbuf(llocalcnt),irootbuf(lglocnt))
         allocate (ltemp_array(lglocnt,ltotvar)) 
         igath_sendcnt = llocalcnt 
         lgathercnts = 0
         lgathercnts(mype) = llocalcnt
         call global_sum(lgathercnts,igathercnts)
         idispls(0) = 0 
         do lproc = 1,numpes-1 
            idispls(lproc) = idispls(lproc-1) + igathercnts(lproc-1)  
         end do 

! gather information from all processor 
         lcount = 1
         do k = 1,dimn 
            call des_gather(des_pos_new(:,k))
            ltemp_array(:,lcount) = drootbuf(:); lcount=lcount+1
         end do  
         do k = 1,dimn 
            call des_gather(des_vel_new(:,k))
            ltemp_array(:,lcount) = drootbuf(:); lcount=lcount+1
         end do  
         if(dimn.eq.2) then 
            call des_gather(omega_new(:,1))
            ltemp_array(:,lcount) = drootbuf(:); lcount=lcount+1
         end if   
         call des_gather(des_radius)
         ltemp_array(:,lcount) = drootbuf(:); lcount=lcount+1
         call des_gather(ro_sol)
         ltemp_array(:,lcount) = drootbuf(:); lcount=lcount+1
         call des_gather(mark_part)
         ltemp_array(:,lcount) = irootbuf(:); lcount=lcount+1

! write the data into file 
         if (mype.eq.pe_io) then 
            if(dimn.eq.3) then
               do l =1,lglocnt
                  write (des_data,'(8(2x,es12.5),I5)') (ltemp_array(l,k),k=1,8),int(ltemp_array(l,9))
               end do 
            else
               do l =1,lglocnt
                  write (des_data,'(7(2x,es12.5),I5)') (ltemp_array(l,k),k=1,7),int(ltemp_array(l,8))
               end do 
            end if 
         end if 
         deallocate (dprocbuf,drootbuf,iprocbuf,irootbuf,ltemp_array)
      end if  
! close the files 
      if (mype.eq.pe_io .or. bdist_io) then 
         close(des_data)
         close(des_ex)
         close(des_eps)
      end if
      return 

  999 write(*,"(/1x,70('*'),//,a,/,a,/1x,70('*'))")&
         ' From: write_des_tecplot ',&
         ' message: error opening des tecplot file. terminating run.'


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
! Close the file and keep
      CLOSE(BH_UNIT, STATUS="KEEP")

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

! Close the file and keep
      CLOSE(GT_UNIT, STATUS='KEEP')

      RETURN

 1001 FORMAT(/1X,70('*')//, ' From: WRITE_DES_THETA',/,&
         ' Message: ', A, ' already exists in the run',/10X,&
         'directory. Run terminated to prevent accidental overwriting',&
         /10X,'of files.',/1X,70('*')/)

      END SUBROUTINE WRITE_DES_THETA
!-----------------------------------------------
