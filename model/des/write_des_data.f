!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_DES_DATA                                         C
!  Purpose: Writing DES output in Paraview format                      C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 26-Jul-06  C
!  Reviewer: Sreekanth Pannala                        Date: 31-Oct-06  C
!  Reviewer: Rahul Garg                               Dare: 01-Aug-07  C
!  Comments: Added one more output file containing averaged bed height C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE WRITE_DES_DATA

      USE param1      
      USE discretelement
      USE run
      IMPLICIT NONE

      DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT 

      INTEGER DES_UNIT, LN, K , j, i
      INTEGER POS_Z, VEL_W
      CHARACTER*5 FILENAME
      CHARACTER*6 IPART
      CHARACTER*115 INUMBER
      INTEGER, SAVE :: unitht, tcount = 1
      DOUBLE PRECISION, PARAMETER :: tmin = 5.d0
      LOGICAL:: filexist, isopen
      DOUBLE PRECISION, DIMENSION(5000), SAVE :: bed_height_time, dt_time
      
      DOUBLE PRECISION :: height_avg, height_rms
!---------------------------------------------------------------------------

      DES_UNIT = 99
      WRITE (FILENAME, 3020) IFI
      WRITE (IPART, 3021) PARTICLES
      IPART = ADJUSTL(IPART)
      IPART = TRIM(IPART)
      INUMBER = '     <Piece NumberOfPoints="'//IPART//'" NumberOfVerts="0" NumberOfLines="0" NumberOfStrips="0" NumberOfPolys="0">'
      OPEN(UNIT=DES_UNIT, FILE=TRIM(RUN_NAME)//'_DES_'//FILENAME//'.vtp', STATUS='NEW')

      POS_Z = 0
      VEL_W = 0


      WRITE(DES_UNIT,3029) '<?xml version="1.0"?>'
!     WRITE(DES_UNIT,*) '<?Time =',S_TIME,'s?>'
      WRITE(DES_UNIT,3030) ' <VTKFile type="PolyData" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">'
      WRITE(DES_UNIT,*) '  <PolyData>'
      WRITE(DES_UNIT,3031) INUMBER
      WRITE(DES_UNIT,*) '      <PointData Scalars="Diameter">'
      WRITE(DES_UNIT,*) '        <DataArray type="Float64" Name="Diameter" format="ascii">'
      DO LN = 1, PARTICLES
         WRITE (DES_UNIT,*) (real(2D0 * DES_RADIUS(LN))) 
      END DO
      WRITE(DES_UNIT,*) '       </DataArray>'
      WRITE(DES_UNIT,3032) '       <DataArray type="Float64" Name="Velocity" NumberOfComponents="3" format="ascii">'
      IF(DIMN.EQ.2) THEN
         DO LN = 1, PARTICLES
            WRITE (DES_UNIT,*) (real(DES_VEL_NEW(LN,K)),K=1,DIMN), VEL_W 
         END DO
      ELSE
         DO LN = 1, PARTICLES
            WRITE (DES_UNIT,*) (real(DES_VEL_NEW(LN,K)),K=1,DIMN) 
         END DO
      END IF
      WRITE(DES_UNIT,*) '       </DataArray>'
      WRITE(DES_UNIT,*) '     </PointData>'
      WRITE(DES_UNIT,*) '     <CellData>'
      WRITE(DES_UNIT,*) '     </CellData>'
      WRITE(DES_UNIT,*) '     <Points>'
      WRITE(DES_UNIT,*) '       <DataArray type="Float32" NumberOfComponents="3" format="ascii">'
      IF(DIMN.EQ.2) THEN
         DO LN = 1, PARTICLES
            WRITE (DES_UNIT,*) (real(DES_POS_NEW(LN,K)),K=1,DIMN), POS_Z 
         END DO
      ELSE
         DO LN = 1, PARTICLES
            WRITE (DES_UNIT,*) (real(DES_POS_NEW(LN,K)),K=1,DIMN) 
         END DO
      END IF
      WRITE(DES_UNIT,*) '       </DataArray>'
      WRITE(DES_UNIT,*) '     </Points>'
      WRITE(DES_UNIT,*) '     <Verts>'
      WRITE(DES_UNIT,*) '       <DataArray type="Int32" Name="connectivity" format="ascii">'
      WRITE(DES_UNIT,*) '       </DataArray>'
      WRITE(DES_UNIT,*) '       <DataArray type="Int32" Name="offsets" format="ascii">'
      WRITE(DES_UNIT,*) '       </DataArray>'
      WRITE(DES_UNIT,*) '     </Verts>'
      WRITE(DES_UNIT,*) '     <Lines>'
      WRITE(DES_UNIT,*) '       <DataArray type="Int32" Name="connectivity" format="ascii">'
      WRITE(DES_UNIT,*) '       </DataArray>'
      WRITE(DES_UNIT,*) '       <DataArray type="Int32" Name="offsets" format="ascii">'
      WRITE(DES_UNIT,*) '       </DataArray>'
      WRITE(DES_UNIT,*) '     </Lines>'
      WRITE(DES_UNIT,*) '     <Strips>'
      WRITE(DES_UNIT,*) '       <DataArray type="Int32" Name="connectivity" format="ascii">'
      WRITE(DES_UNIT,*) '       </DataArray>'
      WRITE(DES_UNIT,*) '       <DataArray type="Int32" Name="offsets" format="ascii">'
      WRITE(DES_UNIT,*) '       </DataArray>'
      WRITE(DES_UNIT,*) '     </Strips>'
      WRITE(DES_UNIT,*) '     <Polys>'
      WRITE(DES_UNIT,*) '       <DataArray type="Int32" Name="connectivity" format="ascii">'
      WRITE(DES_UNIT,*) '       </DataArray>'
      WRITE(DES_UNIT,*) '       <DataArray type="Int32" Name="offsets" format="ascii">'
      WRITE(DES_UNIT,*) '       </DataArray>'
      WRITE(DES_UNIT,*) '     </Polys>'
      WRITE(DES_UNIT,*) '   </Piece>'
      WRITE(DES_UNIT,*) ' </PolyData>'
      WRITE(DES_UNIT,*) ' </VTKFile>'
    
      CLOSE(DES_UNIT)

      IFI = IFI+1

      
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
               end do
               height_avg = height_avg/(tcount-1)
               do i = 1, tcount-1,1
                  height_rms = height_rms + ((bed_height_time(i)&
                       &-height_avg)**2)!*dt_time(i)
               end do
               
               
               height_rms = sqrt(height_rms/(tcount-1))
            end if
         end if
         
      end if
      
      
      INQUIRE(FILE='bed_height.dat',EXIST=filexist,OPENED=isopen)
      IF (.NOT.filexist.OR.(.NOT.isopen)) THEN

         unitht = 400
       OPEN(unit=unitht,file='bed_height.dat',form='formatted', status='replace')
    END IF
      write(*, '(10(2x,e20.12))') time, (bed_height(j) , j = 1, size(bed_height,1))
      !read(*,*)
      write(unitht, '(10(2x,e20.12))') time, (bed_height(j) , j = 1,&
           & size(bed_height,1)), height_avg, height_rms

      
      RETURN

 3029 FORMAT(A21)
 3030 FORMAT(A102)
 3031 FORMAT(A115)
 3032 FORMAT(A87)

 3020 FORMAT(I5.5)
 3021 FORMAT(I6)
      END SUBROUTINE WRITE_DES_DATA 
