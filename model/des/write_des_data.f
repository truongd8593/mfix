!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_DES_DATA                                         C
!>  Purpose: Writing DES output in Paraview format                      
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 26-Jul-06  C
!  Reviewer: Sreekanth Pannala                        Date: 31-Oct-06  C
!  Reviewer: Rahul Garg                               Dare: 01-Aug-07  C
!  Comments: Added one more output file containing averaged bed height C
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
      IMPLICIT NONE

      DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT 
      LOGICAL:: filexist, isopen
      LOGICAL, SAVE :: ONCE_OPEN = .FALSE.
      INTEGER DES_UNIT, LN, K , j, i, M, IJK, NP, IDIM
      INTEGER POS_Z, VEL_W
      CHARACTER*5 FILENAME
      CHARACTER*6 IPART
      CHARACTER*115 INUMBER
      INTEGER, SAVE :: unitht, tcount = 1, ROUTINE_COUNT = 0
      DOUBLE PRECISION, PARAMETER :: tmin = 5.d0
      DOUBLE PRECISION, DIMENSION(5000), SAVE :: bed_height_time, dt_time
      DOUBLE PRECISION :: height_avg, height_rms, AVG_EPS(JMAX2, MMAX), AVG_THETA(JMAX2, MMAX)
            CHARACTER*50 :: FILENAME_EXTRA, FILENAME_EXTRA2, FILENAME_DES
      CHARACTER*100 ::  TMP_CHARLINE2

!---------------------------------------------------------------------------

      INCLUDE 'function.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'ep_s2.inc'
      ROUTINE_COUNT  = ROUTINE_COUNT + 1
      DES_UNIT = 99
      IF(DEM_OUTPUT_DATA_TECPLOT) GOTO 200


      WRITE (FILENAME, 3020) IFI
      WRITE (IPART, 3021) PARTICLES
      IPART = ADJUSTL(IPART)
      IPART = TRIM(IPART)
      INUMBER = '     <Piece NumberOfPoints="'//IPART//'" NumberOfVerts="0" NumberOfLines="0" NumberOfStrips="0" NumberOfPolys="0">'
      OPEN(UNIT=DES_UNIT, FILE=TRIM(RUN_NAME)//'_DES_'//FILENAME//'.vtp', STATUS='NEW')

      POS_Z = 0
      VEL_W = 0

      WRITE(DES_UNIT,3029) '<?xml version="1.0"?>'
      WRITE(DES_UNIT,*) '<?Time =',S_TIME,'s?>'
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
      WRITE(DES_UNIT,*) '       <DataArray type="Float32" NAME="Position" NumberOfComponents="3" format="ascii">'
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

      return
 
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
           


 200   continue 

      FILENAME_EXTRA = TRIM(RUN_NAME)//'_DES_'//'EXTRA'//'.dat'
      FILENAME_EXTRA2 = TRIM(RUN_NAME)//'_AVG_EPS'//'.dat'
      FILENAME_DES = TRIM(RUN_NAME)//'_DES_DATA'//'.dat'


      IF(ROUTINE_COUNT.GT.1) THEN 
         
         INQUIRE(FILE=FILENAME_EXTRA,EXIST=filexist,OPENED=isopen)
         
         WRITE(TMP_CHARLINE2, *) 'VARIABLES = '
         DO M = 1, MMAX 
            WRITE(TMP_CHARLINE2, '(A,A,i2,A)') TRIM(TMP_CHARLINE2),', "`e_S_,_',M,'",'
         ENDDO
         DO M = 1, MMAX 
            WRITE(TMP_CHARLINE2, '(A,A,i2,A)') TRIM(TMP_CHARLINE2),', "`Q_S_,_',M,'",'
         ENDDO
         
         WRITE(TMP_CHARLINE2, '(A,A)') TRIM(TMP_CHARLINE2), ' "y" '

         IF (.NOT.filexist) THEN
            
            OPEN(unit=des_extra_unit,FILE=FILENAME_EXTRA, status='new')
            write(des_extra_unit,*)'VARIABLES= ',' "t" ', ' "MAX_NEIGH" ', ' "MAX_OVERLAP" ', ' "GRAN_ENERGY" ', ' "GRAN_TEMP" '
            
            OPEN(unit=des_volfrac_unit,FILE=FILENAME_EXTRA2, status='new')
            WRITE(des_volfrac_unit,*)  TRIM(TMP_CHARLINE2)

            OPEN(UNIT = DES_UNIT, FILE=FILENAME_DES, status='new')
            IF(DIMN.EQ.3) THEN 
               
               WRITE (DES_UNIT, '(10(A))')  'VARIABLES = ',  ' "x" ',  ' "y" ',&
               ' "z" ',  ' "vx" ',  ' "vy" ',  ' "vz" ',  ' "rad" ', ' "den" ', ' "mark" '
            ELSE 
               WRITE (DES_UNIT, '(10(A))')  'VARIABLES = ',  ' "x" ',  ' "y" ',&
              ' "vx" ',  ' "vy" ', ' "omega" ',  ' "rad" ', ' "den" ', ' "mark" '
            ENDIF
            
            ONCE_OPEN = .TRUE.
         ELSEIF(filexist.AND.(.NOT.isopen)) THEN 
            IF(RUN_TYPE.EQ.'NEW') THEN
               IF(.NOT.ONCE_OPEN) THEN 
                  OPEN(unit=des_extra_unit,FILE=FILENAME_EXTRA, status="replace")

                  write(des_extra_unit,*)'VARIABLES= ',' "t" ', ' "MAX_NEIGH" ',&
              ' "MAX_OVERLAP" ', ' "GRAN_ENERGY" ', ' "GRAN_TEMP" '
               
                  OPEN(unit=des_volfrac_unit,FILE=FILENAME_EXTRA2, status='replace')
               
                  WRITE(des_volfrac_unit,*)  TRIM(TMP_CHARLINE2)

                  
                  OPEN(UNIT = DES_UNIT, FILE=FILENAME_DES, status="replace")
                  IF(DIMN.EQ.3) THEN 
                     
                     WRITE (DES_UNIT, '(10(A))')  'VARIABLES = ',  ' "x" ',  ' "y" '&
             ,  ' "z" ',  ' "vx" ',  ' "vy" ',  ' "vz" ',  ' "rad" ', ' "den" ', ' "mark" '
                  ELSE 
                     WRITE (DES_UNIT, '(10(A))')  'VARIABLES = ',  ' "x" ',  ' "y" '&
             , ' "vx" ',  ' "vy" ', ' "omega" ',  ' "rad" ', ' "den" ', ' "mark" '
                  ENDIF
                  
                  ONCE_OPEN = .TRUE.
               ELSE
                  OPEN(unit=des_extra_unit,FILE=FILENAME_EXTRA,POSITION="append")
                  
                  OPEN(unit=des_volfrac_unit,FILE=FILENAME_EXTRA2, POSITION="append")
                  
                  OPEN(UNIT = DES_UNIT, FILE=FILENAME_DES, POSITION="append")
               ENDIF
            ELSE
               OPEN(unit=des_extra_unit,FILE=FILENAME_EXTRA,POSITION="append")
            
               OPEN(unit=des_volfrac_unit,FILE=FILENAME_EXTRA2, POSITION="append")
               
               OPEN(UNIT = DES_UNIT, FILE=FILENAME_DES, POSITION="append")
            END IF
         ENDIF
         



         WRITE(des_extra_unit,3022) s_time, NEIGH_MAX, OVERLAP_MAX,  &
         SUM(GLOBAL_GRAN_ENERGY(1:DIMN))*HALF, SUM(GLOBAL_GRAN_TEMP(1:DIMN))*1.d0/3.d0
         
         Write(DES_VOLFRAC_UNIT, *)'ZONE T = "', ROUTINE_COUNT, '"',',DATAPACKING=POINT, J=',JMAX, ', SOLUTIONTIME=', s_time
         
         WRITE (DES_UNIT, *) 'ZONE T = "', s_time, '"'!, SOLUTIONTIME=', s_time
         
          DO NP = 1, PARTICLES 
             IF(DIMN.EQ.3) THEN
                WRITE (DES_UNIT, '(10(2x,g12.5))') (DES_POS_NEW(NP, IDIM), IDIM = 1,DIMN), &
             (DES_VEL_NEW(NP, IDIM), IDIM = 1,DIMN), DES_RADIUS(NP), Ro_Sol(NP), mark_part(NP)
             ELSE
                WRITE (DES_UNIT, '(10(2x,g12.5))') (DES_POS_NEW(NP, IDIM), IDIM = 1,DIMN), &
             (DES_VEL_NEW(NP, IDIM), IDIM = 1,DIMN),OMEGA_NEW(NP,1), DES_RADIUS(NP), Ro_Sol(NP), mark_part(NP)
             ENDIF
          ENDDO
          
         DO J = JMIN1, JMAX1
            AVG_EPS(J,:) = ZERO 
            AVG_THETA(J,:) = ZERO
            DO K = KMIN1, KMAX1
               DO I = IMIN1, IMAX1
                  IJK = FUNIJK(I,J,K)
                  DO M = 1, MMAX
                     !PRINT*,'EP_S1=',EP_S(IJK,1)
                     AVG_EPS(J,M) =  AVG_EPS(J,M) + EP_S(IJK,M)
                     AVG_THETA(J,M) =  AVG_THETA(J,M) + DES_THETA(IJK,M)
                  ENDDO
               ENDDO
            ENDDO
            AVG_EPS(J,:) = AVG_EPS(J,:)/(IMAX*KMAX)
            AVG_THETA(J,:) = AVG_THETA(J,:)/(IMAX*KMAX)
            WRITE(DES_VOLFRAC_UNIT,3023) (AVG_EPS(J,M), M = 1, MMAX) ,  (AVG_THETA(J,M), M = 1, MMAX),0.5d0*(YN(J)+YN(J-1))
         ENDDO
         !READ(*,*)
         CLOSE(des_extra_unit, status = "keep")
         CLOSE(des_volfrac_unit, status = "keep")
         CLOSE(des_unit, status = "keep")
      ENDIF
      
      RETURN


 3029 FORMAT(A21)
 3030 FORMAT(A102)
 3031 FORMAT(A115)
 3032 FORMAT(A87)
 3022 FORMAT(g12.5,2x,i4,2x,3(g12.5,2x))
 3023 FORMAT(20(g12.5,2x))
 3020 FORMAT(I5.5)
 3021 FORMAT(I6)
      END SUBROUTINE WRITE_DES_DATA 
