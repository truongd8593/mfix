!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: stl_dbg_des                                            !
!  Author: Rahul Garg                                 Date: 24-Oct-13  !
!                                                                      !
!  Purpose: Random functions for debugging STLs with DES.              !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      MODULE STL_DBG_DES

      IMPLICIT NONE

! Use this module only to define functions and subroutines.
      CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: STL_DBG_DG_REPORT                                       !
!  Author: Rahul Garg                                 Date: 24-Oct-13  !
!                                                                      !
!  Purpose: Reports the total number of facets in each DES grid cell.  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE STL_DBG_DG_REPORT

      use desgrid, only: DG_IJKSTART2, DG_IJKEND2

      use stl, only: FACETS_AT_DG

      use compar, only: myPE, numPEs

      IMPLICIT NONE

      INTEGER :: IJK, TOTAL_FACETS, LC

      CHARACTER(LEN=100) :: FN

      IF(numPEs == 1) THEN
         WRITE(FN,'("FACETS_DG_GRID.DAT")')
      ELSE
         WRITE(FN,'("FACETS_DG_GRID_",I5.5,".DAT")') myPE
      ENDIF

      OPEN(1001, file=TRIM(FN))

      DO IJK=DG_IJKSTART2, DG_IJKEND2
         TOTAL_FACETS = FACETS_AT_DG(IJK)%COUNT  
         IF(TOTAL_FACETS < 1) CYCLE
         WRITE(1001,2000) IJK, TOTAL_FACETS
         DO LC=1, TOTAL_FACETS
            WRITE(1001,'(2x,I10)') FACETS_AT_DG(IJK)%ID(LC)
         ENDDO
      ENDDO

      CLOSE(1001, STATUS = "keep")

 2000 FORMAT(2/2x,'DG CELL: ',I10,3x,'Total STLs: ',I4)

      RETURN
      END SUBROUTINE STL_DBG_DG_REPORT


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: STL_DBG_WRITE_INPUT_FACETS                              !
!  Author: Rahul Garg                                 Date: 24-Oct-13  !
!                                                                      !
!  Purpose: Write back out the STL files read from input files.        !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE STL_DBG_WRITE_INPUT_FACETS

! Number of facets 
      use stl, only: N_FACETS
! Facet Vertices and normal
      use stl, only: VERTEX, NORM_FACE
! Processor rank and rank of IO
      use compar, only: myPE, PE_IO

      IMPLICIT NONE

      INTEGER :: LC

      IF(myPE /= PE_IO) RETURN

      OPEN(UNIT=444, FILE='INPUT_FACETS.stl')

      WRITE(444,*) 'solid vcg'
      DO LC = 1, N_FACETS
         WRITE(444,*) '   facet normal ', NORM_FACE(:,LC)
         WRITE(444,*) '      outer loop'
         WRITE(444,*) '         vertex ', VERTEX(1,1:3,LC)
         WRITE(444,*) '         vertex ', VERTEX(2,1:3,LC)
         WRITE(444,*) '         vertex ', VERTEX(3,1:3,LC)
         WRITE(444,*) '      endloop'
         WRITE(444,*) '   endfacet'
      ENDDO
      WRITE(444,*)'endsolid vcg'

      CLOSE(555)

      RETURN
      END SUBROUTINE STL_DBG_WRITE_INPUT_FACETS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: STL_DBG_WRITE_AUTO_FACETS                               !
!  Author: J.Musser                                   Date: 03-Nov-15  !
!                                                                      !
!  Purpose: Write out auto-generated facets for BCs, ISs, and default  !
!  walls.                                                              !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE STL_DBG_WRITE_AUTO_FACETS

! Number of facets (plus) auto generated facets
      use stl, only: N_FACETS, N_FACETS_DES
! Facet Vertices and normal
      use stl, only: VERTEX, NORM_FACE
! Processor rank and rank of IO
      use compar, only: myPE, PE_IO

      IMPLICIT NONE

      INTEGER :: LC

! Raw STL data is common to all ranks.
      IF(myPE /= PE_IO) RETURN

      OPEN(UNIT=444, FILE='AUTO_FACETS.stl')

      WRITE(444,*) 'solid vcg'
      DO LC = N_FACETS+1, N_FACETS_DES
         WRITE(444,*) '   facet normal ', NORM_FACE(:,LC)
         WRITE(444,*) '      outer loop'
         WRITE(444,*) '         vertex ', VERTEX(1,1:3,LC)
         WRITE(444,*) '         vertex ', VERTEX(2,1:3,LC)
         WRITE(444,*) '         vertex ', VERTEX(3,1:3,LC)
         WRITE(444,*) '      endloop'
         WRITE(444,*) '   endfacet'
      ENDDO
      WRITE(444,*)'endsolid vcg'

      CLOSE(555)

      RETURN
      END SUBROUTINE STL_DBG_WRITE_AUTO_FACETS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: DEBUG_write_stl_from_grid_facet                         !
!  Author: Rahul Garg                                 Date: 24-Oct-13  !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE STL_DBG_WRITE_STL_FROM_DG(WRITE_EACH_CELL)

! DES grid loop bounds
      use desgrid, only: DG_IJKSTART2, DG_IJKEND2
! Max numer of STLs
      use stl, only: DIM_STL
! Facets binned to DES grid
      use stl, only: FACETS_AT_DG
! Facet normal and vertex data
      use stl, only: NORM_FACE, VERTEX
! Processor ID and total proc count
      USE compar, only: myPE, numPEs

      IMPLICIT NONE

      LOGICAL, INTENT(IN), OPTIONAL :: WRITE_EACH_CELL

      INTEGER :: IJK, LC1, LC2
      CHARACTER(LEN=200) :: FN
      LOGICAL :: EACH_CELL

      LOGICAL, ALLOCATABLE :: WRITE_FACET(:)

! Initialize flag.
      EACH_CELL = .FALSE.
      IF(present(WRITE_EACH_CELL)) EACH_CELL = WRITE_EACH_CELL

      ALLOCATE (WRITE_FACET(DIM_STL))
      WRITE_FACET = .TRUE.

      IF(numPEs == 1) THEN
         WRITE(FN,'("DG_FACETS.stl")')
      ELSE
         WRITE(FN,'("DG_FACETS_",I5.5,".stl")') MYPE
      ENDIF

      OPEN(UNIT=444, FILE=trim(FN))

      write(444,*)'solid vcg'
      DO IJK=DG_IJKSTART2,DG_IJKEND2
         IF(FACETS_AT_DG(IJK)%COUNT< 1) CYCLE

         IF(EACH_CELL) CALL WRITE_STLS_THIS_DG(IJK)

         DO LC1 = 1, FACETS_AT_DG(IJK)%COUNT
            LC2 = FACETS_AT_DG(IJK)%ID(LC1)
            IF(WRITE_FACET(LC2)) THEN
               write(444,*) '   facet normal ', NORM_FACE(:,LC2)
               write(444,*) '      outer loop'
               write(444,*) '         vertex ', VERTEX(1,:,LC2)
               write(444,*) '         vertex ', VERTEX(2,:,LC2)
               write(444,*) '         vertex ', VERTEX(3,:,LC2)
               write(444,*) '      endloop'
               write(444,*) '   endfacet'
               WRITE_FACET(LC2) = .FALSE.
            ENDIF
         ENDDO
      ENDDO
      write(444,*)'endsolid vcg'

      close(444)

      DEALLOCATE (WRITE_FACET)

      RETURN
      END SUBROUTINE STL_DBG_WRITE_STL_FROM_DG




!----------------------------------------------------------------------!
!                                                                      !
!                                                                      !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE write_stls_this_dg(dg)

! STL Vertices
      use stl, only: VERTEX
! STL Facet normals
      use stl, only: NORM_FACE
! Facets binned to DES grid
      use stl, only: facets_at_dg


      IMPLICIT NONE
!-----------------------------------------------
      integer, intent(in) :: dg
      integer :: lc, nf

      logical :: EXISTS
      character(len=6) :: IDX

      write(idx,"(I6.6)") dg
      open(unit=555, file='dg_'//idx//'.stl', status='UNKNOWN')

      write(555,*) 'solid vcg'

      DO lc=1, facets_at_dg(dg)%count

         NF = facets_at_dg(dg)%id(lc)

         write(555,*) '   facet normal ', NORM_FACE(:,NF)
         write(555,*) '      outer loop'
         write(555,*) '         vertex ', VERTEX(1,1:3,NF)
         write(555,*) '         vertex ', VERTEX(2,1:3,NF)
         write(555,*) '         vertex ', VERTEX(3,1:3,NF)
         write(555,*) '      endloop'
         write(555,*) '   endfacet'
      enddo
      close(555)

      RETURN
      END SUBROUTINE write_stls_this_dg


!----------------------------------------------------------------------!
!                                                                      !
!                                                                      !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE write_this_stl(this)


! STL Vertices
      use stl, only: VERTEX
! STL Facet normals
      use stl, only: NORM_FACE
      use compar, only: myPE

      IMPLICIT NONE
!-----------------------------------------------
      integer, intent(in) :: this

      logical :: EXISTS
      character(len=4) :: IDX
      character(len=4) :: IPE


      write(idx,"(I4.4)") this
      write(ipe,"(I4.4)") myPE
      open(unit=555, file='geo_'//idx//'_'//IPE//'.stl',&
         status='UNKNOWN')
      write(555,*) 'solid vcg'
      write(555,*) '   facet normal ', NORM_FACE(:,this)
      write(555,*) '      outer loop'
      write(555,*) '         vertex ', VERTEX(1,1:3,this)
      write(555,*) '         vertex ', VERTEX(2,1:3,this)
      write(555,*) '         vertex ', VERTEX(3,1:3,this)
      write(555,*) '      endloop'
      write(555,*) '   endfacet'
      close(555)


      RETURN
      END SUBROUTINE write_this_stl

      END MODULE STL_DBG_DES


