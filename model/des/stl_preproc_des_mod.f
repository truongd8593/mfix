!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: des_stl_functions                                      !
!  Author: Rahul Garg                                 Date: 24-Oct-13  !
!                                                                      !
!  Purpose: This module containd routines for geometric interaction    !
!  required for STL files.                                             !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      MODULE STL_PREPROC_DES

      use des_stl_functions, only: TestTriangleAABB

      use desgrid
      use error_manager

      IMPLICIT NONE

! Use this module only to define functions and subroutines.
      CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: des_stl_preprocessing                                   !
!  Author: Rahul Garg                                 Date: 24-Oct-13  !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DES_STL_PREPROCESSING

      USE cutcell, only: use_stl

      USE discretelement, only: DES_CONVERT_BOX_TO_FACETS
      USE discretelement, only: STL_FACET_TYPE
      USE discretelement, only: FACET_TYPE_NORMAL
      USE discretelement, only: FACET_TYPE_MI
      USE discretelement, only: FACET_TYPE_PO
      USE discretelement, only: COUNT_FACET_TYPE_NORMAL
      USE discretelement, only: COUNT_FACET_TYPE_MI
      USE discretelement, only: COUNT_FACET_TYPE_PO

      USE desgrid

      USE param
      USE stl

      implicit none

      INTEGER :: IJK, COUNT, NF

! Pre-procssing for the des in order to assign facets to grid cells.
      WRITE(ERR_MSG,"('Pre-Processing geometry for DES.')")
      CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)

      CALL ALLOCATE_DES_STL_ARRAYS

      N_FACETS_DES = 0

! Set N_facets_des to add any more facets needed by dem and not to.
! contaminate the Eulerian-Eulerian CG stuff
      IF(USE_STL) N_FACETS_DES = N_FACETS

! Set all the facets generated so far as normal facet types
      STL_FACET_TYPE(1:N_FACETS) = FACET_TYPE_NORMAL


      CALL BIN_FACETS_TO_GRID_DES

      DO IJK = 1, DG_IJKSIZE2
         COUNT = LIST_FACET_AT_DES(IJK)%COUNT_FACETS
         IF(COUNT.eq.0) DEALLOCATE(LIST_FACET_AT_DES(IJK)%FACET_LIST)
      ENDDO

      COUNT_FACET_TYPE_NORMAL = 0
      COUNT_FACET_TYPE_PO     = 0
      COUNT_FACET_TYPE_MI     = 0
      DO NF = 1, N_FACETS_DES
         IF(STL_FACET_TYPE(NF) == FACET_TYPE_NORMAL) &
            COUNT_FACET_TYPE_NORMAL = COUNT_FACET_TYPE_NORMAL + 1

         IF(STL_FACET_TYPE(NF).EQ.FACET_TYPE_MI) &
            COUNT_FACET_TYPE_MI = COUNT_FACET_TYPE_MI + 1

         IF(STL_FACET_TYPE(NF).EQ.FACET_TYPE_PO) &
            COUNT_FACET_TYPE_PO = COUNT_FACET_TYPE_PO + 1
      ENDDO

!      CALL DEBUG_WRITE_GRID_FACEINFO
!      CALL DEBUG_WRITE_STL_FROM_GRID_FACET(WRITE_FACETS_EACH_CELL=.FALSE.)

!      CALL DEBUG_WRITE_ALL_READIN_FACETS

      RETURN
      END SUBROUTINE DES_STL_PREPROCESSING

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: ALLOCATE_DES_STL_ARRAYS                                 !
!  Author: Rahul Garg                                 Date: 24-Oct-13  !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE ALLOCATE_DES_STL_ARRAYS

      USE desgrid
      USE discretelement, only: STL_Facet_type
      USE discretelement, only: CELLNEIGHBOR_FACET
      USE discretelement, only: CELLNEIGHBOR_FACET_MAX
      USE discretelement, only: CELLNEIGHBOR_FACET_NUM
      USE discretelement, only: CELLNEIGHBOR_FACET_NUM
      USE stl, only: LIST_FACET_AT_DES
      USE stl, only: NO_NEIGHBORING_FACET_DES
      USE param
      USE stl, only: DIM_STL, MAX_FACETS_PER_CELL_DES

      IMPLICIT NONE

      integer :: ii,ijk

      ALLOCATE(LIST_FACET_AT_DES(DG_IJKSIZE2))

      Allocate(  CELLNEIGHBOR_FACET (DG_IJKSIZE2) )
      Allocate(  CELLNEIGHBOR_FACET_MAX (DG_IJKSIZE2) )
      Allocate(  CELLNEIGHBOR_FACET_NUM (DG_IJKSIZE2) )
      DO II = 1, DG_IJKSIZE2
         CELLNEIGHBOR_FACET_MAX(II) = 4
         allocate(CELLNEIGHBOR_FACET(II)%P(CELLNEIGHBOR_FACET_MAX(II)))
         allocate(CELLNEIGHBOR_FACET(II)%EXTENTDIR(CELLNEIGHBOR_FACET_MAX(II)))
         allocate(CELLNEIGHBOR_FACET(II)%EXTENTMIN(CELLNEIGHBOR_FACET_MAX(II)))
         allocate(CELLNEIGHBOR_FACET(II)%EXTENTMAX(CELLNEIGHBOR_FACET_MAX(II)))
         CELLNEIGHBOR_FACET_NUM(II) = 0
      ENDDO

      DO IJK = 1, DG_IJKSIZE2
         LIST_FACET_AT_DES(IJK)%COUNT_FACETS = 0
         ALLOCATE(LIST_FACET_AT_DES(IJK)%FACET_LIST(MAX_FACETS_PER_CELL_DES))
      ENDDO

      ALLOCATE(NO_NEIGHBORING_FACET_DES(DG_IJKSIZE2))
      NO_NEIGHBORING_FACET_DES = .false.

      Allocate(stl_facet_type(DIM_STL))

      END SUBROUTINE ALLOCATE_DES_STL_ARRAYS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: bin_facets_to_grid_des                                  !
!  Author: Rahul Garg                                 Date: 24-Oct-13  !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      Subroutine BIN_FACETS_TO_GRID_DES

      USE desgrid
      USE discretelement, only: dimn, xe, yn, zt
      USE functions
      USE geometry
      USE indices
      USE mpi_utility
      USE param1, only: zero
      Use stl, only: tol_stl, n_facets_des
      Use stl, only: vertex, NO_NEIGHBORING_FACET_DES, LIST_FACET_AT_DES

      IMPLICIT NONE

! DES Grid cell index.
      INTEGER :: IJK, IJK2
! Loop counters:
      INTEGER :: I, I1, I2, II  ! X-axis
      INTEGER :: J, J1, J2, JJ  ! Y-axis
      INTEGER :: K, K1, K2, KK  ! Z-axis
      INTEGER :: N              ! STLs
! Generic accumulator
      INTEGER :: COUNT_FAC

! Maximum and minimum extents of the indexed STL
      DOUBLE PRECISION:: X1,Y1,Z1
      DOUBLE PRECISION:: X2,Y2,Z2



      DO N = 1,N_FACETS_DES

         X1 = minval(VERTEX(1:3,1,N))
         X2 = maxval(VERTEX(1:3,1,N))
         Y1 = minval(VERTEX(1:3,2,N))
         Y2 = maxval(VERTEX(1:3,2,N))
         Z1 = minval(VERTEX(1:3,3,N))
         Z2 = maxval(VERTEX(1:3,3,N))

         I1 = DG_IEND2
         I2 = DG_ISTART2
         IF(X2>=ZERO .AND. X1<=XLENGTH+TOL_STL) THEN
            I1 = iofpos(X1)
            I2 = iofpos(1/dg_dxinv+X2)
         ENDIF

         J1 = DG_JEND2
         J2 = DG_JSTART2
         IF(Y2>=ZERO .AND. Y1<=YLENGTH+TOL_STL) THEN
            J1 = jofpos(Y1)
            J2 = jofpos(1/dg_dyinv+Y2)
         ENDIF

         K1 = DG_KEND2
         K2 = DG_KSTART2
         IF(DO_K) THEN
            IF(Z2>=ZERO .AND. Z1<=ZLENGTH+TOL_STL) THEN
               K1 = kofpos(Z1)
               K2 = kofpos(1/dg_dzinv+Z2)
            ENDIF
         ENDIF


         DO K=K1,K2
         DO J=J1,J2
         DO I=I1,I2
            IF(dg_is_ON_myPE_plus1layers(I,J,K)) THEN
               IJK = DG_FUNIJK(I,J,K)
               CALL ADD_FACET_FOR_DES(I,J,K,IJK,N)
            ENDIF
         ENDDO
         ENDDO
         ENDDO

      ENDDO

      DO K = DG_KSTART1, DG_KEND1
      DO J = DG_JSTART1, DG_JEND1
      DO I = DG_ISTART1, DG_IEND1

         IJK = DG_FUNIJK(I, J, K)

         I1 =  MAX( I - 1, DG_ISTART2)
         I2 =  MIN( I + 1, DG_IEND2)

         J1 =  MAX( J - 1, DG_JSTART2)
         J2 =  MIN (J + 1, DG_JEND2)

         K1 =  MAX( K - 1, DG_KSTART2)
         K2 =  MIN (K + 1, DG_KEND2)

         NO_NEIGHBORING_FACET_DES(IJK)  = .FALSE.

         COUNT_FAC = 0
         DO KK = K1, k2
         DO JJ = J1, j2
         DO II = I1, i2
            IJK2 = DG_FUNIJK(II, JJ, KK)
            COUNT_FAC = COUNT_FAC + LIST_FACET_AT_DES(IJK2)%COUNT_FACETS
         ENDDO
         ENDDO
         ENDDO

         IF(COUNT_FAC == 0) NO_NEIGHBORING_FACET_DES(IJK) = .TRUE.

      ENDDO
      ENDDO
      ENDDO

      RETURN
      END SUBROUTINE BIN_FACETS_TO_GRID_DES

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: ADD_FACET_for_des                                      C
!  Purpose: Add facet to list in IJK scalar cell for the               C
!           discrete modules.                                          C
!                                                                      C
!  Author: Rahul Garg                              Date: 24-Oct-13     C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE ADD_FACET_FOR_DES(I,J,K,IJK,N)

      USE param1, only: zero, one
      use geometry, only: DO_K, ZLENGTH

      use run, only: RUN_NAME

      use compar, only: myPE
      use compar, only: NODESi, NODESj, NODESk

      use desgrid, only: dg_xstart, dg_dxinv
      use desgrid, only: dg_ystart, dg_dyinv
      use desgrid, only: dg_zstart, dg_dzinv

      Use stl
      use error_manager

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: I,J,K,IJK, N
      !LOCAL VARIABLES

      INTEGER ::  CURRENT_COUNT, COUNT

      double precision ::   box_origin(3), box_extents(3)
      Logical :: sa_exist
      integer :: sep_axis

      CHARACTER(LEN=100) :: FNAME
      integer :: stl_unit, fid

      stl_unit = 1001

      BOX_ORIGIN(1) = DG_XSTART + (I-DG_ISTART1)/DG_DXINV
      BOX_EXTENTS(1) = 1.0/DG_DXINV

      BOX_ORIGIN(2) = DG_YSTART + (J-DG_JSTART1)/DG_DYINV
      BOX_EXTENTS(2) = 1.0/DG_DYINV

      IF(DO_K)THEN
         BOX_ORIGIN(3) = DG_ZSTART + (K-DG_KSTART1)/DG_DZINV
         BOX_EXTENTS(3) = 1.0/DG_DZINV
      ELSE
         BOX_ORIGIN(3) = 0.0D0
         BOX_EXTENTS(3) = ZLENGTH
      ENDIF

! Use the separating axis test to determine if the cell and facet cannot
! intersect one another.
      CALL TESTTRIANGLEAABB(VERTEX(:,:,N), NORM_FACE(:,N),             &
         BOX_ORIGIN(:), BOX_EXTENTS(:), SA_EXIST, SEP_AXIS,I,J,K )

! A sparating axis exists so the cell and the triangle do not intersect.
      IF(SA_EXIST) RETURN

      CURRENT_COUNT = LIST_FACET_AT_DES(IJK)%COUNT_FACETS

      IF(CURRENT_COUNT .LT. MAX_FACETS_PER_CELL_DES) THEN

         LIST_FACET_AT_DES(IJK)%COUNT_FACETS = CURRENT_COUNT+1
         LIST_FACET_AT_DES(IJK)%FACET_LIST(CURRENT_COUNT+1) = N

      ELSE
         CALL INIT_ERR_MSG("des_stl_functions_mod::add_facets_for_des")

         WRITE(err_msg, 200) MAX_FACETS_PER_CELL_DES, IJK, &
            I, J, K, mype, .FALSE.
         CALL flush_err_msg(footer = .false.)


         write(err_msg, *) "current_list for this cell is"
         CALL flush_err_msg(header = .false., footer = .false.)


         IF(nodesI*nodesJ*nodesK.gt.1) then

            WRITE(fname,'(A,"_TROUBLE_CELL",A, I5.5, 3(A,I5.5), ".stl")') &
            TRIM(RUN_NAME), '_pid', mype, '_I', I, '_J', J, '_K', K
         else

            WRITE(fname,'(A,"_TROUBLE_CELL", 3(A,I5.5), ".stl")') &
            TRIM(RUN_NAME), '_I', I, '_J', J, '_K', K
         endif

         open(stl_unit, file = fname, form='formatted')
         write(stl_unit,*)'solid vcg'

         DO COUNT  = 1, CURRENT_COUNT

             !write(*, '(/,I10)') LIST_FACET_AT(IJK)%FACET_LIST(COUNT)
            FID = LIST_FACET_AT_DES(IJK)%FACET_LIST(COUNT)
            write(err_msg, '(I10)') FID

            CALL flush_err_msg(header = .false., footer = .false.)


             write(stl_unit,*) '   facet normal ', NORM_FACE(1:3,FID)
             write(stl_unit,*) '      outer loop'
             write(stl_unit,*) '         vertex ', VERTEX(1,1:3,FID)
             write(stl_unit,*) '         vertex ', VERTEX(2,1:3,FID)
             write(stl_unit,*) '         vertex ', VERTEX(3,1:3,FID)
             write(stl_unit,*) '      endloop'
             write(stl_unit,*) '   endfacet'

          ENDDO

          write(stl_unit,*)'endsolid vcg'
          close(stl_unit, status = 'keep')


          write(err_msg, *) 'Stopping'
          CALL flush_err_msg(header = .false., abort = .true.)
      ENDIF



 200  FORMAT('ERROR MESSAGE FROM CUT_CELL_PREPROCESSING', /10x,        &
         'INCREASE MAX_FACETS_PER_CELL_DES from the current value of', &
         I3, /10x,'Happening for cell IJK, I, J, K = ', 4(2x, i5),/10X,&
         'mype, Is on myPe? ', I6, L2, /10X,'see the file ',&
         'TROUBLE_CELL for all the current facets in this cell')

      END SUBROUTINE ADD_FACET_FOR_DES



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: DEBUG_WRITE_GRID_FACEINFO                               !
!  Author: Rahul Garg                                 Date: 24-Oct-13  !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE DEBUG_WRITE_GRID_FACEINFO
      USE param
      USE param1
      USE run
      USE geometry
      USE indices
      USE compar
      USE stl
      USE functions

      IMPLICIT NONE

      INTEGER :: CELL_ID, I, J, K, COUNT, COUNT_FACETS

      CHARACTER(LEN=100) :: FN

      IF(NODESI*NODESJ*NODESK == 1) THEN
         WRITE(FN,'("FACETS_DG_GRID.DAT")')
      ELSE
         WRITE(FN,'("FACETS_DG_GRID_",I5.5,".DAT")') myPE
      ENDIF

      OPEN(1001, file = TRIM(FN), form ='formatted')

      DO K=DG_KSTART2, DG_KEND2
      DO J=DG_JSTART2, DG_JEND2
      DO I=DG_ISTART2, DG_IEND2

         CELL_ID = DG_FUNIJK(I,J,K)
         COUNT_FACETS =  LIST_FACET_AT_DES(CELL_ID)%COUNT_FACETS
         IF(COUNT_FACETS.eq.0) cycle

         WRITE(1001,2000) CELL_ID, I, J, K,  COUNT_FACETS

 2000 FORMAT(50('*'),/2X,'CELL IJK, I, J, K =        = ', i20, 2x,     &
         4(2x,i10),/2X,'TOTAL FACETS',18(' '),'= ',3(2x, i10))

         DO COUNT = 1, COUNT_FACETS
            WRITE(1001, '(2x, i20)')                                   &
               LIST_FACET_AT_DES(CELL_ID)%FACET_LIST(COUNT)
         ENDDO

      ENDDO
      ENDDO
      ENDDO

      CLOSE(1001, STATUS = "keep")

      RETURN
      END SUBROUTINE DEBUG_WRITE_GRID_FACEINFO



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: DEBUG_write_all_readin_facets                           !
!  Author: Rahul Garg                                 Date: 24-Oct-13  !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      Subroutine  DEBUG_write_all_readin_facets

      USE stl
      USE compar

      IMPLICIT NONE

      INTEGER :: N
!      CHARACTER(LEN=100) :: FN


      WRITE(ERR_MSG, 2000) N_FACETS, N_FACETS_DES
      CALL FLUSH_ERR_MSG(HEADER=.FALSE.,FOOTER=.FALSE.)

 2000 FORMAT('Saving STL geometry on DES grid: FACETS_TO_DG.stl',      &
         11x,'Facet Count:',I10,/2x,'DES Grid Facet Count:',I10)

      OPEN(UNIT=444, FILE='FACETS_TO_DG.stl')
      write(444,*)'solid vcg'

      DO N = 1, N_FACETS_DES
         write(444,*) '   facet normal ', NORM_FACE(1:3,N)
         write(444,*) '      outer loop'
         write(444,*) '         vertex ', VERTEX(1,1:3,N)
         write(444,*) '         vertex ', VERTEX(2,1:3,N)
         write(444,*) '         vertex ', VERTEX(3,1:3,N)
         write(444,*) '      endloop'
         write(444,*) '   endfacet'
      ENDDO
      write(444,*)'endsolid vcg'
      close(444)



      END Subroutine DEBUG_write_all_readin_facets


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: DEBUG_write_stl_from_grid_facet                         !
!  Author: Rahul Garg                                 Date: 24-Oct-13  !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DEBUG_WRITE_STL_FROM_GRID_FACET(WRITE_FACETS_EACH_CELL)

      use run
      USE stl

      USE functions
      USE geometry
      USE indices
      USE compar
      USE discretelement, only: STL_FACET_TYPE, FACET_TYPE_NORMAL, &
           FACET_TYPE_MI, FACET_TYPE_PO, COUNT_FACET_TYPE_NORMAL, &
           COUNT_FACET_TYPE_MI, COUNT_FACET_TYPE_PO

      IMPLICIT NONE

      LOGICAL, INTENT(IN),optional  :: WRITE_FACETS_EACH_CELL
      INTEGER ::  CELL_ID, N, I, J, K, COUNT, COUNT_FACETS, W_UNIT
      CHARACTER(LEN=200) :: FN, FN_PO
      LOGICAL :: write_each_cell
      LOGICAL, DIMENSION(:), allocatable :: FACET_WRITTEN

      ALLOCATE (FACET_WRITTEN(DIM_STL))

      write_each_cell = .false.
      if(present(WRITE_FACETS_EACH_CELL)) then
         write_each_cell = WRITE_FACETS_EACH_CELL
      endif

      FACET_WRITTEN = .false.

      IF(NODESI*NODESJ*NODESK == 1) THEN
         WRITE(FN,'("GEOM_DG.stl")')
         WRITE(FN_PO,'("GEOM_DG_PO.stl")')
      ELSE
         WRITE(FN,'("GEOM_DG_",I5.5,".stl")') MYPE
         WRITE(FN_PO,'("GEOM_DG_PO_",I5.5,".stl")')  MYPE
      ENDIF

      IF(COUNT_FACET_TYPE_PO.GE.1) THEN
         OPEN(UNIT=443, FILE=TRIM(FN_PO))
         write(443,*)'solid vcg'
      endif

      OPEN(UNIT=444, FILE=trim(FN))
      write(444,*)'solid vcg'

      DO K=DG_KSTART2, DG_KEND2
      DO J=DG_JSTART2, DG_JEND2
      DO I=DG_ISTART2, DG_IEND2

         CELL_ID = DG_FUNIJK(I,J,K)
         COUNT_FACETS =  LIST_FACET_AT_DES(CELL_ID)%COUNT_FACETS

         IF(COUNT_FACETS.EQ.0) CYCLE

         IF(WRITE_EACH_CELL) THEN

            WRITE(FN,'("GEOM_DG",3("_",I3.3),"_",I8.8,".stl")')        &
               I, J, K, CELL_ID


            WRITE(FN_PO,'("GEOM_DG_PO",3("_",I3.3),"_",I8.8,".stl")') &
               I,J,K,CELL_ID

            OPEN(UNIT=446, FILE=FN)
            WRITE(446,*)'solid vcg'

            OPEN(UNIT=445, FILE=FN_PO)
            WRITE(445,*)'solid vcg'
         ENDIF

         DO COUNT = 1, COUNT_FACETS
            N = LIST_FACET_AT_DES(CELL_ID)%FACET_LIST(COUNT)

            IF(WRITE_EACH_CELL) THEN
               IF(STL_FACET_TYPE(N).EQ.FACET_TYPE_NORMAL) THEN
                  W_UNIT  = 446
               ELSE
                  W_UNIT  = 445
               ENDIF

               write(w_unit,*) '   facet normal ', NORM_FACE(:,N)
               write(w_unit,*) '      outer loop'
               write(w_unit,*) '         vertex ', VERTEX(1,1:3,N)
               write(w_unit,*) '         vertex ', VERTEX(2,1:3,N)
               write(w_unit,*) '         vertex ', VERTEX(3,1:3,N)
               write(w_unit,*) '      endloop'
               write(w_unit,*) '   endfacet'

            ENDIF

            IF (FACET_WRITTEN(N)) CYCLE

            IF(STL_FACET_TYPE(N).EQ.FACET_TYPE_NORMAL) THEN
               W_UNIT  = 444
            ELSE
               W_UNIT  = 443
            ENDIF

            write(w_unit,*) '   facet normal ', NORM_FACE(:,N)
            write(w_unit,*) '      outer loop'
            write(w_unit,*) '         vertex ', VERTEX(1,1:3,N)
            write(w_unit,*) '         vertex ', VERTEX(2,1:3,N)
            write(w_unit,*) '         vertex ', VERTEX(3,1:3,N)
            write(w_unit,*) '      endloop'
            write(w_unit,*) '   endfacet'
            facet_written(N) = .true.
         ENDDO

         if(write_each_cell) then
            write(445,*)'endsolid vcg'
            close(445)

            write(446,*)'endsolid vcg'
            close(446)
         ENDIF

         ENDDO
      ENDDO
      ENDDO
      write(444,*)'endsolid vcg'
      close(444)

      if(count_facet_type_po.gt.0) then
         write(443,*)'endsolid vcg'
         close(443)
      ENDIF

      DEALLOCATE (FACET_WRITTEN)

!      IF(MyPE == PE_IO) THEN
!         WRITE(*,*) ' The file geometry_from_grid_facets.stl was sucessfully written.'
!         WRITE(*,*) ' This is the based on facets stored grid wise for DES modules'
!         WRITE(*,*) ' and is provided for convenience and debugging (it is not used).'
!      ENDIF
      END Subroutine  DEBUG_write_stl_from_grid_facet


      END MODULE STL_PREPROC_DES


