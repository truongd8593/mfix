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

      Use discretelement, only: MINIMIZE_DES_FACET_LIST

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

! A separating axis exists so the cell and the triangle do not intersect.
! JFD: If MINIMIZE_DES_FACET_LIST is .FALSE. then the facet is still added to the list.
!      This seems to prevent particle leakage for some des gris setups.
      IF(SA_EXIST.AND.MINIMIZE_DES_FACET_LIST ) RETURN

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

         open(stl_unit, file = fname, form='formatted',convert='big_endian')
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

      OPEN(1001, file = TRIM(FN), form ='formatted',CONVERT='BIG_ENDIAN')

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

      OPEN(UNIT=444, FILE='FACETS_TO_DG.stl',CONVERT='BIG_ENDIAN')
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
         OPEN(UNIT=443, FILE=TRIM(FN_PO),CONVERT='BIG_ENDIAN')
         write(443,*)'solid vcg'
      endif

      OPEN(UNIT=444, FILE=trim(FN),CONVERT='BIG_ENDIAN')
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

            OPEN(UNIT=446, FILE=FN,CONVERT='BIG_ENDIAN')
            WRITE(446,*)'solid vcg'

            OPEN(UNIT=445, FILE=FN_PO,CONVERT='BIG_ENDIAN')
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



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: ADD_FACET                                               !
!                                                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE ADD_FACET(CELL_ID, FACET_ID)

      Use discretelement
      USE stl

      implicit none

      INTEGER, INTENT(IN) :: cell_id, facet_id

      INTEGER, DIMENSION(:), ALLOCATABLE :: int_tmp
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: real_tmp

      INTEGER :: lSIZE2, ii
      DOUBLE PRECISION :: smallest_extent, min_temp, max_temp

      IF(STL_FACET_TYPE(facet_id) /= FACET_TYPE_NORMAL) RETURN

      DO II = 1, CELLNEIGHBOR_FACET_NUM(CELL_ID)
         IF(FACET_ID .EQ. CELLNEIGHBOR_FACET(CELL_ID)%P(II)) RETURN
      ENDDO

      CELLNEIGHBOR_FACET_NUM(CELL_ID) = &
         CELLNEIGHBOR_FACET_NUM(CELL_ID) + 1

      NO_NEIGHBORING_FACET_DES(CELL_ID)  = .FALSE.

      IF(cellneighbor_facet_num(cell_id) > &
         cellneighbor_facet_max(cell_id)) THEN

         cellneighbor_facet_max(cell_id) = &
         2*cellneighbor_facet_max(cell_id)

         lSIZE2 = size(cellneighbor_facet(cell_id)%p)
         allocate(int_tmp(cellneighbor_facet_max(cell_id)))
         int_tmp(1:lSIZE2) = cellneighbor_facet(cell_id)%p(1:lSIZE2)
         call move_alloc(int_tmp,cellneighbor_facet(cell_id)%p)

         lSIZE2 = size(cellneighbor_facet(cell_id)%extentdir)
         allocate(int_tmp(cellneighbor_facet_max(cell_id)))
         int_tmp(1:lSIZE2) = &
            cellneighbor_facet(cell_id)%extentdir(1:lSIZE2)
         call move_alloc(int_tmp,cellneighbor_facet(cell_id)%extentdir)

         lSIZE2 = size(cellneighbor_facet(cell_id)%extentmin)
         allocate(real_tmp(cellneighbor_facet_max(cell_id)))
         real_tmp(1:lSIZE2) = &
            cellneighbor_facet(cell_id)%extentmin(1:lSIZE2)
         call move_alloc(real_tmp,cellneighbor_facet(cell_id)%extentmin)

         lSIZE2 = size(cellneighbor_facet(cell_id)%extentmax)
         allocate(real_tmp(cellneighbor_facet_max(cell_id)))
         real_tmp(1:lSIZE2) = &
            cellneighbor_facet(cell_id)%extentmax(1:lSIZE2)
         call move_alloc(real_tmp,cellneighbor_facet(cell_id)%extentmax)

      ENDIF

      CELLNEIGHBOR_FACET(CELL_ID)%&
         P(CELLNEIGHBOR_FACET_NUM(CELL_ID)) = FACET_ID
      SMALLEST_EXTENT = HUGE(0.0)

      DO II=1,3
         MIN_TEMP = MINVAL(VERTEX(:,II,FACET_ID))
         MAX_TEMP = MAXVAL(VERTEX(:,II,FACET_ID))

         IF(ABS(MAX_TEMP - MIN_TEMP) < SMALLEST_EXTENT ) THEN
             CELLNEIGHBOR_FACET(CELL_ID)%&
                EXTENTDIR(CELLNEIGHBOR_FACET_NUM(CELL_ID)) = II
             CELLNEIGHBOR_FACET(CELL_ID)%&
                EXTENTMIN(CELLNEIGHBOR_FACET_NUM(CELL_ID)) = MIN_TEMP
             CELLNEIGHBOR_FACET(CELL_ID)%&
                EXTENTMAX(CELLNEIGHBOR_FACET_NUM(CELL_ID)) = MAX_TEMP
             SMALLEST_EXTENT = ABS(MAX_TEMP - MIN_TEMP)
         ENDIF
      ENDDO

      END SUBROUTINE ADD_FACET

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  SUBROUTINE: CHECK_IF_PARTICLE_OVERLAPS_STL                           C
!                                                                      C
!  Purpose: This subroutine is special written to check if a particle  C
!          overlaps any of the STL faces. The routine exits on         C
!          detecting an overlap. It is called after initial            C
!          generation of lattice configuration to remove out of domain C
!          particles                                                   C
!                                                                      C
!  Authors: Rahul Garg                               Date: 21-Mar-2014 C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CHECK_IF_PARTICLE_OVERLAPS_STL(POSITION, RADIUS, &
         OVERLAP_EXISTS)

      USE run
      USE param1
      USE discretelement, only: dimn, xe, yn, zt
      USE geometry
      USE constant
      USE cutcell
      USE indices
      USE stl
      USE compar
      USE des_stl_functions
      use desgrid
      USE functions
      Implicit none

      DOUBLE PRECISION, INTENT(IN) :: POSITION(DIMN), RADIUS
      LOGICAL, INTENT(OUT) :: OVERLAP_EXISTS

      INTEGER I, J, K, IJK, NF

      DOUBLE PRECISION :: RADSQ, DISTSQ, DIST(DIMN), CLOSEST_PT(DIMN)
      INTEGER :: COUNT_FAC, COUNT, contact_facet_count, NEIGH_CELLS, &
      NEIGH_CELLS_NONNAT, &
      LIST_OF_CELLS(27), CELL_ID, I_CELL, J_CELL, K_CELL, cell_count , &
      IMINUS1, IPLUS1, JMINUS1, JPLUS1, KMINUS1, KPLUS1, PHASELL, LOC_MIN_PIP, &
      LOC_MAX_PIP!, focus_particle


      OVERLAP_EXISTS = .FALSE.

      I_CELL = MIN(DG_IEND2,MAX(DG_ISTART2,IOFPOS(POSITION(1))))
      J_CELL = MIN(DG_JEND2,MAX(DG_JSTART2,JOFPOS(POSITION(2))))
      K_CELL = 1
      IF(DO_K) K_CELL =MIN(DG_KEND2,MAX(DG_KSTART2,KOFPOS(POSITION(3))))

      CELL_ID = DG_FUNIJK(I_CELL, J_CELL, K_CELL)
      IF (NO_NEIGHBORING_FACET_DES(CELL_ID)) RETURN

      LIST_OF_CELLS(:) = -1
      NEIGH_CELLS = 0
      NEIGH_CELLS_NONNAT  = 0

      COUNT_FAC = LIST_FACET_AT_DES(CELL_ID)%COUNT_FACETS

      RADSQ = RADIUS*RADIUS

! Add the facets in the cell the particle currently resides in
      IF (COUNT_FAC.gt.0)   then
         NEIGH_CELLS = NEIGH_CELLS + 1
         LIST_OF_CELLS(NEIGH_CELLS) = CELL_ID
      ENDIF

      IPLUS1  =  MIN (I_CELL + 1, DG_IEND2)
      IMINUS1 =  MAX (I_CELL - 1, DG_ISTART2)

      JPLUS1  =  MIN (J_CELL + 1, DG_JEND2)
      JMINUS1 =  MAX (J_CELL - 1, DG_JSTART2)

      KPLUS1  =  MIN (K_CELL + 1, DG_KEND2)
      KMINUS1 =  MAX (K_CELL - 1, DG_KSTART2)

      DO K = KMINUS1, KPLUS1
         DO J = JMINUS1, JPLUS1
            DO I = IMINUS1, IPLUS1

               IF(.NOT.dg_is_ON_myPE_plus1layers(I,J,K)) CYCLE

               IJK = DG_FUNIJK(I,J,K)
               COUNT_FAC = LIST_FACET_AT_DES(IJK)%COUNT_FACETS

               IF(COUNT_FAC.EQ.0) CYCLE
               distsq = zero

               IF(POSITION(1) > XE(I)) DISTSQ = DISTSQ &
                  + (POSITION(1)-XE(I))*(POSITION(1)-XE(I))

               IF(POSITION(1) < XE(I) - DX(I)) DISTSQ = DISTSQ &
               + (XE(I) - DX(I) - POSITION(1))*(XE(I) - DX(I) - POSITION(1))

               IF(POSITION(2) > YN(J)) DISTSQ = DISTSQ &
               + (POSITION(2)-YN(J))* (POSITION(2)-YN(J))

               IF(POSITION(2) < YN(J) - DY(J)) DISTSQ = DISTSQ &
               + (YN(J) - DY(J) - POSITION(2))* (YN(J) - DY(J) - POSITION(2))

               IF(POSITION( 3) > ZT(K)) DISTSQ = DISTSQ &
               + (POSITION(3)-ZT(K))*(POSITION(3)-ZT(K))

               IF(POSITION(3) < ZT(K) - DZ(K)) DISTSQ = DISTSQ &
               + (ZT(K) - DZ(K) - POSITION(3))*(ZT(K) - DZ(K) - POSITION(3))
               IF (DISTSQ < RADSQ) then
                  NEIGH_CELLS_NONNAT = NEIGH_CELLS_NONNAT + 1
                  NEIGH_CELLS = NEIGH_CELLS + 1
                  LIST_OF_CELLS(NEIGH_CELLS) = IJK
                                !WRITE(*,'(A10, 4(2x,i5))') 'WCELL  = ', IJK, I,J,K
               ENDIF
            ENDDO
         ENDDO
      ENDDO

      CONTACT_FACET_COUNT = 0

      DO CELL_COUNT = 1, NEIGH_CELLS
         IJK = LIST_OF_CELLS(CELL_COUNT)

         DO COUNT = 1, LIST_FACET_AT_DES(IJK)%COUNT_FACETS
            NF = LIST_FACET_AT_DES(IJK)%FACET_LIST(COUNT)

            CALL ClosestPtPointTriangle(POSITION(:), VERTEX(:,:,NF), CLOSEST_PT(:))

            DIST(:) = POSITION(:) - CLOSEST_PT(:)
            DISTSQ = DOT_PRODUCT(DIST, DIST)

            IF(DISTSQ .GE. RADSQ) CYCLE !No overlap exists, move on to the next facet

            !Overlap detected
            !Set overlap_exists to true and exit
            OVERLAP_EXISTS = .true.
            RETURN

         ENDDO

      end DO

      RETURN

    END SUBROUTINE CHECK_IF_PARTICLE_OVERLAPS_STL


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!                                                                      !
!                                                                      !
!                                                                      !
!                                                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE write_this_facet_and_part(FID, PID)
      USE run
      USE param1
      USE discretelement
      USE geometry
      USE compar
      USE constant
      USE cutcell
      USE funits
      USE indices
      USE physprop
      USE parallel
      USE stl
      USE des_stl_functions
      Implicit none
      !facet id and particle id
      Integer, intent(in) :: fid, pid
      Integer :: stl_unit, vtp_unit , k
      CHARACTER(LEN=100) :: stl_fname, vtp_fname
      real :: temp_array(3)

      stl_unit = 1001
      vtp_unit = 1002

      WRITE(vtp_fname,'(A,"_OFFENDING_PARTICLE",".vtp")') TRIM(RUN_NAME)
      WRITE(stl_fname,'(A,"_STL_FACE",".stl")') TRIM(RUN_NAME)

      open(vtp_unit, file = vtp_fname, form='formatted',convert='big_endian')
      open(stl_unit, file = stl_fname, form='formatted',convert='big_endian')

      write(vtp_unit,"(a)") '<?xml version="1.0"?>'
      write(vtp_unit,"(a,es24.16,a)") '<!-- time =',s_time,'s -->'
      write(vtp_unit,"(a,a)") '<VTKFile type="PolyData"',&
           ' version="0.1" byte_order="LittleEndian">'
      write(vtp_unit,"(3x,a)") '<PolyData>'
      write(vtp_unit,"(6x,a,i10.10,a,a)")&
           '<Piece NumberOfPoints="',1,'" NumberOfVerts="0" ',&
           'NumberOfLines="0" NumberOfStrips="0" NumberOfPolys="0">'
      write(vtp_unit,"(9x,a)")&
           '<PointData Scalars="Diameter" Vectors="Velocity">'
      write(vtp_unit,"(12x,a)")&
           '<DataArray type="Float32" Name="Diameter" format="ascii">'
      write (vtp_unit,"(15x,es13.6)") (2*des_radius(pid))
      write(vtp_unit,"(12x,a)") '</DataArray>'

      temp_array = zero
      temp_array(:) = des_vel_new(:,pid)
      write(vtp_unit,"(12x,a,a)") '<DataArray type="Float32" ',&
           'Name="Velocity" NumberOfComponents="3" format="ascii">'
      write (vtp_unit,"(15x,3(es13.6,3x))")&
           ((temp_array(k)),k=1,3)
      write(vtp_unit,"(12x,a,/9x,a)") '</DataArray>','</PointData>'
      ! skip cell data
      write(vtp_unit,"(9x,a)") '<CellData></CellData>'

      temp_array = zero
      temp_array(1:dimn) = des_pos_new(1:dimn, pid)
      write(vtp_unit,"(9x,a)") '<Points>'
      write(vtp_unit,"(12x,a,a)") '<DataArray type="Float32" ',&
           'Name="Position" NumberOfComponents="3" format="ascii">'
      write (vtp_unit,"(15x,3(es13.6,3x))")&
           ((temp_array(k)),k=1,3)
      write(vtp_unit,"(12x,a,/9x,a)")'</DataArray>','</Points>'
      ! Write tags for data not included (vtp format style)
      write(vtp_unit,"(9x,a,/9x,a,/9x,a,/9x,a)")'<Verts></Verts>',&
           '<Lines></Lines>','<Strips></Strips>','<Polys></Polys>'
      write(vtp_unit,"(6x,a,/3x,a,/a)")&
           '</Piece>','</PolyData>','</VTKFile>'

      !Now write the facet info

      write(stl_unit,*)'solid vcg'

      write(stl_unit,*) '   facet normal ', NORM_FACE(1:3,FID)
      write(stl_unit,*) '      outer loop'
      write(stl_unit,*) '         vertex ', VERTEX(1,1:3,FID)
      write(stl_unit,*) '         vertex ', VERTEX(2,1:3,FID)
      write(stl_unit,*) '         vertex ', VERTEX(3,1:3,FID)
      write(stl_unit,*) '      endloop'
      write(stl_unit,*) '   endfacet'

      write(stl_unit,*)'endsolid vcg'

      close(vtp_unit, status = 'keep')
      close(stl_unit, status = 'keep')

      end SUBROUTINE write_this_facet_and_part

      END MODULE STL_PREPROC_DES


