!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SET_BC_DEM_MO                                           !
!                                                                      !
!                                                                      !
!  Author: J.Musser                                   Date: 23-Nov-09  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE SET_BC_DEM_MO

      use bc, only: BC_PLANE
      use bc, only: BC_I_w, BC_I_e, BC_J_s, BC_J_n, BC_K_b, BC_K_t

      use des_bc, only: DEM_BCMO, DEM_BCMO_MAP, DEM_BCMO_IJK
      use des_bc, only: DEM_BCMO_IJKSTART, DEM_BCMO_IJKEND

      use funits, only: DMP_LOG

      use mpi_utility
      use error_manager

      IMPLICIT NONE

!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER :: BCV, BCV_I      ! BC loop counter

      INTEGER :: LC

      LOGICAL, parameter :: setDBG = .FALSE.
      LOGICAL :: dFlag

      INTEGER :: MAX_CELLS, BND1, BND2

      INTEGER, ALLOCATABLE :: LOC_DEM_BCMO_IJK(:)

      INTEGER :: I,J,K,IJK
      INTEGER :: I_w, I_e, J_s, J_n, K_b, K_t


      include 'function.inc'


      CALL INIT_ERR_MSG("SET_BC_DEM_MO")

      dFlag = (DMP_LOG .AND. setDBG)
      if(dFlag) write(*,"(2/,2x,'DEM outlet count: ',I4)") DEM_BCMO

! Loop over the outflow BCs to get an approximate count of the number
! of fluid cells that are adjacent to the outlet.
      MAX_CELLS = 0
      DO BCV_I = 1, DEM_BCMO
         BCV = DEM_BCMO_MAP(BCV_I)

! Set the search area to the dimensions of the inlet.
         if(dFlag) WRITE(*,"(/2x,'Adding cells for BCV: ',I3)") BCV
         SELECT CASE (BC_PLANE(BCV))
         CASE('N','S')
            BND1 = BC_I_e(BCV) - BC_I_w(BCV)
            BND2 = BC_K_t(BCV) - BC_K_b(BCV)

         CASE('E','W')
            BND1 = BC_J_n(BCV) - BC_J_s(BCV)
            BND2 = BC_K_t(BCV) - BC_K_b(BCV)

         CASE('T','B')
            BND1 = BC_I_e(BCV) - BC_I_w(BCV)
            BND2 = BC_J_n(BCV) - BC_J_s(BCV)
         END SELECT

         MAX_CELLS = MAX_CELLS +                                      &
            2*(BND1+1)*(BND2+1) + 2*(BND1+2) + 2*(BND2+2)

         if(dFlag) WRITE(*,"(4x,'Plane:   ',A)") BC_PLANE(BCV)
         if(dFlag) WRITE(*,"(4x,'Cells: ', I8)") (BND1+1)*(BND2+1)
      ENDDO

      if(dFlag) write(*,"(2x,'Max Cells: ',I8)") MAX_CELLS

! Allocate an array to hold the IJK values. This array should be
! more than enough to store all the IJKs.
      allocate( LOC_DEM_BCMO_IJK(MAX_CELLS) )

! Loop over the IJKs for each BC and store only the IJKs that you
! own as well as the start/end array positions for each BC.
      LC = 1
      DO BCV_I = 1, DEM_BCMO

         DEM_BCMO_IJKSTART(BCV_I) = LC
         BCV = DEM_BCMO_MAP(BCV_I)

         if(dFlag) write(*,"(/2x,'Searching for fluid cells:',I3)") BCV

         I_w = BC_I_w(BCV); I_e = BC_I_e(BCV)
         J_s = BC_J_s(BCV); J_n = BC_J_n(BCV)
         K_b = BC_K_b(BCV); K_t = BC_K_t(BCV)

! Depending on the flow plane, the 'common' index needs shifted to
! reference the fluid cell.
         SELECT CASE (BC_PLANE(BCV))
         CASE('N'); J_s = J_s+1;  J_n = J_s
         CASE('S'); J_s = J_s-1;  J_n = J_s
         CASE('E'); I_w = I_w+1;  I_e = I_w
         CASE('W'); I_w = I_w-1;  I_e = I_w
         CASE('T'); K_b = K_b+1;  K_t = K_b
         CASE('B'); K_b = K_b-1;  K_t = K_b
        END SELECT

         if(dFlag) then
            write(*,"(4x,'Search bounds: ')")
            write(*,"(6x,'I_w/I_e:',2(2x,I6))") I_w, I_e
            write(*,"(6x,'J_s/J_n:',2(2x,I6))") J_s, J_n
            write(*,"(6x,'K_b/K_t:',2(2x,I6))") K_b, K_t
         endif

! Store the IJKs.
         DO K = K_b, K_t
         DO J = J_s, J_n
         DO I = I_w, I_e
! Skip cells that this rank does not own or are considered dead.
            IF(.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
            IF(DEAD_CELL_AT(I,J,K)) CYCLE

            IJK = FUNIJK(I,J,K)
            LOC_DEM_BCMO_IJK(LC) = IJK
            LC = LC+1
         ENDDO
         ENDDO
         ENDDO

         if(dFlag) write(*,"(/2x,'Adding boundary cells:',I3)") BCV

         I_w = BC_I_w(BCV)-1; I_e = BC_I_e(BCV)+1
         J_s = BC_J_s(BCV)-1; J_n = BC_J_n(BCV)+1

         IF(DO_K) THEN
            K_b = BC_K_b(BCV)-1; K_t = BC_K_t(BCV)+1
         ELSE
            K_b = BC_K_b(BCV);   K_t = BC_K_t(BCV)
         ENDIF

! Depending on the flow plane, the 'common' index needs shifted to
! reference the fluid cell.
         SELECT CASE (BC_PLANE(BCV))
         CASE('N','S'); J_s = J_s+1;  J_n = J_n-1
         CASE('E','W'); I_w = I_w+1;  I_e = I_e-1
         CASE('T','B'); K_b = K_b+1;  K_t = K_t-1
         END SELECT

         if(dFlag) then
            write(*,"(4x,'Search bounds: ')")
            write(*,"(6x,'I_w/I_e:',2(2x,I6))") I_w, I_e
            write(*,"(6x,'J_s/J_n:',2(2x,I6))") J_s, J_n
            write(*,"(6x,'K_b/K_t:',2(2x,I6))") K_b, K_t
         endif

! Store the IJKs.
         DO K = K_b, K_t
         DO J = J_s, J_n
         DO I = I_w, I_e
! Skip cells that this rank does not own or are considered dead.
            IF(.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
            IF(DEAD_CELL_AT(I,J,K)) CYCLE

            IJK = FUNIJK(I,J,K)
            LOC_DEM_BCMO_IJK(LC) = IJK
            LC = LC+1
         ENDDO
         ENDDO
         ENDDO

         DEM_BCMO_IJKEND(BCV_I) = LC-1

         if(dFLAG) write(*,1111) BCV, BCV_I,                           &
            DEM_BCMO_IJKSTART(BCV_I),DEM_BCMO_IJKEND(BCV_I)

      ENDDO

 1111 FORMAT(/2x,'DEM Mass Outflow:',/4x,'BC:',I4,3x,'MAP:',I4,&
         /4x,'IJKSTART:',I6,/4x,'IJKEND:  ',I6)

! Allocate the global store arrary array. This changes across MPI ranks.
      IF(LC >=1) THEN
         allocate( DEM_BCMO_IJK(LC-1) )
         DEM_BCMO_IJK(1:LC-1) = LOC_DEM_BCMO_IJK(1:LC-1)
      ELSE
         allocate( DEM_BCMO_IJK(1) )
         DEM_BCMO_IJK(1) = -1
      ENDIF

      deallocate(LOC_DEM_BCMO_IJK)

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE SET_BC_DEM_MO


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: CHECK_DES_LE_BC                                        !
!                                                                      !
!  Purpose: Check/set parameters for DES Lees Edeards BC.              !
!                                                                      !
!  Author: J.Musser                                   Date: 11-DEC-13  !
!                                                                      !
!  Comments: *** DES Lees Edwards BC funcionality has been lost. ***   !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_DES_LE_BC

      use discretelement
      use mpi_utility

      IMPLICIT NONE

! Lees Edwards BC functionality has been lost in current DEM code
      IF(DES_LE_BC) THEN
         IF (DES_CONTINUUM_COUPLED) THEN
            WRITE(UNIT_LOG, 1064)
             CALL MFIX_EXIT(myPE)
         ENDIF
         IF (DES_NEIGHBOR_SEARCH .NE. 4) THEN
            WRITE(UNIT_LOG, 1060)
            CALL MFIX_EXIT(myPE)
         ENDIF
! not all possible shear directions are fully coded         
         IF (DIMN .EQ. 2) THEN
            IF(TRIM(DES_LE_SHEAR_DIR) .NE. 'DUDY' .AND. &
               TRIM(DES_LE_SHEAR_DIR) .NE. 'DVDX') THEN
               WRITE(UNIT_LOG, 1061)
               CALL MFIX_EXIT(myPE)
            ENDIF
         ELSEIF(DIMN.EQ.3) THEN
            IF(TRIM(DES_LE_SHEAR_DIR) .NE. 'DUDY') THEN ! .AND. & 
!               TRIM(DES_LE_SHEAR_DIR) .NE. 'DUDZ' .AND. &
!               TRIM(DES_LE_SHEAR_DIR) .NE. 'DVDX' .AND. &
!               TRIM(DES_LE_SHEAR_DIR) .NE. 'DVDZ' .AND. &
!               TRIM(DES_LE_SHEAR_DIR) .NE. 'DWDX' .AND. &
!               TRIM(DES_LE_SHEAR_DIR) .NE. 'DWDY') THEN
               WRITE(UNIT_LOG, 1062)
               CALL MFIX_EXIT(myPE)
            ENDIF
         ENDIF
         IF (DES_PERIODIC_WALLS) THEN
            DES_PERIODIC_WALLS = .FALSE.
            DES_PERIODIC_WALLS_X = .FALSE.
            DES_PERIODIC_WALLS_Y = .FALSE.
            DES_PERIODIC_WALLS_Z = .FALSE.
            WRITE(UNIT_LOG, 1063)
            WRITE(*,1063)
         ENDIF
      ENDIF

      RETURN

 1060 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'Only the grid based search option is allowed when using',&
         'using',/10X,'Lees & Edwards BC.',/1X,70('*')/)

 1061 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'Invalid option for DES_LE_SHEAR_DIR in mfix.dat. When',/10X,&
         'DIMN=2 shear options are DUDY or DVDX',/1X,70('*')/)

 1062 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'Invalid option for DES_LE_SHEAR_DIR in mfix.dat. When',/10X,&
         'DIMN=3 shear options are DUDY, DUDZ, DVDX, DVDZ, DWDX or',&
         'DWDY.',/1X,70('*')/)

 1063 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'WARNING: DES_PERIODIC_WALLS set to false when DES_LE_BC.',&
         /10X,'DES_LE_BC implies periodic walls, however, the ',&
         'periodicity is handled',/10X, 'independently of ',&
         'DES_PERIODIC_WALLS option and so it is shut off.',&
         /1X,70('*')/)

 1064 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'DES_CONTINUUM_COUPLED cannot be true when using ',&
         'DES_LE_BC.',/1X,70('*')/)

      END SUBROUTINE CHECK_DES_LE_BC
