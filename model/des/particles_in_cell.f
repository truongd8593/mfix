!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Subroutine: PARTICLES_IN_CELL                                       !
!                                                                      !
!  Purpose:                                                            !
!     - For each particle find the computational fluid cell            !
!       containing the particle center.                                !
!     - Calculate the bulk density in each computational fluid         !
!       cell.                                                          !
!     - Calculate the volume average solids velocity in each           !
!       computational fluid cell.                                      !
!     - For parallel processing indices are altered                    !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE PARTICLES_IN_CELL

      use tmp_array, only: PARTICLE_COUNT => ARRAY1

      USE param1
      USE fldvar
      USE geometry
      USE indices
      USE physprop
      USE compar
      USE parallel
      USE sendrecv
      USE discretelement
      use desgrid
      use desmpi
      USE cutcell
      USE mfix_pic
      USE des_rxns
      USE run
      USE error_manager
      USE functions

! Number of particles in the I/J/K direction
      use param, only: DIMENSION_I, DIMENSION_J, DIMENSION_K


      IMPLICIT NONE
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! particle no.
      INTEGER L
! accounted for particles
      INTEGER PC
! solids phase no.
      INTEGER M
! ijk indices
      INTEGER I, J, K, IJK
! variables that count/store the number of particles in i, j, k cell
      INTEGER:: npic, pos

      INTEGER :: RECOVERED
      INTEGER :: DELETED

! following quantities are reset every call to particles_in_cell
      PINC(:) = 0

! Use an incremental approach to determine the new particle location.
!-----------------------------------------------------------------------
!!$omp parallel default(shared) private(L, I, J, K, IJK)
!!$omp do reduction(+:PINC) schedule (guided,50)

      DO L = 1, MAX_PIP
! skipping particles that do not exist
         IF(.NOT.PEA(L,1)) CYCLE
! skipping ghost particles
         IF(PEA(L,4)) CYCLE


         I = PIJK(L,1)
         IF(I <= ISTART2 .OR. I >= IEND2) THEN
            CALL PIC_SEARCH(I, DES_POS_NEW(1,L), XE,                   &
               DIMENSION_I, ISTART2, IEND2)
         ELSE
            IF((DES_POS_NEW(1,L) >= XE(I-1)) .AND.                     &
               (DES_POS_NEW(1,L) <  XE(I))) THEN
               I = I
            ELSEIF((DES_POS_NEW(1,L) >= XE(I)) .AND.                   &
               (DES_POS_NEW(1,L) < XE(I+1))) THEN
              I = I+1
            ELSEIF((DES_POS_NEW(1,L) >= XE(I-2)) .AND.                 &
               (DES_POS_NEW(1,L) < XE(I-1))) THEN
               I = I-1
            ELSE
               CALL PIC_SEARCH(I, DES_POS_NEW(1,L), XE,                &
                  DIMENSION_I, ISTART2, IEND2)
            ENDIF
         ENDIF


         J = PIJK(L,2)
         IF(J <= JSTART2 .OR. J >= JEND2) THEN
            CALL PIC_SEARCH(J, DES_POS_NEW(2,L), YN,                   &
               DIMENSION_J, JSTART2, JEND2)
         ELSE
            IF((DES_POS_NEW(2,L) >= YN(J-1)) .AND.                     &
               (DES_POS_NEW(2,L) < YN(J))) THEN
               J = J
            ELSEIF((DES_POS_NEW(2,L) >= YN(J)) .AND.                   &
               (DES_POS_NEW(2,L) < YN(J+1))) THEN
               J = J+1
            ELSEIF((DES_POS_NEW(2,L) >= YN(J-2)) .AND.                 &
               (DES_POS_NEW(2,L) < YN(J-1)))THEN
               J = J-1
            ELSE
               CALL PIC_SEARCH(J, DES_POS_NEW(2,L), YN,                &
                  DIMENSION_J, JSTART2, JEND2)
            ENDIF
         ENDIF


         IF(NO_K) THEN
            K = 1
         ELSE
            K = PIJK(L,3)
            IF(K <= KSTART2 .OR. K >= KEND2) THEN
               CALL PIC_SEARCH(K, DES_POS_NEW(3,L), ZT,                &
                  DIMENSION_K, KSTART2, KEND2)
            ELSE
               IF((DES_POS_NEW(3,L) >= ZT(K-1)) .AND.                  &
                  (DES_POS_NEW(3,L) < ZT(K))) THEN
                  K = K
                ELSEIF((DES_POS_NEW(3,L) >= ZT(K)) .AND.               &
                  (DES_POS_NEW(3,L) < ZT(K+1))) THEN
                  K = K+1
               ELSEIF((DES_POS_NEW(3,L) >= ZT(K-2)) .AND.              &
                  (DES_POS_NEW(3,L) >= ZT(K-1))) THEN
                  K = K-1
               ELSE
                  CALL PIC_SEARCH(K, DES_POS_NEW(3,L), ZT,             &
                     DIMENSION_K, KSTART2, KEND2)
               ENDIF
            ENDIF
         ENDIF

! Calculate the fluid cell index.
         IJK = FUNIJK(I,J,K)
! Increment the number of particles in cell IJK
         PINC(IJK) = PINC(IJK) + 1
! Assign PIJK(L,1:4)
         PIJK(L,1) = I
         PIJK(L,2) = J
         PIJK(L,3) = K
         PIJK(L,4) = IJK


      ENDDO
!!$omp end parallel

      CALL CHECK_CELL_MOVEMENT(RECOVERED, DELETED)

! Assigning the variable PIC(IJK)%p(:). For each computational fluid
! cell compare the number of current particles in the cell to what was
! in the cell previously. If different reallocate. Store the particle
! ids
! ---------------------------------------------------------------->>>
!!$omp parallel do if(ijkend3 .ge. 2000) default(shared)           &
!!$omp private(ijk,npic) !schedule (guided,50)
      DO IJK = IJKSTART3, IJKEND3

! checking all cells (including ghost cells); updating entering/exiting
! particle regions
         NPIC =  PINC(IJK)
         IF (ASSOCIATED(PIC(IJK)%p)) THEN
            IF (NPIC.NE.SIZE(PIC(IJK)%p)) THEN
               DEALLOCATE(PIC(IJK)%p)
               IF (NPIC.GT.0) ALLOCATE(PIC(IJK)%p(NPIC))
            ENDIF
         ELSE
            IF (NPIC.GT.0) ALLOCATE(PIC(IJK)%p(NPIC))
         ENDIF
      ENDDO
!!$omp end parallel do


      particle_count(:) = 1
      PC = 1
      DO L = 1, MAX_PIP
! exiting loop if reached max number of particles in processor
         IF(PC.GT.PIP) exit
! skipping indices with no particles (non-existent particles)
         IF(.NOT.PEA(L,1)) CYCLE
! incrementing particle account when particle exists
         pc = pc+1
! skipping ghost particles
         IF(PEA(L,4)) CYCLE
         IJK = PIJK(L,4)
         pos = particle_count(IJK)
         pic(IJK)%p(pos) = L
         particle_count(IJK) = particle_count(IJK) + 1
      ENDDO

! Calculate interpolation weights
!      CALL CALC_INTERP_WEIGHTS

! Calculate mean fields using either interpolation or cell averaging.
!      CALL COMP_MEAN_FIELDS

!      IF(MPPIC) CALL REPORT_PIC_STATS(RECOVERED, DELETED)

      RETURN
      END SUBROUTINE PARTICLES_IN_CELL


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: INIT_PARTICLES_IN_CELL                                  !
!                                                                      !
!  Purpose:                                                            !
!     - For each particle find the computational fluid cell            !
!       containing the particle center.                                !
!     - Calculate the bulk density in each computational fluid         !
!       cell.                                                          !
!     - Calculate the volume average solids velocity in each           !
!       computational fluid cell.                                      !
!     - For parallel processing indices are altered                    !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE INIT_PARTICLES_IN_CELL

      USE physprop, only: SMAX

      use discretelement, only: PEA, PIJK, PINC
      USE discretelement, only: DES_POS_NEW
      USE discretelement, only: MAX_PIP
      USE discretelement, only: XE, YN, ZT
      use desmpi, only: des_par_exchange

      USE run, only: RUN_TYPE
      USE run, only: ANY_SPECIES_EQ

! Number of particles in the I/J/K direction
      use param, only: DIMENSION_I, DIMENSION_J, DIMENSION_K

      use mpi_utility
      use sendrecv

      USE error_manager
      USE functions

      IMPLICIT NONE
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! particle no.
      INTEGER :: L
! solids phase no.
      INTEGER :: M
! ijk indices
      INTEGER :: I, J, K, IJK
! particle x,y,z position
      DOUBLE PRECISION :: lPOS
! IER for error reporting
      INTEGER :: IER

      CALL INIT_ERR_MSG("INIT_PARTICLES_IN_CELL")

! following quantities are reset every call to particles_in_cell
      PINC(:) = 0

! Call exchange particles - this will exchange particle crossing
! boundaries as well as updates ghost particles information
      CALL DES_PAR_EXCHANGE

! Assigning PIJK(L,1), PIJK(L,2) and PIJK(L,3) the i, j, k indices
! of particle L (locating particle on fluid grid). Also determine
! composite ijk index. If first_pass, also assigning PIJK(L,5) the
! solids phase index of particle.
! ---------------------------------------------------------------->>>
!!$omp parallel default(shared)                                       &
!!$omp private(L,M,lPOS,I,J,K,IJK)
!!$omp do reduction(+:PINC) schedule (guided,50)
      DO L = 1, MAX_PIP
! skipping particles that do not exist
         IF(.NOT.PEA(L,1)) CYCLE
! skipping ghost particles
         IF(PEA(L,4)) CYCLE

! Use a brute force technique to determine the particle locations in
! the Eulerian fluid grid.

         CALL PIC_SEARCH(I, DES_POS_NEW(1,L), XE,                      &
            DIMENSION_I, ISTART2, IEND2)
         PIJK(L,1) = I

         CALL PIC_SEARCH(J, DES_POS_NEW(2,L), YN,                      &
            DIMENSION_J, JSTART2, JEND2)
         PIJK(L,2) = J

         IF(NO_K) THEN
            K=1
            PIJK(L,3) = 1
         ELSE
            CALL PIC_SEARCH(K, DES_POS_NEW(3,L), ZT,                   &
               DIMENSION_K, KSTART2, KEND2)
            PIJK(L,3) = K
         ENDIF

! Assigning PIJK(L,4) now that particles have been located on the fluid
         IJK = FUNIJK(I,J,K)
         PIJK(L,4) = IJK

         PINC(IJK) = PINC(IJK) + 1
      ENDDO
!!$omp end parallel
! Calling exchange particles - this will exchange particle crossing
! boundaries as well as updates ghost particles information
! unclear why this needs to be called again.
      CALL DES_PAR_EXCHANGE

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE INIT_PARTICLES_IN_CELL



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: PIC_SEARCH                                              !
!                                                                      !
!  Purpose: Identify the I (or J or K) index of the fluid cell that    !
!  contains the particle centroid.                                     !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE PIC_SEARCH(IDX, lPOS, ENT_POS, lDIMN, lSTART, lEND)

      IMPLICIT NONE
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! Index being searched for (I, J, or K)
      INTEGER, INTENT(OUT) :: IDX
! Particle x,y,z position
      DOUBLE PRECISION, INTENT(IN) :: lPOS
! Dimension of ENT_POS array
      INTEGER, INTENT(IN) :: lDIMN
! East, North, or Top cell face location
      DOUBLE PRECISION, INTENT(IN) :: ENT_POS(0:lDIMN)
! Search bounds (by rank)
      INTEGER, INTENT(IN) :: lSTART, lEND

      DO IDX = lSTART,lEND
         IF(lPOS >= ENT_POS(IDX-1) .AND. lPOS < ENT_POS(IDX)) EXIT
      ENDDO

      RETURN
      END SUBROUTINE PIC_SEARCH
