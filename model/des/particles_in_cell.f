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

      USE param
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

!-----------------------------------------------

! following quantities are reset every call to particles_in_cell
      PINC(:) = 0

! Call exchange particles - this will exchange particle crossing
! boundaries as well as updates ghost particles information
      CALL DES_PAR_EXCHANGE

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
         IF(DES_POS_NEW(1,L) < XE(I-1)) THEN
            I = I-1
         ELSEIF(DES_POS_NEW(1,L) >= XE(I)) THEN
            I = I+1
         ENDIF

         J = PIJK(L,2)
         IF(DES_POS_NEW(2,L)< YN(J-1)) THEN
            J = J-1
         ELSEIF(DES_POS_NEW(2,L) >= YN(J))THEN
            J = J+1
         ENDIF

         IF(NO_K) THEN
            K = 1
         ELSE
            K = PIJK(L,3)
            IF(DES_POS_NEW(3,L) < ZT(K-1)) THEN
               K = K-1
            ELSEIF(DES_POS_NEW(3,L) >= ZT(K)) THEN
               K = K+1
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
      CALL CALC_INTERP_WEIGHTS

! Calculate mean fields using either interpolation or cell averaging.
      CALL COMP_MEAN_FIELDS

      IF(MPPIC) CALL REPORT_PIC_STATS(RECOVERED, DELETED)

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
         lPOS = DES_POS_NEW(1,L)
         DO I = ISTART2,IEND2
            IF( lPOS >= XE(I-1) .and. lPOS < XE(I)) THEN
               PIJK(L,1) = I
               EXIT
            ENDIF
         ENDDO

         lPOS = DES_POS_NEW(2,L)
         DO J = jstart2,jend2
            IF(lPOS >= YN(J-1) .and. lPOS < YN(J)) THEN
               PIJK(L,2) = J
               EXIT
            ENDIF
         ENDDO

         IF(NO_K) THEN
            K=1
            PIJK(L,3) = 1
         ELSE
            lPOS = DES_POS_NEW(3,L)
            DO K = KSTART2, KEND2
               IF(lPOS >= ZT(K-1) .and. lPOS < ZT(K)) THEN
                  PIJK(L,3) = K
                  EXIT
               ENDIF
            ENDDO
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

