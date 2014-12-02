!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: DIFFUSE_MEAN_FIELDS                                     !
!  Author: J.Musser                                   Date: 11-NOV-14  !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DIFFUSE_MEAN_FIELDS

      use discretelement, only: DES_MMAX
      use discretelement, only: DES_ROP_S
      use particle_filter, only: DES_DIFFUSE_MEAN_FIELDS

      IMPLICIT NONE

      INTEGER :: M


      IF(.NOT.DES_DIFFUSE_MEAN_FIELDS) THEN
         write(*,*) 'Not diffusing mean fields'
         return
      ENDIF


      write(*,*) '>>> DIFFUSE_MEAN_FIELDS'

      DO M=1, DES_MMAX
         CALL DIFFUSE_MEAN_FIELD(DES_ROP_S(:,M))
      ENDDO

      write(*,*) '<<< DIFFUSE_MEAN_FIELDS'


      RETURN
      END SUBROUTINE DIFFUSE_MEAN_FIELDS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: DIFFUSE_MEAN_FIELDS                                     !
!  Author: J.Musser                                   Date: 11-NOV-14  !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DIFFUSE_MEAN_FIELD(PHI)

      USE param
      USE param1
      USE run
      USE geometry
      USE compar
      USE sendrecv
      Use xsi_array
      USE mpi_utility
      USE indices

      use particle_filter, only: DIF_TSTOP
      use discretelement, only: MAX_RADIUS

      use ambm, only: A_M, B_M
      use leqsol

      use rxns

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(INOUT) :: PHI(DIMENSION_3)

! Integer error flag
      INTEGER :: IER

! Linear equation solver method and iterations
      INTEGER :: LEQM, LEQI, IJK

      DOUBLE PRECISION :: DIF_TIME
      DOUBLE PRECISION :: DIF_DT

      IER = 0

      DIF_TIME = 0.0d0
      DIF_DT = DIF_TSTOP*0.05 ! min(1.0d-3, DIF_TSTOP)

      write(*,*) 'DIF TSTOP:    ', DIF_TSTOP
      write(*,*) 'DIF DT:       ', DIF_DT

      IF(DIF_TSTOP > 0.0) THEN
         write(*,*) 'DIF STEPS:    ', int(DIF_TSTOP/DIF_DT)
      ELSE
         write(*,*) 'DIF STEPS:    ', 0.0
      ENDIF

      DO WHILE(DIF_TIME < DIF_TSTOP)

         CALL INIT_AB_M (A_M, B_M, IJKMAX2, 0, IER)

         CALL DIF_PHI_DES(0, A_M, B_M, IER)

         CALL DIF_PHI_BC_DES(PHI, 0, A_M, B_M, IER)

         CALL DIF_PHI_SOURCE_DES(PHI, 0, A_M, B_M, DIF_DT, IER)

         CALL ADJUST_LEQ(0.0d0, 10, 2, LEQI, LEQM, IER)

         CALL SOLVE_LIN_EQ ('DEM', 9, PHI, A_M, B_M, 0, LEQI, LEQM, &
            LEQ_SWEEP(9), LEQ_TOL(9), LEQ_PC(9), IER)

         DIF_TIME = DIF_TIME + DIF_DT
      ENDDO

      RETURN
      END SUBROUTINE DIFFUSE_MEAN_FIELD
