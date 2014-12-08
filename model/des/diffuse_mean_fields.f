!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: DIFFUSE_MEAN_FIELDS                                     !
!  Author: J.Musser                                   Date: 11-NOV-14  !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DIFFUSE_MEAN_FIELD(PHI, VNAME)

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

      use error_manager

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(INOUT) :: PHI(DIMENSION_3)
      CHARACTER(LEN=*), INTENT(IN) :: VNAME

! Integer error flag
      INTEGER :: IER

! Linear equation solver method and iterations
      INTEGER :: LEQM, LEQI, IJK

      DOUBLE PRECISION :: DIF_TIME
      DOUBLE PRECISION :: DIF_DT

      DOUBLE PRECISION :: WALL_START, WALL_TIME

      LOGICAL, PARAMETER :: setDBG = .FALSE.

      IER = 0

      DIF_TIME = 0.0d0
      DIF_DT = DIF_TSTOP*0.05 ! min(1.0d-3, DIF_TSTOP)


      IF(setDBG) THEN
         WRITE(ERR_MSG, 9000) VNAME, DIF_TSTOP, DIF_DT, &
            1+int(DIF_TSTOP/(DIF_DT+DIF_TSTOP*1.0d-3))
         CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
         WALL_START = WALL_TIME()
      ENDIF

 9000 FORMAT(/3x,'Diffusing Variable: ',A,/5X,'Diffusion time scale: ',&
         g11.4,/5x,'DT: ',g11.4,3x,'Steps:',I4)


      DO WHILE(DIF_TIME < DIF_TSTOP)

         CALL INIT_AB_M (A_M, B_M, IJKMAX2, 0, IER)

         CALL DIF_PHI_DES(0, A_M, B_M, IER)

         CALL DIF_PHI_BC_DES(PHI, 0, A_M, B_M, IER)

         CALL DIF_PHI_SOURCE_DES(PHI, 0, A_M, B_M, DIF_DT, IER)

         CALL ADJUST_LEQ(0.0d0, 10, 2, LEQI, LEQM, IER)

         CALL SOLVE_LIN_EQ (VNAME, 9, PHI, A_M, B_M, 0, LEQI, LEQM, &
            LEQ_SWEEP(9), LEQ_TOL(9), LEQ_PC(9), IER)

         DIF_TIME = DIF_TIME + DIF_DT
      ENDDO


      IF(setDBG) THEN
         WRITE(ERR_MSG, 9001) WALL_TIME() - WALL_START
         CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
      ENDIF

 9001 FORMAT(5x,'Wall Time: ',g11.4)

      RETURN
      END SUBROUTINE DIFFUSE_MEAN_FIELD
