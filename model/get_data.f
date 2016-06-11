! -*- f90 -*-
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  MODULE: GET_DATA                                                    C
!  Purpose: read and verify input data, open files                     C
!                                                                      C
!  Author: P. Nicoletti                               Date: 04-DEC-91  C
!  Reviewer: M.SYAMLAL, W.ROGERS, P.NICOLETTI         Date: 24-JAN-92  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
MODULE read_input

  CONTAINS

      SUBROUTINE GET_DATA

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE compar
      USE cutcell
      USE dashboard
      USE des_rxns
      USE des_thermo
      USE discretelement
      USE funits
      USE leqsol
      USE main, only: mfix_dat
      USE mfix_pic
      use mpi_utility, only: BCAST
      USE parallel
      USE param
      USE param1
      USE run

      IMPLICIT NONE

      LOGICAL :: PRESENT

! This module call routines to initialize the namelist variables.
      CALL INIT_NAMELIST
! Read in the namelist variables from the ascii input file.
      CALL READ_NAMELIST(0,MFIX_DAT)
! Set RUN_TYPE to RESTART_1 when adjusting partition
! and read partition layout in gridmap.dat if it exists
      IF(ADJUST_PARTITION) THEN
         RUN_TYPE = 'RESTART_1'

         INQUIRE(FILE='gridmap.dat',EXIST=PRESENT)
         IF(PRESENT) THEN
            IF(MyPE == PE_IO) THEN
               WRITE(*,*)'Reading partition layout from grimap.dat...'
               OPEN(UNIT=777, FILE='gridmap.dat', STATUS='OLD')

                READ (777, *) NODESI,NODESJ,NODESK

                CLOSE(777)
            ENDIF

            CALL BCAST(NODESI)
            CALL BCAST(NODESJ)
            CALL BCAST(NODESK)
         ENDIF

      ENDIF

      RETURN

    END SUBROUTINE GET_DATA

    END MODULE read_input
