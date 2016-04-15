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
      USE iterate, only: max_nit
      USE leqsol
      USE main, only: mfix_dat
      USE mfix_pic
      USE parallel
      USE param
      USE param1
      USE run
      USE turb, only: L_SCALE0, K_EPSILON, l_scale
      USE visc_g, only: MU_GMAX

      IMPLICIT NONE

! This module call routines to initialize the namelist variables.
      CALL INIT_NAMELIST
! Read in the namelist variables from the ascii input file.
<<<<<<< 1782d64de81a0bb387e58b08b743190186856b10
      print *,"FOUL WASTELAND ", mfix_dat
=======
>>>>>>> specify alternate mfix.dat filename on the command line with -f
      CALL READ_NAMELIST(0,MFIX_DAT)

      RETURN

    END SUBROUTINE GET_DATA

    END MODULE read_input
