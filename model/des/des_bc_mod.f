!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: DES_INLET                                              !
!                                                                      !
!  Purpose: Common elements needed for the des mass inflow boundary    !
!  condition.                                                          !
!                                                                      !
!  Author: J.Musser                                   Date: 13-Jul-09  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      MODULE DES_BC

!     Logical Flag from mfix.dat file
      LOGICAL DES_MI

!     Physical injection location
      DOUBLE PRECISION DES_BC_X_w, DES_BC_X_e
      DOUBLE PRECISION DES_BC_Y_s, DES_BC_Y_n
      DOUBLE PRECISION DES_BC_Z_b, DES_BC_Z_t

!     Specification for inflow into the system 
      DOUBLE PRECISION DES_BC_VOLFLOW_s
      DOUBLE PRECISION DES_BC_MASSFLOW_s

!     Particles Injection Factor
      INTEGER PI_FACTOR

!     Particle Injection Count (injection number)
      INTEGER PI_COUNT

      END MODULE DES_BC

