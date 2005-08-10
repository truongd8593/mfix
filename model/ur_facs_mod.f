      MODULE ur_facs


      Use param
      Use param1


!     ur_facs.inc
!     1   p
!     2   rho, ep
!     3   u
!     4   v
!     5   w
!     6   T
!     7   X
!     8   Th
!     9   S	
!                      Under relaxation factors
      DOUBLE PRECISION UR_FAC(9)
      
!
!
      INTEGER           STEPS_LAST
!

!      Under relaxation factors for coefficient update: 0 - every time step (explicit),
!      1 - every iteration (implicit); value between 0-1 for underrelaxation.
!      Note that these values need to be temporarily set to 1 before the calc_coeff
!      call in time_march.  And after the call reset to their original value.
!                      Under relaxation factor for gas-solids drag coefficient
      DOUBLE PRECISION UR_F_gs

      END MODULE ur_facs                                                                         
