!include user defined modules in this file.  Remember to allocate user-defined
!adjustable arrays in usr0.f

MODULE usr


      Use param
      Use param1
      
!
!                      Sherwood number
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: N_sh
!
      
      
END MODULE usr
