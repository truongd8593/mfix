MODULE ambm

      Use param
      Use param1


! IMPORTANT:  For using these arrays in a subroutine
! lock the module in the beginning of the subroutine
!   call lock_ambm
!
! and unlock the module at the end of the subroutine
!   call unlock_ambm
!

!
!linear equation matrix and vector
      DOUBLE PRECISION, DIMENSION(:, :, :), ALLOCATABLE :: A_m
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: B_m

!HPF$ align A_m(:, *, *) with TT(:)
!HPF$ align B_m(:, *) with TT(:)


      
   LOGICAL :: ambm_locked = .false.
      
   CONTAINS
      SUBROUTINE lock_ambm
        IF(ambm_locked)Then
	  Write(*,*)'Error:  Multiple use of ambm (ambm_mod.f)'
	  Stop
	Else
	  ambm_locked = .true.
	Endif
      END SUBROUTINE lock_ambm
      
      SUBROUTINE unlock_ambm
        ambm_locked = .false.
      END SUBROUTINE unlock_ambm
      
END MODULE ambm                                                          
