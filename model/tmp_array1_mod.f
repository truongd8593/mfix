MODULE tmp_array1

      Use param
      Use param1

! IMPORTANT:  For using these arrays in a subroutine
! lock the module in the beginning of the subroutine
!   call lock_tmp_array1
!
! and unlock the module at the end of the subroutine
!   call unlock_tmp_array1
!

!  temporary storage
   DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: &
      ARRAYm1

!!!HPF$ align Arraym1(:, *) with TT(:)

      
   LOGICAL :: tmp_array1_locked = .false.
      
   CONTAINS
      SUBROUTINE lock_tmp_array1
        IF(tmp_array1_locked)Then
	  Write(*,*)'Error:  Multiple use of tmp_array1 (tmp_array1_mod.f)'
	  Stop
	Else
	  tmp_array1_locked = .true.
	Endif
      END SUBROUTINE lock_tmp_array1
      
      SUBROUTINE unlock_tmp_array1
        tmp_array1_locked = .false.
      END SUBROUTINE unlock_tmp_array1
      
END MODULE tmp_array1
