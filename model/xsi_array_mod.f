MODULE xsi_array

      Use param
      Use param1

! IMPORTANT:  For using these arrays in a subroutine
! lock the module in the beginning of the subroutine
!   call lock_xsi_array
!
! and unlock the module at the end of the subroutine
!   call unlock_xsi_array
!

!  discretization factors xsi
   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: &
      XSI_e, XSI_n, XSI_t

!HPF$ align Xsi_e(:) with TT(:)
!HPF$ align Xsi_n(:) with TT(:)
!HPF$ align Xsi_t(:) with TT(:)


   LOGICAL :: xsi_array_locked = .false.
      
   CONTAINS
      SUBROUTINE lock_xsi_array
        IF(xsi_array_locked)Then
	  Write(*,*)'Error:  Multiple use of xsi_array (xsi_array_mod.f)'
	  Stop
	Else
	  xsi_array_locked = .true.
	Endif
      END SUBROUTINE lock_xsi_array
      
      SUBROUTINE unlock_xsi_array
        xsi_array_locked = .false.
      END SUBROUTINE unlock_xsi_array
      
END MODULE xsi_array
