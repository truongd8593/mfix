MODULE tmp_array

      Use param
      Use param1

! IMPORTANT:  For using these arrays in a subroutine
! lock the module in the beginning of the subroutine
!   call lock_tmp_array
!
! and unlock the module at the end of the subroutine
!   call unlock_tmp_array
!

!  temporary storage
   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: &
      ARRAY1, ARRAY2, ARRAY3,  ARRAY4

!  temporary storage for 4th order scheme
   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: TMP4

   INTEGER, DIMENSION(:), ALLOCATABLE :: &
      ARRAY1I

   character*3,  dimension(:), pointer :: ARRAY1C     
   
!!!HPF$ align Array1(:) with TT(:)
!!!HPF$ align Array2(:) with TT(:)
!!!HPF$ align Array3(:) with TT(:)
!!!HPF$ align Array4(:) with TT(:)
!!!HPF$ align Array1i(:) with TT(:)
!!!HPF$ align Array1c(:) with TT(:)

      
   LOGICAL :: tmp_array_locked = .false.
   LOGICAL :: tmp4_array_locked = .false.
      
   CONTAINS
      SUBROUTINE lock_tmp_array
        IF(tmp_array_locked)Then
	  Write(*,*)'Error:  Multiple use of tmp_array (tmp_array_mod.f)'
	  CALL MFIX_EXIT(myPE)
	Else
	  tmp_array_locked = .true.
	Endif
      END SUBROUTINE lock_tmp_array
      
      SUBROUTINE unlock_tmp_array
        tmp_array_locked = .false.
      END SUBROUTINE unlock_tmp_array
      
      SUBROUTINE lock_tmp4_array
        IF(tmp4_array_locked)Then
	  Write(*,*)'Error:  Multiple use of tmp_array (tmp4_array_mod.f)'
	  CALL MFIX_EXIT(myPE)
	Else
	  tmp4_array_locked = .true.
	Endif
      END SUBROUTINE lock_tmp4_array
      
      SUBROUTINE unlock_tmp4_array
        tmp4_array_locked = .false.
      END SUBROUTINE unlock_tmp4_array
      

END MODULE tmp_array

!// Comments on the modifications for DMP version implementation      
!//PG Allocatable causes PG internal error, Ed's soln: pointers for ARRAY1c


