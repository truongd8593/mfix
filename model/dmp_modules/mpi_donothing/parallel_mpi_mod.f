module parallel_mpi

  !       A module to carry out init, finalize and check for any parallel errors

  use geometry
  use compar
  implicit none

contains

  subroutine parallel_init()
    numPEs = 1
    myPE = 0

    return
  end subroutine parallel_init

  subroutine parallel_fin()

    return
  end subroutine parallel_fin

  subroutine MPI_Check( msg, ierr )
    character(len=*),intent(in) :: msg
    integer, intent(in) :: ierr

    return
  end subroutine MPI_Check


end module parallel_mpi
