	module parallel_mpi

!	A module to carry out init, finalize and check for any parallel errors

	use geometry
	use compar
	implicit none

	contains

	subroutine parallel_init()
	
	integer :: ierr


	call MPI_Init(ierr)
	call MPI_Check( 'parallel_init:MPI_Init ', ierr)

	call MPI_COMM_SIZE( MPI_COMM_WORLD, numPEs, ierr )
	call MPI_Check( 'parallel_init:MPI_Comm_size ', ierr )

	call MPI_COMM_RANK( MPI_COMM_WORLD, myPE, ierr )
	call MPI_Check( 'parallel_init:MPI_Comm_size ', ierr )

	return
	end subroutine parallel_init

	subroutine parallel_fin()

	integer :: ierr

	call MPI_Finalize(ierr)
        call MPI_Check( 'parallel_init:MPI_Finalize ', ierr)

	return
	end subroutine parallel_fin

        subroutine MPI_Check( msg, ierr )
        character(len=*),intent(in) :: msg
        integer, intent(in) :: ierr

        character(len=512) :: errmsg
        integer :: resultlen

        if (ierr .ne. MPI_SUCCESS ) then
                call MPI_Error_string( ierr, errmsg, resultlen )
                print*, 'Error: ', msg
                print*, errmsg(1:resultlen)
                stop '** ERROR ** '
        endif

        return
        end subroutine MPI_Check


	end module parallel_mpi
	
