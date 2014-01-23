!	MPI Modules written at ORNL by Ed and Sreekanth for MFIX
!	under joint effort with FETC - 06/08/99.

	module mpi_utility

!	module to perform most of the mpi functionalities like scatter,
!	gather, bcast, globalsum and so on.

	use geometry
	use compar
	use parallel_mpi
	use debug
	use indices
	implicit none

!	Object-oriented approach to direct to the correct procedure
!	depending on the argument type. i stands for integer, r for real
!	and d for double precision. 0 for scalar, 1 for vector, 2 for
!	2-D array and similarly 3.

!==============================================================================
!  JFD: Interfaces used for vtk file writting (Cartesian grid):
!==============================================================================

	interface allgather
	module  procedure allgather_1i
	end interface 

	interface gatherv
	module  procedure gatherv_1d
	end interface 

!==============================================================================
!  JFD: End of Interfaces used for vtk file writting (Cartesian grid):
!==============================================================================


	interface scatter
	module  procedure scatter_1i, scatter_2i, scatter_3i, &
                          scatter_1r, scatter_2r, scatter_3r, &
                          scatter_1d, scatter_2d, scatter_3d, &
			  scatter_1c,scatter_1l
	end interface 

	interface gather
	module  procedure gather_1i, gather_2i, gather_3i, &
                          gather_1r, gather_2r, gather_3r, &
                          gather_1d, gather_2d, gather_3d, &
			  gather_1c, gather_1l
	end interface 

	interface bcast
	module  procedure bcast_0i, bcast_1i, bcast_2i, bcast_3i, &
                          bcast_0r, bcast_1r, bcast_2r, bcast_3r, &
                          bcast_0d, bcast_1d, bcast_2d, bcast_3d, &
			  bcast_0l, bcast_1l, bcast_0c, bcast_1c
	end interface 

        interface global_sum
        module  procedure global_sum_0i, global_sum_1i, global_sum_2i, global_sum_3i, &
                          global_sum_0r, global_sum_1r, global_sum_2r, global_sum_3r, &
                          global_sum_0d, global_sum_1d, global_sum_2d, global_sum_3d
        end interface

        interface global_all_sum
        module  procedure &
                          global_all_sum_0i, global_all_sum_1i, &
                          global_all_sum_2i, global_all_sum_3i, &
                          global_all_sum_0r, global_all_sum_1r, &
                          global_all_sum_2r, global_all_sum_3r, &
                          global_all_sum_0d, global_all_sum_1d, &
                          global_all_sum_2d, global_all_sum_3d, &
                          global_all_sum_onevar_0i, global_all_sum_onevar_1i, &
                          global_all_sum_onevar_2i, global_all_sum_onevar_3i, &
                          global_all_sum_onevar_0r, global_all_sum_onevar_1r, &
                          global_all_sum_onevar_2r, global_all_sum_onevar_3r, &
                          global_all_sum_onevar_0d, global_all_sum_onevar_1d, &
                          global_all_sum_onevar_2d, global_all_sum_onevar_3d
        end interface

        interface global_min
        module  procedure global_min_0i, global_min_1i, global_min_2i, global_min_3i, &
                          global_min_0r, global_min_1r, global_min_2r, global_min_3r, &
                          global_min_0d, global_min_1d, global_min_2d, global_min_3d
        end interface

        interface global_all_min
        module  procedure &
                          global_all_min_0i, global_all_min_1i, &
                          global_all_min_2i, global_all_min_3i, &
                          global_all_min_0r, global_all_min_1r, &
                          global_all_min_2r, global_all_min_3r, &
                          global_all_min_0d, global_all_min_1d, &
                          global_all_min_2d, global_all_min_3d, &
                          global_all_min_onevar_0i, global_all_min_onevar_1i, &
                          global_all_min_onevar_2i, global_all_min_onevar_3i, &
                          global_all_min_onevar_0r, global_all_min_onevar_1r, &
                          global_all_min_onevar_2r, global_all_min_onevar_3r, &
                          global_all_min_onevar_0d, global_all_min_onevar_1d, &
                          global_all_min_onevar_2d, global_all_min_onevar_3d
        end interface

        interface global_max
        module  procedure global_max_0i, global_max_1i, global_max_2i, global_max_3i, &
                          global_max_0r, global_max_1r, global_max_2r, global_max_3r, &
                          global_max_0d, global_max_1d, global_max_2d, global_max_3d
        end interface

        interface global_all_max
        module  procedure &
                          global_all_max_0i, global_all_max_1i, &
                          global_all_max_2i, global_all_max_3i, &
                          global_all_max_0r, global_all_max_1r, &
                          global_all_max_2r, global_all_max_3r, &
                          global_all_max_0d, global_all_max_1d, &
                          global_all_max_2d, global_all_max_3d, &
                          global_all_max_onevar_0i, global_all_max_onevar_1i, &
                          global_all_max_onevar_2i, global_all_max_onevar_3i, &
                          global_all_max_onevar_0r, global_all_max_onevar_1r, &
                          global_all_max_onevar_2r, global_all_max_onevar_3r, &
                          global_all_max_onevar_0d, global_all_max_onevar_1d, &
                          global_all_max_onevar_2d, global_all_max_onevar_3d
        end interface

        interface global_all_and
        module procedure &
                global_all_and_0d, global_all_and_1d, &
                global_all_and_onevar_0d, global_all_and_onevar_1d
        end interface

        interface global_all_or
        module procedure &
                global_all_or_0d, global_all_or_1d, &
                global_all_or_onevar_0d, global_all_or_onevar_1d
        end interface

        contains

!==============================================================================
!  JFD: Subroutines used for vtk file writting (Cartesian grid):
!==============================================================================

	subroutine allgather_1i( lbuf, gbuf, idebug )
        integer, intent(in) :: lbuf
        integer, intent(out), dimension(:) :: gbuf
	integer, optional, intent(in) :: idebug
!       do nothing
        gbuf = 0
	return
	end subroutine allgather_1i


	subroutine gatherv_1d( lbuf, sendcnt, gbuf, rcount, disp, mroot, idebug )
        double precision, intent(in), dimension(:) :: lbuf
        double precision, intent(out), dimension(:) :: gbuf
        integer, intent(in), dimension(:) :: rcount
        integer, intent(in), dimension(:) :: disp
	integer, optional, intent(in) :: mroot, idebug
	integer :: sendtype,recvtype,sendcnt,recvcnt,lroot,ierr,lidebug
	gbuf = lbuf
	return
	end subroutine gatherv_1d

!==============================================================================
!  JFD: End of Subroutines used for vtk file writting (Cartesian grid):
!==============================================================================


!	Routine to scatter gbuf available on root to all the processors

	subroutine scatter_1i( lbuf, gbuf, mroot, idebug )
        integer, intent(in), dimension(:) :: gbuf
        integer, intent(out), dimension(:) :: lbuf
	integer, optional, intent(in) :: mroot, idebug
	lbuf = gbuf
	return
	end subroutine scatter_1i

        subroutine scatter_2i( lbuf, gbuf, mroot, idebug )
        integer, intent(in), dimension(:,:) :: gbuf       
        integer, intent(out), dimension(:,:) :: lbuf
        integer, optional, intent(in) :: mroot, idebug
	lbuf = gbuf
	return
	end subroutine scatter_2i

        subroutine scatter_3i( lbuf, gbuf, mroot, idebug )
        integer, intent(in), dimension(:,:,:) :: gbuf       
        integer, intent(out), dimension(:,:,:) :: lbuf
        integer, optional, intent(in) :: mroot, idebug
	lbuf = gbuf
	return
	end subroutine scatter_3i

	subroutine scatter_1r( lbuf, gbuf, mroot, idebug )
        real, intent(in), dimension(:) :: gbuf       
        real, intent(out), dimension(:) :: lbuf
	integer, optional, intent(in) :: mroot, idebug
	lbuf = gbuf
	return
	end subroutine scatter_1r

	
        subroutine scatter_2r( lbuf, gbuf, mroot, idebug )
        real, intent(in), dimension(:,:) :: gbuf       
        real, intent(out), dimension(:,:) :: lbuf
        integer, optional, intent(in) :: mroot, idebug
	lbuf = gbuf
	return
	end subroutine scatter_2r

        subroutine scatter_3r( lbuf, gbuf, mroot, idebug )
        real, intent(in), dimension(:,:,:) :: gbuf       
        real, intent(out), dimension(:,:,:) :: lbuf
        integer, optional, intent(in) :: mroot, idebug
	lbuf = gbuf
	return
	end subroutine scatter_3r


	subroutine scatter_1d( lbuf, gbuf, mroot, idebug )
        double precision, intent(in), dimension(:) :: gbuf       
        double precision, intent(out), dimension(:) :: lbuf
	integer, optional, intent(in) :: mroot, idebug
	lbuf = gbuf
	return
	end subroutine scatter_1d

	
        subroutine scatter_2d( lbuf, gbuf, mroot, idebug )
        double precision, intent(in), dimension(:,:) :: gbuf       
        double precision, intent(out), dimension(:,:) :: lbuf
        integer, optional, intent(in) :: mroot, idebug
	lbuf = gbuf
	return
	end subroutine scatter_2d

        subroutine scatter_3d( lbuf, gbuf, mroot, idebug )
        double precision, intent(in), dimension(:,:,:) :: gbuf       
        double precision, intent(out), dimension(:,:,:) :: lbuf
        integer, optional, intent(in) :: mroot, idebug
	lbuf = gbuf
	return
	end subroutine scatter_3d


        subroutine scatter_1c( lbuf, gbuf, mroot, idebug )
        character(len=*), intent(in), dimension(:) :: gbuf
        character(len=*), intent(out), dimension(:) :: lbuf
        integer, optional, intent(in) :: mroot, idebug
	lbuf = gbuf
	return
	end subroutine scatter_1c


        subroutine scatter_1l( lbuf, gbuf, mroot, idebug )
        logical, intent(in), dimension(:) :: gbuf
        logical, intent(out), dimension(:) :: lbuf
        integer, optional, intent(in) :: mroot, idebug
	lbuf = gbuf
	return
	end subroutine scatter_1l


!	Routines to gather lbuf from individual processors and put it on
!	processor root in gbuf
!	Logic is similar to the scatter routines above.

	subroutine gather_1i( lbuf, gbuf, mroot, idebug )
        integer, intent(in), dimension(:) :: lbuf       
        integer, intent(out), dimension(:) :: gbuf
	integer, optional, intent(in) :: mroot, idebug
	gbuf = lbuf
	return
	end subroutine gather_1i

	
        subroutine gather_2i( lbuf, gbuf, mroot, idebug )
        integer, intent(in), dimension(:,:) :: lbuf
        integer, intent(out), dimension(:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug
	gbuf = lbuf
	return
	end subroutine gather_2i

        subroutine gather_3i( lbuf, gbuf, mroot, idebug )
        integer, intent(in), dimension(:,:,:) :: lbuf
        integer, intent(out), dimension(:,:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug
	gbuf = lbuf
	return
	end subroutine gather_3i

	subroutine gather_1r( lbuf, gbuf, mroot, idebug )
        real, intent(in), dimension(:) :: lbuf
        real, intent(out), dimension(:) :: gbuf
	integer, optional, intent(in) :: mroot, idebug
	gbuf = lbuf
	return
	end subroutine gather_1r

	
        subroutine gather_2r( lbuf, gbuf, mroot, idebug )
        real, intent(in), dimension(:,:) :: lbuf
        real, intent(out), dimension(:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug
	gbuf = lbuf
	return
	end subroutine gather_2r

        subroutine gather_3r( lbuf, gbuf, mroot, idebug )
        real, intent(in), dimension(:,:,:) :: lbuf
        real, intent(out), dimension(:,:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug
	gbuf = lbuf
	return
	end subroutine gather_3r


	subroutine gather_1d( lbuf, gbuf, mroot, idebug )
        double precision, intent(in), dimension(:) :: lbuf
        double precision, intent(out), dimension(:) :: gbuf
	integer, optional, intent(in) :: mroot, idebug
	gbuf = lbuf
	return
	end subroutine gather_1d

	
        subroutine gather_2d( lbuf, gbuf, mroot, idebug )
        double precision, intent(in), dimension(:,:) :: lbuf
        double precision, intent(out), dimension(:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug
	gbuf = lbuf
	return
	end subroutine gather_2d

        subroutine gather_3d( lbuf, gbuf, mroot, idebug )
        double precision, intent(in), dimension(:,:,:) :: lbuf
        double precision, intent(out), dimension(:,:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug
	gbuf = lbuf
	return
	end subroutine gather_3d


        subroutine gather_1c( lbuf, gbuf, mroot, idebug )
        character(len=*), intent(in), dimension(:) :: lbuf
        character(len=*), intent(out), dimension(:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug
	gbuf = lbuf
	return
	end subroutine gather_1c


        subroutine gather_1l( lbuf, gbuf, mroot, idebug )
        logical, intent(in), dimension(:) :: lbuf
        logical, intent(out), dimension(:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug
	gbuf = lbuf
	return
	end subroutine gather_1l


!	Routines to broadcast information from processor 0 in buffer to all
!	the processors	

        subroutine bcast_0i( buffer, mroot, idebug )
        integer, intent(inout) :: buffer
        integer, optional, intent(in) :: mroot, idebug
	return
	end subroutine bcast_0i


	subroutine bcast_1i( buffer, mroot, idebug )
        integer, intent(inout), dimension(:) :: buffer       
	integer, optional, intent(in) :: mroot, idebug
	return
	end subroutine bcast_1i

	
        subroutine bcast_2i( buffer, mroot, idebug )
        integer, intent(inout), dimension(:,:) :: buffer
        integer, optional, intent(in) :: mroot, idebug
	return
	end subroutine bcast_2i

        subroutine bcast_3i( buffer, mroot, idebug )
        integer, intent(inout), dimension(:,:,:) :: buffer
        integer, optional, intent(in) :: mroot, idebug
	return
	end subroutine bcast_3i

        subroutine bcast_0r( buffer, mroot, idebug )
        real, intent(inout) :: buffer
        integer, optional, intent(in) :: mroot, idebug

	return
        end subroutine bcast_0r


	subroutine bcast_1r( buffer, mroot, idebug )
        real, intent(inout), dimension(:) :: buffer
	integer, optional, intent(in) :: mroot, idebug

	return
	end subroutine bcast_1r

	
        subroutine bcast_2r( buffer, mroot, idebug )
        real, intent(inout), dimension(:,:) :: buffer
        integer, optional, intent(in) :: mroot, idebug

	return
	end subroutine bcast_2r

        subroutine bcast_3r( buffer, mroot, idebug )
        real, intent(inout), dimension(:,:,:) :: buffer
        integer, optional, intent(in) :: mroot, idebug

	return
        end subroutine bcast_3r

        subroutine bcast_0d( buffer, mroot, idebug )
        double precision, intent(inout) :: buffer
        integer, optional, intent(in) :: mroot, idebug

	return
        end subroutine bcast_0d


	subroutine bcast_1d( buffer, mroot, idebug )
        double precision, intent(inout), dimension(:) :: buffer
	integer, optional, intent(in) :: mroot, idebug

	return
	end subroutine bcast_1d

	
        subroutine bcast_2d( buffer, mroot, idebug )
        double precision, intent(inout), dimension(:,:) :: buffer
        integer, optional, intent(in) :: mroot, idebug

	return
	end subroutine bcast_2d

        subroutine bcast_3d( buffer, mroot, idebug )
        double precision, intent(inout), dimension(:,:,:) :: buffer
        integer, optional, intent(in) :: mroot, idebug

	return
        end subroutine bcast_3d

        subroutine bcast_0c( buffer, mroot, idebug )
        character(len=*), intent(inout) :: buffer
        integer, optional, intent(in) :: mroot, idebug
	character, allocatable, dimension(:) :: buffer1

	return
        end subroutine bcast_0c


        subroutine bcast_1c( buffer, mroot, idebug )
        character(len=*), intent(inout), dimension(:) :: buffer
        integer, optional, intent(in) :: mroot, idebug
        character, allocatable, dimension(:) :: buffer1

	return
        end subroutine bcast_1c

        subroutine bcast_0l( buffer, mroot, idebug )
        logical, intent(inout) :: buffer
        integer, optional, intent(in) :: mroot, idebug

	return
        end subroutine bcast_0l


        subroutine bcast_1l( buffer, mroot, idebug )
        logical, intent(inout), dimension(:) :: buffer
        integer, optional, intent(in) :: mroot, idebug

	return
        end subroutine bcast_1l


!	Procedures to do global operations (Sum, Min, Max). _all_ routines
!	send the information to all the processors otherwise they are
!	kept on processor 0.

        subroutine global_sum_0i( lbuf, gbuf, mroot, idebug )
        integer, intent(in) :: lbuf
        integer, intent(out) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

	gbuf = lbuf

	return
        end subroutine global_sum_0i


	subroutine global_sum_1i( lbuf, gbuf, mroot, idebug )
        integer, intent(in), dimension(:) :: lbuf       
        integer, intent(out), dimension(:) :: gbuf
	integer, optional, intent(in) :: mroot, idebug

        gbuf = lbuf

	return
	end subroutine global_sum_1i

        subroutine global_sum_2i( lbuf, gbuf, mroot, idebug )
        integer, intent(in), dimension(:,:) :: lbuf
        integer, intent(out), dimension(:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        gbuf = lbuf

	return
	end subroutine global_sum_2i

        subroutine global_sum_3i( lbuf, gbuf, mroot, idebug )
        integer, intent(in), dimension(:,:,:) :: lbuf
        integer, intent(out), dimension(:,:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        gbuf = lbuf

	return
        end subroutine global_sum_3i

        subroutine global_sum_0r( lbuf, gbuf, mroot, idebug )
        real, intent(in) :: lbuf
        real, intent(out) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        gbuf = lbuf

	return
        end subroutine global_sum_0r


	subroutine global_sum_1r( lbuf, gbuf, mroot, idebug )
        real, intent(in), dimension(:) :: lbuf       
        real, intent(out), dimension(:) :: gbuf
	integer, optional, intent(in) :: mroot, idebug

        gbuf = lbuf

	return
	end subroutine global_sum_1r

        subroutine global_sum_2r( lbuf, gbuf, mroot, idebug )
        real, intent(in), dimension(:,:) :: lbuf
        real, intent(out), dimension(:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        gbuf = lbuf

	return
	end subroutine global_sum_2r

        subroutine global_sum_3r( lbuf, gbuf, mroot, idebug )
        real, intent(in), dimension(:,:,:) :: lbuf
        real, intent(out), dimension(:,:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        gbuf = lbuf

	return
        end subroutine global_sum_3r

        subroutine global_sum_0d( lbuf, gbuf, mroot, idebug )
        double precision, intent(in) :: lbuf
        double precision, intent(out) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        gbuf = lbuf

	return
        end subroutine global_sum_0d


	subroutine global_sum_1d( lbuf, gbuf, mroot, idebug )
        double precision, intent(in), dimension(:) :: lbuf       
        double precision, intent(out), dimension(:) :: gbuf
	integer, optional, intent(in) :: mroot, idebug

        gbuf = lbuf

	return
	end subroutine global_sum_1d

        subroutine global_sum_2d( lbuf, gbuf, mroot, idebug )
        double precision, intent(in), dimension(:,:) :: lbuf
        double precision, intent(out), dimension(:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        gbuf = lbuf

	return
	end subroutine global_sum_2d

        subroutine global_sum_3d( lbuf, gbuf, mroot, idebug )
        double precision, intent(in), dimension(:,:,:) :: lbuf
        double precision, intent(out), dimension(:,:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        gbuf = lbuf

	return
        end subroutine global_sum_3d

        subroutine global_all_sum_onevar_0d( gbuf )
        doubleprecision, intent(inout) :: gbuf
        doubleprecision :: lbuf

        return
        end subroutine global_all_sum_onevar_0d


        subroutine global_all_sum_onevar_1d( gbuf )
        doubleprecision, dimension(:), intent(inout) :: gbuf
        doubleprecision, dimension(size(gbuf)) :: lbuf

        return
        end subroutine global_all_sum_onevar_1d

        subroutine global_all_sum_onevar_2d( gbuf )
        doubleprecision, dimension(:,:), intent(inout) :: gbuf
        doubleprecision, dimension(size(gbuf,1),size(gbuf,2)) :: lbuf

        return
        end subroutine global_all_sum_onevar_2d


        subroutine global_all_sum_onevar_3d( gbuf )
        doubleprecision, dimension(:,:,:), intent(inout) :: gbuf
        doubleprecision, dimension(size(gbuf,1),size(gbuf,2),size(gbuf,3)) :: lbuf

        return
        end subroutine global_all_sum_onevar_3d



        subroutine global_all_sum_onevar_0i( gbuf )
        integer, intent(inout) :: gbuf
        integer :: lbuf

        return
        end subroutine global_all_sum_onevar_0i

        subroutine global_all_sum_onevar_1i( gbuf )
        integer, dimension(:), intent(inout) :: gbuf
        integer, dimension(size(gbuf)) :: lbuf

        return
        end subroutine global_all_sum_onevar_1i

        subroutine global_all_sum_onevar_2i( gbuf )
        integer, dimension(:,:), intent(inout) :: gbuf
        integer, dimension(size(gbuf,1),size(gbuf,2)) :: lbuf

        return
        end subroutine global_all_sum_onevar_2i


        subroutine global_all_sum_onevar_3i( gbuf )
        integer, dimension(:,:,:), intent(inout) :: gbuf
        integer, dimension(size(gbuf,1),size(gbuf,2),size(gbuf,3)) :: lbuf

        return
        end subroutine global_all_sum_onevar_3i

        subroutine global_all_sum_onevar_0r( gbuf )
        real, intent(inout) :: gbuf
        real :: lbuf

        return
        end subroutine global_all_sum_onevar_0r


        subroutine global_all_sum_onevar_1r( gbuf )
        real, dimension(:), intent(inout) :: gbuf
        real, dimension(size(gbuf)) :: lbuf

        return
        end subroutine global_all_sum_onevar_1r

        subroutine global_all_sum_onevar_2r( gbuf )
        real, dimension(:,:), intent(inout) :: gbuf
        real, dimension(size(gbuf,1),size(gbuf,2)) :: lbuf

        return
        end subroutine global_all_sum_onevar_2r


        subroutine global_all_sum_onevar_3r( gbuf )
        real, dimension(:,:,:), intent(inout) :: gbuf
        real, dimension(size(gbuf,1),size(gbuf,2),size(gbuf,3)) :: lbuf

        return
        end subroutine global_all_sum_onevar_3r


        subroutine global_all_sum_0i( lbuf, gbuf, mroot, idebug )
        integer, intent(in) :: lbuf
        integer, intent(out) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        gbuf = lbuf

	return
        end subroutine global_all_sum_0i


	subroutine global_all_sum_1i( lbuf, gbuf, mroot, idebug )
        integer, intent(in), dimension(:) :: lbuf       
        integer, intent(out), dimension(:) :: gbuf
	integer, optional, intent(in) :: mroot, idebug

        gbuf = lbuf

	return
	end subroutine global_all_sum_1i

        subroutine global_all_sum_2i( lbuf, gbuf, mroot, idebug )
        integer, intent(in), dimension(:,:) :: lbuf
        integer, intent(out), dimension(:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        gbuf = lbuf

	return
	end subroutine global_all_sum_2i

        subroutine global_all_sum_3i( lbuf, gbuf, mroot, idebug )
        integer, intent(in), dimension(:,:,:) :: lbuf
        integer, intent(out), dimension(:,:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        gbuf = lbuf

	return
        end subroutine global_all_sum_3i

        subroutine global_all_sum_0r( lbuf, gbuf, mroot, idebug )
        real, intent(in) :: lbuf
        real, intent(out) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        gbuf = lbuf

	return
        end subroutine global_all_sum_0r


	subroutine global_all_sum_1r( lbuf, gbuf, mroot, idebug )
        real, intent(in), dimension(:) :: lbuf       
        real, intent(out), dimension(:) :: gbuf
	integer, optional, intent(in) :: mroot, idebug

        gbuf = lbuf

	return
	end subroutine global_all_sum_1r

        subroutine global_all_sum_2r( lbuf, gbuf, mroot, idebug )
        real, intent(in), dimension(:,:) :: lbuf
        real, intent(out), dimension(:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        gbuf = lbuf

	return
	end subroutine global_all_sum_2r

        subroutine global_all_sum_3r( lbuf, gbuf, mroot, idebug )
        real, intent(in), dimension(:,:,:) :: lbuf
        real, intent(out), dimension(:,:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        gbuf = lbuf

	return
        end subroutine global_all_sum_3r

        subroutine global_all_sum_0d( lbuf, gbuf, mroot, idebug )
        double precision, intent(in) :: lbuf
        double precision, intent(out) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        gbuf = lbuf

	return
        end subroutine global_all_sum_0d


	subroutine global_all_sum_1d( lbuf, gbuf, mroot, idebug )
        double precision, intent(in), dimension(:) :: lbuf       
        double precision, intent(out), dimension(:) :: gbuf
	integer, optional, intent(in) :: mroot, idebug

        gbuf = lbuf

	return
	end subroutine global_all_sum_1d

        subroutine global_all_sum_2d( lbuf, gbuf, mroot, idebug )
        double precision, intent(in), dimension(:,:) :: lbuf
        double precision, intent(out), dimension(:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        gbuf = lbuf

	return
	end subroutine global_all_sum_2d

        subroutine global_all_sum_3d( lbuf, gbuf, mroot, idebug )
        double precision, intent(in), dimension(:,:,:) :: lbuf
        double precision, intent(out), dimension(:,:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        gbuf = lbuf

	return
        end subroutine global_all_sum_3d

        subroutine global_min_0i( lbuf, gbuf, mroot, idebug )
        integer, intent(in) :: lbuf
        integer, intent(out) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        gbuf = lbuf

	return
        end subroutine global_min_0i


	subroutine global_min_1i( lbuf, gbuf, mroot, idebug )
        integer, intent(in), dimension(:) :: lbuf       
        integer, intent(out), dimension(:) :: gbuf
	integer, optional, intent(in) :: mroot, idebug

        gbuf = lbuf

	return
	end subroutine global_min_1i

        subroutine global_min_2i( lbuf, gbuf, mroot, idebug )
        integer, intent(in), dimension(:,:) :: lbuf
        integer, intent(out), dimension(:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        gbuf = lbuf

	return
	end subroutine global_min_2i

        subroutine global_min_3i( lbuf, gbuf, mroot, idebug )
        integer, intent(in), dimension(:,:,:) :: lbuf
        integer, intent(out), dimension(:,:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        gbuf = lbuf

	return
        end subroutine global_min_3i

        subroutine global_min_0r( lbuf, gbuf, mroot, idebug )
        real, intent(in) :: lbuf
        real, intent(out) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        gbuf = lbuf

	return
        end subroutine global_min_0r


	subroutine global_min_1r( lbuf, gbuf, mroot, idebug )
        real, intent(in), dimension(:) :: lbuf       
        real, intent(out), dimension(:) :: gbuf
	integer, optional, intent(in) :: mroot, idebug

        gbuf = lbuf

	return
	end subroutine global_min_1r

        subroutine global_min_2r( lbuf, gbuf, mroot, idebug )
        real, intent(in), dimension(:,:) :: lbuf
        real, intent(out), dimension(:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        gbuf = lbuf

	return
	end subroutine global_min_2r

        subroutine global_min_3r( lbuf, gbuf, mroot, idebug )
        real, intent(in), dimension(:,:,:) :: lbuf
        real, intent(out), dimension(:,:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        gbuf = lbuf

	return
        end subroutine global_min_3r

        subroutine global_min_0d( lbuf, gbuf, mroot, idebug )
        double precision, intent(in) :: lbuf
        double precision, intent(out) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        gbuf = lbuf

	return
        end subroutine global_min_0d


	subroutine global_min_1d( lbuf, gbuf, mroot, idebug )
        double precision, intent(in), dimension(:) :: lbuf       
        double precision, intent(out), dimension(:) :: gbuf
	integer, optional, intent(in) :: mroot, idebug

        gbuf = lbuf

	return
	end subroutine global_min_1d

        subroutine global_min_2d( lbuf, gbuf, mroot, idebug )
        double precision, intent(in), dimension(:,:) :: lbuf
        double precision, intent(out), dimension(:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        gbuf = lbuf

	return
	end subroutine global_min_2d

        subroutine global_min_3d( lbuf, gbuf, mroot, idebug )
        double precision, intent(in), dimension(:,:,:) :: lbuf
        double precision, intent(out), dimension(:,:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        gbuf = lbuf

	return
        end subroutine global_min_3d

        subroutine global_all_min_onevar_0d( gbuf )
        doubleprecision, intent(inout) :: gbuf
        doubleprecision :: lbuf

        return
        end subroutine global_all_min_onevar_0d


        subroutine global_all_min_onevar_1d( gbuf )
        doubleprecision, dimension(:), intent(inout) :: gbuf
        doubleprecision, dimension(size(gbuf)) :: lbuf

        return
        end subroutine global_all_min_onevar_1d

        subroutine global_all_min_onevar_2d( gbuf )
        doubleprecision, dimension(:,:), intent(inout) :: gbuf
        doubleprecision, dimension(size(gbuf,1),size(gbuf,2)) :: lbuf

        return
        end subroutine global_all_min_onevar_2d


        subroutine global_all_min_onevar_3d( gbuf )
        doubleprecision, dimension(:,:,:), intent(inout) :: gbuf
        doubleprecision, dimension(size(gbuf,1),size(gbuf,2),size(gbuf,3)) :: lbuf

        return
        end subroutine global_all_min_onevar_3d




        subroutine global_all_min_onevar_0i( gbuf )
        integer, intent(inout) :: gbuf
        integer :: lbuf

        return
        end subroutine global_all_min_onevar_0i


        subroutine global_all_min_onevar_1i( gbuf )
        integer, dimension(:), intent(inout) :: gbuf
        integer, dimension(size(gbuf)) :: lbuf

        return
        end subroutine global_all_min_onevar_1i

        subroutine global_all_min_onevar_2i( gbuf )
        integer, dimension(:,:), intent(inout) :: gbuf
        integer, dimension(size(gbuf,1),size(gbuf,2)) :: lbuf

        return
        end subroutine global_all_min_onevar_2i


        subroutine global_all_min_onevar_3i( gbuf )
        integer, dimension(:,:,:), intent(inout) :: gbuf
        integer, dimension(size(gbuf,1),size(gbuf,2),size(gbuf,3)) :: lbuf

        return
        end subroutine global_all_min_onevar_3i

        subroutine global_all_min_onevar_0r( gbuf )
        real, intent(inout) :: gbuf
        real :: lbuf

        return
        end subroutine global_all_min_onevar_0r


        subroutine global_all_min_onevar_1r( gbuf )
        real, dimension(:), intent(inout) :: gbuf
        real, dimension(size(gbuf)) :: lbuf

        return
        end subroutine global_all_min_onevar_1r

        subroutine global_all_min_onevar_2r( gbuf )
        real, dimension(:,:), intent(inout) :: gbuf
        real, dimension(size(gbuf,1),size(gbuf,2)) :: lbuf

        return
        end subroutine global_all_min_onevar_2r



        subroutine global_all_min_onevar_3r( gbuf )
        real, dimension(:,:,:), intent(inout) :: gbuf
        real, dimension(size(gbuf,1),size(gbuf,2),size(gbuf,3)) :: lbuf

        return
        end subroutine global_all_min_onevar_3r


        subroutine global_all_min_0i( lbuf, gbuf, mroot, idebug )
        integer, intent(in) :: lbuf
        integer, intent(out) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        gbuf = lbuf

	return
        end subroutine global_all_min_0i


	subroutine global_all_min_1i( lbuf, gbuf, mroot, idebug )
        integer, intent(in), dimension(:) :: lbuf       
        integer, intent(out), dimension(:) :: gbuf
	integer, optional, intent(in) :: mroot, idebug

        gbuf = lbuf

	return
	end subroutine global_all_min_1i

        subroutine global_all_min_2i( lbuf, gbuf, mroot, idebug )
        integer, intent(in), dimension(:,:) :: lbuf
        integer, intent(out), dimension(:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        gbuf = lbuf

	return
	end subroutine global_all_min_2i

        subroutine global_all_min_3i( lbuf, gbuf, mroot, idebug )
        integer, intent(in), dimension(:,:,:) :: lbuf
        integer, intent(out), dimension(:,:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        gbuf = lbuf

	return
        end subroutine global_all_min_3i

        subroutine global_all_min_0r( lbuf, gbuf, mroot, idebug )
        real, intent(in) :: lbuf
        real, intent(out) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        gbuf = lbuf

	return
        end subroutine global_all_min_0r


	subroutine global_all_min_1r( lbuf, gbuf, mroot, idebug )
        real, intent(in), dimension(:) :: lbuf       
        real, intent(out), dimension(:) :: gbuf
	integer, optional, intent(in) :: mroot, idebug

        gbuf = lbuf

	return
	end subroutine global_all_min_1r

        subroutine global_all_min_2r( lbuf, gbuf, mroot, idebug )
        real, intent(in), dimension(:,:) :: lbuf
        real, intent(out), dimension(:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        gbuf = lbuf

	return
	end subroutine global_all_min_2r

        subroutine global_all_min_3r( lbuf, gbuf, mroot, idebug )
        real, intent(in), dimension(:,:,:) :: lbuf
        real, intent(out), dimension(:,:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        gbuf = lbuf

	return
        end subroutine global_all_min_3r

        subroutine global_all_min_0d( lbuf, gbuf, mroot, idebug )
        double precision, intent(in) :: lbuf
        double precision, intent(out) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        gbuf = lbuf

	return
        end subroutine global_all_min_0d


	subroutine global_all_min_1d( lbuf, gbuf, mroot, idebug )
        double precision, intent(in), dimension(:) :: lbuf       
        double precision, intent(out), dimension(:) :: gbuf
	integer, optional, intent(in) :: mroot, idebug

        gbuf = lbuf

	return
	end subroutine global_all_min_1d

        subroutine global_all_min_2d( lbuf, gbuf, mroot, idebug )
        double precision, intent(in), dimension(:,:) :: lbuf
        double precision, intent(out), dimension(:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        gbuf = lbuf

	return
	end subroutine global_all_min_2d

        subroutine global_all_min_3d( lbuf, gbuf, mroot, idebug )
        double precision, intent(in), dimension(:,:,:) :: lbuf
        double precision, intent(out), dimension(:,:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        gbuf = lbuf

	return
        end subroutine global_all_min_3d

        subroutine global_max_0i( lbuf, gbuf, mroot, idebug )
        integer, intent(in) :: lbuf
        integer, intent(out) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        gbuf = lbuf

	return
        end subroutine global_max_0i


	subroutine global_max_1i( lbuf, gbuf, mroot, idebug )
        integer, intent(in), dimension(:) :: lbuf       
        integer, intent(out), dimension(:) :: gbuf
	integer, optional, intent(in) :: mroot, idebug

        gbuf = lbuf

	return
	end subroutine global_max_1i

        subroutine global_max_2i( lbuf, gbuf, mroot, idebug )
        integer, intent(in), dimension(:,:) :: lbuf
        integer, intent(out), dimension(:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        gbuf = lbuf

	return
	end subroutine global_max_2i

        subroutine global_max_3i( lbuf, gbuf, mroot, idebug )
        integer, intent(in), dimension(:,:,:) :: lbuf
        integer, intent(out), dimension(:,:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        gbuf = lbuf

	return
        end subroutine global_max_3i

        subroutine global_max_0r( lbuf, gbuf, mroot, idebug )
        real, intent(in) :: lbuf
        real, intent(out) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        gbuf = lbuf

	return
        end subroutine global_max_0r


	subroutine global_max_1r( lbuf, gbuf, mroot, idebug )
        real, intent(in), dimension(:) :: lbuf       
        real, intent(out), dimension(:) :: gbuf
	integer, optional, intent(in) :: mroot, idebug

        gbuf = lbuf

	return
	end subroutine global_max_1r

        subroutine global_max_2r( lbuf, gbuf, mroot, idebug )
        real, intent(in), dimension(:,:) :: lbuf
        real, intent(out), dimension(:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        gbuf = lbuf

	return
	end subroutine global_max_2r

        subroutine global_max_3r( lbuf, gbuf, mroot, idebug )
        real, intent(in), dimension(:,:,:) :: lbuf
        real, intent(out), dimension(:,:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        gbuf = lbuf

	return
        end subroutine global_max_3r

        subroutine global_max_0d( lbuf, gbuf, mroot, idebug )
        double precision, intent(in) :: lbuf
        double precision, intent(out) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        gbuf = lbuf

	return
        end subroutine global_max_0d


	subroutine global_max_1d( lbuf, gbuf, mroot, idebug )
        double precision, intent(in), dimension(:) :: lbuf       
        double precision, intent(out), dimension(:) :: gbuf
	integer, optional, intent(in) :: mroot, idebug

        gbuf = lbuf

	return
	end subroutine global_max_1d

        subroutine global_max_2d( lbuf, gbuf, mroot, idebug )
        double precision, intent(in), dimension(:,:) :: lbuf
        double precision, intent(out), dimension(:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        gbuf = lbuf

	return
	end subroutine global_max_2d

        subroutine global_max_3d( lbuf, gbuf, mroot, idebug )
        double precision, intent(in), dimension(:,:,:) :: lbuf
        double precision, intent(out), dimension(:,:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        gbuf = lbuf

	return
        end subroutine global_max_3d

        subroutine global_all_max_onevar_0d( gbuf )
        doubleprecision, intent(inout) :: gbuf
        doubleprecision :: lbuf

        return
        end subroutine global_all_max_onevar_0d


        subroutine global_all_max_onevar_1d( gbuf )
        doubleprecision, dimension(:), intent(inout) :: gbuf
        doubleprecision, dimension(size(gbuf)) :: lbuf

        return
        end subroutine global_all_max_onevar_1d

        subroutine global_all_max_onevar_2d( gbuf )
        doubleprecision, dimension(:,:), intent(inout) :: gbuf
        doubleprecision, dimension(size(gbuf,1),size(gbuf,2)) :: lbuf

        return
        end subroutine global_all_max_onevar_2d


        subroutine global_all_max_onevar_3d( gbuf )
        doubleprecision, dimension(:,:,:), intent(inout) :: gbuf
        doubleprecision, dimension(size(gbuf,1),size(gbuf,2),size(gbuf,3)) :: lbuf

        return
        end subroutine global_all_max_onevar_3d




        subroutine global_all_max_onevar_0i( gbuf )
        integer, intent(inout) :: gbuf
        integer :: lbuf

        return
        end subroutine global_all_max_onevar_0i


        subroutine global_all_max_onevar_1i( gbuf )
        integer, dimension(:), intent(inout) :: gbuf
        integer, dimension(size(gbuf)) :: lbuf

        return
        end subroutine global_all_max_onevar_1i

        subroutine global_all_max_onevar_2i( gbuf )
        integer, dimension(:,:), intent(inout) :: gbuf
        integer, dimension(size(gbuf,1),size(gbuf,2)) :: lbuf

        return
        end subroutine global_all_max_onevar_2i


        subroutine global_all_max_onevar_3i( gbuf )
        integer, dimension(:,:,:), intent(inout) :: gbuf
        integer, dimension(size(gbuf,1),size(gbuf,2),size(gbuf,3)) :: lbuf

        return
        end subroutine global_all_max_onevar_3i

        subroutine global_all_max_onevar_0r( gbuf )
        real, intent(inout) :: gbuf
        real :: lbuf

        return
        end subroutine global_all_max_onevar_0r

        subroutine global_all_max_onevar_1r( gbuf )
        real, dimension(:), intent(inout) :: gbuf
        real, dimension(size(gbuf)) :: lbuf

        return
        end subroutine global_all_max_onevar_1r

        subroutine global_all_max_onevar_2r( gbuf )
        real, dimension(:,:), intent(inout) :: gbuf
        real, dimension(size(gbuf,1),size(gbuf,2)) :: lbuf

        return
        end subroutine global_all_max_onevar_2r


        subroutine global_all_max_onevar_3r( gbuf )
        real, dimension(:,:,:), intent(inout) :: gbuf
        real, dimension(size(gbuf,1),size(gbuf,2),size(gbuf,3)) :: lbuf

        return
        end subroutine global_all_max_onevar_3r



        subroutine global_all_max_0i( lbuf, gbuf, mroot, idebug )
        integer, intent(in) :: lbuf
        integer, intent(out) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        gbuf = lbuf

	return
        end subroutine global_all_max_0i


	subroutine global_all_max_1i( lbuf, gbuf, mroot, idebug )
        integer, intent(in), dimension(:) :: lbuf       
        integer, intent(out), dimension(:) :: gbuf
	integer, optional, intent(in) :: mroot, idebug

        gbuf = lbuf

	return
	end subroutine global_all_max_1i

        subroutine global_all_max_2i( lbuf, gbuf, mroot, idebug )
        integer, intent(in), dimension(:,:) :: lbuf
        integer, intent(out), dimension(:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        gbuf = lbuf

	return
	end subroutine global_all_max_2i

        subroutine global_all_max_3i( lbuf, gbuf, mroot, idebug )
        integer, intent(in), dimension(:,:,:) :: lbuf
        integer, intent(out), dimension(:,:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        gbuf = lbuf

	return
        end subroutine global_all_max_3i

        subroutine global_all_max_0r( lbuf, gbuf, mroot, idebug )
        real, intent(in) :: lbuf
        real, intent(out) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        gbuf = lbuf

	return
        end subroutine global_all_max_0r


	subroutine global_all_max_1r( lbuf, gbuf, mroot, idebug )
        real, intent(in), dimension(:) :: lbuf       
        real, intent(out), dimension(:) :: gbuf
	integer, optional, intent(in) :: mroot, idebug

        gbuf = lbuf

	return
	end subroutine global_all_max_1r

        subroutine global_all_max_2r( lbuf, gbuf, mroot, idebug )
        real, intent(in), dimension(:,:) :: lbuf
        real, intent(out), dimension(:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        gbuf = lbuf

	return
	end subroutine global_all_max_2r

        subroutine global_all_max_3r( lbuf, gbuf, mroot, idebug )
        real, intent(in), dimension(:,:,:) :: lbuf
        real, intent(out), dimension(:,:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        gbuf = lbuf

	return
        end subroutine global_all_max_3r

        subroutine global_all_max_0d( lbuf, gbuf, mroot, idebug )
        double precision, intent(in) :: lbuf
        double precision, intent(out) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        gbuf = lbuf

	return
        end subroutine global_all_max_0d


	subroutine global_all_max_1d( lbuf, gbuf, mroot, idebug )
        double precision, intent(in), dimension(:) :: lbuf       
        double precision, intent(out), dimension(:) :: gbuf
	integer, optional, intent(in) :: mroot, idebug

        gbuf = lbuf

	return
	end subroutine global_all_max_1d

        subroutine global_all_max_2d( lbuf, gbuf, mroot, idebug )
        double precision, intent(in), dimension(:,:) :: lbuf
        double precision, intent(out), dimension(:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        gbuf = lbuf

	return
	end subroutine global_all_max_2d

        subroutine global_all_max_3d( lbuf, gbuf, mroot, idebug )
        double precision, intent(in), dimension(:,:,:) :: lbuf
        double precision, intent(out), dimension(:,:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        gbuf = lbuf

	return
        end subroutine global_all_max_3d



        subroutine global_all_or_onevar_0d( gbuf )
        logical, intent(inout) :: gbuf
        logical :: lbuf

        return
        end subroutine global_all_or_onevar_0d

        subroutine global_all_or_onevar_1d( gbuf )
        logical, dimension(:), intent(inout) :: gbuf
        logical, dimension(size(gbuf)) :: lbuf

        return
        end subroutine global_all_or_onevar_1d

        subroutine global_all_and_onevar_0d( gbuf )
        logical, intent(inout) :: gbuf
        logical :: lbuf

        return
        end subroutine global_all_and_onevar_0d

        subroutine global_all_and_onevar_1d( gbuf )
        logical, dimension(:), intent(inout) :: gbuf
        logical, dimension(size(gbuf)) :: lbuf

        return
        end subroutine global_all_and_onevar_1d


        subroutine global_all_and_0d( lvalue, gvalue, mroot, idebug )
        logical, intent(in) :: lvalue
        logical, intent(out) :: gvalue
        integer, optional, intent(in) :: mroot, idebug

        gvalue = lvalue

	return
	end subroutine  global_all_and_0d


        subroutine global_all_and_1d( lvalue, gvalue, mroot, idebug )
        logical, intent(in), dimension(:) :: lvalue
        logical, intent(out), dimension(:) :: gvalue
        integer, optional, intent(in) :: mroot, idebug

        gvalue = lvalue

	return
	end subroutine global_all_and_1d



        subroutine global_all_or_0d( lvalue, gvalue, mroot, idebug )
        logical, intent(in) :: lvalue
        logical, intent(out) :: gvalue
        integer, optional, intent(in) :: mroot, idebug

        gvalue = lvalue

	return
	end subroutine global_all_or_0d


        subroutine global_all_or_1d( lvalue, gvalue, mroot, idebug )
        logical, intent(in), dimension(:) :: lvalue
        logical, intent(out), dimension(:) :: gvalue
        integer, optional, intent(in) :: mroot, idebug

        gvalue = lvalue

	return
	end subroutine global_all_or_1d


      SUBROUTINE ExitMPI(myid)
      USE funits
      integer, optional, intent(in) :: myid

      RETURN
      END subroutine exitMPI



	subroutine MPI_BARRIER(MPI_COMM_WORLD, mpierr)
	integer, intent(in) :: MPI_COMM_WORLD, mpierr

	return
	end subroutine MPI_BARRIER

	subroutine MPI_ABORT(MPI_COMM_WORLD, mpierr, ierr)
	integer, intent(in) :: MPI_COMM_WORLD, mpierr, ierr

	STOP
	end subroutine MPI_ABORT

	end module mpi_utility


