	module sendrecv
	use parallel_mpi
	use debug
 	use geometry
 	use compar
	use indices
	implicit none

	integer, pointer, dimension(:) :: &
		recvproc1, recvtag1, xrecv1, recvijk1, &
		sendproc1, sendtag1, xsend1, sendijk1, &
		recvproc2, recvtag2, xrecv2, recvijk2, &
		sendproc2, sendtag2, xsend2, sendijk2

	integer :: nrecv1,nsend1, nrecv2,nsend2
	logical,parameter :: localfunc=.true.


	double precision, dimension(:), allocatable :: &
		dsendbuffer, drecvbuffer
	integer, dimension(:), allocatable :: &
		isendbuffer, irecvbuffer
	character, dimension(:), pointer :: &
		csendbuffer, crecvbuffer

	integer :: nrecv,nsend
	integer, pointer, dimension(:) :: &
		recvrequest, sendrequest, &
		xrecv,recvproc, recvijk, recvtag, &
		xsend,sendproc, sendijk, sendtag

	integer :: communicator

!	-----------------
!	generic interface 
!	-----------------
	interface sendrecv_begin
	module procedure &
		sendrecv_begin_1d, &
		sendrecv_begin_1i, &
		sendrecv_begin_1c 
	end interface

	interface sendrecv_end
	module procedure &
		sendrecv_end_1d, &
		sendrecv_end_1i, &
		sendrecv_end_1c  
	end interface

        interface send_recv
        module procedure &
                send_recv_1d, send_recv_2d, send_recv_3d, &
                send_recv_1i, &
                send_recv_1c
        end interface


	contains

	subroutine ijk_of( ijkp, i,j,k )
	integer, intent(in) :: ijkp
	integer, intent(out) :: i,j,k

	integer :: k1,k2, j1,j2, i1,i2, &
		ijk, isize,jsize,ksize, gijk

	character(len=32), parameter :: name = "ijk_of"
	logical :: isok_k, isok_j, isok_i, is_same, isok

	ijk = ijkp

	i1 = istart3_all(myPE)
	i2 = iend3_all(myPE)
	j1 = jstart3_all(myPE)
	j2 = jend3_all(myPE)
	k1 = kstart3_all(myPE)
	k2 = kend3_all(myPE)

	ksize = (k2-k1+1)
	jsize = (j2-j1+1)
	isize = (i2-i1+1)

	
	if (mod(ijk,isize*jsize).ne.0) then
		k = int( ijk/(isize*jsize) ) + k1
	else
		k = int( ijk/(isize*jsize) ) + k1 -1
	endif
	ijk = ijk - (k-k1)*(isize*jsize)

	if (mod(ijk,isize).ne.0) then
		j = int( ijk/isize ) + j1
	else
		j = int( ijk/isize ) + j1 - 1
	endif
	ijk = ijk - (j-j1)*isize

	i = (ijk-1) + i1
!	------------
!	double check
!	------------
	isok_i = (i1 <= i) .and. (i <= i2)
	isok_j = (j1 <= j) .and. (j <= j2)
	isok_k = (k1 <= k) .and. (k <= k2)
	gijk = 1 + (i-i1) + (j-j1)*(i2-i1+1) + &
		(k-k1)*(j2-j1+1)*(i2-i1+1)
	is_same = (gijk .eq. ijkp)
	isok = isok_i .and. isok_j .and. isok_k .and. is_same
	if (.not.isok) then
	    call write_debug( name, 'i,j,k ', i,j,k )
	    call write_debug( name, 'ijkp, gijk ', ijkp, gijk )
	endif


	return
	end subroutine ijk_of


	subroutine ijk_of_gl( ijkp, i,j,k )
	integer, intent(in) :: ijkp
	integer, intent(out) :: i,j,k

	integer :: k1,k2, j1,j2, i1,i2, &
		ijk, isize,jsize,ksize, gijk

	character(len=32), parameter :: name = "ijk_of_gl"
	logical :: isok_k, isok_j, isok_i, is_same, isok

	ijk = ijkp

	k1 = minval( kstart3_all(:) )
	k2 = maxval( kend3_all(:) )

	j1 = minval( jstart3_all(:) )
	j2 = maxval( jend3_all(:) )

	i1 = minval( istart3_all(:) )
	i2 = maxval( iend3_all(:) )

	ksize = (k2-k1+1)
	jsize = (j2-j1+1)
	isize = (i2-i1+1)

	
	if (mod(ijk,isize*jsize).ne.0) then
		k = int( ijk/(isize*jsize) ) + k1
	else
		k = int( ijk/(isize*jsize) ) + k1 -1
	endif
	ijk = ijk - (k-k1)*(isize*jsize)

	if (mod(ijk,isize).ne.0) then
		j = int( ijk/isize ) + j1
	else
		j = int( ijk/isize ) + j1 - 1
	endif
	ijk = ijk - (j-j1)*isize

	i = (ijk-1) + i1
!	------------
!	double check
!	------------
	isok_i = (i1 <= i) .and. (i <= i2)
	isok_j = (j1 <= j) .and. (j <= j2)
	isok_k = (k1 <= k) .and. (k <= k2)
	gijk = 1 + (i-i1) + (j-j1)*(i2-i1+1) + &
		(k-k1)*(j2-j1+1)*(i2-i1+1)
	is_same = (gijk .eq. ijkp)
	isok = isok_i .and. isok_j .and. isok_k .and. is_same
	if (.not.isok) then
	    call write_debug( name, 'i,j,k ', i,j,k )
	    call write_debug( name, 'ijkp, gijk ', ijkp, gijk )
	endif


	return
	end subroutine ijk_of_gl

	subroutine sendrecv_init( 	 &
		comm,			 &
		cyclic_i,cyclic_j,cyclic_k, idebug )
        implicit none

	integer, intent(in) :: comm
	logical,intent(in) :: cyclic_i,cyclic_j,cyclic_k

	integer, intent(in), optional :: idebug

	return
	end subroutine sendrecv_init



	subroutine sendrecv_begin_1d( X, ilayer, idebug )
        implicit none

	integer, intent(in),optional :: ilayer
	double precision, intent(inout), dimension(:) :: X
	integer, intent(in), optional :: idebug

	return
	end subroutine sendrecv_begin_1d


	subroutine sendrecv_begin_1i( X, ilayer, idebug )
        implicit none

	integer, intent(in),optional :: ilayer
	integer, intent(inout), dimension(:) :: X
	integer, intent(in), optional :: idebug

	return
	end subroutine sendrecv_begin_1i


	subroutine sendrecv_begin_1c( X, ilayer, idebug )
        implicit none

	integer, intent(in),optional :: ilayer
	character(len=*), intent(inout), dimension(:) :: X
	integer, intent(in), optional :: idebug

	return
	end subroutine sendrecv_begin_1c


	subroutine sendrecv_end_1d( X, idebug )
	double precision, intent(inout), dimension(:) :: X
	integer, intent(in), optional :: idebug
	return
	end subroutine sendrecv_end_1d


	subroutine sendrecv_end_1c( X, idebug )
	character(len=*), intent(inout), dimension(:) :: X
	integer, intent(in), optional :: idebug
	return
	end subroutine sendrecv_end_1c


	subroutine sendrecv_end_1i( X, idebug )
	integer, intent(inout), dimension(:) :: X
	integer, intent(in), optional :: idebug
	return
	end subroutine sendrecv_end_1i


	subroutine send_recv_1c( X, ilayer, idebug )
	character(len=*),  dimension(:), intent(inout) :: X
	integer, intent(in), optional :: ilayer,idebug
	return
	end subroutine send_recv_1c

	subroutine send_recv_1d( X, ilayer, idebug )
	double precision,  dimension(:), intent(inout) :: X
	integer, intent(in), optional :: ilayer,idebug
	return
	end subroutine send_recv_1d

        subroutine send_recv_2d( X, ilayer, idebug )
        double precision,  dimension(:,:), intent(inout) :: X
        integer, intent(in), optional :: ilayer,idebug
	return
	end subroutine send_recv_2d

        subroutine send_recv_3d( X, ilayer, idebug )
        double precision,  dimension(:,:,:), intent(inout) :: X
        integer, intent(in), optional :: ilayer,idebug
	return
	end subroutine send_recv_3d

	subroutine send_recv_1i( X, ilayer, idebug )
	integer,  dimension(:), intent(inout) :: X
	integer, intent(in), optional :: ilayer,idebug
	return
	end subroutine send_recv_1i

	
	end module sendrecv
