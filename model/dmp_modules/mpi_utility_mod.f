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
 	module  procedure global_all_sum_0i, global_all_sum_1i, &
                          global_all_sum_2i, global_all_sum_3i, &
                          global_all_sum_0r, global_all_sum_1r, &
                          global_all_sum_2r, global_all_sum_3r, &
                          global_all_sum_0d, global_all_sum_1d, &
                          global_all_sum_2d, global_all_sum_3d
 	end interface 

	interface global_min
	module  procedure global_min_0i, global_min_1i, global_min_2i, global_min_3i, &
                          global_min_0r, global_min_1r, global_min_2r, global_min_3r, &
                          global_min_0d, global_min_1d, global_min_2d, global_min_3d
	end interface 

 	interface global_all_min
 	module  procedure global_all_min_0i, global_all_min_1i, &
                          global_all_min_2i, global_all_min_3i, &
                          global_all_min_0r, global_all_min_1r, &
                          global_all_min_2r, global_all_min_3r, &
                          global_all_min_0d, global_all_min_1d, &
                          global_all_min_2d, global_all_min_3d
 	end interface 

	interface global_max
	module  procedure global_max_0i, global_max_1i, global_max_2i, global_max_3i, &
                          global_max_0r, global_max_1r, global_max_2r, global_max_3r, &
                          global_max_0d, global_max_1d, global_max_2d, global_max_3d
	end interface 

 	interface global_all_max
 	module  procedure global_all_max_0i, global_all_max_1i, &
                          global_all_max_2i, global_all_max_3i, &
                          global_all_max_0r, global_all_max_1r, &
                          global_all_max_2r, global_all_max_3r, &
                          global_all_max_0d, global_all_max_1d, &
                          global_all_max_2d, global_all_max_3d
 	end interface 

        interface global_all_and
        module procedure global_all_and_0d, global_all_and_1d
        end interface

        interface global_all_or
        module procedure global_all_or_0d, global_all_or_1d
        end interface

	contains


!	Routine to scatter gbuf available on root to all the processors

	subroutine scatter_1i( lbuf, gbuf, mroot, idebug )
        integer, intent(in), dimension(:) :: gbuf
        integer, intent(out), dimension(:) :: lbuf
	integer, optional, intent(in) :: mroot, idebug

	integer, allocatable, dimension(:) :: gbuf_pack

	integer :: sendtype, recvtype, ijk1, ijk2, recvcnt, ierr,lroot, lidebug
	integer :: i,j,k,ibuffer,iproc
        integer :: ijk
        include 'function.inc'

!	check to see whether there is root

	if (.not. present(mroot)) then
	   lroot = 0
	else
	   lroot = mroot
	endif

        if (.not. present(idebug)) then
           lidebug = 0
        else
           lidebug = idebug
        endif

	if(myPE.eq.lroot) then
	allocate(gbuf_pack(sum(ijksize3_all(:))))
	else
	allocate(gbuf_pack(10))
	endif

	if( myPE.eq.lroot) then
        ibuffer = 0
          do iproc = 0,numPEs-1
            do k = kstart3_all(iproc), kend3_all(iproc)
              do j = jstart3_all(iproc), jend3_all(iproc)
                do i = istart3_all(iproc), iend3_all(iproc)

                  ibuffer = ibuffer + 1
                  gbuf_pack(ibuffer) = gbuf(funijk_gl(i,j,k))
                
                enddo
              enddo
            enddo
          enddo
	endif 

	sendtype = MPI_INTEGER
	recvtype = sendtype

	ijk1 = ijkstart3
	ijk2 = ijkend3

	recvcnt = ijk2-ijk1+1

!	Call MPI routines

	call MPI_Scatterv( gbuf_pack, ijksize3_all, displs, sendtype, &
			  lbuf, recvcnt, recvtype,  &
			  lroot, MPI_COMM_WORLD, ierr )
	call MPI_Check( 'scatter_1i:MPI_Scatterv', ierr )

	deallocate(gbuf_pack)

	return
	end subroutine scatter_1i

        subroutine scatter_2i( lbuf, gbuf, mroot, idebug )
        integer, intent(in), dimension(:,:) :: gbuf       
        integer, intent(out), dimension(:,:) :: lbuf
        integer, optional, intent(in) :: mroot, idebug

	integer :: i,j,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
	   lroot = mroot
        endif

        if (.not. present(idebug)) then 
           lidebug = 0
        else
           lidebug = idebug
        endif

	if(myPE.eq.lroot) then
	call assert( size(lbuf,2).eq.size(gbuf,2),  &
		'** scatter_2i: size(lbuf,2).ne.size(gbuf,2) ', &
		size(lbuf,2), size(gbuf,2) )
	endif

	do j=lbound(lbuf,2),ubound(lbuf,2)
	  call scatter_1i( lbuf(:,j), gbuf(:,j), lroot, lidebug )
	enddo

	return
	end subroutine scatter_2i

        subroutine scatter_3i( lbuf, gbuf, mroot, idebug )
        integer, intent(in), dimension(:,:,:) :: gbuf       
        integer, intent(out), dimension(:,:,:) :: lbuf
        integer, optional, intent(in) :: mroot, idebug

        integer :: j,k,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

	if(myPE.eq.lroot) then
        call assert( size(lbuf,2).eq.size(gbuf,2),  &
                '** scatter_3i: size(lbuf,2).ne.size(gbuf,2) ', &
                size(lbuf,2), size(gbuf,2) )

        call assert( size(lbuf,3).eq.size(gbuf,3),  &
                '** scatter_3i: size(lbuf,3).ne.size(gbuf,3) ', &
                size(lbuf,3), size(gbuf,3) )
	endif

        do k=lbound(lbuf,3),ubound(lbuf,3)
        do j=lbound(lbuf,2),ubound(lbuf,2)
          call scatter_1i( lbuf(:,j,k), gbuf(:,j,k), lroot, lidebug )
        enddo
	enddo

        return
        end subroutine scatter_3i

	subroutine scatter_1r( lbuf, gbuf, mroot, idebug )
        real, intent(in), dimension(:) :: gbuf       
        real, intent(out), dimension(:) :: lbuf
	integer, optional, intent(in) :: mroot, idebug

        real, allocatable, dimension(:) :: gbuf_pack

	integer :: sendtype, recvtype, ijk1, ijk2, recvcnt, ierr,lroot, lidebug
        integer :: i,j,k,ibuffer,iproc
        integer :: ijk
        include 'function.inc'

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

        if(myPE.eq.lroot) then
	allocate(gbuf_pack(sum(ijksize3_all(:))))
	else
	allocate(gbuf_pack(10))
	endif

        if( myPE.eq.lroot) then
        ibuffer = 0
          do iproc = 0,numPEs-1
            do k = kstart3_all(iproc), kend3_all(iproc)
              do j = jstart3_all(iproc), jend3_all(iproc)
                do i = istart3_all(iproc), iend3_all(iproc)

                  ibuffer = ibuffer + 1
                  gbuf_pack(ibuffer) = gbuf(funijk_gl(i,j,k))

                enddo
              enddo
            enddo
          enddo
        endif

	sendtype = MPI_REAL
	recvtype = sendtype

        ijk1 = ijkstart3
        ijk2 = ijkend3

	recvcnt = ijk2-ijk1+1

	call MPI_Scatterv( gbuf_pack, ijksize3_all, displs, sendtype, &
			  lbuf, recvcnt, recvtype,  &
			  lroot, MPI_COMM_WORLD, ierr )
	call MPI_Check( 'scatter_1r:MPI_Scatterv', ierr )

        deallocate(gbuf_pack)

	return
	end subroutine scatter_1r

	
        subroutine scatter_2r( lbuf, gbuf, mroot, idebug )
        real, intent(in), dimension(:,:) :: gbuf       
        real, intent(out), dimension(:,:) :: lbuf
        integer, optional, intent(in) :: mroot, idebug

	integer :: i,j,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

	if(myPE.eq.lroot) then
	call assert( size(lbuf,2).eq.size(gbuf,2),  &
		'** scatter_2r: size(lbuf,2).ne.size(gbuf,2) ', &
		size(lbuf,2), size(gbuf,2) )
	endif

	do j=lbound(lbuf,2),ubound(lbuf,2)
	  call scatter_1r( lbuf(:,j), gbuf(:,j), lroot, lidebug )
	enddo

	return
	end subroutine scatter_2r

        subroutine scatter_3r( lbuf, gbuf, mroot, idebug )
        real, intent(in), dimension(:,:,:) :: gbuf       
        real, intent(out), dimension(:,:,:) :: lbuf
        integer, optional, intent(in) :: mroot, idebug

        integer :: j,k,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

	if(myPE.eq.lroot) then
        call assert( size(lbuf,2).eq.size(gbuf,2),  &
                '** scatter_3r: size(lbuf,2).ne.size(gbuf,2) ', &
                size(lbuf,2), size(gbuf,2) )

        call assert( size(lbuf,3).eq.size(gbuf,3),  &
                '** scatter_3r: size(lbuf,3).ne.size(gbuf,3) ', &
                size(lbuf,3), size(gbuf,3) )
	endif

        do k=lbound(lbuf,3),ubound(lbuf,3)
        do j=lbound(lbuf,2),ubound(lbuf,2)
          call scatter_1r( lbuf(:,j,k), gbuf(:,j,k), lroot, lidebug )
        enddo
	enddo

        return
        end subroutine scatter_3r


	subroutine scatter_1d( lbuf, gbuf, mroot, idebug )
        double precision, intent(in), dimension(:) :: gbuf       
        double precision, intent(out), dimension(:) :: lbuf
	integer, optional, intent(in) :: mroot, idebug

        double precision, allocatable, dimension(:) :: gbuf_pack

	integer :: sendtype, recvtype, ijk1,ijk2,recvcnt, ierr,lroot, lidebug
        integer :: i,j,k,ibuffer,iproc
        integer :: ijk
        include 'function.inc'

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

        if(myPE.eq.lroot) then
	allocate(gbuf_pack(sum(ijksize3_all(:))))
	else
	allocate(gbuf_pack(10))
	endif

        if( myPE.eq.lroot) then
        ibuffer = 0
          do iproc = 0,numPEs-1
            do k = kstart3_all(iproc), kend3_all(iproc)
              do j = jstart3_all(iproc), jend3_all(iproc)
                do i = istart3_all(iproc), iend3_all(iproc)

                  ibuffer = ibuffer + 1
                  gbuf_pack(ibuffer) = gbuf(funijk_gl(i,j,k))

                enddo
              enddo
            enddo
          enddo
        endif

	sendtype = MPI_DOUBLE_PRECISION
	recvtype = sendtype

        ijk1 = ijkstart3
        ijk2 = ijkend3

	recvcnt = ijk2-ijk1+1

	call MPI_Scatterv( gbuf_pack, ijksize3_all, displs, sendtype, &
			  lbuf, recvcnt, recvtype,  &
			  lroot, MPI_COMM_WORLD, ierr )
	call MPI_Check( 'scatter_1d:MPI_Scatterv', ierr )

        deallocate(gbuf_pack)

	return
	end subroutine scatter_1d

	
        subroutine scatter_2d( lbuf, gbuf, mroot, idebug )
        double precision, intent(in), dimension(:,:) :: gbuf       
        double precision, intent(out), dimension(:,:) :: lbuf
        integer, optional, intent(in) :: mroot, idebug

	integer :: i,j,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

	if(myPE.eq.lroot) then
	call assert( size(lbuf,2).eq.size(gbuf,2),  &
		'** scatter_2d: size(lbuf,2).ne.size(gbuf,2) ', &
		size(lbuf,2), size(gbuf,2) )
	endif

	do j=lbound(lbuf,2),ubound(lbuf,2)
	  call scatter_1d( lbuf(:,j), gbuf(:,j), lroot, lidebug )
	enddo

	return
	end subroutine scatter_2d

        subroutine scatter_3d( lbuf, gbuf, mroot, idebug )
        double precision, intent(in), dimension(:,:,:) :: gbuf       
        double precision, intent(out), dimension(:,:,:) :: lbuf
        integer, optional, intent(in) :: mroot, idebug

        integer :: j,k,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

	if(myPE.eq.lroot) then
        call assert( size(lbuf,2).eq.size(gbuf,2),  &
                '** scatter_3d: size(lbuf,2).ne.size(gbuf,2) ', &
                size(lbuf,2), size(gbuf,2) )

        call assert( size(lbuf,3).eq.size(gbuf,3),  &
                '** scatter_3d: size(lbuf,3).ne.size(gbuf,3) ', &
                size(lbuf,3), size(gbuf,3) )
	endif

        do k=lbound(lbuf,3),ubound(lbuf,3)
        do j=lbound(lbuf,2),ubound(lbuf,2)
          call scatter_1d( lbuf(:,j,k), gbuf(:,j,k), lroot, lidebug )
        enddo
	enddo

        return
        end subroutine scatter_3d


        subroutine scatter_1c( lbuf, gbuf, mroot, idebug )
        character(len=*), intent(in), dimension(:) :: gbuf
        character(len=*), intent(out), dimension(:) :: lbuf
        integer, optional, intent(in) :: mroot, idebug

        character(len=len(lbuf(1))), allocatable, dimension(:) :: gbuf_pack
        character, allocatable, dimension(:) :: gbuf_pack1,lbuf1
	character(len=len(lbuf(1))) :: string

        integer :: sendtype, recvtype, ijk1, ijk2, recvcnt, ierr,lroot, lidebug
        integer :: i,j,k,ibuffer,iproc
        integer :: ijk
        integer :: lenchar, icount
        include 'function.inc'

!       check to see whether there is root

        if (.not. present(mroot)) then
           lroot = 0
        else
           lroot = mroot
        endif

        if (.not. present(idebug)) then
           lidebug = 0
        else
           lidebug = idebug
        endif

	lenchar = len(lbuf(1))

        if(myPE.eq.lroot) then
	  allocate(gbuf_pack(sum(ijksize3_all(:))))
          allocate(gbuf_pack1(sum(ijksize3_all(:))*lenchar))
	else
	  allocate(gbuf_pack(10))
	  allocate(gbuf_pack1(10))
	endif

        allocate(lbuf1(lenchar*size(lbuf)))

        if( myPE.eq.lroot) then
        ibuffer = 0
          do iproc = 0,numPEs-1
            do k = kstart3_all(iproc), kend3_all(iproc)
              do j = jstart3_all(iproc), jend3_all(iproc)
                do i = istart3_all(iproc), iend3_all(iproc)

                  ibuffer = ibuffer + 1
                  gbuf_pack(ibuffer)(1:lenchar) = gbuf(funijk_gl(i,j,k))(1:lenchar)

                enddo
              enddo
            enddo
          enddo
        endif

	if(myPE.eq.lroot) then
	icount = 0
	do i = 1,size(gbuf_pack)
	  do j = 1,lenchar

	    icount = icount+1
	    string = gbuf_pack(i)(1:lenchar)
	    gbuf_pack1(icount) = string(j:j)

	  enddo
	enddo
	endif

        sendtype = MPI_CHARACTER
        recvtype = sendtype

        ijk1 = ijkstart3
        ijk2 = ijkend3

        recvcnt = ijk2-ijk1+1

!       Call MPI routines

        call MPI_Scatterv( gbuf_pack1, ijksize3_all*lenchar, displs*lenchar, sendtype, &
                          lbuf1, recvcnt*lenchar, recvtype,  &
                          lroot, MPI_COMM_WORLD, ierr )
        call MPI_Check( 'scatter_1c:MPI_Scatterv', ierr )

        icount = 0
        do i = 1,size(lbuf)
          do j = 1,lenchar

            icount = icount+1
            lbuf(i)(j:j) = lbuf1(icount)

          enddo
        enddo

        deallocate(gbuf_pack)
        deallocate(gbuf_pack1)
        deallocate(lbuf1)

        return
        end subroutine scatter_1c


        subroutine scatter_1l( lbuf, gbuf, mroot, idebug )
        logical, intent(in), dimension(:) :: gbuf
        logical, intent(out), dimension(:) :: lbuf
        integer, optional, intent(in) :: mroot, idebug

        logical, allocatable, dimension(:) :: gbuf_pack

        integer :: sendtype, recvtype, ijk1, ijk2, recvcnt, ierr,lroot, lidebug
        integer :: i,j,k,ibuffer,iproc
        integer :: ijk
        include 'function.inc'

!       check to see whether there is root

        if (.not. present(mroot)) then
           lroot = 0
        else
           lroot = mroot
        endif

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

        if(myPE.eq.lroot) then
	allocate(gbuf_pack(sum(ijksize3_all(:))))
	else
	allocate(gbuf_pack(10))
	endif

        if( myPE.eq.lroot) then
        ibuffer = 0
          do iproc = 0,numPEs-1
            do k = kstart3_all(iproc), kend3_all(iproc)
              do j = jstart3_all(iproc), jend3_all(iproc)
                do i = istart3_all(iproc), iend3_all(iproc)

                  ibuffer = ibuffer + 1
                  gbuf_pack(ibuffer) = gbuf(funijk_gl(i,j,k))

                enddo
              enddo
            enddo
          enddo
        endif

        sendtype = MPI_LOGICAL
        recvtype = sendtype

        ijk1 = ijkstart3
        ijk2 = ijkend3

        recvcnt = ijk2-ijk1+1

!       Call MPI routines

        call MPI_Scatterv( gbuf_pack, ijksize3_all, displs, sendtype, &
                          lbuf, recvcnt, recvtype,  &
                          lroot, MPI_COMM_WORLD, ierr )
        call MPI_Check( 'scatter_1l:MPI_Scatterv', ierr )

        deallocate(gbuf_pack)

        return
        end subroutine scatter_1l


!	Routines to gather lbuf from individual processors and put it on
!	processor root in gbuf
!	Logic is similar to the scatter routines above.

	subroutine gather_1i( lbuf, gbuf, mroot, idebug )
        integer, intent(in), dimension(:) :: lbuf       
        integer, intent(out), dimension(:) :: gbuf
	integer, optional, intent(in) :: mroot, idebug

        integer, allocatable, dimension(:) :: gbuf_pack

	integer :: recvtype, sendtype, ijk1,ijk2,sendcnt, ierr,lroot, lidebug
        integer :: i,j,k,ibuffer,iproc
        integer :: ijk, ijk_gl
	logical :: isok_k,isok_j,isok_i, isinterior 
	logical :: isbc_k,isbc_j,isbc_i, isboundary, need_copy
        include 'function.inc'

!       check to see whether there is root

        if (.not. present(mroot)) then
           lroot = 0
        else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

        if(myPE.eq.lroot) then
	  allocate(gbuf_pack(sum(ijksize3_all(:))))
	else
	  allocate(gbuf_pack(10))
        endif

	recvtype = MPI_INTEGER
	sendtype = recvtype

        ijk1 = ijkstart3
        ijk2 = ijkend3

	sendcnt = ijk2-ijk1+1

        call MPI_Gatherv( lbuf, sendcnt, sendtype,  &
                           gbuf_pack, ijksize3_all, displs, recvtype, &
                          lroot, MPI_COMM_WORLD, ierr )
	call MPI_Check( 'gather_1i:MPI_Gatherv', ierr )

        if( myPE.eq.lroot) then
        ibuffer = 0
          do iproc = 0,numPEs-1
            do k = kstart3_all(iproc), kend3_all(iproc)
              do j = jstart3_all(iproc), jend3_all(iproc)
                do i = istart3_all(iproc), iend3_all(iproc)

                  ibuffer = ibuffer + 1
                  isok_k = (kstart1_all(iproc) <= k) .and. (k <=kend1_all(iproc))
                  isok_j = (jstart1_all(iproc) <= j) .and. (j <=jend1_all(iproc))
                  isok_k = (istart1_all(iproc) <= i) .and. (i <=iend1_all(iproc))

                  isinterior = isok_k .and. isok_j .and. isok_i

                  isbc_k = (k <= kmin1) .or. (k >= kmax1)
                  isbc_j = (j <= jmin1) .or. (j >= jmax1)
                  isbc_i = (i <= imin1) .or. (i >= imax1)

                  isboundary = isbc_k .or. isbc_j .or. isbc_i

                  need_copy = isinterior .or. isboundary
                  if (need_copy) then
                      ijk_gl = funijk_gl(i,j,k)
                      gbuf( ijk_gl ) = gbuf_pack(ibuffer)
                  endif

                enddo
              enddo
            enddo
          enddo
        endif

        deallocate(gbuf_pack)

	return
	end subroutine gather_1i

	
        subroutine gather_2i( lbuf, gbuf, mroot, idebug )
        integer, intent(in), dimension(:,:) :: lbuf
        integer, intent(out), dimension(:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

	integer :: i,j,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

	if(myPE.eq.lroot) then
	call assert( size(lbuf,2).eq.size(gbuf,2),  &
		'** gather_2i: size(lbuf,2).ne.size(gbuf,2) ', &
		size(lbuf,2), size(gbuf,2) )
	endif

	do j=lbound(lbuf,2),ubound(lbuf,2)
	  call gather_1i( lbuf(:,j), gbuf(:,j), lroot, lidebug )
	enddo

	return
	end subroutine gather_2i

        subroutine gather_3i( lbuf, gbuf, mroot, idebug )
        integer, intent(in), dimension(:,:,:) :: lbuf
        integer, intent(out), dimension(:,:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        integer :: j,k,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

	if(myPE.eq.lroot) then
        call assert( size(lbuf,2).eq.size(gbuf,2),  &
                '** gather_3i: size(lbuf,2).ne.size(gbuf,2) ', &
                size(lbuf,2), size(gbuf,2) )

        call assert( size(lbuf,3).eq.size(gbuf,3),  &
                '** gather_3i: size(lbuf,3).ne.size(gbuf,3) ', &
                size(lbuf,3), size(gbuf,3) )
	endif

        do k=lbound(lbuf,3),ubound(lbuf,3)
        do j=lbound(lbuf,2),ubound(lbuf,2)
          call gather_1i( lbuf(:,j,k), gbuf(:,j,k), lroot, lidebug )
        enddo
	enddo

        return
        end subroutine gather_3i

	subroutine gather_1r( lbuf, gbuf, mroot, idebug )
        real, intent(in), dimension(:) :: lbuf
        real, intent(out), dimension(:) :: gbuf
	integer, optional, intent(in) :: mroot, idebug

        real, allocatable, dimension(:) :: gbuf_pack

	integer :: recvtype, sendtype, ijk1,ijk2,sendcnt, ierr,lroot, lidebug
        integer :: i,j,k,ibuffer,iproc
        integer :: ijk, ijk_gl
        logical :: isok_k,isok_j,isok_i, isinterior
        logical :: isbc_k,isbc_j,isbc_i, isboundary, need_copy
        include 'function.inc'

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

        if(myPE.eq.lroot) then
	  allocate(gbuf_pack(sum(ijksize3_all(:))))
	else
	  allocate(gbuf_pack(10))
	endif

	recvtype = MPI_REAL
	sendtype = recvtype

        ijk1 = ijkstart3
        ijk2 = ijkend3

	sendcnt = ijk2-ijk1+1

        call MPI_Gatherv( lbuf, sendcnt, sendtype,  &
                           gbuf_pack, ijksize3_all, displs, recvtype, &
                          lroot, MPI_COMM_WORLD, ierr )
	call MPI_Check( 'gather_1r:MPI_Gatherv', ierr )

        if( myPE.eq.lroot) then
        ibuffer = 0
          do iproc = 0,numPEs-1
            do k = kstart3_all(iproc), kend3_all(iproc)
              do j = jstart3_all(iproc), jend3_all(iproc)
                do i = istart3_all(iproc), iend3_all(iproc)

                  ibuffer = ibuffer + 1
                  isok_k = (kstart1_all(iproc) <= k) .and. (k <=kend1_all(iproc))
                  isok_j = (jstart1_all(iproc) <= j) .and. (j <=jend1_all(iproc))
                  isok_k = (istart1_all(iproc) <= i) .and. (i <=iend1_all(iproc))

                  isinterior = isok_k .and. isok_j .and. isok_i

                  isbc_k = (k <= kmin1) .or. (k >= kmax1)
                  isbc_j = (j <= jmin1) .or. (j >= jmax1)
                  isbc_i = (i <= imin1) .or. (i >= imax1)

                  isboundary = isbc_k .or. isbc_j .or. isbc_i

                  need_copy = isinterior .or. isboundary
                  if (need_copy) then
                      ijk_gl = funijk_gl(i,j,k)
                      gbuf( ijk_gl ) = gbuf_pack(ibuffer)
                  endif

                enddo
              enddo
            enddo
          enddo
        endif

        deallocate(gbuf_pack)

	return
	end subroutine gather_1r

	
        subroutine gather_2r( lbuf, gbuf, mroot, idebug )
        real, intent(in), dimension(:,:) :: lbuf
        real, intent(out), dimension(:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

	integer :: i,j,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

	if(myPE.eq.lroot) then
	call assert( size(lbuf,2).eq.size(gbuf,2),  &
		'** gather_2r: size(lbuf,2).ne.size(gbuf,2) ', &
		size(lbuf,2), size(gbuf,2) )
	endif

	do j=lbound(lbuf,2),ubound(lbuf,2)
	  call gather_1r( lbuf(:,j), gbuf(:,j), lroot, lidebug )
	enddo

	return
	end subroutine gather_2r

        subroutine gather_3r( lbuf, gbuf, mroot, idebug )
        real, intent(in), dimension(:,:,:) :: lbuf
        real, intent(out), dimension(:,:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        integer :: j,k,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

	if(myPE.eq.lroot) then
        call assert( size(lbuf,2).eq.size(gbuf,2),  &
                '** gather_3r: size(lbuf,2).ne.size(gbuf,2) ', &
                size(lbuf,2), size(gbuf,2) )

        call assert( size(lbuf,3).eq.size(gbuf,3),  &
                '** gather_3r: size(lbuf,3).ne.size(gbuf,3) ', &
                size(lbuf,3), size(gbuf,3) )
	endif

        do k=lbound(lbuf,3),ubound(lbuf,3)
        do j=lbound(lbuf,2),ubound(lbuf,2)
          call gather_1r( lbuf(:,j,k), gbuf(:,j,k), lroot, lidebug )
        enddo
	enddo

        return
        end subroutine gather_3r


	subroutine gather_1d( lbuf, gbuf, mroot, idebug )
        double precision, intent(in), dimension(:) :: lbuf
        double precision, intent(out), dimension(:) :: gbuf
	integer, optional, intent(in) :: mroot, idebug

        double precision, allocatable, dimension(:) :: gbuf_pack

	integer :: recvtype, sendtype, ijk1,ijk2,sendcnt, ierr,lroot, lidebug
        integer :: i,j,k,ibuffer,iproc
        integer :: ijk, ijk_gl
        logical :: isok_k,isok_j,isok_i, isinterior
        logical :: isbc_k,isbc_j,isbc_i, isboundary, need_copy
        include 'function.inc'

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

        if(myPE.eq.lroot) then
	  allocate(gbuf_pack(sum(ijksize3_all(:))))
	else
	  allocate(gbuf_pack(10))
	endif

	recvtype = MPI_DOUBLE_PRECISION
	sendtype = recvtype

        ijk1 = ijkstart3
        ijk2 = ijkend3

	sendcnt = ijk2-ijk1+1

        call MPI_Gatherv( lbuf, sendcnt, sendtype,  &
                           gbuf_pack, ijksize3_all, displs, recvtype, &
                          lroot, MPI_COMM_WORLD, ierr )
	call MPI_Check( 'gather_1d:MPI_Gatherv', ierr )

        if( myPE.eq.lroot) then
        ibuffer = 0
          do iproc = 0,numPEs-1
            do k = kstart3_all(iproc), kend3_all(iproc)
              do j = jstart3_all(iproc), jend3_all(iproc)
                do i = istart3_all(iproc), iend3_all(iproc)

                  ibuffer = ibuffer + 1
                  isok_k = (kstart1_all(iproc) <= k) .and. (k <=kend1_all(iproc))
                  isok_j = (jstart1_all(iproc) <= j) .and. (j <=jend1_all(iproc))
                  isok_k = (istart1_all(iproc) <= i) .and. (i <=iend1_all(iproc))

                  isinterior = isok_k .and. isok_j .and. isok_i

                  isbc_k = (k <= kmin1) .or. (k >= kmax1)
                  isbc_j = (j <= jmin1) .or. (j >= jmax1)
                  isbc_i = (i <= imin1) .or. (i >= imax1)

                  isboundary = isbc_k .or. isbc_j .or. isbc_i

                  need_copy = isinterior .or. isboundary
                  if (need_copy) then
                      ijk_gl = funijk_gl(i,j,k)
                      gbuf( ijk_gl ) = gbuf_pack(ibuffer)
                  endif

                enddo
              enddo
            enddo
          enddo
        endif

        deallocate(gbuf_pack)

	return
	end subroutine gather_1d

	
        subroutine gather_2d( lbuf, gbuf, mroot, idebug )
        double precision, intent(in), dimension(:,:) :: lbuf
        double precision, intent(out), dimension(:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

	integer :: i,j,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

	if(myPE.eq.lroot) then
	call assert( size(lbuf,2).eq.size(gbuf,2),  &
		'** gather_2d: size(lbuf,2).ne.size(gbuf,2) ', &
		size(lbuf,2), size(gbuf,2) )
	endif

	do j=lbound(lbuf,2),ubound(lbuf,2)
	  call gather_1d( lbuf(:,j), gbuf(:,j), lroot, lidebug )
	enddo

	return
	end subroutine gather_2d

        subroutine gather_3d( lbuf, gbuf, mroot, idebug )
        double precision, intent(in), dimension(:,:,:) :: lbuf
        double precision, intent(out), dimension(:,:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        integer :: j,k,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

	if(myPE.eq.lroot) then
        call assert( size(lbuf,2).eq.size(gbuf,2),  &
                '** gather_3d: size(lbuf,2).ne.size(gbuf,2) ', &
                size(lbuf,2), size(gbuf,2) )

        call assert( size(lbuf,3).eq.size(gbuf,3),  &
                '** gather_3d: size(lbuf,3).ne.size(gbuf,3) ', &
                size(lbuf,3), size(gbuf,3) )
	endif

        do k=lbound(lbuf,3),ubound(lbuf,3)
        do j=lbound(lbuf,2),ubound(lbuf,2)
          call gather_1d( lbuf(:,j,k), gbuf(:,j,k), lroot, lidebug )
        enddo
	enddo

        return
        end subroutine gather_3d


        subroutine gather_1c( lbuf, gbuf, mroot, idebug )
        character(len=*), intent(in), dimension(:) :: lbuf
        character(len=*), intent(out), dimension(:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        character(len=len(lbuf(1))), allocatable, dimension(:) :: gbuf_pack
        character, allocatable, dimension(:) :: gbuf_pack1,lbuf1
	character(len=len(lbuf(1))) :: string

        integer :: recvtype, sendtype, ijk1,ijk2,sendcnt, ierr,lroot, lidebug
        integer :: i,j,k,ibuffer,iproc
        integer :: ijk, ijk_gl
        integer :: lenchar, icount
        logical :: isok_k,isok_j,isok_i, isinterior
        logical :: isbc_k,isbc_j,isbc_i, isboundary, need_copy
        include 'function.inc'

!       check to see whether there is root

        if (.not. present(mroot)) then
           lroot = 0
        else
           lroot = mroot
        endif

        if (.not. present(idebug)) then
           lidebug = 0
        else
           lidebug = idebug
        endif

	lenchar = len(lbuf(1))

        if(myPE.eq.lroot) then
	allocate(gbuf_pack(sum(ijksize3_all(:))))
        allocate(gbuf_pack1(sum(ijksize3_all(:))*lenchar))
	else
	allocate(gbuf_pack(10))
	allocate(gbuf_pack1(10))
	endif

        allocate(lbuf1(lenchar*size(lbuf)))


        recvtype = MPI_CHARACTER
        sendtype = recvtype

        ijk1 = ijkstart3
        ijk2 = ijkend3

        sendcnt = ijk2-ijk1+1

        icount = 0
        do i = 1,size(lbuf)
	    string = lbuf(i)(1:lenchar)
          do j = 1,lenchar

            icount = icount+1
            lbuf1(icount) = string(j:j)

          enddo
        enddo

        call MPI_Gatherv( lbuf1, sendcnt*lenchar, sendtype,  &
                           gbuf_pack1, ijksize3_all*lenchar, displs*lenchar, recvtype, &
                          lroot, MPI_COMM_WORLD, ierr )
        call MPI_Check( 'gather_1c:MPI_Gatherv', ierr )


	if(myPE.eq.lroot) then
        icount = 0
        do i = 1,size(gbuf_pack)
          do j = 1,lenchar

            icount = icount+1
            string(j:j) = gbuf_pack1(icount)

          enddo
	  gbuf_pack(i)(1:lenchar) = string(1:lenchar)

        enddo
	endif

        if( myPE.eq.lroot) then
        ibuffer = 0
          do iproc = 0,numPEs-1
            do k = kstart3_all(iproc), kend3_all(iproc)
              do j = jstart3_all(iproc), jend3_all(iproc)
                do i = istart3_all(iproc), iend3_all(iproc)

                  ibuffer = ibuffer + 1
                  isok_k = (kstart1_all(iproc) <= k) .and. (k <=kend1_all(iproc))
                  isok_j = (jstart1_all(iproc) <= j) .and. (j <=jend1_all(iproc))
                  isok_k = (istart1_all(iproc) <= i) .and. (i <=iend1_all(iproc))

                  isinterior = isok_k .and. isok_j .and. isok_i
                  isbc_k = (k <= kmin1) .or. (k >= kmax1)
                  isbc_j = (j <= jmin1) .or. (j >= jmax1)
                  isbc_i = (i <= imin1) .or. (i >= imax1)

                  isboundary = isbc_k .or. isbc_j .or. isbc_i
                  need_copy = isinterior .or. isboundary
                  if (need_copy) then
                      ijk_gl = funijk_gl(i,j,k)
		      string = gbuf_pack(ibuffer)(1:lenchar)
                      gbuf( ijk_gl )(1:lenchar) = string
                  endif

                enddo
              enddo
            enddo
          enddo
        endif

        deallocate(gbuf_pack)
        deallocate(gbuf_pack1)
        deallocate(lbuf1)


        return
        end subroutine gather_1c


        subroutine gather_1l( lbuf, gbuf, mroot, idebug )
        logical, intent(in), dimension(:) :: lbuf
        logical, intent(out), dimension(:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        logical, allocatable, dimension(:) :: gbuf_pack

        integer :: recvtype, sendtype, ijk1,ijk2,sendcnt, ierr,lroot, lidebug
        integer :: i,j,k,ibuffer,iproc
        integer :: ijk, ijk_gl
        logical :: isok_k,isok_j,isok_i, isinterior
        logical :: isbc_k,isbc_j,isbc_i, isboundary, need_copy
        include 'function.inc'

!       check to see whether there is root

        if (.not. present(mroot)) then
           lroot = 0
        else
           lroot = mroot
        endif

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

        if(myPE.eq.lroot) then
	allocate(gbuf_pack(sum(ijksize3_all(:))))
	else
	allocate(gbuf_pack(10))
	endif

        recvtype = MPI_LOGICAL
        sendtype = recvtype

        ijk1 = ijkstart3
        ijk2 = ijkend3

        sendcnt = ijk2-ijk1+1

        call MPI_Gatherv( lbuf, sendcnt, sendtype,  &
                           gbuf_pack, ijksize3_all, displs, recvtype, &
                          lroot, MPI_COMM_WORLD, ierr )
        call MPI_Check( 'gather_1l:MPI_Gatherv', ierr )

        if( myPE.eq.lroot) then
        ibuffer = 0
          do iproc = 0,numPEs-1
            do k = kstart3_all(iproc), kend3_all(iproc)
              do j = jstart3_all(iproc), jend3_all(iproc)
                do i = istart3_all(iproc), iend3_all(iproc)

                  ibuffer = ibuffer + 1
                  isok_k = (kstart1_all(iproc) <= k) .and. (k <=kend1_all(iproc))
                  isok_j = (jstart1_all(iproc) <= j) .and. (j <=jend1_all(iproc))
                  isok_k = (istart1_all(iproc) <= i) .and. (i <=iend1_all(iproc))

                  isinterior = isok_k .and. isok_j .and. isok_i

                  isbc_k = (k <= kmin1) .or. (k >= kmax1)
                  isbc_j = (j <= jmin1) .or. (j >= jmax1)
                  isbc_i = (i <= imin1) .or. (i >= imax1)

                  isboundary = isbc_k .or. isbc_j .or. isbc_i

                  need_copy = isinterior .or. isboundary
                  if (need_copy) then
                      ijk_gl = funijk_gl(i,j,k)
                      gbuf( ijk_gl ) = gbuf_pack(ibuffer)
                  endif

                enddo
              enddo
            enddo
          enddo
        endif

        deallocate(gbuf_pack)

        return
        end subroutine gather_1l


!	Routines to broadcast information from processor 0 in buffer to all
!	the processors	

        subroutine bcast_0i( buffer, mroot, idebug )
        integer, intent(inout) :: buffer
        integer, optional, intent(in) :: mroot, idebug

        integer :: datatype, count, ierr,lroot, lidebug

        if (.not. present(mroot)) then
           lroot = 0
        else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

        datatype = MPI_INTEGER

        count = 1

        call MPI_Bcast( buffer, count, datatype, lroot, MPI_COMM_WORLD, ierr)
        call MPI_Check( 'bcast_0i:MPI_Bcast', ierr )

        return
        end subroutine bcast_0i


	subroutine bcast_1i( buffer, mroot, idebug )
        integer, intent(inout), dimension(:) :: buffer       
	integer, optional, intent(in) :: mroot, idebug

	integer :: datatype, count, ierr,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

	datatype = MPI_INTEGER

	count = size(buffer,1)

        call MPI_Bcast( buffer, count, datatype, lroot, MPI_COMM_WORLD, ierr)
	call MPI_Check( 'bcast_1i:MPI_Bcast', ierr )

	return
	end subroutine bcast_1i

	
        subroutine bcast_2i( buffer, mroot, idebug )
        integer, intent(inout), dimension(:,:) :: buffer
        integer, optional, intent(in) :: mroot, idebug

	integer :: i,j,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

	do j=lbound(buffer,2),ubound(buffer,2)
	  call bcast_1i( buffer(:,j), lroot, lidebug )
	enddo

	return
	end subroutine bcast_2i

        subroutine bcast_3i( buffer, mroot, idebug )
        integer, intent(inout), dimension(:,:,:) :: buffer
        integer, optional, intent(in) :: mroot, idebug

        integer :: j,k,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

        do k=lbound(buffer,3),ubound(buffer,3)
        do j=lbound(buffer,2),ubound(buffer,2)
          call bcast_1i( buffer(:,j,k), lroot, lidebug )
        enddo
	enddo

        return
        end subroutine bcast_3i

        subroutine bcast_0r( buffer, mroot, idebug )
        real, intent(inout) :: buffer
        integer, optional, intent(in) :: mroot, idebug

        integer :: datatype, count, ierr,lroot, lidebug

        if (.not. present(mroot)) then
           lroot = 0
        else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

        datatype = MPI_REAL

        count = 1

        call MPI_Bcast( buffer, count, datatype, lroot, MPI_COMM_WORLD, ierr)
        call MPI_Check( 'bcast_0r:MPI_Bcast', ierr )

        return
        end subroutine bcast_0r


	subroutine bcast_1r( buffer, mroot, idebug )
        real, intent(inout), dimension(:) :: buffer
	integer, optional, intent(in) :: mroot, idebug

	integer :: datatype, count, ierr,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

	datatype = MPI_REAL

        count = size(buffer,1)

        call MPI_Bcast( buffer, count, datatype, lroot, MPI_COMM_WORLD, ierr)
	call MPI_Check( 'bcast_1r:MPI_Bcast', ierr )

	return
	end subroutine bcast_1r

	
        subroutine bcast_2r( buffer, mroot, idebug )
        real, intent(inout), dimension(:,:) :: buffer
        integer, optional, intent(in) :: mroot, idebug

	integer :: i,j,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

	do j=lbound(buffer,2),ubound(buffer,2)
	  call bcast_1r( buffer(:,j), lroot, lidebug )
	enddo

	return
	end subroutine bcast_2r

        subroutine bcast_3r( buffer, mroot, idebug )
        real, intent(inout), dimension(:,:,:) :: buffer
        integer, optional, intent(in) :: mroot, idebug

        integer :: j,k,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

        do k=lbound(buffer,3),ubound(buffer,3)
        do j=lbound(buffer,2),ubound(buffer,2)
          call bcast_1r( buffer(:,j,k), lroot, lidebug )
        enddo
	enddo

        return
        end subroutine bcast_3r

        subroutine bcast_0d( buffer, mroot, idebug )
        double precision, intent(inout) :: buffer
        integer, optional, intent(in) :: mroot, idebug

        integer :: datatype, count, ierr,lroot, lidebug

        if (.not. present(mroot)) then
           lroot = 0
        else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

        datatype = MPI_DOUBLE_PRECISION

        count = 1

        call MPI_Bcast( buffer, count, datatype, lroot, MPI_COMM_WORLD, ierr)
        call MPI_Check( 'bcast_0d:MPI_Bcast', ierr )

        return
        end subroutine bcast_0d


	subroutine bcast_1d( buffer, mroot, idebug )
        double precision, intent(inout), dimension(:) :: buffer
	integer, optional, intent(in) :: mroot, idebug

	integer :: datatype, count, ierr,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

	datatype = MPI_DOUBLE_PRECISION

        count = size(buffer,1)

        call MPI_Bcast( buffer, count, datatype, lroot, MPI_COMM_WORLD, ierr)
	call MPI_Check( 'bcast_1d:MPI_Bcast', ierr )

	return
	end subroutine bcast_1d

	
        subroutine bcast_2d( buffer, mroot, idebug )
        double precision, intent(inout), dimension(:,:) :: buffer
        integer, optional, intent(in) :: mroot, idebug

	integer :: i,j,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

	do j=lbound(buffer,2),ubound(buffer,2)
	  call bcast_1d( buffer(:,j), lroot, lidebug )
	enddo

	return
	end subroutine bcast_2d

        subroutine bcast_3d( buffer, mroot, idebug )
        double precision, intent(inout), dimension(:,:,:) :: buffer
        integer, optional, intent(in) :: mroot, idebug

        integer :: j,k,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

        do k=lbound(buffer,3),ubound(buffer,3)
        do j=lbound(buffer,2),ubound(buffer,2)
          call bcast_1d( buffer(:,j,k), lroot, lidebug )
        enddo
	enddo

        return
        end subroutine bcast_3d

        subroutine bcast_0c( buffer, mroot, idebug )
        character(len=*), intent(inout) :: buffer
        integer, optional, intent(in) :: mroot, idebug
	character, allocatable, dimension(:) :: buffer1

        integer :: datatype, count, ierr,lroot, lidebug
        integer :: lenchar,icount, i, j

        if (.not. present(mroot)) then
           lroot = 0
        else
           lroot = mroot
        endif

        if (.not. present(idebug)) then
           lidebug = 0
        else
           lidebug = idebug
        endif

	lenchar = len(buffer)

	allocate(buffer1(lenchar))

        icount = 0
        do j = 1,lenchar

          icount = icount+1
          buffer1(icount) = buffer(j:j)

        enddo

        datatype = MPI_CHARACTER

        count = 1

        call MPI_Bcast( buffer1, count*lenchar, datatype, lroot, MPI_COMM_WORLD, ierr)
        call MPI_Check( 'bcast_0c:MPI_Bcast', ierr )

        icount = 0
        do j = 1,lenchar

          icount = icount+1
          buffer(j:j) = buffer1(icount)

        enddo

        return
        end subroutine bcast_0c


        subroutine bcast_1c( buffer, mroot, idebug )
        character(len=*), intent(inout), dimension(:) :: buffer
        integer, optional, intent(in) :: mroot, idebug
        character, allocatable, dimension(:) :: buffer1

        integer :: datatype, count, ierr,lroot, lidebug
        integer :: lenchar,icount, i, j
	character(len=len(buffer(1))) :: string

        if (.not. present(mroot)) then
           lroot = 0
        else
           lroot = mroot
        endif 

        if (.not. present(idebug)) then
           lidebug = 0
        else
           lidebug = idebug
        endif

        lenchar = len(buffer(1))

        allocate(buffer1(size(buffer)*lenchar))

        icount = 0
        do i = 1,size(buffer)
          string = buffer(i)(1:lenchar)
          do j = 1,lenchar

            icount = icount+1
            buffer1(icount) = string(j:j)

          enddo
        enddo

        datatype = MPI_CHARACTER

        count = size(buffer,1)

        call MPI_Bcast( buffer1, count*lenchar, datatype, lroot, MPI_COMM_WORLD, ierr)
        call MPI_Check( 'bcast_1c:MPI_Bcast', ierr )

        icount = 0
        do i = 1,size(buffer)
          do j = 1,lenchar

            icount = icount+1
            string(j:j) = buffer1(icount)

          enddo
	    buffer(i) = string
        enddo

        return
        end subroutine bcast_1c

        subroutine bcast_0l( buffer, mroot, idebug )
        logical, intent(inout) :: buffer
        integer, optional, intent(in) :: mroot, idebug

        integer :: datatype, count, ierr,lroot, lidebug

        if (.not. present(mroot)) then
           lroot = 0
        else
           lroot = mroot
        endif

        if (.not. present(idebug)) then
           lidebug = 0
        else
           lidebug = idebug
        endif

        datatype = MPI_LOGICAL

        count = 1

        call MPI_Bcast( buffer, count, datatype, lroot, MPI_COMM_WORLD, ierr)
        call MPI_Check( 'bcast_0l:MPI_Bcast', ierr )

        return
        end subroutine bcast_0l


        subroutine bcast_1l( buffer, mroot, idebug )
        logical, intent(inout), dimension(:) :: buffer
        integer, optional, intent(in) :: mroot, idebug

        integer :: datatype, count, ierr,lroot, lidebug

        if (.not. present(mroot)) then
           lroot = 0
        else
           lroot = mroot
        endif 

        if (.not. present(idebug)) then
           lidebug = 0
        else
           lidebug = idebug
        endif

        datatype = MPI_LOGICAL

        count = size(buffer,1)

        call MPI_Bcast( buffer, count, datatype, lroot, MPI_COMM_WORLD, ierr)
        call MPI_Check( 'bcast_1l:MPI_Bcast', ierr )

        return
        end subroutine bcast_1l


!	Procedures to do global operations (Sum, Min, Max). _all_ routines
!	send the information to all the processors otherwise they are
!	kept on processor 0.

        subroutine global_sum_0i( lbuf, gbuf, mroot, idebug )
        integer, intent(in) :: lbuf
        integer, intent(out) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

        if (.not. present(mroot)) then
           lroot = 0
        else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

        recvtype = MPI_INTEGER
        sendtype = recvtype

        sendcnt = 1

        call MPI_Reduce( lbuf, gbuf, sendcnt, sendtype,  MPI_SUM, &
                          lroot, MPI_COMM_WORLD, ierr )
        call MPI_Check( 'global_sum_0i:MPI_Reduce', ierr )

        return
        end subroutine global_sum_0i


	subroutine global_sum_1i( lbuf, gbuf, mroot, idebug )
        integer, intent(in), dimension(:) :: lbuf       
        integer, intent(out), dimension(:) :: gbuf
	integer, optional, intent(in) :: mroot, idebug

	integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

	recvtype = MPI_INTEGER
	sendtype = recvtype

	sendcnt = size(lbuf)

        call MPI_Reduce( lbuf, gbuf, sendcnt, sendtype,  MPI_SUM, &
                          lroot, MPI_COMM_WORLD, ierr )
	call MPI_Check( 'global_sum_1i:MPI_Reduce', ierr )

	return
	end subroutine global_sum_1i

        subroutine global_sum_2i( lbuf, gbuf, mroot, idebug )
        integer, intent(in), dimension(:,:) :: lbuf
        integer, intent(out), dimension(:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

	integer :: i,j,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

	if(myPE.eq.lroot) then
	call assert( size(lbuf,2).eq.size(gbuf,2),  &
		'** global_sum_2i: size(lbuf,2).ne.size(gbuf,2) ', &
		size(lbuf,2), size(gbuf,2) )
	endif

	do j=lbound(lbuf,2),ubound(lbuf,2)
	  call global_sum_1i( lbuf(:,j), gbuf(:,j), lroot, lidebug )
	enddo

	return
	end subroutine global_sum_2i

        subroutine global_sum_3i( lbuf, gbuf, mroot, idebug )
        integer, intent(in), dimension(:,:,:) :: lbuf
        integer, intent(out), dimension(:,:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        integer :: j,k,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

	if(myPE.eq.lroot) then
        call assert( size(lbuf,2).eq.size(gbuf,2),  &
                '** global_sum_3i: size(lbuf,2).ne.size(gbuf,2) ', &
                size(lbuf,2), size(gbuf,2) )

        call assert( size(lbuf,3).eq.size(gbuf,3),  &
                '** global_sum_3i: size(lbuf,3).ne.size(gbuf,3) ', &
                size(lbuf,3), size(gbuf,3) )
	endif

        do k=lbound(lbuf,3),ubound(lbuf,3)
        do j=lbound(lbuf,2),ubound(lbuf,2)
          call global_sum_1i( lbuf(:,j,k), gbuf(:,j,k), lroot, lidebug )
        enddo
	enddo

        return
        end subroutine global_sum_3i

        subroutine global_sum_0r( lbuf, gbuf, mroot, idebug )
        real, intent(in) :: lbuf
        real, intent(out) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

        if (.not. present(mroot)) then
           lroot = 0
        else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

        recvtype = MPI_REAL
        sendtype = recvtype

        sendcnt = 1

        call MPI_Reduce( lbuf, gbuf, sendcnt, sendtype,  MPI_SUM, &
                          lroot, MPI_COMM_WORLD, ierr )
        call MPI_Check( 'global_sum_0r:MPI_Reduce', ierr )

        return
        end subroutine global_sum_0r


	subroutine global_sum_1r( lbuf, gbuf, mroot, idebug )
        real, intent(in), dimension(:) :: lbuf       
        real, intent(out), dimension(:) :: gbuf
	integer, optional, intent(in) :: mroot, idebug

	integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

	recvtype = MPI_REAL
	sendtype = recvtype

	sendcnt = size(lbuf)

        call MPI_Reduce( lbuf, gbuf, sendcnt, sendtype,  MPI_SUM, &
                          lroot, MPI_COMM_WORLD, ierr )
	call MPI_Check( 'global_sum_1r:MPI_Reduce', ierr )

	return
	end subroutine global_sum_1r

        subroutine global_sum_2r( lbuf, gbuf, mroot, idebug )
        real, intent(in), dimension(:,:) :: lbuf
        real, intent(out), dimension(:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

	integer :: i,j,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

	if(myPE.eq.lroot) then
	call assert( size(lbuf,2).eq.size(gbuf,2),  &
		'** global_sum_2r: size(lbuf,2).ne.size(gbuf,2) ', &
		size(lbuf,2), size(gbuf,2) )
	endif

	do j=lbound(lbuf,2),ubound(lbuf,2)
	  call global_sum_1r( lbuf(:,j), gbuf(:,j), lroot, lidebug )
	enddo

	return
	end subroutine global_sum_2r

        subroutine global_sum_3r( lbuf, gbuf, mroot, idebug )
        real, intent(in), dimension(:,:,:) :: lbuf
        real, intent(out), dimension(:,:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        integer :: j,k,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

	if(myPE.eq.lroot) then
        call assert( size(lbuf,2).eq.size(gbuf,2),  &
                '** global_sum_3i: size(lbuf,2).ne.size(gbuf,2) ', &
                size(lbuf,2), size(gbuf,2) )

        call assert( size(lbuf,3).eq.size(gbuf,3),  &
                '** global_sum_3i: size(lbuf,3).ne.size(gbuf,3) ', &
                size(lbuf,3), size(gbuf,3) )
	endif

        do k=lbound(lbuf,3),ubound(lbuf,3)
        do j=lbound(lbuf,2),ubound(lbuf,2)
          call global_sum_1r( lbuf(:,j,k), gbuf(:,j,k), lroot, lidebug )
        enddo
	enddo

        return
        end subroutine global_sum_3r

        subroutine global_sum_0d( lbuf, gbuf, mroot, idebug )
        double precision, intent(in) :: lbuf
        double precision, intent(out) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

        if (.not. present(mroot)) then
           lroot = 0
        else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

        recvtype = MPI_DOUBLE_PRECISION
        sendtype = recvtype

        sendcnt = 1

        call MPI_Reduce( lbuf, gbuf, sendcnt, sendtype,  MPI_SUM, &
                          lroot, MPI_COMM_WORLD, ierr )
        call MPI_Check( 'global_sum_0d:MPI_Reduce', ierr )

        return
        end subroutine global_sum_0d


	subroutine global_sum_1d( lbuf, gbuf, mroot, idebug )
        double precision, intent(in), dimension(:) :: lbuf       
        double precision, intent(out), dimension(:) :: gbuf
	integer, optional, intent(in) :: mroot, idebug

	integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

	recvtype = MPI_DOUBLE_PRECISION
	sendtype = recvtype

	sendcnt = size(lbuf)

        call MPI_Reduce( lbuf, gbuf, sendcnt, sendtype,  MPI_SUM, &
                          lroot, MPI_COMM_WORLD, ierr )
	call MPI_Check( 'global_sum_1d:MPI_Reduce', ierr )

	return
	end subroutine global_sum_1d

        subroutine global_sum_2d( lbuf, gbuf, mroot, idebug )
        double precision, intent(in), dimension(:,:) :: lbuf
        double precision, intent(out), dimension(:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

	integer :: i,j,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

	if(myPE.eq.lroot) then
	call assert( size(lbuf,2).eq.size(gbuf,2),  &
		'** global_sum_2d: size(lbuf,2).ne.size(gbuf,2) ', &
		size(lbuf,2), size(gbuf,2) )
	endif

	do j=lbound(lbuf,2),ubound(lbuf,2)
	  call global_sum_1d( lbuf(:,j), gbuf(:,j), lroot, lidebug )
	enddo

	return
	end subroutine global_sum_2d

        subroutine global_sum_3d( lbuf, gbuf, mroot, idebug )
        double precision, intent(in), dimension(:,:,:) :: lbuf
        double precision, intent(out), dimension(:,:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        integer :: j,k,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

	if(myPE.eq.lroot) then
        call assert( size(lbuf,2).eq.size(gbuf,2),  &
                '** global_sum_3i: size(lbuf,2).ne.size(gbuf,2) ', &
                size(lbuf,2), size(gbuf,2) )

        call assert( size(lbuf,3).eq.size(gbuf,3),  &
                '** global_sum_3i: size(lbuf,3).ne.size(gbuf,3) ', &
                size(lbuf,3), size(gbuf,3) )
	endif

        do k=lbound(lbuf,3),ubound(lbuf,3)
        do j=lbound(lbuf,2),ubound(lbuf,2)
          call global_sum_1d( lbuf(:,j,k), gbuf(:,j,k), lroot, lidebug )
        enddo
	enddo

        return
        end subroutine global_sum_3d

        subroutine global_all_sum_0i( lbuf, gbuf, mroot, idebug )
        integer, intent(in) :: lbuf
        integer, intent(out) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

        if (.not. present(mroot)) then
           lroot = 0
        else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

        recvtype = MPI_INTEGER
        sendtype = recvtype

        sendcnt = 1

        call MPI_Allreduce( lbuf, gbuf, sendcnt, sendtype, MPI_SUM, &
                          MPI_COMM_WORLD, ierr )
        call MPI_Check( 'global_all_sum_0i:MPI_Allreduce', ierr )

        return
        end subroutine global_all_sum_0i


	subroutine global_all_sum_1i( lbuf, gbuf, mroot, idebug )
        integer, intent(in), dimension(:) :: lbuf       
        integer, intent(out), dimension(:) :: gbuf
	integer, optional, intent(in) :: mroot, idebug

	integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

	recvtype = MPI_INTEGER
	sendtype = recvtype

	sendcnt = size(lbuf)

        call MPI_Allreduce( lbuf, gbuf, sendcnt, sendtype, MPI_SUM, &
                          MPI_COMM_WORLD, ierr )
	call MPI_Check( 'global_all_sum_1i:MPI_Allreduce', ierr )

	return
	end subroutine global_all_sum_1i

        subroutine global_all_sum_2i( lbuf, gbuf, mroot, idebug )
        integer, intent(in), dimension(:,:) :: lbuf
        integer, intent(out), dimension(:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

	integer :: i,j,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

	call assert( size(lbuf,2).eq.size(gbuf,2),  &
		'** global_all_sum_2i: size(lbuf,2).ne.size(gbuf,2) ', &
		size(lbuf,2), size(gbuf,2) )

	do j=lbound(lbuf,2),ubound(lbuf,2)
	  call global_all_sum_1i( lbuf(:,j), gbuf(:,j), lroot, lidebug )
	enddo

	return
	end subroutine global_all_sum_2i

        subroutine global_all_sum_3i( lbuf, gbuf, mroot, idebug )
        integer, intent(in), dimension(:,:,:) :: lbuf
        integer, intent(out), dimension(:,:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        integer :: j,k,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

        call assert( size(lbuf,2).eq.size(gbuf,2),  &
                '** global_all_sum_3i: size(lbuf,2).ne.size(gbuf,2) ', &
                size(lbuf,2), size(gbuf,2) )

        call assert( size(lbuf,3).eq.size(gbuf,3),  &
                '** global_all_sum_3i: size(lbuf,3).ne.size(gbuf,3) ', &
                size(lbuf,3), size(gbuf,3) )

        do k=lbound(lbuf,3),ubound(lbuf,3)
        do j=lbound(lbuf,2),ubound(lbuf,2)
          call global_all_sum_1i( lbuf(:,j,k), gbuf(:,j,k), lroot, lidebug )
        enddo
	enddo

        return
        end subroutine global_all_sum_3i

        subroutine global_all_sum_0r( lbuf, gbuf, mroot, idebug )
        real, intent(in) :: lbuf
        real, intent(out) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

        if (.not. present(mroot)) then
           lroot = 0
        else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

        recvtype = MPI_REAL
        sendtype = recvtype

        sendcnt = 1

        call MPI_Allreduce( lbuf, gbuf, sendcnt, sendtype, MPI_SUM, &
                          MPI_COMM_WORLD, ierr )
        call MPI_Check( 'global_all_sum_0r:MPI_Allreduce', ierr )


        return
        end subroutine global_all_sum_0r


	subroutine global_all_sum_1r( lbuf, gbuf, mroot, idebug )
        real, intent(in), dimension(:) :: lbuf       
        real, intent(out), dimension(:) :: gbuf
	integer, optional, intent(in) :: mroot, idebug

	integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

	recvtype = MPI_REAL
	sendtype = recvtype

	sendcnt = size(lbuf)

        call MPI_Allreduce( lbuf, gbuf, sendcnt, sendtype, MPI_SUM, &
                            MPI_COMM_WORLD, ierr )
	call MPI_Check( 'global_all_sum_1r:MPI_Allreduce', ierr )

	return
	end subroutine global_all_sum_1r

        subroutine global_all_sum_2r( lbuf, gbuf, mroot, idebug )
        real, intent(in), dimension(:,:) :: lbuf
        real, intent(out), dimension(:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

	integer :: i,j,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

	call assert( size(lbuf,2).eq.size(gbuf,2),  &
		'** global_all_sum_2r: size(lbuf,2).ne.size(gbuf,2) ', &
		size(lbuf,2), size(gbuf,2) )

	do j=lbound(lbuf,2),ubound(lbuf,2)
	  call global_all_sum_1r( lbuf(:,j), gbuf(:,j), lroot, lidebug )
	enddo

	return
	end subroutine global_all_sum_2r

        subroutine global_all_sum_3r( lbuf, gbuf, mroot, idebug )
        real, intent(in), dimension(:,:,:) :: lbuf
        real, intent(out), dimension(:,:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        integer :: j,k,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

        call assert( size(lbuf,2).eq.size(gbuf,2),  &
                '** global_all_sum_3i: size(lbuf,2).ne.size(gbuf,2) ', &
                size(lbuf,2), size(gbuf,2) )

        call assert( size(lbuf,3).eq.size(gbuf,3),  &
                '** global_all_sum_3i: size(lbuf,3).ne.size(gbuf,3) ', &
                size(lbuf,3), size(gbuf,3) )

        do k=lbound(lbuf,3),ubound(lbuf,3)
        do j=lbound(lbuf,2),ubound(lbuf,2)
          call global_all_sum_1r( lbuf(:,j,k), gbuf(:,j,k), lroot, lidebug )
        enddo
	enddo

        return
        end subroutine global_all_sum_3r

        subroutine global_all_sum_0d( lbuf, gbuf, mroot, idebug )
        double precision, intent(in) :: lbuf
        double precision, intent(out) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

        if (.not. present(mroot)) then
           lroot = 0
        else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

        recvtype = MPI_DOUBLE_PRECISION
        sendtype = recvtype

        sendcnt = 1

        call MPI_Allreduce( lbuf, gbuf, sendcnt, sendtype, MPI_SUM, &
                          MPI_COMM_WORLD, ierr )
        call MPI_Check( 'global_all_sum_0d:MPI_Allreduce', ierr )

        return
        end subroutine global_all_sum_0d


	subroutine global_all_sum_1d( lbuf, gbuf, mroot, idebug )
        double precision, intent(in), dimension(:) :: lbuf       
        double precision, intent(out), dimension(:) :: gbuf
	integer, optional, intent(in) :: mroot, idebug

	integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

	recvtype = MPI_DOUBLE_PRECISION
	sendtype = recvtype

	sendcnt = size(lbuf)

        call MPI_Allreduce( lbuf, gbuf, sendcnt, sendtype,  MPI_SUM, &
                          MPI_COMM_WORLD, ierr )
	call MPI_Check( 'global_all_sum_1d:MPI_Allreduce', ierr )

	return
	end subroutine global_all_sum_1d

        subroutine global_all_sum_2d( lbuf, gbuf, mroot, idebug )
        double precision, intent(in), dimension(:,:) :: lbuf
        double precision, intent(out), dimension(:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

	integer :: i,j,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

	call assert( size(lbuf,2).eq.size(gbuf,2),  &
		'** global_all_sum_2d: size(lbuf,2).ne.size(gbuf,2) ', &
		size(lbuf,2), size(gbuf,2) )

	do j=lbound(lbuf,2),ubound(lbuf,2)
	  call global_all_sum_1d( lbuf(:,j), gbuf(:,j), lroot, lidebug )
	enddo

	return
	end subroutine global_all_sum_2d

        subroutine global_all_sum_3d( lbuf, gbuf, mroot, idebug )
        double precision, intent(in), dimension(:,:,:) :: lbuf
        double precision, intent(out), dimension(:,:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        integer :: j,k,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

        call assert( size(lbuf,2).eq.size(gbuf,2),  &
                '** global_all_sum_3i: size(lbuf,2).ne.size(gbuf,2) ', &
                size(lbuf,2), size(gbuf,2) )

        call assert( size(lbuf,3).eq.size(gbuf,3),  &
                '** global_all_sum_3i: size(lbuf,3).ne.size(gbuf,3) ', &
                size(lbuf,3), size(gbuf,3) )

        do k=lbound(lbuf,3),ubound(lbuf,3)
        do j=lbound(lbuf,2),ubound(lbuf,2)
          call global_all_sum_1d( lbuf(:,j,k), gbuf(:,j,k), lroot, lidebug )
        enddo
	enddo

        return
        end subroutine global_all_sum_3d

        subroutine global_min_0i( lbuf, gbuf, mroot, idebug )
        integer, intent(in) :: lbuf
        integer, intent(out) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

        if (.not. present(mroot)) then
           lroot = 0
        else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

        recvtype = MPI_INTEGER
        sendtype = recvtype

        sendcnt = 1

        call MPI_Reduce( lbuf, gbuf, sendcnt, sendtype,  MPI_MIN, &
                          lroot, MPI_COMM_WORLD, ierr )
        call MPI_Check( 'global_min_0i:MPI_Reduce', ierr )

        return
        end subroutine global_min_0i


	subroutine global_min_1i( lbuf, gbuf, mroot, idebug )
        integer, intent(in), dimension(:) :: lbuf       
        integer, intent(out), dimension(:) :: gbuf
	integer, optional, intent(in) :: mroot, idebug

	integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

	recvtype = MPI_INTEGER
	sendtype = recvtype

	sendcnt = size(lbuf)

        call MPI_Reduce( lbuf, gbuf, sendcnt, sendtype,  MPI_MIN, &
                          lroot, MPI_COMM_WORLD, ierr )
	call MPI_Check( 'global_min_1i:MPI_Reduce', ierr )

	return
	end subroutine global_min_1i

        subroutine global_min_2i( lbuf, gbuf, mroot, idebug )
        integer, intent(in), dimension(:,:) :: lbuf
        integer, intent(out), dimension(:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

	integer :: i,j,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

	if(myPE.eq.lroot) then
	call assert( size(lbuf,2).eq.size(gbuf,2),  &
		'** global_min_2i: size(lbuf,2).ne.size(gbuf,2) ', &
		size(lbuf,2), size(gbuf,2) )
	endif

	do j=lbound(lbuf,2),ubound(lbuf,2)
	  call global_min_1i( lbuf(:,j), gbuf(:,j), lroot, lidebug )
	enddo

	return
	end subroutine global_min_2i

        subroutine global_min_3i( lbuf, gbuf, mroot, idebug )
        integer, intent(in), dimension(:,:,:) :: lbuf
        integer, intent(out), dimension(:,:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        integer :: j,k,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

	if(myPE.eq.lroot) then
        call assert( size(lbuf,2).eq.size(gbuf,2),  &
                '** global_min_3i: size(lbuf,2).ne.size(gbuf,2) ', &
                size(lbuf,2), size(gbuf,2) )

        call assert( size(lbuf,3).eq.size(gbuf,3),  &
                '** global_min_3i: size(lbuf,3).ne.size(gbuf,3) ', &
                size(lbuf,3), size(gbuf,3) )
	endif

        do k=lbound(lbuf,3),ubound(lbuf,3)
        do j=lbound(lbuf,2),ubound(lbuf,2)
          call global_min_1i( lbuf(:,j,k), gbuf(:,j,k), lroot, lidebug )
        enddo
	enddo

        return
        end subroutine global_min_3i

        subroutine global_min_0r( lbuf, gbuf, mroot, idebug )
        real, intent(in) :: lbuf
        real, intent(out) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

        if (.not. present(mroot)) then
           lroot = 0
        else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

        recvtype = MPI_REAL
        sendtype = recvtype

        sendcnt = 1

        call MPI_Reduce( lbuf, gbuf, sendcnt, sendtype,  MPI_MIN, &
                          lroot, MPI_COMM_WORLD, ierr )
        call MPI_Check( 'global_min_0r:MPI_Reduce', ierr )

        return
        end subroutine global_min_0r


	subroutine global_min_1r( lbuf, gbuf, mroot, idebug )
        real, intent(in), dimension(:) :: lbuf       
        real, intent(out), dimension(:) :: gbuf
	integer, optional, intent(in) :: mroot, idebug

	integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

	recvtype = MPI_REAL
	sendtype = recvtype

	sendcnt = size(lbuf)

        call MPI_Reduce( lbuf, gbuf, sendcnt, sendtype,  MPI_MIN, &
                          lroot, MPI_COMM_WORLD, ierr )
	call MPI_Check( 'global_min_1r:MPI_Reduce', ierr )

	return
	end subroutine global_min_1r

        subroutine global_min_2r( lbuf, gbuf, mroot, idebug )
        real, intent(in), dimension(:,:) :: lbuf
        real, intent(out), dimension(:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

	integer :: i,j,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

	if(myPE.eq.lroot) then
	call assert( size(lbuf,2).eq.size(gbuf,2),  &
		'** global_min_2r: size(lbuf,2).ne.size(gbuf,2) ', &
		size(lbuf,2), size(gbuf,2) )
	endif

	do j=lbound(lbuf,2),ubound(lbuf,2)
	  call global_min_1r( lbuf(:,j), gbuf(:,j), lroot, lidebug )
	enddo

	return
	end subroutine global_min_2r

        subroutine global_min_3r( lbuf, gbuf, mroot, idebug )
        real, intent(in), dimension(:,:,:) :: lbuf
        real, intent(out), dimension(:,:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        integer :: j,k,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

	if(myPE.eq.lroot) then
        call assert( size(lbuf,2).eq.size(gbuf,2),  &
                '** global_min_3i: size(lbuf,2).ne.size(gbuf,2) ', &
                size(lbuf,2), size(gbuf,2) )

        call assert( size(lbuf,3).eq.size(gbuf,3),  &
                '** global_min_3i: size(lbuf,3).ne.size(gbuf,3) ', &
                size(lbuf,3), size(gbuf,3) )
	endif

        do k=lbound(lbuf,3),ubound(lbuf,3)
        do j=lbound(lbuf,2),ubound(lbuf,2)
          call global_min_1r( lbuf(:,j,k), gbuf(:,j,k), lroot, lidebug )
        enddo
	enddo

        return
        end subroutine global_min_3r

        subroutine global_min_0d( lbuf, gbuf, mroot, idebug )
        double precision, intent(in) :: lbuf
        double precision, intent(out) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

        if (.not. present(mroot)) then
           lroot = 0
        else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

        recvtype = MPI_DOUBLE_PRECISION
        sendtype = recvtype

        sendcnt = 1

        call MPI_Reduce( lbuf, gbuf, sendcnt, sendtype,  MPI_MIN, &
                          lroot, MPI_COMM_WORLD, ierr )
        call MPI_Check( 'global_min_0d:MPI_Reduce', ierr )

        return
        end subroutine global_min_0d


	subroutine global_min_1d( lbuf, gbuf, mroot, idebug )
        double precision, intent(in), dimension(:) :: lbuf       
        double precision, intent(out), dimension(:) :: gbuf
	integer, optional, intent(in) :: mroot, idebug

	integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

	recvtype = MPI_DOUBLE_PRECISION
	sendtype = recvtype

	sendcnt = size(lbuf)

        call MPI_Reduce( lbuf, gbuf, sendcnt, sendtype,  MPI_MIN, &
                          lroot, MPI_COMM_WORLD, ierr )
	call MPI_Check( 'global_min_1d:MPI_Reduce', ierr )

	return
	end subroutine global_min_1d

        subroutine global_min_2d( lbuf, gbuf, mroot, idebug )
        double precision, intent(in), dimension(:,:) :: lbuf
        double precision, intent(out), dimension(:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

	integer :: i,j,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

	if(myPE.eq.lroot) then
	call assert( size(lbuf,2).eq.size(gbuf,2),  &
		'** global_min_2d: size(lbuf,2).ne.size(gbuf,2) ', &
		size(lbuf,2), size(gbuf,2) )
	endif

	do j=lbound(lbuf,2),ubound(lbuf,2)
	  call global_min_1d( lbuf(:,j), gbuf(:,j), lroot, lidebug )
	enddo

	return
	end subroutine global_min_2d

        subroutine global_min_3d( lbuf, gbuf, mroot, idebug )
        double precision, intent(in), dimension(:,:,:) :: lbuf
        double precision, intent(out), dimension(:,:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        integer :: j,k,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

	if(myPE.eq.lroot) then
        call assert( size(lbuf,2).eq.size(gbuf,2),  &
                '** global_min_3i: size(lbuf,2).ne.size(gbuf,2) ', &
                size(lbuf,2), size(gbuf,2) )

        call assert( size(lbuf,3).eq.size(gbuf,3),  &
                '** global_min_3i: size(lbuf,3).ne.size(gbuf,3) ', &
                size(lbuf,3), size(gbuf,3) )
	endif

        do k=lbound(lbuf,3),ubound(lbuf,3)
        do j=lbound(lbuf,2),ubound(lbuf,2)
          call global_min_1d( lbuf(:,j,k), gbuf(:,j,k), lroot, lidebug )
        enddo
	enddo

        return
        end subroutine global_min_3d

        subroutine global_all_min_0i( lbuf, gbuf, mroot, idebug )
        integer, intent(in) :: lbuf
        integer, intent(out) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

        if (.not. present(mroot)) then
           lroot = 0
        else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

        recvtype = MPI_INTEGER
        sendtype = recvtype

        sendcnt = 1

        call MPI_Allreduce( lbuf, gbuf, sendcnt, sendtype, MPI_MIN, &
                          MPI_COMM_WORLD, ierr )
        call MPI_Check( 'global_all_min_0i:MPI_Allreduce', ierr )

        return
        end subroutine global_all_min_0i


	subroutine global_all_min_1i( lbuf, gbuf, mroot, idebug )
        integer, intent(in), dimension(:) :: lbuf       
        integer, intent(out), dimension(:) :: gbuf
	integer, optional, intent(in) :: mroot, idebug

	integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

	recvtype = MPI_INTEGER
	sendtype = recvtype

	sendcnt = size(lbuf)

        call MPI_Allreduce( lbuf, gbuf, sendcnt, sendtype, MPI_MIN, &
                          MPI_COMM_WORLD, ierr )
	call MPI_Check( 'global_all_min_1i:MPI_Allreduce', ierr )

	return
	end subroutine global_all_min_1i

        subroutine global_all_min_2i( lbuf, gbuf, mroot, idebug )
        integer, intent(in), dimension(:,:) :: lbuf
        integer, intent(out), dimension(:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

	integer :: i,j,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

	call assert( size(lbuf,2).eq.size(gbuf,2),  &
		'** global_all_min_2i: size(lbuf,2).ne.size(gbuf,2) ', &
		size(lbuf,2), size(gbuf,2) )

	do j=lbound(lbuf,2),ubound(lbuf,2)
	  call global_all_min_1i( lbuf(:,j), gbuf(:,j), lroot, lidebug )
	enddo

	return
	end subroutine global_all_min_2i

        subroutine global_all_min_3i( lbuf, gbuf, mroot, idebug )
        integer, intent(in), dimension(:,:,:) :: lbuf
        integer, intent(out), dimension(:,:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        integer :: j,k,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

        call assert( size(lbuf,2).eq.size(gbuf,2),  &
                '** global_all_min_3i: size(lbuf,2).ne.size(gbuf,2) ', &
                size(lbuf,2), size(gbuf,2) )

        call assert( size(lbuf,3).eq.size(gbuf,3),  &
                '** global_all_min_3i: size(lbuf,3).ne.size(gbuf,3) ', &
                size(lbuf,3), size(gbuf,3) )

        do k=lbound(lbuf,3),ubound(lbuf,3)
        do j=lbound(lbuf,2),ubound(lbuf,2)
          call global_all_min_1i( lbuf(:,j,k), gbuf(:,j,k), lroot, lidebug )
        enddo
	enddo

        return
        end subroutine global_all_min_3i

        subroutine global_all_min_0r( lbuf, gbuf, mroot, idebug )
        real, intent(in) :: lbuf
        real, intent(out) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

        if (.not. present(mroot)) then
           lroot = 0
        else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

        recvtype = MPI_REAL
        sendtype = recvtype

        sendcnt = 1

        call MPI_Allreduce( lbuf, gbuf, sendcnt, sendtype, MPI_MIN, &
                          MPI_COMM_WORLD, ierr )
        call MPI_Check( 'global_all_min_0r:MPI_Allreduce', ierr )


        return
        end subroutine global_all_min_0r


	subroutine global_all_min_1r( lbuf, gbuf, mroot, idebug )
        real, intent(in), dimension(:) :: lbuf       
        real, intent(out), dimension(:) :: gbuf
	integer, optional, intent(in) :: mroot, idebug

	integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

	recvtype = MPI_REAL
	sendtype = recvtype

	sendcnt = size(lbuf)

        call MPI_Allreduce( lbuf, gbuf, sendcnt, sendtype, MPI_MIN, &
                            MPI_COMM_WORLD, ierr )
	call MPI_Check( 'global_all_min_1r:MPI_Allreduce', ierr )

	return
	end subroutine global_all_min_1r

        subroutine global_all_min_2r( lbuf, gbuf, mroot, idebug )
        real, intent(in), dimension(:,:) :: lbuf
        real, intent(out), dimension(:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

	integer :: i,j,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

	call assert( size(lbuf,2).eq.size(gbuf,2),  &
		'** global_all_min_2r: size(lbuf,2).ne.size(gbuf,2) ', &
		size(lbuf,2), size(gbuf,2) )

	do j=lbound(lbuf,2),ubound(lbuf,2)
	  call global_all_min_1r( lbuf(:,j), gbuf(:,j), lroot, lidebug )
	enddo

	return
	end subroutine global_all_min_2r

        subroutine global_all_min_3r( lbuf, gbuf, mroot, idebug )
        real, intent(in), dimension(:,:,:) :: lbuf
        real, intent(out), dimension(:,:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        integer :: j,k,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

        call assert( size(lbuf,2).eq.size(gbuf,2),  &
                '** global_all_min_3i: size(lbuf,2).ne.size(gbuf,2) ', &
                size(lbuf,2), size(gbuf,2) )

        call assert( size(lbuf,3).eq.size(gbuf,3),  &
                '** global_all_min_3i: size(lbuf,3).ne.size(gbuf,3) ', &
                size(lbuf,3), size(gbuf,3) )

        do k=lbound(lbuf,3),ubound(lbuf,3)
        do j=lbound(lbuf,2),ubound(lbuf,2)
          call global_all_min_1r( lbuf(:,j,k), gbuf(:,j,k), lroot, lidebug )
        enddo
	enddo

        return
        end subroutine global_all_min_3r

        subroutine global_all_min_0d( lbuf, gbuf, mroot, idebug )
        double precision, intent(in) :: lbuf
        double precision, intent(out) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

        if (.not. present(mroot)) then
           lroot = 0
        else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

        recvtype = MPI_DOUBLE_PRECISION
        sendtype = recvtype

        sendcnt = 1

        call MPI_Allreduce( lbuf, gbuf, sendcnt, sendtype, MPI_MIN, &
                          MPI_COMM_WORLD, ierr )
        call MPI_Check( 'global_all_min_0d:MPI_Allreduce', ierr )

        return
        end subroutine global_all_min_0d


	subroutine global_all_min_1d( lbuf, gbuf, mroot, idebug )
        double precision, intent(in), dimension(:) :: lbuf       
        double precision, intent(out), dimension(:) :: gbuf
	integer, optional, intent(in) :: mroot, idebug

	integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

	recvtype = MPI_DOUBLE_PRECISION
	sendtype = recvtype

	sendcnt = size(lbuf)

        call MPI_Allreduce( lbuf, gbuf, sendcnt, sendtype,  MPI_MIN, &
                          MPI_COMM_WORLD, ierr )
	call MPI_Check( 'global_all_min_1d:MPI_Allreduce', ierr )

	return
	end subroutine global_all_min_1d

        subroutine global_all_min_2d( lbuf, gbuf, mroot, idebug )
        double precision, intent(in), dimension(:,:) :: lbuf
        double precision, intent(out), dimension(:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

	integer :: i,j,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

	call assert( size(lbuf,2).eq.size(gbuf,2),  &
		'** global_all_min_2d: size(lbuf,2).ne.size(gbuf,2) ', &
		size(lbuf,2), size(gbuf,2) )

	do j=lbound(lbuf,2),ubound(lbuf,2)
	  call global_all_min_1d( lbuf(:,j), gbuf(:,j), lroot, lidebug )
	enddo

	return
	end subroutine global_all_min_2d

        subroutine global_all_min_3d( lbuf, gbuf, mroot, idebug )
        double precision, intent(in), dimension(:,:,:) :: lbuf
        double precision, intent(out), dimension(:,:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        integer :: j,k,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

        call assert( size(lbuf,2).eq.size(gbuf,2),  &
                '** global_all_min_3i: size(lbuf,2).ne.size(gbuf,2) ', &
                size(lbuf,2), size(gbuf,2) )

        call assert( size(lbuf,3).eq.size(gbuf,3),  &
                '** global_all_min_3i: size(lbuf,3).ne.size(gbuf,3) ', &
                size(lbuf,3), size(gbuf,3) )

        do k=lbound(lbuf,3),ubound(lbuf,3)
        do j=lbound(lbuf,2),ubound(lbuf,2)
          call global_all_min_1d( lbuf(:,j,k), gbuf(:,j,k), lroot, lidebug )
        enddo
	enddo

        return
        end subroutine global_all_min_3d

        subroutine global_max_0i( lbuf, gbuf, mroot, idebug )
        integer, intent(in) :: lbuf
        integer, intent(out) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

        if (.not. present(mroot)) then
           lroot = 0
        else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

        recvtype = MPI_INTEGER
        sendtype = recvtype

        sendcnt = 1

        call MPI_Reduce( lbuf, gbuf, sendcnt, sendtype,  MPI_MAX, &
                          lroot, MPI_COMM_WORLD, ierr )
        call MPI_Check( 'global_max_0i:MPI_Reduce', ierr )

        return
        end subroutine global_max_0i


	subroutine global_max_1i( lbuf, gbuf, mroot, idebug )
        integer, intent(in), dimension(:) :: lbuf       
        integer, intent(out), dimension(:) :: gbuf
	integer, optional, intent(in) :: mroot, idebug

	integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

	recvtype = MPI_INTEGER
	sendtype = recvtype

	sendcnt = size(lbuf)

        call MPI_Reduce( lbuf, gbuf, sendcnt, sendtype,  MPI_MAX, &
                          lroot, MPI_COMM_WORLD, ierr )
	call MPI_Check( 'global_max_1i:MPI_Reduce', ierr )

	return
	end subroutine global_max_1i

        subroutine global_max_2i( lbuf, gbuf, mroot, idebug )
        integer, intent(in), dimension(:,:) :: lbuf
        integer, intent(out), dimension(:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

	integer :: i,j,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

	if(myPE.eq.lroot) then
	call assert( size(lbuf,2).eq.size(gbuf,2),  &
		'** global_max_2i: size(lbuf,2).ne.size(gbuf,2) ', &
		size(lbuf,2), size(gbuf,2) )
	endif

	do j=lbound(lbuf,2),ubound(lbuf,2)
	  call global_max_1i( lbuf(:,j), gbuf(:,j), lroot, lidebug )
	enddo

	return
	end subroutine global_max_2i

        subroutine global_max_3i( lbuf, gbuf, mroot, idebug )
        integer, intent(in), dimension(:,:,:) :: lbuf
        integer, intent(out), dimension(:,:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        integer :: j,k,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

	if(myPE.eq.lroot) then
        call assert( size(lbuf,2).eq.size(gbuf,2),  &
                '** global_max_3i: size(lbuf,2).ne.size(gbuf,2) ', &
                size(lbuf,2), size(gbuf,2) )

        call assert( size(lbuf,3).eq.size(gbuf,3),  &
                '** global_max_3i: size(lbuf,3).ne.size(gbuf,3) ', &
                size(lbuf,3), size(gbuf,3) )
	endif

        do k=lbound(lbuf,3),ubound(lbuf,3)
        do j=lbound(lbuf,2),ubound(lbuf,2)
          call global_max_1i( lbuf(:,j,k), gbuf(:,j,k), lroot, lidebug )
        enddo
	enddo

        return
        end subroutine global_max_3i

        subroutine global_max_0r( lbuf, gbuf, mroot, idebug )
        real, intent(in) :: lbuf
        real, intent(out) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

        if (.not. present(mroot)) then
           lroot = 0
        else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

        recvtype = MPI_REAL
        sendtype = recvtype

        sendcnt = 1

        call MPI_Reduce( lbuf, gbuf, sendcnt, sendtype,  MPI_MAX, &
                          lroot, MPI_COMM_WORLD, ierr )
        call MPI_Check( 'global_max_0r:MPI_Reduce', ierr )

        return
        end subroutine global_max_0r


	subroutine global_max_1r( lbuf, gbuf, mroot, idebug )
        real, intent(in), dimension(:) :: lbuf       
        real, intent(out), dimension(:) :: gbuf
	integer, optional, intent(in) :: mroot, idebug

	integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

	recvtype = MPI_REAL
	sendtype = recvtype

	sendcnt = size(lbuf)

        call MPI_Reduce( lbuf, gbuf, sendcnt, sendtype,  MPI_MAX, &
                          lroot, MPI_COMM_WORLD, ierr )
	call MPI_Check( 'global_max_1r:MPI_Reduce', ierr )

	return
	end subroutine global_max_1r

        subroutine global_max_2r( lbuf, gbuf, mroot, idebug )
        real, intent(in), dimension(:,:) :: lbuf
        real, intent(out), dimension(:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

	integer :: i,j,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

	if(myPE.eq.lroot) then
	call assert( size(lbuf,2).eq.size(gbuf,2),  &
		'** global_max_2r: size(lbuf,2).ne.size(gbuf,2) ', &
		size(lbuf,2), size(gbuf,2) )
	endif

	do j=lbound(lbuf,2),ubound(lbuf,2)
	  call global_max_1r( lbuf(:,j), gbuf(:,j), lroot, lidebug )
	enddo

	return
	end subroutine global_max_2r

        subroutine global_max_3r( lbuf, gbuf, mroot, idebug )
        real, intent(in), dimension(:,:,:) :: lbuf
        real, intent(out), dimension(:,:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        integer :: j,k,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

	if(myPE.eq.lroot) then
        call assert( size(lbuf,2).eq.size(gbuf,2),  &
                '** global_max_3i: size(lbuf,2).ne.size(gbuf,2) ', &
                size(lbuf,2), size(gbuf,2) )

        call assert( size(lbuf,3).eq.size(gbuf,3),  &
                '** global_max_3i: size(lbuf,3).ne.size(gbuf,3) ', &
                size(lbuf,3), size(gbuf,3) )
	endif

        do k=lbound(lbuf,3),ubound(lbuf,3)
        do j=lbound(lbuf,2),ubound(lbuf,2)
          call global_max_1r( lbuf(:,j,k), gbuf(:,j,k), lroot, lidebug )
        enddo
	enddo

        return
        end subroutine global_max_3r

        subroutine global_max_0d( lbuf, gbuf, mroot, idebug )
        double precision, intent(in) :: lbuf
        double precision, intent(out) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

        if (.not. present(mroot)) then
           lroot = 0
        else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

        recvtype = MPI_DOUBLE_PRECISION
        sendtype = recvtype

        sendcnt = 1

        call MPI_Reduce( lbuf, gbuf, sendcnt, sendtype,  MPI_MAX, &
                          lroot, MPI_COMM_WORLD, ierr )
        call MPI_Check( 'global_max_0d:MPI_Reduce', ierr )

        return
        end subroutine global_max_0d


	subroutine global_max_1d( lbuf, gbuf, mroot, idebug )
        double precision, intent(in), dimension(:) :: lbuf       
        double precision, intent(out), dimension(:) :: gbuf
	integer, optional, intent(in) :: mroot, idebug

	integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

	recvtype = MPI_DOUBLE_PRECISION
	sendtype = recvtype

	sendcnt = size(lbuf)

        call MPI_Reduce( lbuf, gbuf, sendcnt, sendtype,  MPI_MAX, &
                          lroot, MPI_COMM_WORLD, ierr )
	call MPI_Check( 'global_max_1d:MPI_Reduce', ierr )

	return
	end subroutine global_max_1d

        subroutine global_max_2d( lbuf, gbuf, mroot, idebug )
        double precision, intent(in), dimension(:,:) :: lbuf
        double precision, intent(out), dimension(:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

	integer :: i,j,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

	if(myPE.eq.lroot) then
	call assert( size(lbuf,2).eq.size(gbuf,2),  &
		'** global_max_2d: size(lbuf,2).ne.size(gbuf,2) ', &
		size(lbuf,2), size(gbuf,2) )
	endif

	do j=lbound(lbuf,2),ubound(lbuf,2)
	  call global_max_1d( lbuf(:,j), gbuf(:,j), lroot, lidebug )
	enddo

	return
	end subroutine global_max_2d

        subroutine global_max_3d( lbuf, gbuf, mroot, idebug )
        double precision, intent(in), dimension(:,:,:) :: lbuf
        double precision, intent(out), dimension(:,:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        integer :: j,k,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

	if(myPE.eq.lroot) then
        call assert( size(lbuf,2).eq.size(gbuf,2),  &
                '** global_max_3i: size(lbuf,2).ne.size(gbuf,2) ', &
                size(lbuf,2), size(gbuf,2) )

        call assert( size(lbuf,3).eq.size(gbuf,3),  &
                '** global_max_3i: size(lbuf,3).ne.size(gbuf,3) ', &
                size(lbuf,3), size(gbuf,3) )
	endif

        do k=lbound(lbuf,3),ubound(lbuf,3)
        do j=lbound(lbuf,2),ubound(lbuf,2)
          call global_max_1d( lbuf(:,j,k), gbuf(:,j,k), lroot, lidebug )
        enddo
	enddo

        return
        end subroutine global_max_3d

        subroutine global_all_max_0i( lbuf, gbuf, mroot, idebug )
        integer, intent(in) :: lbuf
        integer, intent(out) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

        if (.not. present(mroot)) then
           lroot = 0
        else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

        recvtype = MPI_INTEGER
        sendtype = recvtype

        sendcnt = 1

        call MPI_Allreduce( lbuf, gbuf, sendcnt, sendtype, MPI_MAX, &
                          MPI_COMM_WORLD, ierr )
        call MPI_Check( 'global_all_max_0i:MPI_Allreduce', ierr )

        return
        end subroutine global_all_max_0i


	subroutine global_all_max_1i( lbuf, gbuf, mroot, idebug )
        integer, intent(in), dimension(:) :: lbuf       
        integer, intent(out), dimension(:) :: gbuf
	integer, optional, intent(in) :: mroot, idebug

	integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

	recvtype = MPI_INTEGER
	sendtype = recvtype

	sendcnt = size(lbuf)

        call MPI_Allreduce( lbuf, gbuf, sendcnt, sendtype, MPI_MAX, &
                          MPI_COMM_WORLD, ierr )
	call MPI_Check( 'global_all_max_1i:MPI_Allreduce', ierr )

	return
	end subroutine global_all_max_1i

        subroutine global_all_max_2i( lbuf, gbuf, mroot, idebug )
        integer, intent(in), dimension(:,:) :: lbuf
        integer, intent(out), dimension(:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

	integer :: i,j,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

	call assert( size(lbuf,2).eq.size(gbuf,2),  &
		'** global_all_max_2i: size(lbuf,2).ne.size(gbuf,2) ', &
		size(lbuf,2), size(gbuf,2) )

	do j=lbound(lbuf,2),ubound(lbuf,2)
	  call global_all_max_1i( lbuf(:,j), gbuf(:,j), lroot, lidebug )
	enddo

	return
	end subroutine global_all_max_2i

        subroutine global_all_max_3i( lbuf, gbuf, mroot, idebug )
        integer, intent(in), dimension(:,:,:) :: lbuf
        integer, intent(out), dimension(:,:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        integer :: j,k,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

        call assert( size(lbuf,2).eq.size(gbuf,2),  &
                '** global_all_max_3i: size(lbuf,2).ne.size(gbuf,2) ', &
                size(lbuf,2), size(gbuf,2) )

        call assert( size(lbuf,3).eq.size(gbuf,3),  &
                '** global_all_max_3i: size(lbuf,3).ne.size(gbuf,3) ', &
                size(lbuf,3), size(gbuf,3) )

        do k=lbound(lbuf,3),ubound(lbuf,3)
        do j=lbound(lbuf,2),ubound(lbuf,2)
          call global_all_max_1i( lbuf(:,j,k), gbuf(:,j,k), lroot, lidebug )
        enddo
	enddo

        return
        end subroutine global_all_max_3i

        subroutine global_all_max_0r( lbuf, gbuf, mroot, idebug )
        real, intent(in) :: lbuf
        real, intent(out) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

        if (.not. present(mroot)) then
           lroot = 0
        else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

        recvtype = MPI_REAL
        sendtype = recvtype

        sendcnt = 1

        call MPI_Allreduce( lbuf, gbuf, sendcnt, sendtype, MPI_MAX, &
                          MPI_COMM_WORLD, ierr )
        call MPI_Check( 'global_all_max_0r:MPI_Allreduce', ierr )


        return
        end subroutine global_all_max_0r


	subroutine global_all_max_1r( lbuf, gbuf, mroot, idebug )
        real, intent(in), dimension(:) :: lbuf       
        real, intent(out), dimension(:) :: gbuf
	integer, optional, intent(in) :: mroot, idebug

	integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif  

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

	recvtype = MPI_REAL
	sendtype = recvtype

	sendcnt = size(lbuf)

        call MPI_Allreduce( lbuf, gbuf, sendcnt, sendtype, MPI_MAX, &
                            MPI_COMM_WORLD, ierr )
	call MPI_Check( 'global_all_max_1r:MPI_Allreduce', ierr )

	return
	end subroutine global_all_max_1r

        subroutine global_all_max_2r( lbuf, gbuf, mroot, idebug )
        real, intent(in), dimension(:,:) :: lbuf
        real, intent(out), dimension(:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

	integer :: i,j,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif 

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

	call assert( size(lbuf,2).eq.size(gbuf,2),  &
		'** global_all_max_2r: size(lbuf,2).ne.size(gbuf,2) ', &
		size(lbuf,2), size(gbuf,2) )

	do j=lbound(lbuf,2),ubound(lbuf,2)
	  call global_all_max_1r( lbuf(:,j), gbuf(:,j), lroot, lidebug )
	enddo

	return
	end subroutine global_all_max_2r

        subroutine global_all_max_3r( lbuf, gbuf, mroot, idebug )
        real, intent(in), dimension(:,:,:) :: lbuf
        real, intent(out), dimension(:,:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        integer :: j,k,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif 

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

        call assert( size(lbuf,2).eq.size(gbuf,2),  &
                '** global_all_max_3i: size(lbuf,2).ne.size(gbuf,2) ', &
                size(lbuf,2), size(gbuf,2) )

        call assert( size(lbuf,3).eq.size(gbuf,3),  &
                '** global_all_max_3i: size(lbuf,3).ne.size(gbuf,3) ', &
                size(lbuf,3), size(gbuf,3) )

        do k=lbound(lbuf,3),ubound(lbuf,3)
        do j=lbound(lbuf,2),ubound(lbuf,2)
          call global_all_max_1r( lbuf(:,j,k), gbuf(:,j,k), lroot, lidebug )
        enddo
	enddo

        return
        end subroutine global_all_max_3r

        subroutine global_all_max_0d( lbuf, gbuf, mroot, idebug )
        double precision, intent(in) :: lbuf
        double precision, intent(out) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

        if (.not. present(mroot)) then
           lroot = 0
        else
           lroot = mroot
        endif 

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

        recvtype = MPI_DOUBLE_PRECISION
        sendtype = recvtype

        sendcnt = 1

        call MPI_Allreduce( lbuf, gbuf, sendcnt, sendtype, MPI_MAX, &
                          MPI_COMM_WORLD, ierr )
        call MPI_Check( 'global_all_max_0d:MPI_Allreduce', ierr )

        return
        end subroutine global_all_max_0d


	subroutine global_all_max_1d( lbuf, gbuf, mroot, idebug )
        double precision, intent(in), dimension(:) :: lbuf       
        double precision, intent(out), dimension(:) :: gbuf
	integer, optional, intent(in) :: mroot, idebug

	integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif 

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

	recvtype = MPI_DOUBLE_PRECISION
	sendtype = recvtype

	sendcnt = size(lbuf)

        call MPI_Allreduce( lbuf, gbuf, sendcnt, sendtype,  MPI_MAX, &
                          MPI_COMM_WORLD, ierr )
	call MPI_Check( 'global_all_max_1d:MPI_Allreduce', ierr )

	return
	end subroutine global_all_max_1d

        subroutine global_all_max_2d( lbuf, gbuf, mroot, idebug )
        double precision, intent(in), dimension(:,:) :: lbuf
        double precision, intent(out), dimension(:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

	integer :: i,j,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif 

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

	call assert( size(lbuf,2).eq.size(gbuf,2),  &
		'** global_all_max_2d: size(lbuf,2).ne.size(gbuf,2) ', &
		size(lbuf,2), size(gbuf,2) )

	do j=lbound(lbuf,2),ubound(lbuf,2)
	  call global_all_max_1d( lbuf(:,j), gbuf(:,j), lroot, lidebug )
	enddo

	return
	end subroutine global_all_max_2d

        subroutine global_all_max_3d( lbuf, gbuf, mroot, idebug )
        double precision, intent(in), dimension(:,:,:) :: lbuf
        double precision, intent(out), dimension(:,:,:) :: gbuf
        integer, optional, intent(in) :: mroot, idebug

        integer :: j,k,lroot, lidebug

	if (.not. present(mroot)) then
	   lroot = 0
	else
           lroot = mroot
        endif 

        if (.not. present(idebug)) then  
           lidebug = 0
        else
           lidebug = idebug
        endif

        call assert( size(lbuf,2).eq.size(gbuf,2),  &
                '** global_all_max_3i: size(lbuf,2).ne.size(gbuf,2) ', &
                size(lbuf,2), size(gbuf,2) )

        call assert( size(lbuf,3).eq.size(gbuf,3),  &
                '** global_all_max_3i: size(lbuf,3).ne.size(gbuf,3) ', &
                size(lbuf,3), size(gbuf,3) )

        do k=lbound(lbuf,3),ubound(lbuf,3)
        do j=lbound(lbuf,2),ubound(lbuf,2)
          call global_all_max_1d( lbuf(:,j,k), gbuf(:,j,k), lroot, lidebug )
        enddo
	enddo

        return
        end subroutine global_all_max_3d

        subroutine global_all_and_0d( lvalue, gvalue, mroot, idebug )
        logical, intent(in) :: lvalue
        logical, intent(out) :: gvalue
        integer, optional, intent(in) :: mroot, idebug

!       ---------------
!       local variables
!       ---------------
        integer :: ierror, icount
        integer :: lroot, lidebug

        if (.not. present(mroot)) then
           lroot = 0
        else
           lroot = mroot
        endif

        if (.not. present(idebug)) then
           lidebug = 0
        else
           lidebug = idebug
        endif

        icount = 1

        call MPI_Allreduce( lvalue, gvalue, icount, MPI_LOGICAL, &
                        MPI_LAND, MPI_COMM_WORLD, ierror )

        call MPI_Check( 'global_all_and_0d ', ierror )
        return
        end subroutine  global_all_and_0d


        subroutine global_all_and_1d( lvalue, gvalue, mroot, idebug )
        logical, intent(in), dimension(:) :: lvalue
        logical, intent(out), dimension(:) :: gvalue
        integer, optional, intent(in) :: mroot, idebug

!       ---------------
!       local variables
!       ---------------
        integer :: ierror, icount
        integer :: lroot, lidebug

        if (.not. present(mroot)) then
           lroot = 0
        else
           lroot = mroot
        endif

        if (.not. present(idebug)) then
           lidebug = 0
        else
           lidebug = idebug
        endif


        icount = size( lvalue )

        call MPI_Allreduce( lvalue, gvalue, icount, MPI_LOGICAL, &
                        MPI_LAND, MPI_COMM_WORLD, ierror )

        call MPI_Check( 'global_all_and_1d ', ierror )
        return
        end subroutine global_all_and_1d



        subroutine global_all_or_0d( lvalue, gvalue, mroot, idebug )
        logical, intent(in) :: lvalue
        logical, intent(out) :: gvalue
        integer, optional, intent(in) :: mroot, idebug

!       ---------------
!       local variables
!       ---------------
        integer :: ierror, icount
        integer :: lroot, lidebug

        if (.not. present(mroot)) then
           lroot = 0
        else
           lroot = mroot
        endif

        if (.not. present(idebug)) then
           lidebug = 0
        else
           lidebug = idebug
        endif


        icount = 1

        call MPI_Allreduce( lvalue, gvalue, icount, MPI_LOGICAL, &
                        MPI_LOR, MPI_COMM_WORLD, ierror )

        call MPI_Check( 'global_all_or_0d ', ierror )
        return
        end subroutine global_all_or_0d


        subroutine global_all_or_1d( lvalue, gvalue, mroot, idebug )
        logical, intent(in), dimension(:) :: lvalue
        logical, intent(out), dimension(:) :: gvalue
        integer, optional, intent(in) :: mroot, idebug

!       ---------------
!       local variables
!       ---------------
        integer :: ierror, icount
        integer :: lroot, lidebug

        if (.not. present(mroot)) then
           lroot = 0
        else
           lroot = mroot
        endif

        if (.not. present(idebug)) then
           lidebug = 0
        else
           lidebug = idebug
        endif


        icount = size( lvalue )

        call MPI_Allreduce( lvalue, gvalue, icount, MPI_LOGICAL, &
                        MPI_LOR, MPI_COMM_WORLD, ierror )

        call MPI_Check( 'global_all_or_1d ', ierror )
        return
        end subroutine global_all_or_1d


        subroutine exitMPI(myid)
        USE funits
        integer, optional, intent(in) :: myid

        INTEGER :: mylid

        if (.not. present(myid)) then
           mylid = myPE
        else
           mylid = myid
        endif

        write(*,100) mylid
        write(UNIT_LOG,100) mylid

100   	format(/,'*****************',&
          '********************************************',/, &
         '(PE ',I2,') : A FATAL ERROR OCCURRED ',/,9X, &
         'ABORTING ALL PROCESSES ',/,'*****************', &
         '********************************************',/)
  
        call MPI_Barrier(MPI_COMM_WORLD,mpierr)
        write(*,"('(PE ',I2,') : MPI_Barrier return = ',I2)") &
               mylid,mpierr
  
        call MPI_Finalize(mpierr)

        STOP 'MPI terminated'
        END subroutine exitMPI

	end module mpi_utility


