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
                send_recv_1d, &
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

!	-------------------------------------
!	set up tables and data structures for
!	exchanging ghost regions
!	-------------------------------------

!	---------------
!	local variables
!	---------------
	character(len=80), parameter :: name = 'sendrecv_init'

	character(len=80), allocatable, dimension(:) :: line
	integer :: ip, lmax

	integer :: lidebug
	integer :: isize,jsize,ksize, ijksize

	integer :: iter, i,j,k, ii, jj,kk, &
		ntotal, icount,ipos, &
		isize, ilayer,        i1,i2,  j1,j2, k1,k2,  &
		ijk, ijk2, iproc, jproc, src,dest, &
		ierror, &
	        kstart1,kend1,  jstart1,jend1, istart1,iend1, &
	        kstart2,kend2,  jstart2,jend2, istart2,iend2, &
	        kstart3,kend3,  jstart3,jend3, istart3,iend3

	logical :: isok, isvalid, ismine, is_halobc

	integer, dimension(:,:,:), allocatable :: ijk2proc
        integer, pointer, dimension(:) :: &
		istartx,iendx, jstartx,jendx, kstartx,kendx, &
		ncount, &
		recvproc, recvtag, xrecv, recvijk,  &
		sendproc, sendtag, xsend, sendijk



	integer, parameter :: message_tag_offset = 1000


!	----------------
!	inline functions
!	----------------
	integer :: message_tag

	message_tag(src,dest) = message_tag_offset + (1+src + dest*10*numPEs)

	include 'function.inc'


!	--------------------
!	initialize variables
!	--------------------
	lidebug = 0
	if (present(idebug)) then
	  lidebug = idebug
	endif

	communicator = comm
	call MPI_COMM_SIZE( comm, numPEs, ierror )
	call MPI_Check( 'sendrecv_init:MPI_COMM_SIZE ', ierror )

	call MPI_COMM_RANK( comm, myPE, ierror )
	call MPI_Check( 'sendrecv_init:MPI_COMM_RANK ', ierror )

!	---------------------------
!	determine bounds of domain
!	---------------------------


	kstart1 = minval( kstart1_all(:) )
	kstart2 = minval( kstart2_all(:) )
	kstart3 = minval( kstart3_all(:) )
	kend1 = maxval( kend1_all(:) )
	kend2 = maxval( kend2_all(:) )
	kend3 = maxval( kend3_all(:) )

	jstart1 = minval( jstart1_all(:) )
	jstart2 = minval( jstart2_all(:) )
	jstart3 = minval( jstart3_all(:) )
	jend1 = maxval( jend1_all(:) )
	jend2 = maxval( jend2_all(:) )
	jend3 = maxval( jend3_all(:) )


	istart1 = minval( istart1_all(:) )
	istart2 = minval( istart2_all(:) )
	istart3 = minval( istart3_all(:) )
	iend1 = maxval( iend1_all(:) )
	iend2 = maxval( iend2_all(:) )
	iend3 = maxval( iend3_all(:) )


	call assert( jstart1 .le. jend1, &
		'** sendrecv_init: jstart1,jend1 ', jstart1,jend1 )
	call assert( jstart2 .le. jend2, &
		'** sendrecv_init: jstart2,jend2 ', jstart2,jend2 )
	call assert( jstart3 .le. jend3, &
		'** sendrecv_init: jstart3,jend3 ', jstart3,jend3 )
	
	call assert( kstart1 .le. kend1, &
		'** sendrecv_init: kstart1,kend1 ', kstart1,kend1 )
	call assert( kstart2 .le. kend2, &
		'** sendrecv_init: kstart2,kend2 ', kstart2,kend2 )
	call assert( kstart3 .le. kend3, &
		'** sendrecv_init: kstart3,kend3 ', kstart3,kend3 )
	
	call assert( istart1 .le. iend1, &
		'** sendrecv_init: istart1,iend1 ', istart1,iend1 )
	call assert( istart2 .le. iend2, &
		'** sendrecv_init: istart2,iend2 ', istart2,iend2 )
	call assert( istart3 .le. iend3, &
		'** sendrecv_init: istart3,iend3 ', istart3,iend3 )
	
	
	
	




	k1 = min( kstart1, min(kstart2, kstart3) )
	k2 = max( kend1, max(kend2, kend3) )
	j1 = min( jstart1, min(jstart2, jstart3) )
	j2 = max( jend1, max(jend2, jend3) )
	i1 = min( istart1, min(istart2, istart3) )
	i2 = max( iend1, max(iend2, iend3) )

	allocate( ijk2proc( i1:i2, j1:j2, k1:k2 ) )


	if(localfunc) then
!	--------------------------------------
!	double check ijk_of()
!	--------------------------------------

	do k=kstart3_all(myPE),kend3_all(myPE)
	do j=jstart3_all(myPE),jend3_all(myPE)
	do i=istart3_all(myPE),iend3_all(myPE)
	   ijk = funijk(i,j,k)
	   call ijk_of(ijk, ii,jj,kk)
	   ijk2 = funijk( ii,jj,kk)

	   isvalid = (ii.eq.i).and.(jj.eq.j).and.(kk.eq.k).and.(ijk.eq.ijk2)
	   if (.not.isvalid) then
		call write_debug( name, 'error with ijk_of ')

		call write_debug( name, 'istart3_all(myPE),iend3_all(myPE) ', &
				istart3_all(myPE),iend3_all(myPE) )
		call write_debug( name, 'jstart3_all(myPE),jend3_all(myPE) ', &
				jstart3_all(myPE),jend3_all(myPE) )
		call write_debug( name, 'kstart3_all(myPE),kend3_all(myPE) ', &
				kstart3_all(myPE),kend3_all(myPE) )

		call write_debug( name, 'i,j,k, ijk ', i,j,k, ijk )
	        call write_debug( name, 'ii,jj,kk,  ijk2 ',ii,jj,kk,ijk2 )

	   endif
	enddo
	enddo
	enddo
	endif ! Local Function

	
	if (lidebug.ge.1) then
	   call write_debug( name, 'imap ', imap )
	   call write_debug( name, 'jmap ', jmap )
	   call write_debug( name, 'kmap ', kmap )
	endif
	    

!	----------------------------
!	set up table ijk2proc(:,:,:)
!
!	ijk2proc(i,j,k) maps (i,j,k) index to
!	unique processor that 'owns' that node.
!	----------------------------



	ijk2proc( :,:,: ) = 0

!	--------------------------------------------------
!	double check domain decomposition that
!	each interior node is assigned to UNIQUE processor
!	--------------------------------------------------
	do iproc=0,numPEs-1

	   i1 = istart1_all(iproc)
	   i2 = iend1_all(iproc)
	   j1 = jstart1_all(iproc)
	   j2 = jend1_all(iproc)
	   k1 = kstart1_all(iproc)
	   k2 = kend1_all(iproc)


	   do k=k1,k2
	   do j=j1,j2
	   do i=i1,i2
		ijk2proc(i,j,k) = ijk2proc(i,j,k) + 1
	   enddo
	   enddo
	   enddo

	enddo

	do k=kstart1,kend1
	do j=jstart1,jend1
	do i=istart1,iend1
	  isvalid = (ijk2proc(i,j,k) .eq. 1)
	  if (.not.isvalid) then

		call write_debug(name, ' invalid decomposition ')
		call write_debug(name, 'i,j,k ',i,j,k )
		call write_debug(name, 'ijk2proc(i,j,k) ', ijk2proc(i,j,k))

		stop '** error ** '
	  endif
	enddo
	enddo
	enddo

	ijk2proc(:,:,:) = -1
	do iproc=0,numPEs-1
           i1 = istart1_all(iproc)
           i2 = iend1_all(iproc)
           j1 = jstart1_all(iproc)
           j2 = jend1_all(iproc)
           k1 = kstart1_all(iproc)
           k2 = kend1_all(iproc)


	   do k=k1,k2
	   do j=j1,j2
	   do i=i1,i2
		ijk2proc(i,j,k) = iproc
	   enddo
	   enddo
	   enddo

	enddo

	   

   allocate( ncount(0:numPEs-1) )

   allocate( istartx(0:numPEs-1) )
   allocate( jstartx(0:numPEs-1) )
   allocate( kstartx(0:numPEs-1) )

   allocate( iendx(0:numPEs-1) )
   allocate( jendx(0:numPEs-1) )
   allocate( kendx(0:numPEs-1) )

do ilayer=1,2

   if (ilayer.eq.1) then
	kstartx(:) = kstart2_all(:)
	kendx(:) = kend2_all(:)
	jstartx(:) = jstart2_all(:)
	jendx(:) = jend2_all(:)
	istartx(:) = istart2_all(:)
	iendx(:) = iend2_all(:)
   else
	kstartx(:) = kstart3_all(:)
	kendx(:) = kend3_all(:)
	jstartx(:) = jstart3_all(:)
	jendx(:) = jend3_all(:)
	istartx(:) = istart3_all(:)
	iendx(:) = iend3_all(:)
   endif



   if (lidebug.ge.1) then
	call write_debug(name, 'determine send schedule ', myPE )
   endif

!  -----------------------
!  determine send schedule
!
!  examine all neighboring processors
!  to see if they need my data
!  -----------------------


! -----------------------------------
! first pass to determine array sizes
! -----------------------------------

   ncount(:) = 0

   do iproc=0,numPEs-1
    if (iproc.ne.myPE) then


	k1 = lbound(ijk2proc,3)
	k2 = ubound(ijk2proc,3)
	j1 = lbound(ijk2proc,2)
	j2 = ubound(ijk2proc,2)
	i1 = lbound(ijk2proc,1)
	i2 = ubound(ijk2proc,1)
	

	do k=kstartx(iproc),kendx(iproc)
	do j=jstartx(iproc),jendx(iproc)
	do i=istartx(iproc),iendx(iproc)

	  ii = imap(i)
	  jj = jmap(j)
	  kk = kmap(k)

	  isvalid  = (k1.le.kk).and.(kk.le.k2)
	  call assert( isvalid, '** sendrecv_init: invalid kk ', kk )

	  isvalid  = (j1.le.jj).and.(jj.le.j2)
	  call assert( isvalid, '** sendrecv_init: invalid jj ', jj )

	  isvalid  = (i1.le.ii).and.(ii.le.i2)
	  call assert( isvalid, '** sendrecv_init: invalid ii ', ii )

	  jproc = ijk2proc( ii,jj,kk )

	  ismine = (jproc .eq. myPE) 
	  if (ismine) then
	        ncount(iproc) = ncount(iproc) + 1
	  endif

	enddo
	enddo
	enddo

    endif
   enddo 


!  --------------
!  prepare arrays
!  --------------
   ntotal = 0
   nsend = 0
   do iproc=0,numPEs-1
      ntotal = ntotal + ncount(iproc)
      if (ncount(iproc).ge.1) then
	nsend = nsend + 1
      endif
   enddo

   if (lidebug.ge.1) then
	call write_debug( name, 'ncount = ', ncount )
	call write_debug( name, 'nsend, ntotal ', nsend, ntotal )
   endif


   allocate( xsend(nsend+1) )

   allocate( sendijk( max(1,ntotal) ) )
   allocate( sendproc(max(1,nsend)) )

   nsend = 0
   do iproc=0,numPEs-1
     if (ncount(iproc).ne.0) then
        nsend = nsend + 1
        sendproc(nsend) = iproc
     endif
   enddo

   xsend(1) = 1
   do i=1,nsend
     iproc = sendproc(i)
     xsend(i+1) = xsend(i) + ncount(iproc)
   enddo

   allocate( sendtag( max(1,nsend) ) )
   do ii=1,nsend
     iproc = sendproc(ii)
     src = myPE
     dest = iproc
     sendtag(ii) = message_tag( src, dest )
   enddo

! -----------------------------
! second pass to fill in arrays
! -----------------------------


  ipos = 1
  do iter=1,nsend
     iproc = sendproc(iter)
     icount = 0

        do k=kstartx(iproc),kendx(iproc)
        do j=jstartx(iproc),jendx(iproc)
        do i=istartx(iproc),iendx(iproc)

          
	  ii = imap(i)
	  jj = jmap(j)
	  kk = kmap(k)
	  jproc = ijk2proc(ii,jj,kk)
          ismine = (jproc.eq.myPE)
          if (ismine) then
              icount = icount + 1
	      ijk = funijk(ii,jj,kk)

	      ipos = xsend(iter)-1 + icount
	      sendijk( ipos ) = ijk
          endif

       enddo
       enddo
       enddo

     isvalid = (icount .eq. ncount(iproc))
     call assert( isvalid, '** sendrecv_init: icount != ncount(iproc) ', iproc)

   enddo



   if (lidebug.ge.1) then
	call write_debug(name, 'determine recv schedule ', myPE )
   endif

! ---------------------------
! determine recv schedule
!
! examine nodes in my ghost region and
! see what data is needed from my neighbors
! ---------------------------

! -----------------------------------
! first pass to determine array sizes
! -----------------------------------

  ncount(:) = 0

        k1 = lbound(ijk2proc,3)
        k2 = ubound(ijk2proc,3)
        j1 = lbound(ijk2proc,2)
        j2 = ubound(ijk2proc,2)
        i1 = lbound(ijk2proc,1)
        i2 = ubound(ijk2proc,1)

      
	do k=kstartx(myPE),kendx(myPE)
	do j=jstartx(myPE),jendx(myPE)
	do i=istartx(myPE),iendx(myPE)

	  ii = imap(i)
	  jj = jmap(j)
	  kk = kmap(k)

          isvalid  = (k1.le.kk).and.(kk.le.k2)
          call assert( isvalid, '** sendrecv_init: invalid kk ', kk )

          isvalid  = (j1.le.jj).and.(jj.le.j2)
          call assert( isvalid, '** sendrecv_init: invalid jj ', jj )

          isvalid  = (i1.le.ii).and.(ii.le.i2)
          call assert( isvalid, '** sendrecv_init: invalid ii ', ii )





	  iproc = ijk2proc(ii,jj,kk)
	  is_halobc = (iproc.eq.-1)
          ismine = (iproc.eq.myPE) .or. is_halobc
          if (.not.ismine) then

	      isvalid = (0 .le. iproc) .and. (iproc.le.numPEs-1)
	      call assert( isvalid,'** sendrecv_init: invalid iproc ',iproc)

              ncount(iproc) = ncount(iproc) + 1
          endif
	enddo
	enddo
	enddo

   ncount(myPE) = 0

   ntotal = 0
   do iproc=0,numPEs-1
     ntotal = ntotal + ncount(iproc) 
   enddo

   nrecv = count( ncount(:) .ne. 0)

   allocate( recvproc( max(1,nrecv) ) )

   nrecv = 0
   do iproc=0,numPEs-1
	if (ncount(iproc).ne.0) then
	   nrecv = nrecv + 1
	   recvproc(nrecv) = iproc
	endif
   enddo

   allocate( xrecv(nrecv+1) )
   allocate( recvijk(max(1,ntotal)) )

   xrecv(1) = 1
   do iter=1,nrecv
     iproc = recvproc(iter)
     xrecv(iter+1) = xrecv(iter) + ncount(iproc)
   enddo

   allocate( recvtag( max(1,nrecv) ) )

   do iter=1,nrecv
      iproc = recvproc(iter)
      src = iproc
      dest = myPE 
      recvtag(iter) = message_tag( src, dest )
   enddo
 

! ----------------------------
! second pass to fill in array
! ----------------------------
   if (lidebug.ge.1) then
	call write_debug( name, 'recv second pass ', myPE )
   endif

        ipos = 1

  do iter=1,nrecv
      jproc = recvproc(iter)

        do k=kstartx(myPE),kendx(myPE)
        do j=jstartx(myPE),jendx(myPE)
        do i=istartx(myPE),iendx(myPE)

          ii = imap(i)
          jj = jmap(j)
          kk = kmap(k)

          iproc = ijk2proc(ii,jj,kk)
          is_halobc = (iproc.eq.-1)
          ismine = (iproc.eq.myPE).or.(is_halobc)
          if ((.not.ismine) .and. (iproc.eq.jproc)) then


		ijk = funijk( i,j,k)
		recvijk( ipos ) = ijk
		ipos = ipos + 1
          endif
        enddo
        enddo
        enddo

   enddo

    if (ilayer.eq.1) then

        nsend1 = nsend
        xsend1 => xsend
	sendijk1 => sendijk
        sendproc1 => sendproc
	sendtag1 => sendtag

        nrecv1 = nrecv
        xrecv1 => xrecv
	recvijk1 => recvijk
        recvproc1 => recvproc
	recvtag1 => recvtag

    else

        nsend2 = nsend
        xsend2 => xsend
        sendijk2 => sendijk
        sendproc2 => sendproc
        sendtag2 => sendtag

        nrecv2 = nrecv
        xrecv2 => xrecv
        recvijk2 => recvijk
        recvproc2 => recvproc
        recvtag2 => recvtag

    endif


    nullify( xsend )
    nullify( sendijk )
    nullify( sendproc )
    nullify( sendtag )

    nullify( xrecv )
    nullify( recvijk )
    nullify( recvproc )
    nullify( recvtag )

enddo ! do ilayer

	
	deallocate( ncount )
	deallocate( ijk2proc )
	
	deallocate( istartx )
	deallocate( jstartx )
	deallocate( kstartx )
	deallocate( iendx ) 
	deallocate( jendx ) 
	deallocate( kendx ) 

	

    if (lidebug.ge.1) then

	call write_debug( name, ' allocate message buffers ' )
	call write_debug( name, 'nrecv1 ', nrecv1 )
	call write_debug( name, 'recvproc1 ', recvproc1 )
	call write_debug( name, 'recvtag1 ', recvtag1 )
	call write_debug( name, 'xrecv1 ', xrecv1 )


	lmax = size(recvijk1)
	allocate( line(lmax) )
        line(:) = " "

	ip = 1
	do ii=lbound(recvijk1,1),ubound(recvijk1,1)
 	   ijk = recvijk1(ii)
	   if(localfunc) then
	   call ijk_of(ijk,i,j,k)
	   else
	   i = i_of(ijk)
	   j = j_of(ijk)
	   k = k_of(ijk)
	   endif
	   write(line(ip),9001) ii,ijk, i,j,k
 9001      format('recvijk1( ', i6,') = ', i6, '( ', i6,',',i6,',',i6,') ')
           ip = ip + 1
	enddo
	call write_error( name, line, lmax )
	deallocate( line )

        lmax = size(recvijk2)
        allocate( line(lmax) )
        line(:) = " "

        ip = 1
        do ii=lbound(recvijk2,1),ubound(recvijk2,1)
           ijk = recvijk2(ii)
           if(localfunc) then
           call ijk_of(ijk,i,j,k)
           else
           i = i_of(ijk)
           j = j_of(ijk)
           k = k_of(ijk)
           endif

           write(line(ip),9101) ii,ijk, i,j,k
 9101      format('recvijk2( ', i6,') = ', i6, '( ', i6,',',i6,',',i6,') ')
           ip = ip + 1
        enddo
        call write_error( name, line, lmax )
        deallocate( line )


        call write_debug( name, ' allocate message buffers ' )
        call write_debug( name, 'nsend1 ', nsend1 )
        call write_debug( name, 'sendproc1 ', sendproc1 )
        call write_debug( name, 'sendtag1 ', sendtag1 )
        call write_debug( name, 'xsend1 ', xsend1 )



	lmax = size(sendijk1)
	allocate(line(lmax))
        line(:) = " "

	ip = 1
        do ii=lbound(sendijk1,1),ubound(sendijk1,1)
           ijk = sendijk1(ii)
           if(localfunc) then
           call ijk_of(ijk,i,j,k)
           else
           i = i_of(ijk)
           j = j_of(ijk)
           k = k_of(ijk)
           endif

           write(line(ip),9002) ii,ijk,   i,j,k
 9002      format('sendijk1( ', i6,') = ', i6, '( ', i6,',',i6,',',i6,') ')
	   ip = ip + 1
        enddo

	call write_error( name, line, lmax )
	deallocate( line )


        lmax = size(sendijk2)
        allocate(line(lmax))
        line(:) = " "

        ip = 1
        do ii=lbound(sendijk2,1),ubound(sendijk2,1)
           ijk = sendijk2(ii)
           if(localfunc) then
           call ijk_of(ijk,i,j,k)
           else
           i = i_of(ijk)
           j = j_of(ijk)
           k = k_of(ijk)
           endif

           write(line(ip),9102) ii,ijk,   i,j,k
 9102      format('sendijk2( ', i6,') = ', i6, '( ', i6,',',i6,',',i6,') ')
           ip = ip + 1
        enddo

        call write_error( name, line, lmax )
        deallocate( line )



    endif



! ------------------------
! allocate message buffers
! ------------------------
	isize = max(1, max(nsend1,nsend2))
	allocate( sendrequest( isize ) )

	isize = max(1, max(nrecv1,nrecv2))
	allocate( recvrequest( isize ) )
	



        if (lidebug.ge.1) then
		call write_debug(name, ' end of sendrecv_init ', myPE )
        endif
	
	end subroutine sendrecv_init



	subroutine sendrecv_begin_1d( X, ilayer, idebug )
        implicit none

	integer, intent(in),optional :: ilayer
	double precision, intent(inout), dimension(:) :: X
	integer, intent(in), optional :: idebug

!	interface
!
!	subroutine MPI_ISEND( buffer, count, datatype, dest, tag, &
!                        comm, request, ierror )
!	double precision buffer(*)
!	integer count,datatype,dest,tag,comm,request,ierror
!	end subroutine MPI_ISEND
!
!        subroutine MPI_IRECV( buffer, count, datatype, source, tag, &
!			comm, request, ierror )
!	double precision buffer(*)
!	integer count,datatype,source,tag,comm,request,ierror
!	end subroutine MPI_IRECV
!
!	end interface




!	---------------
!	local variables
!	---------------
	character(len=80), parameter :: name = 'sendrecv_begin_1d'

        integer :: lidebug

	integer ::  layer, datatype, comm, recvsize, sendsize, &
		ijk,jj,j1,j2, request, ii,count,source,dest, tag, ierror

	include 'function.inc'

        lidebug = 0
        if (present(idebug)) then
           lidebug = idebug
        endif

	layer = 1
	if (present(ilayer)) then
	   layer = ilayer
	endif

	if (layer.eq.1) then
	  nrecv = nrecv1
	  recvtag =>recvtag1
	  recvproc => recvproc1
	  recvijk => recvijk1
	  xrecv => xrecv1

	  nsend = nsend1
	  sendtag => sendtag1
	  sendproc => sendproc1
	  sendijk => sendijk1
	  xsend => xsend1
	else
	  nrecv = nrecv2
	  recvtag =>recvtag2
	  recvproc => recvproc2
	  recvijk => recvijk2
	  xrecv => xrecv2

	  nsend = nsend2
	  sendtag => sendtag2
	  sendproc => sendproc2
	  sendijk => sendijk2
	  xsend => xsend2
	endif
	  

!   --------------------------
!   post asynchronous receives
!   --------------------------

    if (lidebug.ge.1) then
	call write_debug(name, 'post asynchronous receives, nrecv = ', nrecv )
    endif

    if (nrecv.ge.1) then
	recvsize = xrecv( nrecv+1)-1

	allocate( drecvbuffer( recvsize ) )

        if (lidebug.ge.1) then
		call write_debug( name, 'recvsize, ubound(drecvbuffer,1) ', &
			recvsize, ubound(drecvbuffer,1) )

		call write_debug( name, 'ubound(xrecv,1) ', &
					ubound(xrecv,1) )
		call write_debug( name, 'ubound(recvproc,1) ', &
					ubound(recvproc,1) )
		call write_debug( name, 'ubound(recvtag,1) ', &
					ubound(recvtag,1) )

        endif

!	-------------
!	post receives
!	-------------
	datatype = MPI_DOUBLE_PRECISION
	comm = communicator

	do ii=1,nrecv
	   j1 = xrecv(ii)
	   j2 = xrecv(ii+1)-1
	   count = j2-j1+1
	   source = recvproc( ii )
	   tag = recvtag( ii )


           if (lidebug.ge.2) then

		call write_debug(name, 'mpi_irecv: ii,j1,j2 ', ii,j1,j2 )
		call write_debug(name, 'count, source, tag ', &
					count,source,tag )
           endif

	   call MPI_IRECV( drecvbuffer(j1), count, datatype, source, tag, &
			comm, request, ierror )

	   call MPI_Check( 'sendrecv_begin_1d:MPI_IRECV ', ierror )

	    recvrequest( ii ) = request
	enddo

   endif

!  -----------------------
!  post asynchronous sends
!  -----------------------

   if (lidebug.ge.1) then

	call write_debug(name, 'post asynchronous sends ')
   endif

   if (nsend.ge.1) then
        sendsize = xsend( nsend+1)-1

	allocate( dsendbuffer( sendsize ) )

        if (lidebug.ge.1) then

                call write_debug( name, 'sendsize, ubound(dsendbuffer,1) ', &
                        sendsize, ubound(dsendbuffer,1) )

                call write_debug( name, 'ubound(xsend,1) ', &
                                        ubound(xsend,1) )
                call write_debug( name, 'ubound(sendproc,1) ', &
                                        ubound(sendproc,1) )
                call write_debug( name, 'ubound(sendtag,1) ', &
                                        ubound(sendtag,1) )

        endif




!	-------------
!	perform sends
!	-------------
	datatype = MPI_DOUBLE_PRECISION
	comm = communicator

	do ii=1,nsend

!	----------------------------
!	perform copy into dsendbuffer
!	----------------------------

	    j1 = xsend(ii)
	    j2 = xsend(ii+1)-1
	    count = j2-j1+1

	    do jj=j1,j2
		ijk = sendijk( jj )
		dsendbuffer(jj) = X(ijk)
	    enddo

	    dest = sendproc( ii )
	    tag = sendtag( ii )

	    if (lidebug.ge.2) then

		call write_debug(name, 'mpi_isend: ii,j1,j2 ', ii,j1,j2)
		call write_debug(name, 'count, dest, tag ', count,dest,tag )
	    endif


	    call MPI_ISEND( dsendbuffer(j1), count, datatype, dest, tag, &
			comm, request, ierror )
	   call MPI_Check( 'sendrecv_begin_1d:MPI_ISEND ', ierror )

	   sendrequest( ii ) = request
	enddo

    endif

	return
	end subroutine sendrecv_begin_1d


	subroutine sendrecv_begin_1i( X, ilayer, idebug )
        implicit none

	integer, intent(in),optional :: ilayer
	integer, intent(inout), dimension(:) :: X
	integer, intent(in), optional :: idebug

!	interface
!
!	subroutine MPI_ISEND( buffer, count, datatype, dest, tag, &
!                        comm, request, ierror )
!	integer buffer(*)
!	integer count,datatype,dest,tag,comm,request,ierror
!	end subroutine MPI_ISEND
!
!        subroutine MPI_IRECV( buffer, count, datatype, source, tag, &
!			comm, request, ierror )
!	integer buffer(*)
!	integer count,datatype,source,tag,comm,request,ierror
!	end subroutine MPI_IRECV
!
!	end interface




!	---------------
!	local variables
!	---------------
	character(len=80), parameter :: name = 'sendrecv_begin_1i'

        integer :: lidebug

	integer ::  layer, datatype, comm, recvsize, sendsize, &
		ijk,jj,j1,j2, request, ii,count,source,dest, tag, ierror

	include 'function.inc'

        lidebug = 0
        if (present(idebug)) then
           lidebug = idebug
        endif

	layer = 1
	if (present(ilayer)) then
	   layer = ilayer
	endif

	if (layer.eq.1) then
	  nrecv = nrecv1
	  recvtag =>recvtag1
	  recvproc => recvproc1
	  recvijk => recvijk1
	  xrecv => xrecv1

	  nsend = nsend1
	  sendtag => sendtag1
	  sendproc => sendproc1
	  sendijk => sendijk1
	  xsend => xsend1
	else
	  nrecv = nrecv2
	  recvtag =>recvtag2
	  recvproc => recvproc2
	  recvijk => recvijk2
	  xrecv => xrecv2

	  nsend = nsend2
	  sendtag => sendtag2
	  sendproc => sendproc2
	  sendijk => sendijk2
	  xsend => xsend2
	endif
	  

!   --------------------------
!   post asynchronous receives
!   --------------------------

    if (lidebug.ge.1) then
	call write_debug(name, 'post asynchronous receives, nrecv = ', nrecv )
    endif

    if (nrecv.ge.1) then
	recvsize = xrecv( nrecv+1)-1
	allocate( irecvbuffer( recvsize ) )

        if (lidebug.ge.1) then
		call write_debug( name, 'recvsize, ubound(irecvbuffer,1) ', &
			recvsize, ubound(irecvbuffer,1) )

		call write_debug( name, 'ubound(xrecv,1) ', &
					ubound(xrecv,1) )
		call write_debug( name, 'ubound(recvproc,1) ', &
					ubound(recvproc,1) )
		call write_debug( name, 'ubound(recvtag,1) ', &
					ubound(recvtag,1) )

        endif

!	-------------
!	post receives
!	-------------
	datatype = MPI_INTEGER
	comm = communicator

	do ii=1,nrecv
	   j1 = xrecv(ii)
	   j2 = xrecv(ii+1)-1
	   count = j2-j1+1
	   source = recvproc( ii )
	   tag = recvtag( ii )


           if (lidebug.ge.2) then

		call write_debug(name, 'mpi_irecv: ii,j1,j2 ', ii,j1,j2 )
		call write_debug(name, 'count, source, tag ', &
					count,source,tag )
           endif

	   call MPI_IRECV( irecvbuffer(j1), count, datatype, source, tag, &
			comm, request, ierror )

	   call MPI_Check( 'sendrecv_begin_1i:MPI_IRECV ', ierror )

	    recvrequest( ii ) = request
	enddo

   endif

!  -----------------------
!  post asynchronous sends
!  -----------------------

   if (lidebug.ge.1) then

	call write_debug(name, 'post asynchronous sends ')
   endif

   if (nsend.ge.1) then
        sendsize = xsend( nsend+1)-1
	allocate( isendbuffer( sendsize ) )

        if (lidebug.ge.1) then

                call write_debug( name, 'sendsize, ubound(isendbuffer,1) ', &
                        sendsize, ubound(isendbuffer,1) )

                call write_debug( name, 'ubound(xsend,1) ', &
                                        ubound(xsend,1) )
                call write_debug( name, 'ubound(sendproc,1) ', &
                                        ubound(sendproc,1) )
                call write_debug( name, 'ubound(sendtag,1) ', &
                                        ubound(sendtag,1) )

        endif




!	-------------
!	perform sends
!	-------------
	datatype = MPI_INTEGER
	comm = communicator

	do ii=1,nsend

!	----------------------------
!	perform copy into sendbuffer
!	----------------------------

	    j1 = xsend(ii)
	    j2 = xsend(ii+1)-1
	    count = j2-j1+1

	    do jj=j1,j2
		ijk = sendijk( jj )
		isendbuffer(jj) = X(ijk)
	    enddo

	    dest = sendproc( ii )
	    tag = sendtag( ii )

	    if (lidebug.ge.2) then

		call write_debug(name, 'mpi_isend: ii,j1,j2 ', ii,j1,j2)
		call write_debug(name, 'count, dest, tag ', count,dest,tag )
	    endif


	    call MPI_ISEND( isendbuffer(j1), count, datatype, dest, tag, &
			comm, request, ierror )
	   call MPI_Check( 'sendrecv_begin_1i:MPI_ISEND ', ierror )

	   sendrequest( ii ) = request
	enddo

    endif

	return
	end subroutine sendrecv_begin_1i


	subroutine sendrecv_begin_1c( X, ilayer, idebug )
        implicit none

	integer, intent(in),optional :: ilayer
	character(len=*), intent(inout), dimension(:) :: X
	integer, intent(in), optional :: idebug

!	interface
!
!	subroutine MPI_ISEND( buffer, count, datatype, dest, tag, &
!                        comm, request, ierror )
!	character(len=*) buffer(*)
!	integer count,datatype,dest,tag,comm,request,ierror
!	end subroutine MPI_ISEND
!
!        subroutine MPI_IRECV( buffer, count, datatype, source, tag, &
!			comm, request, ierror )
!	character(len=*) buffer(*)
!	integer count,datatype,source,tag,comm,request,ierror
!	end subroutine MPI_IRECV
!
!	end interface




!	---------------
!	local variables
!	---------------
	character(len=80), parameter :: name = 'sendrecv_begin_1c'

        integer :: lidebug

	integer ::  layer, datatype, comm, recvsize, sendsize, &
		ijk,jj,j1,j2, request, ii,count,source,dest, tag, ierror

	integer :: ic, clen, jpos

	include 'function.inc'

        lidebug = 0
        if (present(idebug)) then
           lidebug = idebug
        endif

	layer = 1
	if (present(ilayer)) then
	   layer = ilayer
	endif

	jpos = lbound(X,1)
	clen = len( X( jpos ) )

	if (layer.eq.1) then
	  nrecv = nrecv1
	  recvtag =>recvtag1
	  recvproc => recvproc1
	  recvijk => recvijk1
	  xrecv => xrecv1

	  nsend = nsend1
	  sendtag => sendtag1
	  sendproc => sendproc1
	  sendijk => sendijk1
	  xsend => xsend1
	else
	  nrecv = nrecv2
	  recvtag =>recvtag2
	  recvproc => recvproc2
	  recvijk => recvijk2
	  xrecv => xrecv2

	  nsend = nsend2
	  sendtag => sendtag2
	  sendproc => sendproc2
	  sendijk => sendijk2
	  xsend => xsend2
	endif
	  

!   --------------------------
!   post asynchronous receives
!   --------------------------

    if (lidebug.ge.1) then
	call write_debug(name, 'post asynchronous receives, nrecv = ', nrecv )
    endif

    if (nrecv.ge.1) then
	recvsize = xrecv( nrecv+1)-1

	allocate( crecvbuffer( recvsize*clen ) )

        if (lidebug.ge.1) then
		call write_debug( name, 'recvsize, ubound(crecvbuffer,1) ', &
			recvsize, ubound(crecvbuffer,1) )

		call write_debug( name, 'ubound(xrecv,1) ', &
					ubound(xrecv,1) )
		call write_debug( name, 'ubound(recvproc,1) ', &
					ubound(recvproc,1) )
		call write_debug( name, 'ubound(recvtag,1) ', &
					ubound(recvtag,1) )

        endif

!	-------------
!	post receives
!	-------------
	datatype = MPI_CHARACTER
	comm = communicator

	do ii=1,nrecv
	   j1 = xrecv(ii)
	   j2 = xrecv(ii+1)-1

	   count = j2-j1+1
	   count = count*clen

	   source = recvproc( ii )
	   tag = recvtag( ii )


           if (lidebug.ge.2) then

		call write_debug(name, 'mpi_irecv: ii,j1,j2 ', ii,j1,j2 )
		call write_debug(name, 'count, source, tag ', &
					count,source,tag )
           endif

	   jpos = 1 + (j1-1)*clen
	   call MPI_IRECV( crecvbuffer(jpos), count, datatype, source, tag, &
			comm, request, ierror )

	   call MPI_Check( 'sendrecv_begin_1c:MPI_IRECV ', ierror )

	    recvrequest( ii ) = request
	enddo

   endif

!  -----------------------
!  post asynchronous sends
!  -----------------------

   if (lidebug.ge.1) then

	call write_debug(name, 'post asynchronous sends ')
   endif

   if (nsend.ge.1) then
        sendsize = xsend( nsend+1)-1

	allocate( csendbuffer( sendsize*clen ) )

        if (lidebug.ge.1) then

                call write_debug( name, 'sendsize, ubound(csendbuffer,1) ', &
                        sendsize, ubound(csendbuffer,1) )

                call write_debug( name, 'ubound(xsend,1) ', &
                                        ubound(xsend,1) )
                call write_debug( name, 'ubound(sendproc,1) ', &
                                        ubound(sendproc,1) )
                call write_debug( name, 'ubound(sendtag,1) ', &
                                        ubound(sendtag,1) )

        endif




!	-------------
!	perform sends
!	-------------
	datatype = MPI_CHARACTER
	comm = communicator

	do ii=1,nsend

!	----------------------------
!	perform copy into sendbuffer
!	----------------------------

	    j1 = xsend(ii)
	    j2 = xsend(ii+1)-1


	    count = j2-j1+1
	    count = count*clen

	    do jj=j1,j2
		ijk = sendijk( jj )
		do ic=1,clen
		    jpos = (jj-1)*clen + ic
		    csendbuffer(jpos) = X(ijk)(ic:ic)
		enddo
	    enddo

	    dest = sendproc( ii )
	    tag = sendtag( ii )

	    if (lidebug.ge.2) then

		call write_debug(name, 'mpi_isend: ii,j1,j2 ', ii,j1,j2)
		call write_debug(name, 'count, dest, tag ', count,dest,tag )
	    endif


	    jpos = (j1-1)*clen + 1
	    call MPI_ISEND( csendbuffer(jpos), count, datatype, dest, tag, &
			comm, request, ierror )
	   call MPI_Check( 'sendrecv_begin_1c:MPI_ISEND ', ierror )

	   sendrequest( ii ) = request
	enddo

    endif

	return
	end subroutine sendrecv_begin_1c


	subroutine sendrecv_end_1d( X, idebug )
	double precision, intent(inout), dimension(:) :: X
	integer, intent(in), optional :: idebug

	interface

	subroutine MPI_WAITANY(count, array_of_requests, jindex, &
				status, ierror)
	use mpi

	integer count
	integer array_of_requests(*)
	integer jindex
	integer status(MPI_STATUS_SIZE)
	integer ierror
	end subroutine MPI_WAITANY

	subroutine MPI_WAITALL( count, array_of_requests,  &
		array_of_status, ierror )
	use mpi

	integer count
	integer array_of_requests(*)
	integer array_of_status( MPI_STATUS_SIZE,*)
	integer ierror
	end subroutine MPI_WAITALL

	end interface

!	---------------
!	local variables
!	---------------
	character(len=80), parameter :: name = 'sendrecv_end_1d'

	logical, parameter :: use_waitany = .false.

	integer :: lidebug
	integer :: jj,ijk,  jindex, ii,j1,j2, ierror

	integer, dimension(MPI_STATUS_SIZE) :: recv_status
	integer, dimension(:,:), allocatable :: send_status

!	---------------
!	inline function
!	---------------
	include 'function.inc'
	
!	--------------------------
!	wait for sends to complete
!	--------------------------

	lidebug = 0
	if (present(idebug)) then
	   lidebug = idebug
	endif

	if (nsend .ge.1) then

	if (lidebug.ge.1) then

	   call write_debug(name, &
		'waiting for sends to complete, nsend  = ', nsend )
	endif

	allocate( send_status(MPI_STATUS_SIZE,nsend))
	call MPI_WAITALL( nsend, sendrequest, send_status, ierror )
	call MPI_Check( 'sendrecv_end_1d:MPI_WAITALL ', ierror )
	deallocate( send_status )

	deallocate( dsendbuffer )

	endif


!	--------------------------
!	wait for recvs to complete
!	--------------------------
	if (nrecv.ge.1) then

	if (lidebug.ge.1) then

	   call write_debug( name, &
		'waiting for receives to complete, nrecv =  ', nrecv )
	endif

      if (use_waitany) then
	do ii=1,nrecv
	   call MPI_WAITANY( nrecv, recvrequest, jindex, recv_status, ierror )
	   call MPI_Check( 'sendrecv_end_1d:MPI_WAITANY ', ierror )

	   j1 = xrecv( jindex )
	   j2 = xrecv( jindex + 1)-1

	   if (lidebug.ge.2) then
		call write_debug(name, 'jindex, j1,j2 ', jindex,j1,j2 )
	   endif

	   do jj=j1,j2
	     ijk = recvijk( jj )
             X(ijk) = drecvbuffer(jj)
	   enddo
	enddo
      else 
	   call MPI_WAITALL( nrecv, recvrequest, recv_status, ierror )
	   call MPI_Check( 'sendrecv_end_1d:MPI_WAITALL recv ', ierror )

	   j1 = xrecv(1)
	   j2 = xrecv( nrecv +1)-1
	   do jj=j1,j2
	       ijk = recvijk( jj )
	       X(ijk) = drecvbuffer(jj)
	   enddo
	endif

	deallocate( drecvbuffer )

	endif

	return
	end subroutine sendrecv_end_1d


	subroutine sendrecv_end_1c( X, idebug )
	character(len=*), intent(inout), dimension(:) :: X
	integer, intent(in), optional :: idebug

	interface

	subroutine MPI_WAITANY(count, array_of_requests, jindex, &
				status, ierror)
	use mpi

	integer count
	integer array_of_requests(*)
	integer jindex
	integer status(MPI_STATUS_SIZE)
	integer ierror
	end subroutine MPI_WAITANY

	subroutine MPI_WAITALL( count, array_of_requests,  &
		array_of_status, ierror )
	use mpi

	integer count
	integer array_of_requests(*)
	integer array_of_status( MPI_STATUS_SIZE,*)
	integer ierror
	end subroutine MPI_WAITALL

	end interface

!	---------------
!	local variables
!	---------------
	character(len=80), parameter :: name = 'sendrecv_end_1c'

	integer :: ic, clen, jpos

	logical, parameter :: use_waitany = .false.

	integer :: lidebug
	integer :: jj,ijk,  jindex, ii,j1,j2, ierror

	integer, dimension(MPI_STATUS_SIZE) :: recv_status
	integer, dimension(:,:), allocatable :: send_status

!	---------------
!	inline function
!	---------------
	include 'function.inc'
	
!	--------------------------
!	wait for sends to complete
!	--------------------------

	lidebug = 0
	if (present(idebug)) then
	   lidebug = idebug
	endif

	jpos = lbound(X,1)
	clen = len(X(jpos))

	if (nsend .ge.1) then

	if (lidebug.ge.1) then

	   call write_debug(name, &
		'waiting for sends to complete, nsend  = ', nsend )
	endif


	allocate( send_status(MPI_STATUS_SIZE,nsend))
	call MPI_WAITALL( nsend, sendrequest, send_status, ierror )
	call MPI_Check( 'sendrecv_end_1c:MPI_WAITALL ', ierror )
	deallocate( send_status )

	deallocate( csendbuffer )

	endif


!	--------------------------
!	wait for recvs to complete
!	--------------------------
	if (nrecv.ge.1) then

	if (lidebug.ge.1) then

	   call write_debug( name, &
		'waiting for receives to complete, nrecv =  ', nrecv )
	endif

      if (use_waitany) then
	do ii=1,nrecv
	   call MPI_WAITANY( nrecv, recvrequest, jindex, recv_status, ierror )
	   call MPI_Check( 'sendrecv_end_1c:MPI_WAITANY ', ierror )

	   j1 = xrecv( jindex )
	   j2 = xrecv( jindex + 1)-1

	   if (lidebug.ge.2) then
		call write_debug(name, 'jindex, j1,j2 ', jindex,j1,j2 )
	   endif

	   do jj=j1,j2
	     ijk = recvijk( jj )

               do ic=1,clen
                jpos = (jj-1)*clen + ic
                X(ijk)(ic:ic) = crecvbuffer(jpos)
               enddo


	   enddo
	enddo
      else 
	   call MPI_WAITALL( nrecv, recvrequest, recv_status, ierror )
	   call MPI_Check( 'sendrecv_end_1c:MPI_WAITALL recv ', ierror )

	   j1 = xrecv(1)
	   j2 = xrecv( nrecv +1)-1
	   do jj=j1,j2
	       ijk = recvijk( jj )

	       do ic=1,clen
		jpos = (jj-1)*clen + ic
	        X(ijk)(ic:ic) = crecvbuffer(jpos)
	       enddo



	   enddo
	endif

	deallocate( crecvbuffer )

	endif

	return
	end subroutine sendrecv_end_1c


	subroutine sendrecv_end_1i( X, idebug )
	integer, intent(inout), dimension(:) :: X
	integer, intent(in), optional :: idebug

	interface

	subroutine MPI_WAITANY(count, array_of_requests, jindex, &
				status, ierror)
	use mpi

	integer count
	integer array_of_requests(*)
	integer jindex
	integer status(MPI_STATUS_SIZE)
	integer ierror
	end subroutine MPI_WAITANY

	subroutine MPI_WAITALL( count, array_of_requests,  &
		array_of_status, ierror )
	use mpi

	integer count
	integer array_of_requests(*)
	integer array_of_status( MPI_STATUS_SIZE,*)
	integer ierror
	end subroutine MPI_WAITALL

	end interface

!	---------------
!	local variables
!	---------------
	character(len=80), parameter :: name = 'sendrecv_end_1i'

	logical, parameter :: use_waitany = .false.

	integer :: lidebug
	integer :: jj,ijk,  jindex, ii,j1,j2, ierror

	integer, dimension(MPI_STATUS_SIZE) :: recv_status
	integer, dimension(:,:), allocatable :: send_status

!	---------------
!	inline function
!	---------------
	include 'function.inc'
	
!	--------------------------
!	wait for sends to complete
!	--------------------------

	lidebug = 0
	if (present(idebug)) then
	   lidebug = idebug
	endif

	if (nsend .ge.1) then

	if (lidebug.ge.1) then

	   call write_debug(name, &
		'waiting for sends to complete, nsend  = ', nsend )
	endif

	allocate( send_status(MPI_STATUS_SIZE,nsend))
	call MPI_WAITALL( nsend, sendrequest, send_status, ierror )
	call MPI_Check( 'sendrecv_end_1i:MPI_WAITALL ', ierror )
	deallocate( send_status )

	deallocate( isendbuffer )

	endif


!	--------------------------
!	wait for recvs to complete
!	--------------------------
	if (nrecv.ge.1) then

	if (lidebug.ge.1) then

	   call write_debug( name, &
		'waiting for receives to complete, nrecv =  ', nrecv )
	endif

      if (use_waitany) then
	do ii=1,nrecv
	   call MPI_WAITANY( nrecv, recvrequest, jindex, recv_status, ierror )
	   call MPI_Check( 'sendrecv_end_1i:MPI_WAITANY ', ierror )

	   j1 = xrecv( jindex )
	   j2 = xrecv( jindex + 1)-1

	   if (lidebug.ge.2) then
		call write_debug(name, 'jindex, j1,j2 ', jindex,j1,j2 )
	   endif

	   do jj=j1,j2
	     ijk = recvijk( jj )
             X(ijk) = irecvbuffer(jj)
	   enddo
	enddo
      else 
	   call MPI_WAITALL( nrecv, recvrequest, recv_status, ierror )
	   call MPI_Check( 'sendrecv_end_1i:MPI_WAITALL recv ', ierror )

	   j1 = xrecv(1)
	   j2 = xrecv( nrecv +1)-1
	   do jj=j1,j2
	       ijk = recvijk( jj )
	       X(ijk) = irecvbuffer(jj)
	   enddo
	endif

	deallocate( irecvbuffer )

	endif

	return
	end subroutine sendrecv_end_1i


	subroutine send_recv_1c( X, ilayer, idebug )
	character(len=*),  dimension(:), intent(inout) :: X
	integer, intent(in), optional :: ilayer,idebug

	integer :: lidebug, layer

	lidebug = 0
	if (present(idebug)) then
	   lidebug = idebug
	endif

	layer = 1
	if (present(ilayer)) then
	   layer = ilayer
	endif

	call sendrecv_begin(X,layer,lidebug)
	call sendrecv_end( X, lidebug )

	return 
	end subroutine send_recv_1c

	subroutine send_recv_1d( X, ilayer, idebug )
	double precision,  dimension(:), intent(inout) :: X
	integer, intent(in), optional :: ilayer,idebug

	integer :: lidebug, layer

	lidebug = 0
	if (present(idebug)) then
	   lidebug = idebug
	endif

	layer = 1
	if (present(ilayer)) then
	   layer = ilayer
	endif

	call sendrecv_begin(X,layer,lidebug)
	call sendrecv_end( X, lidebug )

	return 
	end subroutine send_recv_1d


	subroutine send_recv_1i( X, ilayer, idebug )
	integer,  dimension(:), intent(inout) :: X
	integer, intent(in), optional :: ilayer,idebug

	integer :: lidebug, layer

	lidebug = 0
	if (present(idebug)) then
	   lidebug = idebug
	endif

	layer = 1
	if (present(ilayer)) then
	   layer = ilayer
	endif

	call sendrecv_begin(X,layer,lidebug)
	call sendrecv_end( X, lidebug )

	return 
	end subroutine send_recv_1i

	
	end module sendrecv
