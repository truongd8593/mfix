	module gridmap
	use parallel_mpi
	use geometry
	use sendrecv
	use compar
	use indices
	implicit none

	contains

!	Routine to partition the grid. It works for 1-d, 2-d decomposition 
!	in the current implementation

	subroutine partition

	integer, dimension(0:nodesi-1) :: isize1_all
	integer, dimension(0:nodesj-1) :: jsize1_all
	integer, dimension(0:nodesk-1) :: ksize1_all

	integer :: ip, iproc, isize, iremain, &
		   kp, kproc, ksize, kremain, ikproc, ierr

        if(numPEs.ne.(nodesi*nodesj*nodesk)) then
	write(*,*) 'incorrect distribution of processors'
	call MPI_abort( MPI_COMM_WORLD, ierr)
	endif

!	   Determine the size in i direction and add the remainder sequentially

	isize = (imax1-imin1+1)/nodesi
	isize1_all(0:nodesi-1) = isize

	iremain = (imax1-imin1+1) - nodesi*isize
	if (iremain.ge.1) then
	    isize1_all( 0:(iremain-1) ) = isize + 1
	endif

	jstart1_all(0:numPEs-1) = jmin1
	jend1_all(0:numPEs-1)   = jmax1

!          Determine the size in k direction and add the remainder sequentially

	ksize = (kmax1-kmin1+1)/nodesk
	ksize1_all(0:nodesk-1) = ksize

	kremain = (kmax1-kmin1+1) - nodesk*ksize
	if (kremain.ge.1) then
	    ksize1_all( 0:(kremain-1) ) = ksize + 1
	endif

!	The following is general for 1-d or 2-d decompostion
!       Determining  istart and kstart for all the processors

	ikproc = 0
	kp = kmin1
	do kproc=0,nodesk-1
	   ip = imin1
	   do iproc=0,nodesi-1

	      istart1_all(ikproc) = ip + sum(isize1_all(0:iproc-1))
	      iend1_all(ikproc) = istart1_all(ikproc) + isize1_all(iproc)-1

	      kstart1_all(ikproc) = kp + sum(ksize1_all(0:kproc-1))
	      kend1_all(ikproc) = kstart1_all(ikproc) + ksize1_all(kproc)-1

	      ikproc = ikproc+1

	   enddo
	enddo

	return
	end subroutine partition

	
! 	Initializing all the variables from the information obtained in the
!	above routine.

	subroutine gridmap_init
	
	integer :: iproc 
	integer :: ijk, ii, jj, kk, idebug, comm
	include 'function.inc'

!       ******************************************************************

        nodesj = 1

        if(no_k) then
        nodesi = numPEs
        nodesk = 1
        else
        nodesi = 1
        nodesk = numPEs
        endif

	allocate( ijksize3_all(0:numPEs-1) )

	allocate( ijkstart3_all(0:numPEs-1) )

	allocate( ijkend3_all(0:numPEs-1) )

	allocate( istart_all(0:numPEs-1) )
	allocate( jstart_all(0:numPEs-1) )
	allocate( kstart_all(0:numPEs-1) )

	allocate( istart1_all(0:numPEs-1) )
	allocate( jstart1_all(0:numPEs-1) )
	allocate( kstart1_all(0:numPEs-1) )

	allocate( istart2_all(0:numPEs-1) )
	allocate( jstart2_all(0:numPEs-1) )
	allocate( kstart2_all(0:numPEs-1) )

	allocate( istart3_all(0:numPEs-1) )
	allocate( jstart3_all(0:numPEs-1) )
	allocate( kstart3_all(0:numPEs-1) )

	allocate( iend_all(0:numPEs-1) )
	allocate( jend_all(0:numPEs-1) )
	allocate( kend_all(0:numPEs-1) )

	allocate( iend1_all(0:numPEs-1) )
	allocate( jend1_all(0:numPEs-1) )
	allocate( kend1_all(0:numPEs-1) )

	allocate( iend2_all(0:numPEs-1) )
	allocate( jend2_all(0:numPEs-1) )
	allocate( kend2_all(0:numPEs-1) )

	allocate( iend3_all(0:numPEs-1) )
	allocate( jend3_all(0:numPEs-1) )
	allocate( kend3_all(0:numPEs-1) )
	 
	allocate( displs(0:numPEs-1) )

	call partition

!	The upper and lower bounds are prescribed such that two ghost 
!	layers are allowed at the physical boundaries - this is consistent
!	with our present approach - need to be generalized if only one ghost 
!	layer is needed

	do iproc=0,numPEs-1

	   istart2_all(iproc) = max(imin1-1,min(imax1+1,istart1_all(iproc)-1))
	   istart3_all(iproc) = max(imin1-2,min(imax1+2,istart2_all(iproc)-1))

	   jstart2_all(iproc) = max(jmin1-1,min(jmax1+1,jstart1_all(iproc)-1))
	   jstart3_all(iproc) = max(jmin1-2,min(jmax1+2,jstart2_all(iproc)-1))

	   if(no_k) then
	   kstart2_all(iproc) = kstart1_all(iproc)
	   kstart3_all(iproc) = kstart1_all(iproc)

	   else

	   kstart2_all(iproc) = max(kmin1-1,min(kmax1+1,kstart1_all(iproc)-1))
	   kstart3_all(iproc) = max(kmin1-2,min(kmax1+2,kstart2_all(iproc)-1))

	   endif


	   iend2_all(iproc) = max(imin1-1,min(imax1+1,iend1_all(iproc)+1))
	   iend3_all(iproc) = max(imin1-2,min(imax1+2,iend2_all(iproc)+1))

	   jend2_all(iproc) = max(jmin1-1,min(jmax1+1,jend1_all(iproc)+1))
	   jend3_all(iproc) = max(jmin1-2,min(jmax1+2,jend2_all(iproc)+1))

           if(no_k) then
           kend2_all(iproc) = kend1_all(iproc)
	   kend3_all(iproc) = kend1_all(iproc)

	   else

	   kend2_all(iproc) = max(kmin1-1,min(kmax1+1,kend1_all(iproc)+1))
	   kend3_all(iproc) = max(kmin1-2,min(kmax1+2,kend2_all(iproc)+1))
	   endif

	enddo



	do iproc=0,numPEs-1

	   ijkstart3_all(iproc) = 1

 	   ijkend3_all(iproc) =  1 + (iend3_all(iproc) - istart3_all(iproc)) &
           + (jend3_all(iproc)-jstart3_all(iproc))*(iend3_all(iproc)-istart3_all(iproc)+1) &
           + (kend3_all(iproc)-kstart3_all(iproc))*(jend3_all(iproc)-jstart3_all(iproc)+1)* &
             (iend3_all(iproc)-istart3_all(iproc)+1)

	enddo

	do iproc=0,numPEs-1

	   ijksize3_all(iproc) = ijkend3_all(iproc) - ijkstart3_all(iproc) + 1

	enddo

	displs(0) = 0
	do iproc=1,numPEs-1
	displs(iproc) = displs(iproc-1)+ijksize3_all(iproc-1)
	enddo

	ijkstart3 = ijkstart3_all(myPE)

	ijkend3   = ijkend3_all(myPE)
	   
	ijksize3  = ijksize3_all(myPE)

	istart1   = istart1_all(myPE)
	iend1     = iend1_all(myPE)
	jstart1   = jstart1_all(myPE)
	jend1     = jend1_all(myPE)
	kstart1   = kstart1_all(myPE)
	kend1     = kend1_all(myPE)

	istart2   = istart2_all(myPE)
	iend2     = iend2_all(myPE)
	jstart2   = jstart2_all(myPE)
	jend2     = jend2_all(myPE)
	kstart2   = kstart2_all(myPE)
	kend2     = kend2_all(myPE)

	istart3   = istart3_all(myPE)
	iend3     = iend3_all(myPE)
	jstart3   = jstart3_all(myPE)
	jend3     = jend3_all(myPE)
	kstart3   = kstart3_all(myPE)
	kend3     = kend3_all(myPE)

!	Setup mapping to take care of cyclic boundary conditions

!       -----------------------------------
!       consider cyclic boundary condition
!       using the imap(:),jmap(:),kmap(:)
!       indirection arrays
!       -----------------------------------

        allocate( imap( imin3:imax3 ) )
        allocate( jmap( jmin3:jmax3 ) )
        allocate( kmap( kmin3:kmax3 ) )

        do kk=kmin3,kmax3
          kmap(kk) = kk
        enddo

        do jj=jmin3,jmax3
          jmap(jj) = jj
        enddo

        do ii=imin3,imax3
          imap(ii) = ii
        enddo

        if (cyclic_z) then
           kmap( kmax2 ) = kmin1
           kmap( kmin2 ) = kmax1
           if (kmax3.gt.kmax2) then
                kmap(kmax3) = kmap(kmax2)+1
           endif
           if (kmin3.lt.kmin2) then
                kmap(kmin3) = kmap(kmin2)-1
           endif
        endif

        if (cyclic_y) then
           jmap( jmax2 ) = jmin1
           jmap( jmin2 ) = jmax1
           if (jmax3.gt.jmax2) then
                jmap(jmax3) = jmap(jmax2)+1
           endif
           if (jmin3.lt.jmin2) then
                jmap(jmin3) = jmap(jmin2)-1
           endif
        endif

        if (cyclic_x) then
           imap( imax2 ) = imin1
           imap( imin2 ) = imax1
           if (imax3.gt.imax2) then
                imap(imax3) = imap(imax2)+1
           endif
           if (imin3.lt.imin2) then
                imap(imin3) = imap(imin2)-1
           endif
        endif

!	Defining new set of varaibles to define upper and lower bound of the indices to include 
!	actual physical boundaries of the problem

	do iproc = 0, numPEs-1

        istart = istart1
        iend = iend1
        jstart = jstart1
        jend = jend1
        kstart = kstart1
        kend = kend1
    
        if(istart3.eq.imin3) istart = istart2
        if(iend3.eq.imax3) iend = iend2
        if(jstart3.eq.jmin3) jstart = jstart2
        if(jend3.eq.jmax3) jend = jend2
        if(kstart3.eq.kmin3) kstart = kstart2
        if(kend3.eq.kmax3) kend = kend2

	istart_all(iproc) = istart
	iend_all(iproc)   = iend
	jstart_all(iproc) = jstart
	jend_all(iproc)   = jend
	kstart_all(iproc) = kstart
	kend_all(iproc)   = kend

	enddo


!	Call to sendrecv_init to set all the communication pattern

	comm = MPI_COMM_WORLD
        call sendrecv_init( comm, &
                cyclic_x, cyclic_y, cyclic_z, idebug=0 )
	

!//AIKEPARDBGSTOP 0922
!      write(*,"('(PE ',I2,'): from gridmap_init ',&
!                 /,9X,'Kmin2 = ',I6,'  Kmax3 = ',I6,'  Kmax = ',I6, &
!                 /,9X,'Jmin2 = ',I6,'  Jmax3 = ',I6,'  Jmax = ',I6,&
!		 /,9X,'Imin2 = ',I6,'  Imax3 = ',I6,'  Imax = ',I6)") &
!                 myPE,Kmin2,Kmax3,Kmax,Jmin2,Jmax3,Jmax,Imin2,Imax3,Imax

	return
	end subroutine gridmap_init

	end module gridmap
