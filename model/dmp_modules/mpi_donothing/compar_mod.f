	Module compar

	use mpi

!       Variables to be declared for parallel information. Added by
!       Ed and Sreekanth on 06/22/99. Removed from geometry_mod and put
!	in COMPAR module with some additional variables used by AEOLUS

!       myPE - my processor id (it varies from 0 to nproc-1)
!
!       numPEs - total number of nodes

        integer :: myPE, numPEs

!	mpierr - used by AEOLUS for error checking

	INTEGER :: mpierr

!//PAR_I/O 0811 specify the rank of the PE to be used for I/O
        INTEGER :: PE_IO = 0

!       nodesi, nodesj and nodesk represent the number of nodes
!       in i, j, k directions respectiviely.
!
!       nodesj = 1 (No decomposition along j-direction
!
!       For 1-D decomposition, nodesk = nproc for a 3d problem and
!       nodesi = nproc for a 2D problem.

        integer :: nodesi, nodesj, nodesk

!       root represents the 'root' processor. For now it is defaulted to
!       zero

        integer :: root

        data root /0/

!       istart1_all contains the starting i value for all the processors
!       excluding the ghost regions. istart2_all is for
!       one extra ghost layer and istart3_all is for two ghost layers.
!
!       Similarly iend1_all, iend2_all and iend3_all contain the ending
!       values. Similarly for j and k, jstart..., kstart.... are prescribed.
!
!       All the variables without the '_all' represent that processor values.
!
!       ijkstart3 denotes the starting value of ijk, which belongs to
!       the processor = funijk(istart3_all(myid),jstart3_all(myid),kstart3_all
!       (myid) for a 1-d decompostion of a 3D problem. For more details see
!       gridmap_mod.f90.
!
!       Similarly the end values are denoted by ijkend3_all
!
!       displs has the necessary shift information to position the buffer
!       in the scatterv and gatherv routines.
!
!       ijksize3 is the size of the element owned by each processor plus the ghost
!       regions.
!
!       '_all' has information about all the processor mapping based on above
!       convention
!
        integer, allocatable,dimension(:) ::  &
                ijkstart3_all,ijkend3_all,    &
                istart_all,istart1_all,istart2_all,istart3_all, &
                jstart_all,jstart1_all,jstart2_all,jstart3_all, &
                kstart_all,kstart1_all,kstart2_all,kstart3_all, &
                iend_all,iend1_all,iend2_all,iend3_all, &
                jend_all,jend1_all,jend2_all,jend3_all, &
                kend_all,kend1_all,kend2_all,kend3_all, &
                ijksize3_all, displs

!	Variables used for mapping i, j, k to ii, jj, kk to take care of cyclic conditions...

	integer, allocatable,dimension(:) :: imap, jmap, kmap
	integer, allocatable,dimension(:) :: imap_c, jmap_c, kmap_c

        integer :: &
                ijksize3, ijkstart3,ijkend3, istart3, iend3, jstart3, jend3, &
                kstart3, kend3, istart2, iend2, jstart2, jend2, &
                kstart2, kend2, istart1, iend1, jstart1, jend1, &
                kstart1, kend1

        integer :: istart, iend, jstart, jend, kstart, kend

!	Variables added at Aytekin's request on 09/01/99 for some code implementation

!// 500 0803 declaration for storing filebasename, e.g. mfix000.dat
    	CHARACTER(len=3) :: fbname
    	INTEGER :: idbg = 1

!       Funijk coefficients

        integer :: c0, c1, c2

!	contains

!       integer function funijk(i,j,k)
!       integer,  intent(in):: i,j,k

!       if(i.lt.istart3_all(myPE).or.i.gt.iend3_all(myPE).or. &
!        j.lt.jstart3_all(myPE).or.j.gt.jend3_all(myPE).or. &
!        k.lt.kstart3_all(myPE).or.k.gt.kend3_all(myPE)) then

!        write(*,"('(PE ',I2,'): i = ',I5,'  j = ',I5,'  k = ',I5&
!          &,' DOES NOT reside on this PE')") myPE,i,j,k
!// 375 1028 added a mechanism which assigns FUNIJK=-99999 if the current cell
!            is not residing on the associated processor's subdomain.
!            the return of this function calls needs to be checked.
!         call mfix_exit(myPE)
!        FUNIJK = -99999
	 
!        else

!       FUNIJK = 1 + (I - istart3_all(myPE)) &
!     + (J-jstart3_all(myPE))*(iend3_all(myPE)-istart3_all(myPE)+1) &
!     + (K-kstart3_all(myPE))*(jend3_all(myPE)-jstart3_all(myPE)+1)* &
!       (iend3_all(myPE)-istart3_all(myPE)+1)
!       endif
	
!       end function funijk


      END MODULE compar

