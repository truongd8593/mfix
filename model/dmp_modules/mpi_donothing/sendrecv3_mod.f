        module sendrecv3
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


        integer,pointer, dimension(:) :: &
            send_persistent_request, recv_persistent_request,     &
            send_persistent_request1, send_persistent_request2,   &
            recv_persistent_request1, recv_persistent_request2

        integer :: nrecv1,nsend1, nrecv2,nsend2

!---------------
!EFD extra layer
!---------------
        integer, pointer, dimension(:) :: &
                recvproc3, recvtag3, xrecv3, recvijk3, &
                sendproc3, sendtag3, xsend3, sendijk3, &
                send_persistent_request3,recv_persistent_request3
        integer :: nrecv3,nsend3
        integer, parameter :: nlayers = 3






        logical,parameter :: localfunc=.false.

        logical,parameter :: use_persistent_message=.true.


        double precision, dimension(:), pointer :: &
                dsendbuffer, drecvbuffer
        integer, dimension(:), pointer :: &
                isendbuffer, irecvbuffer
        character, dimension(:), pointer :: &
                csendbuffer, crecvbuffer

        integer :: nrecv,nsend
        integer, pointer, dimension(:) :: &
                recvrequest, sendrequest, &
                xrecv,recvproc, recvijk, recvtag, &
                xsend,sendproc, sendijk, sendtag

        integer :: &
           kstart_all_myPE, jstart_all_myPE, istart_all_myPE, &
           kend_all_myPE, jend_all_myPE, iend_all_myPE

        integer :: communicator

!       -----------------
!       generic interface
!       -----------------
        interface sendrecv3_begin
        module procedure &
                sendrecv3_begin_1d, &
                sendrecv3_begin_1i, &
                sendrecv3_begin_1c
        end interface

        interface sendrecv3_end
        module procedure &
                sendrecv3_end_1d, &
                sendrecv3_end_1i, &
                sendrecv3_end_1c
        end interface

        interface send_recv3
        module procedure &
                send_recv3_1d, send_recv3_2d, send_recv3_3d, &
                send_recv3_1i, &
                send_recv3_1c
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


!---------------
!EFD extra layer
!---------------
        i1 = istart4_all(myPE)
        i2 = iend4_all(myPE)
        j1 = jstart4_all(myPE)
        j2 = jend4_all(myPE)
        k1 = kstart4_all(myPE)
        k2 = kend4_all(myPE)






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
!       ------------
!       double check
!       ------------
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


!---------------
!EFD extra layer
!---------------
        k1 = minval( kstart4_all(:) )
        k2 = maxval( kend4_all(:) )

        j1 = minval( jstart4_all(:) )
        j2 = maxval( jend4_all(:) )

        i1 = minval( istart4_all(:) )
        i2 = maxval( iend4_all(:) )



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
!       ------------
!       double check
!       ------------
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

        subroutine sendrecv3_init(        &
                comm,                    &
                cyclic_i,cyclic_j,cyclic_k, idebug )
        implicit none

        integer, intent(in) :: comm
        logical,intent(in) :: cyclic_i,cyclic_j,cyclic_k

        integer, intent(in), optional :: idebug

        return
        end subroutine sendrecv3_init



        subroutine sendrecv3_begin_1d( X, ilayer, idebug )
        implicit none

        integer, intent(in),optional :: ilayer
        double precision, intent(inout), dimension(:) :: X
        integer, intent(in), optional :: idebug

        return
        end subroutine sendrecv3_begin_1d


        subroutine sendrecv3_begin_1i( X, ilayer, idebug )
        implicit none

        integer, intent(in),optional :: ilayer
        integer, intent(inout), dimension(:) :: X
        integer, intent(in), optional :: idebug

        return
        end subroutine sendrecv3_begin_1i


        subroutine sendrecv3_begin_1c( X, ilayer, idebug )
        implicit none

        integer, intent(in),optional :: ilayer
        character(len=*), intent(inout), dimension(:) :: X
        integer, intent(in), optional :: idebug

        return
        end subroutine sendrecv3_begin_1c


        subroutine sendrecv3_end_1d( X, idebug )
        implicit none

        double precision, intent(inout), dimension(:) :: X
        integer, intent(in), optional :: idebug

        return
        end subroutine sendrecv3_end_1d


        subroutine sendrecv3_end_1c( X, idebug )
        implicit none

        character(len=*), intent(inout), dimension(:) :: X
        integer, intent(in), optional :: idebug

        return
        end subroutine sendrecv3_end_1c


        subroutine sendrecv3_end_1i( X, idebug )
        implicit none

        integer, intent(inout), dimension(:) :: X
        integer, intent(in), optional :: idebug

        return
        end subroutine sendrecv3_end_1i


        subroutine send_recv3_1c( X, ilayer, idebug )
        implicit none

        character(len=*),  dimension(:), intent(inout) :: X
        integer, intent(in), optional :: ilayer,idebug

        return
        end subroutine send_recv3_1c

        subroutine send_recv3_1d( X, ilayer, idebug )
        implicit none

        double precision,  dimension(:), intent(inout) :: X
        integer, intent(in), optional :: ilayer,idebug

        return
        end subroutine send_recv3_1d

        subroutine send_recv3_2d( X, ilayer, idebug )
        implicit none

        double precision,  dimension(:,:), intent(inout) :: X
        integer, intent(in), optional :: ilayer,idebug

        return
        end subroutine send_recv3_2d

        subroutine send_recv3_3d( X, ilayer, idebug )
        implicit none

        double precision,  dimension(:,:,:), intent(inout) :: X
        integer, intent(in), optional :: ilayer,idebug

        return
        end subroutine send_recv3_3d

        subroutine send_recv3_1i( X, ilayer, idebug )
        implicit none

        integer,  dimension(:), intent(inout) :: X
        integer, intent(in), optional :: ilayer,idebug

        return
        end subroutine send_recv3_1i


        end module sendrecv3
