module pair_manager

  use discretelement, only: MAX_PIP

  integer, parameter :: MAX_NUM_NEIGH = 30
  integer, dimension(:,:), allocatable :: pairs

  integer :: current_row, current_column

contains

  subroutine init_pairs
    implicit none

    current_row = 1
    current_column = 1
    if (.not.allocated(pairs)) allocate(pairs(MAX_PIP,MAX_NUM_NEIGH))
    pairs(:,:) = 0

  end subroutine init_pairs

  subroutine get_pair(pair)
    implicit none

    integer, intent(out) :: pair(2)

    pair(1) = 0
    pair(2) = 0

    do while (current_row <= size(pairs,1))

       do while (current_column <= MAX_NUM_NEIGH)
          !print *,"row,col == ",current_row, current_column
          if (pairs(current_row, current_column).ne.0) then
             pair(1) = current_row
             pair(2) = pairs(current_row, current_column)
             print *,"RETURNING PAIR: ",pair(1),pair(2)
             current_column = current_column + 1
             return
          endif
          current_column = current_column + 1
       enddo

       current_row = current_row + 1
       current_column = 1

    enddo

    print *,"RETURNING PAIR: ",pair(1),pair(2)

  end subroutine get_pair

  subroutine add_pair(ii,jj)
    use discretelement
    implicit none
    integer, intent(in) :: ii,jj
    integer :: nn, tmp, i0, j0

    i0 = min(ii,jj)
    j0 = max(ii,jj)

    if (ii.eq.jj) then
       print *,"tried to add pair with same index: ",ii
       stop __LINE__
    endif

    if (.not. allocated(pairs)) then
       allocate(pairs(MAX_PIP,MAX_NUM_NEIGH))
    else if(size(pairs,1) < MAX_PIP) then
       stop __LINE__
       !integer_grow2_reverse(pairs,MAX_PIP)
    endif

    do nn=1, MAX_NUM_NEIGH
       if (pairs(i0,nn).eq.0 .or. pairs(i0,nn).eq.j0) then
          pairs(i0,nn) = j0
          return
       endif
    end do

    print *,"particle ",i0," had more than ",MAX_NUM_NEIGH," neighbors."
    do nn=1, MAX_NUM_NEIGH
       print *,"PAIRS(",i0,",",nn,") = ",pairs(i0,nn)
    end do

    stop __LINE__

  end subroutine add_pair

  subroutine del_pair(ii,jj)
    implicit none
    integer, intent(in) :: ii,jj
    integer :: nn, tmp, i0, j0

    i0 = min(ii,jj)
    j0 = max(ii,jj)

    if (ii.eq.jj) then
       print *,"tried to add pair with same index: ",ii
       stop __LINE__
    endif

    if (.not. allocated(pairs)) then
       allocate(pairs(MAX_PIP,MAX_NUM_NEIGH))
    else if(size(pairs,1) < MAX_PIP) then
       stop __LINE__
       !integer_grow2_reverse(pairs,MAX_PIP)
    endif

    do nn=1, MAX_NUM_NEIGH
       if (pairs(i0,nn).eq.j0) then
       pairs(i0,nn) = 0
       return
       endif
    end do

    !print *,"pair ",i0,j0," was not in the pair manager"

  end subroutine del_pair

end module pair_manager
