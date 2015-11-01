module pair_manager

  !use discretelement, only :: MAX_PIP
  integer, parameter :: MAX_PIP = 1000

  integer, parameter :: MAX_NUM_NEIGH = 30
  integer, dimension(:,:), allocatable :: pairs

  integer :: current_row, current_column

contains

  subroutine init_pair_iterator
    implicit none

    current_row = 1
    current_column = 1

  end subroutine init_pair_iterator

  subroutine get_pair(pair)
    implicit none

    integer, intent(out) :: pair(2)

    pair(1) = 0
    pair(2) = 0

    do while (current_row <= size(pairs,1))

       do while (current_column <= MAX_NUM_NEIGH)
          if (pairs(current_row, current_column).ne.0) then
             pair(1) = current_row
             pair(2) = pairs(current_row, current_column)
             current_column = current_column + 1
             return
          endif
          current_column = current_column + 1
       enddo

       current_row = current_row + 1
       current_column = 1

    enddo

  end subroutine get_pair

  subroutine add_pair(ii,jj)
    implicit none
    integer, intent(in) :: ii,jj
    integer :: nn

    if (.not. allocated(pairs)) then
       allocate(pairs(MAX_PIP,MAX_NUM_NEIGH))
    else if(size(pairs,1) < MAX_PIP) then
       !integer_grow2_reverse(pairs,MAX_PIP)
    endif

    do nn=1, MAX_NUM_NEIGH
       if (pairs(ii,nn).eq.0) then
       pairs(ii,nn) = jj
       return
    endif
 end do

 print *,"particle ",ii," had more than ",MAX_NUM_NEIGH," neighbors."
 stop 1111

end subroutine add_pair

  subroutine del_pair(ii,jj)
    implicit none
    integer, intent(in) :: ii,jj
    integer :: nn

    if (.not. allocated(pairs)) then
       allocate(pairs(MAX_PIP,MAX_NUM_NEIGH))
    else if(size(pairs,1) < MAX_PIP) then
       !integer_grow2_reverse(pairs,MAX_PIP)
    endif

    do nn=1, MAX_NUM_NEIGH
       if (pairs(ii,nn).eq.jj) then
       pairs(ii,nn) = 0
       return
       endif
    end do

    print *,"pair ",ii,jj," was not in the pair manager"

  end subroutine del_pair

end module pair_manager
