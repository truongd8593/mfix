module pair_manager

  integer :: current_hash

  type pair_t
     integer(kind=4) :: ii
     integer(kind=4) :: jj
      integer(kind=1) :: count
  end type pair_t

  type(pair_t), dimension(:), allocatable :: table

  ! if table_size/size(table) >> 50%, time to resize the hashtable
  integer :: table_size

contains

  subroutine init_pairs
    implicit none

    current_hash = 0
    if (.not.allocated(table)) allocate(table(0:1000))
    table(:)%ii = 0
    table(:)%jj = 0
    table_size = 0

  end subroutine init_pairs

  logical function check_table()
    implicit none
    integer :: nn, blanks, deleted, full

    check_table = .true.
    return

    blanks = 0
    deleted = 0
    full = 0
    do nn=0, size(table)-1
       if (table(nn)%ii > 0 .and. table(nn)%jj > 0) then
          full = full + 1
       else if (table(nn)%ii .eq. 0 .and. table(nn)%jj .eq. 0) then
          blanks = blanks + 1
       else if (table(nn)%ii .eq. 0 .and. table(nn)%jj .eq. 1) then
          deleted = deleted + 1
       else
          print *,"SHOULD NEVER OCCUR"
          check_table = .false.
          return
       endif
    end do

    if (full .ne. table_size) then
       print *,"SIZE = ",size(table)
       print *,"blanks = ",blanks
       print *,"deleted = ",deleted
       print *,"full = ",full
       print *,"table_size = ",table_size
       check_table = .false.
    endif

    if (full+deleted+blanks .ne. size(table)) then
       print *,"SIZE = ",size(table)
       print *,"blanks = ",blanks
       print *,"deleted = ",deleted
       print *,"full = ",full
       print *,"table_size = ",table_size
       check_table = .false.
    endif

  end function check_table

  ! subroutine get_pair(pair)
  !   implicit none

  !   integer, intent(out) :: pair(2)

  !   pair(1) = 0
  !   pair(2) = 0

  !   do while (current_hash < size(table))
  !      if (0.ne.table(current_hash)%ii .and. 0.ne.table(current_hash)%jj) then
  !         pair(1) = table(current_hash)%ii
  !         pair(2) = table(current_hash)%jj
  !         current_hash = current_hash + 1
  !         return
  !      endif
  !      current_hash = current_hash + 1
  !   enddo
  ! end subroutine get_pair

  logical function is_pair(i0,j0)
    implicit none
    integer, intent(in) :: i0, j0
    integer :: ii, jj
    integer(kind=8) :: hash, init_hash

    ii = min(i0,j0)
    jj = max(i0,j0)

    if (ii < 1 .or. jj < 1) then
       print *,"invalid pair: ",ii,jj
       stop __LINE__
    endif

    ! assign ii to hash to convert to 64-bit
    hash = ii
    hash = mod(ishft(hash,32)+jj,size(table))
    init_hash = hash
    !print *,"INIT HASH IS =",hash," TABLE IS ",table_size,"/",size(table)

    do
       if (table(hash)%ii .eq. ii .and. table(hash)%jj .eq. jj) then
          is_pair = (table(hash)%count > 0)
          if (ii.eq.105 .and. 106.eq.jj) print *,"FOUND PAIR:",ii,jj,"   AT LOCATION:",hash,"   IN TABLE OF SIZE:  ",table_size,"/",size(table)
          return
       endif
       if (table(hash)%ii .eq. 0 .and. table(hash)%jj .eq. 0) then
          if (ii.eq.105 .and. 106.eq.jj) print *,"DID NOT FIND PAIR:",ii,jj,"   AT LOCATION:",hash,"   IN TABLE OF SIZE:  ",table_size,"/",size(table)
          is_pair = .false.
          return
       endif
       hash = mod(hash+1,size(table))
       !print *,"HASH IS =",hash," TABLE IS ",table_size,"/",size(table)
       if (hash .eq. init_hash) exit
    enddo

    print *,"loop in hash addressing, this should not occur"
    stop __LINE__

  end function is_pair

  recursive subroutine add_pair(i0,j0)
    use discretelement
    implicit none
    integer, intent(in) :: i0,j0
    integer :: ii, jj, nn, old_size, cc
    integer(kind=8) :: hash, init_hash
    type(pair_t), dimension(:), allocatable :: table_tmp

    if (i0 < 1 .or. j0 < 1) then
       print *,"invalid pair: ",i0,j0
       stop __LINE__
    endif

    if (size(table) < 2*table_size ) then
       old_size = size(table)
       allocate(table_tmp(0:old_size-1))
       if (size(table_tmp).ne.old_size) then
          !print *,"size = ",size(table_tmp)
          !print *,"old_size = ",old_size
          stop __LINE__
       endif
       table_tmp(0:old_size-1) = table(0:old_size-1)

       deallocate(table)
       allocate(table(0:2*old_size))
       table(:)%ii = 0
       table(:)%jj = 0
       do nn=0, old_size-1
          if ( table_tmp(nn)%ii .ne. 0 .and. table_tmp(nn)%jj .ne. 0) then
             do cc=1, table_tmp(nn)%count
                call add_pair(table_tmp(nn)%ii,table_tmp(nn)%jj)
             enddo
          endif
       enddo
       deallocate(table_tmp)
    endif

    ii = min(i0,j0)
    jj = max(i0,j0)

      if (ii.eq. 105 .and. jj.eq.106) then
         print *,"|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||   ADDING ",ii,jj
      endif

    if (ii < 1 .or. jj < 1) then
       print *,"invalid pair: ",ii,jj
       stop __LINE__
    endif

    ! assign ii to hash to convert to 64-bit
    hash = ii
    hash = mod(ishft(hash,32)+jj,size(table))
    init_hash = hash
    !print *,"INIT HASH IS =",hash," TABLE IS ",table_size,"/",size(table)

    do
       if (table(hash)%ii .eq. ii .and. table(hash)%jj .eq. jj) then
          table(hash)%count = table(hash)%count + 1
          if (ii.eq.105 .and. 106.eq.jj) print *,"INCREMENT PAIR:",ii,jj," COUNT ===",table(hash)%count,"   TO LOCATION:",hash,"   IN TABLE OF SIZE:  ",table_size,"/",size(table)
          return
       endif
       if (table(hash)%ii .eq. 0 .or. table(hash)%jj .eq. 0) then
          table(hash)%ii = ii
          table(hash)%jj = jj
          table(hash)%count = 1
          table_size = table_size + 1
          if (ii.eq.105 .and. 106.eq.jj) print *,"ADDED PAIR:",ii,jj,"   TO LOCATION:",hash,"   IN TABLE OF SIZE:  ",table_size,"/",size(table)
          if(.not. check_table()) stop __LINE__
          return
       endif
       hash = mod(hash+1,size(table))
       !!print *,"HASH IS =",hash," TABLE IS ",table_size,"/",size(table)
       if (hash .eq. init_hash) exit
    enddo

    print *,"loop in hash addressing, this should not occur.  maybe hash table is full"
    stop __LINE__

  end subroutine add_pair

  subroutine del_pair(i0,j0)
    use discretelement
    implicit none
    integer, intent(in) :: i0,j0
    integer :: ii, jj
    integer(kind=8) :: hash, init_hash

    ii = min(i0,j0)
    jj = max(i0,j0)

      if (ii.eq. 105 .and. jj.eq.106) then
         print *,"|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||   DELETING ",ii,jj
      endif

      if (ii < 1 .or. jj < 1) then
       print *,"invalid pair: ",ii,jj
       stop __LINE__
    endif

    ! assign ii to hash to convert to 64-bit
    hash = ii
    hash = mod(ishft(hash,32)+jj,size(table))
    init_hash = hash

    do
       if (table(hash)%ii .eq. 0 .and. table(hash)%jj .eq. 0) then
          if(.not. check_table()) stop __LINE__
          return
       endif
       if (table(hash)%ii .eq. ii .and. table(hash)%jj .eq. jj) then
          if (table(hash)%count > 1 ) then
             table(hash)%count = table(hash)%count - 1
             if (ii.eq.105 .and. 106.eq.jj) print *,"DECREMENT PAIR:",ii,jj," COUNT ===",table(hash)%count,"   TO LOCATION:",hash,"   IN TABLE OF SIZE:  ",table_size,"/",size(table)
          else
             ! 0,1 signifies DELETED hash entry
             table(hash)%ii = 0
             table(hash)%jj = 1
             table_size = table_size - 1
             if (ii.eq.105 .and. 106.eq.jj) print *,"REMOVED PAIR:",ii,jj,"   FROM LOCATION:",hash,"   IN TABLE OF SIZE:  ",table_size,"/",size(table)
          endif
          if(.not. check_table()) stop __LINE__
          return
       endif
       hash = mod(hash+1,size(table))
       if (hash .eq. init_hash) exit
    enddo

    if(.not. check_table()) stop __LINE__

    print *,"loop in hash addressing. must be a lot of DELETED entries:  ",table_size,"/",size(table)

  end subroutine del_pair

end module pair_manager
