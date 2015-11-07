module pair_manager

  use discretelement, only: MAX_PIP

  integer :: current_hash

  type pair_t
     integer(kind=4) :: ii
     integer(kind=4) :: jj
  end type pair_t

  type(pair_t), dimension(:), allocatable :: table

  ! if table_size/size(table) >> 50%, time to resize the hashtable
  integer :: table_size

contains

  subroutine init_pairs
    implicit none

    current_hash = 0
    if (.not.allocated(table)) allocate(table(0:2*MAX_PIP))
    table(:)%ii = 0
    table(:)%jj = 0
    table_size = 0

  end subroutine init_pairs

  subroutine get_pair(pair)
    implicit none

    integer, intent(out) :: pair(2)

    pair(1) = 0
    pair(2) = 0

    do while (current_hash < size(table))
       if (0.ne.table(current_hash)%ii .and. 0.ne.table(current_hash)%jj) then
          pair(1) = table(current_hash)%ii
          pair(2) = table(current_hash)%jj
          current_hash = current_hash + 1
          return
       endif
       current_hash = current_hash + 1
    enddo
  end subroutine get_pair

  logical function is_pair(ii,jj)
    implicit none

    integer, intent(in) :: ii,jj
    integer(kind=8) :: hash, init_hash

    if (ii < 1 .or. jj < 1) then
       print *,"invalid pair: ",ii,jj
       stop __LINE__
    endif

    ! assign ii to hash to convert to 64-bit
    hash = ii
    hash = mod(ishft(hash,32)+jj,size(table))
    init_hash = hash

    if (table(hash)%ii .eq. ii .and. table(hash)%jj .eq. jj) then
       is_pair = .true.
       return
    endif
    if (table(hash)%ii .eq. 0 .and. table(hash)%jj .eq. 0) then
       is_pair = .false.
       return
    endif
    hash = mod(hash+1,size(table))

    do while(hash .ne. init_hash)
       if (table(hash)%ii .eq. ii .and. table(hash)%jj .eq. jj) then
          is_pair = .true.
          return
       endif
       if (table(hash)%ii .eq. 0 .and. table(hash)%jj .eq. 0) then
          is_pair = .false.
          return
       endif
       hash = mod(hash+1,size(table))
    enddo

    print *,"loop in hash addressing, this should not occur"
    stop __LINE__

  end function is_pair

  recursive subroutine add_pair(i0,j0)
    use discretelement
    implicit none
    integer, intent(in) :: i0,j0
    integer :: ii, jj, nn, old_size
    integer(kind=8) :: hash, init_hash
    type(pair_t), dimension(:), allocatable :: table_tmp

    if (size(table) < 2*table_size ) then
       old_size = size(table)
       allocate(table_tmp(0:old_size))
       table_tmp(0:old_size-1) = table(0:old_size-1)
       deallocate(table)
       allocate(table(0:2*old_size))
       do nn=1, old_size
          if ( table_tmp(ii)%ii .ne. 0 .and. table_tmp(ii)%jj .ne. 0) then
             call add_pair(table_tmp(ii)%ii,table_tmp(ii)%jj)
          endif
       enddo
       deallocate(table_tmp)
    endif

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

    if (table(hash)%ii .eq. 0 .and. table(hash)%jj .eq. 0) then
       table(hash)%ii = ii
       table(hash)%ii = jj
       table_size = table_size + 1
       return
    endif
    hash = mod(hash+1,size(table))
    !print *,"HASH IS =",hash," TABLE IS ",table_size,"/",size(table)
    
    do while(hash .ne. init_hash)
       if (table(hash)%ii .eq. 0 .and. table(hash)%jj .eq. 0) then
          table(hash)%ii = ii
          table(hash)%ii = jj
          table_size = table_size + 1
          return
       endif
       hash = mod(hash+1,size(table))
       !print *,"HASH IS =",hash," TABLE IS ",table_size,"/",size(table)
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

    if (ii < 1 .or. jj < 1) then
       print *,"invalid pair: ",ii,jj
       stop __LINE__
    endif

    ! assign ii to hash to convert to 64-bit
    hash = ii
    hash = mod(ishft(hash,32)+jj,size(table))
    init_hash = hash

    if (table(hash)%ii .eq. 0 .and. table(hash)%jj .eq. 0) return 
    if (table(hash)%ii .eq. ii .and. table(hash)%jj .eq. jj) then
       ! 0,1 signifies DELETED hash entry
       table(hash)%ii = 0
       table(hash)%ii = 1
       table_size = table_size - 1
       return
    endif
    hash = mod(hash+1,size(table))

    do while(hash .ne. init_hash)
       if (table(hash)%ii .eq. 0 .and. table(hash)%jj .eq. 0) return
       if (table(hash)%ii .eq. ii .and. table(hash)%jj .eq. jj) then
          ! 0,1 signifies DELETED hash entry
          table(hash)%ii = 0
          table(hash)%ii = 1
          table_size = table_size - 1
          return
       endif
       hash = mod(hash+1,size(table))

    enddo

    print *,"loop in hash addressing, this should not occur.  maybe hash table is full"
    stop __LINE__

  end subroutine del_pair

end module pair_manager
