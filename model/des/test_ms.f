      program test_pm

        use multi_sweep_and_prune
      use pair_manager

      implicit none

      integer :: nn,mm
      integer :: pair(2)
      type(multisap_t) :: ms
      type(aabb_t) :: aabb
      type(boxhandlelist_t) :: boxhandles(100)
      real :: mins(3), maxs(3)

      mins(1) = 0
      mins(2) = 0
      mins(3) = 0
      maxs(1) = 10
      maxs(2) = 10
      maxs(3) = 10
      call init_multisap(ms,2,2,2,mins,maxs)

      aabb%minendpoint(1) = 3
      aabb%minendpoint(2) = 3
      aabb%minendpoint(3) = 3
      aabb%maxendpoint(1) = 6
      aabb%maxendpoint(2) = 6
      aabb%maxendpoint(3) = 6
      call multisap_add(ms,aabb,10,boxhandles(10))

      aabb%minendpoint(1) = 4
      aabb%minendpoint(2) = 4
      aabb%minendpoint(3) = 4
      aabb%maxendpoint(1) = 7
      aabb%maxendpoint(2) = 7
      aabb%maxendpoint(3) = 7
      call multisap_add(ms,aabb,20,boxhandles(20))

      call multisap_quicksort(ms)
      call multisap_sweep(ms)
      call multisap_sort(ms)

      if (.not. is_pair(ms%hashtable,10,20)) stop __LINE__

      ! call del_pair(ms,1,1)
      ! call del_pair(ms,4,16)
      ! call del_pair(ms,9,1)
      ! call del_pair(ms,10,10)
      ! call del_pair(ms,10,100)

      ! do
      !    call get_pair(ms,pair)
      !    print *,"pair:  ",pair(1),pair(2)
      !    if (pair(1).eq.0 .and. pair(2).eq.0) then
      !       exit
      !    endif
      ! enddo

print *,"TEST SUCCESSFUL"

      end program
