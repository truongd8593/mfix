      program test_pm

      use pair_manager

      implicit none

      integer :: nn,mm
      integer :: pair(2)
      type(hashtable_t) :: tab

      call init_pairs(tab)

      ! print *,"abs(-82124) == ",abs(-892124)
      do nn=1 , 10000
         call add_pair(tab,1+nn*nn,100000+mod(1+nn*nn*nn,100000))
      enddo

      do nn=1 , 10000
         if (.not. is_pair(tab,1+nn*nn,100000+mod(1+nn*nn*nn,100000))) stop __LINE__
      enddo

      call del_pair(tab,1,1)
      call del_pair(tab,4,16)
      call del_pair(tab,9,1)
      call del_pair(tab,10,10)
      call del_pair(tab,10,100)

      do
         call get_pair(tab,pair)
         print *,"pair:  ",pair(1),pair(2)
         if (pair(1).eq.0 .and. pair(2).eq.0) then
            exit
         endif
      enddo

print *,"TEST SUCCESSFUL"

      end program
