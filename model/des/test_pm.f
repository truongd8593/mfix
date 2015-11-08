      program test_pm

      use pair_manager

      implicit none

      integer :: nn,mm
      integer :: pair(2)

      call init_pairs

      do nn=1 , 10000
         call add_pair(1+abs(nn*nn),1+abs(nn*nn*nn))
      enddo

      do nn=1 , 10000
         if (.not. is_pair(1+abs(nn*nn),1+abs(nn*nn*nn))) stop __LINE__
      enddo

      call del_pair(1,1)
      call del_pair(4,16)
      call del_pair(9,1)
      call del_pair(10,10)
      call del_pair(10,100)

      do
         call get_pair(pair)
         print *,"pair:  ",pair(1),pair(2)
         if (pair(1).eq.0 .and. pair(2).eq.0) then
            exit
         endif
      enddo

      end program
