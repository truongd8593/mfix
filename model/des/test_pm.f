      program test_pm

      use pair_manager

      implicit none

      integer :: nn,mm
      integer :: pair(2)

      do nn=1 , 100
         do mm=1 , 20
      call add_pair(nn,nn*nn)
      enddo
      enddo

      call del_pair(1,1)
      call del_pair(4,16)
      call del_pair(9,1)
      call del_pair(10,10)
      call del_pair(10,100)

      call init_pair_iterator

      do
         call get_pair(pair)
         print *,"pair:  ",pair(1),pair(2)
         if (pair(1).eq.0 .and. pair(2).eq.0) then
            exit
         endif
      enddo

      end program
