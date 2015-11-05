      program test_sap

      use pair_manager, only: pairs, init_pair_iterator, add_pair,get_pair
      use sweep_and_prune

      implicit none

      integer :: nn,mm,box_A,box_B
      integer :: pair(2)
      type(sap_t) sap
      type(aabb_t) aabb

      call init_sap(sap)

      aabb%minendpoint(1) = 10
      aabb%minendpoint(2) = 10
      aabb%minendpoint(3) = 10
      aabb%maxendpoint(1) = 12.5
      aabb%maxendpoint(2) = 12.5
      aabb%maxendpoint(3) = 12.5

      print *,"after init"
      !call print_boxes(sap)

      call add_box(sap,aabb,box_A)

      print *,"after add"
      call print_boxes(sap)

      aabb%minendpoint(1) = 12
      aabb%minendpoint(2) = 12
      aabb%minendpoint(3) = 12
      aabb%maxendpoint(1) = 17
      aabb%maxendpoint(2) = 17
      aabb%maxendpoint(3) = 17

      call add_box(sap,aabb,box_B)

      print *,"after second add:"
      call print_boxes(sap)

      print *,"PRESORT"
      call sort(sap)
      print *,"POSTSORT"

      print *,"after sort:"
      call print_boxes(sap)
      call check_boxes(sap)

      call init_pair_iterator
      call sweep(sap)
      print *,"after sweep:"
      call print_boxes(sap)
      call check_boxes(sap)

      call get_pair(pair)
      if (pair(1).ne. 1 .or. pair(2).ne.2) then
         print *,"SHOULD HAVE FOUND A PAIR HERE"
         stop __LINE__
      endif

      aabb%minendpoint(1) = 100
      aabb%minendpoint(2) = 10
      aabb%minendpoint(3) = 10
      aabb%maxendpoint(1) = 125
      aabb%maxendpoint(2) = 12.5
      aabb%maxendpoint(3) = 12.5

      call update_box(sap,box_A,aabb)

      print *,"after update"
      call print_boxes(sap)

      call sort(sap)
      call init_pair_iterator
      call sweep(sap)

      print *,"after second sweep"
      call print_boxes(sap)

      call get_pair(pair)
      if (pair(1).ne. 0 .or. pair(2).ne.0) then
         print *,"SHOULD NOT HAVE FOUND A PAIR HERE"
         stop __LINE__
      endif

      print *,"TEST SUCCESSFUL"

      end program
