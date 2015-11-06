      program test_sap

        use discretelement
        use pair_manager, only: pairs, init_pairs, add_pair,get_pair
        use sweep_and_prune

      implicit none

      integer :: nn,mm,box_A,box_B
      integer :: pair(2)
      type(aabb_t) aabb

      MAX_PIP = 10

      call init_sap(sap)

      aabb%minendpoint(1) = 5
      aabb%minendpoint(2) = 5
      aabb%minendpoint(3) = 5
      aabb%maxendpoint(1) = 15
      aabb%maxendpoint(2) = 15
      aabb%maxendpoint(3) = 15

      print *,"after init"
      !call print_boxes(sap)

      call add_box(sap,aabb,box_A)

      print *,"after add"
      call print_boxes(sap)

      aabb%minendpoint(1) = 10
      aabb%minendpoint(2) = 10
      aabb%minendpoint(3) = 10
      aabb%maxendpoint(1) = 17
      aabb%maxendpoint(2) = 17
      aabb%maxendpoint(3) = 17

      call add_box(sap,aabb,box_B)

      print *,"after second add:"
      call print_boxes(sap)

      call init_pairs

      print *,"PRESORT"
      call sort(sap)
      print *,"POSTSORT"

      print *,"after sort:"
      call print_boxes(sap)
      call check_boxes(sap)

      !call sweep(sap)
      !print *,"after sweep:"
      !call print_boxes(sap)
      !call check_boxes(sap)

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

      call init_pairs
      call sort(sap)
      !call sweep(sap)

      print *,"after second sweep"
      call print_boxes(sap)

      call get_pair(pair)
      if (pair(1).ne. 0 .or. pair(2).ne.0) then
         print *,"SHOULD NOT HAVE FOUND A PAIR HERE"
         stop __LINE__
      endif

      print *,"TEST SUCCESSFUL"

      end program
