module sweep_and_prune

  ! endpoint data structure should be 64 bits in size
  type endpoint_t
     ! owner of the endpoint
     ! if Min endpoint, id = -box_id
     ! if Max endpoint, id = box_id
     integer(kind=4) :: box_id

     ! actual value of the endpoint; single precision to save space
     real :: value
  end type endpoint_t

  type aabb_t
     real :: maxendpoint(3)
     real :: minendpoint(3)
  end type aabb_t

  type box_t
     integer :: box_id
     integer :: maxendpoint_id(3)
     integer :: minendpoint_id(3)
  end type box_t

  type sap_t
     type(endpoint_t), dimension(:), allocatable :: x_endpoints
     type(endpoint_t), dimension(:), allocatable :: y_endpoints
     type(endpoint_t), dimension(:), allocatable :: z_endpoints
     integer :: x_endpoints_len
     integer :: y_endpoints_len
     integer :: z_endpoints_len
     type(box_t), dimension(:), allocatable :: boxes
     integer :: boxes_len
  end type sap_t

  ! the global SAP
  type(sap_t) sap

  public :: add_box, del_box, update_box, sort, sweep
  private :: partition

  contains

    subroutine check_boxes(this)
      implicit none
      type(sap_t), intent(inout) :: this
      integer :: ii

      do ii=1, this%boxes_len
         if (abs(this%x_endpoints(abs(this%boxes(ii)%minendpoint_id(1)))%box_id) .ne. this%boxes(ii)%box_id) then
            print *,"this%boxes(",ii,")%box_id = ",this%boxes(ii)%box_id
            print *,"this%boxes(",ii,")%minendpoint_id(1)",this%boxes(ii)%minendpoint_id(1)
            print *,"this%x_endpoints(this%boxes(",ii,")%minendpoint_id(1))%box_id",this%x_endpoints(abs(this%boxes(ii)%minendpoint_id(1)))%box_id
            stop __LINE__
         endif
         if (abs(this%y_endpoints(abs(this%boxes(ii)%minendpoint_id(2)))%box_id) .ne. this%boxes(ii)%box_id) stop __LINE__
         if (abs(this%z_endpoints(abs(this%boxes(ii)%minendpoint_id(3)))%box_id) .ne. this%boxes(ii)%box_id) stop __LINE__

         if (abs(this%x_endpoints(abs(this%boxes(ii)%maxendpoint_id(1)))%box_id) .ne. this%boxes(ii)%box_id) then
            print *,"this%boxes(",ii,")%box_id = ",this%boxes(ii)%box_id
            print *,"this%boxes(",ii,")%maxendpoint_id(1)",this%boxes(ii)%maxendpoint_id(1)
            print *,"this%x_endpoints(this%boxes(",ii,")%maxendpoint_id(1))%box_id",this%x_endpoints(abs(this%boxes(ii)%maxendpoint_id(1)))%box_id
            stop __LINE__
         endif
         if (abs(this%y_endpoints(abs(this%boxes(ii)%maxendpoint_id(2)))%box_id) .ne. this%boxes(ii)%box_id) stop __LINE__
         if (abs(this%z_endpoints(abs(this%boxes(ii)%maxendpoint_id(3)))%box_id) .ne. this%boxes(ii)%box_id) stop __LINE__
      enddo

    end subroutine check_boxes

    subroutine print_boxes(this)
      implicit none
      type(sap_t), intent(inout) :: this
      integer :: ii

      do ii=1, size(this%x_endpoints)
         print *,"ENDPOINTX ",this%x_endpoints_len,": ",ii,this%x_endpoints(ii)%box_id," has value ",this%x_endpoints(ii)%value
      enddo

      do ii=1, size(this%boxes)
         print *,"BOXX ",this%boxes_len,": ",this%boxes(ii)%box_id," exists from ",this%boxes(ii)%minendpoint_id(1)," to ",this%boxes(ii)%maxendpoint_id(1)
         print *,"BOXY ",this%boxes_len,": ",this%boxes(ii)%box_id," exists from ",this%boxes(ii)%minendpoint_id(2)," to ",this%boxes(ii)%maxendpoint_id(2)
         print *,"BOXZ ",this%boxes_len,": ",this%boxes(ii)%box_id," exists from ",this%boxes(ii)%minendpoint_id(3)," to ",this%boxes(ii)%maxendpoint_id(3)
      enddo

    end subroutine print_boxes

    subroutine init_sap(this)
      implicit none
      type(sap_t), intent(inout) :: this

      this%x_endpoints_len = 0
      this%y_endpoints_len = 0
      this%z_endpoints_len = 0

      allocate (this%x_endpoints(10))
      allocate (this%y_endpoints(10))
      allocate (this%z_endpoints(10))
      allocate (this%boxes(10))

    end subroutine init_sap

    SUBROUTINE boxes_GROW(box_array,new_size)
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: new_size
      type(box_t), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: box_array
      type(box_t), DIMENSION(:), ALLOCATABLE :: box_tmp
      INTEGER lSIZE

      lSIZE = size(box_array,1)
      allocate(box_tmp(new_size))
      box_tmp(1:lSIZE) = box_array(1:lSIZE)
      call move_alloc(box_tmp,box_array)

    END SUBROUTINE boxes_GROW

    SUBROUTINE endpoints_GROW(endpoint_array,new_size)
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: new_size
      type(endpoint_t), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: endpoint_array
      type(endpoint_t), DIMENSION(:), ALLOCATABLE :: endpoint_tmp
      INTEGER lSIZE

      lSIZE = size(endpoint_array,1)
      allocate(endpoint_tmp(new_size))
      endpoint_tmp(1:lSIZE) = endpoint_array(1:lSIZE)
      call move_alloc(endpoint_tmp,endpoint_array)

    END SUBROUTINE endpoints_GROW

    subroutine add_box(this,aabb,id)

      implicit none
      type(sap_t), intent(inout) :: this
      type(aabb_t), intent(in) :: aabb
      integer, intent(out) :: id

      if (this%x_endpoints_len+2 > size(this%x_endpoints)) then
         call endpoints_GROW(this%x_endpoints,2*(this%x_endpoints_len+2))
      endif

      if (this%y_endpoints_len+2 > size(this%y_endpoints)) then
         call endpoints_GROW(this%y_endpoints,2*(this%y_endpoints_len+2))
      endif

      if (this%z_endpoints_len+2 > size(this%z_endpoints)) then
         call endpoints_GROW(this%z_endpoints,2*(this%z_endpoints_len+2))
      endif

      if (this%boxes_len+1 > size(this%boxes)) then
         call boxes_GROW(this%boxes,2*(this%boxes_len+1))
      endif

      this%boxes_len = this%boxes_len + 1
      id = this%boxes_len

      this%x_endpoints_len = this%x_endpoints_len + 2
      this%y_endpoints_len = this%y_endpoints_len + 2
      this%z_endpoints_len = this%z_endpoints_len + 2

      this%boxes(id)%box_id = id
      this%boxes(id)%maxendpoint_id(1) = this%x_endpoints_len
      this%boxes(id)%maxendpoint_id(2) = this%y_endpoints_len
      this%boxes(id)%maxendpoint_id(3) = this%z_endpoints_len
      this%boxes(id)%minendpoint_id(1) = this%x_endpoints_len-1
      this%boxes(id)%minendpoint_id(2) = this%y_endpoints_len-1
      this%boxes(id)%minendpoint_id(3) = this%z_endpoints_len-1

      this%x_endpoints(this%x_endpoints_len-1)%box_id = -id
      this%x_endpoints(this%x_endpoints_len-1)%value = aabb%minendpoint(1)
      this%y_endpoints(this%y_endpoints_len-1)%box_id = -id
      this%y_endpoints(this%y_endpoints_len-1)%value = aabb%minendpoint(2)
      this%z_endpoints(this%z_endpoints_len-1)%box_id = -id
      this%z_endpoints(this%z_endpoints_len-1)%value = aabb%minendpoint(3)

      this%x_endpoints(this%x_endpoints_len)%box_id = id
      this%x_endpoints(this%x_endpoints_len)%value = aabb%maxendpoint(1)
      this%y_endpoints(this%y_endpoints_len)%box_id = id
      this%y_endpoints(this%y_endpoints_len)%value = aabb%maxendpoint(2)
      this%z_endpoints(this%z_endpoints_len)%box_id = id
      this%z_endpoints(this%z_endpoints_len)%value = aabb%maxendpoint(3)

    end subroutine add_box


    subroutine del_box(this,id)
      implicit none
      type(sap_t), intent(in) :: this
      integer, intent(in) :: id

      ! TODO

    end subroutine del_box

    subroutine update_box(this,id,aabb)
      implicit none
      type(sap_t), intent(inout) :: this
      integer, intent(in) :: id
      type(aabb_t), intent(in) :: aabb

      ! update end points
      this%x_endpoints(this%boxes(id)%maxendpoint_id(1))%value = aabb%maxendpoint(1)
      this%y_endpoints(this%boxes(id)%maxendpoint_id(2))%value = aabb%maxendpoint(2)
      this%z_endpoints(this%boxes(id)%maxendpoint_id(3))%value = aabb%maxendpoint(3)
      this%x_endpoints(this%boxes(id)%minendpoint_id(1))%value = aabb%minendpoint(1)
      this%y_endpoints(this%boxes(id)%minendpoint_id(2))%value = aabb%minendpoint(2)
      this%z_endpoints(this%boxes(id)%minendpoint_id(3))%value = aabb%minendpoint(3)

    end subroutine update_box

    subroutine sort(this)
      implicit none
      type(sap_t), intent(inout) :: this

      ! sort end points
      !print *,"SORTING FROM ",1, " TO ",this%x_endpoints_len

      !print *,"????? ",this%x_endpoints(this%x_endpoints_len),this%y_endpoints(this%x_endpoints_len),this%z_endpoints(this%x_endpoints_len)
      !print *,"????? ",this%x_endpoints(this%x_endpoints_len)%box_id,this%y_endpoints(this%x_endpoints_len)%box_id,this%z_endpoints(this%x_endpoints_len)%box_id
      !print *,"????? ",this%boxes(this%x_endpoints(this%x_endpoints_len)%box_id),this%boxes(this%y_endpoints(this%x_endpoints_len)%box_id),this%boxes(this%z_endpoints(this%x_endpoints_len)%box_id)

      call sort_endpoints(this%x_endpoints(1:this%x_endpoints_len),this,1)
      call sort_endpoints(this%y_endpoints(1:this%y_endpoints_len),this,2)
      call sort_endpoints(this%z_endpoints(1:this%z_endpoints_len),this,3)
    end subroutine sort

    subroutine sweep(this)

      use pair_manager, only: add_pair
      use resize, only: integer_grow
      use geometry

      implicit none
      type(sap_t), intent(inout) :: this

      ! active list of box id's
      integer :: ii,aa,minmax,ai
      logical :: debug

      type active_t
         integer, dimension(:), allocatable :: list
         integer :: list_len
      end type active_t

      type(active_t) :: active

      !stop __LINE__

      call active_init(active)

      do ii=1, this%x_endpoints_len
         minmax = this%x_endpoints(ii)%box_id

         if ( minmax < 0) then
            ! add pairs for new box x ( active pairs )

            debug = .false.
            if (this%boxes(abs(minmax))%box_id.eq.46) then
               ! print *,"DEBUG!  NOW CHECKING: 46 WHICH HAS RADIUS ",des_radius(46)," AND ACTIVE LENGTH ",active_get_length(active)
               ! print *,"MIN1 46 is ",sap%boxes(46)%minendpoint_id(1),sap%x_endpoints(sap%boxes(46)%minendpoint_id(1))%value
               ! print *,"MAX1 46 is ",sap%boxes(46)%maxendpoint_id(1),sap%x_endpoints(sap%boxes(46)%maxendpoint_id(1))%value
               ! print *,"MIN2 46 is ",sap%boxes(46)%minendpoint_id(2),sap%y_endpoints(sap%boxes(46)%minendpoint_id(2))%value
               ! print *,"MAX2 46 is ",sap%boxes(46)%maxendpoint_id(2),sap%y_endpoints(sap%boxes(46)%maxendpoint_id(2))%value
               debug = .true.
            endif

            if (this%boxes(abs(minmax))%box_id.eq.77) then
               ! print *,"DEBUG!  NOW CHECKING: 77 WHICH HAS RADIUS ",des_radius(77)," AND ACTIVE LENGTH ",active_get_length(active)
               ! print *,"MIN1 77 is ",sap%boxes(77)%minendpoint_id(1),sap%x_endpoints(sap%boxes(77)%minendpoint_id(1))%value
               ! print *,"MAX1 77 is ",sap%boxes(77)%maxendpoint_id(1),sap%x_endpoints(sap%boxes(77)%maxendpoint_id(1))%value
               ! print *,"MIN2 77 is ",sap%boxes(77)%minendpoint_id(2),sap%y_endpoints(sap%boxes(77)%minendpoint_id(2))%value
               ! print *,"MAX2 77 is ",sap%boxes(77)%maxendpoint_id(2),sap%y_endpoints(sap%boxes(77)%maxendpoint_id(2))%value
               debug = .true.
            endif

            do ai=1, active_get_length(active)
               aa = active_get(active,ai)
               ! compare other two axes

               !if (debug .and. aa.eq.77) then
                  if (this%boxes(abs(minmax))%box_id.eq.46) then
                     ! print *,"....NOW CHECKING: ",aa," WHICH HAS radius ",des_radius(aa)
                     ! print *,"MIN1 aa is ",sap%boxes(aa)%minendpoint_id(1),sap%x_endpoints(sap%boxes(aa)%minendpoint_id(1))%value
                     ! print *,"MAX1 aa is ",sap%boxes(aa)%maxendpoint_id(1),sap%x_endpoints(sap%boxes(aa)%maxendpoint_id(1))%value
                     ! print *,"MIN2 aa is ",sap%boxes(aa)%minendpoint_id(2),sap%y_endpoints(sap%boxes(aa)%minendpoint_id(2))%value
                     ! print *,"MAX2 aa is ",sap%boxes(aa)%maxendpoint_id(2),sap%y_endpoints(sap%boxes(aa)%maxendpoint_id(2))%value
                  endif
               !endif

               if (max(this%boxes(aa)%minendpoint_id(2),this%boxes(-minmax)%minendpoint_id(2)) <= min(this%boxes(-minmax)%maxendpoint_id(2),this%boxes(aa)%maxendpoint_id(2)) .and. &
                    (NO_K .or. max(this%boxes(aa)%minendpoint_id(3),this%boxes(-minmax)%minendpoint_id(3)) <= min(this%boxes(-minmax)%maxendpoint_id(3),this%boxes(aa)%maxendpoint_id(3)))) then
                  !if (46.eq.-minmax .and. aa.eq.77) print *,"FOUND PAIR! ADDING...",-minmax,aa
                  if (46.eq.-minmax) print *,"FOUND PAIR! ADDING...",-minmax,aa
                  call add_pair(aa,-minmax)
               else
                  if (46.eq.-minmax .and. aa.eq.77) then
                     print *,"COULDNT FIND...",-minmax,aa
                     print *,"NO_K =====  ",no_k
                     print *,"maxofmin2 ===",max(this%boxes(aa)%minendpoint_id(2),this%boxes(-minmax)%minendpoint_id(2)) <= min(this%boxes(-minmax)%maxendpoint_id(2),this%boxes(aa)%maxendpoint_id(2))
                     print *,"maxofmin3 ===",(NO_K .or. max(this%boxes(aa)%minendpoint_id(3),this%boxes(-minmax)%minendpoint_id(3)) <= min(this%boxes(-minmax)%maxendpoint_id(3),this%boxes(aa)%maxendpoint_id(3)))
                     print *,"this%boxes(aa)%minendpoint_id(2) = ",this%boxes(aa)%minendpoint_id(2),this%y_endpoints(this%boxes(aa)%minendpoint_id(2))%value
                     print *,"this%boxes(-minmax)%minendpoint_id(2) = ",this%boxes(-minmax)%minendpoint_id(2),this%y_endpoints(this%boxes(-minmax)%minendpoint_id(2))%value
                  print *,"max(this%boxes(aa)%minendpoint_id(2),this%boxes(-minmax)%minendpoint_id(2)) = ",max(this%boxes(aa)%minendpoint_id(2),this%boxes(-minmax)%minendpoint_id(2))
                  print *,"this%boxes(-minmax)%maxendpoint_id(2) = ",this%boxes(-minmax)%maxendpoint_id(2),this%y_endpoints(this%boxes(-minmax)%maxendpoint_id(2))%value
                  print *,"this%boxes(aa)%maxendpoint_id(2) = ",this%boxes(aa)%maxendpoint_id(2),this%y_endpoints(this%boxes(aa)%maxendpoint_id(2))%value
                  print *,"min(this%boxes(-minmax)%maxendpoint_id(2),this%boxes(aa)%maxendpoint_id(2)) = ",min(this%boxes(-minmax)%maxendpoint_id(2),this%boxes(aa)%maxendpoint_id(2))
                  print *,"this%boxes(aa)%minendpoint_id(3) = ",this%boxes(aa)%minendpoint_id(3)
                  print *,"this%boxes(-minmax)%minendpoint_id(3) = ",this%boxes(-minmax)%minendpoint_id(3)
                  print *,"max(this%boxes(aa)%minendpoint_id(3),this%boxes(-minmax)%minendpoint_id(3)) = ",max(this%boxes(aa)%minendpoint_id(3),this%boxes(-minmax)%minendpoint_id(3))
                  print *,"this%boxes(-minmax)%maxendpoint_id(3) = ",this%boxes(-minmax)%maxendpoint_id(3)
                  print *,"this%boxes(aa)%maxendpoint_id(3) = ",this%boxes(aa)%maxendpoint_id(3)
                  print *,"min(this%boxes(-minmax)%maxendpoint_id(3),this%boxes(aa)%maxendpoint_id(3)) = ",min(this%boxes(-minmax)%maxendpoint_id(3),this%boxes(aa)%maxendpoint_id(3))
               endif
               endif
            enddo

            ! add new box to active pair list
            call active_add(active,-minmax)

         else if ( 0 < minmax ) then
            !print *,"minmax = ",minmax
            ! remove box from active pair list
            call active_del(active,minmax)
         else
            print *,"minmax shouldn't be zero: ",minmax
            stop __LINE__
         endif
      enddo

    contains

      subroutine active_init(this)
        implicit none
        type(active_t), intent(inout) :: this

        allocate(this%list(10))
        this%list = 0
        this%list_len = 0

      end subroutine active_init

      subroutine active_add(this,new_active)
        implicit none
        type(active_t), intent(inout) :: this
        integer, intent(in) :: new_active
        integer :: old_len

        this%list_len = this%list_len + 1
        if (size(this%list) < this%list_len) then
           old_len = size(this%list)
           call integer_GROW(this%list,this%list_len)
           this%list(old_len+1:size(this%list)) = 0
        endif

        if (0.eq.this%list(this%list_len)) then
           this%list(this%list_len) = new_active
        else
           print *,"ACTIVE LIST SHOULD END IN ZERO"
           print *,"new_active = ",new_active
           do old_len = 1, size(this%list)
              print *,"list: ",old_len,this%list(old_len)
           enddo
           stop __LINE__
        endif
      end subroutine active_add

      subroutine active_del(this,active_to_delete)
        implicit none
        type(active_t), intent(inout) :: this
        integer, intent(in) :: active_to_delete
        integer :: ai, aa

        do ai=1, size(this%list)
           aa = this%list(ai)
           if (active_to_delete.eq.aa) then
              this%list(ai) = this%list(this%list_len)
              this%list(this%list_len) = 0
              this%list_len = this%list_len - 1
              exit
           endif
        enddo

      end subroutine active_del

      integer function active_get(this,index)
        implicit none
        type(active_t), intent(in) :: this
        integer :: index

        active_get = this%list(index)

      end function active_get

      integer function active_get_length(this)
        implicit none
        type(active_t), intent(in) :: this

        active_get_length = this%list_len

      end function active_get_length

    end subroutine sweep

     subroutine sort_endpoints(endpoints, sap, axis)
       use pair_manager, only: add_pair, del_pair
       !use discretelement
       implicit none
       type(sap_t), intent(inout) :: sap
       integer, intent(in) :: axis
       type(endpoint_t), dimension(:), intent(inout) :: endpoints
       type(endpoint_t) :: swap
       integer :: ii, jj
       integer :: tmp_ii, tmp_jj
       type(endpoint_t) :: sweeppoint, swappoint
       real :: sweepval
       !real :: zzz(3), asdf

! private void SortAxis(List<SweepPoint> axis)
!  {
!      for (int j = 1; j < axis.Count; j++)
!      {
!          SweepPoint keyelement = axis[j];
!          float key = keyelement.Value;
!          int i = j - 1;
!          while (i >= 0 && axis[i].Value > key)
!          {
!              SweepPoint swapper = axis[i];
!              if (keyelement.Begin && !swapper.Begin)
!              {
!                  if (CheckBoundingBoxes(swapper.Body, keyelement.Body))
!                  {
!                      lock (fullOverlaps) fullOverlaps.Add(new BroadphasePair(swapper.Body, keyelement.Body));
!                  }
!              }
!              if (!keyelement.Begin && swapper.Begin)
!              {
!                  lock (fullOverlaps) fullOverlaps.Remove(new BroadphasePair(swapper.Body, keyelement.Body));
!              }
!              axis[i + 1] = swapper;
!              i = i - 1;
!          }
!          axis[i + 1] = keyelement;
!      }
!  }

       do ii=2, size(endpoints)

          sweeppoint = endpoints(ii)
          sweepval = sweeppoint%value
          jj = ii-1

          do while ( 0 < jj )
             if ( endpoints(jj)%value <= sweepval ) exit

             swappoint = endpoints(jj)

            !print *,"partioning...",i,j
            !call print_boxes(sap)
            !call check_boxes(sap)

            !print *,"NOW COMPARING ENDPOINTS BELONGING TO BOXES:",sweeppoint%box_id,swappoint%box_id

            !print *,"ii,jj = ",ii,jj
            !print *,"endpoints(ii) = ",endpoints(ii)
            !print *,"endpoints(ii)%box_id = ",endpoints(ii)%box_id

            if (sweeppoint%box_id < 0 .and. 0 < swappoint%box_id ) then

               !print *,"MAYBE OVERLAPPPP"

               !zzz = des_pos_NEW(:,-sweeppoint%box_id)-des_pos_NEW(:,swappoint%box_id)
               !asdf = dot_product(zzz,zzz)
               !asdf = sqrt(asdf)
               !print *,"DIST==",ii,jj,asdf

               if (fullcheck(sap,-sweeppoint%box_id,swappoint%box_id,axis)) then
                  !print *,"PAIR FOUND"
                  call add_pair(-sweeppoint%box_id,swappoint%box_id)
               endif
            endif

            !zzz = des_pos_NEW(:,abs(sweeppoint%box_id))-des_pos_NEW(:,abs(swappoint%box_id))
            !asdf = dot_product(zzz,zzz)
            !asdf = sqrt(asdf)

            !if (asdf > 0.5 ) then
               !print *,"SHOULD REMOVE!!!!!!!",ii,jj,asdf
            !endif

            if (0 < sweeppoint%box_id .and. swappoint%box_id < 0 ) then
               !print *,"PAIR RMEOVEd"
               call del_pair(sweeppoint%box_id,-swappoint%box_id)
            endif

            !print *,"partioned! ",i,j
            !call print_boxes(sap)
            !call check_boxes(sap)

             !call swap_ep(endpoints(ii),endpoints(jj),axis)
            endpoints(jj+1) = swappoint

            if (swappoint%box_id < 0) then
               sap%boxes(-swappoint%box_id)%minendpoint_id(axis) = jj+1
            else
               sap%boxes(swappoint%box_id)%maxendpoint_id(axis) = jj+1
            endif

            jj = jj - 1
          enddo

          endpoints(jj+1) = sweeppoint
          if (sweeppoint%box_id < 0) then
             sap%boxes(-sweeppoint%box_id)%minendpoint_id(axis) = jj+1
          else
             sap%boxes(sweeppoint%box_id)%maxendpoint_id(axis) = jj+1
          endif

       enddo
     end subroutine sort_endpoints

     logical function fullcheck(this,id,id2,curr_axis)
       use geometry
       implicit none
       type(sap_t), intent(inout) :: this
       ! box ids to compare
       integer, intent(in) :: id, id2, curr_axis

        !print *,"FULLCHECK FOR ",id,id2

        !print *,"this%boxes(id)%minendpoint_id(1)==",this%boxes(id)%minendpoint_id(1),this%x_endpoints(this%boxes(id)%minendpoint_id(1))%value
        !print *,"this%boxes(id2)%minendpoint_id(1)==",this%boxes(id2)%minendpoint_id(1),this%x_endpoints(this%boxes(id2)%minendpoint_id(1))%value
        !print *,"this%boxes(id)%maxendpoint_id(1)==",this%boxes(id)%maxendpoint_id(1),this%x_endpoints(this%boxes(id)%maxendpoint_id(1))%value
        !print *,"this%boxes(id2)%maxendpoint_id(1)==",this%boxes(id2)%maxendpoint_id(1),this%x_endpoints(this%boxes(id2)%maxendpoint_id(1))%value
        !print *,""
        !print *,"this%boxes(id)%minendpoint_id(2)==",this%boxes(id)%minendpoint_id(2),this%y_endpoints(this%boxes(id)%minendpoint_id(2))%value
        !print *,"this%boxes(id2)%minendpoint_id(2)==",this%boxes(id2)%minendpoint_id(2),this%y_endpoints(this%boxes(id2)%minendpoint_id(2))%value
        !print *,"this%boxes(id)%maxendpoint_id(2)==",this%boxes(id)%maxendpoint_id(2),this%y_endpoints(this%boxes(id)%maxendpoint_id(2))%value
        !print *,"this%boxes(id2)%maxendpoint_id(2)==",this%boxes(id2)%maxendpoint_id(2),this%y_endpoints(this%boxes(id2)%maxendpoint_id(2))%value

       if (max(this%boxes(id)%minendpoint_id(1),this%boxes(id2)%minendpoint_id(1)) <= min(this%boxes(id2)%maxendpoint_id(1),this%boxes(id)%maxendpoint_id(1))) then
          !print *,"OVERLAP ON X AXIS"
          else
             !print *,"NOVERLAP ON X AXIS"
             endif
             if (max(this%boxes(id)%minendpoint_id(2),this%boxes(id2)%minendpoint_id(2)) <= min(this%boxes(id2)%maxendpoint_id(2),this%boxes(id)%maxendpoint_id(2))) then
                !print *,"OVERLAP ON Y AXIS"
             else
                !print *,"NOVERLAP ON Y AXIS"
             endif

             if (NO_K .or. max(this%boxes(id)%minendpoint_id(3),this%boxes(id2)%minendpoint_id(3)) <= min(this%boxes(id2)%maxendpoint_id(3),this%boxes(id)%maxendpoint_id(3))) then
                !print *,"OVERLAP ON Z AXIS"
    else
       !print *,"NOVERLAP ON Z AXIS"
    endif

    fullcheck = ((curr_axis.eq.1 .or. max(this%boxes(id)%minendpoint_id(1),this%boxes(id2)%minendpoint_id(1)) <= min(this%boxes(id2)%maxendpoint_id(1),this%boxes(id)%maxendpoint_id(1))) &
         .and. (curr_axis.eq.2 .or. max(this%boxes(id)%minendpoint_id(2),this%boxes(id2)%minendpoint_id(2)) <= min(this%boxes(id2)%maxendpoint_id(2),this%boxes(id)%maxendpoint_id(2))) &
         .and. (NO_K .or. curr_axis.eq.3 .or. max(this%boxes(id)%minendpoint_id(3),this%boxes(id2)%minendpoint_id(3)) <= min(this%boxes(id2)%maxendpoint_id(3),this%boxes(id)%maxendpoint_id(3))))

    !print *,"FROM FUCLLCHECK RETURNING:   ",fullcheck

     end function fullcheck

    recursive subroutine quicksort_endpoints(A, sap, axis)
      type(sap_t), intent(inout) :: sap
      integer, intent(in) :: axis
      type(endpoint_t), intent(in out), dimension(:) :: A
      integer :: iq

      !call check_boxes(sap)

      if(size(A) > 1) then
         call Partition(A, iq, sap, axis)
         call quicksort_endpoints(A(:iq-1),sap,axis)
         call quicksort_endpoints(A(iq:),sap,axis)
      endif
    end subroutine quicksort_endpoints

    subroutine Partition(A, marker, sap, axis)
      type(sap_t), intent(inout) :: sap
      integer, intent(in) :: axis
      type(endpoint_t), intent(in out), dimension(:) :: A
      integer, intent(out) :: marker
      type(endpoint_t) :: temp
      integer :: tmp_ii, tmp_jj
      integer :: i, j
      real :: x      ! pivot point
      x = A(1)%value
      i= 0
      j= size(A) + 1

      do
         j = j-1
         do
            if (A(j)%value <= x) exit
            j = j-1
         end do
         i = i+1
         do
            if (A(i)%value >= x) exit
            i = i+1
         end do
         if (i < j) then
            ! exchange A(i) and A(j)
            !print *,"partioning...",i,j
            !call print_boxes(sap)
            !call check_boxes(sap)

            !print *,"NOW SWAPPING ENDPOINTS BELONGING TO BOXES:",A(i)%box_id,A(j)%box_id

            !print *,"i = ",i
            !print *,"size(A) = ",size(A)
            !print *,"A(i) = ",A(i)
            !print *,"A(i)%box_id = ",A(i)%box_id

            if (A(i)%box_id < 0) then
               tmp_ii = sap%boxes(-A(i)%box_id)%minendpoint_id(axis)
            else
               tmp_ii = sap%boxes(A(i)%box_id)%maxendpoint_id(axis)
            endif
            if (A(j)%box_id < 0) then
               tmp_jj = sap%boxes(-A(j)%box_id)%minendpoint_id(axis)
            else
               tmp_jj = sap%boxes(A(j)%box_id)%maxendpoint_id(axis)
            endif

            temp = A(i)
            A(i) = A(j)
            !print *,"SET ENDPOINT TO VALUE ",A(i)%value," WHICH BELONG TO BOX ",abs(A(i)%box_id)
            A(j) = temp
            !print *,"SET ENDPOINT TO VALUE ",A(j)%value," WHICH BELONG TO BOX ",abs(A(j)%box_id)

            if (A(i)%box_id < 0) then
               !print *,j,"setting min endpoint of box ",-A(j)%box_id, " ON AXIS ",axis," TO ",i
               sap%boxes(-A(i)%box_id)%minendpoint_id(axis) = tmp_ii
            else
               !print *,j,"setting max endpoint of box ",A(j)%box_id, " ON AXIS ",axis," TO ",i
               sap%boxes(A(i)%box_id)%maxendpoint_id(axis) = tmp_ii
            endif
            if (A(j)%box_id < 0) then
               !print *,j,"setting min endpoint of box ",-A(j)%box_id, " ON AXIS ",axis," TO ",i
               sap%boxes(-A(j)%box_id)%minendpoint_id(axis) = tmp_jj
            else
               !print *,j,"setting max endpoint of box ",A(j)%box_id, " ON AXIS ",axis," TO ",i
               sap%boxes(A(j)%box_id)%maxendpoint_id(axis) = tmp_jj
            endif

            !print *,"partioned! ",i,j
            !call print_boxes(sap)
            !call check_boxes(sap)

         elseif (i == j) then
            marker = i+1
            return
         else
            marker = i
            return
         endif
      end do

    end subroutine Partition

end module sweep_and_prune
