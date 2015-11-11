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
     integer :: maxendpoint_id(3)
     integer :: minendpoint_id(3)
     integer :: particle_id
  end type box_t

  type active_t
     integer, dimension(:), allocatable :: list
     integer :: list_len
  end type active_t

  type sap_t
     type(endpoint_t), dimension(:), allocatable :: x_endpoints
     type(endpoint_t), dimension(:), allocatable :: y_endpoints
     type(endpoint_t), dimension(:), allocatable :: z_endpoints
     integer :: x_endpoints_len
     integer :: y_endpoints_len
     integer :: z_endpoints_len
     integer :: id
     type(box_t), dimension(:), allocatable :: boxes
     integer :: boxes_len
     type(active_t) :: deleted_boxes
  end type sap_t

  type boxhandle_t
     integer :: sap_id
     integer :: box_id
  end type boxhandle_t

  integer, parameter :: MAX_SAPS = 8 ! maybe 4 for 2d, eventually
  type boxhandlelist_t
     type(boxhandle_t) :: list(MAX_SAPS)
  end type boxhandlelist_t

  type(boxhandlelist_t), DIMENSION(:),  ALLOCATABLE :: boxhandle         !(PARTICLES)

  public :: add_box, del_box, update_box, sort, sweep
  private :: partition

  contains

    logical function check_boxes(this)
      use discretelement
      use pair_manager
      use geometry
      implicit none
      type(sap_t), intent(inout) :: this
      integer :: ii, nn, ss, particle_id
      double precision :: diff,asdf

      check_boxes = .true.
      return

      do ii=1, this%boxes_len

         if (0 .ne.this%boxes(ii)%minendpoint_id(1)) cycle

         if (0.ne.this%boxes(ii)%minendpoint_id(1) .and. abs(this%x_endpoints(abs(this%boxes(ii)%minendpoint_id(1)))%box_id) .ne. ii) then
            print *,"SAP_ID=",this%id,"   this%boxes(",ii,")%minendpoint_id(1)",this%boxes(ii)%minendpoint_id(1)
            print *,"SAP_ID=",this%id,"this%x_endpoints(this%boxes(",ii,")%minendpoint_id(1))%box_id",this%x_endpoints(abs(this%boxes(ii)%minendpoint_id(1)))%box_id
            check_boxes = .false.
            return
         endif
         if (abs(this%y_endpoints(abs(this%boxes(ii)%minendpoint_id(2)))%box_id) .ne. ii) stop __LINE__
         if (abs(this%z_endpoints(abs(this%boxes(ii)%minendpoint_id(3)))%box_id) .ne. ii) stop __LINE__
         if (0.ne.this%boxes(ii)%maxendpoint_id(1) .and. abs(this%x_endpoints(abs(this%boxes(ii)%maxendpoint_id(1)))%box_id) .ne. ii) then
            print *,"SAP_ID=",this%id,"this%boxes(",ii,")%maxendpoint_id(1)",this%boxes(ii)%maxendpoint_id(1)
            print *,"SAP_ID=",this%id,"this%x_endpoints(this%boxes(",ii,")%maxendpoint_id(1))%box_id",this%x_endpoints(abs(this%boxes(ii)%maxendpoint_id(1)))%box_id
            check_boxes = .false.
            return
         endif
         if (0.ne.this%boxes(ii)%maxendpoint_id(1) .and. abs(this%y_endpoints(abs(this%boxes(ii)%maxendpoint_id(2)))%box_id) .ne. ii) stop __LINE__
         if (do_k) then
            if (abs(this%z_endpoints(abs(this%boxes(ii)%maxendpoint_id(3)))%box_id) .ne. ii) stop __LINE__
            endif

if (.true.) then
      particle_id = this%boxes(ii)%particle_id

      asdf = 0.1
      diff = (DES_POS_NEW(1,particle_id)-DES_RADIUS(particle_id) - this%x_endpoints(this%boxes(ii)%minendpoint_id(1))%value)
      if (asdf < abs(diff) .and. abs(diff) < 1000000000.0) then
         print *,"SAP_ID=",this%id,"  min xdiff ======== ",diff
         print *,"SAP_ID=",this%id,"DES_POS_NEW(1) = ",DES_POS_NEW(1,particle_id)
         print *,"SAP_ID=",this%id,"DES_RADIUS = ",DES_RADIUS(particle_id)
         print *,"SAP_ID=",this%id,"this%x_endpoints(this%boxes%minendpoint_id(1))%value = ",this%x_endpoints(this%boxes(ii)%minendpoint_id(1))%value
         check_boxes = .false.
         return
      endif

      diff = (DES_POS_NEW(1,particle_id)+DES_RADIUS(particle_id) - this%x_endpoints(this%boxes(ii)%maxendpoint_id(1))%value)
      if (asdf < abs(diff) .and. abs(diff) < 1000000000.0) then
         print *,"SAP_ID=",this%id,"SAP_ID::  ",this%id,"  max xdiff ======== ",diff
         print *,"SAP_ID=",this%id,"DES_POS_NEW(1,",particle_id,") = ",DES_POS_NEW(1,particle_id)
         print *,"SAP_ID=",this%id,"DES_RADIUS(",particle_id,") = ",DES_RADIUS(particle_id)
         print *,"SAP_ID=",this%id,"this%boxes%maxendpoint_id(1) = ",this%boxes(ii)%maxendpoint_id(1)
         print *,"SAP_ID=",this%id,"this%x_endpoints(this%boxes%maxendpoint_id(1))%value = ",this%x_endpoints(this%boxes(ii)%maxendpoint_id(1))%value
         print *,"SAP_ID=",this%id,"this%x_endpoints(",this%boxes(ii)%maxendpoint_id(1),")%value = ",this%x_endpoints(this%boxes(ii)%maxendpoint_id(1))%value
         check_boxes = .false.
         return
      endif

      diff = (DES_POS_NEW(2,particle_id)-DES_RADIUS(particle_id) - this%y_endpoints(this%boxes(ii)%minendpoint_id(2))%value)
      if (asdf < abs(diff) .and. abs(diff) < 1000000000.0) then
         print *,"SAP_ID=",this%id,"  min ydiff ======== ",diff
         print *,"SAP_ID=",this%id,"DES_POS_NEW(2,particle_id) = ",DES_POS_NEW(2,particle_id)
         print *,"SAP_ID=",this%id,"DES_RADIUS(particle_id) = ",DES_RADIUS(particle_id)
         print *,"SAP_ID=",this%id,"this%y_endpoints(this%boxes%minendpoint_id(2))%value = ",this%y_endpoints(this%boxes(ii)%minendpoint_id(2))%value
         check_boxes = .false.
         return
      endif

      diff = (DES_POS_NEW(2,particle_id)+DES_RADIUS(particle_id) - this%y_endpoints(this%boxes(ii)%maxendpoint_id(2))%value)
      if (asdf < abs(diff) .and. abs(diff) < 1000000000.0) then
         print *,"SAP_ID=",this%id,"  max ydiff ======== ",diff
         print *,"SAP_ID=",this%id,"DES_POS_NEW(2,particle_id) = ",DES_POS_NEW(2,particle_id)
         print *,"SAP_ID=",this%id,"DES_RADIUS(particle_id) = ",DES_RADIUS(particle_id)
         print *,"SAP_ID=",this%id,"this%y_endpoints(this%boxes%maxendpoint_id(2))%value = ",this%y_endpoints(this%boxes(ii)%maxendpoint_id(2))%value
         check_boxes = .false.
         return
      endif
endif
   enddo

      check_boxes = .true.
      return

    end function check_boxes

      logical function check_sort(this)
        use geometry, only: do_k
        implicit none
        type(sap_t), intent(inout) :: this
        integer :: ii, nn

         ! CHECK SORT
        do ii=2, this%x_endpoints_len
            if (this%x_endpoints(ii)%value < this%x_endpoints(ii-1)%value) then
               print *,"SAP_ID=",this%id,"*********x endpoints********************************************************************"
               print *,"SAP_ID=",this%id,"ii-1:",ii-1,"  endpoints(ii):",this%x_endpoints(ii-1)%box_id,this%x_endpoints(ii-1)%value
               print *,"SAP_ID=",this%id,"ii:",ii,"  endpoints(ii):",this%x_endpoints(ii)%box_id,this%x_endpoints(ii)%value
               print *,"SAP_ID=",this%id,"****************************************************************************************"
               check_sort = .false.
               return
            endif
         enddo

         ! CHECK SORT
         do ii=2, this%y_endpoints_len
            if (this%y_endpoints(ii)%value < this%y_endpoints(ii-1)%value) then
               print *,"SAP_ID=",this%id,"****************y endpoints************************************************************************"
               print *,"SAP_ID=",this%id,"ii-1:",ii-1,"  endpoints(ii):",this%y_endpoints(ii-1)%box_id,this%y_endpoints(ii-1)%value
               print *,"SAP_ID=",this%id,"ii:",ii,"  endpoints(ii):",this%y_endpoints(ii)%box_id,this%y_endpoints(ii)%value
               print *,"SAP_ID=",this%id,"****************************************************************************************"
               check_sort = .false.
               return
            endif
         enddo

         ! CHECK SORT
         if (do_k) then
         do ii=2, this%z_endpoints_len
            if (this%z_endpoints(ii)%value < this%z_endpoints(ii-1)%value) then
               print *,"SAP_ID=",this%id,"***********z endpoints******************************************************************"
               print *,"SAP_ID=",this%id,"ii-1:",ii-1,"  endpoints(ii):",this%z_endpoints(ii-1)%box_id,this%z_endpoints(ii-1)%value
               print *,"SAP_ID=",this%id,"ii:",ii,"  endpoints(ii):",this%z_endpoints(ii)%box_id,this%z_endpoints(ii)%value
               print *,"SAP_ID=",this%id,"****************************************************************************************"
               check_sort = .false.
               return
            endif
         enddo
         endif

      check_sort = .true.
      return

    end function check_sort

    subroutine print_boxes(this)
      implicit none
      type(sap_t), intent(inout) :: this
      integer :: ii

      ! CHECK SORT
      ! do ii=2, this%x_endpoints_len
      !    if (this%x_endpoints(ii)%value < this%x_endpoints(ii-1)%value) then
      !       print *,"SAP_ID=",this%id,"****************************************************************************************"
      !       print *,"SAP_ID=",this%id,"ii-1:",ii-1,"  endpoints(ii):",this%x_endpoints(ii-1)%box_id,this%x_endpoints(ii-1)%value
      !       print *,"SAP_ID=",this%id,"ii:",ii,"  endpoints(ii):",this%x_endpoints(ii)%box_id,this%x_endpoints(ii)%value
      !       print *,"SAP_ID=",this%id,"****************************************************************************************"
      !       !stop __LINE__
      !    endif
      ! enddo

      do ii=1, this%x_endpoints_len
         print *,"SAP_ID=",this%id,"ENDPOINTX ",this%x_endpoints_len,": ",ii,this%x_endpoints(ii)%box_id," has value ",this%x_endpoints(ii)%value
      enddo

      do ii=1, this%boxes_len 
         print *,"SAP_ID=",this%id,"BOXX ",this%boxes_len,": ",ii," exists from ",this%boxes(ii)%minendpoint_id(1)," to ",this%boxes(ii)%maxendpoint_id(1)
         print *,"SAP_ID=",this%id,"BOXY ",this%boxes_len,": ",ii," exists from ",this%boxes(ii)%minendpoint_id(2)," to ",this%boxes(ii)%maxendpoint_id(2)
         print *,"SAP_ID=",this%id,"BOXZ ",this%boxes_len,": ",ii," exists from ",this%boxes(ii)%minendpoint_id(3)," to ",this%boxes(ii)%maxendpoint_id(3)
      enddo

    end subroutine print_boxes

    subroutine init_sap(this)
      implicit none
      type(sap_t), intent(inout) :: this

      this%x_endpoints_len = 0
      this%y_endpoints_len = 0
      this%z_endpoints_len = 0
      this%boxes_len = 0

      allocate (this%x_endpoints(10))
      allocate (this%y_endpoints(10))
      allocate (this%z_endpoints(10))
      allocate (this%boxes(10))

      call active_init(this%deleted_boxes)


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

    subroutine add_box(this,aabb,particle_id,id)
      use discretelement

      implicit none
      type(sap_t), intent(inout) :: this
      type(aabb_t), intent(in) :: aabb
      integer, intent(in) :: particle_id
      integer, intent(out) :: id

      ! print *,"SAP_ID=",this%id,"SAP:: ",this%id,"ADDING PARTICLE #   ",particle_id

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

      id = active_pop(this%deleted_boxes)

      if ( id < 0 ) then
         this%boxes_len = this%boxes_len + 1
         id = this%boxes_len
      endif

      this%x_endpoints_len = this%x_endpoints_len + 2
      this%y_endpoints_len = this%y_endpoints_len + 2
      this%z_endpoints_len = this%z_endpoints_len + 2

      this%boxes(id)%particle_id = particle_id
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

      ! print *,"SAP_ID=",this%id,"JUST ADDED PARTICLE",particle_id," TO SAP ",this%id
      ! print *,"SAP_ID=",this%id," PARTICLE IS AT ",des_pos_new(:,particle_id)
      ! print *,"SAP_ID=",this%id,"XXX === ",aabb%minendpoint(1)," TO ",aabb%maxendpoint(1)
      ! print *,"SAP_ID=",this%id,"YYY === ",aabb%minendpoint(2)," TO ",aabb%maxendpoint(2)
      ! print *,"SAP_ID=",this%id,"ZZZ === ",aabb%minendpoint(3)," TO ",aabb%maxendpoint(3)
      ! if (.not. check_boxes(this)) stop __LINE__

    end subroutine add_box


    subroutine del_box(this,id)
      implicit none
      type(sap_t), intent(inout) :: this
      integer, intent(in) :: id
      integer :: len

      this%x_endpoints(this%boxes(id)%maxendpoint_id(1))%value = HUGE(0.0)
      this%y_endpoints(this%boxes(id)%maxendpoint_id(2))%value = HUGE(0.0)
      this%z_endpoints(this%boxes(id)%maxendpoint_id(3))%value = HUGE(0.0)
      this%x_endpoints(this%boxes(id)%minendpoint_id(1))%value = HUGE(0.0)
      this%y_endpoints(this%boxes(id)%minendpoint_id(2))%value = HUGE(0.0)
      this%z_endpoints(this%boxes(id)%minendpoint_id(3))%value = HUGE(0.0)

      !call sort(this)

      ! this%x_endpoints_len = this%x_endpoints_len - 2
      ! this%y_endpoints_len = this%y_endpoints_len - 2
      ! this%z_endpoints_len = this%z_endpoints_len - 2

      ! this%x_endpoints(this%boxes(len)%maxendpoint_id(1))%box_id = id
      ! this%y_endpoints(this%boxes(len)%maxendpoint_id(2))%box_id = id
      ! this%z_endpoints(this%boxes(len)%maxendpoint_id(3))%box_id = id
      ! this%x_endpoints(this%boxes(len)%minendpoint_id(1))%box_id = -id
      ! this%y_endpoints(this%boxes(len)%minendpoint_id(2))%box_id = -id
      ! this%z_endpoints(this%boxes(len)%minendpoint_id(3))%box_id = -id

       ! if (abs(this%boxes(len)%maxendpoint_id(1)).eq.289) then
       !    print *,"SAP_ID=",this%id,"289289289289289289289289    this%boxes(",id,")%particle_id = ",this%boxes(id)%particle_id
       ! endif

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

       if (this%boxes(id)%maxendpoint_id(1).eq.289) then
          !print *,"SAP_ID=",this%id,"000000UPDATE    this%boxes(",289,"):   particle_id,value = ",this%boxes(id)%particle_id,aabb%maxendpoint(1)
          !print *,"SAP_ID=",this%id,"SAP_ID: ",this%id," 000000UPDATE    this%x_endpoints(",this%boxes(id)%maxendpoint_id(1),")%value = ",this%x_endpoints(this%boxes(id)%maxendpoint_id(1))%value
       endif
    end subroutine update_box

    subroutine sort(this)
      use geometry, only: do_k
      implicit none
      type(sap_t), intent(inout) :: this
      integer :: ii

      !if (.not. check_boxes(this)) stop __LINE__

      ! sort end points
      call sort_endpoints(this%x_endpoints(1:this%x_endpoints_len),this,1)
      !if (.not. check_boxes(this)) stop __LINE__
      call sort_endpoints(this%y_endpoints(1:this%y_endpoints_len),this,2)
      if (.not. check_boxes(this)) stop __LINE__
      if (do_k) call sort_endpoints(this%z_endpoints(1:this%z_endpoints_len),this,3)
      if (.not. check_boxes(this)) stop __LINE__

! cleanup HUGE endpoints
      ii = this%x_endpoints_len
      do
         if (ii .le. 0) exit
         if (this%x_endpoints(ii)%value .ne. HUGE(0.0)) exit
         print *,"SAPID: ",this%id, " endpoint_id= ",ii," IS HUGE",this%x_endpoints(ii)%box_id
!     mark corresponding box as deleted
         this%boxes(abs(this%x_endpoints(ii)%box_id))%minendpoint_id(:) = 0
         this%boxes(abs(this%x_endpoints(ii)%box_id))%maxendpoint_id(:) = 0
         call active_add(this%deleted_boxes,this%x_endpoints(ii)%box_id)
         ii = ii - 1
      enddo

      this%x_endpoints_len = ii
      this%y_endpoints_len = ii
      this%z_endpoints_len = ii

    end subroutine sort

    subroutine sweep(this,ss)

      use pair_manager, only: add_pair
      use geometry
      use discretelement

      implicit none
      type(sap_t), intent(inout) :: this
      integer, intent(in) :: ss

      ! active list of box id's
      integer :: ii,aa,minmax,ai

      type(active_t) :: active

      !stop __LINE__

      call active_init(active)

      do ii=1, this%x_endpoints_len

         minmax = this%x_endpoints(ii)%box_id

         if ( minmax < 0) then
            ! add pairs for new box x ( active pairs )

            do ai=1, active_get_length(active)
               aa = active_get(active,ai)
               ! compare other two axes

               if (max(this%boxes(aa)%minendpoint_id(2),this%boxes(-minmax)%minendpoint_id(2)) <= min(this%boxes(-minmax)%maxendpoint_id(2),this%boxes(aa)%maxendpoint_id(2)) .and. &
                    (NO_K .or. max(this%boxes(aa)%minendpoint_id(3),this%boxes(-minmax)%minendpoint_id(3)) <= min(this%boxes(-minmax)%maxendpoint_id(3),this%boxes(aa)%maxendpoint_id(3)))) then
                  !print *,"SAP_ID=",this%id,"FOUND PAIR! ADDING...",this%boxes(-minmax)%particle_id,this%boxes(aa)%particle_id
                  call add_pair(this%boxes(-minmax)%particle_id,this%boxes(aa)%particle_id)
               endif
            enddo

            ! add new box to active pair list
            call active_add(active,-minmax)

         else if ( 0 < minmax ) then
            !print *,"SAP_ID=",this%id,"minmax = ",minmax
            ! remove box from active pair list
            call active_del(active,minmax)
         else
            print *,"SAP_ID=",this%id,"minmax shouldn't be zero: ",minmax
            stop __LINE__
         endif
      enddo

    end subroutine sweep

      subroutine active_init(this)
        implicit none
        type(active_t), intent(inout) :: this

        allocate(this%list(10))
        this%list = 0
        this%list_len = 0

      end subroutine active_init

      subroutine active_add(this,new_active)
        use resize, only: integer_grow
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

      integer function active_pop(this)
        implicit none
        type(active_t), intent(inout) :: this

        if (this%list_len .eq. 0) then
           active_pop = -1
           return
        endif

        active_pop = this%list(this%list_len)
        this%list(this%list_len) = 0

        this%list_len = this%list_len - 1

      end function active_pop

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

       ! print *,"SAP_ID=",sap%id,"NOW SORTING AXIS ",axis,"############################################################"

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

          !print *,"SAP_ID=",sap%id,"EVERYTHING BELOW ii=",ii," / ",size(endpoints),"  IS SORTED"

          sweeppoint = endpoints(ii)
          sweepval = sweeppoint%value
          jj = ii-1

          do while ( 0 < jj )
             if ( endpoints(jj)%value <= sweepval ) then
                !print *,"SAP_ID=",sap%id," AHAH, ",endpoints(jj)%value , " IS LESS THAN ",sweepval," . DONE WITH SORTING",ii
                exit
             endif

             swappoint = endpoints(jj)

            !print *,"SAP_ID=",sap%id,"SORTSTART................................................",ii,jj
            !print *,"SAP_ID=",sap%id,"swappoint is jj=",jj,",  "," BOX:",swappoint%box_id," VAL:",swappoint%value
            !print *,"SAP_ID=",sap%id,"NOW COMPARING ENDPOINTS BELONGING TO BOXES:",sweeppoint%box_id,swappoint%box_id
            !print *,"SAP_ID=",sap%id,"endpoints(ii) = ",endpoints(ii)
            !print *,"SAP_ID=",sap%id,"endpoints(ii)%box_id = ",endpoints(ii)%box_id

            if (sweeppoint%box_id < 0 .and. 0 < swappoint%box_id ) then

               if (fullcheck(sap,-sweeppoint%box_id,swappoint%box_id,axis)) then
                  !print *,"SAP_ID=",sap%id,"PAIR FOUND",sap%boxes(-sweeppoint%box_id)%particle_id,sap%boxes(swappoint%box_id)%particle_id

                  if (105.eq. min(sap%boxes(-sweeppoint%box_id)%particle_id,sap%boxes(swappoint%box_id)%particle_id) .and. 106.eq.max(sap%boxes(-sweeppoint%box_id)%particle_id,sap%boxes(swappoint%box_id)%particle_id))then
                     print *,"SAP_ID=",sap%id," ADDING TO SAP:::::::::::::::::::",sap%id,ii,jj
                  endif

                  call add_pair(sap%boxes(-sweeppoint%box_id)%particle_id,sap%boxes(swappoint%box_id)%particle_id)
               else
                  !print *,"SAP_ID=",sap%id,"NOPE, FAILED FULLAXISCHECK"
               endif
            endif

            if (0 < sweeppoint%box_id .and. swappoint%box_id < 0 ) then
               !print *,"SAP_ID=",sap%id,"PAIR RMEOVEd"
               if ((2302 .eq. sap%boxes(sweeppoint%box_id)%particle_id .and. 2334 .eq. sap%boxes(-swappoint%box_id)%particle_id) .or. (2334 .eq. sap%boxes(sweeppoint%box_id)%particle_id .and. 2302 .eq. sap%boxes(-swappoint%box_id)%particle_id)) then
                  print *,"SAP_ID=",sap%id,"GOT SWEEPPOINT FROM ",ii
                  print *,"SAP_ID=",sap%id,"GOT SWAPPOINT FROM ",jj
                  print *,"SAP_ID=",sap%id,"sweeppoint%box_id ==== ",sweeppoint%box_id,"sweeppoint%particle_id ==== ",sap%boxes(sweeppoint%box_id)%particle_id
                  print *,"SAP_ID=",sap%id,"swappoint%box_id ==== ",swappoint%box_id,"swappoint%particle_id ==== ",sap%boxes(-swappoint%box_id)%particle_id
               endif

               if (fullcheck(sap,sweeppoint%box_id,-swappoint%box_id,axis)) then

                  if (105.eq. min(sap%boxes(sweeppoint%box_id)%particle_id,sap%boxes(-swappoint%box_id)%particle_id) .and. 106.eq.max(sap%boxes(sweeppoint%box_id)%particle_id,sap%boxes(-swappoint%box_id)%particle_id))then
                     print *,"SAP_ID=",sap%id," DELETING FROM SAP:::::::::::::::::::",sap%id,ii,jj
                  endif

                  call del_pair(sap%boxes(sweeppoint%box_id)%particle_id,sap%boxes(-swappoint%box_id)%particle_id)
                  endif
            endif

            endpoints(jj+1) = swappoint

             ! if (jj+1.eq.289) then
             !    print *,"SAP_ID=",sap%id,"00000000000000000SWAP    endpoints(",289,")%particle_id = ",sap%boxes(abs(endpoints(jj+1)%box_id))%particle_id
             !    print *,"SAP_ID=",sap%id,"00000000000000000SWAP    endpoints(",289,")%value = ",endpoints(289)%value
             ! endif

            ! if box_id is zero, that is an endpoint corresponding to a deleted box
            if (swappoint%box_id < 0) then
               sap%boxes(-swappoint%box_id)%minendpoint_id(axis) = jj+1
            else if (swappoint%box_id > 0) then
               sap%boxes(swappoint%box_id)%maxendpoint_id(axis) = jj+1
            endif

            jj = jj - 1
          enddo

          endpoints(jj+1) = sweeppoint

          ! if box_id is zero, that is an endpoint corresponding to a deleted box
          if (sweeppoint%box_id < 0) then
             sap%boxes(-sweeppoint%box_id)%minendpoint_id(axis) = jj+1
          else if (sweeppoint%box_id > 0) then
             sap%boxes(sweeppoint%box_id)%maxendpoint_id(axis) = jj+1
          endif

          !print *,"SAP_ID=",sap%id,"DONE SORTING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ",ii,jj
          !call print_boxes(sap)
          !if (.not. check_boxes(sap)) stop __LINE__

       enddo

       if (1.eq.axis) then
          do ii=1, sap%boxes_len
             if (sap%boxes(ii)%minendpoint_id(1) .eq. 0 ) cycle
             if (abs(endpoints(abs(sap%boxes(ii)%minendpoint_id(1)))%box_id) .ne. ii) then
             print *,"SAP_ID=",sap%id,"sap%boxes(",ii,")%minendpoint_id(1)",sap%boxes(ii)%minendpoint_id(1)
             print *,"SAP_ID=",sap%id,"endpoints(sap%boxes(",ii,")%minendpoint_id(1))%box_id",endpoints(abs(sap%boxes(ii)%minendpoint_id(1)))%box_id
             stop __LINE__
             return
          endif
          if (abs(endpoints(abs(sap%boxes(ii)%maxendpoint_id(1)))%box_id) .ne. ii) then
             print *,"SAP_ID=",sap%id,"sap%boxes(",ii,")%maxendpoint_id(1)",sap%boxes(ii)%maxendpoint_id(1)
             print *,"SAP_ID=",sap%id,"endpoints(sap%boxes(",ii,")%maxendpoint_id(1))%box_id",endpoints(abs(sap%boxes(ii)%maxendpoint_id(1)))%box_id
             stop __LINE__
             return
          endif
       enddo
    endif

    if (2.eq.axis) then
       do ii=1, sap%boxes_len
          if (sap%boxes(ii)%minendpoint_id(2) .eq. 0 ) cycle
          if (abs(endpoints(abs(sap%boxes(ii)%minendpoint_id(2)))%box_id) .ne. ii) then
             print *,"SAP_ID=",sap%id,"sap%boxes(",ii,")%minendpoint_id(2)",sap%boxes(ii)%minendpoint_id(2)
             print *,"SAP_ID=",sap%id,"endpoints(sap%boxes(",ii,")%minendpoint_id(2))%box_id",endpoints(abs(sap%boxes(ii)%minendpoint_id(2)))%box_id
             stop __LINE__
             return
          endif
          if (abs(endpoints(abs(sap%boxes(ii)%maxendpoint_id(2)))%box_id) .ne. ii) then
             print *,"SAP_ID=",sap%id,"sap%boxes(",ii,")%maxendpoint_id(2)",sap%boxes(ii)%maxendpoint_id(2)
             print *,"SAP_ID=",sap%id,"endpoints(sap%boxes(",ii,")%maxendpoint_id(2))%box_id",endpoints(abs(sap%boxes(ii)%maxendpoint_id(2)))%box_id
             stop __LINE__
             return
          endif
       enddo
    endif

    if (3.eq.axis) then
       do ii=1, sap%boxes_len
          if (sap%boxes(ii)%minendpoint_id(3) .eq. 0 ) cycle
          if (abs(endpoints(abs(sap%boxes(ii)%minendpoint_id(3)))%box_id) .ne. ii) then
             print *,"SAP_ID=",sap%id,"sap%boxes(",ii,")%minendpoint_id(3)",sap%boxes(ii)%minendpoint_id(3)
             print *,"SAP_ID=",sap%id,"endpoints(sap%boxes(",ii,")%minendpoint_id(3))%box_id",endpoints(abs(sap%boxes(ii)%minendpoint_id(3)))%box_id
             stop __LINE__
             return
          endif
          if (abs(endpoints(abs(sap%boxes(ii)%maxendpoint_id(3)))%box_id) .ne. ii) then
             print *,"SAP_ID=",sap%id,"sap%boxes(",ii,")%maxendpoint_id(3)",sap%boxes(ii)%maxendpoint_id(3)
             print *,"SAP_ID=",sap%id,"endpoints(sap%boxes(",ii,")%maxendpoint_id(3))%box_id",endpoints(abs(sap%boxes(ii)%maxendpoint_id(3)))%box_id
             stop __LINE__
             return
          endif
       enddo
    endif


       ! TODO...delete endpoints at the end of the list of size HUGE, box_id==0, and then reduce endpoints_len
     end subroutine sort_endpoints

     logical function fullcheck(this,id,id2,curr_axis)
       use geometry
       implicit none
       type(sap_t), intent(inout) :: this
       ! box ids to compare
       integer, intent(in) :: id, id2, curr_axis

        !print *,"SAP_ID=",this%id,"FULLCHECK FOR ",id,id2

        !print *,"SAP_ID=",this%id,"this%boxes(id)%minendpoint_id(1)==",this%boxes(id)%minendpoint_id(1),this%x_endpoints(this%boxes(id)%minendpoint_id(1))%value
        !print *,"SAP_ID=",this%id,"this%boxes(id2)%minendpoint_id(1)==",this%boxes(id2)%minendpoint_id(1),this%x_endpoints(this%boxes(id2)%minendpoint_id(1))%value
        !print *,"SAP_ID=",this%id,"this%boxes(id)%maxendpoint_id(1)==",this%boxes(id)%maxendpoint_id(1),this%x_endpoints(this%boxes(id)%maxendpoint_id(1))%value
        !print *,"SAP_ID=",this%id,"this%boxes(id2)%maxendpoint_id(1)==",this%boxes(id2)%maxendpoint_id(1),this%x_endpoints(this%boxes(id2)%maxendpoint_id(1))%value
        !print *,"SAP_ID=",this%id,""
        !print *,"SAP_ID=",this%id,"this%boxes(id)%minendpoint_id(2)==",this%boxes(id)%minendpoint_id(2),this%y_endpoints(this%boxes(id)%minendpoint_id(2))%value
        !print *,"SAP_ID=",this%id,"this%boxes(id2)%minendpoint_id(2)==",this%boxes(id2)%minendpoint_id(2),this%y_endpoints(this%boxes(id2)%minendpoint_id(2))%value
        !print *,"SAP_ID=",this%id,"this%boxes(id)%maxendpoint_id(2)==",this%boxes(id)%maxendpoint_id(2),this%y_endpoints(this%boxes(id)%maxendpoint_id(2))%value
        !print *,"SAP_ID=",this%id,"this%boxes(id2)%maxendpoint_id(2)==",this%boxes(id2)%maxendpoint_id(2),this%y_endpoints(this%boxes(id2)%maxendpoint_id(2))%value

       if (max(this%boxes(id)%minendpoint_id(1),this%boxes(id2)%minendpoint_id(1)) <= min(this%boxes(id2)%maxendpoint_id(1),this%boxes(id)%maxendpoint_id(1))) then
          !print *,"SAP_ID=",this%id,"OVERLAP ON X AXIS"
          else
             !print *,"SAP_ID=",this%id,"NOVERLAP ON X AXIS"
             endif
             if (max(this%boxes(id)%minendpoint_id(2),this%boxes(id2)%minendpoint_id(2)) <= min(this%boxes(id2)%maxendpoint_id(2),this%boxes(id)%maxendpoint_id(2))) then
                !print *,"SAP_ID=",this%id,"OVERLAP ON Y AXIS"
             else
                !print *,"SAP_ID=",this%id,"NOVERLAP ON Y AXIS"
             endif

             if (NO_K .or. max(this%boxes(id)%minendpoint_id(3),this%boxes(id2)%minendpoint_id(3)) <= min(this%boxes(id2)%maxendpoint_id(3),this%boxes(id)%maxendpoint_id(3))) then
                !print *,"SAP_ID=",this%id,"OVERLAP ON Z AXIS"
    else
       !print *,"SAP_ID=",this%id,"NOVERLAP ON Z AXIS"
    endif

    fullcheck = ((curr_axis.eq.1 .or. max(this%boxes(id)%minendpoint_id(1),this%boxes(id2)%minendpoint_id(1)) <= min(this%boxes(id2)%maxendpoint_id(1),this%boxes(id)%maxendpoint_id(1))) &
         .and. (curr_axis.eq.2 .or. max(this%boxes(id)%minendpoint_id(2),this%boxes(id2)%minendpoint_id(2)) <= min(this%boxes(id2)%maxendpoint_id(2),this%boxes(id)%maxendpoint_id(2))) &
         .and. (NO_K .or. curr_axis.eq.3 .or. max(this%boxes(id)%minendpoint_id(3),this%boxes(id2)%minendpoint_id(3)) <= min(this%boxes(id2)%maxendpoint_id(3),this%boxes(id)%maxendpoint_id(3))))

     end function fullcheck

    recursive subroutine quicksort_endpoints(A, sap, axis)
      type(sap_t), intent(inout) :: sap
      integer, intent(in) :: axis
      type(endpoint_t), intent(in out), dimension(:) :: A
      integer :: iq,ii

      !call check_boxes(sap)

      if(size(A) > 1) then
         call Partition(A, iq, sap, axis)
         call quicksort_endpoints(A(:iq-1),sap,axis)
         call quicksort_endpoints(A(iq:),sap,axis)
      endif

      do ii=2, size(A)
         if (A(ii)%value < A(ii-1)%value) then
            print *,"SAP_ID=",sap%id,"********************************************************************************************"
            print *,"SAP_ID=",sap%id,"ii-1:",ii-1,"  endpoints(ii):",A%box_id,A%value
            print *,"SAP_ID=",sap%id,"ii:",ii,"  endpoints(ii):",A%box_id,A%value
            print *,"SAP_ID=",sap%id,"*******************************************************************************************"
            stop __LINE__
            return
         endif
      enddo

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
            !print *,"SAP_ID=",sap%id,"partioning...",i,j
            !call print_boxes(sap)
            !call check_boxes(sap)

            !print *,"SAP_ID=",sap%id,"NOW SWAPPING ENDPOINTS BELONGING TO BOXES:",A(i)%box_id,A(j)%box_id

            !print *,"SAP_ID=",sap%id,"i = ",i
            !print *,"SAP_ID=",sap%id,"size(A) = ",size(A)
            !print *,"SAP_ID=",sap%id,"A(i) = ",A(i)
            !print *,"SAP_ID=",sap%id,"A(i)%box_id = ",A(i)%box_id

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
            !print *,"SAP_ID=",sap%id,"SET ENDPOINT TO VALUE ",A(i)%value," WHICH BELONG TO BOX ",abs(A(i)%box_id)
            A(j) = temp
            !print *,"SAP_ID=",sap%id,"SET ENDPOINT TO VALUE ",A(j)%value," WHICH BELONG TO BOX ",abs(A(j)%box_id)

            if (A(i)%box_id < 0) then
               !print *,"SAP_ID=",sap%id,j,"setting min endpoint of box ",-A(j)%box_id, " ON AXIS ",axis," TO ",i
               sap%boxes(-A(i)%box_id)%minendpoint_id(axis) = tmp_ii
            else
               !print *,"SAP_ID=",sap%id,j,"setting max endpoint of box ",A(j)%box_id, " ON AXIS ",axis," TO ",i
               sap%boxes(A(i)%box_id)%maxendpoint_id(axis) = tmp_ii
            endif
            if (A(j)%box_id < 0) then
               !print *,"SAP_ID=",sap%id,j,"setting min endpoint of box ",-A(j)%box_id, " ON AXIS ",axis," TO ",i
               sap%boxes(-A(j)%box_id)%minendpoint_id(axis) = tmp_jj
            else
               !print *,"SAP_ID=",sap%id,j,"setting max endpoint of box ",A(j)%box_id, " ON AXIS ",axis," TO ",i
               sap%boxes(A(j)%box_id)%maxendpoint_id(axis) = tmp_jj
            endif

            !print *,"SAP_ID=",sap%id,"partioned! ",i,j
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
