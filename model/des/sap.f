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
  end type box_t

  type sap_t
     type(endpoint_t), dimension(:), allocatable :: x_endpoints
     type(endpoint_t), dimension(:), allocatable :: y_endpoints
     type(endpoint_t), dimension(:), allocatable :: z_endpoints
     integer, allocatable :: x_endpoints_len
     integer, allocatable :: y_endpoints_len
     integer, allocatable :: z_endpoints_len
     type(box_t), dimension(:), allocatable :: boxes
     integer :: boxes_len
  end type sap_t

  contains

    subroutine add_box(this,aabb,id)
      implicit none
      type(sap_t), intent(inout) :: this
      type(aabb_t), intent(in) :: aabb
      integer, intent(out) :: id

      if (this%x_endpoints_len+2 <= size(this%x_endpoints)) then
         !integer_grow(x_endpoints,x_endpoints_len)
      endif

      if (this%y_endpoints_len+2 <= size(this%y_endpoints)) then
         !integer_grow(y_endpoints,y_endpoints_len)
      endif

      if (this%z_endpoints_len+2 <= size(this%z_endpoints)) then
         !integer_grow(z_endpoints,z_endpoints_len)
      endif

      if (this%boxes_len+1 <= size(this%boxes)) then
         !integer_grow(boxes,this%boxes)
      endif

      this%x_endpoints(this%x_endpoints_len)%box_id = this%boxes_len
      this%x_endpoints(this%x_endpoints_len)%value = aabb%minendpoint(1)
      this%y_endpoints(this%y_endpoints_len)%box_id = this%boxes_len
      this%y_endpoints(this%y_endpoints_len)%value = aabb%minendpoint(2)
      this%z_endpoints(this%z_endpoints_len)%box_id = this%boxes_len
      this%z_endpoints(this%z_endpoints_len)%value = aabb%minendpoint(3)

      this%x_endpoints(this%x_endpoints_len+1)%box_id = this%boxes_len
      this%x_endpoints(this%x_endpoints_len+1)%value = aabb%maxendpoint(1)
      this%y_endpoints(this%y_endpoints_len+1)%box_id = this%boxes_len
      this%y_endpoints(this%y_endpoints_len+1)%value = aabb%maxendpoint(2)
      this%z_endpoints(this%z_endpoints_len+1)%box_id = this%boxes_len
      this%z_endpoints(this%z_endpoints_len+1)%value = aabb%maxendpoint(3)

      this%x_endpoints_len = this%x_endpoints_len + 2
      this%y_endpoints_len = this%y_endpoints_len + 2
      this%z_endpoints_len = this%z_endpoints_len + 2

      ! boxes(boxes_len)%minendpoint_id(:) = x_endpoints_len
      ! boxes(boxes_len)%maxendpoint_id(:) = x_endpoints_len+1
      this%boxes_len = this%boxes_len + 1

      id = this%boxes_len

      ! sort x_endpoints
      ! sort y_endpoints
      ! sort z_endpoints

    end subroutine add_box

    subroutine del_box(this,id)
      implicit none
      type(sap_t), intent(in) :: this
      integer, intent(in) :: id

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

      ! sort end points
      call sort_endpoints(this%x_endpoints)
      call sort_endpoints(this%y_endpoints)
      call sort_endpoints(this%z_endpoints)

    end subroutine update_box

    subroutine sort_endpoints(endpoints)
      implicit none
      type(endpoint_t), dimension(:) :: endpoints
      type(endpoint_t) :: swap
      integer :: ii, jj
      do ii=1, size(endpoints)
         jj = ii
         do while ((0<jj) .and. (endpoints(jj)%value < endpoints(jj-1)%value))
            swap = endpoints(jj)
            endpoints(jj) = endpoints(ii)
            endpoints(ii) = swap
            jj = jj - 1
         enddo
      enddo
    end subroutine sort_endpoints

end module sweep_and_prune
