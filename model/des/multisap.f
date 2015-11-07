module multi_sweep_and_prune

  use sweep_and_prune

  type multisap_t
     type(sap_t), dimension(:), allocatable :: saps
     integer :: saps_len
     real :: minbounds(3), maxbounds(3)
     real :: one_over_xcell_length, one_over_ycell_length, one_over_zcell_length
     integer :: x_grid, y_grid, z_grid

  end type multisap_t

  type boxhandle_t
     integer :: sap_id(8)
     integer :: box_id(8)
  end type boxhandle_t

  ! the global multisap
  type(multisap_t) multisap

  type(boxhandle_t), DIMENSION(:),  ALLOCATABLE :: boxhandle         !(PARTICLES)

contains

  subroutine init_multisap(this,x_grid,y_grid,z_grid,minbounds,maxbounds)
    implicit none
    type(multisap_t), intent(inout) :: this
    integer, intent(in) :: x_grid,y_grid, z_grid
    real, dimension(3), intent(in) :: minbounds, maxbounds
    integer :: ii,jj,kk

    this%x_grid = x_grid
    this%y_grid = y_grid
    this%z_grid = z_grid

    allocate(this%saps(0:this%x_grid*this%y_grid*this%z_grid))
    do ii=0,x_grid-1
       do jj=0,y_grid-1
          do kk=0,z_grid-1
             call init_sap(this%saps(ii*this%y_grid*this%z_grid+jj*this%z_grid+kk))
          enddo
       enddo
    enddo

    this%one_over_xcell_length = x_grid/(maxbounds(1)-minbounds(1))
    this%one_over_ycell_length = y_grid/(maxbounds(2)-minbounds(2))
    this%one_over_zcell_length = z_grid/(maxbounds(2)-minbounds(2))

    this%minbounds(:) = minbounds(:)
    this%maxbounds(:) = maxbounds(:)

  end subroutine init_multisap

  subroutine multisap_raster(this,aabb,sap_ids)
    implicit none
    type(multisap_t), intent(inout) :: this
    type(aabb_t), intent(in) :: aabb
    integer, intent(out) :: sap_ids(8)
    integer :: xx,yy,zz,x2,y2,z2

    sap_ids(:) = -1

    ! bin to grid
    xx = floor((aabb%minendpoint(1)-this%minbounds(1))*this%one_over_xcell_length)
    yy = floor((aabb%minendpoint(2)-this%minbounds(2))*this%one_over_ycell_length)
    zz = floor((aabb%minendpoint(3)-this%minbounds(3))*this%one_over_zcell_length)

    x2 = floor((aabb%maxendpoint(1)-this%minbounds(1))*this%one_over_xcell_length)
    y2 = floor((aabb%maxendpoint(2)-this%minbounds(2))*this%one_over_ycell_length)
    z2 = floor((aabb%maxendpoint(3)-this%minbounds(3))*this%one_over_zcell_length)

    call add_to_list(xx*this%y_grid*this%z_grid+yy*this%z_grid+zz)
    call add_to_list(xx*this%y_grid*this%z_grid+yy*this%z_grid+z2)
    call add_to_list(xx*this%y_grid*this%z_grid+y2*this%z_grid+zz)
    call add_to_list(xx*this%y_grid*this%z_grid+y2*this%z_grid+z2)
    call add_to_list(x2*this%y_grid*this%z_grid+yy*this%z_grid+zz)
    call add_to_list(x2*this%y_grid*this%z_grid+yy*this%z_grid+z2)
    call add_to_list(x2*this%y_grid*this%z_grid+y2*this%z_grid+zz)
    call add_to_list(x2*this%y_grid*this%z_grid+y2*this%z_grid+z2)

  contains

    subroutine add_to_list(nn)
      integer, intent(in) :: nn
      integer :: first_blank, mm

      first_blank = -1
      do mm=1, size(sap_ids)
         if (sap_ids(mm).eq.-1) then
            first_blank = mm
         endif
         if (sap_ids(mm).eq.nn) then
            ! already in the list
            return
         endif
      enddo
      sap_ids(first_blank) = nn

    end subroutine add_to_list

  end subroutine multisap_raster

  subroutine multisap_add(this,aabb,handle)
    implicit none
    type(multisap_t), intent(inout) :: this
    type(aabb_t), intent(in) :: aabb
    type(boxhandle_t), intent(out) :: handle
    integer :: sap_ids(8)
    integer :: nn,id

    call multisap_raster(this,aabb,sap_ids)

    ! add to each individual SAP
    do nn=1, size(sap_ids)
       handle%sap_id(nn) = sap_ids(nn)
       if (sap_ids(nn) >= 0) then
          call add_box(this%saps(sap_ids(nn)),aabb,handle%box_id(nn))
       endif
    enddo

  end subroutine multisap_add

  subroutine multisap_del(this,aabb,handle)
    implicit none
    type(multisap_t), intent(inout) :: this
    type(aabb_t), intent(in) :: aabb
    type(boxhandle_t), intent(in) :: handle
    integer :: nn

    ! remove from each SAP listed for this id
    do nn=1, size(handle%sap_id)
       if (handle%sap_id(nn) >= 0) then
          call del_box(this%saps(handle%sap_id(nn)),handle%box_id(nn))
       endif
    enddo

  end subroutine multisap_del

  subroutine multisap_update(this,aabb,handle)
    implicit none
    type(multisap_t), intent(inout) :: this
    type(aabb_t), intent(in) :: aabb
    type(boxhandle_t), intent(inout) :: handle
    integer :: sap_ids(8)
    logical :: found
    integer :: mm,nn

    call multisap_raster(this,aabb,sap_ids)

    ! update for each SAP listed for this id
    do nn=1, size(handle%sap_id)
       found = .false.
       do mm=1, size(sap_ids)
          !print *,"size = ",size(this%sap_ids,1),,size(this%sap_ids,2)
          if (sap_ids(mm) .eq. handle%sap_id(nn) .and. sap_ids(mm) > 0) then
             ! update SAPs listed for this id
             !print *,"sap_ids(nn) == ",sap_ids(nn)
             !print *,"this%sap_ids(nn,id) == ",this%sap_ids(mm,id)
             call update_box(this%saps(handle%sap_id(nn)),handle%box_id(nn),aabb)
             found = .true.
             exit
          endif
       enddo
       if (.not.found .and. handle%sap_id(nn) > 0) then
          ! add to SAPs not listed for this id
          call del_box(this%saps(handle%sap_id(nn)),handle%box_id(nn))
       endif
    enddo

    do mm=1, size(sap_ids)
       found = .false.
       do nn=1, size(handle%sap_id)
          if (sap_ids(mm) .eq. handle%sap_id(nn) .and. sap_ids(mm) > 0) then
             found = .true.
             exit
          endif
          if (.not.found .and. handle%sap_id(nn) > 0) then
             ! remove from SAPs not listed for this id
             call del_box(this%saps(handle%sap_id(nn)),handle%box_id(nn))
          endif
       enddo
    enddo

  end subroutine multisap_update

  subroutine multisap_sort(this,aabb)
    implicit none
    type(multisap_t), intent(inout) :: this
    type(aabb_t), intent(in) :: aabb
    integer :: ii, jj, kk

    do ii=1,this%x_grid
       do jj=1,this%y_grid
          do kk=1,this%z_grid
             call sort(this%saps(ii*this%y_grid*this%z_grid+jj*this%z_grid+kk))
          enddo
       enddo
    enddo

  end subroutine multisap_sort

  SUBROUTINE boxhandle_GROW(boxhandle_array,new_size)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: new_size
    type(boxhandle_t), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: boxhandle_array
    type(boxhandle_t), DIMENSION(:), ALLOCATABLE :: boxhandle_tmp
    INTEGER lSIZE

    lSIZE = size(boxhandle_array,1)
    allocate(boxhandle_tmp(new_size))
    boxhandle_tmp(1:lSIZE) = boxhandle_array(1:lSIZE)
    call move_alloc(boxhandle_tmp,boxhandle_array)

  END SUBROUTINE boxhandle_GROW

end module multi_sweep_and_prune
