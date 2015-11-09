module multi_sweep_and_prune

  use sweep_and_prune

  type multisap_t
     type(sap_t), dimension(:), allocatable :: saps
     integer :: saps_len
     real :: minbounds(3), maxbounds(3)
     real :: one_over_cell_length(3)
     integer :: grid(3)

  end type multisap_t

  ! the global multisap
  type(multisap_t) multisap

contains

  subroutine init_multisap(this,x_grid,y_grid,z_grid2,minbounds,maxbounds)
    use geometry
    implicit none
    type(multisap_t), intent(inout) :: this
    integer, intent(in) :: x_grid,y_grid, z_grid2
    real, dimension(3), intent(in) :: minbounds, maxbounds
    integer :: ii,jj,kk, z_grid


    if (NO_K) then
       z_grid = 1
    else
       z_grid = z_grid2
    endif

    this%grid(1) = x_grid
    this%grid(2) = y_grid
    this%grid(3) = z_grid

    allocate(this%saps(0:x_grid*y_grid*z_grid-1))
    do ii=0,x_grid-1
       do jj=0,y_grid-1
          do kk=0,z_grid-1
             call init_sap(this%saps(ii*y_grid*z_grid+jj*z_grid+kk))
             this%saps(ii*y_grid*z_grid+jj*z_grid+kk)%id = (ii*y_grid*z_grid+jj*z_grid+kk)
          enddo
       enddo
    enddo

    this%one_over_cell_length(:) = this%grid(:)/(maxbounds(:)-minbounds(:))

    this%minbounds(:) = minbounds(:)
    this%maxbounds(:) = maxbounds(:)

  end subroutine init_multisap

  subroutine multisap_raster(this,aabb,sap_ids,particle_id)
    implicit none
    type(multisap_t), intent(inout) :: this
    type(aabb_t), intent(in) :: aabb
    integer, intent(in) :: particle_id
    integer, intent(out) :: sap_ids(MAX_SAPS)
    integer :: max_grid(3),min_grid(3)

    sap_ids(:) = -1

    if (particle_id.eq.46) then
       print *,"aabb%minendpoint(:) = ",aabb%minendpoint(:)
       print *,"aabb%maxendpoint(:) = ",aabb%maxendpoint(:)

       print *,"this%minbounds(:) = ",this%minbounds(:)
       print *,"this%maxbounds(:) = ",this%maxbounds(:)

       print *,"this%one_over_cell_length(:)",this%one_over_cell_length(:)
    endif

    ! bin to grid
    min_grid(:) = floor((aabb%minendpoint(:)-this%minbounds(:))*this%one_over_cell_length(:))
    max_grid(:) = floor((aabb%maxendpoint(:)-this%minbounds(:))*this%one_over_cell_length(:))

    if (particle_id.eq.46) then
       print *,"min_grid:  ",min_grid(:)
       print *,"aabb%maxendpoint(:) = ",aabb%maxendpoint(:)
       print *,"(aabb%maxendpoint(:)-this%maxbounds(:) = ",(aabb%maxendpoint(:)-this%minbounds(:))
       print *,"(aabb%maxendpoint(:)-this%maxbounds(:))*this%one_over_cell_length(:) = ",(aabb%maxendpoint(:)-this%minbounds(:))*this%one_over_cell_length(:)
       print *,"max_grid:  ",max_grid(:)
    endif

    min_grid(:) = max(min_grid(:),0)
    max_grid(:) = min(max_grid(:),this%grid(:)-1)

    !if (any(min_grid.ne.0)) stop __LINE__
    !if (any(max_grid.ne.0)) stop __LINE__

    if (particle_id.eq.46) then
       print *,"min_grid:  ",min_grid(:)
       print *,"max_grid:  ",max_grid(:)
    endif

    call add_to_list((min_grid(1)*this%grid(2) + min_grid(2))*this%grid(3) + min_grid(3))
    call add_to_list((min_grid(1)*this%grid(2) + min_grid(2))*this%grid(3) + max_grid(3))
    call add_to_list((min_grid(1)*this%grid(2) + max_grid(2))*this%grid(3) + min_grid(3))
    call add_to_list((min_grid(1)*this%grid(2) + max_grid(2))*this%grid(3) + max_grid(3))
    call add_to_list((max_grid(1)*this%grid(2) + min_grid(2))*this%grid(3) + min_grid(3))
    call add_to_list((max_grid(1)*this%grid(2) + min_grid(2))*this%grid(3) + max_grid(3))
    call add_to_list((max_grid(1)*this%grid(2) + max_grid(2))*this%grid(3) + min_grid(3))
    call add_to_list((max_grid(1)*this%grid(2) + max_grid(2))*this%grid(3) + max_grid(3))

  contains

    subroutine add_to_list(nn)
      integer, intent(in) :: nn
      integer :: mm

      do mm=1, size(sap_ids)
         if (sap_ids(mm).eq.nn) then
            ! already in the list
            return
         endif
         if (sap_ids(mm).eq.-1) then
            sap_ids(mm) = nn
            return
         endif
      enddo

    end subroutine add_to_list

  end subroutine multisap_raster

  subroutine multisap_add(this,aabb,particle_id,handlelist)
    implicit none
    type(multisap_t), intent(inout) :: this
    type(aabb_t), intent(in) :: aabb
    integer, intent(in) :: particle_id
    type(boxhandlelist_t), intent(out) :: handlelist
    integer :: sap_ids(MAX_SAPS)
    integer :: nn,id

    call multisap_raster(this,aabb,sap_ids,particle_id)

    ! add to each individual SAP
    do nn=1, size(sap_ids)
       handlelist%list(nn)%sap_id = sap_ids(nn)
       if (sap_ids(nn) >= 0) then
          call add_box(this%saps(sap_ids(nn)),aabb,particle_id,handlelist%list(nn)%box_id)
       endif
    enddo

  end subroutine multisap_add

  subroutine multisap_del(this,aabb,handlelist)
    implicit none
    type(multisap_t), intent(inout) :: this
    type(aabb_t), intent(in) :: aabb
    type(boxhandlelist_t), intent(in) :: handlelist
    integer :: nn

    ! remove from each SAP listed for this id
    do nn=1, size(handlelist%list)
       if (handlelist%list(nn)%sap_id >= 0) then
          call del_box(this%saps(handlelist%list(nn)%sap_id),handlelist%list(nn)%box_id)
       endif
    enddo

  end subroutine multisap_del

  subroutine multisap_update(this,aabb,handlelist)
    implicit none
    type(multisap_t), intent(inout) :: this
    type(aabb_t), intent(in) :: aabb
    type(boxhandlelist_t), intent(inout) :: handlelist
    integer :: particle_id
    integer, DIMENSION(MAX_SAPS) :: new_sap_ids
    logical :: found
    integer :: mm,nn,first_blank

    particle_id = this%saps(handlelist%list(1)%sap_id)%boxes(handlelist%list(1)%box_id)%particle_id

    call multisap_raster(this,aabb,new_sap_ids,-1)

    ! update for each SAP listed for this id
    do nn=1, size(new_sap_ids)
       if (new_sap_ids(nn) < 0) exit

       found = .false.
       do mm=1, size(handlelist%list)
          if (handlelist%list(mm)%sap_id < 0) first_blank = mm

          if (handlelist%list(mm)%sap_id .eq. new_sap_ids(nn)) then
             ! update existing SAPs
             call update_box(this%saps(new_sap_ids(nn)),handlelist%list(mm)%box_id,aabb)
             found = .true.
             exit
          endif
       enddo
       if (.not.found) then
          ! add SAPs yet not listed for this id
          call add_box(this%saps(new_sap_ids(nn)),aabb,particle_id,handlelist%list(first_blank)%box_id)
          handlelist%list(first_blank)%sap_id = new_sap_ids(nn)
       endif
    enddo

    do mm=1, size(handlelist%list)
       if (handlelist%list(mm)%sap_id < 0) cycle

       found = .false.
       do nn=1, size(new_sap_ids)
          if (new_sap_ids(nn) < 0) exit

          if (handlelist%list(mm)%sap_id .eq. new_sap_ids(nn)) then
             found = .true.
             exit
          endif
       enddo

       if (.not.found) then
          ! remove SAP not found in handle%sap_id
          call del_box(this%saps(handlelist%list(mm)%sap_id),handlelist%list(mm)%box_id)
          handlelist%list(mm)%sap_id = -1
       endif
    enddo

  end subroutine multisap_update

  ! logical function inarray(aa,bb,array)
  !   implicit none
  !   integer, intent(in) :: aa,bb
  !   type(boxhandle_t), intent(in) :: array
  !   integer :: ii

  !   if (size(array%sap_id) .ne. size(array%box_id)) then
  !      stop __LINE__
  !   endif

  !   do ii=1,size(array%sap_id)
  !      if (aa.eq.array%sap_id(ii) .and. bb.eq.array%box_id(ii)) then
  !         inarray = .false.
  !         return
  !      endif
  !   enddo
  !   inarray = .false.
  ! end function inarray

  subroutine multisap_sort(this)
    USE discretelement
    implicit none
    type(multisap_t), intent(inout) :: this
    integer :: ii, jj, kk
    type(sap_t) :: sap
    type(box_t) :: box

    do ii=0,this%grid(1)-1
       do jj=0,this%grid(2)-1
          do kk=0,this%grid(3)-1
             !print *,"NOW GOING TO SORT SAP::",ii,jj,kk,":::   ",ii*this%grid(2)*this%grid(3)+jj*this%grid(3)+kk
             !if (.not. check_boxes(this%saps(ii*this%grid(2)*this%grid(3)+jj*this%grid(3)+kk))) stop __LINE__
             call sort(this%saps(ii*this%grid(2)*this%grid(3)+jj*this%grid(3)+kk))
             !if (.not. check_boxes(this%saps(ii*this%grid(2)*this%grid(3)+jj*this%grid(3)+kk))) stop __LINE__
          enddo
       enddo
    enddo

  end subroutine multisap_sort

  subroutine multisap_quicksort(this)
    implicit none
    type(multisap_t), intent(inout) :: this
    integer :: ii, jj, kk
    integer :: sap_id

    do ii=0,this%grid(1)-1
       do jj=0,this%grid(2)-1
          do kk=0,this%grid(3)-1
             sap_id = ii*this%grid(2)*this%grid(3)+jj*this%grid(3)+kk
             print *,"NOW GOING DO QUICKSORT THE THING:::::   ",sap_id
             call do_quicksort(this,this%saps(ii*this%grid(2)*this%grid(3)+jj*this%grid(3)+kk))

             !if (.not.check_boxes(multisap%saps(sap_id))) stop __LINE__
             if (.not.check_sort(multisap%saps(sap_id))) stop __LINE__
          enddo
       enddo
    enddo
  end subroutine multisap_quicksort

  subroutine do_quicksort(this,sap)
    implicit none
    type(multisap_t), intent(inout) :: this
    type(sap_t), intent(inout) :: sap
    integer :: nn

    call quicksort_endpoints(sap%x_endpoints(1:sap%x_endpoints_len),sap,1)
    !if (.not.check_boxes(sap)) stop __LINE__
    !if (.not.check_sort(sap)) stop __LINE__

    call quicksort_endpoints(sap%y_endpoints(1:sap%y_endpoints_len),sap,2)
    !if (.not.check_boxes(sap)) stop __LINE__
    !if (.not.check_sort(sap)) stop __LINE__

    call quicksort_endpoints(sap%z_endpoints(1:sap%z_endpoints_len),sap,3)
    !if (.not.check_boxes(sap)) stop __LINE__
    if (.not.check_sort(sap)) stop __LINE__

  end subroutine do_quicksort

  subroutine multisap_sweep(this)
    USE discretelement
    use pair_manager
    implicit none
    type(multisap_t), intent(inout) :: this
    integer :: ii, jj, kk
    integer :: ll, i, ss
    integer :: minenx, mineny, minenz, minenx2, mineny2, minenz2
    integer :: maxenx, maxeny, maxenz, maxenx2, maxeny2, maxenz2

    do ii=0,this%grid(1)-1
       do jj=0,this%grid(2)-1
          do kk=0,this%grid(3)-1
             ss = ii*this%grid(2)*this%grid(3)+jj*this%grid(3)+kk
             print *,"NOW GOING TO sweep THE THING:::::   ",ss
             call sweep(this%saps(ii*this%grid(2)*this%grid(3)+jj*this%grid(3)+kk),ss)

             if (.not.check_boxes( this%saps(ii*this%grid(2)*this%grid(3)+jj*this%grid(3)+kk) )) stop __LINE__

          enddo
       enddo
    enddo

  end subroutine multisap_sweep

  SUBROUTINE boxhandle_GROW(boxhandles,new_size)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: new_size
    type(boxhandlelist_t), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: boxhandles
    type(boxhandlelist_t), DIMENSION(:), ALLOCATABLE :: boxhandle_tmp
    INTEGER lSIZE

    lSIZE = size(boxhandles,1)
    allocate(boxhandle_tmp(new_size))
    boxhandle_tmp(1:lSIZE) = boxhandles(1:lSIZE)
    call move_alloc(boxhandle_tmp,boxhandles)

  END SUBROUTINE boxhandle_GROW

end module multi_sweep_and_prune
