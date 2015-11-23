module multi_sweep_and_prune

  use sweep_and_prune

  type multisap_t
     type(sap_t), dimension(:), allocatable :: saps
     integer :: saps_len
     real :: minbounds(3), maxbounds(3)
     real :: one_over_cell_length(3)
     integer :: grid(3)
     type(hashtable_t) :: hashtable
  end type multisap_t

contains

  subroutine init_multisap(this,x_grid,y_grid,z_grid2,minbounds,maxbounds)
    use geometry
    implicit none
    type(multisap_t), intent(inout) :: this
    integer, intent(in) :: x_grid,y_grid, z_grid2
    real, dimension(3), intent(in) :: minbounds, maxbounds
    integer :: ii,jj,kk, z_grid, sap_id


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
             sap_id = (ii*this%grid(2)+jj)*this%grid(3)+kk
             print *,"INITIALIZING...",sap_id
             call init_sap(this%saps(sap_id),sap_id)
          enddo
       enddo
    enddo

    this%one_over_cell_length(:) = this%grid(:)/(maxbounds(:)-minbounds(:))

    this%minbounds(:) = minbounds(:)
    this%maxbounds(:) = maxbounds(:)

  end subroutine init_multisap

  subroutine multisap_raster(this,aabb,sap_ids,debug)
    implicit none
    type(multisap_t), intent(inout) :: this
    type(aabb_t), intent(in) :: aabb
    logical, intent(in) :: debug
    integer, intent(out) :: sap_ids(MAX_SAPS)
    integer :: max_grid(3),min_grid(3)

    sap_ids(:) = -1

    if (debug) then
       print *,"aabb%minendpoint(:) = ",aabb%minendpoint(:)
       print *,"aabb%maxendpoint(:) = ",aabb%maxendpoint(:)

       print *,"this%minbounds(:) = ",this%minbounds(:)
       print *,"this%maxbounds(:) = ",this%maxbounds(:)

       print *,"this%one_over_cell_length(:)",this%one_over_cell_length(:)
    endif

    ! bin to grid
    min_grid(:) = floor((aabb%minendpoint(:)-this%minbounds(:))*this%one_over_cell_length(:))
    max_grid(:) = floor((aabb%maxendpoint(:)-this%minbounds(:))*this%one_over_cell_length(:))

    if (debug) then
       print *,"min_grid:  ",min_grid(:)
       print *,"aabb%maxendpoint(:) = ",aabb%maxendpoint(:)
       print *,"(aabb%maxendpoint(:)-this%maxbounds(:) = ",(aabb%maxendpoint(:)-this%minbounds(:))
       print *,"(aabb%maxendpoint(:)-this%maxbounds(:))*this%one_over_cell_length(:) = ",(aabb%maxendpoint(:)-this%minbounds(:))*this%one_over_cell_length(:)
       print *,"max_grid:  ",max_grid(:)
    endif

    min_grid(:) = max(min_grid(:),0)
    max_grid(:) = min(max_grid(:),this%grid(:)-1)

    !if (any(min_grid.ne.0)) error stop __LINE__
    !if (any(max_grid.ne.0)) error stop __LINE__

    if (debug) then
       print *,"NORMALIZED min_grid:  ",min_grid(:)
       print *,"NORMALIZED max_grid:  ",max_grid(:)
    endif

    call add_to_list((min_grid(1)*this%grid(2) + min_grid(2))*this%grid(3) + min_grid(3))
    call add_to_list((min_grid(1)*this%grid(2) + min_grid(2))*this%grid(3) + max_grid(3))
    call add_to_list((min_grid(1)*this%grid(2) + max_grid(2))*this%grid(3) + min_grid(3))
    call add_to_list((min_grid(1)*this%grid(2) + max_grid(2))*this%grid(3) + max_grid(3))
    call add_to_list((max_grid(1)*this%grid(2) + min_grid(2))*this%grid(3) + min_grid(3))
    call add_to_list((max_grid(1)*this%grid(2) + min_grid(2))*this%grid(3) + max_grid(3))
    call add_to_list((max_grid(1)*this%grid(2) + max_grid(2))*this%grid(3) + min_grid(3))
    call add_to_list((max_grid(1)*this%grid(2) + max_grid(2))*this%grid(3) + max_grid(3))

    if (debug) then
       print *,"AND THUS IT COMES TO:::::    ",sap_ids(:)
    endif

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

    call multisap_raster(this,aabb,sap_ids,.false.)

    handlelist%particle_id = particle_id

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
          print *,"DELETED FROM SAP......",handlelist%list(nn)%sap_id
          call del_box(this%saps(handlelist%list(nn)%sap_id),handlelist%list(nn)%box_id)
       endif
    enddo

  end subroutine multisap_del

  subroutine multisap_update(this,aabb,handlelist)
    implicit none
    type(multisap_t), intent(inout) :: this
    type(aabb_t), intent(in) :: aabb
    type(boxhandlelist_t), intent(inout) :: handlelist
    ! integer :: particle_id
    integer, DIMENSION(MAX_SAPS) :: new_sap_ids
    logical :: found
    integer :: mm,nn,first_blank, ii, box_id
    real ::  asdf,diff
    logical :: debug

    ! particle_id = -1
    ! do nn=1, size(handlelist%list)
    !    if (handlelist%list(nn)%sap_id < 0) cycle
    !    particle_id = this%saps(handlelist%list(nn)%sap_id)%boxes(handlelist%list(nn)%box_id)%particle_id
    ! enddo
    ! if (particle_id < 0 ) error stop __LINE__

    debug = (handlelist%particle_id.eq.114 .or. handlelist%particle_id.eq. 115)
    ! if (debug) print *,""
    ! if (debug) print *," PARTICLE ",handlelist%particle_id
    call multisap_raster(this,aabb,new_sap_ids,.false.)

    ! update for each SAP listed for this id
    do nn=1, size(new_sap_ids)
       if (new_sap_ids(nn) < 0) exit

       found = .false.
       first_blank = -1
       do mm=1, size(handlelist%list)
          if (handlelist%list(mm)%sap_id < 0) first_blank = mm

          if (handlelist%list(mm)%sap_id .eq. new_sap_ids(nn)) then
             ! update existing SAPs

             box_id = handlelist%list(mm)%box_id

              ! if (particle_id.eq.43) then
              !    print *,"PREUPDATE  ",this%saps(new_sap_ids(nn))%x_endpoints(this%saps(new_sap_ids(nn))%boxes(box_id)%minendpoint_id(1))%value, &
              !         this%saps(new_sap_ids(nn))%x_endpoints(this%saps(new_sap_ids(nn))%boxes(box_id)%maxendpoint_id(1))%value, &
              !         "boxid=",box_id," endpoint_ids =",this%saps(new_sap_ids(nn))%boxes(box_id)%minendpoint_id(1),this%saps(new_sap_ids(nn))%boxes(box_id)%maxendpoint_id(1)
              ! endif

              ! asdf = 0.1
              ! diff = (DES_POS_NEW(1,particle_id)-DES_RADIUS(particle_id) - this%saps(new_sap_ids(nn))%x_endpoints(this%saps(new_sap_ids(nn))%boxes(ii)%minendpoint_id(1))%value)
              ! if (asdf < abs(diff) .and. abs(diff) < 1000000000.0) then
              !    print *,"  min xdiff ======== ",diff
              !    print *,"DES_POS_NEW(1) = ",DES_POS_NEW(1,particle_id)
              !    print *,"DES_RADIUS = ",DES_RADIUS(particle_id)
              !    print *,"value = ",this%saps(new_sap_ids(nn))%x_endpoints(this%saps(new_sap_ids(nn))%boxes(ii)%minendpoint_id(1))%value
              !    error stop __LINE__
              ! endif

              ! asdf = 0.1
              ! diff = (DES_POS_NEW(1,particle_id)-DES_RADIUS(particle_id) - this%saps(new_sap_ids(nn))%x_endpoints(this%saps(new_sap_ids(nn))%boxes(ii)%maxendpoint_id(1))%value)
              ! if (asdf < abs(diff) .and. abs(diff) < 1000000000.0) then
              !    print *,"  max xdiff ======== ",diff
              !    print *,"DES_POS_NEW(1) = ",DES_POS_NEW(1,particle_id)
              !    print *,"DES_RADIUS = ",DES_RADIUS(particle_id)
              !    print *,"value = ",this%saps(new_sap_ids(nn))%x_endpoints(this%saps(new_sap_ids(nn))%boxes(ii)%maxendpoint_id(1))%value
              !    error stop __LINE__
              ! endif


             call update_box(this%saps(new_sap_ids(nn)),handlelist%list(mm)%box_id,aabb)

             ! if (particle_id.eq.43) then
             !    print *,"POSTUPDATE  ",this%saps(new_sap_ids(nn))%x_endpoints(this%saps(new_sap_ids(nn))%boxes(box_id)%minendpoint_id(1))%value, &
             !         this%saps(new_sap_ids(nn))%x_endpoints(this%saps(new_sap_ids(nn))%boxes(box_id)%maxendpoint_id(1))%value, &
             !         "boxid=",box_id," endpoint_ids =",this%saps(new_sap_ids(nn))%boxes(box_id)%minendpoint_id(1),this%saps(new_sap_ids(nn))%boxes(box_id)%maxendpoint_id(1)
             ! endif

             ! ii = handlelist%list(mm)%box_id
             !  asdf = 0.1
             !  diff = (DES_POS_NEW(1,particle_id)-DES_RADIUS(particle_id) - this%saps(new_sap_ids(nn))%x_endpoints(this%saps(new_sap_ids(nn))%boxes(ii)%minendpoint_id(1))%value)
             !  if (asdf < abs(diff) .and. abs(diff) < 1000000000.0) then
             !     print *,"  min xdiff ======== ",diff
             !     print *,"DES_POS_NEW(1) = ",DES_POS_NEW(1,particle_id)
             !     print *,"DES_RADIUS = ",DES_RADIUS(particle_id)
             !     print *,"value = ",this%saps(new_sap_ids(nn))%x_endpoints(this%saps(new_sap_ids(nn))%boxes(ii)%minendpoint_id(1))%value
             !     error stop __LINE__
             !  endif

              ! asdf = 0.1
              ! diff = (DES_POS_NEW(1,particle_id)+DES_RADIUS(particle_id) - this%saps(new_sap_ids(nn))%x_endpoints(this%saps(new_sap_ids(nn))%boxes(ii)%maxendpoint_id(1))%value)
              ! if (asdf < abs(diff) .and. abs(diff) < 1000000000.0) then
              !    print *,"  max xdiff ======== ",diff
              !    print *,"DES_POS_NEW(1) = ",DES_POS_NEW(1,particle_id)
              !    print *,"DES_RADIUS = ",DES_RADIUS(particle_id)
              !    print *,"value = ",this%saps(new_sap_ids(nn))%x_endpoints(this%saps(new_sap_ids(nn))%boxes(ii)%maxendpoint_id(1))%value
              !    error stop __LINE__
              ! endif

             found = .true.
             exit
          endif
       enddo
       if (.not.found) then
          !if (particle_id.eq.43) then
             !print *,"ADDDDDDD  ",this%saps(new_sap_ids(nn))%x_endpoints(this%saps(new_sap_ids(nn))%boxes(handlelist%list(mm)%box_id)%maxendpoint_id(1))%value
          !endif

          ! add SAPs yet not listed for this id
          if (first_blank .eq. -1) then
             print *,"FAIL",first_blank,handlelist%list
             error stop __LINE__
          endif
          call add_box(this%saps(new_sap_ids(nn)),aabb,handlelist%particle_id,handlelist%list(first_blank)%box_id)
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
          ! if (particle_id.eq.43) then
             ! print *,"DELETED FROM SAP......",handlelist%list(mm)%sap_id
          ! endif
          call del_box(this%saps(handlelist%list(mm)%sap_id),handlelist%list(mm)%box_id)
          handlelist%list(mm)%sap_id = -1
       endif

    enddo

  end subroutine multisap_update

  subroutine multisap_sort(this)
    implicit none
    type(multisap_t), intent(inout) :: this
    integer :: ii, jj, kk, sap_id
    type(sap_t) :: sap
    type(box_t) :: box
    integer :: pair(2)
    logical :: asdf

    call init_pairs(this%hashtable)

!$omp parallel default(none) private(ii,jj,kk,sap_id,pair) shared(this)
!$omp do collapse(3)
    do ii=0,this%grid(1)-1
       do jj=0,this%grid(2)-1
          do kk=0,this%grid(3)-1
             sap_id = ii*this%grid(2)*this%grid(3)+jj*this%grid(3)+kk
             !print *,"NOW GOING TO SORT SAP::",ii,jj,kk,":::   ",ii*this%grid(2)*this%grid(3)+jj*this%grid(3)+kk
             if (.not. check_boxes(this%saps(sap_id))) error stop __LINE__
             call sort(this%saps(sap_id))
             if (.not.check_sort(this%saps(sap_id))) error stop __LINE__
             if (.not. check_boxes(this%saps(sap_id))) error stop __LINE__
          enddo
       enddo
    enddo
!$omp end parallel

    do ii=0,this%grid(1)-1
       do jj=0,this%grid(2)-1
          do kk=0,this%grid(3)-1
             sap_id = ii*this%grid(2)*this%grid(3)+jj*this%grid(3)+kk
             call reset_pairs(this%saps(sap_id)%hashtable)
             do
                call get_pair(this%saps(sap_id)%hashtable,pair)

                ! if (pair(1).eq. 114 .and. pair(2).eq.115) then
                !    print *,"GOTTA PAIR::::",pair, " IN SAP ",sap_id
                !    error stop __LINE__
                ! endif

                if (pair(1).eq.0 .and. pair(2).eq.0) exit
                call add_pair(this%hashtable,pair(1),pair(2))
             enddo
          enddo
       enddo
    enddo

  end subroutine multisap_sort

  subroutine multisap_quicksort(this)
    implicit none
    type(multisap_t), intent(inout) :: this
    integer :: ii, jj, kk
    integer :: sap_id, boxcount

    boxcount = 0
    print *,"grid ",this%grid(:)

    do ii=0,this%grid(1)-1
       do jj=0,this%grid(2)-1
          do kk=0,this%grid(3)-1
             sap_id = ii*this%grid(2)*this%grid(3)+jj*this%grid(3)+kk
             ! print *,"NOW GOING DO QUICKSORT THE THING:::::   ",sap_id
             call do_quicksort(this,this%saps(sap_id))
             boxcount = boxcount + this%saps(sap_id)%boxes_len
             !if (.not.check_boxes(multisap%saps(sap_id))) error stop __LINE__
             ! print *,"NOW GOING TO CHECK",this%saps(sap_id)%id
             if (.not.check_sort(this%saps(sap_id))) error stop __LINE__
          enddo
       enddo
    enddo

    print *,"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx TOTAL BOXES::::::::::::",boxcount

  end subroutine multisap_quicksort

  subroutine do_quicksort(this,sap)
    implicit none
    type(multisap_t), intent(inout) :: this
    type(sap_t), intent(inout) :: sap
    integer :: nn

    call quicksort_endpoints(sap%x_endpoints(1:sap%x_endpoints_len),sap,1)
    !if (.not.check_boxes(sap)) error stop __LINE__
    !if (.not.check_sort(sap)) error stop __LINE__

    call quicksort_endpoints(sap%y_endpoints(1:sap%y_endpoints_len),sap,2)
    !if (.not.check_boxes(sap)) error stop __LINE__
    !if (.not.check_sort(sap)) error stop __LINE__

    call quicksort_endpoints(sap%z_endpoints(1:sap%z_endpoints_len),sap,3)
    !if (.not.check_boxes(sap)) error stop __LINE__
    if (.not.check_sort(sap)) error stop __LINE__

  end subroutine do_quicksort

  subroutine multisap_sweep(this)
    ! use pair_manager
    implicit none
    type(multisap_t), intent(inout) :: this
    integer :: ii, jj, kk
    integer :: ll, i, sap_id
    integer :: minenx, mineny, minenz, minenx2, mineny2, minenz2
    integer :: maxenx, maxeny, maxenz, maxenx2, maxeny2, maxenz2

    do ii=0,this%grid(1)-1
       do jj=0,this%grid(2)-1
          do kk=0,this%grid(3)-1
             sap_id = (ii*this%grid(2)+jj)*this%grid(3)+kk
             ! print *,"NOW GOING TO sweep THE THING:::::   ",sap_id
             call sweep(this%saps(sap_id),sap_id)

             !if (.not.check_boxes( this%saps(ii*this%grid(2)*this%grid(3)+jj*this%grid(3)+kk) )) error stop __LINE__

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
