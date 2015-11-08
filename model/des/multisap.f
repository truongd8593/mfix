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

  subroutine multisap_add(this,aabb,particle_id,handle)
    implicit none
    type(multisap_t), intent(inout) :: this
    type(aabb_t), intent(in) :: aabb
    integer, intent(in) :: particle_id
    type(boxhandle_t), intent(out) :: handle
    integer :: sap_ids(MAX_SAPS)
    integer :: nn,id

    call multisap_raster(this,aabb,sap_ids,particle_id)

    if (particle_id.eq.46) then
       print *,"RASTERED 46 TO:  ",sap_ids
    endif

    if (particle_id.eq.77) then
       print *,"RASTERED 77 TO:  ",sap_ids
    endif


    ! add to each individual SAP
    do nn=1, size(sap_ids)
       handle%sap_id(nn) = sap_ids(nn)
       if (sap_ids(nn) >= 0) then
          call add_box(this%saps(sap_ids(nn)),aabb,particle_id,handle%box_id(nn))
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
    integer :: sap_ids(MAX_SAPS)
    logical :: found
    integer :: mm,nn

    call multisap_raster(this,aabb,sap_ids,-1)

    ! update for each SAP listed for this id
    do nn=1, size(handle%sap_id)
       if (handle%sap_id(nn) < 0) exit

       found = .false.
       do mm=1, size(sap_ids)
          if (sap_ids(mm) < 0) exit

          !print *,"size = ",size(this%sap_ids,1),,size(this%sap_ids,2)
          if (sap_ids(mm) .eq. handle%sap_id(nn)) then
             ! update SAPs listed for this id
             !print *,"sap_ids(nn) == ",sap_ids(nn)
             !print *,"this%sap_ids(nn,id) == ",this%sap_ids(mm,id)
             call update_box(this%saps(handle%sap_id(nn)),handle%box_id(nn),aabb)
             found = .true.
             exit
          endif
       enddo
       if (.not.found) then
          ! add to SAPs not listed for this id
          call del_box(this%saps(handle%sap_id(nn)),handle%box_id(nn))
       endif
    enddo

    do mm=1, size(sap_ids)
       if (sap_ids(mm) < 0) exit

       found = .false.
       do nn=1, size(handle%sap_id)
          if (handle%sap_id(nn) < 0) exit

          if (sap_ids(mm) .eq. handle%sap_id(nn)) then
             found = .true.
             exit
          endif
       enddo
       if (.not.found) then
          ! remove from SAPs not listed for this id
          ! print *,"HANDLE: ",handle%sap_id,handle%box_id
          ! print *," IS NOT FOUND IN SAP ",mm
          call del_box(this%saps(handle%sap_id(nn)),handle%box_id(nn))
       endif
    enddo

  end subroutine multisap_update

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
             !print *,"NOW GOING TO SORT THE THING::",ii,jj,kk,":::   ",ii*this%grid(2)*this%grid(3)+jj*this%grid(3)+kk
             call sort(this%saps(ii*this%grid(2)*this%grid(3)+jj*this%grid(3)+kk))
             if (.not. check_boxes(this%saps(ii*this%grid(2)*this%grid(3)+jj*this%grid(3)+kk))) stop __LINE__
          enddo
       enddo
    enddo

    ! do ii = 1, MAX_PIP
    !    do jj = 1, MAX_PIP
    !       sap = multisap%saps(boxhandle(ii)%sap_id(jj))
    !       box = sap%boxes(boxhandle(ii)%box_id(jj))
    !       if (sap%x_endpoints(box%minendpoint_id(1))%box_id .ne. boxhandle(ii)%box_id(jj)) stop __LINE__
    !       if (sap%x_endpoints(box%minendpoint_id(1))%value > 100) stop __LINE__
    !    enddo
    ! enddo

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

    ll = 46
    i = 77

    do ss = 10000, MAX_SAPS

       print *,"SWEEP SS = ",ss

               minenx = multisap%saps(boxhandle(ll)%sap_id(ss))%boxes(boxhandle(ll)%box_id(1))%minendpoint_id(1)
               maxenx = multisap%saps(boxhandle(ll)%sap_id(ss))%boxes(boxhandle(ll)%box_id(1))%maxendpoint_id(1)
               print *,"PARTICLE (",ll,"):  ",des_pos_new(:,ll)," at x endpoints ",minenx," AND ",maxenx,multisap%saps(boxhandle(ll)%sap_id(ss))%x_endpoints(minenx),multisap%saps(boxhandle(ll)%sap_id(ss))%x_endpoints(maxenx)

               minenx2 = multisap%saps(boxhandle(i)%sap_id(ss))%boxes(boxhandle(i)%box_id(1))%minendpoint_id(1)
               maxenx2 = multisap%saps(boxhandle(i)%sap_id(ss))%boxes(boxhandle(i)%box_id(1))%maxendpoint_id(1)
               print *,"PARTICLE (",i,"):  ",des_pos_new(:,i)," at x endpoints ",minenx2," AND ",maxenx2,multisap%saps(boxhandle(i)%sap_id(ss))%x_endpoints(minenx2),multisap%saps(boxhandle(i)%sap_id(ss))%x_endpoints(maxenx2)

               mineny = multisap%saps(boxhandle(ll)%sap_id(ss))%boxes(boxhandle(ll)%box_id(2))%minendpoint_id(2)
               maxeny = multisap%saps(boxhandle(ll)%sap_id(ss))%boxes(boxhandle(ll)%box_id(2))%maxendpoint_id(2)
               print *,"PARTICLE (",ll,"):  ",des_pos_new(:,ll)," at y endpoints ",mineny," AND ",maxeny,multisap%saps(boxhandle(ll)%sap_id(ss))%y_endpoints(mineny),multisap%saps(boxhandle(ll)%sap_id(ss))%y_endpoints(maxeny)

               mineny2 = multisap%saps(boxhandle(i)%sap_id(ss))%boxes(boxhandle(i)%box_id(2))%minendpoint_id(2)
               maxeny2 = multisap%saps(boxhandle(i)%sap_id(ss))%boxes(boxhandle(i)%box_id(2))%maxendpoint_id(2)
               print *,"PARTICLE (",i,"):  ",des_pos_new(:,i)," at y endpoints ",mineny2," AND ",maxeny2,multisap%saps(boxhandle(i)%sap_id(ss))%y_endpoints(mineny2),multisap%saps(boxhandle(i)%sap_id(ss))%y_endpoints(maxeny2)

               minenz = multisap%saps(boxhandle(ll)%sap_id(ss))%boxes(boxhandle(ll)%box_id(3))%minendpoint_id(3)
               maxenz = multisap%saps(boxhandle(ll)%sap_id(ss))%boxes(boxhandle(ll)%box_id(3))%maxendpoint_id(3)
               print *,"PARTICLE (",ll,"):  ",des_pos_new(:,ll)," at z endpoints ",minenz," AND ",maxenz,multisap%saps(boxhandle(ll)%sap_id(ss))%z_endpoints(minenz),multisap%saps(boxhandle(ll)%sap_id(ss))%z_endpoints(maxenz)

               minenz2 = multisap%saps(boxhandle(i)%sap_id(ss))%boxes(boxhandle(i)%box_id(3))%minendpoint_id(3)
               maxenz2 = multisap%saps(boxhandle(i)%sap_id(ss))%boxes(boxhandle(i)%box_id(3))%maxendpoint_id(3)
               print *,"PARTICLE (",i,"):  ",des_pos_new(:,i), " at z endpoints ",minenz2," AND ",maxenz2,multisap%saps(boxhandle(i)%sap_id(ss))%z_endpoints(minenz2),multisap%saps(boxhandle(i)%sap_id(ss))%z_endpoints(maxenz2)

               if (max(minenx,minenx2) < min(maxenx,maxenx2)) then
                  print *,"X INTERSECTION "
                  !stop __LINE__
                  if (max(mineny,mineny2) < min(maxenx,maxeny2)) then
                     print *,"X AND Y INTERSECTION "

                     print *,"Z INTERSECT? ",max(minenz,minenz2),min(maxenz,maxenz2)
                     print *,"Z INTERSECT? ",(max(minenz,minenz2) < min(maxenz,maxenz2))
                     if (.not.is_pair(ll,i)) then
                        print *,"FAIL"

                        print *,"LL sap_id = ",boxhandle(LL)%sap_id(1)
                        print *,"I sap_id = ",boxhandle(I)%sap_id(1)

                        stop __LINE__
                     endif

                     stop __LINE__
                     if (max(minenz,minenz2) < min(maxenx,maxenz2)) then
                        print *,"INTERSECTION...."
                        if (.not.is_pair(ll,i)) then
                           print *,"FAIL"
                           stop __LINE__
                        endif
                        if (.not.is_pair(i,ll)) then
                           print *,"fail"
                           stop __LINE__
                        endif
                        stop __LINE__
                     endif
                  endif
               endif

            enddo

  end subroutine multisap_sweep

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
