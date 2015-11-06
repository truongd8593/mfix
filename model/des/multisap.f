module multisap

  type multisap_t
     type(sap_t), dimension(:,:), allocatable :: saps
     integer :: saps_len
     real(3) :: minbounds, maxbounds
     real :: one_over_ycell_length, one_over_zcell_length
  end type multisap_t

contains

  subroutine init_multisap(this,y_grid,z_grid,minbounds,maxbounds)
    implicit none
    type(multisap_t), intent(inout) :: this
    integer, intent(in) :: y_grid, z_grid
    real, dimension(3), intent(in) :: minbounds, maxbounds

    allocate(saps(y_grid,z_grid))
    do ii=1,y_grid
       do jj=1,z_grid
          init_sap(saps(ii,jj))
       enddo
    enddo

    one_over_ycell_length = y_grid/(maxbounds(2)-minbounds(2))
    one_over_zcell_length = z_grid/(maxbounds(2)-minbounds(2))

  end subroutine init_multisap

  subroutine multisap_add(this,aabb,id)
    implicit none
    type(multisap_t), intent(inout) :: this
    type(aabb_t), intent(in) :: aabb
    integer, intent(out) :: id

    ! bin to y/z grid
    yy = floor((aabb%minendpoint(2)-minbounds(2))*one_over_ycell_length)
    zz = floor((aabb%minendpoint(3)-minbounds(3))*one_over_ycell_length)

    y2 = floor((aabb%maxendpoint(2)-minbounds(2))*one_over_ycell_length)
    z2 = floor((aabb%maxendpoint(3)-minbounds(3))*one_over_ycell_length)

    call sap_add(saps(yy,zz),aabb,id)

    if (y2.ne.yy) then
       call sap_add(saps(y2,zz),aabb,id)
    endif

    if (z2.ne.zz) then
       call sap_add(saps(yy,z2),aabb,id)
    endif

    yy = floor((aabb%maxendpoint(2)-minbounds(2))*one_over_ycell_length)
    zz = floor((aabb%minendpoint(3)-minbounds(3))*one_over_ycell_length)

    yy = floor((aabb%minendpoint(2)-minbounds(2))*one_over_ycell_length)
    zz = floor((aabb%maxendpoint(3)-minbounds(3))*one_over_ycell_length)

    ! add to each individual SAP

  end subroutine multisap_add

  subroutine multisap_del(this,aabb,id)
    implicit none
    type(multisap_t), intent(inout) :: this
    type(aabb_t), intent(in) :: aabb
    integer, intent(out) :: id

    ! bin to y/z grid
    ! remove from each SAP listed for this id
  end subroutine multisap_del

  subroutine multisap_update(this,aabb,id)
    implicit none
    type(multisap_t), intent(inout) :: this
    type(aabb_t), intent(in) :: aabb
    integer, intent(out) :: id

    ! bin to y/z grid
    ! add to SAPs not listed for this id
    ! update SAPs listed for this id
    ! remove from SAPs not listed for this id
  end subroutine multisap_update

end module multisap
