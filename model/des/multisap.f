module multisap

  type multisap_t
     type(sap_t), dimension(:,:), allocatable :: saps
     integer :: saps_len
  end type multisap_t

contains

  subroutine init_multisap(this,y_grid,z_grid,y_min,y_max,z_min,z_max)
    implicit none
    type(multisap_t), intent(inout) :: this
    integer, intent(in) :: y_grid, z_grid
    real, intent(in) :: y_min, z_min, y_max, z_max

    allocate(saps(y_grid,z_grid))
    do ii=1,y_grid
       do jj=1,z_grid
          init_sap(saps(ii,jj))
       enddo
    enddo

  end subroutine init_multisap

end module multisap
