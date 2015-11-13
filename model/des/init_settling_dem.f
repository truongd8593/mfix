!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Module name: MAKE_ARRAYS_DES                                        !
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  !
!                                                                      !
!  Purpose: DES - allocating DES arrays
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE INIT_SETTLING_DEM

      USE desgrid, ONLY: desgrid_pic
      USE discretelement
      USE error_manager
      USE mpi_funs_des, ONLY: DES_PAR_EXCHANGE
      USE run
      use multi_sweep_and_prune
      use functions, only: is_nonexistent

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER :: FACTOR, nn

      type(aabb_t) :: aabb

!-----------------------------------------------
! Include statement functions
!-----------------------------------------------


! Skip this routine if there are no particles.
      IF(PARTICLES == 0) RETURN
! Skip this routine if not a new run.
      IF(RUN_TYPE /= 'NEW') RETURN

! Skip if not coupled.
      IF(.NOT.DES_CONTINUUM_COUPLED) RETURN

! Write the initial configuration before settling
      IF(PRINT_DES_DATA .AND. NFACTOR>0) CALL WRITE_DES_DATA

      WRITE(ERR_MSG, 1100) trim(iVal(NFACTOR))
      CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
 1100 FORMAT('Beginning DEM settling period: ',A,' steps.')

! Disable the coupling flag.
      DES_CONTINUUM_COUPLED = .FALSE.

      ! initialize SAP
      do nn=1, MAX_PIP
         if(is_nonexistent(nn)) cycle
         aabb%minendpoint(:) = DES_POS_NEW(:,nn)-DES_RADIUS(nn)-0.001
         aabb%maxendpoint(:) = DES_POS_NEW(:,nn)+DES_RADIUS(nn)+0.001

         if ( any(DES_RADIUS(nn)*multisap%one_over_cell_length(:) > 0.5 ) ) then
            print *,"BAD RAD...grid too fine, need to have radius=",des_radius(nn),"  less than ",0.5/multisap%one_over_cell_length(:)
            stop __LINE__
         endif

         ! if (0.eq.mod(nn,1000)) print *,"PARTICLE #  ",nn
         ! print *,""
         ! print *,"aabb%minendpoint(:) = ",aabb%minendpoint(:)
         ! print *,"aabb%maxendpoint(:) = ",aabb%maxendpoint(:)
         ! print *,"DES_POS_NEW(:,nn) = ",DES_POS_NEW(:,nn)
         ! print *,""

         call multisap_add(multisap,aabb,nn,boxhandle(nn))
      enddo

      call multisap_quicksort(multisap)

      do nn=0,size(multisap%saps)-1
         if (.not.check_sort(multisap%saps(nn))) stop __LINE__
      enddo

      call multisap_sweep(multisap)

      do nn=0,size(multisap%saps)-1
         if (.not.check_boxes(multisap%saps(nn))) stop __LINE__
         if (.not.check_sort(multisap%saps(nn))) stop __LINE__
      enddo

      DO FACTOR = 1, NFACTOR
         print *,"FACTOR  =  ",FACTOR
! calculate forces

         do nn=0,size(multisap%saps)-1
            if (.not.check_boxes(multisap%saps(nn))) stop __LINE__
            if (.not.check_sort(multisap%saps(nn))) stop __LINE__
         enddo

         CALL CALC_FORCE_DEM
! update particle position/velocity

         do nn=0,size(multisap%saps)-1
            if (.not.check_boxes(multisap%saps(nn))) stop __LINE__
            if (.not.check_sort(multisap%saps(nn))) stop __LINE__
         enddo

         CALL CFNEWVALUES
! set the flag do_nsearch before calling particle in cell (for mpi)
         DO_NSEARCH = (MOD(FACTOR,NEIGHBOR_SEARCH_N)==0)

         do nn=0,size(multisap%saps)-1
            if (.not.check_boxes(multisap%saps(nn))) stop __LINE__
            if (.not.check_sort(multisap%saps(nn))) stop __LINE__
         enddo

! Bin the particles to the DES grid.
         CALL DESGRID_PIC(.TRUE.)
! exchange particle crossing boundaries and updates ghost particles
         CALL DES_PAR_EXCHANGE
! find particles on grid
         CALL PARTICLES_IN_CELL
! perform neighbor search
         IF(DO_NSEARCH) CALL NEIGHBOUR
      ENDDO

! Reset the comoupling flag.
      DES_CONTINUUM_COUPLED = .TRUE.

      WRITE(ERR_MSG, 1200)
      CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
 1200 FORMAT('DEM settling period complete.')

! this write_des_data is needed to properly show the initial state of
! the simulation (granular or coupled). In the coupled case, the
! particles may have 'settled' according to above.  In the granular
! case, the initial state won't be written until after the particles
! have moved without this call.
!      IF(PRINT_DES_DATA) CALL WRITE_DES_DATA

      RETURN
      END SUBROUTINE INIT_SETTLING_DEM

