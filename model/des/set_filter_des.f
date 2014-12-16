!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SET_FILTER_DES                                   !
!  Author: J.Musser                                   Date: 25-Nov-14  !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE SET_FILTER_DES

! Runtime Flag: Utilize cutcell geometry.
      use cutcell, only: CARTESIAN_GRID
! Runtime Flag: Invoke gas/solids coupled simulation.
      use discretelement, only: DES_CONTINUUM_COUPLED

      use geometry, only: DX, IMIN1, IMAX1
      use geometry, only: DY, JMIN1, JMAX1
      use geometry, only: DZ, KMIN1, KMAX1, DO_K

      use particle_filter, only: DES_INTERP_SCHEME_ENUM
      use particle_filter, only: DES_INTERP_NONE
      use particle_filter, only: DES_INTERP_GARG
      use particle_filter, only: DES_INTERP_DPVM
      use particle_filter, only: DES_INTERP_GAUSS
      use particle_filter, only: DES_INTERP_SCHEME

      use particle_filter, only: DIF_TSTOP
      use particle_filter, only: FILTER_WIDTH
      use particle_filter, only: FILTER_WIDTH_INTERP
      use particle_filter, only: DES_DIFFUSE_MEAN_FIELDS
      use particle_filter, only: DES_INTERP_MEAN_FIELDS
      use particle_filter, only: DES_INTERP_ON

      use particle_filter, only: OoFILTER_VOL
      use particle_filter, only: FILTER_WIDTH_INTERPx3

      use param1, only: UNDEFINED, UNDEFINED_C

      use error_manager

      IMPLICIT NONE

      DOUBLE PRECISION :: DXYZ_MIN


!......................................................................!

! Initialize the error manager.
      CALL INIT_ERR_MSG("SET_FILTER_DES")

! Verify that the interpolation scheme doesn't exceed the grid.
      IF(FILTER_WIDTH /= UNDEFINED) THEN

         DXYZ_MIN = min(minval(DX(IMIN1:IMAX1)),minval(DY(JMIN1:JMAX1)))
         IF(DO_K) DXYZ_MIN = min(DXYZ_MIN,minval(DZ(KMIN1:KMAX1)))

         IF(0.5d0*FILTER_WIDTH < DXYZ_MIN) THEN
             FILTER_WIDTH_INTERP = 0.500d0*FILTER_WIDTH
         ELSEIF(0.5d0*FILTER_WIDTH == DXYZ_MIN) THEN
             FILTER_WIDTH_INTERP = 0.999d0*DXYZ_MIN
         ELSEIF(DES_DIFFUSE_MEAN_FIELDS) THEN
             FILTER_WIDTH_INTERP = 0.999d0*DXYZ_MIN
         ELSE
            WRITE(ERR_MSG,2130) DXYZ_MIN, FILTER_WIDTH
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

 2130 FORMAT('Error 2130: The specified FILTER_WIDTH is too large to ',&
         'interpolate',/'to the fuild grid as it spans further than ', &
         'adjacent neighbors. Either',/'discrease the filter width, ', &
         'or enable diffusive filtering.',2/3x,'Minimum Cell ',        &
         'dimension: ',g12.4,/3x,'Filter Width: ',g12.4)


         IF(DES_DIFFUSE_MEAN_FIELDS) THEN

            DES_INTERP_MEAN_FIELDS= .TRUE.

            DIF_TSTOP = MAX((0.5*FILTER_WIDTH)**2-DXYZ_MIN**2,0.0d0) / &
               16.0d0*log(2.0)

            IF(DIF_TSTOP == 0.0) THEN
               WRITE(ERR_MSG,2131)
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF

         ENDIF

 2131 FORMAT('Error 2131: The FILTER_WIDTH is too small for ',       &
         'diffusive filtering ',/'on the current mesh and will ',    &
         'have no effect. Either increase the',/'FILTER_WIDTH or ',  &
         'set DES_DIFFUSE_MEAN_FIELDS=.FALSE.')

      ENDIF

! Calculate reused quanties
      SELECT CASE(DES_INTERP_SCHEME_ENUM)
      CASE(DES_INTERP_DPVM, DES_INTERP_GAUSS)
         OoFILTER_VOL = 0.25d0/(FILTER_WIDTH_INTERP**3)
         FILTER_WIDTH_INTERPx3 = FILTER_WIDTH_INTERP*3
      END SELECT

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE SET_FILTER_DES
