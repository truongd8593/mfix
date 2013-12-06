!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CARTESIAN_GRID_NAMELIST.INC                            C
!  Purpose: Variables needed for CARTESIAN_GRID input                  C
!                                                                      C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 26-Aug-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
       NAMELIST / CARTESIAN_GRID_INPUT_DATA /&

      CARTESIAN_GRID , N_QUADRIC,USE_STL,USE_MSH,USE_POLYGON,N_USR_DEF,&
      quadric_form, quadric_scale, &
      lambda_x, lambda_y, lambda_z, dquadric, &
      theta_x, theta_y, theta_z, &
      Radius, Half_angle,n_x,n_y,n_z, &
      t_x, t_y, t_z, &
      clip_xmin, clip_xmax, clip_ymin, clip_ymax, clip_zmin, clip_zmax, &
      piece_xmin, piece_xmax, piece_ymin, piece_ymax, piece_zmin, piece_zmax, &
      FLUID_IN_CLIPPED_REGION, &
      BC_ID_Q,&
      N_GROUP, GROUP_SIZE,GROUP_Q,GROUP_RELATION,RELATION_WITH_PREVIOUS, &
      TOL_SNAP,TOL_DELH, TOL_SMALL_CELL,TOL_MERGE,&
      TOL_SMALL_AREA, ALPHA_MAX, &
      TOL_F, TOL_POLY,ITERMAX_INT, &
      TOL_STL,STL_SMALL_ANGLE,RAY_DIR,OUT_STL_VALUE,STL_BC_ID,CAD_PROPAGATE_ORDER,&
      TX_STL,TY_STL,TZ_STL,SCALE_STL,&
      TOL_MSH,MSH_SMALL_ANGLE,OUT_MSH_VALUE,&
      TX_MSH,TY_MSH,TZ_MSH,SCALE_MSH,&
      SET_CORNER_CELLS,&
      FAC_DIM_MAX_CUT_CELL, &
      WRITE_VTK_FILES,TIME_DEPENDENT_FILENAME, VTK_DT, VTK_VAR,FRAME,VTU_DIR,&
      PG_OPTION,&
      CG_SAFE_MODE,PRINT_WARNINGS,&
      CG_UR_FAC,&
      PRINT_PROGRESS_BAR,BAR_WIDTH, BAR_CHAR, BAR_RESOLUTION,  &
      WRITE_DASHBOARD, F_DASHBOARD, &
      CPX,CPY,CPZ,NCX,NCY,NCZ,ERX,ERY,ERZ, &
  FIRST_DX,LAST_DX,FIRST_DY,LAST_DY,FIRST_DZ,LAST_DZ, &
  MAX_FACETS_PER_CELL_DES


  


