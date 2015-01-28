!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!     Module name: DEPRECATED_OR_UNKNOWN                               !
!     Author: J.Musser                                Date:  5-SEPT-14 !
!                                                                      !
!     Purpose: This routine is called when a keyword was not matched   !
!     to any of the keywords in the namelist files. This routine       !
!     reports if the keyword was deprecated or incorrect.              !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DEPRECATED_OR_UNKNOWN(LINE_NO, INPUT)

      use param
      use param1
      use compar, only: myPE, PE_IO
      use error_manager

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: LINE_NO
      CHARACTER(len=*), INTENT(IN) :: INPUT

      CHARACTER(len=256) :: STRING

! Old keyword for solids density :: Replaced by RO_s0
      DOUBLE PRECISION :: RO_s
      LOGICAL :: BC_APPLY_TO_MPPIC(DIMENSION_BC)
      INTEGER :: COHESION_DEBUG, DES_MMAX,DES_NMAX_s(DIM_M),DIMN
      CHARACTER(LEN=16) :: DES_BC_TYPE(DIMENSION_BC)
      DOUBLE PRECISION :: DES_BC_MASSFLOW_s(DIMENSION_BC, DIM_M)
      DOUBLE PRECISION :: DES_BC_ROP_s (DIMENSION_BC, DIM_M)
      DOUBLE PRECISION :: DES_BC_T_s (DIMENSION_IC, DIM_M)
      DOUBLE PRECISION :: DES_BC_VOLFLOW_s(DIMENSION_BC, DIM_M)
      DOUBLE PRECISION :: DES_BC_X_e(DIMENSION_BC)
      DOUBLE PRECISION :: DES_BC_X_w(DIMENSION_BC)
      DOUBLE PRECISION :: DES_BC_Y_n(DIMENSION_BC)
      DOUBLE PRECISION :: DES_BC_Y_s(DIMENSION_BC)
      DOUBLE PRECISION :: DES_BC_Z_b(DIMENSION_BC)
      DOUBLE PRECISION :: DES_BC_Z_t(DIMENSION_BC)
      DOUBLE PRECISION :: DES_C_ps0(DIM_M), DES_K_s0(DIM_M)
      LOGICAL :: DES_CALC_BEDHEIGHT, DES_CONV_EQ, DES_ENERGY_EQ
      LOGICAL :: DES_COND_EQ,DES_COND_EQ_PFP,DES_COND_EQ_PP, DES_RADI_EQ

      DOUBLE PRECISION :: DES_D_P0 (DIM_M), DES_F, DES_GAMMA
      DOUBLE PRECISION :: DES_EPS_XSTART,DES_EPS_YSTART,DES_EPS_ZSTART

      DOUBLE PRECISION :: DES_IC_X_e(DIMENSION_IC)
      DOUBLE PRECISION :: DES_IC_X_w(DIMENSION_IC)
      DOUBLE PRECISION :: DES_IC_Y_n(DIMENSION_IC)
      DOUBLE PRECISION :: DES_IC_Y_s(DIMENSION_IC)
      DOUBLE PRECISION :: DES_IC_Z_b(DIMENSION_IC)
      DOUBLE PRECISION :: DES_IC_Z_t(DIMENSION_IC)
      DOUBLE PRECISION :: DES_IC_T_s (DIMENSION_IC, DIM_M)
      DOUBLE PRECISION :: DES_MW_s(DIM_M, DIM_N_s), DES_RO_s (DIM_M)
      DOUBLE PRECISION :: DTSOLID_FACTOR, LID_VEL
      DOUBLE PRECISION :: MASTER_WELL_DEPTH, MASTER_WALL_WELL_DEPTH
      CHARACTER(len=18) DES_SPECIES_s(DIM_M, DIM_N_s)
      CHARACTER(len=32)  DES_SPECIES_ALIAS_s(DIM_M, DIM_N_s)
      LOGICAL :: DES_PERIODIC_WALLS, DES_PERIODIC_WALLS_X,DES_SPECIES_EQ
      LOGICAL :: DES_PERIODIC_WALLS_Y, DES_PERIODIC_WALLS_Z,TSUJI_DRAG

      INTEGER :: MAX_DES_BC_CELL, NPC_PIC(DIM_M)
      INTEGER :: QLM, QLN, INIT_QUAD_COUNT
      LOGICAL :: MPPIC_CONSTANTNPC, MPPIC_CONSTANTWT, SQUARE_WELL
      LOGICAL :: WALLDTSPLIT,WALLREFLECT, CALL_DI, CALL_GROW, CALL_ISAT
      DOUBLE PRECISION :: MQUAD_FACTOR,STATWT_PIC(DIM_M), ISATdt
      DOUBLE PRECISION :: RADIUS_RATIO,WALL_RADIUS_RATIO
      DOUBLE PRECISION :: pvel_mean, PVEL_StDev,VOL_FRAC(DIM_M)
      CHARACTER(LEN=64) :: REACTION_MODEL

      LOGICAL :: DISCRETE_ELEMENT, MPPIC, DES_CONTINUUM_HYBRID

! These are deprecated keywords. Some are no longer in the code and
! others are used, but no longer specifed in the mfix.dat file.
      NAMELIST / DEPRECATED / RO_s, BC_APPLY_TO_MPPIC, COHESION_DEBUG,&
         DES_BC_MASSFLOW_s, DES_BC_ROP_s,DES_BC_T_s,DES_BC_TYPE,       &
         DES_BC_VOLFLOW_s,DES_BC_X_e,DES_BC_X_w, DES_BC_Y_n,DES_BC_Y_s,&
         DES_BC_Z_b,DES_BC_Z_t,DES_CALC_BEDHEIGHT,QLM,QLN,DES_COND_EQ, &
         DES_COND_EQ_PFP,DES_COND_EQ_PP,DES_CONV_EQ,DES_C_ps0,DES_D_p0,&
         DES_ENERGY_EQ,DES_EPS_XSTART,DES_F,DES_EPS_YSTART,DIMN,       &
         DES_EPS_ZSTART,DES_GAMMA,DES_IC_T_s,DES_IC_X_e,DES_IC_X_w,    &
         DES_IC_Y_n,DES_IC_Y_s,DES_IC_Z_b,DES_IC_Z_t,DES_K_s0,DES_MMAX,&
         DES_MW_s,DES_NMAX_s,DES_PERIODIC_WALLS,DES_SPECIES_s,LID_VEL, &
         DES_PERIODIC_WALLS_X,DES_PERIODIC_WALLS_Y,WALLREFLECT,        &
         DES_PERIODIC_WALLS_Z,DES_RADI_EQ,DES_RO_s,DES_SPECIES_ALIAS_s,&
         DES_SPECIES_EQ,DTSOLID_FACTOR,INIT_QUAD_COUNT,MAX_DES_BC_CELL,&
         MASTER_WALL_WELL_DEPTH,MASTER_WELL_DEPTH,MPPIC_CONSTANTNPC,   &
         MPPIC_CONSTANTWT,MQUAD_FACTOR,NPC_PIC,pvel_mean,pvel_stdev,   &
         RADIUS_RATIO,REACTION_MODEL,SQUARE_WELL,STATWT_PIC,STATWT_PIC,&
         TSUJI_DRAG,VOL_FRAC,WALLDTSPLIT,WALL_RADIUS_RATIO, CALL_DI,   &
         CALL_GROW, CALL_ISAT, ISATdt, DISCRETE_ELEMENT, MPPIC,        &
         DES_CONTINUUM_HYBRID


      STRING = '&DEPRECATED '//trim(adjustl(INPUT))//'/'
      READ(STRING,NML=DEPRECATED,ERR=999, END=999)

      IF(myPE == 0) WRITE(*,2000) trim(iVAL(LINE_NO)), trim(INPUT)
      CALL MFIX_EXIT(myPE)

 2000 FORMAT(//1X,70('*')/' From DEPRECATED_OR_UNKNOWN:',/1x,          &
         'Error 2000: A keyword pair on line ',A,' of the mfix.dat ',  &
         'file was',/' identified as being deprecated.',//3x,A,//1x,   &
         'Please see the user documentation and update the mfix.dat ', &
         'file.',/1X,70('*')//)

  999 if(myPE == 0)WRITE (*, 2001) trim(iVAL(LINE_NO)), trim(INPUT)
      CALL MFIX_EXIT(myPE)


 2001 FORMAT(//1X,70('*')/' From: DEPRECATED_OR_UNKNOWN',/1x,          &
         'Error 2001: Unable to process line ',A,' of the mfix.dat ',  &
         'file.',2/3x,A,2/1x,'Possible causes are',/3x,'*  Incorrect', &
         ' or illegal keyword format',/3x,'*  Unknown or mistyped name',&
         /3x,'*  The dimensioned item is too small (array overflow).', &
         2/1x,'Please see the user documentation and update the ',     &
         'mfix.dat file. ',/1X,70('*')//)

      END SUBROUTINE DEPRECATED_OR_UNKNOWN
