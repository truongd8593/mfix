      INTEGER       N_PLOT_MENU
      PARAMETER     (N_PLOT_MENU = 7)
C
      INTEGER       N_VAR_MENU
      PARAMETER     (N_VAR_MENU = 13)
C
C                   Flag to tell whether an SPX file is open
      LOGICAL       SPX_OPEN(N_SPX)
C
C                   Last record in an SPX file, number of files
      INTEGER       LAST_REC(N_SPX), NUM_REC(N_SPX)
C
C                   index of variable within the file to plot
      INTEGER       VAR_INDEX
      INTEGER       LOC_X,LOC_Y,LOC_Z
      INTEGER       PLOT_TYPE,PLOT_TYPE2
      INTEGER       SELECTION
      INTEGER       COUNT 
      CHARACTER*60  PLOT_MENU(N_PLOT_MENU)
      CHARACTER*60  VAR_MENU(N_VAR_MENU)
      REAL          VERSION_NUMBER
      REAL          XDIST_SC(DIMENSION_I) , XDIST_VEC(DIMENSION_I)
      REAL          YDIST_SC(DIMENSION_J) , YDIST_VEC(DIMENSION_J)
      REAL          ZDIST_SC(DIMENSION_K) , ZDIST_VEC(DIMENSION_K)
      REAL          VOL_SC(DIMENSION_3) , VOL_U(DIMENSION_3)
      REAL          VOL_V(DIMENSION_3)  , VOL_W(DIMENSION_3)
C
      COMMON /POST3D_I/ VAR_INDEX, LOC_X, LOC_Y, LOC_Z, PLOT_TYPE,
     &                 SELECTION, COUNT, PLOT_TYPE2, LAST_REC, NUM_REC
C
      COMMON /POST3D_C/ PLOT_MENU, VAR_MENU
C
      COMMON /POST3D_L/ SPX_OPEN
C
      COMMON /POST3D_R/ XDIST_SC,XDIST_VEC,YDIST_SC,YDIST_VEC,
     &                  ZDIST_SC,ZDIST_VEC,VOL_SC,VOL_U,
     &                  VOL_V,VOL_W, VERSION_NUMBER
