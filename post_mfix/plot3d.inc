C                maximum number of files in the plotting menu
      INTEGER    MAX_FILES
      PARAMETER  (MAX_FILES=30)
C
C                names of the menu files
      CHARACTER  MENU_FILES(MAX_FILES)*60
C
C                number of files in plot3d.control
      INTEGER    NUM_FILES
C
C                flag .. if (.true.) ... the menu file's format is OK
      LOGICAL    FORMAT_OK(MAX_FILES)
C
C                index of the file to plot
      INTEGER    FILE1_INDEX , FILE2_INDEX
C
C                index of variable within the file to plot
      INTEGER    VAR1_INDEX , VAR2_INDEX
C
      INTEGER    KCON
      INTEGER    IMIN_PLOT , IMAX_PLOT , JMIN_PLOT , JMAX_PLOT
      INTEGER    KMIN_PLOT , KMAX_PLOT
C
      INTEGER    GID(25) , NSPEED1 , NSPEED2 , MENU_NUM
C
      INTEGER    MENUX1 , MENUX2 , MENUY1 , MENUY2
C
      REAL       PLOT_MIN,PLOT_MAX,var1_min_color,var1_max_color
      REAL       var1_factor, ZOOM_X, ZOOM_Y, VAR1MIN,VAR1MAX
      REAL       VERT1_LINE
      REAL       VAR2_MIN_COLOR,VAR2_MAX_COLOR,VAR2_FACTOR
      REAL       VAR2MIN , VAR2MAX , VERT2_LINE
C
      LOGICAL DO_3D , FREEZE , DO_READ , DOING2 , CHANGE1
C
      COMMON /PLOT3D_C/ MENU_FILES
C
      COMMON /PLOT3D_I/ NUM_FILES, FILE1_INDEX, VAR1_INDEX, VAR2_INDEX,
     &                  GID , MENUX1 , MENUX2 , MENUY1 , MENUY2,
     &                  NSPEED1 , KCON , IMIN_PLOT , IMAX_PLOT ,
     &                  JMIN_PLOT , JMAX_PLOT , KMIN_PLOT ,
     &                  KMAX_PLOT , MENU_NUM , FILE2_INDEX ,
     &                  NSPEED2
C
      COMMON /PLOT3D_L/ FORMAT_OK , DO_3D , FREEZE , DO_READ ,
     &                  DOING2 , CHANGE1
C
      COMMON /PLOT3D_R/ PLOT_MIN,PLOT_MAX,VAR1_MIN_COLOR,
     &                  VAR1_MAX_COLOR, VAR1_FACTOR, ZOOM_X,ZOOM_Y,
     &                  VAR1MIN,VAR1MAX,VERT1_LINE,
     &                  VAR2_MIN_COLOR,VAR2_MAX_COLOR,VAR2_FACTOR,
     &                  VAR2MIN , VAR2MAX , VERT2_LINE
