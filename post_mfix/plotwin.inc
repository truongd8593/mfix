CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name:                                                        C
C  Purpose:                                                            C
C                                                                      C
C  Author: P. Nicoletti                               Date: 12-NOV-92  C
C  Reviewer:                                                           C
C                                                                      C
C  Revision Number:                                                    C
C  Purpose:                                                            C
C  Author:                                            Date: dd-mmm-yy  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Literature/Document References:                                     C
C                                                                      C
C  Variables referenced:                                               C
C  Variables modified:                                                 C
C                                                                      C
C  Local variables:                                                    C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
C
      INTEGER  MAX_WINDOWS
      PARAMETER (MAX_WINDOWS=30)
C
C
      INTEGER   GID(MAX_WINDOWS) , XS1(MAX_WINDOWS) , XS2(MAX_WINDOWS)
      INTEGER   YS1(MAX_WINDOWS) , YS2(MAX_WINDOWS) , XXMAX , YYMAX
      INTEGER   MENU_WINDOW_INDEX
C
      REAL      XPIXELS(MAX_WINDOWS) , YPIXELS(MAX_WINDOWS)
      REAL      XRATIO , YRATIO , vcr_scale
C
      LOGICAL   CAN_DOUBLE_BUFFER , VCR_SIZE
C
      COMMON /PLOT_WINDOWS_I/ GID,XS1,XS2,YS1,YS2,XXMAX,YYMAX,
     &       MENU_WINDOW_INDEX
C
      COMMON /PLOT_WINDOWS_R/ XPIXELS,YPIXELS,XRATIO,YRATIO,vcr_scale
C
      COMMON /PLOT_WINDOWS_L/ CAN_DOUBLE_BUFFER,VCR_SIZE
