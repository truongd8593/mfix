CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: output_control.inc                                     C
C  Purpose: Common block containing data for output control            C
C                                                                      C
C  Author: M. Syamlal                                 Date: dd-mmm-yy  C
C  Reviewer:                                          Date: dd-mmm-yy  C
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
C                      Interval at which restart (.RES) file data is
C                      updated.
      DOUBLE PRECISION RES_DT
C
C                      Interval at which REAL restart (.SPx) files data
C                      are updated.
      DOUBLE PRECISION SPX_DT(N_SPX)
C
C                      Interval at which standard output (.OUT) file data
C                      is updated.
      DOUBLE PRECISION OUT_DT
C
C
C                      x coordinate of the west face of a region where
C                      output is to be done
      DOUBLE PRECISION USR_X_w (DIMENSION_USR) 
C
C                      x coordinate of the east face of a region where
C                      output is to be done
      DOUBLE PRECISION USR_X_e (DIMENSION_USR)
C
C                      y coordinate of the south face of a region where
C                      output is to be done
      DOUBLE PRECISION USR_Y_s (DIMENSION_USR)
C
C                      y coordinate of the north face of a region where
C                      output is to be done
      DOUBLE PRECISION USR_Y_n (DIMENSION_USR)
C
C                      z coordinate of the bottom face of a region where
C                      output is to be done
      DOUBLE PRECISION USR_Z_b (DIMENSION_USR)
C
C                      z coordinate of the top face of a region where
C                      output is to be done
      DOUBLE PRECISION USR_Z_t (DIMENSION_USR)
C
C                      i index of the west face of a region where
C                      output is to be done
      INTEGER          USR_I_w (DIMENSION_USR) 
C
C                      i index of the east face of a region where
C                      output is to be done
      INTEGER          USR_I_e (DIMENSION_USR)
C
C                      j index of the south face of a region where
C                      output is to be done
      INTEGER          USR_J_s (DIMENSION_USR)
C
C                      j index of the north face of a region where
C                      output is to be done
      INTEGER          USR_J_n (DIMENSION_USR)
C
C                      k index of the bottom face of a region where
C                      output is to be done
      INTEGER          USR_K_b (DIMENSION_USR)
C
C                      k index of the top face of a region where
C                      output is to be done
      INTEGER          USR_K_t (DIMENSION_USR)
C
C                      Interval at which user-defined output files are updated.
      DOUBLE PRECISION USR_DT (DIMENSION_USR)
C
C                      Type of user-defined output: BINARY or ASCII.
      CHARACTER*6      USR_TYPE (DIMENSION_USR)
C
C                      Variables to be written in the user-defined output file.
      CHARACTER*60     USR_VAR (DIMENSION_USR)
C
C                      Format for writing user-defined (ASCII) output file.
      CHARACTER*60     USR_FORMAT (DIMENSION_USR)
C
C                      Extension for the user-defined output file.
      CHARACTER*16     USR_EXT (DIMENSION_USR)
C
C                      Interval in number of time steps at which LOG file
C                      is written
      INTEGER          NLOG
C
C                      If true, display the residuals on the screen and messages
C                      about convergence on the screen and in the .LOG file.
      LOGICAL          FULL_LOG
C
      INTEGER          NLOG_HEADER
      LOGICAL          LOG_HEADER
C
      COMMON / OUTPUT_CONTROL_DP /
     &       RES_DT, SPX_DT , OUT_DT, 
     &       USR_DT, USR_X_w , USR_X_e , USR_Y_s , USR_Y_n ,
     &       USR_Z_b , USR_Z_t
C
      COMMON / OUTPUT_CONTROL_I /
     &       USR_I_w , USR_I_e , USR_J_s , USR_J_n ,
     &       USR_K_b , USR_K_t,
     &       NLOG, NLOG_HEADER
C
      COMMON / OUTPUT_CONTROL_C /
     &       USR_TYPE, USR_VAR, USR_FORMAT, USR_EXT
C
      COMMON /OUTPUT_CONTROL_L /
     &       LOG_HEADER, FULL_LOG
