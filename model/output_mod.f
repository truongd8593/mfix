!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: output_control.inc                                     C
!  Purpose: Common block containing data for output control            C
!                                                                      C
!  Author: M. Syamlal                                 Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C


      MODULE output


      Use param
      Use param1


!
!                      Interval at which restart (.RES) file data is
!                      updated.
      DOUBLE PRECISION RES_DT
!
!                      Interval at which REAL restart (.SPx) files data
!                      are updated.
      DOUBLE PRECISION SPX_DT(N_SPX)
!
!                      Interval at which standard output (.OUT) file data
!                      is updated.
      DOUBLE PRECISION OUT_DT
!
!
!                      x coordinate of the west face of a region where
!                      output is to be done
      DOUBLE PRECISION USR_X_w (DIMENSION_USR)
!
!                      x coordinate of the east face of a region where
!                      output is to be done
      DOUBLE PRECISION USR_X_e (DIMENSION_USR)
!
!                      y coordinate of the south face of a region where
!                      output is to be done
      DOUBLE PRECISION USR_Y_s (DIMENSION_USR)
!
!                      y coordinate of the north face of a region where
!                      output is to be done
      DOUBLE PRECISION USR_Y_n (DIMENSION_USR)
!
!                      z coordinate of the bottom face of a region where
!                      output is to be done
      DOUBLE PRECISION USR_Z_b (DIMENSION_USR)
!
!                      z coordinate of the top face of a region where
!                      output is to be done
      DOUBLE PRECISION USR_Z_t (DIMENSION_USR)
!
!                      i index of the west face of a region where
!                      output is to be done
      INTEGER          USR_I_w (DIMENSION_USR)
!
!                      i index of the east face of a region where
!                      output is to be done
      INTEGER          USR_I_e (DIMENSION_USR)
!
!                      j index of the south face of a region where
!                      output is to be done
      INTEGER          USR_J_s (DIMENSION_USR)
!
!                      j index of the north face of a region where
!                      output is to be done
      INTEGER          USR_J_n (DIMENSION_USR)
!
!                      k index of the bottom face of a region where
!                      output is to be done
      INTEGER          USR_K_b (DIMENSION_USR)
!
!                      k index of the top face of a region where
!                      output is to be done
      INTEGER          USR_K_t (DIMENSION_USR)
!
!                      Interval at which user-defined output files are updated.
      DOUBLE PRECISION USR_DT (DIMENSION_USR)
!
!                      Type of user-defined output: BINARY or ASCII.
      CHARACTER*6      USR_TYPE (DIMENSION_USR)
!
!                      Variables to be written in the user-defined output file.
      CHARACTER*60     USR_VAR (DIMENSION_USR)
!
!                      Format for writing user-defined (ASCII) output file.
      CHARACTER*60     USR_FORMAT (DIMENSION_USR)
!
!                      Extension for the user-defined output file.
      CHARACTER*16     USR_EXT (DIMENSION_USR)
!
!                      Interval in number of time steps at which LOG file
!                      is written
      INTEGER          NLOG
!
!                      If true, display the residuals on the screen and messages
!                      about convergence on the screen and in the .LOG file.
      LOGICAL          FULL_LOG
!
      INTEGER          NLOG_HEADER
      LOGICAL          LOG_HEADER
!
!
!              User defined flag to enable log file writing from each processor in DMP mode
!              If .false. only PE_IO writes the log file
!AEOLUS added to mfix.dat namelist as input parameter
!     LOGICAL, PARAMETER          :: ENABLE_DMP_LOG = .FALSE.
      LOGICAL :: ENABLE_DMP_LOG


!                      Interval at which mass balance is checked and reported in .LOG file.
!                      This value is used in check_mass_balance and is by default undefined.
      DOUBLE PRECISION report_mass_balance_dt

!AEOLUS Flag variable which controls debug print of whole index layout to be used in determining ijk<=>i,j,k and other similar
! debugging tasks
      LOGICAL          DBGPRN_LAYOUT

      END MODULE output                                                                          
