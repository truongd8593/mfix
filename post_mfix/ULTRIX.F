CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: MACHINE_CONS                                           C
C  Purpose: set the machine constants    ( SGI ONLY )                  C
C                                                                      C
C  Author: P. Nicoletti                               Date: 28-JAN-92  C
C  Reviewer: P. Nicoletti, W. Rogers, M. Syamlal      Date:            C
C                                                                      C
C  Revision Number:                                                    C
C  Purpose:                                                            C
C  Author:                                            Date: dd-mmm-yy  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Literature/Document References:                                     C
C                                                                      C
C  Variables referenced: None                                          C
C  Variables modified: OPEN_N1, NWORDS_DP, NWORDS_R, N_WORDS_I         C
C                                                                      C
C  Local variables: None                                               C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
C
      SUBROUTINE MACHINE_CONS
C
      IMPLICIT NONE
C
      INCLUDE 'machine.inc'
C
      OPEN_N1   = 128
      NWORDS_DP =  64
      NWORDS_R  = 128
      NWORDS_I  = 128
C
      RETURN
      END
CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: GET_RUN_ID                                             C
C  Purpose: get the run id for this run  (**** SGI ONLY ****)          C
C                                                                      C
C  Author: P. Nicoletti                               Date: 16-DEC-91  C
C  Reviewer: P. Nicoletti, W. Rogers, M. Syamlal      Date:            C
C                                                                      C
C  Revision Number: 1                                                  C
C  Purpose: add ndoe name                                              C
C  Author: P.Nicoletti                                Date: 07-FEB-92  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Literature/Document References:                                     C
C                                                                      C
C  Variables referenced: None                                          C
C  Variables modified: ID_MONTH, ID_DAY, ID_YEAR, ID_HOUR, ID_MINUTE   C
C                      ID_SECOND, ID_NODE                              C
C                                                                      C
C  Local variables: TIME_ARRAY                                         C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
C
      SUBROUTINE GET_RUN_ID
C
      IMPLICIT NONE
      INCLUDE 'param.inc'
      INCLUDE 'run.inc'
C
C             temporary array to hold time data
      INTEGER TIME_ARRAY(3)
C
      CALL IDATE (ID_MONTH,ID_DAY,ID_YEAR)
      CALL ITIME (TIME_ARRAY)
      ID_HOUR   = TIME_ARRAY(1)
      ID_MINUTE = TIME_ARRAY(2)
      ID_SECOND = TIME_ARRAY(3)
C
C      CALL GETHOSTNAME(ID_NODE,64)
      ID_NODE = 'DEC ULTRIX Machine'
C
      RETURN
      END
CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: CPU_TIME (CPU)                                         C
C  Purpose: get the CPU time for the run  (**** SGI ONLY ****)         C
C                                                                      C
C  Author: P. Nicoletti                               Date: 10-JAN-92  C
C  Reviewer: P. Nicoletti, W. Rogers, M. Syamlal      Date:            C
C                                                                      C
C  Revision Number:                                                    C
C  Purpose:                                                            C
C  Author:                                            Date: dd-mmm-yy  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Literature/Document References:                                     C
C                                                                      C
C  Variables referenced: None                                          C
C  Variables modified: None                                            C
C                                                                      C
C  Local variables: TA, XT                                             C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
C
      SUBROUTINE CPU_TIME(CPU)
C
      IMPLICIT NONE
C
C passed arguments
C
C                      cpu time since start of run
      DOUBLE PRECISION CPU
C
C local variables
C
C                      TA(1) = user cpu time   TA(2) = system cpu time
      REAL             TA(2) 
C
C                      XT = TA(1) + TA(2)
      REAL             XT
C
C                      ETIME is an SGI system function which returns
C                      the elasped CPU time
      REAL             ETIME
C
      XT  = ETIME(TA)
      CPU = XT
C
      RETURN
      END
CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: START_LOG                                              C
C  Purpose: does nothing ... for VAX/VMS compatibility (SGI ONLY)      C
C                                                                      C
C  Author: P. Nicoletti                               Date: 28-JAN-92  C
C  Reviewer: P. Nicoletti, W. Rogers, M. Syamlal      Date:            C
C                                                                      C
C  Revision Number:                                                    C
C  Purpose:                                                            C
C  Author:                                            Date: dd-mmm-yy  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Literature/Document References:                                     C
C                                                                      C
C  Variables referenced: None                                          C
C  Variables modified: None                                            C
C                                                                      C
C  Local variables: None                                               C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
C
      SUBROUTINE START_LOG
      IMPLICIT NONE
      RETURN
      END
CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: END_LOG                                                C
C  Purpose: flushes the log file                                       C
C                                                                      C
C  Author: P. Nicoletti                               Date: 28-JAN-92  C
C  Reviewer: P. Nicoletti, W. Rogers, M. Syamlal      Date:            C
C                                                                      C
C  Revision Number:                                                    C
C  Purpose:                                                            C
C  Author:                                            Date: dd-mmm-yy  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Literature/Document References:                                     C
C                                                                      C
C  Variables referenced: None                                          C
C  Variables modified: None                                            C
C                                                                      C
C  Local variables: None                                               C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
C
      SUBROUTINE END_LOG
      IMPLICIT NONE
      INCLUDE 'funits.inc'
      CALL FLUSH (UNIT_LOG)
      RETURN
      END
c
      subroutine slumber
      return
      end
c
      subroutine pc_quickwin
      return
      end
