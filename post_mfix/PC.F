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
      CALL GETHOSTNAME(ID_NODE,64)
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
c
      IMPLICIT NONE
      INCLUDE 'funits.inc'
      CALL FLUSH (UNIT_LOG)
      RETURN
      END
c
c
      subroutine flush(n)
	use dflib
	integer n,i
	i = commitqq(n)
	return
	end
c
      double precision function dnrm2(n,x,i)
	integer n,i
	double precision x
	dnrm2 = 1.0
	return
	end
c
      double precision function ddot(n,x1,i,x2,k)
	integer n,i,k
	double precision x1,x2
	ddot = 1.0
	return
	end
c
      subroutine daxpy(n,x1,x2,i,x3,k)
	integer n,i,k
	double precision x1,x2,x3
	return
	end
c
      subroutine dcopy(n,x1,i,x2,k)
	integer n,i,k
	double precision x1,x2
	return
	end
c
      subroutine dscal(n,x1,x2,i)
	integer n,i
	double precision x1,x2
	return
	end
c
      subroutine itime(it)
	integer it(3)
	it(1) = 1
	it(2) = 1
	it(3) = 1
	return
	end
c
      subroutine gethostname(id_node,n)
	character*64 id_node
	integer n
	id_node = 'pc'
	return
	end
c
      real function etime(ta)
	real ta(2)
	ta(1) = 1.0
	ta(2) = 3.0
	etime = ta(1) + ta(2)
	return
	end
c
c
      subroutine pc_quickwin
c
      USE DFLIB
c
      implicit none
c
      include 'param.inc'
      include 'run.inc'
C
C set the version.release of the software
C
      TYPE (windowconfig) wc
      LOGICAL status /.FALSE./

      LOGICAL(4) result
      INTEGER(2) numfonts, fontnum
      TYPE (qwinfo) winfo
      TYPE (xycoord) pos
c
      character title*81
c
      integer i , nb , nc
c
C
c Set the x & y pixels to 800X600
c

      do i = 1,len(run_name)
	   if (run_name(i:i).eq.' ') then
	      nb = i
	      goto 111
	   end if
	end do
	nb = len(run_name)
 111  continue
      title = run_name
	title(nb:nb) = ':'
	nc = 80-nb
	if (nc.ge.len(description)) then
	   title(nb+1:80) = description(1:len(description))
	else
	   title(nb+1:80) = description(1:80-nb)
	end if
	title(81:81) = char(0)
c
      wc.numxpixels  = 800
      wc.numypixels  = 600
      wc.numtextcols = -1
      wc.numtextrows = -1
      wc.numcolors   = -1
c      wc.title= "This is a test"C
      wc.title= title
c      wc.fontsize = #0008000C
      wc.fontsize = -1
      if (.NOT.status) status = SETWINDOWCONFIG(wc)
c	status = setbkcolorrgb(#330000)
c
      winfo.TYPE = QWIN$MAX
c      result = SETWSIZEQQ(QWIN$FRAMEWINDOW, winfo)
      ! Maximize child window
      result=SETWSIZEQQ(0, winfo)
c
      return
      end
c
c
      subroutine slumber
c
      use dfport
c
      do while (.true.)
         call sleep(1000)
      end do
c
      return
      end
