!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: MACHINE_CONS                                           C
!  Purpose: set the machine constants    ( SGI ONLY )                  C
!                                                                      C
!  Author: P. Nicoletti                               Date: 28-JAN-92  C
!  Reviewer: P. Nicoletti, W. Rogers, M. Syamlal      Date:            C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: None                                          C
!  Variables modified: OPEN_N1, NWORDS_DP, NWORDS_R, N_WORDS_I         C
!                                                                      C
!  Local variables: None                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE MACHINE_CONS
!
!
      Use machine
      IMPLICIT NONE
!
      OPEN_N1   = 128
      NWORDS_DP =  64
      NWORDS_R  = 128
      NWORDS_I  = 128
!
      RETURN
      END
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: GET_RUN_ID                                             C
!  Purpose: get the run id for this run                                C
!                                                                      C
!  Author: P. Nicoletti                               Date: 16-DEC-91  C
!  Reviewer: P. Nicoletti, W. Rogers, M. Syamlal      Date:            C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: add ndoe name                                              C
!  Author: P.Nicoletti                                Date: 07-FEB-92  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: None                                          C
!  Variables modified: ID_MONTH, ID_DAY, ID_YEAR, ID_HOUR, ID_MINUTE   C
!                      ID_SECOND, ID_NODE                              C
!                                                                      C
!  Local variables: TIME_ARRAY                                         C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE GET_RUN_ID
!
      USE param
      USE run
      IMPLICIT NONE
!
!             temporary array to hold time data
      INTEGER DAT(8)
      CHARACTER*10 DATE, TIM, ZONE

      CALL DATE_AND_TIME(DATE, TIM, ZONE, DAT)
      ID_YEAR   = DAT(1)
      ID_MONTH  = DAT(2)
      ID_DAY    = DAT(3)
      ID_HOUR   = DAT(5)
      ID_MINUTE = DAT(6)
      ID_SECOND = DAT(7)
      
!     For SGI only
      CALL GETHOSTNAME(ID_NODE,64)
!
      RETURN
      END
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CPU_TIME (CPU)                                         C
!  Purpose: get the CPU time for the run                               C
!                                                                      C
!  Author: P. Nicoletti                               Date: 10-JAN-92  C
!  Reviewer: P. Nicoletti, W. Rogers, M. Syamlal      Date:            C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: None                                          C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables: TA, XT                                             C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CPU_TIME(CPU)
!
      IMPLICIT NONE
!
! passed arguments
!
!                      cpu time since start of run
      DOUBLE PRECISION CPU
!
! local variables
!

!                       clock cycle
      INTEGER           COUNT

!                       number of cycles per second
      INTEGER           COUNT_RATE
      
!                       max number of cycles, after which count is reset to 0
      INTEGER           COUNT_MAX

      CALL SYSTEM_CLOCK(COUNT, COUNT_RATE, COUNT_MAX)
      
      CPU           = DBLE(COUNT)/DBLE(COUNT_RATE)
!
      RETURN
      END
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: START_LOG                                              C
!  Purpose: does nothing ... for VAX/VMS compatibility (SGI ONLY)      C
!                                                                      C
!  Author: P. Nicoletti                               Date: 28-JAN-92  C
!  Reviewer: P. Nicoletti, W. Rogers, M. Syamlal      Date:            C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: None                                          C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables: None                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE START_LOG
      IMPLICIT NONE
      RETURN
      END
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: END_LOG                                                C
!  Purpose: flushes the log file                                       C
!                                                                      C
!  Author: P. Nicoletti                               Date: 28-JAN-92  C
!  Reviewer: P. Nicoletti, W. Rogers, M. Syamlal      Date:            C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: None                                          C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables: None                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE END_LOG
      USE funits
      IMPLICIT NONE
      CALL FLUSH (UNIT_LOG)
      RETURN
      END
!
!
      subroutine flush(n)
	use dflib
	integer n,i
	i = commitqq(n)
	return
	end
!
      double precision function dnrm2(n,x,i)
	integer n,i
	double precision x
	dnrm2 = 1.0
	return
	end
!
      double precision function ddot(n,x1,i,x2,k)
	integer n,i,k
	double precision x1,x2
	ddot = 1.0
	return
	end
!
      subroutine daxpy(n,x1,x2,i,x3,k)
	integer n,i,k
	double precision x1,x2,x3
	return
	end
!
      subroutine dcopy(n,x1,i,x2,k)
	integer n,i,k
	double precision x1,x2
	return
	end
!
      subroutine dscal(n,x1,x2,i)
	integer n,i
	double precision x1,x2
	return
	end
!
      subroutine itime(it)
	integer it(3)
	it(1) = 1
	it(2) = 1
	it(3) = 1
	return
	end
!
      subroutine gethostname(id_node,n)
	character*64 id_node
	integer n
	id_node = 'pc'
	return
	end
!
      real function etime(ta)
	real ta(2)
	ta(1) = 1.0
	ta(2) = 3.0
	etime = ta(1) + ta(2)
	return
	end
!
!
      subroutine pc_quickwin
!
      USE DFLIB
!
!
      Use param
      Use run
      implicit none
!
! set the version.release of the software
!
      TYPE (windowconfig) wc
      LOGICAL status /.FALSE./

      LOGICAL(4) result
      INTEGER(2) numfonts, fontnum
      TYPE (qwinfo) winfo
      TYPE (xycoord) pos
!
      character title*81
!
      integer i , nb , nc
!
!
! Set the x & y pixels to 800X600
!

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
!
      wc.numxpixels  = 800
      wc.numypixels  = 600
      wc.numtextcols = -1
      wc.numtextrows = -1
      wc.numcolors   = -1
!      wc.title= "This is a test"C
      wc.title= title
!      wc.fontsize = #0008000C
      wc.fontsize = -1
      if (.NOT.status) status = SETWINDOWCONFIG(wc)
!	status = setbkcolorrgb(#330000)
!
      winfo.TYPE = QWIN$MAX
!      result = SETWSIZEQQ(QWIN$FRAMEWINDOW, winfo)
      ! Maximize child window
      result=SETWSIZEQQ(0, winfo)
!
      return
      end
!
!
      subroutine slumber
!
      use dfport
!
      do while (.true.)
         call sleep(1000)
      end do
!
      return
      end
!
      subroutine flush_bin(iunit,code)
      implicit none
      integer :: iunit , code
      call flush(iunit)
      return
      end
