      MODULE time_cpu


      Use param
      Use param1


!                      cpu time/second
      DOUBLE PRECISION CPUos
!                      old cpu time and time for calculating CPUos
      DOUBLE PRECISION CPU_NLOG, TIME_NLOG
!                       Initial value of CPU time.
      DOUBLE PRECISION CPU0
!
!			Time for IO
	DOUBLE PRECISION CPU0_IO, CPU1_IO, CPU_IO


!AEOLUS STOP Initial value of CPU time at the begin of MFIX, prior any I/O
      DOUBLE PRECISION CPU00

      END MODULE time_cpu                                                                        
