C                      cpu time/second
      DOUBLE PRECISION CPUos
C                      old cpu time and time for calculating CPUos
      DOUBLE PRECISION CPU_NLOG, TIME_NLOG
C                       Initial value of CPU time.
      DOUBLE PRECISION CPU0
C
      COMMON /TIME_CPU/ CPU0, CPUos, CPU_NLOG, TIME_NLOG
