C
C  First, define the user-defined namelist variables as DOUBLE PRECISION,
C  INTEGER, LOGICAL, or CHARACTER.  Second, include the variable in the
C  appropriate common.  Third, add the variable name in the NAMELIST section.
C  A dummy variable is defined, so that the code can be compiled with the
C  template include file.
C
      DOUBLE PRECISION DUMMY_DP
C
C     Double Precision variables
      COMMON / USR_DATA_DP /
     &        DUMMY_DP
C
C     Integer variables
C      COMMON / USR_DATA_I /
C
C     Logical variables
C      COMMON / USR_DATA_L /
C
C     Character variables
C      COMMON / USR_DATA_C /
C
C
      NAMELIST / USR_INPUT_DATA /
     &         DUMMY_DP
