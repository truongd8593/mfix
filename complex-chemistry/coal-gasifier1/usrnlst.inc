!
!  First, define the user-defined namelist variables as DOUBLE PRECISION,
!  INTEGER, LOGICAL, or CHARACTER.  Second, include the variable in the
!  appropriate common.  Third, add the variable name in the NAMELIST section.
!  A dummy variable is defined, so that the code can be compiled with the
!  template include file.
!
! Definitions for Devolatilization Reactions
!
!                      Higher heating value of coal
      DOUBLE PRECISION HHVC
!
!                      Higher heating value of tar
      DOUBLE PRECISION HHVT
!
!                      Proximate Analysis
      DOUBLE PRECISION PAFC, PAVM, PAA, PAM
!
!                      Ultimate Analysis
      DOUBLE PRECISION UAC, UAH, UAO, UAN, UAS
!
!     Double Precision variables
      COMMON / USR_DATA_DP /&
             PAFC, PAVM, PAA, PAM, UAC, UAH, UAO, UAN, UAS,&
             HHVC, HHVT
!
!                      Coal #
      INTEGER          COAL
!
!     Integer variables
      COMMON / USR_DATA_I /&
             COAL
!
!                      Flag to show whether devolatilization is considered
      LOGICAL         USE_DEVOL
!     Logical variables
      COMMON / USR_DATA_L /&
             USE_DEVOL
!
!     Character variables
!      COMMON / USR_DATA_C /
!
!
      NAMELIST / USR_INPUT_DATA /&
             COAL, PAFC, PAVM, PAA, PAM, UAC, UAH, UAO, UAN, UAS,&
             HHVC, HHVT, USE_DEVOL
