!
!  First, define the user-defined namelist variables as DOUBLE PRECISION,
!  INTEGER, LOGICAL, or CHARACTER.  Second, include the variable in the
!  appropriate common.  Third, add the variable name in the NAMELIST section.
!  A dummy variable is defined, so that the code can be compiled with the
!  template include file.
!
!
!                      Sherwood number
      DOUBLE PRECISION N_Sh(DIMENSION_3, DIMENSION_M)
!
!                      Proximate Analysis
      DOUBLE PRECISION PAFC, PAA
!
!                      Function of ash-layer void fraction
      DOUBLE PRECISION f_EP_A
!
!     Double Precision variables
      COMMON / USR_DATA_DP / &
             N_Sh, PAFC, PAA, f_EP_A
!
!     Integer variables
!      COMMON / USR_DATA_I /
!
!     Logical variables
!      COMMON / USR_DATA_L /
!
!     Character variables
!      COMMON / USR_DATA_C /
!
!
      NAMELIST / USR_INPUT_DATA / &
            PAFC, PAA
