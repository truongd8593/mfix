CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: functions.inc                                          C
C  Purpose: This is an include file containing all in-line functions   C
C                                                                      C
C  Author: M. Syamlal                                 Date: dd-mmm-yy  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Revision Number: 1                                                  C
C  Purpose: Add volume-weighted averaging statement functions for      C
C           variable grid capability                                   C
C  Author:  W. Rogers                                 Date: 15-JUL-92  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Literature/Document References: None.                               C
C                                                                      C
C  Variables referenced: None.                                         C
C  Variables modified: None.                                           C
C                                                                      C
C  Local variables: None.                                              C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
C
C                      Function for generating the 3-D array index IJK
C                      from the 1-D indices I, J, and K.
      INTEGER          FUNIJK
C
C                      FUNIJK function with bounds on I, J, and K
      INTEGER          BOUND_FUNIJK
C
C                      Function for generating the index for the entries
C                      to the upper triangle (excluding the diagonal) of
C                      an (L,M) matrix.
      INTEGER          FUNLM
C
C                      Dummy arguments
      DOUBLE PRECISION XXX, XXX1, XXX2, XXX3, XXX4, XXX5
C
C                      Dummy index for L and M
      INTEGER          L1, L2

C                      Dummy indices for I, J, K
      INTEGER          LI, LJ, LK
CC
C                      Function for calculating IJKW
      INTEGER          WEST_OF
C
C                      Function for calculating IJKE
      INTEGER          EAST_OF
C
C                      Function for calculating IJKS
      INTEGER          SOUTH_OF
C
C                      Function for calculating IJKN
      INTEGER          NORTH_OF
C
C                      Function for calculating IJKB
      INTEGER          BOTTOM_OF
C
C                      Function for calculating IJKT
      INTEGER          TOP_OF
C
C                      Function for calculating IMJK
      INTEGER          IM_OF
C
C                      Function for calculating IPJK
      INTEGER          IP_OF
C
C                      Function for calculating IJMK
      INTEGER          JM_OF
C
C                      Function for calculating IJPK
      INTEGER          JP_OF
C
C                      Function for calculating IJKM
      INTEGER          KM_OF
C
C                      Function for calculating IJKP
      INTEGER          KP_OF
C
C                      Logical function to identify a fluid cell
      LOGICAL          FLUID_AT
C
C                      cyclic condition at east boundary
      LOGICAL          CYCLIC_AT_E
C
C                      cyclic condition at north boundary
      LOGICAL          CYCLIC_AT_N
C
C                      cyclic condition at top boundary
      LOGICAL          CYCLIC_AT_T
C
C                      identify flow at east boundary
      LOGICAL          FLOW_AT_E
C
C                      identify flow at north boundary
      LOGICAL          FLOW_AT_N
C
C                      identify flow at top boundary
      LOGICAL          FLOW_AT_T
C
C                      identify const. pressure flow east boundary
      LOGICAL          PFLOW_AT_E
C
C                      identify const. pressure flow north boundary
      LOGICAL          PFLOW_AT_N
C
C                      identify const. pressure flow top boundary
      LOGICAL          PFLOW_AT_T
C
C                      identify specified flow east boundary
      LOGICAL          MFLOW_AT_E
C
C                      identify specified flow north boundary
      LOGICAL          MFLOW_AT_N
C
C                      identify specified flow top boundary
      LOGICAL          MFLOW_AT_T
C
C                      FLUID_AT or P_FLOW_AT (simplified check)
      LOGICAL          FLUIDorP_FLOW_AT
C
C                      P_FLOW_AT 
      LOGICAL          P_FLOW_AT
C
C                      Logical function to identify default walls
      LOGICAL          DEFAULT_WALL_AT
C
C                      Logical function to identify wall ICBC_FLAG
      LOGICAL          WALL_ICBC_FLAG
C
C                      Logical function to identify a No-slip wall cell
      LOGICAL          NS_WALL_AT
C
C                      Logical function to identify a Free-slip wall cell
      LOGICAL          FS_WALL_AT
C
C                      Logical function to identify a Partial-slip wall cell
      LOGICAL          PS_WALL_AT
C
C                      Logical function to identify a wall cell
      LOGICAL          WALL_AT
C
C                      Logical function to identify a cyclic cell
      LOGICAL          CYCLIC_AT
C
C                      Logical function to identify IP at East of the cell
      LOGICAL          IP_AT_E
C
C                      Logical function to identify IP at North of the cell
      LOGICAL          IP_AT_N
C
C                      Logical function to identify IP at Top of the cell
      LOGICAL          IP_AT_T
C
C                      Logical function to identify IP or SP at East of the cell
      LOGICAL          SIP_AT_E
C
C                      Logical function to identify IP or SP at North of the cell
      LOGICAL          SIP_AT_N
C
C                      Logical function to identify IP or SP at Top of the cell
      LOGICAL          SIP_AT_T
C
C                      Logical function to identify SP at East of the cell
      LOGICAL          SP_AT_E
C
C                      Logical function to identify SP at North of the cell
      LOGICAL          SP_AT_N
C
C                      Logical function to identify SP at Top of the cell
      LOGICAL          SP_AT_T
C
C                      Internal surface ID for east face
      INTEGER          IS_ID_AT_E
C
C                      Internal surface ID for north face
      INTEGER          IS_ID_AT_N
C
C                      Internal surface ID for top face
      INTEGER          IS_ID_AT_T
C
C                      Logical function to identify IS at East of the cell
      LOGICAL          IS_AT_E
C
C                      Logical function to identify IS at North of the cell
      LOGICAL          IS_AT_N
C
C                      Logical function to identify IS at Top of the cell
      LOGICAL          IS_AT_T
C
C                      Logical function to identify No IS at East of the cell
      LOGICAL          NO_IS_AT_E
C
C                      Logical function to identify No IS at North of the cell
      LOGICAL          NO_IS_AT_N
C
C                      Logical function to identify No IS at Top of the cell
      LOGICAL          NO_IS_AT_T
C
C                      Maximum of zero or input
      DOUBLE PRECISION ZMAX
C
C 
      FUNIJK (LI, LJ, LK) = LI + (LJ-1)*IMAX2 + (LK-1)*IJMAX2
      BOUND_FUNIJK (LI, LJ, LK) = FUNIJK ( MIN( IMAX2, MAX (1, LI) ),
     &                                  MIN( JMAX2, MAX (1, LJ) ),
     &                                  MIN( KMAX2, MAX (1, LK) ) )
      WEST_OF  (IJK)   = IJK + INCREMENT_FOR_w (CELL_CLASS(IJK))
      EAST_OF  (IJK)   = IJK + INCREMENT_FOR_e (CELL_CLASS(IJK))
      SOUTH_OF (IJK)   = IJK + INCREMENT_FOR_s (CELL_CLASS(IJK))
      NORTH_OF (IJK)   = IJK + INCREMENT_FOR_n (CELL_CLASS(IJK))
      BOTTOM_OF(IJK)   = IJK + INCREMENT_FOR_b (CELL_CLASS(IJK))
      TOP_OF   (IJK)   = IJK + INCREMENT_FOR_t (CELL_CLASS(IJK))
      IM_OF  (IJK)     = IJK + INCREMENT_FOR_im(CELL_CLASS(IJK))
      IP_OF  (IJK)     = IJK + INCREMENT_FOR_ip(CELL_CLASS(IJK))
      JM_OF (IJK)      = IJK + INCREMENT_FOR_jm(CELL_CLASS(IJK))
      JP_OF (IJK)      = IJK + INCREMENT_FOR_jp(CELL_CLASS(IJK))
      KM_OF(IJK)       = IJK + INCREMENT_FOR_km(CELL_CLASS(IJK))
      KP_OF   (IJK)    = IJK + INCREMENT_FOR_kp(CELL_CLASS(IJK))
      FUNLM (L1, L2)   = STORE_LM (L1, L2)
      FLUID_AT(IJK)    = FLAG(IJK) .EQ. 1

      CYCLIC_AT_E(IJK)   = FLAG_E(IJK) .EQ. 2000
      CYCLIC_AT_N(IJK)   = FLAG_N(IJK) .EQ. 2000
      CYCLIC_AT_T(IJK)   = FLAG_T(IJK) .EQ. 2000

      FLOW_AT_E(IJK)   = FLAG_E(IJK) .GE. 2000 .AND.
     &                   FLAG_E(IJK) .LE. 2011
      FLOW_AT_N(IJK)   = FLAG_N(IJK) .GE. 2000 .AND.
     &                   FLAG_N(IJK) .LE. 2011
      FLOW_AT_T(IJK)   = FLAG_T(IJK) .GE. 2000 .AND.
     &                   FLAG_T(IJK) .LE. 2011
      PFLOW_AT_E(IJK)  = FLAG_E(IJK) .EQ. 2010 .OR.
     &                   FLAG_E(IJK) .EQ. 2011
      PFLOW_AT_N(IJK)  = FLAG_N(IJK) .EQ. 2010 .OR.
     &                   FLAG_N(IJK) .EQ. 2011
      PFLOW_AT_T(IJK)  = FLAG_T(IJK) .EQ. 2010 .OR.
     &                   FLAG_T(IJK) .EQ. 2011
      MFLOW_AT_E(IJK)  = FLAG_E(IJK) .EQ. 2020 .OR.
     &                   FLAG_E(IJK) .EQ. 2021
      MFLOW_AT_N(IJK)  = FLAG_N(IJK) .EQ. 2020 .OR.
     &                   FLAG_N(IJK) .EQ. 2021
      MFLOW_AT_T(IJK)  = FLAG_T(IJK) .EQ. 2020 .OR.
     &                   FLAG_T(IJK) .EQ. 2021
      FLUIDorP_FLOW_AT(IJK)= FLAG(IJK) .LE. 11
      P_FLOW_AT(IJK)   = FLAG(IJK) .EQ. 10 .OR.
     &                   FLAG(IJK) .EQ. 11 
      NS_WALL_AT(IJK)  = FLAG(IJK) .EQ. 100
      FS_WALL_AT(IJK)  = FLAG(IJK) .EQ. 101
      PS_WALL_AT(IJK)  = FLAG(IJK) .EQ. 102
      CYCLIC_AT(IJK)   = FLAG(IJK) .EQ. 106 .OR.
     &                   FLAG(IJK) .EQ. 107
      WALL_AT(IJK)     = FLAG(IJK) .GE. 100
      IP_AT_E(IJK)     = FLAG_E(IJK) .LT. 1000
      IP_AT_N(IJK)     = FLAG_N(IJK) .LT. 1000
      IP_AT_T(IJK)     = FLAG_T(IJK) .LT. 1000
      SIP_AT_E(IJK)    = (FLAG_E(IJK) .LT. 2000)
      SIP_AT_N(IJK)    = (FLAG_N(IJK) .LT. 2000)
      SIP_AT_T(IJK)    = (FLAG_T(IJK) .LT. 2000)
      SP_AT_E(IJK)     = (FLAG_E(IJK) .LT. 2000) .AND.
     &                   (FLAG_E(IJK) .GE. 1000)
      SP_AT_N(IJK)     = (FLAG_N(IJK) .LT. 2000) .AND.
     &                   (FLAG_N(IJK) .GE. 1000)
      SP_AT_T(IJK)     = (FLAG_T(IJK) .LT. 2000) .AND.
     &                   (FLAG_T(IJK) .GE. 1000)
      IS_ID_AT_E(IJK)  = FLAG_E(IJK) - 1000
      IS_ID_AT_N(IJK)  = FLAG_N(IJK) - 1000
      IS_ID_AT_T(IJK)  = FLAG_T(IJK) - 1000
      IS_AT_E(IJK)     = FLAG_E(IJK) .LT. 2000
      IS_AT_N(IJK)     = FLAG_N(IJK) .LT. 2000
      IS_AT_T(IJK)     = FLAG_T(IJK) .LT. 2000
      NO_IS_AT_E(IJK)  = FLAG_E(IJK) .GE. 2000
      NO_IS_AT_N(IJK)  = FLAG_N(IJK) .GE. 2000
      NO_IS_AT_T(IJK)  = FLAG_T(IJK) .GE. 2000
      WALL_ICBC_FLAG(IJK) = ICBC_FLAG(IJK)(1:1) .EQ. 'W' .OR.
     &                      ICBC_FLAG(IJK)(1:1) .EQ. 'S' .OR.
     &                      ICBC_FLAG(IJK)(1:1) .EQ. 's' .OR.
     &                      ICBC_FLAG(IJK)(1:1) .EQ. 'c' .OR.
     &                      ICBC_FLAG(IJK)(1:1) .EQ. 'C'
      DEFAULT_WALL_AT(IJK) = ICBC_FLAG(IJK)(2:3) .EQ. '--' 
      ZMAX (XXX)       = MAX(XXX, ZERO)
