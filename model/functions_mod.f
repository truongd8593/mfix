MODULE functions

  USE compar
  USE geometry
  USE indices

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  !
  !                      Function for generating the LOCAL 3-D array index IJK
  !                      from the 1-D indices I, J, and K.
  !//FUNIJK is moved to compar for debugging purposes - Sreekanth-10/26/99
  !      INTEGER          FUNIJK,FUNIJK_0
  !      INTEGER          FUNIJK_0
  !                      Function for generating the LOCAL 3-D array index IJK
  !                      from the 1-D indices I, J, K and IPROC.
  !      INTEGER          FUNIJK_PROC
  !                      Function for generating the GLOBAL 3-D array index IJK
  !                      from the 1-D indices I, J, and K.
  !      INTEGER          FUNIJK_GL
  !
  !                      Function for generating the index for the entries
  !                      to the upper triangle (excluding the diagonal) of
  !                      an (L,M) matrix.
  !      INTEGER          FUNLM
  !
  !                      Function for generating the 3-D array index IJK
  !                      from the 1-D indices I, J, and K in IO format
  !      INTEGER          FUNIJK_IO
  !
  !                      Dummy arguments
  !      DOUBLE PRECISION XXX, YYY, XXX1, XXX2, XXX3, XXX4, XXX5
  !
  !                      Dummy index for L and M
  !      INTEGER          L1, L2

  !                      Dummy indices for I, J, K
  !      INTEGER          LI, LJ, LK
  !                      Dummy indices for IPROC
  !      INTEGER          LIPROC
  !
  !                      Function for calculating IJKW
  !      INTEGER          WEST_OF, WEST_OF_0
  !
  !                      Function for calculating IJKE
  !      INTEGER          EAST_OF, EAST_OF_0
  !
  !                      Function for calculating IJKS
  !      INTEGER          SOUTH_OF, SOUTH_OF_0
  !
  !                      Function for calculating IJKN
  !      INTEGER          NORTH_OF, NORTH_OF_0
  !
  !                      Function for calculating IJKB
  !      INTEGER          BOTTOM_OF, BOTTOM_OF_0
  !
  !                      Function for calculating IJKT
  !      INTEGER          TOP_OF, TOP_OF_0
  !
  !                      Function for calculating IMJK
  !      INTEGER          IM_OF, IM_OF_0
  !
  !                      Function for calculating IPJK
  !      INTEGER          IP_OF, IP_OF_0
  !
  !                      Function for calculating IJMK
  !      INTEGER          JM_OF, JM_OF_0
  !
  !                      Function for calculating IJPK
  !      INTEGER          JP_OF, JP_OF_0
  !
  !                      Function for calculating IJKM
  !      INTEGER          KM_OF, KM_OF_0
  !
  !                      Function for calculating IJKP
  !      INTEGER          KP_OF, KP_OF_0
  !
  !                      Logical function to identify a fluid cell
  !      LOGICAL          FLUID_AT
  !
  !                      cyclic condition at east boundary
  !      LOGICAL          CYCLIC_AT_E
  !
  !                      cyclic condition at north boundary
  !      LOGICAL          CYCLIC_AT_N
  !
  !                      cyclic condition at top boundary
  !      LOGICAL          CYCLIC_AT_T
  !
  !                      identify flow at east boundary
  !      LOGICAL          FLOW_AT_E
  !
  !                      identify flow at north boundary
  !      LOGICAL          FLOW_AT_N
  !
  !                      identify flow at top boundary
  !      LOGICAL          FLOW_AT_T
  !
  !                      identify const. pressure flow east boundary
  !      LOGICAL          PFLOW_AT_E
  !
  !                      identify const. pressure flow north boundary
  !      LOGICAL          PFLOW_AT_N
  !
  !                      identify const. pressure flow top boundary
  !      LOGICAL          PFLOW_AT_T
  !
  !                      identify specified flow east boundary
  !      LOGICAL          MFLOW_AT_E
  !
  !                      identify specified flow north boundary
  !      LOGICAL          MFLOW_AT_N
  !
  !                      identify specified flow top boundary
  !      LOGICAL          MFLOW_AT_T
  !
  !                      FLUID_AT or P_FLOW_AT (simplified check)
  !      LOGICAL          FLUIDorP_FLOW_AT
  !
  !                      P_FLOW_AT
  !      LOGICAL          P_FLOW_AT
  !
  !                      P_OUTFLOW_AT
  !      LOGICAL          P_OUTFLOW_AT
  !
  !                      MASS_OUTFLOW_AT
  !      LOGICAL          MASS_OUTFLOW_AT
  !
  !                      OUTFLOW_AT
  !      LOGICAL          OUTFLOW_AT
  !
  !                      FLOW_AT
  !      LOGICAL          FLOW_AT
  !
  !                      Logical function to identify default walls
  !      LOGICAL          DEFAULT_WALL_AT
  !
  !                      Logical function to identify wall ICBC_FLAG
  !      LOGICAL          WALL_ICBC_FLAG
  !
  !                      Logical function to identify a No-slip wall cell
  !      LOGICAL          NS_WALL_AT
  !
  !                      Logical function to identify a Free-slip wall cell
  !      LOGICAL          FS_WALL_AT
  !
  !                      Logical function to identify a Partial-slip wall cell
  !      LOGICAL          PS_WALL_AT
  !
  !                      Logical function to identify a wall cell
  !      LOGICAL          WALL_AT
  !
  !                      Logical function to identify a cyclic cell
  !      LOGICAL          CYCLIC_AT
  !
  !                      Logical function to identify IP at East of the cell
  !      LOGICAL          IP_AT_E
  !
  !                      Logical function to identify IP at North of the cell
  !      LOGICAL          IP_AT_N
  !
  !                      Logical function to identify IP at Top of the cell
  !      LOGICAL          IP_AT_T
  !
  !                      Logical function to identify IP or SP at East of the cell
  !      LOGICAL          SIP_AT_E
  !
  !                      Logical function to identify IP or SP at North of the cell
  !      LOGICAL          SIP_AT_N
  !
  !                      Logical function to identify IP or SP at Top of the cell
  !      LOGICAL          SIP_AT_T
  !
  !                      Logical function to identify SP at East of the cell
  !      LOGICAL          SP_AT_E
  !
  !                      Logical function to identify SP at North of the cell
  !      LOGICAL          SP_AT_N
  !
  !                      Logical function to identify SP at Top of the cell
  !      LOGICAL          SP_AT_T
  !
  !                      Internal surface ID for east face
  !      INTEGER          IS_ID_AT_E
  !
  !                      Internal surface ID for north face
  !      INTEGER          IS_ID_AT_N
  !
  !                      Internal surface ID for top face
  !      INTEGER          IS_ID_AT_T
  !
  !                      Logical function to identify IS at East of the cell
  !      LOGICAL          IS_AT_E
  !
  !                      Logical function to identify IS at North of the cell
  !      LOGICAL          IS_AT_N
  !
  !                      Logical function to identify IS at Top of the cell
  !      LOGICAL          IS_AT_T
  !
  !                      Logical function to identify No IS at East of the cell
  !      LOGICAL          NO_IS_AT_E
  !
  !                      Logical function to identify No IS at North of the cell
  !      LOGICAL          NO_IS_AT_N
  !
  !                      Logical function to identify No IS at Top of the cell
  !      LOGICAL          NO_IS_AT_T
  !\\   Added to determine whether I am on my PE
  !\\   Logical function to determine whether I am on my PE's domain
  !      LOGICAL          IS_ON_myPE_owns
  !\\   Logical function to determine whether I am on my PE's domain + one &
  !                                                       extra ghost layer
  !      LOGICAL          IS_ON_myPE_plus1layer
  !\\   Logical function to determine whether I am on my PE's domain + two &
  !                                                       extra ghost layers
  !      LOGICAL          IS_ON_myPE_plus2layers
  !\\   Logical function to determine whether I am on my PE's domain - Boundaries
  !      LOGICAL          IS_ON_myPE_wobnd
  !
  !                      Maximum of zero or input
  !      DOUBLE PRECISION ZMAX
  !
  !
CONTAINS

  !//FUNIJK is moved to compar for debugging purposes - Sreekanth-10/26/99
  !     FUNIJK (LI, LJ, LK) = c0 + LI + (LJ-jstart3_all(myPE))*c1 + (LK-kstart3_all(myPE))* c2
  !      funijk(li,lj,lk) = lj + c0 + li*c1 + lk*c2
  INTEGER FUNCTION funijk_0(li,lj,lk)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: LI, LJ, LK
    funijk_0 = lj + c0 + li*c1 + lk*c2
  END FUNCTION funijk_0

  INTEGER FUNCTION funijk(li,lj,lk)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: LI, LJ, LK
    funijk = IJK_ARRAY_OF(li,lj,lk)
  END FUNCTION funijk

  !//SP
  !     FUNIJK_PROC(LI, LJ, LK, LIPROC) = 1 + (LI - istart3_all(LIPROC))+ &
  !     (LJ-jstart3_all(LIPROC))*(iend3_all(LIPROC)-istart3_all(LIPROC)+1) &
  !     + (LK-kstart3_all(LIPROC))*(jend3_all(LIPROC)-jstart3_all(LIPROC)+1)* &
  !     (iend3_all(LIPROC)-istart3_all(LIPROC)+1)

  INTEGER FUNCTION FUNIJK_PROC(LI, LJ, LK, LIPROC)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: LI, LJ, LK, LIPROC
    FUNIJK_PROC = 1 + (LJ - jstart3_all(LIPROC))+ &
         (LI-Istart3_all(LIPROC))*(jend3_all(LIPROC)-jstart3_all(LIPROC)+1) &
         + (LK-kstart3_all(LIPROC))*(jend3_all(LIPROC)-jstart3_all(LIPROC)+1)* &
         (iend3_all(LIPROC)-istart3_all(LIPROC)+1)
  END FUNCTION FUNIJK_PROC

  !     FUNIJK_GL (LI, LJ, LK) = 1 + (LI - imin3) + (LJ-jmin3)*(imax3-imin3+1) &
  !     + (LK-kmin3)*(jmax3-jmin3+1)*(imax3-imin3+1)
  !
  INTEGER FUNCTION FUNIJK_GL (LI, LJ, LK)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: LI, LJ, LK
    FUNIJK_GL = 1 + (LJ - jmin3) + (LI-imin3)*(jmax3-jmin3+1) &
         + (LK-kmin3)*(jmax3-jmin3+1)*(imax3-imin3+1)
  END FUNCTION FUNIJK_GL

  INTEGER FUNCTION FUNIJK_IO (LI, LJ, LK)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: LI, LJ, LK
    FUNIJK_IO = 1 + (LI - imin2) + (LJ-jmin2)*(imax2-imin2+1) &
         + (LK-kmin2)*(jmax2-jmin2+1)*(imax2-imin2+1)
  END FUNCTION FUNIJK_IO

  LOGICAL FUNCTION IS_ON_myPE_plus2layers (LI, LJ, LK)
    IS_ON_myPE_plus2layers = LI.ge.istart3.and.LI.le.iend3.and. &
         LJ.ge.jstart3.and.LJ.le.jend3.and. &
         LK.ge.kstart3.and.LK.le.kend3
  END FUNCTION IS_ON_myPE_plus2layers

  !      IS_ON_myPE_plus2layers (LI, LJ, LK) = LI.ge.istart3.and.LI.le.iend3.and. &
  !                               LJ.ge.jstart3.and.LJ.le.jend3.and. &
  !                                LK.ge.kstart3.and.LK.le.kend3.and. &
  !                                (.NOT.DEAD_CELL_AT(LI,LJ,LK))


  LOGICAL FUNCTION IS_ON_myPE_plus1layer (LI, LJ, LK)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: LI, LJ, LK

    IS_ON_myPE_plus1layer = LI.ge.istart2.and.LI.le.iend2.and. &
         LJ.ge.jstart2.and.LJ.le.jend2.and. &
         LK.ge.kstart2.and.LK.le.kend2
  END FUNCTION IS_ON_myPE_plus1layer

  LOGICAL FUNCTION IS_ON_myPE_owns (LI, LJ, LK)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: LI, LJ, LK

    IS_ON_myPE_owns = LI.ge.istart.and.LI.le.iend.and. &
         LJ.ge.jstart.and.LJ.le.jend.and. &
         LK.ge.kstart.and.LK.le.kend
  END FUNCTION IS_ON_myPE_owns

  !      IS_ON_myPE_owns (LI, LJ, LK) = LI.ge.istart.and.LI.le.iend.and. &
  !                               LJ.ge.jstart.and.LJ.le.jend.and. &
  !                                LK.ge.kstart.and.LK.le.kend.and. &
  !                                (.NOT.DEAD_CELL_AT(LI,LJ,LK))


  LOGICAL FUNCTION IS_ON_myPE_wobnd (LI, LJ, LK)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: LI, LJ, LK

    IS_ON_myPE_wobnd = LI.ge.istart1.and.LI.le.iend1.and. &
         LJ.ge.jstart1.and.LJ.le.jend1.and. &
         LK.ge.kstart1.and.LK.le.kend1
  END FUNCTION IS_ON_myPE_wobnd

  !      WEST_OF  (IJK)   = IJK + INCREMENT_FOR_w (CELL_CLASS(IJK))
  !      EAST_OF  (IJK)   = IJK + INCREMENT_FOR_e (CELL_CLASS(IJK))
  !      SOUTH_OF (IJK)   = IJK + INCREMENT_FOR_s (CELL_CLASS(IJK))
  !      NORTH_OF (IJK)   = IJK + INCREMENT_FOR_n (CELL_CLASS(IJK))
  !      BOTTOM_OF(IJK)   = IJK + INCREMENT_FOR_b (CELL_CLASS(IJK))
  !      TOP_OF   (IJK)   = IJK + INCREMENT_FOR_t (CELL_CLASS(IJK))
  !      IM_OF  (IJK)     = IJK + INCREMENT_FOR_im(CELL_CLASS(IJK))
  !      IP_OF  (IJK)     = IJK + INCREMENT_FOR_ip(CELL_CLASS(IJK))
  !      JM_OF (IJK)      = IJK + INCREMENT_FOR_jm(CELL_CLASS(IJK))
  !      JP_OF (IJK)      = IJK + INCREMENT_FOR_jp(CELL_CLASS(IJK))
  !      KM_OF(IJK)       = IJK + INCREMENT_FOR_km(CELL_CLASS(IJK))
  !      KP_OF   (IJK)    = IJK + INCREMENT_FOR_kp(CELL_CLASS(IJK))

  INTEGER FUNCTION EAST_OF  (IJK)
    IMPLICIT NONE
    INTEGER IJK

    EAST_OF = IJK + INCREMENT_FOR_NB (1,CELL_CLASS(IJK))
  END FUNCTION EAST_OF

  INTEGER FUNCTION WEST_OF  (IJK)
    IMPLICIT NONE
    INTEGER IJK
    WEST_OF = IJK + INCREMENT_FOR_NB (2,CELL_CLASS(IJK))
  END FUNCTION WEST_OF

  INTEGER FUNCTION SOUTH_OF (IJK)
    IMPLICIT NONE
    INTEGER IJK
    SOUTH_OF = IJK + INCREMENT_FOR_NB (3,CELL_CLASS(IJK))
  END FUNCTION SOUTH_OF

  INTEGER FUNCTION NORTH_OF (IJK)
    IMPLICIT NONE
    INTEGER IJK
    NORTH_OF = IJK + INCREMENT_FOR_NB (4,CELL_CLASS(IJK))
  END FUNCTION NORTH_OF

  INTEGER FUNCTION BOTTOM_OF(IJK)
    IMPLICIT NONE
    INTEGER IJK
    BOTTOM_OF = IJK + INCREMENT_FOR_NB (5,CELL_CLASS(IJK))
  END FUNCTION BOTTOM_OF

  INTEGER FUNCTION TOP_OF   (IJK)
    IMPLICIT NONE
    INTEGER IJK
    TOP_OF = IJK + INCREMENT_FOR_NB (6,CELL_CLASS(IJK))
  END FUNCTION TOP_OF


  INTEGER FUNCTION IM_OF  (IJK)
    IMPLICIT NONE
    INTEGER IJK
    IM_OF = IJK + INCREMENT_FOR_MP(1,CELL_CLASS(IJK))
  END FUNCTION IM_OF

  INTEGER FUNCTION IP_OF  (IJK)
    IMPLICIT NONE
    INTEGER IJK
    IP_OF = IJK + INCREMENT_FOR_MP(2,CELL_CLASS(IJK))
  END FUNCTION IP_OF

  INTEGER FUNCTION JM_OF (IJK)
    IMPLICIT NONE
    INTEGER IJK
    JM_OF = IJK + INCREMENT_FOR_MP(3,CELL_CLASS(IJK))
  END FUNCTION JM_OF

  INTEGER FUNCTION JP_OF (IJK)
    IMPLICIT NONE
    INTEGER IJK

    JP_OF = IJK + INCREMENT_FOR_MP(4,CELL_CLASS(IJK))
  END FUNCTION JP_OF

  INTEGER FUNCTION KM_OF(IJK)
    IMPLICIT NONE
    INTEGER IJK
    KM_OF = IJK + INCREMENT_FOR_MP(5,CELL_CLASS(IJK))
  END FUNCTION KM_OF

  INTEGER FUNCTION KP_OF   (IJK)
    IMPLICIT NONE
    INTEGER IJK
    KP_OF = IJK + INCREMENT_FOR_MP(6,CELL_CLASS(IJK))
  END FUNCTION KP_OF


  INTEGER FUNCTION WEST_OF_0  (IJK)
    IMPLICIT NONE
    INTEGER IJK
    WEST_OF_0 = IJK + INCREMENT_FOR_w (CELL_CLASS(IJK))
  END FUNCTION WEST_OF_0

  INTEGER FUNCTION EAST_OF_0  (IJK)
    IMPLICIT NONE
    INTEGER IJK
    EAST_OF_0 = IJK + INCREMENT_FOR_e (CELL_CLASS(IJK))
  END FUNCTION EAST_OF_0

  INTEGER FUNCTION SOUTH_OF_0 (IJK)
    IMPLICIT NONE
    INTEGER IJK
    SOUTH_OF_0 = IJK + INCREMENT_FOR_s (CELL_CLASS(IJK))
  END FUNCTION SOUTH_OF_0

  INTEGER FUNCTION NORTH_OF_0 (IJK)
    IMPLICIT NONE
    INTEGER IJK
    NORTH_OF_0 = IJK + INCREMENT_FOR_n (CELL_CLASS(IJK))
  END FUNCTION NORTH_OF_0

  INTEGER FUNCTION BOTTOM_OF_0(IJK)
    IMPLICIT NONE
    INTEGER IJK
    BOTTOM_OF_0 = IJK + INCREMENT_FOR_b (CELL_CLASS(IJK))
  END FUNCTION BOTTOM_OF_0

  INTEGER FUNCTION TOP_OF_0   (IJK)
    IMPLICIT NONE
    INTEGER IJK
    TOP_OF_0 = IJK + INCREMENT_FOR_t (CELL_CLASS(IJK))
  END FUNCTION TOP_OF_0

  INTEGER FUNCTION IM_OF_0  (IJK)
    IMPLICIT NONE
    INTEGER IJK

    IM_OF_0 = IJK + INCREMENT_FOR_im(CELL_CLASS(IJK))
  END FUNCTION IM_OF_0

  INTEGER FUNCTION IP_OF_0  (IJK)
    IMPLICIT NONE
    INTEGER IJK
    IP_OF_0 = IJK + INCREMENT_FOR_ip(CELL_CLASS(IJK))
  END FUNCTION IP_OF_0

  INTEGER FUNCTION JM_OF_0 (IJK)
    IMPLICIT NONE
    INTEGER IJK
    JM_OF_0 = IJK + INCREMENT_FOR_jm(CELL_CLASS(IJK))
  END FUNCTION JM_OF_0

  INTEGER FUNCTION JP_OF_0 (IJK)
    IMPLICIT NONE
    INTEGER IJK
    JP_OF_0 = IJK + INCREMENT_FOR_jp(CELL_CLASS(IJK))
  END FUNCTION JP_OF_0

  INTEGER FUNCTION KM_OF_0(IJK)
    IMPLICIT NONE
    INTEGER IJK
    KM_OF_0 = IJK + INCREMENT_FOR_km(CELL_CLASS(IJK))
  END FUNCTION KM_OF_0

  INTEGER FUNCTION KP_OF_0   (IJK)
    IMPLICIT NONE
    INTEGER IJK
    KP_OF_0 = IJK + INCREMENT_FOR_kp(CELL_CLASS(IJK))
  END FUNCTION KP_OF_0


  !      WEST_OF  (IJK)   = WEST_ARRAY_OF(IJK)
  !      EAST_OF  (IJK)   = EAST_ARRAY_OF(IJK)
  !      SOUTH_OF  (IJK)   = SOUTH_ARRAY_OF(IJK)
  !      NORTH_OF  (IJK)   = NORTH_ARRAY_OF(IJK)
  !      BOTTOM_OF  (IJK)   = BOTTOM_ARRAY_OF(IJK)
  !      TOP_OF  (IJK)   = TOP_ARRAY_OF(IJK)
  !      IM_OF  (IJK)   = IM_ARRAY_OF(IJK)
  !      IP_OF  (IJK)   = IP_ARRAY_OF(IJK)
  !      JM_OF  (IJK)   = JM_ARRAY_OF(IJK)
  !      JP_OF  (IJK)   = JP_ARRAY_OF(IJK)
  !      KM_OF  (IJK)   = KM_ARRAY_OF(IJK)
  !      KP_OF  (IJK)   = KP_ARRAY_OF(IJK)

  INTEGER FUNCTION FUNLM (L1, L2)
    IMPLICIT NONE
    INTEGER L1, L2
    FUNLM = STORE_LM (L1, L2)
  END FUNCTION FUNLM

  LOGICAL FUNCTION FLUID_AT(IJK)
    INTEGER IJK
    FLUID_AT    = FLAG(IJK) .EQ. 1
  END FUNCTION FLUID_AT


  LOGICAL FUNCTION CYCLIC_AT_E(IJK)
    INTEGER IJK
    CYCLIC_AT_E   = FLAG_E(IJK) .EQ. 2000
  END FUNCTION CYCLIC_AT_E

  LOGICAL FUNCTION CYCLIC_AT_N(IJK)
    INTEGER IJK
    CYCLIC_AT_N   = FLAG_N(IJK) .EQ. 2000
  END FUNCTION CYCLIC_AT_N

  LOGICAL FUNCTION CYCLIC_AT_T(IJK)
    INTEGER IJK
    CYCLIC_AT_T   = FLAG_T(IJK) .EQ. 2000
  END FUNCTION CYCLIC_AT_T

  LOGICAL FUNCTION FLOW_AT_E(IJK)
    INTEGER IJK
    FLOW_AT_E   = FLAG_E(IJK) .GE. 2000 .AND.&
         FLAG_E(IJK) .LE. 2011
  END FUNCTION FLOW_AT_E

  LOGICAL FUNCTION FLOW_AT_N(IJK)
    INTEGER IJK
    FLOW_AT_N   = FLAG_N(IJK) .GE. 2000 .AND.&
         FLAG_N(IJK) .LE. 2011
  END FUNCTION FLOW_AT_N

  LOGICAL FUNCTION FLOW_AT_T(IJK)
    INTEGER IJK
    FLOW_AT_T   = FLAG_T(IJK) .GE. 2000 .AND.&
         FLAG_T(IJK) .LE. 2011
  END FUNCTION FLOW_AT_T

  LOGICAL FUNCTION PFLOW_AT_E(IJK)
    INTEGER IJK
    PFLOW_AT_E  = FLAG_E(IJK) .EQ. 2010 .OR.&
         FLAG_E(IJK) .EQ. 2011
  END FUNCTION PFLOW_AT_E

  LOGICAL FUNCTION PFLOW_AT_N(IJK)
    INTEGER IJK
    PFLOW_AT_N  = FLAG_N(IJK) .EQ. 2010 .OR.&
         FLAG_N(IJK) .EQ. 2011
  END FUNCTION PFLOW_AT_N

  LOGICAL FUNCTION PFLOW_AT_T(IJK)
    INTEGER IJK
    PFLOW_AT_T  = FLAG_T(IJK) .EQ. 2010 .OR.&
         FLAG_T(IJK) .EQ. 2011
  END FUNCTION PFLOW_AT_T

  LOGICAL FUNCTION MFLOW_AT_E(IJK)
    INTEGER IJK
    MFLOW_AT_E  = FLAG_E(IJK) .EQ. 2020 .OR. &
         FLAG_E(IJK) .EQ. 2021 .OR. &
         FLAG_E(IJK) .EQ. 2031
  END FUNCTION MFLOW_AT_E

  LOGICAL FUNCTION MFLOW_AT_N(IJK)
    INTEGER IJK
    MFLOW_AT_N  = FLAG_N(IJK) .EQ. 2020 .OR. &
         FLAG_N(IJK) .EQ. 2021 .OR. &
         FLAG_N(IJK) .EQ. 2031
  END FUNCTION MFLOW_AT_N

  LOGICAL FUNCTION MFLOW_AT_T(IJK)
    INTEGER IJK
    MFLOW_AT_T  = FLAG_T(IJK) .EQ. 2020 .OR. &
         FLAG_T(IJK) .EQ. 2021 .OR. &
         FLAG_T(IJK) .EQ. 2031
  END FUNCTION MFLOW_AT_T

  LOGICAL FUNCTION FLUIDorP_FLOW_AT(IJK)
    INTEGER IJK
    FLUIDorP_FLOW_AT = FLAG(IJK) .LE. 11
  END FUNCTION FLUIDorP_FLOW_AT

  LOGICAL FUNCTION P_FLOW_AT(IJK)
    INTEGER IJK
    P_FLOW_AT = FLAG(IJK) .EQ. 10 .OR. &
         FLAG(IJK) .EQ. 11
  END FUNCTION P_FLOW_AT

  LOGICAL FUNCTION P_OUTFLOW_AT(IJK)
    INTEGER IJK
    P_OUTFLOW_AT= FLAG(IJK) .EQ. 11
  END FUNCTION P_OUTFLOW_AT

  LOGICAL FUNCTION MASS_OUTFLOW_AT(IJK)
    INTEGER IJK
    MASS_OUTFLOW_AT= FLAG(IJK) .EQ. 21
  END FUNCTION MASS_OUTFLOW_AT

  LOGICAL FUNCTION OUTFLOW_AT(IJK)
    INTEGER IJK
    OUTFLOW_AT  = FLAG(IJK) .EQ. 31
  END FUNCTION OUTFLOW_AT

  LOGICAL FUNCTION FLOW_AT(IJK)
    INTEGER IJK
    FLOW_AT     = FLAG(IJK) .GE. 10 .AND. FLAG(IJK) .LE. 31
  END FUNCTION FLOW_AT

  LOGICAL FUNCTION NS_WALL_AT(IJK)
    INTEGER IJK
    NS_WALL_AT  = FLAG(IJK) .EQ. 100
  END FUNCTION NS_WALL_AT

  LOGICAL FUNCTION FS_WALL_AT(IJK)
    INTEGER IJK
    FS_WALL_AT  = FLAG(IJK) .EQ. 101
  END FUNCTION FS_WALL_AT

  LOGICAL FUNCTION PS_WALL_AT(IJK)
    INTEGER IJK
    PS_WALL_AT  = FLAG(IJK) .EQ. 102
  END FUNCTION PS_WALL_AT

  LOGICAL FUNCTION CYCLIC_AT(IJK)
    INTEGER IJK
    CYCLIC_AT = FLAG(IJK) .EQ. 106 .OR. &
         FLAG(IJK) .EQ. 107
  END FUNCTION CYCLIC_AT

  LOGICAL FUNCTION WALL_AT(IJK)
    INTEGER IJK
    WALL_AT     = FLAG(IJK) .GE. 100
  END FUNCTION WALL_AT

  LOGICAL FUNCTION IP_AT_E(IJK)
    INTEGER IJK
    IP_AT_E     = FLAG_E(IJK) .LT. 1000
  END FUNCTION IP_AT_E

  LOGICAL FUNCTION IP_AT_N(IJK)
    INTEGER IJK
    IP_AT_N     = FLAG_N(IJK) .LT. 1000
  END FUNCTION IP_AT_N

  LOGICAL FUNCTION IP_AT_T(IJK)
    INTEGER IJK
    IP_AT_T     = FLAG_T(IJK) .LT. 1000
  END FUNCTION IP_AT_T

  LOGICAL FUNCTION SIP_AT_E(IJK)
    INTEGER IJK
    SIP_AT_E    = (FLAG_E(IJK) .LT. 2000)
  END FUNCTION SIP_AT_E

  LOGICAL FUNCTION SIP_AT_N(IJK)
    INTEGER IJK
    SIP_AT_N    = (FLAG_N(IJK) .LT. 2000)
  END FUNCTION SIP_AT_N

  LOGICAL FUNCTION SIP_AT_T(IJK)
    INTEGER IJK
    SIP_AT_T    = (FLAG_T(IJK) .LT. 2000)
  END FUNCTION SIP_AT_T

  LOGICAL FUNCTION SP_AT_E(IJK)
    INTEGER IJK
    SP_AT_E     = (FLAG_E(IJK) .LT. 2000) .AND. &
         (FLAG_E(IJK) .GE. 1000)
  END FUNCTION SP_AT_E

  LOGICAL FUNCTION SP_AT_N(IJK)
    INTEGER IJK
    SP_AT_N     = (FLAG_N(IJK) .LT. 2000) .AND. &
         (FLAG_N(IJK) .GE. 1000)
  END FUNCTION SP_AT_N

  LOGICAL FUNCTION SP_AT_T(IJK)
    INTEGER IJK
    SP_AT_T     = (FLAG_T(IJK) .LT. 2000) .AND. &
         (FLAG_T(IJK) .GE. 1000)
  END FUNCTION SP_AT_T

  INTEGER FUNCTION IS_ID_AT_E(IJK)
    IMPLICIT NONE
    INTEGER IJK
    IS_ID_AT_E = FLAG_E(IJK) - 1000
  END FUNCTION IS_ID_AT_E

  INTEGER FUNCTION IS_ID_AT_N(IJK)
    IMPLICIT NONE
    INTEGER IJK
    IS_ID_AT_N  = FLAG_N(IJK) - 1000
  END FUNCTION IS_ID_AT_N

  INTEGER FUNCTION IS_ID_AT_T(IJK)
    IMPLICIT NONE
    INTEGER IJK
    IS_ID_AT_T  = FLAG_T(IJK) - 1000
  END FUNCTION IS_ID_AT_T

  LOGICAL FUNCTION IS_AT_E(IJK)
    IMPLICIT NONE
    INTEGER IJK
    IS_AT_E     = FLAG_E(IJK) .LT. 2000
  END FUNCTION IS_AT_E

  LOGICAL FUNCTION IS_AT_N(IJK)
    IMPLICIT NONE
    INTEGER IJK
    IS_AT_N     = FLAG_N(IJK) .LT. 2000
  END FUNCTION IS_AT_N

  LOGICAL FUNCTION IS_AT_T(IJK)
    IMPLICIT NONE
    INTEGER IJK
    IS_AT_T     = FLAG_T(IJK) .LT. 2000
  END FUNCTION IS_AT_T

  LOGICAL FUNCTION NO_IS_AT_E(IJK)
    IMPLICIT NONE
    INTEGER IJK
    NO_IS_AT_E  = FLAG_E(IJK) .GE. 2000
  END FUNCTION NO_IS_AT_E

  LOGICAL FUNCTION NO_IS_AT_N(IJK)
    IMPLICIT NONE
    INTEGER IJK
    NO_IS_AT_N  = FLAG_N(IJK) .GE. 2000
  END FUNCTION NO_IS_AT_N

  LOGICAL FUNCTION NO_IS_AT_T(IJK)
    IMPLICIT NONE
    INTEGER IJK
    NO_IS_AT_T  = FLAG_T(IJK) .GE. 2000
  END FUNCTION NO_IS_AT_T

  LOGICAL FUNCTION WALL_ICBC_FLAG(IJK)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: IJK
    WALL_ICBC_FLAG = ICBC_FLAG(IJK)(1:1) .EQ. 'W' .OR. &
         ICBC_FLAG(IJK)(1:1) .EQ. 'S' .OR. &
         ICBC_FLAG(IJK)(1:1) .EQ. 's' .OR. &
         ICBC_FLAG(IJK)(1:1) .EQ. 'c' .OR. &
         ICBC_FLAG(IJK)(1:1) .EQ. 'C'
  END FUNCTION WALL_ICBC_FLAG

  LOGICAL FUNCTION DEFAULT_WALL_AT(IJK)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: IJK
    DEFAULT_WALL_AT = ICBC_FLAG(IJK)(2:3) .EQ. '--' .AND. &
         (ICBC_FLAG(IJK)(1:1) .NE. 'c'  .AND. &
         ICBC_FLAG(IJK)(1:1) .NE. 'C')
  END FUNCTION DEFAULT_WALL_AT

  DOUBLE PRECISION FUNCTION ZMAX(XXX)
    IMPLICIT NONE
    DOUBLE PRECISION XXX
    ZMAX       = MAX(XXX, ZERO)
  END FUNCTION ZMAX

END MODULE functions
