!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                         C 
!     Module name: DES_INIT_NAMELIST                                      C
!     Purpose: DES - initialize the des-namelist                          C
!                                                                         C
!                                                                         C
!     Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!     Reviewer:                                          Date:            C
!                                                                         C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!     
      SUBROUTINE DES_INIT_NAMELIST 

      USE discretelement
      
      IMPLICIT NONE
!-----------------------------------------------
!     G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!     L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!     L o c a l   V a r i a b l e s
!-----------------------------------------------
!     loop counters
      INTEGER :: LC, LCM, M, N, K, KKK
!     
!     Coefficient of restitution (old symbol)
      DOUBLE PRECISION :: E
!-----------------------------------------------
!     
!     
      
      INCLUDE 'des/desnamelist.inc'
      
      
                  DISCRETE_ELEMENT = .FALSE.
                  DO_NBS = .FALSE.
                  DO_QUADTREE = .FALSE.
                  DO_OCTREE = .FALSE.
                  DO_NSQUARE = .FALSE.
                  WALLFIXEDOVERLAP = .FALSE.
                  WALLDTSPLIT = .FALSE.
                  WALLREFLECT = .FALSE.
                  DES_PERIODIC_WALLS = .FALSE.
                  DES_PERIODIC_WALLS_X = .FALSE.
                  DES_PERIODIC_WALLS_Y = .FALSE.
                  DES_PERIODIC_WALLS_Z = .FALSE.
                  EQUIVALENT_RADIUS = .FALSE.
                  INLET_OUTLET = .FALSE.
                  TIME_ADJUST  = .TRUE.
                  DES_CONTINUUM_COUPLED = .FALSE.
                  COUPLED_FLOW = .FALSE.
                  TSUJI_DRAG = .FALSE.
                  KUIPERS_DRAG = .FALSE.

                  IF(DISCRETE_ELEMENT) THEN
                    TIME_ADJUST = .FALSE.
                  END IF
                  
                  DO LC = 1, NPARTICLES 
                     DES_RADIUS(LC) = 0D0
                     PMASS(LC) = 0D0
                     MOI(LC) = 0D0
                     RO_Sol(LC) = 0.D0
                     DO K = 1, NDIM
                        DES_POS_OLD(K,LC) = 0D0
                        DES_POS_NEW(K,LC) = 0D0
                        DES_VEL_OLD(K,LC) = 0D0
                        DES_VEL_NEW(K,LC) = 0D0
                        FC(K,LC) = 0D0
                        FN(K,LC) = 0D0
                        FT(K,LC) = 0D0
                        TOW(K,LC) = 0D0
                        OMEGA_OLD(K,LC) = 0D0
                        OMEGA_NEW(K,LC) = 0D0
                     END DO
                  END DO

                  DO LC = 1, MAXQUADS
                     DO K = 1, NMQD 
                        LQUAD(K,LC) = UNDEFINED_I
                     END DO
                     DO K = 1, NWALLS
                        CQUAD(K,LC) = UNDEFINED
                     END DO
                  END DO

                  DO K= 1, NPARTICLES
                     DO KKK = 2, MAXNEIGHBORS
                        NEIGHBOURS(KKK,K) = -1
                        PN(KKK,K) = -1
                        PV(KKK,K) = 1
                        DO LC = 1, NDIM
                           PFN(LC,KKK,K) = 0.0
                           PFT(LC,KKK,K) = 0.0
                        END DO
                     END DO
                     NEIGHBOURS(1,K) = 0
                     PN(1,K) = 0
                     PV(1,K) = 1
                  END DO

                  DO K = 1, NWALLS
                     DO KKK = 1,NDIM
                        DES_WALL_POS(KKK,K) = UNDEFINED
                        DES_WALL_VEL(KKK,K) = UNDEFINED
                        WALL_NORMAL(KKK,K) = UNDEFINED
                     END DO
                  END DO

                  DO K = 1,NDIM
                     FNS1(K) = UNDEFINED
                     FTS1(K) = UNDEFINED
                     GRAV(K) = UNDEFINED
                  END DO

                  KN = UNDEFINED
                  KT = UNDEFINED
                  KN_W = UNDEFINED
                  KT_W = UNDEFINED
                  ETA_DES_N = UNDEFINED
                  ETA_DES_T = UNDEFINED
                  ETA_N_W = UNDEFINED
                  ETA_T_W = UNDEFINED
                  MEW = UNDEFINED
                  MEW_W = UNDEFINED
                  E_RESTITUTION = UNDEFINED
                  DES_GAMMA = UNDEFINED
                  DES_F = UNDEFINED
                  DES_KE = UNDEFINED
                  DES_PE = UNDEFINED
                  CALLED = 0D0
                  PCOUNT = 0D0
                  PARTICLES = UNDEFINED_I
                  OUTOFBOX = UNDEFINED_I
                  DTSOLID = UNDEFINED
                  FACTOR = UNDEFINED_I
                  ROs = UNDEFINED
                  ZONES = 1
                  ZN1 = 1
                  ZN2 = 0
                  ZN3 = 0
                  ZN4 = 0

                  N2CT = UNDEFINED
                  NBSCT = UNDEFINED
                  QUADCT = UNDEFINED
                  OCTCT = UNDEFINED

                  WX1 = UNDEFINED
                  EX2 = UNDEFINED
                  BY1 = UNDEFINED
                  TY2 = UNDEFINED
                  SZ1 = UNDEFINED
                  NZ2 = UNDEFINED
                  RADIUS_EQ = UNDEFINED
                  NQUAD = UNDEFINED_I
                  DIMN = UNDEFINED_I

                  RETURN
                  END SUBROUTINE DES_INIT_NAMELIST
