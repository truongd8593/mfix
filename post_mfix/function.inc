!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
!                      Function for generating the LOCAL 3-D array index IJK
!                      from the 1-D indices I, J, and K.
!//FUNIJK is moved to compar for debugging purposes - Sreekanth-10/26/99
      INTEGER          FUNIJK
!                      Function for generating the LOCAL 3-D array index IJK
!                      from the 1-D indices I, J, K and IPROC.
      INTEGER          FUNIJK_PROC
!                      Function for generating the GLOBAL 3-D array index IJK
!                      from the 1-D indices I, J, and K.
      INTEGER          FUNIJK_GL
!
!                      Function for generating the index for the entries
!                      to the upper triangle (excluding the diagonal) of
!                      an (L,M) matrix.
      INTEGER          FUNLM
!
!                      Function for generating the 3-D array index IJK
!                      from the 1-D indices I, J, and K in IO format
      INTEGER          FUNIJK_IO
!
!                      Dummy arguments
      DOUBLE PRECISION XXX, YYY, XXX1, XXX2, XXX3, XXX4, XXX5
!
!                      Dummy index for L and M
      INTEGER          L1, L2

!                      Dummy indices for I, J, K
      INTEGER          LI, LJ, LK
!                      Dummy indices for IPROC
      INTEGER          LIPROC
!
!                      Function for calculating IJKW
      INTEGER          WEST_OF
!
!                      Function for calculating IJKE
      INTEGER          EAST_OF
!
!                      Function for calculating IJKS
      INTEGER          SOUTH_OF
!
!                      Function for calculating IJKN
      INTEGER          NORTH_OF
!
!                      Function for calculating IJKB
      INTEGER          BOTTOM_OF
!
!                      Function for calculating IJKT
      INTEGER          TOP_OF
!
!                      Function for calculating IMJK
      INTEGER          IM_OF
!
!                      Function for calculating IPJK
      INTEGER          IP_OF
!
!                      Function for calculating IJMK
      INTEGER          JM_OF
!
!                      Function for calculating IJPK
      INTEGER          JP_OF
!
!                      Function for calculating IJKM
      INTEGER          KM_OF
!
!                      Function for calculating IJKP
      INTEGER          KP_OF
!
!                      Logical function to identify a fluid cell
      LOGICAL          FLUID_AT
!
!                      cyclic condition at east boundary
      LOGICAL          CYCLIC_AT_E
!
!                      cyclic condition at north boundary
      LOGICAL          CYCLIC_AT_N
!
!                      cyclic condition at top boundary
      LOGICAL          CYCLIC_AT_T
!
!                      identify flow at east boundary
      LOGICAL          FLOW_AT_E
!
!                      identify flow at north boundary
      LOGICAL          FLOW_AT_N
!
!                      identify flow at top boundary
      LOGICAL          FLOW_AT_T
!
!                      identify const. pressure flow east boundary
      LOGICAL          PFLOW_AT_E
!
!                      identify const. pressure flow north boundary
      LOGICAL          PFLOW_AT_N
!
!                      identify const. pressure flow top boundary
      LOGICAL          PFLOW_AT_T
!
!                      identify specified flow east boundary
      LOGICAL          MFLOW_AT_E
!
!                      identify specified flow north boundary
      LOGICAL          MFLOW_AT_N
!
!                      identify specified flow top boundary
      LOGICAL          MFLOW_AT_T
!
!                      FLUID_AT or P_FLOW_AT (simplified check)
      LOGICAL          FLUIDorP_FLOW_AT
!
!                      P_FLOW_AT 
      LOGICAL          P_FLOW_AT
!
!                      P_OUTFLOW_AT 
      LOGICAL          P_OUTFLOW_AT
!
!                      MASS_OUTFLOW_AT 
      LOGICAL          MASS_OUTFLOW_AT
!
!                      OUTFLOW_AT 
      LOGICAL          OUTFLOW_AT
!
!                      FLOW_AT 
      LOGICAL          FLOW_AT
!
!                      Logical function to identify default walls
      LOGICAL          DEFAULT_WALL_AT
!
!                      Logical function to identify wall ICBC_FLAG
      LOGICAL          WALL_ICBC_FLAG
!
!                      Logical function to identify a No-slip wall cell
      LOGICAL          NS_WALL_AT
!
!                      Logical function to identify a Free-slip wall cell
      LOGICAL          FS_WALL_AT
!
!                      Logical function to identify a Partial-slip wall cell
      LOGICAL          PS_WALL_AT
!
!                      Logical function to identify a wall cell
      LOGICAL          WALL_AT
!
!                      Logical function to identify a cyclic cell
      LOGICAL          CYCLIC_AT
!
!                      Logical function to identify IP at East of the cell
      LOGICAL          IP_AT_E
!
!                      Logical function to identify IP at North of the cell
      LOGICAL          IP_AT_N
!
!                      Logical function to identify IP at Top of the cell
      LOGICAL          IP_AT_T
!
!                      Logical function to identify IP or SP at East of the cell
      LOGICAL          SIP_AT_E
!
!                      Logical function to identify IP or SP at North of the cell
      LOGICAL          SIP_AT_N
!
!                      Logical function to identify IP or SP at Top of the cell
      LOGICAL          SIP_AT_T
!
!                      Logical function to identify SP at East of the cell
      LOGICAL          SP_AT_E
!
!                      Logical function to identify SP at North of the cell
      LOGICAL          SP_AT_N
!
!                      Logical function to identify SP at Top of the cell
      LOGICAL          SP_AT_T
!
!                      Internal surface ID for east face
      INTEGER          IS_ID_AT_E
!
!                      Internal surface ID for north face
      INTEGER          IS_ID_AT_N
!
!                      Internal surface ID for top face
      INTEGER          IS_ID_AT_T
!
!                      Logical function to identify IS at East of the cell
      LOGICAL          IS_AT_E
!
!                      Logical function to identify IS at North of the cell
      LOGICAL          IS_AT_N
!
!                      Logical function to identify IS at Top of the cell
      LOGICAL          IS_AT_T
!
!                      Logical function to identify No IS at East of the cell
      LOGICAL          NO_IS_AT_E
!
!                      Logical function to identify No IS at North of the cell
      LOGICAL          NO_IS_AT_N
!
!                      Logical function to identify No IS at Top of the cell
      LOGICAL          NO_IS_AT_T
!\\   Added to determine whether I am on my PE  
!\\   Logical function to determine whether I am on my PE's domain 
      LOGICAL          IS_ON_myPE_owns
!\\   Logical function to determine whether I am on my PE's domain + one &
!                                                       extra ghost layer
      LOGICAL          IS_ON_myPE_plus1layer
!\\   Logical function to determine whether I am on my PE's domain + two &
!                                                       extra ghost layers
      LOGICAL          IS_ON_myPE_plus2layers
!\\   Logical function to determine whether I am on my PE's domain - Boundaries
      LOGICAL          IS_ON_myPE_wobnd
!
!                      Maximum of zero or input
      DOUBLE PRECISION ZMAX
!
! 
!//FUNIJK is moved to compar for debugging purposes - Sreekanth-10/26/99
!     FUNIJK (LI, LJ, LK) = c0 + LI + (LJ-jstart3_all(myPE))*c1 + (LK-kstart3_all(myPE))* c2
      funijk(li,lj,lk) = lj + c0 + li*c1 + lk*c2

!//SP
!     FUNIJK_PROC(LI, LJ, LK, LIPROC) = 1 + (LI - istart3_all(LIPROC))+ &
!     (LJ-jstart3_all(LIPROC))*(iend3_all(LIPROC)-istart3_all(LIPROC)+1) &
!     + (LK-kstart3_all(LIPROC))*(jend3_all(LIPROC)-jstart3_all(LIPROC)+1)* &
!     (iend3_all(LIPROC)-istart3_all(LIPROC)+1) 

      FUNIJK_PROC(LI, LJ, LK, LIPROC) = 1 + (LJ - jstart3_all(LIPROC))+ &
      (LI-Istart3_all(LIPROC))*(jend3_all(LIPROC)-jstart3_all(LIPROC)+1) &
      + (LK-kstart3_all(LIPROC))*(jend3_all(LIPROC)-jstart3_all(LIPROC)+1)* &
      (iend3_all(LIPROC)-istart3_all(LIPROC)+1) 

!     FUNIJK_GL (LI, LJ, LK) = 1 + (LI - imin3) + (LJ-jmin3)*(imax3-imin3+1) &
!     + (LK-kmin3)*(jmax3-jmin3+1)*(imax3-imin3+1)
!
      FUNIJK_GL (LI, LJ, LK) = 1 + (LJ - jmin3) + (LI-imin3)*(jmax3-jmin3+1) &
      + (LK-kmin3)*(jmax3-jmin3+1)*(imax3-imin3+1)

      FUNIJK_IO (LI, LJ, LK) = 1 + (LI - imin2) + (LJ-jmin2)*(imax2-imin2+1) &
      + (LK-kmin2)*(jmax2-jmin2+1)*(imax2-imin2+1)

      IS_ON_myPE_plus2layers (LI, LJ, LK) = LI.ge.istart3.and.LI.le.iend3.and. &
				LJ.ge.jstart3.and.LJ.le.jend3.and. &
                                LK.ge.kstart3.and.LK.le.kend3

      IS_ON_myPE_plus1layer (LI, LJ, LK) = LI.ge.istart2.and.LI.le.iend2.and. & 
      					   LJ.ge.jstart2.and.LJ.le.jend2.and. & 
      					   LK.ge.kstart2.and.LK.le.kend2

      IS_ON_myPE_owns (LI, LJ, LK) = LI.ge.istart.and.LI.le.iend.and. &
				LJ.ge.jstart.and.LJ.le.jend.and. &
                                LK.ge.kstart.and.LK.le.kend

      IS_ON_myPE_wobnd (LI, LJ, LK) = LI.ge.istart1.and.LI.le.iend1.and. &
                                LJ.ge.jstart1.and.LJ.le.jend1.and. &
                                LK.ge.kstart1.and.LK.le.kend1

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

      FLOW_AT_E(IJK)   = FLAG_E(IJK) .GE. 2000 .AND.&
                        FLAG_E(IJK) .LE. 2011
      FLOW_AT_N(IJK)   = FLAG_N(IJK) .GE. 2000 .AND.&
                        FLAG_N(IJK) .LE. 2011
      FLOW_AT_T(IJK)   = FLAG_T(IJK) .GE. 2000 .AND.&
                        FLAG_T(IJK) .LE. 2011
      PFLOW_AT_E(IJK)  = FLAG_E(IJK) .EQ. 2010 .OR.&
                        FLAG_E(IJK) .EQ. 2011
      PFLOW_AT_N(IJK)  = FLAG_N(IJK) .EQ. 2010 .OR.&
                        FLAG_N(IJK) .EQ. 2011
      PFLOW_AT_T(IJK)  = FLAG_T(IJK) .EQ. 2010 .OR.&
                        FLAG_T(IJK) .EQ. 2011
      MFLOW_AT_E(IJK)  = FLAG_E(IJK) .EQ. 2020 .OR. &
                        FLAG_E(IJK) .EQ. 2021 .OR. &
                        FLAG_E(IJK) .EQ. 2031
      MFLOW_AT_N(IJK)  = FLAG_N(IJK) .EQ. 2020 .OR. &
                        FLAG_N(IJK) .EQ. 2021 .OR. &
                        FLAG_N(IJK) .EQ. 2031
      MFLOW_AT_T(IJK)  = FLAG_T(IJK) .EQ. 2020 .OR. &
                        FLAG_T(IJK) .EQ. 2021 .OR. &
                        FLAG_T(IJK) .EQ. 2031
      FLUIDorP_FLOW_AT(IJK)= FLAG(IJK) .LE. 11
      P_FLOW_AT(IJK)   = FLAG(IJK) .EQ. 10 .OR. &
                        FLAG(IJK) .EQ. 11
      P_OUTFLOW_AT(IJK)= FLAG(IJK) .EQ. 11
      MASS_OUTFLOW_AT(IJK)= FLAG(IJK) .EQ. 21
      OUTFLOW_AT(IJK)  = FLAG(IJK) .EQ. 31
      FLOW_AT(IJK)     = FLAG(IJK) .GE. 10 .AND. FLAG(IJK) .LE. 31
      NS_WALL_AT(IJK)  = FLAG(IJK) .EQ. 100
      FS_WALL_AT(IJK)  = FLAG(IJK) .EQ. 101
      PS_WALL_AT(IJK)  = FLAG(IJK) .EQ. 102
      CYCLIC_AT(IJK)   = FLAG(IJK) .EQ. 106 .OR. &
                        FLAG(IJK) .EQ. 107
      WALL_AT(IJK)     = FLAG(IJK) .GE. 100 
      IP_AT_E(IJK)     = FLAG_E(IJK) .LT. 1000
      IP_AT_N(IJK)     = FLAG_N(IJK) .LT. 1000
      IP_AT_T(IJK)     = FLAG_T(IJK) .LT. 1000
      SIP_AT_E(IJK)    = (FLAG_E(IJK) .LT. 2000)
      SIP_AT_N(IJK)    = (FLAG_N(IJK) .LT. 2000)
      SIP_AT_T(IJK)    = (FLAG_T(IJK) .LT. 2000)
      SP_AT_E(IJK)     = (FLAG_E(IJK) .LT. 2000) .AND. &
                        (FLAG_E(IJK) .GE. 1000)
      SP_AT_N(IJK)     = (FLAG_N(IJK) .LT. 2000) .AND. &
                        (FLAG_N(IJK) .GE. 1000)
      SP_AT_T(IJK)     = (FLAG_T(IJK) .LT. 2000) .AND. &
                        (FLAG_T(IJK) .GE. 1000)
      IS_ID_AT_E(IJK)  = FLAG_E(IJK) - 1000
      IS_ID_AT_N(IJK)  = FLAG_N(IJK) - 1000
      IS_ID_AT_T(IJK)  = FLAG_T(IJK) - 1000
      IS_AT_E(IJK)     = FLAG_E(IJK) .LT. 2000
      IS_AT_N(IJK)     = FLAG_N(IJK) .LT. 2000
      IS_AT_T(IJK)     = FLAG_T(IJK) .LT. 2000
      NO_IS_AT_E(IJK)  = FLAG_E(IJK) .GE. 2000
      NO_IS_AT_N(IJK)  = FLAG_N(IJK) .GE. 2000
      NO_IS_AT_T(IJK)  = FLAG_T(IJK) .GE. 2000
      WALL_ICBC_FLAG(IJK) = ICBC_FLAG(IJK)(1:1) .EQ. 'W' .OR. &
                           ICBC_FLAG(IJK)(1:1) .EQ. 'S' .OR. &
                           ICBC_FLAG(IJK)(1:1) .EQ. 's' .OR. &
                           ICBC_FLAG(IJK)(1:1) .EQ. 'c' .OR. &
                           ICBC_FLAG(IJK)(1:1) .EQ. 'C'
      DEFAULT_WALL_AT(IJK) = ICBC_FLAG(IJK)(2:3) .EQ. '--' .AND. &
                             (ICBC_FLAG(IJK)(1:1) .NE. 'c'  .AND. &
                             ICBC_FLAG(IJK)(1:1) .NE. 'C')
      
      ZMAX (XXX)       = MAX(XXX, ZERO)
