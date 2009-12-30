!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFWALLPOSVEL
!  Purpose:  DES -Calculate the position and velocity of wall particle 
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CFWALLPOSVEL(L, IW)

      USE discretelement
      USE des_bc
      USE param
      USE param1
      USE parallel
      USE fldvar
      USE run
      USE geometry
      USE matrix
      USE indices
      USE physprop
      USE drag
      USE constant
      USE compar
      IMPLICIT NONE

!-----------------------------------------------
! Local variables
!-----------------------------------------------
! index particle no. 
      INTEGER L
! index of wall no. (1-6)
      INTEGER IW
! local quantities for a vibrating wall 
      DOUBLE PRECISION A, OMEGA_W, OOMEGAW2, COSOMEGAT, SINOMEGAT
! local value for particle radius
      DOUBLE PRECISION DES_R

!-----------------------------------------------      
     
! initialize wall velocity
      DES_WALL_VEL(IW,1) = ZERO
      DES_WALL_VEL(IW,2) = ZERO
      IF(DIMN.EQ.3) DES_WALL_VEL(IW,3) = ZERO
     
! initialize local values
      A = ZERO
      OMEGA_W = ZERO
      SINOMEGAT = ZERO
      COSOMEGAT = ZERO

! if moving wall assign values accordingly
      IF(DES_F.NE.ZERO) THEN
         OMEGA_W = 2.0d0*PI*DES_F
         OOMEGAW2 = ONE/(OMEGA_W**2)
         A = DES_GAMMA*GRAV(2)*OOMEGAW2
         SINOMEGAT = SIN(OMEGA_W*S_TIME)
         COSOMEGAT = COS(OMEGA_W*S_TIME)
      END IF

      DES_R = DES_RADIUS(L)
      IF(WALLREFLECT) THEN
         DES_R = ZERO
      ENDIF

! Assigning wall position and velocity

! west wall; in 3D on yz plane (x=wx1->0)
      IF(IW.EQ.1) THEN
         DES_WALL_POS(IW,1) = WX1 - DES_R
         DES_WALL_POS(IW,2) = DES_POS_NEW(L,2)
         IF(DIMN.EQ.3) DES_WALL_POS(IW,3) = DES_POS_NEW(L,3)

         IF(DES_F.NE.ZERO) THEN
         ELSE                 
            DES_WALL_VEL(IW,1) = ZERO
            DES_WALL_VEL(IW,2) = DES_BC_Vw_s(IW,1)
            IF(DIMN.EQ.3) DES_WALL_VEL(IW,3) = DES_BC_Ww_s(IW,1)
         ENDIF

         WALL_NORMAL(1,1) = -ONE
         WALL_NORMAL(1,2) = ZERO
         IF(DIMN.EQ.3) WALL_NORMAL(1,3) = ZERO

! east wall; in 3D on yz plane (x=ex2->xlength)
      ELSEIF(IW.EQ.2) THEN
         DES_WALL_POS(IW,1) = EX2 + DES_R
         DES_WALL_POS(IW,2) = DES_POS_NEW(L,2)
         IF(DIMN.EQ.3) DES_WALL_POS(IW,3) = DES_POS_NEW(L,3)

         IF(DES_F.NE.ZERO) THEN
         ELSE                 
            DES_WALL_VEL(IW,1) = ZERO
            DES_WALL_VEL(IW,2) = DES_BC_Vw_s(IW,1)
            IF(DIMN.EQ.3) DES_WALL_VEL(IW,3) = DES_BC_Ww_s(IW,1)
         ENDIF

         WALL_NORMAL(2,1) = ONE
         WALL_NORMAL(2,2) = ZERO
         IF(DIMN.EQ.3) WALL_NORMAL(2,3) = ZERO

! bottom wall; in 3D on xz plane (y=by1->0)
      ELSEIF(IW.EQ.3) THEN
         DES_WALL_POS(IW,1) = DES_POS_NEW(L,1)
         DES_WALL_POS(IW,2) = BY1 - DES_R + (A*SINOMEGAT)
         IF(DIMN.EQ.3) DES_WALL_POS(IW,3) = DES_POS_NEW(L,3)

         IF(DES_F.NE.ZERO) THEN
            DES_WALL_VEL(IW,1) = lid_vel
            DES_WALL_VEL(IW,2) =  A*OMEGA_W*COSOMEGAT
         ELSE
            DES_WALL_VEL(IW,1) = DES_BC_Uw_s(IW,1)
            DES_WALL_VEL(IW,2) = ZERO
            IF(DIMN.EQ.3) DES_WALL_VEL(IW,3) = DES_BC_Ww_s(IW,1)
         ENDIF

         WALL_NORMAL(3,1) = ZERO
         WALL_NORMAL(3,2) = -ONE
         IF(DIMN.EQ.3) WALL_NORMAL(3,3) = ZERO

! top wall; in 3D on xz plane (y=ty2->ylength)
      ELSEIF(IW.EQ.4) THEN
         DES_WALL_POS(IW,1) = DES_POS_NEW(L,1)
         DES_WALL_POS(IW,2) = TY2 + DES_R
         IF(DIMN.EQ.3) DES_WALL_POS(IW,3) = DES_POS_NEW(L,3)

         IF(DES_F.NE.ZERO) THEN
            DES_WALL_VEL(IW,1) = -lid_vel
         ELSE
            DES_WALL_VEL(IW,1) = DES_BC_Uw_s(IW,1)
            DES_WALL_VEL(IW,2) = ZERO
            IF(DIMN.EQ.3) DES_WALL_VEL(IW,3) = DES_BC_Ww_s(IW,1)
         ENDIF

         WALL_NORMAL(4,1) = ZERO
         WALL_NORMAL(4,2) = ONE
         IF(DIMN.EQ.3) WALL_NORMAL(4,3) = ZERO

! south wall; in 3D on xy plane (z=sz1->0)
      ELSEIF(IW.EQ.5) THEN
         DES_WALL_POS(IW,1) = DES_POS_NEW(L,1)
         DES_WALL_POS(IW,2) = DES_POS_NEW(L,2)
         DES_WALL_POS(IW,3) = SZ1 - DES_R

         IF(DES_F.NE.ZERO) THEN
         ELSE
            DES_WALL_VEL(IW,1) = DES_BC_Uw_s(IW,1)
            DES_WALL_VEL(IW,2) = DES_BC_Vw_s(IW,1)
            DES_WALL_VEL(IW,3) = ZERO 
         ENDIF

         WALL_NORMAL(5,1) = ZERO
         WALL_NORMAL(5,2) = ZERO
         WALL_NORMAL(5,3) = -ONE

! north wall; in 3D on xy plane (z=nz2->zlength)
      ELSEIF(IW.EQ.6) THEN
         DES_WALL_POS(IW,1) = DES_POS_NEW(L,1)
         DES_WALL_POS(IW,2) = DES_POS_NEW(L,2)
         DES_WALL_POS(IW,3) = NZ2 + DES_R

         IF(DES_F.NE.ZERO) THEN
         ELSE
            DES_WALL_VEL(IW,1) = DES_BC_Uw_s(IW,1)
            DES_WALL_VEL(IW,2) = DES_BC_Vw_s(IW,1)
            DES_WALL_VEL(IW,3) = ZERO 
         ENDIF

         WALL_NORMAL(6,1) = ZERO
         WALL_NORMAL(6,2) = ZERO
         WALL_NORMAL(6,3) = ONE

      ENDIF
      
      RETURN
      END SUBROUTINE CFWALLPOSVEL



