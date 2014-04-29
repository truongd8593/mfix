!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!      
!  Subroutine: CFWALLCONTACT(WALL, L, WALLCONTACT) 	 
!  Purpose: Check if particle L is in contact with WALL.  If so, set 	 
!           WALLCONTACT to 1, else WALLCONTACT is 0
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

      SUBROUTINE CFWALLCONTACT(WALL, L, WALLCONTACT)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param1
      USE constant
      USE parallel      
      USE compar
      Use discretelement      
      USE des_bc
      use geometry, only: DO_K
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments      
!-----------------------------------------------
! Given wall ID number (1=west, 2=east, 3=south, 4=north, 5=bottom,
! 6=top)
      INTEGER, INTENT (IN) :: WALL
! Given particle ID number      
      INTEGER, INTENT (IN) :: L
! Flag to indicate whether given particle is in contact with given wall
! (1=contact, 0 = no contact)
      INTEGER, INTENT (INOUT) :: WALLCONTACT
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! local variables for x, y, z position of the particle
      DOUBLE PRECISION :: XPOS, YPOS, ZPOS      
! local variables to define system dimensions
      DOUBLE PRECISION :: LXE, LXW, LYN, LYS, LZT, LZB
! local variables: distance between particle surface and wall
      DOUBLE PRECISION :: DistApart
!-----------------------------------------------

! assign temporary local variables for quick reference
      LXE = EX2
      LXW = WX1
      LYN = TY2
      LYS = BY1
      LZT = NZ2
      LZB = SZ1

! assign temporary local variables for manipulation/use
      XPOS = DES_POS_NEW(L,1)
      YPOS = DES_POS_NEW(L,2)
      IF(DO_K) ZPOS = DES_POS_NEW(L,3)
 

! initialize
      WALLCONTACT = 0

      IF (DES_LE_BC) THEN
! Current implementation of Lees & Edwards boundaries implies all other
! boundaries are periodic (i.e. no walls in system)
         RETURN
      ELSEIF (DES_PERIODIC_WALLS) THEN
! Check if current wall corresponds to a periodic boundary (i.e. no wall) 
         IF( (DES_PERIODIC_WALLS_X .AND. (WALL.EQ.1.OR.WALL.EQ.2)).OR.&
             (DES_PERIODIC_WALLS_Y .AND. (WALL.EQ.3.OR.WALL.EQ.4)).OR.&
             (DO_K.AND.DES_PERIODIC_WALLS_Z .AND. &
             (WALL.EQ.5.OR.WALL.EQ.6)) ) THEN
            RETURN
         ENDIF         
      ENDIF

! Note that if no cohesion is used WALL_VDW_OUTER_CUTOFF = zero.
! west wall (X)
      IF(WALL.EQ.1) THEN
         DistApart = XPOS-LXW-DES_RADIUS(L)
! consider wall contact calculations for particle wall distances less 
! than the cutoff
         IF( DistApart <= WALL_VDW_OUTER_CUTOFF ) THEN
            IF(DES_MO_X)THEN
               CALL DES_MASS_OUTLET(L,'XW',WALLCONTACT,DistApart)
            ELSE
               WALLCONTACT = 1
            ENDIF
         ENDIF

! east wall (X)
      ELSEIF(WALL.EQ.2) THEN
         DistApart = LXE-XPOS-DES_RADIUS(L)
         IF( DistApart <= WALL_VDW_OUTER_CUTOFF ) THEN
            IF(DES_MO_X)THEN
               CALL DES_MASS_OUTLET(L,'XE',WALLCONTACT,DistApart)
            ELSE
               WALLCONTACT = 1
            ENDIF
         ENDIF

! south wall (Y)
      ELSEIF(WALL.EQ.3) THEN
         DistApart = YPOS-(LYS)-DES_RADIUS(L)
         IF( DistApart <= WALL_VDW_OUTER_CUTOFF ) THEN
            IF(DES_MO_Y)THEN
               CALL DES_MASS_OUTLET(L,'YS',WALLCONTACT,DistApart)
            ELSE
               WALLCONTACT = 1
            ENDIF
         ENDIF

! north wall (Y)
      ELSEIF(WALL.EQ.4) THEN
         DistApart = LYN-YPOS-DES_RADIUS(L)
         IF( DistApart <= WALL_VDW_OUTER_CUTOFF ) THEN
            IF(DES_MO_Y)THEN
               CALL DES_MASS_OUTLET(L,'YN',WALLCONTACT,DistApart)
            ELSE
               WALLCONTACT = 1
            ENDIF
         ENDIF

! bottom wall (Z)
      ELSEIF(WALL.EQ.5) THEN
         DistApart = ZPOS-LZB-DES_RADIUS(L)
         IF( DistApart <= WALL_VDW_OUTER_CUTOFF ) THEN
            IF(DES_MO_Z)THEN
               CALL DES_MASS_OUTLET(L,'ZB',WALLCONTACT,DistApart)
            ELSE
               WALLCONTACT = 1
            ENDIF
         ENDIF

! top wall (Z)
      ELSEIF(WALL.EQ.6) THEN
         DistApart = LZT-ZPOS-DES_RADIUS(L)
         IF( DistApart <= WALL_VDW_OUTER_CUTOFF ) THEN
            IF(DES_MO_Z)THEN
               CALL DES_MASS_OUTLET(L,'ZT',WALLCONTACT,DistApart)
            ELSE
               WALLCONTACT = 1
            ENDIF
         ENDIF
      ENDIF

      RETURN
      END SUBROUTINE CFWALLCONTACT


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: DES_MASS_OUTLET(NP, BCC, WC, dist)                      !
!                                                                      !
!  Purpose:                                                            !
!  Check if the particle's position overlaps the outlet position 
!  by more than its radius (i.e., the particle is more than just
!  touching the boundary)
!                                                                      !
!  Author: J.Musser                                   Date:  5-Oct-09  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE DES_MASS_OUTLET(NP, BCC, WC, dist)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE des_bc
      USE discretelement
      USE param1
      use geometry, only: NO_K

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Particle ID Number
      INTEGER, INTENT (IN) :: NP        
! Wall contact value (0: No Contact, 1: Contact)      
      INTEGER, INTENT (INOUT) :: WC
! distance between particle surface and walls      
      DOUBLE PRECISION, INTENT (IN) :: dist
! Boundary condition class (Xw,Xe,Ys,Yn,Zb,Zt)      
      CHARACTER*2, INTENT (IN) :: BCC   
!-----------------------------------------------  
! Local variables
!-----------------------------------------------
      INTEGER BCV_I, BCV
      DOUBLE PRECISION XPOS, YPOS, ZPOS ! particle x,y,z position
      DOUBLE PRECISION RAD ! particle radius

!-----------------------------------------------  

      RAD = DES_RADIUS(NP)

      XPOS = DES_POS_NEW(NP,1) 
      YPOS = DES_POS_NEW(NP,2)
      ZPOS = merge(ONE, DES_POS_NEW(NP,3), NO_K)

! Set the flag identifying a wall contact to 1 (contact exists).  If the
! particle is in contact with a region of the wall that is listed as an
! outlet for DES particles, then the wall contact is ignored by setting 
! the flag back to zero.
      WC = 1
      
! Loop over all outflow boundary conditions
      DO BCV_I = 1, DES_BCMO
         BCV = DES_BC_MO_ID(BCV_I)
! Find which outlet boundary(s) index coincide with the passed boundary
! condition class.
         IF(DES_MO_CLASS(BCV_I) == BCC)THEN 

            IF(NO_K)THEN   ! 2D domain

               IF(BCC == 'XW' .OR. BCC == 'XE')THEN
                  IF(DES_BC_Y_s(BCV) < (YPOS - RAD) .AND. &
                     DES_BC_Y_n(BCV) > (YPOS + RAD) )THEN
                     WC = 0
! When using Van der Waals cohesion model this subroutine may be invoked for
! particle wall distances greater than zero.  In such an event, the particle
! should not be considered in contact with the wall/outlet.
                     IF(dist <= ZERO) PEA(NP,3) = .TRUE.
                  ENDIF
               ENDIF
               IF(BCC == 'YS' .OR. BCC == 'YN' )THEN
                  IF(DES_BC_X_w(BCV) < (XPOS - RAD) .AND. &
                     DES_BC_X_e(BCV) > (XPOS + RAD) )THEN
                     WC = 0
                     IF(dist <= ZERO) PEA(NP,3) = .TRUE.
                  ENDIF
               ENDIF

            ELSE  ! 3D domain

                IF(BCC == 'XW' .OR. BCC == 'XE')THEN
                  IF(DES_BC_Y_s(BCV) < (YPOS - RAD) .AND. &
                     DES_BC_Y_n(BCV) > (YPOS + RAD) .AND. &
                     DES_BC_Z_b(BCV) < (ZPOS - RAD) .AND. &
                     DES_BC_Z_t(BCV) > (ZPOS + RAD) )THEN
                     WC = 0
                     IF(dist <= ZERO) PEA(NP,3) = .TRUE.
                  ENDIF
               ENDIF
               IF(BCC == 'YS' .OR. BCC == 'YN')THEN
                  IF(DES_BC_X_w(BCV) < (XPOS - RAD) .AND. &
                     DES_BC_X_e(BCV) > (XPOS + RAD) .AND. &
                     DES_BC_Z_b(BCV) < (ZPOS - RAD) .AND. &
                     DES_BC_Z_t(BCV) > (ZPOS + RAD) )THEN
                     WC = 0
                     IF(dist <= ZERO) PEA(NP,3) = .TRUE.
                  ENDIF
               ENDIF
               IF(BCC == 'ZB' .OR. BCC == 'ZT')THEN
                  IF(DES_BC_X_w(BCV) < (XPOS - RAD) .AND. &
                     DES_BC_X_e(BCV) > (XPOS + RAD) .AND. &
                     DES_BC_Y_s(BCV) < (YPOS - RAD) .AND. &
                     DES_BC_Y_n(BCV) > (YPOS + RAD) )THEN
                     WC = 0
                     IF(dist <= ZERO) PEA(NP,3) = .TRUE.
                  ENDIF
               ENDIF

            ENDIF   ! endif NO_K
         ENDIF   ! endif DES_MO_CLASS(BCV_I) == BCC

      ENDDO  ! loop over BCV_I the no. of outlet boundaries

 9000 FORMAT(F9.4,2X,'PEA(3)=T ',I4,1X,6(F9.4,1X))

      END SUBROUTINE DES_MASS_OUTLET
