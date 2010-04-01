!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFWALLCONTACT(WALL, L, WALLCONTACTI)                   C
!  Purpose: DES - Checking for contact with walls
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CFWALLCONTACT(WALL, L, WALLCONTACTI)

      Use discretelement
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
      USE des_bc
      IMPLICIT NONE

!-----------------------------------------------
! Local variables
!-----------------------------------------------        
      INTEGER L, I, K, WALL, WALLCONTACTI, WC
      DOUBLE PRECISION A, OMEGA, OOMEGA2, ASINOMEGAT 

!-----------------------------------------------        
     

      A = ZERO
      OMEGA = ZERO
      ASINOMEGAT = ZERO
      IF(DES_F.NE.ZERO) THEN
         OMEGA = 2.0D0*PI*DES_F
         OOMEGA2 = ONE/(OMEGA**2)
         A = DES_GAMMA*GRAV(2)*OOMEGA2
         ASINOMEGAT = A*SIN(OMEGA*S_TIME)
      ENDIF

      WALLCONTACTI = 0

! west wall (X)
      IF(WALL.EQ.1.AND.(.NOT.DES_PERIODIC_WALLS_X)) THEN
         IF((DES_POS_NEW(L,1)-WX1).LE.DES_RADIUS(L)) THEN
            IF(DES_MO_X)THEN
               CALL DES_MASS_OUTLET(L,'XW',WALLCONTACTI)
            ELSE
               WALLCONTACTI = 1
            ENDIF
         ENDIF

! east wall (X)
      ELSE IF(WALL.EQ.2.AND.(.NOT.DES_PERIODIC_WALLS_X)) THEN
         IF((EX2-DES_POS_NEW(L,1)).LE.DES_RADIUS(L)) THEN
            IF(DES_MO_X)THEN
               CALL DES_MASS_OUTLET(L,'XE',WALLCONTACTI)
            ELSE
               WALLCONTACTI = 1
            ENDIF
         ENDIF

! south wall (Y)
      ELSE IF(WALL.EQ.3.AND.(.NOT.DES_PERIODIC_WALLS_Y)) THEN
         IF((DES_POS_NEW(L,2)-(BY1+ASINOMEGAT)).LE.DES_RADIUS(L)) THEN
            IF(DES_MO_Y)THEN
               CALL DES_MASS_OUTLET(L,'YS',WALLCONTACTI)
            ELSE
               WALLCONTACTI = 1
            ENDIF
         ENDIF

! north wall (Y)
      ELSE IF(WALL.EQ.4.AND.(.NOT.DES_PERIODIC_WALLS_Y)) THEN
         IF((TY2-DES_POS_NEW(L,2)).LE.DES_RADIUS(L)) THEN
            IF(DES_MO_Y)THEN
               CALL DES_MASS_OUTLET(L,'YN',WALLCONTACTI)
            ELSE
               WALLCONTACTI = 1
            ENDIF
         ENDIF

! bottom wall (Z)
      ELSE IF(WALL.EQ.5.AND.(.NOT.DES_PERIODIC_WALLS_Z)) THEN
         IF((DES_POS_NEW(L,3)-SZ1).LE.DES_RADIUS(L)) THEN
            IF(DES_MO_Z)THEN
               CALL DES_MASS_OUTLET(L,'ZB',WALLCONTACTI)
            ELSE
               WALLCONTACTI = 1
            ENDIF
         ENDIF

! top wall (Z)
      ELSE IF(WALL.EQ.6.AND.(.NOT.DES_PERIODIC_WALLS_Z)) THEN
         IF((NZ2-DES_POS_NEW(L,3)).LE.DES_RADIUS(L)) THEN
            IF(DES_MO_Z)THEN
               CALL DES_MASS_OUTLET(L,'ZT',WALLCONTACTI)
            ELSE
               WALLCONTACTI = 1
            ENDIF
         ENDIF
      ENDIF

      RETURN
      END SUBROUTINE CFWALLCONTACT


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: DES_MASS_OUTLET(NP, BCC, WC)                           !
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

      SUBROUTINE DES_MASS_OUTLET(NP, BCC, WC)

      USE des_bc
      USE discretelement
      USE param1

      IMPLICIT NONE

!-----------------------------------------------  
! Local variables
!-----------------------------------------------
      INTEGER NP        ! Particle ID Number
      INTEGER WC        ! Wall contact value (0: No Contact, 1: Contact)
      INTEGER BCV_I, BCV

      DOUBLE PRECISION XPOS, YPOS, ZPOS ! particle x,y,z position
      DOUBLE PRECISION RAD ! particle radius
 
      CHARACTER*2 BCC   ! Boundary condition class (Xw,Xe,Ys,Yn,Zb,Zt)

!-----------------------------------------------  

      RAD = DES_RADIUS(NP)

      XPOS = DES_POS_NEW(NP,1) 
      YPOS = DES_POS_NEW(NP,2)
      IF (DIMN .EQ. 2) THEN
         ZPOS = ONE
      ELSE
         ZPOS = DES_POS_NEW(NP,3)
      ENDIF

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

            IF(DIMN == 2)THEN   ! 2D domain

               IF(BCC == 'XW' .OR. BCC == 'XE')THEN
                  IF(DES_BC_Y_s(BCV) < (YPOS - RAD) .AND. &
                     DES_BC_Y_n(BCV) > (YPOS + RAD) )THEN
                     WC = 0
                     PEA(NP,3) = .TRUE.
                  ENDIF
               ENDIF
               IF(BCC == 'YS' .OR. BCC == 'YN' )THEN
                  IF(DES_BC_X_w(BCV) < (XPOS - RAD) .AND. &
                     DES_BC_X_e(BCV) > (XPOS + RAD) )THEN
                     WC = 0
                     PEA(NP,3) = .TRUE.
                  ENDIF
               ENDIF

            ELSE  ! 3D domain

                IF(BCC == 'XW' .OR. BCC == 'XE')THEN
                  IF(DES_BC_Y_s(BCV) < (YPOS - RAD) .AND. &
                     DES_BC_Y_n(BCV) > (YPOS + RAD) .AND. &
                     DES_BC_Z_b(BCV) < (ZPOS - RAD) .AND. &
                     DES_BC_Z_t(BCV) > (ZPOS + RAD) )THEN
                     WC = 0
                     PEA(NP,3) = .TRUE.
                  ENDIF
               ENDIF
               IF(BCC == 'YS' .OR. BCC == 'YN')THEN
                  IF(DES_BC_X_w(BCV) < (XPOS - RAD) .AND. &
                     DES_BC_X_e(BCV) > (XPOS + RAD) .AND. &
                     DES_BC_Z_b(BCV) < (ZPOS - RAD) .AND. &
                     DES_BC_Z_t(BCV) > (ZPOS + RAD) )THEN
                     WC = 0
                     PEA(NP,3) = .TRUE.
                  ENDIF
               ENDIF
               IF(BCC == 'ZB' .OR. BCC == 'ZT')THEN
                  IF(DES_BC_X_w(BCV) < (XPOS - RAD) .AND. &
                     DES_BC_X_e(BCV) > (XPOS + RAD) .AND. &
                     DES_BC_Y_s(BCV) < (YPOS - RAD) .AND. &
                     DES_BC_Y_n(BCV) > (YPOS + RAD) )THEN
                     WC = 0
                     PEA(NP,3) = .TRUE.
                  ENDIF
               ENDIF

            ENDIF   ! endif dimn == 2
         ENDIF   ! endif DES_MO_CLASS(BCV_I) == BCC

      ENDDO  ! loop over BCV_I the no. of outlet boundaries

 9000 FORMAT(F9.4,2X,'PEA(3)=T ',I4,1X,6(F9.4,1X))

      END SUBROUTINE DES_MASS_OUTLET
