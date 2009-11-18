!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: DES_MASS_INLET(BCV_I)                                    !
!                                                                      !
!  Purpose:  This routine fills in the necessary information for new   !
!  particles entereing the system.                                     !
!                                                                      !
!                                                                      !
!  Author: J.Musser                                   Date: 13-Jul-09  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE DES_MASS_INLET(BCV_I)

      USE compar
      USE constant
      USE des_bc
      USE discretelement
      USE funits
      USE geometry
      USE indices
      USE param1
      USE physprop

      IMPLICIT NONE

!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER IJK             ! Necessary for function.inc
      INTEGER I, J, K         ! Loop Counter
      INTEGER IS, IE          ! Loop Counters Start/End
      INTEGER JS, JE          ! Loop Counters Start/End
      INTEGER KS, KE          ! Loop Counters Start/End
      INTEGER IP, NP, LL      ! Counter, Index
      INTEGER LS              ! Loop Counter start value
      INTEGER TMPI, NPG
      INTEGER BCV_I, BCV      ! Boundary Condition ID
! this variable temporarily stores the number ids of all particles 
! in the ijk location of interest
      INTEGER, DIMENSION(:), ALLOCATABLE :: HOLDER
!-----------------------------------------------

      INCLUDE 'function.inc'

      LS = 1
      BCV = DES_BC_MI_ID(BCV_I)

      DO IP = 1, PI_COUNT(BCV_I) !Loop over the particles being injected

! Check to see if MAX_PIS has been exceeded, if so, STOP
         IF(PIS .GE. MAX_PIS) THEN
            WRITE(UNIT_LOG, 1000)
             WRITE(*,1000)
             CALL MFIX_EXIT(myPE)
         ENDIF 

! Find the first free space in the particle existance array
         DO NP = LS, MAX_PIS
            IF(.NOT.PEA(NP,1)) THEN
               LS = NP
               EXIT
            ENDIF
         ENDDO

! Set the flag in the particle existance array
         PEA(NP,1) = .TRUE. 

! Set the flag that the particle is new.  This allows it to be ignored
! by various subroutines
         PEA(NP,2) = .TRUE.
      
! Increment the particle in system value by one
         PIS = PIS + 1

! Set the initial velocity values
         DES_VEL_OLD(NP,1) = DES_BC_U_s(BCV)
         DES_VEL_OLD(NP,2) = DES_BC_V_s(BCV)
         IF(DIMN == 3) DES_VEL_OLD(NP,3) = DES_BC_W_s(BCV)

         DES_VEL_NEW(NP,:) = DES_VEL_OLD(NP,:)

! Set the initial angular velocity values
         OMEGA_OLD(NP,:) = 0
         OMEGA_NEW(NP,:) = 0

! Set the particle radius value
         DES_RADIUS(NP) = (D_P0(1) * HALF)

! Set the particle density value
         RO_Sol(NP) = RO_S(1)

! Set the particle mass phase
         PIJK(NP,5) = 1

! Calculate the new particle's Volume, Mass, OMOI
         PVOL(NP) = PI * DES_RADIUS(NP)**3
         PMASS(NP) = PVOL(NP) * RO_Sol(NP)
         OMOI(NP) = 5.d0 / (2.d0 * PMASS(NP) * DES_RADIUS(NP)**2) 

! Set the initial position values based on mass inlet class
         CALL DES_PLACE_NEW_PARTICLE(NP, BCV_I)
         DES_POS_NEW(NP,:) = DES_POS_OLD(NP,:)

! Determine the I,J,K indices of the cell containing the new particle(s)
! by checking the cells near the mass inlet using GS_ARRAY; note the
! indices will place the particle in a ghost cell
         DO I = GS_ARRAY(BCV_I,1), GS_ARRAY(BCV_I,2)
            IF(DES_POS_NEW(NP,1) .LT. XE(1))THEN 
               PIJK(NP,1) = 1
               EXIT
            ELSEIF(DES_POS_NEW(NP,1) .GT. XLENGTH)THEN 
               PIJK(NP,1) = IMAX2
               EXIT
            ELSEIF((DES_POS_NEW(NP,1) .GE. XE(I-1)) .AND. &
            (DES_POS_NEW(NP,1) .LT. XE(I))) THEN
               PIJK(NP,1) = I
               EXIT
            ENDIF
         ENDDO

         DO J = GS_ARRAY(BCV_I,3), GS_ARRAY(BCV_I,4)
            IF(DES_POS_NEW(NP,2) .LT. YN(1))THEN
               PIJK(NP,2) = 1
               EXIT
            ELSEIF(DES_POS_NEW(NP,2) .GT. YLENGTH)THEN 
               PIJK(NP,2) = JMAX2
               EXIT
            ELSEIF((DES_POS_NEW(NP,2) .GE. YN(J-1)) .AND. &
            (DES_POS_NEW(NP,2) .LT. YN(J))) THEN
               PIJK(NP,2) = J
               EXIT
            ENDIF
         ENDDO

         IF(DIMN.EQ.2) THEN
            PIJK(NP,3)  = 1
         ELSE
            DO K = GS_ARRAY(BCV_I,5), GS_ARRAY(BCV_I,6)
               IF(DES_POS_NEW(NP,3) .LT. ZT(1))THEN
                  PIJK(NP,3) = 1
                  EXIT
               ELSEIF(DES_POS_NEW(NP,3) .GT. ZLENGTH)THEN 
                  PIJK(NP,3) = KMAX2
                  EXIT               
               ELSEIF((DES_POS_NEW(NP,3) .GT. ZT(K-1)) .AND. &
               (DES_POS_NEW(NP,3) .LE. ZT(K))) THEN 
                  PIJK(NP,3) = K
                  EXIT
               ENDIF
            ENDDO
         ENDIF


! update the PIC array for new particles so that 
! 1) any subsequent particles that are to be injected will be checked to 
!    prevent overlap with previously injected particles and
         I = PIJK(NP,1)
         J = PIJK(NP,2)
         K = PIJK(NP,3)
         IJK = FUNIJK(I,J,K)
         PIJK(NP,4) = IJK

         IF (ASSOCIATED(PIC(I,J,K)%P)) THEN
            NPG = SIZE(PIC(I,J,K)%P)
            ALLOCATE( HOLDER (NPG) )

! store the particle no. id of all particles at the ijk location            
            TMPI = 1
            DO LL = 1, NPG
               HOLDER(TMPI) = PIC(I,J,K)%P(LL)
               TMPI = TMPI + 1
            ENDDO
            DEALLOCATE(PIC(I,J,K)%P)

! essentially increasing the no. of particles at the ijk location by 1
            ALLOCATE(PIC(I,J,K)%P(TMPI))
            DO LL = 1, TMPI - 1
               PIC(I,J,K)%P(LL) = HOLDER(LL)
            ENDDO
! storing the new particle no. id in the list
            PIC(I,J,K)%P(TMPI) = NP
            DEALLOCATE(HOLDER)
         ELSE
! no other particles were at this ijk location
            TMPI = 1
            ALLOCATE(PIC(I,J,K)%P(TMPI))
            PIC(I,J,K)%P(TMPI) = NP
         ENDIF 
         PINC(IJK) = TMPI

      ENDDO   ! end loop over the no. of injected particles (IP)

 1000 FORMAT(/1X,70('*')//&
         ' From: DES_MASS_INLET -',/&
         ' Message: Maximum number of particles in the system MAX_PIS',/&
         ' has been exceeded.  Increase the value in mfix.dat',/&
         1X,70('*')/)
 
      RETURN
      END SUBROUTINE DES_MASS_INLET



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: DES_PLACE_NEW_PARTICLE                                 !
!                                                                      !
!  Purpose:  This routine uses the classification information to place !
!  a new particle in the proper location.                              !
!                                                                      !
!  Author: J.Musser                                   Date: 14-Aug-09  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE DES_PLACE_NEW_PARTICLE(NP, BCV_I)

      USE compar
      USE des_bc
      USE discretelement
      USE funits
      USE geometry
      USE param1
      USE physprop

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! the current particle no.
      INTEGER NP   
! the associated bc no.      
      INTEGER BCV_I, BCV  
! a random number between 0-1
      DOUBLE PRECISION RAND1, RAND2, TMPDP
! 'lower' x,y,z location of current window
      DOUBLE PRECISION MI_WIN
! a random index number for I_OF_MI or J_OF_MI
      INTEGER MI_ORD
!      
      LOGICAL TOUCHING 
!-----------------------------------------------

      TOUCHING = .TRUE.
      BCV = DES_BC_MI_ID(BCV_I)

      SELECT CASE(DES_MI_CLASS(BCV_I))   

! 2D domain with verticle mass inlet: (west face)
! ----------------------------------------
         CASE ('XW')
! Offset Y axis by 1 particle diameter into ghost cell
              DES_POS_OLD(NP,1) = -D_P0(1)

            IF(PARTICLE_PLCMNT(BCV_I) == 'RAND')THEN
               DO WHILE (TOUCHING)
! Place randomly between DES_BC_Y_s(BCV) and DES_BC_Y_n(BCV)
                  CALL RANDOM_NUMBER(RAND1)
                  DES_POS_OLD(NP,2) = DES_BC_Y_s(BCV) + RAND1*&
                     (DES_BC_Y_n(BCV) - DES_BC_Y_s(BCV))
! Test that no new particles are touching
                  CALL DES_NEW_PARTICLE_TEST(NP, TOUCHING, BCV_I)
               ENDDO
            ENDIF
            IF(PARTICLE_PLCMNT(BCV_I) == 'ORDR')THEN
! Obtain a random number from (0,1]
               CALL RANDOM_NUMBER(RAND1)
               MI_WIN = dble(MI_ORDER(BCV_I)%VALUE(MI_FACTOR(BCV_I)))*&
                  MI_WINDOW(BCV_I)   
               DES_POS_OLD(NP,2) = RAND1 *&
                  (MI_WINDOW(BCV_I) - D_P0(1)) +&   !play room
                  MI_WIN + D_P0(1)*HALF +&   !shift loc. into window 
                  DES_BC_Y_s(BCV)   !shift relative to BC location
               IF(MI_FACTOR(BCV_I) < SIZE(MI_ORDER(BCV_I)%VALUE))THEN
                  MI_FACTOR(BCV_I) = MI_FACTOR(BCV_I) + 1
               ELSE
                  MI_FACTOR(BCV_I) = 1
               ENDIF
            ENDIF

! 2D domain with verticle mass inlet: (east face)
! ----------------------------------------
         CASE ('XE')
! Offset Y axis by 1 particle diameter into ghost cell
            DES_POS_OLD(NP,1) = XLENGTH + D_P0(1)

            IF(PARTICLE_PLCMNT(BCV_I) == 'RAND')THEN
               DO WHILE (TOUCHING)
! Obtain a random number from (0,1]
                  CALL RANDOM_NUMBER(RAND1)
! Place randomly between DES_BC_Y_s(BCV) and DES_BC_Y_n(BCV)
                  DES_POS_OLD(NP,2) = DES_BC_Y_s(BCV) + RAND1*&
                     (DES_BC_Y_n(BCV)-DES_BC_Y_s(BCV))
! Test that no new particles are touching
                  CALL DES_NEW_PARTICLE_TEST(NP, TOUCHING, BCV_I)
               ENDDO
            ENDIF
            IF(PARTICLE_PLCMNT(BCV_I) == 'ORDR')THEN
! Obtain a random number from (0,1]
               CALL RANDOM_NUMBER(RAND1)
               MI_WIN = dble(MI_ORDER(BCV_I)%VALUE(MI_FACTOR(BCV_I)))*&
                  MI_WINDOW(BCV_I)
               DES_POS_OLD(NP,2) = RAND1 *&
                  (MI_WINDOW(BCV_I) - D_P0(1)) +&   !play room
                  MI_WIN + D_P0(1)*HALF +&   !shift loc. into window
                  DES_BC_Y_s(BCV)   !shift relative to BC location
               IF(MI_FACTOR(BCV_I) < SIZE(MI_ORDER(BCV_I)%VALUE))THEN
                  MI_FACTOR(BCV_I) = MI_FACTOR(BCV_I) + 1
               ELSE
                  MI_FACTOR(BCV_I) = 1
               ENDIF
            ENDIF

! 2D domain with horizontal mass inlet: (south face)
! ----------------------------------------
         CASE ('YS')
! Offset X axis by 1 particle diameter into ghost cell
            DES_POS_OLD(NP,2) = -D_P0(1)

            IF(PARTICLE_PLCMNT(BCV_I) == 'RAND')THEN
               DO WHILE (TOUCHING)
! Obtain a random number from (0,1]
                  CALL RANDOM_NUMBER(RAND1)
! Place randomly between DES_BC_X_w(BCV) and DES_BC_X_e(BCV)
                  DES_POS_OLD(NP,1) = DES_BC_X_w(BCV) + RAND1*&
                     (DES_BC_X_e(BCV) - DES_BC_X_w(BCV))
! Test that no new particles are touching
                  CALL DES_NEW_PARTICLE_TEST(NP, TOUCHING, BCV_I)
                ENDDO
            ENDIF
            IF(PARTICLE_PLCMNT(BCV_I) == 'ORDR')THEN
! Obtain a random number from (0,1]
               CALL RANDOM_NUMBER(RAND1)
               MI_WIN = dble(MI_ORDER(BCV_I)%VALUE(MI_FACTOR(BCV_I)))*&
                  MI_WINDOW(BCV_I)
               DES_POS_OLD(NP,1) = RAND1 *&
                  (MI_WINDOW(BCV_I) - D_P0(1)) +&   !play room
                  MI_WIN + D_P0(1)*HALF +&   !shift loc. into window
                  DES_BC_X_w(BCV)   !shift relative to BC location
               IF(MI_FACTOR(BCV_I) < SIZE(MI_ORDER(BCV_I)%VALUE))THEN
                  MI_FACTOR(BCV_I) = MI_FACTOR(BCV_I) + 1
               ELSE
                  MI_FACTOR(BCV_I) = 1
               ENDIF
            ENDIF

! 2D domain with horizontal mass inlet: (north face)
! ----------------------------------------
         CASE ('YN')
! Offset X axis by 1 particle diameter into ghost cell
            DES_POS_OLD(NP,2) = YLENGTH + D_P0(1)

            IF(PARTICLE_PLCMNT(BCV_I) == 'RAND')THEN
               DO WHILE (TOUCHING)
! Obtain a random number from (0,1]
                  CALL RANDOM_NUMBER(RAND1)
! Place randomly between DES_BC_X_w(BCV) and DES_BC_X_e(BCV)
                  DES_POS_OLD(NP,1) = DES_BC_X_w(BCV) + RAND1*&
                     (DES_BC_X_e(BCV) - DES_BC_X_w(BCV))
! Test that no new particles are touching
                  CALL DES_NEW_PARTICLE_TEST(NP, TOUCHING, BCV_I)
               ENDDO
            ENDIF
            IF(PARTICLE_PLCMNT(BCV_I) == 'ORDR')THEN
! Obtain a random number from (0,1]
               CALL RANDOM_NUMBER(RAND1)
               MI_WIN = dble(MI_ORDER(BCV_I)%VALUE(MI_FACTOR(BCV_I)))*&
                  MI_WINDOW(BCV_I)
               DES_POS_OLD(NP,1) = RAND1 *&
                  (MI_WINDOW(BCV_I) - D_P0(1)) + &   !play room
                  MI_WIN + D_P0(1)*HALF +&   !shift loc. into window
                  DES_BC_X_w(BCV)   !shift relative to BC location
               IF(MI_FACTOR(BCV_I) < SIZE(MI_ORDER(BCV_I)%VALUE))THEN
                  MI_FACTOR(BCV_I) = MI_FACTOR(BCV_I) + 1
               ELSE
                  MI_FACTOR(BCV_I) = 1
               ENDIF
            ENDIF


! 3D domain with mass inlet on XZ plane: (south face)
! ----------------------------------------
         CASE ('XZs')
! Offset XZ plane by 1 particle diameter into ghost cell
            DES_POS_OLD(NP,2) = -D_P0(1)
            IF(PARTICLE_PLCMNT(BCV_I) == 'RAND')THEN
               DO WHILE (TOUCHING)
! Obtain a random number from (0,1]
                  CALL RANDOM_NUMBER(RAND1)
                  CALL RANDOM_NUMBER(RAND2)
! Place randomly between DES_BC_X_w(BCV) and DES_BC_X_e(BCV)
                  DES_POS_OLD(NP,1) = DES_BC_X_w(BCV)+ RAND1*&
                     (DES_BC_X_e(BCV) - DES_BC_X_w(BCV))
! Place randomly between DES_BC_Z_b(BCV) and DES_BC_Z_t(BCV)
                  DES_POS_OLD(NP,3) = DES_BC_Z_b(BCV) + RAND2*&
                     (DES_BC_Z_t(BCV) - DES_BC_Z_b(BCV))
! Test that no new particles are touching
                  CALL DES_NEW_PARTICLE_TEST(NP, TOUCHING, BCV_I)
               ENDDO
            ENDIF
            IF(PARTICLE_PLCMNT(BCV_I) == 'ORDR')THEN
! Obtain a random number from (0,1]
               CALL RANDOM_NUMBER(RAND1)
               CALL RANDOM_NUMBER(RAND2)
! Determine x-coordinate from ordered grid position I
               MI_ORD = MI_ORDER(BCV_I)%VALUE(MI_FACTOR(BCV_I))
               MI_WIN = dble(I_OF_MI(BCV_I)%VALUE(MI_ORD))*&
                  MI_WINDOW(BCV_I)
               DES_POS_OLD(NP,1) = RAND1 *&
                  (MI_WINDOW(BCV_I) - D_P0(1)) +&   !play room
                  MI_WIN + D_P0(1)*HALF +&   !shift loc. into window
                  DES_BC_X_w(BCV)   !shift relative to BC location
! Determine z-coordinate from ordered grid position J
               MI_WIN = dble(J_OF_MI(BCV_I)%VALUE(MI_ORD))*&
                  MI_WINDOW(BCV_I)
               DES_POS_OLD(NP,3) = RAND2 *&
                  (MI_WINDOW(BCV_I) - D_P0(1)) +&   !play room
                  MI_WIN + D_P0(1)*HALF +&   !shift loc. into window
                  DES_BC_Z_b(BCV)   !shift relative to BC location
! Increment MI_FACTOR by one or loop back to one
               IF(MI_FACTOR(BCV_I) < SIZE(MI_ORDER(BCV_I)%VALUE))THEN
                  MI_FACTOR(BCV_I) = MI_FACTOR(BCV_I) + 1
               ELSE
                  MI_FACTOR(BCV_I) = 1
               ENDIF
            ENDIF

! 3D domain with mass inlet on XZ plane: (north face)
! ----------------------------------------
         CASE ('XZn')
! Offset XZ plane by 1 particle diameter into ghost cell
            DES_POS_OLD(NP,2) = YLENGTH + D_P0(1)
            IF(PARTICLE_PLCMNT(BCV_I) == 'RAND')THEN
               DO WHILE (TOUCHING)
! Obtain a random number from (0,1]
                  CALL RANDOM_NUMBER(RAND1)
                  CALL RANDOM_NUMBER(RAND2)
! Place randomly between DES_BC_X_w(BCV) and DES_BC_X_e(BCV)
                  DES_POS_OLD(NP,1) = DES_BC_X_w(BCV) + RAND1*&
                     (DES_BC_X_e(BCV) - DES_BC_X_w(BCV))
! Place randomly between DES_BC_Z_b(BCV) and DES_BC_Z_t(BCV)
                  DES_POS_OLD(NP,3) = DES_BC_Z_b(BCV) + RAND2*&
                     (DES_BC_Z_t(BCV) - DES_BC_Z_b(BCV))
! Test that no new particles are touching
                  CALL DES_NEW_PARTICLE_TEST(NP, TOUCHING, BCV_I)
               ENDDO
            ENDIF
            IF(PARTICLE_PLCMNT(BCV_I) == 'ORDR')THEN
! Obtain a random number from (0,1]
               CALL RANDOM_NUMBER(RAND1)
               CALL RANDOM_NUMBER(RAND2)
               MI_ORD = MI_ORDER(BCV_I)%VALUE(MI_FACTOR(BCV_I))
! Determine x-coordinate from ordered grid position I
               MI_WIN = dble(I_OF_MI(BCV_I)%VALUE(MI_ORD))*&
                  MI_WINDOW(BCV_I)
               DES_POS_OLD(NP,1) = RAND1 *&
                  (MI_WINDOW(BCV_I) - D_P0(1)) +&   !play room
                  MI_WIN + D_P0(1)*HALF +&   !shift loc. into window
                  DES_BC_X_w(BCV)   !shift relative to BC location
! Determine z-coordinate from ordered grid position J
               MI_WIN = dble(J_OF_MI(BCV_I)%VALUE(MI_ORD))*&
                  MI_WINDOW(BCV_I)                  
               DES_POS_OLD(NP,3) = RAND2 *&
                  (MI_WINDOW(BCV_I) - D_P0(1)) +&   !play room
                  MI_WIN + D_P0(1)*HALF +&   !shift loc. into window
                  DES_BC_Z_b(BCV)   !shift relative to BC location
! Increment MI_FACTOR by one or loop back to one
               IF(MI_FACTOR(BCV_I) < SIZE(MI_ORDER(BCV_I)%VALUE))THEN
                  MI_FACTOR(BCV_I) = MI_FACTOR(BCV_I) + 1
               ELSE
                  MI_FACTOR(BCV_I) = 1
               ENDIF
            ENDIF

! 3D domain with mass inlet on XY plane: (bottom face)
! ----------------------------------------
         CASE ('XYb')
! Offset XY plane by 1 particle diameter into ghost cell
            DES_POS_OLD(NP,3) = -D_P0(1)

            IF(PARTICLE_PLCMNT(BCV_I) == 'RAND')THEN
               DO WHILE (TOUCHING)
! Obtain a random number from (0,1]
                  CALL RANDOM_NUMBER(RAND1)
                  CALL RANDOM_NUMBER(RAND2)
! Place randomly between DES_BC_X_w(BCV) and DES_BC_X_e(BCV)
                  DES_POS_OLD(NP,1) = DES_BC_X_w(BCV) + RAND1*&
                     (DES_BC_X_e(BCV) - DES_BC_X_w(BCV))
! Place randomly between DES_BC_Y_s(BCV) and DES_BC_Y_n(BCV)
                  DES_POS_OLD(NP,2) = DES_BC_Y_s(BCV) + RAND2*&
                     (DES_BC_Y_n(BCV) - DES_BC_Y_s(BCV))
! Test that no new particles are touching
                  CALL DES_NEW_PARTICLE_TEST(NP, TOUCHING, BCV_I)
               ENDDO
            ENDIF
            IF(PARTICLE_PLCMNT(BCV_I) == 'ORDR')THEN
! Obtain a random number from (0,1]
               CALL RANDOM_NUMBER(RAND1)
               CALL RANDOM_NUMBER(RAND2)
               MI_ORD = MI_ORDER(BCV_I)%VALUE(MI_FACTOR(BCV_I))
! Determine x-coordinate from ordered grid position I
               MI_WIN = dble(I_OF_MI(BCV_I)%VALUE(MI_ORD))*&
                  MI_WINDOW(BCV_I)
               DES_POS_OLD(NP,1) = RAND1 *&
                  (MI_WINDOW(BCV_I) - D_P0(1)) +&   !play room
                  MI_WIN + D_P0(1)*HALF +&   !shift loc. into window
                  DES_BC_X_w(BCV)   !shift relative to BC location
! Determine y-coordinate from ordered grid position J
               MI_WIN = dble(J_OF_MI(BCV_I)%VALUE(MI_ORD))*&
                  MI_WINDOW(BCV_I)
               DES_POS_OLD(NP,2) = RAND2 *&
                  (MI_WINDOW(BCV_I) - D_P0(1)) +&   !play room
                  MI_WIN + D_P0(1)*HALF +&   !shift loc. into window
                  DES_BC_Y_s(BCV)   !shift relative to BC location
! Increment MI_FACTOR by one or loop back to one
                 IF(MI_FACTOR(BCV_I) < SIZE(MI_ORDER(BCV_I)%VALUE))THEN
                    MI_FACTOR(BCV_I) = MI_FACTOR(BCV_I) + 1
                 ELSE
                    MI_FACTOR(BCV_I) = 1
                 ENDIF
              ENDIF

! 3D domain with mass inlet on XY plane: (top face)
! ----------------------------------------
         CASE ('XYt')
! Offset XY plane by 1 particle diameter into ghost cell
            DES_POS_OLD(NP,3) = ZLENGTH + D_P0(1)

            IF(PARTICLE_PLCMNT(BCV_I) == 'RAND')THEN
               DO WHILE (TOUCHING)
! Obtain a random number from (0,1]
                  CALL RANDOM_NUMBER(RAND1)
                  CALL RANDOM_NUMBER(RAND2)
! Place randomly between DES_BC_X_w(BCV) and DES_BC_X_e(BCV)
                  DES_POS_OLD(NP,1) = DES_BC_X_w(BCV) + RAND1*&
                     (DES_BC_X_e(BCV) - DES_BC_X_w(BCV))
! Place randomly between DES_BC_Y_s(BCV) and DES_BC_Y_n(BCV)
                  DES_POS_OLD(NP,2) = DES_BC_Y_s(BCV) + RAND2*&
                     (DES_BC_Y_n(BCV) - DES_BC_Y_s(BCV))
! Test that no new particles are touching
                  CALL DES_NEW_PARTICLE_TEST(NP, TOUCHING, BCV_I)
               ENDDO
            ENDIF
            IF(PARTICLE_PLCMNT(BCV_I) == 'ORDR')THEN
! Obtain a random number from (0,1]
               CALL RANDOM_NUMBER(RAND1)
               CALL RANDOM_NUMBER(RAND2)
               MI_ORD = MI_ORDER(BCV_I)%VALUE(MI_FACTOR(BCV_I))
! Determine x-coordinate from ordered grid position I
               MI_WIN = dble(I_OF_MI(BCV_I)%VALUE(MI_ORD))*&
                  MI_WINDOW(BCV_I)
               DES_POS_OLD(NP,1) = RAND1 *&
                  (MI_WINDOW(BCV_I) - D_P0(1)) +&   !play room
                  MI_WIN + D_P0(1)*HALF +&   !shift loc. into window
                  DES_BC_X_w(BCV)   !shift relative to BC location
! Determine y-coordinate from ordered grid position J
               MI_WIN = dble(J_OF_MI(BCV_I)%VALUE(MI_ORD))*&
                  MI_WINDOW(BCV_I)
               DES_POS_OLD(NP,2) = RAND2 *&
                  (MI_WINDOW(BCV_I) - D_P0(1)) +&   !play room
                  MI_WIN + D_P0(1)*HALF +&   !shift loc. into window
                  DES_BC_Y_s(BCV)   !shift relative to BC location
! Increment MI_FACTOR by one or loop back to one
               IF(MI_FACTOR(BCV_I) < SIZE(MI_ORDER(BCV_I)%VALUE))THEN
                  MI_FACTOR(BCV_I) = MI_FACTOR(BCV_I) + 1
               ELSE
                  MI_FACTOR(BCV_I) = 1
               ENDIF
            ENDIF

! 3D domain with mass inlet on YZ plane: (west face)
! ----------------------------------------
         CASE ('YZw')
! Offset YZ plane by 1 particle diameter into ghost cell
            DES_POS_OLD(NP,1) = -D_P0(1)

            IF(PARTICLE_PLCMNT(BCV_I) == 'RAND')THEN
               DO WHILE (TOUCHING)
! Obtain a random number from (0,1]
                  CALL RANDOM_NUMBER(RAND1)
                  CALL RANDOM_NUMBER(RAND2)
! Place randomly between DES_BC_Y_s(BCV) and DES_BC_Y_n(BCV)
                  DES_POS_OLD(NP,2) = DES_BC_Y_s(BCV) + RAND1*&
                     (DES_BC_Y_n(BCV) - DES_BC_Y_s(BCV))
! Place randomly between DES_BC_Z_b(BCV) and DES_BC_Z_t(BCV)
                  DES_POS_OLD(NP,3) = DES_BC_Z_b(BCV) + RAND2*&
                     (DES_BC_Z_t(BCV) - DES_BC_Z_b(BCV))
! Test that no new particles are touching
                  CALL DES_NEW_PARTICLE_TEST(NP, TOUCHING, BCV_I)
               ENDDO
            ENDIF
            IF(PARTICLE_PLCMNT(BCV_I) == 'ORDR')THEN
! Obtain a random number from (0,1]
               CALL RANDOM_NUMBER(RAND1)
               CALL RANDOM_NUMBER(RAND2)
               MI_ORD = MI_ORDER(BCV_I)%VALUE(MI_FACTOR(BCV_I))
! Determine y-coordinate from ordered grid position I
               MI_WIN = dble(I_OF_MI(BCV_I)%VALUE(MI_ORD))*&
                  MI_WINDOW(BCV_I)
               DES_POS_OLD(NP,2) = RAND1 *&
                  (MI_WINDOW(BCV_I) - D_P0(1)) +&   !play room
                  MI_WIN + D_P0(1)*HALF +&   !shift loc. into window
                  DES_BC_Y_s(BCV)   !shift relative to BC location
! Determine z-coordinate from ordered grid position J
               MI_WIN = dble(J_OF_MI(BCV_I)%VALUE(MI_ORD))*&
                  MI_WINDOW(BCV_I)
               DES_POS_OLD(NP,3) = RAND2 *&
                  (MI_WINDOW(BCV_I) - D_P0(1)) + &   !play room
                  MI_WIN + D_P0(1)*HALF +&   !shift loc. into window
                  DES_BC_Z_b(BCV)   !shift relative to BC location
! Increment MI_FACTOR by one or loop back to one
               IF(MI_FACTOR(BCV_I) < SIZE(MI_ORDER(BCV_I)%VALUE))THEN
                  MI_FACTOR(BCV_I) = MI_FACTOR(BCV_I) + 1
               ELSE
                  MI_FACTOR(BCV_I) = 1
               ENDIF
            ENDIF

! 3D domain with mass inlet on YZ plane: (east face)
! ----------------------------------------
         CASE ('YZe')
! Offset YZ plane by 1 particle diamter into ghost cell
            DES_POS_OLD(NP,1) = XLENGTH + D_P0(1)

            IF(PARTICLE_PLCMNT(BCV_I) == 'RAND')THEN
               DO WHILE (TOUCHING)
!  Obtain a random number from (0,1]
                  CALL RANDOM_NUMBER(RAND1)
                  CALL RANDOM_NUMBER(RAND2)
! Place randomly between DES_BC_Y_s(BCV) and DES_BC_Y_n(BCV)
                  DES_POS_OLD(NP,2) = DES_BC_Y_s(BCV) + RAND1*&
                     (DES_BC_Y_n(BCV) - DES_BC_Y_s(BCV))
! Place randomly between DES_BC_Z_b(BCV) and DES_BC_Z_t(BCV)
                  DES_POS_OLD(NP,3) = DES_BC_Z_b(BCV) + RAND2*&
                     (DES_BC_Z_t(BCV) - DES_BC_Z_b(BCV)) 
! Test that no new particles are touching
                  CALL DES_NEW_PARTICLE_TEST(NP, TOUCHING, BCV_I)
               ENDDO
            ENDIF
            IF(PARTICLE_PLCMNT(BCV_I) == 'ORDR')THEN
! Obtain a random number from (0,1]
               CALL RANDOM_NUMBER(RAND1)
               CALL RANDOM_NUMBER(RAND2)
               MI_ORD = MI_ORDER(BCV_I)%VALUE(MI_FACTOR(BCV_I))
! Determine y-coordinate from ordered grid position I
               MI_WIN = dble(I_OF_MI(BCV_I)%VALUE(MI_ORD))*&
                  MI_WINDOW(BCV_I)
               DES_POS_OLD(NP,2) = RAND1 *&
                  (MI_WINDOW(BCV_I) - D_P0(1)) +&   !play room
                  MI_WIN + D_P0(1)*HALF +&   !shift loc. into window
                  DES_BC_Y_s(BCV)   !shift relative to BC location
! Determine z-coordinate from ordered grid position J
               MI_WIN = dble(J_OF_MI(BCV_I)%VALUE(MI_ORD))*&
                  MI_WINDOW(BCV_I)
               DES_POS_OLD(NP,3) = RAND2 *&
                  (MI_WINDOW(BCV_I) - D_P0(1)) +&   !play room
                  MI_WIN + D_P0(1)*HALF +&   !shift loc. into window
                  DES_BC_Z_b(BCV)   !shift relative to BC location
! Increment MI_FACTOR by one or loop back to one
                 IF(MI_FACTOR(BCV_I) < SIZE(MI_ORDER(BCV_I)%VALUE))THEN
                    MI_FACTOR(BCV_I) = MI_FACTOR(BCV_I) + 1
                 ELSE
                    MI_FACTOR(BCV_I) = 1
                 ENDIF
              ENDIF

         CASE DEFAULT
            PRINT*,'INVALID DES MASS INLET CLASSIFICATION'
            CALL MFIX_EXIT(myPE)
      END SELECT

      RETURN
      END SUBROUTINE DES_PLACE_NEW_PARTICLE


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name:  DES_NEW_PARTICLE_TEST                                 !
!                                                                      !
!  Purpose:  This routine checks if a new particle placed using the    !
!  random inlet was placed in contact with an existing particle.  If   !
!  so a flag is set indicating contact, and the new particle is        !
!  repositioned within the inlet domain.                               !
!                                                                      !
!  Author: J.Musser                                   Date: 14-Aug-09  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE DES_NEW_PARTICLE_TEST(NP, TOUCHING, BCV_I)

      USE compar
      USE constant
      USE des_bc
      USE discretelement
      USE funits
      USE geometry
      USE indices
      USE param1
      USE physprop

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! particle number id of interest
      INTEGER NP
! particle number id of a potential overlapping/contacting particle
      INTEGER NP2      
! index of boundary condition 
      INTEGER BCV_I
! total number of particles in current ijk cell and loop counter
      INTEGER NPG, LL
! i, j, k indices along boundary used for loop counters
      INTEGER I, J, K, IJK

      LOGICAL TOUCHING

      DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT 
      DOUBLE PRECISION  DISTVEC(DIMN), DIST, R_LM
!-----------------------------------------------

      INCLUDE 'function.inc'

      TOUCHING = .FALSE.

      DO K = GS_ARRAY(BCV_I,5), GS_ARRAY(BCV_I,6)
         DO J = GS_ARRAY(BCV_I,3), GS_ARRAY(BCV_I,4)
           DO I =  GS_ARRAY(BCV_I,1), GS_ARRAY(BCV_I,2)
             IJK = FUNIJK(I,J,K)
             IF(ASSOCIATED(PIC(I,J,K)%P)) THEN
               NPG =  SIZE(PIC(I,J,K)%P)                     
               DO LL = 1, NPG
                  NP2 = PIC(I,J,K)%P(LL)
                  DISTVEC(:) = DES_POS_OLD(NP,:) - DES_POS_NEW(NP2,:)
                  DIST = SQRT(DES_DOTPRDCT(DISTVEC,DISTVEC))
                  R_LM = DES_RADIUS(NP) + DES_RADIUS(NP2)
                  IF(DIST .LE. R_LM) TOUCHING = .TRUE.
               ENDDO
             ENDIF
           ENDDO
         ENDDO
       ENDDO

      RETURN
      END SUBROUTINE DES_NEW_PARTICLE_TEST
