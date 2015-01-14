!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: MASS_OUTFLOW_DEM                                        !
!  Author: J.Musser                                   Date: 13-Jul-09  !
!                                                                      !
!  Purpose:  This routine fills in the necessary information for new   !
!  particles entereing the system.                                     !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE MASS_OUTFLOW_DEM(FORCE_NSEARCH)

      use discretelement
      use des_bc
      use bc

      use mpi_utility, only: GLOBAL_ALL_OR

      implicit none

      LOGICAL, INTENT(INOUT) :: FORCE_NSEARCH

      INTEGER :: IJK
      INTEGER :: LC, LP, NP
      INTEGER :: BCV, BCV_I

      DOUBLE PRECISION :: DIST

      DOUBLE PRECISION :: FREEZE(3)

      DO BCV_I = 1, DEM_BCMO

         BCV = DEM_BCMO_MAP(BCV_I)

         SELECT CASE (BC_PLANE(BCV))
         CASE('E','W'); FREEZE = (/1.0d0, 0.0d0, 0.0d0/)
         CASE('N','S'); FREEZE = (/0.0d0, 1.0d0, 0.0d0/)
         CASE('T','B'); FREEZE = (/0.0d0, 0.0d0, 1.0d0/)
         END SELECT

         DO LC=DEM_BCMO_IJKSTART(BCV_I), DEM_BCMO_IJKEND(BCV_I)
            IJK = DEM_BCMO_IJK(LC)
            DO LP= 1,PINC(IJK)

               NP = PIC(IJK)%p(LP)

               IF(PEA(NP,4)) CYCLE

               SELECT CASE (BC_PLANE(BCV))
               CASE('S'); DIST = YN(BC_J_s(BCV)-1) - DES_POS_NEW(2,NP)
               CASE('N'); DIST = DES_POS_NEW(2,NP) - YN(BC_J_s(BCV))
               CASE('W'); DIST = XE(BC_I_w(BCV)-1) - DES_POS_NEW(1,NP)
               CASE('E'); DIST = DES_POS_NEW(1,NP) - XE(BC_I_w(BCV))
               CASE('B'); DIST = ZT(BC_K_b(BCV)-1) - DES_POS_NEW(3,NP)
               CASE('T'); DIST = DES_POS_NEW(3,NP) - ZT(BC_K_b(BCV))
               END SELECT

! The particle is still inside the domain
               IF(DIST > DES_RADIUS(NP)) THEN
                  PEA(NP,3) = .FALSE.

! Check if the particle is crossing over the outlet plane. The velocity
! is 'frozen' normal to the outlet and the particle exiting flags are
! set. This implementation is strict as complex BCs (via STLs) can let
! particles pop through the wall along the outlet.
               ELSEIF(DIST > ZERO) THEN
                  DES_VEL_NEW(:,NP) = DES_VEL_NEW(:,NP)*FREEZE(:)
                  PEA(NP,2) = .TRUE.
                  PEA(NP,3) = .TRUE.

! Ladies and gentlemen, the particle has left the building.
               ELSE
                  CALL DELETE_PARTICLE(NP)
                  FORCE_NSEARCH = .TRUE.
               ENDIF

            ENDDO
         ENDDO
      ENDDO

! Sync the search flag across all processes.
      CALL GLOBAL_ALL_OR(FORCE_NSEARCH)

      RETURN
      END SUBROUTINE MASS_OUTFLOW_DEM




!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: DELETE_PARTICLE                                         !
!  Author: J.Musser                                   Date: 13-Jul-09  !
!                                                                      !
!  Purpose:  This routine is used to check if a new particle has fully !
!  entered the domain.  If so, the flag classifying the particle as new!
!  is removed, allowing the particle to respond to contact forces from !
!  walls and other particles.                                          !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DELETE_PARTICLE(NP)

      USE compar
      USE constant
      USE des_bc
      USE discretelement
      USE funits
      USE geometry
      USE indices
      USE param1
      USE physprop
      USE functions

      IMPLICIT NONE


      INTEGER, INTENT(IN) :: NP

!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER I,J,K,IJK       ! Unused indices at this time
      INTEGER PC              ! Loop counter
      INTEGER NPT             ! Temp Particle Index
! Tmp holder of particle position
      DOUBLE PRECISION XPOS, YPOS, ZPOS

! Loop indices used for clearing forces associated with exiting particles
      INTEGER NLIMNP, NLIM, NEIGHNP, NLIM_NEIGHNP

! Logical for local debug warnings
      LOGICAL DES_LOC_DEBUG
!-----------------------------------------------

      PEA(NP,:) = .FALSE.
      DES_POS_NEW(:,NP) = ZERO
      DES_VEL_NEW(:,NP) = ZERO
      OMEGA_NEW(:,NP) = ZERO

      IF (DO_OLD) THEN
         DES_POS_OLD(:,NP) = ZERO
         DES_VEL_OLD(:,NP) = ZERO
         OMEGA_OLD(:,NP) = ZERO
      ENDIF

      DES_RADIUS(NP) = ZERO
      PMASS(NP) = ZERO
      PVOL(NP) = ZERO
      RO_Sol(NP) = ZERO
      OMOI(NP) = ZERO

      FC(:,NP) = ZERO
      TOW(:,NP) = ZERO

      PPOS(:,NP) = ZERO

      PIP = PIP - 1

      RETURN
      END SUBROUTINE DELETE_PARTICLE
