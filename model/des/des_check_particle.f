!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: DES_CHECK_PARTICLE                                     !
!                                                                      !
!  Purpose:  This routine is used to check if a new particle has fully 
!  entered the domain.  If so, the flag classifying the particle as new
!  is removed, allowing the particle to respond to contact forces from 
!  walls and other particles.                                          
!                                                                      !
!                                                                      !
!  Author: J.Musser                                   Date: 13-Jul-09  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE DES_CHECK_PARTICLE

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
      INTEGER I,J,K,IJK       ! Unused indices at this time
      INTEGER PC              ! Loop counter
      INTEGER NP              ! Particle Index
! Logical for local debug warnings
      LOGICAL DES_LOC_DEBUG
!-----------------------------------------------
      INCLUDE 'function.inc'


      DES_LOC_DEBUG = .FALSE.

      PC = 1
      DO NP = 1, MAX_PIS
         IF(PC .GT. PIS) EXIT
         IF(.NOT.PEA(NP,1)) CYCLE

! check the status of a new injected particle. if the new particle
! fully entered the domain change its 'new' status PEA(NP,2) to false
         IF(PEA(NP,2))THEN
            IF(DIMN == 2 .AND. &
               (DES_POS_NEW(NP,1) - DES_RADIUS(NP)) .GT. 0 .AND.&
               (DES_POS_NEW(NP,1) + DES_RADIUS(NP)) .LT. XLENGTH.AND.&
               (DES_POS_NEW(NP,2) - DES_RADIUS(NP)) .GT. 0 .AND.&
               (DES_POS_NEW(NP,2) + DES_RADIUS(NP)) .LT. YLENGTH &
            .OR. &
               DIMN == 3 .AND. &
               (DES_POS_NEW(NP,1) - DES_RADIUS(NP)) .GT. 0 .AND.&
               (DES_POS_NEW(NP,1) + DES_RADIUS(NP)) .LT. XLENGTH.AND.&
               (DES_POS_NEW(NP,2) - DES_RADIUS(NP)) .GT. 0 .AND.&
               (DES_POS_NEW(NP,2) + DES_RADIUS(NP)) .LT. YLENGTH .AND. &
               (DES_POS_NEW(NP,3) - DES_RADIUS(NP)) .GT. 0 .AND.&
               (DES_POS_NEW(NP,3) + DES_RADIUS(NP)) .LT. ZLENGTH)THEN

               PEA(NP,2) = .FALSE.
               PC = PC + 1

            ENDIF  ! end if particle has fully entered the system

! check the status of an exiting particle (a particle that has come into
! contact with a mass outlet).
         ELSEIF(PEA(NP,3))THEN
            IF(DIMN == 2 .AND. &
               (DES_POS_NEW(NP,1) - DES_RADIUS(NP)) .GT. 0 .AND.&
               (DES_POS_NEW(NP,1) + DES_RADIUS(NP)) .LT. XLENGTH.AND.&
               (DES_POS_NEW(NP,2) - DES_RADIUS(NP)) .GT. 0 .AND.&
               (DES_POS_NEW(NP,2) + DES_RADIUS(NP)) .LT. YLENGTH &
            .OR. &
               DIMN == 3 .AND. &
               (DES_POS_NEW(NP,1) - DES_RADIUS(NP)) .GT. 0 .AND.&
               (DES_POS_NEW(NP,1) + DES_RADIUS(NP)) .LT. XLENGTH.AND.&
               (DES_POS_NEW(NP,2) - DES_RADIUS(NP)) .GT. 0 .AND.&
               (DES_POS_NEW(NP,2) + DES_RADIUS(NP)) .LT. YLENGTH .AND. &
               (DES_POS_NEW(NP,3) - DES_RADIUS(NP)) .GT. 0 .AND.&
               (DES_POS_NEW(NP,3) + DES_RADIUS(NP)) .LT. ZLENGTH)THEN
! if an exiting particle has altered course so that it is moving back
! into the domain change (no longer on track to exit the domain) and 
! it has fully entered the domain then change its 'exiting' status 
! PEA(NP,3) to false

               PEA(NP,3) = .FALSE.
               PC = PC + 1

            ELSEIF(DIMN == 2 .AND. &
               (DES_POS_NEW(NP,1) + DES_RADIUS(NP)) .LE. 0 .OR.&
               (DES_POS_NEW(NP,1) - DES_RADIUS(NP)) .GE. XLENGTH .OR.&
               (DES_POS_NEW(NP,2) + DES_RADIUS(NP)) .LE. 0 .OR.&
               (DES_POS_NEW(NP,2) - DES_RADIUS(NP)) .GT. YLENGTH &
            .OR. &
               DIMN == 3 .AND. &
               (DES_POS_NEW(NP,1) + DES_RADIUS(NP)) .LE. 0 .OR.&
               (DES_POS_NEW(NP,1) - DES_RADIUS(NP)) .GE. XLENGTH.OR.&
               (DES_POS_NEW(NP,2) + DES_RADIUS(NP)) .LE. 0 .OR.&
               (DES_POS_NEW(NP,2) - DES_RADIUS(NP)) .GE. YLENGTH .OR. &
               (DES_POS_NEW(NP,3) + DES_RADIUS(NP)) .LE. 0 .OR.&
               (DES_POS_NEW(NP,3) - DES_RADIUS(NP)) .GE. ZLENGTH)THEN
! if an exiting particle has fully exited the domain, effectively
! remove the particle (reset PEA(NP,:) to false)       

               PEA(NP,:) = .FALSE.
              
               DES_POS_OLD(NP,:) = ZERO
               DES_POS_NEW(NP,:) = ZERO
               DES_VEL_OLD(NP,:) = ZERO
               DES_VEL_NEW(NP,:) = ZERO
               OMEGA_OLD(NP,:) = ZERO
               OMEGA_NEW(NP,:) = ZERO
               DES_RADIUS(NP) = ZERO
               PMASS(NP) = ZERO
               PVOL(NP) = ZERO
               RO_Sol(NP) = ZERO
               OMOI(NP) = ZERO

               FC(NP,:) = ZERO
               FN(NP,:) = ZERO
               FT(NP,:) = ZERO
               TOW(NP,:) = ZERO
               PN(NP,:) = -1
               PN(NP,1) = 0
               PV(NP,:) = 1
               PFN(NP,:,:) = ZERO
               PFT(NP,:,:) = ZERO

! Note that if particle NP has any neighbors then the particle NP will
! still exist in the neighbor's neighbours list, however, this will
! eventually be cleared in the second call to calc_force_des following
! the removal of NP
               NEIGHBOURS(NP,:) = -1
               NEIGHBOURS(NP,1) = 0 
               
               PPOS(NP,:) = ZERO

               PIS = PIS - 1
! in this case do not increment pc

            ENDIF   ! endif particle has reentered/exited of the domain

         ELSE
            PC = PC + 1 
         ENDIF   ! endif PEA(NP,3) or PEA(NP,2)

      ENDDO   ! loop over NP=1,MAX_PIS

      IF (DES_LOC_DEBUG) WRITE(*,1001)

 1000 FORMAT(3X,'---------- START DES_CHECK_PARTICLE ---------->')
 1001 FORMAT(3X,'<---------- END DES_CHECK_PARTICLE ----------') 

 9000 FORMAT(F9.4,2X,'NOT NEW  ',I4,1X,6(F9.4,1X))
 9001 FORMAT(F9.4,2X,'RE-ENTER ',I4,1X,6(F9.4,1X))
 9002 FORMAT(F9.4,2X,'EXITED   ',I4)

      RETURN
      END SUBROUTINE DES_CHECK_PARTICLE
