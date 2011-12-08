!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: THERMO_NBR                                             !
!                                                                      !
!  Purpose: Check the data provided for the des                        !
!                                                                      !
!                                                                      !
!  Author: J.Musser                                   Date: 17-Feb-11  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE THERMO_NBR(L,LL,DIST)

      Use des_thermo
      Use discretelement

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: L, LL

      DOUBLE PRECISION, INTENT(IN) :: DIST

      INTEGER NEIGH_L, NEIGH_LL

      LOGICAL ALREADY_NEIGHBORS


! Check to see if the particles are already neighbors.
      ALREADY_NEIGHBORS = .FALSE.
      DO NEIGH_L = 2, THERMO_NBRHD(LL,1)+1
         IF(L .EQ. THERMO_NBRHD(LL,NEIGH_L)) ALREADY_NEIGHBORS=.TRUE.
      ENDDO

      IF((DIST .LE. NBRHD_SZ) .AND. .NOT.ALREADY_NEIGHBORS) THEN

         THERMO_NBRHD(L,1) = THERMO_NBRHD(L,1) + 1
         THERMO_NBRHD(LL,1) = THERMO_NBRHD(LL,1) + 1

         NEIGH_L = THERMO_NBRHD(L,1)
         NEIGH_LL = THERMO_NBRHD(LL,1)

         IF(NEIGH_L .LE. MN) THEN
            THERMO_NBRHD(L,NEIGH_L+1) = LL
         ELSE
            WRITE(*,1002)L,MN
         ENDIF

         IF(NEIGH_LL .LE. MN) THEN
            THERMO_NBRHD(LL,NEIGH_LL+1) = L
         ELSE
            WRITE(*,1002)LL,MN
         ENDIF

      ENDIF

      RETURN

 1002 FORMAT(/1X,70('*')/,1X,'From: THERMO_NBR -',/&
         ' Message: The number of particles in the thermodynamic',&
         ' neighborhood',/' of particle ',I6,' has exceeded',I3,&
         /1X,70('*'))

      END SUBROUTINE THERMO_NBR
