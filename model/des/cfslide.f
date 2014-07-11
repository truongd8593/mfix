!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!
!  Subroutine: CFSLIDE(LL, V_TANG, PARTICLE_SLIDE, MU)
!  Purpose:  Check for Coulombs friction law - calculate sliding
!            friction
!
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04
!  Reviewer:                                          Date:
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

      SUBROUTINE CFSLIDE(LL, V_TANG, PARTICLE_SLIDE, MU)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE discretelement
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! particle index no.
      INTEGER, INTENT(IN) :: LL
! tangent to the plane of contact
      DOUBLE PRECISION, INTENT(IN) :: V_TANG(3)
! logic set to T when a sliding contact occurs
      LOGICAL, INTENT(INOUT) :: PARTICLE_SLIDE
! Coefficient of friction
      DOUBLE PRECISION, INTENT(IN) :: MU
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! squared magnitude of tangential and normal forces
      DOUBLE PRECISION FTMD, FNMD
!-----------------------------------------------

      FTMD = dot_product(FT(:,LL),FT(:,LL))
      FNMD = dot_product(FN(:,LL),FN(:,LL))

      IF (FTMD.GT.(MU*MU*FNMD)) THEN
! tangential force based on sliding friction
         PARTICLE_SLIDE = .TRUE.
         IF(ALL(V_TANG.EQ.0)) THEN
            FT(:,LL) =  MU * FT(:,LL) * SQRT(FNMD/FTMD)
         ELSE
            FT(:,LL) = -MU * V_TANG(:) * SQRT(FNMD/dot_product(V_TANG,V_TANG))
         ENDIF
      ENDIF

      IF(DEBUG_DES .AND. PARTICLE_SLIDE) THEN
         WRITE(*,'(7X,A)') &
            'FROM CFSLIDE.F ---------->'
         WRITE(*,'(9X,A)') 'PARTICLE_SLIDE = T'
         WRITE(*,'(9X,A,2(ES15.7,X))')&
            'FTMD, mu*FNMD = ', FTMD, MU*FNMD
         WRITE(*,'(7X,A)') '<----------END CFSLIDE.F'
      ENDIF

      RETURN
      END SUBROUTINE CFSLIDE
