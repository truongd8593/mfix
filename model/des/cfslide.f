!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!
!  Subroutine: CFSLIDE(LL, TANGNT_VREL)
!  Purpose:  Check for Coulombs friction law - calculate sliding
!            friction
!
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04
!  Reviewer:                                          Date:
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

      SUBROUTINE CFSLIDE(LL, TANGNT, PARTICLE_SLIDE, MU)

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
      DOUBLE PRECISION, INTENT(IN) :: TANGNT(3)
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
! Functions
!-----------------------------------------------
      DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT
!-----------------------------------------------

      FTMD = DES_DOTPRDCT(FT(:,LL),FT(:,LL))
      FNMD = DES_DOTPRDCT(FN(:,LL),FN(:,LL))

      IF (FTMD.GT.(MU*MU*FNMD)) THEN
! tangential force based on sliding friction
         PARTICLE_SLIDE = .TRUE.
         IF(DES_DOTPRDCT(TANGNT,TANGNT).EQ.0) THEN
            FT(:,LL) =  MU * FT(:,LL) * SQRT(FNMD/FTMD)
         ELSE
            FT(:,LL) = -MU * SQRT(FNMD) * TANGNT(:)
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
