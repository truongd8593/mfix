!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!
!  Subroutine: CFSLIDE(L, TANGNT_VREL)
!  Purpose:  Check for Coulombs friction law - calculate sliding 
!            friction
!
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04
!  Reviewer:                                          Date:
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

      SUBROUTINE CFSLIDE(L, TANGNT, PARTICLE_SLIDE)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param1
      USE discretelement
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------      
! particle index no.      
      INTEGER, INTENT(IN) :: L
! tangent to the plane of contact      
      DOUBLE PRECISION, INTENT(IN) :: TANGNT(DIMN)      
! logic set to T when a sliding contact occurs      
      LOGICAL, INTENT(INOUT) :: PARTICLE_SLIDE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! local variables for tangential and normal forces at point of contact
      DOUBLE PRECISION TMP_FT(DIMN), TMP_FN(DIMN)
! magnitude of tangential and normal forces
      DOUBLE PRECISION FTMD, FNMD
!-----------------------------------------------      
! Functions
!-----------------------------------------------     
      DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT
!-----------------------------------------------     

      TMP_FN(:) = FN(L, :)
      TMP_FT(:) = FT(L, :)

      FTMD = SQRT(DES_DOTPRDCT(TMP_FT,TMP_FT))
      FNMD = SQRT(DES_DOTPRDCT(TMP_FN,TMP_FN))

      IF (FTMD.GT.(MEW*FNMD)) THEN
! tangential force based on sliding friction              
         PARTICLE_SLIDE = .TRUE.
         IF(DES_DOTPRDCT(TANGNT,TANGNT).EQ.0) THEN
            FT(L,:) =  MEW * FNMD * TMP_FT(:)/FTMD
         ELSE
            FT(L,:) = -MEW * FNMD * TANGNT(:)
         ENDIF
      ELSE
! no sliding friction so tangential force is not changed              
         FT(L, :) = TMP_FT(:)
      ENDIF

      IF(DEBUG_DES .AND. PARTICLE_SLIDE) THEN
         WRITE(*,'(7X,A)') &
            'FROM CFSLIDE.F ---------->'
         WRITE(*,'(9X,A)') 'PARTICLE_SLIDE = T'
         WRITE(*,'(9X,A,2(ES15.7,X))')&
            'FTMD, mu*FNMD = ', FTMD, MEW*FNMD
         WRITE(*,'(7X,A)') '<----------END CFSLIDE.F'
      ENDIF

      RETURN
      END SUBROUTINE CFSLIDE


