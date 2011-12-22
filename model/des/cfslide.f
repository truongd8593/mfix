!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFSLIDE(L, TANGNT_VREL)
!
!  Purpose: DES - calculate sliding between particles                  
!
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CFSLIDE(L, TANGNT,PARTICLE_SLIDE)

      USE param1
      USE discretelement
      IMPLICIT NONE

!-----------------------------------------------
! Local Variables
!-----------------------------------------------
      INTEGER L
      DOUBLE PRECISION FTMD, FNMD, TANGNT(DIMN)
      DOUBLE PRECISION TMP_FT(DIMN), TMP_FN(DIMN)
      logical PARTICLE_SLIDE

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


