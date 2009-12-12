!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFSLIDEWALL(L, TANGNT)                                 C
!>
!!  Purpose: DES - Calculate slide between particles and walls 
!<
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CFSLIDEWALL(L, TANGNT, TMP_FT)
      
      USE discretelement
      USE param1
      IMPLICIT NONE
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
      INTEGER L, K
      DOUBLE PRECISION FTMD, FNMD, TANGNT(DIMN)
      DOUBLE PRECISION TMP_FT(DIMN), TMP_FN(DIMN)

!-----------------------------------------------      
! Functions
!-----------------------------------------------         
      DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT
      
!----------------------------------------------- 

      TMP_FN(:) = FN(L, :)

      FTMD = SQRT(DES_DOTPRDCT(TMP_FT,TMP_FT))
      FNMD = SQRT(DES_DOTPRDCT(TMP_FN,TMP_FN))

      IF (FTMD.GT.(MEW_W*FNMD)) THEN
         PARTICLE_SLIDE = .TRUE.
         IF(DES_DOTPRDCT(TANGNT,TANGNT).EQ.0) THEN
            FT(L,:) =  MEW_W * FNMD * TMP_FT(:)/FTMD
         ELSE
            FT(L,:) = -MEW_W * FNMD * TANGNT(:)
         ENDIF
      ELSE
         FT(L,:) = TMP_FT(:)
      ENDIF

      IF(DEBUG_DES .AND. PARTICLE_SLIDE) THEN
         WRITE(*,'(7X,A)') &
            'FROM CFSLIDEWALL.F ---------->'
         WRITE(*,'(9X,A)') 'PARTICLE_SLIDE = T'
         WRITE(*,'(9X,A,2(ES15.7,X))') &
         'FTMD, mu_w*FNMD = ', FTMD, MEW_W*FNMD
         WRITE(*,'(7X,A)') '<----------END CFSLIDEWALL.F'
      ENDIF

      RETURN
      END SUBROUTINE CFSLIDEWALL


