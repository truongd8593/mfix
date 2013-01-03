!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: ADJUST_DT(IER, NIT)                                    C
!  Purpose: Automatically adjust time step                             C
!                                                                      C
!  Author: M. Syamlal                                 Date: FEB-10-97  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!                                                                      C
!  Variables modified:
!                                                                      C
!  Local variables: NONE                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      LOGICAL FUNCTION ADJUST_DT (IER, NIT) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
!//d      USE funits 
      USE run
      USE output
      USE compar      
      USE mpi_utility

      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER IER, NIT 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: STEPS_MIN = 5 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      DOUBLE PRECISION :: NITOS_NEW 
      CHARACTER, DIMENSION(1) :: LINE*70 
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
!-----------------------------------------------
!
!
!//  only do this routine if the root processor
!//SP
!     if (myPE.ne.PE_IO) return
!QX
!     restart with new dt
      RESTART_REACTION = .FALSE.
!end
!
      ADJUST_DT = .FALSE.                     !No need to iterate again
      IF (DT==UNDEFINED .OR. DT<ZERO) RETURN 
      

!     Initialize
      IF (IER == 100) THEN 
         DT_DIR = -1 
         STEPS_TOT = 0 
         NITOS = 0. 
         NIT_TOT = 0 
!
!     Iterate gave converged results
      ELSE IF (IER == 0) THEN 

!AE TIME 041601 Set back the timestep to original size which was
!               halved previously for 2nd order accurate time implementation
!         IF (CN_ON.AND.NSTEP>1) DT = 2.D0*DT      
         IF ((CN_ON.AND.NSTEP>1.AND.RUN_TYPE == 'NEW') .OR. & 
          (CN_ON.AND.RUN_TYPE /= 'NEW' .AND. NSTEP >= (NSTEPRST+1))) &
	     DT = 2.D0*DT      

         IF (STEPS_TOT >= STEPS_MIN) THEN 
            NITOS_NEW = DBLE(NIT_TOT)/(STEPS_TOT*DT) 
            IF (NITOS_NEW > NITOS) DT_DIR = DT_DIR*(-1) 
!
            STEPS_TOT = 0 
            NITOS = NITOS_NEW 
            NIT_TOT = 0 
!
            IF (DT_DIR > 0) THEN 
               DT = MIN(DT_MAX,DT/DT_FAC) 
            ELSE 
               DT = DT*DT_FAC 
            ENDIF 
!//SP
            IF (FULL_LOG.and.myPE.eq.PE_IO) &
            WRITE (*, *) 'DT = ', DT, '  NIT/s = ', NINT(NITOS) 
         ELSE 
            STEPS_TOT = STEPS_TOT + 1 
            NIT_TOT = NIT_TOT + NIT 
         ENDIF 
         ADJUST_DT = .FALSE.                     !No need to iterate again 

!AE TIME 041601 Cut the timestep into half for 2nd order accurate time implementation
!         IF (CN_ON.AND.NSTEP>1) DT = 0.5D0*DT      
         IF ((CN_ON.AND.NSTEP>1.AND.RUN_TYPE == 'NEW') .OR. & 
          (CN_ON.AND.RUN_TYPE /= 'NEW' .AND. NSTEP >= (NSTEPRST+1))) &
	      DT = 0.5D0*DT      


!
!     No convergence in iterate
      ELSE 
!
!AE TIME 041601 Set back the timestep to original size which was
!               halved previously for 2nd order accurate time implementation
!         IF (CN_ON.AND.NSTEP>1) DT = 2.D0*DT      
         IF ((CN_ON.AND.NSTEP>1.AND.RUN_TYPE == 'NEW') .OR. & 
          (CN_ON.AND.RUN_TYPE /= 'NEW' .AND. NSTEP >= (NSTEPRST+1))) &
              DT = 2.D0*DT 

         IF (DT < DT_MIN) THEN 
!//SP
            IF(AUTO_RESTART) THEN
              LINE(1) = 'DT < DT_MIN.  Recovery not possible! Trying Automatic Restart'
	    ELSE
	      LINE(1) = 'DT < DT_MIN.  Recovery not possible!'
	    ENDIF
            IF (FULL_LOG.and.myPE.eq.PE_IO) WRITE (*, *) LINE(1) 
            CALL WRITE_ERROR ('ADJUST_DT', LINE, 1) 
!           CALL MFIX_EXIT(myPE)
         ELSE IF (DT_FAC >= ONE) THEN 
            LINE(1) = 'DT_FAC >= 1.  Recovery not possible!' 
!//SP
            IF (FULL_LOG.and.myPE.eq.PE_IO) WRITE (*, *) LINE(1) 
            CALL WRITE_ERROR ('ADJUST_DT', LINE, 1) 
            CALL MFIX_EXIT(myPE)
         ELSE 
!
            DT_DIR = -1 
            STEPS_TOT = 0 
            NITOS = 0. 
            NIT_TOT = 0 
!
            DT = DT*DT_FAC 
!
!QX smaller dt and reiterating
          if(DT > 0.8d0 * DT_OLD) then

            IF (FULL_LOG) THEN 
!//SP
	      IF(myPE.eq.PE_IO) then
               WRITE (*, '(12X,A,G11.5,9X,A)') ' Dt=', DT, &
                  ' Recovered            :-)' 
	      ENDIF
!
               CALL START_LOG 
               IF(DMP_LOG)WRITE (UNIT_LOG, '(12X,A,G11.5,9X,A)') ' Dt=', DT, &
                  ' Recovered            :-)' 
               CALL END_LOG 
            ENDIF 
!
            CALL RESET_NEW(0) 
            ADJUST_DT = .TRUE.                   !Iterate again with new dt 

         elseif (DT .le. 0.8d0 * DT_OLD) then
            ADJUST_DT = .FALSE.                   !NOT Iterate again with new dt 
            RESTART_REACTION = .TRUE.           ! re-iterate reaction with new dt

            IF (FULL_LOG) THEN 
!//SP
	      IF(myPE.eq.PE_IO) then
               WRITE (*, '(12X,A,G11.5,9X,A)') ' Dt=', DT, &
                  ' Re-iterate reaction            :-(' 
	      ENDIF
!
               CALL START_LOG 
               IF(DMP_LOG)WRITE (UNIT_LOG, '(12X,A,G11.5,9X,A)') ' Dt=', DT, &
                  ' Re-iterate reaction            :-(' 
               CALL END_LOG 
            ENDIF 
         endif
!end

!AE TIME 041601 Cut the timestep into half for 2nd order accurate time implementation
!            IF (CN_ON.AND.NSTEP>1) DT = 0.5D0*DT      
            IF ((CN_ON.AND.NSTEP>1.AND.RUN_TYPE == 'NEW') .OR. & 
              (CN_ON.AND.RUN_TYPE /= 'NEW' .AND. NSTEP >= (NSTEPRST+1))) &
	          DT = 0.5D0*DT      
         ENDIF 
!
!
!       get courant numbers and max p_star for the diverged iteration
!        call get_stats(IER)
!
      ENDIF 

!AIKEDEBUG Get the min. Courant number displayed.
!      call get_stats(IER)

!
      ODT = ONE/DT 

!
      call bcast (dt,root)
      call bcast (odt,root)
!
      RETURN  
      END FUNCTION ADJUST_DT 
      
!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization 
!// 400 Added mpi_utility module and other global reduction (bcast) calls
