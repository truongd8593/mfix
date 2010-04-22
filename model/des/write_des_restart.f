!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_DES_RESTART                                      C
!  Purpose: Writing DES data for restart
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 26-Jul-06  C
!  Reviewer:Sreekanth Pannala                         Date: 31-Oct-06  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE WRITE_DES_RESTART

      USE param1      
      USE discretelement
      USE run
      USE des_bc

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER BCV     ! Loop counter for no. of DES_BCMI
      LOGICAL ASSOC   ! If a derived data type is associated

!-----------------------------------------------


      OPEN (UNIT=901,FILE=TRIM(RUN_NAME)//'_DES.RES',FORM='Unformatted',STATUS='unknown')

      REWIND (901)

      WRITE (901) PARTICLES
      WRITE (901) VTP_FINDEX
      WRITE (901) TECPLOT_FINDEX      
      WRITE (901) DTSOLID

      WRITE (901) DES_POS_NEW
      WRITE (901) DES_VEL_NEW
      WRITE (901) OMEGA_NEW

      WRITE (901) DES_RADIUS
      WRITE (901) RO_Sol

      WRITE (901) NEIGHBOURS

! Needed for particle contact history in tangential direction
      WRITE (901) PFT
      WRITE (901) PN
      WRITE (901) PV

! J. Musser DES boundary condition data
      WRITE (901) PIS
      WRITE (901) PEA

! These arrays are allocated only if inlet exists
      IF(DES_MI)THEN
         WRITE (901) DES_MI_TIME
         WRITE (901) MI_FACTOR
         WRITE (901) MI_WINDOW

         DO BCV =1, DES_BCMI

            WRITE (901) PARTICLE_PLCMNT(BCV)

            IF(PARTICLE_PLCMNT(BCV) == 'ORDR')THEN

               ASSOC = ASSOCIATED(MI_ORDER(BCV)%VALUE)
               WRITE (901) ASSOC
               IF(ASSOC)THEN
                  WRITE (901) SIZE(MI_ORDER(BCV)%VALUE)
                  WRITE (901) MI_ORDER(BCV)%VALUE
               ENDIF

               ASSOC = ASSOCIATED(I_OF_MI(BCV)%VALUE)
               WRITE (901) ASSOC
               IF(ASSOC)THEN
                  WRITE (901) SIZE(I_OF_MI(BCV)%VALUE)
                  WRITE (901) I_OF_MI(BCV)%VALUE
               ENDIF

               ASSOC = ASSOCIATED(J_OF_MI(BCV)%VALUE)
               WRITE (901) ASSOC
               IF(ASSOC)THEN
                  WRITE (901) SIZE(J_OF_MI(BCV)%VALUE)
                  WRITE (901) J_OF_MI(BCV)%VALUE
               ENDIF

            ENDIF  
         ENDDO
      ENDIF

      END SUBROUTINE WRITE_DES_RESTART 

