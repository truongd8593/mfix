!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: READ_DES_RESTART                                       C
!  Purpose: Reading DES data for restart                               
!                                                                      C
!                                                                      C
!  Author: Sreekanth Pannala                          Date: 09-Nov-06  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE READ_DES_RESTART

      USE param1      
      USE run
      USE discretelement
      USE des_bc

      IMPLICIT NONE

!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER BCV    ! Loop counter for no. of DES_BCMI
      INTEGER ALL_SZ ! Size used to allocate derived data types

      LOGICAL ASSOC  ! If a derived data type is associated
!-----------------------------------------------

      OPEN (UNIT=901,FILE=TRIM(RUN_NAME)//'_DES.RES',FORM='Unformatted',STATUS='unknown')

      REWIND (901)

      READ (901) PARTICLES
      READ (901) VTP_FINDEX
      READ (901) TECPLOT_FINDEX
      READ (901) DTSOLID

      READ (901) DES_POS_OLD
      READ (901) DES_VEL_OLD
      READ (901) OMEGA_OLD

      READ (901) DES_RADIUS
      READ (901) RO_Sol

      READ (901) NEIGHBOURS

! Needed for particle contact history in tangential direction
      READ (901) PFT
      READ (901) PN
      READ (901) PV

! J. Musser : DES boundary condition data
      READ (901) PIS
      READ (901) PEA

! These arrays are allocated only if inlet exists
      IF(DES_BCMI)THEN
         READ (901) DES_MI_TIME
         READ (901) MI_FACTOR
         READ (901) MI_WINDOW

         DO BCV =1, DES_BCMI

            READ (901) PARTICLE_PLCMNT(BCV)

            IF(PARTICLE_PLCMNT(BCV) == 'ORDR')THEN

               READ (901) ASSOC
               IF(ASSOC)THEN
                  READ (901) ALL_SZ
                  ALLOCATE( MI_ORDER(BCV)%VALUE( ALL_SZ ) )
                  READ (901) MI_ORDER(BCV)%VALUE
               ENDIF

               READ (901) ASSOC
               IF(ASSOC)THEN
                  READ (901) ALL_SZ
                  ALLOCATE( I_OF_MI(BCV)%VALUE( ALL_SZ ) )
                  READ (901) I_OF_MI(BCV)%VALUE
               ENDIF

               READ (901) ASSOC
               IF(ASSOC)THEN
                  READ (901) ALL_SZ
                  ALLOCATE( J_OF_MI(BCV)%VALUE( ALL_SZ ) )
                  READ (901) J_OF_MI(BCV)%VALUE
               ENDIF

            ENDIF  
         ENDDO
      ENDIF

      END SUBROUTINE READ_DES_RESTART 


