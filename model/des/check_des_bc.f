!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: CHECK_DES_BC                                           !
!                                                                      !
!  Purpose: Check the data provided for the des mass inflow boundary   !
!  condition and flag errors if the data is improper.  This module is  !
!  also used to convert the proveded information into the format       !
!  necessary for the dependent subrountines to function properly.      !
!                                                                      !
!                                                                      !
!  Author: J.Musser                                   Date: 13-Jul-09  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE CHECK_DES_BC

      USE compar
      USE constant
      USE des_bc
      USE discretelement 
      USE funits  
      USE geometry
      USE indices
      USE param
      USE param1
      USE physprop
      USE run

      IMPLICIT NONE

!-----------------------------------------------
! Local parameters
!-----------------------------------------------
      INTEGER, PARAMETER :: DIM_BCTYPE = 4

!-----------------------------------------------
! Local variables
!-----------------------------------------------

      LOGICAL DES_BC_DEFINED(DIMENSION_BC)

      INTEGER BCV, I, J, K         ! Loop Counter
      INTEGER BCV_I

      INTEGER BC_MI ! Total DES mass inlets
      INTEGER BC_MO ! Total DES mass outlets

! the number of particles injected in a solids time step
      DOUBLE PRECISION NPpDT 

!     valid boundary condition types
      CHARACTER*16, DIMENSION(1:DIM_BCTYPE) ::VALID_BC_TYPE = (/&
           'MASS_INFLOW     ', 'MI              ',&
           'MASS_OUTFLOW    ', 'MO              '/)
!-----------------------------------------------

! Initialize

      BC_MI = 0; DES_MI = .FALSE.
      BC_MO = 0
      WRITE(*,'(3X,A)') '---------- START CHECK_DES_BC ---------->'

! Check for des inlet/outlet information:
      CHECK_BC: DO BCV = 1, DIMENSION_BC 
         DES_BC_DEFINED(BCV) = .FALSE. 
         IF (DES_BC_X_w(BCV) /= UNDEFINED) DES_BC_DEFINED(BCV) = .TRUE. 
         IF (DES_BC_X_e(BCV) /= UNDEFINED) DES_BC_DEFINED(BCV) = .TRUE. 
         IF (DES_BC_Y_s(BCV) /= UNDEFINED) DES_BC_DEFINED(BCV) = .TRUE. 
         IF (DES_BC_Y_n(BCV) /= UNDEFINED) DES_BC_DEFINED(BCV) = .TRUE. 
         IF (DES_BC_Z_b(BCV) /= UNDEFINED) DES_BC_DEFINED(BCV) = .TRUE. 
         IF (DES_BC_Z_t(BCV) /= UNDEFINED) DES_BC_DEFINED(BCV) = .TRUE. 

! If a boundary location is specified, verify necessary data
         IF (DES_BC_DEFINED(BCV)) THEN 
            IF(DES_BC_X_w(BCV) == UNDEFINED .OR. &
               DES_BC_X_e(BCV) == UNDEFINED .OR. &
               DES_BC_Y_s(BCV) == UNDEFINED .OR. &
               DES_BC_Y_n(BCV) == UNDEFINED) THEN
               WRITE (UNIT_LOG, 1000) BCV
               WRITE (*, 1000) BCV
               CALL MFIX_EXIT(myPE)
            ENDIF 
            IF(DIMN == 3)THEN
              IF(DES_BC_Z_b(BCV) == UNDEFINED .OR. &
                 DES_BC_Z_t(BCV) == UNDEFINED) THEN
               WRITE (UNIT_LOG, 1000) BCV
               WRITE (*, 1000) BCV
               CALL MFIX_EXIT(myPE)
              ENDIF
            ENDIF
            IF(DES_BC_X_w(BCV) .LT. ZERO .OR. &
               DES_BC_Y_s(BCV) .LT. ZERO .OR. &
               DES_BC_X_e(BCV) .GT. XLENGTH .OR. &
               DES_BC_Y_n(BCV) .GT. YLENGTH .OR. &
               DES_BC_X_w(BCV) .GT. DES_BC_X_e(BCV) .OR. &
               DES_BC_Y_s(BCV) .GT. DES_BC_Y_n(BCV))THEN
               WRITE (UNIT_LOG, 1001) BCV
               WRITE (*, 1001) BCV
               CALL MFIX_EXIT(myPE)
            ENDIF
            IF(DIMN == 3)THEN
               IF(DES_BC_Z_b(BCV) .LT. 0 .OR. &
                  DES_BC_Z_t(BCV) .GT. ZLENGTH .OR. &
                  DES_BC_Z_b(BCV) .GT. DES_BC_Z_t(BCV))THEN
                  WRITE (UNIT_LOG, 1001) BCV
                  WRITE (*, 1001)BCV
                  CALL MFIX_EXIT(myPE)
               ENDIF
            ENDIF

! Require the inlet to be on a boundary/wall of the system
            CALL DES_CHECK_MIO_LOCATION(BCV)

            DO I = 1, DIM_BCTYPE 
               VALID_BC_TYPE(I) = TRIM(VALID_BC_TYPE(I))
               IF (VALID_BC_TYPE(I) == DES_BC_TYPE(BCV)) THEN 
                  IF (MOD(I,2) == 0) DES_BC_TYPE(BCV) = VALID_BC_TYPE(I-1)

! Check if solids phase velocity at the BC plane is specified, if not,
! set solids velocity to zero and flag log file
                  IF(DES_BC_TYPE(BCV) == 'MASS_INFLOW')THEN

                     BC_MI = BC_MI + 1

                     IF(DES_BC_U_s(BCV) == UNDEFINED)THEN
                       DES_BC_U_s(BCV) = ZERO
                       WRITE (UNIT_LOG, 1004)'DES_BC_U_s',BCV
                     ENDIF
                     IF(DES_BC_V_s(BCV) == UNDEFINED)THEN
                       DES_BC_V_s(BCV) = ZERO
                       WRITE (UNIT_LOG, 1004)'DES_BC_V_s',BCV
                     ENDIF
                     IF(DES_BC_W_s(BCV) == UNDEFINED .AND. &
                        DIMN == 3)THEN
                       DES_BC_W_s(BCV) = ZERO
                       WRITE (UNIT_LOG, 1004)'DES_BC_W_s',BCV
                     ENDIF

! Check to verify that one and only one form of mass inlet value is 
! specified. Flag error and stop otherwise
                     IF(DES_BC_VOLFLOW_s(BCV) == UNDEFINED .AND. &
                       DES_BC_MASSFLOW_s(BCV) == UNDEFINED) THEN
                       WRITE (UNIT_LOG, 1006) BCV
                       WRITE (*, 1006) BCV
                       CALL MFIX_EXIT(myPE)
                     ENDIF
                     IF(DES_BC_VOLFLOW_s(BCV) /= UNDEFINED .AND. &
                       DES_BC_MASSFLOW_s(BCV) /= UNDEFINED) THEN
                       WRITE (UNIT_LOG, 1007) BCV
                       WRITE (*, 1007) BCV
                       CALL MFIX_EXIT(myPE)
                     ENDIF
                  ELSE   ! if des_bc_type is not mass_inflow
                     BC_MO = BC_MO + 1
                  ENDIF   ! end if des_bc_type is MI, etc

! gives way to exit loop over des_bc_defined(bcv) when valid_bc_type is t                  
                  CYCLE CHECK_BC

               ENDIF   ! endif des_bc_type is valid
            ENDDO   ! end loop over dim_bc_type

! exit if des_bc_defined but bc_type is not valid            
            WRITE (UNIT_LOG, 1002) BCV, DES_BC_TYPE(BCV)
            WRITE (*, 1002) BCV, DES_BC_TYPE(BCV) 
            WRITE (UNIT_LOG, 1003) VALID_BC_TYPE
            WRITE (*, 1003) VALID_BC_TYPE
            WRITE (*, 1008)
            CALL MFIX_EXIT(myPE)  

         ENDIF   ! end if des_bc_defined(bcv)

      ENDDO CHECK_BC

      WRITE(*,1012) BC_MI, BC_MO

      IF(BC_MI /= 0 .OR. BC_MO /=0)THEN
         DES_MI = .TRUE.
! Allocate necessary arrays for discrete mass inlets
         CALL ALLOCATE_DES_MIO(BC_MI, BC_MO)

! Verify that either the nsqare or grid based neighbor searches are
! used, otherwise flage and exit
         IF((DES_NEIGHBOR_SEARCH == 2) .OR. &
            (DES_NEIGHBOR_SEARCH == 3)) THEN
            WRITE (UNIT_LOG, 1005)
            WRITE (*, 1005)
            CALL MFIX_EXIT(myPE)
         ENDIF
      ENDIF


! If one or more discrete mass inlets exist, calculate necessary data 
      IF(BC_MI/=0)THEN

! Check each discrete mass inlet for necessary data
         BCV_I = 1
         DO BCV = 1, DIMENSION_BC 

           IF(DES_BC_DEFINED(BCV) .AND. &
              (DES_BC_TYPE(BCV) == 'MASS_INFLOW'))THEN

               DES_BC_MI_ID(BCV_I) = BCV

! convert mass flow rate or volumetric flow rate to particles per
! soldis time step               
               IF(DES_BC_MASSFLOW_s(BCV) /= UNDEFINED)THEN
                  NPpDT = (DES_BC_MASSFLOW_s(BCV) * DTSOLID) / &
                       (PI/6.d0 * D_P0(1)**3 * RO_S(1))
               ELSEIF(DES_BC_VOLFLOW_s(BCV) /= UNDEFINED)THEN
                  NPpDT = (DES_BC_VOLFLOW_s(BCV) * DTSOLID) / &
                       (PI/6.d0 * D_P0(1)**3)
               ENDIF

               IF(NPpDT .LT. 1)THEN
                  PI_FACTOR(BCV_I) = FLOOR(real(1.d0/NPpDT))
                  PI_COUNT(BCV_I) = 1
               ELSE
                  PI_FACTOR(BCV_I) = 1
                  PI_COUNT(BCV_I) = CEILING(real(NPpDT))
               ENDIF
! Calculate des mass inlet time; time between injection
               DES_MI_TIME(BCV_I) = TIME +&
                  dble(PI_FACTOR(BCV_I)) * DTSOLID 

               WRITE(*,1013) BCV, NPpDT, PI_FACTOR(BCV_I),&
                  PI_COUNT(BCV_I), DES_MI_TIME(BCV_I)

! Classify boundary condition and verify appropriate solids velocity
               CALL DES_MI_CLASSIFY(BCV_I, BCV)

! Determine the computational cells near the inlet
               CALL DES_MI_CELLS(BCV, BCV_I)

               BCV_I = BCV_I + 1 
            ENDIF
         ENDDO
      ENDIF

! If the system is started without any particles and an inlet is not
! specified, the run is aborted.
      IF(BC_MI == 0 .AND. PARTICLES == 0)THEN
         WRITE(UNIT_LOG, 1009)
         WRITE(*, 1009)
         CALL MFIX_EXIT(myPE)
      ELSEIF(PARTICLES == 0)THEN
         WRITE(*,'(5X,A)') &
            'Run initiated with no particles in the system'
      ENDIF

! Check MAX_PIS requirements
      IF(BC_MI == 0 .AND. MAX_PIS == UNDEFINED_I)THEN
         WRITE(*,'(5X,A)')'Setting MAX_PIS = PARTICLES'
         MAX_PIS = PARTICLES
      ELSEIF(BC_MI /= 0 .AND. MAX_PIS == UNDEFINED_I)THEN
         WRITE(UNIT_LOG, 1010)
         WRITE(*, 1010)
         CALL MFIX_EXIT(myPE)
      ENDIF
      IF(BC_MO/=0)THEN


! Check each discrete mass outlet for necessary data
         BCV_I = 1
         DO BCV = 1, DIMENSION_BC 

           IF(DES_BC_DEFINED(BCV) .AND. &
              (DES_BC_TYPE(BCV) == 'MASS_OUTFLOW'))THEN

               DES_BC_MO_ID(BCV_I) = BCV

               CALL DES_MO_CLASSIFY(BCV_I, BCV)

               BCV_I = BCV_I + 1 
            ENDIF
         ENDDO
      ENDIF




      WRITE(*,'(3X,A)') '<---------- END CHECK_DES_BC ----------'



 1000 FORMAT(/1X,70('*')//&
         ' From: CHECK_DES_BC -',/&
         ' Message: Insufficient DES boundary condition infomation',/&
         ' Check boundary condition number: ',I3,/&
         1X,70('*')/)

 1001 FORMAT(/1X,70('*')//&
         ' From: CHECK_DES_BC -',/&
         ' Message: Improper DES boundary condition information',/&
         ' Check boundary condition number: ',I3,/&
         1X,70('*')/)

 1002 FORMAT(/1X,70('*')//' From: CHECK_DES_BC -',/&
         ' Message: Illegal BC_TYPE for boundary condition ',I3,/&
         ' BC_TYPE = ',A,' and the valid types are:') 
 1003 FORMAT(5X,A16)

 1004 FORMAT(/1X,70('*')//&
         ' From: CHECK_DES_BC -',/&
         ' Message: ',A,' not specified in mfix.dat file for BC ',I3,/&
         ' Setting value to ZERO',/&
         1X,70('*')/)

 1005 FORMAT(/1X,70('*')//&
         ' From: CHECK_DES_BC -',/&
         ' Message: Currently, DES_MI only supports the NSQARE and',/&
         ' GRID_BASED_NEIGHBOR search methods.',/&
         1X,70('*')/)

 1006 FORMAT(/1X,70('*')//&
         ' From: CHECK_DES_BC -',/&
         ' Message: No mass inlet value set for the DES_MI for BC ',I3,/&
         ' DES_BC_VOLFLOW_s or DES_BC_MASSFLOW_s must be specified',/&
         1X,70('*')/)

 1007 FORMAT(/1X,70('*')//&
       ' From: CHECK_DES_BC -',/&
       ' Message: To many mass inlet values set for the DES_MI for BC ',I3,/&
       ' DES_BC_VOLFLOW_s or DES_BC_MASSFLOW_s must be specified',/&
         1X,70('*')/)

 1008 FORMAT(/1X,70('*')/) 

 1009 FORMAT(/1X,70('*')//&
       ' From: CHECK_DES_BC -',/&
       ' Message: The system was initiated with no particles',&
       ' and no DES inlet.',/&
       ' The run is being terminated.',/&
         1X,70('*')/)

 1010 FORMAT(/1X,70('*')//&
       ' From: CHECK_DES_BC -',/&
       ' Message: The maximum number of particles permitted',&
       ' in the system (MAX_PIS)',/&
       ' must be set in the mfix.dat file if an inlet is specified.',/&
         1X,70('*')/)

 1012 FORMAT(5X,'No. of mass inlet BC = ', I4,/,&
         5X,'No. of mass outlet BC = ', I4)

 1013 FORMAT(5X,'For mass inlet BC: ', I3,/,&
         7X,'No. particles injected per solids time step = ', ES15.8,/,&
         7X,'PI_FACTOR = ', I,' PI_COUNT = ', I5,/,&
         7X,'start DES_MI_TIME = ', ES15.8)

      RETURN
      END SUBROUTINE CHECK_DES_BC



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: DES_CHECK_MIO_LOCATION                                  !
!                                                                      !
!  Purpose:  This subroutine verifies that the location of the DES Mass!
!  Inlet was specified to on a wall of the system.  If an inlet is     !
!  located elsewhere, an error is flagged and the prgram terminated.   !
!                                                                      !
!  Author: J.Musser                                   Date: 17-Aug-09  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE DES_CHECK_MIO_LOCATION(BCV)

      USE compar
      USE constant
      USE des_bc
      USE discretelement
      USE funits
      USE geometry
      USE param1
      USE physprop
      USE run
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER BCV         ! Boundary Condition 
!-----------------------------------------------

      IF(DIMN == 2) THEN
! Check verticle mass inlet
         IF(DES_BC_X_w(BCV) == DES_BC_X_e(BCV)) THEN
            IF(DES_BC_X_w(BCV) /= ZERO .AND. &
               DES_BC_X_w(BCV) /= XLENGTH)THEN
               WRITE (UNIT_LOG, 1100)BCV; WRITE (*, 1100)BCV
               CALL MFIX_EXIT(myPE)
            ENDIF
         ENDIF 
! Check horizontal mass inlet
         IF(DES_BC_Y_s(BCV) == DES_BC_Y_n(BCV)) THEN
            IF(DES_BC_Y_s(BCV) /= ZERO .AND. &
               DES_BC_Y_s(BCV) /= YLENGTH)THEN
               WRITE (UNIT_LOG, 1100)BCV; WRITE (*, 1100)BCV
               CALL MFIX_EXIT(myPE)
            ENDIF
         ENDIF

! Require the inlet be an edge in 2D not a point
         IF(DES_BC_X_w(BCV) == DES_BC_X_e(BCV) .AND. &
            DES_BC_Y_s(BCV) == DES_BC_Y_n(BCV)) THEN
            WRITE (UNIT_LOG, 1101)&
            'DES mass inlet must larger than a single point.',BCV
            WRITE (*, 1101)&
            'DES mass inlet must larger than a single point.',BCV
            CALL MFIX_EXIT(myPE)
         ENDIF
! Require the inlet be an edge in 2D not a plane
         IF(DES_BC_X_w(BCV) /= DES_BC_X_e(BCV) .AND. &
            DES_BC_Y_s(BCV) /= DES_BC_Y_n(BCV)) THEN
            WRITE (UNIT_LOG, 1101)&
            'For DIMN=2, DES mass inlet cannot be an area.',BCV
            WRITE (*, 1101)&
            'For DIMN=2, DES mass inlet cannot be an area.',BCV
            CALL MFIX_EXIT(myPE)
         ENDIF

      ELSE   ! if dimn != 2

         IF(DES_BC_X_w(BCV) == DES_BC_X_e(BCV)) THEN
            IF(DES_BC_X_w(BCV) /= ZERO .AND. DES_BC_X_w(BCV) /= XLENGTH)THEN
               IF(DES_BC_Y_s(BCV) == DES_BC_Y_n(BCV)) THEN
!  Require the inlet be a slit/edge on the XZ-face
                  IF(DES_BC_Y_s(BCV) /= ZERO .AND. &
                     DES_BC_Y_s(BCV) /= YLENGTH)THEN
                     WRITE (UNIT_LOG, 1100)BCV; WRITE (*, 1100)BCV
                     CALL MFIX_EXIT(myPE)
                  ENDIF
               ENDIF
               IF(DES_BC_Z_b(BCV) == DES_BC_Z_t(BCV)) THEN
! Require the inlet be a slit/edge on the XY-face
                  IF(DES_BC_Z_b(BCV) /= ZERO .AND. &
                     DES_BC_Z_b(BCV) /= ZLENGTH)THEN
                     WRITE (UNIT_LOG, 1100)BCV; WRITE (*, 1100)BCV
                     CALL MFIX_EXIT(myPE)
                  ENDIF
               ENDIF
! Require an area inlet to be on the YZ-face and the face must be at
! one of the x-boundary edges
               IF(DES_BC_Y_s(BCV) /= DES_BC_Y_n(BCV) .AND. & 
                  DES_BC_Z_b(BCV) /= DES_BC_Z_t(BCV)) THEN
                  WRITE (UNIT_LOG, 1100)BCV; WRITE (*, 1100)BCV
                  CALL MFIX_EXIT(myPE)
               ENDIF
            ENDIF
         ENDIF

         IF(DES_BC_Y_s(BCV) == DES_BC_Y_n(BCV)) THEN
            IF(DES_BC_Y_s(BCV) /= ZERO .AND. DES_BC_Y_s(BCV) /= YLENGTH)THEN
               IF(DES_BC_X_w(BCV) == DES_BC_X_e(BCV)) THEN
! Require the inlet be a slit/edge on the YZ-face
                  IF(DES_BC_X_w(BCV) /= ZERO .AND. &
                     DES_BC_X_w(BCV) /= XLENGTH)THEN
                     WRITE (UNIT_LOG, 1100)BCV; WRITE (*, 1100)BCV
                     CALL MFIX_EXIT(myPE)
                  ENDIF
               ENDIF
               IF(DES_BC_Z_b(BCV) == DES_BC_Z_t(BCV)) THEN
! Require the inlet be a slit/edge on the XY-face
                  IF(DES_BC_Z_b(BCV) /= ZERO .AND. &
                     DES_BC_Z_b(BCV) /= ZLENGTH)THEN
                     WRITE (UNIT_LOG, 1100)BCV; WRITE (*, 1100)BCV
                     CALL MFIX_EXIT(myPE)
                  ENDIF
               ENDIF
               IF(DES_BC_X_w(BCV) /= DES_BC_X_e(BCV) .AND. & 
                  DES_BC_Z_b(BCV) /= DES_BC_Z_t(BCV)) THEN
! Require an area inlet to be on the XZ-face and the face must be at
! one of the y-boundary edges
                  WRITE (UNIT_LOG, 1100)BCV; WRITE (*, 1100)BCV
                  CALL MFIX_EXIT(myPE)
               ENDIF
            ENDIF
         ENDIF

         IF(DES_BC_Z_b(BCV) == DES_BC_Z_t(BCV)) THEN
            IF(DES_BC_Z_b(BCV) /= ZERO .AND. DES_BC_Z_b(BCV) /= ZLENGTH)THEN
               IF(DES_BC_X_w(BCV) == DES_BC_X_e(BCV)) THEN
! Require the inlet be a slit/edge on the YZ-face
                  IF(DES_BC_X_w(BCV) /= ZERO .AND. &
                     DES_BC_X_w(BCV) /= XLENGTH)THEN
                     WRITE (UNIT_LOG, 1100)BCV; WRITE (*, 1100)BCV
                     CALL MFIX_EXIT(myPE)
                  ENDIF
               ENDIF
               IF(DES_BC_Y_s(BCV) == DES_BC_Y_n(BCV)) THEN
! Require the inlet be a slit/edge on the XZ-face
                  IF(DES_BC_Y_s(BCV) /= ZERO .AND. &
                     DES_BC_Y_s(BCV) /= YLENGTH)THEN
                     WRITE (UNIT_LOG, 1100)BCV; WRITE (*, 1100)BCV
                     CALL MFIX_EXIT(myPE)
                  ENDIF
                  IF(DES_BC_X_w(BCV) /= DES_BC_X_e(BCV) .AND. & 
                     DES_BC_Y_s(BCV) /= DES_BC_Y_n(BCV)) THEN
! Require an area inlet to be on the XY-face and the face must be at
! one of the z-boundary edges
                     WRITE (UNIT_LOG, 1100)BCV; WRITE (*, 1100)BCV
                     CALL MFIX_EXIT(myPE)
                  ENDIF
               ENDIF
            ENDIF
         ENDIF

      ENDIF ! end if dimn == 2 


! The following checks target potential bc problems in 3D systems   
! ----------------------------------------
      IF(DES_BC_X_w(BCV) == DES_BC_X_e(BCV) .AND. &
         DES_BC_Y_s(BCV) == DES_BC_Y_n(BCV) .AND. &
         DES_BC_Z_b(BCV) == DES_BC_Z_t(BCV))THEN
! Prohibit point-based mass inlet
         WRITE (UNIT_LOG, 1101)&
         'DES mass inlet must larger than a single point.',BCV
         WRITE (*, 1101)&
         'DES mass inlet must larger than a single point.',BCV
         CALL MFIX_EXIT(myPE)
      ENDIF
      IF(DES_BC_X_w(BCV) /= DES_BC_X_e(BCV) .AND. &
         DES_BC_Y_s(BCV) /= DES_BC_Y_n(BCV) .AND. &
         DES_BC_Z_b(BCV) /= DES_BC_Z_t(BCV))THEN
! Prohibit volume-based mass inlet
         WRITE (UNIT_LOG, 1101)&
         'For DIMN=3, DES mass inlet cannot be a volume.',BCV
         WRITE (*, 1101)&
         'For DIMN=3, DES mass inlet cannot be a volume.',BCV
         CALL MFIX_EXIT(myPE)
      ENDIF

      IF((DES_BC_X_w(BCV) == DES_BC_X_e(BCV) .AND. &
          DES_BC_Y_s(BCV) == DES_BC_Y_n(BCV) .AND. &
             ((DES_BC_X_w(BCV) == ZERO .AND. DES_BC_Y_s(BCV) == ZERO) &
            .OR. (DES_BC_X_w(BCV) == ZERO .AND. DES_BC_Y_s(BCV) == YLENGTH) &
            .OR. (DES_BC_X_w(BCV) == XLENGTH .AND. DES_BC_Y_s(BCV) == ZERO) &
            .OR. (DES_BC_X_w(BCV) == XLENGTH .AND. DES_BC_Y_s(BCV) == YLENGTH))) .OR. &
         (DES_BC_X_w(BCV) == DES_BC_X_e(BCV) .AND. &
          DES_BC_Z_b(BCV) == DES_BC_Z_t(BCV) .AND. &
            ((DES_BC_X_w(BCV) == ZERO .AND. DES_BC_Z_b(BCV) == ZERO) &
            .OR. (DES_BC_X_w(BCV) == ZERO .AND. DES_BC_Z_b(BCV) == ZLENGTH) &
            .OR. (DES_BC_X_w(BCV) == XLENGTH .AND. DES_BC_Z_b(BCV) == ZERO) &
            .OR. (DES_BC_X_w(BCV) == XLENGTH .AND. DES_BC_Z_b(BCV) == ZLENGTH))) .OR. &
         (DES_BC_Y_s(BCV) == DES_BC_Y_n(BCV) .AND. &
          DES_BC_Z_b(BCV) == DES_BC_Z_t(BCV) .AND. &
             ((DES_BC_Y_s(BCV) == ZERO .AND. DES_BC_Z_b(BCV) == ZERO) &
            .OR. (DES_BC_Y_s(BCV) == ZERO .AND. DES_BC_Z_b(BCV) == ZLENGTH) &
            .OR. (DES_BC_Y_s(BCV) == YLENGTH .AND. DES_BC_Z_b(BCV) == ZERO) &
            .OR. (DES_BC_Y_s(BCV) == YLENGTH .AND. DES_BC_Z_b(BCV) == ZLENGTH))))THEN
! Prohibit DES mass inlet along domain edges
         WRITE (UNIT_LOG, 1101)&
         'DES mass inlet cannot be positioned in a corner.',BCV
         WRITE (*, 1101)&
         'DES mass inlet cannot be positioned in a corner.',BCV
         CALL MFIX_EXIT(myPE)
      ENDIF

 1100 FORMAT(/1X,70('*')//&
         ' From: DES_CHECK_MIO_LOCATION -',/&
         ' Message: DES boundary condition ',I4,&
         ' must be specified on a wall.',/1X,70('*'))

 1101 FORMAT(/1X,70('*')//&
         ' From: DES_CHECK_MIO_LOCATION -',/&
         ' Message: ',A,/&
         ' Check boundary condition ',I4,&
         ' in the mfix.dat file.',/1X,70('*'))

      RETURN
      END SUBROUTINE DES_CHECK_MIO_LOCATION



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: ALLOCATE_DES_MIO                                       !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!  Author: J.Musser                                   Date: 17-Aug-09  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE ALLOCATE_DES_MIO(BC_MI, BC_MO)

      USE des_bc
      USE discretelement

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER BC_MI ! Number of valid discrete mass inlets
      INTEGER BC_MO ! Number of valid discrete mass outlets
      INTEGER I     ! Loop counter for no. of BC_MI
!-----------------------------------------------

      IF(BC_MI /= 0)THEN

! Boundary condition ID array
         Allocate( DES_BC_MI_ID (BC_MI) )

! Particle injection factor
         Allocate( PI_FACTOR (BC_MI) )

! Particle injection count (injection number)
         Allocate( PI_COUNT (BC_MI) )

! Particle injection time scale
         Allocate( DES_MI_TIME (BC_MI) )

! Boundary classification
         Allocate( DES_MI_CLASS (BC_MI) )
         Allocate( PARTICLE_PLCMNT (BC_MI) )

! Order inlet condition variables
! (only needed if particle_plcmt is assigned 'ordr')
         Allocate( MI_FACTOR (BC_MI) )
         Allocate( MI_WINDOW (BC_MI) )
         Allocate( MI_ORDER (BC_MI) )   ! type dmi
         Allocate( I_OF_MI ( BC_MI) )   ! type dmi
         Allocate( J_OF_MI ( BC_MI) )   ! type dmi

! Initializiation
         DO I = 1,BC_MI
            NULLIFY( MI_ORDER(I)%VALUE )
            NULLIFY( I_OF_MI(I)%VALUE )
            NULLIFY( J_OF_MI(I)%VALUE )
         ENDDO

! Grid search loop counter array; 6 = no. of faces
         Allocate(  GS_ARRAY (BC_MI, 6) )
      ENDIF


      IF(BC_MO /= 0)THEN
         DES_MO_X = .FALSE.
         DES_MO_Y = .FALSE.
         DES_MO_Z = .FALSE.

! Boundary Condition ID array
         Allocate( DES_BC_MO_ID (BC_MO) )

         Allocate( DES_MO_CLASS (BC_MO) )
      ENDIF

      RETURN
      END SUBROUTINE ALLOCATE_DES_MIO



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: DES_MI_CLASSIFY                                        !
!                                                                      !
!  Purpose:  This subroutine is used to give a classification to the   !
!  inlet.  The classification is used in the placement of new particles!
!  in ghost cells behind the inlet.                                    !
!                                                                      !
!  Author: J.Musser                                   Date: 14-Aug-09  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE DES_MI_CLASSIFY(BCV_I, BCV)

      USE compar
      USE constant
      USE des_bc
      USE discretelement 
      USE funits  
      USE geometry
      USE param1
      USE physprop
      USE run

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! passed arguments giving index information of the boundary
      INTEGER BCV_I, BCV
! the length of each side of the inlet boundary
      DOUBLE PRECISION LEN1, LEN2
!-----------------------------------------------


      IF(DIMN == 2)THEN   ! 2D domain

! Check verticle mass inlet: 2D
! ----------------------------------------
         IF(DES_BC_X_w(BCV) == DES_BC_X_e(BCV))THEN
! Verify that a velocity has been specifed for DES_BC_U_s(BCV)
            IF(DES_BC_U_s(BCV) == ZERO)THEN
               WRITE (UNIT_LOG, 1200)BCV,'U_s'; WRITE (*, 1200)BCV,'U_s'
               CALL MFIX_EXIT(myPE)
            ENDIF
            IF(DES_BC_X_w(BCV) == ZERO)THEN
! Correct velocity direction, if wrong
               IF(DES_BC_U_s(BCV) .LT. ZERO) THEN 
                  DES_BC_U_s(BCV) = -DES_BC_U_s(BCV)
                  WRITE(*,1201) 'XW'
               ENDIF
! Classify the boundary condition
               DES_MI_CLASS(BCV_I) = 'XW'
            ENDIF
            IF(DES_BC_X_w(BCV) == XLENGTH)THEN
! Correct velocity direction, if wrong
               IF(DES_BC_U_s(BCV) .GT. ZERO) THEN
                  DES_BC_U_s(BCV) = -DES_BC_U_s(BCV)
                  WRITE(*,1201) 'XE'
               ENDIF
! Classify the boundary condition
               DES_MI_CLASS(BCV_I) = 'XE'
            ENDIF
            LEN1 = ABS(DES_BC_Y_n(BCV) - DES_BC_Y_s(BCV))

! This subroutine determines the pattern that the particles will need to
! enter the system, if any.
            CALL DES_MI_LAYOUT(LEN1, ZERO, DES_BC_U_s(BCV), BCV_I, BCV)
         ENDIF

! Check horizontal mass inlet: 2D 
! ----------------------------------------
         IF(DES_BC_Y_s(BCV) == DES_BC_Y_n(BCV))THEN
! Verify that a velocity has been specifed for DES_BC_V_s(BCV)
            IF(DES_BC_V_s(BCV) == ZERO)THEN
               WRITE (UNIT_LOG, 1200)BCV,'V_s'; WRITE (*, 1200)BCV,'V_s'
               CALL MFIX_EXIT(myPE)
            ENDIF
            IF(DES_BC_Y_s(BCV) == ZERO)THEN
! Correct velocity direction, if wrong
               IF(DES_BC_V_s(BCV) .LT. ZERO) THEN
                  DES_BC_V_s(BCV) = -DES_BC_V_s(BCV)
                  WRITE(*,1201) 'YS'
               ENDIF
! Classify the boundary condition
               DES_MI_CLASS(BCV_I) = 'YS'
            ENDIF
            IF(DES_BC_Y_s(BCV) == YLENGTH)THEN
! Correct velocity direction, if wrong
               IF(DES_BC_V_s(BCV) .GT. ZERO) THEN
                  DES_BC_V_s(BCV) = -DES_BC_V_s(BCV)
                  WRITE(*,1201) 'YN'
               ENDIF
! Classify the boundary condition
               DES_MI_CLASS(BCV_I) = 'YN'
            ENDIF
            LEN1 = DES_BC_X_e(BCV) - DES_BC_X_w(BCV)

! This subroutine determines the pattern that the particles will need to
! enter the system, if any.
            CALL DES_MI_LAYOUT(LEN1, ZERO, DES_BC_V_s(BCV), BCV_I, BCV)
         ENDIF

      ELSE   ! 3D domain

! Check mass inlet on XZ face: 3D 
! ----------------------------------------
! see comments following the yz face section for elucidation on this
! series of if statements      
         IF(DES_BC_Y_s(BCV) == DES_BC_Y_n(BCV) .AND. &
            DES_BC_X_w(BCV) /= DES_BC_X_e(BCV) .AND. &
            DES_BC_Z_b(BCV) /= DES_BC_Z_t(BCV) &
         .OR. &
            DES_BC_Y_s(BCV) == DES_BC_Y_n(BCV) .AND. &
            DES_BC_X_w(BCV) == DES_BC_X_e(BCV) .AND. &
            DES_BC_X_w(BCV) /= ZERO .AND. &
            DES_BC_X_w(BCV) /= XLENGTH &
         .OR. &
            DES_BC_Y_s(BCV) == DES_BC_Y_n(BCV) .AND. &
            DES_BC_Z_b(BCV) == DES_BC_Z_t(BCV) .AND. &
            DES_BC_Z_b(BCV) /= ZERO .AND. &
            DES_BC_Z_b(BCV) /= ZLENGTH)THEN

! Verify that a velocity has been specifed for DES_BC_V_s(BCV)
            IF(DES_BC_V_s(BCV) == ZERO)THEN
               WRITE (UNIT_LOG, 1200)BCV,'V_s'; WRITE (*, 1200)BCV,'V_s'
               CALL MFIX_EXIT(myPE)
            ENDIF
            IF(DES_BC_Y_s(BCV) == ZERO)THEN
! Correct velocity direction, if wrong
               IF(DES_BC_V_s(BCV) .LT. ZERO) THEN
                  DES_BC_V_s(BCV) = -DES_BC_V_s(BCV)
                  WRITE(*,1201) 'XZs'
               ENDIF
! Classify the boundary condition
               DES_MI_CLASS(BCV_I) = 'XZs'
            ENDIF
            IF(DES_BC_Y_s(BCV) == YLENGTH)THEN
! Correct velocity direction, if wrong
               IF(DES_BC_V_s(BCV) .GT. ZERO) THEN
                  DES_BC_V_s(BCV) = -DES_BC_V_s(BCV)
                  WRITE(*,1201) 'XZn'
               ENDIF
! Classify the boundary condition
               DES_MI_CLASS(BCV_I) = 'XZn'
            ENDIF

! If a slit flow inlet is specified near an edge, expand the slit to 
! the width of one particle diameter :
! see comments in the yz face section for further explanation   
            IF(DES_BC_X_w(BCV) == DES_BC_X_e(BCV))THEN
               IF(((DES_BC_X_w(BCV) - D_P0(1)*HALF) .GE. ZERO) .AND. &
                  ((DES_BC_X_w(BCV) + D_P0(1)*HALF) .LE. XLENGTH))THEN
                  DES_BC_X_w(BCV) = DES_BC_X_w(BCV) - D_P0(1)*HALF
                  DES_BC_X_e(BCV) = DES_BC_X_e(BCV) + D_P0(1)*HALF
               ELSEIF((DES_BC_X_w(BCV) - D_P0(1)*HALF) .LT. ZERO)THEN
                  DES_BC_X_w(BCV) = ZERO
                  DES_BC_X_e(BCV) = D_P0(1)
               ELSEIF((DES_BC_X_e(BCV) + D_P0(1)*HALF) .GT. XLENGTH)THEN
                  DES_BC_X_w(BCV) = XLENGTH - D_P0(1)
                  DES_BC_X_e(BCV) = XLENGTH 
               ENDIF
            ENDIF
            IF(DES_BC_Z_b(BCV) == DES_BC_Z_t(BCV))THEN
               IF(((DES_BC_Z_b(BCV) - D_P0(1)*HALF) .GE. ZERO) .AND. &
                  ((DES_BC_Z_b(BCV) + D_P0(1)*HALF) .LE. ZLENGTH))THEN
                  DES_BC_Z_b(BCV) = DES_BC_Z_b(BCV) - D_P0(1)*HALF
                  DES_BC_Z_t(BCV) = DES_BC_Z_t(BCV) + D_P0(1)*HALF
               ELSEIF((DES_BC_Z_b(BCV) - D_P0(1)*HALF) .LT. ZERO)THEN
                  DES_BC_Z_b(BCV) = ZERO
                  DES_BC_Z_t(BCV) = D_P0(1)
               ELSEIF((DES_BC_Z_t(BCV) + D_P0(1)*HALF) .GT. ZLENGTH)THEN
                  DES_BC_Z_b(BCV) = ZLENGTH - D_P0(1)
                  DES_BC_Z_t(BCV) = ZLENGTH 
               ENDIF
            ENDIF
! This subroutine determines the pattern that the particles will need to
! enter the system, if any.
            LEN1 = DES_BC_X_e(BCV) - DES_BC_X_w(BCV)
            LEN2 = DES_BC_Z_t(BCV) - DES_BC_Z_b(BCV)
            CALL DES_MI_LAYOUT(LEN1, LEN2, DES_BC_V_s(BCV), BCV_I, BCV)
         ENDIF   
! End check mass inlet on XZ face: 3D 

! Check mass inlet on XY face: 3D 
! ----------------------------------------
! see comments following the yz face section for elucidation on this
! series of if statements 
         IF(DES_BC_Z_b(BCV) == DES_BC_Z_t(BCV) .AND. &
            DES_BC_X_w(BCV) /= DES_BC_X_e(BCV) .AND. &
            DES_BC_Y_s(BCV) /= DES_BC_Y_n(BCV) &
         .OR. &
            DES_BC_Z_b(BCV) == DES_BC_Z_t(BCV) .AND. &
            DES_BC_X_w(BCV) == DES_BC_X_e(BCV) .AND. &
            DES_BC_X_w(BCV) /= ZERO .AND. &
            DES_BC_X_w(BCV) /= XLENGTH &
         .OR. &
            DES_BC_Z_b(BCV) == DES_BC_Z_t(BCV) .AND. &
            DES_BC_Y_s(BCV) == DES_BC_Y_n(BCV) .AND. &
            DES_BC_Y_s(BCV) /= ZERO .AND. &
            DES_BC_Y_s(BCV) /= YLENGTH)THEN

! Verify that a velocity has been specifed for DES_BC_W_s(BCV)
            IF(DES_BC_W_s(BCV) == ZERO)THEN
               WRITE (UNIT_LOG, 1200)BCV,'W_s'; WRITE (*, 1200)BCV,'W_s'
               CALL MFIX_EXIT(myPE)
            ENDIF
            IF(DES_BC_Z_b(BCV) == ZERO)THEN
! Correct velocity direction, if wrong
               IF(DES_BC_W_s(BCV) .LT. ZERO) THEN
                  DES_BC_W_s(BCV) = -DES_BC_W_s(BCV)
                  WRITE(*,1201) 'XYb'
               ENDIF                  
! Classify the boundary condition
               DES_MI_CLASS(BCV_I) = 'XYb'
            ENDIF
            IF(DES_BC_Z_b(BCV) == ZLENGTH)THEN
! Correct velocity direction, if wrong
               IF(DES_BC_W_s(BCV) .GT. ZERO) THEN
                  DES_BC_W_s(BCV) = -DES_BC_W_s(BCV)
                  WRITE(*,1201) 'XYt'
               ENDIF                  
! Classify the boundary condition
               DES_MI_CLASS(BCV_I) = 'XYt'
            ENDIF

! If a slit flow inlet is specified near an edge, expand the slit to 
! the width of one particle diameter :
! see comments in the yz face section for further explanation            
            IF(DES_BC_X_w(BCV) == DES_BC_X_e(BCV))THEN
               IF(((DES_BC_X_w(BCV) - D_P0(1)*HALF) .GE. ZERO) .AND. &
                  ((DES_BC_X_w(BCV) + D_P0(1)*HALF) .LE. XLENGTH))THEN
                  DES_BC_X_w(BCV) = DES_BC_X_w(BCV) - D_P0(1)*HALF
                  DES_BC_X_e(BCV) = DES_BC_X_e(BCV) + D_P0(1)*HALF
               ELSEIF((DES_BC_X_w(BCV) - D_P0(1)*HALF) .LT. ZERO)THEN
                  DES_BC_X_w(BCV) = ZERO
                  DES_BC_X_e(BCV) = D_P0(1)
               ELSEIF((DES_BC_X_e(BCV) + D_P0(1)*HALF) .GT. XLENGTH)THEN
                  DES_BC_X_w(BCV) = XLENGTH - D_P0(1)
                  DES_BC_X_e(BCV) = XLENGTH 
               ENDIF
            ENDIF
            IF(DES_BC_Y_s(BCV) == DES_BC_Y_n(BCV))THEN
               IF(((DES_BC_Y_s(BCV) - D_P0(1)*HALF) .GE. ZERO) .AND. &
                  ((DES_BC_Y_s(BCV) + D_P0(1)*HALF) .LE. YLENGTH))THEN
                  DES_BC_Y_s(BCV) = DES_BC_Y_s(BCV) - D_P0(1)*HALF
                  DES_BC_Y_n(BCV) = DES_BC_Y_n(BCV) + D_P0(1)*HALF
               ELSEIF((DES_BC_Y_s(BCV) - D_P0(1)*HALF) .LT. ZERO)THEN
                  DES_BC_Y_s(BCV) = ZERO
                  DES_BC_Y_n(BCV) = D_P0(1)
               ELSEIF((DES_BC_Y_n(BCV) + D_P0(1)*HALF) .GT. YLENGTH)THEN
                  DES_BC_Y_s(BCV) = YLENGTH - D_P0(1)
                  DES_BC_Y_n(BCV) = YLENGTH 
               ENDIF
            ENDIF
! This subroutine determines the pattern that the particles will need to
! enter the system, if any.
            LEN1 = (DES_BC_X_e(BCV) - DES_BC_X_w(BCV))
            LEN2 = (DES_BC_Y_n(BCV) - DES_BC_Y_s(BCV))
            CALL DES_MI_LAYOUT(LEN1, LEN2, DES_BC_W_s(BCV), BCV_I, BCV)
         ENDIF
! End check mass inlet on XY face: 3D 

! Check mass inlet on YZ face: 3D 
! ----------------------------------------
! from previous checks in des_check_mio_location the combination of xe==xw,
! yn==ys, ys!=zero and ys!=ylength implies two items:
!   that xw is either zero or xlength, since if xw was zero or xlength then 
!     ys could not be zero or ylength (or the boundary would be on the
!     domain edge which is prohibited)
!   that zt!=zb, since zt==zb would result in a prohibited point boundary
! similarly the combination of xe==xw, zt==zb, zb!=zero and zb!=zlength
! implies the following:
!   that xw is either zero or xlength, since if xw was zero or xlength then 
!     zb could not be zero or zlength (or the boundary would be on the
!     domain edge which is prohibited)
!   that yn!=ys, since yn==ys would result in a prohibited point boundary
         IF(DES_BC_X_w(BCV) == DES_BC_X_e(BCV) .AND. &
            DES_BC_Y_s(BCV) /= DES_BC_Y_n(BCV) .AND. &
            DES_BC_Z_b(BCV) /= DES_BC_Z_t(BCV) &
         .OR. &
            DES_BC_X_w(BCV) == DES_BC_X_e(BCV) .AND. &
            DES_BC_Y_s(BCV) == DES_BC_Y_n(BCV) .AND. &
            DES_BC_Y_s(BCV) /= ZERO .AND. &
            DES_BC_Y_s(BCV) /= YLENGTH &
         .OR. &
            DES_BC_X_w(BCV) == DES_BC_X_e(BCV) .AND. &
            DES_BC_Z_b(BCV) == DES_BC_Z_t(BCV) .AND. &
            DES_BC_Z_b(BCV) /= ZERO .AND. &
            DES_BC_Z_b(BCV) /= ZLENGTH)THEN

! Verify that a velocity has been specifed for DES_BC_U_s(BCV)
            IF(DES_BC_U_s(BCV) == ZERO)THEN
               WRITE (UNIT_LOG, 1200)BCV,'U_s'; WRITE (*, 1200)BCV,'U_s'
               CALL MFIX_EXIT(myPE)
            ENDIF
            IF(DES_BC_X_w(BCV) == ZERO)THEN
! Correct velocity direction, if wrong
               IF(DES_BC_U_s(BCV) .LT. ZERO) THEN
                  DES_BC_U_s(BCV) = -DES_BC_U_s(BCV)
                  WRITE(*,1201) 'YZw'
               ENDIF                  
! Classify the boundary condition
               DES_MI_CLASS(BCV_I) = 'YZw'
            ENDIF
            IF(DES_BC_X_w(BCV) == XLENGTH)THEN
! Correct velocity direction, if wrong
               IF(DES_BC_U_s(BCV) .GT. ZERO) THEN
                  DES_BC_U_s(BCV) = -DES_BC_U_s(BCV)
                  WRITE(*,1201) 'YZe'
               ENDIF                  
! Classify the boundary condition
               DES_MI_CLASS(BCV_I) = 'YZe'
            ENDIF

! If a slit flow inlet is specified near an edge, expand the slit to 
! the width of one particle diameter :
! 1) check that the inlet slit +/- a particle radius is not near either
!    domain boundary edge (south/north or bottom/top) and if not shift
!    the south/bottom location of the inlet down by a particle radius 
!    and the north/top location of the inlet up by a particle radius
! 2) if the inlet slit is near the south/bottom domain edge then move
!    the south/bottom location of the inlet to the domain boundary edge 
!    and the north/top location one particle diameter distance away
! 3) if the inlet slit is near the north/top domain edge then move the
!    north/top location of the inlet to the domain boundary edge and the
!    south/bottom location one particle diameter distance away
            IF(DES_BC_Y_s(BCV) == DES_BC_Y_n(BCV))THEN
               IF(((DES_BC_Y_s(BCV) - D_P0(1)*HALF) .GE. ZERO) .AND. &
                  ((DES_BC_Y_s(BCV) + D_P0(1)*HALF) .LE. YLENGTH))THEN
                  DES_BC_Y_s(BCV) = DES_BC_Y_s(BCV) - D_P0(1)*HALF
                  DES_BC_Y_n(BCV) = DES_BC_Y_n(BCV) + D_P0(1)*HALF
               ELSEIF((DES_BC_Y_s(BCV) - D_P0(1)*HALF) .LT. ZERO)THEN
                  DES_BC_Y_s(BCV) = ZERO
                  DES_BC_Y_n(BCV) = D_P0(1) 
               ELSEIF((DES_BC_Y_n(BCV) + D_P0(1)*HALF) .GT. YLENGTH)THEN
                  DES_BC_Y_s(BCV) = YLENGTH - D_P0(1)
                  DES_BC_Y_n(BCV) = YLENGTH 
               ENDIF
            ENDIF
            IF(DES_BC_Z_b(BCV) == DES_BC_Z_t(BCV))THEN
               IF(((DES_BC_Z_b(BCV) - D_P0(1)*HALF) .GE. ZERO) .AND. &
                  ((DES_BC_Z_b(BCV) + D_P0(1)*HALF) .LE. ZLENGTH))THEN
                  DES_BC_Z_b(BCV) = DES_BC_Z_b(BCV) - D_P0(1)*HALF
                  DES_BC_Z_t(BCV) = DES_BC_Z_t(BCV) + D_P0(1)*HALF
               ELSEIF((DES_BC_Z_b(BCV) - D_P0(1)*HALF) .LT. ZERO)THEN
                  DES_BC_Z_b(BCV) = ZERO
                  DES_BC_Z_t(BCV) = D_P0(1)
               ELSEIF((DES_BC_Z_t(BCV) + D_P0(1)*HALF) .GT. ZLENGTH)THEN
                  DES_BC_Z_b(BCV) = ZLENGTH - D_P0(1)
                  DES_BC_Z_t(BCV) = ZLENGTH 
               ENDIF
            ENDIF
! This subroutine determines the pattern that the particles will need to
! enter the system, if any.
            LEN1 = (DES_BC_Y_n(BCV) - DES_BC_Y_s(BCV))
            LEN2 = (DES_BC_Z_t(BCV) - DES_BC_Z_b(BCV))
            CALL DES_MI_LAYOUT(LEN1, LEN2, DES_BC_U_s(BCV), BCV_I, BCV)
         ENDIF
! End check mass inlet on YZ face: 3D  

      ENDIF   ! endif dimn == 2

 1200 FORMAT(/1X,70('*')//&
         ' From: DES_MI_CLASSIFY -',/&
         ' Message: Boundary condition ',I4,' requires DES_BC_',A,&
         ' be nonzero.',/1X,70('*')/)
 1201 FORMAT(/1X,70('*')//&
         ' From: DES_MI_CLASSIFY -',/&
         ' Message: Direction of boundary condition velocity at ',A,&
         '-face was reversed.',/1X,70('*')/)

      RETURN
      END SUBROUTINE DES_MI_CLASSIFY



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: DES_MI_LAYOUT                                          !
!                                                                      !
!  Purpose:  This routine determines the layout of the mass inlet as   !
!  either ordered or random based upon the inlet conditions.  This     !
!  routine also verifies that the specified inlet conditions for mass  !
!  or volumetric flow rates along with inlet size (length or area) and !
!  particle inlet velocity will work.  If not an error is flagged and  !
!  the program is exited.                                              !
!                                                                      !
!                                                                      !
!  Author: J.Musser                                   Date: 14-Aug-09  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE DES_MI_LAYOUT(LEN1, LEN2, BC_VEL, BCV_I, BCV)

      USE compar
      USE constant
      USE des_bc
      USE discretelement 
      USE funits  
      USE geometry
      USE param1
      USE physprop
      USE run

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER TMP_FACTOR
      INTEGER LL, LC, I, J, IJ
! max number of particle diameters that fit along length of inlet
      INTEGER TMP_LEN1, TMP_LEN2
! passed arguments giving index information of the boundary
      INTEGER BCV_I, BCV 
! passed arguments giving the length of the boundary edges
      DOUBLE PRECISION LEN1, LEN2
! a random number between 0-1
      DOUBLE PRECISION TMP_DP
! a random integer between 0 and tmp_factor-1      
      INTEGER TMP
! passued argument giving the assigned velocity at the boundary      
      DOUBLE PRECISION BC_VEL
! lower threshold velocities that will work given the flow boundary
! specifications:
! lower maximum inlet particle velocity (MAXIPV); above which particles
! enter fast enough not to require additional handling
! lowest minimum inlet particle velocity (MINIPV) below which particles
! do not enter fast enough to meet the massflow rate and an error is
! flagged; above minipv and below maxipv the particle inlet conditions
! are adequate but must be carefully controlled
      DOUBLE PRECISION MAXIPV, MINIPV 
! used to clarify log messages in the event of an error      
      CHARACTER*17 SMIC
!-----------------------------------------------

      IF(DES_BC_VOLFLOW_s(BCV)  /= UNDEFINED)SMIC = 'DES_BC_VOLFLOW_s'
      IF(DES_BC_MASSFLOW_s(BCV) /= UNDEFINED)SMIC = 'DES_BC_MASSFLOW_s'
      
      WRITE(*,'(7X,A)') '---------- START DES_MI_LAYOUT ---------->'

      IF(LEN2 == ZERO)THEN   ! 2D domain 
         TMP_LEN1 = FLOOR(real(LEN1/D_P0(1)))
         TMP_LEN2 = ZERO
! notes :
!   dtsolid*pi_factor(:)        = time elapsed between particle injections
!                               = des_mi_time(:)
!   d_p0/(dtsolid*pi_factor(:)) = approx velocity needed to move one particle
!                                 diameter at the specified mass flow rate
!   ceiling(tmp_len1/2)         = the minimum no. of particles that can
!                                 be arranged along the inlet so that an 
!                                 additional particle cannot fit

         MAXIPV = D_P0(1)/( DTSOLID*dble(PI_FACTOR(BCV_I))*&
                  dble( CEILING(real(TMP_LEN1)/2.0)) ) 
         MINIPV = D_P0(1)/( DTSOLID*dble(PI_FACTOR(BCV_I))*&
                  dble(TMP_LEN1) )
         IF (MINIPV .LT. SMALL_NUMBER) MINIPV = ZERO

! The bc velocity is lowered by a small number to prevent possible
! issues with the mantissa truncation
         IF(MAXIPV .LE. (ABS(BC_VEL) - SMALL_NUMBER))THEN
! The inlet velocity is sufficient to permit random placement of the new
! particles without risk of overlap
            PARTICLE_PLCMNT(BCV_I) = 'RAND'
         ELSEIF(MINIPV .LE. ABS(BC_VEL) - SMALL_NUMBER .AND. &
         ABS(BC_VEL) .LT. MAXIPV + SMALL_NUMBER)THEN
! Then inlet velocity will require that the new particles be placed with
! order to prevent overlap.
            PARTICLE_PLCMNT(BCV_I) = 'ORDR'
         ELSE
! The inlet velocity is too low for the other inlet conditions, flag
! error and exit.
            WRITE(UNIT_LOG,1302) BCV, SMIC
            WRITE(*,1302) BCV, SMIC
            WRITE(*,1303) (MINIPV+SMALL_NUMBER)
            CALL MFIX_EXIT(myPE)
         ENDIF

         IF(PARTICLE_PLCMNT(BCV_I) == 'ORDR')THEN
! In 2D calculate the approx. particle line density (e.g., the no. of 
! particles along the inlet) needed to achieve the specified particle 
! mass flow rate and particle velocity;
            TMP_FACTOR = CEILING(real(D_P0(1) / &
               (dble(PI_FACTOR(BCV_I)) * DTSOLID * ABS(BC_VEL))))
            ALLOCATE( MI_ORDER(BCV_I)%VALUE( TMP_FACTOR ) )
            ALLOCATE( I_OF_MI(BCV_I)%VALUE( TMP_FACTOR ) )

! Initialize            
            MI_ORDER(BCV_I)%VALUE(:) = -1
            MI_FACTOR(BCV_I) = 1
! Dimension of grid cell; this may be larger than than the particle
! diameter but not smaller              
            MI_WINDOW(BCV_I) = LEN1/dble(TMP_FACTOR)            
            DO I = 1, TMP_FACTOR 
               I_OF_MI(BCV_I)%VALUE(I) = I - 1
            ENDDO

! Construct an array of integers from 1 to TMP_FACTOR in a random
! order. This is used when placing new particles. 
            LL = 1
            DO WHILE (MI_ORDER(BCV_I)%VALUE(TMP_FACTOR) .EQ. -1)
               CALL RANDOM_NUMBER(TMP_DP)
               TMP = CEILING(real(TMP_DP*dble(TMP_FACTOR)))
               DO LC = 1, LL
                 IF(TMP .EQ. MI_ORDER(BCV_I)%VALUE(LC) )EXIT
                 IF(LC .EQ. LL)THEN
                    MI_ORDER(BCV_I)%VALUE(LC) = TMP
                    LL = LL + 1
                 ENDIF
               ENDDO
            ENDDO           
         ENDIF     ! endif particle_plcmnt(bcv_i) == 'ordr'


      ELSEIF(LEN2 /= ZERO) THEN   ! 3D domain
         TMP_LEN1 = FLOOR(real(LEN1/D_P0(1)))
         TMP_LEN2 = FLOOR(real(LEN2/D_P0(1)))

! In the 3D case the calculation for MAXIPV is conservative.  That is,
! the actual bc velocity could be somewhat lower than the calculated 
! value of MAXIPV and still allow for random particle placement
         MAXIPV = D_P0(1)/( DTSOLID*dble(PI_FACTOR(BCV_I)) * &
                     dble( CEILING(real(TMP_LEN1*TMP_LEN2)/2.0)) ) 
! The cutoff is associated with square packing of disks on a plane
! A lower velocity would be possible with hexagonal packing
         MINIPV = D_P0(1)/( DTSOLID*dble(PI_FACTOR(BCV_I)) * &
                     dble(TMP_LEN1*TMP_LEN2) )

         IF (MINIPV .LT. SMALL_NUMBER) MINIPV = ZERO

         IF(MAXIPV .LE. (ABS(BC_VEL) - SMALL_NUMBER))THEN
! The inlet velocity is sufficient to permit random placement of the new
! particles without risk of overlap
            PARTICLE_PLCMNT(BCV_I) = 'RAND'
         ELSE IF(MINIPV .LE. ABS(BC_VEL) - SMALL_NUMBER .AND. &
                 ABS(BC_VEL) .LT. MAXIPV + SMALL_NUMBER)THEN
! Then inlet velocity will require that the new particles be placed with
! order to prevent overlap.
            PARTICLE_PLCMNT(BCV_I) = 'ORDR'
         ELSE
! The inlet velocity is too low for the other inlet conditions, flag
! error and exit.
            WRITE(UNIT_LOG,1302) BCV, SMIC;
            WRITE(*,1302) BCV, SMIC
            WRITE(*,1303) (MINIPV+SMALL_NUMBER)
            CALL MFIX_EXIT(myPE)
         ENDIF

         IF(PARTICLE_PLCMNT(BCV_I) == 'ORDR')THEN
! In 3D calculate the the number of grid cells within the inlet face
! with length dimension of 1 particle diameter
            TMP_FACTOR = TMP_LEN1 * TMP_LEN2
            ALLOCATE( MI_ORDER(BCV_I)%VALUE( TMP_FACTOR ) )
            ALLOCATE( I_OF_MI(BCV_I)%VALUE( TMP_FACTOR ) )
            ALLOCATE( J_OF_MI(BCV_I)%VALUE( TMP_FACTOR ) )

! Initialize            
            MI_ORDER(BCV_I)%VALUE(:) = -1
            MI_FACTOR(BCV_I) = 1
! Dimension of grid cell; this may be larger than than the particle
! diameter but not smaller: 
! if len1/tmp_len1 < len2/tmp_len2, len1/tmp_len1*tmp_len2 < len2 or
! if len1/tmp_len1 > len2/tmp_len2, len2/tmp_len2*tmp_len1 < len1
            MI_WINDOW(BCV_I) = MIN(LEN1/TMP_LEN1, LEN2/TMP_LEN2)

            DO I = 1, TMP_LEN1
               DO J = 1, TMP_LEN2
                  IJ = J + (I-1)*TMP_LEN2
                  I_OF_MI(BCV_I)%VALUE(IJ) = I - 1
                  J_OF_MI(BCV_I)%VALUE(IJ) = J - 1
               ENDDO
            ENDDO

! Construct an array of integers from 1 to TMP_FACTOR in a random
! order. This is used when placing new particles.
            LL = 1
            DO WHILE (MI_ORDER(BCV_I)%VALUE(TMP_FACTOR) .EQ. -1)
               CALL RANDOM_NUMBER(TMP_DP)
               TMP = CEILING(real(TMP_DP*dble(TMP_FACTOR)))
               DO LC = 1, LL
                 IF(TMP .EQ. MI_ORDER(BCV_I)%VALUE(LC) )EXIT
                 IF(LC .EQ. LL)THEN
                    MI_ORDER(BCV_I)%VALUE(LC) = TMP
                    LL = LL + 1
                 ENDIF
               ENDDO
            ENDDO
         ENDIF   ! endif particle_plcmnt(bcv_i) == 'ordr'

      ENDIF   ! endif len2 == zero


      WRITE(*,1304) BCV, LEN1, LEN2, TMP_LEN1, TMP_LEN2,&
         MAXIPV, MINIPV, PARTICLE_PLCMNT(BCV_I)

      WRITE(*,'(7X,A)') '<---------- END DES_MI_LAYOUT ----------'

 1302 FORMAT(/1X,70('*')//&
         ' From: DES_MI_LAYOUT -',/&
         ' Message: Inlet velocity for BC ',I3,&
         ' is too low '/3X,'for the specified ',A,' flow rate.')
 1303 FORMAT(3X,&
         'The particle inlet velocity should have a ',/,3X,&
         'magnitude greater than ',ES17.5,/1X,70('*')/)

 1304 FORMAT(9X,'For mass inlet BC: ', I3,&
            /12X,'LEN1 = ', ES17.5, ' LEN2 = ', ES17.5,&
            ' TMP_LEN1 = ', I6, ' TMP_LEN2 = ', I6,&
            /12X,'MAXIPV = ', ES17.5,' MINIPV = ', ES17.5,&
            ' and PLCMNT = ',A)

      RETURN
      END SUBROUTINE DES_MI_LAYOUT



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: DES_MI_CELLS                                           !
!                                                                      !
!  Purpose:                                                            !
!    Locate the i, j, k position of the cells near the DES_MI so that 
!    when particles are injected into the system, a complete grid 
!    search is not necessary to prevent/determine:
!      1) the injected particle from overlapping with an existing
!         particle nearby when particle_plcmnt = 'rand'      
!      2) the i, j, k locations of the newly injected particle      
!                                                                      !
!  Author: J.Musser                                   Date:  5-Oct-09  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE DES_MI_CELLS(BCV, BCV_I)

      USE des_bc
      USE discretelement
      USE geometry

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! passed arguments giving index information of the boundary
      INTEGER BCV_I, BCV

      INTEGER I,J,K
! dummy variable for length calculations
      DOUBLE PRECISION LOCATION
!-----------------------------------------------

     
      IF(DES_BC_X_w(BCV) == DES_BC_X_e(BCV))THEN
         IF(DES_BC_X_w(BCV) == ZERO) GS_ARRAY(BCV_I,1:2) = 1
         IF(DES_BC_X_w(BCV) == XLENGTH) GS_ARRAY(BCV_I,1:2) = IMAX2
      ELSE
         I=2; LOCATION = ZERO
         DO WHILE (LOCATION .LE. XLENGTH)
            IF((DES_BC_X_w(BCV) .GT. LOCATION) .AND. &
            (DES_BC_X_w(BCV) .LE. LOCATION + DX(I))) THEN
               GS_ARRAY(BCV_I,1) = I
            ENDIF
            IF((DES_BC_X_e(BCV) .GE. LOCATION) .AND. &
            (DES_BC_X_e(BCV) .LT. LOCATION + DX(I))) THEN
               GS_ARRAY(BCV_I,2) = I
               EXIT
            ENDIF
            LOCATION = LOCATION + DX(I)
            I=I+1
         ENDDO
      ENDIF

      IF(DES_BC_Y_s(BCV) == DES_BC_Y_n(BCV))THEN
         IF(DES_BC_Y_s(BCV) == ZERO) GS_ARRAY(BCV_I,3:4) = 1
         IF(DES_BC_Y_s(BCV) == YLENGTH) GS_ARRAY(BCV_I,3:4) = JMAX2
      ELSE
         J=2; LOCATION = ZERO
         DO WHILE (LOCATION .LE. YLENGTH)
            IF((DES_BC_Y_s(BCV) .GT. LOCATION) .AND. &
            (DES_BC_Y_s(BCV) .LE. LOCATION + DY(J))) THEN
               GS_ARRAY(BCV_I,3) = J
            ENDIF
            IF((DES_BC_Y_n(BCV) .GE. LOCATION) .AND. &
            (DES_BC_Y_n(BCV) .LT. LOCATION + DY(J))) THEN
               GS_ARRAY(BCV_I,4) = J            
               EXIT
            ENDIF
            LOCATION = LOCATION + DY(J)
            J=J+1
         ENDDO
      ENDIF

      IF(DIMN == 2) THEN
         GS_ARRAY(BCV_I,5:6) = 1
      ELSE
         IF(DES_BC_Z_b(BCV) == DES_BC_Z_t(BCV))THEN
            IF(DES_BC_Z_b(BCV) == ZERO) GS_ARRAY(BCV_I,5:6) = 1
            IF(DES_BC_Z_b(BCV) == ZLENGTH) GS_ARRAY(BCV_I,5:6) = KMAX2
         ELSE
            K=2; LOCATION = ZERO
            DO WHILE (LOCATION .LE. ZLENGTH)
               IF((DES_BC_Z_b(BCV) .GT. LOCATION) .AND. &
               (DES_BC_Z_b(BCV) .LE. LOCATION + DZ(K))) THEN
                  GS_ARRAY(BCV_I,5) = K
               ENDIF
               IF((DES_BC_Z_t(BCV) .GE. LOCATION) .AND. &
               (DES_BC_Z_t(BCV) .LT. LOCATION + DZ(K))) THEN
                  GS_ARRAY(BCV_I,6) = K
                  EXIT
               ENDIF
               LOCATION = LOCATION + DZ(K)
               K=K+1
            ENDDO
         ENDIF
      ENDIF

      RETURN
      END SUBROUTINE DES_MI_CELLS



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: DES_MO_CLASSIFY                                        !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!  Author: J.Musser                                   Date:  5-Oct-09  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE DES_MO_CLASSIFY(BCV_I, BCV)

      USE compar
      USE constant
      USE des_bc
      USE discretelement 
      USE funits  
      USE geometry
      USE param1
      USE physprop
      USE run

      IMPLICIT NONE

!-----------------------------------------------
! Local variables
!-----------------------------------------------
! passed arguments giving index information of the boundary
      INTEGER BCV_I, BCV
!-----------------------------------------------


      IF(DIMN == 2)THEN   ! 2D domain

! Check verticle mass outlet: 2D
! ----------------------------------------
         IF(DES_BC_X_w(BCV) == DES_BC_X_e(BCV))THEN
            DES_MO_X = .TRUE.
            DES_MO_CLASS(BCV_I) = 'Xwe'   ! 'XW' or 'XE'
         ENDIF

! Check horizontal mass outlet: 2D 
! ----------------------------------------
         IF(DES_BC_Y_s(BCV) == DES_BC_Y_n(BCV))THEN
            DES_MO_Y = .TRUE.
            DES_MO_CLASS(BCV_I) = 'Ysn'   ! 'YN' or 'YS'
         ENDIF

     ELSE   !  3D domain

! Check mass outlet on YZ face: 3D 
! ----------------------------------------
! see comments following the yz face section of des_mi_classify
! for elucidation on this, and subsequent series of if statements
         IF(DES_BC_X_w(BCV) == DES_BC_X_e(BCV) .AND. &
            DES_BC_Y_s(BCV) /= DES_BC_Y_n(BCV) .AND. &
            DES_BC_Z_b(BCV) /= DES_BC_Z_t(BCV) &
         .OR. &
            DES_BC_X_w(BCV) == DES_BC_X_e(BCV) .AND. &
            DES_BC_Y_s(BCV) == DES_BC_Y_n(BCV) .AND. &
            DES_BC_Y_s(BCV) /= ZERO .AND. &
            DES_BC_Y_s(BCV) /= YLENGTH &
         .OR. &
            DES_BC_X_w(BCV) == DES_BC_X_e(BCV) .AND. &
            DES_BC_Z_b(BCV) == DES_BC_Z_t(BCV) .AND. &
            DES_BC_Z_b(BCV) /= ZERO .AND. &
            DES_BC_Z_b(BCV) /= ZLENGTH)THEN

            DES_MO_X = .TRUE.
            DES_MO_CLASS(BCV_I) = 'Xwe'  ! 'YZw' or 'YZe'

         ENDIF
! End check mass outlet on YZ face: 3D 


! Check mass inlet on XZ face: 3D 
! ----------------------------------------
         IF(DES_BC_Y_s(BCV) == DES_BC_Y_n(BCV) .AND. &
            DES_BC_X_w(BCV) /= DES_BC_X_e(BCV) .AND. &
            DES_BC_Z_b(BCV) /= DES_BC_Z_t(BCV) &
         .OR. &
            DES_BC_Y_s(BCV) == DES_BC_Y_n(BCV) .AND. &
            DES_BC_X_w(BCV) == DES_BC_X_e(BCV) .AND. &
            DES_BC_X_w(BCV) /= ZERO .AND. &
            DES_BC_X_w(BCV) /= XLENGTH &
         .OR. &
            DES_BC_Y_s(BCV) == DES_BC_Y_n(BCV) .AND. &
            DES_BC_Z_b(BCV) == DES_BC_Z_t(BCV) .AND. &
            DES_BC_Z_b(BCV) /= ZERO .AND. &
            DES_BC_Z_b(BCV) /= ZLENGTH)THEN

            DES_MO_Y = .TRUE.
            DES_MO_CLASS(BCV_I) = 'Ysn'   ! 'XZs' or 'XZn'

         ENDIF
! End check mass outlet on XZ face: 3D 

! Check mass inlet on XY face: 3D 
! ----------------------------------------
         IF(DES_BC_Z_b(BCV) == DES_BC_Z_t(BCV) .AND. &
            DES_BC_X_w(BCV) /= DES_BC_X_e(BCV) .AND. &
            DES_BC_Y_s(BCV) /= DES_BC_Y_n(BCV) &
         .OR. &
            DES_BC_Z_b(BCV) == DES_BC_Z_t(BCV) .AND. &
            DES_BC_X_w(BCV) == DES_BC_X_e(BCV) .AND. &
            DES_BC_X_w(BCV) /= ZERO .AND. &
            DES_BC_X_w(BCV) /= XLENGTH &
         .OR. &
            DES_BC_Z_b(BCV) == DES_BC_Z_t(BCV) .AND. &
            DES_BC_Y_s(BCV) == DES_BC_Y_n(BCV) .AND. &
            DES_BC_Y_s(BCV) /= ZERO .AND. &
            DES_BC_Y_s(BCV) /= YLENGTH)THEN

            DES_MO_Z = .TRUE.
            DES_MO_CLASS(BCV_I) = 'Zbt'   ! 'XYt' or 'XYb'

         ENDIF
! End check mass outlet on XY face: 3D 

      ENDIF   ! endif dimn == 2

      RETURN
      END SUBROUTINE DES_MO_CLASSIFY

