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
      USE mfix_pic
      IMPLICIT NONE

!-----------------------------------------------
! Local parameters
!-----------------------------------------------
      INTEGER, PARAMETER :: DIM_BCTYPE = 4

!-----------------------------------------------
! Local variables
!-----------------------------------------------

! Loop counters
      INTEGER BCV, I
! Solids phase index
      INTEGER M
! tmp variable to calculate solids volume fraction at inlet   
      DOUBLE PRECISION EPs_tmp
! valid boundary condition types
      CHARACTER*16, DIMENSION(1:DIM_BCTYPE) ::VALID_BC_TYPE = (/&
           'MASS_INFLOW     ', 'MI              ',&
           'MASS_OUTFLOW    ', 'MO              '/)

!-----------------------------------------------
!   External functions
!-----------------------------------------------
      LOGICAL , EXTERNAL :: COMPARE 

!-----------------------------------------------           

! Initialize
      DES_BCMI = 0; DES_MI = .FALSE.
      DES_BCMO = 0
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
               IF(DMP_LOG) WRITE (UNIT_LOG, 1000) BCV
               WRITE (*, 1000) BCV
               CALL MFIX_EXIT(myPE)
            ENDIF 
            IF(DIMN == 3)THEN
              IF(DES_BC_Z_b(BCV) == UNDEFINED .OR. &
                 DES_BC_Z_t(BCV) == UNDEFINED) THEN
               IF(DMP_LOG) WRITE (UNIT_LOG, 1000) BCV
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
               IF(DMP_LOG) WRITE (UNIT_LOG, 1001) BCV
               WRITE (*, 1001) BCV
               CALL MFIX_EXIT(myPE)
            ENDIF
            IF(DIMN == 3)THEN
               IF(DES_BC_Z_b(BCV) .LT. 0 .OR. &
                  DES_BC_Z_t(BCV) .GT. ZLENGTH .OR. &
                  DES_BC_Z_b(BCV) .GT. DES_BC_Z_t(BCV))THEN
                  IF(DMP_LOG) WRITE (UNIT_LOG, 1001) BCV
                  WRITE (*, 1001)BCV
                  CALL MFIX_EXIT(myPE)
               ENDIF
            ENDIF

! Require the inlet/outlet to be on a boundary/wall of the system
            CALL DES_CHECK_MIO_LOCATION(BCV)

            DO I = 1, DIM_BCTYPE 
               VALID_BC_TYPE(I) = TRIM(VALID_BC_TYPE(I))
               IF (VALID_BC_TYPE(I) == DES_BC_TYPE(BCV)) THEN 
! If des_bc_type set with abbreviation form, then change to corresponding full form
                  IF (MOD(I,2) == 0) DES_BC_TYPE(BCV) = VALID_BC_TYPE(I-1)

! Check if solids phase velocity at the BC plane is specified, if not,
! set solids velocity to zero and flag log file
                  IF(DES_BC_TYPE(BCV) == 'MASS_INFLOW')THEN
                     DES_MI = .TRUE.
                     DES_BCMI = DES_BCMI + 1

                     EPs_tmp = ZERO
                     DO M=1, MMAX
! If DES_BC_ROP_s is well defined (defined and not zero), check that either
! DES_BC_MASSFLOW or DES_BC_VOLFLOW is also well defined.  If not, flag an
! error and exit.
                        IF(DES_BC_ROP_s(BCV,M) /= UNDEFINED .AND. &
                          .NOT. COMPARE(DES_BC_ROP_s(BCV,M),ZERO) ) THEN

! Check that density information is not missing, if so exit.  RO_s is 
! checked for unrealistic values in check_data_04.f                              
                           IF(COMPARE(RO_s(M),ZERO))THEN
                              IF(DMP_LOG) WRITE(UNIT_LOG,1006) M, BCV
                              WRITE(*,1006) M, BCV
                              CALL MFIX_EXIT(myPE)
                           ENDIF

! Check the values given for MASSFLOW and/or VOLFLOW
                           IF(DES_BC_MASSFLOW_s(BCV,M) /= UNDEFINED) THEN
                              IF(DES_BC_VOLFLOW_s(BCV,M) /=UNDEFINED) THEN
! Both volumetric and mass flow rates have been defined. 
! Verify that the values match.
                                 IF(.NOT.COMPARE(DES_BC_VOLFLOW_s(BCV,M),&
                                   DES_BC_MASSFLOW_s(BCV,M)/RO_s(M))) THEN
                                   IF(DMP_LOG) WRITE(UNIT_LOG,1006)BCV,M
                                    WRITE(*,1006)BCV,M
                                    CALL MFIX_EXIT(myPE)
                                 ENDIF
                              ELSE
! Only MASSFLOW was given. Calculate VOLFLOW.
                                 DES_BC_VOLFLOW_s(BCV,M) = &
                                    DES_BC_MASSFLOW_s(BCV,M)/RO_s(M)
                              ENDIF
                           ELSE   ! no mass flow rate is specified
                              IF(DES_BC_VOLFLOW_s(BCV,M) == UNDEFINED) THEN
! If neither a volumetric or mass flow rate is specified, exit                                      
                                 IF(DMP_LOG) WRITE(UNIT_LOG,1011)BCV,M
                                 WRITE(*,1011)BCV,M
                                 CALL MFIX_EXIT(myPE)
                              ENDIF
                           ENDIF   
! Add solids volume to total solds volume to check for overflow
                           EPs_tmp = EPs_tmp + &
                              (DES_BC_ROP_s(BCV,M)/RO_s(M))
                        ENDIF  ! endif des_bc_rop_s is well defined

! Back check that if either DES_BC_MASSFLOW_s or DES_BC_VOLFLOW_s are
! well defined (defined and not zero), then ROP_s is also well defined.
                        IF((DES_BC_MASSFLOW_s(BCV,M) /= UNDEFINED .AND. &
                           .NOT.COMPARE(DES_BC_MASSFLOW_s(BCV,M),ZERO)) &
                        .OR. &
                           (DES_BC_VOLFLOW_s(BCV,M) /= UNDEFINED .AND. &
                            .NOT.COMPARE(DES_BC_VOLFLOW_s(BCV,M),ZERO))) THEN

                           IF(DES_BC_ROP_s(BCV,M) == UNDEFINED .OR. &
                              COMPARE(DES_BC_ROP_s(BCV,M),ZERO))THEN
! A nonzero mass or volumetric flow rate is defined for BCV on mass
! phase M, and either ROP_s is zero or undefined.  Flag error and exit.
                              IF(DMP_LOG) WRITE(UNIT_LOG,1015)BCV,M
                              WRITE(*,1015)BCV, M
                              CALL MFIX_EXIT(myPE)
                           ENDIF
                        ENDIF
                     ENDDO   ! end loop of M=1,MMAX

                     IF(COMPARE(EPs_tmp,ZERO))THEN
! A des inlet has been defined, but there is no specifed flow information
! (no solids are entering)                             
                        IF(DMP_LOG) WRITE(UNIT_LOG,1013)BCV
                        WRITE(*,1013)BCV
                        CALL MFIX_EXIT(myPE)
                     ELSEIF(EPs_tmp > ONE) THEN
! The total solids volume exceeds one
                        IF(DMP_LOG) WRITE(UNIT_LOG,1014)BCV
                        WRITE(*,1014)BCV
                        CALL MFIX_EXIT(myPE)
                     ENDIF


                  ELSE   ! if des_bc_type is not mass_inflow
                     DES_BCMO = DES_BCMO + 1
                  ENDIF   ! end if des_bc_type is MI, etc

! gives way to exit loop over des_bc_defined(bcv) when valid_bc_type is T                 
                  CYCLE CHECK_BC

               ENDIF   ! endif des_bc_type is valid
            ENDDO   ! end loop over dim_bc_type

! exit if des_bc_defined but bc_type is not valid            
            IF(DMP_LOG) WRITE (UNIT_LOG, 1002) BCV, DES_BC_TYPE(BCV)
            WRITE (*, 1002) BCV, DES_BC_TYPE(BCV) 
            IF(DMP_LOG) WRITE (UNIT_LOG, 1003) VALID_BC_TYPE
            WRITE (*, 1003) VALID_BC_TYPE
            WRITE (*, 1008)
            CALL MFIX_EXIT(myPE)  

         ENDIF   ! end if des_bc_defined(bcv)

      ENDDO CHECK_BC

      WRITE(*,1012) DES_BCMI, DES_BCMO

      IF(DES_BCMI /= 0 .OR. DES_BCMO /=0)THEN
         DES_MIO = .TRUE.

! Verify that either the nsquare or grid based neighbor searches are
! used, otherwise flag and exit
         IF((DES_NEIGHBOR_SEARCH == 2) .OR. &
            (DES_NEIGHBOR_SEARCH == 3)) THEN
            IF(DMP_LOG) WRITE (UNIT_LOG, 1005)
            WRITE (*, 1005)
            CALL MFIX_EXIT(myPE)
         ENDIF
      ENDIF



! If the particle count is not defined, but MAX_PIS, the maximum number
! of particles permitted in the system at any given time, is set, assume
! the system is starting empty then check whether a mass inlet has been
! specified.

! The variable PARTICLES should already be set by this point if using
! gener_part_config option      
      IF(PARTICLES == UNDEFINED_I .AND. MAX_PIS /= UNDEFINED_I)THEN
         PARTICLES = 0
      ELSEIF(PARTICLES == UNDEFINED_I .AND. MAX_PIS == UNDEFINED_I)THEN
         IF(DMP_LOG) WRITE(UNIT_LOG,'(3X,A)')&
            'Either PARTICLES or MAX_PIS must specified in mfix.dat'
         CALL MFIX_EXIT(myPE)
      ELSEIF(PARTICLES == 0 .AND. MAX_PIS == UNDEFINED_I) THEN
         IF(DMP_LOG) WRITE(UNIT_LOG,'(3X,A,/3X,A,A)') &
            'If starting with 0 PARTICLES, MAX_PIS must be ', &
            'specified in mfix.dat'
         CALL MFIX_EXIT(myPE)         
      ENDIF 

      IF (.NOT.GENER_PART_CONFIG) THEN
         IF(DMP_LOG) WRITE(UNIT_LOG,'(3X,A,I10)') &
            'Total number of particles read from input file = ', PARTICLES
         IF(DMP_LOG) WRITE(*,'(3X,A,I10)') &
            'Total number of particles read from input file = ', PARTICLES
      ENDIF

! If the system is started without any particles and an inlet is not
! specified, the run is likely aborted.
! Inlet/outlet for MPPIC are based off the regular mfix declarations, 
! and so des_bcmi could still be zero.
      IF (PARTICLES == 0) THEN
         IF(DES_BCMI == 0 .AND. (.NOT.MPPIC))THEN
            IF(DMP_LOG) WRITE(UNIT_LOG, 1009)
            WRITE(*, 1009)
           CALL MFIX_EXIT(myPE)
         ENDIF
         WRITE(*,'(5X,A)') &
            'Run initiated with no particles in the system'
      ENDIF

! Check MAX_PIS requirements
      IF(DES_BCMI == 0 .AND. MAX_PIS == UNDEFINED_I)THEN
         WRITE(*,'(5X,A)')'Setting MAX_PIS = PARTICLES'
         MAX_PIS = PARTICLES
      ELSEIF(DES_BCMI /= 0 .AND. MAX_PIS == UNDEFINED_I)THEN
         IF(DMP_LOG) WRITE(UNIT_LOG, 1010)
         WRITE(*, 1010)
         CALL MFIX_EXIT(myPE)
      ENDIF

      WRITE(*,'(3X,A)') '<---------- END CHECK_DES_BC ----------'



 1000 FORMAT(/1X,70('*')//, ' From: CHECK_DES_BC',/, ' Message: ',&
         'Insufficient DEM boundary condition infomation',/10X,&
         'Check boundary condition number: ',I3,/1X,70('*')/)

 1001 FORMAT(/1X,70('*')//, ' From: CHECK_DES_BC',/, ' Messsage: ',&
         'Improper DEM boundary condition information',/10X,&
         'Check boundary condition number: ',I3,/1X,70('*')/)

 1002 FORMAT(/1X,70('*')//, ' From: CHECK_DES_BC',/, ' Message: ',&
         'Illegal BC_TYPE for boundary condition ',I3,/10X,&
         'BC_TYPE = ',A,' and the valid types are: ') 
 
 1003 FORMAT(5X,A16)


 1005 FORMAT(/1X,70('*')//, ' From: CHECK_DES_BC',/, ' Message: ',&
         'Currently, DEM inlet/outlet only supports the NSQUARE',/10X,&
         'and GRID_BASED_NEIGHBOR search methods.',/1X,70('*')/)

 1006 FORMAT(/1X,70('*')//, ' From: CHECK_DES_BC',/, ' Message: ',&
         'Solids phase ',I3,' is used at DEM boundary condition ',I3,/10X,&
         'but its density RO_s is defined as zero.',/1X,70('*')/)

 1008 FORMAT(/1X,70('*')/) 

 1009 FORMAT(/1X,70('*')//, ' From: CHECK_DES_BC',/, ' Message: ',&
         'The system is initiated with no particles and no DEM',/10X,&
         'inlet. The run is being terminated.',/1X,70('*')/)

 1010 FORMAT(/1X,70('*')//, ' From: CHECK_DES_BC',/, ' Message: ',&
         'If a DEM inlet is specified then the maximum number of',/10X,&
         'particles permitted in the system (MAX_PIS) must be set',/10X,&
         'in mfix.dat',/1X,70('*')/)

 1011 FORMAT(/1X,70('*')//, ' From: CHECK_DES_BC',/, ' Message: ',&
         'DES_BC_ROP_s is defined (and not zero) for boundary ',/10X,&
         'condition ',I3, ' and solids phase ',I2,' without defining',/10X,&
         'either DES_BC_MASSFLOW_s or DES_BC_VOLFOW.'/,1X,70('*')/)

 1012 FORMAT(3X,'No. of mass inlet BC = ', I4,/,&
             3X,'No. of mass outlet BC = ', I4)

 1013 FORMAT(/1X,70('*')//, ' From: CHECK_DES_BC',/, ' Message: ',&
         'Boundary condition ',I3,' is identifed as an DEM inlet',/10X,&
         'but has no associated inflow information. Check mfix.dat',&
         /1X,70('*')/)
         
 1014 FORMAT(/1X,70('*')//, ' From: CHECK_DES_BC',/, ' Message: ',&
         'The total solids volume fraction for boundary condition',/10X,&
         I3,' exceeds one. Check the mfix.dat file.',/1X,70('*')/)

 1015 FORMAT(/1X,70('*')//, ' From: CHECK_DES_BC',/, ' Message: ',&
        'A nonzero mass or volumetric flow rate has been specified',&
        'for',/10X,'boundary condition ', I3, ' and solids phase ',I3,/10X,& 
        'and the associated DES_BC_ROP_s is either undefined or zero',/10X,&
        'Check mfix.dat',/,1X,70('*')/)


      RETURN
      END SUBROUTINE CHECK_DES_BC



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: DES_CHECK_MIO_LOCATION                                 !
!                                                                      !
!  Purpose: This subroutine verifies that the location of the DES mass !
!  inlet/outlet was specified to on a wall of the system.  If an inlet !
!  or outlet is located elsewhere, an error is flagged and the program !
!  terminated.                                                         !
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
               IF(DMP_LOG) WRITE (UNIT_LOG, 1100)BCV; WRITE (*, 1100)BCV
               CALL MFIX_EXIT(myPE)
            ENDIF
         ENDIF 
! Check horizontal mass inlet
         IF(DES_BC_Y_s(BCV) == DES_BC_Y_n(BCV)) THEN
            IF(DES_BC_Y_s(BCV) /= ZERO .AND. &
               DES_BC_Y_s(BCV) /= YLENGTH)THEN
               IF(DMP_LOG) WRITE (UNIT_LOG, 1100)BCV; WRITE (*, 1100)BCV
               CALL MFIX_EXIT(myPE)
            ENDIF
         ENDIF

! Require the inlet be an edge in 2D not a point
         IF(DES_BC_X_w(BCV) == DES_BC_X_e(BCV) .AND. &
            DES_BC_Y_s(BCV) == DES_BC_Y_n(BCV)) THEN
            IF(DMP_LOG) WRITE (UNIT_LOG, 1101)&
            'DES mass inlet/outlet must larger than a single point.',BCV
            WRITE (*, 1101)&
            'DES mass inlet/outlet must larger than a single point.',BCV
            CALL MFIX_EXIT(myPE)
         ENDIF
! Require the inlet be an edge in 2D not a plane
         IF(DES_BC_X_w(BCV) /= DES_BC_X_e(BCV) .AND. &
            DES_BC_Y_s(BCV) /= DES_BC_Y_n(BCV)) THEN
            IF(DMP_LOG) WRITE (UNIT_LOG, 1101)&
            'For DIMN=2, DES mass inlet/outlet cannot be an area.',BCV
            WRITE (*, 1101)&
            'For DIMN=2, DES mass inlet/outlet cannot be an area.',BCV
            CALL MFIX_EXIT(myPE)
         ENDIF

      ELSE   ! if dimn != 2

         IF(DES_BC_X_w(BCV) == DES_BC_X_e(BCV)) THEN
            IF(DES_BC_X_w(BCV) /= ZERO .AND. DES_BC_X_w(BCV) /= XLENGTH)THEN
               IF(DES_BC_Y_s(BCV) == DES_BC_Y_n(BCV)) THEN
!  Require the inlet be a slit/edge on the XZ-face
                  IF(DES_BC_Y_s(BCV) /= ZERO .AND. &
                     DES_BC_Y_s(BCV) /= YLENGTH)THEN
                     IF(DMP_LOG) WRITE (UNIT_LOG, 1100)BCV; WRITE (*, 1100)BCV
                     CALL MFIX_EXIT(myPE)
                  ENDIF
               ENDIF
               IF(DES_BC_Z_b(BCV) == DES_BC_Z_t(BCV)) THEN
! Require the inlet be a slit/edge on the XY-face
                  IF(DES_BC_Z_b(BCV) /= ZERO .AND. &
                     DES_BC_Z_b(BCV) /= ZLENGTH)THEN
                     IF(DMP_LOG) WRITE (UNIT_LOG, 1100)BCV; WRITE (*, 1100)BCV
                     CALL MFIX_EXIT(myPE)
                  ENDIF
               ENDIF
! Require an area inlet to be on the YZ-face and the face must be at
! one of the x-boundary edges
               IF(DES_BC_Y_s(BCV) /= DES_BC_Y_n(BCV) .AND. & 
                  DES_BC_Z_b(BCV) /= DES_BC_Z_t(BCV)) THEN
                  IF(DMP_LOG) WRITE (UNIT_LOG, 1100)BCV; WRITE (*, 1100)BCV
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
                     IF(DMP_LOG) WRITE (UNIT_LOG, 1100)BCV; WRITE (*, 1100)BCV
                     CALL MFIX_EXIT(myPE)
                  ENDIF
               ENDIF
               IF(DES_BC_Z_b(BCV) == DES_BC_Z_t(BCV)) THEN
! Require the inlet be a slit/edge on the XY-face
                  IF(DES_BC_Z_b(BCV) /= ZERO .AND. &
                     DES_BC_Z_b(BCV) /= ZLENGTH)THEN
                     IF(DMP_LOG) WRITE (UNIT_LOG, 1100)BCV; WRITE (*, 1100)BCV
                     CALL MFIX_EXIT(myPE)
                  ENDIF
               ENDIF
               IF(DES_BC_X_w(BCV) /= DES_BC_X_e(BCV) .AND. & 
                  DES_BC_Z_b(BCV) /= DES_BC_Z_t(BCV)) THEN
! Require an area inlet to be on the XZ-face and the face must be at
! one of the y-boundary edges
                  IF(DMP_LOG) WRITE (UNIT_LOG, 1100)BCV; WRITE (*, 1100)BCV
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
                     IF(DMP_LOG) WRITE (UNIT_LOG, 1100)BCV; WRITE (*, 1100)BCV
                     CALL MFIX_EXIT(myPE)
                  ENDIF
               ENDIF
               IF(DES_BC_Y_s(BCV) == DES_BC_Y_n(BCV)) THEN
! Require the inlet be a slit/edge on the XZ-face
                  IF(DES_BC_Y_s(BCV) /= ZERO .AND. &
                     DES_BC_Y_s(BCV) /= YLENGTH)THEN
                     IF(DMP_LOG) WRITE (UNIT_LOG, 1100)BCV; WRITE (*, 1100)BCV
                     CALL MFIX_EXIT(myPE)
                  ENDIF
                  IF(DES_BC_X_w(BCV) /= DES_BC_X_e(BCV) .AND. & 
                     DES_BC_Y_s(BCV) /= DES_BC_Y_n(BCV)) THEN
! Require an area inlet to be on the XY-face and the face must be at
! one of the z-boundary edges
                     IF(DMP_LOG) WRITE (UNIT_LOG, 1100)BCV; WRITE (*, 1100)BCV
                     CALL MFIX_EXIT(myPE)
                  ENDIF
               ENDIF
            ENDIF
         ENDIF

      ENDIF ! end if dimn == 2 / else


! The following checks target potential bc problems in 3D systems   
! ----------------------------------------
      IF(DES_BC_X_w(BCV) == DES_BC_X_e(BCV) .AND. &
         DES_BC_Y_s(BCV) == DES_BC_Y_n(BCV) .AND. &
         DES_BC_Z_b(BCV) == DES_BC_Z_t(BCV))THEN
! Prohibit point-based mass inlet
         IF(DMP_LOG) WRITE (UNIT_LOG, 1101)&
         'DES mass inlet/outlet must larger than a single point.',BCV
         WRITE (*, 1101)&
         'DES mass inlet/outlet must larger than a single point.',BCV
         CALL MFIX_EXIT(myPE)
      ENDIF
      IF(DES_BC_X_w(BCV) /= DES_BC_X_e(BCV) .AND. &
         DES_BC_Y_s(BCV) /= DES_BC_Y_n(BCV) .AND. &
         DES_BC_Z_b(BCV) /= DES_BC_Z_t(BCV))THEN
! Prohibit volume-based mass inlet
         IF(DMP_LOG) WRITE (UNIT_LOG, 1101)&
         'For DIMN=3, DES mass inlet/outlet cannot be a volume.',BCV
         WRITE (*, 1101)&
         'For DIMN=3, DES mass inlet/outlet cannot be a volume.',BCV
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
            IF(DMP_LOG) WRITE (UNIT_LOG, 1101)&
         'DES mass inlet/outlet cannot be positioned in a corner.',BCV
         WRITE (*, 1101)&
         'DES mass inlet/outlet cannot be positioned in a corner.',BCV
         CALL MFIX_EXIT(myPE)
      ENDIF

 1100 FORMAT(/1X,70('*')//, ' From: DES_CHECK_MIO_LOCATION',/10X,&
         ' Message: DES boundary condition ',I4, ' must be ',/10X,&
         ' specified on a wall.',/1X,70('*'))

 1101 FORMAT(/1X,70('*')//, ' From: DES_CHECK_MIO_LOCATION',/10X,&
         ' Message: ',A, /10X,'Check boundary condition ',I4,&
         ' in mfix.dat.',/1X,70('*'))

      RETURN
      END SUBROUTINE DES_CHECK_MIO_LOCATION








