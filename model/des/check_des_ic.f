!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: CHECK_DES_IC                                           !
!                                                                      !
!  Purpose: Check the data provided for des initial conditions.        !
!                                                                      !
!                                                                      !
!  Author: J.Musser                                   Date: 17-Feb-11  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_DES_IC

      USE des_ic
      USE discretelement 
      USE param
      USE param1
      USE des_thermo
      USE des_rxns
      USE compar
      USE constant
      USE funits  
      USE geometry
      USE indices
      USE physprop
      USE run

      IMPLICIT NONE


      INTEGER ICV

!-----------------------------------------------
!   External functions
!-----------------------------------------------
      LOGICAL , EXTERNAL :: COMPARE 

!-----------------------------------------------           

      IF(DMP_LOG) WRITE(*,'(1X,A)') &
         '---------- START CHECK_DES_IC ---------->'

! Initialize the initial condition logical.
      DES_IC_EXIST = .FALSE.

! Check for des inlet/outlet information:
      CHECK_IC: DO ICV = 1, DIMENSION_IC 
         DES_IC_DEFINED(ICV) = .FALSE. 
         IF (DES_IC_X_w(ICV) /= UNDEFINED) DES_IC_DEFINED(ICV) = .TRUE. 
         IF (DES_IC_X_e(ICV) /= UNDEFINED) DES_IC_DEFINED(ICV) = .TRUE. 
         IF (DES_IC_Y_s(ICV) /= UNDEFINED) DES_IC_DEFINED(ICV) = .TRUE. 
         IF (DES_IC_Y_n(ICV) /= UNDEFINED) DES_IC_DEFINED(ICV) = .TRUE. 
         IF (DES_IC_Z_b(ICV) /= UNDEFINED) DES_IC_DEFINED(ICV) = .TRUE. 
         IF (DES_IC_Z_t(ICV) /= UNDEFINED) DES_IC_DEFINED(ICV) = .TRUE. 

! If an initial condition is specified, verify necessary data
! Check that all dimensions have been defined
         IF (DES_IC_DEFINED(ICV)) THEN
            IF(DES_IC_X_w(ICV) == UNDEFINED .OR. &
               DES_IC_X_e(ICV) == UNDEFINED .OR. &
               DES_IC_Y_s(ICV) == UNDEFINED .OR. &
               DES_IC_Y_n(ICV) == UNDEFINED) THEN
               IF(DMP_LOG) THEN
                  WRITE (UNIT_LOG, 1000) ICV
                  WRITE (*, 1000) ICV
               ENDIF
               CALL MFIX_EXIT(myPE)
            ENDIF 
            IF(DIMN == 3)THEN
               IF(DES_IC_Z_b(ICV) == UNDEFINED .OR. &
                 DES_IC_Z_t(ICV) == UNDEFINED) THEN
               IF(DMP_LOG) THEN
                  WRITE (UNIT_LOG, 1000) ICV
                  WRITE (*, 1000) ICV
               ENDIF
               CALL MFIX_EXIT(myPE)
              ENDIF
            ENDIF
! Check that the region does not overlap
            IF(DES_IC_X_w(ICV) .LT. ZERO .OR. &
               DES_IC_Y_s(ICV) .LT. ZERO .OR. &
               DES_IC_X_e(ICV) .GT. XLENGTH .OR. &
               DES_IC_Y_n(ICV) .GT. YLENGTH .OR. &
               DES_IC_X_w(ICV) .GT. DES_IC_X_e(ICV) .OR. &
               DES_IC_Y_s(ICV) .GT. DES_IC_Y_n(ICV))THEN
               IF(DMP_LOG) THEN
                  WRITE (UNIT_LOG, 1001) ICV
                  WRITE (*, 1001) ICV
               ENDIF
               CALL MFIX_EXIT(myPE)
            ENDIF

            IF(DIMN == 3)THEN
               IF(DES_IC_Z_b(ICV) .LT. 0 .OR. &
                  DES_IC_Z_t(ICV) .GT. ZLENGTH .OR. &
                  DES_IC_Z_b(ICV) .GT. DES_IC_Z_t(ICV))THEN
                  IF(DMP_LOG) THEN
                     WRITE (UNIT_LOG, 1001) ICV
                     WRITE (*, 1001)ICV
                  ENDIF
                  CALL MFIX_EXIT(myPE)
               ENDIF
            ENDIF

! This is used to invoke the call to DES_SET_IC which is only done
! if there are valid initial conditions specified.
            DES_IC_EXIST = .TRUE.

         ENDIF   ! end if DES_IC_DEFINED(ICV)

      ENDDO CHECK_IC


! If the DES energy equations or DES species equation are being solve
! then initial conditions must be specified to provide the intial 
! sate of the particles.
      IF(.NOT.DES_IC_EXIST) THEN
         IF(DES_ENERGY_EQ) THEN
            IF(DMP_LOG) THEN
               WRITE(*,1002)'temperature','energy equation'
               WRITE(UNIT_LOG,1002)
            ENDIF
            CALL MFIX_EXIT(myPE)
         ELSEIF(ANY_DES_SPECIES_EQ) THEN
            IF(DMP_LOG) THEN
               WRITE(*,1002)'species mass fractions','species equation'
               WRITE(UNIT_LOG,1002)
            ENDIF
            CALL MFIX_EXIT(myPE)
         ENDIF
      ENDIF

      WRITE(*,'(1X,A)') '<---------- END CHECK_DES_IC ----------'

 1000 FORMAT(/1X,70('*')/, ' From: CHECK_DES_IC',/, ' Message: ',&
         'Insufficient DEM initial condition infomation',/10X,&
         'Check initial condition number: ',I3,/1X,70('*')/)

 1001 FORMAT(/1X,70('*')/, ' From: CHECK_DES_IC',/, ' Messsage: ',&
         'Improper DEM initial condition information',/10X,&
         'Check initial condition number: ',I3,/1X,70('*')/)

 1002 FORMAT(/1X,70('*')/,' From: CHECK_DES_IC',/'Message: Initial',   &
         ' condition data specifying the paritcle',/A,' is requied',   &
         ' when solving the ',A,'.',/' Check mfix.dat.',/1X,70('*')/)


      RETURN
      END SUBROUTINE CHECK_DES_IC
