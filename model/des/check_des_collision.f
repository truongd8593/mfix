!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_DES_COLLISION                                     !
!  Purpose: Check user input data for DES collision calculations.      !
!                                                                      !
!  Author: J.Musser                                   Date: 11-Dec-13  !
!  Reviewer:                                                           !
!  Comments: Most all user's input data are checked here               !
!  Revision: Some of the checks made in des_allocate_arrays are        !
!            moved here. In addition write statments are now made to   !
!            log file                                                  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_DES_COLLISION

! Global Variables:
!---------------------------------------------------------------------//

! User specified collision model
      USE discretelement, only: DES_COLL_MODEL, DES_COLL_MODEL_ENUM, LSD, HERTZIAN
! Number of discrete solids phases
      USE discretelement, only: DES_MMAX
! Particle and wall friction coeff.
      USE discretelement, only: MEW, MEW_W
! Particle and wall normal spring constants
      USE discretelement, only: KN, KN_W
! Particle and wall tangential spring factor := KN/KT
      USE discretelement, only: KT_FAC, KT_W_FAC
! Coefficients of restitution: Normal and Tangential
      USE discretelement, only: DES_EN_INPUT, DES_EN_WALL_INPUT
      USE discretelement, only: DES_ET_INPUT, DES_ET_WALL_INPUT
! Tangential damping factors := ET/EN
      USE discretelement, only: DES_ETAT_FAC, DES_ETAT_W_FAC
! Particle and wall Young's modulus
      USE discretelement, only: E_YOUNG, EW_YOUNG
! Particle and wall Pooisson ratio
      USE discretelement, only: V_POISSON, VW_POISSON
! Parameter constatns.
      USE param1, only: ONE, ZERO, UNDEFINED
! File unit for .LOG file.
      USE funits, only: UNIT_LOG

      USE mpi_utility


      IMPLICIT NONE


! Local Variables:
!---------------------------------------------------------------------//
! Loop index.
      INTEGER :: M
! Number of phase interactions
      INTEGER :: MxM
! Flag to warn user.
      LOGICAL :: FLAG_WARN
! Integer error flag.
      INTEGER :: IER
! Message for formated output.
      CHARACTER(len=64) :: MSG

! Initialize the error flag.
      IER = 0

! Check settings for collision models

! check particle-particle normal restitution coefficient
      MxM = DES_MMAX+DES_MMAX*(DES_MMAX-1)/2
      DO M = 1, MxM
         IF(DES_EN_INPUT(M) == UNDEFINED) THEN
            IF(myPE == PE_IO) WRITE(*, 1003)                           &
               'particle normal', 'DES_EN_INPUT', MxM
            IF(DMP_LOG) WRITE(UNIT_LOG,1003)                           &
               'particle normal', 'DES_EN_INPUT', MxM
            IER = 1
         ENDIF
      ENDDO
      DO M = 1, MxM
         IF(DES_EN_INPUT(M) > ONE .OR. DES_EN_INPUT(M) < ZERO) THEN
            MSG=''; WRITE(MSG,"('DES_EN_INPUT(',I2,')')") M
            IF(myPE == PE_IO) WRITE(*, 1002) trim(MSG)
            IF(DMP_LOG) WRITE(UNIT_LOG,1002) trim(MSG)
            IER = 1
         ENDIF
      ENDDO

! Check particle-wall normal restitution coefficient.
      DO M = 1, DES_MMAX
         IF(DES_EN_WALL_INPUT(M) == UNDEFINED) THEN
            IF(myPE == PE_IO) WRITE(*, 1003)                           &
               'wall normal', 'DES_EN_WALL_INPUT', DES_MMAX
            IF(DMP_LOG) WRITE(UNIT_LOG,1003)                           &
               'wall normal', 'DES_EN_WALL_INPUT', DES_MMAX
            IER = 1
         ENDIF
      ENDDO

      DO M = 1, DES_MMAX
         IF(DES_EN_WALL_INPUT(M) > ONE .OR.                            &
            DES_EN_WALL_INPUT(M) < ZERO) THEN
            MSG=''; WRITE(MSG,"('DES_EN_WALL_INPUT(',I2,')')") M
            IF(myPE == PE_IO) WRITE(*, 1002) trim(MSG)
            IF(DMP_LOG) WRITE(UNIT_LOG,1002) trim(MSG)
            IER = 1
         ENDIF
      ENDDO


! Check coefficient friction
      IF(MEW == UNDEFINED) THEN
         MSG='Particle friction coefficient.'
         IF(myPE == PE_IO) WRITE(*, 1001) 'MEW', MSG
         IF(DMP_LOG) WRITE(UNIT_LOG,1001) 'MEW', MSG
         IER = 1
      ELSEIF (MEW < ZERO .OR. MEW_W > ONE) THEN
         IF(myPE == PE_IO) WRITE(*, 1002) 'MEW'
         IF(DMP_LOG) WRITE(UNIT_LOG,1002) 'MEW'
         IER = 1
      ENDIF
      IF(MEW_W == UNDEFINED) THEN
         MSG='Wall friction coefficient.'
         IF(myPE == PE_IO) WRITE(*, 1001) 'MEW_W', MSG
         IF(DMP_LOG) WRITE(UNIT_LOG,1001) 'MEW_W', MSG
         IER = 1
      ELSEIF(MEW_w < ZERO .OR. MEW_W > ONE) THEN
         IF(myPE == PE_IO) WRITE(*, 1002) 'MEW_W'
         IF(DMP_LOG) WRITE(UNIT_LOG,1002) 'MEW_W'
         IER = 1
      ENDIF





! Check collision model specific model parameters.
      SELECT CASE (DES_COLL_MODEL_ENUM)

!**********************************************************************!
!*                                                                    *!
!*                     linear spring-dashpot model                    *!
!*                                                                    *!
!**********************************************************************!
      CASE (LSD)

         IF(KN == UNDEFINED) THEN
            IF(myPE == PE_IO) WRITE(*, 1101) 'KN'
            IF(DMP_LOG) WRITE(UNIT_LOG,1101) 'KN'
            IER = 1
         ENDIF
         IF(KN_W == UNDEFINED) THEN
            IF(myPE == PE_IO) WRITE(*, 1101) 'KN_W'
            IF(DMP_LOG) WRITE(UNIT_LOG,1101) 'KN_W'
            IER = 1
         ENDIF

! Check for tangential spring constant factors
         IF(KT_FAC == UNDEFINED) THEN
            IF(DMP_LOG) WRITE (UNIT_LOG, 1102) 'KT_FAC'
         ELSEIF(KT_FAC > ONE .OR. KT_FAC < ZERO) THEN
            IF(myPE == PE_IO) WRITE(*, 1103) 'KT_FAC'
            IF(DMP_LOG) WRITE(UNIT_LOG,1103) 'KT_FAC'
            IER = 1
         ENDIF
         IF(KT_W_FAC == UNDEFINED) THEN
            IF(DMP_LOG) WRITE (UNIT_LOG, 1102) 'KT_W_FAC'
         ELSEIF(KT_W_FAC > ONE .OR. KT_W_FAC < ZERO) THEN
            IF(myPE == PE_IO) WRITE(*, 1103) 'KT_W_FAC'
            IF(DMP_LOG) WRITE(UNIT_LOG,1103) 'KT_W_FAC'
            IER = 1
         ENDIF

! Check for tangential damping factor
         IF(DES_ETAT_FAC == UNDEFINED) THEN
            IF(DMP_LOG) WRITE (UNIT_LOG, 1104) 'DES_ETAT_FAC'
         ELSEIF(DES_ETAT_FAC > ONE .OR. DES_ETAT_FAC < ZERO) THEN
            IF(myPE == PE_IO) WRITE(*, 1103) 'DES_ETAT_FAC'
            IF(DMP_LOG) WRITE(UNIT_LOG,1103) 'DES_ETAT_FAC'
            IER = 1
         ENDIF
         IF(DES_ETAT_W_FAC == UNDEFINED) THEN
            IF(DMP_LOG) WRITE (UNIT_LOG, 1104) 'DES_ETAT_W_FAC'
         ELSEIF(DES_ETAT_W_FAC > ONE .OR. DES_ETAT_W_FAC < ZERO) THEN
            IF(myPE == PE_IO) WRITE(*, 1103) 'DES_ETAT_W_FAC'
            IF(DMP_LOG) WRITE(UNIT_LOG,1103) 'DES_ETAT_W_FAC'
            IER = 1
         ENDIF

! if following are assigned warn user they are discarded
         FLAG_WARN = .FALSE.
         DO M = 1, DES_MMAX+DES_MMAX*(DES_MMAX-1)/2
            IF(DES_ET_INPUT(M) .NE. UNDEFINED) FLAG_WARN = .TRUE.
         ENDDO
         IF (FLAG_WARN) WRITE(UNIT_LOG,1105) 'DES_ET_INPUT'
         FLAG_WARN = .FALSE.
         DO M = 1, DES_MMAX
            IF(DES_ET_WALL_INPUT(M) .NE. UNDEFINED) FLAG_WARN = .TRUE.
         ENDDO
         IF (FLAG_WARN) WRITE(UNIT_LOG,1105) 'DES_ET_WALL_INPUT'
         FLAG_WARN = .FALSE.



!**********************************************************************!
!*                                                                    *!
!*                           Hertzian model                           *!
!*                                                                    *!
!**********************************************************************!
      CASE (HERTZIAN)

! check young's modulus and poisson ratio
         IF(EW_YOUNG == UNDEFINED ) THEN
            MSG='Wall value for Youngs modulus'
            IF(myPE == PE_IO) WRITE(*, 1201) 'EW_YOUNG', MSG
            IF(DMP_LOG) WRITE(UNIT_LOG,1201) 'EW_YOUNG', MSG
            IER = 1
         ENDIF

         IF(VW_POISSON == UNDEFINED) THEN
            MSG='Wall value for Poissons ratio'
            IF(myPE == PE_IO) WRITE(*, 1201) 'VW_POISSON', MSG
            IF(DMP_LOG) WRITE(UNIT_LOG,1201) 'VW_POISSON', MSG
            IER = 1
         ELSEIF (VW_POISSON > 0.5d0 .OR. VW_POISSON <= -ONE) THEN
            IF(myPE == PE_IO) WRITE(*, 1202) 'VW_POISSON'
            IF(DMP_LOG) WRITE(UNIT_LOG,1202) 'VW_POISSON'
            IER = 1
         ENDIF

         DO M = 1, DES_MMAX
            IF(E_YOUNG(M) == UNDEFINED) THEN
               MSG=''; WRITE(MSG,"('Phase ',I2,' Youngs modulus')") M
               IF(myPE == PE_IO) WRITE(*, 1201) 'E_YOUNG', MSG
               IF(DMP_LOG) WRITE(UNIT_LOG,1201) 'E_YOUNG', MSG
               IER = 1
            ENDIF
            IF(V_POISSON(M) == UNDEFINED) THEN
               MSG=''; WRITE(MSG,"('Phase ',I2,' Poissons ratio')") M
               IF(myPE == PE_IO) WRITE(*, 1201) 'V_POISSON', MSG
               IF(DMP_LOG) WRITE(UNIT_LOG,1201) 'V_POISSON', MSG
               IER = 1
            ELSEIF(V_POISSON(M) > 0.5d0 .OR. V_POISSON(M) <= -ONE) THEN
               IF(myPE == PE_IO) WRITE(*, 1202) 'V_POISSON'
               IF(DMP_LOG) WRITE(UNIT_LOG,1202) 'V_POISSON'
               IER = 1
            ENDIF
         ENDDO

! check particle-particle tangential restitution coefficient
         MxM = DES_MMAX+DES_MMAX*(DES_MMAX-1)/2
         DO M = 1, MxM
            IF(DES_ET_INPUT(M) == UNDEFINED) THEN
               IF(myPE == PE_IO) WRITE(*, 1203)                        &
                  'particle', 'DES_ET_INPUT', MxM
               IF(DMP_LOG) WRITE(UNIT_LOG,1203)                        &
                  'particle', 'DES_ET_INPUT', MxM
               IER = 1
            ENDIF
         ENDDO
         DO M = 1, MxM
            IF(DES_ET_INPUT(M) > ONE .OR. DES_ET_INPUT(M) < ZERO) THEN
               MSG=''; WRITE(MSG,"('DES_ET_INPUT(',I2,')')") M
               IF(myPE == PE_IO) WRITE(*, 1204) trim(MSG)
               IF(DMP_LOG) WRITE(UNIT_LOG,1204) trim(MSG)
               IER = 1
            ENDIF
         ENDDO

! check particle-wall tangential restitution coefficient
         DO M = 1, DES_MMAX
            IF(DES_ET_WALL_INPUT(M) == UNDEFINED) THEN
               IF(myPE == PE_IO) WRITE(*, 1203)                        &
                  'wall', 'DES_ET_WALL_INPUT', DES_MMAX
               IF(DMP_LOG) WRITE(UNIT_LOG,1203)                        &
                  'wall', 'DES_ET_WALL_INPUT', DES_MMAX
               IER = 1
            ENDIF
         ENDDO
         DO M = 1, DES_MMAX
            IF(DES_ET_WALL_INPUT(M) > ONE .OR.                         &
               DES_ET_WALL_INPUT(M) < ZERO) THEN
               MSG=''; WRITE(MSG,"('DES_ET_WALL_INPUT(',I2,')')") M
               IF(myPE == PE_IO) WRITE(*, 1204) trim(MSG)
               IF(DMP_LOG) WRITE(UNIT_LOG,1204) trim(MSG)
               IER = 1
            ENDIF
         ENDDO

! If following are assigned warn user they are discarded.
         IF(DMP_LOG)THEN
            IF(KN .NE. UNDEFINED)                                      &
               WRITE(UNIT_LOG, 1205) 'KN'
            IF(KN_W .NE. UNDEFINED)                                    &
               WRITE(UNIT_LOG, 1205) 'KN_W'
            IF(KT_FAC .NE. UNDEFINED)                                  &
               WRITE(UNIT_LOG, 1205) 'KT_FAC'
            IF(KT_W_FAC .NE. UNDEFINED)                                &
               WRITE(UNIT_LOG, 1205) 'KT_W_FAC'
            IF(DES_ETAT_FAC .NE. UNDEFINED)                            &
               WRITE(UNIT_LOG, 1205) 'DES_ETAT_FAC'
            IF(DES_ETAT_W_FAC .NE. UNDEFINED)                          &
               WRITE(UNIT_LOG, 1205) 'DES_ETAT_W_FAC'
         ENDIF


! Unknown collision model.
      CASE DEFAULT

         IF(myPE == PE_IO) WRITE(*, 1000)TRIM(DES_COLL_MODEL)
         IF(DMP_LOG) WRITE(UNIT_LOG,1000) TRIM(DES_COLL_MODEL)
         IER =1

      END SELECT


     IF(IER /=0) CALL MFIX_EXIT(myPE)



      RETURN


 1000 FORMAT(/1X,70('*')/' From: CHECK_DES_COLLISION',/' Error 1000:', &
         ' Invalid collision model specified: ',A,'.',/                &
         ' Available DEM collision modles include:',/                  &
         '   > LSD ',/'   > HERTZIAN',/1X,70('*')/)

 1001 FORMAT(/1X,70('*')/' From: CHECK_DES_COLLISION',/' Error 1001:', &
         ' Undefined collision parameter: ',A,/' Description: ',A,/    &
         1X,70('*')/)

 1002 FORMAT(/1X,70('*')/' From: CHECK_DES_COLLISION',/' Error 1002:', &
         ' Unphysical value specified for ',A,/' Valid parameter',     &
         ' range is between zero and one.'/1X,70('*')/)

 1003 FORMAT(/1X,70('*')/' From: CHECK_DES_COLLISION',/' Error 1003:', &
         ' Particle-',A,' restitution coefficient',/1x,A,              &
         ' must be specified for all ',I3,' interactions.',/1X,70('*')/)


 1101 FORMAT(/1X,70('*')/' From: CHECK_DES_COLLISION',/' Error 1101:', &
         ' Spring constants ',A,' not specified in mfix.dat',/1X,70('*'))



 1102 FORMAT(/1X,70('*')/' From: CHECK_DES_COLLISION',/' Warning 1102:',&
         ' Tangential spring factor ',A,' not specified in mfix.dat.',/&
         ' It will be defined in cfassign.f as 2/7.',/1X,70('*')/)

 1103 FORMAT(/1X,70('*')/' From: CHECK_DES_COLLISION',/' Error 1103:', &
         ' Unphysical value specified for ',A,/' Valid parameter',     &
         ' range is between zero and one.'/1X,70('*')/)

 1104 FORMAT(/1X,70('*')/' From: CHECK_DES_COLLISION',/' Warning 1104:',&
         ' Tangential damping factor ',A,' not specified.',/           &
         ' The default of 1/2 will be set in cfassign.f. ',/1X,70('*')/)

 1105 FORMAT(/1X,70('*')/' From: CHECK_DES_COLLISION',/' Warning 1105:',&
         ' Specified values for ',A,' are not used',/' with the',       &
         ' default collision model.',/1X,70('*')/)




 1201 FORMAT(/1X,70('*')/' From: CHECK_DES_COLLISION',/' Error 1201:', &
         ' Undefined collision parameter: ',A,/' Description: ',A,/    &
         1X,70('*')/)

 1202 FORMAT(/1X,70('*')/' From: CHECK_DES_COLLISION',/' Error 1202:', &
         ' Unphysical value specified for ',A,/' Valid parameter',     &
         ' range is >0.50 or =< 1.0'/1X,70('*')/)

 1203 FORMAT(/1X,70('*')/' From: CHECK_DES_COLLISION',/' Error 1203:', &
         'Particle-',A,' tangential restitution coefficient',/1x,A,    &
         ' must be specified for all ',I3,' interactions.',/1X,70('*')/)

 1204 FORMAT(/1X,70('*')/' From: CHECK_DES_COLLISION',/' Error 1204:', &
         ' Unphysical value specified for ',A,/' Valid parameter',     &
         ' range is between zero and one.'/1X,70('*')/)

 1205 FORMAT(/1X,70('*')/' From: CHECK_DES_COLLISION',/' Warning 1205:',&
         'Specified parameter not used in Hertz model',A,/1X,70('*')/)


         END SUBROUTINE CHECK_DES_COLLISION
