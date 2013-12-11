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

      USE param1
      USE geometry
      USE funits
      USE compar      
      USE discretelement
      USE run
      USE constant
      USE physprop
      use desgrid 
      USE indices
      USE fldvar
      USE toleranc 
      USE mpi_utility
      USE output 
      USE mfix_pic
      USE cutcell
      USE qmom_kinetic_equation
      USE ic 

      IMPLICIT NONE
!-----------------------------------------------
! Local Variables
!-----------------------------------------------      
      INTEGER :: M
      INTEGER :: I, J, K, IJK
      LOGICAL :: FLAG_WARN


! Check settings for collision models
      IF (DES_COLL_MODEL /= UNDEFINED_C) THEN      
         IF (TRIM(DES_COLL_MODEL) .NE. 'HERTZIAN') THEN
            IF(DMP_LOG) WRITE(UNIT_LOG, 1006)
            CALL MFIX_EXIT(myPE)
         ENDIF
      ENDIF

      IF (TRIM(DES_COLL_MODEL) .EQ. 'HERTZIAN') THEN
! check young's modulus and poisson ratio
         IF(EW_YOUNG == UNDEFINED ) THEN
            IF(DMP_LOG) WRITE (UNIT_LOG, 1020)
            CALL MFIX_EXIT(myPE)
         ENDIF
         IF(VW_POISSON == UNDEFINED) THEN
            IF(DMP_LOG) WRITE (UNIT_LOG, 1021)
            CALL MFIX_EXIT(myPE)
         ELSE
            IF (VW_POISSON > 0.5d0 .OR. VW_POISSON <= -ONE) THEN
               IF(DMP_LOG) WRITE (UNIT_LOG, 1033)
            ENDIF
         ENDIF
         DO M = 1, DES_MMAX
            IF(E_YOUNG(M) == UNDEFINED) THEN
               IF(DMP_LOG) WRITE (UNIT_LOG, 1022)
               CALL MFIX_EXIT(myPE)
            ENDIF
            IF(V_POISSON(M) == UNDEFINED) THEN
               IF(DMP_LOG) WRITE (UNIT_LOG, 1023)
               CALL MFIX_EXIT(myPE)
            ELSE
               IF(V_POISSON(M) > 0.5d0 .OR. V_POISSON(M) <= -ONE) THEN
                  IF(DMP_LOG) WRITE (UNIT_LOG, 1033)
               ENDIF
            ENDIF
         ENDDO
! check particle-particle tangential restitution coefficient      
         DO M = 1, DES_MMAX+DES_MMAX*(DES_MMAX-1)/2
            IF(DES_ET_INPUT(M) == UNDEFINED) THEN
               IF(DMP_LOG) WRITE (UNIT_LOG, 1024) DES_MMAX+DES_MMAX*(DES_MMAX-1)/2
               CALL MFIX_EXIT(myPE)
            ENDIF
         ENDDO
         DO M = 1, DES_MMAX+DES_MMAX*(DES_MMAX-1)/2
            IF(DES_ET_INPUT(M) > ONE .OR. DES_ET_INPUT(M) < ZERO) THEN
               IF(DMP_LOG) WRITE (UNIT_LOG, 1025)
               CALL MFIX_EXIT(myPE)
            ENDIF
         ENDDO
! check particle-wall tangential restitution coefficient
         DO M = 1, DES_MMAX
            IF(DES_ET_WALL_INPUT(M) == UNDEFINED) THEN
               IF(DMP_LOG) WRITE (UNIT_LOG, 1026)
               CALL MFIX_EXIT(myPE)
            ENDIF
         ENDDO
         DO M = 1, DES_MMAX
            IF(DES_ET_WALL_INPUT(M) > ONE .OR. DES_ET_WALL_INPUT(M) < ZERO) THEN
               IF(DMP_LOG) WRITE (UNIT_LOG, 1027)
               CALL MFIX_EXIT(myPE)
            ENDIF
         ENDDO
! if following are assigned warn user they are discarded
         IF(KN .NE. UNDEFINED .OR. KN_W .NE. UNDEFINED) THEN 
            IF(DMP_LOG) WRITE (UNIT_LOG, 1028)
         ENDIF
         IF(KT_FAC .NE. UNDEFINED .OR. KT_W_FAC .NE. UNDEFINED) then 
            IF(DMP_LOG) WRITE (UNIT_LOG, 1029)
         ENDIF
         IF(DES_ETAT_FAC .NE. UNDEFINED .OR. &
            DES_ETAT_W_FAC .NE. UNDEFINED) then 
            IF(DMP_LOG) WRITE (UNIT_LOG, 1030)
         ENDIF


      ELSE  ! default linear spring-dashpot model

! check normal spring constants 
         IF(KN == UNDEFINED .OR. KN_W == UNDEFINED) THEN
            IF(DMP_LOG) WRITE (UNIT_LOG, 1004)
            CALL MFIX_EXIT(myPE)
         ENDIF
! Check for tangential spring constant factors         
         IF(KT_FAC == UNDEFINED .OR. KT_W_FAC == UNDEFINED) THEN
            IF(DMP_LOG) WRITE (UNIT_LOG, 1005)
         ENDIF
         IF(KT_FAC .NE. UNDEFINED) THEN
            IF(KT_FAC > ONE .OR. KT_FAC < ZERO) THEN
               IF(DMP_LOG) WRITE (UNIT_LOG, 1016)
               CALL MFIX_EXIT(myPE)
            ENDIF
         ENDIF
         IF(KT_W_FAC .NE. UNDEFINED) THEN
            IF(KT_W_FAC > ONE .OR. KT_W_FAC < ZERO) THEN
               IF(DMP_LOG) WRITE (UNIT_LOG, 1016)
               CALL MFIX_EXIT(myPE)
            ENDIF
         ENDIF
! check for tangential damping factor
         IF(DES_ETAT_FAC == UNDEFINED) THEN
            IF(DMP_LOG) WRITE (UNIT_LOG, 1010)
         ELSEIF(DES_ETAT_FAC > ONE .OR. DES_ETAT_FAC < ZERO) THEN
            IF(DMP_LOG) WRITE (UNIT_LOG, 1015)
            CALL MFIX_EXIT(myPE)
         ENDIF
         IF(DES_ETAT_W_FAC == UNDEFINED) THEN
            IF(DMP_LOG) WRITE (UNIT_LOG, 1011)
         ELSEIF(DES_ETAT_W_FAC > ONE .OR. DES_ETAT_W_FAC < ZERO) THEN
            IF(DMP_LOG) WRITE (UNIT_LOG, 1015)
            CALL MFIX_EXIT(myPE)
         ENDIF
! if following are assigned warn user they are discarded
         FLAG_WARN = .FALSE.
         DO M = 1, DES_MMAX+DES_MMAX*(DES_MMAX-1)/2
            IF(DES_ET_INPUT(M) .NE. UNDEFINED) FLAG_WARN = .TRUE.
         ENDDO
         IF (FLAG_WARN) WRITE(UNIT_LOG,1031)
         FLAG_WARN = .FALSE.
         DO M = 1, DES_MMAX
            IF(DES_ET_WALL_INPUT(M) .NE. UNDEFINED) FLAG_WARN = .TRUE.
         ENDDO
         IF (FLAG_WARN) WRITE(UNIT_LOG,1032)
         FLAG_WARN = .FALSE.

      ENDIF   ! end if/else(des_coll_model.eq.hertzian)

! check particle-particle normal restitution coefficient      
      DO M = 1, DES_MMAX+DES_MMAX*(DES_MMAX-1)/2
         IF(DES_EN_INPUT(M) == UNDEFINED) THEN
            IF(DMP_LOG) WRITE (UNIT_LOG, 1008) DES_MMAX+DES_MMAX*(DES_MMAX-1)/2
            CALL MFIX_EXIT(myPE)
         ENDIF
      ENDDO
      DO M = 1, DES_MMAX+DES_MMAX*(DES_MMAX-1)/2
         IF(DES_EN_INPUT(M) > ONE .OR. DES_EN_INPUT(M) < ZERO) THEN
            IF(DMP_LOG) WRITE (UNIT_LOG, 1012)
            CALL MFIX_EXIT(myPE)
         ENDIF
      ENDDO

! check particle-wall normal restitution coefficient (needed by all
! current collision models)
! MPPIC has separate key words now for the particle-wall interaction
      DO M = 1, DES_MMAX
         IF(DES_EN_WALL_INPUT(M) == UNDEFINED) THEN
            IF(DMP_LOG) WRITE (UNIT_LOG, 1009)
            CALL MFIX_EXIT(myPE)
         ENDIF
      ENDDO
      DO M = 1, DES_MMAX
         IF(DES_EN_WALL_INPUT(M) > ONE .OR. DES_EN_WALL_INPUT(M) < ZERO) THEN
            IF(DMP_LOG) WRITE (UNIT_LOG, 1013)
            CALL MFIX_EXIT(myPE)
         ENDIF
      ENDDO

! check coefficient friction 
      IF(MEW /= UNDEFINED .AND. MEW_W /= UNDEFINED) THEN
         IF (MEW> ONE .OR. MEW_W > ONE .OR. &
             MEW < ZERO .OR. MEW_W < ZERO) THEN
            IF(DMP_LOG) WRITE (UNIT_LOG, 1014)
            CALL MFIX_EXIT(myPE)
         ENDIF
      ELSE
         WRITE(UNIT_LOG, 1035)
         CALL MFIX_EXIT(myPE)
      ENDIF


      RETURN


 1004 FORMAT(/1X,70('*')//' From: CHECK_DES_COLLISION',/' Message: ',&
         'Spring constants KN or KN_W not specified in mfix.dat',&
         /1X,70('*')/)

 1005 FORMAT(/1X,70('*')//' From: CHECK_DES_COLLISION',/' Message: ',&
         'WARNING: Tangential spring factors KT_FAC or KT_W_FAC not',/10X,&
         'specified in mfix.dat.  These factors will be defined in ',&
         'cfassign.f',/10X,'as 2/7.  See subroutine for references.',&
         /1X,70('*')/)

 1006 FORMAT(/1X,70('*')//' From: CHECK_DES_COLLISION',/' Message: ',&
         'Change the value DES_COLL_MODEL in mfix.dat. Options are',/10X,&
         'leave it undefined to use the (default) linear spring-',/10X,&
         'dashpot model or set it to HERTZIAN for hertz model.',&
          /1X,70('*')/)

 1008 FORMAT(/1X,70('*')//' From: CHECK_DES_COLLISION',/' Message: ',&
         'Particle-particle restitution coefficient DES_EN_INPUT(M)',/10X,&
         'not specified in mfix.dat for interactions M = 1 to ',I5,&
          /1X,70('*')/)

 1009 FORMAT(/1X,70('*')//' From: CHECK_DES_COLLISION',/' Message: ',&
         'Particle-wall restitution coefficients DES_EN_WALL_INPUT(M)',&
         /10X,'not specified in mfix.dat for interactions M= 1 to MMAX',&
         /1X,70('*')/)

 1010 FORMAT(/1X,70('*')//' From: CHECK_DES_COLLISION',/' Message: ',&
         'WARNING: Tangential damping factors DES_ETAT_FAC not ',&
         'specified in',/10X, 'mfix.dat. This factor will be set in ',&
         'cfassign.f as 1/2. See subroutine for references.',&
         /1X,70('*')/)

 1011 FORMAT(/1X,70('*')//' From: CHECK_DES_COLLISION',/' Message: ',&
         'WARNING: Tangential damping factors DES_ETAT_W_FAC not ',&
         'specified in mfix.dat',/10X,'This factor will be set in ',&
         'cfassign.f as 1/2. See subroutine for references.',&
         /1X,70('*')/)

 1012 FORMAT(/1X,70('*')//' From: CHECK_DES_COLLISION',/' Message: ',&
         'Unphysical ( > 1 or < 0) values of DES_EN_INPUT(M)',&
         /1X,70('*')/)

 1013 FORMAT(/1X,70('*')//' From: CHECK_DES_COLLISION',/' Message: ',&
         'Unphysical ( > 1 or < 0) values of DES_EN_WALL_INPUT(M)',&
         /1X,70('*')/)

 1014 FORMAT(/1X,70('*')//' From: CHECK_DES_COLLISION',/' Message: ',&
         'Unphysical ( > 1 or < 0) values of friction coefficients',&
         /1X,70('*')/)

 1015 FORMAT(/1X,70('*')//' From: CHECK_DES_COLLISION',/' Message: ',&
         'Values of DES_ETAT_FAC or DES_ETAT_W_FAC unphysical ',/10X,&
         '(< 0 or > 1) defined in mfix.dat',/1X,70('*')/)

 1016 FORMAT(/1X,70('*')//' From: CHECK_DES_COLLISION',/' Message: ',&
         'Values of KT_FAC or KT_W_FAC unphysical (< 0 or > 1) ',/10X,&
         'defined in mfix.dat',/1X,70('*')/)


 1020 FORMAT(/1X,70('*')//' From: CHECK_DES_COLLISION',/' Message: ',&
         'Wall value for Youngs modulus (EW_YOUNG) must be,'/10X,&
         'specified in mfix.dat',/1X,70('*')/)

 1021 FORMAT(/1X,70('*')//' From: CHECK_DES_COLLISION',/' Message: ',&
         'Wall value for Poissons ratio (VW_POISSON) must be',/10X,&
         'specified in mfix.dat',/1X,70('*')/)

 1022 FORMAT(/1X,70('*')//' From: CHECK_DES_COLLISION',/' Message: ',&
         'Youngs modulus (E_YOUNG) must be specified in mfix.dat',/10X,&
         'for all particle types M = 1 to MMAX',/1X,70('*')/)

 1023 FORMAT(/1X,70('*')//' From: CHECK_DES_COLLISION',/' Message: ',&
         'Poissons ratio (V_POISSON) must be specified in mfix.dat',/10X,&
         'for all particle types M = 1 to MMAX',/1X,70('*')/)

 1024 FORMAT(/1X,70('*')//' From: CHECK_DES_COLLISION',/' Message: ',&
         'Particle-particle tangential restitution coefficient',/10X,&
         'DES_ET_INPUT(M) must be specified in mfix.dat for all',/10X,&
         'interactions M = 1 to ',I5,/1X,70('*')/)

 1025 FORMAT(/1X,70('*')//' From: CHECK_DES_COLLISION',/' Message: ',&
         'Unphysical ( > 1 or < 0) values of DES_ET_INPUT(M)',/1X,70('*')/)

 1026 FORMAT(/1X,70('*')//' From: CHECK_DES_COLLISION',/' Message: ',&
         'Particle-wall tangential restitution coefficient',/10X,&
         'DES_ET_WALL_INPUT(M) must be specified in mfix.dat for all',/10X,&
         'wall interactions M = 1 to MMAX',/1X,70('*')/)

 1027 FORMAT(/1X,70('*')//' From: CHECK_DES_COLLISION',/' Message: ',&
         'Unphysical ( > 1 or < 0) values of DES_ET_WALL_INPUT(M)',/1X,70('*')/)

 1028 FORMAT(/1X,70('*')//' From: CHECK_DES_COLLISION',/' Message: ',&
         'WARNING: specified values of KN and KN_W are not used with',/10X,&
         'Hertz collision model',/1X,70('*')/)

 1029 FORMAT(/1X,70('*')//' From: CHECK_DES_COLLISION',/' Message: ',&
         'WARNING: specified values of KT_FAC and KT_W_WALL are not',/10X,&
         'used with Hertz collision model',/1X,70('*')/)

 1030 FORMAT(/1X,70('*')//' From: CHECK_DES_COLLISION',/' Message: ',&
         'WARNING: specified values of DES_ETAT_FAC and',/10X,&
         'DES_ETAT_W_FAC are not used with Hertz collision model',/1X,70('*')/)

 1031 FORMAT(/1X,70('*')//' From: CHECK_DES_COLLISION',/' Message: ',&
         'WARNING: specified values of DES_ET_INPUT(M) are not used',/10X,&
         'with default collision model',/1X,70('*')/)

 1032 FORMAT(/1X,70('*')//' From: CHECK_DES_COLLISION',/' Message: ',&
         'WARNING: specified values of DES_ET_WALL_INPUT(M) are not',/10X,&
         'used with default collision model',/1X,70('*')/)

 1033 FORMAT(/1X,70('*')//' From: CHECK_DES_COLLISION',/' Message: ',&
         'WARNING: Values of VW_POISSON OR V_POISSON unphysical',/10X,&
         '(> 0.50 or =< -1.d0) defined in mfix.dat',/1X,70('*')/)


 1035 FORMAT(/1X,70('*')//' From: CHECK_DES_COLLISION',/' Message: ',&
         'Friction coefficients (MEW and MEW_w) must be ',&
         'specified in mfix.dat',/1X,70('*')/)

         END SUBROUTINE CHECK_DES_COLLISION
