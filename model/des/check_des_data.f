!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CHECK_DES_DATA                                         C
!>  Purpose: DES - Check user input Data   
!                                                                      C
!                                                                      C
!  Author: S. Benyahia                                Date: 10-Nov-08  C
!  Reviewer:                                                           C
!  Comments: All user's input data are checked here                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CHECK_DES_DATA
      
      USE param1
      USE geometry
      USE funits
      USE compar      
      USE discretelement
      USE run
      USE constant
      USE physprop

      IMPLICIT NONE
      
      INTEGER CHECK_MPI, M
      
      IF(COORDINATES == 'CYLINDRICAL') THEN
         WRITE (UNIT_LOG, 1000)
         CALL MFIX_EXIT(myPE)
      END IF

      CHECK_MPI = NODESI * NODESJ * NODESK
      IF(CHECK_MPI.NE.1) THEN
         WRITE (UNIT_LOG, 1001)
         CALL MFIX_EXIT(myPE)
      END IF
!
      IF(DES_PERIODIC_WALLS) THEN
         IF(.NOT.DES_PERIODIC_WALLS_X .AND. .NOT.DES_PERIODIC_WALLS_Y .AND. &
	    .NOT.DES_PERIODIC_WALLS_Z) THEN
	    WRITE (UNIT_LOG, 1002)
            CALL MFIX_EXIT(myPE)
         ENDIF
      END IF
!
      IF(INLET_OUTLET) THEN
         IF(.NOT.INLET_OUTLET_X .AND. .NOT.INLET_OUTLET_Y .AND. &
	    .NOT.INLET_OUTLET_Z) THEN
	    WRITE (UNIT_LOG, 1003)
            CALL MFIX_EXIT(myPE)
         ENDIF
      END IF
!
      IF(KN == UNDEFINED .OR. KT == UNDEFINED) THEN
         WRITE (UNIT_LOG, 1004)
         CALL MFIX_EXIT(myPE)
      END IF
!
      IF(KN_W == UNDEFINED .OR. KT_W == UNDEFINED) THEN
         WRITE (UNIT_LOG, 1005)
         CALL MFIX_EXIT(myPE)
      END IF
!
      IF(ETA_DES_N .NE. UNDEFINED .OR. ETA_DES_T .NE. UNDEFINED .OR. &
         ETA_N_W   .NE. UNDEFINED .OR. ETA_T_W   .NE. UNDEFINED ) THEN
         WRITE (UNIT_LOG, 1006)
         CALL MFIX_EXIT(myPE)
      END IF
!
      DO M = 1, MMAX+MMAX*(MMAX-1)/2
         IF(DES_EN_INPUT(M) == UNDEFINED) THEN
            WRITE (UNIT_LOG, 1008)
	    write(UNIT_LOG,*) 'Must be specified in mfix.dat for interactions M = 1, to', MMAX+MMAX*(MMAX-1)/2
            WRITE (UNIT_LOG, 1009)
            CALL MFIX_EXIT(myPE)
         END IF
      ENDDO
!
      DO M = 1, MMAX
         IF(DES_EN_WALL_INPUT(M) == UNDEFINED) THEN
            WRITE (UNIT_LOG, 1010)
            CALL MFIX_EXIT(myPE)
         END IF
      ENDDO
!
      DO M = 1, MMAX+MMAX*(MMAX-1)/2
         IF(DES_ET_INPUT(M) == UNDEFINED) THEN
            WRITE (UNIT_LOG, 1011)
	    write(UNIT_LOG,*) 'Must be specified in mfix.dat for interactions M = 1, to', MMAX+MMAX*(MMAX-1)/2
            WRITE (UNIT_LOG, 1009)
            CALL MFIX_EXIT(myPE)
         END IF
      ENDDO
!
      DO M = 1, MMAX
         IF(DES_ET_WALL_INPUT(M) == UNDEFINED) THEN
            WRITE (UNIT_LOG, 1012)
            CALL MFIX_EXIT(myPE)
         END IF
      ENDDO
!
      DO M = 1, MMAX+MMAX*(MMAX-1)/2
         IF(DES_EN_INPUT(M) > ONE .OR. DES_ET_INPUT(M) > ONE .OR. DES_EN_INPUT(M) < ZERO .OR. &
	    DES_ET_INPUT(M) < ZERO) THEN
            WRITE (UNIT_LOG, 1013)
            CALL MFIX_EXIT(myPE)
         END IF
      ENDDO
!
      DO M = 1, MMAX
         IF(DES_EN_WALL_INPUT(M) > ONE .OR. DES_ET_WALL_INPUT(M) > ONE .OR. &
	    DES_EN_WALL_INPUT(M) < ZERO .OR. DES_ET_WALL_INPUT(M) < ZERO ) THEN
            WRITE (UNIT_LOG, 1014)
            CALL MFIX_EXIT(myPE)
         END IF
      ENDDO
!
      IF(MEW > ONE .OR. MEW_W > ONE .OR. MEW < ZERO .OR. MEW_W < ZERO) THEN
         WRITE (UNIT_LOG, 1015)
         CALL MFIX_EXIT(myPE)
      END IF
      
      RETURN
 1000 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA.',/' Message: ',&
         'DES should only be run using CARTESIAN coordinates',/1X,70('*')/)
 1001 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA.',/' Message: ',&
         'DES being run on multiple processors. Only serial runs allowed',/1X,70('*')/)
 1002 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA.',/' Message: ',&
         'Direction of periodicity not defined in mfix.dat',/1X,70('*')/)
 1003 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA.',/' Message: ',&
         'Direction of inlet_outlet not defined in mfix.dat',/1X,70('*')/)
 1004 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA.',/' Message: ',&
         'Spring constants KN or KT not specified in mfix.dat',/1X,70('*')/)
 1005 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA.',/' Message: ',&
         'Spring constants KN_W or KT_W not specified in mfix.dat',/1X,70('*')/)
 1006 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA.',/' Message: ',&
         'Damping constants MUST NOT be specified in mfix.dat,'/10X,&
         'Specify restitution coefficients instead',/1X,70('*')/)
 1007 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA.',/' Message: ',&
         'Friction coefficients MEW or MEW_W not specified in mfix.dat',/1X,70('*')/)
 1008 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA.',/' Message: ',&
         'Particle-particle restitution coefficient DES_EN_INPUT(M)')
 1009 FORMAT(/1X,70('*')/)
 1010 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA.',/' Message: ',&
         'Particle-wall restitution coefficients DES_EN_WALL_INPUT(M),'/10X,&
         'Must be specified in mfix.dat for interactions M = 1 to MMAX',/1X,70('*')/)
 1011 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA.',/' Message: ',&
         'Particle-particle Tangential restitution coefficient DES_ET_INPUT(M)')
 1012 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA.',/' Message: ',&
         'Particle-wall Tangential restitution coefficients DES_ET_WALL_INPUT(M),'/10X,&
         'Must be specified in mfix.dat for interactions M = 1 to MMAX',/1X,70('*')/)
 1013 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA.',/' Message: ',&
         'Unphysical ( > 1 or < 0) values of part-part coefficients of restitution',/1X,70('*')/)
 1014 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA.',/' Message: ',&
         'Unphysical ( > 1 or < 0) values of part-wall coefficients of restitution',/1X,70('*')/)
 1015 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA.',/' Message: ',&
         'Unphysical ( > 1 or < 0) values of friction coefficients',/1X,70('*')/)
      END SUBROUTINE CHECK_DES_DATA
