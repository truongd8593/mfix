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

      IF(DES_PERIODIC_WALLS) THEN
         IF(.NOT.DES_PERIODIC_WALLS_X .AND. .NOT.DES_PERIODIC_WALLS_Y .AND. &
            .NOT.DES_PERIODIC_WALLS_Z) THEN
            WRITE (UNIT_LOG, 1002)
            CALL MFIX_EXIT(myPE)
         ENDIF
         IF (DES_NEIGHBOR_SEARCH .EQ. 1) THEN
            WRITE (UNIT_LOG, 1018)
            WRITE (*,1018)
         ENDIF
         IF (DES_NEIGHBOR_SEARCH .EQ. 2 .OR. &
          DES_NEIGHBOR_SEARCH .EQ. 3) THEN
            WRITE (UNIT_LOG, 1019)
            CALL MFIX_EXIT(myPE)
         ENDIF
      END IF
! Code commented since inlet_outlet option is not working in the current code. Dec 04 2008 (sof)
!      IF(INLET_OUTLET) THEN
!         IF(.NOT.INLET_OUTLET_X .AND. .NOT.INLET_OUTLET_Y .AND. &
!	    .NOT.INLET_OUTLET_Z) THEN
!	    WRITE (UNIT_LOG, 1003)
!            CALL MFIX_EXIT(myPE)
!         ENDIF
!      END IF
!
      IF(INLET_OUTLET .OR. INLET_OUTLET_X .OR. INLET_OUTLET_Y .OR. &
            INLET_OUTLET_Z) THEN
         WRITE (UNIT_LOG, 1003)
         CALL MFIX_EXIT(myPE)
      END IF

      IF(KN == UNDEFINED .OR. KN_W == UNDEFINED) THEN
         WRITE (UNIT_LOG, 1004)
         CALL MFIX_EXIT(myPE)
      END IF

      IF(KT_FAC == UNDEFINED .OR. KT_W_FAC == UNDEFINED) THEN
         WRITE (UNIT_LOG, 1005)
      END IF

      IF(KT_FAC .NE. UNDEFINED) THEN
         IF(KT_FAC > ONE .OR. KT_FAC < ZERO) THEN
            WRITE (UNIT_LOG, 1017)
            CALL MFIX_EXIT(myPE)
         ENDIF
      ENDIF

      IF(KT_W_FAC .NE. UNDEFINED) THEN
         IF(KT_W_FAC > ONE .OR. KT_W_FAC < ZERO) THEN
            WRITE (UNIT_LOG, 1017)
            CALL MFIX_EXIT(myPE)
         ENDIF
      ENDIF

      DO M = 1, MMAX+MMAX*(MMAX-1)/2
         IF(DES_EN_INPUT(M) == UNDEFINED) THEN
            WRITE (UNIT_LOG, 1008) MMAX+MMAX*(MMAX-1)/2
            CALL MFIX_EXIT(myPE)
         END IF
      ENDDO

      DO M = 1, MMAX
         IF(DES_EN_WALL_INPUT(M) == UNDEFINED) THEN
            WRITE (UNIT_LOG, 1010)
            CALL MFIX_EXIT(myPE)
         END IF
      ENDDO

      IF(DES_ETAT_FAC == UNDEFINED) THEN
         WRITE (UNIT_LOG, 1011)
      ELSEIF(DES_ETAT_FAC > ONE .OR. DES_ETAT_FAC < ZERO) THEN
         WRITE (UNIT_LOG, 1016)
         CALL MFIX_EXIT(myPE)
      END IF

      IF(DES_ETAT_W_FAC == UNDEFINED) THEN
         WRITE (UNIT_LOG, 1012)
      ELSEIF(DES_ETAT_W_FAC > ONE .OR. DES_ETAT_W_FAC < ZERO) THEN
         WRITE (UNIT_LOG, 1016)
         CALL MFIX_EXIT(myPE)
      END IF

      DO M = 1, MMAX+MMAX*(MMAX-1)/2
         IF(DES_EN_INPUT(M) > ONE .OR. DES_EN_INPUT(M) < ZERO) THEN
            WRITE (UNIT_LOG, 1013)
            CALL MFIX_EXIT(myPE)
         END IF
      ENDDO

      DO M = 1, MMAX
         IF(DES_EN_WALL_INPUT(M) > ONE .OR. DES_EN_WALL_INPUT(M) < ZERO) THEN
            WRITE (UNIT_LOG, 1014)
            CALL MFIX_EXIT(myPE)
         END IF
      ENDDO

      IF(MEW > ONE .OR. MEW_W > ONE .OR. MEW < ZERO .OR. MEW_W < ZERO) THEN
         WRITE (UNIT_LOG, 1015)
         CALL MFIX_EXIT(myPE)
      END IF

! Overwrite user's input in case of DEM (no fluid)
      IF(.NOT.DES_CONTINUUM_COUPLED) DES_INTERP_ON = .FALSE.
      
      RETURN
 1000 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'DES should only be run using CARTESIAN coordinates',/1X,70('*')/)
 1001 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'DES being run on multiple processors. Only serial runs allowed',/1X,70('*')/)
 1002 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'Direction of periodicity not defined in mfix.dat',/1X,70('*')/)
 1003 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA.',/' Message: ',&
         'INLET_OUTLET option is disabled: does not work properly in this MFIX version',/1X,70('*')/)
 1004 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'Spring constants KN or KN_W not specified in mfix.dat',/1X,70('*')/)
 1005 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'WARNING: Tangential spring factors KT_FAC or KT_W_FAC not specified in mfix.dat'/,&
         ' These factors will be defined in cfassign.f as 2/7 based on:'/,&
         ' Silbert et al, 2001, Physical Review E, vol. 64-5, see page 051302-5',/1X,70('*')/)
 1006 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'Damping constants MUST NOT be specified in mfix.dat,'/10X,&
         'Specify restitution coefficients instead',/1X,70('*')/)
 1007 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'Friction coefficients MEW or MEW_W not specified in mfix.dat',/1X,70('*')/)
 1008 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'Particle-particle restitution coefficient DES_EN_INPUT(M)',/10X,&
         'Must be specified in mfix.dat for interactions M = 1 to ',I5,&
          /1X,70('*')/)
 1010 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'Particle-wall restitution coefficients DES_EN_WALL_INPUT(M),'/10X,&
         'Must be specified in mfix.dat for interactions M = 1 to MMAX',/1X,70('*')/)
 1011 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'WARNING: Tangential damping factors DES_ETAT_FAC not specified in mfix.dat'/,&
         ' This factor will be set in cfassign.f as 1/2 based on:'/,&
         ' Silbert et al, 2003, Physics of Fluids, vol. 15-1, see page 3',/1X,70('*')/)
 1012 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'WARNING: Tangential damping factors DES_ETAT_W_FAC not specified in mfix.dat'/,&
         ' This factors will be set in cfassign.f as 1/2 based on:'/,&
         ' Silbert et al, 2003, Physics of Fluids, vol. 15-1, see page 3',/1X,70('*')/)
 1013 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'Unphysical ( > 1 or < 0) values of DES_EN_INPUT(M)',/1X,70('*')/)
 1014 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'Unphysical ( > 1 or < 0) values of DES_EN_WALL_INPUT(M)',/1X,70('*')/)
 1015 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'Unphysical ( > 1 or < 0) values of friction coefficients',/1X,70('*')/)
 1016 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'Values of DES_ETAT_FAC or DES_ETAT_W_FAC unphysical (< 0 or > 1) defined in mfix.dat',/1X,70('*')/)
 1017 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'Values of KT_FAC or KT_W_FAC unphysical (< 0 or > 1) defined in mfix.dat',/1X,70('*')/)
 1018 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message:',&
         ' WARNING: nsquare neighbor search may be slow with periodic',/&
         ' boundaries.  Grid based search is recommended.',/1X,70('*')/)
 1019 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message:',&
         ' octree and quadtree neighbor search methods do not',/&
         ' currently work with periodic boundaries.  Change',&
         ' value of',/&
         ' des_neighbor_search in mfix.dat',/1X,70('*')/)

      END SUBROUTINE CHECK_DES_DATA
