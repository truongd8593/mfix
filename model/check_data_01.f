!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CHECK_DATA_01                                           C
!  Purpose: check the run control namelist section                     C
!                                                                      C
!  Author: P. Nicoletti                               Date: 27-NOV-91  C
!  Reviewer: M.SYAMLAL, W.ROGERS, P.NICOLETTI         Date: 24-JAN-92  C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: allow Si unit                                              C
!  Author: S. Dartevelle                              Date: 01-Jul-02  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CHECK_DATA_01 

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param 
      USE param1 
      USE constant 
      USE run
      USE physprop
      USE indices
      USE scalars
      USE funits
      USE rxns
      USE cutcell, only : cartesian_grid
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER :: M, N, L
      Character*80  Line(1)
      CHARACTER*85 LONG_STRING
!-----------------------------------------------


      IF (DESCRIPTION == UNDEFINED_C) DESCRIPTION = ' ' 

      IF ( (UNITS /= 'CGS') .AND. (UNITS /= 'SI') )&
         CALL ERROR_ROUTINE ('check_data_01', &
         'UNITS not specified or illegal in mfix.dat',1,1) 

      IF (DT /= UNDEFINED) THEN 
         IF (TIME==UNDEFINED .OR. TIME<ZERO) &
            CALL ERROR_ROUTINE ('check_data_01', &
            'TIME not specified OR < 0.0 in mfix.dat',1,1) 
      ELSE 
         TIME = ZERO 
      ENDIF 

      IF (DT /= UNDEFINED) THEN 
         IF (TSTOP==UNDEFINED .OR. TSTOP<TIME) &
            CALL ERROR_ROUTINE ('check_data_01',&
            'TSTOP not specified OR TSTOP < TIME in mfix.dat',1,1) 
      ENDIF 

      IF (DT==UNDEFINED .OR. DT<ZERO) THEN 
         ODT = ZERO 
      ELSE 
         ODT = ONE/DT 
      ENDIF 

      IF(.NOT.(RUN_TYPE=='NEW' .OR. RUN_TYPE=='RESTART_1' &
         .OR. RUN_TYPE=='RESTART_2') )&
         CALL ERROR_ROUTINE ('check_data_01', &
         'RUN_TYPE not specified or illegal in mfix.dat', 1, 1) 

! Solids phase quantities
      IF (FRICTION) THEN 
! Sof: Schaeffer formulation is not used when friction is used 
         SCHAEFFER = .FALSE.    
         IF (.NOT.GRANULAR_ENERGY) &
            CALL ERROR_ROUTINE ('check_data_01', &
            'For FRICTION=.T., GRANULAR_ENERGY must be turned on',1,1)
         IF (SAVAGE>2 .OR. SAVAGE<0) &
            CALL ERROR_ROUTINE ('check_data_01', &
            'Value of SAVAGE should be 0, 1, or 2',1,1) 
         IF (BLENDING_STRESS) &
            CALL ERROR_ROUTINE ('check_data_01', &
            'Use blending_stress with Schaeffer not with Friction',1,1)
      ENDIF

      IF (JENKINS .AND. .NOT.GRANULAR_ENERGY) &
         CALL ERROR_ROUTINE ('check_data_01', &
         'GRANULAR_ENERGY must = .T. if JENKINS = .T.', 1, 1)

! GHD Theory requires some do loops over MMAX - 1, which is the real
! number of solids phases
      SMAX = MMAX
      IF(TRIM(KT_TYPE) == 'GHD') SMAX = MMAX - 1
      IF(TRIM(KT_TYPE) == 'GHD' .AND. SMAX > 2) THEN
        IF(DMP_LOG)WRITE (UNIT_LOG, 1504)
      ENDIF

! Sof: cannot use both Ahmadi and Simonin models at the same time.
      IF (AHMADI .AND. SIMONIN) &
         CALL ERROR_ROUTINE ('check_data_01', &
         'Cannot set both AHMADI = .T. and SIMONIN = .T.', 1, 1)  

      IF (.NOT.GRANULAR_ENERGY .AND. (AHMADI .OR. SIMONIN)) &
         CALL ERROR_ROUTINE ('check_data_01', &
         'GRANULAR_ENERGY must = .T. if AHMADI OR SIMONIN = .T.', 1, 1)

      IF (.NOT.K_EPSILON .AND. (AHMADI .OR. SIMONIN)) &
         CALL ERROR_ROUTINE ('check_data_01', &
         'K_EPSILON must = .T. if AHMADI OR SIMONIN = .T.', 1, 1)

! Sof: cannot use both L_scale0 and K-epsilon models at the same time.
      IF (K_Epsilon .AND. L_SCALE0/=ZERO) &
         CALL ERROR_ROUTINE ('check_data_01', &
         'Cannot set both K_Epsilon = .T. and L_SCALE0 /= ZERO',1,1)


! Sebastien May 2013, added SUBGRID conditions to be used in MFIX
! SUBGRID/LES model has some conditions with it
      IF (SUBGRID_TYPE /= UNDEFINED_C) THEN   
         IF(SUBGRID_TYPE /= 'IGCI' .AND. SUBGRID_TYPE /= 'MILIOLI') &
            CALL ERROR_ROUTINE ('check_data_01', &
            'The only option for SUBGRID_TYPE is IGCI or MILIOLI',1,1)

         LONG_STRING = 'GRANULAR_ENERGY and K_EPSILON &
            &must be FALSE when invoking a SUBGRID model'
         IF (GRANULAR_ENERGY .OR. K_Epsilon) &
            CALL ERROR_ROUTINE ('check_data_01', TRIM(LONG_STRING),1,1)

         LONG_STRING = 'Currently only WEN_YU DRAG_TYPE allowed &
            &when invoking a SUBGRID model'
         IF (DRAG_TYPE /= 'WEN_YU') &
            CALL ERROR_ROUTINE ('check_data_01', TRIM(LONG_STRING),1,1)

         LONG_STRING = 'FILTER_SIZE_RATIO unphysical'
         IF (filter_size_ratio .LE. ZERO) &
            CALL ERROR_ROUTINE ('check_data_01',TRIM(LONG_STRING),1,1)

         LONG_STRING = 'BLENDING_STRESS must be FALSE when invoking a &
            &SUBGRID model'
         IF (BLENDING_STRESS) CALL ERROR_ROUTINE ('check_data_01', &
            TRIM(LONG_STRING),1,1)

         LONG_STRING = 'FRICTION must be FALSE when invoking a &
            &SUBGRID model'
         IF (FRICTION) CALL ERROR_ROUTINE ('check_data_01', &
            TRIM(LONG_STRING),1,1)

         LONG_STRING = 'SHEAR must be FALSE when invoking a &
            &SUBGRID model'           
         IF (SHEAR) CALL ERROR_ROUTINE ('check_data_01',&
            TRIM(LONG_STRING),1,1)

         LONG_STRING = 'PHIP must be ZERO when SUBGRID_WALL is TRUE'
         IF ( SUBGRID_Wall .AND. (PHIP .NE. 0.0d0) ) &     
            CALL ERROR_ROUTINE ('check_data_01',TRIM(LONG_STRING),1,1)

         LONG_STRING = 'IF SUBGRID_WALL, then CARTESIAN_GRID must&
            &also be TRUE'
         IF (SUBGRID_Wall .AND. .NOT.CARTESIAN_GRID) &
            CALL ERROR_ROUTINE ('check_data_01', TRIM(LONG_STRING),1,1)
      ELSE   ! subgrid_type is undefined
         LONG_STRING = 'SUBGRID_WALL must be FALSE if SUBGRID_TYPE is &
            &undefined'
         IF(SUBGRID_WALL) &
            CALL ERROR_ROUTINE ('check_data_01',TRIM(LONG_STRING),1,1)
      ENDIF


! Check that phase number where added mass applies is properly defined.
      IF (ADDED_MASS) THEN
         LONG_STRING = 'Must set disperse phase number M_AM where &
            &virtual mass applies.'
         IF(M_AM == UNDEFINED_I) &
            CALL ERROR_ROUTINE ('check_data_01', TRIM(LONG_STRING),1,1)
         IF(M_AM == 0 .OR. M_AM > MMAX) &
            CALL ERROR_ROUTINE ('check_data_01', &
            'M_AM must be > 0 and <= MMAX', 1, 1)

         LONG_STRING = 'Added mass force cannot be implemented with &
            &GHD theory that solves for mixture eqs.'
         IF(KT_TYPE == 'GHD') &
            CALL ERROR_ROUTINE ('check_data_01', TRIM(LONG_STRING),1,1)
      ENDIF

! check for valid options for KT type
      IF (KT_TYPE /= UNDEFINED_C) THEN
         LONG_STRING = 'The only options for KT_TYPE are &
            &IA_NONEP, GD_99, or GHD'
         IF(KT_TYPE /= 'IA_NONEP' .AND. KT_TYPE /= 'GD_99' .AND. &
            KT_TYPE /= 'GHD') &
            CALL ERROR_ROUTINE ('check_data_01',TRIM(LONG_STRING),1,1)
         IF(KT_TYPE == 'GHD') THEN
            LONG_STRING = 'When KT_TYPE is GHD the only valid &
               &options for DRAG_TYPE is WEN_YU or HYS'
            IF(DRAG_TYPE /= 'WEN_YU' .AND. DRAG_TYPE /= 'HYS') &
               CALL ERROR_ROUTINE ('check_data_01',TRIM(LONG_STRING),1,1)
         ENDIF
      ENDIF

! sof: Check name of radial distribution function
      IF (RDF_TYPE /= 'LEBOWITZ') THEN
         IF(RDF_TYPE /= 'MODIFIED_LEBOWITZ' .AND. &
            RDF_TYPE /= 'MANSOORI' .AND. &
            RDF_TYPE /= 'MODIFIED_MANSOORI') &
            CALL ERROR_ROUTINE ('check_data_01','Unknown RDF_TYPE',1,1)
      ENDIF

! Set variable ANY_SPECIES_EQ
! modify to only check for real number of solid phases in case KT_TYPE='GHD'
      ANY_SPECIES_EQ = .FALSE. 
      M = 0 
      IF (SMAX + 1 > 0) THEN 
         ANY_SPECIES_EQ = ANY_SPECIES_EQ .OR. ANY(SPECIES_EQ(:SMAX)) 
         M = SMAX + 1 
      ENDIF 
      
! Check phase specification for Scalars
      DO N = 1, NScalar
         IF(Phase4Scalar(N) < 0 .OR. Phase4Scalar(N) > MMAX) THEN
            WRITE (Line,'(A, I3, A, I4, A)')&
               'Phase4Scalar( ', N, ') = ', Phase4Scalar(N), &
               ' is invalid'
            CALL ERROR_ROUTINE ('check_data_01', Line, 1, 1) 
         ENDIF
      ENDDO
      

! CHECK THE DISCRETIZATION METHODS
      IF (FPFOI) THEN
         DO L = 1,8
            IF(DISCRETIZE(L).LE.1) DISCRETIZE(L) = 2
         ENDDO
         IF(DMP_LOG)WRITE (UNIT_LOG, 1502) 
      ENDIF
      
      IF(Chi_scheme)THEN
        IF(DISCRETIZE(7).NE.3 .and. DISCRETIZE(7).NE.6) THEN ! works with smart & muscl
           IF(DMP_LOG)WRITE (UNIT_LOG, 1503)
           CALL Mfix_exit(0) 
        ENDIF
      ENDIF

! k4phi, phip0 for variable specularity coefficient      
      IF(BC_JJ_M) THEN
         IF(phi_w .eq. UNDEFINED) &
            CALL ERROR_ROUTINE ('check_data_01',&
            'Need to specify phi_w when BC_JJ_M is TRUE',1,1)
         WRITE (UNIT_LOG, 1505) e_w
         WRITE (UNIT_LOG, 1506) tan(phi_w*Pi/180.d0)

         k4phi = 7.d0/2.d0*tan(phi_w*Pi/180.d0)*(1.d0+e_w)
         IF (phip0 .eq. UNDEFINED) THEN
            phip0 = -0.0012596340709032689 + &
               0.10645510095633175*k4phi - &
               0.04281476447854031*k4phi**2 + &
               0.009759402181229842*k4phi**3 - &
               0.0012508257938705263*k4phi**4 + &
               0.00008369829630479206*k4phi**5 - &
               0.000002269550565981776*k4phi**6
! if k4phi is less than 0.2, the analytical expression for phi is used
! to estimate the phi at r->0
            IF (k4phi .le. 0.2d0) THEN
               phip0=0.09094568176225006*k4phi
            ENDIF
            WRITE (UNIT_LOG, 1507) phip0
         ENDIF
 
         IF (phip0 < 0) THEN
            CALL ERROR_ROUTINE ('check_data_01','phip0 less than zero',1,1)
         ENDIF

         IF (PHIP_OUT_JJ) THEN
            IF (nRR < 1) CALL ERROR_ROUTINE ('check_data_01',&
                'nRR should be at least 1 for storing specularity',1,1)
            WRITE(UNIT_LOG, 1508)
         ENDIF

      ENDIF
      

      RETURN  
 1500 FORMAT(/1X,70('*')//' From: CHECK_DATA_01',/' Message: ',&
         'Continuity equations can be discretized only with FOUP.',/&
         ' DISCRETIZE (1 or 2 or both) reset to 0.',G12.5,/1X,70('*')/) 
 1501 FORMAT(/1X,70('*')//' From: CHECK_DATA_01',/' Message: ',&
         'Only deferred correction method can be used.',/&
         ' DEF_COR reset to .true.',G12.5,/1X,70('*')/) 
 1502 FORMAT(/1X,70('*')//' From: CHECK_DATA_01',/' Message: ',&
         'Discretize>=2 has to be used along with',/&
         ' fourth order scheme',G12.5,/1X,70('*')/) 
 1503 FORMAT(/1X,70('*')//' From: CHECK_DATA_01',/' Message: ',&
         'Only schemes 3 and 5 are Chi-scheme enabled for',/&
         'species equations [DISCRETIZE(7)].',/1X,70('*')/) 
 1504 FORMAT(/1X,70('*')//' From: CHECK_DATA_01',/' Message: ',&
         'GHD theory may not be valid for more than two solids phases',/&
         ' it requires further development.',/1X,70('*')/) 
 1505 FORMAT(/1X,70('*')//' From: CHECK_DATA_01',/' Message: ',&
         'BC_JJ_M is TRUE, particle-wall restitution coefficient is' &
          ,G12.5,/1X,70('*')/) 
 1506 FORMAT(/1X,70('*')//' From: CHECK_DATA_01',/' Message: ',&
         'BC_JJ_M is TRUE, particle-wall friction coefficient is' &
          ,G12.5,/1X,70('*')/)           
 1507 FORMAT(/1X,70('*')//' From: CHECK_DATA_01',/' Message: ',&
         'No input for phip0 available, working expression is used' &
          ,G12.5,/1X,70('*')/)          
 1508 FORMAT(/1X,70('*')//' From: CHECK_DATA_01',/' Message: ',&
         'Specularity will be stored as the first element of ',/&
         ' ReactionRates array in WALL CELLs. Please avoid ', &
         'overwriting it',/1X,'when reacting flow is simulated',&
         /1X,70('*')/)

      END SUBROUTINE CHECK_DATA_01 
