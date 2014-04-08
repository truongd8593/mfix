!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_CONTINUUM_SOLIDS                                  !
!  Purpose: Check kinetic the run control namelist section             !
!                                                                      !
!  Author: P. Nicoletti                               Date: 27-NOV-91  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_SOLIDS_CONTINUUM

! Global Variables:
!---------------------------------------------------------------------//
      USE constant 
      USE run
      USE physprop
      USE indices
      USE scalars
      USE funits
      USE rxns
      USE cutcell, only : cartesian_grid

! Global Parameters:
!---------------------------------------------------------------------//
      USE param 
      USE param1 

! Global Module procedures:
!---------------------------------------------------------------------//
      use error_manager

      IMPLICIT NONE

! Local Variables:
!---------------------------------------------------------------------//
      INTEGER :: M, N, L, LC


!......................................................................!


! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_SOLIDS_CONTINUUM")

! Check EP_star. This is used to populate ep_star_array which is what
! should be used elsewhere in the code. (see set_constprop)
      IF(EP_STAR == UNDEFINED) THEN
         WRITE(ERR_MSG,1000) 'EP_STAR'
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ELSEIF(EP_STAR < ZERO .OR. EP_STAR > ONE) THEN
         WRITE(ERR_MSG, 1001)'EP_STAR', iVal(EP_STAR)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

! CHECK MU_s0
      IF (MU_S0 < ZERO) THEN 
         WRITE(ERR_MSG, 1001) 'MU_s0', MU_s0
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF 

! CHECK DIF_s0
      IF (DIF_S0 < ZERO) THEN 
         WRITE(ERR_MSG, 1001) 'DIF_s0', DIF_s0
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF 

! Check kinetic theory specifications.
      CALL CHECK_KT_TYPE

! Fedors_Landel correlation checks.
      IF(FEDORS_LANDEL .AND. SMAX /= 2) THEN
         WRITE(ERR_MSG, 1200)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF
 1200 FORMAT('Error 1200: FEDORS_LANDEL correlation is for binary ',   &
         'mixtures (MMAX=2).',/'Please correct the mfix.dat file.')

! Yu_Standish correlation checks
      IF(YU_STANDISH .AND. SMAX < 2) THEN
         WRITE(ERR_MSG, 1201)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF
 1201 FORMAT('Error 1201: YU_STANDISH correlation is for polydisperse',&
         ' mixtures',/'(MMAX >= 2). Please correct the mfix.dat file.')

      IF(YU_STANDISH .AND. FEDORS_LANDEL) THEN
         WRITE(ERR_MSG, 1202)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF
 1202 FORMAT('Error 1202: FEDORS_LANDEL and YU_STANDISH correlations ',&
         'cannot be',/'used at the same time. Please correct the ',    &
         'mfix.dat file.')

! If frictional stress modeling check various dependent/conflicting 
! settings
      IF (FRICTION) THEN 
! Turn off the default when friction is set.
         SCHAEFFER = .FALSE.

! Check that the granular energy PDE is solved.
         IF (.NOT.GRANULAR_ENERGY) THEN
            WRITE(ERR_MSG,1203)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
 1203 FORMAT('Error 1203: The FRICTION solids stress model requires ', &
         'setting',/'GRANULAR_ENERGY=.TRUE. Please correct the ',      &
         'mfix.dat file.')

! Check the value specified for SAVAGE.
         IF(SAVAGE>2 .OR. SAVAGE<0) THEN
            WRITE(ERR_MSG, 1001)'SAVAGE', iVal(SAVAGE)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

! Verify that stress blending is not turned on.
         IF(BLENDING_STRESS) THEN
            WRITE(ERR_MSG, 1204)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
 1204 FORMAT('Error 1204: The BLENDING_STRESS is used with SCHAEFFER ',&
         'not FRICTION.',/'Please correct the mfix.dat file.')

      ENDIF


      IF (JENKINS .AND. .NOT.GRANULAR_ENERGY) THEN
         WRITE(ERR_MSG, 1205)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF
 1205 FORMAT('Error 1205: The JENKINS small frication boundary ',      &
         'condition is only',/'valid when solving GRANULAR_ENERGY',/    &
         'Please correct the mfix.dat file.')


! if constant solids viscosity is employed, then the calculation for
! many of the solids phase transport coefficients that are needed for
! the granular energy equation will have be skipped in calc_mu_s. so
! do not evaluate granular energy when the solids viscosity is set to
! a constant!
      IF (MU_S0 .NE. UNDEFINED .AND. GRANULAR_ENERGY) THEN
         WRITE(ERR_MSG,1206)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF
 1206 FORMAT('Error 1206: A constant viscosity is specified but the ', &
         'setting',/'GRANULAR_ENERGY=.TRUE. Please correct the ',      &
         'mfix.dat file.')                 


! Set the flags for blending stresses.
      IF(BLENDING_STRESS) THEN
! Turn off the default if SIGM_BLEND is set.
         IF(SIGM_BLEND)  TANH_BLEND = .FALSE.
      ELSE
         TANH_BLEND  = .FALSE.
         SIGM_BLEND  = .FALSE.
      ENDIF


! this version of the restitution coefficient is needed by most KT_TYPE
! models. it is also needed in default solids-solids drag model      
      IF (C_E == UNDEFINED .AND. (MU_S0 == UNDEFINED .OR. &
          (SMAX>=2 .AND. KT_TYPE .EQ. UNDEFINED_C))) THEN
         WRITE(ERR_MSG,1300) 
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF
 1300 FORMAT('Error 1300: Coefficient of restitution (C_E) not ',      &
         'specified.',/'Please correct the mfix.dat file.')


      IF(C_F == UNDEFINED .AND. SMAX>=2 .AND.                          &
         KT_TYPE .EQ. UNDEFINED_C) THEN
         WRITE(ERR_MSG, 1301)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF
 1301 FORMAT('Error 1301: Coefficient of friction (C_F) not ',         &
         'specified.',/'Please correct the mfix.dat file.')


      IF((FRICTION .OR. SCHAEFFER) .AND. (PHI == UNDEFINED)) THEN
         WRITE(ERR_MSG, 1302)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF
 1302 FORMAT('Error 1302: Angle of internal friction (PHI) not ',      &
         'specified.',/'Please correct the mfix.dat file.')


      IF((FRICTION .OR. JENKINS) .AND. (PHI_W == UNDEFINED)) THEN
         WRITE(ERR_MSG, 1303)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF
 1303 FORMAT('Error 1303: Angle of particle-wall friction (PHI_W) not',&
         ' specified.',/'Please correct the mfix.dat file.')


      IF(MODEL_B) THEN
         DO LC = 1, MMAX 
            IF(.NOT.CLOSE_PACKED(LC)) THEN 
               WRITE(ERR_MSG, 1304), LC
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF 
         ENDDO 
      ENDIF 
 1304 FORMAT('Error 1304: Solids phase ',I2,' is not CLOSE_PACKED.',/, &
         'All solids phases must be CLOSE_PACKED with MODEL_B=.TURE.',/ &
         'Please correct the mfix.dat file.')


! Check that phase number where added mass applies is properly defined.
      IF (ADDED_MASS) THEN

         IF(M_AM == UNDEFINED_I)THEN
            WRITE(ERR_MSG, 1305)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
 1305 FORMAT('Error 1305: Must specify a disperse phase, M_AM, where ',&
         'the',/'virtual mass applies (ADDED_MASS).',/'Please correct',&
         ' the mfix.dat file.')

         ELSEIF(M_AM == 0 .OR. M_AM > MMAX) THEN
            WRITE(ERR_MSG,1306)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
 1306 FORMAT('Error 1306: M_AM is out of range. [1,MMAX]',/'Please ',  &
         'correct the mfix.dat file.')
         ENDIF

      ENDIF

! Check name of radial distribution function
      IF (RDF_TYPE /= 'LEBOWITZ') THEN
         IF(RDF_TYPE /= 'MODIFIED_LEBOWITZ' .AND. &
            RDF_TYPE /= 'MANSOORI' .AND. &
            RDF_TYPE /= 'MODIFIED_MANSOORI') &
            CALL ERROR_ROUTINE ('CHECK_CONTINUUM_SOLIDS',&
               'Unknown RDF_TYPE',1,1)
      ENDIF



 1000 FORMAT('Error 1000: Required input not specified: ',A,/'Please ',&
            'correct the mfix.dat file.')

 1001 FORMAT('Error 1001: Illegal or unphysical input: ',A,' = ',A,/   &
         'Please correct the mfix.dat file.')


! Plastic/frictional regime stress (jenkins is small frictional bc)
      IF (FRICTION .OR. JENKINS .OR. SCHAEFFER) THEN
! PHI and PHI_W are given in degree but calculated in radian within 
! the fortran codes 
         TAN_PHI_W = TAN(PHI_W*PI/180.D0)   ! jenkins or friction
         SIN_PHI = SIN(PHI*PI/180.D0)   ! for friction/schaeffer
         SIN2_PHI = SIN_PHI*SIN_PHI    ! for schaeffer
         F_PHI = (3.0D0 - 2.0D0*SIN2_PHI)/3.0D0    ! for commented
                                                   ! schaeffer section
      ELSE
         TAN_PHI_W = UNDEFINED
         SIN_PHI = UNDEFINED
         SIN2_PHI = UNDEFINED
         F_PHI = UNDEFINED
      ENDIF


! This variable is only used for GHD at this point...
! Define restitution coefficient matrix
      DO LC = 1, SMAX 
         DO N = 1, SMAX
            IF(r_p(LC,N) == UNDEFINED) r_p(LC,N) = C_e
! just need to define r_p(1,2) and r_p(2,1) will be set.
            r_p(N,LC) = r_p(LC,N)
         ENDDO
      ENDDO


! These quantities will only be used if invoking jj_bc_PS which implies
! granular_energy is true.
! k4phi, phip0 for variable specularity coefficient
      k4phi = UNDEFINED
      IF(BC_JJ_M) THEN
         IF(phi_w .eq. UNDEFINED) &
            CALL ERROR_ROUTINE ('CHECK_CONTINUUM_SOLIDS',&
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
            CALL ERROR_ROUTINE ('CHECK_CONTINUUM_SOLIDS',&
               'phip0 less than zero',1,1)
         ENDIF

         IF (PHIP_OUT_JJ) THEN
            IF (nRR < 1) CALL ERROR_ROUTINE ('CHECK_CONTINUUM_SOLIDS',&
                'nRR should be at least 1 for storing specularity',1,1)
            WRITE(UNIT_LOG, 1508)
         ENDIF

      ENDIF


      CALL FINL_ERR_MSG

      

      RETURN  
 1505 FORMAT(/1X,70('*')//' From: CHECK_CONTINUUM_SOLIDS',/' Message: ',&
         'BC_JJ_M is TRUE, particle-wall restitution coefficient is' &
          ,G12.5,/1X,70('*')/) 
 1506 FORMAT(/1X,70('*')//' From: CHECK_CONTINUUM_SOLIDS',/' Message: ',&
         'BC_JJ_M is TRUE, particle-wall friction coefficient is' &
          ,G12.5,/1X,70('*')/)           
 1507 FORMAT(/1X,70('*')//' From: CHECK_CONTINUUM_SOLIDS',/' Message: ',&
         'No input for phip0 available, working expression is used' &
          ,G12.5,/1X,70('*')/)          
 1508 FORMAT(/1X,70('*')//' From: CHECK_CONTINUUM_SOLIDS',/' Message: ',&
         'Specularity will be stored as the first element of ',/&
         ' ReactionRates array in WALL CELLs. Please avoid ', &
         'overwriting it',/1X,'when reacting flow is simulated',&
         /1X,70('*')/)

      END SUBROUTINE CHECK_SOLIDS_CONTINUUM



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_KT_TYPE                                           !
!  Purpose: Check kinetic theory input specifications. These checks    !
!  are almost all related to the KT_TYPE keyword.                      !
!                                                                      !
!  Author: J.Musser                                   Date: 04-FEB-14  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_KT_TYPE


! Global Variables:
!---------------------------------------------------------------------//
      USE constant 
      USE run
      USE physprop
      USE indices
      USE scalars
      USE funits
      USE rxns
      USE cutcell, only : cartesian_grid

! Global Parameters:
!---------------------------------------------------------------------//
      USE param 
      USE param1 

! Global Module procedures:
!---------------------------------------------------------------------//
      use error_manager

      IMPLICIT NONE

! Local Variables:
!---------------------------------------------------------------------//
! NONE

!......................................................................!


! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_KT_TYPE")


! These are some checks to satisfy legacy input:
      IF (AHMADI .AND. SIMONIN) THEN
         WRITE(ERR_MSG, 9001)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

 9001 FORMAT('Error 9001: Cannot specify AHMADI and SIMONIN together.',&
         /'Please correct the mfix.dat file.')

      ELSEIF(AHMADI) THEN
         IF(KT_TYPE(1:6) /= 'AHMADI' .AND.                             &
            KT_TYPE(1:8) /= 'LUN_1984')THEN
            WRITE(ERR_MSG,9002)trim(KT_TYPE)
            CALL FLUSH_ERR_MSG(ABORT = .TRUE.)

 9002 FORMAT('Error 9002: Cannot specify AHMADI and KT_TYPE = ',A,'.', &
         /'Please correct the mfix.dat file.')

         ELSE
            KT_TYPE='AHMADI'
         ENDIF

      ELSEIF(SIMONIN) THEN
         IF(KT_TYPE(1:7) /= 'SIMONIN' .AND.                            &
            KT_TYPE(1:8) /= 'LUN_1984')THEN
            WRITE(ERR_MSG,9003)trim(KT_TYPE)
            CALL FLUSH_ERR_MSG(ABORT = .TRUE.)

 9003 FORMAT('Error 9003: Cannot specify SIMONIN and KT_TYPE = ',A,'.',&
         /'Please correct the mfix.dat file.')

         ELSE
            KT_TYPE='SIMONIN'
         ENDIF
      ENDIF



! Check for valid options for kinetic theory models (KT_TYPE)
      SELECT CASE(trim(adjustl(KT_TYPE)))

!``````````````````````````````````````````````````````````````````````
      CASE ('IA_NONEP')
      CASE ('GD_99')
      CASE ('GTSH')


!``````````````````````````````````````````````````````````````````````
      CASE ('GHD')

         IF(DRAG_TYPE_ENUM /= WEN_YU .AND. DRAG_TYPE_ENUM /= HYS) THEN
            WRITE(ERR_MSG, 1101)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

         ELSEIF(ADDED_MASS) THEN
            WRITE(ERR_MSG,1102)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

         IF(SMAX > 2) THEN
            WRITE(ERR_MSG, 1103)
            CALL FLUSH_ERR_MSG
         ENDIF

! Automatically set SPECIES_EQ(MMAX) = .FALSE. to avoid potential
! checks/loops over the mmax species type eqn which has no meaning         
         SPECIES_EQ(MMAX) = .FALSE.
         NMAX_s(MMAX) = 1

! currently set to avoid an overflow error in write_res0 
! legacy variable?         
         NMAX(MMAX) = 1

 1101 FORMAT('Error 1101: KT_TYPE = "GHD" is restricted to DRAG_TYPE', &
         'values of WEN_YU and HYS.',/'Please correct the mfix.dat ',  &
         'file.')

 1102 FORMAT('Error 1102: ADDED_MASS force cannot be applied with ',   &
         'GHD theory that',/'solves for mixture equations.',/'Please', &
         'correct the mifx.dat file.')

 1103 FORMAT('Warning 1103: GHD theory may not be valid for more ',    &
         'than two solids phases',/'it requires further development.') 


!``````````````````````````````````````````````````````````````````````
      CASE ('AHMADI')
         KT_TYPE = UNDEFINED_C
         AHMADI = .TRUE.

         IF(.NOT.GRANULAR_ENERGY)THEN
            WRITE(ERR_MSG,1111) 'GRANULAR_ENERGY = .TRUE.'
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

         ELSEIF(.NOT.K_EPSILON) THEN
            WRITE(ERR_MSG,1111) 'K_EPSILON = .TRUE.'
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

 1111 FORMAT('Error 1111: KT_TYPE = "AHMADI" requires ',A,/            &
         'Please correct the mfix.dat file.')  


!``````````````````````````````````````````````````````````````````````
      CASE ('SIMONIN')
         KT_TYPE = UNDEFINED_C
         SIMONIN = .TRUE.

         IF(.NOT.GRANULAR_ENERGY)THEN
            WRITE(ERR_MSG,1121) 'GRANULAR_ENERGY = .TRUE.'
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

         ELSEIF(.NOT.K_EPSILON) THEN
            WRITE(ERR_MSG,1121) 'K_EPSILON = .TRUE.'
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

 1121 FORMAT('Error 1121: KT_TYPE = "SIMONIN" requires ',A,/           &
         'Please correct the mfix.dat file.')  


! Lun is the defalut implementation.
!``````````````````````````````````````````````````````````````````````
      CASE ('LUN_1984')
         KT_TYPE = UNDEFINED_C

      CASE DEFAULT
         WRITE(ERR_MSG,1100) trim(adjustl(KT_TYPE))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
 1100 FORMAT('Error 1100: Invalid or unknown KT_TYPE: ',A,/            &
         'Please correct the mfix.dat file.')

      END SELECT

      RETURN
      END SUBROUTINE CHECK_KT_TYPE
