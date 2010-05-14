!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: DES_INIT_BC                                            !
!                                                                      !
!                                                                      !
!                                                                      !
!  Author: J.Musser                                   Date: 23-Nov-09  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE DES_INIT_BC

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
! Local variables
!-----------------------------------------------
      INTEGER BCV, BCV_I      ! BC loop counter
      INTEGER M, MM           ! Mass phase loop counter
      INTEGER HOLD, I         ! Dummy values
      INTEGER RANGE_TOP, RANGE_BOT ! Dummy values
      INTEGER PHASE_CNT        ! Number of solid phases at bc
      INTEGER PHASE_LIST(MMAX) ! List of phases used in current bc

! the number of particles injected in a solids time step
      DOUBLE PRECISION NPMpSEC(MMAX) ! For solid phase m
      DOUBLE PRECISION NPpSEC
      DOUBLE PRECISION NPpDT        ! Total for BC
      DOUBLE PRECISION SCALED_VAL
      DOUBLE PRECISION MAX_DIA ! Max diameter of incoming particles at bc

!-----------------------------------------------
!   External functions
!-----------------------------------------------
      LOGICAL , EXTERNAL :: COMPARE 
!-----------------------------------------------


      WRITE(*,'(3X,A)') '---------- START DES_INIT_BC ---------->'

! If one or more discrete mass inlets exist, calculate necessary
! data 
      IF(DES_BCMI/=0)THEN
         BCV_I = 1
         DO BCV = 1, DIMENSION_BC 
           IF((DES_BC_DEFINED(BCV)) .AND. &         
              (DES_BC_TYPE(BCV) == 'MASS_INFLOW'))THEN

               DES_BC_MI_ID(BCV_I) = BCV

! Initialize temp variables
!    Tthe number of mass phases at this inlet.  While a system may be
!    polydisperse, the inlet could consist of a single mass phase
               PHASE_CNT = 0
!    The mass phase indices of incoming particles at this inlet
               PHASE_LIST(:) = -1
!    The max diameter of incoming particles at this inlet
               MAX_DIA = ZERO

! Determine if the inlet is mono or polydisperse               
               DO M=1, MMAX
                  IF(DES_BC_ROP_s(BCV,M) /= UNDEFINED .AND. &
                     .NOT.COMPARE(DES_BC_ROP_s(BCV,M),ZERO)) THEN
                     PHASE_CNT = PHASE_CNT + 1
                     PHASE_LIST(PHASE_CNT) = M
                     MAX_DIA = MAX(MAX_DIA,D_P0(M))
                  ENDIF
               ENDDO

! Flag that the inlet is polydisperse if true
               IF(PHASE_CNT > 1)DES_BC_POLY(BCV_I) = .TRUE.

! Set the value of the boundary condtion offset value used in the
! placement of new particles.
               DES_BC_OFFSET(BCV_I) = MAX_DIA

! Initialize temp variables
!    Number of phase m particles at BCV injected per second
               NPMpSEC(:) = ZERO
!    Total number of particles at BCV injected per second
               NPpSEC = ZERO
!    Total number of partices injected per solids time step
               NPpDT = ZERO

! calculate the number of particles of mass phase M are injected per 
! second for each solid phase present at the boundary
               DO MM=1,PHASE_CNT
                  M = PHASE_LIST(MM)
                  NPMpSEC(M) = (DES_BC_VOLFLOW_s(BCV,M) / &
                     (PI/6.d0 * D_P0(M)**3))
                  WRITE(*,"(5X,A,I2,A,F9.5)") &
                     'NPMpSEC(',M,'): ',NPMpSEC(M)
! calculate the total number of particles per second at the inlet
                     NPpSEC = NPpSEC + NPMpSEC(M)
               ENDDO

               WRITE(*,"(5X,A,F12.6)") 'NPpSEC: ',NPpSEC

! For polydisperse inlets, construct the DES_POLY_LAYOUT array
               IF(DES_BC_POLY(BCV_I)) THEN
!                  HOLD = 1
                  RANGE_BOT = 1
                  DO MM=1,PHASE_CNT - 1
                     M = PHASE_LIST(MM)
                     SCALED_VAL = dble(NUMFRAC_LIMIT)*(NPMpSEC(M)/NPpSEC)
                     RANGE_TOP = FLOOR(SCALED_VAL) + (RANGE_BOT-1)
!                     DO I=HOLD,FLOOR(SCALED_VAL)+(HOLD-1)
!                        DES_POLY_LAYOUT(BCV_I,I) = M
!                     ENDDO
                      DES_BC_POLY_LAYOUT(BCV_I,RANGE_BOT:RANGE_TOP) = M 
                      RANGE_BOT = RANGE_TOP+1
!                      HOLD = I
                  ENDDO
                  M = PHASE_LIST(PHASE_CNT)
!                  DO I=HOLD,NUMFRAC_LIMIT
!                     DES_POLY_LAYOUT(BCV_I,I) = M
!                  ENDDO
                  DES_BC_POLY_LAYOUT(BCV_I,RANGE_BOT:NUMFRAC_LIMIT) = M
! For monodisperse inlets, store the single mass phase used
               ELSE
                  DES_BC_POLY_LAYOUT(BCV_I,:) = PHASE_LIST(1)
               ENDIF

! The number of total particles per solid time step DTSOLID
               NPpDT = NPpSEC * DTSOLID

               IF(NPpDT .LT. 1)THEN
! The number of solid time steps between injection of a particle
                  PI_FACTOR(BCV_I) = FLOOR(real(1.d0/NPpDT))
                  PI_COUNT(BCV_I) = 1
               ELSE
                  PI_FACTOR(BCV_I) = 1
! The number of particles injected in a solid time step                  
                  PI_COUNT(BCV_I) = CEILING(real(NPpDT))
               ENDIF

               WRITE(*,"(5X,A,F12.6)") 'NPpDT: ',NPpDT

! Calculate des mass inlet time; time between injection.  If the run
! type is RESTART_1, DES_MI_TIME will be picked up from the restart file
! with an updated value.
               IF(RUN_TYPE == 'NEW')THEN
                  DES_MI_TIME(BCV_I) = TIME +&
                     dble(PI_FACTOR(BCV_I)) * DTSOLID 
               ENDIF
               WRITE(*,1000) BCV, NPpDT, PI_FACTOR(BCV_I),&
                  PI_COUNT(BCV_I), DES_MI_TIME(BCV_I)

! Classify boundary condition and verify appropriate solids velocity
               CALL DES_MI_CLASSIFY(BCV, BCV_I, MAX_DIA, PHASE_CNT,&
                  PHASE_LIST)

! Determine the computational cells near the inlet
               CALL DES_MI_CELLS(BCV, BCV_I)

! Verify that the inlet is not 'too' close to a periodic boundary condition
               IF(DES_PERIODIC_WALLS)THEN
                  CALL DES_MIO_PERIODIC(BCV, BCV_I,'MI',MAX_DIA)
               ENDIF

               BCV_I = BCV_I + 1 
            ENDIF
         ENDDO
      ENDIF

! Check each discrete mass outlet for necessary data
      IF(DES_BCMO/=0)THEN
         BCV_I = 1
         DO BCV = 1, DIMENSION_BC 
           IF((DES_BC_DEFINED(BCV)) .AND. &
              (DES_BC_TYPE(BCV) == 'MASS_OUTFLOW'))THEN

               DES_BC_MO_ID(BCV_I) = BCV

               CALL DES_MO_CLASSIFY(BCV_I, BCV)

! Verify that the outlet is not too close to a periodic boundary condition
               IF(DES_PERIODIC_WALLS)THEN
                  CALL DES_MIO_PERIODIC(BCV,BCV_I,'MO',MAX_DIA)
               ENDIF

               BCV_I = BCV_I + 1 
            ENDIF
         ENDDO
      ENDIF

      WRITE(*,'(3X,A)') '<---------- END DES_INIT_BC ----------'

 1000 FORMAT(/5X,'For mass inlet BC: ', I3,/,&
         7X,'No. particles injected per solids time step = ', ES15.8,/,&
         7X,'PI_FACTOR = ', I10,' PI_COUNT = ', I5,/,&
         7X,'start DES_MI_TIME = ', ES15.8)

      RETURN
      END SUBROUTINE DES_INIT_BC


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

      SUBROUTINE DES_MI_CLASSIFY(BCV, BCV_I, MAX_DIA, PHASE_CNT, &
         PHASE_LIST)

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
      INTEGER M, MM
! passed arguments giving index information of the boundary
      INTEGER BCV_I, BCV
! Number of solid phases at bc
      INTEGER PHASE_CNT
! List of phases used in current bc
      INTEGER PHASE_LIST(MMAX)
! the length of each side of the inlet boundary
      DOUBLE PRECISION LEN1, LEN2

! Max diameter of incoming particles at bc
      DOUBLE PRECISION MAX_DIA
! temp value of DES_MI_CLASS for comparison tests
      CHARACTER*4 DMC
!-----------------------------------------------


      IF(DIMN == 2)THEN   ! 2D domain

! Check verticle mass inlet: 2D
! ----------------------------------------
         IF(DES_BC_X_w(BCV) == DES_BC_X_e(BCV))THEN

! Classify the boundary condition
            IF(DES_BC_X_w(BCV) == ZERO) DES_MI_CLASS(BCV_I) = 'XW'
            IF(DES_BC_X_w(BCV) == XLENGTH) DES_MI_CLASS(BCV_I) = 'XE'

! Inlet length
            LEN1 = ABS(DES_BC_Y_n(BCV) - DES_BC_Y_s(BCV))

         ENDIF

! Check horizontal mass inlet: 2D 
! ----------------------------------------
         IF(DES_BC_Y_s(BCV) == DES_BC_Y_n(BCV))THEN

! Classify the boundary condition
            IF(DES_BC_Y_s(BCV) == ZERO) DES_MI_CLASS(BCV_I) = 'YS'
            IF(DES_BC_Y_s(BCV) == YLENGTH) DES_MI_CLASS(BCV_I) = 'YN'

! Inlet length
            LEN1 = DES_BC_X_e(BCV) - DES_BC_X_w(BCV)
         ENDIF

! max_dia is used for the 'depth/width' of the inlet 
! (see des_bc_vel_assign for details)         
         CALL DES_BC_VEL_ASSIGN(BCV,BCV_I,LEN1,MAX_DIA,PHASE_CNT,PHASE_LIST)

! This subroutine determines the pattern that the particles will need to
! enter the system, if any. This routine only needs to be called if a
! run is new.  If a run is a RESTART_1, all of the setup information
! provided by this subroutine is will be obtained from the *_DES.RES file.
! This is done due to this routine's strong dependence on the 
! RANDOM_NUMBER() subroutine.
            IF(RUN_TYPE == 'NEW')THEN
               CALL DES_MI_LAYOUT(LEN1, ZERO, BCV, BCV_I, MAX_DIA,&
                  PHASE_CNT, PHASE_LIST)
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

! If a slit flow inlet is specified near an edge, expand the slit to 
! the width of one particle diameter :
! see comments in the yz face section for further explanation   
            IF(DES_BC_X_w(BCV) == DES_BC_X_e(BCV))THEN
               IF(((DES_BC_X_w(BCV) - MAX_DIA*HALF) .GE. ZERO) .AND. &
                  ((DES_BC_X_w(BCV) + MAX_DIA*HALF) .LE. XLENGTH))THEN
                  DES_BC_X_w(BCV) = DES_BC_X_w(BCV) - MAX_DIA*HALF
                  DES_BC_X_e(BCV) = DES_BC_X_e(BCV) + MAX_DIA*HALF
               ELSEIF((DES_BC_X_w(BCV) - MAX_DIA*HALF) .LT. ZERO)THEN
                  DES_BC_X_w(BCV) = ZERO
                  DES_BC_X_e(BCV) = MAX_DIA
               ELSEIF((DES_BC_X_e(BCV) + MAX_DIA*HALF) .GT. XLENGTH)THEN
                  DES_BC_X_w(BCV) = XLENGTH - MAX_DIA
                  DES_BC_X_e(BCV) = XLENGTH 
               ENDIF
            ENDIF
            IF(DES_BC_Z_b(BCV) == DES_BC_Z_t(BCV))THEN
               IF(((DES_BC_Z_b(BCV) - MAX_DIA*HALF) .GE. ZERO) .AND. &
                  ((DES_BC_Z_b(BCV) + MAX_DIA*HALF) .LE. ZLENGTH))THEN
                  DES_BC_Z_b(BCV) = DES_BC_Z_b(BCV) - MAX_DIA*HALF
                  DES_BC_Z_t(BCV) = DES_BC_Z_t(BCV) + MAX_DIA*HALF
               ELSEIF((DES_BC_Z_b(BCV) - MAX_DIA*HALF) .LT. ZERO)THEN
                  DES_BC_Z_b(BCV) = ZERO
                  DES_BC_Z_t(BCV) = MAX_DIA
               ELSEIF((DES_BC_Z_t(BCV) + MAX_DIA*HALF) .GT. ZLENGTH)THEN
                  DES_BC_Z_b(BCV) = ZLENGTH - MAX_DIA
                  DES_BC_Z_t(BCV) = ZLENGTH 
               ENDIF
            ENDIF

! Classify the boundary condition
            IF(DES_BC_Y_s(BCV) == ZERO) DES_MI_CLASS(BCV_I) = 'XZs'
            IF(DES_BC_Y_s(BCV) == YLENGTH) DES_MI_CLASS(BCV_I) = 'XZn'

! Inlet dimensions
            LEN1 = DES_BC_X_e(BCV) - DES_BC_X_w(BCV)
            LEN2 = DES_BC_Z_t(BCV) - DES_BC_Z_b(BCV)

         ENDIF ! End check mass inlet on XZ face: 3D 

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

! If a slit flow inlet is specified near an edge, expand the slit to 
! the width of one particle diameter :
! see comments in the yz face section for further explanation            
            IF(DES_BC_X_w(BCV) == DES_BC_X_e(BCV))THEN
               IF(((DES_BC_X_w(BCV) - MAX_DIA*HALF) .GE. ZERO) .AND. &
                  ((DES_BC_X_w(BCV) + MAX_DIA*HALF) .LE. XLENGTH))THEN
                  DES_BC_X_w(BCV) = DES_BC_X_w(BCV) - MAX_DIA*HALF
                  DES_BC_X_e(BCV) = DES_BC_X_e(BCV) + MAX_DIA*HALF
               ELSEIF((DES_BC_X_w(BCV) - MAX_DIA*HALF) .LT. ZERO)THEN
                  DES_BC_X_w(BCV) = ZERO
                  DES_BC_X_e(BCV) = MAX_DIA
               ELSEIF((DES_BC_X_e(BCV) + MAX_DIA*HALF) .GT. XLENGTH)THEN
                  DES_BC_X_w(BCV) = XLENGTH - MAX_DIA
                  DES_BC_X_e(BCV) = XLENGTH 
               ENDIF
            ENDIF
            IF(DES_BC_Y_s(BCV) == DES_BC_Y_n(BCV))THEN
               IF(((DES_BC_Y_s(BCV) - MAX_DIA*HALF) .GE. ZERO) .AND. &
                  ((DES_BC_Y_s(BCV) + MAX_DIA*HALF) .LE. YLENGTH))THEN
                  DES_BC_Y_s(BCV) = DES_BC_Y_s(BCV) - MAX_DIA*HALF
                  DES_BC_Y_n(BCV) = DES_BC_Y_n(BCV) + MAX_DIA*HALF
               ELSEIF((DES_BC_Y_s(BCV) - MAX_DIA*HALF) .LT. ZERO)THEN
                  DES_BC_Y_s(BCV) = ZERO
                  DES_BC_Y_n(BCV) = MAX_DIA
               ELSEIF((DES_BC_Y_n(BCV) + MAX_DIA*HALF) .GT. YLENGTH)THEN
                  DES_BC_Y_s(BCV) = YLENGTH - MAX_DIA
                  DES_BC_Y_n(BCV) = YLENGTH 
               ENDIF
            ENDIF

! Classify the boundary condition
            IF(DES_BC_Z_b(BCV) == ZERO) DES_MI_CLASS(BCV_I) = 'XYb'
            IF(DES_BC_Z_b(BCV) == ZLENGTH) DES_MI_CLASS(BCV_I) = 'XYt'

! Inlet dimensions
            LEN1 = (DES_BC_X_e(BCV) - DES_BC_X_w(BCV))
            LEN2 = (DES_BC_Y_n(BCV) - DES_BC_Y_s(BCV))

         ENDIF ! End check mass inlet on XY face: 3D 

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
               IF(((DES_BC_Y_s(BCV) - MAX_DIA*HALF) .GE. ZERO) .AND. &
                  ((DES_BC_Y_s(BCV) + MAX_DIA*HALF) .LE. YLENGTH))THEN
                  DES_BC_Y_s(BCV) = DES_BC_Y_s(BCV) - MAX_DIA*HALF
                  DES_BC_Y_n(BCV) = DES_BC_Y_n(BCV) + MAX_DIA*HALF
               ELSEIF((DES_BC_Y_s(BCV) - MAX_DIA*HALF) .LT. ZERO)THEN
                  DES_BC_Y_s(BCV) = ZERO
                  DES_BC_Y_n(BCV) = MAX_DIA 
               ELSEIF((DES_BC_Y_n(BCV) + MAX_DIA*HALF) .GT. YLENGTH)THEN
                  DES_BC_Y_s(BCV) = YLENGTH - MAX_DIA
                  DES_BC_Y_n(BCV) = YLENGTH 
               ENDIF
            ENDIF
            IF(DES_BC_Z_b(BCV) == DES_BC_Z_t(BCV))THEN
               IF(((DES_BC_Z_b(BCV) - MAX_DIA*HALF) .GE. ZERO) .AND. &
                  ((DES_BC_Z_b(BCV) + MAX_DIA*HALF) .LE. ZLENGTH))THEN
                  DES_BC_Z_b(BCV) = DES_BC_Z_b(BCV) - MAX_DIA*HALF
                  DES_BC_Z_t(BCV) = DES_BC_Z_t(BCV) + MAX_DIA*HALF
               ELSEIF((DES_BC_Z_b(BCV) - MAX_DIA*HALF) .LT. ZERO)THEN
                  DES_BC_Z_b(BCV) = ZERO
                  DES_BC_Z_t(BCV) = MAX_DIA
               ELSEIF((DES_BC_Z_t(BCV) + MAX_DIA*HALF) .GT. ZLENGTH)THEN
                  DES_BC_Z_b(BCV) = ZLENGTH - MAX_DIA
                  DES_BC_Z_t(BCV) = ZLENGTH 
               ENDIF
            ENDIF

! Classify the boundary condition
            IF(DES_BC_X_w(BCV) == ZERO) DES_MI_CLASS(BCV_I) = 'YZw'
            IF(DES_BC_X_w(BCV) == XLENGTH) DES_MI_CLASS(BCV_I) = 'YZe'

! Inlet dimensions
            LEN1 = (DES_BC_Y_n(BCV) - DES_BC_Y_s(BCV))
            LEN2 = (DES_BC_Z_t(BCV) - DES_BC_Z_b(BCV))

         ENDIF ! End check mass inlet on YZ face: 3D

         CALL DES_BC_VEL_ASSIGN(BCV,BCV_I,LEN1,LEN2,PHASE_CNT,PHASE_LIST)

! This subroutine determines the pattern that the particles will need to
! enter the system, if any. This routine only needs to be called if a
! run is new.  If a run is a RESTART_1, all of the setup information
! provided by this subroutine is will be obtained from the *_DES.RES file.
! This is done due to this routine's strong dependence on the 
! RANDOM_NUMBER() subroutine.
         IF(RUN_TYPE == 'NEW')THEN
            CALL DES_MI_LAYOUT(LEN1, LEN2, BCV, BCV_I, MAX_DIA,&
               PHASE_CNT, PHASE_LIST)
         ENDIF

      ENDIF   ! end if dimn == 2/else

! Verify that an inlet is not on a face that is connected to a periodic
! boundary condition.  If so, write error message and exit.
      DMC = DES_MI_CLASS(BCV_I)
      IF (DES_PERIODIC_WALLS) THEN
! No XW, XE, YZw or YZe inlet with X direction periodic walls
         IF((DMC == 'XW'  .OR. DMC == 'XE'  .OR. &
             DMC == 'YZw' .OR. DMC == 'YZe') .AND. DES_PERIODIC_WALLS_X) THEN
            WRITE(UNIT_LOG, 1202) BCV, 'DES_PERIODIC_WALLS_X'
            WRITE(*, 1202) BCV, 'DES_PERIODIC_WALLS_X'
            CALL MFIX_EXIT(myPE)
         ENDIF
! No YS, YN, XZs or XZn inlet with Y direction periodic walls
         IF((DMC == 'YS'  .OR. DMC == 'YN' .OR. &
             DMC == 'XZs' .OR. DMC == 'XZn') .AND. DES_PERIODIC_WALLS_Y) THEN
            WRITE(UNIT_LOG, 1202 )BCV, 'DES_PERIODIC_WALLS_Y'
            WRITE(*, 1202) BCV, 'DES_PERIODIC_WALLS_Y'
            CALL MFIX_EXIT(myPE)
         ENDIF
! No XYb or XYt inlet with Z direction periodic walls
         IF((DMC == 'XYb' .OR. DMC == 'XYt') .AND. DES_PERIODIC_WALLS_Z) THEN
            WRITE (UNIT_LOG, 1202) BCV, 'DES_PERIODIC_WALLS_Z'
            WRITE (*, 1202) BCV, 'DES_PERIODIC_WALLS_Z'
            CALL MFIX_EXIT(myPE)
         ENDIF
      ENDIF

! Set flags for individual mass inlets for used in the grid based search
! Flag that a X-face inlet exits
      IF(DMC == 'XW' .OR. DMC == 'XE' .OR. &
         DMC == 'YZw' .OR. DMC == 'YZe') THEN
         DES_MI_X = .TRUE.
      ENDIF
! Flag that a Y-face inlet exits
      IF(DMC == 'YS' .OR. DMC == 'YN' .OR. &
         DMC == 'XZs' .OR. DMC == 'XZn') THEN
         DES_MI_Y = .TRUE.
      ENDIF
! Flag that a Z-face inlet exits
      IF(DMC == 'XYb' .OR. DMC == 'XYt') THEN
         DES_MI_Z = .TRUE.
      ENDIF

 1202 FORMAT(/1X,70('*')//, ' From: DES_MI_CLASSIFY -',/,&
         ' Message : DEM inlet can not be placed on a periodic',&
         ' boundary.',/10X,'Check DEM boundary condtion ',I3,&
         ' and ',A,'.',/1X,70('*')/)

      RETURN
      END SUBROUTINE DES_MI_CLASSIFY

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: DES_BC_VEL_ASSIGN                                      !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!  Author: J.Musser                                   Date: 21-Jan-10  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE DES_BC_VEL_ASSIGN(BCV,BCV_I,LEN1,LEN2,PHASE_CNT,&
         PHASE_LIST)

      USE compar
      USE des_bc
      USE discretelement
      USE funits
      USE geometry
      USE param
      USE physprop

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! passed arguments giving index information of the boundary
      INTEGER BCV_I, BCV
! Number of solid phases at bc
      INTEGER PHASE_CNT
! List of phases used in current bc
      INTEGER PHASE_LIST(MMAX)
! the length of each side of the inlet boundary
      DOUBLE PRECISION LEN1, LEN2, BCV_AREA
! the number of trial velocities to test for UI_VEL in a given iteration
      INTEGER, PARAMETER :: NUM_VEL = 100
! Dummy index
      INTEGER I
! Index of solid phases at bc
      INTEGER M, MM
! Solids volume fraction of solids phase M
      DOUBLE PRECISION EP_sM
! Temp inlet velocity for solids phase M
      DOUBLE PRECISION TMP_VEL(DIM_M)
! Minimum/maximum solids velocity at inlet.  Also used in the iterative
! steps as the starting and ending velocities 
      DOUBLE PRECISION  MIN_VEL, MAX_VEL
! Temp variable for DES_MI_CLASS(BCV_I)
      CHARACTER*3 DMCL
! Uniform inlet velocity
      DOUBLE PRECISION UI_VEL
! Error in solids phase bulk density for phase M corresponding 
! to UI_VEL compared to the user set bulk density for phase M
      DOUBLE PRECISION ERR_ROPSM
! Bulk density of phase M calculated using STP_VEL or UI_VEL
      DOUBLE PRECISION CAL_ROPSM

!---  Iterative process variables
! Counts the number of iterations performed and exits if >100
      INTEGER LC_EXIT
! Velocity minimizing the error for the PReVious and CURrent iterations
      DOUBLE PRECISION CUR_VEL, PRV_VEL
! Incremental velocity step size (delta velocity)
      DOUBLE PRECISION DELTA_VEL
! Incremental step velocity
      DOUBLE PRECISION STP_VEL
! Total error for all solid phase generated by STP_VEL
      DOUBLE PRECISION ERR_TOT
! Minimum error generated on current iteration
      DOUBLE PRECISION ERR_MIN
! Convergence tolerance
      DOUBLE PRECISION, PARAMETER :: TOL = 1.0e-4
!-----------------------------------------------
! Initialize
      MIN_VEL = LARGE_NUMBER
      MAX_VEL = ZERO

! BC inlet area
! To calculate velocity based on specified mass/volumetric inflow and
! bulk density the inlet needs to have dimension of length^2 (area).  
! In the 2D case the depth of the inlet is taken to be the maximum
! particle diameter which is consistent with how the criteria for
! maximum and minimum velocity for the boundary are defined and how
! particles are seeded for a 2D inlet (particles are seeded along a 
! line that has a depth equal to the maximum particle diameter).  
! Alternatively, the  depth of the inlet could be taken as zlength, 
! but zlength is not guaranteed to be equal to the maximum particle 
! diameter.
      BCV_AREA = LEN1*LEN2

! Calculate the individual velocities for each solid phase
      DO MM = 1, PHASE_CNT
         M = PHASE_LIST(MM)
! Solids volume fraction for phase M
         EP_sM = DES_BC_ROP_s(BCV,M)/RO_s(M)
! Inlet velocity for solids phase M
         TMP_VEL(M) = ( DES_BC_VOLFLOW_s(BCV,M) / &
            BCV_AREA ) / EP_sM
! Check for min/max inlet velocity
         MIN_VEL = MIN(ABS(TMP_VEL(M)), MIN_VEL)
         MAX_VEL = MAX(ABS(TMP_VEL(M)), MAX_VEL)
      ENDDO

! Determine necessary uniform inlet velocity
      IF(MIN_VEL == MAX_VEL) THEN
         UI_VEL = MIN_VEL
      ELSE

! Initialize Values
         LC_EXIT = 0
         PRV_VEL = UNDEFINED
         CUR_VEL = ZERO

         DO WHILE (ABS(CUR_VEL - PRV_VEL)>TOL)
! Initialize Values
!   Store current velocity value in previous velocity
            PRV_VEL = CUR_VEL 
!   Set minimum error value to a large number
            ERR_MIN = LARGE_NUMBER
!   Determine step size between tested velocity values
            DELTA_VEL = (MAX_VEL-MIN_VEL)/DBLE(NUM_VEL)

! Check i=NUM_VEL velocites between MAX_VEL and MIN_VEL for one that
! minimizes the total error between calculated bulk density values 
! and the values provided in the mfix.dat file.  The first loop (i=1)
! acts to initialize the values of cur_vel and err_min. 
            DO I = 1, (NUM_VEL+1)
! Initialize values
               STP_VEL = MIN_VEL + DELTA_VEL*(I-1)
               ERR_TOT = ZERO
! Calculate total error over all solid phases with respect to STP_VEL
               DO MM = 1, PHASE_CNT
                  M = PHASE_LIST(MM)
! Calculate bulk density value based upon velocity 
                  CAL_ROPSM = (RO_s(M) * DES_BC_VOLFLOW_s(BCV,M)) / &
                     (STP_VEL * BCV_AREA)
                  ERR_TOT = ERR_TOT+ABS(DES_BC_ROP_s(BCV,M)-CAL_ROPSM)
               ENDDO
! Compare to determine if ERR_TOT is minimum error over all STP_VEL
               IF(ERR_TOT == MIN(ERR_TOT, ERR_MIN)) THEN
                  ERR_MIN = ERR_TOT
                  CUR_VEL = STP_VEL
               ENDIF
! Narrow search range for next iterative set.  These values only change
! if cur_vel is updated which only occurs when the current stp_vel gives
! a smaller total error than previously
               MIN_VEL = CUR_VEL - DELTA_VEL
               MAX_VEL = CUR_VEL + DELTA_VEL
            ENDDO

! Loop control to prevent hang up. * This should not be necessary *
            LC_EXIT = LC_EXIT +1
            IF(LC_EXIT > 100) THEN
               WRITE(*,1256) BCV
               EXIT
            ENDIF
         ENDDO
         UI_VEL = CUR_VEL
      ENDIF

! Assign the uniform inlet velocity to each solid phase
      DMCL = DES_MI_CLASS(BCV_I)
      DO MM = 1, PHASE_CNT
         M = PHASE_LIST(MM)
         IF(DMCL == 'XW' .OR. DMCL == 'YZw')THEN
            DES_BC_U_s(BCV) =  UI_VEL
         ELSEIF(DMCL == 'XE' .OR. DMCL == 'YZe')THEN
            DES_BC_U_s(BCV) = -UI_VEL
         ELSEIF(DMCL == 'YS' .OR. DMCL == 'XZs')THEN
            DES_BC_V_s(BCV) =  UI_VEL
         ELSEIF(DMCL == 'YN' .OR. DMCL == 'XZn')THEN
            DES_BC_V_s(BCV) = -UI_VEL
         ELSEIF(DMCL == 'XYb')THEN
            DES_BC_W_s(BCV) =  UI_VEL
         ELSEIF(DMCL == 'XYt')THEN
            DES_BC_W_s(BCV) = -UI_VEL
         ELSE
            WRITE(*,1257) BCV, DMCL 
            CALL MFIX_EXIT(myPE)
         ENDIF
      ENDDO

! Output the new data values (screen/log file)
      WRITE(*,1249)
! Table Header
      WRITE(*,1250) BCV
      WRITE(*,1251) UI_VEL
      WRITE(*,1252)
      WRITE(UNIT_LOG,1250)BCV
      WRITE(UNIT_LOG,1251)UI_VEL
      WRITE(UNIT_LOG,1252)

! Column Labels
      WRITE(*,1253)
      WRITE(*,1254)
      WRITE(*,1253)
      WRITE(UNIT_LOG,1253)
      WRITE(UNIT_LOG,1254)
      WRITE(UNIT_LOG,1253)

! Fill Table Rows
      DO MM = 1, PHASE_CNT
         M = PHASE_LIST(MM)
! Calculate bulk density value based upon velocity 
         CAL_ROPSM = (RO_s(M) * DES_BC_VOLFLOW_s(BCV,M)) / &
            (UI_VEL * BCV_AREA)
         ERR_ROPSM =  ABS(DES_BC_ROP_s(BCV,M) - CAL_ROPSM)
         WRITE(*,1255) M, DES_BC_ROP_s(BCV,M), CAL_ROPSM, ERR_ROPSM
         WRITE(*,1253)
         WRITE(UNIT_LOG,1255) M, DES_BC_ROP_s(BCV,M), CAL_ROPSM, &
            ERR_ROPSM
         WRITE(UNIT_LOG,1253)
      ENDDO
      WRITE(*,"(//)")
      WRITE(UNIT_LOG,"(//)")

 1249 FORMAT(//,5X,'From: DES_BC_VEL_ASSIGN - ')
 1250 FORMAT(5X,'|<--- Boundary Condition ',I2,1X,26('-'),'>|')
 1251 FORMAT(5X,'| Uniform Inlet Velocity: ',ES11.4,18(' '),'|')
 1252 FORMAT(5X,'| Adjusted DES_BC_ROP_s Values',25(' '),'|')
 1253 FORMAT(5X,'|',54('-'),'|')
 1254 FORMAT(5X,'|',3X,'Phase',3X,'|',2X,'Specified',2X,'|',2X,&
         'Calculated',2X,'|',2X,'ABS Error',2X,'|')
 1255 FORMAT(5X,'|',4X,I2,5X,'|',1X,ES11.4,1X,'|',1X,ES11.4,2X,'|',&
         1X,ES11.4,1X,'|')

 1256 FORMAT(/1X,70('*')//,' From: DES_BC_VEL_ASSIGN -',/,&
         ' Message : Tolerance not met on uniform inlet velocity for',&
         /10X,'boundary ', I3, /1X,70('*')/)

 1257 FORMAT(/1X,70('*')//,' From: DES_BC_VEL_ASSIGN -',/,&
         ' Message : INVALID BOUNDARY CLASSIFICATION FOR BOUNDARY',I3,&
         /10X, 'with classification ',A,/1X,70('*')/)
         

      RETURN
      END SUBROUTINE DES_BC_VEL_ASSIGN


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

      SUBROUTINE DES_MI_LAYOUT(LEN1, LEN2, BCV, BCV_I, MAX_DIA,&
                               PHASE_CNT, PHASE_LIST)

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
! passed arguments giving the length of the boundary edges
      DOUBLE PRECISION LEN1, LEN2
! passed arguments giving index information of the boundary
      INTEGER BCV_I, BCV 
! Max diameter of incoming particles at bc
      DOUBLE PRECISION MAX_DIA
! Number of solid phases at bc
      INTEGER PHASE_CNT
! List of phases used in current bc
      INTEGER PHASE_LIST(MMAX)

! max number of particle diameters that fit along length of inlet
      INTEGER TMP_LEN1, TMP_LEN2

      DOUBLE PRECISION MAX_ROPs
! tmp variable to store the inlet velocity at the boundary      
      DOUBLE PRECISION BC_VEL
      INTEGER TMP_FACTOR
! indices      
      INTEGER LL, LC, I, J, IJ, M, MM
! a random number between 0 and 1
      DOUBLE PRECISION TMP_DP
! a random integer between 1 and tmp_factor      
      INTEGER TMP_INT
! Temp variable for DES_MI_CLASS(BCV_I)
      CHARACTER*3 DMCL

! lower threshold velocities that will work given the flow boundary
! specifications:
! lower maximum inlet particle velocity (MAXIPV); above which particles
! enter fast enough not to require additional handling
! lowest minimum inlet particle velocity (MINIPV) below which particles
! do not enter fast enough to meet the massflow rate and an error is
! flagged; above minipv and below maxipv the particle inlet conditions
! are adequate but must be carefully controlled
      DOUBLE PRECISION MAXIPV, MINIPV 

!-----------------------------------------------

      IF(LEN2 == ZERO)THEN   ! 2D domain
!----------------------------------------------- 
         TMP_LEN1 = FLOOR(real(LEN1/MAX_DIA))
         TMP_LEN2 = ZERO
! notes :
!   dtsolid*pi_factor(:)        = time elapsed between particle injections
!                               = des_mi_time(:)
!   d_p0/(dtsolid*pi_factor(:)) = approx velocity needed to move one particle
!                                 diameter at the specified mass flow rate
!   ceiling(tmp_len1/2)         = the minimum no. of particles that can
!                                 be arranged along the inlet so that an 
!                                 additional particle cannot fit

         MAXIPV = MAX_DIA/( DTSOLID*dble(PI_FACTOR(BCV_I))*&
                  dble( CEILING(real(TMP_LEN1)/2.0)) ) 
         MINIPV = MAX_DIA/( DTSOLID*dble(PI_FACTOR(BCV_I))*&
                  dble(TMP_LEN1) )
         IF (MINIPV .LT. SMALL_NUMBER) MINIPV = ZERO

         DMCL = DES_MI_CLASS(BCV_I)

! Check each solid phase to determine how the particles are to be
! placed and if the mass inflow conditions can be achevied.
         IF(DMCL == 'XW' .OR. DMCL == "XE")THEN
            BC_VEL = ABS(DES_BC_U_s(BCV))
         ELSEIF(DMCL == 'YS' .OR. DMCL == "YN")THEN
            BC_VEL = ABS(DES_BC_V_s(BCV))
         ENDIF

! The bc velocity is lowered by a small number to prevent possible
! issues with the mantissa truncation
         IF(MAXIPV .LE. (BC_VEL - SMALL_NUMBER))THEN
! The inlet velocity is sufficient to permit random placement of the new
! particles without risk of overlap
            PARTICLE_PLCMNT(BCV_I) = 'RAND'
! Override the random placement for an ordered inlet
            IF(FORCE_ORD_BC) PARTICLE_PLCMNT(BCV_I) = 'ORDR'
         ELSEIF(MINIPV .LE. (BC_VEL - SMALL_NUMBER) .AND. &
            (BC_VEL .LT. MAXIPV + SMALL_NUMBER)) THEN
! Then inlet velocity will require that the new particles be placed with
! order to prevent overlap.
                PARTICLE_PLCMNT(BCV_I) = 'ORDR'
            ELSE
! The inlet velocity is too low to satisfy the other inlet conditions.
! Determine maximum possible ROP_s values that will satisfy the inlet
! conditions. Flag error, prompt with new values, and exit.
               WRITE(UNIT_LOG,1300) BCV
               WRITE(*,1300) BCV
               DO MM = 1, PHASE_CNT
                  M = PHASE_LIST(MM)
! Even though the system is 2D, an area of the inlet is needed for
! the following calculation.  This 'depth' is taken as max_dia.
                  MAX_ROPs = (DES_BC_VOLFLOW_s(BCV,M)*RO_s(M)) / &
                     (MINIPV * LEN1 *  MAX_DIA)
                  WRITE(UNIT_LOG,1301) BCV, M, MAX_ROPs
                  WRITE(*,1301) BCV, M, MAX_ROPs
               ENDDO
               WRITE(UNIT_LOG,1302)
               WRITE(*,1302)
               CALL MFIX_EXIT(myPE)
            ENDIF

         IF(PARTICLE_PLCMNT(BCV_I) == 'ORDR')THEN
! In 2D calculate the approx. particle line density (e.g., the no. of 
! particles along the inlet) needed to achieve the specified particle 
! mass flow rate and particle velocity;
            TMP_FACTOR = CEILING(real(MAX_DIA / &
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
               TMP_INT = CEILING(real(TMP_DP*dble(TMP_FACTOR)))
               DO LC = 1, LL
                 IF(TMP_INT .EQ. MI_ORDER(BCV_I)%VALUE(LC) )EXIT
                 IF(LC .EQ. LL)THEN
                    MI_ORDER(BCV_I)%VALUE(LC) = TMP_INT
                    LL = LL + 1
                 ENDIF
               ENDDO
            ENDDO           
         ENDIF     ! endif particle_plcmnt(bcv_i) == 'ordr'


      ELSEIF(LEN2 /= ZERO) THEN   ! 3D domain
!-----------------------------------------------      
         TMP_LEN1 = FLOOR(real(LEN1/MAX_DIA))
         TMP_LEN2 = FLOOR(real(LEN2/MAX_DIA))

! In the 3D case the calculation for MAXIPV is conservative.  That is,
! the actual bc velocity could be somewhat lower than the calculated 
! value of MAXIPV and still allow for random particle placement
         MAXIPV = MAX_DIA/( DTSOLID*dble(PI_FACTOR(BCV_I)) * &
                     dble( CEILING(real(TMP_LEN1*TMP_LEN2)/2.0)) ) 
! The cutoff is associated with square packing of disks on a plane
! A lower velocity would be possible with hexagonal packing
         MINIPV = MAX_DIA/( DTSOLID*dble(PI_FACTOR(BCV_I)) * &
                     dble(TMP_LEN1*TMP_LEN2) )

         IF (MINIPV .LT. SMALL_NUMBER) MINIPV = ZERO

         DMCL = DES_MI_CLASS(BCV_I)
! Check each solid phase to determine how the particles are to be
! placed or if the mass inflow conditions can be achevied.
         IF(DMCL == 'YZe' .OR. DMCL == "YZw")THEN
            BC_VEL = ABS(DES_BC_U_s(BCV))
         ELSEIF(DMCL == 'XZs' .OR. DMCL == "XZn")THEN
            BC_VEL = ABS(DES_BC_V_s(BCV))
         ELSEIF(DMCL == 'XYb' .OR. DMCL == "XYt")THEN
            BC_VEL = ABS(DES_BC_W_s(BCV))
         ENDIF
! The bc velocity is lowered by a small number to prevent possible
! issues with the mantissa truncation
         IF(MAXIPV .LE. (ABS(BC_VEL) - SMALL_NUMBER))THEN
! The inlet velocity is sufficient to permit random placement of the new
! particles without risk of overlap
            PARTICLE_PLCMNT(BCV_I) = 'RAND'
! Override the random placement for an ordered inlet
            IF(FORCE_ORD_BC) PARTICLE_PLCMNT(BCV_I) = 'ORDR'
         ELSEIF(MINIPV .LE. (BC_VEL - SMALL_NUMBER) .AND. &
            (BC_VEL .LT. MAXIPV + SMALL_NUMBER)) THEN
! Then inlet velocity will require that the new particles be placed with
! order to prevent overlap.
             PARTICLE_PLCMNT(BCV_I) = 'ORDR'
         ELSE
! The inlet velocity is too low to satisfy the other inlet conditions.
! Determine maximum possible ROP_s values that will satisfy the inlet
! conditions. Flag error, propt with new values, and exit.
            WRITE(UNIT_LOG,1300)BCV; WRITE(*,1300)BCV
            DO MM = 1, PHASE_CNT
               M = PHASE_LIST(MM)
               MAX_ROPs = (DES_BC_VOLFLOW_s(BCV,M)*RO_s(M)) / &
                  (MINIPV * LEN1 * LEN2)
               WRITE(UNIT_LOG,1301)BCV,M,MAX_ROPs
               WRITE(*,1301)BCV,M,MAX_ROPs
            ENDDO
            WRITE(UNIT_LOG,1302)
            WRITE(*,1302)
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
               TMP_INT = CEILING(real(TMP_DP*dble(TMP_FACTOR)))
               DO LC = 1, LL
                 IF(TMP_INT .EQ. MI_ORDER(BCV_I)%VALUE(LC) )EXIT
                 IF(LC .EQ. LL)THEN
                    MI_ORDER(BCV_I)%VALUE(LC) = TMP_INT
                    LL = LL + 1
                 ENDIF
               ENDDO
            ENDDO
         ENDIF   ! endif particle_plcmnt(bcv_i) == 'ordr'

      ENDIF   ! endif len2 == zero


      WRITE(*,1303) BCV, LEN1, LEN2, TMP_LEN1, TMP_LEN2,&
         MAXIPV, MINIPV, PARTICLE_PLCMNT(BCV_I)


 1300 FORMAT(/1X,70('*')//,' From: DES_MI_LAYOUT -',/,&
         ' Message: DES_BC_ROP_s values for BC ',I3,' are set too ',&
         'high.',/,10X,'The MAXIMUM values for the supplied inlet ',&
         'conditions are:')

 1301 FORMAT(10X,'DES_BC_ROP_s(',I3,',',I3,') < ',F9.4)

 1302 FORMAT(1X,70('*'))

 1303 FORMAT(5X,'From: DES_MI_LAYOUT - BC: ', I3, /7X,&
         'LEN1 = ', ES17.5,4X, ' LEN2 = ', ES17.5,/7X,&
         'TMP_LEN1 = ', I6,11X, ' TMP_LEN2 = ', I6,/7X,&
         'MAXIPV = ', ES17.5,2X,' MINIPV = ', ES17.5,/7X,&
         'PLCMNT = ',3X,A,//)

      RETURN
      END SUBROUTINE DES_MI_LAYOUT



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: DES_MI_CELLS                                           !
!                                                                      !
!  Purpose:                                                            !
!    Locate the i, j, k position of the cells encompassing the DES
!    mass inlet BC so that when particles are injected into the system,
!    a complete grid search is not necessary to prevent/determine:
!      1) the injected particle from overlapping with an existing
!         particle nearby when particle_plcmnt = 'rand'      
!      2) the i, j, k indices of the newly injected particle on the 
!         eulerian grid defined by imax, jmax, kmax
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
! note dx(1) = length of ghost cell (x < zero)
         I=2; LOCATION = ZERO
         DO WHILE (LOCATION <= XLENGTH)
            IF((DES_BC_X_w(BCV) >= LOCATION) .AND. &
            (DES_BC_X_w(BCV) < LOCATION + DX(I))) THEN
               GS_ARRAY(BCV_I,1) = I
            ENDIF
            IF((DES_BC_X_e(BCV) >= LOCATION) .AND. &
            (DES_BC_X_e(BCV) < LOCATION + DX(I))) THEN
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
! note dy(1) = length of ghost cell (y < zero)
         J=2; LOCATION = ZERO
         DO WHILE (LOCATION <= YLENGTH)
            IF((DES_BC_Y_s(BCV) >= LOCATION) .AND. &
            (DES_BC_Y_s(BCV) < LOCATION + DY(J))) THEN
               GS_ARRAY(BCV_I,3) = J
            ENDIF
            IF((DES_BC_Y_n(BCV) >= LOCATION) .AND. &
            (DES_BC_Y_n(BCV) < LOCATION + DY(J))) THEN
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
            DO WHILE (LOCATION <= ZLENGTH)
               IF((DES_BC_Z_b(BCV) >= LOCATION) .AND. &
               (DES_BC_Z_b(BCV) < LOCATION + DZ(K))) THEN
                  GS_ARRAY(BCV_I,5) = K
               ENDIF
               IF((DES_BC_Z_t(BCV) >= LOCATION) .AND. &
               (DES_BC_Z_t(BCV) < LOCATION + DZ(K))) THEN
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
            IF (DES_BC_X_w(BCV) == ZERO) DES_MO_CLASS(BCV_I) = 'XW'
            IF (DES_BC_X_w(BCV) == XLENGTH) DES_MO_CLASS(BCV_I) = 'XE'
         ENDIF

! Check horizontal mass outlet: 2D 
! ----------------------------------------
         IF(DES_BC_Y_s(BCV) == DES_BC_Y_n(BCV))THEN
            DES_MO_Y = .TRUE.
            IF (DES_BC_Y_s(BCV) == ZERO) DES_MO_CLASS(BCV_I) = 'YS'
            IF (DES_BC_Y_s(BCV) == YLENGTH) DES_MO_CLASS(BCV_I) = 'YN'
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
            IF (DES_BC_X_w(BCV) == ZERO) DES_MO_CLASS(BCV_I) = 'XW'
            IF (DES_BC_X_w(BCV) == XLENGTH) DES_MO_CLASS(BCV_I) = 'XE'

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
            IF (DES_BC_Y_s(BCV) == ZERO) DES_MO_CLASS(BCV_I) = 'YS'
            IF (DES_BC_Y_s(BCV) == YLENGTH) DES_MO_CLASS(BCV_I) = 'YN'

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
            IF (DES_BC_Z_b(BCV) == ZERO) DES_MO_CLASS(BCV_I) = 'ZB'
            IF (DES_BC_Z_b(BCV) == ZLENGTH) DES_MO_CLASS(BCV_I) = 'ZT'

         ENDIF
! End check mass outlet on XY face: 3D 

      ENDIF   ! endif dimn == 2

! Verify that an inlet is not on a face that is connected to a periodic
! boundary condition.  If so, write error message and exit.
! No Xew outlet with X direction periodic walls
      IF((DES_MO_CLASS(BCV_I) == 'XW' .OR. DES_MO_CLASS(BCV_I) == 'XE') &
          .AND. DES_PERIODIC_WALLS_X) THEN
         WRITE (UNIT_LOG, 1500) BCV, 'DES_PERIODIC_WALLS_X'
         WRITE (*, 1500) BCV, 'DES_PERIODIC_WALLS_X'
         CALL MFIX_EXIT(myPE)
      ENDIF
! No Ysn outlet with Y direction periodic walls
      IF((DES_MO_CLASS(BCV_I) == 'YS' .OR. DES_MO_CLASS(BCV_I) == 'YN') &
          .AND. DES_PERIODIC_WALLS_Y) THEN
         WRITE (UNIT_LOG, 1500) BCV, 'DES_PERIODIC_WALLS_Y'
         WRITE (*, 1500) BCV, 'DES_PERIODIC_WALLS_Y'
         CALL MFIX_EXIT(myPE)
      ENDIF
! No Zbt outlet with Z direction periodic walls
      IF((DES_MO_CLASS(BCV_I) == 'ZB' .OR. DES_MO_CLASS(BCV_I) == 'ZT') &
          .AND. DES_PERIODIC_WALLS_Z) THEN
         WRITE (UNIT_LOG, 1500) BCV, 'DES_PERIODIC_WALLS_Z'
         WRITE (*, 1500) BCV, 'DES_PERIODIC_WALLS_Z'
         CALL MFIX_EXIT(myPE)
      ENDIF

 1500 FORMAT(/1X,70('*')//,' From: DES_MI_CLASSIFY -',/10X,&
         'DEM outlets can not be placed on a periodic boundary.',/10X,&
         'Check DEM boundary condtion ',I2,' and ',A,'.',/1X,70('*')/)

      RETURN
      END SUBROUTINE DES_MO_CLASSIFY


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: DES_MIO_PERIODIC                                       !
!                                                                      !
!  Purpose: ! Verify that any inlet or outlet is placed at least       !
!  "SPACER" distace away from a periodic boundary condition. This      !
!  will prevent any particle that crossing a periodic boundary from    !
!  being moved into contact with an entering or exiting particle.      !
!  This is necessary since particles are not incrementally moved       !
!  across a periodic boundary.                                         !
!                                                                      !
!  Author: J.Musser                                   Date: 22-Dec-09  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE DES_MIO_PERIODIC(BCV,BCV_I,CLASS, MAX_DIA)

      USE compar
      USE des_bc
      USE discretelement
      USE funits
      USE geometry
      USE physprop

      IMPLICIT NONE

!-----------------------------------------------
! Local variables
!-----------------------------------------------

! boundary classification of inlet or outlet (MI or MO)
      CHARACTER*2 CLASS

! boundary classification of inlet or outlet (full text for messages)
      CHARACTER*6 IO_ID

! value of DES_MI/MO_CLASS (boundary plane that MI or MO is 
! located on)
      CHARACTER*3 DMC

      INTEGER BCV   ! absolute boundary condition number
      INTEGER BCV_I ! index boundary condtion number

! Max diameter of incoming particles at bc
      DOUBLE PRECISION MAX_DIA      

! buffer space around inlets and outlets with boardering periodic
! boundary conditions
      DOUBLE PRECISION SPACER

!-----------------------------------------------

      SPACER = 1.05d0 * MAX_DIA

      IF(CLASS == 'MI')THEN
         DMC = DES_MI_CLASS(BCV_I)
         IO_ID = "inlet"
      ELSE
         DMC = DES_MO_CLASS(BCV_I)
         IO_ID = "outlet"
      ENDIF

      IF(DMC == 'XW'  .OR. DMC == 'XE'  .OR. & ! 2D inlet \ 2D-3D outlet
         DMC == 'YZw' .OR. DMC == 'YZe') THEN  ! 3D inlet

         IF(DES_PERIODIC_WALLS_Y) THEN
            IF(DES_BC_Y_s(BCV) < SPACER .OR. &
               DES_BC_Y_n(BCV) > (YLENGTH - SPACER)) THEN
               WRITE(UNIT_LOG, 1600) IO_ID, BCV, 'Y'
               WRITE(*, 1600) IO_ID, BCV, 'Y'
               CALL MFIX_EXIT(myPE)
            ENDIF
         ENDIF
         IF(DES_PERIODIC_WALLS_Z) THEN
            IF(DES_BC_Z_b(BCV) < SPACER .OR. &
               DES_BC_Z_t(BCV) > (ZLENGTH - SPACER)) THEN
               WRITE(UNIT_LOG, 1600) IO_ID, BCV, 'Z'
               WRITE(*, 1600) IO_ID, BCV, 'Z'
               CALL MFIX_EXIT(myPE)
            ENDIF
         ENDIF
      ENDIF

      IF(DMC == 'YS'  .OR. DMC == 'YN'  .OR. & ! 2D inlet \ 2D-3D outlet
         DMC == 'XZs' .OR. DMC == 'XZn') THEN  ! 3D inlet
         IF(DES_PERIODIC_WALLS_X) THEN
            IF(DES_BC_X_w(BCV) < SPACER .OR. &
               DES_BC_X_e(BCV) > (XLENGTH - SPACER)) THEN
               WRITE(UNIT_LOG, 1600) IO_ID, BCV, 'X'
               WRITE(*, 1600) IO_ID, BCV, 'X'
               CALL MFIX_EXIT(myPE)
            ENDIF
         ENDIF
         IF(DES_PERIODIC_WALLS_Z) THEN
            IF(DES_BC_Z_b(BCV) < SPACER .OR. &
               DES_BC_Z_t(BCV) > (ZLENGTH - SPACER)) THEN
               WRITE(UNIT_LOG, 1600) IO_ID, BCV, 'Z'
               WRITE(*, 1600) IO_ID, BCV, 'Z'
               CALL MFIX_EXIT(myPE)
            ENDIF
         ENDIF
      ENDIF

      IF(DMC == 'XYb' .OR. DMC == 'XYt' .OR. & ! 3D inlet
         DMC == 'ZB' .OR. DMC == 'ZT') THEN    ! 3D outlet
         IF(DES_PERIODIC_WALLS_X) THEN
            IF(DES_BC_X_w(BCV) < SPACER .OR. &
               DES_BC_X_e(BCV) > (XLENGTH - SPACER)) THEN
               WRITE(UNIT_LOG, 1600) IO_ID, BCV, 'X'
               WRITE(*, 1600) IO_ID, BCV, 'X'
               CALL MFIX_EXIT(myPE)
            ENDIF
         ENDIF
         IF(DES_PERIODIC_WALLS_Y) THEN
            IF(DES_BC_Y_s(BCV) < SPACER .OR. &
               DES_BC_Y_n(BCV) > (YLENGTH - SPACER)) THEN
               WRITE(UNIT_LOG, 1600) IO_ID, BCV, 'Y'
               WRITE(*, 1600) IO_ID, BCV, 'Y'
               CALL MFIX_EXIT(myPE)
            ENDIF
         ENDIF
      ENDIF

 1600 FORMAT(/1X,70('*')//, ' From: CHECK_DES_BC -',/,&
         ' Message: DES ',A,' boundary ',I3,' is too close to',&
         ' DES_PERIODIC_WALLS_',A,'.',/1X,70('*')/) 

      RETURN
      END SUBROUTINE DES_MIO_PERIODIC

