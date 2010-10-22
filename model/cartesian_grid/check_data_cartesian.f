!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CHECK_DATA_CARTESIAN                                   C
!  Purpose: check the data related to cartesian grid implementation    C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CHECK_DATA_CARTESIAN
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE constant 
      USE run
      USE physprop
      USE indices
      USE scalars
      USE funits
      USE leqsol
      USE compar             
      USE mpi_utility        
      USE bc
      USE DISCRETELEMENT

      USE cutcell
      USE quadric
      USE vtk
      USE polygon
      USE dashboard
      USE stl


      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: G,I,J,IJK,Q
      Character*80  Line(1)
      DOUBLE PRECISION :: norm, tan_half_angle
      CHARACTER(LEN=9) :: GR
!-----------------------------------------------
!

      IF(.NOT.CARTESIAN_GRID) RETURN

      IF(GRANULAR_ENERGY) THEN
         WRITE(*,*)'INPUT ERROR: CARTESIAN GRID OPTION NOT CURRENTLY'
         WRITE(*,*)'AVALAIBLE WHEN SOLVING GRANULAR ENERGY EQUATION.'
         WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
         CALL MFIX_EXIT(MYPE)
      ENDIF

      IF(DISCRETE_ELEMENT) THEN
         WRITE(*,*)'INPUT ERROR: CARTESIAN GRID OPTION NOT CURRENTLY'
         WRITE(*,*)'AVALAIBLE WITH DISCRETE ELEMENT MODEL.'
         WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
         CALL MFIX_EXIT(MYPE)
      ENDIF

      IF(COORDINATES=='CYLINDRICAL') THEN
         WRITE(*,*)'INPUT ERROR: CARTESIAN GRID OPTION NOT AVAILABLE'
         WRITE(*,*)'WITH CYLINDRICAL COORDINATE SYSTEM.'
         WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
         CALL MFIX_EXIT(MYPE)
      ENDIF

      IF(USE_STL) THEN
         IF(DO_K) THEN 
            CALL GET_STL_DATA
         ELSE
            WRITE(*,*) 'ERROR: STL METHOD VALID ONLY IN 3D.'
            CALL MFIX_EXIT(MYPE) 
         ENDIF
         IF(N_QUADRIC > 0) THEN
            WRITE(*,*) 'ERROR: BOTH QUADRIC(S) AND STL INPUT ARE SPECIFIED.'
            WRITE(*,*) 'MFIX HANDLES ONLY ONE TYPE OF SURFACE INPUT.'
            CALL MFIX_EXIT(MYPE) 
         ENDIF
         IF(STL_BC_ID == UNDEFINED_I) THEN
            WRITE(*,*) 'ERROR: STL_BC_ID NOT DEFINED.'
            CALL MFIX_EXIT(MYPE) 
         ENDIF
      ENDIF

      IF(USE_POLYGON) THEN
         IF(DO_K) THEN 
            WRITE(*,*) 'ERROR: POLYGON METHOD VALID ONLY IN 2D.'
            CALL MFIX_EXIT(MYPE) 
         ELSE
            CALL GET_POLY_DATA
         ENDIF
      ENDIF

      IF(N_QUADRIC > 0) THEN
         IF(N_POLYGON > 0) THEN 
            WRITE(*,*) 'ERROR: BOTH QUADRIC(S) AND POLYGON(S) DEFINED.'
            WRITE(*,*) 'MFIX HANDLES ONLY ONE TYPE OF SURFACE INPUT.'
            CALL MFIX_EXIT(MYPE) 
         ENDIF
         IF(N_USR_DEF > 0) THEN 
            WRITE(*,*) 'ERROR: BOTH QUADRIC(S) AND USER-DEFINED FUNTION DEFINED.'
            WRITE(*,*) 'MFIX HANDLES ONLY ONE TYPE OF SURFACE.'
            CALL MFIX_EXIT(MYPE) 
         ENDIF
      ELSE
         IF((N_POLYGON > 0).AND.(N_USR_DEF > 0)) THEN 
            WRITE(*,*) 'ERROR: POLYGON(S) AND USER-DEFINED FUNTION DEFINED.'
            WRITE(*,*) 'MFIX HANDLES ONLY ONE TYPE OF SURFACE.'
            CALL MFIX_EXIT(MYPE) 
         ENDIF
      ENDIF

      
      IF(N_QUADRIC > DIM_QUADRIC) THEN
         WRITE(*,*)'INPUT ERROR: INVALID VALUE OF N_QUADRIC =', N_QUADRIC
         WRITE(*,*)'MAXIMUM ACCEPTABLE VALUE IS DIM_QUADRIC =', DIM_QUADRIC
         WRITE(*,*)'CHANGE MAXIMUM VALUE IN QUADRIC_MOD.F IF NECESSARY.'
         WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
         CALL MFIX_EXIT(MYPE)
      ENDIF


      DO Q = 1, N_QUADRIC



         SELECT CASE (TRIM(QUADRIC_FORM(Q)))
    
            CASE ('NORMAL')       

               norm = dsqrt(lambda_x(Q)**2 + lambda_y(Q)**2 + lambda_z(Q)**2)

               IF(norm < TOL_F) THEN
                  WRITE(*,*)'INPUT ERROR: QUADRIC:', Q, ' HAS ZERO COEFFICIENTS.'
                  WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
                  CALL MFIX_EXIT(MYPE)             
               ENDIF

            CASE ('PLANE')   ! The quadric is predefined as a plane

               lambda_x(Q) = n_x(Q)
               lambda_y(Q) = n_y(Q)
               lambda_z(Q) = n_z(Q)

               norm = dsqrt(lambda_x(Q)**2 + lambda_y(Q)**2 + lambda_z(Q)**2)

               IF( norm > TOL_F) THEN
                  lambda_x(Q) = lambda_x(Q) / norm
                  lambda_y(Q) = lambda_y(Q) / norm
                  lambda_z(Q) = lambda_z(Q) / norm
               ELSE
                  WRITE(*,*)'INPUT ERROR: PLANE:', Q, ' HAS ZERO NORMAL VECTOR.'
                  WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
                  CALL MFIX_EXIT(MYPE)             
               ENDIF

              dquadric(Q) = - (lambda_x(Q)*t_x(Q) + lambda_y(Q)*t_y(Q) + lambda_z(Q)*t_z(Q))

            CASE ('X_CYL_INT')   ! The quadric is predefined as a cylinder, along x-axis
                                 ! Internal flow

               IF( Radius(Q) <= ZERO) THEN
                  WRITE(*,*)'INPUT ERROR: CYLINDER:', Q, ' HAS ZERO RADIUS.'
                  WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
                  CALL MFIX_EXIT(MYPE)             
               ELSE
                  lambda_x(Q) = ZERO
                  lambda_y(Q) = ONE
                  lambda_z(Q) = ONE
                  dquadric(Q) = -Radius(Q)**2
               ENDIF

            CASE ('Y_CYL_INT')   ! The quadric is predefined as a cylinder, along y-axis
                                 ! Internal flow

               IF( Radius(Q) <= ZERO) THEN
                  WRITE(*,*)'INPUT ERROR: CYLINDER:', Q, ' HAS ZERO RADIUS.'
                  WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
                  CALL MFIX_EXIT(MYPE)             
               ELSE
                  lambda_x(Q) = ONE
                  lambda_y(Q) = ZERO
                  lambda_z(Q) = ONE
                  dquadric(Q) = -Radius(Q)**2
               ENDIF

            CASE ('Z_CYL_INT')   ! The quadric is predefined as a cylinder, along z-axis
                                 ! Internal flow

               IF( Radius(Q) <= ZERO) THEN
                  WRITE(*,*)'INPUT ERROR: CYLINDER:', Q, ' HAS ZERO RADIUS.'
                  WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
                  CALL MFIX_EXIT(MYPE)             
               ELSE
                  lambda_x(Q) = ONE
                  lambda_y(Q) = ONE
                  lambda_z(Q) = ZERO
                  dquadric(Q) = -Radius(Q)**2
               ENDIF


            CASE ('X_CYL_EXT')   ! The quadric is predefined as a cylinder, along x-axis
                                 ! External flow

               IF( Radius(Q) <= ZERO) THEN
                  WRITE(*,*)'INPUT ERROR: CYLINDER:', Q, ' HAS ZERO RADIUS.'
                  WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
                  CALL MFIX_EXIT(MYPE)             
               ELSE
                  lambda_x(Q) = ZERO
                  lambda_y(Q) = -ONE
                  lambda_z(Q) = -ONE
                  dquadric(Q) = Radius(Q)**2
               ENDIF

            CASE ('Y_CYL_EXT')   ! The quadric is predefined as a cylinder, along y-axis
                                 ! External flow

               IF( Radius(Q) <= ZERO) THEN
                  WRITE(*,*)'INPUT ERROR: CYLINDER:', Q, ' HAS ZERO RADIUS.'
                  WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
                  CALL MFIX_EXIT(MYPE)             
               ELSE
                  lambda_x(Q) = -ONE
                  lambda_y(Q) = ZERO
                  lambda_z(Q) = -ONE
                  dquadric(Q) = Radius(Q)**2
               ENDIF

            CASE ('Z_CYL_EXT')   ! The quadric is predefined as a cylinder, along z-axis
                                 ! External flow

               IF( Radius(Q) <= ZERO) THEN
                  WRITE(*,*)'INPUT ERROR: CYLINDER:', Q, ' HAS ZERO RADIUS.'
                  WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
                  CALL MFIX_EXIT(MYPE)             
               ELSE
                  lambda_x(Q) = -ONE
                  lambda_y(Q) = -ONE
                  lambda_z(Q) = ZERO
                  dquadric(Q) = Radius(Q)**2
               ENDIF

            CASE ('X_CONE')    ! The quadric is predefined as a cone, along x-axis
                               ! Internal flow

            IF(HALF_ANGLE(Q) <= ZERO .OR. HALF_ANGLE(Q) >= 90.0) THEN
                  WRITE(*,*)'INPUT ERROR: CONE:', Q, ' HAS INCORRECT HALF-ANGLE.'
                  WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
                  CALL MFIX_EXIT(MYPE)             
               ELSE
                  tan_half_angle = DTAN(HALF_ANGLE(Q)/180.0*PI)
                  lambda_x(Q) = -ONE
                  lambda_y(Q) = ONE/(tan_half_angle)**2
                  lambda_z(Q) = ONE/(tan_half_angle)**2
                  dquadric(Q) = ZERO
               ENDIF

            CASE ('Y_CONE')    ! The quadric is predefined as a cone, along y-axis
                               ! Internal flow

            IF(HALF_ANGLE(Q) <= ZERO .OR. HALF_ANGLE(Q) >= 90.0) THEN
                  WRITE(*,*)'INPUT ERROR: CONE:', Q, ' HAS INCORRECT HALF-ANGLE.'
                  WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
                  CALL MFIX_EXIT(MYPE)             
               ELSE
                  tan_half_angle = DTAN(HALF_ANGLE(Q)/180.0*PI)
                  lambda_x(Q) = ONE/(tan_half_angle)**2
                  lambda_y(Q) = -ONE
                  lambda_z(Q) = ONE/(tan_half_angle)**2
                  dquadric(Q) = ZERO
               ENDIF

            CASE ('Z_CONE')    ! The quadric is predefined as a cone, along z-axis
                               ! Internal flow

            IF(HALF_ANGLE(Q) <= ZERO .OR. HALF_ANGLE(Q) >= 90.0) THEN
                  WRITE(*,*)'INPUT ERROR: CONE:', Q, ' HAS INCORRECT HALF-ANGLE.'
                  WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
                  CALL MFIX_EXIT(MYPE)             
               ELSE
                  tan_half_angle = DTAN(HALF_ANGLE(Q)/180.0*PI)
                  lambda_x(Q) = ONE/(tan_half_angle)**2
                  lambda_y(Q) = ONE/(tan_half_angle)**2
                  lambda_z(Q) = -ONE
                  dquadric(Q) = ZERO
               ENDIF


            CASE DEFAULT
               WRITE(*,*)'INPUT ERROR: QUADRIC:', Q, ' HAS INCORRECT FORM: ',quadric_form(Q)
               WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
               CALL MFIX_EXIT(MYPE)             

         END SELECT

         IF(BC_ID_Q(Q) == UNDEFINED_I) THEN
            WRITE(*,*)'INPUT ERROR: QUADRIC:', Q, ' HAS NO ASSIGNED BC ID.'
            WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
            CALL MFIX_EXIT(MYPE)
         ENDIF

      ENDDO
 

      IF(N_QUADRIC>0) THEN


         IF(N_GROUP > DIM_GROUP) THEN
            WRITE(*,*)'INPUT ERROR: INVALID VALUE OF N_GROUP =', N_GROUP
            WRITE(*,*)'MAXIMUM ACCEPTABLE VALUE IS DIM_GROUP =', DIM_GROUP
            WRITE(*,*)'CHANGE MAXIMUM VALUE IN QUADRIC_MOD.F IF NECESSARY.'
            WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
            CALL MFIX_EXIT(MYPE)
         ENDIF


         DO I = 1,N_GROUP

            IF(GROUP_SIZE(I) < 1 .OR. GROUP_SIZE(I) > N_QUADRIC) THEN
               WRITE(*,*)'INPUT ERROR: GROUP:', I, ' HAS INCORRECT SIZE:', GROUP_SIZE(I)
               WRITE(*,*)'VALID GROUP SIZE RANGE IS:', 1, ' TO ', N_QUADRIC
               WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
               CALL MFIX_EXIT(MYPE)
            ENDIF

            DO J = 1,GROUP_SIZE(I)
               IF(GROUP_Q(I,J) < 1 .OR. GROUP_Q(I,J) > N_QUADRIC) THEN
                  WRITE(*,*)'INPUT ERROR: GROUP_Q(', I,',',J, ') HAS INCORRECT VALUE:', GROUP_Q(I,J)
                  WRITE(*,*)'VALID GROUP_Q RANGE IS:', 1, ' TO ', N_QUADRIC
                  WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
                  CALL MFIX_EXIT(MYPE)
               ENDIF
            ENDDO
   
            GR = TRIM(GROUP_RELATION(I)) 

            IF(GR/='OR'.AND.GR/='AND'.AND.GR/='PIECEWISE') THEN
               WRITE(*,*)'INPUT ERROR: GROUP:', I, ' HAS INCORRECT GROUP RELATION: ', GR
               WRITE(*,*)'VALID GROUP RELATIONS ARE ''OR'',''AND'', AND ''PIECEWISE''. '
               WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
               CALL MFIX_EXIT(MYPE)
            ENDIF

         ENDDO

         DO I = 2,N_GROUP

            GR = TRIM(RELATION_WITH_PREVIOUS(I)) 

            IF(GR/='OR'.AND.GR/='AND') THEN
               WRITE(*,*)'INPUT ERROR: GROUP:', I, ' HAS INCORRECT RELATION WITH PREVIOUS: ', GR
               WRITE(*,*)'VALID GROUP RELATIONS ARE ''OR'', AND ''AND''. '
               WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
               CALL MFIX_EXIT(MYPE)
            ENDIF
         
         ENDDO

      ENDIF


      IF(TOL_SNAP(1)<ZERO.OR.TOL_SNAP(1)>HALF) THEN
         WRITE(*,*)'INPUT ERROR: INVALID VALUE OF TOL_SNAP IN X-DIRECTION =', TOL_SNAP(1)
         WRITE(*,*)'ACCEPTABLE VALUES ARE BETWEEN 0.0 AND 0.5.'
         WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
         CALL MFIX_EXIT(MYPE)
      ENDIF   
                
      IF(TOL_SNAP(2)<ZERO.OR.TOL_SNAP(2)>HALF) THEN
         WRITE(*,*)'INPUT ERROR: INVALID VALUE OF TOL_SNAP IN Y-DIRECTION =', TOL_SNAP(2)
         WRITE(*,*)'ACCEPTABLE VALUES ARE BETWEEN 0.0 AND 0.5.'
         WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
         CALL MFIX_EXIT(MYPE)
      ENDIF   
                
      IF(TOL_SNAP(3)<ZERO.OR.TOL_SNAP(3)>HALF) THEN
         WRITE(*,*)'INPUT ERROR: INVALID VALUE OF TOL_SNAP IN Z-DIRECTION =', TOL_SNAP(3)
         WRITE(*,*)'ACCEPTABLE VALUES ARE BETWEEN 0.0 AND 0.5.'
         WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
         CALL MFIX_EXIT(MYPE)
      ENDIF   


      IF(TOL_DELH<ZERO.OR.TOL_DELH>ONE) THEN
         WRITE(*,*)'INPUT ERROR: INVALID VALUE OF TOL_DELH =', TOL_DELH
         WRITE(*,*)'ACCEPTABLE VALUES ARE BETWEEN 0.0 AND 1.0.'
         WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
         CALL MFIX_EXIT(MYPE)
      ENDIF

      IF(TOL_SMALL_CELL<ZERO.OR.TOL_SMALL_CELL>ONE) THEN
         WRITE(*,*)'INPUT ERROR: INVALID VALUE OF TOL_SMALL_CELL =', TOL_SMALL_CELL
         WRITE(*,*)'ACCEPTABLE VALUES ARE BETWEEN 0.0 AND 1.0.'
         WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
         CALL MFIX_EXIT(MYPE)
      ENDIF

      IF(TOL_SMALL_AREA<ZERO.OR.TOL_SMALL_AREA>ONE) THEN
         WRITE(*,*)'INPUT ERROR: INVALID VALUE OF TOL_SMALL_AREA =', TOL_SMALL_AREA
         WRITE(*,*)'ACCEPTABLE VALUES ARE BETWEEN 0.0 AND 1.0.'
         WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
         CALL MFIX_EXIT(MYPE)
      ENDIF

      IF(ALPHA_MAX<ZERO) THEN
         WRITE(*,*)'INPUT ERROR: NEGATIVE VALUE OF ALPHA_MAX =', ALPHA_MAX
         WRITE(*,*)'ACCEPTABLE VALUES ARE POSITIVE NUMBERS (E.G. 1.0).'
         WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
         CALL MFIX_EXIT(MYPE)
      ENDIF


      IF(TOL_F<ZERO) THEN
         WRITE(*,*)'INPUT ERROR: NEGATIVE VALUE OF TOL_F =', TOL_F
         WRITE(*,*)'ACCEPTABLE VALUES ARE SMALL POSITIVE NUMBERS (E.G. 1.0E-9).'
         WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
         CALL MFIX_EXIT(MYPE)
      ENDIF

      IF(TOL_POLY<ZERO) THEN
         WRITE(*,*)'INPUT ERROR: NEGATIVE VALUE OF TOL_POLY =', TOL_POLY
         WRITE(*,*)'ACCEPTABLE VALUES ARE SMALL POSITIVE NUMBERS (E.G. 1.0E-9).'
         WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
         CALL MFIX_EXIT(MYPE)
      ENDIF

      IF(ITERMAX_INT<0) THEN
         WRITE(*,*)'INPUT ERROR: NEGATIVE VALUE OF ITERMAX_INT =', ITERMAX_INT
         WRITE(*,*)'ACCEPTABLE VALUES ARE LARGE POSITIVE INTEGERS (E.G. 10000).'
         WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
         CALL MFIX_EXIT(MYPE)
      ENDIF

      IF(FAC_DIM_MAX_CUT_CELL<0.05.OR.FAC_DIM_MAX_CUT_CELL>ONE) THEN
         WRITE(*,*)'INPUT ERROR: NEGATIVE VALUE OF FAC_DIM_MAX_CUT_CELL =', FAC_DIM_MAX_CUT_CELL
         WRITE(*,*)'ACCEPTABLE VALUES ARE BETWEEN 0.05 AND 1.0.'
         WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
         CALL MFIX_EXIT(MYPE)
      ENDIF

      IF(.NOT.CARTESIAN_GRID) THEN
         IF(WRITE_VTK_FILES) THEN
            WRITE(*,*)'INPUT ERROR: VTK FILES CAN BE WRITTEN ONLY WHEN CARTESIAN GRID IS ACTIVATED.'
            WRITE(*,*)'PLEASE SET WRITE_VTK_FILES = .FALSE. IN MFIX.DAT AND TRY AGAIN.'
            CALL MFIX_EXIT(MYPE) 
         ENDIF
      ENDIF


      IF(VTK_DT<ZERO) THEN
         WRITE(*,*)'INPUT ERROR: NEGATIVE VALUE OF VTK_DT =', VTK_DT
         WRITE(*,*)'ACCEPTABLE VALUES ARE POSITIVE NUMBERS (E.G. 0.1).'
         WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
         CALL MFIX_EXIT(MYPE)
      ENDIF

      IF(FRAME<-1) THEN
         WRITE(*,*)'INPUT ERROR: INALID VALUE OF FRAME =', FRAME
         WRITE(*,*)'ACCEPTABLE VALUES ARE INTEGERS >= -1.'
         WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
         CALL MFIX_EXIT(MYPE)
      ENDIF

      IF((CG_SAFE_MODE(1)==1).AND.(PG_OPTION/=0)) THEN
         PG_OPTION = 0
         WRITE(*,*)'WARNING: SAFE_MODE ACTIVATED FOR GAS PRESSURE, REVERTING TO PG_OPTION = 0'
      ENDIF

      IF(PG_OPTION <0 .OR. PG_OPTION>2) THEN
         WRITE(*,*)'INPUT ERROR: INALID VALUE OF PG_OPTION =', PG_OPTION
         WRITE(*,*)'ACCEPTABLE VALUES ARE 0,1,AND 2.'
         WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
         CALL MFIX_EXIT(MYPE)
      ENDIF

      IF(CG_UR_FAC(2)<ZERO.OR.CG_UR_FAC(2)>ONE) THEN
         WRITE(*,*)'INPUT ERROR: INVALID VALUE OF CG_UR_FAC(2) =', CG_UR_FAC(2)
         WRITE(*,*)'ACCEPTABLE VALUES ARE BETWEEN 0.0 AND 1.0.'
         WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
         CALL MFIX_EXIT(MYPE)
      ENDIF

      IF(BAR_WIDTH<10.OR.BAR_WIDTH>80) THEN
         WRITE(*,*)'INPUT ERROR: INVALID VALUE OF BAR_WIDTH =', BAR_WIDTH
         WRITE(*,*)'ACCEPTABLE VALUES ARE BETWEEN 10 AND 80.'
         WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
         CALL MFIX_EXIT(MYPE)
      ENDIF

      IF(BAR_RESOLUTION<ONE.OR.BAR_RESOLUTION>100.0) THEN
         WRITE(*,*)'INPUT ERROR: INVALID VALUE OF BAR_RESOLUTION =', BAR_RESOLUTION
         WRITE(*,*)'ACCEPTABLE VALUES ARE BETWEEN 0.0 AND 100.0.'
         WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
         CALL MFIX_EXIT(MYPE)
      ENDIF

      IF(F_DASHBOARD<1) THEN
         WRITE(*,*)'INPUT ERROR: INALID VALUE OF F_DASHBOARD =', F_DASHBOARD
         WRITE(*,*)'ACCEPTABLE VALUES ARE INTEGERS >= 1.'
         WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
         CALL MFIX_EXIT(MYPE)
      ENDIF


!======================================================================
! Data initialization for Dashboard
!======================================================================
      INIT_TIME = TIME
      SMMIN =  LARGE_NUMBER
      SMMAX = -LARGE_NUMBER

      DTMIN =  LARGE_NUMBER
      DTMAX = -LARGE_NUMBER

      NIT_MIN = MAX_NIT
      NIT_MAX = 0

      N_DASHBOARD = 0

      RETURN  
      END SUBROUTINE CHECK_DATA_CARTESIAN



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CHECK_BC_FLAGS                                         C
!  Purpose: check the boundary conditions flags                        C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CHECK_BC_FLAGS
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE constant 
      USE run
      USE physprop
      USE indices
      USE scalars
      USE funits
      USE leqsol
      USE compar             
      USE mpi_utility        
      USE bc
      
      USE fldvar
      USE cutcell
      USE quadric
      USE vtk
      USE polygon
      USE dashboard


      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I,J,IJK,IJKW,IJKS,IJKB,M,N
      INTEGER :: IJKWW,IJKSS,IJKBB
      INTEGER :: BCV,BCV_U,BCV_V,BCV_W
      Character*80  Line(1)
!-----------------------------------------------
!
      DOUBLE PRECISION SUM, SUM_EP
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      LOGICAL , EXTERNAL :: COMPARE 
!-----------------------------------------------
      INCLUDE 'function.inc'

!======================================================================
! Boundary conditions
!======================================================================

      DO IJK = ijkstart3, ijkend3
         BCV = BC_ID(IJK)
         IF(BCV>0) THEN

            IF(BC_TYPE(BCV)  == 'CG_MI') THEN
 
               FLAG(IJK) = 20
               FLAG_E(IJK) = UNDEFINED_I
               FLAG_N(IJK) = UNDEFINED_I
               FLAG_T(IJK) = UNDEFINED_I

            ELSEIF(BC_TYPE(BCV)  == 'CG_PO') THEN
 
               FLAG(IJK) = 11
               FLAG_E(IJK) = UNDEFINED_I
               FLAG_N(IJK) = UNDEFINED_I
               FLAG_T(IJK) = UNDEFINED_I

               IJKW = WEST_OF(IJK)
               BCV_U = BC_U_ID(IJKW)
               IF(BCV_U>0) THEN
                  IF(BC_TYPE(BCV_U)  == 'CG_PO') THEN
                     FLAG(IJKW) = 11
                     FLAG_E(IJKW) = UNDEFINED_I
                     FLAG_N(IJKW) = UNDEFINED_I
                     FLAG_T(IJKW) = UNDEFINED_I
                  ENDIF
               ENDIF

               IJKS = SOUTH_OF(IJK)
               BCV_V = BC_V_ID(IJKS)
               IF(BCV_V>0) THEN
                  IF(BC_TYPE(BCV_V)  == 'CG_PO') THEN
                     FLAG(IJKS) = 11
                     FLAG_E(IJKS) = UNDEFINED_I
                     FLAG_N(IJKS) = UNDEFINED_I
                     FLAG_T(IJKS) = UNDEFINED_I
                  ENDIF
               ENDIF

               IF(DO_K) THEN
                  IJKB = BOTTOM_OF(IJK)
                  BCV_W = BC_W_ID(IJKB)
                  IF(BCV_W>0) THEN
                     IF(BC_TYPE(BCV_W)  == 'CG_PO') THEN
                        FLAG(IJKB) = 11
                        FLAG_E(IJKB) = UNDEFINED_I
                        FLAG_N(IJKB) = UNDEFINED_I
                        FLAG_T(IJKB) = UNDEFINED_I
                     ENDIF
                  ENDIF
               ENDIF

            ENDIF
         ENDIF
      ENDDO


      DO IJK = ijkstart3, ijkend3
         BCV = BC_ID(IJK)
         IF(BCV>0) THEN
            IF(BC_TYPE(BCV)  == 'CG_MI') THEN

               IJKW = WEST_OF(IJK)
               IF(FLUID_AT(IJKW)) THEN
                  FLAG_E(IJKW) = 2020
               ENDIF           

               IJKS = SOUTH_OF(IJK)
               IF(FLUID_AT(IJKS)) THEN
                  FLAG_N(IJKS) = 2020
               ENDIF           

               IJKB = BOTTOM_OF(IJK)
               IF(FLUID_AT(IJKB)) THEN
                  FLAG_T(IJKB) = 2020
               ENDIF           

               IF (BC_U_G(BCV) == UNDEFINED) THEN 
                   IF (NO_I) THEN 
                       BC_U_G(BCV) = ZERO 
                   ELSE 
                       IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 'BC_U_g', BCV 
                       call mfix_exit(myPE)
                   ENDIF 
               ENDIF 
               IF (BC_V_G(BCV) == UNDEFINED) THEN 
                   IF (NO_J) THEN 
                       BC_V_G(BCV) = ZERO 
                   ELSE 
                       IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 'BC_V_g', BCV 
                       call mfix_exit(myPE)
                   ENDIF 
               ENDIF 
               IF (BC_W_G(BCV) == UNDEFINED) THEN 
                   IF (NO_K) THEN 
                       BC_W_G(BCV) = ZERO 
                   ELSE 
                       IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 'BC_W_g', BCV 
                       call mfix_exit(myPE)
                   ENDIF 
               ENDIF  
               IF (K_Epsilon .AND. BC_K_Turb_G(BCV) == UNDEFINED) THEN
                   IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 'BC_K_Turb_G', BCV 
                   call mfix_exit(myPE) 
               ENDIF   
               IF (K_Epsilon .AND. BC_E_Turb_G(BCV) == UNDEFINED) THEN
                   IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 'BC_E_Turb_G', BCV 
                   call mfix_exit(myPE) 
               ENDIF 

!               Check whether the bc velocity components have the correct sign
!               SELECT CASE (BC_PLANE(BCV))  
!               CASE ('W')  
!                   IF (BC_U_G(BCV) > ZERO) THEN 
!                       IF(DMP_LOG)WRITE (UNIT_LOG, 1050) BCV, 'BC_U_g', '<' 
!                       CALL MFIX_EXIT(myPE) 
!                   ENDIF 
!               CASE ('E')  
!                   IF (BC_U_G(BCV) < ZERO) THEN 
!                       IF(DMP_LOG)WRITE (UNIT_LOG, 1050) BCV, 'BC_U_g', '>' 
!                       CALL MFIX_EXIT(myPE) 
!                   ENDIF 
!               CASE ('S')  
!                   IF (BC_V_G(BCV) > ZERO) THEN 
!                       IF(DMP_LOG)WRITE (UNIT_LOG, 1050) BCV, 'BC_V_g', '<' 
!                       CALL MFIX_EXIT(myPE) 
!                   ENDIF 
!               CASE ('N')  
!                   IF (BC_V_G(BCV) < ZERO) THEN 
!                       IF(DMP_LOG)WRITE (UNIT_LOG, 1050) BCV, 'BC_V_g', '>' 
!                       CALL MFIX_EXIT(myPE) 
!                   ENDIF 
!               CASE ('B')  
!                   IF (BC_W_G(BCV) > ZERO) THEN 
!                       IF(DMP_LOG)WRITE (UNIT_LOG, 1050) BCV, 'BC_W_g', '<' 
!                       CALL MFIX_EXIT(myPE) 
!                   ENDIF 
!               CASE ('T')  
!                   IF (BC_W_G(BCV) < ZERO) THEN 
!                       IF(DMP_LOG)WRITE (UNIT_LOG, 1050) BCV, 'BC_W_g', '>' 
!                       CALL MFIX_EXIT(myPE) 
!                   ENDIF 
!               END SELECT 

               SUM_EP = BC_EP_G(BCV) 
               DO M = 1, MMAX 
                  IF (BC_ROP_S(BCV,M) == UNDEFINED) THEN 
                     IF (BC_EP_G(BCV) == ONE) THEN 
                        BC_ROP_S(BCV,M) = ZERO 
                     ELSEIF (MMAX == 1) THEN 
                         BC_ROP_S(BCV,M) = (ONE - BC_EP_G(BCV))*RO_S(M) 
                     ELSE 
                         IF(DMP_LOG)WRITE (UNIT_LOG, 1100) 'BC_ROP_s', BCV, M 
                         call mfix_exit(myPE)
                     ENDIF 
                  ENDIF 

                  SUM_EP = SUM_EP + BC_ROP_S(BCV,M)/RO_S(M) 
                  IF (SPECIES_EQ(M)) THEN 
                     SUM = ZERO 
                        DO N = 1, NMAX(M) 
                           IF(BC_X_S(BCV,M,N)/=UNDEFINED)SUM=SUM+BC_X_S(BCV,M,N) 
                        ENDDO 

                     IF (BC_ROP_S(BCV,M)==ZERO .AND. SUM==ZERO) THEN 
                        BC_X_S(BCV,M,1) = ONE 
                        SUM = ONE 
                     ENDIF 

                     DO N = 1, NMAX(M) 
                        IF (BC_X_S(BCV,M,N) == UNDEFINED) THEN 
                           IF(.NOT.COMPARE(ONE,SUM) .AND. DMP_LOG)WRITE (UNIT_LOG,1110)BCV,M,N 
                              BC_X_S(BCV,M,N) = ZERO 
                        ENDIF 
                     ENDDO 

                     IF (.NOT.COMPARE(ONE,SUM)) THEN 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1120) BCV, M 
                           call mfix_exit(myPE)
                     ENDIF 
                  ENDIF 

                  IF (BC_U_S(BCV,M) == UNDEFINED) THEN 
                     IF (BC_ROP_S(BCV,M)==ZERO .OR. NO_I) THEN 
                        BC_U_S(BCV,M) = ZERO 
                     ELSE 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1100) 'BC_U_s', BCV, M 
                            call mfix_exit(myPE)
                     ENDIF 
                  ENDIF 
                  
                  IF (BC_V_S(BCV,M) == UNDEFINED) THEN 
                     IF (BC_ROP_S(BCV,M)==ZERO .OR. NO_J) THEN 
                        BC_V_S(BCV,M) = ZERO 
                     ELSE 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1100) 'BC_V_s', BCV, M 
                            call mfix_exit(myPE)
                     ENDIF 
                  ENDIF 
                  
                  IF (BC_W_S(BCV,M) == UNDEFINED) THEN 
                     IF (BC_ROP_S(BCV,M)==ZERO .OR. NO_K) THEN 
                        BC_W_S(BCV,M) = ZERO 
                     ELSE 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1100) 'BC_W_s', BCV, M 
                           call mfix_exit(myPE)
                     ENDIF 
                  ENDIF 

                  IF (ENERGY_EQ .AND. BC_T_S(BCV,M)==UNDEFINED) THEN 
                     IF (BC_ROP_S(BCV,M) == ZERO) THEN 
                        BC_T_S(BCV,M) = BC_T_G(BCV) 
                     ELSE 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1100) 'BC_T_s', BCV, M 
                           call mfix_exit(myPE)
                     ENDIF 
                  ENDIF 

                  IF (GRANULAR_ENERGY .AND. BC_THETA_M(BCV,M)==UNDEFINED) THEN 
                     IF (BC_ROP_S(BCV,M) == ZERO) THEN 
                        BC_THETA_M(BCV,M) = ZERO 
                     ELSE 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1100) 'BC_Theta_m', BCV, M 
                          call mfix_exit(myPE)
                     ENDIF 
                  ENDIF 

!                   Check whether the bc velocity components have the correct sign
!                    SELECT CASE (TRIM(BC_PLANE(BCV)))  
!                    CASE ('W')  
!                        IF (BC_U_S(BCV,M) > ZERO) THEN 
!                            IF(DMP_LOG)WRITE (UNIT_LOG, 1150) BCV, 'BC_U_s', M, '<' 
!                            CALL MFIX_EXIT(myPE) 
!                        ENDIF 
!                    CASE ('E')  
!                        IF (BC_U_S(BCV,M) < ZERO) THEN 
!                            IF(DMP_LOG)WRITE (UNIT_LOG, 1150) BCV, 'BC_U_s', M, '>' 
!                            CALL MFIX_EXIT(myPE) 
!                        ENDIF 
!                    CASE ('S')  
!                        IF (BC_V_S(BCV,M) > ZERO) THEN 
!                            IF(DMP_LOG)WRITE (UNIT_LOG, 1150) BCV, 'BC_V_s', M, '<' 
!                            CALL MFIX_EXIT(myPE) 
!                        ENDIF 
!                    CASE ('N')  
!                        IF (BC_V_S(BCV,M) < ZERO) THEN 
!                            IF(DMP_LOG)WRITE (UNIT_LOG, 1150) BCV, 'BC_V_s', M, '>' 
!                            CALL MFIX_EXIT(myPE) 
!                        ENDIF 
!                    CASE ('B')  
!                        IF (BC_W_S(BCV,M) > ZERO) THEN 
!                            IF(DMP_LOG)WRITE (UNIT_LOG, 1150) BCV, 'BC_W_s', M, '<' 
!                            CALL MFIX_EXIT(myPE) 
!                        ENDIF 
!                    CASE ('T')  
!                        IF (BC_W_S(BCV,M) < ZERO) THEN 
!                            IF(DMP_LOG)WRITE (UNIT_LOG, 1150) BCV, 'BC_W_s', M, '>' 
!                            CALL MFIX_EXIT(myPE) 
!                        ENDIF 
!                    END SELECT 


               ENDDO 

               IF (.NOT.COMPARE(ONE,SUM_EP)) THEN 
                  IF(DMP_LOG)WRITE (UNIT_LOG, 1125) BCV 
                     call mfix_exit(myPE)  
               ENDIF 
       
               DO N = 1, NScalar
                  IF (BC_Scalar(BCV,N) == UNDEFINED) THEN 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1004) 'BC_Scalar', BCV, N 
                        CALL MFIX_EXIT(myPE)
                  ENDIF 
               ENDDO


            ELSEIF(BC_TYPE(BCV)  == 'CG_PO') THEN

               IJKW = WEST_OF(IJK)
               IF(FLUID_AT(IJKW)) THEN
                  FLAG_E(IJKW) = 2011
               ENDIF           

               BCV_U = BC_U_ID(IJKW)
               IF(BCV_U>0) THEN
                  IF(BC_TYPE(BCV_U)  == 'CG_PO') THEN
                    IJKWW = WEST_OF(IJKW)
                    IF(FLUID_AT(IJKWW)) THEN
                       FLAG_E(IJKWW) = 2011
                    ENDIF           
                  ENDIF
               ENDIF

               IJKS = SOUTH_OF(IJK)
               IF(FLUID_AT(IJKS)) THEN
                  FLAG_N(IJKS) = 2011
               ENDIF           

               BCV_V = BC_V_ID(IJKS)
               IF(BCV_V>0) THEN
                  IF(BC_TYPE(BCV_V)  == 'CG_PO') THEN
                    IJKSS = SOUTH_OF(IJKS)
                    IF(FLUID_AT(IJKSS)) THEN
                       FLAG_N(IJKSS) = 2011
                    ENDIF           
                  ENDIF
               ENDIF


               IF(DO_K) THEN
                  IJKB = BOTTOM_OF(IJK)
                  IF(FLUID_AT(IJKB)) THEN
                     FLAG_T(IJKB) = 2011
                  ENDIF 

                  BCV_W = BC_W_ID(IJKB)
                  IF(BCV_W>0) THEN
                     IF(BC_TYPE(BCV_W)  == 'CG_PO') THEN
                       IJKBB = BOTTOM_OF(IJKB)
                       IF(FLUID_AT(IJKBB)) THEN
                          FLAG_T(IJKBB) = 2011
                       ENDIF           
                     ENDIF
                  ENDIF

               ENDIF

               IF (BC_P_G(BCV) == UNDEFINED) THEN 
                   IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 'BC_P_g', BCV 
                   call mfix_exit(myPE)  
               ELSEIF (BC_P_G(BCV)<=ZERO .AND. RO_G0==UNDEFINED) THEN 
                   IF(DMP_LOG)WRITE (UNIT_LOG, 1010) BCV, BC_P_G(BCV) 
                   call mfix_exit(myPE)  
               ENDIF 

            ENDIF

         ENDIF

      ENDDO

      RETURN  

 1000 FORMAT(/1X,70('*')//' From: CHECK_BC_FLAGS',/' Message: ',A,'(',I2,&
         ') not specified',/1X,70('*')/) 
 1001 FORMAT(/1X,70('*')//' From: CHECK_BC_FLAGS',/&
         ' Message: Illegal BC_TYPE for BC # = ',I2,/'   BC_TYPE = ',A,/&
         '  Valid BC_TYPE are: ') 
 1002 FORMAT(5X,A16) 
 1003 FORMAT(/1X,70('*')//' From: CHECK_BC_FLAGS',/' Message: ',A,'(',I2,&
         ') value is unphysical',/1X,70('*')/) 
 1004 FORMAT(/1X,70('*')//' From: CHECK_BC_FLAGS',/' Message: ',A,'(',I2,',',I2,&
         ') not specified',/1X,70('*')/) 
 1005 FORMAT(/1X,70('*')//' From: CHECK_BC_FLAGS',/' Message: ',A,'(',I2,',',I2,&
         ') value is unphysical',/1X,70('*')/) 
 1010 FORMAT(/1X,70('*')//' From: CHECK_BC_FLAGS',/' Message: BC_P_g( ',I2,&
         ') = ',G12.5,/&
         ' Pressure should be greater than zero for compressible flow',/1X,70(&
         '*')/) 
 1050 FORMAT(/1X,70('*')//' From: CHECK_BC_FLAGS',/' Message: BC number:',I2,&
         ' - ',A,' should be ',A,' zero.',/1X,70('*')/) 
 1060 FORMAT(/1X,70('*')//' From: CHECK_BC_FLAGS',/' Message: BC_X_g(',I2,',',I2&
         ,') not specified',/1X,70('*')/) 
 1065 FORMAT(/1X,70('*')//' From: CHECK_BC_FLAGS',/' Message: BC number:',I2,&
         ' - Sum of gas mass fractions is NOT equal to one',/1X,70('*')/) 
 1100 FORMAT(/1X,70('*')//' From: CHECK_BC_FLAGS',/' Message: ',A,'(',I2,',',I1,&
         ') not specified',/1X,70('*')/) 
 1103 FORMAT(/1X,70('*')//' From: CHECK_BC_FLAGS',/' Message: ',A,'(',I2,',',I1,&
         ') value is unphysical',/1X,70('*')/) 
 1104 FORMAT(/1X,70('*')//' From: CHECK_BC_FLAGS',/' Message: ',A,'(',I2,',',I2,&
         ',',I2,') not specified',/1X,70('*')/) 
 1105 FORMAT(/1X,70('*')//' From: CHECK_BC_FLAGS',/' Message: ',A,'(',I2,',',I2,&
         ',',I2,') value is unphysical',/1X,70('*')/) 
 1110 FORMAT(/1X,70('*')//' From: CHECK_BC_FLAGS',/' Message: BC_X_s(',I2,',',I2&
         ,',',I2,') not specified',/1X,70('*')/) 
 1120 FORMAT(/1X,70('*')//' From: CHECK_BC_FLAGS',/' Message: BC number:',I2,&
         ' - Sum of solids-',I1,' mass fractions is NOT equal to one',/1X,70(&
         '*')/) 
 1125 FORMAT(/1X,70('*')//' From: CHECK_BC_FLAGS',/' Message: BC number:',I2,&
         ' - Sum of volume fractions is NOT equal to one',/1X,70('*')/) 
 1150 FORMAT(/1X,70('*')//' From: CHECK_BC_FLAGS',/' Message: BC number:',I2,&
         ' - ',A,I1,' should be ',A,' zero.',/1X,70('*')/) 
 1160 FORMAT(/1X,70('*')//' From: CHECK_BC_FLAGS',/&
         ' Message: Boundary condition no', &
         I2,' is a second outflow condition.',/1X,&
         '  Only one outflow is allowed.  Consider using P_OUTFLOW.',/1X, 70('*')/) 
 1200 FORMAT(/1X,70('*')//' From: CHECK_BC_FLAGS',/' Message: ',A,'(',I2,&
         ') specified',' for an undefined BC location',/1X,70('*')/) 
 1300 FORMAT(/1X,70('*')//' From: CHECK_BC_FLAGS',/' Message: ',A,'(',I2,',',I1,&
         ') specified',' for an undefined BC location',/1X,70('*')/) 
 1400 FORMAT(/1X,70('*')//' From: CHECK_BC_FLAGS',/&
         ' Message: No initial or boundary condition specified',/&
         '    I       J       K') 
 1410 FORMAT(I5,3X,I5,3X,I5) 
 1420 FORMAT(/1X,70('*')/) 

      END SUBROUTINE CHECK_BC_FLAGS
