!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CHECK_DATA_07                                           C
!  Purpose: Check boundary condition specifications                    C
!     - convert physical locations to i, j, k's (GET_FLOW_BC)          C
!     - compute area of boundary surfaces (GET_BC_AREA)                C
!     - convert mass and volumetric flows to velocities (FLOW_TO_VEL)  C
!     - check specification of physical quantities                     C 
!  Comments:                                                           C
!                                                                      C
!  Author: P. Nicoletti                               Date: 10-DEC-91  C
!  Reviewer: M.SYAMLAL, W.ROGERS, P.NICOLETTI         Date: 27-JAN-92  C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: Check specification of physical quantities.  Call routine  C
!           to convert mass and volumetric flows to velocities.        C
!  Author: M. Syamlal                                 Date: 24-JUL-92  C
!  Reviewer: W. Rogers                                Date: 11-DEC-92  C
!                                                                      C
!  Revision Number: 2                                                  C
!  Purpose: Include Cartesian grid boundary condition types.           C
!  Author: Jeff Dietiker                              Date: 01-JUL-09  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CHECK_DATA_07 

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param 
      USE param1 
      USE geometry
      USE fldvar
      USE physprop
      USE run
      USE bc
      USE indices
      USE funits 
      USE scalars
      USE compar
      USE sendrecv
      USE discretelement
      USE mfix_pic
      USE cutcell
      IMPLICIT NONE
!-----------------------------------------------
! Local Parameters
!-----------------------------------------------
! JFD: for cartesian grid implementation
      INTEGER, PARAMETER :: DIM_BCTYPE = 21 
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! error flag
      LOGICAL :: ERROR, RECOGNIZED_BC_TYPE
! loop/variable indices
      INTEGER :: BCV
      INTEGER :: I, J ,K, IJK
      INTEGER :: M, N
! Solids phase density in BC region.
      DOUBLE PRECISION :: BC_ROs
! Index of inert species
      INTEGER :: INERT
! valid boundary condition types
      CHARACTER*16, DIMENSION(1:DIM_BCTYPE) ::VALID_BC_TYPE = (/&
           'MASS_INFLOW     ', 'MI              ',&
           'MASS_OUTFLOW    ', 'MO              ',&
           'P_INFLOW        ', 'PI              ',&
           'P_OUTFLOW       ', 'PO              ',&
           'FREE_SLIP_WALL  ', 'FSW             ',&
           'NO_SLIP_WALL    ', 'NSW             ',&
           'PAR_SLIP_WALL   ', 'PSW             ',&
           'OUTFLOW         ', 'OF              ',&
           'CG_NSW          ', 'CG_FSW          ',&
           'CG_PSW          ', 'CG_MI           ',&
           'CG_PO           '                     &
            /)

      DOUBLE PRECISION SUM, SUM_EP

!-----------------------------------------------
! External functions
!-----------------------------------------------
      LOGICAL , EXTERNAL :: COMPARE 
!-----------------------------------------------
! Include statement functions
!-----------------------------------------------
      INCLUDE 'function.inc'
!-----------------------------------------------


! DETERMINE WHICH BOUNDARY CONDITION INDICES HAVE VALUES
      L50: DO BCV = 1, DIMENSION_BC 
         BC_DEFINED(BCV) = .FALSE. 
         IF (BC_X_W(BCV) /= UNDEFINED) BC_DEFINED(BCV) = .TRUE. 
         IF (BC_X_E(BCV) /= UNDEFINED) BC_DEFINED(BCV) = .TRUE. 
         IF (BC_Y_S(BCV) /= UNDEFINED) BC_DEFINED(BCV) = .TRUE. 
         IF (BC_Y_N(BCV) /= UNDEFINED) BC_DEFINED(BCV) = .TRUE. 
         IF (BC_Z_B(BCV) /= UNDEFINED) BC_DEFINED(BCV) = .TRUE. 
         IF (BC_Z_T(BCV) /= UNDEFINED) BC_DEFINED(BCV) = .TRUE. 
         IF (BC_I_W(BCV) /= UNDEFINED_I) BC_DEFINED(BCV) = .TRUE. 
         IF (BC_I_E(BCV) /= UNDEFINED_I) BC_DEFINED(BCV) = .TRUE. 
         IF (BC_J_S(BCV) /= UNDEFINED_I) BC_DEFINED(BCV) = .TRUE. 
         IF (BC_J_N(BCV) /= UNDEFINED_I) BC_DEFINED(BCV) = .TRUE. 
         IF (BC_K_B(BCV) /= UNDEFINED_I) BC_DEFINED(BCV) = .TRUE. 
         IF (BC_K_T(BCV) /= UNDEFINED_I) BC_DEFINED(BCV) = .TRUE. 
         IF (BC_TYPE(BCV) == 'DUMMY') BC_DEFINED(BCV) = .FALSE. 
         IF (BC_TYPE(BCV) == 'CG_NSW') BC_DEFINED(BCV) = .TRUE. 
         IF (BC_TYPE(BCV) == 'CG_FSW') BC_DEFINED(BCV) = .TRUE. 
         IF (BC_TYPE(BCV) == 'CG_PSW') BC_DEFINED(BCV) = .TRUE. 
         IF (BC_TYPE(BCV) == 'CG_MI') BC_DEFINED(BCV) = .TRUE. 
         IF (BC_TYPE(BCV) == 'CG_PO') BC_DEFINED(BCV) = .TRUE. 

        IF(BC_TYPE(BCV) /= UNDEFINED_C.AND.BC_TYPE(BCV) /= 'DUMMY') THEN
            RECOGNIZED_BC_TYPE = .FALSE.
            DO I = 1, DIM_BCTYPE 
                VALID_BC_TYPE(I) = TRIM(VALID_BC_TYPE(I))
                IF (VALID_BC_TYPE(I) == BC_TYPE(BCV)) THEN 
                   RECOGNIZED_BC_TYPE = .TRUE.
                   EXIT
                ENDIF 
            ENDDO 
            
            IF(.NOT.RECOGNIZED_BC_TYPE) THEN
               IF(DMP_LOG)WRITE (UNIT_LOG, 1001) BCV, BC_TYPE(BCV) 
               IF(DMP_LOG)WRITE (UNIT_LOG, 1002) VALID_BC_TYPE 
               call mfix_exit(myPE)
            ENDIF
         ENDIF

! Validate the BC postion for all non cut-cell boundaries:
         IF (BC_DEFINED(BCV).AND.BC_TYPE(BCV)(1:2)/='CG') THEN 

            IF (BC_X_W(BCV)==UNDEFINED .AND. BC_I_W(BCV)==UNDEFINED_I) THEN 
               IF (NO_I) THEN 
                  BC_X_W(BCV) = ZERO
               ELSE 
                  IF(DMP_LOG)WRITE (UNIT_LOG, 1000)&
                     'BC_X_w and BC_I_w ', BCV
                  call mfix_exit(myPE)
               ENDIF 
            ENDIF 
            IF (BC_X_E(BCV)==UNDEFINED .AND. BC_I_E(BCV)==UNDEFINED_I) THEN 
               IF (NO_I) THEN 
                  BC_X_E(BCV) = XLENGTH 
               ELSE 
                  IF(DMP_LOG)WRITE (UNIT_LOG, 1000) &
                     'BC_X_e and BC_I_e ', BCV
                  call mfix_exit(myPE)
               ENDIF 
            ENDIF 
            IF (BC_Y_S(BCV)==UNDEFINED .AND. BC_J_S(BCV)==UNDEFINED_I) THEN 
               IF (NO_J) THEN 
                  BC_Y_S(BCV) = ZERO 
               ELSE 
                  IF(DMP_LOG)WRITE (UNIT_LOG, 1000) &
                     'BC_Y_s and BC_J_s ', BCV
                  call mfix_exit(myPE)
               ENDIF 
            ENDIF 
            IF (BC_Y_N(BCV)==UNDEFINED .AND. BC_J_N(BCV)==UNDEFINED_I) THEN 
               IF (NO_J) THEN 
                  BC_Y_N(BCV) = YLENGTH 
               ELSE 
                  IF(DMP_LOG)WRITE (UNIT_LOG, 1000) &
                     'BC_Y_n and BC_J_n ', BCV 
                  call mfix_exit(myPE)
               ENDIF 
            ENDIF 
            IF (BC_Z_B(BCV)==UNDEFINED .AND. BC_K_B(BCV)==UNDEFINED_I) THEN 
               IF (NO_K) THEN 
                  BC_Z_B(BCV) = ZERO 
               ELSE 
                  IF(DMP_LOG)WRITE (UNIT_LOG, 1000) &
                     'BC_Z_b and BC_K_b ', BCV 
                  call mfix_exit(myPE)
               ENDIF 
            ENDIF 
            IF (BC_Z_T(BCV)==UNDEFINED .AND. BC_K_T(BCV)==UNDEFINED_I) THEN 
               IF (NO_K) THEN 
                  BC_Z_T(BCV) = ZLENGTH 
               ELSE 
                  IF(DMP_LOG)WRITE (UNIT_LOG, 1000) &
                     'BC_Z_t and BC_K_t ', BCV 
                  call mfix_exit(myPE)
               ENDIF 
            ENDIF 

            DO I = 1, DIM_BCTYPE 
               VALID_BC_TYPE(I) = TRIM(VALID_BC_TYPE(I))
               IF (VALID_BC_TYPE(I) == BC_TYPE(BCV)) THEN 
                  IF (MOD(I,2) == 0) BC_TYPE(BCV) = VALID_BC_TYPE(I-1)
                  CYCLE  L50 
               ENDIF 
            ENDDO 

            IF(DMP_LOG)WRITE (UNIT_LOG, 1001) BCV, BC_TYPE(BCV) 
            IF(DMP_LOG)WRITE (UNIT_LOG, 1002) VALID_BC_TYPE 
            call mfix_exit(myPE)  

         ENDIF    ! end if (bc_devined(bcv))
      ENDDO L50   ! end loop over (bcv=1,dimension_bc)


      CALL GET_WALLS_BC 


! Find and validate i, j, k locations of flow BC's
      CALL GET_FLOW_BC 


! Compute area of boundary surfaces
      CALL GET_BC_AREA 


! Check the specification of physical quantities needed in FLOW_TO_VEL
      DO BCV = 1, DIMENSION_BC 
         IF (BC_DEFINED(BCV)) THEN 
            IF (BC_TYPE(BCV) == 'MASS_INFLOW') THEN 
               IF (BC_EP_G(BCV) == UNDEFINED) THEN 
                  IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 'BC_EP_g', BCV 
                  call mfix_exit(myPE)
               ENDIF 
               IF (BC_P_G(BCV)==UNDEFINED .AND. RO_G0==UNDEFINED) THEN 
                  IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 'BC_P_g', BCV 
                  call mfix_exit(myPE)
               ELSEIF (BC_P_G(BCV)<=ZERO .AND. RO_G0==UNDEFINED) THEN 
                  IF(DMP_LOG)WRITE (UNIT_LOG, 1010) BCV, BC_P_G(BCV) 
                  call mfix_exit(myPE)
               ENDIF 
               IF ((ENERGY_EQ .OR. RO_G0==UNDEFINED .OR. &
                    MU_G0==UNDEFINED) .AND. BC_T_G(BCV)==UNDEFINED) THEN
                  IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 'BC_T_g', BCV 
                  call mfix_exit(myPE)
               ENDIF 

               IF (SPECIES_EQ(0) .OR. RO_G0==UNDEFINED .AND. &
                   MW_AVG==UNDEFINED) THEN 
! if gas phase species equations are solved or the gas is compressible
! and the average molecular weight is undefined then check sum of 
! the species mass fractions
                  SUM = ZERO 
! sum together those gas phase species mass fractions that are defined
                  DO N = 1, NMAX(0) 
                     IF (BC_X_G(BCV,N) /= UNDEFINED) &
                        SUM = SUM + BC_X_G(BCV,N) 
                  ENDDO 
! set any undefined species mass fractions to zero. Warn the user if 
! the sum of species mass fractions is not one and indicate which 
! species is undefined (mfix will exit in next check).             
                  DO N = 1, NMAX(0) 
                     IF (BC_X_G(BCV,N) == UNDEFINED) THEN 
                        IF (.NOT.COMPARE(ONE,SUM) .AND. DMP_LOG)&
                           WRITE (UNIT_LOG, 1060) BCV, N 
                        BC_X_G(BCV,N) = ZERO 
                     ENDIF 
                  ENDDO 
! if sum of the gas phase species mass fraction not 1...
                  IF (.NOT.COMPARE(ONE,SUM)) THEN 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1065) BCV 
                     call mfix_exit(myPE)
                  ENDIF 
               ENDIF 

! Verify that species mass fractions are defined for mass flow BCs whe
! using variable solids density. Needed to calculation RO_s
               DO M = 1, SMAX
                  IF(SOLVE_ROs(M)) THEN
                     IF (BC_MASSFLOW_S(BCV,M) /= UNDEFINED .OR.        &
                         BC_VOLFLOW_S(BCV,M)  /= UNDEFINED) THEN
! Sum together those species mass fractions that are defined.
                        SUM = ZERO 
                        DO N = 1, NMAX(M) 
                           IF(BC_X_S(BCV,M,N) /= UNDEFINED) &
                              SUM = SUM + BC_X_S(BCV,M,N)
                        ENDDO 
! Set any undefined species mass fractions to zero.
                        DO N = 1, NMAX(M)
                           IF (BC_X_S(BCV,M,N) == UNDEFINED) THEN
! Only warn the user of this action when the sum does not equal one.
                              IF (.NOT.COMPARE(ONE,SUM) .AND. DMP_LOG)&
                                 WRITE (UNIT_LOG, 1110) BCV, M, N
                              BC_X_S(BCV,M,N) = ZERO 
                           ENDIF 
                        ENDDO
! Exit if sum of the  species mass fraction not one.
                        IF (.NOT.COMPARE(ONE,SUM)) THEN 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1120) BCV 
                           CALL MFIX_EXIT(myPE)
                        ENDIF
                     ENDIF   ! solids mass or vol flow
                  ENDIF   ! variable solids density
               ENDDO   ! M=1,MMAX
            ENDIF   ! end if (bc_type(bcv)=='mass_inflow')
         ENDIF   ! end if (bc_defined(bcv))
      ENDDO   ! end loop over (bcv=1,dimension_bc)


! Convert mass and volumetric flows to velocities 
      CALL FLOW_TO_VEL(.TRUE.) 


! Following section to check quantities for flow type boundaries
! ---------------------------------------------------------------->>>
      ERROR = .FALSE. 
      DO BCV = 1, DIMENSION_BC 
         IF (BC_DEFINED(BCV)) THEN 
            SELECT CASE (TRIM(BC_TYPE(BCV)))

            CASE ('MASS_INFLOW')  
! -------------------------------------------->>>
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
! Comments:                                                            C
!     The velocities at the inflow face are fixed and the momentum     C
!     equations are not solved in the inflow cells. Since the flow is  C
!     into the domain all other scalars that are used need to be       C
!     specified (e.g., mass fractions, void fraction, etc.,)           C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

               IF (K_Epsilon .AND. BC_K_Turb_G(BCV) == UNDEFINED) THEN
                  IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 'BC_K_Turb_G', BCV 
                  call mfix_exit(myPE) 
               ENDIF   
               IF (K_Epsilon .AND. BC_E_Turb_G(BCV) == UNDEFINED) THEN
                  IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 'BC_E_Turb_G', BCV 
                  call mfix_exit(myPE) 
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

! If non-zero, check whether the velocity component through the boundary
! plane has the correct sign.
               SELECT CASE (BC_PLANE(BCV))  
               CASE ('W')  
                  IF (BC_U_G(BCV) > ZERO) THEN 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1050) BCV, &
                        'BC_U_g', '<'
                     CALL MFIX_EXIT(myPE) 
                  ENDIF 
               CASE ('E')  
                  IF (BC_U_G(BCV) < ZERO) THEN 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1050) BCV, &
                        'BC_U_g', '>'
                     CALL MFIX_EXIT(myPE) 
                  ENDIF 
               CASE ('S')  
                  IF (BC_V_G(BCV) > ZERO) THEN 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1050) BCV, &
                        'BC_V_g', '<'
                     CALL MFIX_EXIT(myPE) 
                  ENDIF 
               CASE ('N')  
                  IF (BC_V_G(BCV) < ZERO) THEN 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1050) BCV, &
                        'BC_V_g', '>' 
                     CALL MFIX_EXIT(myPE) 
                  ENDIF 
               CASE ('B')  
                  IF (BC_W_G(BCV) > ZERO) THEN 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1050) BCV, &
                        'BC_W_g', '<'
                     CALL MFIX_EXIT(myPE) 
                  ENDIF 
               CASE ('T')  
                  IF (BC_W_G(BCV) < ZERO) THEN 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1050) BCV, &
                        'BC_W_g', '>'
                     CALL MFIX_EXIT(myPE) 
                  ENDIF 
               END SELECT ! end select case (bc_plane(bcv))


! DEM simulations do not employ variables for continuum solids. So do
! not perform checks on unnecessary data.  Note, MPPIC initialization
! and BC's are based on continuum quantities so that these checks are
! necessary when invoking MPPIC.
               IF(.NOT.DISCRETE_ELEMENT .OR. &
                  DES_CONTINUUM_HYBRID  .OR. MPPIC) THEN

! at this point bc_ep_g must be defined
                  SUM_EP = BC_EP_G(BCV) 
                  DO M = 1, SMAX

! the following check on bc_rop_s is redundant with code in flow_to_vel,
! except that the latter has no calls for MFIX to exit
                     IF (BC_ROP_S(BCV,M) == UNDEFINED) THEN 
                        IF (BC_EP_G(BCV) == ONE) THEN 
                           BC_ROP_S(BCV,M) = ZERO 
                        ELSEIF (SMAX > 1 .OR. DES_CONTINUUM_HYBRID) THEN
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1100) &
                              'BC_ROP_s', BCV, M 
                           call mfix_exit(myPE)
                        ENDIF
                     ENDIF

! Verify that species mass fractions have been specified.
                     SUM = ZERO
                     DO N = 1, NMAX(M)
! Sum together solids phase M species mass fractions that are defined.
                        IF(BC_X_S(BCV,M,N)/=UNDEFINED) &
                           SUM=SUM+BC_X_S(BCV,M,N) 
                     ENDDO
! If no solids phase M present and no species are defined for phase M 
! then set mass fraction of solids phase M species 1 to one
                     IF (SUM==ZERO .AND. BC_ROP_S(BCV,M)==ZERO) THEN 
                        BC_X_S(BCV,M,:) = ZERO
                        BC_X_S(BCV,M,1) = ONE 
                        SUM = ONE 
                     ENDIF

! Set the default if needed.
                     IF(NMAX(M)==1 .AND. BC_X_S(BCV,M,1)==UNDEFINED)THEN
                        BC_X_S(BCV,M,:) = ZERO
                        BC_X_S(BCV,M,1) = ONE
                        SUM = ONE
! Set any undefined species mass fractions to zero. Warn the user if 
! the sum of species mass fractions is not one and indicate which
! species is undefined (mfix will exit in next check). 
                     ELSE
                        DO N = 1, NMAX(M)
                           IF (BC_X_S(BCV,M,N) == UNDEFINED) THEN 
                              IF(.NOT.COMPARE(ONE,SUM) .AND. DMP_LOG)  &
                                 WRITE (UNIT_LOG, 1110)BCV, M, N 
                              BC_X_S(BCV,M,N) = ZERO 
                           ENDIF
                        ENDDO
                     ENDIF
! Verify that the sum of the mass fractions sum to one. Otherwise, exit.
                     IF(.NOT.COMPARE(ONE,SUM)) THEN 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1120) BCV, M 
                        CALL MFIX_EXIT(myPE)
                     ENDIF 

! Calculate the solid density.
                     BC_ROs = UNDEFINED
                     IF(SOLVE_ROs(M))THEN
                        INERT = INERT_SPECIES(M)
! Verify that the species mass fraction for the inert material is not
! zero in the BC region when the solids is present.
                        IF(BC_X_S(BCV,M,INERT) == ZERO) THEN
                           IF(BC_ROP_S(BCV,M) /= ZERO) THEN
                              IF(DMP_LOG) THEN
                                 WRITE(*,1400) M, BCV
                                 WRITE(UNIT_LOG,1400) M, BCV
                              ENDIF
                              CALL MFIX_EXIT(myPE)
                           ELSE
! If the solids isn't present, give it the baseline density.
                              BC_ROs = RO_S0(M)
                           ENDIF
                        ELSE
! Calculate the solids density.
                           BC_ROs = RO_S0(M) * X_s0(M,INERT) / &
                              BC_X_S(BCV,M,INERT)
                        ENDIF
                     ELSE
                        BC_ROs = RO_S0(M)
                     ENDIF

! Sanity check on solids phase density.
                     IF(BC_ROs <= ZERO .OR. BC_ROs==UNDEFINED) THEN
                        IF(DMP_LOG)THEN
                           WRITE(*,1500) M, BCV
                           WRITE(UNIT_LOG,1500) M, BCV
                        ENDIF
                        CALL MFIX_EXIT(myPE)
                     ENDIF

! BC_ROP_S may still be undefined if there is only one solids phase
! and the hybrid model is not in use. Back out the bulk density with
! the assigned solids density.
                     IF(BC_ROP_S(BCV,M) == UNDEFINED) &
                        BC_ROP_S(BCV,M) = (ONE - BC_EP_G(BCV))*BC_ROs

! sum of void fraction and solids volume fractions
                     SUM_EP = SUM_EP + BC_ROP_S(BCV,M)/BC_ROs


                     IF (ENERGY_EQ .AND. BC_T_S(BCV,M)==UNDEFINED) THEN 
                        IF (BC_ROP_S(BCV,M) == ZERO) THEN 
                           BC_T_S(BCV,M) = BC_T_G(BCV) 
                        ELSE 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1100) &
                              'BC_T_s', BCV, M 
                           call mfix_exit(myPE)
                        ENDIF 
                     ENDIF 

                     IF (GRANULAR_ENERGY .AND. &
                         BC_THETA_M(BCV,M)==UNDEFINED) THEN 
                        IF (BC_ROP_S(BCV,M) == ZERO) THEN 
                           BC_THETA_M(BCV,M) = ZERO 
                        ELSE 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1100) &
                              'BC_Theta_m', BCV, M 
                           call mfix_exit(myPE)
                        ENDIF 
                     ENDIF 
                                    
                     IF (BC_U_S(BCV,M) == UNDEFINED) THEN 
                        IF (BC_ROP_S(BCV,M)==ZERO .OR. NO_I) THEN 
                           BC_U_S(BCV,M) = ZERO 
                        ELSE 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1100) &
                              'BC_U_s', BCV, M 
                           call mfix_exit(myPE)
                        ENDIF 
                     ENDIF 
                     IF (BC_V_S(BCV,M) == UNDEFINED) THEN 
                        IF (BC_ROP_S(BCV,M)==ZERO .OR. NO_J) THEN 
                           BC_V_S(BCV,M) = ZERO 
                        ELSE 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1100) &
                              'BC_V_s', BCV, M 
                           call mfix_exit(myPE)
                        ENDIF 
                     ENDIF 
                     IF (BC_W_S(BCV,M) == UNDEFINED) THEN 
                        IF (BC_ROP_S(BCV,M)==ZERO .OR. NO_K) THEN 
                           BC_W_S(BCV,M) = ZERO 
                        ELSE 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1100) &
                              'BC_W_s', BCV, M 
                           call mfix_exit(myPE)
                        ENDIF 
                     ENDIF 

! If non-zero, check whether the velocity component through the boundary
! plane has the correct sign.
                     SELECT CASE (TRIM(BC_PLANE(BCV)))  
                     CASE ('W')  
                        IF (BC_U_S(BCV,M) > ZERO) THEN 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1150) BCV, &
                              'BC_U_s', M, '<' 
                           CALL MFIX_EXIT(myPE) 
                        ENDIF 
                     CASE ('E')  
                        IF (BC_U_S(BCV,M) < ZERO) THEN 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1150) BCV, &
                              'BC_U_s', M, '>' 
                           CALL MFIX_EXIT(myPE) 
                        ENDIF 
                     CASE ('S')  
                        IF (BC_V_S(BCV,M) > ZERO) THEN 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1150) BCV, &
                              'BC_V_s', M, '<' 
                           CALL MFIX_EXIT(myPE) 
                        ENDIF 
                     CASE ('N')  
                        IF (BC_V_S(BCV,M) < ZERO) THEN 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1150) BCV, &
                              'BC_V_s', M, '>' 
                           CALL MFIX_EXIT(myPE) 
                        ENDIF 
                     CASE ('B')  
                     IF (BC_W_S(BCV,M) > ZERO) THEN 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1150) BCV, &
                           'BC_W_s', M, '<' 
                        CALL MFIX_EXIT(myPE) 
                     ENDIF 
                     CASE ('T')  
                        IF (BC_W_S(BCV,M) < ZERO) THEN 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1150) BCV, &
                              'BC_W_s', M, '>' 
                           CALL MFIX_EXIT(myPE) 
                        ENDIF 
                     END SELECT   ! end select case (bc_plane(bcv))

                  ENDDO   ! end do loop m=1,smax


! check sum of gas void fraction and all solids volume fractions
                  IF (.NOT.DES_CONTINUUM_HYBRID) THEN
                     IF (.NOT.COMPARE(ONE,SUM_EP)) THEN 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1125) BCV 
                        call mfix_exit(myPE)  
                     ENDIF 
                  ELSE
! sum_ep not necessarily one at this point since discrete particles
! present in hybrid model
                     IF (SUM_EP>ONE .OR. SUM_EP<ZERO) THEN 
                        IF(DMP_LOG) WRITE (UNIT_LOG, 1130) BCV 
                        call mfix_exit(myPE)
                     ENDIF
                  ENDIF
               ELSE   ! discrete_element and .not.des_continuum_hybrid
                      ! and .not.mppic
                  SUM_EP = BC_EP_G(BCV)
                  IF (SUM_EP>ONE .OR. SUM_EP<ZERO) THEN 
                     IF(DMP_LOG) WRITE (UNIT_LOG, 1130) BCV 
                     call mfix_exit(myPE)
                  ENDIF
               ENDIF   ! end if/else (.not.discrete_element .or.
                       !               des_continuum_hybrid .or. mppic)

       
               DO N = 1, NScalar
                  IF (BC_Scalar(BCV,N) == UNDEFINED) THEN 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1004) 'BC_Scalar', BCV, N 
                     CALL MFIX_EXIT(myPE)
                  ENDIF 
               ENDDO
! case mass_inflow
! --------------------------------------------<<<


            CASE ('MASS_OUTFLOW')  
! -------------------------------------------->>>
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
! Comments:                                                            C
!     The velocities at the outflow face are fixed and the momentum    C
!     equations are not solved in the outflow cells. Since the flow    C
!     is out of the domain none of the other scalars should need to    C
!     be specified (e.g., mass fractions, void fraction, etc.,).       C
!     Such values will become defined according to their adjacent      C
!     fluid cell                                                       C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
               IF (BC_DT_0(BCV) == UNDEFINED) THEN 
                  IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 'BC_DT_0', BCV 
                   call mfix_exit(myPE)  
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

! If non-zero, check whether the velocity component through the boundary
! plane has the correct sign.
               SELECT CASE (TRIM(BC_PLANE(BCV)))  
               CASE ('W')  
                  IF (BC_U_G(BCV) < ZERO) THEN 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1050) BCV, &
                        'BC_U_g', '>'
                     CALL MFIX_EXIT(myPE) 
                  ENDIF 
               CASE ('E')  
                  IF (BC_U_G(BCV) > ZERO) THEN 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1050) BCV, &
                        'BC_U_g', '<' 
                     CALL MFIX_EXIT(myPE) 
                  ENDIF 
               CASE ('S')  
                  IF (BC_V_G(BCV) < ZERO) THEN 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1050) BCV, &
                        'BC_V_g', '>' 
                     CALL MFIX_EXIT(myPE) 
                  ENDIF 
               CASE ('N')  
                  IF (BC_V_G(BCV) > ZERO) THEN 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1050) BCV, &
                        'BC_V_g', '<' 
                     CALL MFIX_EXIT(myPE) 
                  ENDIF 
               CASE ('B')  
                  IF (BC_W_G(BCV) < ZERO) THEN 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1050) BCV, &
                        'BC_W_g', '>' 
                     CALL MFIX_EXIT(myPE) 
                  ENDIF 
               CASE ('T')  
                  IF (BC_W_G(BCV) > ZERO) THEN 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1050) BCV, &
                        'BC_W_g', '<' 
                     CALL MFIX_EXIT(myPE) 
                  ENDIF 
               END SELECT   ! end select case (bc_plane(bcv)) 

! DEM simulations do not employ variables for continuum solids. So do
! not perform checks on unnecessary data. Note, MPPIC initialization
! and BC's are based on continuum quantities so that these checks are
! necessary when invoking MPPIC.
               IF (.NOT.DISCRETE_ELEMENT .OR. DES_CONTINUUM_HYBRID &
                   .OR.MPPIC) THEN
                  DO M = 1, SMAX 
                     IF (BC_U_S(BCV,M) == UNDEFINED) THEN 
                        IF (BC_ROP_S(BCV,M)==ZERO .OR. NO_I) THEN 
                           BC_U_S(BCV,M) = ZERO 
                        ELSE 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1100) &
                              'BC_U_s', BCV, M 
                           call mfix_exit(myPE)  
                        ENDIF 
                     ENDIF 
                     IF (BC_V_S(BCV,M) == UNDEFINED) THEN 
                        IF (BC_ROP_S(BCV,M)==ZERO .OR. NO_J) THEN 
                           BC_V_S(BCV,M) = ZERO 
                        ELSE 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1100) &
                              'BC_V_s', BCV, M 
                           call mfix_exit(myPE)  
                        ENDIF 
                     ENDIF 
                     IF (BC_W_S(BCV,M) == UNDEFINED) THEN 
                        IF (BC_ROP_S(BCV,M)==ZERO .OR. NO_K) THEN 
                           BC_W_S(BCV,M) = ZERO 
                        ELSE 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1100) &
                              'BC_W_s', BCV, M 
                           call mfix_exit(myPE)  
                        ENDIF 
                     ENDIF 

! If non-zero, check whether the velocity component through the boundary
! plane has the correct sign.
                     SELECT CASE (TRIM(BC_PLANE(BCV)))  
                     CASE ('W')  
                        IF (BC_U_S(BCV,M) < ZERO) THEN 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1150) BCV, &
                              'BC_U_s', M, '>' 
                           CALL MFIX_EXIT(myPE) 
                        ENDIF 
                     CASE ('E')  
                        IF (BC_U_S(BCV,M) > ZERO) THEN 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1150) BCV, &
                              'BC_U_s', M, '<' 
                           CALL MFIX_EXIT(myPE) 
                        ENDIF 
                     CASE ('S')  
                        IF (BC_V_S(BCV,M) < ZERO) THEN 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1150) BCV, &
                              'BC_V_s', M, '>' 
                           CALL MFIX_EXIT(myPE) 
                        ENDIF 
                     CASE ('N')  
                        IF (BC_V_S(BCV,M) > ZERO) THEN 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1150) BCV, &
                              'BC_V_s', M, '<' 
                           CALL MFIX_EXIT(myPE) 
                        ENDIF 
                     CASE ('B')  
                        IF (BC_W_S(BCV,M) < ZERO) THEN 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1150) BCV, &
                              'BC_W_s', M, '>' 
                           CALL MFIX_EXIT(myPE) 
                        ENDIF 
                     CASE ('T')  
                        IF (BC_W_S(BCV,M) > ZERO) THEN 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1150) BCV, &
                              'BC_W_s', M, '<' 
                           CALL MFIX_EXIT(myPE) 
                        ENDIF 
                     END SELECT   ! end select case (bc_plane(bcv))
                  ENDDO   ! end loop over (m=1,smax)
               ENDIF   ! end if (.not.discrete_element .or.
                       !         des_continuum_hybrid .or. mppic)
! case mass_outflow
! --------------------------------------------<<<


            CASE ('P_INFLOW')  
! -------------------------------------------->>>
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
! Comments:                                                            C
!     Unlike the MI boundary, for the PI boundary the velocities at    C
!     the inflow face are calculated by solving the momentum eqns      C
!     and are not fixed. In this way, the PI is similar to the PO      C
!     except that the flow is into the domain and hence all other      C
!     scalars (e.g., mass fractions, void fraction, temperature,       C
!     etc.,) at the inflow cells need to be specified. To satisfy      C
!     the error routines at the start of the simulation, both the      C
!     tangential and normal components at the inflow also need to      C
!     be specified. The velocities values essentially serve as IC.     C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
               IF (BC_EP_G(BCV) == UNDEFINED) THEN 
                  IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 'BC_EP_g', BCV 
                  call mfix_exit(myPE)  
               ENDIF 
               IF (BC_P_G(BCV) == UNDEFINED) THEN 
                  IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 'BC_P_g', BCV 
                  call mfix_exit(myPE)  
               ELSEIF (BC_P_G(BCV)<=ZERO .AND. RO_G0==UNDEFINED) THEN 
                  IF(DMP_LOG)WRITE (UNIT_LOG, 1010) BCV, BC_P_G(BCV) 
                  call mfix_exit(myPE)  
               ENDIF 
               IF ((ENERGY_EQ .OR. RO_G0==UNDEFINED .OR. &
                    MU_G0==UNDEFINED) .AND. BC_T_G(BCV)==UNDEFINED) THEN
                  IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 'BC_T_g', BCV 
                  call mfix_exit(myPE)  
               ENDIF 

! if gas species equations are solved then check sum of the gas species
! mass fractions
               IF (SPECIES_EQ(0)) THEN 
                  SUM = ZERO 
! sum together those gas phase species mass fractions that are defined
                  DO N = 1, NMAX(0) 
                     IF (BC_X_G(BCV,N) /= UNDEFINED) &
                        SUM = SUM + BC_X_G(BCV,N) 
                  ENDDO 
! set any undefined species mass fractions to zero. Warn the user if 
! the sum of species mass fractions is not one and indicate which 
! species is undefined (mfix will exit in next check).           
                  DO N = 1, NMAX(0) 
                     IF (BC_X_G(BCV,N) == UNDEFINED) THEN 
                        IF (.NOT.COMPARE(ONE,SUM) .AND. DMP_LOG)&
                           WRITE (UNIT_LOG, 1060) BCV, N 
                        BC_X_G(BCV,N) = ZERO 
                     ENDIF 
                  ENDDO 
! if sum of the gas phase species mass fraction not 1...                  
                  IF (.NOT.COMPARE(ONE,SUM)) THEN 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1065) BCV 
                     call mfix_exit(myPE)  
                  ENDIF 
               ENDIF 

! DEM simulations do not employ variables for continuum solids. So do
! not perform checks on unnecessary data. Note, MPPIC initialization
! and BC's are based on continuum quantities so that these checks are
! necessary when invoking MPPIC.
               IF (.NOT.DISCRETE_ELEMENT .OR. DES_CONTINUUM_HYBRID &
                   .OR.MPPIC) THEN

! at this point bc_ep_g must be defined
                  SUM_EP = BC_EP_G(BCV) 
                  DO M = 1, SMAX 

! the following check on bc_rop_s is redundant with code in flow_to_vel,
! except that the latter has no calls for MFIX to exit
                     IF (BC_ROP_S(BCV,M) == UNDEFINED) THEN 
                        IF (BC_EP_G(BCV) == ONE) THEN 
                           BC_ROP_S(BCV,M) = ZERO 
                        ELSEIF (SMAX > 1 .OR. DES_CONTINUUM_HYBRID) THEN
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1100)  &
                              'BC_ROP_s', BCV, M 
                           call mfix_exit(myPE)
                        ENDIF
                     ENDIF

! Verify that species mass fractions have been specified.
                     SUM = ZERO
                     DO N = 1, NMAX(M)
! Sum together solids phase M species mass fractions that are defined.
                        IF(BC_X_S(BCV,M,N)/=UNDEFINED) &
                           SUM=SUM+BC_X_S(BCV,M,N) 
                     ENDDO
! If no solids phase M present and no species are defined for phase M 
! then set mass fraction of solids phase M species 1 to one
                     IF (SUM==ZERO .AND. BC_ROP_S(BCV,M)==ZERO) THEN 
                        BC_X_S(BCV,M,:) = ZERO
                        BC_X_S(BCV,M,1) = ONE 
                        SUM = ONE 
                     ENDIF

! Set the default if needed.
                     IF(NMAX(M)==1 .AND. BC_X_S(BCV,M,1)==UNDEFINED)THEN
                        BC_X_S(BCV,M,:) = ZERO
                        BC_X_S(BCV,M,1) = ONE
                        SUM = ONE
! Set any undefined species mass fractions to zero. Warn the user if 
! the sum of species mass fractions is not one and indicate which
! species is undefined (mfix will exit in next check). 
                     ELSE
                        DO N = 1, NMAX(M)
                           IF (BC_X_S(BCV,M,N) == UNDEFINED) THEN 
                              IF(.NOT.COMPARE(ONE,SUM) .AND. DMP_LOG)  &
                                 WRITE (UNIT_LOG, 1110)BCV, M, N 
                              BC_X_S(BCV,M,N) = ZERO 
                           ENDIF
                        ENDDO
                     ENDIF
! Verify that the sum of the mass fractions sum to one. Otherwise, exit.
                     IF(.NOT.COMPARE(ONE,SUM)) THEN 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1120) BCV, M 
                        CALL MFIX_EXIT(myPE)
                     ENDIF 

! Calculate the solid density.
                     BC_ROs = UNDEFINED
                     IF(SOLVE_ROs(M))THEN
                        INERT = INERT_SPECIES(M)
! Verify that the species mass fraction for the inert material is not
! zero in the BC region when the solids is present.
                        IF(BC_X_S(BCV,M,INERT) == ZERO) THEN
                           IF(BC_ROP_S(BCV,M) /= ZERO) THEN
                              IF(DMP_LOG) THEN
                                 WRITE(*,1400) M, BCV
                                 WRITE(UNIT_LOG,1400) M, BCV
                              ENDIF
                              CALL MFIX_EXIT(myPE)
                           ELSE
! If the solids isn't present, give it the baseline density.
                              BC_ROs = RO_S0(M)
                           ENDIF
                        ELSE
! Calculate the solids density.
                           BC_ROs = RO_S0(M) * X_s0(M,INERT) / &
                              BC_X_S(BCV,M,INERT)
                        ENDIF
                     ELSE
                        BC_ROs = RO_S0(M)
                     ENDIF

! Sanity check on solids phase density.
                     IF(BC_ROs <= ZERO .OR. BC_ROs==UNDEFINED) THEN
                        IF(DMP_LOG)THEN
                           WRITE(*,1500) M, BCV
                           WRITE(UNIT_LOG,1500) M, BCV
                        ENDIF
                        CALL MFIX_EXIT(myPE)
                     ENDIF

! at this point bc_rop_s must be defined
! sum of void fraction and solids volume fractions
                     SUM_EP = SUM_EP + BC_ROP_S(BCV,M)/BC_ROs

                     IF (ENERGY_EQ .AND. BC_T_S(BCV,M)==UNDEFINED) THEN
                        IF (BC_ROP_S(BCV,M) == ZERO) THEN 
                           BC_T_S(BCV,M) = BC_T_G(BCV) 
                        ELSE 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1100) &
                              'BC_T_s', BCV, M 
                           call mfix_exit(myPE)  
                        ENDIF 
                     ENDIF 


                  ENDDO   ! end loop over (m=1,smax)

! check sum of gas void fraction and all solids volume fractions
                  IF (.NOT.DES_CONTINUUM_HYBRID) THEN
                     IF (.NOT.COMPARE(ONE,SUM_EP)) THEN 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1125) BCV 
                        call mfix_exit(myPE)  
                     ENDIF 
                  ELSE
! sum_ep not necessarily one at this point since discrete particles
! present in hybrid model
                     IF (SUM_EP>ONE .OR. SUM_EP<ZERO) THEN 
                        IF(DMP_LOG) WRITE (UNIT_LOG, 1130) BCV 
                        call mfix_exit(myPE)
                     ENDIF
                  ENDIF
               ELSE   ! discrete_element and .not.des_continuum_hybrid
                      ! and .not.mppic
                  SUM_EP = BC_EP_G(BCV)
                  IF (SUM_EP>ONE .OR. SUM_EP<ZERO) THEN 
                     IF(DMP_LOG) WRITE (UNIT_LOG, 1130) BCV 
                     call mfix_exit(myPE)
                  ENDIF
               ENDIF   ! end if/else (.not.discrete_element .or.
                       !               des_continuum_hybrid .or. mppic)

               DO N = 1, NScalar
                  IF (BC_Scalar(BCV,N) == UNDEFINED) THEN 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1004) &
                        'BC_Scalar', BCV, N 
                     CALL MFIX_EXIT(myPE)
                  ENDIF 
               ENDDO


! Check that velocities are also specified. These are essentially used
! as initial conditions for the boundary region. If they are not
! specified then a deafult value is set here otherwise check_data_20 
! will complain and cause MFIX to exit.
               IF (BC_U_G(BCV) == UNDEFINED) THEN 
                  BC_U_G(BCV) = ZERO 
                  IF (.NOT.NO_I) THEN                   
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1011) 'BC_U_g', BCV 
                  ENDIF
               ENDIF 
               IF (BC_V_G(BCV) == UNDEFINED) THEN 
                  BC_V_G(BCV) = ZERO 
                  IF (.NOT.NO_J) THEN                   
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1011) 'BC_V_g', BCV 
                  ENDIF
               ENDIF 
               IF (BC_W_G(BCV) == UNDEFINED) THEN 
                  BC_W_G(BCV) = ZERO 
                  IF (.NOT.NO_K) THEN
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1011) 'BC_W_g', BCV 
                  ENDIF
               ENDIF  

! DEM simulations do not employ variables for continuum solids. So do
! not perform checks on unnecessary data. Note, MPPIC initialization
! and BC's are based on continuum quantities so that these checks are
! necessary when invoking MPPIC.
               IF (.NOT.DISCRETE_ELEMENT .OR. DES_CONTINUUM_HYBRID &
                   .OR.MPPIC) THEN
                  DO M = 1, SMAX 
                     IF (BC_U_S(BCV,M) == UNDEFINED) THEN 
                        BC_U_S(BCV,M) = ZERO
                        IF(BC_ROP_S(BCV,M) == UNDEFINED) THEN
! do nothing                                
                        ELSEIF (BC_ROP_S(BCV,M) >= ZERO .AND. .NOT.NO_I) THEN
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1111) &
                              'BC_U_s', BCV, M 
                        ENDIF 
                     ENDIF 
                     IF (BC_V_S(BCV,M) == UNDEFINED) THEN 
                        BC_V_S(BCV,M) = ZERO
                        IF (BC_ROP_S(BCV,M) == UNDEFINED) THEN 
                        ELSEIF (BC_ROP_S(BCV,M)>=ZERO .AND. .NOT.NO_J) THEN 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1111) &
                              'BC_V_s', BCV, M 
                        ENDIF 
                     ENDIF 
                     IF (BC_W_S(BCV,M) == UNDEFINED) THEN 
                        BC_W_S(BCV,M) = ZERO 
                        IF (BC_ROP_S(BCV,M) == UNDEFINED) THEN
                        ELSEIF (BC_ROP_S(BCV,M)>=ZERO .AND. .NOT.NO_K) THEN 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1111) &
                              'BC_W_s', BCV, M 
                        ENDIF 
                     ENDIF 
                  ENDDO   ! end loop over (m=1,smax)
               ENDIF   ! end if (.not.discrete_element .or.
                       !         des_continuum_hybrid .or. mppic)

! case p_inflow               
! --------------------------------------------<<<


            CASE ('P_OUTFLOW')  
! -------------------------------------------->>>

               IF (BC_P_G(BCV) == UNDEFINED) THEN 
                  IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 'BC_P_g', BCV 
                  call mfix_exit(myPE)  
               ELSEIF (BC_P_G(BCV)<=ZERO .AND. RO_G0==UNDEFINED) THEN 
                  IF(DMP_LOG)WRITE (UNIT_LOG, 1010) BCV, BC_P_G(BCV) 
                  call mfix_exit(myPE)  
               ENDIF 
! case p_outflow
! --------------------------------------------<<<


            CASE ('OUTFLOW')  
! ------------------------------------------->>>

               IF (.NOT.ERROR) THEN 
                  ERROR = .TRUE. 
               ELSE 
                  IF(DMP_LOG)WRITE (UNIT_LOG, 1160) BCV 
                  call mfix_exit(myPE)  
               ENDIF 
! case outflow
! -------------------------------------------<<<

            END SELECT   ! end select case (bc_type(bcv))

! it is unclear whether specifying bc_ep_g or bc_rop_s at the boundary
! is physical? Should it be prevented in the first place? Do their
! values matter (factor into the solution) or are they inconsequential? 
            IF (BC_TYPE(BCV) == 'OUTFLOW' .OR. &
                BC_TYPE(BCV) == 'MASS_OUTFLOW' .OR. &
                BC_TYPE(BCV) == 'P_OUTFLOW') THEN

! if bc_ep_g is defined at a PO, MO or O boundary, then the sum of ep_g
! and ep_s at the boundary may not equal one given the following code
! in the subroutine set_outflow (see code for details). therefore if
! bc_ep_g is is defined and any solids are present, provide the user
! with a warning.
               IF (BC_EP_G(BCV) /= UNDEFINED) THEN
                  SUM_EP = BC_EP_G(BCV)

                  write(*,*) 'sum_ep: ', sum_ep
! Unclear how the discrete solids volume fraction can be dictated at 
! the boundary, so it is currently prevented!
                  IF (DISCRETE_ELEMENT) THEN
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1101) &
                        BC_TYPE(BCV), BCV
                     call mfix_exit(myPE)
                  ELSE
! by this point the code has checked that no discrete solids are present
! otherwise the routine will have exited
                     DO M = 1, SMAX 
                        IF (BC_ROP_S(BCV,M) == UNDEFINED) THEN 
                           IF (BC_EP_G(BCV) == ONE) THEN 
! what does it mean to force the bulk density to zero at the
! boundary? (again does this value matter anyway)
                              BC_ROP_S(BCV,M) = ZERO 
                           ELSEIF (SMAX == 1 ) THEN 
! no discrete solids are present so a bulk density can be defined from 
! 1-bc_ep_g even for hybrid model
                              BC_ROP_S(BCV,M) = (ONE - BC_EP_G(BCV))*RO_S0(M)
                           ELSE
! bc_ep_g is defined but some bc_rop_s(m) are undefined.
! in this ep_p in the outflow boundary will be based on the user defined
! value of bc_ep_g, while rop_s would become based on the adjacent fluid
! cell. consequently, no check ensures the result is consistent with
! a requirement for ep_g+ep_s=1.
                              IF(DMP_LOG)WRITE (UNIT_LOG, 1102) &
                                 'BC_ROP_s', BCV, M 
                              call mfix_exit(myPE)
                           ENDIF
                        ENDIF  ! end if(bc_rop_s(bcv,m) == undefined)
! by this point bc_rop_s should either be defined or mfix exited
! therefore we can check that sum of void fraction and solids volume 
! fractions
                        SUM_EP = SUM_EP + BC_ROP_S(BCV,M)/RO_S0(M) 
                        write(*,*) 'M, sum_ep: ', M, sum_ep
                     ENDDO

                     IF (.NOT.COMPARE(ONE,SUM_EP)) THEN 
! if both variables are completely defined they should sum to 1.
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1125) BCV 
                           call mfix_exit(myPE)  
                     ENDIF

                  ENDIF   ! end if/else (.not.discrete_element)

               ENDIF   ! end if (bc_ep_g(bcv) /= undefined)

            ENDIF   ! end if (bc_type(bcv)== 'mass_outflow', 'outflow'
                    ! or 'p_outflow')


! elseif branch from if (bc_defined(bcv))
         ELSEIF ( (BC_TYPE(BCV) /= 'DUMMY') .AND. & 
                  (BC_TYPE(BCV)(1:2) /= 'CG') ) THEN    

! Check whether BC values are specified for undefined BC locations
            IF (BC_U_G(BCV) /= UNDEFINED) THEN 
               IF(DMP_LOG)WRITE (UNIT_LOG, 1200) 'BC_U_g', BCV 
               call mfix_exit(myPE)  
            ENDIF 
            IF (BC_V_G(BCV) /= UNDEFINED) THEN 
               IF(DMP_LOG)WRITE (UNIT_LOG, 1200) 'BC_V_g', BCV 
               call mfix_exit(myPE)  
            ENDIF 
            IF (BC_W_G(BCV) /= UNDEFINED) THEN 
               IF(DMP_LOG)WRITE (UNIT_LOG, 1200) 'BC_W_g', BCV 
               call mfix_exit(myPE)  
            ENDIF 
            IF (BC_EP_G(BCV) /= UNDEFINED) THEN 
               IF(DMP_LOG)WRITE (UNIT_LOG, 1200) 'BC_EP_g', BCV 
               call mfix_exit(myPE)  
            ENDIF 
            IF (BC_P_G(BCV) /= UNDEFINED) THEN 
               IF(DMP_LOG)WRITE (UNIT_LOG, 1200) 'BC_P_g', BCV 
               call mfix_exit(myPE)  
            ENDIF 
            IF (BC_T_G(BCV) /= UNDEFINED) THEN 
               IF(DMP_LOG)WRITE (UNIT_LOG, 1200) 'BC_T_g', BCV 
               call mfix_exit(myPE)  
            ENDIF 
            DO N = 1, DIMENSION_N_G 
               IF (BC_X_G(BCV,N) /= UNDEFINED) THEN 
                  IF(DMP_LOG)WRITE (UNIT_LOG, 1200) 'X_g', BCV 
                  call mfix_exit(myPE)  
               ENDIF 
            ENDDO 

            IF (.NOT.DISCRETE_ELEMENT .OR. DES_CONTINUUM_HYBRID &
                .OR. MPPIC) THEN
               DO M = 1, SMAX 
                  IF (BC_ROP_S(BCV,M) /= UNDEFINED) THEN 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1300) &
                        'BC_ROP_s', BCV, M 
                     call mfix_exit(myPE)  
                  ENDIF 
                  DO N = 1, DIMENSION_N_S 
                     IF (BC_X_S(BCV,M,N) /= UNDEFINED) THEN 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1300) &
                           'BC_X_s', BCV, M 
                        call mfix_exit(myPE)  
                     ENDIF 
                  ENDDO 
                  IF (BC_U_S(BCV,M) /= UNDEFINED) THEN 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1300) 'BC_U_s', BCV, M 
                     call mfix_exit(myPE)  
                  ENDIF 
                  IF (BC_V_S(BCV,M) /= UNDEFINED) THEN 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1300) 'BC_V_s', BCV, M 
                     call mfix_exit(myPE)  
                  ENDIF 
                  IF (BC_W_S(BCV,M) /= UNDEFINED) THEN 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1300) 'BC_W_s', BCV, M 
                     call mfix_exit(myPE)  
                  ENDIF 
                  IF (BC_T_S(BCV,M) /= UNDEFINED) THEN 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1300) 'BC_T_s', BCV, M 
                     call mfix_exit(myPE)  
                  ENDIF 
               ENDDO 
            ENDIF   ! end if (.not.discrete_element)
    
            DO N = 1, NScalar
               IF (BC_Scalar(BCV,N) /= UNDEFINED) THEN 
                  IF(DMP_LOG)WRITE (UNIT_LOG, 1200) 'BC_Scalar', BCV
                  CALL MFIX_EXIT(myPE)
               ENDIF 
            ENDDO 
  
         ENDIF  ! end if (bc_defined(bcv)) 
                !     elseif (bc_type(bcv) /= 'dummy' .and.
                !             bc_type(bcv)(1:2) /= 'cg')
      ENDDO    ! end loop over (bcv=1,dimension_bc)

! END Check quantities for flow type boundaries
! ----------------------------------------------------------------<<<



! Following section to check quantities for wall type boundaries
! ---------------------------------------------------------------->>>      
      DO BCV = 1, DIMENSION_BC 
         IF (BC_DEFINED(BCV)) THEN 

! Default specification of Johnson-Jackson bc
            IF (GRANULAR_ENERGY) THEN 
               IF (BC_JJ_PS(BCV) == UNDEFINED_I) BC_JJ_PS(BCV) = 1 
            ELSE 
               IF (BC_JJ_PS(BCV) == UNDEFINED_I) BC_JJ_PS(BCV) = 0 
            ENDIF 

! Set Jenkins default specification, modify BC_JJ accordingly
            IF (GRANULAR_ENERGY .AND. JENKINS) BC_JJ_PS(BCV) = 1

            IF (BC_TYPE(BCV)=='FREE_SLIP_WALL' .OR. &
                BC_TYPE(BCV)=='NO_SLIP_WALL' .OR. &
                BC_TYPE(BCV)=='PAR_SLIP_WALL') THEN 

! The wall velocities are not needed for no-slip or free-slip                
               IF (BC_TYPE(BCV) == 'PAR_SLIP_WALL') THEN 
                  IF (BC_UW_G(BCV) == UNDEFINED) THEN 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 'BC_Uw_g', BCV 
                     call mfix_exit(myPE)  
                  ENDIF 
                  IF (BC_VW_G(BCV) == UNDEFINED) THEN 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 'BC_Vw_g', BCV 
                     call mfix_exit(myPE)  
                  ENDIF 
                  IF (.NOT.NO_K) THEN 
                     IF (BC_WW_G(BCV) == UNDEFINED) THEN 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 'BC_Ww_g', BCV
                        call mfix_exit(myPE)  
                     ENDIF 
                  ELSE 
                     BC_WW_G(BCV) = ZERO 
                  ENDIF 
               ENDIF   ! end if (bc_type(bcv)=='par_slip_wall')

               IF (.NOT.DISCRETE_ELEMENT .OR. DES_CONTINUUM_HYBRID &
                   .OR.MPPIC) THEN
                  IF (BC_TYPE(BCV)=='PAR_SLIP_WALL' .OR. BC_JJ_PS(BCV)==1) THEN
                     DO M = 1, MMAX
                        IF (BC_UW_S(BCV,M) == UNDEFINED) THEN 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1100) 'BC_Uw_s', BCV, M 
                           call mfix_exit(myPE)  
                        ENDIF 
                        IF (BC_VW_S(BCV,M) == UNDEFINED) THEN 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1100) 'BC_Vw_s', BCV, M 
                           call mfix_exit(myPE)  
                        ENDIF 
                        IF (.NOT.NO_K) THEN 
                           IF (BC_WW_S(BCV,M) == UNDEFINED) THEN 
                             IF(DMP_LOG)WRITE (UNIT_LOG, 1100) 'BC_Ww_s', BCV, M 
                             call mfix_exit(myPE)  
                           ENDIF 
                        ELSE 
                           BC_WW_S(BCV,M) = ZERO 
                        ENDIF 
                     ENDDO 
                  ENDIF   ! end if (bc_type(bcv)=='par_slip_wall' or bc_jj_ps(bcv)==1)
               ENDIF   ! end if (.not.discrete_element .or. des_continuum_hybrid .or.
                       !          mppic)

               IF (ENERGY_EQ) THEN 
                  IF (BC_HW_T_G(BCV) < ZERO) THEN 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1003) 'BC_hw_T_g', BCV 
                     call mfix_exit(myPE)  
                  ENDIF 
                  IF (BC_HW_T_G(BCV)/=ZERO .AND. &
                      BC_TW_G(BCV)==UNDEFINED) THEN 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 'BC_Tw_g', BCV 
                     call mfix_exit(myPE)  
                  ENDIF 
                  IF (BC_HW_T_G(BCV)/=UNDEFINED .AND. &
                      BC_C_T_G(BCV)==UNDEFINED) THEN 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 'BC_C_T_g', BCV 
                     call mfix_exit(myPE)  
                  ENDIF 

                  IF (.NOT.DISCRETE_ELEMENT .OR. DES_CONTINUUM_HYBRID &
                      .OR. MPPIC) THEN
                     DO M = 1, SMAX 
                        IF (BC_HW_T_S(BCV,M) < ZERO) THEN 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1103) 'BC_hw_T_s', BCV, M 
                           call mfix_exit(myPE)  
                        ENDIF 
                        IF (BC_HW_T_S(BCV,M)/=ZERO .AND.&
                            BC_TW_S(BCV,M)==UNDEFINED) THEN 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1100) 'BC_Tw_s', BCV, M 
                           call mfix_exit(myPE)  
                        ENDIF 
                        IF (BC_HW_T_S(BCV,M)/=UNDEFINED .AND. &
                            BC_C_T_S(BCV,M)==UNDEFINED) THEN 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1100) 'BC_C_T_s', BCV, M 
                           call mfix_exit(myPE)  
                        ENDIF 
                     ENDDO 
                  ENDIF   ! end if (.not.discrete_element .or. des_continuum_hybrid)
               ENDIF   ! end if (energy_eq)

               IF (.NOT. DISCRETE_ELEMENT .OR. DES_CONTINUUM_HYBRID  &
                   .OR. MPPIC) THEN
                  IF (GRANULAR_ENERGY .AND. BC_JJ_PS(BCV)==0) THEN 
                     DO M = 1, MMAX 
                        IF (BC_HW_THETA_M(BCV,M) < ZERO) THEN 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1103) &
                              'BC_hw_Theta_m', BCV, M 
                           call mfix_exit(myPE)  
                        ENDIF 
                        IF (BC_HW_THETA_M(BCV,M)/=ZERO .AND. &
                            BC_THETAW_M(BCV,M)==UNDEFINED) THEN 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1100) &
                              'BC_Thetaw_m', BCV, M 
                           call mfix_exit(myPE)  
                        ENDIF 
                        IF (BC_HW_THETA_M(BCV,M)/=UNDEFINED .AND. &
                            BC_C_THETA_M(BCV,M)==UNDEFINED) THEN 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1100) &
                              'BC_C_Theta_m', BCV, M 
                           call mfix_exit(myPE)  
                        ENDIF 
                     ENDDO 
                  ENDIF   ! end if (granular_energy .and. bc_jj_ps(bcv)==0)
               ENDIF    ! end if (.not.discrete_element .or. des_continuum_hybrid 
                        !         .or. mppic)

               IF (SPECIES_EQ(0)) THEN 
                  DO N = 1, NMAX(0) 
                     IF (BC_HW_X_G(BCV,N) < ZERO) THEN 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1005) 'BC_hw_X_g', BCV, N 
                        call mfix_exit(myPE)  
                     ENDIF 
                     IF (BC_HW_X_G(BCV,N)/=ZERO .AND. &
                         BC_XW_G(BCV,N)==UNDEFINED) THEN 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1004) 'BC_Xw_g', BCV, N 
                        call mfix_exit(myPE)  
                     ENDIF 
                     IF (BC_HW_X_G(BCV,N)/=UNDEFINED .AND. &
                         BC_C_X_G(BCV,N)==UNDEFINED) THEN 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1004) 'BC_C_X_g', BCV, N 
                        call mfix_exit(myPE)  
                     ENDIF 
                  ENDDO 
               ENDIF   ! end if (species_eq(0))

               IF (.NOT.DISCRETE_ELEMENT .OR. DES_CONTINUUM_HYBRID .OR.&
                   MPPIC) THEN
                  DO M = 1, SMAX 
                     IF (SPECIES_EQ(M)) THEN 
                        DO N = 1, NMAX(M) 
                           IF (BC_HW_X_S(BCV,M,N) < ZERO) THEN 
                              IF(DMP_LOG)WRITE (UNIT_LOG, 1105) &
                                 'BC_hw_X_s', BCV, M, N 
                              call mfix_exit(myPE)  
                           ENDIF 
                           IF (BC_HW_X_S(BCV,M,N)/=ZERO .AND. &
                               BC_XW_S(BCV,M,N)==UNDEFINED) THEN 
                              IF(DMP_LOG)WRITE (UNIT_LOG, 1104) &
                                 'BC_Xw_s', BCV, M, N 
                              call mfix_exit(myPE)  
                           ENDIF 
                           IF (BC_HW_X_S(BCV,M,N)/=UNDEFINED .AND. &
                               BC_C_X_S(BCV,M,N)==UNDEFINED) THEN 
                              IF(DMP_LOG)WRITE (UNIT_LOG, 1104) &
                                 'BC_C_X_s', BCV, M, N 
                              call mfix_exit(myPE)  
                           ENDIF 
                        ENDDO 
                     ENDIF    ! end if (species_eq(m))
                  ENDDO    ! end loop over (m=1,smax)
               ENDIF   ! end if (.not.discrete_element .or. des_continuum_hybrid)
                       !         .or. mppic)

               DO N = 1, NScalar
                  IF (BC_HW_Scalar(BCV,N) < ZERO) THEN 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1005) 'BC_hw_Scalar', BCV, N 
                     CALL MFIX_EXIT(myPE)
                  ENDIF 
                  IF (BC_HW_Scalar(BCV,N)/=ZERO .AND. &
                      BC_Scalarw(BCV,N)==UNDEFINED) THEN 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1004) 'BC_ScalarW', BCV, N 
                     CALL MFIX_EXIT(myPE)
                  ENDIF 
                  IF (BC_HW_Scalar(BCV,N)/=UNDEFINED .AND. &
                      BC_C_Scalar(BCV,N)==UNDEFINED) THEN 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1004) 'BC_C_Scalar', BCV, N 
                     CALL MFIX_EXIT(myPE)
                  ENDIF 
               ENDDO    ! end loop over (n=1,nscalar)

            ENDIF  ! end if (bc_type(bcv)=='free_slip_wall' .or. 
                   !        bc_type(bcv)=='no_slip_wall' .or.
                   !        bc_type(bcv)=='par_slip_wall')

         ENDIF    ! end if (bc_defined(bcv))
      ENDDO    ! end loop over (bcv = 1,dimension_bc)

! END Check quantities for wall type boundaries
! ----------------------------------------------------------------<<<


      IF (RUN_TYPE(1:3) /= 'NEW') RETURN  
      ERROR = .FALSE. 

      DO K = kstart2, kend2 
         DO J = jstart2, Jend2 
            DO I = istart2, Iend2
               IJK = FUNIJK(I,J,K) 
               IF (ICBC_FLAG(IJK) == '   ') THEN 
                  IF (.NOT.ERROR .AND. DMP_LOG)WRITE (UNIT_LOG, 1400) 
                  IF(DMP_LOG)WRITE (UNIT_LOG, 1410) I, J, K 
                  ERROR = .TRUE. 
               ENDIF 
            ENDDO 
         ENDDO 
      ENDDO 

      call send_recv(icbc_flag,2)

      IF (ERROR) THEN 
         IF(DMP_LOG)WRITE (UNIT_LOG, 1420) 
         call mfix_exit(myPE)  
      ENDIF 

      RETURN  

 1000 FORMAT(/1X,70('*')//' From: CHECK_DATA_07',/' Message: ',A,'(',I2,&
         ') not specified',/1X,70('*')/) 
 1001 FORMAT(/1X,70('*')//' From: CHECK_DATA_07',/&
         ' Message: Illegal BC_TYPE for BC # = ',I2,/'   BC_TYPE = ',A,/&
         '  Valid BC_TYPE are: ') 
 1002 FORMAT(5X,A16) 
 1003 FORMAT(/1X,70('*')//' From: CHECK_DATA_07',/' Message: ',A,'(',I2,&
         ') value is unphysical',/1X,70('*')/) 
 1004 FORMAT(/1X,70('*')//' From: CHECK_DATA_07',/' Message: ',A,'(',I2,',',I2,&
         ') not specified',/1X,70('*')/) 
 1005 FORMAT(/1X,70('*')//' From: CHECK_DATA_07',/' Message: ',A,'(',I2,',',I2,&
         ') value is unphysical',/1X,70('*')/) 
 1010 FORMAT(/1X,70('*')//' From: CHECK_DATA_07',/' Message: BC_P_g( ',I2,&
         ') = ',G12.5,/&
         ' Pressure should be greater than zero for compressible flow',/1X,70(&
         '*')/) 
 1011 FORMAT(/1X,70('*')//' From: CHECK_DATA_07',/' Message: ',A,'(',I2,&
         ') not specified.',/1X,'These serve as initial values in ',&
         'the boundary region. Set to 0 by default',/1X,70('*')/) 
 1050 FORMAT(/1X,70('*')//' From: CHECK_DATA_07',/' Message: BC number:',I2,&
         ' - ',A,' should be ',A,' zero.',/1X,70('*')/) 
 1060 FORMAT(/1X,70('*')//' From: CHECK_DATA_07',/' Message: BC_X_g(',I2,',',I2&
         ,') not specified',/1X,70('*')/) 
 1065 FORMAT(/1X,70('*')//' From: CHECK_DATA_07',/' Message: BC number:',I2,&
         ' - Sum of gas mass fractions is NOT equal to one',/1X,70('*')/) 
 1100 FORMAT(/1X,70('*')//' From: CHECK_DATA_07',/' Message: ',A,'(',I2,',',I1,&
         ') not specified',/1X,70('*')/) 
 1101 FORMAT(/1X,70('*')//' From: CHECK_DATA_07',/' Message: For ', &
         'BC_TYPE= ', A, ' BC_EP_G(',I2,') should not be defined',/1X,&
         'or the sum of the volume fractions may not equal one',&
         /1X,70('*')/)
 1102 FORMAT(/1X,70('*')//' From: CHECK_DATA_07',/' Message: For ',& 
         'BC_TYPE= ', A, ' Since BC_EP_G('I2,') is defined,',/1X,&
         'BC_ROP_S should also be defined to ensure that the ',&
         'sum of their',/1X, 'volume fractions will equal one',&
         /1X,70('*')/)

 1103 FORMAT(/1X,70('*')//' From: CHECK_DATA_07',/' Message: ',A,'(',I2,',',I1,&
         ') value is unphysical',/1X,70('*')/) 
 1104 FORMAT(/1X,70('*')//' From: CHECK_DATA_07',/' Message: ',A,'(',I2,',',I2,&
         ',',I2,') not specified',/1X,70('*')/) 
 1105 FORMAT(/1X,70('*')//' From: CHECK_DATA_07',/' Message: ',A,'(',I2,',',I2,&
         ',',I2,') value is unphysical',/1X,70('*')/) 
 1110 FORMAT(/1X,70('*')//' From: CHECK_DATA_07',/' Message: BC_X_s(',I2,',',I2&
         ,',',I2,') not specified',/1X,70('*')/) 
 1111 FORMAT(/1X,70('*')//' From: CHECK_DATA_07',/' Message: ',A,'(',I2,',',I1,&
         ') not specified.',/1X,'These serve as initial values in ',&
         'the boundary region. Set to 0 by default',/1X,70('*')/) 
 1120 FORMAT(/1X,70('*')//' From: CHECK_DATA_07',/' Message: BC number:',I2,&
         ' - Sum of solids-',I1,' mass fractions is NOT equal to one',/1X,70(&
         '*')/) 
 1125 FORMAT(/1X,70('*')//' From: CHECK_DATA_07',/' Message: BC number:',I2,&
         ' - Sum of volume fractions is NOT equal to one',/1X,70('*')/) 
 1130 FORMAT(/1X,70('*')//' From: CHECK_DATA_07',/' Message: BC number:',I2,&
         ' - void fraction is unphysical (>1 or <0)',/1X,70('*')/) 
 1150 FORMAT(/1X,70('*')//' From: CHECK_DATA_07',/' Message: BC number:',I2,&
         ' - ',A,I1,' should be ',A,' zero.',/1X,70('*')/) 
 1160 FORMAT(/1X,70('*')//' From: CHECK_DATA_07',/&
         ' Message: Boundary condition no', &
         I2,' is a second outflow condition.',/1X,&
         '  Only one outflow is allowed.  Consider using P_OUTFLOW.',/1X, 70('*')/) 
 1200 FORMAT(/1X,70('*')//' From: CHECK_DATA_07',/' Message: ',A,'(',I2,&
         ') specified',' for an undefined BC location',/1X,70('*')/) 
 1300 FORMAT(/1X,70('*')//' From: CHECK_DATA_07',/' Message: ',A,'(',I2,',',I1,&
         ') specified',' for an undefined BC location',/1X,70('*')/) 
 1400 FORMAT(/1X,70('*')//' From: CHECK_DATA_07',/&
         ' Message: No initial or boundary condition specified',/&
         '    I       J       K') 
 1410 FORMAT(I5,3X,I5,3X,I5) 
 1420 FORMAT(/1X,70('*')/) 

 1500 FORMAT(//1X,70('*')/' From: CHECK_DATA_07',/,' Error 1500:'      &
         ' Solids phase ',I2,' failed sanity check in IC region ',I3,  &
         '. ',/' Please check mfix.dat file.',/1X,70('*')//)


      END SUBROUTINE CHECK_DATA_07 


