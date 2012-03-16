!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CHECK_DATA_07                                          C
!  Purpose: Check boundary condition specifications, convert           C
!           physical locations to i, j, k's, compute BC_AREA           C
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
!  Variables referenced: BC_X_w, BC_X_e, BC_Y_s, BC_Y_n, BC_Z_b        C
!                        BC_Z_t, BC_I_w, BC_I_e, BC_J_s, BC_J_n        C
!                        BC_K_b, BC_K_t, IMAX2, JMAX2, KMAX2           C
!  Variables modified: BC_DEFINED                                      C
!                                                                      C
!  Local variables: BCV, I, J, K, IJK, VALID_BC_TYPE                  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CHECK_DATA_07 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
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
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
      INTEGER, PARAMETER :: DIM_BCTYPE = 21 
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!
!     error flag
      LOGICAL ERROR,RECOGNIZED_BC_TYPE
!
!     loop/variable indices
      INTEGER BCV , I , J , K , IJK, M, N
!
!     valid boundary condition types
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
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
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================

      DOUBLE PRECISION SUM, SUM_EP

!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      LOGICAL , EXTERNAL :: COMPARE 
!-----------------------------------------------
      INCLUDE 'function.inc'


!     DETERMINE WHICH BOUNDARY CONDITION INDICES HAVE VALUES
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
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
        IF (BC_TYPE(BCV) == 'CG_NSW') BC_DEFINED(BCV) = .FALSE. 
        IF (BC_TYPE(BCV) == 'CG_FSW') BC_DEFINED(BCV) = .FALSE. 
        IF (BC_TYPE(BCV) == 'CG_PSW') BC_DEFINED(BCV) = .FALSE. 
        IF (BC_TYPE(BCV) == 'CG_MI') BC_DEFINED(BCV) = .FALSE. 
        IF (BC_TYPE(BCV) == 'CG_PO') BC_DEFINED(BCV) = .FALSE. 

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
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
        IF (BC_DEFINED(BCV)) THEN 
            IF (BC_X_W(BCV)==UNDEFINED .AND. BC_I_W(BCV)==UNDEFINED_I) THEN 
                IF (NO_I) THEN 
                    BC_X_W(BCV) = ZERO
                ELSE 
                    IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 'BC_X_w and BC_I_w ', BCV
                    call mfix_exit(myPE)
                ENDIF 
            ENDIF 
            IF (BC_X_E(BCV)==UNDEFINED .AND. BC_I_E(BCV)==UNDEFINED_I) THEN 
                IF (NO_I) THEN 
                    BC_X_E(BCV) = XLENGTH 
                ELSE 
                    IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 'BC_X_e and BC_I_e ', BCV
                    call mfix_exit(myPE)
                ENDIF 
            ENDIF 
            IF (BC_Y_S(BCV)==UNDEFINED .AND. BC_J_S(BCV)==UNDEFINED_I) THEN 
                IF (NO_J) THEN 
                    BC_Y_S(BCV) = ZERO 
                ELSE 
                    IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 'BC_Y_s and BC_J_s ', BCV
                    call mfix_exit(myPE)
                ENDIF 
            ENDIF 
            IF (BC_Y_N(BCV)==UNDEFINED .AND. BC_J_N(BCV)==UNDEFINED_I) THEN 
                IF (NO_J) THEN 
                    BC_Y_N(BCV) = YLENGTH 
                ELSE 
                    IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 'BC_Y_n and BC_J_n ', BCV 
                    call mfix_exit(myPE)
                ENDIF 
            ENDIF 
            IF (BC_Z_B(BCV)==UNDEFINED .AND. BC_K_B(BCV)==UNDEFINED_I) THEN 
                IF (NO_K) THEN 
                    BC_Z_B(BCV) = ZERO 
                ELSE 
                    IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 'BC_Z_b and BC_K_b ', BCV 
                    call mfix_exit(myPE)
                ENDIF 
            ENDIF 
            IF (BC_Z_T(BCV)==UNDEFINED .AND. BC_K_T(BCV)==UNDEFINED_I) THEN 
                IF (NO_K) THEN 
                    BC_Z_T(BCV) = ZLENGTH 
                ELSE 
                   IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 'BC_Z_t and BC_K_t ', BCV 
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
        ENDIF 
      ENDDO L50 


      CALL GET_WALLS_BC 

!     Find and validate i, j, k locations of flow BC's
      CALL GET_FLOW_BC 

!     Compute area of boundary surfaces
      CALL GET_BC_AREA 

      DO BCV = 1, DIMENSION_BC 
        IF (BC_DEFINED(BCV)) THEN 
!           Check the specification of physical quantities needed in FLOW_TO_VEL
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
                IF ((ENERGY_EQ .OR. RO_G0==UNDEFINED .OR. MU_G0==UNDEFINED)&
                .AND. BC_T_G(BCV)==UNDEFINED) THEN 
                    IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 'BC_T_g', BCV 
                    call mfix_exit(myPE)
                ENDIF 
                IF (SPECIES_EQ(0) .OR. RO_G0==UNDEFINED .AND. MW_AVG==UNDEFINED&
                ) THEN 
                    SUM = ZERO 
                    DO N = 1, NMAX(0) 
                        IF (BC_X_G(BCV,N) /= UNDEFINED) SUM = SUM + BC_X_G(BCV,N) 
                    ENDDO 
                    DO N = 1, NMAX(0) 
                        IF (BC_X_G(BCV,N) == UNDEFINED) THEN 
                            IF (.NOT.COMPARE(ONE,SUM) .AND. DMP_LOG)WRITE (UNIT_LOG, 1060) BCV, N 
                            BC_X_G(BCV,N) = ZERO 
                        ENDIF 
                    ENDDO 
                    IF (.NOT.COMPARE(ONE,SUM)) THEN 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1065) BCV 
                        call mfix_exit(myPE)
                    ENDIF 
                ENDIF 
            ENDIF 
        ENDIF 
      ENDDO 

      CALL FLOW_TO_VEL 

      ERROR = .FALSE. 
      DO BCV = 1, DIMENSION_BC 
        IF (BC_DEFINED(BCV)) THEN 

!           Check the specification of physical quantities for inflow cells.
            SELECT CASE (TRIM(BC_TYPE(BCV)))
            CASE ('MASS_INFLOW')  
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
                SELECT CASE (BC_PLANE(BCV))  
                CASE ('W')  
                    IF (BC_U_G(BCV) > ZERO) THEN 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1050) BCV, 'BC_U_g', '<' 
                        CALL MFIX_EXIT(myPE) 
                    ENDIF 
                CASE ('E')  
                    IF (BC_U_G(BCV) < ZERO) THEN 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1050) BCV, 'BC_U_g', '>' 
                        CALL MFIX_EXIT(myPE) 
                    ENDIF 
                CASE ('S')  
                    IF (BC_V_G(BCV) > ZERO) THEN 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1050) BCV, 'BC_V_g', '<' 
                        CALL MFIX_EXIT(myPE) 
                    ENDIF 
                CASE ('N')  
                    IF (BC_V_G(BCV) < ZERO) THEN 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1050) BCV, 'BC_V_g', '>' 
                        CALL MFIX_EXIT(myPE) 
                    ENDIF 
                CASE ('B')  
                    IF (BC_W_G(BCV) > ZERO) THEN 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1050) BCV, 'BC_W_g', '<' 
                        CALL MFIX_EXIT(myPE) 
                    ENDIF 
                CASE ('T')  
                    IF (BC_W_G(BCV) < ZERO) THEN 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1050) BCV, 'BC_W_g', '>' 
                        CALL MFIX_EXIT(myPE) 
                    ENDIF 
                END SELECT 

                SUM_EP = BC_EP_G(BCV) 
                DO M = 1, SMAX 
                    IF (BC_ROP_S(BCV,M) == UNDEFINED) THEN 
                        IF (BC_EP_G(BCV) == ONE) THEN 
                            BC_ROP_S(BCV,M) = ZERO 
                        ELSEIF (SMAX == 1) THEN 
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
                    SELECT CASE (TRIM(BC_PLANE(BCV)))  
                    CASE ('W')  
                        IF (BC_U_S(BCV,M) > ZERO) THEN 
                            IF(DMP_LOG)WRITE (UNIT_LOG, 1150) BCV, 'BC_U_s', M, '<' 
                            CALL MFIX_EXIT(myPE) 
                        ENDIF 
                    CASE ('E')  
                        IF (BC_U_S(BCV,M) < ZERO) THEN 
                            IF(DMP_LOG)WRITE (UNIT_LOG, 1150) BCV, 'BC_U_s', M, '>' 
                            CALL MFIX_EXIT(myPE) 
                        ENDIF 
                    CASE ('S')  
                        IF (BC_V_S(BCV,M) > ZERO) THEN 
                            IF(DMP_LOG)WRITE (UNIT_LOG, 1150) BCV, 'BC_V_s', M, '<' 
                            CALL MFIX_EXIT(myPE) 
                        ENDIF 
                    CASE ('N')  
                        IF (BC_V_S(BCV,M) < ZERO) THEN 
                            IF(DMP_LOG)WRITE (UNIT_LOG, 1150) BCV, 'BC_V_s', M, '>' 
                            CALL MFIX_EXIT(myPE) 
                        ENDIF 
                    CASE ('B')  
                        IF (BC_W_S(BCV,M) > ZERO) THEN 
                            IF(DMP_LOG)WRITE (UNIT_LOG, 1150) BCV, 'BC_W_s', M, '<' 
                            CALL MFIX_EXIT(myPE) 
                        ENDIF 
                    CASE ('T')  
                        IF (BC_W_S(BCV,M) < ZERO) THEN 
                            IF(DMP_LOG)WRITE (UNIT_LOG, 1150) BCV, 'BC_W_s', M, '>' 
                            CALL MFIX_EXIT(myPE) 
                        ENDIF 
                    END SELECT 
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

            CASE ('MASS_OUTFLOW')  
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

!               Check whether the bc velocity components have the correct sign
                SELECT CASE (TRIM(BC_PLANE(BCV)))  
                CASE ('W')  
                    IF (BC_U_G(BCV) < ZERO) THEN 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1050) BCV, 'BC_U_g', '>' 
                        CALL MFIX_EXIT(myPE) 
                    ENDIF 
                CASE ('E')  
                    IF (BC_U_G(BCV) > ZERO) THEN 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1050) BCV, 'BC_U_g', '<' 
                        CALL MFIX_EXIT(myPE) 
                    ENDIF 
                CASE ('S')  
                    IF (BC_V_G(BCV) < ZERO) THEN 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1050) BCV, 'BC_V_g', '>' 
                        CALL MFIX_EXIT(myPE) 
                    ENDIF 
                CASE ('N')  
                    IF (BC_V_G(BCV) > ZERO) THEN 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1050) BCV, 'BC_V_g', '<' 
                        CALL MFIX_EXIT(myPE) 
                    ENDIF 
                CASE ('B')  
                    IF (BC_W_G(BCV) < ZERO) THEN 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1050) BCV, 'BC_W_g', '>' 
                        CALL MFIX_EXIT(myPE) 
                    ENDIF 
                CASE ('T')  
                    IF (BC_W_G(BCV) > ZERO) THEN 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1050) BCV, 'BC_W_g', '<' 
                        CALL MFIX_EXIT(myPE) 
                    ENDIF 
                END SELECT 

                DO M = 1, SMAX 
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

!                   Check whether the bc velocity components have the correct sign
                    SELECT CASE (TRIM(BC_PLANE(BCV)))  
                    CASE ('W')  
                        IF (BC_U_S(BCV,M) < ZERO) THEN 
                            IF(DMP_LOG)WRITE (UNIT_LOG, 1150) BCV, 'BC_U_s', M, '>' 
                            CALL MFIX_EXIT(myPE) 
                        ENDIF 
                    CASE ('E')  
                        IF (BC_U_S(BCV,M) > ZERO) THEN 
                            IF(DMP_LOG)WRITE (UNIT_LOG, 1150) BCV, 'BC_U_s', M, '<' 
                            CALL MFIX_EXIT(myPE) 
                        ENDIF 
                    CASE ('S')  
                        IF (BC_V_S(BCV,M) < ZERO) THEN 
                            IF(DMP_LOG)WRITE (UNIT_LOG, 1150) BCV, 'BC_V_s', M, '>' 
                            CALL MFIX_EXIT(myPE) 
                        ENDIF 
                    CASE ('N')  
                        IF (BC_V_S(BCV,M) > ZERO) THEN 
                            IF(DMP_LOG)WRITE (UNIT_LOG, 1150) BCV, 'BC_V_s', M, '<' 
                            CALL MFIX_EXIT(myPE) 
                        ENDIF 
                    CASE ('B')  
                        IF (BC_W_S(BCV,M) < ZERO) THEN 
                            IF(DMP_LOG)WRITE (UNIT_LOG, 1150) BCV, 'BC_W_s', M, '>' 
                            CALL MFIX_EXIT(myPE) 
                        ENDIF 
                    CASE ('T')  
                        IF (BC_W_S(BCV,M) > ZERO) THEN 
                            IF(DMP_LOG)WRITE (UNIT_LOG, 1150) BCV, 'BC_W_s', M, '<' 
                            CALL MFIX_EXIT(myPE) 
                        ENDIF 
                    END SELECT 
                ENDDO 

            CASE ('P_INFLOW')  
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
                IF ((ENERGY_EQ .OR. RO_G0==UNDEFINED .OR. MU_G0==UNDEFINED)&
                .AND. BC_T_G(BCV)==UNDEFINED) THEN 
                    IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 'BC_T_g', BCV 
                    call mfix_exit(myPE)  
                ENDIF 

                IF (SPECIES_EQ(0)) THEN 
                    SUM = ZERO 
                    DO N = 1, NMAX(0) 
                        SUM = SUM + BC_X_G(BCV,N) 
                        IF (BC_X_G(BCV,N) == UNDEFINED) THEN 
                            IF(DMP_LOG)WRITE (UNIT_LOG, 1060) BCV, N 
                            call mfix_exit(myPE)  
                        ENDIF 
                    ENDDO 
                    IF (.NOT.COMPARE(ONE,SUM)) THEN 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1065) BCV 
                        call mfix_exit(myPE)  
                    ENDIF 
                ENDIF 

                DO M = 1, SMAX 
                    IF (BC_ROP_S(BCV,M) == UNDEFINED) THEN 
                        IF (BC_EP_G(BCV) == ONE) THEN 
                            BC_ROP_S(BCV,M) = ZERO 
                        ELSEIF (SMAX == 1) THEN 
                            BC_ROP_S(BCV,M) = (ONE - BC_EP_G(BCV))*RO_S(M) 
                        ELSE 
                            IF(DMP_LOG)WRITE (UNIT_LOG, 1100) 'BC_ROP_s', BCV, M 
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
                    IF (SPECIES_EQ(M)) THEN 
                        SUM = ZERO 
                        DO N = 1, NMAX(M) 
                            IF (BC_X_S(BCV,M,N) == UNDEFINED) THEN 
                                IF (BC_ROP_S(BCV,M) == ZERO) THEN 
                                    BC_X_S(BCV,M,N) = ZERO 
                                ELSE 
                                    IF(DMP_LOG)WRITE (UNIT_LOG, 1110) BCV, M, N 
                                    call mfix_exit(myPE)  
                                ENDIF 
                            ENDIF 
                            SUM = SUM + BC_X_S(BCV,M,N) 
                        ENDDO 
                        IF (.NOT.COMPARE(ONE,SUM)) THEN 
                            IF (SUM /= ZERO) THEN 
                                IF(DMP_LOG)WRITE (UNIT_LOG, 1120) BCV, M 
                                call mfix_exit(myPE)  
                            ENDIF 
                        ENDIF 
                    ENDIF 
                ENDDO 
 

                DO N = 1, NScalar
                    IF (BC_Scalar(BCV,N) == UNDEFINED) THEN 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1004) 'BC_Scalar', BCV, N 
                        CALL MFIX_EXIT(myPE)
                    ENDIF 
                ENDDO


            CASE ('P_OUTFLOW')  
                IF (BC_P_G(BCV) == UNDEFINED) THEN 
                    IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 'BC_P_g', BCV 
                    call mfix_exit(myPE)  
                ELSEIF (BC_P_G(BCV)<=ZERO .AND. RO_G0==UNDEFINED) THEN 
                    IF(DMP_LOG)WRITE (UNIT_LOG, 1010) BCV, BC_P_G(BCV) 
                    call mfix_exit(myPE)  
                ENDIF 


            CASE ('OUTFLOW')  
                IF (.NOT.ERROR) THEN 
                    ERROR = .TRUE. 
                ELSE 
                    IF(DMP_LOG)WRITE (UNIT_LOG, 1160) BCV 
                    call mfix_exit(myPE)  
                ENDIF 

            END SELECT 

!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
        ELSEIF ((BC_TYPE(BCV) /= 'DUMMY').AND.(BC_TYPE(BCV)(1:2) /= 'CG')) THEN 
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
!           Check whether BC values are specified for undefined BC locations
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

            DO M = 1, SMAX 
                IF (BC_ROP_S(BCV,M) /= UNDEFINED) THEN 
                    IF(DMP_LOG)WRITE (UNIT_LOG, 1300) 'BC_ROP_s', BCV, M 
                    call mfix_exit(myPE)  
                ENDIF 
                DO N = 1, DIMENSION_N_S 
                    IF (BC_X_S(BCV,M,N) /= UNDEFINED) THEN 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1300) 'BC_X_s', BCV, M 
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
    
            DO N = 1, NScalar
                IF (BC_Scalar(BCV,N) /= UNDEFINED) THEN 
                    IF(DMP_LOG)WRITE (UNIT_LOG, 1200) 'BC_Scalar', BCV
                    CALL MFIX_EXIT(myPE)
                ENDIF 
            ENDDO 
  
        ENDIF 
      ENDDO 


      DO BCV = 1, DIMENSION_BC 
        IF (BC_DEFINED(BCV)) THEN 

!           Default specification of Johnson-Jackson bc
            IF (GRANULAR_ENERGY) THEN 
                IF (BC_JJ_PS(BCV) == UNDEFINED_I) BC_JJ_PS(BCV) = 1 
            ELSE 
                IF (BC_JJ_PS(BCV) == UNDEFINED_I) BC_JJ_PS(BCV) = 0 
            ENDIF 

!           Set Jenkins default specification, modify BC_JJ accordingly
            IF (GRANULAR_ENERGY .AND. JENKINS) BC_JJ_PS(BCV) = 1

            IF (BC_TYPE(BCV)=='FREE_SLIP_WALL' .OR. BC_TYPE(BCV)=='NO_SLIP_WALL'&
            .OR. BC_TYPE(BCV)=='PAR_SLIP_WALL') THEN 
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
                ENDIF 

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
                ENDIF 

                IF (ENERGY_EQ) THEN 
                    IF (BC_HW_T_G(BCV) < ZERO) THEN 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1003) 'BC_hw_T_g', BCV 
                        call mfix_exit(myPE)  
                    ENDIF 

                    IF (BC_HW_T_G(BCV)/=ZERO .AND. BC_TW_G(BCV)==UNDEFINED) THEN 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 'BC_Tw_g', BCV 
                        call mfix_exit(myPE)  
                    ENDIF 

                    IF (BC_HW_T_G(BCV)/=UNDEFINED .AND. BC_C_T_G(BCV)==UNDEFINED) &
                    THEN 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 'BC_C_T_g', BCV 
                        call mfix_exit(myPE)  
                    ENDIF 

                    DO M = 1, SMAX 
                        IF (BC_HW_T_S(BCV,M) < ZERO) THEN 
                            IF(DMP_LOG)WRITE (UNIT_LOG, 1103) 'BC_hw_T_s', BCV, M 
                            call mfix_exit(myPE)  
                        ENDIF 

                        IF (BC_HW_T_S(BCV,M)/=ZERO .AND. BC_TW_S(BCV,M)==UNDEFINED) &
                        THEN 
                            IF(DMP_LOG)WRITE (UNIT_LOG, 1100) 'BC_Tw_s', BCV, M 
                            call mfix_exit(myPE)  
                        ENDIF 

                        IF (BC_HW_T_S(BCV,M)/=UNDEFINED .AND. BC_C_T_S(BCV,M)==&
                        UNDEFINED) THEN 
                            IF(DMP_LOG)WRITE (UNIT_LOG, 1100) 'BC_C_T_s', BCV, M 
                            call mfix_exit(myPE)  
                        ENDIF 
                    ENDDO 
                ENDIF 

                IF (GRANULAR_ENERGY .AND. BC_JJ_PS(BCV) .EQ. ZERO) THEN 
                    DO M = 1, MMAX 
                        IF (BC_HW_THETA_M(BCV,M) < ZERO) THEN 
                            IF(DMP_LOG)WRITE (UNIT_LOG, 1103) 'BC_hw_Theta_m', BCV, M 
                            call mfix_exit(myPE)  
                        ENDIF 

                        IF (BC_HW_THETA_M(BCV,M)/=ZERO .AND. BC_THETAW_M(BCV,M)==&
                        UNDEFINED) THEN 
                            IF(DMP_LOG)WRITE (UNIT_LOG, 1100) 'BC_Thetaw_m', BCV, M 
                            call mfix_exit(myPE)  
                        ENDIF 

                        IF (BC_HW_THETA_M(BCV,M)/=UNDEFINED .AND. BC_C_THETA_M(BCV,M&
                        )==UNDEFINED) THEN 
                            IF(DMP_LOG)WRITE (UNIT_LOG, 1100) 'BC_C_Theta_m', BCV, M 
                            call mfix_exit(myPE)  
                        ENDIF 
                    ENDDO 
                ENDIF 

                IF (SPECIES_EQ(0)) THEN 
                    DO N = 1, NMAX(0) 
                        IF (BC_HW_X_G(BCV,N) < ZERO) THEN 
                            IF(DMP_LOG)WRITE (UNIT_LOG, 1005) 'BC_hw_X_g', BCV, N 
                            call mfix_exit(myPE)  
                        ENDIF 

                        IF (BC_HW_X_G(BCV,N)/=ZERO .AND. BC_XW_G(BCV,N)==UNDEFINED) &
                        THEN 
                            IF(DMP_LOG)WRITE (UNIT_LOG, 1004) 'BC_Xw_g', BCV, N 
                            call mfix_exit(myPE)  
                        ENDIF 

                        IF (BC_HW_X_G(BCV,N)/=UNDEFINED .AND. BC_C_X_G(BCV,N)==&
                        UNDEFINED) THEN 
                            IF(DMP_LOG)WRITE (UNIT_LOG, 1004) 'BC_C_X_g', BCV, N 
                            call mfix_exit(myPE)  
                        ENDIF 
                    ENDDO 
                ENDIF 

                DO M = 1, SMAX 
                    IF (SPECIES_EQ(M)) THEN 
                        DO N = 1, NMAX(M) 
                            IF (BC_HW_X_S(BCV,M,N) < ZERO) THEN 
                                IF(DMP_LOG)WRITE (UNIT_LOG, 1105) 'BC_hw_X_s', BCV, M, N 
                                call mfix_exit(myPE)  
                            ENDIF 

                            IF (BC_HW_X_S(BCV,M,N)/=ZERO .AND. BC_XW_S(BCV,M,N)==&
                            UNDEFINED) THEN 
                                IF(DMP_LOG)WRITE (UNIT_LOG, 1104) 'BC_Xw_s', BCV, M, N 
                                call mfix_exit(myPE)  
                            ENDIF 

                            IF (BC_HW_X_S(BCV,M,N)/=UNDEFINED .AND. BC_C_X_S(BCV,M,N)&
                            ==UNDEFINED) THEN 
                                IF(DMP_LOG)WRITE (UNIT_LOG, 1104) 'BC_C_X_s', BCV, M, N 
                                call mfix_exit(myPE)  
                            ENDIF 
                        ENDDO 
                    ENDIF 
                ENDDO 
  
                DO N = 1, NScalar
                    IF (BC_HW_Scalar(BCV,N) < ZERO) THEN 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1005) 'BC_hw_Scalar', BCV, N 
                        CALL MFIX_EXIT(myPE)
                    ENDIF 

                    IF (BC_HW_Scalar(BCV,N)/=ZERO .AND. BC_Scalarw(BCV,N)==UNDEFINED) &
                    THEN 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1004) 'BC_ScalarW', BCV, N 
                        CALL MFIX_EXIT(myPE)
                    ENDIF 

                    IF (BC_HW_Scalar(BCV,N)/=UNDEFINED .AND. BC_C_Scalar(BCV,N)==&
                    UNDEFINED) THEN 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1004) 'BC_C_Scalar', BCV, N 
                        CALL MFIX_EXIT(myPE)
                    ENDIF 
                ENDDO 
            ENDIF 
        ENDIF 
      ENDDO 
!
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

!//SP Send Receive
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
 1050 FORMAT(/1X,70('*')//' From: CHECK_DATA_07',/' Message: BC number:',I2,&
         ' - ',A,' should be ',A,' zero.',/1X,70('*')/) 
 1060 FORMAT(/1X,70('*')//' From: CHECK_DATA_07',/' Message: BC_X_g(',I2,',',I2&
         ,') not specified',/1X,70('*')/) 
 1065 FORMAT(/1X,70('*')//' From: CHECK_DATA_07',/' Message: BC number:',I2,&
         ' - Sum of gas mass fractions is NOT equal to one',/1X,70('*')/) 
 1100 FORMAT(/1X,70('*')//' From: CHECK_DATA_07',/' Message: ',A,'(',I2,',',I1,&
         ') not specified',/1X,70('*')/) 
 1103 FORMAT(/1X,70('*')//' From: CHECK_DATA_07',/' Message: ',A,'(',I2,',',I1,&
         ') value is unphysical',/1X,70('*')/) 
 1104 FORMAT(/1X,70('*')//' From: CHECK_DATA_07',/' Message: ',A,'(',I2,',',I2,&
         ',',I2,') not specified',/1X,70('*')/) 
 1105 FORMAT(/1X,70('*')//' From: CHECK_DATA_07',/' Message: ',A,'(',I2,',',I2,&
         ',',I2,') value is unphysical',/1X,70('*')/) 
 1110 FORMAT(/1X,70('*')//' From: CHECK_DATA_07',/' Message: BC_X_s(',I2,',',I2&
         ,',',I2,') not specified',/1X,70('*')/) 
 1120 FORMAT(/1X,70('*')//' From: CHECK_DATA_07',/' Message: BC number:',I2,&
         ' - Sum of solids-',I1,' mass fractions is NOT equal to one',/1X,70(&
         '*')/) 
 1125 FORMAT(/1X,70('*')//' From: CHECK_DATA_07',/' Message: BC number:',I2,&
         ' - Sum of volume fractions is NOT equal to one',/1X,70('*')/) 
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
      END SUBROUTINE CHECK_DATA_07 

!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization
!// 220 Use local FUNIJK for triple DO i,j,k loop
!// 350 1206 change do loop limits: 1,kmax2->kmin3,kmax3      
!// 990 Replace STOP with exitMPI to terminate all processors
