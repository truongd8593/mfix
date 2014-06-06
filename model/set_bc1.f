!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SET_BC1                                                 C
!  Purpose: Set transient flow boundary conditions                     C
!                                                                      C
!  Author: M. Syamlal                                 Date: 29-JAN-92  C
!  Reviewer: P. Nicoletti, W. Rogers, S. Venkatesan   Date: 29-JAN-92  C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: Add calculations for mass outflow boundary condition       C
!  Author: M. Syamlal                                 Date: 23-OCT-92  C
!  Reviewer: M. Syamlal                               Date: 11-DEC-92  C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: BC_DEFINED, BC_I_w, BC_I_e, BC_J_s, BC_J_n,   C
!                        BC_K_b, BC_K_t, BC_TYPE, BC_TIME, BC_PLANE,   C
!                        TIME, TSTOP, DT, BC_TIME, BC_DT_0, SMAX,      C
!                        BC_JET_G, BC_JET_GH, BC_JET_GL, BC_DT_H,      C
!                        BC_DT_L, KT_TYPE                             C
!                                                                      C
!  Variables modified:                                                 C
!     For MI: BC_TIME, BC_JET_G, U_g, V_g, W_g                         C
!     FOR PI: nothing                                                  C
!     For PO, O: BC_TIME, BC_MOUT_G, BC_VOUT_G, BC_OUT_N               C
!                BC_MOUT_S, BC_VOUT_S                                  C
!     For MO: BC_TIME, BC_MOUT_G, BC_VOUT_G, BC_OUT_N                  C
!             BC_MOUT_S, BC_VOUT_S,                                    C
!             BC_U_g, BC_V_g, BC_W_g, BC_U_S, BC_V_S, BC_W_S,          C
!             U_g, V_g, W_g, U_s, V_s, W_s                             C
!     For PO, O & MO: set_outflow is called which modifies many        C
!                     quantities                                       C
!                                                                      C
!  Local variables: IJK, IJK2, I, J, K I1, I2, J1, J2, K1, K2, L, M    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE SET_BC1 

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param 
      USE param1 
      USE bc
      USE fldvar
      USE geometry
      USE indices
      USE physprop
      USE run
      USE funits 
      USE compar 
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices 
      INTEGER :: I, J, K, IJK 
! IJK index for setting velocity bc 
      INTEGER :: IJK2 
! Solids phase index
      INTEGER :: M 
! Local index for boundary condition 
      INTEGER :: L 
! Starting and ending I index 
      INTEGER :: I1, I2 
! Starting and ending J index 
      INTEGER :: J1, J2
! Starting and ending K index 
      INTEGER :: K1, K2
! local solids velocity for mixture (for ghd)
      DOUBLE PRECISION :: lvel_s
!-----------------------------------------------
! Include statement functions
!-----------------------------------------------
      INCLUDE 'function.inc'
!-----------------------------------------------

! Set the boundary conditions
      DO L = 1, DIMENSION_BC 
         IF (BC_DEFINED(L)) THEN 

! The range of boundary cells
            I1 = BC_I_W(L) 
            I2 = BC_I_E(L) 
            J1 = BC_J_S(L) 
            J2 = BC_J_N(L) 
            K1 = BC_K_B(L) 
            K2 = BC_K_T(L) 

! ---------------------------------------------------------------->>>
            IF (BC_TYPE(L) == 'MASS_OUTFLOW') THEN 

               CALL SET_OUTFLOW (L, I1, I2, J1, J2, K1, K2) 

               CALL CALC_OUTFLOW (L)

! Calculate and accumulate the actual mass and volume outflow               
               IF (TIME + 0.1d0*DT>=BC_TIME(L) .OR. TIME+0.1d0*DT>=TSTOP) THEN 
                  BC_TIME(L) = TIME + BC_DT_0(L) 

! Average and print out the flow rates
                  BC_MOUT_G(L) = ABS(BC_MOUT_G(L))/BC_OUT_N(L) 
                  BC_VOUT_G(L) = ABS(BC_VOUT_G(L))/BC_OUT_N(L) 
                  CALL START_LOG 
                  IF(DMP_LOG)WRITE (UNIT_LOG, 1000) L, TIME 
                  IF(DMP_LOG)WRITE (UNIT_LOG, 1100) BC_MOUT_G(L), BC_VOUT_G(L) 
                  DO M = 1, SMAX 
                     BC_MOUT_S(L,M) = ABS(BC_MOUT_S(L,M))/BC_OUT_N(L) 
                     BC_VOUT_S(L,M) = ABS(BC_VOUT_S(L,M))/BC_OUT_N(L) 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1200) M, BC_MOUT_S(L,M), BC_VOUT_S(L,M) 
                  END DO 
                  CALL END_LOG 
                  BC_OUT_N(L) = 0 

! Adjust the velocities if needed
                  IF (BC_MASSFLOW_G(L) /= UNDEFINED) THEN 
                     IF (BC_MOUT_G(L) > SMALL_NUMBER) THEN 
                        SELECT CASE (TRIM(BC_PLANE(L)))  
                        CASE ('W')  
                           BC_U_G(L) = BC_U_G(L)*BC_MASSFLOW_G(L)/BC_MOUT_G(L) 
                        CASE ('E')  
                           BC_U_G(L) = BC_U_G(L)*BC_MASSFLOW_G(L)/BC_MOUT_G(L) 
                        CASE ('S')  
                           BC_V_G(L) = BC_V_G(L)*BC_MASSFLOW_G(L)/BC_MOUT_G(L) 
                        CASE ('N')  
                           BC_V_G(L) = BC_V_G(L)*BC_MASSFLOW_G(L)/BC_MOUT_G(L) 
                        CASE ('B')  
                           BC_W_G(L) = BC_W_G(L)*BC_MASSFLOW_G(L)/BC_MOUT_G(L) 
                        CASE ('T')  
                           BC_W_G(L) = BC_W_G(L)*BC_MASSFLOW_G(L)/BC_MOUT_G(L) 
                        END SELECT 
                     ENDIF 
                  ELSEIF (BC_VOLFLOW_G(L) /= UNDEFINED) THEN 
                     IF (BC_VOUT_G(L) > SMALL_NUMBER) THEN 
                        SELECT CASE (TRIM(BC_PLANE(L)))  
                        CASE ('W')  
                           BC_U_G(L) = BC_U_G(L)*BC_VOLFLOW_G(L)/BC_VOUT_G(L) 
                        CASE ('E')  
                           BC_U_G(L) = BC_U_G(L)*BC_VOLFLOW_G(L)/BC_VOUT_G(L) 
                        CASE ('S')  
                           BC_V_G(L) = BC_V_G(L)*BC_VOLFLOW_G(L)/BC_VOUT_G(L) 
                        CASE ('N')  
                           BC_V_G(L) = BC_V_G(L)*BC_VOLFLOW_G(L)/BC_VOUT_G(L) 
                        CASE ('B')  
                           BC_W_G(L) = BC_W_G(L)*BC_VOLFLOW_G(L)/BC_VOUT_G(L) 
                        CASE ('T')  
                           BC_W_G(L) = BC_W_G(L)*BC_VOLFLOW_G(L)/BC_VOUT_G(L) 
                        END SELECT 
                     ENDIF 
                  ENDIF 

                  BC_MOUT_G(L) = ZERO 
                  BC_VOUT_G(L) = ZERO 
                  DO M = 1, SMAX 
                     IF (BC_MASSFLOW_S(L,M) /= UNDEFINED) THEN 
                        IF (BC_MOUT_S(L,M) > SMALL_NUMBER) THEN 
                           SELECT CASE (TRIM(BC_PLANE(L)))  
                           CASE ('W')  
                              BC_U_S(L,M) = BC_U_S(L,M)*BC_MASSFLOW_S(L,M)/&
                                 BC_MOUT_S(L,M) 
                           CASE ('E')  
                              BC_U_S(L,M) = BC_U_S(L,M)*BC_MASSFLOW_S(L,M)/&
                                 BC_MOUT_S(L,M) 
                           CASE ('S')  
                              BC_V_S(L,M) = BC_V_S(L,M)*BC_MASSFLOW_S(L,M)/&
                                 BC_MOUT_S(L,M) 
                           CASE ('N')  
                              BC_V_S(L,M) = BC_V_S(L,M)*BC_MASSFLOW_S(L,M)/&
                                 BC_MOUT_S(L,M) 
                           CASE ('B')  
                              BC_W_S(L,M) = BC_W_S(L,M)*BC_MASSFLOW_S(L,M)/&
                                 BC_MOUT_S(L,M) 
                           CASE ('T')  
                              BC_W_S(L,M) = BC_W_S(L,M)*BC_MASSFLOW_S(L,M)/&
                                 BC_MOUT_S(L,M) 
                           END SELECT 
                        ENDIF 
                     ELSEIF (BC_VOLFLOW_S(L,M) /= UNDEFINED) THEN 
                        IF (BC_VOUT_S(L,M) > SMALL_NUMBER) THEN 
                           SELECT CASE (TRIM(BC_PLANE(L)))  
                           CASE ('W')  
                              BC_U_S(L,M) = BC_U_S(L,M)*BC_VOLFLOW_S(L,M)/&
                                 BC_VOUT_S(L,M) 
                           CASE ('E')  
                              BC_U_S(L,M) = BC_U_S(L,M)*BC_VOLFLOW_S(L,M)/&
                                 BC_VOUT_S(L,M) 
                           CASE ('S')  
                              BC_V_S(L,M) = BC_V_S(L,M)*BC_VOLFLOW_S(L,M)/&
                                 BC_VOUT_S(L,M) 
                           CASE ('N')  
                              BC_V_S(L,M) = BC_V_S(L,M)*BC_VOLFLOW_S(L,M)/&
                                 BC_VOUT_S(L,M) 
                           CASE ('B')  
                              BC_W_S(L,M) = BC_W_S(L,M)*BC_VOLFLOW_S(L,M)/&
                                 BC_VOUT_S(L,M) 
                           CASE ('T')  
                              BC_W_S(L,M) = BC_W_S(L,M)*BC_VOLFLOW_S(L,M)/&
                                 BC_VOUT_S(L,M) 
                           END SELECT 
                        ENDIF 
                     ENDIF 
                     BC_MOUT_S(L,M) = ZERO 
                     BC_VOUT_S(L,M) = ZERO 
                  ENDDO

! Apply the boundary velocities - Defining the field variables at the 
! boundaries according to user specifications with any modifications
! from the above calculations.  If a W, S, or B plane (i.e., fluid cell
! is west, south or bottom of boundary cell) then define the velocity of
! the fluid cell according to the boundary velocity.
                  DO K = BC_K_B(L), BC_K_T(L) 
                     DO J = BC_J_S(L), BC_J_N(L) 
                        DO I = BC_I_W(L), BC_I_E(L) 
                         IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                        IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
                           IJK = FUNIJK(I,J,K) 
                           SELECT CASE (TRIM(BC_PLANE(L)))  
                           CASE ('W')  
                              IJK2 = IM_OF(IJK) 
                              U_G(IJK2) = BC_U_G(L) 
                           CASE ('E')  
                              U_G(IJK) = BC_U_G(L) 
                           CASE ('S')  
                              IJK2 = JM_OF(IJK) 
                              V_G(IJK2) = BC_V_G(L) 
                           CASE ('N')  
                              V_G(IJK) = BC_V_G(L) 
                           CASE ('B')  
                              IJK2 = KM_OF(IJK) 
                              W_G(IJK2) = BC_W_G(L) 
                           CASE ('T')  
                              W_G(IJK) = BC_W_G(L) 
                           END SELECT 
                           DO M = 1, SMAX 
                              SELECT CASE (TRIM(BC_PLANE(L)))  
                              CASE ('W')  
                                 IJK2 = IM_OF(IJK) 
                                 U_S(IJK2,M) = BC_U_S(L,M) 
                              CASE ('E')  
                                 U_S(IJK,M) = BC_U_S(L,M) 
                              CASE ('S')  
                                 IJK2 = JM_OF(IJK) 
                                 V_S(IJK2,M) = BC_V_S(L,M) 
                              CASE ('N')  
                                 V_S(IJK,M) = BC_V_S(L,M) 
                              CASE ('B')  
                                 IJK2 = KM_OF(IJK) 
                                 W_S(IJK2,M) = BC_W_S(L,M) 
                              CASE ('T')  
                                 W_S(IJK,M) = BC_W_S(L,M) 
                              END SELECT 
                           ENDDO  

! compute mixutre velocity BC for GHD theory
                           IF(KT_TYPE_ENUM == GHD_2007) THEN
                              lvel_s = zero

! bulk density is set by set_outflow and is set in bc according to
! neighboring fluid appropriately. so we don't need to use an ijk2
! index here (ijk value will = ijk2 value).
                              DO M = 1, SMAX 
                                 SELECT CASE (TRIM(BC_PLANE(L)))
                                 CASE ('W', 'E')  
                                    lvel_s = lvel_s + BC_U_S(L,M)*ROP_S(IJK,M)
                                 CASE ('S', 'N')  
                                    lvel_s = lvel_s + BC_V_S(L,M)*ROP_S(IJK,M)
                                 CASE ('B', 'T')  
                                    lvel_s = lvel_s + BC_W_S(L,M)*ROP_S(IJK,M)
                                 END SELECT 
                              ENDDO 

                              IF (ROP_S(IJK,MMAX) > 0) THEN
                                 lvel_s = lvel_s /ROP_S(IJK,MMAX)
                              ELSE
                                 lvel_s = ZERO
                              ENDIF
                              
                              SELECT CASE (TRIM(BC_PLANE(L)))
                              CASE ('W') 
                                 IJK2 = IM_OF(IJK)
                                 U_S(IJK2,MMAX) =  lvel_s 
                              CASE ('E')
                                 U_S(IJK,MMAX) = lvel_s 
                              CASE ('S')
                                 IJK2 = JM_OF(IJK)
                                 V_S(IJK2,MMAX) = lvel_s 
                              CASE ('N')
                                 V_S(IJK,MMAX) = lvel_s 
                              CASE ('B')
                                 IJK2 = KM_OF(IJK)
                                 W_S(IJK2,MMAX) = lvel_s 
                              CASE ('T')  
                                 W_S(IJK,MMAX) = lvel_s 
                              END SELECT 

                           ENDIF   ! end if (kt_type_enum==ghd_2007)

                        ENDDO 
                     ENDDO 
                  ENDDO 
               ENDIF 
! end setting 'mass_outflow'
! ---------------------------------------------------------------->>>


! ---------------------------------------------------------------->>>
            ELSEIF (BC_TYPE(L) == 'MASS_INFLOW') THEN
! update transient jet conditions

               IF (TIME + 0.1d0*DT>=BC_TIME(L) .AND. BC_JET_G(L)/=UNDEFINED) THEN 
                  IF (BC_JET_G(L) == BC_JET_GH(L)) THEN 
                     BC_JET_G(L) = BC_JET_GL(L) 
                     BC_TIME(L) = TIME + BC_DT_L(L) 
                  ELSE IF (BC_JET_G(L) == BC_JET_GL(L)) THEN 
                     BC_JET_G(L) = BC_JET_GH(L) 
                     BC_TIME(L) = TIME + BC_DT_H(L) 
                  ELSE 
                     BC_JET_G(L) = BC_JET_GH(L) 
                     BC_TIME(L) = TIME + BC_DT_H(L) 
                  ENDIF 
                  DO K = BC_K_B(L), BC_K_T(L) 
                     DO J = BC_J_S(L), BC_J_N(L) 
                        DO I = BC_I_W(L), BC_I_E(L)
                          IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                        IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
                          IJK = FUNIJK(I,J,K) 
                           SELECT CASE (TRIM(BC_PLANE(L)))  
                           CASE ('W')  
                              IJK2 = IM_OF(IJK) 
                              U_G(IJK2) = BC_JET_G(L) 
                           CASE ('E')  
                              U_G(IJK) = BC_JET_G(L) 
                           CASE ('S')  
                              IJK2 = JM_OF(IJK) 
                              V_G(IJK2) = BC_JET_G(L) 
                           CASE ('N')  
                              V_G(IJK) = BC_JET_G(L) 
                           CASE ('B')  
                              IJK2 = KM_OF(IJK) 
                              W_G(IJK2) = BC_JET_G(L) 
                           CASE ('T')  
                              W_G(IJK) = BC_JET_G(L) 
                           END SELECT 
                        END DO 
                     END DO 
                  END DO 
               ENDIF 
! end setting 'mass_inflow'
! ----------------------------------------------------------------<<<

! ---------------------------------------------------------------->>>
            ELSEIF (BC_TYPE(L) == 'P_INFLOW') THEN 
! No need to do anything
! ----------------------------------------------------------------<<<


! ---------------------------------------------------------------->>>
            ELSEIF (BC_TYPE(L)=='P_OUTFLOW' .OR. BC_TYPE(L)=='OUTFLOW') THEN 

               CALL SET_OUTFLOW (L, I1, I2, J1, J2, K1, K2) 

               IF (BC_DT_0(L) /= UNDEFINED) THEN 
! Calculate and accumulate the actual mass and volume outflow
                  CALL CALC_OUTFLOW (L) 
                  IF (TIME + 0.1d0*DT>=BC_TIME(L) .OR. TIME+0.1d0*DT>=TSTOP) THEN 
                     BC_TIME(L) = TIME + BC_DT_0(L) 

! Average and Print out the flow rates
                     BC_MOUT_G(L) = ABS(BC_MOUT_G(L))/BC_OUT_N(L) 
                     BC_VOUT_G(L) = ABS(BC_VOUT_G(L))/BC_OUT_N(L) 
                     CALL START_LOG 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1000) L, TIME 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1100) BC_MOUT_G(L), BC_VOUT_G(L) 
                     BC_MOUT_G(L) = ZERO 
                     BC_VOUT_G(L) = ZERO 
                     DO M = 1, SMAX 
                        BC_MOUT_S(L,M) = ABS(BC_MOUT_S(L,M))/BC_OUT_N(L) 
                        BC_VOUT_S(L,M) = ABS(BC_VOUT_S(L,M))/BC_OUT_N(L) 
                        IF(DMP_LOG)WRITE(UNIT_LOG,1200)M,BC_MOUT_S(L,M),BC_VOUT_S(L,M) 
                        BC_MOUT_S(L,M) = ZERO 
                        BC_VOUT_S(L,M) = ZERO 
                     END DO 
                     CALL END_LOG 
                     BC_OUT_N(L) = 0 
                  ENDIF 
               ENDIF 
! end setting 'p_outflow' or 'outflow'
! ----------------------------------------------------------------<<<

            ENDIF    ! end if/else branch if('mass_outflow'; 
                     ! mass_inflow'; p_inflow'; 'p_outflow' or 
                     ! 'outflow')
         ENDIF   ! end if (bc_defined(l))

      ENDDO    ! end do loop (l=1,dimension_bc)


      RETURN  
 1000 FORMAT(/,1X,'Average outflow rates at BC No. ',I2,'  At Time = ',G12.5) 
 1100 FORMAT(3X,'Gas : Mass flow = ',G12.5,'     Volumetric flow = ',G12.5) 
 1200 FORMAT(3X,'Solids-',I1,' : Mass flow = ',G12.5,&
         '     Volumetric flow = ',G12.5) 

      END SUBROUTINE SET_BC1 


