!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: FLOW_TO_VEL                                             C
!  Purpose: Convert volumetric and mass flow rates to velocities       C
!     A specified mass flow rate is first converted to volumetric      C
!     flow rate. The volumetric flow rate is then converted to a       C
!     velocity.                                                        C
!                                                                      C
!  Author: M. Syamlal                                 Date: 28-JUL-92  C
!  Reviewer: W. Rogers                                Date: 11-DEC-92  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE FLOW_TO_VEL(DO_VEL_CHECK)

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
      USE scales
      USE indices
      USE funits 
      USE compar 
      USE discretelement
      USE mfix_pic

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! loop/variable indices 
      INTEGER :: BCV, M 
! Volumetric flow rate computed from mass flow rate 
      DOUBLE PRECISION :: VOLFLOW 
! Velocity computed from volumetric flow rate 
      DOUBLE PRECISION :: VEL 
! Solids phase volume fraction 
      DOUBLE PRECISION :: EPS 
! Whether any volumetric flow conversion was done 
      LOGICAL :: CONVERTED,DO_VEL_CHECK 
! Average molecular weight 
      DOUBLE PRECISION :: MW 
!-----------------------------------------------
! External functions
!-----------------------------------------------
      DOUBLE PRECISION, EXTERNAL :: EOSG, CALC_MW 
      LOGICAL, EXTERNAL :: COMPARE 
!-----------------------------------------------

! When both flow rates and velocities are specified, a consistency check is done
! The first time flow_to_vel is called in 
! by setting the logical DO_VEL_CHECK to .TRUE.
! If cut-cells are not used, flow_to_vel is only called once.
! When cut-cells are used, flow_to_vel is called another time after
! the cut-cell preprocessing stage. During, the second call, the velocity check 
! should not be performed, because the velocity assigned suring the first call
! will not match the flow rate. Therfore, when called from cut_cell_preprocessing.f
! DO_VEL_CHECK is set to .FALSE.




! initialize
      VOLFLOW = UNDEFINED

      CONVERTED = .FALSE. 
      DO BCV = 1, DIMENSION_BC 
         IF (BC_DEFINED(BCV)) THEN 
            IF (BC_TYPE(BCV)=='MASS_INFLOW' .OR.&
                BC_TYPE(BCV)=='MASS_OUTFLOW') THEN 

! If gas mass flow is defined convert it to volumetric flow
! ---------------------------------------------------------------->>>
               IF (BC_MASSFLOW_G(BCV) /= UNDEFINED) THEN 
                  IF (RO_G0 /= UNDEFINED) THEN 
                     VOLFLOW = BC_MASSFLOW_G(BCV)/RO_G0 
                  ELSE 
                     IF (BC_P_G(BCV)/=UNDEFINED .AND. &
                         BC_T_G(BCV)/=UNDEFINED) THEN
                        IF (MW_AVG == UNDEFINED) THEN 
                           MW = CALC_MW(BC_X_G,DIMENSION_BC,BCV,NMAX(0),MW_G) 
                        ELSE 
                           MW = MW_AVG 
                        ENDIF 
                        VOLFLOW = BC_MASSFLOW_G(BCV)/&
                           EOSG(MW,(BC_P_G(BCV)-P_REF),BC_T_G(BCV))
                          
                     ELSE
! for mass_inflow, check_data_07 has already required that either ro_g0
! be defined or bc_p_g and bc_t_g be defined. So this branch will never
! be entered when mass_inflow. if mass_outflow, ro_g0 and either bc_p_g, 
! or bc_t_g, or both, must be undefined.
! if the code comes through this branch without exiting and massflow is
! non-zero, then volflow will remain undefined.

                        IF (BC_TYPE(BCV) == 'MASS_OUTFLOW') THEN ! this check seems unnecessary

! if no mass flow through the boundary, the volume flow is zero.
! otherwise check that the value of velocity component through the
! boundary plane is defined, and is non-zero (otherwise would be caught
! by bc_massflow_g == zero branch)
                           IF (BC_MASSFLOW_G(BCV) == ZERO) THEN 
                              VOLFLOW = ZERO 
                           ELSEIF (BC_PLANE(BCV)=='W' .OR. &
                                   BC_PLANE(BCV)=='E') THEN 
                              IF (BC_U_G(BCV)==UNDEFINED .OR. &
                                  BC_U_G(BCV)/=ZERO) THEN 
                                 IF(DMP_LOG)WRITE (UNIT_LOG, 1010)&
                                    BCV, 'BC_U_g' 
                                 call mfix_exit(myPE)  
                              ENDIF 
                           ELSEIF (BC_PLANE(BCV)=='N' .OR. &
                                   BC_PLANE(BCV)=='S') THEN 
                              IF (BC_V_G(BCV)==UNDEFINED .OR. &
                                  BC_V_G(BCV)/=ZERO) THEN 
                                 IF(DMP_LOG)WRITE (UNIT_LOG, 1010) &
                                    BCV, 'BC_V_g' 
                                 call mfix_exit(myPE)  
                              ENDIF 
                           ELSEIF (BC_PLANE(BCV)=='T' .OR. &
                                   BC_PLANE(BCV)=='B') THEN 
                              IF (BC_W_G(BCV)==UNDEFINED .OR. &
                                  BC_W_G(BCV)/=ZERO) THEN 
                                 IF(DMP_LOG)WRITE (UNIT_LOG, 1010) &
                                    BCV, 'BC_W_g' 
                                 call mfix_exit(myPE)  
                              ENDIF 
                           ENDIF
                        ELSE   ! else branch if(bc_type='mass_outflow')
! not sure how this branch will be reached by mass_inflow
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1020) BCV 
                           call mfix_exit(myPE)
                        ENDIF   ! end if (bc_type(bcv)=='mass_outflow')
                     ENDIF 
                  ENDIF   ! end if/else (ro_g0 /=undefined)

! If volumetric flow is also specified compare both
                  IF (BC_VOLFLOW_G(BCV) /= UNDEFINED) THEN 
! volflow may be undefined for mass_outflow boundaries wherein ro_g0 and
! either bc_p_g, or bc_t_g, or both, were undefined.                       
                     IF (.NOT.COMPARE(VOLFLOW,BC_VOLFLOW_G(BCV))) THEN 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1000) BCV, &
                           VOLFLOW, BC_VOLFLOW_G(BCV) 
                        call mfix_exit(myPE)  
                     ENDIF 
                  ELSE 
                     BC_VOLFLOW_G(BCV) = VOLFLOW 
                  ENDIF 
               ENDIF   ! end if (bc_massflow_g(bcv) /= undefined)
! end gas mass flow conversion to volumetric flow               
! ----------------------------------------------------------------<<<



! If gas volumetric flow is defined convert it to velocity
! ---------------------------------------------------------------->>>
               IF (BC_VOLFLOW_G(BCV) /= UNDEFINED) THEN 

                  IF (BC_EP_G(BCV) /= UNDEFINED) THEN
! volumetric flow rate and void fraction at the boundary are specified
! (known) so that the corresponding gas velocity through the boundary
! plane may be calculated. 
                     VEL = BC_VOLFLOW_G(BCV)/(BC_AREA(BCV)*BC_EP_G(BCV))
                  ELSE  
! for mass_inflow, check_data_07 has already required that bc_ep_g be
! defined (i.e., this section will not happen for MI). For mass_outflow
! the routine will exit here if bc_ep_g is not defined.  However, for
! this MO the boundary velocities must also be set or mfix will exit due
! to later checks in check_data_07.
                     RETURN          
                  ENDIF 

! if the user also defined the boundary velocity through the plane, then
! check that the calculated value agrees with the specified value. if
! the user did not define the boundary velocity through the plane, then
! if mass_inflow set the value of the boundary velocity to the
! calculated value. otherwise do nothing.
                  CONVERTED = .TRUE. 
                  SELECT CASE (TRIM(BC_PLANE(BCV)))  
                  CASE ('W')  
                     IF (BC_U_G(BCV) /= UNDEFINED.AND.DO_VEL_CHECK) THEN 
                        IF (BC_TYPE(BCV)=='MASS_INFLOW' .AND. &
                            .NOT.COMPARE((-VEL),BC_U_G(BCV))) THEN 
                           IF(DMP_LOG) WRITE (UNIT_LOG, 1100) BCV,&
                              (-VEL), 'BC_U_g', BC_U_G(BCV) 
                           call mfix_exit(myPE)  
                        ENDIF 
                        IF (BC_TYPE(BCV)=='MASS_OUTFLOW' .AND. &
                            .NOT.COMPARE(VEL,BC_U_G(BCV))) THEN 
                           IF(DMP_LOG) WRITE (UNIT_LOG, 1100) BCV, &
                              VEL, 'BC_U_g', BC_U_G(BCV) 
                           call mfix_exit(myPE)  
                        ENDIF 
                     ELSE 
                        IF (BC_TYPE(BCV) == 'MASS_INFLOW') THEN 
                           BC_U_G(BCV) = -VEL 
                           BC_V_G(BCV) = ZERO 
                           BC_W_G(BCV) = ZERO 
                        ELSE 
                           BC_U_G(BCV) = VEL 
                        ENDIF 
                     ENDIF 
                  CASE ('E')  
                     IF (BC_U_G(BCV) /= UNDEFINED.AND.DO_VEL_CHECK) THEN 
                        IF (BC_TYPE(BCV)=='MASS_INFLOW' .AND. &
                            .NOT.COMPARE(VEL,BC_U_G(BCV))) THEN 
                           IF(DMP_LOG) WRITE (UNIT_LOG, 1100) BCV, VEL,&
                              'BC_U_g', BC_U_G(BCV) 
                           call mfix_exit(myPE)  
                        ENDIF 
                        IF (BC_TYPE(BCV)=='MASS_OUTFLOW' .AND. &
                            .NOT.COMPARE((-VEL),BC_U_G(BCV))) THEN 
                           IF(DMP_LOG) WRITE (UNIT_LOG, 1100) BCV, &
                              (-VEL), 'BC_U_g', BC_U_G(BCV) 
                           call mfix_exit(myPE)  
                        ENDIF 
                     ELSE 
                        IF (BC_TYPE(BCV) == 'MASS_INFLOW') THEN 
                           BC_U_G(BCV) = VEL 
                           BC_V_G(BCV) = ZERO
                           BC_W_G(BCV) = ZERO
                        ELSE 
                           BC_U_G(BCV) = -VEL 
                        ENDIF 
                     ENDIF 
                  CASE ('S')  
                     IF (BC_V_G(BCV) /= UNDEFINED.AND.DO_VEL_CHECK) THEN 
                        IF (BC_TYPE(BCV)=='MASS_INFLOW' .AND. &
                            .NOT.COMPARE((-VEL),BC_V_G(BCV))) THEN 
                           IF(DMP_LOG) WRITE (UNIT_LOG, 1100) BCV,&
                              (-VEL), 'BC_V_g', BC_V_G(BCV) 
                           call mfix_exit(myPE)  
                        ENDIF 
                        IF (BC_TYPE(BCV)=='MASS_OUTFLOW' .AND. &
                            .NOT.COMPARE(VEL,BC_V_G(BCV))) THEN 
                           IF(DMP_LOG) WRITE (UNIT_LOG, 1100) BCV, VEL,&
                              'BC_V_g', BC_V_G(BCV) 
                           call mfix_exit(myPE)  
                        ENDIF 
                     ELSE 
                        IF (BC_TYPE(BCV) == 'MASS_INFLOW') THEN 
                           BC_U_G(BCV) = ZERO 
                           BC_V_G(BCV) = -VEL 
                           BC_W_G(BCV) = ZERO 
                        ELSE 
                           BC_V_G(BCV) = VEL 
                        ENDIF 
                     ENDIF 
                  CASE ('N')  
                     IF (BC_V_G(BCV) /= UNDEFINED.AND.DO_VEL_CHECK) THEN 
                        IF (BC_TYPE(BCV)=='MASS_INFLOW' .AND. &
                            .NOT.COMPARE(VEL,BC_V_G(BCV))) THEN 
                           IF(DMP_LOG) WRITE (UNIT_LOG, 1100) BCV, VEL,&
                              'BC_V_g', BC_V_G(BCV) 
                           call mfix_exit(myPE)  
                        ENDIF 
                        IF (BC_TYPE(BCV)=='MASS_OUTFLOW' .AND. &
                            .NOT.COMPARE((-VEL),BC_V_G(BCV))) THEN 
                           IF(DMP_LOG) WRITE (UNIT_LOG, 1100) BCV, &
                              (-VEL), 'BC_V_g', BC_V_G(BCV) 
                           call mfix_exit(myPE)  
                        ENDIF 
                     ELSE 
                        IF (BC_TYPE(BCV) == 'MASS_INFLOW') THEN 
                           BC_U_G(BCV) = ZERO
                           BC_V_G(BCV) = VEL 
                           BC_W_G(BCV) = ZERO
                        ELSE 
                           BC_V_G(BCV) = -VEL 
                        ENDIF 
                     ENDIF 
                  CASE ('B')  
                     IF (BC_W_G(BCV) /= UNDEFINED.AND.DO_VEL_CHECK) THEN 
                        IF (BC_TYPE(BCV)=='MASS_INFLOW' .AND.  &
                            .NOT.COMPARE((-VEL),BC_W_G(BCV))) THEN 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1100) BCV,&
                              (-VEL), 'BC_W_g', BC_W_G(BCV) 
                           call mfix_exit(myPE)  
                        ENDIF 
                        IF (BC_TYPE(BCV)=='MASS_OUTFLOW' .AND. &
                            .NOT.COMPARE(VEL,BC_W_G(BCV))) THEN 
                           IF(DMP_LOG) WRITE (UNIT_LOG, 1100) BCV, VEL,&
                              'BC_W_g', BC_W_G(BCV) 
                           call mfix_exit(myPE)  
                        ENDIF 
                     ELSE 
                        IF (BC_TYPE(BCV) == 'MASS_INFLOW') THEN 
                           BC_U_G(BCV) = ZERO 
                           BC_V_G(BCV) = ZERO 
                           BC_W_G(BCV) = -VEL 
                        ELSE 
                           BC_W_G(BCV) = VEL 
                        ENDIF 
                     ENDIF 
                  CASE ('T')  
                     IF (BC_W_G(BCV) /= UNDEFINED.AND.DO_VEL_CHECK) THEN 
                        IF (BC_TYPE(BCV)=='MASS_INFLOW' .AND. &
                            .NOT.COMPARE(VEL,BC_W_G(BCV))) THEN 
                           IF(DMP_LOG) WRITE (UNIT_LOG, 1100) BCV, VEL,&
                              'BC_W_g', BC_W_G(BCV) 
                           call mfix_exit(myPE)  
                        ENDIF 
                        IF (BC_TYPE(BCV)=='MASS_OUTFLOW' .AND. &
                            .NOT.COMPARE((-VEL),BC_W_G(BCV))) THEN 
                           IF(DMP_LOG) WRITE (UNIT_LOG, 1100) BCV,&
                              (-VEL), 'BC_W_g', BC_W_G(BCV) 
                           call mfix_exit(myPE)  
                        ENDIF 
                     ELSE 
                        IF (BC_TYPE(BCV) == 'MASS_INFLOW') THEN 
                           BC_U_G(BCV) = ZERO
                           BC_V_G(BCV) = ZERO
                           BC_W_G(BCV) = VEL 
                        ELSE 
                           BC_W_G(BCV) = -VEL 
                        ENDIF 
                     ENDIF 
                  END SELECT    ! end select (trim(bc_plane(bcv))
               ENDIF   ! end if (bc_volflow_g(bcv)/=undefined)
! end gas volumetric flow conversion to velocity
! ----------------------------------------------------------------<<<               
               
               IF (.NOT.DISCRETE_ELEMENT .OR. (DISCRETE_ELEMENT &
                   .AND. DES_CONTINUUM_HYBRID).OR. (DISCRETE_ELEMENT &
                   .AND. MPPIC)) THEN
! The following quantities should not be required for DEM simulations
! To ensure this is the case leave them undefined in mfix.dat
! MPPIC BC's are based on Two Fluid Model based specification. 
! So call the below for MPPIC. 

! Do flow conversions for solids phases
               DO M = 1, SMAX 

! initialize
                  VOLFLOW = UNDEFINED

! If solids mass flow is defined convert it to volumetric flow
! ---------------------------------------------------------------->>>
                  IF (BC_MASSFLOW_S(BCV,M) /= UNDEFINED) THEN 

                     IF (RO_S(M) /= UNDEFINED) THEN 
! RO_S must be defined for solids phases (see check_data_04).
                        VOLFLOW = BC_MASSFLOW_S(BCV,M)/RO_S(M) 
                     ELSE 
! this section should never happen (redundant).
                        RETURN    
                     ENDIF 

! If volumetric flow is also specified compare both
                     IF (BC_VOLFLOW_S(BCV,M) /= UNDEFINED) THEN 
                        IF (.NOT.COMPARE(VOLFLOW,BC_VOLFLOW_S(BCV,M))) THEN 
                           IF(DMP_LOG) WRITE(UNIT_LOG,1200) &
                              BCV,VOLFLOW,M,BC_VOLFLOW_S(BCV,M) 
                           call mfix_exit(myPE)  
                        ENDIF 
                     ELSE 
                        BC_VOLFLOW_S(BCV,M) = VOLFLOW 
                     ENDIF 
                  ENDIF   ! end if (bc_massflow_s(bcv,m)/=undefined)
! end solids mass flow conversion to volumetric flow
! ----------------------------------------------------------------<<<

                  
! if possible, define bulk density based on ep_g. note if mass_outflow,
! bc_ep_g may still be undefined at this point wherein bc_rop_s would 
! become set as 1-undefined. to avoid this issue, a check for bc_ep_g 
! defined was added. note if mass_inflow, check_data_07 later performs
! the same calculations as below but with additional checks.
                  IF (BC_ROP_S(BCV,M) == UNDEFINED .AND. &
                     BC_EP_G(BCV) /= UNDEFINED ) THEN 
                     IF (BC_EP_G(BCV) == ONE) THEN 
                         BC_ROP_S(BCV,M) = ZERO 
                     ELSEIF (SMAX == 1 .AND. &
                              .NOT.DES_CONTINUUM_HYBRID) THEN 
! bulk density must be explicitly defined for hybrid model and cannot be
! defined from 1-bc_ep_g
                         BC_ROP_S(BCV,M) = (ONE - BC_EP_G(BCV))*RO_S(M)
                     ENDIF
                  ENDIF
! note bc_rop_s may still be undefined at this point


! If solids volumetric flow is defined convert it to velocity
! ---------------------------------------------------------------->>>
                  IF (BC_VOLFLOW_S(BCV,M) /= UNDEFINED) THEN 

                     IF (BC_ROP_S(BCV,M) /= UNDEFINED) THEN 
                        EPS = BC_ROP_S(BCV,M)/RO_S(M)
! volumetric flow rate and solids volume fraction at the boundary are
! specified (known) so that the corresponding solids velocity through
! the boundary plane may be calculated. 
                        IF (EPS /= ZERO) THEN 
                           VEL = BC_VOLFLOW_S(BCV,M)/(BC_AREA(BCV)*EPS) 
                        ELSE 
                           IF (BC_VOLFLOW_S(BCV,M) == ZERO) THEN 
                              VEL = ZERO 
                           ELSE 
                              IF(DMP_LOG)WRITE (UNIT_LOG, 1250) BCV, M 
                              call mfix_exit(myPE)  
                           ENDIF 
                        ENDIF
                     ELSE   ! bc_rop_s is undefined
                        IF (BC_VOLFLOW_S(BCV,M) == ZERO) THEN 
                           VEL = ZERO 
                        ELSE 
! if bc_rop_s is undefined and bc_volflow_s is not zero exit MFIX with
! error                                
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1260) BCV, M 
                           call mfix_exit(myPE)  
                        ENDIF 
                     ENDIF   

! if the user also defined the boundary velocity through the plane, then
! check that the calculated value agrees with the specified value. if
! the user did not define the boundary velocity through the plane, then
! if mass_inflow set the value of the boundary velocity to the
! calculated value. otherwise do nothing.
                     CONVERTED = .TRUE. 
                     SELECT CASE (TRIM(BC_PLANE(BCV)))  
                     CASE ('W')  
                        IF (BC_U_S(BCV,M) /= UNDEFINED.AND.DO_VEL_CHECK) THEN 
                           IF (BC_TYPE(BCV)=='MASS_INFLOW' .AND. &
                               .NOT.COMPARE((-VEL),BC_U_S(BCV,M))) THEN 
                              IF(DMP_LOG)WRITE (UNIT_LOG, 1300) BCV, &
                                 (-VEL), 'BC_U_s', M, BC_U_S(BCV,M) 
                              call mfix_exit(myPE)  
                           ENDIF 
                           IF (BC_TYPE(BCV)=='MASS_OUTFLOW' .AND. &
                               .NOT.COMPARE(VEL,BC_U_S(BCV,M))) THEN 
                              IF(DMP_LOG)WRITE (UNIT_LOG, 1300) BCV, &
                                 VEL, 'BC_U_s', M, BC_U_S(BCV,M) 
                              call mfix_exit(myPE)  
                           ENDIF 
                        ELSE 
                           IF (BC_TYPE(BCV) == 'MASS_INFLOW') THEN 
                              BC_U_S(BCV,M) = -VEL 
                              BC_V_S(BCV,M) = ZERO 
                              BC_W_S(BCV,M) = ZERO 
                           ELSE 
                              BC_U_S(BCV,M) = VEL 
                           ENDIF 
                        ENDIF 
                     CASE ('E')  
                        IF (BC_U_S(BCV,M) /= UNDEFINED.AND.DO_VEL_CHECK) THEN 
                           IF (BC_TYPE(BCV)=='MASS_INFLOW' .AND. &
                               .NOT.COMPARE(VEL,BC_U_S(BCV,M))) THEN 
                              IF(DMP_LOG)WRITE (UNIT_LOG, 1300) BCV, &
                                 VEL, 'BC_U_s', M, BC_U_S(BCV,M) 
                              call mfix_exit(myPE)  
                           ENDIF 
                           IF (BC_TYPE(BCV)=='MASS_OUTFLOW' .AND. &
                               .NOT.COMPARE((-VEL),BC_U_S(BCV,M))) THEN 
                              IF(DMP_LOG)WRITE (UNIT_LOG, 1300) BCV, &
                                 (-VEL), 'BC_U_s', M, BC_U_S(BCV,M) 
                              call mfix_exit(myPE)  
                           ENDIF 
                        ELSE 
                           IF (BC_TYPE(BCV) == 'MASS_INFLOW') THEN 
                              BC_U_S(BCV,M) = VEL 
                              BC_V_S(BCV,M) = ZERO 
                              BC_W_S(BCV,M) = ZERO 
                           ELSE 
                              BC_U_S(BCV,M) = -VEL 
                           ENDIF 
                        ENDIF 
                     CASE ('S')  
                        IF (BC_V_S(BCV,M) /= UNDEFINED.AND.DO_VEL_CHECK) THEN 
                           IF (BC_TYPE(BCV)=='MASS_INFLOW' .AND. &
                               .NOT.COMPARE((-VEL),BC_V_S(BCV,M))) THEN 
                              IF(DMP_LOG)WRITE (UNIT_LOG, 1300) BCV, &
                                 (-VEL), 'BC_V_s', M, BC_V_S(BCV,M) 
                              call mfix_exit(myPE)  
                           ENDIF 
                           IF (BC_TYPE(BCV)=='MASS_OUTFLOW' .AND. &
                               .NOT.COMPARE(VEL,BC_V_S(BCV,M))) THEN 
                              IF(DMP_LOG)WRITE (UNIT_LOG, 1300) BCV, &
                                 VEL, 'BC_V_s', M, BC_V_S(BCV,M) 
                              call mfix_exit(myPE)  
                           ENDIF 
                        ELSE 
                           IF (BC_TYPE(BCV) == 'MASS_INFLOW') THEN 
                              BC_U_S(BCV,M) = ZERO 
                              BC_V_S(BCV,M) = -VEL 
                              BC_W_S(BCV,M) = ZERO 
                           ELSE 
                              BC_V_S(BCV,M) = VEL 
                           ENDIF 
                        ENDIF 
                     CASE ('N')  
                        IF (BC_V_S(BCV,M) /= UNDEFINED.AND.DO_VEL_CHECK) THEN 
                           IF (BC_TYPE(BCV)=='MASS_INFLOW' .AND. &
                               .NOT.COMPARE(VEL,BC_V_S(BCV,M))) THEN 
                              IF(DMP_LOG)WRITE (UNIT_LOG, 1300) BCV, &
                                 VEL, 'BC_V_s', M, BC_V_S(BCV,M) 
                              call mfix_exit(myPE)  
                           ENDIF 
                           IF (BC_TYPE(BCV)=='MASS_OUTFLOW' .AND. &
                               .NOT.COMPARE((-VEL),BC_V_S(BCV,M))) THEN 
                              IF(DMP_LOG)WRITE (UNIT_LOG, 1300) BCV, &
                                 (-VEL), 'BC_V_s', M, BC_V_S(BCV,M) 
                              call mfix_exit(myPE)  
                           ENDIF 
                        ELSE 
                           IF (BC_TYPE(BCV) == 'MASS_INFLOW') THEN 
                              BC_U_S(BCV,M) = ZERO
                              BC_V_S(BCV,M) = VEL 
                              BC_W_S(BCV,M) = ZERO
                           ELSE 
                              BC_V_S(BCV,M) = -VEL 
                           ENDIF 
                        ENDIF 
                     CASE ('B')  
                        IF (BC_W_S(BCV,M) /= UNDEFINED.AND.DO_VEL_CHECK) THEN 
                           IF (BC_TYPE(BCV)=='MASS_INFLOW' .AND. &
                               .NOT.COMPARE((-VEL),BC_W_S(BCV,M))) THEN 
                              IF(DMP_LOG)WRITE (UNIT_LOG, 1300) BCV, &
                                 (-VEL), 'BC_W_s', M, BC_W_S(BCV,M) 
                              call mfix_exit(myPE)  
                           ENDIF 
                           IF (BC_TYPE(BCV)=='MASS_OUTFLOW' .AND. &
                               .NOT.COMPARE(VEL,BC_W_S(BCV,M))) THEN 
                              IF(DMP_LOG)WRITE (UNIT_LOG, 1300) BCV, &
                                 VEL, 'BC_W_s', M, BC_W_S(BCV,M) 
                              call mfix_exit(myPE)  
                           ENDIF 
                        ELSE 
                           IF (BC_TYPE(BCV) == 'MASS_INFLOW') THEN 
                              BC_U_S(BCV,M) = ZERO 
                              BC_V_S(BCV,M) = ZERO 
                              BC_W_S(BCV,M) = -VEL 
                           ELSE 
                              BC_W_S(BCV,M) = VEL 
                           ENDIF 
                        ENDIF 
                     CASE ('T')  
                        IF (BC_W_S(BCV,M) /= UNDEFINED.AND.DO_VEL_CHECK) THEN 
                           IF (BC_TYPE(BCV)=='MASS_INFLOW' .AND. &
                               .NOT.COMPARE(VEL,BC_W_S(BCV,M))) THEN 
                              IF(DMP_LOG)WRITE (UNIT_LOG, 1300) BCV, &
                                 VEL, 'BC_W_s', M, BC_W_S(BCV,M) 
                              call mfix_exit(myPE)  
                           ENDIF 
                           IF (BC_TYPE(BCV)=='MASS_OUTFLOW' .AND. &
                              .NOT.COMPARE((-VEL),BC_W_S(BCV,M))) THEN 
                              IF(DMP_LOG)WRITE (UNIT_LOG, 1300) BCV, &
                                 (-VEL), 'BC_W_s', M, BC_W_S(BCV,M) 
                              call mfix_exit(myPE)  
                           ENDIF 
                        ELSE 
                           IF (BC_TYPE(BCV) == 'MASS_INFLOW') THEN 
                              BC_U_S(BCV,M) = ZERO 
                              BC_V_S(BCV,M) = ZERO 
                              BC_W_S(BCV,M) = VEL 
                           ELSE 
                              BC_W_S(BCV,M) = -VEL 
                           ENDIF 
                        ENDIF 
                     END SELECT    ! end select (trim(bc_plane(bcv))
                  ENDIF   ! end if (bc_volflow_s(bcv,m) /= undefined)
! end solids volumetric flow conversion to velocity
! ----------------------------------------------------------------<<<

               ENDDO   ! end do m = 1,smax
               ENDIF   ! end if (.not.discrete_element)

            ENDIF   ! end if (bc_type(bcv)=='mass_inflow' or 'mass_outflow')
         ENDIF   ! end if (bc_defined(bcv)) 
      ENDDO   ! end do (bcv =1,dimension_bc)


      IF (CONVERTED .AND. (NO_I .OR. NO_J .OR. NO_K) &
          .AND. DMP_LOG)WRITE (UNIT_LOG, 1500) 

      RETURN  

 1000 FORMAT(/1X,70('*')//' From: FLOW_TO_VEL',/' Message: BC No:',I2,/,&
         ' Computed volumetric flow is not equal to specified value',/,&
         ' Value computed from mass flow  = ',G14.7,/,&
         ' Specified value (BC_VOLFLOW_g) = ',G14.7,/1X,70('*')/) 
 1010 FORMAT(/1X,70('*')//' From: FLOW_TO_VEL',/' Message: BC No:',I2,&
         '  BC_P_g, BC_T_g, and BC_X_g or',/' a nonzero value for ',A,&
         ' should be specified',/1X,70('*')/) 
 1020 FORMAT(/1X,70('*')//' From: FLOW_TO_VEL',/' Message: BC No:',I2,&
         '  BC_P_g, BC_T_g, and BC_X_g',/' should be specified',/1X,70('*')/)
 1100 FORMAT(/1X,70('*')//' From: FLOW_TO_VEL',/' Message: BC No:',I2,/,&
         ' Computed velocity is not equal to specified value',/,&
         ' Value computed from vol. or mass flow  = ',G14.7,/,&
         ' Specified value (',A,') = ',G14.7,/1X,70('*')/) 
 1200 FORMAT(/1X,70('*')//' From: FLOW_TO_VEL',/' Message: BC No:',I2,/,&
         ' Computed volumetric flow is not equal to specified value',/,&
         ' Value computed from mass flow  = ',G14.7,/,&
         ' Specified value (BC_VOLFLOW_s',I1,') = ',G14.7,/1X,70('*')/) 
 1250 FORMAT(/1X,70('*')//' From: FLOW_TO_VEL',/' Message: BC No:',I2,/,&
         ' Non-zero vol. or mass flow specified with BC_ROP_s',&
         I1,' = 0.',/1X,70('*')/) 
 1260 FORMAT(/1X,70('*')//' From: FLOW_TO_VEL',/' Message: BC No:',I2,/,&
         ' BC_ROP_s',I1,' not specified',/1X,70('*')/) 
 1300 FORMAT(/1X,70('*')//' From: FLOW_TO_VEL',/' Message: BC No:',I2,/,&
         ' Computed velocity is not equal to specified value',/,&
         ' Value computed from vol. or mass flow  = ',G14.7,/,&
         ' Specified value (',A,I1,') = ',G14.7,/1X,70('*')/) 
 1500 FORMAT(/1X,70('*')//' From: FLOW_TO_VEL',/&
         ' Message: Some volumetric or mass flow rates have been',/&
         '   converted to velocity values.  In 2D simulations ensure',/&
         '   that the third (unused) dimension is correctly specified;',/&
         '   e.g. in axisymmetric cylindrical coordinates ZLENGTH = 2*Pi'/1X,70&
         ('*')/) 

      END SUBROUTINE FLOW_TO_VEL 

