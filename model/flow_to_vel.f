!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: FLOW_TO_VEL                                            C
!  Purpose: Convert volumetric and mass flow rates to velocities       C
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
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE FLOW_TO_VEL 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
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
      USE compar   !//
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
! 
!             loop/variable indices 
      INTEGER BCV, M 
! 
!                      Volumetric flow rate computed from mass flow rate 
      DOUBLE PRECISION VOLFLOW 
! 
!                      Velocity computed from volumetric flow rate 
      DOUBLE PRECISION VEL 
! 
!                      Solids phase volume fraction 
      DOUBLE PRECISION EPS 
! 
!                      Whether any volumetric flow conversion was done 
      LOGICAL          CONVERTED 
! 
!                      Average molecular weight 
      DOUBLE PRECISION MW 
! 
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      DOUBLE PRECISION , EXTERNAL :: EOSG, CALC_MW 
      LOGICAL , EXTERNAL :: COMPARE 
!-----------------------------------------------
!
!
      CONVERTED = .FALSE. 
      DO BCV = 1, DIMENSION_BC 
         IF (BC_DEFINED(BCV)) THEN 
            IF (BC_TYPE(BCV)=='MASS_INFLOW' .OR. BC_TYPE(BCV)=='MASS_OUTFLOW') &
               THEN 
!
!           If gas mass flow is defined convert it to volumetric flow
!
               IF (BC_MASSFLOW_G(BCV) /= UNDEFINED) THEN 
                  IF (RO_G0 /= UNDEFINED) THEN 
                     VOLFLOW = BC_MASSFLOW_G(BCV)/RO_G0 
                  ELSE 
                     IF (BC_P_G(BCV)/=UNDEFINED .AND. BC_T_G(BCV)/=UNDEFINED) &
                        THEN 
                        IF (MW_AVG == UNDEFINED) THEN 
                           MW = CALC_MW(BC_X_G,DIMENSION_BC,BCV,NMAX(0),MW_G) 
                        ELSE 
                           MW = MW_AVG 
                        ENDIF 
                        VOLFLOW = BC_MASSFLOW_G(BCV)/EOSG(MW,BC_P_G(BCV),BC_T_G(&
                           BCV)) 
                     ELSE 
                        IF (BC_TYPE(BCV) == 'MASS_OUTFLOW') THEN 
                           IF (BC_MASSFLOW_G(BCV) == ZERO) THEN 
                              VOLFLOW = ZERO 
                           ELSE IF (BC_PLANE(BCV)=='W' .OR. BC_PLANE(BCV)=='E') &
                                 THEN 
                              IF (BC_U_G(BCV)==UNDEFINED .OR. BC_U_G(BCV)==ZERO) &
                                 THEN 
                                 WRITE (UNIT_LOG, 1010) BCV, 'BC_U_g' 
                                 call mfix_exit(myPE)  
                              ENDIF 
                           ELSE IF (BC_PLANE(BCV)=='N' .OR. BC_PLANE(BCV)=='S') &
                                 THEN 
                              IF (BC_V_G(BCV)==UNDEFINED .OR. BC_V_G(BCV)==ZERO) &
                                 THEN 
                                 WRITE (UNIT_LOG, 1010) BCV, 'BC_V_g' 
                                 call mfix_exit(myPE)  
                              ENDIF 
                           ELSE IF (BC_PLANE(BCV)=='T' .OR. BC_PLANE(BCV)=='B') &
                                 THEN 
                              IF (BC_W_G(BCV)==UNDEFINED .OR. BC_W_G(BCV)==ZERO) &
                                 THEN 
                                 WRITE (UNIT_LOG, 1010) BCV, 'BC_W_g' 
                                 call mfix_exit(myPE)  
                              ENDIF 
                           ENDIF 
                        ELSE 
                           WRITE (UNIT_LOG, 1020) BCV 
                           call mfix_exit(myPE)  
                        ENDIF 
                     ENDIF 
                  ENDIF 
!
!             If volumetric flow is also specified compare both
!
                  IF (BC_VOLFLOW_G(BCV) /= UNDEFINED) THEN 
                     IF (.NOT.COMPARE(VOLFLOW,BC_VOLFLOW_G(BCV))) THEN 
                        WRITE (UNIT_LOG, 1000) BCV, VOLFLOW, BC_VOLFLOW_G(BCV) 
                        call mfix_exit(myPE)  
                     ENDIF 
                  ELSE 
                     BC_VOLFLOW_G(BCV) = VOLFLOW 
                  ENDIF 
               ENDIF 
!
!           If gas volumetric flow is defined convert it to velocity
!
               IF (BC_VOLFLOW_G(BCV) /= UNDEFINED) THEN 
                  IF (BC_EP_G(BCV) /= UNDEFINED) THEN 
                     VEL = BC_VOLFLOW_G(BCV)/(BC_AREA(BCV)*BC_EP_G(BCV)) 
                  ELSE 
                     RETURN                      !Error caught in Check_data_07 
                  ENDIF 
                  CONVERTED = .TRUE. 
                  SELECT CASE (BC_PLANE(BCV))  
                  CASE ('W')  
                     IF (BC_U_G(BCV) /= UNDEFINED) THEN 
                        IF (BC_TYPE(BCV)=='MASS_INFLOW' .AND.  .NOT.COMPARE((-&
                           VEL),BC_U_G(BCV))) THEN 
                           WRITE (UNIT_LOG, 1100) BCV, (-VEL), 'BC_U_g', BC_U_G(&
                              BCV) 
                           call mfix_exit(myPE)  
                        ENDIF 
                        IF (BC_TYPE(BCV)=='MASS_OUTFLOW' .AND.  .NOT.COMPARE(VEL&
                           ,BC_U_G(BCV))) THEN 
                           WRITE (UNIT_LOG, 1100) BCV, VEL, 'BC_U_g', BC_U_G(BCV) 
                           call mfix_exit(myPE)  
                        ENDIF 
                     ELSE 
                        IF (BC_TYPE(BCV) == 'MASS_INFLOW') THEN 
                           BC_U_G(BCV) = -VEL 
                        ELSE 
                           BC_U_G(BCV) = VEL 
                        ENDIF 
                     ENDIF 
                  CASE ('E')  
                     IF (BC_U_G(BCV) /= UNDEFINED) THEN 
                        IF (BC_TYPE(BCV)=='MASS_INFLOW' .AND.  .NOT.COMPARE(VEL,&
                           BC_U_G(BCV))) THEN 
                           WRITE (UNIT_LOG, 1100) BCV, VEL, 'BC_U_g', BC_U_G(BCV) 
                           call mfix_exit(myPE)  
                        ENDIF 
                        IF (BC_TYPE(BCV)=='MASS_OUTFLOW' .AND.  .NOT.COMPARE((-&
                           VEL),BC_U_G(BCV))) THEN 
                           WRITE (UNIT_LOG, 1100) BCV, (-VEL), 'BC_U_g', BC_U_G(&
                              BCV) 
                           call mfix_exit(myPE)  
                        ENDIF 
                     ELSE 
                        IF (BC_TYPE(BCV) == 'MASS_INFLOW') THEN 
                           BC_U_G(BCV) = VEL 
                        ELSE 
                           BC_U_G(BCV) = -VEL 
                        ENDIF 
                     ENDIF 
                  CASE ('S')  
                     IF (BC_V_G(BCV) /= UNDEFINED) THEN 
                        IF (BC_TYPE(BCV)=='MASS_INFLOW' .AND.  .NOT.COMPARE((-&
                           VEL),BC_V_G(BCV))) THEN 
                           WRITE (UNIT_LOG, 1100) BCV, (-VEL), 'BC_V_g', BC_V_G(&
                              BCV) 
                           call mfix_exit(myPE)  
                        ENDIF 
                        IF (BC_TYPE(BCV)=='MASS_OUTFLOW' .AND.  .NOT.COMPARE(VEL&
                           ,BC_V_G(BCV))) THEN 
                           WRITE (UNIT_LOG, 1100) BCV, VEL, 'BC_V_g', BC_V_G(BCV) 
                           call mfix_exit(myPE)  
                        ENDIF 
                     ELSE 
                        IF (BC_TYPE(BCV) == 'MASS_INFLOW') THEN 
                           BC_V_G(BCV) = -VEL 
                        ELSE 
                           BC_V_G(BCV) = VEL 
                        ENDIF 
                     ENDIF 
                  CASE ('N')  
                     IF (BC_V_G(BCV) /= UNDEFINED) THEN 
                        IF (BC_TYPE(BCV)=='MASS_INFLOW' .AND.  .NOT.COMPARE(VEL,&
                           BC_V_G(BCV))) THEN 
                           WRITE (UNIT_LOG, 1100) BCV, VEL, 'BC_V_g', BC_V_G(BCV) 
                           call mfix_exit(myPE)  
                        ENDIF 
                        IF (BC_TYPE(BCV)=='MASS_OUTFLOW' .AND.  .NOT.COMPARE((-&
                           VEL),BC_V_G(BCV))) THEN 
                           WRITE (UNIT_LOG, 1100) BCV, (-VEL), 'BC_V_g', BC_V_G(&
                              BCV) 
                           call mfix_exit(myPE)  
                        ENDIF 
                     ELSE 
                        IF (BC_TYPE(BCV) == 'MASS_INFLOW') THEN 
                           BC_V_G(BCV) = VEL 
                        ELSE 
                           BC_V_G(BCV) = -VEL 
                        ENDIF 
                     ENDIF 
                  CASE ('B')  
                     IF (BC_W_G(BCV) /= UNDEFINED) THEN 
                        IF (BC_TYPE(BCV)=='MASS_INFLOW' .AND.  .NOT.COMPARE((-&
                           VEL),BC_W_G(BCV))) THEN 
                           WRITE (UNIT_LOG, 1100) BCV, (-VEL), 'BC_W_g', BC_W_G(&
                              BCV) 
                           call mfix_exit(myPE)  
                        ENDIF 
                        IF (BC_TYPE(BCV)=='MASS_OUTFLOW' .AND.  .NOT.COMPARE(VEL&
                           ,BC_W_G(BCV))) THEN 
                           WRITE (UNIT_LOG, 1100) BCV, VEL, 'BC_W_g', BC_W_G(BCV) 
                           call mfix_exit(myPE)  
                        ENDIF 
                     ELSE 
                        IF (BC_TYPE(BCV) == 'MASS_INFLOW') THEN 
                           BC_W_G(BCV) = -VEL 
                        ELSE 
                           BC_W_G(BCV) = VEL 
                        ENDIF 
                     ENDIF 
                  CASE ('T')  
                     IF (BC_W_G(BCV) /= UNDEFINED) THEN 
                        IF (BC_TYPE(BCV)=='MASS_INFLOW' .AND.  .NOT.COMPARE(VEL,&
                           BC_W_G(BCV))) THEN 
                           WRITE (UNIT_LOG, 1100) BCV, VEL, 'BC_W_g', BC_W_G(BCV) 
                           call mfix_exit(myPE)  
                        ENDIF 
                        IF (BC_TYPE(BCV)=='MASS_OUTFLOW' .AND.  .NOT.COMPARE((-&
                           VEL),BC_W_G(BCV))) THEN 
                           WRITE (UNIT_LOG, 1100) BCV, (-VEL), 'BC_W_g', BC_W_G(&
                              BCV) 
                           call mfix_exit(myPE)  
                        ENDIF 
                     ELSE 
                        IF (BC_TYPE(BCV) == 'MASS_INFLOW') THEN 
                           BC_W_G(BCV) = VEL 
                        ELSE 
                           BC_W_G(BCV) = -VEL 
                        ENDIF 
                     ENDIF 
                  END SELECT 
               ENDIF 
!
!  Do flow conversions for solids phases
!
               DO M = 1, MMAX 
!
!             If solids mass flow is defined convert it to volumetric flow
!
                  IF (BC_MASSFLOW_S(BCV,M) /= UNDEFINED) THEN 
                     IF (RO_S(M) /= UNDEFINED) THEN 
                        VOLFLOW = BC_MASSFLOW_S(BCV,M)/RO_S(M) 
                     ELSE 
                        RETURN                   !  This error will be caught in a previous routine 
                     ENDIF 
!
!               If volumetric flow is also specified compare both
!
                     IF (BC_VOLFLOW_S(BCV,M) /= UNDEFINED) THEN 
                        IF (.NOT.COMPARE(VOLFLOW,BC_VOLFLOW_S(BCV,M))) THEN 
                           WRITE(UNIT_LOG,1200)BCV,VOLFLOW,M,BC_VOLFLOW_S(BCV,M) 
                           call mfix_exit(myPE)  
                        ENDIF 
                     ELSE 
                        BC_VOLFLOW_S(BCV,M) = VOLFLOW 
                     ENDIF 
                  ENDIF 
                  IF (BC_ROP_S(BCV,M)==UNDEFINED .AND. MMAX==1) BC_ROP_S(BCV,M)&
                      = (ONE - BC_EP_G(BCV))*RO_S(M) 
                  IF (BC_VOLFLOW_S(BCV,M) /= UNDEFINED) THEN 
                     IF (BC_ROP_S(BCV,M) /= UNDEFINED) THEN 
                        EPS = BC_ROP_S(BCV,M)/RO_S(M) 
                        IF (EPS /= ZERO) THEN 
                           VEL = BC_VOLFLOW_S(BCV,M)/(BC_AREA(BCV)*EPS) 
                        ELSE 
                           IF (BC_VOLFLOW_S(BCV,M) == ZERO) THEN 
                              VEL = ZERO 
                           ELSE 
                              WRITE (UNIT_LOG, 1250) BCV, M 
                              call mfix_exit(myPE)  
                           ENDIF 
                        ENDIF 
                     ELSE 
                        IF (BC_VOLFLOW_S(BCV,M) == ZERO) THEN 
                           VEL = ZERO 
                        ELSE 
                           WRITE (UNIT_LOG, 1260) BCV, M 
                           call mfix_exit(myPE)  
                        ENDIF 
                     ENDIF 
                     CONVERTED = .TRUE. 
                     SELECT CASE (BC_PLANE(BCV))  
                     CASE ('W')  
                        IF (BC_U_S(BCV,M) /= UNDEFINED) THEN 
                           IF (BC_TYPE(BCV)=='MASS_INFLOW' .AND.  .NOT.COMPARE((&
                              -VEL),BC_U_S(BCV,M))) THEN 
                              WRITE (UNIT_LOG, 1300) BCV, (-VEL), 'BC_U_s', M, &
                                 BC_U_S(BCV,M) 
                              call mfix_exit(myPE)  
                           ENDIF 
                           IF (BC_TYPE(BCV)=='MASS_OUTFLOW' .AND.  .NOT.COMPARE(&
                              VEL,BC_U_S(BCV,M))) THEN 
                              WRITE (UNIT_LOG, 1300) BCV, VEL, 'BC_U_s', M, &
                                 BC_U_S(BCV,M) 
                              call mfix_exit(myPE)  
                           ENDIF 
                        ELSE 
                           IF (BC_TYPE(BCV) == 'MASS_INFLOW') THEN 
                              BC_U_S(BCV,M) = -VEL 
                           ELSE 
                              BC_U_S(BCV,M) = VEL 
                           ENDIF 
                        ENDIF 
                     CASE ('E')  
                        IF (BC_U_S(BCV,M) /= UNDEFINED) THEN 
                           IF (BC_TYPE(BCV)=='MASS_INFLOW' .AND.  .NOT.COMPARE(&
                              VEL,BC_U_S(BCV,M))) THEN 
                              WRITE (UNIT_LOG, 1300) BCV, VEL, 'BC_U_s', M, &
                                 BC_U_S(BCV,M) 
                              call mfix_exit(myPE)  
                           ENDIF 
                           IF (BC_TYPE(BCV)=='MASS_OUTFLOW' .AND.  .NOT.COMPARE(&
                              (-VEL),BC_U_S(BCV,M))) THEN 
                              WRITE (UNIT_LOG, 1300) BCV, (-VEL), 'BC_U_s', M, &
                                 BC_U_S(BCV,M) 
                              call mfix_exit(myPE)  
                           ENDIF 
                        ELSE 
                           IF (BC_TYPE(BCV) == 'MASS_INFLOW') THEN 
                              BC_U_S(BCV,M) = VEL 
                           ELSE 
                              BC_U_S(BCV,M) = -VEL 
                           ENDIF 
                        ENDIF 
                     CASE ('S')  
                        IF (BC_V_S(BCV,M) /= UNDEFINED) THEN 
                           IF (BC_TYPE(BCV)=='MASS_INFLOW' .AND.  .NOT.COMPARE((&
                              -VEL),BC_V_S(BCV,M))) THEN 
                              WRITE (UNIT_LOG, 1300) BCV, (-VEL), 'BC_V_s', M, &
                                 BC_V_S(BCV,M) 
                              call mfix_exit(myPE)  
                           ENDIF 
                           IF (BC_TYPE(BCV)=='MASS_OUTFLOW' .AND.  .NOT.COMPARE(&
                              VEL,BC_V_S(BCV,M))) THEN 
                              WRITE (UNIT_LOG, 1300) BCV, VEL, 'BC_V_s', M, &
                                 BC_V_S(BCV,M) 
                              call mfix_exit(myPE)  
                           ENDIF 
                        ELSE 
                           IF (BC_TYPE(BCV) == 'MASS_INFLOW') THEN 
                              BC_V_S(BCV,M) = -VEL 
                           ELSE 
                              BC_V_S(BCV,M) = VEL 
                           ENDIF 
                        ENDIF 
                     CASE ('N')  
                        IF (BC_V_S(BCV,M) /= UNDEFINED) THEN 
                           IF (BC_TYPE(BCV)=='MASS_INFLOW' .AND.  .NOT.COMPARE(&
                              VEL,BC_V_S(BCV,M))) THEN 
                              WRITE (UNIT_LOG, 1300) BCV, VEL, 'BC_V_s', M, &
                                 BC_V_S(BCV,M) 
                              call mfix_exit(myPE)  
                           ENDIF 
                           IF (BC_TYPE(BCV)=='MASS_OUTFLOW' .AND.  .NOT.COMPARE(&
                              (-VEL),BC_V_S(BCV,M))) THEN 
                              WRITE (UNIT_LOG, 1300) BCV, (-VEL), 'BC_V_s', M, &
                                 BC_V_S(BCV,M) 
                              call mfix_exit(myPE)  
                           ENDIF 
                        ELSE 
                           IF (BC_TYPE(BCV) == 'MASS_INFLOW') THEN 
                              BC_V_S(BCV,M) = VEL 
                           ELSE 
                              BC_V_S(BCV,M) = -VEL 
                           ENDIF 
                        ENDIF 
                     CASE ('B')  
                        IF (BC_W_S(BCV,M) /= UNDEFINED) THEN 
                           IF (BC_TYPE(BCV)=='MASS_INFLOW' .AND.  .NOT.COMPARE((&
                              -VEL),BC_W_S(BCV,M))) THEN 
                              WRITE (UNIT_LOG, 1300) BCV, (-VEL), 'BC_W_s', M, &
                                 BC_W_S(BCV,M) 
                              call mfix_exit(myPE)  
                           ENDIF 
                           IF (BC_TYPE(BCV)=='MASS_OUTFLOW' .AND.  .NOT.COMPARE(&
                              VEL,BC_W_S(BCV,M))) THEN 
                              WRITE (UNIT_LOG, 1300) BCV, VEL, 'BC_W_s', M, &
                                 BC_W_S(BCV,M) 
                              call mfix_exit(myPE)  
                           ENDIF 
                        ELSE 
                           IF (BC_TYPE(BCV) == 'MASS_INFLOW') THEN 
                              BC_W_S(BCV,M) = -VEL 
                           ELSE 
                              BC_W_S(BCV,M) = VEL 
                           ENDIF 
                        ENDIF 
                     CASE ('T')  
                        IF (BC_W_S(BCV,M) /= UNDEFINED) THEN 
                           IF (BC_TYPE(BCV)=='MASS_INFLOW' .AND.  .NOT.COMPARE(&
                              VEL,BC_W_S(BCV,M))) THEN 
                              WRITE (UNIT_LOG, 1300) BCV, VEL, 'BC_W_s', M, &
                                 BC_W_S(BCV,M) 
                              call mfix_exit(myPE)  
                           ENDIF 
                           IF (BC_TYPE(BCV)=='MASS_OUTFLOW' .AND.  .NOT.COMPARE(&
                              (-VEL),BC_W_S(BCV,M))) THEN 
                              WRITE (UNIT_LOG, 1300) BCV, (-VEL), 'BC_W_s', M, &
                                 BC_W_S(BCV,M) 
                              call mfix_exit(myPE)  
                           ENDIF 
                        ELSE 
                           IF (BC_TYPE(BCV) == 'MASS_INFLOW') THEN 
                              BC_W_S(BCV,M) = VEL 
                           ELSE 
                              BC_W_S(BCV,M) = -VEL 
                           ENDIF 
                        ENDIF 
                     END SELECT 
                  ENDIF 
               END DO 
            ENDIF 
         ENDIF 
      END DO 
      IF (CONVERTED .AND. (NO_I .OR. NO_J .OR. NO_K)) WRITE (UNIT_LOG, 1500) 
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
         ' Non-zero vol. or mass flow specified with BC_ROP_s',I1,' = 0.',/1X,&
         70('*')/) 
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
