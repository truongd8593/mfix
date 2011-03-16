!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CHECK_DATA_04                                          C
!  Purpose: check the solid phase input section                        C
!                                                                      C
!  Author: P. Nicoletti                               Date: 02-DEC-91  C
!  Reviewer: M.SYAMLAL, W.ROGERS, P.NICOLETTI         Date: 24-JAN-92  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: MMAX, D_p, RO_s, EP_star                      C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables: LC                                                 C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CHECK_DATA_04 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE run
      USE indices
      USE physprop
      USE constant
      USE discretelement
      USE funits 
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
      INTEGER :: LC, N
      CHARACTER*85 LONG_STRING      
!-----------------------------------------------


! Check MMAX
      IF (MMAX<0 .OR. MMAX>DIMENSION_M) THEN 
         CALL ERROR_ROUTINE ('check_data_04', &
            'MMAX not specified or unphysical', 0, 2) 
            IF(DMP_LOG)WRITE (UNIT_LOG, 1000) MMAX, DIMENSION_M 
         CALL ERROR_ROUTINE (' ', ' ', 1, 3) 
      ENDIF 

      IF (MMAX > 0) THEN 
         IF (C_E == UNDEFINED) CALL ERROR_ROUTINE ('check_data_04', &
            'Coefficient of restitution (C_e) not specified', 1, 1) 
         IF (C_F==UNDEFINED .AND. MMAX>=2 .AND. KT_TYPE .EQ. UNDEFINED_C) &
            CALL ERROR_ROUTINE ('check_data_04',&
               'Coefficient of friction (C_f) not specified',1,1) 

         IF ((FRICTION .OR. SCHAEFFER) .AND. (PHI == UNDEFINED)) &
            CALL ERROR_ROUTINE ('check_data_04', &
               'Angle of internal friction (Phi) not specified',1,1)
         LONG_STRING = 'Angle of wall-particle friction (Phi_w) &
            &not specified'
         IF ((FRICTION .OR. JENKINS) .AND. (PHI_W == UNDEFINED)) &
             CALL ERROR_ROUTINE ('check_data_04',TRIM(LONG_STRING),1,1)
      ENDIF 

! Check D_p0
! Only need to check for real phases (for GHD theory)      
      DO LC = 1, SMAX 
         IF (D_P0(LC)<ZERO .OR. D_P0(LC)==UNDEFINED) THEN 
            CALL ERROR_ROUTINE ('check_data_04', &
               'D_p0 not specified or unphysical', 0, 2) 
               IF(DMP_LOG)WRITE (UNIT_LOG, 1100) LC, D_P0(LC) 
            CALL ERROR_ROUTINE (' ', ' ', 1, 3) 
         ENDIF 
      ENDDO 
      DO LC = SMAX + 1, DIMENSION_M 
         IF (D_P0(LC) /= UNDEFINED) THEN 
            CALL ERROR_ROUTINE ('check_data_04', &
               'too many D_p0 values specified', 0, 2) 
               IF(DMP_LOG)WRITE (UNIT_LOG, 1200) LC, D_P0(LC), MMAX
            CALL ERROR_ROUTINE (' ', ' ', 1, 3) 
         ENDIF 
      ENDDO 

! Check RO_s      
! Only need to check for real phases (for GHD theory) 
      DO LC = 1, SMAX 
         IF (RO_S(LC)<ZERO .OR. RO_S(LC)==UNDEFINED) THEN 
            CALL ERROR_ROUTINE ('check_data_04', &
               'RO_s not specified or unphysical', 0, 2) 
               IF(DMP_LOG)WRITE (UNIT_LOG, 1300) LC, RO_S(LC) 
            CALL ERROR_ROUTINE (' ', ' ', 1, 3) 
         ENDIF 
      ENDDO 
      DO LC = SMAX + 1, DIMENSION_M 
         IF (RO_S(LC) /= UNDEFINED) THEN 
            CALL ERROR_ROUTINE ('check_data_04', &
               'too many RO_s values specified', 0, 2) 
               IF(DMP_LOG)WRITE (UNIT_LOG, 1400) LC, RO_S(LC), MMAX 
            CALL ERROR_ROUTINE (' ', ' ', 1, 3) 
         ENDIF 
      ENDDO 

      IF (MODEL_B) THEN 
         DO LC = 1, MMAX 
            IF (.NOT.CLOSE_PACKED(LC)) THEN 
               CALL ERROR_ROUTINE ('check_data_04', ' ', 0, 2) 
                  IF(DMP_LOG)WRITE (UNIT_LOG, 1420) LC
               CALL ERROR_ROUTINE (' ', ' ', 1, 3) 
            ENDIF 
         ENDDO 
      ENDIF 

! Check EP_star
      IF (MMAX > 0) THEN 
         IF (EP_STAR<ZERO .OR. EP_STAR>ONE) &
            CALL ERROR_ROUTINE ('check_data_04',&
            'Value of EP_star is unphysical', 1, 1) 
      ENDIF 

! Yu_Standish and Fedors_Landel correlations are used with more than one solids phase
      LONG_STRING = 'MMAX must be >= 2 for Yu_Standish or &
         &Fedors_Landel correlations'
      IF (SMAX < 2 .AND. (YU_STANDISH .OR. FEDORS_LANDEL)) &
         CALL ERROR_ROUTINE ('check_data_04',TRIM(LONG_STRING),1,1)

! Fedors_Landel correlations is limited to a binary mixture of powders
      IF (SMAX > 2 .AND. FEDORS_LANDEL) &
         CALL ERROR_ROUTINE ('check_data_04', &
            'Fedors_Landel requires MMAX = 2 ', 1, 1)

! Must choose between Yu_Standish and Fedors_Landel correlations, can't use both.
      LONG_STRING = 'Cannot use both Yu_Standish and Fedors_Landel &
         &correlations'
      IF (YU_STANDISH .AND. FEDORS_LANDEL) &
         CALL ERROR_ROUTINE ('check_data_04',TRIM(LONG_STRING),1,1)
      IF (SIGM_BLEND) THEN
         TANH_BLEND = .FALSE. ! Setting tanh blending to be false
      ENDIF

! Define restitution coefficient matrix      
      DO LC = 1, SMAX 
         DO N = 1, SMAX
            IF(r_p(LC,N) == UNDEFINED) r_p(LC,N) = C_e
            r_p(N,LC) = r_p(LC,N) ! just need to define r_p(1,2) and r_p(2,1) will be set.
         ENDDO
      ENDDO


! Check NMAX - nmax(M) must be defined if species_eq(M) true
! Only need to check species for real phases (for GHD theory) 
      DO LC = 1, MMAX 
         IF (NMAX(LC) == UNDEFINED_I) THEN 
            IF (SPECIES_EQ(LC)) THEN 
               CALL ERROR_ROUTINE ('check_data_04', &
                  'Number of species not specified', 0, 2) 
                  IF(DMP_LOG)WRITE (UNIT_LOG, 1045) LC 
               CALL ERROR_ROUTINE (' ', ' ', 1, 3) 
            ELSE 
               NMAX(LC) = 1 
            ENDIF 
         ENDIF 
         IF (NMAX(LC) > DIMENSION_N_S) THEN 
            CALL ERROR_ROUTINE ('check_data_04', 'NMAX is too large', 0, 2) 
               IF(DMP_LOG)WRITE (UNIT_LOG, 1050) LC, NMAX, DIMENSION_N_S 
            CALL ERROR_ROUTINE (' ', ' ', 1, 3) 
         ENDIF 
      ENDDO 

! Check MW_s if solids species are present    
      DO LC = 1, SMAX 
         IF (SPECIES_EQ(LC)) THEN 
            DO N = 1, NMAX(LC) 
               IF (MW_S(LC,N) == UNDEFINED) THEN 
                  CALL ERROR_ROUTINE ('check_data_04', &
                     'Species molecular weight undefined', 0, 2) 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1410) LC, N 
                  !CALL ERROR_ROUTINE (' ', ' ', 1, 3) no need to abort as they will be read from database
               ENDIF 
            ENDDO 
            DO N = NMAX(LC) + 1, DIMENSION_N_S 
               IF (MW_S(LC,N) /= UNDEFINED) THEN 
                  CALL ERROR_ROUTINE ('check_data_04', &
                     'MW_s defined for N > NMAX(m)', 0, 2) 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1410) LC, N 
                  CALL ERROR_ROUTINE (' ', ' ', 1, 3) 
               ENDIF 
            ENDDO 
         ENDIF 
      ENDDO 

! CHECK MU_s0
      IF (MU_S0 < ZERO) THEN 
         CALL ERROR_ROUTINE ('CHECK_DATA_04',&
            'MU_s0 value is unphysical', 0, 2) 
         IF(DMP_LOG)WRITE (UNIT_LOG, 1500) MU_S0 
         CALL ERROR_ROUTINE (' ', ' ', 1, 3) 
      ENDIF 

! CHECK K_s0
      IF (K_S0 < ZERO) THEN 
         CALL ERROR_ROUTINE ('CHECK_DATA_04',&
            'K_s0 value is unphysical', 0, 2) 
         IF(DMP_LOG)WRITE (UNIT_LOG, 1510) K_S0 
         CALL ERROR_ROUTINE (' ', ' ', 1, 3) 
      ENDIF 
      
! CHECK C_ps0
      IF (C_PS0 < ZERO) THEN 
         CALL ERROR_ROUTINE ('CHECK_DATA_04', &
            'C_ps0 value is unphysical', 0, 2) 
         IF(DMP_LOG)WRITE (UNIT_LOG, 1520) C_PS0 
         CALL ERROR_ROUTINE (' ', ' ', 1, 3) 
      ENDIF 

! CHECK DIF_s0
      IF (DIF_S0 < ZERO) THEN 
         CALL ERROR_ROUTINE ('CHECK_DATA_04',&
            'DIF_s0 value is unphysical', 0, 2) 
         IF(DMP_LOG)WRITE (UNIT_LOG, 1530) DIF_S0 
         CALL ERROR_ROUTINE (' ', ' ', 1, 3) 
      ENDIF 


      RETURN  

 1000 FORMAT(1X,/,1X,'MMAX        in  mfix.dat = ',I6,/,1X,&
         'DIMENSION_M in  param.inc  = ',I6,/) 
 1045 FORMAT(1X,/,1X,'NMAX is not specified for solids phase',I2) 
 1050 FORMAT(1X,/,1X,'NMAX(',I2,')   in  mfix.dat = ',I6,/,1X,&
         'DIMENSION_N_s in  param.inc  = ',I6,/) 
 1100 FORMAT(1X,/,1X,'D_p0(',I2,') in mfix.dat = ',G12.5) 
 1200 FORMAT(1X,/,1X,'D_p0(',I2,') = ',G12.5,/,1X,'MMAX in mfix = ',I2,/) 
 1300 FORMAT(1X,/,1X,'RO_s(',I2,') in mfix.dat = ',G12.5) 
 1400 FORMAT(1X,/,1X,'RO_s(',I2,') = ',G12.5,/,1X,'MMAX in mfix = ',I2,/) 
 1410 FORMAT(1X,/,1X,'Solids phase = ',I2,'   Species = ',I3) 
 1420 FORMAT(1X,/,1X,'Solids phase = ',I2,' is not Close_Packed.',/,&
         ' With Model B all solids phases should have that property')

 1500 FORMAT(1X,/,1X,'MU_s0   in mfix.dat = ',G12.5)
 1510 FORMAT(1X,/,1X,'K_s0   in mfix.dat = ',G12.5)
 1520 FORMAT(1X,/,1X,'C_ps0   in mfix.dat = ',G12.5)
 1530 FORMAT(1X,/,1X,'DIF_s0   in mfix.dat = ',G12.5)

      END SUBROUTINE CHECK_DATA_04 
