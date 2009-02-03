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
!-----------------------------------------------
!
!
! Overwrite any user's DES logicals set with discrete_element = .false.
!
      IF(.NOT.DISCRETE_ELEMENT) THEN
         DES_CONTINUUM_COUPLED = .FALSE.
         DES_INTERP_ON = .FALSE.
         TSUJI_DRAG = .FALSE.
         WALLDTSPLIT = .FALSE.
         PRINT_DES_DATA = .FALSE.
      ENDIF
!
!
! CHECK MMAX
!
      IF (MMAX<0 .OR. MMAX>DIMENSION_M) THEN 
         CALL ERROR_ROUTINE ('check_data_04', &
            'MMAX not specified or unphysical', 0, 2) 
            IF(DMP_LOG)WRITE (UNIT_LOG, 1000) MMAX, DIMENSION_M 
         CALL ERROR_ROUTINE (' ', ' ', 1, 3) 
      ENDIF 
!
! CHECK NMAX
!
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
      END DO 
      IF (MMAX > 0) THEN 
         IF (C_E == UNDEFINED) CALL ERROR_ROUTINE ('check_data_04', &
            'Coefficient of restitution (C_e) not specified', 1, 1) 
         IF (C_F==UNDEFINED .AND. MMAX>=2) CALL ERROR_ROUTINE ('check_data_04'&
            , 'Coefficient of friction (C_f) not specified', 1, 1) 
         IF (PHI == UNDEFINED) CALL ERROR_ROUTINE ('check_data_04', &
            'Angle of internal friction (Phi) not specified', 1, 1) 
         IF (PHI_W == UNDEFINED) CALL ERROR_ROUTINE ('check_data_04', &
            'Angle of wall-particle friction (Phi_w) not specified', 1, 1) 
      ENDIF 
!
! CHECK D_p
!
      DO LC = 1, MMAX 
         IF (D_P0(LC)<ZERO .OR. D_P0(LC)==UNDEFINED) THEN 
            CALL ERROR_ROUTINE ('check_data_04', &
               'D_p0 not specified or unphysical', 0, 2) 
               IF(DMP_LOG)WRITE (UNIT_LOG, 1100) LC, D_P0(LC) 
            CALL ERROR_ROUTINE (' ', ' ', 1, 3) 
         ENDIF 
      END DO 
      DO LC = MMAX + 1, DIMENSION_M 
         IF (D_P0(LC) /= UNDEFINED) THEN 
            CALL ERROR_ROUTINE ('check_data_04', &
               'too many D_p0 values specified', 0, 2) 
               IF(DMP_LOG)WRITE (UNIT_LOG, 1200) LC, D_P0(LC), MMAX
            CALL ERROR_ROUTINE (' ', ' ', 1, 3) 
         ENDIF 
      END DO 
      DO LC = 1, MMAX 
         IF (RO_S(LC)<ZERO .OR. RO_S(LC)==UNDEFINED) THEN 
            CALL ERROR_ROUTINE ('check_data_04', &
               'RO_s not specified or unphysical', 0, 2) 
               IF(DMP_LOG)WRITE (UNIT_LOG, 1300) LC, RO_S(LC) 
            CALL ERROR_ROUTINE (' ', ' ', 1, 3) 
         ENDIF 
      END DO 
      DO LC = MMAX + 1, DIMENSION_M 
         IF (RO_S(LC) /= UNDEFINED) THEN 
            CALL ERROR_ROUTINE ('check_data_04', &
               'too many RO_s values specified', 0, 2) 
               IF(DMP_LOG)WRITE (UNIT_LOG, 1400) LC, RO_S(LC), MMAX 
            CALL ERROR_ROUTINE (' ', ' ', 1, 3) 
         ENDIF 
      END DO 
      DO LC = 1, MMAX 
         IF (SPECIES_EQ(LC)) THEN 
            DO N = 1, NMAX(LC) 
               IF (MW_S(LC,N) == UNDEFINED) THEN 
                  CALL ERROR_ROUTINE ('check_data_04', &
                     'Species molecular weight undefined', 0, 2) 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1410) LC, N 
                  !CALL ERROR_ROUTINE (' ', ' ', 1, 3) no need to abort as they will be read from database
               ENDIF 
            END DO 
            DO N = NMAX(LC) + 1, DIMENSION_N_S 
               IF (MW_S(LC,N) /= UNDEFINED) THEN 
                  CALL ERROR_ROUTINE ('check_data_04', &
                     'MW_s defined for N > NMAX(m)', 0, 2) 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1410) LC, N 
                  CALL ERROR_ROUTINE (' ', ' ', 1, 3) 
               ENDIF 
            END DO 
         ENDIF 
      END DO 
      IF (MODEL_B) THEN 
!
         DO LC = 1, MMAX 
!
            IF (.NOT.CLOSE_PACKED(LC)) THEN 
               CALL ERROR_ROUTINE ('check_data_04', ' ', 0, 2) 
                  IF(DMP_LOG)WRITE (UNIT_LOG, 1420) LC
               CALL ERROR_ROUTINE (' ', ' ', 1, 3) 
!
            ENDIF 
         END DO 
      ENDIF 
!
!
! check EP_star
!
      IF (MMAX > 0) THEN 
         IF (EP_STAR<ZERO .OR. EP_STAR>ONE) CALL ERROR_ROUTINE ('check_data_04'&
            , 'Value of EP_star is unphysical', 1, 1) 
      ENDIF 
!
! Yu_Standish and Fedors_Landel correlations are used with more than one solids phase
!
      IF (MMAX < 2 .AND. (YU_STANDISH .OR. FEDORS_LANDEL)) &
         CALL ERROR_ROUTINE ('check_data_04', &
       'MMAX must be >= 2 for Yu_Standish or Fedors_Landel correlations', 1, 1)
!
! Fedors_Landel correlations is limited to a binary mixture of powders
!
      IF (MMAX > 2 .AND. FEDORS_LANDEL) &
         CALL ERROR_ROUTINE ('check_data_04' &
         , 'Fedors_Landel requires MMAX = 2 ', 1, 1)
!
! Must choose between Yu_Standish and Fedors_Landel correlations, can't use both.
!
      IF (YU_STANDISH .AND. FEDORS_LANDEL) &
         CALL ERROR_ROUTINE ('check_data_04', &
         'Cannot use both Yu_Standish and Fedors_Landel correlations', 1, 1)
!
      RETURN  
!
      IF (SIGM_BLEND) THEN
         TANH_BLEND = .FALSE. ! Setting tanh blending to be false
      ENDIF
!
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
!
      END SUBROUTINE CHECK_DATA_04 
