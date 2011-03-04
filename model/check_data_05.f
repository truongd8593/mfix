!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CHECK_DATA_05                                          C
!  Purpose: check the gas phase input section                          C
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
!  Variables referenced: MU_g0, MW_AVG, RO_g0                          C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables: None                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CHECK_DATA_05 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE physprop
      USE funits 
      USE run
      USE indices
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
      INTEGER :: N 
!-----------------------------------------------


! CHECK NMAX(0) 
      IF (NMAX(0) == UNDEFINED_I) THEN 
         IF (SPECIES_EQ(0)) THEN 
            CALL ERROR_ROUTINE ('CHECK_DATA_05', &
               'Number of gas species (NMAX(0)) not specified',0,2) 
            CALL ERROR_ROUTINE (' ', ' ', 1, 3) 
         ELSE 
            NMAX(0) = 1 
         ENDIF 
      ELSEIF (NMAX(0) > DIMENSION_N_G) THEN 
         CALL ERROR_ROUTINE ('CHECK_DATA_05',&
            'NMAX(0) is too large',0,2)
         IF(DMP_LOG)WRITE (UNIT_LOG, 1000) NMAX(0), DIMENSION_N_G 
         CALL ERROR_ROUTINE (' ', ' ', 1, 3) 
      ENDIF 

! CHECK MW_AVG
      IF (SPECIES_EQ(0)) THEN 
         IF (MW_AVG /= UNDEFINED) THEN 
            IF(DMP_LOG)WRITE (UNIT_LOG, 1410) 
            MW_AVG = UNDEFINED 
         ENDIF 
      ELSE 
         IF (RO_G0 == UNDEFINED) THEN 
            IF (MW_AVG <= ZERO) THEN 
               CALL ERROR_ROUTINE ('CHECK_DATA_05',&
                  'MW_AVG is unphysical',0,2) 
               IF(DMP_LOG)WRITE (UNIT_LOG, 1200) MW_AVG 
               CALL ERROR_ROUTINE (' ', ' ', 1, 3) 
            ENDIF 
         ELSE 
            IF (RO_G0 < ZERO) THEN 
               CALL ERROR_ROUTINE ('CHECK_DATA_05', &
                  'Value of RO_g0 is unphysical', 0, 2) 
               IF(DMP_LOG)WRITE (UNIT_LOG, 1300) RO_G0 
               CALL ERROR_ROUTINE (' ', ' ', 1, 3) 
            ENDIF 
            IF (MW_AVG /= UNDEFINED .AND. DMP_LOG)WRITE (UNIT_LOG, 1400) 
         ENDIF 
      ENDIF 
      
! CHECK MW_g
      IF (SPECIES_EQ(0) .OR. MW_AVG==UNDEFINED .AND. RO_G0==UNDEFINED) THEN 
         DO N = 1, NMAX(0) 
            IF (MW_G(N) == UNDEFINED) THEN 
               CALL ERROR_ROUTINE ('CHECK_DATA_05', &
                  'Value of MW_g not specified', 0, 2) 
               IF(DMP_LOG)WRITE (UNIT_LOG, 1500) N 
! No need to abort since MW will be read from database
!               CALL ERROR_ROUTINE (' ', ' ', 1, 3)                
            ENDIF 
         ENDDO 
         DO N = NMAX(0) + 1, DIMENSION_N_G 
            IF (MW_G(N) /= UNDEFINED) THEN 
               CALL ERROR_ROUTINE ('CHECK_DATA_05', &
                  'MW_g defined for N > NMAX(0)', 0, 2) 
               IF(DMP_LOG)WRITE (UNIT_LOG, 1501) N 
               CALL ERROR_ROUTINE (' ', ' ', 1, 3) 
            ENDIF 
         ENDDO 
      ENDIF 

! CHECK MU_g0
      IF (MU_G0 <= ZERO) THEN 
         CALL ERROR_ROUTINE ('CHECK_DATA_05', &
            'MU_g0 value is unphysical', 0, 2) 
         IF(DMP_LOG)WRITE (UNIT_LOG, 1100) MU_G0 
         CALL ERROR_ROUTINE (' ', ' ', 1, 3) 
      ENDIF 

! CHECK K_g0
      IF (K_G0 < ZERO) THEN 
         CALL ERROR_ROUTINE ('CHECK_DATA_05', &
            'K_g0 value is unphysical', 0, 2) 
         IF(DMP_LOG)WRITE (UNIT_LOG, 1110) K_G0 
         CALL ERROR_ROUTINE (' ', ' ', 1, 3) 
      ENDIF 

! CHECK C_pg0
      IF (C_PG0 < ZERO) THEN 
         CALL ERROR_ROUTINE ('CHECK_DATA_05',&
            'C_pg0 value is unphysical', 0, 2) 
         IF(DMP_LOG)WRITE (UNIT_LOG, 1120) C_PG0 
         CALL ERROR_ROUTINE (' ', ' ', 1, 3) 
      ENDIF 

! CHECK DIF_g0
      IF (DIF_G0 < ZERO) THEN 
         CALL ERROR_ROUTINE ('CHECK_DATA_05',&
            'DIF_g0 value is unphysical', 0, 2) 
         IF(DMP_LOG)WRITE (UNIT_LOG, 1130) DIF_G0 
         CALL ERROR_ROUTINE (' ', ' ', 1, 3) 
      ENDIF 


      RETURN  

 1000 FORMAT(1X,/,1X,'NMAX(0) in mfix.dat = ',I3,/,1X,&
         'DIMENSION_N_g in param.inc = ',I3) 
 1100 FORMAT(1X,/,1X,'MU_g0   in mfix.dat = ',G12.5) 
 1110 FORMAT(1X,/,1X,'K_g0   in mfix.dat = ',G12.5) 
 1120 FORMAT(1X,/,1X,'C_pg0   in mfix.dat = ',G12.5) 
 1130 FORMAT(1X,/,1X,'DIF_g0   in mfix.dat = ',G12.5)          
 1200 FORMAT(1X,/,1X,'MW_AVG in mfix.dat = ',G12.5) 
 1300 FORMAT(1X,/,1X,'RO_g0   in mfix.dat = ',G12.5) 
 1400 FORMAT(/1X,70('*')//' From: CHECK_DATA_05',/&
         ' Message: Since RO_g0 is specified MW_AVG will be ignored',/, 1X&
         ,70('*')/) 
 1410 FORMAT(/1X,70('*')//' From: CHECK_DATA_05',/&
         ' Message: Since gas phase rxns are specified MW_AVG',&
         ' will be ignored',/1X,70('*')/) 
 1500 FORMAT(1X,/,1X,'MW_g for gas species ',I3,' not specified') 
 1501 FORMAT(1X,/,1X,'MW_g for gas species ',I3,' specified') 

      END SUBROUTINE CHECK_DATA_05 
