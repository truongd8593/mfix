!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CHECK_DATA_02                                          C
!  Purpose: check the output control namelist section                  C
!                                                                      C
!  Author: P. Nicoletti                               Date: 27-NOV-91  C
!  Reviewer: M.SYAMLAL, W.ROGERS, P.NICOLETTI         Date: 27-JAN-92  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: RES_DT, SPX_DT, OUT_DT, NLOG                  C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables: LC                                                 C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CHECK_DATA_02 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE output
      USE leqsol 
      USE geometry
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
!             loop counter
      INTEGER :: LC 
      LOGICAL :: MSG 
      CHARACTER :: LINE*80 
!-----------------------------------------------
!
!
      IF (RES_DT==UNDEFINED .OR. RES_DT<=ZERO) CALL ERROR_ROUTINE (&
         'check_data_02', 'RES_DT not specified OR RES_DT <= 0 in mfix.dat', 1&
         , 1) 
!
      DO LC = 1, N_SPX 
         IF (SPX_DT(LC)==UNDEFINED .OR. SPX_DT(LC)<=ZERO) CALL ERROR_ROUTINE (&
            'check_data_02', 'SPX_DT not specified OR SPX_DT <= 0 in mfix.dat'&
            , 1, 1) 
      END DO 
      IF (NLOG == 0) CALL ERROR_ROUTINE ('check_data_02', &
         'NLOG = 0 in mfix.dat', 1, 1) 
!
!  If cyclic conditions are used the GMRES should be used as the linear
!  equation solver
!
!      IF (COORDINATES=='CYLINDRICAL' .AND.  .NOT.NO_K .OR. CYCLIC_X .OR. &
!         CYCLIC_X_PD .OR. CYCLIC_Y .OR. CYCLIC_Y_PD .OR. CYCLIC_Z .OR. &
!         CYCLIC_Z_PD) THEN 
!
!         MSG = .FALSE. 
!         DO LC = 1, 8 
!            IF (LEQ_METHOD(LC) == 2) THEN 
!               LEQ_METHOD(LC) = 3 
!               MSG = .TRUE. 
!            ENDIF 
!         END DO 
!         IF (MSG) THEN 
!            WRITE (LINE, '(A)') 'Warning: LEQ_METHOD changed from 2 to 3' 
!            CALL WRITE_ERROR ('CHECK_DATA_02', LINE, 1) 
!         ENDIF 
!      ENDIF 
!
      RETURN  
      END SUBROUTINE CHECK_DATA_02 
