!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_OUT1                                             C
!  Purpose: write out the field variables to standard output           C
!                                                                      C
!  Author: P. Nicoletti                               Date: 03-DEC-91  C
!  Reviewer: W. Rogers, M. Syamlal, S. Venkatesan     Date: 31-JAN-92  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: TIME, P_g, EP_g, RO_g, ROP_g, MMAX, ROP_s     C
!                        T_g, T_s, U_g, V_g, W_g, U_s, V_s, W_s C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables: LC, N                                              C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE WRITE_OUT1 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE physprop
      USE fldvar
      USE run
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
!             form feed character = CHAR(12)
!
      WRITE (UNIT_OUT, 1000) CHAR(12), TIME 
      CALL OUT_ARRAY (P_G, 'P_g') 
!
      WRITE (UNIT_OUT, 1050) CHAR(12), TIME 
      CALL OUT_ARRAY (P_STAR, 'P_star') 
!
      WRITE (UNIT_OUT, 1100) CHAR(12), TIME 
      CALL OUT_ARRAY (EP_G, 'EP_g') 
!
      WRITE (UNIT_OUT, 1200) CHAR(12), TIME 
      CALL OUT_ARRAY (RO_G, 'RO_g') 
!
      DO LC = 1, MMAX 
         WRITE (UNIT_OUT, 1400) CHAR(12), LC, TIME 
         CALL OUT_ARRAY (ROP_S(1,LC), 'ROP_s') 
      END DO 
      WRITE (UNIT_OUT, 1500) CHAR(12), TIME 
      CALL OUT_ARRAY (T_G, 'T_g') 
!
      DO LC = 1, MMAX 
         WRITE (UNIT_OUT, 1600) CHAR(12), LC, TIME 
         CALL OUT_ARRAY (T_S(1,LC), 'T_s') 
      END DO 
      IF (SPECIES_EQ(0)) THEN 
         DO N = 1, NMAX(0) 
            WRITE (UNIT_OUT, 1710) CHAR(12), N, TIME 
            CALL OUT_ARRAY (X_G(1,N), 'X_g') 
         END DO 
      ENDIF 
!
      DO LC = 1, MMAX 
         IF (SPECIES_EQ(LC)) THEN 
            DO N = 1, NMAX(LC) 
               WRITE (UNIT_OUT, 1720) CHAR(12), LC, N, TIME 
               CALL OUT_ARRAY (X_S(1,LC,N), 'X_s') 
            END DO 
         ENDIF 
      END DO 
      WRITE (UNIT_OUT, 1800) CHAR(12), TIME 
      CALL OUT_ARRAY (U_G, 'U_g') 
!
      WRITE (UNIT_OUT, 1900) CHAR(12), TIME 
      CALL OUT_ARRAY (V_G, 'V_g') 
!
      WRITE (UNIT_OUT, 2000) CHAR(12), TIME 
      CALL OUT_ARRAY (W_G, 'W_g') 
!
      DO LC = 1, MMAX 
         WRITE (UNIT_OUT, 2100) CHAR(12), LC, TIME 
         CALL OUT_ARRAY (U_S(1,LC), 'U_s') 
!
         WRITE (UNIT_OUT, 2200) CHAR(12), LC, TIME 
         CALL OUT_ARRAY (V_S(1,LC), 'V_s') 
!
         WRITE (UNIT_OUT, 2300) CHAR(12), LC, TIME 
         CALL OUT_ARRAY (W_S(1,LC), 'W_s') 
!
!         IF(GRANULAR_ENERGY)THEN
         WRITE (UNIT_OUT, 2400) CHAR(12), LC, TIME 
         CALL OUT_ARRAY (THETA_M(1,LC), 'Theta_m') 
      END DO 
      WRITE (UNIT_OUT, '(/1X,1A1)') CHAR(12) 
      IF (CALL_USR) CALL USR_WRITE_OUT1 
      RETURN  
 1000 FORMAT(1X,A1,/5X,'--- Gas pressure (P_g) at time ',G12.5,' ---',2/) 
 1050 FORMAT(1X,A1,/5X,'--- Solids pressure (P_star) at time ',G12.5,' ---',2/) 
 1100 FORMAT(1X,A1,/5X,'--- Void fraction (EP_g) at time ',G12.5,' ---',2/) 
 1200 FORMAT(1X,A1,/5X,'--- Gas density (RO_g) at time ',G12.5,' ---',2/) 
 1400 FORMAT(1X,A1,/5X,'--- Solids Phase-',I1,' density x volume',&
         ' fraction (ROP_s) at time ',G12.5,' ---',2/) 
 1500 FORMAT(1X,A1,/5X,'--- Gas temperature (T_g) at time ',G12.5,' ---',2/) 
 1600 FORMAT(1X,A1,/5X,'--- Solids Phase-',I1,' temperature (T_s)',' at time ',&
         G12.5,' ---',2/) 
 1710 FORMAT(1X,A1,/5X,'--- Mass fraction of gas species (X_g) ',I2,' at time '&
         ,G12.5,' ---',2/) 
 1720 FORMAT(1X,A1,/5X,'--- Mass fraction of solids-',I1,' species (X_s)',I2,&
         ' at time ',G12.5,' ---',2/) 
 1800 FORMAT(1X,A1,/5X,'--- X-component of gas velocity (U_g) at time ',G12.5,&
         ' ---',2/) 
 1900 FORMAT(1X,A1,/5X,'--- Y-component of gas velocity (V_g) at time ',G12.5,&
         ' ---',2/) 
 2000 FORMAT(1X,A1,/5X,'--- Z-component of gas velocity (W_g) at time ',G12.5,&
         ' ---',2/) 
 2100 FORMAT(1X,A1,/5X,'--- X-component of Solids Phase-',I1,&
         ' velocity (U_s) at time ',G12.5,' ---',2/) 
 2200 FORMAT(1X,A1,/5X,'--- Y-component of Solids Phase-',I1,&
         ' velocity (V_s) at time ',G12.5,' ---',2/) 
 2300 FORMAT(1X,A1,/5X,'--- Z-component of Solids Phase-',I1,&
         ' velocity (W_s) at time ',G12.5,' ---',2/) 
 2400 FORMAT(1X,A1,/5X,'--- Granular temperature of Solids Phase-',I1,&
         ' velocity (W_s) at time ',G12.5,' ---',2/) 
      END SUBROUTINE WRITE_OUT1 
