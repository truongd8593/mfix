CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: WRITE_OUT1                                             C
C  Purpose: write out the field variables to standard output           C
C                                                                      C
C  Author: P. Nicoletti                               Date: 03-DEC-91  C
C  Reviewer: W. Rogers, M. Syamlal, S. Venkatesan     Date: 31-JAN-92  C
C                                                                      C
C  Revision Number:                                                    C
C  Purpose:                                                            C
C  Author:                                            Date: dd-mmm-yy  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Literature/Document References:                                     C
C                                                                      C
C  Variables referenced: TIME, P_g, EP_g, RO_g, ROP_g, MMAX, ROP_s     C
C                        T_g, T_s, U_g, V_g, W_g, U_s, V_s, W_s C
C  Variables modified: None                                            C
C                                                                      C
C  Local variables: LC, N                                              C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
C
      SUBROUTINE WRITE_OUT1
C
      IMPLICIT NONE
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE 'physprop.inc'
      INCLUDE 'fldvar.inc'
      INCLUDE 'run.inc'
      INCLUDE 'funits.inc'
C
C             loop counters
      INTEGER LC, N
C
C             form feed character = CHAR(12)
C 
      WRITE (UNIT_OUT,1000) CHAR(12), TIME
      CALL OUT_ARRAY(P_g,'P_g')

      WRITE (UNIT_OUT,1050) CHAR(12), TIME
      CALL OUT_ARRAY(P_star,'P_star')

      WRITE (UNIT_OUT,1100) CHAR(12), TIME
      CALL OUT_ARRAY(EP_g,'EP_g')

      WRITE (UNIT_OUT,1200) CHAR(12), TIME
      CALL OUT_ARRAY(RO_g,'RO_g')

      DO 100 LC = 1,MMAX
         WRITE (UNIT_OUT,1400) CHAR(12), LC, TIME
         CALL OUT_ARRAY(ROP_s(1,LC),'ROP_s')
100   CONTINUE

      WRITE (UNIT_OUT,1500) CHAR(12), TIME
      CALL OUT_ARRAY(T_g,'T_g')

      DO 110 LC = 1,MMAX
        WRITE (UNIT_OUT,1600) CHAR(12), LC,  TIME
        CALL OUT_ARRAY(T_s(1,LC),'T_s')
110   CONTINUE

      IF(SPECIES_EQ(0))THEN
        DO 150 N = 1, NMAX(0)
          WRITE(UNIT_OUT,1710) CHAR(12), N, TIME
          CALL OUT_ARRAY(X_g(1,N),'X_g')
150     CONTINUE
      ENDIF

      DO 170 LC = 1, MMAX
        IF(SPECIES_EQ(LC))THEN
          DO 160 N = 1, NMAX(LC)
            WRITE(UNIT_OUT,1720) CHAR(12), LC, N, TIME
            CALL OUT_ARRAY(X_s(1, LC, N), 'X_s')
160       CONTINUE
        ENDIF
170   CONTINUE

      WRITE (UNIT_OUT,1800) CHAR(12), TIME
      CALL OUT_ARRAY(U_g,'U_g')

      WRITE (UNIT_OUT,1900) CHAR(12), TIME
      CALL OUT_ARRAY(V_g,'V_g')

      WRITE (UNIT_OUT,2000) CHAR(12), TIME
      CALL OUT_ARRAY(W_g,'W_g')

      DO 200 LC = 1,MMAX
         WRITE (UNIT_OUT,2100) CHAR(12), LC, TIME
         CALL OUT_ARRAY(U_s(1,LC),'U_s')

         WRITE (UNIT_OUT,2200) CHAR(12), LC, TIME
         CALL OUT_ARRAY(V_s(1,LC),'V_s')

         WRITE (UNIT_OUT,2300) CHAR(12), LC, TIME
         CALL OUT_ARRAY(W_s(1,LC),'W_s')

c         IF(GRANULAR_ENERGY)THEN
           WRITE (UNIT_OUT,2400) CHAR(12), LC, TIME
           CALL OUT_ARRAY(Theta_m(1,LC),'Theta_m')
c         ENDIF
200   CONTINUE
C
C  Print user-defined output
C
      WRITE (UNIT_OUT,'(/1X,1A1)') CHAR(12)
      IF(CALL_USR) CALL USR_WRITE_OUT1
      RETURN
1000  FORMAT(1X,A1,/5X,'--- Gas pressure (P_g) at time ',G12.5,' ---'
     &       ,//)
1050  FORMAT(1X,A1,/5X,'--- Solids pressure (P_star) at time ',G12.5,
     &       ' ---' ,//)
1100  FORMAT(1X,A1,/5X,'--- Void fraction (EP_g) at time ',G12.5,
     &       ' ---',//)
1200  FORMAT(1X,A1,/5X,'--- Gas density (RO_g) at time ',G12.5,
     &       ' ---',//)
1400  FORMAT(1X,A1,/5X,'--- Solids Phase-', I1, ' density x volume',
     &       ' fraction (ROP_s) at time ',G12.5,' ---',//)
1500  FORMAT(1X,A1,/5X,'--- Gas temperature (T_g) at time ',G12.5,
     &       ' ---', //)
1600  FORMAT(1X,A1,/5X,'--- Solids Phase-', I1, ' temperature (T_s)',
     &      ' at time ', G12.5,' ---',//)
1710  FORMAT(1X,A1,/5X,'--- Mass fraction of gas species (X_g) ', I2, 
     &       ' at time ', G12.5,' ---',//)
1720  FORMAT(1X,A1,/5X,'--- Mass fraction of solids-',I1, 
     &  ' species (X_s)',I2, ' at time ', G12.5,' ---',//)
1800  FORMAT(1X,A1,/5X,'--- X-component of gas velocity (U_g) at time ',
     &       G12.5,' ---', //)
1900  FORMAT(1X,A1,/5X,'--- Y-component of gas velocity (V_g) at time ',
     &       G12.5,' ---', //)
2000  FORMAT(1X,A1,/5X,'--- Z-component of gas velocity (W_g) at time ',
     &       G12.5,' ---', //)
2100  FORMAT(1X,A1,/5X,'--- X-component of Solids Phase-',I1,
     &       ' velocity (U_s) at time ', G12.5,' ---', //)
2200  FORMAT(1X,A1,/5X,'--- Y-component of Solids Phase-',I1,
     &       ' velocity (V_s) at time ', G12.5,' ---', //)
2300  FORMAT(1X,A1,/5X,'--- Z-component of Solids Phase-',I1,
     &       ' velocity (W_s) at time ', G12.5,' ---', //)
2400  FORMAT(1X,A1,/5X,'--- Granular temperature of Solids Phase-',I1,
     &       ' velocity (W_s) at time ', G12.5,' ---', //)
      END
