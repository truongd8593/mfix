!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: TIME_MARCH                                             C
!  Purpose: Controlling module for time marching and finding the       C
!           solution of equations from TIME to TSTOP at intervals of   C
!           DT, updating the b.c.'s, and creating output.              C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-JAN-92  C
!  Reviewer:M. Syamlal, S. Venkatesan, P. Nicoletti,  Date: 29-JAN-92  C
!           W. Rogers                                                  C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: Change subroutine name from SET_GASFLUX to SET_FLUXES      C
!           Add a CALC_THETA call for granular stresses                C
!  Author: M. Syamlal                                 Date: 20-FEB-92  C
!  Reviewer: S. Venkatesan                            Date: 11-DEC-92  C
!  Revision Number: 2                                                  C
!  Purpose: Changes for MFIX 2.0                                       C
!  Author: M. Syamlal                                 Date: 12-APR-96  C
!  Reviewer:                                 Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: RUN_TYPE, TIME, DT, NSTEP, TSTOP, OUT_DT,     C
!                        RES_DT, USR_DT, SPX_DT                        C
!                                                                      C
!  Variables modified: TIME, NSTEP                                     C
!                                                                      C
!  Local variables: OUT_TIME, RES_TIME, USR_TIME, SPX_TIME, USR_TIME,  C
!                   L, FINISH                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE TIME_MARCH 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE run
      USE output
      USE physprop
      USE fldvar
      USE geometry
      USE pgcor
      USE pscor
      USE cont
      USE coeff
      USE tau_g
      USE tau_s
      USE visc_g
      USE visc_s
      USE funits 
      USE compar               !//d
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      DOUBLE PRECISION, PARAMETER :: ONEMEG = 1048576 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
! 
!                      Flag to indicate one pass through iterate for steady 
!                      state conditions. 
      LOGICAL          FINISH 
! 
!                      Time at which standard output is to be written 
      DOUBLE PRECISION OUT_TIME 
! 
!                      Time at which restart file is to be written 
      DOUBLE PRECISION RES_TIME 
! 
!                      Time at which REAL restart file is to be written 
      DOUBLE PRECISION SPX_TIME(N_SPX) 
! 
!                      DIsk space needed for one variable and each 
!                      SPx file 
      DOUBLE PRECISION DISK_ONE, DISK(N_SPX) 
! 
!                      Total DIsk space 
      DOUBLE PRECISION DISK_TOT 
! 
!                      number SPX writes 
      INTEGER          ISPX 
  
      LOGICAL          RES_MSG, SPX_MSG 
! 
!                      Time at which special output is to be written 
      DOUBLE PRECISION USR_TIME (DIMENSION_USR) 
! 
!                      Loop indices 
      INTEGER          L, M 
! 
!                      Error index 
      INTEGER          IER 
  
! 
!                      Number of iterations 
      INTEGER          NIT 
! 
!                      used for activating check_data_30 
      INTEGER          NCHECK, DNCHECK 
! 
!                      dummy logical variable for initializing adjust_dt 
      LOGICAL          dummy 
! 
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      LOGICAL , EXTERNAL :: ADJUST_DT 
!-----------------------------------------------
!
!
      FINISH = .FALSE. 
      NCHECK = NSTEP 
      DNCHECK = 1 
!
!  Initialize times for writing outputs
!
      OUT_TIME = TIME 
      LOG_HEADER = .TRUE. 
      RES_MSG = .TRUE. 
      SPX_MSG = .TRUE. 
!
!  Initialize disk space calculations
!
      DISK_TOT = ZERO 
      DISK_ONE = 4.*IJKMAX2/ONEMEG 
      DISK(1) = DISK_ONE 
      DISK(2) = 2.*DISK_ONE 
      DISK(3) = 3.*DISK_ONE 
      DISK(4) = MMAX*3.*DISK_ONE 
      DISK(5) = MMAX*DISK_ONE 
      DISK(6) = 3.*DISK_ONE 
      DISK(7) = NMAX(0)*DISK_ONE 
      M = 1 
      IF (MMAX > 0) THEN 
         DISK(7) = DISK(7) + SUM(NMAX(1:MMAX)*DISK_ONE) 
         M = MMAX + 1 
      ENDIF 
      DISK(8) = MMAX*DISK_ONE 
!
!
!
      if (myPE.eq.PE_IO) then        !//
       IF (RUN_TYPE == 'NEW') THEN 
         RES_TIME = TIME 
         SPX_TIME(:N_SPX) = TIME 
         L = N_SPX + 1 
       ELSE 
         IF (DT /= UNDEFINED) THEN 
            RES_TIME = (INT((TIME + 0.1*DT)/RES_DT) + 1)*RES_DT 
            SPX_TIME(:N_SPX) = (INT((TIME + 0.1*DT)/SPX_DT(:N_SPX))+1)*SPX_DT(:&
               N_SPX) 
            L = N_SPX + 1 
         ENDIF 
       ENDIF 
!
       DO L = 1, DIMENSION_USR 
         USR_TIME(L) = UNDEFINED 
         IF (USR_DT(L) /= UNDEFINED) THEN 
            IF (RUN_TYPE == 'NEW') THEN 
               USR_TIME(L) = TIME 
            ELSE 
               USR_TIME(L) = (INT((TIME + 0.1*DT)/USR_DT(L))+1)*USR_DT(L) 
            ENDIF 
         ENDIF 
       END DO
    end if             !//
!
!   Parse residual strings
!
      CALL PARSE_RESID_STRING (IER) 
!
!  Initialization for the linear equation solver igcg
!
      CALL IGCG_INIT 
!
!  Call user-defined subroutine to set constants, check data, etc.
!
      IF (CALL_USR) CALL USR0 
!
!  Calculate all the coefficients once before entering the time loop
!
      RRATE = .TRUE. 
      WALL_TR = .TRUE. 
      M = 0 
      IF (MMAX + 1 > 0) THEN 
         DENSITY(:MMAX) = .TRUE. 
         SIZE(:MMAX) = .TRUE. 
         IF (ENERGY_EQ) THEN 
            SP_HEAT(:MMAX) = .TRUE. 
            COND(:MMAX) = .TRUE. 
         ELSE 
            SP_HEAT(:MMAX) = .FALSE. 
            COND(:MMAX) = .FALSE. 
         ENDIF 
         VISC(:MMAX) = .TRUE. 
         DIFF(:MMAX) = .TRUE. 
         L = 0 
         DRAGCOEF(:MMAX,:MMAX) = .TRUE. 
         IF (ENERGY_EQ) THEN 
            HEAT_TR(:MMAX,:MMAX) = .TRUE. 
         ELSE 
            HEAT_TR(:MMAX,:MMAX) = .FALSE. 
         ENDIF 
         L = MMAX + 1 
         M = MMAX + 1 
      ENDIF 
      IF (RO_G0 /= UNDEFINED) DENSITY(0) = .FALSE. 
      IF (MU_S0 /= UNDEFINED) VISC(1) = .FALSE. 
      CALL CALC_COEFF (DENSITY, SIZE, SP_HEAT, VISC, COND, DIFF, RRATE, DRAGCOEF, &
         HEAT_TR, WALL_TR, IER) 
!
!  Remove undefined values at wall cells for scalars
!
      CALL UNDEF_2_0 (ROP_G, IER) 
      DO M = 1, MMAX 
         CALL UNDEF_2_0 (ROP_S(1,M), IER) 
      END DO 

!
!  Initialize d's and e's to zero   !//? pnicol : ??????
!
      DO M = 0, MMAX 
         CALL ZERO_ARRAY (D_E(1,M), IJKMAX2, IER) 
         CALL ZERO_ARRAY (D_N(1,M), IJKMAX2, IER) 
         CALL ZERO_ARRAY (D_T(1,M), IJKMAX2, IER) 
      END DO 
      CALL ZERO_ARRAY (E_E, IJKMAX2, IER) 
      CALL ZERO_ARRAY (E_N, IJKMAX2, IER) 
      CALL ZERO_ARRAY (E_T, IJKMAX2, IER) 
!
!   Initialize adjust_ur
!
      dummy = ADJUST_DT(100, 0)

!
!  The TIME loop begins here.............................................
!
  100 CONTINUE 
      IF (CALL_USR) CALL USR1 
!
!  Remove solids from cells containing very small quantities of solids
!
      CALL ADJUST_EPS 
!
!  Mark the phase whose continuity will be used for forming Pp_g and Pp_s eqs.
!
      CALL MARK_PHASE_4_COR (PHASE_4_P_G, PHASE_4_P_S, DO_CONT, MCP, DO_P_S, &
         SWITCH_4_P_G, SWITCH_4_P_S, IER) 
!
!  Set wall boundary conditions and transient flow b.c.'s
!
      CALL SET_BC1 
!
! Write standard output, if needed
!
      if (myPE.eq.PE_IO) then       !//
         IF (OUT_DT /= UNDEFINED) THEN 
            IF (DT == UNDEFINED) THEN 
               CALL WRITE_OUT1 
            ELSE IF (TIME + 0.1*DT>=OUT_TIME .OR. TIME+0.1*DT>=TSTOP) THEN 
               OUT_TIME = (INT((TIME + 0.1*DT)/OUT_DT) + 1)*OUT_DT 
               CALL WRITE_OUT1 
            ENDIF 
         ENDIF 
      end if                          !//
!
! Write restart file, if needed
!
      CALL START_LOG 
      IF (DT == UNDEFINED) THEN 
         IF (FINISH) THEN 
            CALL WRITE_RES1 
            RES_MSG = .FALSE. 
            WRITE (UNIT_LOG, '('' t='',F10.4, ''  Wrote RES;'')', ADVANCE='NO') TIME 
            IF (FULL_LOG .and. myPE.eq.PE_IO) WRITE (*, 1000,  ADVANCE='NO') TIME  !//
         ENDIF 
      ELSE IF (TIME + 0.1*DT>=RES_TIME .OR. TIME+0.1*DT>=TSTOP) THEN 
         RES_TIME = (INT((TIME + 0.1*DT)/RES_DT) + 1)*RES_DT 
         CALL WRITE_RES1 
         RES_MSG = .FALSE. 
         WRITE (UNIT_LOG, 1000,  ADVANCE='NO') TIME 
         IF (FULL_LOG .and. myPE.eq.PE_IO) WRITE (*, 1000,  ADVANCE='NO') TIME   !//
      ENDIF 
!
! Write SPx files, if needed
!
      if (myPE.eq.PE_IO) then         !//
         ISPX = 0 
         DO L = 1, N_SPX 
            IF (DT == UNDEFINED) THEN 
               IF (FINISH) THEN 
                  CALL WRITE_SPX1 (L) 
                  DISK_TOT = DISK_TOT + DISK(L) 
                  ISPX = ISPX + 1 
!
                  IF (SPX_MSG) THEN 
                     IF (RES_MSG) THEN 
                        WRITE (UNIT_LOG, 1001,  ADVANCE='NO') TIME 
                        IF (FULL_LOG) WRITE (*, 1001,  ADVANCE='NO') TIME 
                     ELSE 
                        WRITE (UNIT_LOG, 1002,  ADVANCE='NO') 
                        IF (FULL_LOG) WRITE (*, 1002,  ADVANCE='NO') 
                     ENDIF 
                     SPX_MSG = .FALSE. 
                  ENDIF 
                  WRITE (UNIT_LOG, 1010,  ADVANCE='NO') L 
                  IF (FULL_LOG) WRITE (*, 1010,  ADVANCE='NO') L 
               ENDIF 
            ELSE IF (TIME + 0.1*DT>=SPX_TIME(L) .OR. TIME+0.1*DT>=TSTOP) THEN 
               SPX_TIME(L) = (INT((TIME + 0.1*DT)/SPX_DT(L))+1)*SPX_DT(L) 
               CALL WRITE_SPX1 (L) 
               DISK_TOT = DISK_TOT + DISK(L) 
               ISPX = ISPX + 1 
!
               IF (SPX_MSG) THEN 
                  IF (RES_MSG) THEN 
                     WRITE (UNIT_LOG, 1001,  ADVANCE='NO') TIME 
                     IF (FULL_LOG) WRITE (*, 1001,  ADVANCE='NO') TIME 
                  ELSE 
                     WRITE (UNIT_LOG, 1002,  ADVANCE='NO') 
                     IF (FULL_LOG) WRITE (*, 1002,  ADVANCE='NO') 
                  ENDIF 
                  SPX_MSG = .FALSE. 
               ENDIF 
               WRITE (UNIT_LOG, 1010,  ADVANCE='NO') L 
               IF (FULL_LOG) WRITE (*, 1010,  ADVANCE='NO') L 
            ENDIF 
         END DO 
      end if                         !//
      IF (.NOT.SPX_MSG) THEN 
         DO L = 1, N_SPX - ISPX 
            WRITE (UNIT_LOG, '(A,$)') '   ' 
            IF (FULL_LOG .and. myPE.eq.PE_IO) WRITE (*, '(A,$)') '   '  !//
         END DO 
         WRITE (UNIT_LOG, 1015) DISK_TOT 
         IF (FULL_LOG.and.myPE.eq.PE_IO) WRITE (*, 1015) DISK_TOT    !//
      ELSE IF (.NOT.RES_MSG) THEN 
         WRITE (UNIT_LOG, *) 
         IF (FULL_LOG .and. myPE.eq.PE_IO) WRITE (*, *)    !//
      ENDIF 
!
      RES_MSG = .TRUE. 
      SPX_MSG = .TRUE. 
      CALL END_LOG 
!
!  Write special output, if needed
!
      if (myPE.eq.PE_IO) then           !//
        DO L = 1, DIMENSION_USR 
         IF (DT == UNDEFINED) THEN 
            IF (FINISH) CALL WRITE_USR1 (L) 
         ELSE IF (USR_TIME(L)/=UNDEFINED .AND. TIME+0.1*DT>=USR_TIME(L)) THEN 
            USR_TIME(L) = (INT((TIME + 0.1*DT)/USR_DT(L))+1)*USR_DT(L) 
            CALL WRITE_USR1 (L) 
         ENDIF 
        END DO 
        IF (DT == UNDEFINED) THEN 
         IF (FINISH) THEN 
            RETURN  
         ELSE 
            FINISH = .TRUE. 
         ENDIF 
        ELSE IF (TIME + 0.1*DT >= TSTOP) THEN 
         RETURN  
        ENDIF 
      end if                !//
!
!  Update previous-time-step values of field variables
!
      CALL UPDATE_OLD 
!
!     Calculate coefficients.  Explicitly set flags for all the quantities
!     that need to be calculated before calling CALC_COEFF.
!
      M = 0 
      IF (MMAX + 1 > 0) THEN 
         IF (ENERGY_EQ) THEN 
            SP_HEAT(:MMAX) = .TRUE. 
            COND(:MMAX) = .TRUE. 
            DIFF(:MMAX) = .TRUE. 
         ENDIF 
         VISC(:MMAX) = .TRUE. 
         L = 0 
         DRAGCOEF(:MMAX,:MMAX) = .TRUE. 
         IF (ENERGY_EQ) HEAT_TR(:MMAX,:MMAX) = .TRUE. 
         L = MMAX + 1 
         M = MMAX + 1 
      ENDIF 
      IF (RO_G0 /= UNDEFINED) DENSITY(0) = .FALSE. 
      IF (MU_S0 /= UNDEFINED) VISC(1) = .FALSE. 
      CALL CALC_COEFF (DENSITY, SIZE, SP_HEAT, VISC, COND, DIFF, RRATE, DRAGCOEF, &
         HEAT_TR, WALL_TR, IER) 
!
!     Calculate the cross terms of the stress tensor
!
      CALL CALC_TRD_G (TRD_G, IER) 
      CALL CALC_TRD_S (TRD_S, IER) 
      CALL CALC_TAU_U_G (TAU_U_G, IER) 
      CALL CALC_TAU_V_G (TAU_V_G, IER) 
      CALL CALC_TAU_W_G (TAU_W_G, IER) 
      CALL CALC_TAU_U_S (TAU_U_S, IER) 
      CALL CALC_TAU_V_S (TAU_V_S, IER) 
      CALL CALC_TAU_W_S (TAU_W_S, IER) 
!
!  Check rates and sums of mass fractions every NLOG time steps
!
      IF (NSTEP == NCHECK) THEN 
         IF (DNCHECK < 256) DNCHECK = DNCHECK*2 
         NCHECK = NCHECK + DNCHECK 
         CALL CHECK_DATA_30 
      ENDIF 
!
!  Advance the solution in time by iteratively solving the equations
!
      CALL ITERATE (IER, NIT) 
!
!
!  Adjust time step and reiterate if necessary
!
      DO WHILE(ADJUST_DT(IER,NIT)) 
         CALL ITERATE (IER, NIT) 
      END DO 
!
!  Advance the time step and continue
!
      IF (DT /= UNDEFINED) THEN 
         TIME = TIME + DT 
         NSTEP = NSTEP + 1 
      ENDIF 
      CALL FLUSH (6) 
      GO TO 100 
!
!  The TIME loop ends here....................................................
!
!1000  FORMAT(' t=',F10.4, '  Wrote RES;',$)
!1001  FORMAT(' t=',F10.4, '  Wrote      SPx:',$)
!1002  FORMAT(' SPx:',$)
!1010  FORMAT(I2,',',$)
!1015  FORMAT(14X, 'Disk=', F7.2,' Mb')
 1000 FORMAT(' t=',F10.4,'  Wrote RES;') 
 1001 FORMAT(' t=',F10.4,'  Wrote      SPx:') 
 1002 FORMAT(' SPx:') 
 1010 FORMAT(I2,',') 
 1015 FORMAT(14X,'Disk=',F7.2,' Mb') 
      END SUBROUTINE TIME_MARCH 
!
