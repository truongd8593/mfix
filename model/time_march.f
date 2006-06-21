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
!  Reviewer:                                          Date:            C
!                                                                      C
!  Reavision Number: 2                                                 C
!  Purpose: To call DES related routines when doing DES                C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!                                                                      C
!  Reavision Number: 3                                                 C
!  Purpose: To call ISAT to calculate chemical rxns                    C
!  Author: Nan Xie                                    Date: 02-Aug-04  C 
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
!...  Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...  Switches: -xf
!-----------------------------------------------
!     M o d u l e s 
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
      USE vshear
      USE scalars
      USE toleranc
      USE drag
      USE rxns
      USE compar     
      USE time_cpu  
      USE discretelement  
      USE mchem
      IMPLICIT NONE
!-----------------------------------------------
!     G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!     L o c a l   P a r a m e t e r s
!-----------------------------------------------
      DOUBLE PRECISION, PARAMETER :: ONEMEG = 1048576
  
!-----------------------------------------------
!     L o c a l   V a r i a b l e s
!-----------------------------------------------
!     
!     Flag to indicate one pass through iterate for steady 
!     state conditions. 
      LOGICAL          FINISH 
!     
!     Time at which standard output is to be written 
      DOUBLE PRECISION OUT_TIME 
!     
!     Time at which restart file is to be written 
      DOUBLE PRECISION RES_TIME 
!     
!     Time at which REAL restart file is to be written 
      DOUBLE PRECISION SPX_TIME(N_SPX) 
!     
!     DIsk space needed for one variable and each 
!     SPx file 
      DOUBLE PRECISION DISK_ONE, DISK(N_SPX) 
!     
!     Total DIsk space 
      DOUBLE PRECISION DISK_TOT 
!     
!     number SPX writes 
      INTEGER          ISPX 
      
      LOGICAL          RES_MSG, SPX_MSG 
!     
!     Time at which special output is to be written 
      DOUBLE PRECISION USR_TIME (DIMENSION_USR) 
!     
!     Loop indices 
      INTEGER          L, M 
!     
!     Error index 
      INTEGER          IER 
      
!     
!     Number of iterations 
      INTEGER          NIT 
!     
!     used for activating check_data_30 
      INTEGER          NCHECK, DNCHECK,ijk 
!     
!     dummy logical variable for initializing adjust_dt 
      LOGICAL          dummy 

      CHARACTER        EXT_END*35 
!     
!     use function vavg_v_g to catch NaN's 
      DOUBLE PRECISION VAVG_U_G, VAVG_V_G, VAVG_W_G, X_vavg     
!
!     use function MAX_VEL_INLET to compute max. velocity at inlet
      DOUBLE PRECISION MAX_VEL_INLET
!     
!-----------------------------------------------
!     E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      LOGICAL , EXTERNAL :: ADJUST_DT 
!-----------------------------------------------
!     
      IF(AUTOMATIC_RESTART) RETURN
!     
      EXT_END = '123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'
!     
      FINISH  = .FALSE. 
      NCHECK  = NSTEP 
      DNCHECK = 1 
      CPU_IO  = ZERO 
!     
!     Initialize times for writing outputs
!     
      OUT_TIME = TIME 
      LOG_HEADER = .TRUE. 
      RES_MSG = .TRUE. 
      SPX_MSG = .TRUE. 
!     
!     Initialize disk space calculations
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
      DISK(9) = NScalar*DISK_ONE
      DISK(10) = nRR*DISK_ONE
!     
!     
!     

      IF (RUN_TYPE == 'NEW') THEN 
         RES_TIME = TIME 
         SPX_TIME(:N_SPX) = TIME 
         L = N_SPX + 1 
      ELSE 
         IF (DT /= UNDEFINED) THEN 
            RES_TIME = (INT((TIME + 0.1d0*DT)/RES_DT) + 1)*RES_DT 
            SPX_TIME(:N_SPX) = (INT((TIME + 0.1d0*DT)/SPX_DT(:N_SPX))+1)*SPX_DT(:&
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
               USR_TIME(L) = (INT((TIME + 0.1d0*DT)/USR_DT(L))+1)*USR_DT(L) 
            ENDIF 
         ENDIF 
      END DO
      
!     
!     Parse residual strings
!     
      CALL PARSE_RESID_STRING (IER) 

!     
!     Call user-defined subroutine to set constants, check data, etc.
!     
      IF (CALL_USR) CALL USR0 
!     
!     Calculate all the coefficients once before entering the time loop
!     
!

 
!  CHEM & ISAT begin (Nan Xie)
      IF (CALL_DI) THEN 
         CALL MCHEM_INIT
         CALL MCHEM_ODEPACK_INIT
      END IF

      IF (CALL_ISAT) THEN 
         CALL MCHEM_INIT
         CALL MCHEM_ODEPACK_INIT
         CALL MISAT_TABLE_INIT
      END IF
!  CHEM & ISAT end (Nan Xie)

      CALL RRATES_INIT(IER)

      DO M=1, MMAX 
         CALL ZERO_ARRAY (F_gs(1,M), IER)
      END DO

!
!  Calculate all the coefficients once before entering the time loop
!

      CALL CALC_COEFF_ALL (0, IER) 

!     
!     Remove undefined values at wall cells for scalars
!     
      CALL UNDEF_2_0 (ROP_G, IER) 
      DO M = 1, MMAX 
         CALL UNDEF_2_0 (ROP_S(1,M), IER) 
      END DO 


!     
!     Initialize d's and e's to zero   
!     
      DO M = 0, MMAX 
         CALL ZERO_ARRAY (D_E(1,M), IER) 	 
         CALL ZERO_ARRAY (D_N(1,M), IER) 
         CALL ZERO_ARRAY (D_T(1,M), IER) 
      END DO 
      CALL ZERO_ARRAY (E_E, IER) 
      CALL ZERO_ARRAY (E_N, IER) 
      CALL ZERO_ARRAY (E_T, IER) 

!     
!     Initialize adjust_ur
      
      dummy = ADJUST_DT(100, 0)

!     calculate shear velocities if periodic shear BCs are used
      IF (SHEAR) THEN
         call CAL_D(V_sh)
      END IF

!     Initialize check_mass_balance.  This routine is not active by default.  Specify a
!     reporting interval (hard-wired in the routine) to activate the routine.
      Call check_mass_balance (0)

!     DES 
!     Jay Boyalakuntla
      IF(DISCRETE_ELEMENT) THEN
         CALL MAKE_ARRAYS_DES(PARTICLES)
      END IF
!     DES end 
!     
!     The TIME loop begins here.............................................
!     
!     
 100  CONTINUE
      
      IF (CALL_USR) CALL USR1 

!     
!     Remove solids from cells containing very small quantities of solids
!     
      CALL ADJUST_EPS 

!     
!     Mark the phase whose continuity will be used for forming Pp_g and Pp_s eqs.
!     
      CALL MARK_PHASE_4_COR (PHASE_4_P_G, PHASE_4_P_S, DO_CONT, MCP, DO_P_S, &
      SWITCH_4_P_G, SWITCH_4_P_S, IER) 

!     
!     Set wall boundary conditions and transient flow b.c.'s
!     
      CALL SET_BC1 

      CALL CPU_TIME(CPU0_IO)
!     
!     Write standard output, if needed
!     
      IF (OUT_DT /= UNDEFINED) THEN 
         IF (DT == UNDEFINED) THEN 
            CALL WRITE_OUT1 
         ELSE IF (TIME + 0.1d0*DT>=OUT_TIME .OR. TIME+0.1d0*DT>=TSTOP) THEN 
            OUT_TIME = (INT((TIME + 0.1d0*DT)/OUT_DT) + 1)*OUT_DT 
            CALL WRITE_OUT1 
         ENDIF 
      ENDIF 


!     
!     Write restart file, if needed
!     
      CALL START_LOG 
      IF (DT == UNDEFINED) THEN 
         IF (FINISH) THEN 
            CALL WRITE_RES1 
            RES_MSG = .FALSE. 
            IF(DMP_LOG)WRITE (UNIT_LOG, '(" t=",F10.4, "  Wrote RES;")', ADVANCE='NO') TIME 
            IF (FULL_LOG .and. myPE.eq.PE_IO) THEN
               WRITE (*, 1000,  ADVANCE="NO") TIME 
            ENDIF
         ENDIF 
      ELSE IF (TIME + 0.1d0*DT>=RES_TIME .OR. TIME+0.1d0*DT>=TSTOP) THEN 
         RES_TIME = (INT((TIME + 0.1d0*DT)/RES_DT) + 1)*RES_DT 
         CALL WRITE_RES1 
         RES_MSG = .FALSE. 
         IF(DMP_LOG)WRITE (UNIT_LOG, 1000,  ADVANCE='NO') TIME 
         IF (FULL_LOG .and. myPE.eq.PE_IO) WRITE (*, 1000,  ADVANCE='NO') TIME 
      ENDIF 
!     

!     
!     Write SPx files, if needed
!     
      ISPX = 0 
      DO L = 1, N_SPX 
         IF (DT == UNDEFINED) THEN 
            IF (FINISH) THEN 
               CALL WRITE_SPX1 (L, 0) 
               DISK_TOT = DISK_TOT + DISK(L) 
               ISPX = ISPX + 1 
!     
               IF (SPX_MSG) THEN 
                  IF (RES_MSG) THEN 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1001,  ADVANCE='NO') TIME 
                     IF (FULL_LOG .and. myPE.eq.PE_IO) WRITE (*, 1001,  ADVANCE='NO') TIME 
                  ELSE 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1002,  ADVANCE='NO') 
                     IF (FULL_LOG .and. myPE.eq.PE_IO) WRITE (*, 1002,  ADVANCE='NO') 
                  ENDIF 
                  SPX_MSG = .FALSE. 
               ENDIF 
               IF(DMP_LOG)WRITE (UNIT_LOG, 1011,  ADVANCE='NO') EXT_END(L:L)
               IF (FULL_LOG .and. myPE.eq.PE_IO) WRITE (*, 1011,  ADVANCE='NO') EXT_END(L:L)
            ENDIF 
         ELSE IF (TIME + 0.1d0*DT>=SPX_TIME(L) .OR. TIME+0.1d0*DT>=TSTOP) THEN 
            SPX_TIME(L) = (INT((TIME + 0.1d0*DT)/SPX_DT(L))+1)*SPX_DT(L) 
            CALL WRITE_SPX1 (L, 0) 
            DISK_TOT = DISK_TOT + DISK(L) 
            ISPX = ISPX + 1 
!     
            IF (SPX_MSG) THEN 
               IF (RES_MSG) THEN 
                  IF(DMP_LOG)WRITE (UNIT_LOG, 1001,  ADVANCE='NO') TIME 
                  IF (FULL_LOG .and. myPE.eq.PE_IO) WRITE (*, 1001,  ADVANCE='NO') TIME 
               ELSE 
                  IF(DMP_LOG)WRITE (UNIT_LOG, 1002,  ADVANCE='NO') 
                  IF (FULL_LOG .and. myPE.eq.PE_IO) WRITE (*, 1002,  ADVANCE='NO') 
               ENDIF 
               SPX_MSG = .FALSE. 
            ENDIF 
            IF(DMP_LOG)WRITE (UNIT_LOG, 1011,  ADVANCE='NO') EXT_END(L:L)
            IF (FULL_LOG .and. myPE.eq.PE_IO) WRITE (*, 1011,  ADVANCE='NO') EXT_END(L:L)
         ENDIF 
      END DO 

      IF (.NOT.SPX_MSG) THEN 
         DO L = 1, N_SPX - ISPX 
            IF(DMP_LOG)WRITE (UNIT_LOG, '(A,$)') '   ' 
            IF (FULL_LOG .and. myPE.eq.PE_IO) WRITE (*, '(A,$)') '   ' !//
         END DO 
         IF(DMP_LOG)WRITE (UNIT_LOG, 1015) DISK_TOT 
         IF (FULL_LOG.and.myPE.eq.PE_IO) WRITE (*, 1015) DISK_TOT !//
      ELSE IF (.NOT.RES_MSG) THEN 
         IF(DMP_LOG)WRITE (UNIT_LOG, *) 
         IF (FULL_LOG .and. myPE.eq.PE_IO) WRITE (*, *) !//
      ENDIF 

!     
      RES_MSG = .TRUE. 
      SPX_MSG = .TRUE. 
      CALL END_LOG 
!     
      CALL CPU_TIME(CPU1_IO)
      CPU_IO = CPU_IO + (CPU1_IO-CPU0_IO)
!     
!     Write special output, if needed
!     
      DO L = 1, DIMENSION_USR 
         IF (DT == UNDEFINED) THEN 
!//   
            IF (FINISH.and.myPE.eq.PE_IO) CALL WRITE_USR1 (L) 
         ELSE IF (USR_TIME(L)/=UNDEFINED .AND. TIME+0.1d0*DT>=USR_TIME(L)) THEN 
            USR_TIME(L) = (INT((TIME + 0.1d0*DT)/USR_DT(L))+1)*USR_DT(L) 
!//   
            if (myPE.eq.PE_IO) CALL WRITE_USR1 (L) 
         ENDIF 
      END DO 
      IF (DT == UNDEFINED) THEN 
         IF (FINISH) THEN 
            RETURN  
         ELSE 
            FINISH = .TRUE. 
         ENDIF 
      ELSE IF (TIME + 0.1d0*DT >= TSTOP) THEN 
         RETURN  
      ENDIF 

!     
!     Update previous-time-step values of field variables
!     
      CALL UPDATE_OLD 

!     
!     Calculate coefficients

      CALL CALC_COEFF_ALL (0, IER) 
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

!     Check rates and sums of mass fractions every NLOG time steps
!     
      IF (NSTEP == NCHECK) THEN 
         IF (DNCHECK < 256) DNCHECK = DNCHECK*2 
         NCHECK = NCHECK + DNCHECK 
	 
!     Upate the reaction rates for checking
         IF (ANY_SPECIES_EQ) RRATE = .TRUE. 
         CALL CALC_RRATE(RRATE)

         CALL CHECK_DATA_30 
      ENDIF 

!     AE TIME 041601 Double the timestep for 2nd order accurate time implementation
!     IF ((CN_ON.AND.NSTEP>1)) THEN
      IF ((CN_ON.AND.NSTEP>1.AND.RUN_TYPE == 'NEW') .OR. & 
      (CN_ON.AND.RUN_TYPE /= 'NEW' .AND. NSTEP >= (NSTEPRST+1))) THEN
         DT = 0.5d0*DT
         ODT = ODT * 2.0d0
      ENDIF
!      
! Check for maximum velocity at inlet to avoid convergence problems 
!
      MAX_INLET_VEL = 100.0d0*MAX_VEL_INLET ()
!
!     if no inlet velocity is specified, use an upper limit defined in toleranc_mod.f
      IF(MAX_INLET_VEL == ZERO) THEN
	MAX_INLET_VEL = MAX_ALLOWED_VEL
	IF (UNITS == 'SI') MAX_INLET_VEL = 1D-2 * MAX_ALLOWED_VEL
      ENDIF

!     Scale the value using a user defined scale factor      
      MAX_INLET_VEL = MAX_INLET_VEL * MAX_INLET_VEL_FAC
      
!     
!     Advance the solution in time by iteratively solving the equations 
!     
      CALL ITERATE (IER, NIT)
      IF(AUTOMATIC_RESTART) RETURN
      
! Just to Check for NaN's, Uncomment the following lines and also lines  
! of code in  VAVG_U_G, VAVG_V_G, VAVG_W_G to use.
!      X_vavg = VAVG_U_G ()
!      X_vavg = VAVG_V_G ()
!      X_vavg = VAVG_W_G ()
!      IF(AUTOMATIC_RESTART) RETURN

!     
!     Adjust time step and reiterate if necessary
!     
         DO WHILE (ADJUST_DT(IER,NIT))
            CALL ITERATE (IER, NIT) 
         END DO
         IF(DT.LT.DT_MIN) THEN
            IF(TIME.LE.RES_DT) THEN
               IF (AUTO_RESTART .AND. DMP_LOG)WRITE(UNIT_LOG,*) &
	                       'Automatic restart not possible as Total Time < RES_DT'
               CALL MFIX_EXIT(MyPE)
            ENDIF
            IF(AUTO_RESTART) AUTOMATIC_RESTART = .TRUE.
            RETURN
         ENDIF

      
!     Check over mass and elemental balances.  This routine is not active by default.
!     Edit the routine and specify a reporting interval to activate it.
      Call check_mass_balance (1)

!     DES begin
!     Jay Boyalakuntla 
      IF (DISCRETE_ELEMENT) CALL DES_TIME_MARCH
!     DES end
!     
!     
!     Advance the time step and continue
!     
! CHEM & ISAT begin (Nan Xie)
      IF (CALL_ISAT .OR. CALL_DI) THEN 
         CALL MCHEM_TIME_MARCH
      END IF
! CHEM & ISAT end (Nan Xie)
!
!
      IF (CALL_DQMOM) CALL USR_DQMOM
!
!  Advance the time step and continue
!
!     AE TIME 041601 Double the timestep for 2nd order accurate time implementation
!     IF (CN_ON.AND.NSTEP>1) then 
      IF ((CN_ON.AND.NSTEP>1.AND.RUN_TYPE == 'NEW') .OR. & 
      (CN_ON.AND.RUN_TYPE /= 'NEW' .AND. NSTEP >= (NSTEPRST+1))) THEN
         DT = 2.d0*DT      
         ODT = ODT * 0.5d0
      ENDIF      

!     AE TIME 043001 Perform the explicit extrapolation for CN implementation
!     IF (CN_ON.AND.NSTEP>1) then
      IF ((CN_ON.AND.NSTEP>1.AND.RUN_TYPE == 'NEW') .OR. & 
      (CN_ON.AND.RUN_TYPE /= 'NEW' .AND. NSTEP >= (NSTEPRST+1))) THEN     
         CALL CN_EXTRAPOL
      ENDIF

      IF (DT /= UNDEFINED) THEN 
         TIME = TIME + DT 
         NSTEP = NSTEP + 1 
      ENDIF 

!     AIKEDEBUG 081101
!     write (*,"('Compute the Courant number')") 
!     call get_stats(IER)
      
      CALL FLUSH (6) 
      GO TO 100 
!     
!     The TIME loop ends here....................................................
!     
!     1000  FORMAT(' t=',F10.4, '  Wrote RES;',$)
!     1001  FORMAT(' t=',F10.4, '  Wrote      SPx:',$)
!     1002  FORMAT(' SPx:',$)
!     1010  FORMAT(I2,',',$)
!     1015  FORMAT(14X, 'Disk=', F7.2,' Mb')
 1000 FORMAT(' t=',F10.4,'  Wrote RES;') 
 1001 FORMAT(' t=',F10.4,'  Wrote      SPx:') 
 1002 FORMAT(' SPx:') 
 1010 FORMAT(I2,',') 
 1011 FORMAT(A2,',') 
 1015 FORMAT(14X,'Disk=',F7.2,' Mb') 
      END SUBROUTINE TIME_MARCH 
!     
!//   Comments on the modifications for DMP version implementation      
!//   001 Include header file and common declarations for parallelization
!//   Additional I/O checks done by root processor
