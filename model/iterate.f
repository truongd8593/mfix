!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: ITERATE(IER, NIT)                                      C
!  Purpose: This module controls the iterations for solving equations  C
!           Version 2.0                                                C
!                                                                      C
!  Author: M. Syamlal                                 Date: 12-APR-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
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
      SUBROUTINE ITERATE(IER, NIT) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE toleranc 
      USE run
      USE physprop
      USE geometry
      USE fldvar
      USE output
      USE indices
      USE funits 
      USE time_cpu 
      USE pscor
      USE coeff
      USE leqsol 
      USE visc_g
      USE pgcor
      USE cont
      USE compar               !//d
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
! 
!                      Error index 
      INTEGER          IER 
! 
!                      Number of iterations 
      INTEGER          NIT 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
! 
!                      current cpu time used 
      DOUBLE PRECISION CPU_NOW 
! 
!                      MUSTIT = 0 implies complete convergence. 
      INTEGER          MUSTIT 
! 
!                      Sum of solids densities 
      DOUBLE PRECISION SUM 
! 
!                      Weight of solids in the reactor 
      DOUBLE PRECISION SMASS 
! 
!                      Heat loss from the reactor 
      DOUBLE PRECISION HLOSS 
! 
!                      phase index 
      INTEGER          M 
! 
!                      Normalization factor for gas pressure residual 
      DOUBLE PRECISION NORMg 
! 
!                      Normalization factor for solids pressure residual 
      DOUBLE PRECISION NORMs 
! 
!                      Set normalization factor for gas pressure residual 
      LOGICAL          SETg 
! 
!                      Set normalization factor for solids pressure residual 
      LOGICAL          SETs 
! 
!                      gas pressure residual 
      DOUBLE PRECISION RESg 
! 
!                      solids pressure residual 
      DOUBLE PRECISION RESs 
  
      DOUBLE PRECISION TLEFT 
      CHARACTER*4 TUNIT 
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      DOUBLE PRECISION , EXTERNAL :: VAVG_U_G, VAVG_V_G, VAVG_W_G, VAVG_U_S, &
         VAVG_V_S, VAVG_W_S 
!-----------------------------------------------
!
!
!//AIKEPARDBG
!    write(*,"('(PE ',I2,'): beginning of iterate')") myPE  !//AIKEPARDBG
!    call mfix_exit(myPE)     !//AIKEPARDBG

      NIT = 0 
      RESG = ZERO 
      RESS = ZERO 
!
      IF (NORM_G == UNDEFINED) THEN 
         NORMG = ONE 
         SETG = .FALSE. 
      ELSE 
         NORMG = NORM_G 
         SETG = .TRUE. 
      ENDIF 
!
      IF (NORM_S == UNDEFINED) THEN 
         NORMS = ONE 
         SETS = .FALSE. 
      ELSE 
         NORMS = NORM_S 
         SETS = .TRUE. 
      ENDIF 
!
      LEQ_ADJUST = .FALSE. 

!//AIKEPARDBG
!    write(*,"('(PE ',I2,'): bef init_resid in iterate')") myPE  !//AIKEPARDBG
!    call mfix_exit(myPE)     !//AIKEPARDBG

!
!     Initialize residuals
!
      CALL INIT_RESID (IER) 

!//AIKEPARDBG
!    write(*,"('(PE ',I2,'): aft init_resid in iterate')") myPE  !//AIKEPARDBG
!    call mfix_exit(myPE)     !//AIKEPARDBG

!
!
!     CPU time left
!
!//SP
      IF (FULL_LOG) THEN      !//
         TLEFT = (TSTOP - TIME)*CPUOS 
         CALL GET_TUNIT (TLEFT, TUNIT) 
!
         IF (DT == UNDEFINED) THEN 
            CALL GET_SMASS (SMASS) 
!//SP
	    IF(myPE.eq.PE_IO) THEN
              WRITE (*, '(/A,G10.5, A,F9.3,1X,A)') ' Starting solids mass = ', &
                 SMASS, '    CPU time left = ', TLEFT, TUNIT 
	    ENDIF
         ELSE 
!//SP
            IF(myPE.eq.PE_IO) THEN
              WRITE (*, '(/A,G12.5, A,G12.5, A,F9.3,1X,A)') ' Time = ', TIME, &
                 '  Dt = ', DT, '    CPU time left = ', TLEFT, TUNIT 
	    ENDIF
!
         ENDIF 
      ENDIF 
!
!
!//AIKEPARDBG
!    write(*,"('(PE ',I2,'): bef iteration loop in iterate')") myPE  !//AIKEPARDBG
!    call mfix_exit(myPE)     !//AIKEPARDBG

!     Begin iterations
!
!------------------------------------------------------------------------------
!
!
   50 CONTINUE 
      MUSTIT = 0 
      NIT = NIT + 1 
!
      IF (.NOT.SETG) THEN 
         IF (RESG > SMALL_NUMBER) THEN 
            NORMG = RESG 
            SETG = .TRUE. 
         ENDIF 
      ENDIF 
!
      IF (.NOT.SETS) THEN 
         IF (RESS > SMALL_NUMBER) THEN 
            NORMS = RESS 
            SETS = .TRUE. 
         ENDIF 
      ENDIF 
!
!
!     Call user-defined subroutine to set quantities that need to be updated
!     every iteration
!
      IF (CALL_USR) CALL USR2 
!//SP
!    write(*,"('(PE ',I2,'): aft USR2 in iterate')") myPE  !//AIKEPARDBG
!    call mfix_exit(myPE)     !//AIKEPARDBG
      
!
!
!     Calculate coefficients.  Explicitly set flags for all the quantities
!     that need to be calculated before calling CALC_COEFF.
!
      IF (RO_G0 == UNDEFINED) DENSITY(0) = .TRUE. 
      IF (ANY_SPECIES_EQ) RRATE = .TRUE. 
!
      VISC(0) = RECALC_VISC_G 
      IF (GRANULAR_ENERGY) THEN 
         M = 1 
         IF (MMAX > 0) THEN 
            VISC(1:MMAX) = .TRUE. 
            M = MMAX + 1 
         ENDIF 
      ENDIF 
!
!//AIKEPARDBG
!    write(*,"('(PE ',I2,'): bef calc_coeff in iterate')") myPE  !//AIKEPARDBG
!    call mfix_exit(myPE)     !//AIKEPARDBG

      CALL CALC_COEFF (DENSITY, SIZE, SP_HEAT, VISC, COND, DIFF, RRATE, DRAGCOEF, &
         HEAT_TR, WALL_TR, IER) 
!
!//AIKEPARDBG
!    write(*,"('(PE ',I2,'): aft calc_coeff in iterate')") myPE  !//AIKEPARDBG
!    call mfix_exit(myPE)     !//AIKEPARDBG
      
!
!     Solve strarred velocitiy components
!
!//SP
      CALL SOLVE_VEL_STAR (IER) 
!
!//AIKEPARDBG
!    write(*,"('(PE ',I2,'): aft solve_vel_star in iterate')") myPE  !//AIKEPARDBG
!    call mfix_exit(myPE)     !//AIKEPARDBG

      write(*,*) RO_G0, myPE
!
!     Solve fluid pressure correction equation
!
      IF (RO_G0 /= ZERO) CALL SOLVE_PP_G (NORMG, RESG, IER) 
!//SP
!    write(*,"('(PE ',I2,'): aft solve_pp_g in iterate')") myPE  !//AIKEPARDBG      
!    call mfix_exit(myPE)     !//AIKEPARDBG
!
!
!     Correct pressure and velocities
!
      IF (RO_G0 /= ZERO) CALL CORRECT_0 (IER) 
!//SP
!      write(*,*) 'after CORRECT_0', myPE, MMAX
!//AIKEPARDBG
!    write(*,"('(PE ',I2,'): aft correct_0 in iterate')") myPE  !//AIKEPARDBG
!    call mfix_exit(myPE)     !//AIKEPARDBG
      
!
!     Solve solids volume fraction correction equation for close-packed
!     solids phases
!
      IF (MMAX > 0) THEN 
        CALL CALC_K_CP (K_CP, IER) 
        CALL SOLVE_EPP (NORMS, RESS, IER) 
        CALL CORRECT_1 (IER) 
!
! IER = 0
        CALL CALC_VOL_FR (P_STAR, RO_G, ROP_G, EP_G, ROP_S, IER) 
        IF (IER == 1) THEN 
           MUSTIT = 2                           !indicates divergence 
           IF(DT/=UNDEFINED)GO TO 1000 
        ENDIF 
!
      ENDIF 
!//AIKEPARDBG
!    write(*,"('(PE ',I2,'): aft calc_vol_fr in iterate')") myPE  !//AIKEPARDBG
!    call mfix_exit(myPE)     !//AIKEPARDBG


!
!  Update wall velocities
 
      CALL SET_WALL_BC (IER) 
!
!//AIKEPARDBG
!     write(*,"('(PE ',I2,'): aft set_wall_bc in iterate')") myPE  !//AIKEPARDBG
!     call mfix_exit(myPE)     !//AIKEPARDBG
!     Calculate P_star in cells where solids continuity equation is
!     solved
!
      IF (MMAX > 0) CALL CALC_P_STAR (EP_G, P_STAR, IER) 
!
!
!     Solve energy equations
!
      IF (ENERGY_EQ) CALL SOLVE_ENERGY_EQ (IER) 
!
!     Solve granular energy equation
!
      IER = 0
      IF (GRANULAR_ENERGY) CALL SOLVE_GRANULAR_ENERGY (IER) 
      IF (IER == 1) THEN 
         MUSTIT = 2                              !indicates divergence 
         IF(DT/=UNDEFINED)GO TO 1000 
!
      ENDIF 

!//AIKEPARDBG
!    write(*,"('(PE ',I2,'): bef solve_species_eq in iterate')") myPE  !//AIKEPARDBG
!    call mfix_exit(myPE)     !//AIKEPARDBG
      
!
!     Solve species mass balance equations
!
      CALL SOLVE_SPECIES_EQ (IER) 

!//AIKEPARDBG
!    write(*,"('(PE ',I2,'): aft solve_species_eq in iterate')") myPE  !//AIKEPARDBG
!    call mfix_exit(myPE)     !//AIKEPARDBG
      
!
!    User-defined linear equation solver parameters may be adjusted after
!    the first iteration
!
      IF (.NOT.CYCLIC) LEQ_ADJUST = .TRUE. 
!
!
!
!     Check for convergence
!
      CALL CHECK_CONVERGENCE (NIT, MUSTIT, IER) 
!//SP
      write(*,*) 'after CHECK_CONVERGENCE, etc.,', myPE
!
!      If not converged continue iterations; else exit subroutine.
!
 1000 CONTINUE 
!
!     Display residuals
!
      IF (FULL_LOG) CALL DISPLAY_RESID (NIT, IER) 
!//SP
      write(*,*) 'after DISPLAY_RESID, etc.,', myPE
      
      IF (MUSTIT == 0) THEN 
         IF (DT==UNDEFINED .AND. NIT==1) GO TO 50!Iterations converged 
         IF (MOD(NSTEP,NLOG) == 0) THEN 
            CALL CPU_TIME (CPU_NOW) 
!
            CPUOS = (CPU_NOW - CPU_NLOG)/(TIME - TIME_NLOG) 
            CPU_NLOG = CPU_NOW 
            TIME_NLOG = TIME 
!
            CPU_NOW = CPU_NOW - CPU0 
            CALL GET_SMASS (SMASS) 
            IF (ENERGY_EQ) CALL GET_HLOSS (HLOSS) 
!
!
            IF (ENERGY_EQ) THEN 
               WRITE (UNIT_LOG, 5000) TIME, DT, NIT, SMASS, HLOSS, CPU_NOW 
               IF(FULL_LOG.and.myPE.eq.PE_IO) &          
                       WRITE(*,5000)TIME,DT,NIT,SMASS,HLOSS,CPU_NOW        !//
            ELSE 
               WRITE (UNIT_LOG, 5001) TIME, DT, NIT, SMASS, CPU_NOW 
               IF (FULL_LOG .and. myPE.eq.PE_IO) &
                       WRITE (*, 5001) TIME, DT, NIT, SMASS, CPU_NOW       !//
            ENDIF 
            CALL START_LOG 
            IF (.NOT.FULL_LOG) THEN 
               TLEFT = (TSTOP - TIME)*CPUOS 
               CALL GET_TUNIT (TLEFT, TUNIT) 
               WRITE (UNIT_LOG, '(46X, A, F9.3, 1X, A)') '    CPU time left = '&
                  , TLEFT, TUNIT 
            ENDIF 
!
            IF (CYCLIC_X .OR. CYCLIC_Y .OR. CYCLIC_Z) THEN 
               IF (DO_I) WRITE (UNIT_LOG, 5050) 'U_g = ', VAVG_U_G() 
               IF (DO_J) WRITE (UNIT_LOG, 5050) 'V_g = ', VAVG_V_G() 
               IF (DO_K) WRITE (UNIT_LOG, 5050) 'W_g = ', VAVG_W_G() 
               DO M = 1, MMAX 
                  IF (DO_I) WRITE (UNIT_LOG, 5060) 'U_s(', M, ') = ', VAVG_U_S(&
                     M) 
                  IF (DO_J) WRITE (UNIT_LOG, 5060) 'V_s(', M, ') = ', VAVG_V_S(&
                     M) 
                  IF (DO_K) WRITE (UNIT_LOG, 5060) 'W_s(', M, ') = ', VAVG_W_S(&
                     M) 
               END DO 
            ENDIF 
!
            CALL END_LOG 
         ENDIF 
         IER = 0 
         RETURN  
!                                                ! diverged or
      ELSE IF (MUSTIT==2 .AND. DT/=UNDEFINED) THEN 
         IF (FULL_LOG) THEN 
            CALL START_LOG 
            WRITE (UNIT_LOG, 5200) TIME, DT, NIT 
            CALL END_LOG 

!//? pnicol : NIT below is only for PE_IO ... can other
!//?          processors have different NIT.  Ans if so,
!//?          what do we want to print out (MAX ?)
            if (myPE.eq.PE_IO) WRITE (*, 5200) TIME, DT, NIT   !//
         ENDIF 
         IER = 1 
         RETURN  
      ENDIF 
!
!//SP
      write(*,*) 'after misc., etc.,', myPE
!     call mfix_exit
!
      IF (NIT < MAX_NIT) THEN 
         MUSTIT = 0 
         GO TO 50 
      ENDIF 
!
!------------------------------------------------------------------------------
!     End iterations
!
      CALL GET_SMASS (SMASS) 
!//? pnicol : NIT below is only for PE_IO ... can other
!//?          processors have different NIT.  Ans if so,
!//?          what do we want to print out (MAX ?)
      if (myPE.eq.PE_IO) WRITE (UNIT_OUT, 5100) TIME, DT, NIT, SMASS    !//
      CALL START_LOG 
      WRITE (UNIT_LOG, 5100) TIME, DT, NIT, SMASS 
      CALL END_LOG 
!
      IER = 0 
      RETURN  
 5000 FORMAT(1X,'t=',F10.4,' Dt=',G10.4,' NIT=',I3,' Sm=',G10.5,' Hl=',G12.5,&
         T67,'CPU=',F8.0,' s') 
 5001 FORMAT(1X,'t=',F10.4,' Dt=',G10.4,' NIT=',I3,' Sm=',G10.5,T67,'CPU=',F8.0&
         ,' s') 
 5050 FORMAT(5X,'Average ',A,G12.5) 
 5060 FORMAT(5X,'Average ',A,I2,A,G12.5) 
 5100 FORMAT(1X,'t=',F10.4,' Dt=',G10.4,' NIT>',I3,' Sm= ',G10.5) 
 5200 FORMAT(1X,'t=',F10.4,' Dt=',G10.4,' NIT=',I3,': Run diverged/stalled :-('&
         ) 
      END SUBROUTINE ITERATE 
!
!
      SUBROUTINE GET_TUNIT(TLEFT, TUNIT) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE toleranc 
      USE run 
      USE physprop 
      USE geometry 
      USE fldvar 
      USE output 
      USE indices 
      USE funits 
      USE time_cpu 
      USE pscor 
      USE coeff 
      USE leqsol 
      USE visc_g 
      USE pgcor 
      USE cont 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      DOUBLE PRECISION TLEFT 
      CHARACTER TUNIT*4 
!-----------------------------------------------
!
!
      IF (TLEFT < 3600.) THEN 
         TUNIT = 's' 
      ELSE 
         TLEFT = TLEFT/3600. 
         TUNIT = 'h' 
         IF (TLEFT >= 24.) THEN 
            TLEFT = TLEFT/24. 
            TUNIT = 'days' 
         ENDIF 
      ENDIF 
!
      RETURN  
      END SUBROUTINE GET_TUNIT 
