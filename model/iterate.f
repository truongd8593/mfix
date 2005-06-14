!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: ITERATE(IER, NIT)                                      C
!  Purpose: This module controls the iterations for solving equations  C
!           Version 2.0                                                C
!                                                                      C
!  Author: M. Syamlal                                 Date: 12-APR-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: To incorporate the "Kinetic" flag so taht the solids       C
!  calculations are not done using kinetic when doing DES              C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
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
      USE scalars
      USE compar   
      USE mpi_utility 
      USE discretelement

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

      DOUBLE PRECISION errorpercent(0:MMAX)
       
      LOGICAL          ABORT_IER
      CHARACTER*4 TUNIT 
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      DOUBLE PRECISION , EXTERNAL :: VAVG_U_G, VAVG_V_G, VAVG_W_G, VAVG_U_S, &
         VAVG_V_S, VAVG_W_S 
	 external cpu_time   !use the subroutine from machine.f
!-----------------------------------------------
!
!
      DT_prev = DT
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

!
!     Initialize residuals
!
      CALL INIT_RESID (IER) 

!	Initialize the routine for holding gas mass flux constant with cyclic bc
      IF(CYCLIC) CALL GoalSeekMassFlux(0, 0, .false.)
!
!
!     CPU time left
!
      IF (FULL_LOG) THEN      !//
         TLEFT = (TSTOP - TIME)*CPUOS 
         CALL GET_TUNIT (TLEFT, TUNIT) 
!
         IF (DT == UNDEFINED) THEN 
            CALL GET_SMASS (SMASS) 

	    IF(myPE.eq.PE_IO) THEN
              WRITE (*, '(/A,G10.5, A,F9.3,1X,A)') ' Starting solids mass = ', &
                 SMASS, '    CPU time left = ', TLEFT, TUNIT 
	    ENDIF
         ELSE 
            IF(myPE.eq.PE_IO) THEN
!AE TIME 050801 add CN_ON check to print original timestep size
!             IF (CN_ON) THEN
             IF ((CN_ON.AND.NSTEP>1.AND.RUN_TYPE == 'NEW') .OR. & 
               (CN_ON.AND.RUN_TYPE /= 'NEW' .AND. NSTEP >= (NSTEPRST+1))) THEN
	        WRITE (*, '(/A,G12.5, A,G12.5, A,F9.3,1X,A)') ' Time = ', TIME, &
                 '  Dt = ', 2.*DT, '    CPU time left = ', TLEFT, TUNIT 
             ELSE
              WRITE (*, '(/A,G12.5, A,G12.5, A,F9.3,1X,A)') ' Time = ', TIME, &
                 '  Dt = ', DT, '    CPU time left = ', TLEFT, TUNIT 
             ENDIF
	    ENDIF
!
         ENDIF 
      ENDIF 
!
!

      call CALC_RESID_MB(0, errorpercent)

!
!     Calculate the face values of densities and mass fluxes for the first solve_vel_star call.
      CALL CONV_ROP(IER)
      CALL CALC_MFLUX (IER)
      
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
!
!
!     Calculate coefficients.  Explicitly set flags for all the quantities
!     that need to be calculated before calling CALC_COEFF.
!
      IF (Call_DQMOM) PSIZE(1:MMAX)=.TRUE.
      IF (RO_G0 == UNDEFINED) DENSITY(0) = .TRUE. 
      IF (ANY_SPECIES_EQ) RRATE = .TRUE. 
!
      VISC(0) = RECALC_VISC_G 
! 	The IF (GRANULAR_ENERGY) statement was commented out to allow calling calc_mu_s even
! 	when the algebraic granular equation is solved (not only the PDE form).
! 	This may enhance convergence. sof, March-10-2005.
       VISC(1:MMAX) = .TRUE. 
!
      CALL CALC_COEFF (DENSITY, PSIZE, SP_HEAT, VISC, COND, DIFF, RRATE, DRAGCOEF, &
         HEAT_TR, WALL_TR, IER) 

!
!     DIffusion coefficient and source terms for user-defined scalars
      IF(NScalar /= 0)CALL SCALAR_PROP(IER) 

!
!     DIffusion coefficient and source terms for K & Epsilon Eq.
      IF(K_Epsilon) CALL K_Epsilon_PROP(IER)

!      
!
!     Solve strarred velocitiy components
!

      CALL SOLVE_VEL_STAR (IER) 
!
!     Calculate density and reaction rates. Do not change density or reaction rate before the call to
!     solve_vel_star.
      IF (RO_G0 == UNDEFINED) DENSITY(0) = .TRUE. 
      IF (ANY_SPECIES_EQ) RRATE = .TRUE. 
      CALL CALC_COEFF (DENSITY, PSIZE, SP_HEAT, VISC, COND, DIFF, RRATE, DRAGCOEF, &
         HEAT_TR, WALL_TR, IER) 

!
!     Solve solids volume fraction correction equation for close-packed
!     solids phases
!
      IF(.NOT.DISCRETE_ELEMENT) THEN
        IF (MMAX > 0) THEN
          IF(MMAX == 1)THEN
            CALL CALC_K_CP (K_CP, IER)
	    CALL SOLVE_EPP (NORMS, RESS, IER)
            CALL CORRECT_1 (IER) 
          ELSE
            DO M=1,MMAX 
!   	      IF (M .EQ. MCp) THEN !Volume fraction correction technique for multiparticle types is 
   	      IF (.FALSE.) THEN    !not implemented.  This will only slow down convergence.
                CALL CALC_K_CP (K_CP, IER)
	        CALL SOLVE_EPP (NORMS, RESS, IER)
                CALL CORRECT_1 (IER) 

      	      ELSE
	        CALL SOLVE_CONTINUITY(M,IER)
	             
	      ENDIF
	    END DO
          ENDIF

          CALL CALC_VOL_FR (P_STAR, RO_G, ROP_G, EP_G, ROP_S, IER) 

	  abort_ier = ier.eq.1
	  call global_all_or(abort_ier)
          IF (abort_ier) THEN 
	      ier = 1
              MUSTIT = 2                           !indicates divergence 
              IF(DT/=UNDEFINED)GO TO 1000 
          ENDIF 
 
        ENDIF
      END IF
!
!
!     Calculate P_star in cells where solids continuity equation is
!     solved
!
      IF(.NOT.DISCRETE_ELEMENT) THEN
        IF (MMAX > 0) CALL CALC_P_STAR (EP_G, P_STAR, IER) 
      END IF
!
!     Calculate the face values of densities.
      CALL CONV_ROP(IER)
!
!     Solve fluid pressure correction equation
!
      IF (RO_G0 /= ZERO) CALL SOLVE_PP_G (NORMG, RESG, IER) 
!
!
!     Correct pressure, velocities, and density
!
      IF (RO_G0 /= ZERO) CALL CORRECT_0 (IER) 
     
      IF (RO_G0 == UNDEFINED) THEN
        DENSITY(0) = .TRUE. 
        CALL CALC_COEFF (DENSITY, PSIZE, SP_HEAT, VISC, COND, DIFF, RRATE, DRAGCOEF, &
         HEAT_TR, WALL_TR, IER) 
      ENDIF 

!
!  Update wall velocities
!  modified by sof to force wall functions even when NSW or FSW are declared
!  default wall BC will still be treated as NSW and no wall functions will be used
 
      IF(.NOT. K_EPSILON) CALL SET_WALL_BC (IER) 
!
!     Calculate the face values of mass fluxes.
      CALL CALC_MFLUX (IER)
!
!     Solve energy equations
!
      IF (ENERGY_EQ) CALL SOLVE_ENERGY_EQ (IER) 
!
!     Solve granular energy equation
!
      IER = 0
      IF (GRANULAR_ENERGY) CALL SOLVE_GRANULAR_ENERGY (IER) 

      abort_ier = ier.eq.1
      call global_all_or(abort_ier)
      IF (abort_ier) THEN
         ier = 1
         MUSTIT = 2                              !indicates divergence 
         IF(DT/=UNDEFINED)GO TO 1000 
!
      ENDIF 

      
!
!     Solve species mass balance equations
!
      CALL SOLVE_SPECIES_EQ (IER) 
!
!     Solve other scalar transport equations
!
      IF(NScalar /= 0) CALL SOLVE_Scalar_EQ (IER) 
!
!     Solve K & Epsilon transport equations
!
      IF(K_Epsilon) CALL SOLVE_K_Epsilon_EQ (IER) 
      
!
!    User-defined linear equation solver parameters may be adjusted after
!    the first iteration
!
      IF (.NOT.CYCLIC) LEQ_ADJUST = .TRUE. 
!
!
!     Check for convergence
!
      call CALC_RESID_MB(1, errorpercent)
      CALL CHECK_CONVERGENCE (NIT, errorpercent(0), MUSTIT, IER) 
      
      IF(CYCLIC)THEN
        IF(MUSTIT==0 .OR. NIT >= MAX_NIT) CALL GoalSeekMassFlux(NIT, MUSTIT, .true.)
        IF(AUTOMATIC_RESTART) RETURN
      ENDIF

!
!      If not converged continue iterations; else exit subroutine.
!
 1000 CONTINUE 
!
!     Display residuals
!
      IF (FULL_LOG) CALL DISPLAY_RESID (NIT, IER) 
      
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
            call CALC_RESID_MB(1, errorpercent)
            CALL GET_SMASS (SMASS) 
            IF (ENERGY_EQ) CALL GET_HLOSS (HLOSS) 
!
!
            IF (ENERGY_EQ) THEN 
               WRITE (UNIT_LOG, 5000) TIME, DT, NIT, SMASS, errorpercent(0), HLOSS, CPU_NOW 
               IF(FULL_LOG.and.myPE.eq.PE_IO) &          
                       WRITE(*,5000)TIME,DT,NIT,SMASS,errorpercent(0), HLOSS,CPU_NOW        !//
            ELSE 
               WRITE (UNIT_LOG, 5001) TIME, DT, NIT, SMASS, errorpercent(0), CPU_NOW 
               IF (FULL_LOG .and. myPE.eq.PE_IO) &
                       WRITE (*, 5001) TIME, DT, NIT, SMASS, errorpercent(0), CPU_NOW       !//
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
            call CALC_RESID_MB(1, errorpercent)
            WRITE (UNIT_LOG, 5200) TIME, DT, NIT, errorpercent(0) 
            CALL END_LOG 

            if (myPE.eq.PE_IO) WRITE (*, 5200) TIME, DT, NIT, errorpercent(0)   !//
         ENDIF 
         IER = 1 
         RETURN  
      ENDIF 
!
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
      if (myPE.eq.PE_IO) WRITE (UNIT_OUT, 5100) TIME, DT, NIT, SMASS    !//
      CALL START_LOG 
      WRITE (UNIT_LOG, 5100) TIME, DT, NIT, SMASS 
      CALL END_LOG 
!
      IER = 0 
      RETURN  
 5000 FORMAT(1X,'t=',F10.4,' Dt=',G10.4,' NIT=',I3,' Sm=',G10.5,'MbErr%=', G10.4,' Hl=',G12.5,&
         T84,'CPU=',F8.0,' s') 
 5001 FORMAT(1X,'t=',F10.4,' Dt=',G10.4,' NIT=',I3,' Sm=',G10.5,'MbErr%=', G10.4,T84,'CPU=',F8.0&
         ,' s') 
 5050 FORMAT(5X,'Average ',A,G12.5) 
 5060 FORMAT(5X,'Average ',A,I2,A,G12.5) 
 5100 FORMAT(1X,'t=',F10.4,' Dt=',G10.4,' NIT>',I3,' Sm= ',G10.5, 'MbErr%=', G10.4) 
 5200 FORMAT(1X,'t=',F10.4,' Dt=',G10.4,' NIT=',&
      I3,'MbErr%=', G10.4, ': Run diverged/stalled :-(') 

      END SUBROUTINE ITERATE 
!
!
      SUBROUTINE GET_TUNIT(TLEFT, TUNIT) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      DOUBLE PRECISION TLEFT 
      CHARACTER TUNIT*4 
!-----------------------------------------------
!
!
      IF (TLEFT < 3600.0d0) THEN 
         TUNIT = 's' 
      ELSE 
         TLEFT = TLEFT/3600.0d0 
         TUNIT = 'h' 
         IF (TLEFT >= 24.) THEN 
            TLEFT = TLEFT/24.0d0 
            TUNIT = 'days' 
         ENDIF 
      ENDIF 
!
      RETURN  
      END SUBROUTINE GET_TUNIT 



!  Purpose:  In the following module the mass flux across a periodic domain
!            with pressure drop is held constant at a user-specified value.  
!            This module is activated only if the user specifies a value for
!            the keyword Flux_g in the datafile.
      subroutine GoalSeekMassFlux(NIT, MUSTIT, doit)

      USE bc
      USE geometry
      USE constant
      USE compar 
      USE run
      USE time_cpu 
      IMPLICIT NONE

      INTEGER, PARAMETER :: MAXOUTIT = 500
      DOUBLE PRECISION, PARAMETER          :: omega = 0.9
      DOUBLE PRECISION, PARAMETER          :: TOL = 1E-03
     
      
      INTEGER :: NIT, MUSTIT
      INTEGER, SAVE :: OUTIT
      
      LOGICAl :: doit
      LOGICAl, SAVE :: firstPass = .true.
      
      DOUBLE PRECISION, Save  :: mdot_n, mdot_nm1, delp_n, delp_nm1, err
      DOUBLE PRECISION          :: mdot_0, delp_xyz 
      
      CHARACTER, Save  :: Direction
      
      DOUBLE PRECISION , EXTERNAL :: VAVG_Flux_U_G, VAVG_Flux_V_G, VAVG_Flux_W_G
      LOGICAL , EXTERNAL :: IsNan
!
      
      IF(CYCLIC_X_MF)THEN
        delp_n = delp_x
      ELSEIF(CYCLIC_Y_MF)THEN
        delp_n = delp_y
      ELSEIF(CYCLIC_Z_MF)THEN
        delp_n = delp_z
      ELSE
        return
      ENDIF
     
      if(.not.doit) then
        OUTIT = 0
	return
      endif
      
      OUTIT = OUTIT + 1
      if(OUTIT > MAXOUTIT) then
        Write(*,5400) MAXOUTIT
        call mfix_exit(0)
      endif
      
      mdot_0 = Flux_g
      
      
      ! calculate the average gas mass flux and error
      IF(CYCLIC_X_MF)THEN
        mdot_n = VAVG_Flux_U_G()
      ELSEIF(CYCLIC_Y_MF)THEN
        mdot_n = VAVG_Flux_V_G()
      ELSEIF(CYCLIC_Z_MF)THEN
        mdot_n = VAVG_Flux_W_G()
      ENDIF
      
      If (isNan(mdot_n) .or. isNan(delp_n)) then
        write(*,*) mdot_n, delp_n, ' NaN being caught in GoalSeekMassFlux '
        AUTOMATIC_RESTART = .TRUE.
	RETURN
      ENDIF
! 
      err = abs((mdot_n - mdot_0)/mdot_0)
      if( err < TOL)then
        MUSTIT = 0
      else
        MUSTIT = 1
	NIT = 1
      endif
      
      ! correct delp
      if(.not.firstPass)then
!        delp_xyz = delp_n - omega * (delp_n - delp_nm1) * (mdot_n - mdot_0) &
!	                   / (mdot_n - mdot_nm1)
! Fail-Safe Newton's method (below) works better than the regular Newton method (above)
!
        delp_xyz = delp_n - omega * (delp_n - delp_nm1) * ((mdot_n - mdot_0)/(mdot_nm1 - mdot_0)) &
	                   / ((mdot_n - mdot_0)/(mdot_nm1 - mdot_0) - ONE)
      else
        firstPass=.false.
        delp_xyz = delp_n*0.99
      endif
      IF(myPE.eq.PE_IO) Write(*,5500) OUTIT, err, delp_xyz, mdot_n
     
      mdot_nm1 = mdot_n
      delp_nm1 = delp_n
      
      IF(CYCLIC_X_MF)THEN
        delp_x = delp_xyz
      ELSEIF(CYCLIC_Y_MF)THEN
        delp_y = delp_xyz
      ELSEIF(CYCLIC_Z_MF)THEN
        delp_z = delp_xyz
      ENDIF
!
      
      return
5400 FORMAT(/1X,70('*')//' From: GoalSeekMassFlux',/&
      ' Message: Number of outer iterations exceeded ', I4,/1X,70('*')/) 
5500  Format('  MassFluxIteration:', I4, ' Err=', G12.5, ' DelP=', G12.5, ' Flux=', G12.5)
     
    
      end subroutine GoalSeekMassFlux
 
