!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: ITERATE                                                 C
!  Purpose: This module controls the iterations for solving equations  C
!                                                                      C
!  Author: M. Syamlal                                 Date: 12-APR-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: To incorporate DES flags so that the solids                C
!  calculations are not done using Continuum when doing DES            C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!                                                                      C
!  Revision Number: 2                                                  C
!  Purpose: To incorporate Cartesian grid modifications                C
!  and utilization of the dashboard                                    C
!  Author: Jeff Dietiker                              Date: 01-Jul-09  C
!      
!  Revision Number: 3                                                  C
!  Purpose: Incorporation of QMOM for the solution of the particle     C
!  kinetic equation                                                    C
!  Author: Alberto Passalacqua - Fox Research Group   Date: 02-Dec-09  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE ITERATE(IER, NIT) 

!-----------------------------------------------
! Modules
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
      USE residual
      USE cutcell
      USE vtk
      USE dashboard
      USE qmom_kinetic_equation
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Error index 
      INTEGER, INTENT(INOUT) :: IER 
! Number of iterations 
      INTEGER, INTENT(OUT) :: NIT 
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! current cpu time used 
      DOUBLE PRECISION :: CPU_NOW 
! cpu time left      
      DOUBLE PRECISION :: TLEFT
! flag indicating convergence status with MUSTIT = 0,1,2 implying
! complete convergence, non-covergence and divergence respectively
      INTEGER :: MUSTIT 
! Normalization factor for gas & solids pressure residual 
      DOUBLE PRECISION :: NORMg, NORMs
! Set normalization factor for gas and solids pressure residual
      LOGICAL :: SETg, SETs 
! gas & solids pressure residual 
      DOUBLE PRECISION :: RESg, RESs
! Weight of solids in the reactor 
      DOUBLE PRECISION :: SMASS 
! Heat loss from the reactor 
      DOUBLE PRECISION :: HLOSS 
! phase index 
      INTEGER :: M
! average velocity
      DOUBLE PRECISION :: Vavg 
      DOUBLE PRECISION :: errorpercent(0:MMAX)
      LOGICAL :: ABORT_IER
      CHARACTER*4 TUNIT 
!-----------------------------------------------
! External functions
!-----------------------------------------------
      DOUBLE PRECISION, EXTERNAL :: VAVG_U_G, VAVG_V_G, VAVG_W_G, &
                                    VAVG_U_S, VAVG_V_S, VAVG_W_S
! use the subroutine from machine.f
      EXTERNAL CPU_TIME   
!-----------------------------------------------
! Include statement functions
!-----------------------------------------------
!-----------------------------------------------

! initializations
      DT_prev = DT
      NIT = 0 
      RESG = ZERO 
      RESS = ZERO 

      IF (NORM_G == UNDEFINED) THEN 
         NORMG = ONE 
         SETG = .FALSE. 
      ELSE 
         NORMG = NORM_G 
         SETG = .TRUE. 
      ENDIF 

      IF (NORM_S == UNDEFINED) THEN 
         NORMS = ONE 
         SETS = .FALSE. 
      ELSE 
         NORMS = NORM_S 
         SETS = .TRUE. 
      ENDIF 

      LEQ_ADJUST = .FALSE. 


! Initialize residuals
      CALL INIT_RESID (IER) 


! Initialize the routine for holding gas mass flux constant with cyclic bc
      IF(CYCLIC) CALL GoalSeekMassFlux(0, 0, .false.)


! CPU time left
      IF (FULL_LOG) THEN  
         TLEFT = (TSTOP - TIME)*CPUOS 
         CALL GET_TUNIT (TLEFT, TUNIT) 

         IF (DT == UNDEFINED) THEN 
            CALL GET_SMASS (SMASS) 
            IF(myPE.eq.PE_IO) THEN
               WRITE (*, '(/A,G10.5, A,F9.3,1X,A)') &
                  ' Starting solids mass = ', SMASS, &
                  '    CPU time left = ', TLEFT, TUNIT 
            ENDIF
         ELSE 
            IF(myPE.eq.PE_IO) THEN
               IF ((CN_ON.AND.NSTEP>1.AND.RUN_TYPE == 'NEW') .OR. & 
                   (CN_ON.AND.RUN_TYPE /= 'NEW' .AND.&
                    NSTEP >= (NSTEPRST+1))) THEN
                  WRITE (*, '(/A,G12.5, A,G12.5, A,F9.3,1X,A)')&
                     ' Time = ', TIME, '  Dt = ', 2.D0*DT,&
                     '    CPU time left = ', TLEFT, TUNIT 
               ELSE
                  WRITE (*, '(/A,G12.5, A,G12.5, A,F9.3,1X,A)') &
                     ' Time = ', TIME, '  Dt = ', DT, &
                     '    CPU time left = ', TLEFT, TUNIT 
               ENDIF
            ENDIF
         ENDIF   ! if/else(dt==undefined)
      ENDIF   ! if(full_log) 

      CALL CALC_RESID_MB(0, errorpercent)

! Calculate the face values of densities and mass fluxes for the first
! solve_vel_star call.
      CALL CONV_ROP(IER)
      CALL CALC_MFLUX (IER)
      CALL SET_BC1

! JFD: modification for cartesian grid implementation
      IF(CARTESIAN_GRID) CALL CG_SET_OUTFLOW


! Begin iterations
!-----------------------------------------------------------------
   50 CONTINUE 
      MUSTIT = 0 
      NIT = NIT + 1 
      PHIP_OUT_ITER=NIT ! To record the output of phip
! mechanism to set the normalization factor for the correction
! after the first iteration to the corresponding residual found 
! in the first iteration
      IF (.NOT.SETG) THEN
         IF (RESG > SMALL_NUMBER) THEN 
            NORMG = RESG 
            SETG = .TRUE. 
         ENDIF 
      ENDIF 
      IF (.NOT.SETS) THEN 
         IF (RESS > SMALL_NUMBER) THEN 
            NORMS = RESS 
            SETS = .TRUE. 
         ENDIF 
      ENDIF 

! Call user-defined subroutine to set quantities that need to be updated
! every iteration
      IF (CALL_USR) CALL USR2 


! Calculate coefficients.  First turn off all flags.  Then explicitly set
! flags for all the quantities that need to be calculated.
      CALL TurnOffCOEFF(DENSITY, PSIZE, SP_HEAT, VISC, COND, DIFF, &
               GRAN_DISS, RRATE, DRAGCOEF, HEAT_TR, WALL_TR, IER)

      IF (Call_DQMOM) PSIZE(1:SMAX)=.TRUE.
!      IF (RO_G0 == UNDEFINED) DENSITY(0) = .TRUE. 
      WALL_TR = .TRUE. 
      IF (ENERGY_EQ) THEN 
         SP_HEAT(:SMAX) = .TRUE. 
         COND(:SMAX) = .TRUE. 
         HEAT_TR(:SMAX,:SMAX) = .TRUE. 
      ENDIF 
      IF(ANY_SPECIES_EQ) DIFF(:SMAX) = .TRUE.
      DRAGCOEF(:MMAX,:MMAX) = .TRUE. 
      VISC(0) = RECALC_VISC_G 
      VISC(1:MMAX) = .TRUE. 

      IF (TRIM(KT_TYPE) .EQ. 'IA_NONEP' .OR. &
          TRIM(KT_TYPE) .EQ. 'GD_99') THEN
         GRAN_DISS(:MMAX) = .TRUE.
      ENDIF

      CALL PHYSICAL_PROP (DENSITY, PSIZE, SP_HEAT, IER)
      CALL TRANSPORT_PROP (VISC, COND, DIFF, GRAN_DISS, IER) 
      CALL EXCHANGE (DRAGCOEF, HEAT_TR, WALL_TR, IER) 

! Diffusion coefficient and source terms for user-defined scalars
      IF(NScalar /= 0) CALL SCALAR_PROP(IER) 

! Diffusion coefficient and source terms for K & Epsilon Eq.
      IF(K_Epsilon) CALL K_Epsilon_PROP(IER)

! Solve starred velocity components
      CALL SOLVE_VEL_STAR (IER) 

! Calculate density and reaction rates. Do not change reaction rate anywhere
! else within this iteration loop.  Fluid density can be changed after the 
! pressure correction step.
      CALL TurnOffCOEFF(DENSITY, PSIZE, SP_HEAT, VISC, COND, DIFF, &
               GRAN_DISS, RRATE, DRAGCOEF, HEAT_TR, WALL_TR, IER)
      IF (RO_G0 == UNDEFINED) DENSITY(0) = .TRUE. 
      IF (ANY_SPECIES_EQ) RRATE = .TRUE.
       
      CALL PHYSICAL_PROP (DENSITY, PSIZE, SP_HEAT, IER)
      IF (Neg_RHO_G) THEN
         MUSTIT = 2                              !indicates divergence 
         Neg_RHO_G = .FALSE.
         IF(FULL_LOG.and.myPE.eq.PE_IO) &
            WRITE(*,6000)'Negative gas density detected.'   ! JFD
         IF(DT/=UNDEFINED)GO TO 1000  
      ENDIF
      CALL CALC_RRATE(RRATE)


! Solve solids volume fraction correction equation for close-packed
! solids phases
      IF(.NOT.(DISCRETE_ELEMENT .OR. QMOMK) .OR. &
         DES_CONTINUUM_HYBRID) THEN
        IF (MMAX > 0) THEN

          IF(MMAX == 1 .AND. MCP /= UNDEFINED_I)THEN 
! if second phase (m=1) can overpack (e.g., bubbles) then solve its
! continuity equation
             CALL CALC_K_CP (K_CP, IER)
             CALL SOLVE_EPP (NORMS, RESS, IER)
             CALL CORRECT_1 (IER) 
          ELSE
             DO M=1,SMAX ! mmax -> smax for GHD theory
! if one chooses to revert back to old mark_phase_4_cor wherein the
! continuity of the gas phase can get marked to be solved then this 
! loop should start at 0.

!   	        IF (M .EQ. MCP) THEN
! Volume fraction correction technique for one of the solids phase
! is not implemented.  This will only slow down convergence.
               IF (.FALSE.) THEN    
                  CALL CALC_K_CP (K_CP, IER)
                  CALL SOLVE_EPP (NORMS, RESS, IER)
                  CALL CORRECT_1 (IER) 
               ELSE
                  CALL SOLVE_CONTINUITY(M,IER)
               ENDIF

             ENDDO
          ENDIF   ! end if/else (mmax==1 .and. mcp /= undefined)

          IF(TRIM(KT_TYPE) .eq. 'GHD') CALL ADJUST_EPS_GHD
          CALL CALC_VOL_FR (P_STAR, RO_G, ROP_G, EP_G, ROP_S, IER) 

          abort_ier = ier.eq.1
          CALL global_all_or(abort_ier)
          IF (abort_ier) THEN 
              ier = 1
              MUSTIT = 2                           !indicates divergence 
              IF(FULL_LOG.and.myPE.eq.PE_IO) &
                 WRITE(*,6000)'Negative void fraction detected.'  ! JFD
              IF(DT/=UNDEFINED) GOTO 1000 
          ENDIF 

        ENDIF  ! endif (mmax >0)
      ENDIF  ! end if (.not.discrete_element)


! Calculate P_star in cells where solids continuity equation is
! solved
      IF(.NOT.(DISCRETE_ELEMENT .OR. QMOMK) .OR. &
         DES_CONTINUUM_HYBRID) THEN
        IF (MMAX > 0 .AND. .NOT.FRICTION) &
           CALL CALC_P_STAR (EP_G, P_STAR, IER)
      ENDIF

! Calculate the face values of densities.
      CALL CONV_ROP(IER)

      IF (RO_G0 /= ZERO) THEN
! Solve fluid pressure correction equation
         CALL SOLVE_PP_G (NORMG, RESG, IER) 
! Correct pressure, velocities, and density
         CALL CORRECT_0 (IER) 
      ENDIF
     
      IF (RO_G0 == UNDEFINED) THEN
        PSIZE(1:MMAX)=.FALSE.
        SP_HEAT(:MMAX) = .FALSE. 
        DENSITY(0) = .TRUE. 
        CALL PHYSICAL_PROP (DENSITY, PSIZE, SP_HEAT, IER)
        IF (Neg_RHO_G) THEN
           MUSTIT = 2                              !indicates divergence 
           Neg_RHO_G = .FALSE.
           IF(FULL_LOG.and.DMP_LOG) &
              WRITE(*,6000)'Negative gas density detected.' ! JFD
           IF(DT/=UNDEFINED) GOTO 1000  
        ENDIF 
      ENDIF 


! Update wall velocities: 
! modified by sof to force wall functions so even when NSW or FSW are 
! declared, default wall BC will still be treated as NSW and no wall 
! functions will be used
      IF(.NOT. K_EPSILON) CALL SET_WALL_BC (IER) 

! Calculate the face values of mass fluxes
      CALL CALC_MFLUX (IER)
      CALL SET_BC1

! JFD: modification for cartesian grid implementation
      IF(CARTESIAN_GRID) CALL CG_SET_OUTFLOW

! Solve energy equations
      IF (ENERGY_EQ) CALL SOLVE_ENERGY_EQ (IER) 

! Solve granular energy equation
      IF(.NOT.DISCRETE_ELEMENT .OR. &
         DES_CONTINUUM_HYBRID) THEN
         IER = 0
         IF (GRANULAR_ENERGY) CALL SOLVE_GRANULAR_ENERGY (IER) 
         ABORT_IER = IER.eq.1
         CALL GLOBAL_ALL_OR(ABORT_IER)
         IF (ABORT_IER) THEN
            IER = 1
            MUSTIT = 2                              !indicates divergence 
            IF(DT/=UNDEFINED) GOTO 1000 
         ENDIF
      ENDIF
      
! Solve species mass balance equations
      CALL SOLVE_SPECIES_EQ (IER) 

! Solve other scalar transport equations
      IF(NScalar /= 0) CALL SOLVE_Scalar_EQ (IER) 

! Solve K & Epsilon transport equations
      IF(K_Epsilon) CALL SOLVE_K_Epsilon_EQ (IER) 
      

! User-defined linear equation solver parameters may be adjusted after
! the first iteration
      IF (.NOT.CYCLIC) LEQ_ADJUST = .TRUE. 


! Check for convergence
      CALL ACCUM_RESID ! Accumulating residuals from all the processors
      RESG = RESID(RESID_P,0)
      RESS = RESID(RESID_P,1)
      CALL CALC_RESID_MB(1, errorpercent)
      CALL CHECK_CONVERGENCE (NIT, errorpercent(0), MUSTIT, IER) 
      
      IF(CYCLIC)THEN
        IF(MUSTIT==0 .OR. NIT >= MAX_NIT) &
           CALL GoalSeekMassFlux(NIT, MUSTIT, .true.)
        IF(AUTOMATIC_RESTART) RETURN
      ENDIF


!  If not converged continue iterations; else exit subroutine.
 1000 CONTINUE 
!-----------------------------------------------------------------

! Display residuals
      IF (FULL_LOG) CALL DISPLAY_RESID (NIT, IER) 

! Determine course of simulation: converge, non-converge, diverge?
      IF (MUSTIT == 0) THEN
! ---------------------------------------------------------------->>>
         IF (DT==UNDEFINED .AND. NIT==1) GOTO 50   !Iterations converged 

! Perform checks and dump to screen every NLOG time steps
         IF (MOD(NSTEP,NLOG) == 0) THEN 
            CALL CPU_TIME (CPU_NOW)
            CPUOS = (CPU_NOW - CPU_NLOG)/(TIME - TIME_NLOG) 
            CPU_NLOG = CPU_NOW 
            TIME_NLOG = TIME 
            CPU_NOW = CPU_NOW - CPU0 

            CALL CALC_RESID_MB(1, errorpercent)
            CALL GET_SMASS (SMASS) 
            IF (ENERGY_EQ) CALL GET_HLOSS (HLOSS) 

            CALL START_LOG 
            IF (ENERGY_EQ) THEN 
               IF(DMP_LOG)WRITE (UNIT_LOG, 5000) TIME, DT, NIT, SMASS,&
                  HLOSS, CPU_NOW 
               IF(FULL_LOG.and.myPE.eq.PE_IO) &          
                  WRITE(*,5000)TIME,DT,NIT,SMASS, HLOSS,CPU_NOW        
            ELSE 
               IF(DMP_LOG)WRITE (UNIT_LOG, 5001) TIME, DT, NIT, &
                  SMASS, CPU_NOW 
               IF (FULL_LOG .and. myPE.eq.PE_IO) &
                  WRITE (*, 5001) TIME, DT, NIT, SMASS, CPU_NOW
            ENDIF 

            IF(DMP_LOG)WRITE (UNIT_LOG, 5002) (errorpercent(M), M=0,MMAX) 
            IF (FULL_LOG .and. myPE.eq.PE_IO) &
               WRITE (*, 5002) (errorpercent(M), M=0,MMAX)      
            IF (.NOT.FULL_LOG) THEN 
               TLEFT = (TSTOP - TIME)*CPUOS 
               CALL GET_TUNIT (TLEFT, TUNIT) 
               IF(DMP_LOG)WRITE (UNIT_LOG, '(46X,A,F9.3,1X,A)') &
                  '    CPU time left = ', TLEFT, TUNIT 
            ENDIF 

            IF (CYCLIC_X .OR. CYCLIC_Y .OR. CYCLIC_Z) THEN 
               IF (DO_I) THEN
                 Vavg = VAVG_U_G()
                 IF(DMP_LOG)WRITE (UNIT_LOG, 5050) 'U_g = ', Vavg
               ENDIF
               IF (DO_J) THEN
                 Vavg = VAVG_V_G()
                 IF(DMP_LOG)WRITE (UNIT_LOG, 5050) 'V_g = ',  Vavg
               ENDIF
               IF (DO_K) THEN
                 Vavg = VAVG_W_G()
                 IF(DMP_LOG)WRITE (UNIT_LOG, 5050) 'W_g = ', Vavg 
               ENDIF
               DO M = 1, SMAX 
                  IF (DO_I) Then
                    Vavg = VAVG_U_S(M)
                    IF(DMP_LOG)WRITE (UNIT_LOG, 5060) 'U_s(', M, ') = ', Vavg
                  ENDIF
                  IF (DO_J) Then
                    Vavg = VAVG_V_S(M)
                    IF(DMP_LOG)WRITE (UNIT_LOG, 5060) 'V_s(', M, ') = ', Vavg
                  ENDIF
                  IF (DO_K) Then
                    Vavg = VAVG_W_S(M) 
                    IF(DMP_LOG)WRITE (UNIT_LOG, 5060) 'W_s(', M, ') = ', Vavg
                  ENDIF
               ENDDO 
            ENDIF   ! end if cyclic_x, cyclic_y or cyclic_z

            CALL END_LOG 
         ENDIF   ! end IF (MOD(NSTEP,NLOG) == 0)

! JFD: modification for cartesian grid implementation
         IF(WRITE_DASHBOARD) THEN
            RUN_STATUS = 'In Progress...'
            N_DASHBOARD = N_DASHBOARD + 1 
            IF(MOD(N_DASHBOARD,F_DASHBOARD)==0) THEN
               TLEFT = (TSTOP - TIME)*CPUOS 
               CALL GET_TUNIT (TLEFT, TUNIT) 
               CALL UPDATE_DASHBOARD(NIT,TLEFT,TUNIT)
            ENDIF
         ENDIF

         IER = 0
         RETURN   ! for if mustit =0 (converged)
! end converged: go back to time_march
! ----------------------------------------------------------------<<<

! diverged
      ELSEIF (MUSTIT==2 .AND. DT/=UNDEFINED) THEN
! ---------------------------------------------------------------->>>
         IF (FULL_LOG) THEN 
            CALL START_LOG 
            CALL CALC_RESID_MB(1, errorpercent)
            IF(DMP_LOG) WRITE(UNIT_LOG,5200) TIME, DT, NIT, errorpercent(0) 
            CALL END_LOG 

            IF (myPE.EQ.PE_IO) WRITE(*,5200) TIME, DT, NIT, errorpercent(0)
         ENDIF 

! JFD: modification for cartesian grid implementation
         IF(WRITE_DASHBOARD) THEN
            RUN_STATUS = 'Diverged/stalled...'
            N_DASHBOARD = N_DASHBOARD + 1 
            IF(MOD(N_DASHBOARD,F_DASHBOARD)==0) THEN
               TLEFT = (TSTOP - TIME)*CPUOS 
               CALL GET_TUNIT (TLEFT, TUNIT) 
               CALL UPDATE_DASHBOARD(NIT,TLEFT,TUNIT)
            ENDIF
         ENDIF

         IER = 1 
         RETURN  ! for if mustit =2 (diverged)
      ENDIF
! end diverged: go back to time_march, decrease time step, try again      
! ----------------------------------------------------------------<<<

! not converged (mustit = 1, !=0,2 )
! ---------------------------------------------------------------->>>
      IF (NIT < MAX_NIT) THEN 
         MUSTIT = 0 
         GOTO 50 
      ENDIF ! continue iterate
! ----------------------------------------------------------------<<<




      CALL GET_SMASS (SMASS) 
      IF (myPE.eq.PE_IO) WRITE(UNIT_OUT, 5100) TIME, DT, NIT, SMASS
      CALL START_LOG
      IF(DMP_LOG) WRITE(UNIT_LOG, 5100) TIME, DT, NIT, SMASS
      CALL END_LOG 

! SOF: MFIX will not go the next time step if MAX_NIT is reached, 
! instead it will decrease the time step. (IER changed from 0 to 1)
      IER = 1
      RETURN  


 5000 FORMAT(1X,'t=',F10.4,' Dt=',G10.4,' NIT=',I3,' Sm=',G10.5,' Hl=',G12.5,&
         T84,'CPU=',F8.0,' s') 
 5001 FORMAT(1X,'t=',F10.4,' Dt=',G10.4,' NIT=',I3,' Sm=',G10.5, T84,'CPU=',F8.0,' s') 
 5002 FORMAT(3X,'MbError%(0,MMAX):', 5(1X,G10.4)) 
 5050 FORMAT(5X,'Average ',A,G12.5) 
 5060 FORMAT(5X,'Average ',A,I2,A,G12.5) 
 5100 FORMAT(1X,'t=',F10.4,' Dt=',G10.4,' NIT>',I3,' Sm= ',G10.5, 'MbErr%=', G10.4) 
 5200 FORMAT(1X,'t=',F10.4,' Dt=',G10.4,' NIT=',&
      I3,'MbErr%=', G10.4, ': Run diverged/stalled :-(') 
 6000 FORMAT(1X,A) 

      END SUBROUTINE ITERATE 



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!  Purpose:
!      
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE GET_TUNIT(TLEFT, TUNIT) 

!-----------------------------------------------
! Modules
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
      DOUBLE PRECISION, INTENT(INOUT) :: TLEFT 
      CHARACTER TUNIT*4 
!-----------------------------------------------

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

      RETURN  
      END SUBROUTINE GET_TUNIT 



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!  Purpose:  In the following subroutine the mass flux across a periodic
!            domain with pressure drop is held constant at a
!            user-specified value.  This module is activated only if 
!            the user specifies a value for the keyword flux_g in the 
!            mfix.dat file.
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

      subroutine GoalSeekMassFlux(NIT, MUSTIT, doit)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE bc
      USE geometry
      USE constant
      USE compar 
      USE run
      USE time_cpu 
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
      INTEGER, INTENT(INOUT) :: NIT, MUSTIT
      LOGICAL, INTENT(IN) :: doit
!-----------------------------------------------
! Local Variables      
!-----------------------------------------------
      INTEGER, PARAMETER :: MAXOUTIT = 500
      DOUBLE PRECISION, PARAMETER          :: omega = 0.9
      DOUBLE PRECISION, PARAMETER          :: TOL = 1E-03
      INTEGER, SAVE :: OUTIT
      LOGICAL, SAVE :: firstPass = .true.

      DOUBLE PRECISION, SAVE  :: mdot_n, mdot_nm1, delp_n, delp_nm1, err
      DOUBLE PRECISION        :: mdot_0, delp_xyz 

      CHARACTER, SAVE :: Direction
!-----------------------------------------------
! Functions      
!-----------------------------------------------
      DOUBLE PRECISION, EXTERNAL :: VAVG_Flux_U_G, VAVG_Flux_V_G, &
                                    VAVG_Flux_W_G
      LOGICAL, EXTERNAL :: IsNan
!-----------------------------------------------     

      IF(CYCLIC_X_MF)THEN
         delp_n = delp_x
      ELSEIF(CYCLIC_Y_MF)THEN
         delp_n = delp_y
      ELSEIF(CYCLIC_Z_MF)THEN
         delp_n = delp_z
      ELSE
         RETURN
      ENDIF
     
      IF(.NOT.doit) THEN
         OUTIT = 0
         RETURN
      ENDIF
      
      OUTIT = OUTIT + 1
      IF(OUTIT > MAXOUTIT) THEN
         IF (myPE.EQ.PE_IO) write(*,5400) MAXOUTIT
         CALL mfix_exit(0)
      ENDIF
      
      mdot_0 = Flux_g
      
      
      ! calculate the average gas mass flux and error
      IF(CYCLIC_X_MF)THEN
        mdot_n = VAVG_Flux_U_G()
      ELSEIF(CYCLIC_Y_MF)THEN
        mdot_n = VAVG_Flux_V_G()
      ELSEIF(CYCLIC_Z_MF)THEN
        mdot_n = VAVG_Flux_W_G()
      ENDIF
      
      IF (isNan(mdot_n) .OR. isNan(delp_n)) THEN
         IF (myPE.eq.PE_IO) write(*,*) mdot_n, delp_n, &
            ' NaN being caught in GoalSeekMassFlux '
         AUTOMATIC_RESTART = .TRUE.
         RETURN
      ENDIF
 
      err = abs((mdot_n - mdot_0)/mdot_0)
      IF( err < TOL) THEN
         MUSTIT = 0
      ELSE
        MUSTIT = 1
        NIT = 1
      ENDIF
      
! correct delp
      if(.not.firstPass)then
!        delp_xyz = delp_n - omega * (delp_n - delp_nm1) * (mdot_n - mdot_0) &
!	                   / (mdot_n - mdot_nm1)
! Fail-Safe Newton's method (below) works better than the regular
! Newton method (above)

         delp_xyz = delp_n - omega * (delp_n - delp_nm1) * &
                     ((mdot_n - mdot_0)/(mdot_nm1 - mdot_0)) / &
                     ((mdot_n - mdot_0)/(mdot_nm1 - mdot_0) - ONE)
      else
         firstPass=.false.
         delp_xyz = delp_n*0.99
      endif

      IF(MUSTIT == 0) then
        IF(myPE.eq.PE_IO) Write(*,5500) TIME, OUTIT, delp_xyz, mdot_n
      ENDIF
     
      mdot_nm1 = mdot_n
      delp_nm1 = delp_n
      
      IF(CYCLIC_X_MF)THEN
        delp_x = delp_xyz
      ELSEIF(CYCLIC_Y_MF)THEN
        delp_y = delp_xyz
      ELSEIF(CYCLIC_Z_MF)THEN
        delp_z = delp_xyz
      ENDIF

      
      RETURN

5400 FORMAT(/1X,70('*')//' From: GoalSeekMassFlux',/&
      ' Message: Number of outer iterations exceeded ', I4,/1X,70('*')/) 
5500  Format('  Time=', G12.5, ' MassFluxIterations=', I4, ' DelP=', G12.5, ' Gas Flux=', G12.5)
    
      END SUBROUTINE GoalSeekMassFlux
 
