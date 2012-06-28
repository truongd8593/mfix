!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!     Subroutine: DES_TIME_MARCH                                       !
!                                                                      !
!     Purpose: Main DEM driver routine                                 !
!                                                                      !
!     Author: Jay Boyalakuntla                        Date: 21-Jun-04  !
!     Reviewer: Sreekanth Pannala                     Date: 09-Nov-06  !
!     Reviewer: Rahul Garg                            Date: 01-Aug-07  !
!                                                                      !
!     Comments:                                                        !
!        Called in model/time_march.f to do DEM calculations.          !
!        do_nsearch has to be set for calling neighbour. Since         !
!        this flag is used during exchange it has to be set before     !
!        calling particle_in_cell (Pradeep G.)                         !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE DES_TIME_MARCH
     
!------------------------------------------------
! Modules
!------------------------------------------------
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
      USE drag
      USE rxns
      USE compar     
      USE time_cpu  
      USE discretelement   
      USE constant
      USE sendrecv
      USE des_bc
      USE cutcell 
      USE mppic_wallbc
      USE mfix_pic
      Use des_thermo
      Use des_rxns
      Use interpolation
      
      IMPLICIT NONE
!------------------------------------------------
! Local variables
!------------------------------------------------
! indices for grid 
      INTEGER :: I, J, K, IJK
! particle loop counter index
      INTEGER :: NP
! time step loop counter index
      INTEGER :: NN      
! loop counter index for any initial particle settling incoupled cases
      INTEGER :: FACTOR
! boundary condition loop counter index
      INTEGER :: BCV_I      
! accounted for particles
      INTEGER :: PC      

! Local variables to keep track of time when dem restart and des
! write data need to be written when des_continuum_coupled is F
      DOUBLE PRECISION DES_RES_TIME, DES_SPX_TIME

! Temporary variables when des_continuum_coupled is T to track
! changes in solid time step 
      DOUBLE PRECISION TMP_DTS, DTSOLID_TMP 
      CHARACTER*5 FILENAME

! Temporary variable used to track reporting frequency for the
! maximum overlap and maximum no. neighbors for given NN loop
      DOUBLE PRECISION DES_TMP_TIME

! Logical to see whether this is the first entry to this routine
      LOGICAL,SAVE:: FIRST_PASS = .TRUE.

! Variables needed for calculating new interpolation quantities for
! species and energy equations
      INTEGER INTERP_IJK(2**DIMN)
      DOUBLE PRECISION INTERP_WEIGHTS(2**DIMN)
      
! Identifies that the indicated particle is of interest for debugging
      LOGICAL FOCUS
!-----------------------------------------------
! Include statement functions
!-----------------------------------------------
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'fun_avg2.inc'
!-----------------------------------------------



! if first_pass
!----------------------------------------------------------------->>>
      IF(FIRST_PASS.AND..NOT.MPPIC) THEN 
         IF(DMP_LOG) WRITE(UNIT_LOG,'(1X,A)')&
            '---------- FIRST PASS DES_TIME_MARCH ---------->'
         S_TIME = ZERO
         FIRST_PASS = .FALSE.

! When no particles are present, skip the startup routines that loop over particles.
! That is, account for a setup that starts with no particles in the system.
         IF(PARTICLES /= 0) THEN
            IF(DO_NSEARCH) CALL NEIGHBOUR         
         
! To do only des in the 1st time step only for a new run so the 
! particles settle down before the coupling is turned on
            IF(RUN_TYPE == 'NEW') THEN
               IF(DES_CONTINUUM_COUPLED.AND.(.NOT.USE_COHESION)) THEN
                  DES_CONTINUUM_COUPLED = .FALSE.
                  DO FACTOR = 1, NFACTOR
                     IF (FACTOR .EQ. 1) THEN
                        IF (DMP_LOG) &
                           WRITE(UNIT_LOG,'(3X,A,/,5X,A,X,I10,X,A)') &
                           'FIRST PASS in DES_TIME_MARCH for new runs',&
                           'DEM settling period performed', NFACTOR, &
                           'times'
                     ENDIF
! calculate forces
                     CALL CALC_FORCE_DES
! update particle position/velocity
                     CALL CFNEWVALUES
! set the flag do_nsearch before calling particle in cell (for mpi)
                     IF(MOD(FACTOR,NEIGHBOR_SEARCH_N).EQ.0) &
                        DO_NSEARCH = .TRUE.
! find particles on grid                        
                     CALL PARTICLES_IN_CELL
! perform neighbor search                     
                     IF(DO_NSEARCH) CALL NEIGHBOUR

                  ENDDO
                  DES_CONTINUUM_COUPLED = .TRUE.
                  IF(DMP_LOG) WRITE(UNIT_LOG,'(3X,A)') &
                     'END DEM settling period'
               ENDIF   ! end if coupled and no cohesion

! this write_des_data is needed to properly show the initial state of
! the simulation (granular or coupled). In the coupled case, the
! particles may have 'settled' according to above.  In the granular
! case, the initial state won't be written until after the particles
! have moved without this call.
               CALL WRITE_DES_DATA

               IF(DMP_LOG) WRITE(UNIT_LOG,'(3X,A,X,ES15.5)') &
                  'DES data file written at time =', S_TIME
            ENDIF   ! end if on new run type

            IF(DMP_LOG) WRITE(UNIT_LOG,'(1X,A)')&
               '<---------- END FIRST PASS DES_TIME_MARCH ----------'

         ENDIF   ! end if particles /= 0         
      ENDIF    ! end if first pass
! end if first_pass
!-----------------------------------------------------------------<<<


! In case of restarts assign S_TIME from MFIX TIME 
      S_TIME = TIME
      TMP_DTS = ZERO
      DTSOLID_TMP = ZERO

      IF(DES_CONTINUUM_COUPLED) THEN
         IF(DT.GE.DTSOLID) THEN
            FACTOR = CEILING(real(DT/DTSOLID))
         ELSE
            FACTOR = 1
            DTSOLID_TMP = DTSOLID
            DTSOLID = DT
         ENDIF
         IF(DMP_LOG) WRITE(UNIT_LOG, 1999) factor, s_time, dt, dtsolid, pip
         IF(DMP_LOG) WRITE(*, 1999) factor, s_time, dt, dtsolid, pip
         
      ELSE
         FACTOR = CEILING(real((TSTOP-TIME)/DTSOLID)) 
         IF(DMP_LOG) WRITE(*,'(1X,A)')&
            '---------- START DES_TIME_MARCH ---------->'
         IF(DMP_LOG) WRITE(*,'(3X,A,X,I10,X,A)') &
            'DEM SIMULATION will be called', FACTOR, 'times'
! Initialization for des_spx_time, des_res_time         
         IF(RUN_TYPE .EQ. 'NEW') THEN
            DES_SPX_TIME =  S_TIME
            DES_RES_TIME =  S_TIME
         ELSE
            DES_SPX_TIME = ( INT((S_TIME+0.1d0*DTSOLID)/DES_SPX_DT) +&
               1 ) * DES_SPX_DT
            DES_RES_TIME = ( INT((S_TIME+0.1d0*DTSOLID)/DES_RES_DT) +&
               1 ) * DES_RES_DT
         ENDIF
      ENDIF   ! end if/else (des_continuum_coupled)

      IF (DES_CONTINUUM_COUPLED) DES_SPX_DT = SPX_DT(1)
      IF (RUN_TYPE .EQ. 'NEW') THEN
         DES_TMP_TIME = S_TIME
      ELSE
         DES_TMP_TIME = ( INT((S_TIME+0.1d0*DTSOLID)/DES_SPX_DT) +1 ) *&
            DES_SPX_DT
      ENDIF

      IF(MPPIC) THEN 
! compute the gas-phase pressure gradient at the beginning of the 
! des loop as the gas-phase pressure field will not change during 
! des calls
         IF(DES_CONTINUUM_COUPLED) CALL COMPUTE_PG_GRAD
      ELSE
         IF(DES_CONTINUUM_COUPLED) CALL COMPUTE_PG_GRAD
      ENDIF 

! Main DEM time loop
!----------------------------------------------------------------->>>
      DO NN = 1, FACTOR 

         IF(DES_CONTINUUM_COUPLED) THEN
            IF(S_TIME.GE.(TIME+DT)) EXIT
! If the current time in the discrete loop exceeds the current time in
! the continuum simulation, exit the discrete loop
            IF((S_TIME+DTSOLID).GT.(TIME+DT)) THEN 
! If next time step in the discrete loop will exceed the current time 
! in the continuum simulation, modify the discrete time step so final
! time will match 
               TMP_DTS = DTSOLID
               DTSOLID = TIME + DT - S_TIME
            ENDIF 
            IF(DEBUG_DES) THEN
              IF(DMP_LOG) write(*,'(3X,A,X,I10,X,A,X,ES15.7)') &
                  'DES-COUPLED LOOP NO. =', NN, ' S_TIME =', S_TIME 
              IF(DMP_LOG) write(*,'(3X,A,X,ES15.7)') &
                  'DTSOLID =', DTSOLID
            ENDIF
         ELSE   ! else if (des_continuum_coupled)
            IF(DEBUG_DES) THEN
               IF(DMP_LOG) WRITE(*,'(3X,A,X,I10,X,A,X,ES15.7)') &
               'DEM LOOP NO. =', NN, ' S_TIME =', S_TIME 
            ENDIF             
         ENDIF   ! end if/else (des_continuum_coupled) 
         
         
! communication between processors have to take place all the time;
! regardless of number of particles 
         IF(MPPIC) THEN 
            CALL MPPIC_COMPUTE_PS_GRAD            
            IF(DES_CONTINUUM_COUPLED) THEN
               CALL CALC_DES_DRAG_GS
               CALL CALC_DES_ROP_S
            ENDIF
            CALL CFUPDATEOLD
         ELSE 
            CALL CALC_FORCE_DES
         ENDIF
    
! Loop over all particles ---------------------------------------->>>
         PC = 1
         DO NP = 1, MAX_PIP
            IF(PC .GT. PIP) EXIT
            IF(.NOT.PEA(NP,1)) CYCLE
            IF(PEA(NP,4)) CYCLE

! Reset the debug flag
            FOCUS = .FALSE.
! Set the debugging flag
            IF(DEBUG_DES .AND. NP.EQ.FOCUS_PARTICLE) FOCUS = .TRUE.

! Calculate time dependent physical properties
            CALL DES_PHYSICAL_PROP(NP, FOCUS)
! Calculate cell-center interpolation weights and determine the 
! associated IJK values for the cells accounting for boundary
! conditions.
            IF(DES_INTERP_ON .AND. &
               (ANY_DES_SPECIES_EQ .OR. DES_CONV_EQ)) THEN
               INTERP_IJK(:) = -1
               INTERP_WEIGHTS(:) = ZERO
               CALL INTERPOLATE_CC(NP, INTERP_IJK, INTERP_WEIGHTS, &
                  FOCUS)
            ENDIF
! Calculate thermodynamic energy exchange
            IF(DES_ENERGY_EQ) CALL CALC_THERMO_DES(NP, &
               INTERP_IJK, INTERP_WEIGHTS, FOCUS)
! Calculate reaction rates and interphase mass transfer
            IF(ANY_DES_SPECIES_EQ) CALL DES_RRATES(NP, &
               INTERP_IJK, INTERP_WEIGHTS, FOCUS, 'SOLIDS')
! Increment the particle counter
            PC = PC + 1
         ENDDO
! End loop over all particles ------------------------------------<<<    
            
         CALL CFNEWVALUES
    
! Loop over all particles ---------------------------------------->>>
         PC = 1
         DO NP = 1, MAX_PIP
            IF(PC .GT. PIP) EXIT
            IF(.NOT.PEA(NP,1)) CYCLE
            IF(PEA(NP,4)) CYCLE            

! Reset the debug flag
            FOCUS = .FALSE.
! Set the debugging flag
            IF(DEBUG_DES .AND. NP.EQ.FOCUS_PARTICLE) FOCUS = .TRUE.

! Update particle temperature
            IF(DES_ENERGY_EQ) &
               CALL DES_THERMO_NEWVALUES(NP, FOCUS)
! Update particle from reactive chemistry process.
            IF(DES_SPECIES_EQ(PIJK(NP,5))) &
               CALL DES_REACTION_MODEL(NP, FOCUS)
! Increment the particle counter
            PC = PC + 1
         ENDDO
! End Loop over all particles ------------------------------------<<<


! Impose the wall-particle boundary condition for mp-pic case 
         IF(MPPIC) CALL MPPIC_APPLY_WALLBC 

! For systems with inlets/outlets check to determine if a particle has
! fully entered or exited the domain.  If the former, remove the status
! of 'new' and if the latter, remove the particle.
         IF (DES_MIO) CALL DES_CHECK_PARTICLE

! set do_nsearch before calling particle_in_cell
         IF(NN.EQ.1 .OR. MOD(NN,NEIGHBOR_SEARCH_N).EQ.0) &
            DO_NSEARCH = .TRUE.
         CALL PARTICLES_IN_CELL
            
         IF (DO_NSEARCH.AND.(.NOT.MPPIC)) THEN
            IF(DEBUG_DES) THEN 
               IF(DMP_LOG) WRITE(UNIT_LOG,'(3X,A,I10,/,5X,A,I10)') &
                  'Calling NEIGHBOUR: during iteration NN =', NN
            ENDIF
            CALL NEIGHBOUR
         ENDIF

! Update time to reflect changes 
         S_TIME = S_TIME + DTSOLID

! When coupled the granular temperature subroutine is only calculated at end 
! of the current DEM simulation 
         IF(DES_CONTINUUM_COUPLED .AND. NN.EQ.FACTOR) &
            CALL DES_GRANULAR_TEMPERATURE

! When coupled, all write calls are made in time_march (the continuum 
! portion) according to user settings for spx_time and res_time.
! The following section targets data writes for DEM only cases:
         IF(.NOT.DES_CONTINUUM_COUPLED) THEN    
! Keep track of TIME for DEM simulations
            TIME = S_TIME
 
! Write data using des_spx_time and des_res_time; note the time will
! reflect current position of particles  
            IF(PRINT_DES_DATA) THEN
               IF ( (S_TIME+0.1d0*DTSOLID >= DES_SPX_TIME) .OR. &
                    (S_TIME+0.1d0*DTSOLID >= TSTOP) .OR. &
                    (NN == FACTOR) ) THEN
                  DES_SPX_TIME = &
                     ( INT((S_TIME+0.1d0*DTSOLID)/DES_SPX_DT) &
                     + 1 )*DES_SPX_DT
! Granular temperature subroutine should be called/calculated when
! writing DES data 
                  CALL DES_GRANULAR_TEMPERATURE
                  CALL WRITE_DES_DATA
                  IF(DMP_LOG) WRITE(UNIT_LOG,'(3X,A,X,ES15.5)') &
                     'DES data file written at time =', S_TIME
               ENDIF
            ENDIF

            IF ( (S_TIME+0.1d0*DTSOLID >= DES_RES_TIME) .OR. &
                 (S_TIME+0.1d0*DTSOLID >= TSTOP) .OR. &
                 (NN == FACTOR) ) THEN
               DES_RES_TIME = &
                  ( INT((S_TIME+0.1d0*DTSOLID)/DES_RES_DT) &
                  + 1 )*DES_RES_DT
                  call des_write_restart
! Write RES1 here since it won't be called in time_march.  This will
! also keep track of TIME
               CALL WRITE_RES1 
               IF(DMP_LOG) WRITE(UNIT_LOG,'(3X,A,X,ES15.5)') &
               'DES.RES and .RES files written at time =', S_TIME
            ENDIF
         ENDIF  ! end if (.not.des_continuum_coupled)


! J.Musser : mass inlet/outlet -> particles entering the system
         IF(DES_MI)THEN 
            DO BCV_I = 1, DES_BCMI
               IF(PI_FACTOR(BCV_I) .GT. 1)THEN
                  IF(DES_MI_TIME(BCV_I) .LE. S_TIME) THEN   !Verify start time
                     CALL DES_MASS_INLET(BCV_I)
                     DES_MI_TIME(BCV_I) = S_TIME + PI_FACTOR(BCV_I) *&
                        DTSOLID
                  ENDIF
               ELSE
                  CALL DES_MASS_INLET(BCV_I)
               ENDIF
            ENDDO
         ENDIF

! Report some statistics on overlap and neighbors to screen log
         IF ( (S_TIME+0.1d0*DTSOLID >= DES_TMP_TIME) .OR. &
              ( (S_TIME+0.1d0*DTSOLID >= TSTOP) .AND. &
               (.NOT.DES_CONTINUUM_COUPLED) ) .OR. &          
              (NN .EQ. FACTOR) ) THEN
            DES_TMP_TIME = ( INT((S_TIME+0.1d0*DTSOLID)/DES_SPX_DT) &
               + 1 )*DES_SPX_DT
         ENDIF

      ENDDO     ! end do NN = 1, FACTOR
! END DEM time loop
!-----------------------------------------------------------------<<<      


! When coupled, and if needed, reset the discrete time step accordingly
      IF(DT.LT.DTSOLID_TMP) THEN
         DTSOLID = DTSOLID_TMP
      ENDIF

      IF(TMP_DTS.NE.ZERO) THEN
         DTSOLID = TMP_DTS
         TMP_DTS = ZERO
      ENDIF

      IF(MPPIC) THEN 
         IF(DMP_LOG) WRITE(UNIT_LOG, 2000) DTSOLID, &
            DTPIC_CFL, DTPIC_TAUP, MIN(DTPIC_CFL, DTPIC_TAUP)
         IF(myPE.eq.pe_IO) WRITE(*, 2000) DTSOLID, & 
            DTPIC_CFL, DTPIC_TAUP, MIN(DTPIC_CFL, DTPIC_TAUP)

         DTPIC_MAX = MIN(DTPIC_CFL, DTPIC_TAUP)

         IF(DTSOLID.GT.DTPIC_MAX) THEN 
            IF(DMP_LOG) WRITE(UNIT_LOG, 2001) DTSOLID
            IF(myPE.eq.pe_IO) WRITE(*, 2001) DTPIC_MAX
            DTSOLID = DTPIC_MAX
         ELSEIF(DTSOLID.LT.DTPIC_MAX) THEN 

            IF(DMP_LOG) WRITE(UNIT_LOG, 2002) DTSOLID
            IF(myPE.eq.pe_IO) WRITE(*, 2002) DTPIC_MAX
            DTSOLID = DTPIC_MAX
         ELSE
            IF(DMP_LOG) WRITE(UNIT_LOG, 2003) DTSOLID 
            IF(mype.eq.pe_IO) WRITE(*,2003) DTSOLID
         ENDIF
      ENDIF

      
      IF(.NOT.DES_CONTINUUM_COUPLED)THEN
         IF(DMP_LOG) WRITE(UNIT_LOG,'(1X,A)')&
         '<---------- END DES_TIME_MARCH ----------'
      ELSE 
         call send_recv(ep_g,2)
         call send_recv(rop_g,2)
         call send_recv(des_u_s,2)
         call send_recv(des_v_s,2) 
         if(dimn.eq.3) call send_recv(des_w_s,2) 
         call send_recv(rop_s,2)
      ENDIF

 1999    FORMAT(/1X,70('- '),//,10X,  & 
         & 'DEM SIMULATION CALLED ', 2x, i5, 2x, 'times this fluid step', /10x &
         & 'S_TIME, DT, DTSOLID and PIP = ', 3(2x,g17.8),2x, i10)

 2000 FORMAT(/10x, & 
      & 'ADAPTING DTSOLID FOR MPPIC', /10x, &
      & 'DTSOLID CURRENT  = ', g17.8, /10x,  &
      & 'DTPIC_CFL :  = ', g17.8, /10x,  &
      & 'DTPIC TAUP:  = ', g17.8, /10x, &
      & 'DTPIC_MAX :  = ', g17.8)

 2001 FORMAT(/10x, & 
      & 'REDUCING CURRENT DTSOLID TO', g17.8)
      
 2002 FORMAT(/10x, & 
      & 'INCREASING CURRENT DTSOLID TO', g17.8)

 2003 FORMAT(/10x, & 
      & 'DTSOLID REMAINS UNCHANGED AT = ', g17.8)

      RETURN

      END SUBROUTINE DES_TIME_MARCH

