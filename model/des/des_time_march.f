!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                         C  
!     Module name: DES_TIME_MARCH                                         C
!
!     Purpose: Called in model/time_march.f to do DES calcs
!     Main DEM driver routine
!
!                                                                         C
!     Author: Jay Boyalakuntla                           Date: 21-Jun-04  C
!     Reviewer: Sreekanth Pannala                        Date: 09-Nov-06  C
!     Reviewer: Rahul Garg                               Date: 01-Aug-07  C
!     Comments: Changed the calling rules to neighbor search routines     C
!                                                                         C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE DES_TIME_MARCH
!     
!-------------------
!     M o d u l e s
!-------------------
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

      IMPLICIT NONE
!-----------------------------------------------
!     L o c a l   V a r i a b l e s
!-----------------------------------------------
     
      INTEGER NN, FACTOR, NSN, NP, IJK, I, J, K, BCV_I

!     Temporary variables when des_continuum_coupled is F to track
!     reporting time 
      DOUBLE PRECISION PTC, DESRESDT 

!     Temporary variables when des_continuum_coupled is T to track
!     changes in solid time step 
      DOUBLE PRECISION TMP_DTS, DTSOLID_TMP 
      CHARACTER*5 FILENAME

!     Logical to see whether this is the first entry to this routine
      LOGICAL,SAVE:: FIRST_PASS = .TRUE.

!     index to track accounted for particles  
      INTEGER PC    

!     logical that tells whether to perform calculations in subroutine
      LOGICAL FLAGTEMP

!-----------------------------------------------

      INCLUDE 'function.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'fun_avg2.inc'
      
      
      IF(FIRST_PASS) THEN 
         
         FIRST_PASS = .FALSE.
         DESRESDT = 0.0d0
         NSN = 0
         S_TIME = ZERO
         TMP_DTS = ZERO
         PTC = ZERO
         INQC = INIT_QUAD_COUNT
         IF(RUN_TYPE == 'NEW') THEN
            DES_SPX_TIME =  TIME
            DES_RES_TIME =  TIME
         ELSE
            DES_SPX_TIME = (INT((TIME+DT+0.1d0*DT)/SPX_DT(1))+1)*SPX_DT(1)
            DES_RES_TIME = (INT((TIME+DT+0.1d0*DT)/RES_DT)   +1)*RES_DT
         ENDIF

         !PRINT *,'SPX TIME', SPX_DT(1), DES_SPX_TIME

         CALL NEIGHBOUR
         
! COHESION INITIALIZE
         IF(USE_COHESION)THEN
            CALL INITIALIZE_COHESION_PARAMETERS
            CALL INITIALIZE_COH_INT_SEARCH
         END IF
! COHESION

! To do only des in the 1st time step only for a new run so the particles settle down
! before the coupling is turned on

         IF(RUN_TYPE == 'NEW') THEN
            IF(DES_CONTINUUM_COUPLED.AND.(.NOT.USE_COHESION)) THEN
               DES_CONTINUUM_COUPLED = .FALSE.
               DO FACTOR = 1, NFACTOR
                  IF (FACTOR .EQ. 1)&
                     WRITE(*,'(/,2X,A,/,4X,A,X,I,X,A)') &
                     'FIRST PASS IN DES_TIME_MARCH BEFORE COUPLING',&
                     'DEM settling period performed', NFACTOR, 'times'

                  ! Force calculation         
                  CALL CALC_FORCE_DES
                  
                  PC = 1
                  DO NP = 1, MAX_PIS
                     IF(PC .GT. PIS) EXIT
                     IF(.NOT.PEA(NP,1)) CYCLE
                     CALL CFNEWVALUES(NP)
                     PC = PC + 1
                  ENDDO

                  CALL PARTICLES_IN_CELL
                  NSN = NSN + 1    

                  ! Neighbor search                      
                  IF(MOD(FACTOR,INT(NEIGHBOR_SEARCH_N)).EQ.0) THEN 
                     CALL NEIGHBOUR
                  ELSEIF(DO_NSEARCH) THEN 
                     CALL NEIGHBOUR
                     DO_NSEARCH = .FALSE.
                     NSN = 0
                  ENDIF
               ENDDO
               DES_CONTINUUM_COUPLED = .TRUE.
               WRITE(*,'(2X,A,/)') 'END DEM settling period'
            ENDIF
            IF(DES_INTERP_ON) THEN 
               CALC_FC = .FALSE.
               CALLFROMDES = .FALSE.
            ENDIF
            CALL PARTICLES_IN_CELL
            CALL WRITE_DES_DATA
         ENDIF ! end if on new run type
      ENDIF    ! end if first pass


! get the right s_time from MFIX time in case of restarts
      S_TIME = TIME

      IF(DES_CONTINUUM_COUPLED) THEN
         IF(DT.GE.DTSOLID) THEN
            FACTOR = CEILING(real(DT/DTSOLID))
         ELSE
            FACTOR = 1
            DTSOLID_TMP = DTSOLID
            DTSOLID = DT
         ENDIF
         write(*,'(2X,A,X,I,X,A)') &
            'DEM SIMULATION will be called', &
            FACTOR, 'times this fluid time step' 
         write(*,'(4X,A,X,ES)') 'dt =', dt
         write(*,'(4X,A,X,ES)') 'dtsolid =', dtsolid
         write(*,'(4X,A,X,I)') 'int(dt/dtsolid) =', nint(dt/dtsolid)
      ELSE
! added TIME for restart & +1 removed
         FACTOR = CEILING(real((TSTOP-TIME)/DTSOLID)) 
         IF(RUN_TYPE .NE. 'NEW') THEN
! for restarts these two counters should not start from zero.
           PTC = DTSOLID 
           DESRESDT = DTSOLID
         ENDIF
         write(*,'(1X,A,X,I,X,A)') &
            'DEM SIMULATION will be called', FACTOR, 'times'
      ENDIF


      IF(NEIGHBOR_SEARCH_N.GT.FACTOR) THEN 
         !NEIGHBOR_SEARCH_N = FACTOR
         NSN = NEIGHBOR_SEARCH_N - 1
      ENDIF 


      DO NN = 1, FACTOR         !  do NN = 1, FACTOR

         IF(DES_CONTINUUM_COUPLED) THEN
            IF(S_TIME.GE.(TIME+DT)) EXIT

            IF((S_TIME+DTSOLID).GT.(TIME+DT)) THEN 
               TMP_DTS = DTSOLID
               DTSOLID = TIME + DT - S_TIME
            ENDIF 

            !PRINT *,"DES-MFIX COUPLED, ITER, TIMESTEP, TIME", NN, DTSOLID, S_TIME
            !PRINT *,"DES-MFIX COUPLED, ITER, TIME, DES_INTERP_ON ?", NN, S_TIME, DES_INTERP_ON
            CALC_FC = .TRUE.
            IF(DES_INTERP_ON .AND. NN.EQ.FACTOR) THEN 
               ! calculate the mean fields only at the last DEM time step
               CALLFROMDES = .FALSE.
            ELSE 
               CALLFROMDES = .TRUE.
            ENDIF
         ELSE
            IF(DEBUG_DES) PRINT *,"DES UNCOUPLED", NN, S_TIME 
         ENDIF 

         ! Force calculation         
         CALL CALC_FORCE_DES

         PC = 1
         DO NP = 1, MAX_PIS
            IF(PC .GT. PIS) EXIT
            IF(.NOT.PEA(NP,1)) CYCLE
            CALL CFNEWVALUES(NP)
            PC = PC + 1
         ENDDO

! For systems with inlets/outlets check to determine if a particle has
! fully entered or exited the domain.  If the former, remove the status
! of 'new' and if the latter, remove the particle.
         IF (DES_MI) CALL DES_CHECK_PARTICLE

         CALL PARTICLES_IN_CELL
         NSN = NSN + 1    
         
         IF(NN.EQ.1.OR.MOD(NN,INT(NEIGHBOR_SEARCH_N)).EQ.0) THEN 
            !WRITE(*,'(4X,A,I)') 'CALLING NEIGHBOR BECAUSE NN = ', NN
            CALL NEIGHBOUR
         ELSEIF(DO_NSEARCH) THEN 
            CALL NEIGHBOUR
            !PRINT*, 'CALLING NEIGHBOR BECASUE DO_NSEARCH = ', DO_NSEARCH
            DO_NSEARCH = .FALSE.
            NSN = 0
         ENDIF

! Update time before any write statements to reflect that the particle has
! already been moved for the current time step                  
         S_TIME = S_TIME + DTSOLID

! So that section of granular_temperature routine is only calculated at
! end of current DEM simulation         
         FLAGTEMP = .FALSE.
         IF(DES_CONTINUUM_COUPLED .AND. NN.EQ.FACTOR) FLAGTEMP = .TRUE.
         CALL DES_GRANULAR_TEMPERATURE(FLAGTEMP)
        
         IF(.NOT.DES_CONTINUUM_COUPLED) THEN    
! Keep track of TIME for DEM simulations
            TIME = S_TIME

! Write DES data for DEM only case 
            IF(PRINT_DES_DATA) THEN
                 PTC = PTC + DTSOLID
! Additional check was added to make sure DEM data are written at exactly NN = FACTOR.
               IF((PTC.GE.P_TIME .AND. NN .NE. (FACTOR-1)) .OR. NN == FACTOR) THEN
                  CALL DES_GRANULAR_TEMPERATURE(FLAGTEMP)
                  CALL WRITE_DES_DATA
                  WRITE(*,*) 'DES_SPX file written at Time= ', S_TIME
                  WRITE(UNIT_LOG,*) 'DES_SPX file written at Time= ', S_TIME
                  PTC = PTC - P_TIME ! this should not be set to zero but to the residual time difference.
               ENDIF
            ENDIF

! Write Restart for DEM only case
            DESRESDT = DESRESDT + DTSOLID
            IF((DESRESDT.GE.RES_DT .AND. NN .NE. (FACTOR-1)) &
            .OR. NN == FACTOR) THEN ! same as PTC
               CALL WRITE_DES_RESTART
! Write RES1 here instead of time_march, this will also keep track of TIME.
               CALL WRITE_RES1 
               WRITE(*,*) 'DES_RES file written at Time= ', S_TIME
               WRITE(UNIT_LOG,*) 'DES_RES file written at Time= ', S_TIME
               DESRESDT = DESRESDT - RES_DT
            ENDIF

         ENDIF  ! end if (.not.des_continuum_coupled)


! J.Musser : mass inlet/outlet -> particles entering the system
         IF(DES_MI)THEN 
            DO BCV_I = 1, SIZE(DES_BC_MI_ID)
               IF(PI_FACTOR(BCV_I) .GT. 1)THEN
                  IF(DES_MI_TIME(BCV_I) .LE. S_TIME) THEN   !Verify start time
                     CALL DES_MASS_INLET(BCV_I)
                     DES_MI_TIME(BCV_I) = S_TIME + PI_FACTOR(BCV_I) * DTSOLID
                  ENDIF
               ELSE
                  CALL DES_MASS_INLET(BCV_I)
               ENDIF
            ENDDO
         ENDIF

         IF (NN .EQ. FACTOR) &
            WRITE(*,'(4X,A,I5,2X,ES15.7)') &
               'MAX no. neigh & % overlap = ', NEIGH_MAX, OVERLAP_MAX

      ENDDO     ! end do NN = 1, FACTOR


! Write Restart
      IF(((TIME+DT+0.1d0*DT)>=DES_RES_TIME).OR.((TIME+DT+0.1d0*DT)>=TSTOP)) THEN
         CALL WRITE_DES_RESTART
! Write RES1 here instead of time_march
         CALL WRITE_RES1 
         WRITE(*,*) 'DES_RES file written at Time= ', TIME
         WRITE(UNIT_LOG,*) 'DES_RES file written at Time= ', TIME
         DES_RES_TIME = (INT((TIME+DT+0.1d0*DT)/RES_DT) + 1)*RES_DT
      ENDIF

      IF(DT.LT.DTSOLID_TMP) THEN
         DTSOLID = DTSOLID_TMP
      ENDIF

      IF(TMP_DTS.NE.ZERO) THEN
         DTSOLID = TMP_DTS
         TMP_DTS = ZERO
      ENDIF


      END SUBROUTINE DES_TIME_MARCH

