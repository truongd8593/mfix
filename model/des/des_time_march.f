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

      IMPLICIT NONE
!------------------------------------------------
! Local variables
!------------------------------------------------
     
      INTEGER NN, FACTOR, NP, IJK, I, J, K, BCV_I

!     Local variables to keep track of time when dem restart and des
!     write data need to be written when des_continuum_coupled is F
      DOUBLE PRECISION DES_RES_TIME, DES_SPX_TIME

!     Temporary variables when des_continuum_coupled is T to track
!     changes in solid time step 
      DOUBLE PRECISION TMP_DTS, DTSOLID_TMP 
      CHARACTER*5 FILENAME

!     Temporary variable used to track reporting frequency for the
!     maximum overlap and maximum no. neighbors for given NN loop
      DOUBLE PRECISION DES_TMP_TIME

!     Logical to see whether this is the first entry to this routine
      LOGICAL,SAVE:: FIRST_PASS = .TRUE.

!     index to track accounted for particles  
      INTEGER PC    

!-----------------------------------------------

      INCLUDE 'function.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'fun_avg2.inc'

! initialization      
      TMP_DTS = ZERO

      IF(FIRST_PASS) THEN 

      WRITE(*,'(1X,A)')&
         '---------- FIRST PASS DES_TIME_MARCH ---------->'
         S_TIME = ZERO
         FIRST_PASS = .FALSE.
         INQC = INIT_QUAD_COUNT

! When no particles are present, skip the startup routines that loop over particles.
! That is, account for a setup that starts with no particles in the system.
         IF(PARTICLES /= 0) THEN
            CALL NEIGHBOUR         

         
! COHESION INITIALIZE
            IF(USE_COHESION)THEN
               CALL INITIALIZE_COHESION_PARAMETERS
               CALL INITIALIZE_COH_INT_SEARCH
            ENDIF

! To do only des in the 1st time step only for a new run so the particles settle down
! before the coupling is turned on

            IF(RUN_TYPE == 'NEW') THEN
               IF(DES_CONTINUUM_COUPLED.AND.(.NOT.USE_COHESION)) THEN
                  DES_CONTINUUM_COUPLED = .FALSE.
                  DO FACTOR = 1, NFACTOR
                     IF (FACTOR .EQ. 1)&
                        WRITE(*,'(3X,A,/,5X,A,X,I,X,A)') &
                        'FIRST PASS in DES_TIME_MARCH for new runs',&
                        'DEM settling period performed', NFACTOR, &
                        'times'

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

                     ! Neighbor search                      
                     IF(MOD(FACTOR,NEIGHBOR_SEARCH_N).EQ.0) THEN 
                        CALL NEIGHBOUR
                     ELSEIF(DO_NSEARCH) THEN 
                        CALL NEIGHBOUR
                        DO_NSEARCH = .FALSE.
                     ENDIF
                  ENDDO
                  DES_CONTINUUM_COUPLED = .TRUE.
                  WRITE(*,'(3X,A)') 'END DEM settling period'
               ENDIF   ! end if coupled and no cohesion
               IF(DES_INTERP_ON) THEN 
                  CALC_FC = .FALSE.
                  CALLFROMDES = .FALSE.
               ENDIF
               CALL PARTICLES_IN_CELL
               CALL WRITE_DES_DATA
               WRITE(*,'(3X,A,X,ES)') &
                  'DES data file written at time =', S_TIME
               WRITE(UNIT_LOG,*) &
                  'DES data file written at time = ', S_TIME
            ENDIF   ! end if on new run type

            WRITE(*,'(1X,A)')&
               '<---------- END FIRST PASS DES_TIME_MARCH ----------'

         ENDIF   ! end if particles /= 0         
      ENDIF    ! end if first pass


! In case of restarts assign S_TIME from MFIX TIME 
      S_TIME = TIME

      IF(DES_CONTINUUM_COUPLED) THEN
         IF(DT.GE.DTSOLID) THEN
            FACTOR = CEILING(real(DT/DTSOLID))
         ELSE
            FACTOR = 1
            DTSOLID_TMP = DTSOLID
            DTSOLID = DT
         ENDIF
         WRITE(*,'(1X,A,X,I,X,A)') &
            'DEM SIMULATION will be called', &
            FACTOR, 'times this fluid time step' 
         WRITE(*,'(3X,2(A,X,ES15.7,2X))') 'dt =', dt,&
            'dtsolid =', dtsolid
         WRITE(*,'(3X,A,X,I)') 'int(dt/dtsolid) =', nint(dt/dtsolid)
      ELSE
         FACTOR = CEILING(real((TSTOP-TIME)/DTSOLID)) 
         WRITE(*,'(1X,A)')&
            '---------- START DES_TIME_MARCH ---------->'
         WRITE(*,'(3X,A,X,I,X,A)') &
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

      IF(NEIGHBOR_SEARCH_N.GT.FACTOR) THEN 
         !NEIGHBOR_SEARCH_N = FACTOR
      ENDIF 


      DO NN = 1, FACTOR         !  do NN = 1, FACTOR

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
               WRITE(*,'(3X,A,X,I,X,A,X,ES15.7)') &
                  'DES-COUPLED LOOP NO. =', NN, ' S_TIME =', S_TIME 
               WRITE(*,'(3X,A,X,ES15.7)') &
                  'DTSOLID =', DTSOLID
            ENDIF

            CALC_FC = .TRUE.
            IF(DES_INTERP_ON .AND. NN.EQ.FACTOR) THEN 
! Toggle flag so the mean fields are calculated only at the last DEM
! time step
               CALLFROMDES = .FALSE.
            ELSE 
               CALLFROMDES = .TRUE.
            ENDIF
         ELSE   ! else if (des_continuum_coupled)
            IF(DEBUG_DES) WRITE(*,'(3X,A,X,I,X,A,X,ES15.7)') &
               'DEM LOOP NO. =', NN, ' S_TIME =', S_TIME 
         ENDIF   ! end if/else (des_continuum_coupled) 
         

! If system is empty, skip force calculation calls
         IF (PIS /= 0) THEN
            CALL CALC_FORCE_DES

            PC = 1
            DO NP = 1, MAX_PIS
               IF(PC .GT. PIS) EXIT
               IF(.NOT.PEA(NP,1)) CYCLE
! Update particle position, velocity            
               CALL CFNEWVALUES(NP)
               PC = PC + 1
            ENDDO

! For systems with inlets/outlets check to determine if a particle has
! fully entered or exited the domain.  If the former, remove the status
! of 'new' and if the latter, remove the particle.
            IF (DES_MI) CALL DES_CHECK_PARTICLE

            CALL PARTICLES_IN_CELL

            IF(NN.EQ.1 .OR. MOD(NN,NEIGHBOR_SEARCH_N).EQ.0) THEN 
               IF(DEBUG_DES) WRITE(*,'(3X,A,A,/,5X,A,I)') &
                  'Calling NEIGHBOUR: NN=1 or ',&
                  'MOD(NN,NEIGHBOR_SEARCH_N)=0',&
                  'NEIGHBOR_SEARCH_N = ', NEIGHBOR_SEARCH_N
               CALL NEIGHBOUR
            ELSEIF(DO_NSEARCH) THEN 
               IF(DEBUG_DES) WRITE(*,'(3X,A,A,/,5X,A,A,L)') &
                  'Calling NEIGHBOUR: a particle moved ',&
                  'more than its radius since', 'last time NEIGHBOUR ',&
                  'was called; DO_NSEARCH = ', DO_NSEARCH
               CALL NEIGHBOUR
               DO_NSEARCH = .FALSE.
            ENDIF
         ENDIF   ! end if particles /= 0


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
                  WRITE(*,'(3X,A,X,ES)') &
                     'DES data file written at time =', S_TIME
                  WRITE(UNIT_LOG,*) &
                     'DES data file written at time = ', S_TIME
               ENDIF
            ENDIF

            IF ( (S_TIME+0.1d0*DTSOLID >= DES_RES_TIME) .OR. &
                 (S_TIME+0.1d0*DTSOLID >= TSTOP) .OR. &
                 (NN == FACTOR) ) THEN
               DES_RES_TIME = &
                  ( INT((S_TIME+0.1d0*DTSOLID)/DES_RES_DT) &
                  + 1 )*DES_RES_DT
                  CALL WRITE_DES_RESTART
! Write RES1 here since it won't be called in time_march.  This will
! also keep track of TIME
               CALL WRITE_RES1 
               WRITE(*,'(3X,A,X,ES)') &
                  'DES.RES and .RES files written at time =', S_TIME
               WRITE(UNIT_LOG,*) &
                  'DES.RES and .RES files written at time = ', S_TIME
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

         IF ( (S_TIME+0.1d0*DTSOLID >= DES_TMP_TIME) .OR. &
              ( (S_TIME+0.1d0*DTSOLID >= TSTOP) .AND. &
               (.NOT.DES_CONTINUUM_COUPLED) ) .OR. &          
              (NN .EQ. FACTOR) ) THEN
            DES_TMP_TIME = ( INT((S_TIME+0.1d0*DTSOLID)/DES_SPX_DT) &
               + 1 )*DES_RES_DT
            WRITE(*,'(3X,A,I,A,/,5X,A,X,I5,2X,A,X,ES15.7)') &
               'For loop ', NN, ' :',&
               'MAX no. of neighbors =',NEIGH_MAX,&
               'and MAX % overlap =', OVERLAP_MAX
         ENDIF

      ENDDO     ! end do NN = 1, FACTOR

! When coupled, and if needed, reset the discrete time step accordingly
      IF(DT.LT.DTSOLID_TMP) THEN
         DTSOLID = DTSOLID_TMP
      ENDIF

      IF(TMP_DTS.NE.ZERO) THEN
         DTSOLID = TMP_DTS
         TMP_DTS = ZERO
      ENDIF

     IF(.NOT.DES_CONTINUUM_COUPLED) WRITE(*,'(1X,A)')&
        '<---------- END DES_TIME_MARCH ----------'

      END SUBROUTINE DES_TIME_MARCH

