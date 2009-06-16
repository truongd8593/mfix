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

      IMPLICIT NONE
!-----------------------------------------------
!     L o c a l   V a r i a b l e s
!-----------------------------------------------
     
      INTEGER NN, FACTOR, NSN, NP, IJK, I, J, K

!     Temporary variables when des_continuum_coupled is F to track
!     reporting time 
      DOUBLE PRECISION PTC, DESRESDT 

!     Temporary variables when des_continuum_coupled is T to track
!     changes in solid time step 
      DOUBLE PRECISION TMP_DTS, DTSOLID_TMP 
      CHARACTER*5 FILENAME

!     Logical to see whether this is the first entry to this routine
      LOGICAL,SAVE:: FIRST_PASS = .TRUE.


      INCLUDE 'function.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'fun_avg2.inc'
      
      
      !write(*,*) ' dt, dtsolid = ', dt, dtsolid
      !read(*,*)
      
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
                  PRINT *,'DES FIRST PASS IN TIME MARCH', FACTOR, DTSOLID, INIT_QUAD_COUNT

                  ! Force calculation         
                  CALL CALC_FORCE_DES
                  
                  DO NP = 1, PARTICLES
                     CALL CFNEWVALUES(NP)
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
            PRINT *,"DT_SOLID greater than DT_FLUID. DEM is called once"
         ENDIF
      ELSE
! added TIME for restart & +1 removed
         FACTOR = CEILING(real((TSTOP-TIME)/DTSOLID)) 
         IF(RUN_TYPE .NE. 'NEW') THEN
! for restarts these two counters should not start from zero.
           PTC = DTSOLID 
           DESRESDT = DTSOLID
         ENDIF
      ENDIF

      write(*,*) '              dt = ', dt
      write(*,*) '         dtsolid = ', dtsolid
      write(*,*),' int(dt/dtsolid) = ', nint(dt/dtsolid)      
      PRINT *,"Discrete Element Simulation is being called"&
             , FACTOR," times"

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
            !PRINT *,"DES-MFIX Co , ITER, TIME, DES_INTERP_ON ?", NN, S_TIME, DES_INTERP_ON
            CALC_FC = .TRUE.
            IF(DES_INTERP_ON .AND. NN.EQ.FACTOR) THEN 
               ! calculate the mean fields only at the last DEM time step
               CALLFROMDES = .FALSE.
            ELSE 
               CALLFROMDES = .TRUE.
            ENDIF
         ELSE
            PRINT *,"DES UNCOUPLED", NN, S_TIME
         ENDIF 
         
         ! Force calculation         
         CALL CALC_FORCE_DES

         DO NP = 1, PARTICLES
            CALL CFNEWVALUES(NP)
         ENDDO

         CALL PARTICLES_IN_CELL
         NSN = NSN + 1    
         
         IF(NN.EQ.1.OR.MOD(NN,INT(NEIGHBOR_SEARCH_N)).EQ.0) THEN 
            !PRINT*, 'CALLING NEIGHBOR BECASUE NN = ', NN
            CALL NEIGHBOUR
         ELSEIF(DO_NSEARCH) THEN 
            CALL NEIGHBOUR
            !PRINT*, 'CALLING NEIGHBOR BECASUE DO_NSEARCH = ', DO_NSEARCH
            DO_NSEARCH = .FALSE.
            NSN = 0
         ENDIF

         IF(DES_CONTINUUM_COUPLED) THEN
            S_TIME = S_TIME + DTSOLID
            IF(NN.EQ.FACTOR) CALL DES_GRANULAR_TEMPERATURE
         ENDIF


         IF(.NOT.DES_CONTINUUM_COUPLED) THEN    
! Update time before any write statements to reflect that the particle has
! already been moved for the current time step                  
            S_TIME = S_TIME + DTSOLID
! Keep track of TIME for DEM simulations
            TIME = S_TIME

! Write DES data for DEM only case 
            IF(PRINT_DES_DATA) THEN
                 PTC = PTC + DTSOLID
! Additional check was added to make sure DEM data are written at exactly NN = FACTOR.
               IF((PTC.GE.P_TIME .AND. NN .NE. (FACTOR-1)) .OR. NN == FACTOR) THEN
                  CALL DES_GRANULAR_TEMPERATURE
                  CALL WRITE_DES_DATA
                  WRITE(*,*) 'DES_SPX file written at Time= ', S_TIME
                  WRITE(UNIT_LOG,*) 'DES_SPX file written at Time= ', S_TIME
                  PTC = PTC - P_TIME ! this should not be set to zero but to the residual time difference.
               ENDIF
            ENDIF

! Write Restart for DEM only case
            DESRESDT = DESRESDT + DTSOLID
            IF((DESRESDT.GE.RES_DT .AND. NN .NE. (FACTOR-1)) .OR. NN == FACTOR) THEN ! same as PTC
               CALL WRITE_DES_RESTART
! Write RES1 here instead of time_march, this will also keep track of TIME.
               CALL WRITE_RES1 
               WRITE(*,*) 'DES_RES file written at Time= ', S_TIME
               WRITE(UNIT_LOG,*) 'DES_RES file written at Time= ', S_TIME
               DESRESDT = DESRESDT - RES_DT
            ENDIF
         ENDIF  ! end if (.not.des_continuum_coupled)


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

      WRITE(*,*) 'MAX neigh & overlap = ', NEIGH_MAX, OVERLAP_MAX
      !write(*,*) 'loop end dt, dtsolid = ', dt, dtsolid     


      END SUBROUTINE DES_TIME_MARCH
