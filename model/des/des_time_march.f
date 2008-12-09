!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                         C  
!     Module name: DES_TIME_MARCH                                         C

!>
!>    Purpose: Called in model/time_march.f to do DES calcs
!>    Main DEM driver routine
!>

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
!     
      INTEGER NN, LN, I, J, K, FACTOR, NSN, NP, L, IJK, GTC_COUNT, GTC_FACTOR
      DOUBLE PRECISION TEMP_DTS, DTSOLIDTEMP 
      CHARACTER*5 FILENAME
!     Logical to see whether this is the first entry to this routine
      LOGICAL,SAVE:: FIRST_PASS = .TRUE.
      LOGICAL DES_CONTINUUM_COUPLED_FT
      DOUBLE PRECISION  pgrad_tmp(1:DIMN), GRAV_TMP(1:DIMN)
      DOUBLE PRECISION:: xmintmp, xmax, ymin, ymax, Ax, bx, ay, by, az, bz , mean_free_path, t_mfp
                                !     
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
         TEMP_DTS = ZERO
         PTC = ZERO
         INQC = INIT_QUAD_COUNT
         IF(RUN_TYPE == 'NEW') THEN
            DES_SPX_TIME =  TIME
            DES_RES_TIME =  TIME
         ELSE
            DES_SPX_TIME = (INT((TIME+DT+0.1d0*DT)/SPX_DT(1))+1)*SPX_DT(1)
            DES_RES_TIME = (INT((TIME+DT+0.1d0*DT)/RES_DT)   +1)*RES_DT
         ENDIF

!        PRINT *,'SPX TIME', SPX_DT(1), DES_SPX_TIME
         
         CALL NEIGHBOUR
         
!     !COHESION INITIALIZE
         IF(USE_COHESION)THEN
            CALL INITIALIZE_COHESION_PARAMETERS
            CALL INITIALIZE_COH_INT_SEARCH
         END IF
!     !COHESION

!     To do only des in the 1st time step only for a new run so the particles settle down
!     before the coupling is turned on

         IF(RUN_TYPE == 'NEW') THEN
            IF(DES_CONTINUUM_COUPLED.AND.(.NOT.USE_COHESION)) THEN
               DES_CONTINUUM_COUPLED = .FALSE.
               DO FACTOR = 1, NFACTOR
                  PRINT *,'DES FIRST PASS IN TIME MARCH', FACTOR, DTSOLID, INIT_QUAD_COUNT
!     New values
                                       
                  !     Force calculation         
                  CALL CALC_FORCE_DES
                  
                  DO LN = 1, PARTICLES
                     CALL CFNEWVALUES(LN)
                  END DO
                  !     Neighbor search     
                  CALL PARTICLES_IN_CELL
                  
                  NSN = NSN + 1    
                  
                  IF(MOD(FACTOR,INT(NEIGHBOR_SEARCH_N)).EQ.0) THEN 
                     CALL NEIGHBOUR
                     
                  ELSE IF(DO_NSEARCH) THEN 
                     CALL NEIGHBOUR
                     DO_NSEARCH = .FALSE.
                     NSN = 0
                  END IF
               END DO
               DES_CONTINUUM_COUPLED = .TRUE.
            END IF
            IF(DES_INTERP_ON) THEN 
               CALC_FC = .FALSE.
               CALLFROMDES = .FALSE.
            end IF
            CALL PARTICLES_IN_CELL
            
            CALL WRITE_DES_DATA
         END IF ! end if on new run type

      END IF
      
      IF(DES_CONTINUUM_COUPLED) THEN
         S_TIME = TIME
         IF(DT.GE.DTSOLID) THEN
            FACTOR = CEILING(real(DT/DTSOLID)) + 1
         ELSE
            FACTOR = 1
            DTSOLIDTEMP = DTSOLID
            DTSOLID = DT
            PRINT *,"DT_SOLID greater than DT_FLUID. DEM is called once"
         END IF
      ELSE
         write(*,*) ' dt, dtsolid = ', dt, dtsolid, nint(dt/dtsolid)
         FACTOR = CEILING(real((TSTOP-TIME)/DTSOLID)) ! added TIME for restart & +1 removed
         S_TIME = TIME  ! get the right s_time from MFIX time in case of DEM restarts
	 IF(RUN_TYPE .NE. 'NEW') THEN
	   PTC = DTSOLID ! for restarts these two counters shoud not start from zero.
	   DESRESDT = DTSOLID
	 ENDIF
      END IF
      
      PRINT *,"Discrete Element Simulation is being called"&
             , FACTOR," times.", dt, dtsolid

      IF(NEIGHBOR_SEARCH_N.GT.FACTOR) THEN 
         !NEIGHBOR_SEARCH_N = FACTOR
         
         NSN = NEIGHBOR_SEARCH_N - 1
      ENDIF 

      DO NN = 1, FACTOR         !  do NN = 1, FACTOR

         IF(DES_CONTINUUM_COUPLED) THEN
!           PRINT *,"DES-MFIX COUPLED, ITER, TIMESTEP, TIME", NN, DTSOLID, S_TIME
            !PRINT *,"DES-MFIX Co , ITER, TIME, DES_INTERP_ON ?", NN, S_TIME, DES_INTERP_ON
            IF(DES_INTERP_ON) THEN 
               CALC_FC = .TRUE.
               IF(NN.EQ.FACTOR) THEN
                  !calculate the mean fields only at the last DEM time step
                  CALLFROMDES = .FALSE.
               ELSE 
                  CALLFROMDES = .TRUE.
               ENDIF
            ELSE
               CALC_FC = .TRUE.
               CALLFROMDES = .TRUE.
            ENDIF

         ELSE
            PRINT *,"DES UNCOUPLED", NN, S_TIME

         END IF 
         

!     Force calculation         
         CALL CALC_FORCE_DES
         
         DO LN = 1, PARTICLES
            CALL CFNEWVALUES(LN)
         END DO

         CALL PARTICLES_IN_CELL

         NSN = NSN + 1    
         
                                !     New values
         IF(NN.EQ.1.OR.MOD(NN,INT(NEIGHBOR_SEARCH_N)).EQ.0) THEN 
            
            !PRINT*, 'CALLING NEIGHBOR BECASUE NN = ', NN
            CALL NEIGHBOUR
         ELSE IF(DO_NSEARCH) THEN 
            CALL NEIGHBOUR
            
            !PRINT*, 'CALLING NEIGHBOR BECASUE DO_NSEARCH = ', DO_NSEARCH
            DO_NSEARCH = .FALSE.
            NSN = 0
         END IF


         
         IF(PRINT_DES_DATA) THEN    
            IF(.NOT.DES_CONTINUUM_COUPLED) THEN
               PTC = PTC + DTSOLID
         ! Additional check was added to make sure DEM data are written at exactly NN = FACTOR.
	       IF((PTC.GE.P_TIME .AND. NN .NE. (FACTOR-1)) .OR. NN == FACTOR) THEN
                  
                  CALL DES_GRANULAR_TEMPERATURE
                  CALL WRITE_DES_DATA
         
                  WRITE(*,*) 'DES_SPX file written at Time= ', S_TIME
                  WRITE(UNIT_LOG,*) 'DES_SPX file written at Time= ', S_TIME
                  PTC = PTC - P_TIME ! this should not be set to zero but to the residual time difference.
               END IF
            END IF
         END IF

!     Write Restart for DEM only case
         IF(.NOT.DES_CONTINUUM_COUPLED) THEN
            DESRESDT = DESRESDT + DTSOLID
            IF((DESRESDT.GE.RES_DT .AND. NN .NE. (FACTOR-1)) .OR. NN == FACTOR) THEN ! same as PTC
               CALL WRITE_DES_RESTART
	       CALL WRITE_RES1 ! write RES1 here instead of time_march, this will also keep track of TIME.
               WRITE(*,*) 'DES_RES file written at Time= ', S_TIME

               WRITE(UNIT_LOG,*) 'DES_RES file written at Time= ', S_TIME
               DESRESDT = DESRESDT - RES_DT
            END IF
         END IF

         IF(DES_CONTINUUM_COUPLED.AND.S_TIME.GE.(TIME+DT)) EXIT

         IF(DES_CONTINUUM_COUPLED) THEN
            IF((S_TIME+DTSOLID).GT.(TIME+DT)) THEN 
               TEMP_DTS = DTSOLID
               DTSOLID = TIME + DT - S_TIME
            END IF
            S_TIME = S_TIME + DTSOLID
         ELSE
            S_TIME = S_TIME + DTSOLID
	    TIME = S_TIME ! keep track of TIME for DEM simulations.
         END IF 
         IF(DES_CONTINUUM_COUPLED.AND.NN.eq.factor) CALL DES_GRANULAR_TEMPERATURE
      END DO                    ! end do NN = 1, FACTOR

      WRITE(*,*) 'max neigh = ',NEIGH_MAX, OVERLAP_MAX
      
!     Write Restart
      IF(((TIME+DT+0.1d0*DT)>=DES_RES_TIME).OR.((TIME+DT+0.1d0*DT)>=TSTOP)) THEN
         CALL WRITE_DES_RESTART
	 CALL WRITE_RES1 ! write RES1 here instead of time_march
         WRITE(*,*) 'DES_RES file written at Time= ', TIME
         WRITE(UNIT_LOG,*) 'DES_RES file written at Time= ', TIME
         DES_RES_TIME = (INT((TIME+DT+0.1d0*DT)/RES_DT) + 1)*RES_DT
      END IF

100   continue
      IF(.NOT.DES_CONTINUUM_COUPLED) TSTOP = DT

      IF(DT.LT.DTSOLIDTEMP) THEN
         DTSOLID = DTSOLIDTEMP
      END IF
      !write(*,*) 'loop end dt, dtsolid = ', dt, dtsolid
      IF(TEMP_DTS.NE.ZERO) THEN
         DTSOLID = TEMP_DTS
         TEMP_DTS = ZERO
      END IF

      !write(*,*) 'end dt, dtsolid = ', dt, dtsolid      
 3020 FORMAT(I5.5)

      END SUBROUTINE DES_TIME_MARCH
