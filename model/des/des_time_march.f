!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                         C  
!     Module name: DES_TIME_MARCH                                         C
!     Purpose: Called in time_march.f to do DES calcs                     C
!                                                                         C
!     Author: Jay Boyalakuntla                           Date: 21-Jun-04  C
!     Reviewer: Sreekanth Pannala                        Date: 09-Nov-06  C
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
      
      IMPLICIT NONE
!-----------------------------------------------
!     G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!     L o c a l   V a r i a b l e s
!-----------------------------------------------
!     
      INTEGER NN, LN, I, J, K, FACTOR, NSN
      DOUBLE PRECISION TEMP_DTS, DTSOLIDTEMP 
      CHARACTER*5 FILENAME
!     Logical to see whether this is the first entry to this routine
      LOGICAL,SAVE:: FIRST_PASS = .TRUE.

      IF(FIRST_PASS) THEN 
         
         FIRST_PASS = .FALSE.
         DESRESDT = 0.0d0
         NSN = 0
         S_TIME = ZERO
         TEMP_DTS = ZERO
         PTC = ZERO
         INQC = INIT_QUAD_COUNT

!        IF(RUN_TYPE == 'NEW') THEN
!           DES_SPX_TIME =  TIME
!           DES_RES_TIME =  TIME
! ELSE
!           DES_SPX_TIME = (INT((TIME+DT+0.1d0*DT)/SPX_DT(1))+1)*SPX_DT(1)
!           DES_RES_TIME = (INT((TIME+DT+0.1d0*DT)/RES_DT)   +1)*RES_DT
!        ENDIF

!        PRINT *,'SPX TIME', SPX_DT(1), DES_SPX_TIME

         CALL CFASSIGN

         CALL PARTICLES_IN_CELL

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
                  PRINT *,'DES', FACTOR 
!     New values
                  DO LN = 1, PARTICLES
                     CALL CFNEWVALUES(LN)
                  END DO
!     Neighbor search     
                  NSN = NSN + 1    
                  IF(DO_NSEARCH.OR.(NSN.EQ.NEIGHBOR_SEARCH_N)) THEN 
                     CALL NEIGHBOUR
                     DO_NSEARCH = .FALSE.
                     NSN = 0
                  END IF
!     Force calculation         
                  IF(DES_PERIODIC_WALLS) THEN
                     CALL PERIODIC_WALL_CALC_FORCE_DES
                  ELSE IF(INLET_OUTLET) THEN
                     CALL DES_INLET_OUTLET
                  ELSE
                     CALL CALC_FORCE_DES
                  END IF 
               END DO
               DES_CONTINUUM_COUPLED = .TRUE.
               CALL PARTICLES_IN_CELL
            END IF
         END IF ! end if on new run type

      END IF

      IF(DES_CONTINUUM_COUPLED) THEN
         NSN = NEIGHBOR_SEARCH_N - 1
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
         FACTOR = CEILING(real(TSTOP/DTSOLID)) + 1
      END IF
      
      PRINT *,"Discrete Element Simulation is being called"&
             , FACTOR," times."

      IF(NEIGHBOR_SEARCH_N.GT.FACTOR) NEIGHBOR_SEARCH_N = FACTOR

      DO NN = 1, FACTOR         !  do NN = 1, FACTOR

         IF(DES_CONTINUUM_COUPLED) THEN
!           PRINT *,"DES-MFIX COUPLED, ITER, TIMESTEP, TIME", NN, DTSOLID, S_TIME
            PRINT *,"DES-MFIX COUPLED, ITER, TIME", NN, S_TIME
         ELSE
            PRINT *,"DES", NN, S_TIME
         END IF 

!     New values
         DO LN = 1, PARTICLES
            CALL CFNEWVALUES(LN)
!           write(*,*) ln, des_pos_new(ln,:), des_vel_new(ln,:)
         END DO
         
!     Neighbor search     
         NSN = NSN + 1    
         IF(DO_NSEARCH.OR.(NSN.EQ.NEIGHBOR_SEARCH_N)) THEN 
            CALL NEIGHBOUR
            DO_NSEARCH = .FALSE.
            NSN = 0
         END IF

!     Force calculation         
         IF(DES_PERIODIC_WALLS) THEN
            CALL PERIODIC_WALL_CALC_FORCE_DES
         ELSE IF(INLET_OUTLET) THEN
            CALL DES_INLET_OUTLET
         ELSE
            CALL CALC_FORCE_DES
         END IF
         
!     CALL DES_GRANULAR_TEMPERATURE(NN, FACTOR)
         
         IF(PRINT_DES_DATA) THEN    
            IF(.NOT.DES_CONTINUUM_COUPLED) THEN
               PTC = PTC + DTSOLID
               IF(PTC.GE.P_TIME) THEN 
!                 WRITE (FILENAME, 3020) IFI
!                 OPEN(UNIT=99, FILE=FILENAME, STATUS='NEW')
                  CALL WRITE_DES_DATA
                  WRITE(*,*) 'DES_SPX file written at Time= ', S_TIME
                  WRITE(UNIT_LOG,*) 'DES_SPX file written at Time= ', S_TIME
                  CLOSE(99)
!                 IFI = IFI + 1
                  PTC = ZERO
               END IF
            END IF
         END IF

!     Write Restart for DEM only case
         IF(.NOT.DES_CONTINUUM_COUPLED) THEN
            DESRESDT = DESRESDT + DTSOLID
            IF(DESRESDT.GE.RES_DT) THEN
               CALL WRITE_DES_RESTART
               WRITE(*,*) 'DES_RES file written at Time= ', S_TIME
               WRITE(UNIT_LOG,*) 'DES_RES file written at Time= ', S_TIME
               DESRESDT = 0.0d0
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
         END IF 

      END DO                    ! end do NN = 1, FACTOR
      
!     IF(PRINT_DES_DATA) THEN
!        IF(((TIME+DT+0.1*DT).GE.DES_SPX_TIME) .OR. ((TIME+DT+0.1*DT).GE.TSTOP)) THEN
!           WRITE (FILENAME, 3020) IFI
!           OPEN(UNIT=99, FILE=TRIM(RUN_NAME)//'_DES_'//FILENAME//'.vtp', STATUS='NEW')
!           CALL WRITE_DES_DATA(99)
!           CLOSE(99)
!           DES_SPX_TIME =  (INT((TIME+DT+0.1*DT)/SPX_DT(1))+1)*SPX_DT(1)
!           IFI = IFI + 1
!           WRITE(*,*) 'DES_SPX file written at Time= ', TIME+DT
!           WRITE(UNIT_LOG,*) 'DES_SPX file written at Time= ', TIME+DT
!        END IF
!     END IF
      
!     Write Restart
!     IF(((TIME+DT+0.1d0*DT)>=DES_RES_TIME).OR.((TIME+DT+0.1d0*DT)>=TSTOP)) THEN
!        CALL WRITE_DES_RESTART
!        WRITE(*,*) 'DES_RES file written at Time= ', TIME+DT
!        WRITE(*,*) 'DES_RES Debug', TIME+1.1*DT, DES_RES_TIME, RES_DT
!        WRITE(UNIT_LOG,*) 'DES_RES file written at Time= ', TIME+DT
!        DES_RES_TIME = (INT((TIME+DT+0.1d0*DT)/RES_DT) + 1)*RES_DT
!     END IF

      IF(.NOT.DES_CONTINUUM_COUPLED) TSTOP = DT

      IF(DT.LT.DTSOLIDTEMP) THEN
         DTSOLID = DTSOLIDTEMP
      END IF

      IF(TEMP_DTS.NE.ZERO) THEN
         DTSOLID = TEMP_DTS
         TEMP_DTS = ZERO
      END IF

 3020 FORMAT(I5.5)

      END SUBROUTINE DES_TIME_MARCH
