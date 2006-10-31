!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                         C  
!     Module name: DES_TIME_MARCH                                         C
!     Purpose: Called in time_march.f to do DES calcs                     C
!                                                                         C
!     Author: Jay Boyalakuntla                           Date: 21-Jun-04  C
!     Reviewer:                                          Date:            C
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
!     Time at which REAL restart file is to be written
      LOGICAL ALREADY_EXISTS
      CHARACTER*20 FILENAME

      IF(TIME.EQ.ZERO) THEN 

         DESRESDT = 0.0
         NSN = 0
         PTC = 0D0 
         S_TIME = TIME
         TEMP_DTS = ZERO
         PTC = ZERO
         INQC = INIT_QUAD_COUNT

         DES_SPX_TIME =  (INT((TIME + 0.1*DT)/SPX_DT(1))+1)*SPX_DT(1)
         PRINT *,'SPX TIME', SPX_DT(1), DES_SPX_TIME
         CALL CFASSIGN(PARTICLES)

!     !COHESION INITIALIZE
         IF(USE_COHESION)THEN
            CALL INITIALIZE_COHESION_PARAMETERS
            CALL INITIALIZE_COH_INT_SEARCH
         END IF
!     !COHESION

!     To do only des in the 1st time step so the particles settle
!     and then start coupling

         IF(DES_CONTINUUM_COUPLED.AND.(.NOT.USE_COHESION)) THEN
            DES_CONTINUUM_COUPLED = .FALSE.
            DO FACTOR = 1, 500
               PRINT *,'DES', FACTOR 
!     New values
               DO LN = 1, PARTICLES
                  CALL CFNEWVALUES(LN)
               END DO
!     Neighbor search     
               NSN = NSN + 1    
               IF(DO_NSEARCH.OR.(NSN.EQ.NEIGHBOR_SEARCH_N)) THEN 
                  CALL NEIGHBOUR(PARTICLES)
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
            CALL PARTICLES_IN_CELL(PARTICLES)
            RETURN              ! exit subroutine
         END IF
      END IF

      IF(DES_CONTINUUM_COUPLED) THEN
         NSN = NEIGHBOR_SEARCH_N - 1
         S_TIME = TIME
         IF(DT.GE.DTSOLID) THEN
            FACTOR = CEILING(real(DT/DTSOLID)) 
         ELSE
            FACTOR = 1
            DTSOLIDTEMP = DTSOLID
            DTSOLID = DT
            PRINT *,"DT_SOLID greater than DT_FLUID. DEM is called once"
         END IF
      ELSE
         FACTOR = CEILING(real(TSTOP/DTSOLID))  
      END IF
      
      PRINT *,"Discrete Element Simulation is being called", FACTOR," times."

      IF(NEIGHBOR_SEARCH_N.GT.FACTOR) NEIGHBOR_SEARCH_N = FACTOR

      DO NN = 1, FACTOR         !  do NN = 1, FACTOR
         
         IF(DES_CONTINUUM_COUPLED) THEN
            PRINT *,"DES-MFIX COUPLED", NN, S_TIME
         ELSE
            PRINT *,"DES", NN, S_TIME
         END IF 

!     New values
         DO LN = 1, PARTICLES
            CALL CFNEWVALUES(LN)
         END DO
         
!     Neighbor search     
         NSN = NSN + 1    
         IF(DO_NSEARCH.OR.(NSN.EQ.NEIGHBOR_SEARCH_N)) THEN 
            CALL NEIGHBOUR(PARTICLES)
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
            
!     Write Restart
            DESRESDT = DESRESDT + DTSOLID
            IF(DESRESDT.EQ.RES_DT) THEN
               CALL WRITE_DES_RESTART(PARTICLES)
               DESRESDT = 0.0d0
            END IF

!     CALL DES_GRANULAR_TEMPERATURE(NN, FACTOR)
            
            IF(PRINT_DES_DATA) THEN    
               IF(.NOT.DES_CONTINUUM_COUPLED) THEN
                  PTC = PTC + DTSOLID
                  IF(PTC.GE.P_TIME) THEN 
                     WRITE (FILENAME, 3020) IFI
                     OPEN(UNIT=99, FILE=FILENAME, STATUS='NEW')
                     CALL WRITE_DES_DATA(99)
                     CLOSE(99)
                     IFI = IFI + 1
                     PTC = ZERO
                  END IF
               END IF
            END IF

            IF(DES_CONTINUUM_COUPLED) THEN
               IF((S_TIME+DTSOLID).GT.(TIME+DT)) THEN 
                  TEMP_DTS = DTSOLID
                  DTSOLID = TIME + DT - S_TIME
               END IF
               S_TIME = S_TIME + DTSOLID
            ELSE
               S_TIME = S_TIME + DTSOLID
            END IF 

         END DO                 ! end do NN = 1, FACTOR
         
         IF(PRINT_DES_DATA) THEN
            IF(((TIME+0.1*DT).GE.DES_SPX_TIME) .OR. ((TIME+DT).GE.TSTOP)) THEN
               WRITE (FILENAME, 3020) IFI
               OPEN(UNIT=99, FILE=FILENAME, STATUS='NEW')
               CALL WRITE_DES_DATA(99)
               CLOSE(99)
               DES_SPX_TIME =  (INT((TIME + 0.1*DT)/SPX_DT(1))+1)*SPX_DT(1)
               IFI = IFI + 1
            END IF
         END IF
         
         IF(.NOT.DES_CONTINUUM_COUPLED) TSTOP = DT

         IF(DT.LT.DTSOLIDTEMP) THEN
            DTSOLID = DTSOLIDTEMP
         END IF

         IF(TEMP_DTS.NE.ZERO) THEN
            DTSOLID = TEMP_DTS
            TEMP_DTS = ZERO
         END IF

 3020    FORMAT('DES_DATA',I4.4,'.vtp')

         END SUBROUTINE DES_TIME_MARCH
