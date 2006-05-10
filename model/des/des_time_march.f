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
!  M o d u l e s
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
      INTEGER NN, LNN, I, J, K, FACTOR, NSN
      DOUBLE PRECISION TEMP_DTS, PTC 
      LOGICAL ALREADY_EXISTS

      IF(TIME.EQ.ZERO) THEN 

        NSN = 0
        S_TIME = TIME
        TEMP_DTS = ZERO
        PTC = ZERO

        CALL CFASSIGN(PARTICLES)

!!COHESION INITIALIZE
         IF(USE_COHESION)THEN
            CALL INITIALIZE_COHESION_PARAMETERS
            CALL INITIALIZE_COH_INT_SEARCH
         END IF
!!COHESION

        IF(DES_CONTINUUM_COUPLED.AND.(.NOT.USE_COHESION)) THEN

           DES_CONTINUUM_COUPLED = .FALSE.

           DO FACTOR = 1, 500
             PRINT *,'DES', FACTOR 
! New values
             DO LNN = 1, PARTICLES
                CALL CFNEWVALUES(LNN)
             END DO
! Neighbor search     
             NSN = NSN + 1    
             IF(NSN.EQ.NEIGHBOR_SEARCH_N) THEN 
                CALL NEIGHBOUR(PARTICLES)
                NSN = 0
             END IF
! Force calculation         
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
           GO TO 100
        END IF
      END IF

      IF(DES_CONTINUUM_COUPLED) THEN
         IF(DT.GE.DTSOLID) THEN
            FACTOR = DT/DTSOLID + 1
         ELSE
            PRINT *,'DT_SOLID greater than DT_FLUID. DEM not called'
         END IF
      ELSE
         FACTOR = TSTOP/DTSOLID + 1 
      END IF
      
      PRINT *,'Discrete Element Simulation is being called', FACTOR,' times.'

      IF(NEIGHBOR_SEARCH_N.GT.FACTOR) NEIGHBOR_SEARCH_N = FACTOR

      DO NN = 1, FACTOR ! do1

         IF(DES_CONTINUUM_COUPLED) THEN
            PRINT *,'DES-MFIX COUPLED', NN, S_TIME
         ELSE
            PRINT *,'DES', NN, S_TIME
         END IF 

! New values
         DO LNN = 1, PARTICLES
            CALL CFNEWVALUES(LNN)
         END DO
         
! Neighbor search     
         NSN = NSN + 1    
         IF(NSN.EQ.NEIGHBOR_SEARCH_N) THEN 
           CALL NEIGHBOUR(PARTICLES)
           NSN = 0
         END IF

! Force calculation         
         IF(DES_PERIODIC_WALLS) THEN
            CALL PERIODIC_WALL_CALC_FORCE_DES
         ELSE IF(INLET_OUTLET) THEN
            CALL DES_INLET_OUTLET
         ELSE
            CALL CALC_FORCE_DES
         END IF
  
!         CALL DES_GRANULAR_TEMPERATURE(NN, FACTOR)

         PTC = PTC + DTSOLID

             IF(PTC.GE.P_TIME) THEN 
               OPEN (UNIT=99, FILE='des_2-particles.out', STATUS='REPLACE')
               WRITE (99,*) real(s_time),real(DES_POS_NEW(1,1)),real(DES_POS_NEW(1,2)),  real(DES_VEL_NEW(2,1)),real(DES_VEl_NEW(2,2))

               IF(S_TIME.LE.0.2*TSTOP) THEN
                  OPEN (UNIT=9, FILE='des_all-particles-1.out', STATUS='REPLACE')
                  WRITE (9,*)' '
                  WRITE (9,*) 'Time=',S_TIME,'s'
                  DO LNN = 1, PARTICLES
                     WRITE (9,*) (DES_POS_NEW(LNN,K),K=1,DIMN),&
                     (DES_VEL_NEW(LNN,K),K=1,DIMN), DES_RADIUS(LNN), RO_Sol(LNN)
                  END DO
               ELSE IF(S_TIME.LE.0.4*TSTOP) THEN
                  OPEN (UNIT=9, FILE='des_all-particles-2.out', STATUS='REPLACE')
                  WRITE (9,*)' '
                  WRITE (9,*) 'Time=',S_TIME,'s'
                  DO LNN = 1, PARTICLES
                     WRITE (9,*) (DES_POS_NEW(LNN,K),K=1,DIMN),&
                     (DES_VEL_NEW(LNN,K),K=1,DIMN), DES_RADIUS(LNN), RO_Sol(LNN)
                  END DO
               ELSE IF(S_TIME.LE.0.6*TSTOP) THEN
                  OPEN (UNIT=9, FILE='des_all-particles-3.out', STATUS='REPLACE')
                  WRITE (9,*)' '
                  WRITE (9,*) 'Time=',S_TIME,'s'
                  DO LNN = 1, PARTICLES
                     WRITE (9,*) (DES_POS_NEW(LNN,K),K=1,DIMN),&
                     (DES_VEL_NEW(LNN,K),K=1,DIMN), DES_RADIUS(LNN), RO_Sol(LNN)
                  END DO
               ELSE IF(S_TIME.LE.0.8*TSTOP) THEN
                  OPEN (UNIT=9, FILE='des_all-particles-4.out', STATUS='REPLACE')
                  WRITE (9,*)' '
                  WRITE (9,*) 'Time=',S_TIME,'s'
                  DO LNN = 1, PARTICLES
                     WRITE (9,*) (DES_POS_NEW(LNN,K),K=1,DIMN),&
                     (DES_VEL_NEW(LNN,K),K=1,DIMN), DES_RADIUS(LNN), RO_Sol(LNN)
                  END DO
               ELSE IF(S_TIME.LE.TIME) THEN
                  OPEN (UNIT=9, FILE='des_all-particles-5.out', STATUS='REPLACE')
                  WRITE (9,*)' '
                  WRITE (9,*) 'Time=',S_TIME,'s'
                  DO LNN = 1, PARTICLES
                     WRITE (9,*) (DES_POS_NEW(LNN,K),K=1,DIMN),&
                     (DES_VEL_NEW(LNN,K),K=1,DIMN), DES_RADIUS(LNN), RO_Sol(LNN)
                  END DO
               END IF

               PTC = ZERO
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

         END DO ! enddo1

         IF(.NOT.DES_CONTINUUM_COUPLED) STOP
         
         IF(TEMP_DTS.NE.ZERO) THEN
            DTSOLID = TEMP_DTS
            TEMP_DTS = ZERO
         END IF

 100   CONTINUE

      END SUBROUTINE DES_TIME_MARCH
