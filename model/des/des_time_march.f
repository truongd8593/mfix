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
!     L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: CALLED_MAX = 50000000  
!-----------------------------------------------
!     L o c a l   V a r i a b l e s
!-----------------------------------------------
!     
      INTEGER          K, PC, PC1, NN, LNN, I, J
      LOGICAL ALREADY_EXISTS


!     
!      CALLED = 0

!      PRINT *, 'COUPLED_FLOW',COUPLED_FLOW

      DTSOLID = DT/FACTOR

      IF(COUPLED_FLOW) THEN            
         IF(TIME.GE.(DTSOLID)) THEN
            FACTOR = 20 
            DTSOLID = DT/FACTOR
            DES_CONTINUUM_COUPLED = .TRUE.
         END IF 
      END IF
      
      PRINT *,'Discrete Element Simulation is called', FACTOR,' times.'

      DO NN = 1, FACTOR

!!COHESION
         IF(CALLED.eq.0)THEN
            CALL INITIALIZE_COHESION_PARAMETERS
            IF(USE_COHESION)THEN
                CALL INITIALIZE_COH_INT_SEARCH
            END IF
         END IF
!!COHESION
      
         IF(DES_CONTINUUM_COUPLED) THEN
            PC = 50*FACTOR
!            OPEN (UNIT=1, FILE='called.out', STATUS='REPLACE')
!            WRITE (1,*) CALLED
            PRINT *,'D.E.S. COUPLED', CALLED
         ELSE 
            PC = 1150 
            PRINT *,'D.E.S', CALLED
         END IF
         
         IF(CALLED.LT.4) THEN
            CALL CFASSIGN(PARTICLES)
         END IF

         IF(CALLED.GT.3) THEN
            IF(DES_PERIODIC_WALLS) THEN
               CALL PERIODIC_WALL_CALC_FORCE_DES(CALLED)
            ELSE IF(INLET_OUTLET) THEN
               CALL DES_INLET_OUTLET(CALLED)
            ELSE
               CALL CALC_FORCE_DES(CALLED)
            END IF
         END IF


         IF(CALLED.GT.3) THEN	

            IF(COUPLED_FLOW) THEN
               IF(PCOUNT.EQ.PC) THEN
                  CALL PRESSURE_DROP(PARTICLES)
               END IF
            END IF
            
            IF(PCOUNT.EQ.PC) THEN

!               OPEN (UNIT=91, FILE='des_XY1.OUT', STATUS='REPLACE')
!               WRITE (91,*) (DES_POS_NEW(K,20)/RADIUS_EQ,K=1,DIMN)
!               OPEN (UNIT=19, FILE='des_XY2.OUT', STATUS='REPLACE')
!               WRITE (19,*) (DES_POS_NEW(K,40)/RADIUS_EQ,K=1,DIMN)
!               OPEN (UNIT=29, FILE='des_XY3.OUT', STATUS='REPLACE')
!               WRITE (29,*) (DES_POS_NEW(K,60)/RADIUS_EQ,K=1,DIMN)
!               OPEN (UNIT=39, FILE='des_XY4.OUT', STATUS='REPLACE')
!               WRITE (39,*) (DES_POS_NEW(K,80)/RADIUS_EQ,K=1,DIMN)
!               OPEN (UNIT=49, FILE='des_XY5.OUT', STATUS='REPLACE') 
!               WRITE (49,*) (DES_POS_NEW(K,100)/RADIUS_EQ,K=1,DIMN)
!               OPEN (UNIT=59, FILE='des_XY6.OUT', STATUS='REPLACE')
!               WRITE (59,*) (DES_POS_NEW(K,130)/RADIUS_EQ,K=1,DIMN)
!               OPEN (UNIT=69, FILE='des_XY7.OUT', STATUS='REPLACE')
!               WRITE (69,*) (DES_POS_NEW(K,160)/RADIUS_EQ,K=1,DIMN)
!               OPEN (UNIT=79, FILE='des_XY8.OUT', STATUS='REPLACE')
!               WRITE (79,*) (DES_POS_NEW(K,200)/RADIUS_EQ,K=1,DIMN)
!               OPEN (UNIT=89, FILE='des_XY9.OUT', STATUS='REPLACE')
!               WRITE (89,*) (DES_POS_NEW(K,240)/RADIUS_EQ,K=1,DIMN)
!               OPEN (UNIT=99, FILE='des_XY10.OUT', STATUS='REPLACE') 
!               WRITE (99,*) (DES_POS_NEW(K,290)/RADIUS_EQ,K=1,DIMN)

               IF(CALLED.LT.(FACTOR/DTSOLID)) THEN
                  OPEN (UNIT=9, FILE='des_all-particles-1.out', STATUS='REPLACE')
                  WRITE (9,*)' '
                  WRITE (9,*)'ZONE T="',CALLED*DTSOLID,'"'
                  DO LNN = 1, PARTICLES
                     WRITE (9,*) ((DES_POS_NEW(K,LNN)/RADIUS_EQ),K=1,DIMN),&
                              (DES_VEL_NEW(K,LNN),K=1,DIMN), PR(LNN), RO_Sol(LNN)
                  END DO
               ELSE IF(CALLED.LT.(2*FACTOR/DTSOLID)) THEN
                  OPEN (UNIT=9, FILE='des_all-particles-2.out', STATUS='REPLACE')
                  WRITE (9,*)' '
                  WRITE (9,*)'ZONE T="',CALLED*DTSOLID,'"'
                  DO LNN = 1, PARTICLES
                     WRITE (9,*) ((DES_POS_NEW(K,LNN)/RADIUS_EQ),K=1,DIMN),&
                              (DES_VEL_NEW(K,LNN),K=1,DIMN), PR(LNN), RO_Sol(LNN)
                  END DO
               ELSE IF(CALLED.LT.(3*FACTOR/DTSOLID)) THEN
                  OPEN (UNIT=9, FILE='des_all-particles-3.out', STATUS='REPLACE')
                  WRITE (9,*)' '
                  WRITE (9,*)'ZONE T="',CALLED*DTSOLID,'"'
                  DO LNN = 1, PARTICLES
                     WRITE (9,*) ((DES_POS_NEW(K,LNN)/RADIUS_EQ),K=1,DIMN),&
                             (DES_VEL_NEW(K,LNN),K=1,DIMN), PR(LNN), RO_Sol(LNN)
                  END DO
               ELSE IF(CALLED.LT.(4*FACTOR/DT)) THEN
                  OPEN (UNIT=9, FILE='des_all-particles-4.out', STATUS='REPLACE')
                  WRITE (9,*)' '
                  WRITE (9,*)'ZONE T="',CALLED*DTSOLID,'"'
                  DO LNN = 1, PARTICLES
                     WRITE (9,*) ((DES_POS_NEW(K,LNN)/RADIUS_EQ),K=1,DIMN),&
                              (DES_VEL_NEW(K,LNN),K=1,DIMN), PR(LNN), RO_Sol(LNN)
                  END DO
               ELSE IF(CALLED.LE.(5*FACTOR/DT)) THEN
                  OPEN (UNIT=9, FILE='des_all-particles-5.out', STATUS='REPLACE')
                  WRITE (9,*)' '
                  WRITE (9,*)'ZONE T="',CALLED*DTSOLID,'"'
                  DO LNN = 1, PARTICLES
                     WRITE (9,*) ((DES_POS_NEW(K,LNN)/RADIUS_EQ),K=1,DIMN),&
                              (DES_VEL_NEW(K,LNN),K=1,DIMN), PR(LNN), RO_Sol(LNN)
                  END DO
               END IF
               PCOUNT = 0
            END IF
            PCOUNT = PCOUNT + 1
         END IF
         CALLED = CALLED + 1


      END DO

      IF(CALLED .GT. CALLED_MAX) THEN
         PRINT *,'TIME =', TIME
         STOP
      END IF

      END SUBROUTINE DES_TIME_MARCH
