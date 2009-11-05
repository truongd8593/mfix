!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: MAKE_ARRAYS_DES                                        C
!  Purpose: DES - allocating DES arrays                                
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer: Rahul Garg                               Date: 01-Aug-07  C
!  Comments: Added some calls that are necessary if INTERPOLATION IS ON C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE MAKE_ARRAYS_DES
      
      USE param1
      USE geometry
      USE funits
      USE compar      
      USE discretelement
      USE des_bc
      USE run
      USE constant
      USE physprop

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------      
      INTEGER K, L

!-----------------------------------------------


      WRITE(*,'(1X,A)')&
         '---------- START MAKE_ARRAYS_DES ---------->'

      IF(DES_NEIGHBOR_SEARCH.EQ.1) THEN
         WRITE(*,'(3X,A)') 'dem using n-square search'
      ELSE IF(DES_NEIGHBOR_SEARCH.EQ.2) THEN
         WRITE(*,'(3X,A)') 'dem using quadtree search'
      ELSE IF(DES_NEIGHBOR_SEARCH.EQ.3) THEN
         WRITE(*,'(3X,A)') 'dem using octree search'
      ELSE IF(DES_NEIGHBOR_SEARCH.EQ.4) THEN
         WRITE(*,'(3X,A)') 'dem using linked cell search'
      ENDIF
      
      IF(RUN_TYPE == 'NEW') THEN ! Fresh run
         
         IF(.NOT.GENER_PART_CONFIG) THEN 
            OPEN(UNIT=10, FILE='particle_input.dat', STATUS='OLD') 
                     
            WRITE(*,'(3X,A,/,5X,A)') &
               'reading particle configuration from the supplied ',& 
               'particle_input.dat file'

            DO L = 1, PARTICLES
               READ (10,*) (DES_POS_OLD(L,K),K=1,DIMN),DES_RADIUS(L),RO_Sol(L) ,(DES_VEL_OLD(L,K),K=1,DIMN)
               OMEGA_OLD(L,:) = ZERO
               DES_POS_NEW(L,:) = DES_POS_OLD(L,:)
               DES_VEL_NEW(L,:) = DES_VEL_OLD(L,:)
               OMEGA_NEW(L,:) = OMEGA_OLD(L,:)
            END DO
         ELSE
            call generate_particle_config
         ENDIF
            
      ELSE IF(RUN_TYPE == 'RESTART_1') THEN !  Read Restart
         CALL READ_DES_RESTART
         DES_POS_NEW(:,:) = DES_POS_OLD(:,:)
         DES_VEL_NEW(:,:) = DES_VEL_OLD(:,:)
         OMEGA_NEW(:,:) = OMEGA_OLD(:,:)
         WRITE(*,'(3X,A,G17.8)') 'DES_RES file read at Time= ', TIME
         WRITE(UNIT_LOG,*) 'DES_RES file read at Time= ', TIME
         IF(USE_COHESION) THEN
            WRITE(UNIT_LOG,*) &
               'Restart 1 is not implemented with DES-COHESION'
            WRITE(*,'(3X,A)') &
               'Restart 1 is not implemented with DES-COHESION'
            CALL MFIX_EXIT(myPE)
         END IF
         
      ELSE IF (RUN_TYPE == 'RESTART_2') THEN 
         WRITE(UNIT_LOG,*) 'Restart 2 is not implemented with DES'
         WRITE(*,'(3X,A)') 'Restart 2 is not implemented with DES'
         CALL MFIX_EXIT(myPE)
      END IF
      !des_radius(2) = 2.d0*des_radius(1)
      CALL CFASSIGN

! J.Musser      
! Check data for des mass inflow boundary condtion 
! dtsolid is needed so call is made after cfassign.f       
      CALL CHECK_DES_BC

      CALL PARTICLES_IN_CELL
      IF(PVEL_StDev.GT.ZERO)       CALL init_particles_jn
      !IF(RUN_TYPE == 'NEW'.and.DES_INTERP_ON.AND.DES_CONTINUUM_COUPLED) CALL SET_INITIAL_VELOCITY
       
      CALL writeic

      WRITE(*,'(1X,A)')&
         '<---------- END MAKE_ARRAYS_DES ----------'

      RETURN

      END SUBROUTINE MAKE_ARRAYS_DES 
