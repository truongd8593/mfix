!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: MAKE_ARRAYS_DES                                        C
!  Purpose: DES - allocating DES arrays                                C
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
      USE run
      USE constant
      USE physprop

      IMPLICIT NONE
      
      INTEGER LN, K, M, NP
!      INTEGER CHECK_MPI
      INTEGER L, I, II, PART_COUNT
      DOUBLE PRECISION DIST, R_LM, DOML(DIMN)

      IF(DES_NEIGHBOR_SEARCH.EQ.UNDEFINED_I) THEN
         DES_NEIGHBOR_SEARCH = 1
         WRITE(*,*) 'Default N-Square search will be used'
      END IF

      IF(DES_NEIGHBOR_SEARCH.EQ.1) THEN
         DO_NSQUARE = .TRUE.
         WRITE(*,*) 'DEM USING N-SQUARE Search'
      ELSE IF(DES_NEIGHBOR_SEARCH.EQ.2) THEN
         DO_QUADTREE = .TRUE.
         WRITE(*,*) 'DEM USING QUADTREE Search'
      ELSE IF(DES_NEIGHBOR_SEARCH.EQ.3) THEN
         DO_OCTREE = .TRUE.
         WRITE(*,*) 'DEM USING OCTREE Search'
      ELSE IF(DES_NEIGHBOR_SEARCH.EQ.4) THEN
         DO_GRID_BASED_SEARCH = .TRUE.
         WRITE(*,*) 'DEM USING CELL LINKED SEARCH'
      END IF
      
      IF(RUN_TYPE == 'NEW') THEN ! Fresh run
         PRINT*,'PARTICLES  = ', particles
         PRINT*,'DIMN  = ', dimn
         
         IF(.NOT.GENER_PART_CONFIG) THEN 
            OPEN(UNIT=10, FILE='particle_input.dat', STATUS='OLD') 
                     
            WRITE(*,*) 'READING PARTICLE CONFIGURATION FROM THE supplied particle_input.dat file'

            DO LN = 1, PARTICLES
               READ (10, *) (DES_POS_OLD(LN,K),K=1,DIMN),DES_RADIUS(LN),RO_Sol(LN) ,(DES_VEL_OLD(LN,K),K=1,DIMN)
               OMEGA_OLD(LN,:) = ZERO
               DES_POS_NEW(LN,:) = DES_POS_OLD(LN,:)
               DES_VEL_NEW(LN,:) = DES_VEL_OLD(LN,:)
               OMEGA_NEW(LN,:) = OMEGA_OLD(LN,:)
            END DO
         ELSE
            call generate_particle_config
         ENDIF
            
      ELSE IF(RUN_TYPE == 'RESTART_1') THEN !  Read Restart
         CALL READ_DES_RESTART
         DES_POS_NEW(:,:) = DES_POS_OLD(:,:)
         DES_VEL_NEW(:,:) = DES_VEL_OLD(:,:)
         OMEGA_NEW(:,:) = OMEGA_OLD(:,:)
         DESRESDT = 0.0d0
         WRITE(*,*) 'DES_RES file read at Time= ', TIME
         WRITE(UNIT_LOG,*) 'DES_RES file read at Time= ', TIME
         IF(USE_COHESION) THEN
            WRITE(UNIT_LOG,*) 'Restart 1 is not implemented with DES-COHESION'
            WRITE(*,*) 'Restart 1 is not implemented with DES-COHESION'
            CALL MFIX_EXIT(myPE)
         END IF
         
      ELSE IF (RUN_TYPE == 'RESTART_2') THEN 
         WRITE(UNIT_LOG,*) 'Restart 2 is not implemented with DES'
         WRITE(*,*) 'Restart 2 is not implemented with DES'
         CALL MFIX_EXIT(myPE)
      END IF
      !des_radius(2) = 2.d0*des_radius(1)
      CALL CFASSIGN
      
      CALL PARTICLES_IN_CELL
      IF(PVEL_StDev.GT.ZERO)       CALL init_particles_jn
      !IF(RUN_TYPE == 'NEW'.and.DES_INTERP_ON.AND.DES_CONTINUUM_COUPLED) CALL SET_INITIAL_VELOCITY
       
      CALL writeic
      
      RETURN
      END SUBROUTINE MAKE_ARRAYS_DES 
