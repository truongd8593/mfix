!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFNEWVALUES(L)                                         C
!
!  Purpose: DES - Calculate the new values of particle velocity,       
!           position, angular velocity etc                            
!
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C 
!
!  Comments: Implements Eqns 1, 2, 3, 4 & 5  from the following paper  
!  Tsuji Y., Kawaguchi T., and Tanak T., "Lagrangian numerical         
!  simulation of plug glow of cohesionless particles in a              
!  horizontal pipe", Powder technology, 71, 239-250, 1992             
!
!                                                                      C 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CFNEWVALUES(L)

      USE param
      USE param1
      USE parallel
      USE matrix
      USE scales
      USE constant
      USE compar
      USE sendrecv
      USE physprop
      USE fldvar
      USE visc_g
      USE rxns
      USE run
      USE geometry
      USE indices
      USE drag
      USE discretelement
      USE des_bc

      IMPLICIT NONE
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
      INTEGER L, I
      DOUBLE PRECISION D(DIMN), DIST, &
                       NEIGHBOR_SEARCH_DIST, &
                       SQR_VEL
! Logical for local debug warnings
      LOGICAL DES_LOC_DEBUG

!-----------------------------------------------
! Functions 
!-----------------------------------------------
      DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT

!-----------------------------------------------      

      DES_LOC_DEBUG = .FALSE.

      IF(L.EQ.1) THEN 
         DES_KE = ZERO
         DES_PE = ZERO 
         DES_VEL_AVG = ZERO
      ENDIF

! If a particle is classified as new, then forces are ignored. 
! Classification from new to existing is performed in routine
! des_check_new_particle.f
      IF(.NOT.PEA(L,2))THEN 
         FC(L, :) = FC(L,:)/PMASS(L) + GRAV(:)
      ELSE 
         FC(L,:) = ZERO
         TOW(L,:) = ZERO         
      ENDIF

      DES_VEL_NEW(L,:) = FC(L,:)
      DES_VEL_NEW(L,:) = DES_VEL_OLD(L,:) + DES_VEL_NEW(L,:)*DTSOLID
      DES_POS_NEW(L,:) = DES_POS_OLD(L,:) + DES_VEL_NEW(L,:)*DTSOLID 
      OMEGA_NEW(L,:) = OMEGA_OLD(L,:) + TOW(L,:)*OMOI(L)*DTSOLID


! Check if the particle has moved a distance greater than or equal to
! its radius since the last time a neighbor search was called; if so, 
! make sure that neighbor is called in des_time_march
      IF(.NOT.DO_NSEARCH) THEN
         D(:) = DES_POS_NEW(L,:) - PPOS(L,:)
         DIST = SQRT(DES_DOTPRDCT(D,D))
         NEIGHBOR_SEARCH_DIST = NEIGHBOR_SEARCH_RAD_RATIO*DES_RADIUS(L)
         IF(DIST.GE.NEIGHBOR_SEARCH_DIST) DO_NSEARCH = .TRUE.
      ENDIF


! Check if the particle has moved a distance greater than or equal to 
! its radius during one solids time step; if so, call stop
      D(:) = DES_POS_NEW(L,:) - DES_POS_OLD(L,:)
      DIST = SQRT(DES_DOTPRDCT(D,D))
      IF(DIST.GE.DES_RADIUS(L)) THEN
         WRITE(*,1002) L
         STOP
      ENDIF


! Warning message for particles moving into ghost cells:
! Note that if this occurs then the particles_in_cell s
! subroutine should call stop
      IF((DES_POS_NEW(L,1) < ZERO .OR. DES_POS_NEW(L,1) > XLENGTH) .AND.&
      .NOT.DES_PERIODIC_WALLS_X .AND. &
      .NOT.PEA(L,2) .AND. .NOT.PEA(L,3)) THEN
! A new or exiting particle may exist in ghost cells
         IF (.NOT.DES_LOC_DEBUG) THEN
            DES_LOC_DEBUG = .TRUE.
            WRITE(*,1000) 
         ENDIF         
         WRITE(*,'(5X,A,I)') &
            'X position outside domain for particle: ', L
         WRITE(*,'(7X,A,(ES15.9))')&
            'particle pos : ', DES_POS_NEW(L,:)
         WRITE(*,'(7X,A,(ES15.9))')&
            'particle vel : ', DES_VEL_NEW(L,:)
      ENDIF 

      IF((DES_POS_NEW(L,2) < ZERO .OR. DES_POS_NEW(L,2) > YLENGTH) .AND.&
      .NOT.DES_PERIODIC_WALLS_Y .AND. &
      .NOT.PEA(L,2) .AND. .NOT.PEA(L,3)) THEN
! A new or exiting particle may exist in ghost cells
         IF (.NOT.DES_LOC_DEBUG) THEN
            DES_LOC_DEBUG = .TRUE.
            WRITE(*,1000) 
         ENDIF         
         WRITE(*,'(5X,A,I)') &
            'Y position outside domain for particle: ', L
         WRITE(*,'(7X,A,(ES15.9))')&
            'particle pos : ', DES_POS_NEW(L,:)
         WRITE(*,'(7X,A,(ES15.9))')&
            'particle vel : ', DES_VEL_NEW(L,:)
      ENDIF 

      IF (DIMN > 2) THEN
         IF((DES_POS_NEW(L,3) < ZERO .OR. DES_POS_NEW(L,3) > ZLENGTH) .AND.&
         .NOT.DES_PERIODIC_WALLS_Z .AND. &
         .NOT.PEA(L,2) .AND. .NOT.PEA(L,3)) THEN
! A new or exiting particle may exist in ghost cells
            IF (.NOT.DES_LOC_DEBUG) THEN
               DES_LOC_DEBUG = .TRUE.
               WRITE(*,1000) 
            ENDIF         
            WRITE(*,'(5X,A,I)') &
               'Z position outside domain for particle: ', L
            WRITE(*,'(7X,A,(ES15.9))')&
               'particle pos : ', DES_POS_NEW(L,:)
            WRITE(*,'(7X,A,(ES15.9))')&
               'particle vel : ', DES_VEL_NEW(L,:)
         ENDIF
      ENDIF 


! Periodic treatment
      IF(DES_PERIODIC_WALLS) THEN
         IF(DES_PERIODIC_WALLS_X) THEN
            IF(DES_POS_NEW(L,1).GT.EX2) THEN
               DES_POS_NEW(L,1) = DES_POS_NEW(L,1) - (EX2 - WX1)
               PIJK(L,1) = 2
            ELSE IF(DES_POS_NEW(L,1).LT.WX1) THEN
               DES_POS_NEW(L,1) = DES_POS_NEW(L,1) + (EX2 - WX1)
               PIJK(L,1) = IMAX1
            ENDIF
         ENDIF
         IF(DES_PERIODIC_WALLS_Y) THEN
            IF(DES_POS_NEW(L,2).GT.TY2) THEN
               DES_POS_NEW(L,2) = DES_POS_NEW(L,2) - (TY2 - BY1)
               PIJK(L,2) = 2
            ELSEIF(DES_POS_NEW(L,2).LT.BY1) THEN
               DES_POS_NEW(L,2) = DES_POS_NEW(L,2) + (TY2 - BY1)
               PIJK(L,2) = JMAX1
            ENDIF
         ENDIF
         IF(DES_PERIODIC_WALLS_Z) THEN
            IF(DES_POS_NEW(L,3).GT.NZ2) THEN
               DES_POS_NEW(L,3) = DES_POS_NEW(L,3) - (NZ2 - SZ1)
               PIJK(L,3) = 2
            ELSEIF(DES_POS_NEW(L,3).LT.SZ1) THEN
               DES_POS_NEW(L,3) = DES_POS_NEW(L,3) + (NZ2 - SZ1)
               PIJK(L,3) = KMAX1
            ENDIF
         ENDIF
      ENDIF


! Calculation of total kinetic and potential energy and average velocity
      SQR_VEL = ZERO 
      DO I = 1, DIMN
         SQR_VEL = SQR_VEL + DES_VEL_NEW(L,I)**2
      ENDDO 
      DES_VEL_AVG(:) =  DES_VEL_AVG(:) + DES_VEL_NEW(L,:)
      DES_KE = DES_KE + HALF*SQR_VEL*PMASS(L) 
      DES_PE = DES_PE + PMASS(L)*DES_POS_NEW(L,2)*SQRT(DES_DOTPRDCT(GRAV,GRAV))


! Reset total contact force and torque      
      FC(L,:) = ZERO
      TOW(L,:) = ZERO


      IF (DES_LOC_DEBUG) WRITE(*,1001)

 1000 FORMAT(3X,'---------- FROM CFNEWVALUES ---------->')
 1001 FORMAT(3X,'<---------- END CFNEWVALUES ----------')  

 1002 FORMAT(/1X,70('*')//&
         ' From: CFNEWVALUES -',/&
         ' Message: Movement undesired for particle: ',I,/&
         1X,70('*')/)         

      RETURN
      END SUBROUTINE CFNEWVALUES


