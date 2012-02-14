!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!                                                                    
!  Subroutine: DES_ALLOCATE_ARRAYS                                     
!  Purpose: allocate arrays for DES
!  
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

      SUBROUTINE DES_ALLOCATE_ARRAYS 
                                                                   
!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param 
      USE param1
      USE constant
      USE discretelement
      Use indices
      Use geometry
      Use compar
      Use physprop
      Use des_bc
      use funits
      use desgrid 
      use desmpi
      USE mfix_pic
      Use des_thermo
      Use des_rxns      
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------      
! indices      
      INTEGER :: I, J, K, IJK, M 
! the number of particles in the system
      INTEGER :: NPARTICLES
!-----------------------------------------------
! Include statment functions
!-----------------------------------------------
      INCLUDE 'function.inc'
!-----------------------------------------------


      IF(DMP_LOG.AND.DEBUG_DES) WRITE(UNIT_LOG,'(1X,A)') &
         '---------- START DES_ALLOCATE_ARRAYS ---------->'

      IF(DMP_LOG) WRITE(UNIT_LOG,'(3X,A,I10)') &
         'Total number of particles = ', PARTICLES      
      IF(DMP_LOG) WRITE(UNIT_LOG,'(3X,A,I5)') 'Dimension = ', DIMN


      NWALLS = 2*DIMN
      MAXNEIGHBORS = MN + 1 + NWALLS

! defining the local variables nparticles  
!----------------------------------------------------------------->>>      

! +nwalls is included since calc_force_des temporarily uses the variables 
! pos, vel, etc at elements beyond the array index given by particles. 
! unclear why additional array space is needed via particles_factor
      NPARTICLES = PARTICLES * PARTICLES_FACTOR + NWALLS

! J.Musser : Dynamic Particle Info
! no real need to print this information out       
      IF(MAX_PIS /= UNDEFINED_I .AND. MAX_PIS .GT. NPARTICLES) THEN
         IF(DMP_LOG) WRITE(UNIT_LOG, 1002) MAX_PIS, NPARTICLES, MAX_PIS
         IF(DMP_LOG) WRITE(*, 1002) MAX_PIS, NPARTICLES, MAX_PIS
         NPARTICLES = MAX_PIS
      ENDIF

! For parallel processing the array size required should be either 
! specified by the user or could be determined from total particles 
! with some factor
      IF (numPEs.GT.1) THEN  
        NPARTICLES = (NPARTICLES/numPEs)
      ENDIF

! minimum size for nparticles      
      IF(NPARTICLES.LT.1000) NPARTICLES = 1000

! max_pip adjusted to accomodate temporary variables used for walls 
      MAX_PIP = NPARTICLES - 2*NWALLS - 3

      IF(DMP_LOG) WRITE(unit_log, 1003)  NPARTICLES, MAX_PIP
      IF(DMP_LOG) WRITE(*, 1003)  NPARTICLES, MAX_PIP 

! check if max_pip is less than maximum number of computational
! particles for the MPPIC case if so, stop the code after printing
! the error message
      IF(GENER_PART_CONFIG.AND.MPPIC) THEN 
         IF(MAX_PIP.LT. SUM(PART_MPHASE(1:MMAX))) THEN 
            IF(DMP_LOG) WRITE(UNIT_LOG,1001) MAX_PIP, SUM(PART_MPHASE(1:MMAX))
            CALL MFIX_EXIT(myPE)
         ENDIF
      ENDIF
!-----------------------------------------------------------------<<<




! DES Allocatable arrays
!-----------------------------------------------
! Dynamic particle info including another index for parallel
! processing for ghost
      ALLOCATE( PEA (NPARTICLES, 4) )

      ALLOCATE (iglobal_id(nparticles))

! J.Musser: Allocate necessary arrays for discrete mass inlets      
      IF(DES_BCMI /= 0 .OR. DES_BCMO /=0) CALL ALLOCATE_DES_MIO


! T. Li: Hertzian collision model
      allocate(hert_kn(MMAX,MMAX))
      allocate(hert_kt(MMAX,MMAX))
      allocate(hert_kwn(MMAX))
      allocate(hert_kwt(MMAX)) 
      allocate(g_mod(MMAX))
      
! Coefficients of normal restitution
      ALLOCATE(REAL_EN(MMAX,MMAX)) 
      ALLOCATE(REAL_EN_WALL(MMAX))

! Coefficients of tangential restitution (needed for hertzian model)      
      ALLOCATE(REAL_ET(MMAX,MMAX)) 
      ALLOCATE(REAL_ET_WALL(MMAX)) 

! normal and tangential dampening factors
      ALLOCATE(DES_ETAN(MMAX,MMAX))
      ALLOCATE(DES_ETAT(MMAX,MMAX))
      ALLOCATE(DES_ETAN_WALL(MMAX), DES_ETAT_WALL(MMAX))
      
! Particle attributes
! Radius, density, mass, moment of inertia           
      Allocate(  DES_RADIUS (NPARTICLES) )
      Allocate(  RO_Sol (NPARTICLES) )
      Allocate(  PVOL (NPARTICLES) )
      Allocate(  PMASS (NPARTICLES) )
      Allocate(  OMOI (NPARTICLES) )

! Old and new particle positions, velocities (translational and
! rotational)       
      Allocate(  DES_POS_OLD (NPARTICLES,DIMN) )
      Allocate(  DES_POS_NEW (NPARTICLES,DIMN) )
      Allocate(  DES_VEL_OLD (NPARTICLES,DIMN) )
      Allocate(  DES_VEL_NEW (NPARTICLES,DIMN) )
      Allocate(  DES_VEL_OOLD(NPARTICLES,DIMN) )
      Allocate(  DES_ACC_OLD (NPARTICLES,DIMN) )

      IF(DIMN.GT.2) THEN
         Allocate(  OMEGA_OLD (NPARTICLES,DIMN) )
         Allocate(  OMEGA_NEW (NPARTICLES,DIMN) )
         ALLOCATE(  ROT_ACC_OLD (NPARTICLES,DIMN))
      ELSE
         Allocate(  OMEGA_OLD (NPARTICLES,1) )
         Allocate(  OMEGA_NEW (NPARTICLES,1) )
         ALLOCATE(  ROT_ACC_OLD (NPARTICLES,1))
      ENDIF        

! Particle positions at the last call neighbor search algorithm call 
      Allocate(  PPOS (NPARTICLES,DIMN) )
    
! Total, normal and tangetial forces      
      Allocate(  FC (NPARTICLES,DIMN) )
      Allocate(  FN (NPARTICLES,DIMN) )
      Allocate(  FT (NPARTICLES,DIMN) )

! Torque     
      IF(DIMN.EQ.3) THEN 
         Allocate(  TOW (NPARTICLES,DIMN) )
      ELSE
         Allocate(  TOW (NPARTICLES,1) )
      ENDIF
     
! Accumulated spring force      
      Allocate(  PFT (NPARTICLES,MAXNEIGHBORS,DIMN) )

! Save the normal direction at previous time step
      Allocate(  PFN (NPARTICLES,MAXNEIGHBORS,DIMN) )

! Tracking variables for particle contact history
      Allocate(  PN (NPARTICLES, MAXNEIGHBORS) )
      Allocate(  PV (NPARTICLES, MAXNEIGHBORS) )

! Temporary variables to store wall position, velocity and normal vector
      Allocate(  WALL_NORMAL  (NWALLS,DIMN) ) 

! Cell faces
! Added 0 to store IMIN3 values 
      Allocate(  XE (0:DIMENSION_I) )
      Allocate(  YN (0:DIMENSION_J) )
      Allocate(  ZT (0:DIMENSION_K) )

! Gravity vector      
      Allocate(  GRAV (DIMN) )

! Neighbor search
      Allocate(  NEIGHBOURS (NPARTICLES, MAXNEIGHBORS) )
      IF (DES_NEIGHBOR_SEARCH .EQ. 2 .OR. DES_NEIGHBOR_SEARCH .EQ. 3) THEN
         MAXQUADS = 5*PARTICLES*MQUAD_FACTOR
         IF(MAXQUADS.LE.80000) MAXQUADS = 80000
         IF(DIMN.EQ.3) THEN
            NMQD = 11
         ELSE
            NMQD = 7
         ENDIF
         Allocate(  LQUAD (MAXQUADS, NMQD) )
         Allocate(  PQUAD (NPARTICLES) )
         Allocate(  CQUAD (MAXQUADS, NWALLS) )
      ENDIF

! Variable that stores the particle in cell information (ID) on the
! computational fluid grid defined by imax, jmax and kmax in mfix.dat
      ALLOCATE(PIC(DIMENSION_3))
      DO IJK=1,DIMENSION_3
        NULLIFY(pic(ijk)%p) 
      ENDDO 

! Particles in a computational fluid cell (for volume fraction)
      Allocate(  PINC (DIMENSION_3) )

! For each particle track its i,j,k location on computational fluid grid
! defined by imax, jmax and kmax in mfix.dat and phase no.  
      Allocate(  PIJK (NPARTICLES,5) )

      IF(DES_INTERP_ON) THEN
         ALLOCATE(DRAG_AM(DIMENSION_3, MMAX))
         ALLOCATE(DRAG_BM(DIMENSION_3, DIMN, MMAX))
         ALLOCATE(WTBAR(DIMENSION_3,  MMAX))
         ALLOCATE(VEL_FP(NPARTICLES,3))
         ALLOCATE(F_gp(NPARTICLES ))  
         F_gp(1:NPARTICLES)  = ZERO
      ENDIF

! force due to gas-pressure gradient 
      ALLOCATE(P_FORCE(DIMENSION_3,DIMN))
      
! Drag exerted by the gas on solids
      Allocate(  SOLID_DRAG (DIMENSION_3, MMAX, DIMN) )
     
! Volume averaged solids volume in a computational fluid cell      
      Allocate(  DES_U_s (DIMENSION_3, MMAX) )
      Allocate(  DES_V_s (DIMENSION_3, MMAX) )
      Allocate(  DES_W_s (DIMENSION_3, MMAX) )

! Volume of nodes 
      ALLOCATE(DES_VOL_NODE(DIMENSION_3))
! ratio of actual volume of nodes to volume of nodes not corrected for
! on the wall or being outside the domain 
      ALLOCATE(DES_VOL_NODE_RATIO(DIMENSION_3))


! MP-PIC related 
      IF(MPPIC) THEN 
! mppic requries the following but if des_interp_on is false
! then these quantities will not be allocated               
         IF (.NOT.DES_INTERP_ON) THEN
            ALLOCATE(DRAG_AM(DIMENSION_3, MMAX))
            ALLOCATE(DRAG_BM(DIMENSION_3, DIMN, MMAX))
            ALLOCATE(WTBAR(DIMENSION_3,  MMAX))
            ALLOCATE(VEL_FP(NPARTICLES,3))
            ALLOCATE(F_gp(NPARTICLES ))  
            F_gp(1:NPARTICLES)  = ZERO
         ENDIF

         Allocate(PS_FORCE_PIC(DIMENSION_3, DIMN))
         ALLOCATE(DES_STAT_WT(NPARTICLES))         
         ALLOCATE(DES_VEL_MAX(DIMN))
         ALLOCATE(PS_GRAD(NPARTICLES, DIMN))         
         ALLOCATE(AVGSOLVEL_P(NPARTICLES, DIMN))
         ALLOCATE(EPG_P(NPARTICLES))
      ENDIF 


! Granular temperature in a computational fluid cell
      Allocate(DES_THETA (DIMENSION_3, MMAX) )      

! Averaged velocity obtained by averaging over all the particles
      ALLOCATE(DES_VEL_AVG(DIMN) )

! Global Granular Energy
      ALLOCATE(GLOBAL_GRAN_ENERGY(DIMN) )
      ALLOCATE(GLOBAL_GRAN_TEMP(DIMN) )
    
! Somewhat free variable used to aid in manipulation
      ALLOCATE(MARK_PART(NPARTICLES))

! variable for bed height of solids phase M      
      ALLOCATE(BED_HEIGHT(MMAX))



     
! ---------------------------------------------------------------->>>
! BEGIN COHESION       
! Square-well potential parameters
      IF(USE_COHESION) THEN 
         IF(SQUARE_WELL)THEN 
           Allocate(  WELL_WIDTH (NPARTICLES) )
           Allocate(  WELL_DEPTH (NPARTICLES) )
! Does particle have at least one linked partner
           Allocate(  IS_LINKED (NPARTICLES) ) 
! Array of linked partners
           Allocate(  LINKS (NPARTICLES, MAXNEIGHBORS) )
         ENDIF

! Matrix location of particle  (should be allocated in case user wishes
! to invoke routines in /cohesion subdirectory
         Allocate(  PART_GRID (NPARTICLES,4) )

         Allocate(  FCohesive (NPARTICLES,DIMN) )
         Allocate(  PostCohesive (NPARTICLES) )
         
      ENDIF
! END COHESION      
! ----------------------------------------------------------------<<<


! ---------------------------------------------------------------->>>
! BEGIN Thermodynamic Allocation      
      IF(DES_ENERGY_EQ)THEN

! Particle temperature
         Allocate( DES_T_s_OLD( NPARTICLES ) )
         Allocate( DES_T_s_NEW( NPARTICLES ) )
! Specific heat
         Allocate( DES_C_PS( NPARTICLES ) )
! Species mass fractions comprising a particle. This array may not be
! needed for all thermo problems.
         Allocate( DES_X_s( NPARTICLES, MAX_DES_NMAX) )
! Convection Specific arrays
         IF(DES_CONV_EQ) &
            Allocate( Qcv( NPARTICLES ) )
! Conduction Specific arrays
         IF(DES_COND_EQ_PP) &
            Allocate( Qpp( NPARTICLES ) )
! Conduction Specific arrays
         IF(DES_COND_EQ_PFP) &
            Allocate( Qpfp( NPARTICLES ) )
! Radiation Specific arrays
         IF(DES_RADI_EQ) &
            Allocate( Qrd( NPARTICLES ) )
! Stores number of neighbors based on neighbor search
         IF(FIND_THERMO_NBRHD) &
            Allocate( THERMO_NBRHD( NPARTICLES, MAXNEIGHBORS ) )
! Allocate the history variables for Adams-Bashforth integration
         IF (TRIM(DES_INTG_METHOD) .EQ. 'ADAMS_BASHFORTH') &
            Allocate( Qtotal_OLD( NPARTICLES ) )
      ENDIF
! End Thermodynamic Allocation
! ----------------------------------------------------------------<<<
      

! ---------------------------------------------------------------->>>
! BEGIN Species Allocation 
      IF(ANY_DES_SPECIES_EQ)THEN
! Rate of solids phase production for each species
         Allocate( DES_R_sp( NPARTICLES, MAX_DES_NMAX) )
! Rate of solids phase consumption for each species
         Allocate( DES_R_sc( NPARTICLES, MAX_DES_NMAX) )
! Allocate the history variables for Adams-Bashforth integration
         IF (TRIM(DES_INTG_METHOD) .EQ. 'ADAMS_BASHFORTH') THEN
! Rate of chnage of particle mass
            Allocate( dMdt_OLD( NPARTICLES ) )
! Rate of chnage of particle mass percent species
            Allocate( dXdt_OLD( NPARTICLES, MAX_DES_NMAX) )
         ENDIF

! Energy generation from reaction (cal/sec)
         Allocate( Qint( NPARTICLES ) )
         IF( TRIM(REACTION_MODEL) == 'SHRINKING_CORE')THEN
! Radius of unreacted core
            Allocate( CORE_RAD( NPARTICLES ) )
! Density of unreacted core
            Allocate( CORE_RHO( NPARTICLES ) )
! Allocate the history variables for Adams-Bashforth integration
            IF (TRIM(DES_INTG_METHOD) .EQ. 'ADAMS_BASHFORTH') &
! Rate of chnage of radius of unreacted core
               Allocate( dRdt_OLD( NPARTICLES ) )
         ENDIF
      ENDIF
! End Species Allocation
! ----------------------------------------------------------------<<<

      IF(DMP_LOG.AND.DEBUG_DES) WRITE(UNIT_LOG,'(1X,A)')&
         '<---------- END DES_ALLOCATE_ARRAYS ----------'

 1001 FORMAT(/1X,70('*')//' From: DES_ALLOCATE_ARRAYS',/,&
         ' Message: MAX_PIP',4X,I4,4X,&
         'is smaller than number of particles ',4x,i4,4x,/,& 
         'INCREASE PARTICLE FACTOR',/1X,70('*')/)

 1002 FORMAT(/2X,'From: DES_ALLOCATE_ARRAYS',/2X,&
         'Message: User supplied MAX_PIS =',I10,/2X,&
         '> NPARTICLES (cummulative size of particle arrays) =',&
         I10,/2X,'Therefore, setting NPARTICLES to ',I10)
         
 1003 FORMAT(/2X,'From: DES_ALLOCATE_ARRAYS',/2X,&
         'Message: particle array size on each proc = ',I10,/2X,&
         'Maximum physical particles (MAX_PIP) one each proc = ',&
         I10,/2X,'Note that this value of MAX_PIP is only ',&
         'relevant for a new run',/2X,'For restarts, max_pip ',&
         'be set later on')

      RETURN
      END SUBROUTINE DES_ALLOCATE_ARRAYS 

      
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: ALLOCATE_DES_MIO                                        !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!  Author: J.Musser                                   Date: 17-Aug-09  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE ALLOCATE_DES_MIO

!-----------------------------------------------
! Modules
!-----------------------------------------------      
      USE des_bc
      USE discretelement
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER :: I     ! Loop counter for no. of DES_BCMI
!-----------------------------------------------

! Allocate/Initialize for inlets      
      IF(DES_BCMI /= 0)THEN

! Boundary condition ID array
         Allocate( DES_BC_MI_ID (DES_BCMI) )

! Distance offset of incoming particles in ghost cell
         Allocate( DES_BC_OFFSET (DES_BCMI) )

! Particle injection factor
         Allocate( PI_FACTOR (DES_BCMI) )

! Particle injection count (injection number)
         Allocate( PI_COUNT (DES_BCMI) )

! Particle injection time scale
         Allocate( DES_MI_TIME (DES_BCMI) )

! Boundary classification
         Allocate( DES_MI_CLASS (DES_BCMI) )
         Allocate( PARTICLE_PLCMNT (DES_BCMI) )

! Order inlet condition variables
! (only needed if particle_plcmt is assigned 'ordr')
         Allocate( MI_FACTOR (DES_BCMI) )
         Allocate( MI_WINDOW (DES_BCMI) )
         Allocate( MI_ORDER (DES_BCMI) )   ! type dmi
         Allocate( I_OF_MI ( DES_BCMI) )   ! type dmi
         Allocate( J_OF_MI ( DES_BCMI) )   ! type dmi

! Grid search loop counter array; 6 = no. of faces
         Allocate(  GS_ARRAY (DES_BCMI, 6) )

! Logical array stating if a bounday condition is polydisperse
         Allocate( DES_BC_POLY( DES_BCMI ) )

! Array used for polydisperse inlets: stores the particle number
! distribution of an inlet scaled with numfrac_limit
         Allocate( DES_BC_POLY_LAYOUT( DES_BCMI, NUMFRAC_LIMIT ) )

! Initializiation
! Logical for whether inlet is polydisperse         
         DES_BC_POLY(:) = .FALSE.
! Logical for inlet existance on IJK face         
         DES_MI_X = .FALSE.
         DES_MI_Y = .FALSE.
         DES_MI_Z = .FALSE.          
! Integer arrays
         DES_BC_MI_ID(:) = -1
         PI_FACTOR(:) = -1
         PI_COUNT(:) = -1
         MI_FACTOR(:) = -1
         MI_WINDOW(:) = -1
         GS_ARRAY(:,:) = -1
         DES_BC_POLY_LAYOUT(:,:) = -1
! Double precision arrays
         DES_MI_TIME(:) = UNDEFINED
! Character precision arrays
         DES_MI_CLASS(:) = UNDEFINED_C
         PARTICLE_PLCMNT(:) = UNDEFINED_C
! Derived data types
         DO I = 1,DES_BCMI
            NULLIFY( MI_ORDER(I)%VALUE )
            NULLIFY( I_OF_MI(I)%VALUE )
            NULLIFY( J_OF_MI(I)%VALUE )
         ENDDO

      ENDIF  ! end if des_bcmi /= 0


! Allocate/Initialize for outlets
      IF(DES_BCMO /= 0)THEN
 
! Boundary Condition ID array
         Allocate( DES_BC_MO_ID (DES_BCMO) )

! Boundary classification
         Allocate( DES_MO_CLASS (DES_BCMO) )

! Initializiation
! Integer arrays         
         DES_BC_MO_ID(:) = -1
! Character arrays         
         DES_MO_CLASS(:) = UNDEFINED_C         
! Logical for outlet existance on IJK face         
         DES_MO_X = .FALSE.
         DES_MO_Y = .FALSE.
         DES_MO_Z = .FALSE.         

      ENDIF   ! end if des_bcmo /= 0


      RETURN
      END SUBROUTINE ALLOCATE_DES_MIO

