!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!                                                                    
!  Module name: DES_ALLOCATE_ARRAYS                                     
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
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------      
! indices      
      INTEGER I, J, K, IJK, M 
! domain volume      
      DOUBLE PRECISION :: VOL_DOMAIN
! the maximum and minimum specified particle diameter 
      DOUBLE PRECISION MAX_DIAM, MIN_DIAM
! for gener_part_config, the total solids volume fraction
      DOUBLE PRECISION TOT_VOL_FRAC
! the number of particles in the system
      INTEGER NPARTICLES
!-----------------------------------------------
      INCLUDE 'function.inc'


      WRITE(*,'(1X,A)')&
         '---------- START DES_ALLOCATE_ARRAYS ---------->'

! Valid D_p0(M) are needed here if gener_part_config.  A second check 
! for realistic d_p0 values is made in check_data_04 but this routine 
! is called after des_allocate_arrays (i.e. would be too late).
! Valid D_P0(M) are also needed to identify which solids phase each
! particle belongs in (although the sorting is performed after 
! check_data_04 is called), and to determine the maximum particle size
! in the system (MAX_RADIUS), which in turn is used for various tasks
       MAX_DIAM = ZERO
       MIN_DIAM = LARGE_NUMBER
       DO M = 1,MMAX   
          IF (D_P0(M)<ZERO .OR. D_P0(M)==UNDEFINED) THEN    
             WRITE (*,'(3X,A,A)') &   
                'D_P0 must be defined and >0 in mfix.dat ',&   
                'for M = 1,MMAX'   
                CALL MFIX_EXIT(myPE)   
          ENDIF   
          MAX_DIAM = MAX(MAX_DIAM, D_P0(M))
          MIN_DIAM = MIN(MIN_DIAM, D_P0(M))
       ENDDO   
       DO M = MMAX+1, DIMENSION_M   
          IF (D_P0(M) /= UNDEFINED) THEN   
             WRITE (*,'(3X,A,A)') &   
                'Too many D_P0 are defined for given MMAX'   
             CALL MFIX_EXIT(myPE)   
          ENDIF   
       ENDDO
       MAX_RADIUS = 0.5d0*MAX_DIAM
       MIN_RADIUS = 0.5d0*MIN_DIAM


! If gener_part_config ensure various quantities are defined and valid
! ------------------------------------------------------------

      IF(GENER_PART_CONFIG) THEN 
         TOT_VOL_FRAC = ZERO
         WRITE(*,'(3X,A)') 'Checking usr info for gener_part_config'
         WRITE(*,'(5X,A,G15.8)') 'DES_EPS_XSTART = ', DES_EPS_XSTART
         WRITE(*,'(5X,A,G15.8)') 'DES_EPS_YSTART = ', DES_EPS_YSTART
         IF (DIMN .EQ. 3) &
            WRITE(*,'(5X,A,G15.8)') 'DES_EPS_ZSTART = ', DES_EPS_ZSTART

! Check if des_eps_xstart and related variables are correctly initialized in mfix.dat file
         IF(DIMN.EQ.2) THEN 
            IF(DES_EPS_XSTART.EQ.UNDEFINED.OR.&
               DES_EPS_YSTART.EQ.UNDEFINED) THEN
               WRITE (*,'(/,5X,A,A,A,/)') 'ERROR MESSAGE: ', &
                  'DES_EPS_XSTART or DES_EPS_YSTART ',&
                  'not defined in mfix.dat'
               CALL MFIX_EXIT(myPE)
            ENDIF
         ELSE
            IF(DES_EPS_XSTART.EQ.UNDEFINED.OR.&
               DES_EPS_YSTART.EQ.UNDEFINED.OR.&
               DES_EPS_ZSTART.EQ.UNDEFINED) THEN
               WRITE (*,'(/,5X,A,A,A,/)') 'ERROR MESSAGE: ', &
                  'DES_EPS_XSTART, DES_EPS_YSTART, or DES_EPS_ZSTART ',&
                  'not defined in mfix.dat'
               CALL MFIX_EXIT(myPE)
            ENDIF
         ENDIF

! perform a quick series of checks for quantities immediately needed in
! calculations; a more comprehensive series of checks is performed later
! (e.g., check_data_03, check_des_data which are called from get_data)
         DO M = 1, MMAX
            IF(VOL_FRAC(M) == UNDEFINED) THEN
               WRITE (*,'(/,5X,A,A,/)') &
                  'VOL_FRAC(M) must be defined in ',&
                  'mfix.dat for M = 1, MMAX'
               CALL MFIX_EXIT(myPE)
            ENDIF
            IF(VOL_FRAC(M) < ZERO .OR. VOL_FRAC(M) > (ONE-EP_STAR)) THEN
               WRITE (*,'(/,5X,A,A,/)') &
                  'Unphysical ( > 1-EP_STAR or < 0) ',&
                  'values of VOL_FRAC(M) set in mfix.dat'
               CALL MFIX_EXIT(myPE)
            ENDIF
            TOT_VOL_FRAC = TOT_VOL_FRAC + VOL_FRAC(M)
         ENDDO

         IF(TOT_VOL_FRAC > (ONE-EP_STAR)) THEN
            WRITE (*,'(/,5X,A,A,/7X,A,ES15.7,2X,A,ES15.7,/)') &
               'Total solids volume fraction should not exceed ', &
               '1-EP_STAR where', 'EP_STAR = ', EP_STAR, &
               ' and TOT_VOL_FRAC = ', TOT_VOL_FRAC
            CALL MFIX_EXIT(myPE)
         ENDIF

         WRITE(*,'(/5X,A,A,/)') 'Particle configuration will ', &
            'automatically be generated'

         PARTICLES = 0
         IF(DIMN.EQ.2) THEN 
! Values of DZ(1) or zlength are not guaranteed at this point; however,
! some value is needed to calculate the number of particles
            IF (DZ(1) == UNDEFINED .AND. ZLENGTH == UNDEFINED) THEN
               IF (MAX_DIAM .EQ. ZERO) MAX_DIAM = ONE ! for calculations
               WRITE(*,'(5X,A,A,/11X,A,A,ES15.7,/11X,A,A,/11X,A,/)') &
                  'NOTE: neither zlength or dz(1) were specified ',&
                  'so ZLENGTH is being set','to the maximum particle ',&
                  'diameter of ', MAX_DIAM, 'to provide a basis ',&
                  'for calculating the number of particles',&
                  'in the system'
! set zlength to ensure consistency with calculations later on especially 
! when conducting coupled simulations when zlength/dz(1) are not set in 
! the mfix.dat file
               ZLENGTH = MAX_DIAM 
               VOL_DOMAIN  = DES_EPS_XSTART*DES_EPS_YSTART*MAX_DIAM
            ELSEIF (DZ(1) == UNDEFINED) THEN
               WRITE(*,'(5X,A,G15.8,A,/11X,A,A,/11X,A,/)') &
                  'NOTE: specified zlength of ', ZLENGTH,&
                  ' is used','to provide a basis for calculating ',&
                  'the number of','particles in the system'
               VOL_DOMAIN  = DES_EPS_XSTART*DES_EPS_YSTART*ZLENGTH
            ELSE
               WRITE(*,'(5X,A,G15.8,A,/11X,A,A,/11X,A,/)') &
                  'NOTE: specified dz(1) of ', DZ(1),&
                  ' is used','to provide a basis for calculating ',&
                  'the number of','particles in the system'
               VOL_DOMAIN  = DES_EPS_XSTART*DES_EPS_YSTART*DZ(1)
            ENDIF
         ELSE 
            VOL_DOMAIN  = DES_EPS_XSTART*DES_EPS_YSTART*DES_EPS_ZSTART
         ENDIF

         DO M = 1, MMAX
            PART_MPHASE(M) = FLOOR((6.D0*VOL_FRAC(M)*VOL_DOMAIN)/&
               (PI*(D_P0(M)**3)))
         ENDDO
         WRITE(*,'(5X,A,I5,2X,A,G15.7)') 'MMAX = ', MMAX, &
           ' VOL_DOMAIN = ', VOL_DOMAIN
         WRITE(*,'(5X,A,/7X,(ES15.7,2X,$))') 'D_P0(M) = ', &
            D_P0(1:MMAX)
         WRITE(*,'')
         WRITE(*,'(5X,A,/7X,(G15.8,2X,$))') &
            'VOL_FRAC(M) (solids volume fraction of phase M) = ', &
            VOL_FRAC(1:MMAX)
         WRITE(*,'')
         WRITE(*,'(5X,A,/7X,(I,2X,$))') &
            'PART_MPHASE(M) (number particles in phase M) = ', &
            PART_MPHASE(1:MMAX)
         WRITE(*,'')
         PARTICLES = SUM(PART_MPHASE(1:MMAX))
      ENDIF !  end if gener_part_config

! End checks if gener_part_config    
! ------------------------------------------------------------


! J.Musser 
! If the particle count is not defined, but MAX_PIS, the maximum number
! of particles permitted in the system at any given time, is set, assume
! that a mass inlet has been specified and that the system is starting
! empty.  Further checks are conducted in check_des_bc which is called
! from get_data to verify this assumption.
      IF (.NOT.GENER_PART_CONFIG) THEN
         IF(PARTICLES == UNDEFINED_I .AND. MAX_PIS /= UNDEFINED_I)THEN
            PARTICLES = 0
         ELSEIF(PARTICLES == UNDEFINED_I .AND. MAX_PIS == UNDEFINED_I)THEN
            WRITE(*,'(3X,A)')&
               'Either PARTICLES or MAX_PIS must specified in mfix.dat'
            CALL MFIX_EXIT(myPE)
         ELSEIF(PARTICLES == 0 .AND. MAX_PIS == UNDEFINED_I) THEN
            WRITE(*,'(3X,A,A)')&
               'If starting with 0 PARTICLES, MAX_PIS must be ', &
               'specified in mfix.dat'
            CALL MFIX_EXIT(myPE)         
         ENDIF 
      ENDIF   ! if .not.gener_part_config


      WRITE(*,'(3X,A,I)') &
         'Total number of particles = ', PARTICLES      
      WRITE(*,'(3X,A,I5)') 'Dimension = ', DIMN
      NWALLS = 2*DIMN

      IF(.NOT.NON_RECT_BC) THEN
! +nwalls is included since calc_force_des temporarily uses the variables 
! pos, vel, etc at elements beyond the array index given by particles. 
! unclear why additional array space is needed via particles_factor
         NPARTICLES = PARTICLES * PARTICLES_FACTOR + NWALLS
      ELSE
! Tingwen 19/01/2008 
! +2 to include the contact with edge and node
! only one edge contact or one node contact with wall is allowed for a particle
         NPARTICLES = PARTICLES * PARTICLES_FACTOR + NWALLS + 2 + NWALLS +1
      ENDIF

! J.Musser : Dynamic Particle Info
      IF(MAX_PIS /= UNDEFINED_I .AND. &
         MAX_PIS .GT. NPARTICLES) NPARTICLES = MAX_PIS

      MAXQUADS = 5*PARTICLES*MQUAD_FACTOR
      IF(MAXQUADS.LE.80000) MAXQUADS = 80000
      MAXNEIGHBORS = MN + 1 + NWALLS

      IF(DIMN.EQ.3) THEN
         NMQD = 11
      ELSE
         NMQD = 7
      ENDIF


! DES Allocatable arrays
!-----------------------------------------------
! J.Musser : Dynamic Particle Info
      ALLOCATE( PEA (NPARTICLES, 3) )

! T. Li: Hertzian collision model
      allocate(hert_kn(MMAX,MMAX))
      allocate(hert_kt(MMAX,MMAX))
      allocate(hert_kwn(MMAX))
      allocate(hert_kwt(MMAX)) 
      allocate(g_mod(MMAX))
      
! COEFF OF RESITUTIONS 
      ALLOCATE(REAL_EN(MMAX,MMAX)) 
      ALLOCATE(REAL_EN_WALL(MMAX))
! for hertzian model need real_et, otherwise specify eta_t_fact 
      ALLOCATE(REAL_ET(MMAX,MMAX)) 
! for hertzian model need real_et_wall, otherwise specifiy eta_t_w_fact
      ALLOCATE(REAL_ET_WALL(MMAX)) 

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
      Allocate(  PPOS (NPARTICLES,DIMN) )
     
! Total, normal and tangetial forces      
      Allocate(  FC (NPARTICLES,DIMN) )
      Allocate(  FN (NPARTICLES,DIMN) )
      Allocate(  FT (NPARTICLES,DIMN) )
      Allocate(  GRAV (DIMN) )

! Torque     
      IF(DIMN.EQ.3) THEN 
         Allocate(  TOW (NPARTICLES,DIMN) )
      ELSE
         Allocate(  TOW (NPARTICLES,1) )
      ENDIF
     
! Accumulated spring force      
      Allocate(  PFT (NPARTICLES,MAXNEIGHBORS,DIMN) )
! Tracking variables for particle contact history
      Allocate(  PN (NPARTICLES, MAXNEIGHBORS) )
      Allocate(  PV (NPARTICLES, MAXNEIGHBORS) )
    
! Temporary variables to store wall position, velocity and normal vector
      Allocate(  DES_WALL_POS (NWALLS,DIMN) )
      Allocate(  DES_WALL_VEL (NWALLS,DIMN) )
      Allocate(  WALL_NORMAL  (NWALLS,DIMN) )

! Variable that stores the particle in cell information (ID) on the
! computational grid defined by imax, jmax and kmax in mfix.dat
      ALLOCATE(PIC(DIMENSION_I,DIMENSION_J,DIMENSION_K))
      DO K = 1,KMAX3
         DO J = 1,JMAX3
            DO I = 1,IMAX3
               NULLIFY(pic(i,j,k)%p) 
            ENDDO 
          ENDDO 
       ENDDO 

! Particles in a computational cell (for volume fraction)
      Allocate(  PINC (DIMENSION_3) )
! For each particle track its i,j,k location on computational grid
! defined by imax, jmax and kmax in mfix.dat and phase no.         
      Allocate(  PIJK (NPARTICLES,5) )

! For each particle track its i, j, k index according to the grid
! based search mesh when des_neighbor_search=4
      IF (DES_NEIGHBOR_SEARCH .EQ. 4) THEN      
         ALLOCATE( DESGRIDSEARCH_PIJK (NPARTICLES,3) )
      ENDIF

      IF(DES_INTERP_ON) THEN
         ALLOCATE(DRAG_AM(DIMENSION_I-1, DIMENSION_J-1, MAX(1,DIMENSION_K-1), MMAX))
         ALLOCATE(DRAG_BM(DIMENSION_I-1, DIMENSION_J-1,  MAX(1,DIMENSION_K-1), DIMN, MMAX))
         ALLOCATE(WTBAR(DIMENSION_I-1, DIMENSION_J-1,  MAX(1,DIMENSION_K-1),  MMAX))
         ALLOCATE(VEL_FP(NPARTICLES,3))
         ALLOCATE(F_gp(NPARTICLES ))  
         F_gp(1:NPARTICLES)  = ZERO
      ENDIF 
    
! Drag exerted by the gas on solids
      Allocate(  SOLID_DRAG (DIMENSION_3, DIMENSION_M, DIMN) )
     
! Neighbor search
      Allocate(  NEIGHBOURS (NPARTICLES, MAXNEIGHBORS) )

      IF (DES_NEIGHBOR_SEARCH .EQ. 2 .OR. DES_NEIGHBOR_SEARCH .EQ. 3) THEN
         Allocate(  LQUAD (MAXQUADS, NMQD) )
         Allocate(  PQUAD (NPARTICLES) )
         Allocate(  CQUAD (MAXQUADS, NWALLS) )
      ENDIF

      ALLOCATE(MARK_PART(NPARTICLES))
      ALLOCATE(BED_HEIGHT(MMAX))

! Volume averaged solids volume in a cell      
      Allocate(  DES_U_s (DIMENSION_3, DIMENSION_M) )
      Allocate(  DES_V_s (DIMENSION_3, DIMENSION_M) )
      Allocate(  DES_W_s (DIMENSION_3, DIMENSION_M) )

! Averaged velocity obtained by avraging over all the particles
      ALLOCATE(  DES_VEL_AVG(DIMN) )

! Granular temperature
      Allocate(  DES_THETA (DIMENSION_3, DIMENSION_M) )

! Global Granular Energy
      ALLOCATE(  GLOBAL_GRAN_ENERGY(DIMN) )
      ALLOCATE(  GLOBAL_GRAN_TEMP(DIMN) )
    
! Cell faces
      Allocate(  XE (DIMENSION_I) )
      Allocate(  YN (DIMENSION_J) )
      Allocate(  ZT (DIMENSION_K) )

     
! COHESION      
! Square-well potential parameters
      Allocate(  WELL_WIDTH (NPARTICLES) )
      Allocate(  WELL_DEPTH (NPARTICLES) )
! Does particle have at least one linked partner
      Allocate(  IS_LINKED (NPARTICLES) )
! Does particle have at least one aggloerated partner
      Allocate(  IS_AGGLOMERATED (NPARTICLES) )
! Array of linked partners
      Allocate(  LINKS (NPARTICLES, MAXNEIGHBORS) )
! Array of agglomerated partners
      Allocate(  AGGS (NPARTICLES, MAXNEIGHBORS) )
! Matrix location of particle 
      Allocate(  PART_GRID (NPARTICLES,4) )

      WRITE(*,'(1X,A)')&
         '<---------- END DES_ALLOCATE_ARRAYS ----------'

      RETURN
      END SUBROUTINE DES_ALLOCATE_ARRAYS 
