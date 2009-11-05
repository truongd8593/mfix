
      SUBROUTINE DES_ALLOCATE_ARRAYS 
      
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!                                                                    
!  Module name: DES_ALLOCATE_ARRAYS                                     
!  Purpose: allocate arrays for DES
!                                  
!                                                                   
!-----------------------------------------------
!   M o d u l e s 
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
! the number of particles in the system      
      INTEGER NPARTICLES
!-----------------------------------------------
      INCLUDE 'function.inc'



      WRITE(*,'(1X,A)')&
         '---------- START DES_ALLOCATE_ARRAYS ---------->'

      IF(GENER_PART_CONFIG) THEN 
         WRITE(*,'(3X,A)') 'Checking usr info for gener_part_config'
         WRITE(*,'(3X,A,3(G15.8))') &
            'des_eps_xstart, ystart, zstart = ',&
            DES_EPS_XSTART, DES_EPS_YSTART, DES_EPS_ZSTART

! Check if des_eps_xstart and related variables are correctly initialized in mfix.dat file
         IF(DIMN.EQ.2) THEN 
            IF(DES_EPS_XSTART.EQ.UNDEFINED.OR.&
               DES_EPS_YSTART.EQ.UNDEFINED) THEN
               WRITE (*,'(3X,A,/,5X,A,/,5X,A)') 'ERROR MESSAGE: ', &
                  'des_eps_xstart or des_eps_ystart ',&
                  'not properly defined in the input file'
               CALL MFIX_EXIT(myPE)
            ENDIF
         ELSE
            IF(DES_EPS_XSTART.EQ.UNDEFINED.OR.&
               DES_EPS_YSTART.EQ.UNDEFINED.OR.&
               DES_EPS_ZSTART.EQ.UNDEFINED) THEN
               WRITE (*,'(3X,A,/,5X,A,/,5X,A)') 'ERROR MESSAGE: ', &
                  'des_eps_xstart, des_eps_ystart or des_eps_zstart ',&
                  'not properly defined in the input file'
               CALL MFIX_EXIT(myPE)
            ENDIF
         ENDIF

         DO M = 1, MMAX
            IF(VOL_FRAC(M) == UNDEFINED) THEN
               WRITE (*,'(3X,A,A)') &
                  'VOL_FRAC(M) must be defined in ',&
                  'mfix.dat for M = 1, MMAX'
               CALL MFIX_EXIT(myPE)
            ENDIF
            IF(VOL_FRAC(M) < ZERO .OR. VOL_FRAC(M) > (ONE-EP_STAR)) THEN
               WRITE (*,'(3X,A,A)') &
                  'Unphysical ( > 1-EP_STAR or < 0) ',&
                  'values of VOL_FRAC(M) set in mfix.dat'
               CALL MFIX_EXIT(myPE)
            ENDIF
         ENDDO
         
         WRITE(*,'(3X,A,A)') 'particle configuration will be ', &
            'automatically generated'
         PARTICLES = 0
         IF(DIMN.EQ.2) THEN 
            VOL_DOMAIN  = DES_EPS_XSTART*DES_EPS_YSTART*ZLENGTH 
            ! DZ(1) is not yet defined here.
         ELSE 
            VOL_DOMAIN  = DES_EPS_XSTART*DES_EPS_YSTART*DES_EPS_ZSTART
         ENDIF
         DO M = 1, MMAX
            PART_MPHASE(M) = FLOOR((6.D0*VOL_FRAC(M)*VOL_DOMAIN)/&
               (PI*(D_P0(M)**3)))
         ENDDO
         WRITE(*,*) '     MMAX = ', MMAX, 'PART_MPHASE = ',&
            PART_MPHASE, VOL_DOMAIN, D_P0(1), pi, &
            (6.D0*VOL_FRAC(1)*VOL_DOMAIN)/(PI*(D_P0(1)**3))
         PARTICLES = SUM(PART_MPHASE(1:MMAX))
         
      ENDIF !  end if gener_part_config

      WRITE(*,'(3X,A,I)') &
         'total number of particles = ', PARTICLES      
      WRITE(*,'(3X,A,I)') 'dimension = ', DIMN
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
! J.Musser : Dynamic Particle Info
      ALLOCATE( PEA (NPARTICLES, 3) )

      
! COEFF OF RESITUTIONS 
      ALLOCATE(REAL_EN(MMAX,MMAX)) !REAL_ET(MMAX,MMAX)) specify eta_t_fact instead of e_t (sof, dec-04-2008)
      ALLOCATE(REAL_EN_WALL(MMAX)) !REAL_ET_WALL(MMAX)) specify eta_t_w_fact instead of e_t_w (sof, dec-04-2008)

      ALLOCATE(DES_ETAN(MMAX,MMAX))
      ALLOCATE(DES_ETAT(MMAX,MMAX))
      ALLOCATE(DES_ETAN_WALL(MMAX), DES_ETAT_WALL(MMAX))

      IF(DES_INTERP_ON) THEN
         ALLOCATE(DRAG_AM(DIMENSION_I-1, DIMENSION_J-1, MAX(1,DIMENSION_K-1), MMAX))
         ALLOCATE(DRAG_BM(DIMENSION_I-1, DIMENSION_J-1,  MAX(1,DIMENSION_K-1), DIMN, MMAX))
         ALLOCATE(WTBAR(DIMENSION_I-1, DIMENSION_J-1,  MAX(1,DIMENSION_K-1),  MMAX))
         ALLOCATE(VEL_FP(NPARTICLES,3))
         ALLOCATE(F_gp(NPARTICLES ))  
         F_gp(1:NPARTICLES)  = ZERO
      ENDIF 

      ALLOCATE(MARK_PART(NPARTICLES))
      ALLOCATE(BED_HEIGHT(MMAX))
      ALLOCATE(AVE_VEL_X(DIMENSION_3,MMAX), AVE_VEL_Y(DIMENSION_3,MMAX), AVE_VEL_Z(DIMENSION_3,MMAX))


      ALLOCATE(PIC(DIMENSION_I,DIMENSION_J,DIMENSION_K))
      DO IJK = ijkstart3, ijkend3
         I = I_OF(IJK)
         J = J_OF(IJK)
         K = K_OF(IJK)
         NULLIFY(pic(i,j,k)%p)
      ENDDO

      
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
      IF(DIMN.GT.2) THEN
         Allocate(  OMEGA_OLD (NPARTICLES,DIMN) )
         Allocate(  OMEGA_NEW (NPARTICLES,DIMN) )
      ELSE
         Allocate(  OMEGA_OLD (NPARTICLES,1) )
         Allocate(  OMEGA_NEW (NPARTICLES,1) )
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
     
! Accumulated spring forces      
      Allocate(  PFN (NPARTICLES,MAXNEIGHBORS,DIMN) )
      Allocate(  PFT (NPARTICLES,MAXNEIGHBORS,DIMN) )
     
! Wall position, velocity and normal vector
      Allocate(  DES_WALL_POS (NWALLS,DIMN) )
      Allocate(  DES_WALL_VEL (NWALLS,DIMN) )
      Allocate(  WALL_NORMAL (NWALLS,DIMN) )
      Allocate(  PN (NPARTICLES, MAXNEIGHBORS) )
      Allocate(  PV (NPARTICLES, MAXNEIGHBORS) )
     
! Particles in a computational cell (for volume fraction) )
      Allocate(  PINC (DIMENSION_3) )
      Allocate(  PIJK (NPARTICLES,5) )
     
! Volume averaged solids volume in a cell      
      Allocate(  DES_U_s (DIMENSION_3, DIMENSION_M) )
      Allocate(  DES_V_s (DIMENSION_3, DIMENSION_M) )
      Allocate(  DES_W_s (DIMENSION_3, DIMENSION_M) )

! Averaged velocity obtained by avraging over all the particles
      ALLOCATE(DES_VEL_AVG(DIMN))

! Global Granular Energy
      ALLOCATE(GLOBAL_GRAN_ENERGY(DIMN))
      ALLOCATE(GLOBAL_GRAN_TEMP(DIMN))
    
! Drag exerted by the gas o solids
      Allocate(  SOLID_DRAG (DIMENSION_3, DIMENSION_M, DIMN) )
     
! Neighbor search
      Allocate(  NEIGHBOURS (NPARTICLES, MAXNEIGHBORS) )

      IF (DES_NEIGHBOR_SEARCH .EQ. 2 .OR. DES_NEIGHBOR_SEARCH .EQ. 3) THEN
         Allocate(  LQUAD (MAXQUADS, NMQD) )
         Allocate(  PQUAD (NPARTICLES) )
         Allocate(  CQUAD (MAXQUADS, NWALLS) )
      ENDIF
     
! Granular temperature
      Allocate(  DES_THETA (DIMENSION_3, DIMENSION_M) )
     
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
