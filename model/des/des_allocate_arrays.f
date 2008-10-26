
      SUBROUTINE DES_ALLOCATE_ARRAYS 
      
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!                                                                    
!  Module name: DES_ALLOCATE_ARRAYS                                     
!  Purpose: allocate arrays for DES
!                                                                    C
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
      IMPLICIT NONE
      
      INTEGER NPARTICLES, I,J,K, NP, M 
      DOUBLE PRECISION :: VOL_DOMAIN
      DIMENSION_I   = IMAX3
      DIMENSION_J   = JMAX3
      DIMENSION_K   = KMAX3
      DIMENSION_3   = (kend3-kstart3+1)*(jend3-jstart3+1)*(iend3-istart3+1)
      DIMENSION_3G   = IJKMAX3            
      DIMENSION_3L  = ijksize3_all(myPE)      
      DIMENSION_M   = MAX(1, MMAX)
      DIMENSION_4   = (kend4-kstart4+1)*(jend4-jstart4+1)*(iend4-istart4+1)
      
      IF(GENER_PART_CONFIG) THEN 
         WRITE(*,*) 'GENER_PART_CONFIG = ', GENER_PART_CONFIG
         WRITE(*,*) 'DES_EPS_XSTART, YSTART, ZSTART = ', DES_EPS_XSTART, DES_EPS_YSTART, DES_EPS_ZSTART

!     Check if des_eps_xstart and related variables are correctly initialized in mfix.dat file
         IF(DIMN.EQ.2) THEN 
            IF(DES_EPS_XSTART.EQ.UNDEFINED.OR.DES_EPS_YSTART.EQ.UNDEFINED) THEN
               WRITE(*,*)'ERROR MESSAGE: DES_EPS_XSTART OR DES_EPS_YSTART NOT PROPERLY DEFINED IN THE INPUT FILE'
               STOP
            ENDIF
         ELSE
            IF(DES_EPS_XSTART.EQ.UNDEFINED.OR.DES_EPS_YSTART.EQ.UNDEFINED.OR.DES_EPS_ZSTART.EQ.UNDEFINED) THEN
               WRITE(*,*)'ERROR MESSAGE: DES_EPS_XSTART OR DES_EPS_YSTART OR DES_EPS_ZSTART NOT PROPERLY DEFINED IN THE INPUT FILE'
               STOP
            ENDIF
         ENDIF
         
         WRITE(*,*) 'PARTILCE CONFIGURATION WILL BE AUTOMATICALLY GENERATED'
         PARTICLES = 0
         IF(DIMN.EQ.2) THEN 
            VOL_DOMAIN  = DES_EPS_XSTART*DES_EPS_YSTART*DZ(1)
         ELSE 
            VOL_DOMAIN  = DES_EPS_XSTART*DES_EPS_YSTART*DES_EPS_ZSTART
         ENDIF
         DO M = 1, MMAX
            PART_MPHASE(M) = FLOOR((6.D0*VOL_FRAC(M)*VOL_DOMAIN)/(PI*(D_P0(M)**3.D0)))
         ENDDO
         WRITE(*,*) 'MMAX = ', MMAX, 'PART_MPHASE = ',  PART_MPHASE, VOL_DOMAIN,D_P0(1), pi,  (6.D0*VOL_FRAC(1)*VOL_DOMAIN)/(PI*(D_P0(1)**3.D0))
         PARTICLES = SUM(PART_MPHASE(1:MMAX))
         
      ENDIF
      
      WRITE(*,*) 'TOTAL NUMBER OF PARTICLES  = ', PARTICLES      
      NWALLS = 2*DIMN

      IF(.NOT.NON_RECT_BC) THEN
                                !particles = npc*(imax)*(jmax)*kmax
         NPARTICLES = PARTICLES * PARTICLES_FACTOR + NWALLS
      ELSE
!     By Tingwen 19/01/2008 6:51:31 PM
!     +2 to include the contact with edge and node
!     only one edge contact or one node contact with wall is allowed for a particle
         NPARTICLES = PARTICLES * PARTICLES_FACTOR + NWALLS + 2 + NWALLS +1
!     by Tingwen 19/01/2008 6:51:31 PM
      END IF
      
      MAXQUADS = 5*PARTICLES*MQUAD_FACTOR
      IF(MAXQUADS.LE.80000) MAXQUADS = 80000
      MAXNEIGHBORS = MN + 1 + NWALLS

      IF(DIMN.EQ.2) THEN
         PBP = PARTICLES/7
      ELSE
         PBP = NPARTICLES/5
      END IF

      IF(DIMN.EQ.3) THEN
         NMQD = 11
      ELSE
         NMQD = 7
      END IF

!     COEFF OF RESITUTIONS 
      ALLOCATE(REAL_EN(MMAX,MMAX),REAL_ET(MMAX,MMAX))
      ALLOCATE(REAL_EN_WALL(MMAX),REAL_ET_WALL(MMAX))
!     
      ALLOCATE(DES_ETAN(MMAX,MMAX))
      ALLOCATE(DES_ETAT(MMAX,MMAX))
      ALLOCATE(DES_ETAN_WALL(MMAX), DES_ETAT_WALL(MMAX))
!     DES Allocatable arrays
!     
!     Particle attributes
!     Radius, density, mass, moment of inertia      
      IF(DES_INTERP_ON) THEN
         ALLOCATE(DRAG_AM(DIMENSION_I-1, DIMENSION_J-1, MAX(1,DIMENSION_K-1), MMAX))
         ALLOCATE(DRAG_BM(DIMENSION_I-1, DIMENSION_J-1,  MAX(1,DIMENSION_K-1), DIMN, MMAX))
         ALLOCATE(WTBAR(DIMENSION_I-1, DIMENSION_J-1,  MAX(1,DIMENSION_K-1),  MMAX))
         ALLOCATE(VEL_FP(NPARTICLES,3))
         ALLOCATE(F_gp(NPARTICLES ))  
         F_gp(1:NPARTICLES)  = ZERO
      ENDIF 
!=======END OF INTERPOLATION RELATED ARRAYS==================!

      ALLOCATE(MARK_PART(NPARTICLES))
      ALLOCATE(BED_HEIGHT(MMAX))
      
      ALLOCATE(AVE_VEL_X(DIMENSION_3,MMAX), AVE_VEL_Y(DIMENSION_3,MMAX), AVE_VEL_Z(DIMENSION_3,MMAX))
      
      IF(DES_NEIGHBOR_SEARCH.EQ.4) THEN
         ALLOCATE(PIC(IMAX2,JMAX2,KMAX2))
         
         DO  k  = 1,KMAX2       !MAX(KMAX1-1,1)
            DO  j  = 1,JMAX2
               DO  i  = 1,IMAX2
                  NULLIFY(pic(i,j,k)%p)
               ENDDO
            ENDDO
         ENDDO
      ENDIF
      
      
      Allocate(  DES_RADIUS (NPARTICLES) )
      Allocate(  RO_Sol (NPARTICLES) )
      
      Allocate(  PVOL (NPARTICLES) )
      Allocate(  PMASS (NPARTICLES) )
      Allocate(  OMOI (NPARTICLES) )
!     
!     Old and new particle positions, velocities (translational and
!     rotational) )      
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
      END IF        
      Allocate(  PPOS (NPARTICLES,DIMN) )
!     
!     Total, normal and tangetial forces      
      Allocate(  FC (NPARTICLES,DIMN) )
      Allocate(  FN (NPARTICLES,DIMN) )
      Allocate(  FT (NPARTICLES,DIMN) )
      Allocate(  FNS2 (DIMN) )
      Allocate(  FTS2 (DIMN) )
      Allocate(  FNS1 (DIMN) )
      Allocate(  FTS1 (DIMN) )
      Allocate(  GRAV (DIMN) )

!     
!     Torque     
      IF(DIMN.EQ.3) THEN 
         Allocate(  TOW (NPARTICLES,DIMN) )
      ELSE
         Allocate(  TOW (NPARTICLES,1) )
      END IF
!     
!     Accumulated spring forces      
      Allocate(  PFN (NPARTICLES,MAXNEIGHBORS,DIMN) )
      Allocate(  PFT (NPARTICLES,MAXNEIGHBORS,DIMN) )
!     
!     Wall position, velocity and normal vector
      Allocate(  DES_WALL_POS (NWALLS,DIMN) )
      Allocate(  DES_WALL_VEL (NWALLS,DIMN) )
      Allocate(  WALL_NORMAL (NWALLS,DIMN) )
      Allocate(  PN (NPARTICLES, MAXNEIGHBORS) )
      Allocate(  PV (NPARTICLES, MAXNEIGHBORS) )
!     
!     Periodic walls
!     
      Allocate(  WWALL (PBP) )
      Allocate(  EWALL (PBP) )
      Allocate(  BWALL (PBP) )
      Allocate(  TWALL (PBP) )
      Allocate(  SWALL (PBP) )
      Allocate(  NWALL (PBP) )
!     
!     Particles in a computational cell (for volume fraction) )
      Allocate(  PINC (DIMENSION_3) )
      Allocate(  PIJK (PARTICLES,5) )
!     
!     Volume averaged solids volume in a cell      
      Allocate(  DES_U_s (DIMENSION_3, DIMENSION_M) )
      Allocate(  DES_V_s (DIMENSION_3, DIMENSION_M) )
      Allocate(  DES_W_s (DIMENSION_3, DIMENSION_M) )

!     Averaged velocity obtained by avraging over all the particles
      ALLOCATE(DES_VEL_AVG(DIMN))

!     Global Granular Energy
      ALLOCATE(GLOBAL_GRAN_ENERGY(DIMN))
      ALLOCATE(GLOBAL_GRAN_TEMP(DIMN))
!     
!     Drag exerted by the gas o solids
      Allocate(  SOLID_DRAG (DIMENSION_3, DIMENSION_M, DIMN) )
!     
!     Neighbor search
      Allocate(  NEIGHBOURS (NPARTICLES, MAXNEIGHBORS) )
      Allocate(  LQUAD (MAXQUADS, NMQD) )
      Allocate(  PQUAD (PARTICLES) )
      Allocate(  CQUAD (MAXQUADS, NWALLS) )
!     
!     Neighbor distances
      Allocate(  PN_DIST (NPARTICLES, MAXNEIGHBORS) )
      Allocate(  PN_RLM (NPARTICLES, MAXNEIGHBORS) )
!     
!     Granular temperature
      Allocate(  DES_THETA (DIMENSION_3, DIMENSION_M) )
!     
!     Cell faces
      Allocate(  XE (DIMENSION_I) )
      Allocate(  YN (DIMENSION_J) )
      Allocate(  ZT (DIMENSION_K) )

!     
!     COHESION      
!     
!     Square-well potential parameters
      Allocate(  WELL_WIDTH (NPARTICLES) )
      Allocate(  WELL_DEPTH (NPARTICLES) )
!     
!     Does particle have at least one linked partner
      Allocate(  IS_LINKED (NPARTICLES) )
!     
!     Does particle have at least one aggloerated partner
      Allocate(  IS_AGGLOMERATED (NPARTICLES) )
!     
!     Array of linked partners
      Allocate(  LINKS (NPARTICLES, MAXNEIGHBORS) )
!     
!     Array of agglomerated partners
      Allocate(  AGGS (NPARTICLES, MAXNEIGHBORS) )
!     
!     Matrix location of particle 
      Allocate(  PART_GRID (NPARTICLES,4) )


      RETURN
      END SUBROUTINE DES_ALLOCATE_ARRAYS 
