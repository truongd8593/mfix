!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFTOTALOVERLAPS(L, II, J, VRN, VRT, N_OVERLAP, T_OVERLAP,CHECK_CON)                                                                   C
!  Purpose:  DES - Calculate the total overlap between particles       C
!                                                                      C
!                                                                      C
!  Author:                                            Date:            C
!  Reviewer:Rahul Garg                                DATE: 01=Aug-07  C
!  Comments: Added allocation of interpolation based arrays            C
!            MAXNEIGHBORS  definition now includes walls also          C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C


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
      USE discretelement
      Use indices
      Use geometry
      Use compar
      Use physprop
      IMPLICIT NONE
      
      INTEGER NPARTICLES, I,J,K

      DIMENSION_I   = IMAX3
      DIMENSION_J   = JMAX3
      DIMENSION_K   = KMAX3
      DIMENSION_3   = (kend3-kstart3+1)*(jend3-jstart3+1)*(iend3-istart3+1)
      DIMENSION_3G   = IJKMAX3            
      DIMENSION_3L  = ijksize3_all(myPE)      
      DIMENSION_M   = MAX(1, MMAX)
      DIMENSION_4   = (kend4-kstart4+1)*(jend4-jstart4+1)*(iend4-istart4+1)
      
      !dimensionint_I = iend3-istart3
      !dimensionint_J = jend3-jstart3
      !dimensionint_K = kend3-kstart3
      NWALLS = 2*DIMN
      !particles = npc*(imax)*(jmax)*kmax
      NPARTICLES = PARTICLES * PARTICLES_FACTOR + NWALLS
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
!
!   DES Allocatable arrays
!
!   Particle attributes
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

      ALLOCATE(bed_height(MMAX))

      
      IF(DES_NEIGHBOR_SEARCH.EQ.4) THEN
         ALLOCATE(pic(IMAX2,JMAX2,KMAX2))
         
      DO 10 k  = 1,KMAX2!MAX(KMAX1-1,1)
      DO 10 j  = 1,JMAX2
      DO 10 i  = 1,IMAX2
        NULLIFY(pic(i,j,k)%p)
      10 CONTINUE
      ENDIF
         
      write(*,*) 'Nparticles = ', nparticles, particles, particles_factor, nwalls

      Allocate(  DES_RADIUS (NPARTICLES) )
      Allocate(  RO_Sol (NPARTICLES) )
      Allocate(  PVOL (NPARTICLES) )
      Allocate(  PMASS (NPARTICLES) )
      Allocate(  OMOI (NPARTICLES) )
!
!   Old and new particle positions, velocities (translational and
!                                                             rotational) )      
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
!   Total, normal and tangetial forces      
      Allocate(  FC (NPARTICLES,DIMN) )
      Allocate(  FN (NPARTICLES,DIMN) )
      Allocate(  FT (NPARTICLES,DIMN) )
      Allocate(  FNS2 (DIMN) )
      Allocate(  FTS2 (DIMN) )
      Allocate(  FNS1 (DIMN) )
      Allocate(  FTS1 (DIMN) )
      Allocate(  GRAV (DIMN) )

!
!   Torque     
      IF(DIMN.EQ.3) THEN 
         Allocate(  TOW (NPARTICLES,DIMN) )
      ELSE
         Allocate(  TOW (NPARTICLES,1) )
      END IF
!
!   Accumulated spring forces      
      Allocate(  PFN (NPARTICLES,MAXNEIGHBORS,DIMN) )
      Allocate(  PFT (NPARTICLES,MAXNEIGHBORS,DIMN) )
!
!   Wall position, velocity and normal vector
      Allocate(  DES_WALL_POS (NWALLS,DIMN) )
      Allocate(  DES_WALL_VEL (NWALLS,DIMN) )
      Allocate(  WALL_NORMAL (NWALLS,DIMN) )
      Allocate(  PN (NPARTICLES, MAXNEIGHBORS) )
      Allocate(  PV (NPARTICLES, MAXNEIGHBORS) )
!
!   Periodic walls
!
      Allocate(  WWALL (PBP) )
      Allocate(  EWALL (PBP) )
      Allocate(  BWALL (PBP) )
      Allocate(  TWALL (PBP) )
      Allocate(  SWALL (PBP) )
      Allocate(  NWALL (PBP) )
!
!   Particles in a computational cell (for volume fraction) )
      Allocate(  PINC (DIMENSION_3) )
      Allocate(  PIJK (PARTICLES,5) )
!
!   Volume averaged solids volume in a cell      
      Allocate(  DES_U_s (DIMENSION_3, DIMENSION_M) )
      Allocate(  DES_V_s (DIMENSION_3, DIMENSION_M) )
      Allocate(  DES_W_s (DIMENSION_3, DIMENSION_M) )
!
!   Drag exerted by the gas o solids
      Allocate(  SOLID_DRAG (DIMENSION_3, DIMENSION_M, DIMN) )
!
!   Neighbor search
      Allocate(  NEIGHBOURS (NPARTICLES, MAXNEIGHBORS) )
      Allocate(  LQUAD (MAXQUADS, NMQD) )
      Allocate(  PQUAD (PARTICLES) )
      Allocate(  CQUAD (MAXQUADS, NWALLS) )
!
!   Neighbor distances
      Allocate(  PN_DIST (NPARTICLES, MAXNEIGHBORS) )
      Allocate(  PN_RLM (NPARTICLES, MAXNEIGHBORS) )
!
!   Granular temperature
      Allocate(  DES_THETA (DIMENSION_3, DIMENSION_M) )
!
!   Cell faces
      Allocate(  XE (DIMENSION_I) )
      Allocate(  YN (DIMENSION_J) )
      Allocate(  ZT (DIMENSION_K) )

!
!  COHESION      
!
!  Square-well potential parameters
      Allocate(  WELL_WIDTH (NPARTICLES) )
      Allocate(  WELL_DEPTH (NPARTICLES) )
!
!  Does particle have at least one linked partner
      Allocate(  IS_LINKED (NPARTICLES) )
!
!  Does particle have at least one aggloerated partner
      Allocate(  IS_AGGLOMERATED (NPARTICLES) )
!
!  Array of linked partners
      Allocate(  LINKS (NPARTICLES, MAXNEIGHBORS) )
!
!  Array of agglomerated partners
      Allocate(  AGGS (NPARTICLES, MAXNEIGHBORS) )
!
!  Matrix location of particle 
      Allocate(  PART_GRID (NPARTICLES,4) )


      RETURN
      END SUBROUTINE DES_ALLOCATE_ARRAYS 
