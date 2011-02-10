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
! the number of particles in the system
      INTEGER NPARTICLES
!-----------------------------------------------
      INCLUDE 'function.inc'


      WRITE(*,'(1X,A)')&
         '---------- START DES_ALLOCATE_ARRAYS ---------->'


! Further initialization         
!-----------------------------------------------         
      NWALLS = 2*DIMN

      IF(.NOT.NON_RECT_BC) THEN
! +nwalls is included since calc_force_des temporarily uses the variables 
! pos, vel, etc at elements beyond the array index given by particles. 
! unclear why additional array space is needed via particles_factor
         NPARTICLES = PARTICLES * PARTICLES_FACTOR + NWALLS
      ELSE
! T. Li : 19/01/2008 
! +2 to include the contact with edge and node
! only one edge contact or one node contact with wall is allowed for a particle
         NPARTICLES = PARTICLES * PARTICLES_FACTOR + NWALLS + 2 + NWALLS +1
      ENDIF

! J.Musser : Dynamic Particle Info
      IF(MAX_PIS /= UNDEFINED_I .AND. &
         MAX_PIS .GT. NPARTICLES) NPARTICLES = MAX_PIS

      MAXNEIGHBORS = MN + 1 + NWALLS


      IF (DES_NEIGHBOR_SEARCH.EQ.2 .OR. DES_NEIGHBOR_SEARCH.EQ.3) THEN
         MAXQUADS = 5*PARTICLES*MQUAD_FACTOR   
         IF(MAXQUADS.LE.80000) MAXQUADS = 80000 
         IF(DIMN.EQ.3) THEN   
            NMQD = 11   
         ELSE   
            NMQD = 7   
         ENDIF 
      ENDIF 
!-----------------------------------------------




! Start allocate DES arrays
!-----------------------------------------------
! J.Musser: BC and dynamic particle info
      ALLOCATE( PEA (NPARTICLES, 3) )

! J.Musser: Allocate necessary arrays for discrete mass inlets      
      IF(DES_BCMI /= 0 .OR. DES_BCMO /=0) CALL ALLOCATE_DES_MIO
 
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
! Variable that stores the particle in cell information (ID) on the
! computational grid defined by cell/grid based search.  Similar to the
! variable PIC but tailored for the grid based neighbor search option
         ALLOCATE(DESGRIDSEARCH_PIC(DESGS_IMAX2,DESGS_JMAX2,DESGS_KMAX2))
         DO K = 1,DESGS_KMAX2
            DO J = 1,DESGS_JMAX2
               DO I = 1,DESGS_IMAX2
                  NULLIFY(DESGRIDSEARCH_PIC(I,J,K)%p) 
               ENDDO 
             ENDDO 
          ENDDO                      
      ENDIF   ! end if des_neighbor_search == 4


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



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: ALLOCATE_DES_MIO                                       !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!  Author: J.Musser                                   Date: 17-Aug-09  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE ALLOCATE_DES_MIO

      USE des_bc
      USE discretelement

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER I     ! Loop counter for no. of DES_BCMI

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

