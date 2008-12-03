!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
! Module name: CFASSIGN                                                C
!>
!! Purpose: Assign the necessary values for DEM                         
!! computation. For example, assigning DEM                           
!! boundaries from the values entered for MFIX                        
!! input in mfix.dat. Assigning DEM gravity vector                     
!! from MFIX input. Calculating RADIUS_EQ.                          
!! Calculating DTSOLID based on rotational and translational       
!! constraints                                                    
!<                                                               
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer: Sreekanth Pannala                        Date: 09-Nov-06  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CFASSIGN

      USE discretelement
      USE param
      USE param1
      USE parallel
      USE fldvar
      USE run
      USE geometry
      USE matrix
      USE indices
      USE physprop
      USE drag
      USE constant
      USE compar
      USE sendrecv

      IMPLICIT NONE
      LOGICAL:: filexist, isopen
      
      INTEGER L, IJK, M, I,J, K, COUNT_E
      DOUBLE PRECISION FOUR_BY_THREE, RAD2, MINMASS, MASS_I, MASS_J, MASS_EFF
      DOUBLE PRECISION :: TCOLL, TCOLL_TMP, AVG_MASS, MAXMASS
      
!     
!---------------------------------------------------------------------
!     Assignments
!---------------------------------------------------------------------
!     

      INCLUDE 'b_force1.inc'
      INCLUDE 'b_force2.inc'
      
      FOUR_BY_THREE = 4.0d0/3.0d0
      MINMASS = LARGE_NUMBER
      MAXMASS = SMALL_NUMBER
      MAX_RADIUS = ZERO
      MIN_RADIUS = LARGE_NUMBER
      TCOLL = LARGE_NUMBER
      RMS_RAD = ZERO
      DO L = 1, PARTICLES
         RAD2 = DES_RADIUS(L)**2
         PVOL(L) = FOUR_BY_THREE*Pi*RAD2*DES_RADIUS(L)
         PMASS(L) = PVOL(L)*RO_Sol(L) 
         OMOI(L) = 2.5d0/(PMASS(L)*RAD2) !one over MOI
         MAX_RADIUS = MAX(MAX_RADIUS, DES_RADIUS(L))
         MIN_RADIUS = MIN(MIN_RADIUS, DES_RADIUS(L))
         RMS_RAD = RMS_RAD + DES_RADIUS(L)**2.d0
         IF(PMASS(L).LT.MINMASS) MINMASS = PMASS(L) 
         MAXMASS = MAX(PMASS(L), MAXMASS)
         MARK_PART(L) = 1
         IF(DES_POS_NEW(L,2).LE.YLENGTH/2.d0) MARK_PART(L) = 0
         
      END DO
      RMS_RAD = SQRT(RMS_RAD/PARTICLES)
      AVG_MASS = SUM(PMASS(1:PARTICLES))/PARTICLES
      AVG_RAD = SUM(DES_RADIUS(1:PARTICLES))/PARTICLES
      KT = (2.d0/7.d0)*KN
      KT_W = (2.d0/7.d0)*KN_W

      IF(.NOT.DES_PERIODIC_WALLS) THEN
         DES_PERIODIC_WALLS_X = .FALSE.
         DES_PERIODIC_WALLS_Y = .FALSE.
         DES_PERIODIC_WALLS_Z = .FALSE.
      ENDIF 
      INTX_PER =  DES_PERIODIC_WALLS_X
      INTY_PER =  DES_PERIODIC_WALLS_Y
      INTZ_PER =  DES_PERIODIC_WALLS_Z
      
      RADIUS_EQ = DES_RADIUS(1)*1.05D0
      NEIGHBOR_SEARCH_DIST = NEIGHBOR_SEARCH_RAD_RATIO * RADIUS_EQ
      
!     DTSOLID = DTSOLID_FACTOR*2.0D0*PI*SQRT((MINMASS)/(15*KN)) ! DTs - Rotational Constraint
!     DTSOLID = DTSOLID_FACTOR*2D0*PI*SQRT(MINMASS/(6*KN)) ! DTs - Translational Constraint
!     DTSOLID = pi*SQRT(one/(KN/PMASS(1) - (ETA_DES_N**2)/4.d0))
!     DTSOLID = DTSOLID/50
      
      Print*,'DTSOLID = ', dtsolid
      Print*,'MAX_RADIUS = ', MAX_RADIUS
      WRITE(*,*) 'KN AND KT = ', KN, KT
!     read(*,*)
      WX1 = ZERO 
      EX2 = XLENGTH 
      BY1 = ZERO
      TY2 = YLENGTH 
      IF((DIMN.EQ.2).AND.(COORDINATES == 'CARTESIAN')) THEN
         SZ1 = ZERO 
         NZ2 = ZLENGTH
      ELSE
         SZ1 = ZERO
         NZ2 = ZLENGTH
!     WRITE(*,*) 'XLENGHT =', XLENGTH, YLENGTH, ZLENGTH
      ENDIF

!     IF((DIMN.EQ.2).AND.(COORDINATES == 'CARTESIAN')) THEN
!     DZ(:) = 2D0*RADIUS_EQ
!     END IF

      GRAV(1) = BFX_s(1,1)
      GRAV(2) = BFY_s(1,1)
      IF(DIMN.EQ.3) GRAV(3) = BFZ_s(1,1)

!     ARRANGE THE COEFF OF RESTITUTION MATRIX FROM INPUT EN VALUES
      count_e = 0
      DO I = 1, MMAX
         DO J = I, MMAX
            COUNT_E = COUNT_E + 1
            REAL_EN(I,J) = DES_EN_INPUT(COUNT_E)
            REAL_ET(I,J) = DES_ET_INPUT(COUNT_E)
            MASS_I = (PI*(D_P0(I)**3.d0)*RO_S(I))/6.d0
            MASS_J = (PI*(D_P0(J)**3.d0)*RO_S(J))/6.d0
            MASS_EFF = (MASS_I*MASS_J)/(MASS_I + MASS_J)
            DES_ETAN(I,J) = 2.D0*SQRT(KN*MASS_EFF)*ABS(LOG(REAL_EN(I,J)))
            DES_ETAN(I,J) = DES_ETAN(I,J)/SQRT(PI*PI + (LOG(REAL_EN(I,J)))**2.0)
            DES_ETAT(I,J) = HALF*DES_ETAN(I,J)

            TCOLL_TMP = PI/SQRT(KN/MASS_EFF - ((DES_ETAN(I,J)/MASS_EFF)**2.d0)/4.d0)
            
            TCOLL = MIN(TCOLL_TMP, TCOLL)
         ENDDO
      ENDDO
      
      COUNT_E = 0 
      DO I = 1, MMAX
         COUNT_E = COUNT_E + 1  
         REAL_EN_WALL(I) = DES_EN_WALL_INPUT(COUNT_E)
         REAL_ET_WALL(I) = DES_ET_WALL_INPUT(COUNT_E)
         MASS_I = (PI*(D_P0(I)**3.d0)*RO_S(I))/6.d0
         MASS_J = MASS_I
!     MASS_EFF = (MASS_I*MASS_J)/(MASS_I + MASS_J)
         
         MASS_EFF = MASS_I
         DES_ETAN_WALL(I) = 2.d0*SQRT(KN_W*MASS_EFF)*ABS(LOG(REAL_EN_WALL(I)))
         DES_ETAN_WALL(I) = DES_ETAN_WALL(I)/SQRT(PI*PI + (LOG(REAL_EN_WALL(I)))**2.0)
         DES_ETAT_WALL(I) = HALF*DES_ETAN_WALL(I)

         TCOLL_TMP = PI/SQRT(KN_W/MASS_EFF - ((DES_ETAN_WALL(I)/MASS_EFF)**2.d0)/4.d0)
         
!     TCOLL = MIN(TCOLL_TMP, TCOLL)
         
      ENDDO

      DO I = 1, MMAX
         DO J = I, MMAX
            REAL_EN(J, I) = REAL_EN(I,J)
            REAL_ET(J, I) = REAL_ET(I,J)
            DES_ETAN(J,I) = DES_ETAN(I,J)
            DES_ETAT(J,I) = DES_ETAT(I,J)
         ENDDO
      ENDDO

      DO I = 1, MMAX
         DO J = 1, MMAX
            WRITE(*,*) 'I AND J = ', I, J
            WRITE(*,'(A,i2,1x,i2,A,2(g17.8))') 'ETA_N AND ETA_T  FOR PAIR',I,J,' = ', DES_ETAN(I,J), DES_ETAT(I,J)
         ENDDO
      ENDDO
      
      DTSOLID = TCOLL/50.d0
      
      WRITE(*,*) 'MIN TCOLL AND DTSOLID = ', TCOLL, DTSOLID
      
      RETURN
      END SUBROUTINE CFASSIGN
