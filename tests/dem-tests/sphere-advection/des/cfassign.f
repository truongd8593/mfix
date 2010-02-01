!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
! Module name: CFASSIGN                                                C
!
! Purpose:
! Assign the necessary values for DEM  computation. For example:
! - assigning DEM boundaries from the values entered for MFIX input
!   in mfix.datat
! - assigning DEM gravity vector from MFIX input. 
! - calculating DTSOLID based on particle properties: spring 
!   coefficient, damping factor & mass 
!      
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
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
      LOGICAL:: filexist, isopen
      INTEGER I, J, K, IJK, L, M
      INTEGER COUNT_E
      DOUBLE PRECISION MIN_MASS, MASS_I, MASS_J, &
                       MASS_EFF, RED_MASS_EFF
      DOUBLE PRECISION TCOLL, TCOLL_TMP, MAXMASS
! local variables for calculation of hertzian contact parameters
      DOUBLE PRECISION R_EFF, E_EFF, G_MOD_EFF     
! local variable to determine minimum grid size
      DOUBLE PRECISION MIN_GRID
!-----------------------------------------------
      INCLUDE 'function.inc'
      INCLUDE 'b_force1.inc'
      INCLUDE 'b_force2.inc'


      WRITE(*,'(3X,A)') '---------- START CFASSIGN ---------->'

      PIS = PARTICLES  ! J.Musser 
      MIN_MASS = LARGE_NUMBER
      MAXMASS = SMALL_NUMBER
      MAX_RADIUS = ZERO
      MIN_RADIUS = LARGE_NUMBER
      TCOLL = LARGE_NUMBER
      DO L = 1, PARTICLES
         PVOL(L) = (4.0d0/3.0d0)*Pi*DES_RADIUS(L)**3
         PMASS(L) = PVOL(L)*RO_SOL(L) 
         OMOI(L) = 2.5d0/(PMASS(L)*DES_RADIUS(L)**2) !one over MOI
         MAX_RADIUS = MAX(MAX_RADIUS, DES_RADIUS(L))
         MIN_RADIUS = MIN(MIN_RADIUS, DES_RADIUS(L))
         IF(PMASS(L).LT.MIN_MASS) MIN_MASS = PMASS(L) 
         MAXMASS = MAX(PMASS(L), MAXMASS)
         MARK_PART(L) = 1
         IF(DES_POS_NEW(L,2).LE.YLENGTH/2.d0) MARK_PART(L) = 0
      ENDDO

      RADIUS_EQ = MAX_RADIUS*1.05d0
      WRITE(*,'(5X,A,ES15.8)') '1.05*MAX_RADIUS = ', MAX_RADIUS     

      IF (DES_NEIGHBOR_SEARCH .EQ. 4) THEN
         MIN_GRID = LARGE_NUMBER
         DO IJK = ijkstart3, ijkend3
            I = I_OF(IJK)
            J = J_OF(IJK)
            K = K_OF(IJK)
            MIN_GRID = MIN(MIN_GRID, DX(I))
            MIN_GRID = MIN(MIN_GRID, DY(J))
            IF (DIMN.EQ.3) MIN_GRID = MIN(MIN_GRID, DZ(K))
         ENDDO
         IF (MIN_GRID <= MIN_RADIUS) THEN
            WRITE(*,'(/,5X,A,A,/14X,A,A,/14X,A,A,/)') &
               'WARNING: for grid based search the grid size should ',&
               'be greater than', 'the radius of the smallest ',&
               'particle or neighbor contacts','may be missed ',&
               'giving bad simulation results' 
         ENDIF
      ENDIF


! Set boundary edges 
! In some instances wx1,ex2, etc have been used and in others
! xlength,zero, etc are used.  todo: code should be modified for
! consistency throughout      
      WX1 = ZERO 
      EX2 = XLENGTH 
      BY1 = ZERO
      TY2 = YLENGTH 
      SZ1 = ZERO 
      NZ2 = ZLENGTH

      
      GRAV(1) = BFX_s(1,1)
      GRAV(2) = BFY_s(1,1)
      IF(DIMN.EQ.3) GRAV(3) = BFZ_s(1,1)

! Calculate collision parameters
!--------------------------------------------------------

      IF (TRIM(DES_COLL_MODEL) == 'HERTZIAN') THEN 

         write(*,'(5X,A)') 'COLLISION MODEL: Hertzian'

! particle-particle contact --------------------
         DO I=1,MMAX
            G_MOD(I) = 0.5d0*e_young(I)/(1.d0+v_poisson(I))   ! shear modulus 
            WRITE(*,'(5X,A,I5,X,A,X,2(ES15.7))') &
               'E_YOUNG AND V_POISSON FOR M = ', I, '=',&
               E_YOUNG(I), V_POISSON(I) 
         ENDDO

         COUNT_E = 0
         DO I=1,MMAX
            DO J=I,MMAX
! Arrange the coefficient of restitution matrix from en_input values
! use coef of rest to determine damping coefficient 
               COUNT_E = COUNT_E + 1
               REAL_EN(I,J) = DES_EN_INPUT(COUNT_E)
               REAL_ET(I,J) = DES_ET_INPUT(COUNT_E)            
               MASS_I = (PI/6.d0)*(D_P0(I)**3)*RO_S(I)
               MASS_J = (PI/6.d0)*(D_P0(J)**3)*RO_S(J)
               MASS_EFF = (MASS_I*MASS_J)/(MASS_I+MASS_J)
! In the Hertzian model introduce a factor of 2/7 to the effective mass 
! for tangential direction to get a reduced mass.  Reference: Van der Hoef
! et al., Multi-scale modeling of gas-fluidized beds, Advances in Chemical
! Engineering, 2006
               RED_MASS_EFF = (2.d0/7.d0)*MASS_EFF               
               R_EFF = 0.5d0*(D_P0(I)*D_P0(J)/(D_P0(I)+D_P0(J)))
               E_EFF = e_young(I)*e_young(J)/ &
                  (e_young(I)*(1.d0-v_poisson(J)**2)+&
                   e_young(J)*(1.d0-v_poisson(I)**2))
               G_MOD_EFF = G_MOD(I)*G_MOD(J)/&
                  (G_MOD(I)*(2.d0-v_poisson(J))+&
                   G_MOD(J)*(2.d0-v_poisson(I)))

               hert_kn(I,J)=(4.d0/3.d0)*E_EFF*SQRT(R_EFF)
               hert_kt(I,J)=(16.d0/3.d0)*G_MOD_EFF*SQRT(R_EFF)

               DES_ETAN(I,J) = 2.d0*SQRT(hert_kn(I,J)*MASS_EFF)*&
                  ABS(LOG(REAL_EN(I,J)))
               DES_ETAN(I,J) = DES_ETAN(I,J)/&
                  SQRT(PI*PI + (LOG(REAL_EN(I,J)))**2)
               DES_ETAT(I,J) = 2.d0*SQRT(hert_kt(I,J)*RED_MASS_EFF)*&
                  ABS(LOG(REAL_ET(I,J)))
               DES_ETAT(I,J) = DES_ETAT(I,J)/&
                  SQRT(PI*PI + (LOG(REAL_ET(I,J)))**2) 

               hert_kn(J,I) = hert_kn(I,J)
               hert_kt(J,I) = hert_kt(I,J)

               TCOLL_TMP = PI/SQRT(hert_kn(I,J)/MASS_EFF - ((DES_ETAN(I,J)/MASS_EFF)**2)/4.d0)
               TCOLL = MIN(TCOLL_TMP, TCOLL) 

               WRITE(*,'(5X,A,I5,X,I5,X,A,X,2(ES15.7))') &
                  'KN AND KT FOR PAIR ',&
                   I, J, '=', hert_kn(I,J), hert_kt(I,J)
            ENDDO
         ENDDO

! particle-wall contact --------------------          
         COUNT_E = 0
         DO I = 1, MMAX
            COUNT_E = COUNT_E + 1  
            REAL_EN_WALL(I) = DES_EN_WALL_INPUT(COUNT_E)
            REAL_ET_WALL(I) = DES_ET_WALL_INPUT(COUNT_E)
            MASS_I = (PI/6.d0)*(D_P0(I)**3)*RO_S(I)
            MASS_J = MASS_I
            MASS_EFF = MASS_I
            RED_MASS_EFF = (2.d0/7.d0)*MASS_I
            R_EFF = 0.5d0*D_P0(I)
            E_EFF = e_young(I)*ew_young/ &
               (e_young(I)*(1.d0-vw_poisson**2)+&
                ew_young*(1.d0-v_poisson(I)**2))
            G_MOD_EFF = G_MOD(I)/(2.d0-v_poisson(I))

            hert_kwn(I) =(4.d0/3.d0)*E_EFF*SQRT(R_EFF)
            hert_kwt(I) = (16.d0/3.d0)*G_MOD_EFF*SQRT(R_EFF)    

            DES_ETAN_WALL(I) = 2.d0*SQRT(hert_kwn(I)*MASS_EFF)*&
               ABS(LOG(REAL_EN_WALL(I)))
            DES_ETAN_WALL(I) = DES_ETAN_WALL(I)/&
               SQRT(PI*PI + (LOG(REAL_EN_WALL(I)))**2)
            DES_ETAT_WALL(I) = 2.d0*SQRT(hert_kwt(I)*RED_MASS_EFF)*&
               ABS(LOG(REAL_ET_WALL(I)))
            DES_ETAT_WALL(I) = DES_ETAT_WALL(I)/&
               SQRT(PI*PI + (LOG(REAL_ET_WALL(I)))**2) 

            TCOLL_TMP = PI/SQRT(hert_kwn(i)/MASS_EFF - ((DES_ETAN_WALL(I)/MASS_EFF)**2)/4.d0) 
         ENDDO

      ELSE   ! Linear spring-dashpot model

         write(*,'(5X,A)') &
            'COLLISION MODEL: Linear Spring-Dashpot (default)'

! User's input for KT_FAC and KT_W_FAC will be used, otherwise these values are
! estimated using: Silbert et al, 2001, Physical Review E, vol. 64-5, see page 051302-5
         IF(KT_FAC == UNDEFINED) THEN
            KT = (2.d0/7.d0)*KN
         ELSE
            KT = KT_FAC*KN
         ENDIF
         IF(KT_W_FAC == UNDEFINED) THEN
            KT_W = (2.d0/7.d0)*KN_W
         ELSE
            KT_W = KT_W_FAC*KN_W
         ENDIF
         WRITE(*,'(5X,A,ES17.10,2X,ES15.7)') 'KN AND KT = ', KN, KT

! particle-particle contact --------------------
         COUNT_E = 0
         DO I = 1, MMAX
            DO J = I, MMAX
! Arrange the coefficient of restitution matrix from en_input values
! use coef of rest to determine damping coefficient 
               COUNT_E = COUNT_E + 1
               REAL_EN(I,J) = DES_EN_INPUT(COUNT_E)
               MASS_I = (PI/6.d0)*(D_P0(I)**3.d0)*RO_S(I)
               MASS_J = (PI/6.d0)*(D_P0(J)**3.d0)*RO_S(J)
               MASS_EFF = (MASS_I*MASS_J)/(MASS_I + MASS_J)
               DES_ETAN(I,J) = 2.D0*SQRT(KN*MASS_EFF)*&
                  ABS(LOG(REAL_EN(I,J)))
               DES_ETAN(I,J) = DES_ETAN(I,J)/&
                  SQRT(PI*PI + (LOG(REAL_EN(I,J)))**2)
 
! User's input for DES_ETAT_FAC will be used, otherwise these values are
! estimated using: Silbert et al, 2003, Physics of Fluids, vol. 15-1, see page 3
               IF(DES_ETAT_FAC == UNDEFINED) THEN
                  DES_ETAT(I,J) = HALF*DES_ETAN(I,J)
               ELSE
                  DES_ETAT(I,J) = DES_ETAT_FAC*DES_ETAN(I,J)
               ENDIF

               TCOLL_TMP = PI/SQRT(KN/MASS_EFF - ((DES_ETAN(I,J)/MASS_EFF)**2)/4.d0)
               TCOLL = MIN(TCOLL_TMP, TCOLL)
            ENDDO
         ENDDO

! particle-wall contact --------------------     
         COUNT_E = 0 
         DO I = 1, MMAX
            COUNT_E = COUNT_E + 1  
            REAL_EN_WALL(I) = DES_EN_WALL_INPUT(COUNT_E)
            MASS_I = (PI*(D_P0(I)**3)*RO_S(I))/6.d0
            MASS_J = MASS_I
            MASS_EFF = MASS_I
            DES_ETAN_WALL(I) = 2.d0*SQRT(KN_W*MASS_EFF)*&
               ABS(LOG(REAL_EN_WALL(I)))
            DES_ETAN_WALL(I) = DES_ETAN_WALL(I)/&
               SQRT(PI*PI + (LOG(REAL_EN_WALL(I)))**2)
 
            IF(DES_ETAT_W_FAC == UNDEFINED) THEN
               DES_ETAT_WALL(I) = HALF*DES_ETAN_WALL(I)
            ELSE
               DES_ETAT_WALL(I) = DES_ETAT_W_FAC*DES_ETAN_WALL(I)
            ENDIF

            TCOLL_TMP = PI/SQRT(KN_W/MASS_EFF - ((DES_ETAN_WALL(I)/MASS_EFF)**2.d0)/4.d0)
            !TCOLL = MIN(TCOLL_TMP, TCOLL)
         ENDDO

      ENDIF   ! endif des_coll_model = 'hertzian'
!--------------------------------------------------------


      DO I = 1, MMAX
         DO J = I, MMAX
            REAL_EN(J, I) = REAL_EN(I,J)
            REAL_ET(J, I) = REAL_ET(J,I)
            DES_ETAN(J,I) = DES_ETAN(I,J)
            DES_ETAT(J,I) = DES_ETAT(I,J)
         ENDDO
      ENDDO

      DO I = 1, MMAX
         DO J = I, MMAX
            WRITE(*,'(5X,A,I,2X,I,A,2(ES15.7))') &
               'ETAN AND ETAT FOR PAIR ',&
               I, J, ' = ', DES_ETAN(I,J), DES_ETAT(I,J)
         ENDDO
      ENDDO

      DTSOLID = DT ! TCOLL/50.d0 - Hardwired for the circle convection case
      
      WRITE(*,'(5X,A,E17.10,2X,E17.10)') 'MIN TCOLL AND DTSOLID = ',&
         TCOLL, DTSOLID

      WRITE(*,'(3X,A)') '<---------- END CFASSIGN ----------'

      RETURN
      END SUBROUTINE CFASSIGN
