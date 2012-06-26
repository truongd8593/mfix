!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: DES_GRANULAR_TEMPERATURE                                C
!  Purpose: Calculate the DES granular temperature                     C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!  Revision: For parallel processing modifications are made to         C
!            accomodate ghost cells  (Pradeep G.)                      C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE DES_GRANULAR_TEMPERATURE

!-----------------------------------------------
! Modules
!-----------------------------------------------       
      USE discretelement
      USE param
      USE param1
      USE run
      USE geometry
      USE indices
      USE compar
      USE sendrecv
      USE physprop
      IMPLICIT NONE
!-----------------------------------------------
! Local Variables
!-----------------------------------------------      
! indices
      INTEGER :: I, J, K, IJK
! 
      INTEGER :: M, LL
! counter for no. of particles in phase m in cell ijk 
      INTEGER :: NP_PHASE(DIMENSION_3, DES_MMAX)
! temporary variable for mth phase granular temperature in cell ijk
      DOUBLE PRECISION :: TEMP(DIMENSION_3, DES_MMAX)
! accounted for particles
      INTEGER :: PC             
! squared particle velocity v.v
      DOUBLE PRECISION :: SQR_VEL
!-----------------------------------------------
! Include statement functions
!----------------------------------------------- 
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'fun_avg2.inc'
!-----------------------------------------------

      IF (DES_CALC_BEDHEIGHT) CALL CALC_DES_BEDHEIGHT      

! Calculate a local species granular temperature for current instant of
! time.  Note that the following calculation of species granular
! temperature employs a fluctuating particle velocity that is defined
! as the difference between a particles velocity and the corresponding
! local mean velocity of that particles species evaluated at the same
! instant of time.
!----------------------------------------------------------------->>> 
! The following calculations are performed on the 'fluid' grid
      TEMP(:,:) = ZERO
      NP_PHASE(:,:) = ZERO
      PC = 0
      DO LL = 1, MAX_PIP
! skipping particles that do not exist
         IF(.NOT.PEA(LL,1)) CYCLE
         PC = PC + 1
! skipping ghost particles         
         IF(PEA(LL,4)) CYCLE    

         I = PIJK(LL,1)
         J = PIJK(LL,2)
         K = PIJK(LL,3)  
         IJK = PIJK(LL,4)

         M = PIJK(LL,5)
         NP_PHASE(IJK,M) = NP_PHASE(IJK,M) + 1

         TEMP(IJK,M) = TEMP(IJK,M) + &
            (DES_VEL_NEW(LL,1)-DES_U_s(IJK,M))**2 
         TEMP(IJK,M) = TEMP(IJK,M) + &
            (DES_VEL_NEW(LL,2)-DES_V_s(IJK,M))**2
         IF(DIMN.EQ.3) THEN 
            TEMP(IJK,M) = TEMP(IJK,M) + &
               (DES_VEL_NEW(LL,3)-DES_W_s(IJK,M))**2 
         ENDIF

         IF(PC .EQ. PIP) EXIT
      ENDDO

! loop over all fluid cells      
      DO IJK = IJKSTART3, IJKEND3
         IF(FLUID_AT(IJK)) THEN

            DO M = 1,DES_MMAX
               IF (NP_PHASE(IJK,M) > 0 ) THEN
                  DES_THETA(IJK,M) = TEMP(IJK,M)/&
                     DBLE(DIMN*NP_PHASE(IJK,M))
               ELSE
                  DES_THETA(IJK,M) = ZERO
               ENDIF 
            ENDDO
         ENDIF
      ENDDO
!-----------------------------------------------------------------<<<


! Calculate global quantities: granular temperature, kinetic energy,
! potential energy and average velocity at the current instant of time
!----------------------------------------------------------------->>>

! initialization for calculations
      DES_KE = ZERO
      DES_PE = ZERO 
      DES_VEL_AVG(:) = ZERO

! Calculate global average velocity, kinetic energy &
! potential energy
      PC = 0
      DO LL = 1, MAX_PIP
         IF(PEA(LL,1)) PC = PC + 1
! skipping ghost particles and particles that don't exist
         IF(.NOT.PEA(LL,1) .OR. PEA(LL,4)) CYCLE    

         SQR_VEL = ZERO
         DO I = 1, DIMN
            SQR_VEL = SQR_VEL + DES_VEL_NEW(LL,I)**2
         ENDDO

         DES_KE = DES_KE + PMASS(LL)/2.d0 * SQR_VEL 
         DES_PE = DES_PE + PMASS(LL)*DBLE(ABS(GRAV(2)))*&
            DES_POS_NEW(LL,2)
         DES_VEL_AVG(:) =  DES_VEL_AVG(:) + DES_VEL_NEW(LL,:)

         IF(PC .EQ. PIP) EXIT
      ENDDO

!J.Musser changed PARTICLES TO PIP 
      IF(PIP > 0) DES_VEL_AVG(:) = DES_VEL_AVG(:)/DBLE(PIP)

! The following quantities are primarily used for debugging/developing
! and allow a quick check of the energy conservation in the system.
! In their current form they are best applied to monodisperse cases. 
! Calculate x,y,z components of global energy & granular temperature
      GLOBAL_GRAN_ENERGY(:) = ZERO
      GLOBAL_GRAN_TEMP(:)  = ZERO
      PC = 0
      DO LL = 1, MAX_PIP
         IF(PEA(LL,1)) PC = PC + 1
! skipping ghost particles and particles that don't exist
         IF(.NOT.PEA(LL,1) .OR. PEA(LL,4)) CYCLE    

         GLOBAL_GRAN_ENERGY(:) = GLOBAL_GRAN_ENERGY(:) + &
            0.5d0*PMASS(LL)*(DES_VEL_NEW(LL,:)-DES_VEL_AVG(:))**2
         GLOBAL_GRAN_TEMP(:) = GLOBAL_GRAN_TEMP(:) + &
            (DES_VEL_NEW(LL,:)-DES_VEL_AVG(:))**2

         IF(PC .EQ. PIP) EXIT
      ENDDO

      IF(PIP > 0) GLOBAL_GRAN_ENERGY(:) =  GLOBAL_GRAN_ENERGY(:)/DBLE(PIP)
      IF(PIP > 0) GLOBAL_GRAN_TEMP(:) =  GLOBAL_GRAN_TEMP(:)/DBLE(PIP)
!-----------------------------------------------------------------<<<

      RETURN

      END SUBROUTINE DES_GRANULAR_TEMPERATURE



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CALC_DES_BEDHEIGHT                                      !
!  Purpose: Calculate an average bed height for each solids phase      !
!           present                                                    !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE CALC_DES_BEDHEIGHT

!-----------------------------------------------
! Modules
!-----------------------------------------------      
      USE indices
      USE geometry      
      USE compar
      USE discretelement
      USE des_bc
      USE physprop
      USE fldvar
      IMPLICIT NONE
!-----------------------------------------------
! Local Variables
!-----------------------------------------------  
! particle index
      INTEGER :: L
! accounted for particles
      INTEGER :: PC
! solids phase no.    
      INTEGER :: M      
! ijk indices      
      INTEGER :: I, J, K, IJK
! average height of fluid cell 
      DOUBLE PRECISION :: hcell
! height of particle (y-position)
      DOUBLE PRECISION :: hpart      
! volume fraction of phase M in fluid cell
      DOUBLE PRECISION :: EP_SM      
! tmp variables for calculations
      DOUBLE PRECISION :: tmp_num(DES_MMAX), tmp_den(DES_MMAX)
!-----------------------------------------------
! Include statement functions      
!-----------------------------------------------
      INCLUDE 'function.inc'
!-----------------------------------------------

! calculation of bed height following the formulation of Goldschmidt et
! al., Powder Technology 138, 2003, 135-159
      tmp_num(:) = ZERO 
      tmp_den(:) = ZERO 
      DO IJK = IJKSTART3, IJKEND3
         J = J_OF(IJK)      
         DO M = 1, DES_MMAX      
            IF(DES_ROP_S(IJK,M) > ZERO) THEN
               hcell = 0.5d0*(YN(J)+YN(J-1))
               EP_SM = DES_ROP_S(IJK,M)/DES_RO_S(M)
               tmp_num(M) = tmp_num(M) + EP_SM*hcell*VOL(IJK)
               tmp_den(M) = tmp_den(M) + EP_SM*VOL(IJK)
            ENDIF
         ENDDO
      ENDDO
! calculate avg height for each phase
      IF (PIP >0) bed_height(:) = tmp_num(:)/tmp_den(:)


! alternative method to calculating bed height (turned off atm)
      IF(.FALSE.) THEN 
      tmp_num(:) = ZERO 
      tmp_den(:) = ZERO 

      PC = 0
      DO L = 1, MAX_PIP
         IF(PEA(L,1)) PC = PC + 1
! skipping ghost particles and particles that don't exist
         IF(.NOT.PEA(L,1) .OR. PEA(L,4)) CYCLE
         M = PIJK(L,5)
         hpart = DES_POS_NEW(L,2)
         tmp_num(M) = tmp_num(M) + hpart
         tmp_den(M) = tmp_den(M) + 1
         IF(PC .EQ. PIP) EXIT
      ENDDO
      ! calculate avg height for each phase
      IF (PIP >0) bed_height(:) = tmp_num(:)/tmp_den(:)
      ENDIF

      RETURN

      END SUBROUTINE CALC_DES_BEDHEIGHT      

