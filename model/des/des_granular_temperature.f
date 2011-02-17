!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: DES_GRANULAR_TEMPERATURE
!  Purpose: DES - Calculate the DES granular temperature               
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE DES_GRANULAR_TEMPERATURE

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
! indices
      INTEGER I, J, K, IJK
! 
      INTEGER M, LL
! counter for no. of particles in phase m in cell ijk 
      INTEGER NP_PHASE(DIMENSION_3, MMAX)
! temporary variable for mth phase granular temperature in cell ijk
      DOUBLE PRECISION TEMP(DIMENSION_3, MMAX)
! accounted for particles
      INTEGER PC             
! squared particle velocity v.v
      DOUBLE PRECISION SQR_VEL
!-----------------------------------------------      

      INCLUDE 'function.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'fun_avg2.inc'


! Calculate a local species granular temperature for current instant of
! time.  Note that the following calculation of species granular
! temperature employs a fluctuating particle velocity that is defined
! as the difference between a particles velocity and the corresponding
! local mean velocity of that particles species evaluated at the same
! instant of time.
!----------------------------------------------- 
! The following calculations are performed on the 'fluid' grid
      TEMP(:,:) = ZERO
      NP_PHASE(:,:) = ZERO
      PC = 1
      DO LL = 1, MAX_PIS
         IF(PC .GT. PIS) EXIT
         IF(.NOT.PEA(LL,1)) CYCLE    

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

         PC = PC + 1
      ENDDO

! loop over all fluid cells      
      DO IJK = IJKSTART3, IJKEND3
         IF(FLUID_AT(IJK)) THEN

            DO M = 1,MMAX
               IF (NP_PHASE(IJK,M) > 0 ) THEN
                  DES_THETA(IJK,M) = TEMP(IJK,M)/&
                     DBLE(DIMN*NP_PHASE(IJK,M))
               ELSE
                  DES_THETA(IJK,M) = ZERO
               ENDIF 
            ENDDO
         ENDIF
      ENDDO


! Calculate global quantities: granular temperature, kinetic energy,
! potential energy and average velocity at the current instant of time
!-----------------------------------------------

! initialization for calculations
      DES_KE = ZERO
      DES_PE = ZERO 
      DES_VEL_AVG(:) = ZERO

! Calculate global average velocity, kinetic energy &
! potential energy
      PC = 1
      DO LL = 1, MAX_PIS
         IF(PC .GT. PIS) EXIT
         IF(.NOT.PEA(LL,1)) CYCLE

         SQR_VEL = ZERO
         DO I = 1, DIMN
            SQR_VEL = SQR_VEL + DES_VEL_NEW(LL,I)**2
         ENDDO

         DES_KE = DES_KE + PMASS(LL)/2.d0 * SQR_VEL 
         DES_PE = DES_PE + PMASS(LL)*DBLE(ABS(GRAV(2)))*&
            DES_POS_NEW(LL,2)
         DES_VEL_AVG(:) =  DES_VEL_AVG(:) + DES_VEL_NEW(LL,:)

         PC = PC + 1
      ENDDO

!J.Musser changed PARTICLES TO PIS
      IF (PIS > 0) DES_VEL_AVG(:) = DES_VEL_AVG(:)/DBLE(PIS)

! The following quantities are primarily used for debugging/developing
! and allow a quick check of the energy conservation in the system.
! In their current form they are best applied to monodisperse cases. 
! Calculate x,y,z components of global energy & granular temperature
      GLOBAL_GRAN_ENERGY(:) = ZERO
      GLOBAL_GRAN_TEMP(:)  = ZERO
      PC = 1
      DO LL = 1, MAX_PIS
         IF(PC .GT. PIS) EXIT
         IF(.NOT.PEA(LL,1)) CYCLE

         GLOBAL_GRAN_ENERGY(:) = GLOBAL_GRAN_ENERGY(:) + &
            0.5d0*PMASS(LL)*(DES_VEL_NEW(LL,:)-DES_VEL_AVG(:))**2
         GLOBAL_GRAN_TEMP(:) = GLOBAL_GRAN_TEMP(:) + &
            (DES_VEL_NEW(LL,:)-DES_VEL_AVG(:))**2

         PC = PC + 1
      ENDDO

      IF (PIS > 0) GLOBAL_GRAN_ENERGY(:) =  GLOBAL_GRAN_ENERGY(:)/DBLE(PIS)
      IF (PIS > 0) GLOBAL_GRAN_TEMP(:) =  GLOBAL_GRAN_TEMP(:)/DBLE(PIS)

      RETURN

      END SUBROUTINE DES_GRANULAR_TEMPERATURE

