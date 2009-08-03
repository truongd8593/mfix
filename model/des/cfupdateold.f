!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFUPDATEOLD                                            C
!>
!!  Purpose: DES - Update arrays                                        
!<
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CFUPDATEOLD
      
      USE discretelement
      USE run

      IMPLICIT NONE

!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER LL
! Accounted for particles
      INTEGER PC             

!-----------------------------------------------


      PC = 1
      DO LL = 1, MAX_PIS
         IF(PC .GT. PIS) EXIT
         IF(.NOT.PEA(LL)) CYCLE

         DES_POS_OLD(LL,:) = DES_POS_NEW(LL,:)
         DES_VEL_OLD(LL,:) = DES_VEL_NEW(LL,:)
         OMEGA_OLD(LL,:) = OMEGA_NEW(LL,:)

         PC = PC + 1
      END DO

!J.Musser changed PARTICLES TO PIS     
      DES_VEL_AVG(:) = DES_VEL_AVG(:)/PIS       
      GLOBAL_GRAN_ENERGY = ZERO
      GLOBAL_GRAN_TEMP  = ZERO

      PC = 1
      DO LL = 1, MAX_PIS
         IF(PC .GT. PIS) EXIT
         IF(.NOT.PEA(LL)) CYCLE

         GLOBAL_GRAN_ENERGY(:) = GLOBAL_GRAN_ENERGY(:) + &
            PMASS(LL)*(DES_VEL_NEW(LL,:)-DES_VEL_AVG(:))**2.d0
         GLOBAL_GRAN_TEMP(:) = GLOBAL_GRAN_TEMP(:) + &
            (DES_VEL_NEW(LL,:)-DES_VEL_AVG(:))**2.d0

         PC = PC + 1
      ENDDO

!J.Musser changed PARTICLES TO PIS
      GLOBAL_GRAN_ENERGY(:) =  GLOBAL_GRAN_ENERGY(:)/PIS 
      GLOBAL_GRAN_TEMP(:) =  GLOBAL_GRAN_TEMP(:)/PIS     
!     WRITE(des_extra_unit,3022) S_TIME, NEIGH_MAX, OVERLAP_MAX,  SUM(GLOBAL_GRAN_ENERGY(1:DIMN))*HALF

      RETURN
 3022 FORMAT(g12.5,2x,i4,2x,2(g12.5,2x))

      END SUBROUTINE CFUPDATEOLD


