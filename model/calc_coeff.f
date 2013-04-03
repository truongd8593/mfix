!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_COEFF_ALL                                          C
!  Purpose: Sets logical values of local flags that in turn tells the  C
!           code whether to perform the indicated type of calculation  C
!           if the value of the flag is true.  Specifically this       C
!           routine directs calculation of all physical and transport  C
!           properties, reaction rates and exchange rates.             C
!                                                                      C
!  Author: M. Syamlal                                 Date: 25-AUG-05  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CALC_COEFF_ALL(FLAG, IER) 

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param 
      USE param1 
      USE physprop
      USE rxns
      USE funits 
      USE compar
      USE ur_facs 
      USE run

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! FLAG = 0, overwrite the coeff arrays, as for example at
!           the beginning of a time step
! FLAG = 1, do not overwrite
      INTEGER, INTENT(IN) :: FLAG
! Error index
      INTEGER, INTENT(INOUT) :: IER
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Flags to tell whether to calculate or not
      LOGICAL :: DENSITY(0:DIMENSION_M), PSIZE(0:DIMENSION_M),&
                 SP_HEAT(0:DIMENSION_M)
! Flags to tell whether to calculate or not
      LOGICAL :: VISC(0:DIMENSION_M), COND(0:DIMENSION_M),&
                 DIFF(0:DIMENSION_M)
! Flag for exchange functions
      LOGICAL :: DRAGCOEF(0:DIMENSION_M, 0:DIMENSION_M),&
                 HEAT_TR(0:DIMENSION_M, 0:DIMENSION_M),&
                 WALL_TR
! Flags to tell whether to calculate or not:
! a granular energy dissipation term
      LOGICAL :: GRAN_DISS(0:DIMENSION_M)
! Temporary storage 
      DOUBLE PRECISION :: UR_F_gstmp, UR_kth_smltmp
!-----------------------------------------------

! First set any under relaxation factors for coefficient calculations
! to 1 and change them to user defined values after the call to
! calc_coefficient
      IF(FLAG == 0) THEN
        UR_F_gstmp = UR_F_gs
        UR_F_gs = ONE
        UR_Kth_smltmp = UR_Kth_sml
        UR_Kth_sml = ONE        
      ENDIF


      CALL TurnOffCOEFF(DENSITY, PSIZE, SP_HEAT, VISC, COND, &
           DIFF, GRAN_DISS, DRAGCOEF, HEAT_TR, WALL_TR, IER) 

      IF(CALL_DQMOM)                          &! Used by PHYSICAL_PROP
         PSIZE(:MMAX) = .TRUE.

      WALL_TR = .TRUE.                         ! Used by EXCHANGE

      IF (ENERGY_EQ) THEN 
         SP_HEAT(:MMAX) = .TRUE.               ! Used by PHYSICAL_PROP
         COND(:MMAX) = .TRUE.                  ! Used by TRANSPORT_PROP
         HEAT_TR(:MMAX,:MMAX) = .TRUE.         ! Used by EXCHANGE
      ENDIF 

      IF (ANY_SPECIES_EQ)                     &! Used by TRANSPORT_PROP
         DIFF(:MMAX) = .TRUE. 

      DRAGCOEF(:MMAX,:MMAX) = .TRUE.           ! Used by EXCHANGE

      IF(RO_G0 == UNDEFINED)                  &! Used in PHYSICAL_PROP
         DENSITY(0) = .TRUE. 

      IF (TRIM(KT_TYPE) .EQ. 'IA_NONEP'       &! Used by TRANSPORT_PROP 
         .OR. TRIM(KT_TYPE) .EQ. 'GD_99')     &
         GRAN_DISS(:MMAX) = .TRUE.

! Viscosity for gas phase is set differently in iterate.
      VISC(0) = .TRUE.                         ! Used by TRANSPORT_PROP
      IF (MU_S0 == UNDEFINED)                 &! Used by TRANSPORT_PROP
         VISC(1:MMAX) = .TRUE. 

      CALL CALC_COEFF (DENSITY, PSIZE, SP_HEAT, VISC, COND, DIFF, &
          GRAN_DISS, DRAGCOEF, HEAT_TR, WALL_TR, IER)

! Now restore all underrelaxation factors for coefficient
! calculations to their original value
      IF(FLAG == 0) THEN
        UR_F_gs = UR_F_gstmp
        UR_Kth_sml = UR_Kth_smltmp
      ENDIF
     
      RETURN  
      END SUBROUTINE CALC_COEFF_ALL



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_COEFF                                              C
!  Purpose: Calculate physical and transport properties, reaction      C
!           rates and exchange rates as directed by the respective     C
!           flags                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 23-APR-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CALC_COEFF(DENSITY, PSIZE, SP_HEAT, VISC, COND, DIFF, &
            GRAN_DISS, DRAG, HEAT_TR, WALL_TR, IER) 

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param 
      USE param1 
      USE physprop
      USE rxns
      USE funits 
      USE compar
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! density, particle size, specific heat
      LOGICAL, INTENT(INOUT) :: DENSITY(0:DIMENSION_M),&
                                PSIZE(0:DIMENSION_M),&
                                SP_HEAT(0:DIMENSION_M)
! vicosity, conductivity, diffusivity
      LOGICAL, INTENT(INOUT) :: VISC(0:DIMENSION_M), &
                                COND(0:DIMENSION_M),&
                                DIFF(0:DIMENSION_M)
! exchange functions
      LOGICAL, INTENT(INOUT) :: DRAG(0:DIMENSION_M, 0:DIMENSION_M),&
                                HEAT_TR(0:DIMENSION_M, 0:DIMENSION_M),&
                                WALL_TR
! a granular energy dissipation term
      LOGICAL, INTENT(INOUT) :: GRAN_DISS(0:DIMENSION_M)
! error index
      INTEGER, INTENT(INOUT) :: IER
!-----------------------------------------------
! Local variables      
!-----------------------------------------------
! Loop indices
      INTEGER :: L, M
!-----------------------------------------------

! Calculate physical properties
      CALL PHYSICAL_PROP (DENSITY, PSIZE, SP_HEAT, IER) 

! Calculate Transport properties
      CALL TRANSPORT_PROP (VISC, COND, DIFF, GRAN_DISS, IER) 

! Calculate reaction rates and interphase mass transfer
      CALL CALC_RRATE

! Calculate interphase momentum, and energy transfers
      CALL EXCHANGE (DRAG, HEAT_TR, WALL_TR, IER) 

! Reset all flags.  The flags need to be set every time this routine is
! called.
      CALL TurnOffCOEFF(DENSITY, PSIZE, SP_HEAT, VISC, COND, DIFF, &
           GRAN_DISS, DRAG, HEAT_TR, WALL_TR, IER)

      RETURN  
      END SUBROUTINE CALC_COEFF 



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_RRATE                                              C
!  Purpose: if rrate then calculate reaction rates and interphase      C
!           mass transfer. if present, calculate discrete reactions    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_RRATE

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE rxns,           only : RRATE, USE_RRATES
      USE funits,         only : DMP_LOG, UNIT_LOG 
      USE compar,         only : myPE
      USE discretelement, only : DISCRETE_ELEMENT
      USE des_rxns,       only : ANY_DES_SPECIES_EQ

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! error index
      INTEGER :: IER

! For DEM simulations that do not have a homogeneous gas phase reaction,
! the gas phase arrays phase change arrays need to be cleared in 
! CALC_RRATE_DES. Other
      LOGICAL CLEAR_ARRAYS

!-----------------------------------------------



      CLEAR_ARRAYS = .TRUE.
! Calculate reaction rates and interphase mass transfer
      IF(RRATE) THEN
! Legacy hook: Calculate reactions from rrates.f.
         IF(USE_RRATES) THEN
            CALL RRATES (IER)
            IF(IER .EQ. 1) THEN
               CALL START_LOG
               IF(DMP_LOG) WRITE (UNIT_LOG, 1000)
               CALL END_LOG 
               CALL MFIX_EXIT(myPE)  
            ENDIF
         ELSE
            CALL RRATES0 (IER)
            CLEAR_ARRAYS = .FALSE.
         ENDIF
      ENDIF 

      IF(DISCRETE_ELEMENT .AND. ANY_DES_SPECIES_EQ) &
         CALL CALC_RRATE_DES(CLEAR_ARRAYS)

      RETURN  
 1000 FORMAT(/1X,70('*')//' From: CALC_COEFF',/,&
         ' Species balance equations are being solved; but chemical',/,&
         ' reactions are not specified in mfix.dat or in rrates.f.',/,&
         ' Copy the file mfix/model/rrates.f into the run directory ',/,&
         ' and remove the initial section that returns IER=1.'&
         ,/1X,70('*')/) 

      END SUBROUTINE CALC_RRATE 



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: TurnoffCoeff                                            C
!  Purpose: Reset all flags to false                                   C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE TurnOffCOEFF(DENSITY, PSIZE, SP_HEAT, VISC, COND, &
           DIFF, GRAN_DISS, DRAG, HEAT_TR, WALL_TR, IER) 

!----------------------------------------------- 
! Modules 
!-----------------------------------------------
      USE param 
      USE param1 
      USE physprop
      IMPLICIT NONE
!----------------------------------------------- 
! Dummy arguments
!-----------------------------------------------
! density, particle size distribution, specific heat      
      LOGICAL, INTENT(INOUT) :: DENSITY(0:DIMENSION_M),&
                                PSIZE(0:DIMENSION_M),&
                                SP_HEAT(0:DIMENSION_M)
! viscosity, conductivity, diffusivity                                
      LOGICAL, INTENT(INOUT) :: VISC(0:DIMENSION_M), &
                                COND(0:DIMENSION_M),&
                                DIFF(0:DIMENSION_M)
! exchange functions
      LOGICAL, INTENT(INOUT) :: DRAG(0:DIMENSION_M, 0:DIMENSION_M),&
                                HEAT_TR(0:DIMENSION_M, 0:DIMENSION_M),&
                                WALL_TR
! a granular energy dissipation term
      LOGICAL, INTENT(INOUT) :: GRAN_DISS(0:DIMENSION_M)
! Error index
      INTEGER, INTENT(INOUT) :: IER
!----------------------------------------------- 
! Variables/Arguments
!----------------------------------------------- 
! Loop indices
      INTEGER :: L, M
!-----------------------------------------------

! Reset all flags
      WALL_TR = .FALSE. 
      DENSITY(:MMAX) = .FALSE. 
      PSIZE(:MMAX) = .FALSE. 
      SP_HEAT(:MMAX) = .FALSE. 
      VISC(:MMAX) = .FALSE. 
      COND(:MMAX) = .FALSE. 
      DIFF(:MMAX) = .FALSE.
      GRAN_DISS(:MMAX) = .FALSE.  
      DRAG(:MMAX,:MMAX) = .FALSE. 
      HEAT_TR(:MMAX,:MMAX) = .FALSE. 

      RETURN  
      END SUBROUTINE TurnOffCOEFF 
