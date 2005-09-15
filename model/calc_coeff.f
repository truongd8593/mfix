!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_COEFF_ALL(FLAG, IER)                              C
!  Purpose: Calculate all physical and transport properties, 
!           reaction rates and exchange rates.                         C
!  Author: M. Syamlal                                 Date: 25-AUG-05  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CALC_COEFF_ALL(FLAG, IER) 
!-----------------------------------------------
!   M o d u l e s 
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
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      FLAG = 0, overwrite the coeff arrays, as for example at
!                              the beginning of a time step
!                      FLAG = 1, do not overwrite
      INTEGER          FLAG
!
!                      Error index
      INTEGER          IER
!
!                      Flags to tell whether to calculate or not
      LOGICAL          DENSITY(0:DIMENSION_M), PSIZE(0:DIMENSION_M),&
                       SP_HEAT(0:DIMENSION_M)
!
!                      Flags to tell whether to calculate or not
      LOGICAL          VISC(0:DIMENSION_M), COND(0:DIMENSION_M),&
                       DIFF(0:DIMENSION_M)
!
!                      Flag for Reaction rates
      LOGICAL          RRATE
!
!                      Flag for exchange functions
      LOGICAL          DRAGCOEF(0:DIMENSION_M, 0:DIMENSION_M),&
                       HEAT_TR(0:DIMENSION_M, 0:DIMENSION_M),&
                       WALL_TR
!     
!     Temporary storage 
      DOUBLE PRECISION UR_F_gstmp
!-----------------------------------------------
      
      
      
      IF(FLAG == 0) then
        ! First set any under relaxation factors for coefficient calculations
	! to 1 and change them to user defined values after the call to
	! calc_coefficient

        UR_F_gstmp = UR_F_gs
        UR_F_gs = ONE
      ENDIF
      
      RRATE = .FALSE.      
      IF (ANY_SPECIES_EQ) RRATE = .TRUE. 
      DENSITY(:MMAX) = .TRUE. 
	 
! add by rong
      IF (Call_DQMOM) THEN
         PSIZE(:MMAX) = .TRUE. 
      ELSE
        PSIZE(:MMAX) = .FALSE.
      ENDIF
! addby rong

      WALL_TR = .TRUE. 
      IF (ENERGY_EQ) THEN 
         SP_HEAT(:MMAX) = .TRUE. 
         COND(:MMAX) = .TRUE. 
         HEAT_TR(:MMAX,:MMAX) = .TRUE. 
      ELSE 
         SP_HEAT(:MMAX) = .FALSE. 
         COND(:MMAX) = .FALSE. 
         HEAT_TR(:MMAX,:MMAX) = .FALSE. 
      ENDIF 
      VISC(:MMAX) = .TRUE. 
      DIFF(:MMAX) = .TRUE. 
      DRAGCOEF(:MMAX,:MMAX) = .TRUE. 
      IF (RO_G0 /= UNDEFINED) DENSITY(0) = .FALSE. 
      IF (MU_S0 /= UNDEFINED) VISC(1) = .FALSE. 

      CALL CALC_COEFF (DENSITY, PSIZE, SP_HEAT, VISC, COND, DIFF, RRATE, DRAGCOEF, &
      HEAT_TR, WALL_TR, IER) 

      
      IF(FLAG == 0) then
        !Now restore all underrelaxation factors for coefficient
	! calculations to their original value
        UR_F_gs = UR_F_gstmp
      ENDIF
!     
      RETURN  
      END SUBROUTINE CALC_COEFF_ALL

!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_COEFF(DENSITY, PSIZE, SP_HEAT, VISC, COND, DIFF,   C
!       RRATE, DRAG, HEAT_TR, WALL_TR, IER)                            C
!  Purpose: Calculate physical and transport properties, reaction ratesC
!           and exchange rates.                                        C
!                                                                      C
!  Author: M. Syamlal                                 Date: 23-APR-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CALC_COEFF(DENSITY, PSIZE, SP_HEAT, VISC, COND, DIFF, RRATE, &
         DRAG, HEAT_TR, WALL_TR, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE physprop
      USE rxns
      USE funits 
      USE compar
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Error index
      INTEGER          IER
!
!                      Loop indices
      INTEGER          L, M
!
!                      Flags to tell whether to calculate or not
      LOGICAL          DENSITY(0:DIMENSION_M), PSIZE(0:DIMENSION_M),&
                       SP_HEAT(0:DIMENSION_M)
!
!                      Flags to tell whether to calculate or not
      LOGICAL          VISC(0:DIMENSION_M), COND(0:DIMENSION_M),&
                       DIFF(0:DIMENSION_M)
!
!                      Flag for Reaction rates
      LOGICAL          RRATE
!
!                      Flag for exchange functions
      LOGICAL          DRAG(0:DIMENSION_M, 0:DIMENSION_M),&
                       HEAT_TR(0:DIMENSION_M, 0:DIMENSION_M),&
                       WALL_TR
!-----------------------------------------------
!
!     Calculate physical properties
!
      CALL PHYSICAL_PROP (DENSITY, PSIZE, SP_HEAT, IER) 

!
!     Calculate Transport properties
!
      CALL TRANSPORT_PROP (VISC, COND, DIFF, IER) 

!
!     Calculate reaction rates and interphase mass transfer
!
      CALL CALC_RRATE(RRATE)
!
!     Calculate interphase momentum, and energy transfers
!
      CALL EXCHANGE (DRAG, HEAT_TR, WALL_TR, IER) 
!
!     Reset all flags.  The flags need to be set every time this routine is
!     called.
!
      CALL TurnOffCOEFF(DENSITY, PSIZE, SP_HEAT, VISC, COND, DIFF, RRATE, &
         DRAG, HEAT_TR, WALL_TR, IER)

      RETURN  
      END SUBROUTINE CALC_COEFF 
      

      SUBROUTINE CALC_RRATE(RRATE) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE rxns
      USE funits 
      USE compar
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Error index
      INTEGER          IER
!
!                      Flag for Reaction rates
      LOGICAL          RRATE
!-----------------------------------------------
!
!     Calculate reaction rates and interphase mass transfer
!
      IF (RRATE) THEN 
         IF (NO_OF_RXNS > 0) THEN 
            CALL RRATES0 (IER)                   !rxns defined in mfix.dat and rrates0.f 
         ELSE
	    IER = 0 
            CALL RRATES (IER)                    !rxns defined in rrates.f 

	    IF(IER .EQ. 1) THEN		!error: rrates.f has not been modified
	      CALL START_LOG 
              IF(DMP_LOG)WRITE (UNIT_LOG, 1000)
              CALL END_LOG 
              call mfix_exit(myPE)  
	    ENDIF

         ENDIF
      ELSE
        !In case mass exchage w/o chemical rxn (e.g., evaporation) occur
        CALL RRATES (IER)           
      ENDIF 

      RETURN  
 1000 FORMAT(/1X,70('*')//' From: CALC_COEFF',/,&
         ' Species balance equations are being solved; but chemical',/,&
	 ' reactions are not specified in mfix.dat or in rrates.f.',/,&
	 ' Copy the file mfix/model/rrates.f into the run directory ',/,&
	 ' and remove the initial section that returns IER=1.'&
         ,/1X,70('*')/) 
      END SUBROUTINE CALC_RRATE 


 
      SUBROUTINE TurnOffCOEFF(DENSITY, PSIZE, SP_HEAT, VISC, COND, DIFF, RRATE, &
         DRAG, HEAT_TR, WALL_TR, IER) 
      USE param 
      USE param1 
      USE physprop
      IMPLICIT NONE
!
!                      Error index
      INTEGER          IER
!
!                      Loop indices
      INTEGER          L, M
!
!                      Flags to tell whether to calculate or not
      LOGICAL          DENSITY(0:DIMENSION_M), PSIZE(0:DIMENSION_M),&
                       SP_HEAT(0:DIMENSION_M)
!
!                      Flags to tell whether to calculate or not
      LOGICAL          VISC(0:DIMENSION_M), COND(0:DIMENSION_M),&
                       DIFF(0:DIMENSION_M)
!
!                      Flag for Reaction rates
      LOGICAL          RRATE
!
!                      Flag for exchange functions
      LOGICAL          DRAG(0:DIMENSION_M, 0:DIMENSION_M),&
                       HEAT_TR(0:DIMENSION_M, 0:DIMENSION_M),&
                       WALL_TR
		       
!     Reset all flags
!
      RRATE = .FALSE. 
      WALL_TR = .FALSE. 
      DENSITY(:MMAX) = .FALSE. 
      PSIZE(:MMAX) = .FALSE. 
      SP_HEAT(:MMAX) = .FALSE. 
      VISC(:MMAX) = .FALSE. 
      COND(:MMAX) = .FALSE. 
      DIFF(:MMAX) = .FALSE. 
      DRAG(:MMAX,:MMAX) = .FALSE. 
      HEAT_TR(:MMAX,:MMAX) = .FALSE. 
      RETURN  
      END SUBROUTINE TurnOffCOEFF 
