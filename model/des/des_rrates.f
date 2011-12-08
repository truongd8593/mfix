!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: RRATES(IER)                                            !
!  Purpose: Calculate reaction rates for various reactions in cell ijk !
!                                                                      !
!  Author:                                            Date:            !
!  Reviewer:                                          Date:            !
!                                                                      !
!  Revision Number:                                                    !
!  Purpose:                                                            !
!  Author:                                            Date: dd-mmm-yy  !
!  Reviewer:                                          Date: dd-mmm-yy  !
!                                                                      !
!  Literature/Document References:                                     !
!                                                                      !
!  Variables referenced: MMAX, IJK, T_g, T_s1, D_p, X_g, X_s, EP_g,    !
!            P_g, HOR_g, HOR_s                                         !
!                                                                      !
!                                                                      !
!  Variables modified: M, N, R_gp, R_sp, RoX_gc, RoX_sc, SUM_R_g,      !
!                      SUM_R_s                                         !
!                                                                      !
!  Local variables:                                                    !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DES_RRATES(NP, INTERP_IJK, INTERP_WEIGHTS, FOCUS, &
                    CALLER)

      USE param 
      USE param1 
      USE parallel 
      USE fldvar
      USE rxns
      USE energy
      USE geometry
      USE run
      USE indices
      USE physprop
      USE constant
      USE funits 
      USE toleranc
      USE compar        !//d
      USE sendrecv      !// 400
      USE usr

      USE des_thermo
      USE des_rxns
      USE discretelement
      USE interpolation

      IMPLICIT NONE
     

! Pass Variables
!-----------------------------------------------
! Particle index number
      INTEGER, INTENT(IN) :: NP
! IJK indicies of fluid cells involved in interpolation
      INTEGER, INTENT(IN) :: INTERP_IJK(2**DIMN)
! Weights associated with interpolation
      DOUBLE PRECISION, INTENT(IN) :: INTERP_WEIGHTS(2**DIMN)
! Indicates that debugging information for the particle
      LOGICAL, INTENT(IN) :: FOCUS
! Indicates which phase is calling for the calculations so that only
! the necessary data is calculated
      CHARACTER(len=*)  :: CALLER

! Local Routine Variables - REQUIRED DO NOT EDIT
!-----------------------------------------------
! Phase and species indices
      INTEGER L, LM, PARTICLE_PHASE, N
! Fluid cell index
      INTEGER IJK

! Interpolated field quantites
      DOUBLE PRECISION Xg(NMAX(0))
! Gas temperature at particle's location (K)
      DOUBLE PRECISION Tg
! Rate of gas phase species production (g/s)
      DOUBLE PRECISION DES_R_gp(NMAX(0))
! Rate of gas phase species consumption per species concentration (g/s)
      DOUBLE PRECISION DES_RoX_gc(NMAX(0))
! Rate of enthalpy transfer
      DOUBLE PRECISION RxH_xfr
! Heat of reaction
      DOUBLE PRECISION HOR
! Mass generated
      DOUBLE PRECISION MASS_GEN

      DOUBLE PRECISION R_tmp(0:MMAX, 0:MMAX)
      DOUBLE PRECISION X_tmp(0:MMAX, 0:MMAX, Dimension_n_all)


      DOUBLE PRECISION, EXTERNAL :: calc_h
      DOUBLE PRECISION, EXTERNAL :: DES_CALC_H

!......................................................................!
! Reaction specific variables:                                         !
!``````````````````````````````````````````````````````````````````````!
! Species indices of gas phase (for readability)
!     Gas Species
      INTEGER, PARAMETER :: N2  = 1   ! Nitrogen
      INTEGER, PARAMETER :: O2  = 2   ! Oxygen
      INTEGER, PARAMETER :: CO2 = 3   ! Carbon dioxide

!     Solids Species
!     - Phase 1
      INTEGER, PARAMETER :: Carbon = 1
      INTEGER, PARAMETER :: Ash = 2
!     - Phase 2
      INTEGER, PARAMETER :: Sand = 1

! Concentration of gas phase species at particle (g/cm^3)
      DOUBLE PRECISION Cg(NMAX(0))

! Particle temperature
      DOUBLE PRECISION Tp

! Reaction Rates:
!     a) Combustion:  C + O2 --> CO2
      DOUBLE PRECISION RXNA, RXNAF, RXNAB

! Binary diffusion coefficient of oxygen in air
      DOUBLE PRECISION D_OXY

! Surface area of particle
      DOUBLE PRECISION A_S

! Rate constats
      DOUBLE PRECISION K_DIFF, K_A, K_EFF

! Sherwood Number of a single particle
      DOUBLE PRECISION Sh

! Logical 
      LOGICAL,SAVE:: FIRST_PASS = .TRUE.
! Time to write output.
      DOUBLE PRECISION, SAVE :: OUTPUT_TIME


      LOGICAL NOISY

!  Function subroutines
!-----------------------------------------------
      LOGICAL COMPARE

      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'usrnlst.inc'

! Set the mass phase of the particle
      PARTICLE_PHASE = PIJK(NP,5)

! Return if the species equations for the particle's phase are not
! being solved.
      IF(.NOT.DES_SPECIES_EQ(PARTICLE_PHASE)) RETURN

! Return if the particle is entering.
      IF(.NOT.PEA(NP,1) .OR. PEA(NP,2)) RETURN

! Initialize variables
      R_tmp = UNDEFINED
! Set the fluid cell of the particle
      IJK = PIJK(NP,4)

      DES_R_gp(:) = ZERO
      DES_RoX_gc(:) = ZERO
      DES_R_sp(NP,:) = ZERO
      DES_R_sc(NP,:) = ZERO

! Set the particle temperature
      Tp = DES_T_s_NEW(NP)

! Debug flag
      NOISY = .FALSE.
!      IF(PARTICLE_PHASE == 1 ) NOISY = .TRUE.
      IF(NOISY)THEN
         WRITE(*,"(//3X,A,I4,/3X,25('*'))")'NP: ',NP
         WRITE(*,"(6X,A,A)")'CALLER: ', CALLER
         WRITE(*,"(6X,A,I4)")'IJK: ',IJK
         WRITE(*,"(6X,A,F14.8)")'Tp: ',Tp
      ENDIF


!  User input is required in sections 1,2,3.

!AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA!
!                                                                      !
! A. Set gas phase quantities.                                         !
!                                                                      !
!    The gas phase species information is set to the array Xg(:) and   !
!    the gas phase temperature is set to the variable Tg.              !
!                                                                      !
!    Additional quantites may also be set by using the function        !
!    TO_PARTICLE. This function will set the correct values depending  !
!    if DES_INTERP_ON is true or false.                                !
!                                                                      !
!``````````````````````````````````````````````````````````````````````!

! Set the gas phase species values. (REQUIRED - DO NOT DELETE)
      DO N=1,NMAX(0)
         Xg(N) = TO_PARTICLE(X_g(:,N))
      ENDDO

! Set the gas phase temperature. (REQUIRED - DO NOT DELETE)
      Tg = TO_PARTICLE(T_G)

! Set the gas phase species values. (REQUIRED - DO NOT DELETE)
      DO N=1,NMAX(0)
         Cg(N) = TO_PARTICLE_2(X_g(:,N), ROP_g(:))
      ENDDO


!1111111111111111111111111111111111111111111111111111111111111111111111!
!                                                                      !
! 1. Write the rates of reactions:                                     !
!                                                                      !
!    Instructions:                                                     !
!                                                                      !
!    - Write the reaction rates for each of the reactions as RXNxF and !
!      RXNxB (both quantities >= 0), where x identifies the reaction,  !
!      F stands for forward rate, and B stands for the backward rate.  !
!      The rates MUST BE in units g/s.                                 !
!                                                                      !
!    - For the sake of clarity, give the reaction scheme and the units !
!      in a comment statement above the rate expression.               !
!                                                                      !
!``````````````````````````````````````````````````````````````````````!

! A) Combustion: 2C + O2 --> 2CO;  g/s
      IF(PARTICLE_PHASE == 1 .AND. &
        (DES_X_s(NP,Carbon) > ZERO) .AND. (Xg(O2) > ZERO)) THEN
! Diffusion coefficient of O2 in air. (300K 1ATM) Mills, 1998 pp. 947
         D_OXY = 0.188*10.0d0**(4) ! cm^2/sec
! Surface area of unreacted core
         A_S = 4.0d0*PI*CORE_RAD(NP)**2 ! cm^2
! Calculate the Sherwood Number of the particle.
         CALL N_Sh(Sh)
! diffustion rate coefficient
         k_DIFF = (Sh * D_OXY) / (2.0d0 * DES_RADIUS(NP))
! rate of chemical reaction A :: Ea/R = 17950 K
         k_A = 59500.0d0*Tp*EXP(-17950/Tp)
! overall rate constant
         k_EFF = 1.0d0 / (1.0d0/k_A + 1.0d0/k_DIFF)
! Rate of forward reaction (first order)
         RXNAF = A_S * k_EFF * Cg(O2) ! g/sec
! No backward reaction
         RXNAB = ZERO

         IF(NOISY)THEN
            WRITE(*,"(6X,A,F14.8)")'DES_RADIUS: ',DES_RADIUS(NP)
            WRITE(*,"(6X,A,F14.8)")'N_Sh: ',Sh
            WRITE(*,"(6X,A,F14.8)")'A_S: ',A_S
            WRITE(*,"(6X,A,F14.8)")'k_DIFF: ',k_DIFF
            WRITE(*,"(6X,A,F14.8)")'k_A: ',k_A
            WRITE(*,"(6X,A,F14.8)")'k_EFF ',k_EFF
         ENDIF

      ELSE
! Rate of forward reaction (first order)
         RXNAF = ZERO
! No backward reaction
         RXNAB = ZERO
      ENDIF

      IF(NOISY)THEN
         WRITE(*,"(6X,A,F14.8)")'RXNAF: ',RXNAF
         WRITE(*,"(6X,A,F14.8)")'RXNAB: ',RXNAB
      ENDIF


!2222222222222222222222222222222222222222222222222222222222222222222222!
!                                                                      !
! 2. Write the formation and consumption rates of various species:     !
!                                                                      !
!    Instructions:                                                     !
!                                                                      !
!    - Obtain the rates of formation and consumption of the various    !
!      species in g/s from the rate expressions RXNxF and RXNxB in the !
!      previous section.                                               !
!                                                                      !
!    - Pay attention to the units of RXNxF and RXNxB.                  !
!      * The formation rates for gas species N are added to get        !
!        R_gp(IJK,N).                                                  !
!      * The consumption rates are added and then divided by X_g(IJK,N)!
!        to get RoX_gc(IJK, n).                                        !
!      * This is done automatically, therefore the user does not need  !
!        to specify these quantities.                                  !
!                                                                      !
!    - The rates of formation and consumption of gas and solids phase  !
!      species are split between separate calls to this routine. When  !
!      the gas phase invokes des_rrates.f (CALLER = 'GAS'), only the   !
!      rates of formation and consumption of the gas phase species are !
!      calculated. Conversely, when the solids phase invokes this      !
!      routine only the rates of formation and consumption for the     !
!      solids phase species are calculated.                            !
!                                                                      !
!    - If X_g(IJK, n) is zero and  species N is likely to be consumed  !
!      in a reaction then it is recommended that RoX_gc(IJK,N) be      !
!      initialized to the derivative of the consumption rate with      !
!      respect to X_g at X_g=0. If the slope is not known analytically,!
!      then a small value (such as 1.0e-9) may be used.                !
!                                                                      !
!    - If you want to clear continuum variables, (e.g. R_gp(IJK,:) =   !
!      ZERO) you may consider doing this in the model/rrates.f file.   !
!      Thus, if there is a heterogeneous gas phase reaction, these     !
!      quantities are not accidentally overwritten.                    !
!                                                                      !
!``````````````````````````````````````````````````````````````````````!


! Solid species
!```````````````````````````````````````````````````````````````````````
      IF(CALLER == 'SOLIDS') THEN

! CARBON
!--------->
         IF(DES_X_s(NP, CARBON) .GT. ZERO) DES_R_sc(NP, CARBON) = RXNAF

         IF(NOISY) WRITE(*,"(6X,A,F14.8)")'DES_R_s: ',DES_R_sc(NP,CARBON)
           

! Gas species
!```````````````````````````````````````````````````````````````````````
      ELSEIF(CALLER == 'GAS') THEN

! O2
!--------->
         DES_R_gp(O2) = ZERO
         IF(Xg(O2) .GT. ZERO) THEN
            DES_RoX_gc(O2) = MW_g(O2) * &
               (RXNAF / DES_MW_s(PARTICLE_PHASE,CARBON))
         ELSE
            DES_RoX_gc(O2) = 1.0D-9
         ENDIF
! CO2
!---------->
         DES_R_gp(CO2) = MW_g(CO2) * &
            (RXNAF / DES_MW_s(PARTICLE_PHASE,CARBON))
         DES_RoX_gc(CO2) = ZERO

         IF(NOISY) THEN
            WRITE(*,"(6X,A,F14.8)")'DES_R_gp(O2): ',DES_R_gp(O2)
            WRITE(*,"(6X,A,F14.8)")'DES_RoX_gc(O2): ',DES_RoX_gc(O2)
            WRITE(*,"(6X,A,F14.8)")'DES_R_gp(CO2): ',DES_R_gp(CO2)
            WRITE(*,"(6X,A,F14.8)")'DES_RoX_gc(CO2): ',DES_RoX_gc(CO2)
         ENDIF

      ENDIF

!3333333333333333333333333333333333333333333333333333333333333333333333!
!                                                                      !
! 3. Determine the rate of mass transfer in g/s transferred between    !
!    the gas phase and the particle.                                   !
!                                                                      !
!   ~-~-~-~-~-~-~-~-~-~-~-~-    WARNING    -~-~-~-~-~-~-~-~-~-~-~-~-   !
!   !                                                              !   !
!   ! The current version of MFIX-DEM des not support mass         !   !
!   ! precipitating from a two fluid reaction or a heterogeneous   !   !
!   ! solids phase reaction!                                       !   !
!   !                                                              !   !
!   ~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-   !
!                                                                      !
!    Instructions:                                                     !
!                                                                      !
!    R_tmp(To phase #, From phase #)                                   !
!                                                                      !
!    Ex: R_tmp(0,PARTICLE_PHASE) specifies the rate of mass generation !
!        of the gas phase from the solids phase containing the         !
!        particle passed in the call to this routine.                  !
!                                                                      !
!      * Ex: Gas is generated from solids phase 1:                     !
!        IF(PARTICLE_PHASE == 1) THEN                                  !
!           R_tmp(0,1) > 0                                             !
!        ENDIF                                                         !
!                                                                      !
!    - The R-phase matrix (R_tmp) is skew-symmetric and the diagonal   !
!      elements are not needed. Furthermore, only one of the two       !
!      skew-symmetric entries must be specified.                       !
!                      R_tmp(0,m) = -R_tmp(m,0)                        !
!                                                                      !
!    - X_tmp(L,M,N) is the mass fraction of species N in the material  !
!      being transferred between phase-L and phase-M. N is the index   !
!      of the species being transferred in the destination phase.      !
!                                                                      !
!      * Ex: In the reaction C + (1/2)O2 -> CO, the destination phase  !
!            is the gas phase. Therefore, N is the index of CO in the  !
!            gas phase.                                                !
!                                                                      !
!``````````````````````````````````````````````````````````````````````!

      IF(PARTICLE_PHASE == 1) THEN
! Carbon from the particle combines with oxygen in the gas phase to
! produce CO2. Thus, the rate of mass transfer from the particle to the
! gas phase is equal to the rate of carbon consumption.
         R_tmp(0, PARTICLE_PHASE) = RXNAF
! All of the products of the reaction are transfered to the carbon
! dioxide species in the gas phase.
         X_tmp(0, PARTICLE_PHASE,:) = ZERO
      	  X_tmp(0, PARTICLE_PHASE, CO2) = 1.0d0
      ENDIF

!======================================================================!
! <---------   No user input is required below this line   --------->  !
!======================================================================!

      IF(CALLER == 'GAS')THEN
! Store the transfer quanties in the continuum variables
!-----------------------------------------------------------------------
         DO N=1,NMAX(0)
            CALL TO_FIELD(R_gp(:,N), DES_R_gp(N))
            CALL TO_FIELD(RoX_gc(:,N), DES_RoX_gc(N))
         ENDDO

! Calculate the rate of mass generate in g/(s.cm^3) for the gas phase.
!-----------------------------------------------------------------------
         IF (SPECIES_EQ(0)) THEN 
! If the species equations are being solved for the gas phase, then
! the total rate of mass generate is the sum of the production and
! consumption terms for each species.
            DO N=1,NMAX(0)
               MASS_GEN = DES_R_gp(N) - DES_RoX_gc(N)*Xg(N)
            ENDDO
	        ELSE
! If the species equations are not being solved, use the values provided
! in the R_tmp array.
   	        IF(R_tmp(0,PARTICLE_PHASE) .NE. UNDEFINED)THEN
           		  MASS_GEN = R_tmp(0,PARTICLE_PHASE)
          		ELSEIF(R_tmp(PARTICLE_PHASE,0) .NE. UNDEFINED)THEN
        	   	  MASS_GEN = R_tmp(PARTICLE_PHASE,0)
           	ENDIF
         ENDIF
         CALL TO_FIELD(SUM_R_G, MASS_GEN)
      ENDIF


! Calculate the rate enthalpy of transfer from material exchange
!-----------------------------------------------------------------------
      IF(ENERGY_EQ) THEN
         RxH_xfr = ZERO
! Mass transfer is occuring.
         IF(R_tmp(0,PARTICLE_PHASE) .NE. UNDEFINED) THEN
! solids phase species >> gas phase species
            IF(R_tmp(0, PARTICLE_PHASE) > ZERO) THEN
               DO N = 1, NMAX(0)
                  RxH_xfr =  RxH_xfr + R_tmp(0,PARTICLE_PHASE) * &
                     X_tmp(0,PARTICLE_PHASE,N) * CALC_H(IJK,0,N)
               ENDDO
! gas phase species >> solids phase species
            ELSE
               DO N = 1, DES_NMAX(PARTICLE_PHASE)
                  RxH_xfr = RxH_xfr + R_tmp(0,PARTICLE_PHASE) * &
                     X_tmp(0,PARTICLE_PHASE,N) * DES_CALC_H(NP,N)
               ENDDO
            ENDIF
       		ELSEIF(R_tmp(PARTICLE_PHASE,0) .NE. UNDEFINED) THEN
! gas phase species >> solids phase species
		          IF(R_tmp(PARTICLE_PHASE,0) > ZERO) THEN
               DO N = 1, DES_NMAX(PARTICLE_PHASE)
                  RxH_xfr = RxH_xfr - R_tmp(PARTICLE_PHASE,0) * &
                     X_tmp(PARTICLE_PHASE,0,N) * DES_CALC_H(NP,N)
               ENDDO
! solids phase species >> gas phase species
            ELSE
               DO N = 1, NMAX(0)
                  RxH_xfr = RxH_xfr - R_tmp(PARTICLE_PHASE,0) * &
                     X_tmp(PARTICLE_PHASE,0,N) * CALC_H(IJK,0,N)
               ENDDO
      		    ENDIF
         ENDIF

! Calculate heats of reactions
!-----------------------------------------------------------------------
         IF(CALLER == 'GAS') THEN
! Gas phase
            HOR = ZERO
            DO N = 1, NMAX(0)
               HOR = HOR + (DES_R_gp(N) - DES_RoX_gc(N) * Xg(N)) * &
                  CALC_H(IJK, 0, N)
            END DO 
            HOR = HOR - RxH_xfr
            IF (UNITS == 'SI') HOR = 4183.925D0*HOR  ! in J/kg K
            CALL TO_FIELD(HOR_G, HOR)

            ELSEIF(CALLER == 'SOLIDS') THEN
! Solids phase
            HOR = ZERO
            DO N = 1, DES_NMAX(PARTICLE_PHASE)
               HOR = HOR + (DES_R_sp(NP,N) - DES_R_sc(NP,N)) * &
                   DES_X_s(NP,N) * DES_CALC_H(NP,N)
            ENDDO 
            HOR = HOR + RxH_xfr
            IF (UNITS == 'SI') HOR = 4183.925D0*HOR ! in J/kg K
! Calculate the rate of heat generation for the particle resulting from
! the reactions. (cal/sec)
            Qint(NP) = HOR
         ENDIF
      ENDIF


! Store R_tmp values in an array.  Only store the upper triangle without
! the diagonal of R_tmp array.
!-----------------------------------------------------------------------

      LM = L + 1 + (PARTICLE_PHASE - 1)*PARTICLE_PHASE/2
      IF (R_TMP(0,PARTICLE_PHASE) /= UNDEFINED) THEN 
         R_PHASE(IJK,LM) = R_TMP(0,PARTICLE_PHASE) 
      ELSE IF (R_TMP(PARTICLE_PHASE,0) /= UNDEFINED) THEN 
         R_PHASE(IJK,LM) = -R_TMP(PARTICLE_PHASE,0) 
      ELSE 
         CALL START_LOG 
         WRITE (UNIT_LOG, 1001) 0, PARTICLE_PHASE 
         CALL END_LOG 
         call mfix_exit(myPE)  
      ENDIF


      IF(CALLER == 'GAS') THEN

         IF(FIRST_PASS) THEN
            OUTPUT_TIME = 0.0d0
            CALL WRITE_DES_RX
            FIRST_PASS = .FALSE.
         ENDIF

         IF(OUTPUT_TIME .LE. S_TIME) THEN
            CALL WRITE_DES_RX
            OUTPUT_TIME = OUTPUT_TIME + 0.10d0
         ENDIF
      ENDIF

      RETURN


   
 1000 FORMAT(/1X,70('*')/' From: CALC_COEFF',/,' Message:',            &
         ' Discrete species equations are being solve; but no',        &
         ' reactions',/' are specified in the des_rrates.f file. For', &
         ' discrete gas-solid',/' reactions:',/4X,'1) Copy',           &
         ' mfix/model/des/des_rrates.f into the run directory',/4X,    &
         '2) Remove the initial secion that returns IER=1',/4X,        &
         '3) Follow the instructions to specify the desired reaction', &
         /1X,70('*'))

 1001 FORMAT(/1X,70('*')/' From: DES_RRATES',/&
         ' Message: Mass transfer between phases ',I2,' and ',I2,&
         ' (R_tmp) not specified',/1X,70('*')/) 

      CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Module name: WRITE_DES_RX                                           !
!                                                                      !
! 	THIS ROUTINE IS NOT DESIGNED TO BE CHECKED INTO CVS. IT IS STRICTLY !
!  HERE FOR PRINTING MESSAGES USED FOR DEBUGGING AND V&V WORK.         !
!                                                                      !
!  Author: J.Musser                                                    !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE WRITE_DES_RX

      USE run

      IMPLICIT NONE

!-----------------------------------------------
! Local variables
!-----------------------------------------------
! file name
      CHARACTER*64 :: FNAME
! logical used for testing is the data file already exists
      LOGICAL :: F_EXISTS
! file unit for heat transfer data
      INTEGER, PARAMETER :: RX_UNIT = 2035
      INTEGER, PARAMETER :: RX2_UNIT = 2036
      INTEGER, PARAMETER :: RX3_UNIT = 2037


      FNAME = TRIM(RUN_NAME)//'_DES_RX.dat'
      INQUIRE(FILE=FNAME,EXIST=F_EXISTS)
      IF (.NOT.F_EXISTS) THEN
         OPEN(UNIT=RX_UNIT,FILE=FNAME,STATUS='NEW')
      ELSE
         OPEN(UNIT=RX_UNIT,FILE=FNAME,&
            POSITION="APPEND",STATUS='OLD')
      ENDIF
      WRITE(RX_UNIT,*)S_TIME, SH, A_S, K_EFF
      CLOSE(RX_UNIT)


      FNAME = TRIM(RUN_NAME)//'_DES_RX2.dat'
      INQUIRE(FILE=FNAME,EXIST=F_EXISTS)
      IF (.NOT.F_EXISTS) THEN
         OPEN(UNIT=RX2_UNIT,FILE=FNAME,STATUS='NEW')
      ELSE
         OPEN(UNIT=RX2_UNIT,FILE=FNAME,&
            POSITION="APPEND",STATUS='OLD')
      ENDIF
      WRITE(RX2_UNIT,*)DES_RoX_gc(O2), DES_R_gp(CO2)
      CLOSE(RX2_UNIT)


      FNAME = TRIM(RUN_NAME)//'_DES_RX3.dat'
      INQUIRE(FILE=FNAME,EXIST=F_EXISTS)
      IF (.NOT.F_EXISTS) THEN
         OPEN(UNIT=RX3_UNIT,FILE=FNAME,STATUS='NEW')
      ELSE
         OPEN(UNIT=RX3_UNIT,FILE=FNAME,&
            POSITION="APPEND",STATUS='OLD')
      ENDIF
      WRITE(RX3_UNIT,*)S_TIME, HOR, RXNAF
      CLOSE(RX3_UNIT)


      END SUBROUTINE WRITE_DES_RX


!----------------------------------------------------------------------!
! Function: N_SH                                                       !
!                                                                      !
! Purpose: Calculate the Sherwood Number of a single particle.         !
!''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''!
      SUBROUTINE N_Sh(Sh)

      DOUBLE PRECISION, INTENT(OUT) :: Sh

      INTEGER I, IMJK, IJMK, IJKM
      DOUBLE PRECISION Sc1o3
      DOUBLE PRECISION UGC, VGC, WGC
      DOUBLE PRECISION VREL, N_Re, N_Sc

! Slip velocity between the particle and fluid
      DOUBLE PRECISION SLIP_VEL(3)

! Functions
!-----------------------------------------------      
      DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT 

! Initialize fluid cell variables
      I =  PIJK(NP,1)
      IJK   = PIJK(NP,4)
      IMJK  = IM_OF(IJK)
      IJMK  = JM_OF(IJK)
      IJKM  = KM_OF(IJK)

! Gas velocity in fluid cell
      UGC = AVG_X_E(U_g(IMJK), U_g(IJK), I)
      VGC = AVG_Y_N(V_g(IJMK), V_g(IJK))
      WGC = AVG_Z_T(W_g(IJKM), W_g(IJK))

! Calculate the slip velocity between the particle and fluid
      SLIP_VEL(1) = UGC - DES_VEL_NEW(NP,1)
      SLIP_VEL(2) = VGC - DES_VEL_NEW(NP,2)
      IF(DIMN == 3) SLIP_VEL(3) = WGC - DES_VEL_NEW(NP,3)
! Calculate the magnitude of the slip velocity
      VREL = SQRT(DES_DOTPRDCT(SLIP_VEL,SLIP_VEL))

! Schmidt Number
      N_Sc = MU_g(IJK)/(RO_g(IJK)*D_OXY)
! Reynods Number
      N_Re = 2.0d0 * DES_RADIUS(NP) * VREL * RO_g(IJK) / MU_g(IJK)
! Sherwood Number: Ranzy and Marshall, 1952
      Sh = 2.0d0 + 0.60d0*(N_Re**(1.0d0/2.0d0))*(N_Sc**(1.0d0/3.0d0))

      END SUBROUTINE N_Sh


!----------------------------------------------------------------------!
! Function: TO_PARTICLE_2                                              !
!                                                                      !
! Purpose: Interpolate scalar field variables to a particle's position.!
!''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''!
      DOUBLE PRECISION FUNCTION TO_PARTICLE_2(FIELD1, FIELD2) RESULT(VALUE)
! Field to be interpolated
      DOUBLE PRECISION, INTENT(IN) :: FIELD1(DIMENSION_3)
      DOUBLE PRECISION, INTENT(IN) :: FIELD2(DIMENSION_3)

! Loop counter
      INTEGER LC
!     Initialize the result
      VALUE = ZERO
      IF(DES_INTERP_ON)THEN
         DO LC=1, 2**DIMN
            VALUE = VALUE + INTERP_WEIGHTS(LC)* &
               (FIELD1(INTERP_IJK(LC))*FIELD2(INTERP_IJK(LC)))
         ENDDO
      ELSE
! As a safe-guard, if DES_INTERP_ON = .FALSE., return the value of the
! field variable associated with the fluid cell containing the
! particle's center.
         VALUE = FIELD1(PIJK(NP,4))*FIELD2(PIJK(NP,4))
      ENDIF
      RETURN
      END FUNCTION TO_PARTICLE_2


!----------------------------------------------------------------------!
! Function: TO_PARTICLE - DO NOT EDIT                                  !
!                                                                      !
! Purpose: Interpolate scalar field variables to a particle's position.!
!''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''!
      DOUBLE PRECISION FUNCTION TO_PARTICLE(FIELD) RESULT(VALUE)
! Field to be interpolated
      DOUBLE PRECISION, INTENT(IN) :: FIELD(DIMENSION_3)
! Loop counter
      INTEGER LC
!     Initialize the result
      VALUE = ZERO
      IF(DES_INTERP_ON)THEN
         DO LC=1, 2**DIMN
            VALUE = VALUE + INTERP_WEIGHTS(LC)*FIELD(INTERP_IJK(LC))
         ENDDO
      ELSE
! As a safe-guard, if DES_INTERP_ON = .FALSE., return the value of the
! field variable associated with the fluid cell containing the
! particle's center.
         VALUE = FIELD(PIJK(NP,4))
      ENDIF
      RETURN
      END FUNCTION TO_PARTICLE

!----------------------------------------------------------------------!
! Function: TO_FIELD - DO NOT EDIT                                     !
!                                                                      !
! Purpose: Distribute particle effects to Eulerian grid. The provided  !
! scalar values should have the units g/sec. These must be mapped to   !
! the grid by dividing by the cell volume to get g/(cm^3.sec).         !
!''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''!
      SUBROUTINE TO_FIELD(FIELD, SCALAR)
! Field to be distributed to from particle
      DOUBLE PRECISION, INTENT(INOUT) :: FIELD(DIMENSION_3)
! Value to be distributed
      DOUBLE PRECISION, INTENT(IN) :: SCALAR
! Loop counter
      INTEGER LC

      IF(DES_INTERP_ON) THEN
         DO LC = 1, 2**DIMN
! Set variable values
! Portion the source term to cell IJK
            FIELD(INTERP_IJK(LC)) = FIELD(INTERP_IJK(LC)) + &
               INTERP_WEIGHTS(LC) * SCALAR / VOL(INTERP_IJK(LC))
         ENDDO ! LC-loop
      ELSE
         FIELD(PIJK(NP,4)) = FIELD(PIJK(NP,4)) + SCALAR/VOL(PIJK(NP,4))
      ENDIF
      RETURN
      END SUBROUTINE TO_FIELD


      END SUBROUTINE DES_RRATES
