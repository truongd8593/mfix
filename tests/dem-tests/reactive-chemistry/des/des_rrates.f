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
      INTEGER, PARAMETER :: A = 1
      INTEGER, PARAMETER :: CC = 2

!     Solids Species
      INTEGER, PARAMETER :: B = 1
      INTEGER, PARAMETER :: D = 2
      INTEGER, PARAMETER :: I = 3

!     gas constant
      DOUBLE PRECISION, PARAMETER :: RGAS = 1.987207D0  !cal/mol.K


! Particle temperature
      DOUBLE PRECISION Tp

! Reaction Rates:
!     a) Combustion:  A(g) + 2B(s) --> CC(g) + D(s)
      DOUBLE PRECISION RXNxF, RXNxB

! Surface area of particle
      DOUBLE PRECISION A_S




      DOUBLE PRECISION H_CC, H_A, H_B
      DOUBLE PRECISION r_C, Q_int
      DOUBLE PRECISION H_s, H_fr

      DOUBLE PRECISION rrB, rrD




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

! A)    A(g) + 2B(s) -> CC(g) + D(s);  g/s
         IF(DES_X_s(NP, B) .GT. ZERO) THEN
            RXNxF = 0.0011781 ! (grams of B)/sec
         ELSE
            RXNxF = ZERO
         ENDIF
           RXNxB = ZERO


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

! B
!--------->
         DES_R_sp(NP,B) = ZERO
         IF(DES_X_s(NP, B) .GT. ZERO) THEN
            DES_R_sc(NP,B) = RXNxF
         ELSE
            DES_R_sc(NP,B) = ZERO
         ENDIF

! D
!--------->
         IF(DES_X_s(NP, B) .GT. ZERO) THEN
            DES_R_sp(NP,D) = (RXNxF/2.0d0)*(DES_MW_s(PARTICLE_PHASE,D)/&
              DES_MW_s(PARTICLE_PHASE,B))
         ELSE
            DES_R_sp(NP,D) = ZERO
         ENDIF
         DES_R_sc(NP,D) = ZERO

! I
!--------->
         DES_R_sp(NP,I) = ZERO
         DES_R_sc(NP,I) = ZERO


         IF(NOISY) WRITE(*,"(6X,A,F14.8)")'DES_R_sp: ',DES_R_sp(NP,B)
         IF(NOISY) WRITE(*,"(6X,A,F14.8)")'DES_R_sc: ',DES_R_sc(NP,B)
            

! Gas species
!```````````````````````````````````````````````````````````````````````
      ELSEIF(CALLER == 'GAS') THEN
! A
!--------->
         DES_R_gp(A) = ZERO
         IF(Xg(A) .GT. ZERO) THEN
            DES_RoX_gc(A) = (RXNxF/2.0d0)*(MW_g(A) / &
               DES_MW_s(PARTICLE_PHASE,B))/ Xg(A)
         ELSE
            DES_RoX_gc(A) = 1.0D-9
         ENDIF
! CC
!---------->
         DES_R_gp(CC) = (RXNXF/2.0d0)*(MW_g(CC) / &
            DES_MW_s(PARTICLE_PHASE,B))
         DES_RoX_gc(CC) = ZERO


         IF(NOISY) THEN
            WRITE(*,"(6X,A,F14.8)")'DES_R_gp(A): ',DES_R_gp(A)
            WRITE(*,"(6X,A,F14.8)")'DES_RoX_gc(A): ',DES_RoX_gc(A)
            WRITE(*,"(6X,A,F14.8)")'DES_R_gp(CC): ',DES_R_gp(CC)
            WRITE(*,"(6X,A,F14.8)")'DES_RoX_gc(CC): ',DES_RoX_gc(CC)
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

         R_tmp(0, PARTICLE_PHASE) = RXNxF * (1.0d0 - ((1.0d0/2.0d0) * &
            (DES_MW_s(PARTICLE_PHASE,D)/DES_MW_s(PARTICLE_PHASE,B))))
         X_tmp(0, PARTICLE_PHASE,:) = ZERO
      	  X_tmp(0, PARTICLE_PHASE, CC) = 1.0d0
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
!----------------------------------------------------------------------->>> This is correct
!            DO N = 1, NMAX(0)
!               RxH_xfr =  RxH_xfr + R_tmp(0,PARTICLE_PHASE) * &
!                  X_tmp(0,PARTICLE_PHASE,N) * CALC_H(IJK,0,N)
!            ENDDO
!-----------------------------------------------------------------------<<< Do not delet


!--------------------------------------------------------------------------->>> 
!  This is not correct. This should be removed. This is only here so
!  that the reactive chemistry code can be tested with a constat
!  enthalpy in the gas phase.
            H_CC = -110.40038888889
            RxH_xfr =  R_tmp(0,PARTICLE_PHASE) * H_CC
!---------------------------------------------------------------------------<<<

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
               HOR = HOR + (DES_R_gp(N) - DES_RoX_gc(N) * Xg(N)) *  CALC_H(IJK, 0, N)
            END DO 
            HOR = HOR - RxH_xfr
            IF (UNITS == 'SI') HOR = 4183.925D0*HOR  ! in J/kg K
            CALL TO_FIELD(HOR_G, HOR)
            ELSEIF(CALLER == 'SOLIDS') THEN
! Solids phase
            HOR = ZERO
            DO N = 1, DES_NMAX(PARTICLE_PHASE)
               HOR = HOR + (DES_R_sp(NP,N) - DES_R_sc(NP,N)) * &
                   DES_X_s(NP, N) * DES_CALC_H(NP,N)
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
