!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: RRATES0(IER)                                           C
!  Purpose: Calculate reaction rates for various reactions in cell ijk C
!           using information from the data file                       C
!                                                                      C
!  Author: M. Syamlal                                 Date: 3-10-98    C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose:Replaced routines with new proceedures for automated        C
!          reaction rate calculations.                                 C
!  Author: J. Musser                                  Date: 10-Oct-12  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: MMAX, IJK, T_g, T_s1, D_p, X_g, X_s, EP_g,    C
!            P_g, HOR_g, HOR_s                                         C
!                                                                      C
!                                                                      C
!  Variables modified: M, N, R_gp, R_sp, RoX_gc, RoX_sc, SUM_R_g,      C
!                      SUM_R_s                                         C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE RRATES0(IER) 

!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
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
      USE compar 
      USE sendrecv 
      Use parse

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

      INTEGER IJK  ! fluid cell index
      INTEGER H    ! Reaction loop counter
      INTEGER L, M ! Global Phase index loop counters
      INTEGER N    ! Global species index
      INTEGER lN   ! Local reaction speices index/loop counter
      INTEGER LM   ! 

      INTEGER mXfr ! Global phase index for mass transfer

! User-defined reaction rates returned from USR_RATES
      DOUBLE PRECISION RATES(NO_OF_RXNS)

      DOUBLE PRECISION lRate

      DOUBLE PRECISION RxH(0:MMAX, 0:MMAX)
      DOUBLE PRECISION lHoRg, LHoRs(1:MMAX)

! External functions for calculating enthalpy (cal/gram)
      DOUBLE PRECISION, EXTERNAL ::CALC_H
      DOUBLE PRECISION, EXTERNAL ::CALC_H0

! External Function for comparing two numbers.
      LOGICAL, EXTERNAL :: COMPARE


!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'ep_s2.inc'


! Initialize global storage arrays to zero
!---------------------------------------------------------------------//
      SUM_R_G(:) = ZERO
      HOR_G(:) = ZERO
      R_GP(:,:) = ZERO
      ROX_GC(:,:) = ZERO
      SUM_R_S(:,:) = ZERO
      HOR_S(:,:) = ZERO
      R_SP(:,:,:) = ZERO
      ROX_SC(:,:,:) = ZERO
      R_PHASE(:,:) = ZERO

! Loop over each fluid cell.
      DO IJK = ijkstart3, ijkend3 
      IF (FLUID_AT(IJK)) THEN 

      RATES(:) = ZERO

! Calculate user defined reaction rates.
      CALL USR_RATES(IJK, RATES)

! Loop over reactions.
      RXN_LP: DO H = 1, NO_OF_RXNS

! Skip empty reactions
         IF(Reaction(H)%nSpecies == 0) CYCLE RXN_LP
         IF(COMPARE(RATES(H),ZERO)) CYCLE RXN_LP

! Initialize local loop arrays
         lHoRg = ZERO
         lHoRs(:) = ZERO
         RxH(:,:) = ZERO

! Calculate the rate of formation/consumption for each species.
!---------------------------------------------------------------------//
         DO lN = 1, Reaction(H)%nSpecies
! Global phase index.
            M = Reaction(H)%Species(lN)%pMap
! Global species index.
            N = Reaction(H)%Species(lN)%sMap
! Index for interphase mass transfer. For a gas/solid reaction, the 
! index is stored with the gas phase. For solid/solid mass transfer
! the index is stored with the source phase.
            mXfr = Reaction(H)%Species(lN)%mXfr
            lRate = RATES(H) * Reaction(H)%Species(lN)%MWxStoich
! Gas Phase:
            IF(M == 0) THEN
! Consumption of gas phase species.
               IF(lRate < ZERO) THEN
                  IF(X_g(IJK,N) > SMALL_NUMBER) THEN
                     RoX_gc(IJK,N) = RoX_gc(IJK,N) - lRate/X_g(IJK,N)
                  ELSE
                     RoX_gc(IJK,N) = 1.0d-9
                  ENDIF
! Enthalpy transfer associated with mass transfer. (gas/solid)
                  IF(M /= mXfr) RxH(M,mXfr) =  RxH(M,mXfr) + &
                     lRate * CALC_H0(T_G(IJK),0,N)
               ELSE
! Formation of gas phase species.
                  R_gp(IJK,N) = R_gp(IJK,N) + lRate
! Enthalpy transfer associated with mass transfer. (gas/solid)
                  IF(M /= mXfr) RxH(M,mXfr) =  RxH(M,mXfr) + &
                     lRate * CALC_H0(T_s(IJK,mXfr),0,N)
               ENDIF
! Solids Phase M:
            ELSE
! Consumption of solids phase species.
               IF(lRate < ZERO) THEN
                  IF(X_s(IJK,M,N) > SMALL_NUMBER) THEN
                     RoX_sc(IJK,M,N) = &
                        RoX_sc(IJK,M,N) - lRate/X_s(IJK,M,N)
                  ELSE
                     RoX_sc(IJK,M,N) = 1.0d-9
                  ENDIF
! Enthalpy transfer associated with mass transfer. (solid/solid) This
! is only calculated from the source (reactant) material.
                  IF(M /= mXfr) THEN
                     IF(M < mXfr) THEN
                        RxH(M,mXfr) =  RxH(M,mXfr) + lRate * &
                          Reaction(H)%Species(lN)%xXfr * CALC_H(IJK,M,N)
                     ELSE
                        RxH(mXfr,M) =  RxH(mXfr,M) - lRate * &
                          Reaction(H)%Species(lN)%xXfr * CALC_H(IJK,M,N)
                     ENDIF
                  ENDIF
               ELSE
! Formation of solids phase species.
                  R_sp(IJK,M,N) = R_sp(IJK,M,N) + lRate
               ENDIF
            ENDIF
         ENDDO ! Loop of species


! Calculate and store the heat of reaction.
!---------------------------------------------------------------------//
         IF(ENERGY_EQ) THEN
! Automated heat of reaction calculations
            IF(Reaction(H)%Calc_DH) THEN
! Loop over reaction species.
               DO lN = 1, Reaction(H)%nSpecies
! Global phase index.
                  M = Reaction(H)%Species(lN)%pMap
! Global species index.
                  N = Reaction(H)%Species(lN)%sMap
! Rate of formation/consumption for speices N
                  lRate = RATES(H) * Reaction(H)%Species(lN)%MWxStoich
! Gas phase enthalpy chnage from energy equation derivation.
                  IF(M == 0) THEN
                     lHORg = lHORg + CALC_H(IJK,0,N) * lRate
! Solid phase enthalpy change from energy equation derivation.
                  ELSE
                     lHORs(M) = lHORs(M) + CALC_H(IJK,M,N) * lRate
                  ENDIF
               ENDDO

! Complete the skew-symettric for enthalpy transfer with mass transfer
               DO M=1, MMAX
                   DO L=0, M-1
                    RxH(M,L) = - RxH(L,M)
                  ENDDO
               ENDDO
! Apply enthalpy transfer associated with mass transfer to get the
! complete heat of reaction of heat phse for Reaction H.
               DO L=0, MMAX
                  DO M = 0, MMAX
                     IF(L == M) CYCLE
                     IF(L == 0) THEN
                        lHORg = lHORg - RxH(L,M)
                     ELSE
                        lHORs(L) = lHORs(L) - RxH(L,M)
                     ENDIF
                  ENDDO
               ENDDO

! Convert the heat of reaction to the appropriate units (if SI), and 
! store in the global array.
               IF(UNITS == 'SI') THEN
                  HOR_g(IJK) = HOR_g(IJK) + 4.183925d3*lHORg
                  DO M=1,MMAX
                     HOR_s(IJK,M) = HOR_s(IJK,M) + 4.183925d3*lHORs(M)
                  ENDDO
               ELSE
                  HOR_g(IJK) = HOR_g(IJK) + lHORg
                  DO M=1,MMAX
                     HOR_s(IJK,M) = HOR_s(IJK,M) + lHORs(M)
                  ENDDO
               ENDIF
            ELSE
! User-defined heat of reaction.
               HOR_g(IJK) = Reaction(H)%HoR(0) * RATES(H)
               DO M=1, MMAX
                  HOR_s(IJK,M) = Reaction(H)%HoR(M) * RATES(H)
               ENDDO
            ENDIF
         ENDIF

! Update rate of interphase mass transfer.
!---------------------------------------------------------------------//
          DO LM=1, (DIMENSION_LM+DIMENSION_M-1)
             R_PHASE(IJK,LM) = R_PHASE(IJK,LM) + &
                RATES(H) * Reaction(H)%rPHASE(LM)
          ENDDO
      ENDDO RXN_LP ! Loop over reactions.


! Calculate the toal rate of formation and consumption for each species.
!---------------------------------------------------------------------//
      IF(SPECIES_EQ(0)) THEN
         SUM_R_G(IJK) = SUM( &
            R_gp(IJK,:NMAX(0)) - &
            ROX_gc(IJK,:NMAX(0))*X_g(IJK,:NMAX(0)))
      ELSE
         DO H=1, NO_OF_RXNS
            IF(Reaction(H)%nPhases <= 0) CYCLE
            DO M=1, MMAX
               LM = 1 + ((M-1)*M)/2
               SUM_R_G(IJK) = SUM_R_G(IJK) + &
                  RATES(H) * Reaction(H)%rPHASE(LM)
            ENDDO
         ENDDO
      ENDIF

      DO M=1, MMAX
         IF(SPECIES_EQ(M)) THEN
            SUM_R_S(IJK,M) = SUM( &
               R_sp(IJK,M,:NMAX(M)) - &
               RoX_sc(IJK,M,:NMAX(M))*X_s(IJK,M,:NMAX(M)))
         ELSE
            DO H=1, NO_OF_RXNS
               IF(Reaction(H)%nPhases <= 0) CYCLE
               DO L=0, M-1
                  LM = 1 + L + ((M-1)*M)/2
                  SUM_R_S(IJK,M) = SUM_R_S(IJK,M) - &
                     RATES(H) * Reaction(H)%rPHASE(LM)
               ENDDO
               DO L=M+1, MMAX
                  LM = 1 + M + ((L-1)*L)/2
                  SUM_R_S(IJK,M) = SUM_R_S(IJK,M) + &
                     RATES(H) * Reaction(H)%rPHASE(LM)
               ENDDO
            ENDDO

         ENDIF
      ENDDO


      ENDIF  ! Fluid_At(IJK)
      END DO ! IJK

      RETURN

      END SUBROUTINE RRATES0 

!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization
!// 350 Changed do loop limits: 1,ijkmax2-> ijkstart3, ijkend3
