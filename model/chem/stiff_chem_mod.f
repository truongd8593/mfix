      MODULE STIFF_CHEM

! External Routines.
!---------------------------------------------------------------------//
! Routine used to calculate the reaction rates and populate the
! fluid variable ODEs for ODEPACK.
      external STIFF_CHEM_RRATES
! Routine used to compare to values.
      LOGICAL, external :: COMPARE
! Routine used to calculate species enthalpies.
      DOUBLE PRECISION, external :: CALC_H0


! Runtime Flags:
!---------------------------------------------------------------------//
! Flag to invoke stiff chemistry solver.
      LOGICAL :: STIFF_CHEMISTRY
! Flag to invoke the variable solids diameter model.
      LOGICAL :: CALL_GROW


! ODEPACK Controlling parameters:
!---------------------------------------------------------------------//
! Dimension of ODEs solved in ISAT or DI
      INTEGER :: ODE_DIMN
! Indicates type of Error control.
      INTEGER :: ODE_ITOL
! Relative error tolerance paramter.
      DOUBLE PRECISION, DIMENSION(1) :: ODE_RTOL
! Absolue error tolerance parameter. (Dimension (ODE_DIMN))
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: ODE_ATOL
! Declared length of RWORK.
      INTEGER :: ODE_LRW
! Declared length of IWORK.
      INTEGER :: ODE_LIW
! Jacobian type indicator.
      INTEGER :: ODE_JT


! Explicit interface for ODEPACK
!---------------------------------------------------------------------//
      INTERFACE
         SUBROUTINE DLSODA (F, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, &
            ITASK,ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, JT)
            external F
            INTEGER :: ITOL, ITASK, ISTATE, IOPT, LRW, LIW, JT
            INTEGER, dimension(2) :: NEQ
            INTEGER, dimension(LIW) :: IWORK
            DOUBLE PRECISION :: T, TOUT
            DOUBLE PRECISION :: JAC
            DOUBLE PRECISION, dimension(1) :: RTOL
            DOUBLE PRECISION, dimension(LRW) :: RWORK
            DOUBLE PRECISION, dimension(NEQ(1)) :: Y, ATOL
         END SUBROUTINE DLSODA
      END INTERFACE


! Legacy Variables:
!---------------------------------------------------------------------//
! Former keyword for invoking stiff solver.
      LOGICAL :: CALL_DI
! Keyword for using ISAT tables with stiff solver. (disabled)
      LOGICAL :: CALL_ISAT
! Time step for isat calculation. (disabled)
      DOUBLE PRECISION :: ISATdt

      contains

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                         C  
!     Module name: MCHEM_TIME_MARCH                                       C
!     Purpose: Called in time_march.f to do rxns calcs                    C
!                                                                         C
!     Author: Nan Xie                                   Date: 02-Aug-04   C
!     Reviewer:                                         Date:             C
!                                                                         C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE STIFF_CHEM_SOLVER(ODE_DT, iErr)

! Global Variables:
!---------------------------------------------------------------------//
      use funits,   only : DMP_LOG
      use output,   only : FULL_LOG
      use run,      only : TIME

      USE indices

      USE mpi_utility   

      implicit none

! Passed Variables:
!----------------------------------------------------------------------!
! Time integral length.
      DOUBLE PRECISION, intent(IN) :: ODE_DT
! Error Flag
      INTEGER, intent(OUT) :: iErr

! Local Variables:
!----------------------------------------------------------------------!

! Fluid Cell Index
      INTEGER :: IJK
! ODEs solved in ISAT or DI                             
      DOUBLE PRECISION, dimension(ODE_DIMN) :: ODE_VARS

! (1) :: Number of ODEs
! (2) :: Fluid cell index (IJK) passed into ODEPACK
      INTEGER, dimension(2) :: lNEQ
! Start time for integration
      DOUBLE PRECISION :: lT
! Stop time for integration
      DOUBLE PRECISION :: lTOUT
! Indicates type of Error control.
      INTEGER :: lITOL
! Relative error tolerance paramter.
      DOUBLE PRECISION :: lRTOL(1)
! Absolue error tolerance parameter. (Dimension (ODE_DIMN))
      DOUBLE PRECISION :: lATOL(ODE_DIMN)
! Index specifying the ODEPACK task.
      INTEGER :: lITASK
! Specifies the state of ODEPACK
      INTEGER :: lISTATE
! Flag indicating optional inputs are used.
      INTEGER :: lIOPT
! Array for REAL* work
      DOUBLE PRECISION :: RWORK(ODE_LRW)
! Declared length of RWORK.
      INTEGER :: lLRW
! Array for Integer work
      INTEGER :: IWORK(ODE_LIW)
! Declared length of IWORK.
      INTEGER :: lLIW
! Jacobain Routine (not used)
      DOUBLE PRECISION :: lJAC
! Jacobian type indicator.
      INTEGER :: lJT

      INCLUDE 'function.inc'

      lNEQ(1) = ODE_DIMN

      IJK_LP: DO IJK = IJKSTART3, IJKEND3
         IF(FLUID_AT(IJK)) THEN


! Map MFIX variables to ODE variables.
            CALL mapMFIXtoODE(IJK, ODE_VARS)

! Set the Fluid Cell IJK value
            lNEQ(2) = IJK

! Forced restset to initial values.
            lT    = 0.0d0
            lTOUT = ODE_DT
            lITOL = ODE_ITOL
            lRTOL = ODE_RTOL
            lATOL = ODE_ATOL
            lLRW  = ODE_LRW
            lLIW  = ODE_LIW
            lJT   = ODE_JT

! Fixed parameters
            lITASK  = 1
            lISTATE = 1
            lIOPT   = 0

            IF(CALC_REACTIONS(IJK)) THEN
! Integrate flow field variables to incorporate reactions.
 100           CALL DLSODA(STIFF_CHEM_RRATES, lNEQ, ODE_VARS, lT,      &
                  lTOUT, lITOL, lRTOL, lATOL, lITASK, lISTATE, lIOPT,  &
                  RWORK, lLRW, IWORK, lLIW, lJAC, lJT)

               IF(lISTATE == -1 ) THEN
                  write(*,"(/,3x,A,4(3x,I4))")                         &
                     'ODEPACK: Additional Iterations Needed --> ',     &
                     lISTATE, IJK, I_OF(IJK), J_OF(IJK)
                  write(*,"(3x,'lT: ',g15.8,3x,'lTOUT: ',g15.8)")      &
                     lT, lTOUT
                  lISTATE = 2
                  GoTo 100
               ELSEIF(lISTATE /= 2) THEN
!                  iErr = 2
!                  EXIT IJK_LP
!                  WRITE(*,1000) TIME, TIME+lTOUT, ODE_DT
                  write(*,"(/,3x,'ODEPACK: Failure --> ',4(3x,I4))")   &
                     lISTATE, IJK, I_OF(IJK), J_OF(IJK)
                  write(*,"(3x,'lT: ',g15.8,3x,'lTOUT: ',g15.8)")      &
                     lT, lTOUT
               ENDIF
            ENDIF
! Map ODE variables to MFIX variables.
            iErr = 0
            CALL mapODEtoMFIX(IJK, ODE_VARS)
         ENDIF
      END DO IJK_LP ! End Loop over fluiod Cells, IJK


! Error Handeling:
!---------------------------------------------------------------------//
!      IF(iErr == 2) THEN
!         IF (FULL_LOG) THEN
!            IF (myPE.EQ.PE_IO) WRITE(*,1000) TIME, TIME+lTOUT, ODE_DT
!         ENDIF 
!         WRITE(*,1000) TIME, TIME+lTOUT, ODE_DT
!         WRITE(*,"(' PE: ',I3,5x,'IJK: ',I6)") myPE, IJK
!
!         IF(WRITE_DASHBOARD) THEN
!            RUN_STATUS = 'Diverged/stalled...'
!            N_DASHBOARD = N_DASHBOARD + 1 
!            IF(MOD(N_DASHBOARD,F_DASHBOARD)==0) THEN
!               TLEFT = (TSTOP - TIME)*CPUOS 
!               CALL GET_TUNIT (TLEFT, TUNIT) 
!               CALL UPDATE_DASHBOARD(NIT,TLEFT,TUNIT)
!            ENDIF
!         ENDIF
!      ENDIF


      RETURN  

 1000 FORMAT(/,' Chemical Reaction diverged/stalled   :-(',/' Time =', &
         F10.4,3x,'ODE Time = ',F10.4, 3x, 'DT = ',G10.4,/)

      END SUBROUTINE STIFF_CHEM_SOLVER



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: mapMFIXtoODE                                           !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!  Author: J.Musser                                   Date: 07-Feb-13  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE mapMFIXtoODE(IJK, VARS)


! Global Variables:
!---------------------------------------------------------------------//
      use fldvar,   only : EP_g, RO_g, T_g, X_g, P_g
      use fldvar,   only : ROP_S, T_s, X_s, D_p
      use param1,   only : ONE
      use physprop, only : NMAX, C_pg, MW_g, MW_MIX_g
      use physprop, only : MMAX, C_ps, MW_s, RO_sv
      use run,      only : SPECIES_EQ

      implicit none

! Passed Variables:
!----------------------------------------------------------------------!
! Fluid cell index
      INTEGER, intent(in) :: IJK
! ODE Variables.
      DOUBLE PRECISION, dimension(ODE_DIMN), intent(out) :: VARS

! Local Variables:
!----------------------------------------------------------------------!
! Loop indicies:
      INTEGER :: M  ! phase
      INTEGER :: N  ! species
! ODE Equation Counter
      INTEGER :: Node

! Wrapper functions for solids phase volume fraction.
      INCLUDE 'ep_s1.inc'
      INCLUDE 'ep_s2.inc'

! Gas phase density.
      VARS(1) = RO_G(IJK)
! Gas phase temperature.
      VARS(2) = T_G(IJK)
! Gas phase species mass fractions.
      DO N=1,NMAX(0)
         Node = 2 + N
         VARS(Node) = X_G(IJK,N)
      ENDDO

      Node = 3+NMAX(0)
      DO M = 1, MMAX
         IF(SPECIES_EQ(M)) THEN
! Solids volume fraction.
            VARS(Node) = EP_S(IJK,M);      Node = Node + 1
! Solids temperature.
            VARS(Node) = T_S(IJK,M);       Node = Node + 1
! Solids phase species mass fractions.
            DO N=1,NMAX(M)
               VARS(Node) = X_S(IJK,M,N);  Node = Node + 1
            ENDDO
! Solids phase diameter.
            IF(CALL_GROW) THEN
               VARS(Node) = D_P(IJK,M);    Node = Node + 1
            ENDIF
         ENDIF
      ENDDO   

      RETURN
      END SUBROUTINE mapMFIXtoODE


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: mapODEtoMFIX                                           !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!  Author: J.Musser                                   Date: 07-Feb-13  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE mapODEtoMFIX(IJK, VARS)

      use constant, only : GAS_CONST
      use fldvar,   only : EP_g, RO_g, T_g, X_g, P_g
      use fldvar,   only : ROP_S, T_s, X_s, D_p
      use param1,   only : ONE
      use physprop, only : NMAX, C_pg, MW_g, MW_MIX_g
      use physprop, only : MMAX, C_ps, MW_s, RO_s, RO_sv
      use run,      only : SPECIES_EQ

      implicit none

      INTEGER, intent(in) :: IJK
      DOUBLE PRECISION, dimension(ODE_DIMN), intent(in) :: VARS

      INTEGER :: L, M, N, Node

! Gas phase density.
      RO_G(IJK) = VARS(1)
! Gas phase temperature.
      T_G(IJK) = VARS(2)
! Gas phase species mass fractions.
      DO N=1,NMAX(0)
         Node = 2 + N
         X_G(IJK,N) = VARS(Node)
      ENDDO

      Node = 3+NMAX(0)
      DO M = 1, MMAX
         IF(SPECIES_EQ(M)) THEN
! Solids volume fraction.
            ROP_S(IJK,M) = VARS(Node) * RO_S(M);  Node = Node + 1
! Solids temperature.
            T_S(IJK,M) = VARS(Node);              Node = Node + 1
! Solids phase species mass fractions.
            DO N=1,NMAX(M)
               X_S(IJK,M,N) = VARS(Node);         Node = Node + 1
            ENDDO
! Solids phase diameter.
            IF(CALL_GROW) THEN
               D_P(IJK,M)  = VARS(Node);          Node = Node + 1
            ENDIF
         ENDIF
      ENDDO   

! Calculate the gas volume fraction from solids volume fractions.
      EP_G(IJK) = ONE - sum(ROP_S(IJK,1:MMAX)/RO_S(1:MMAX))
! Calculate the mixture molecular weight.
      MW_MIX_G(IJK) = sum(X_G(IJK,1:NMAX(0))/MW_g(1:NMAX(0)))
      MW_MIX_G(IJK) = ONE/MW_MIX_G(IJK)
! Calculate the gas phase pressure.
      P_G(IJK) = (RO_G(IJK)*GAS_CONST*T_G(IJK))/MW_MIX_G(IJK)

      RETURN
      END SUBROUTINE mapODEtoMFIX


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: mapODEtoMFIX                                           !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!  Author: J.Musser                                   Date: 07-Feb-13  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      LOGICAL FUNCTION CALC_REACTIONS(IJK)

      use param1, only : SMALL_NUMBER
      use rxns,   only : NO_OF_RXNS

      implicit none

      INTEGER, intent(in) :: IJK

      DOUBLE PRECISION :: RATES(NO_OF_RXNS)

! Initialize
      RATES = 0.0d0

! Calculate user defined reaction rates.
      CALL USR_RATES(IJK, RATES)

! If there is little to no reaction in the cell, then set the ODE
! Time to zero to avoid calling the stiff solver.
      CALC_REACTIONS = .TRUE.
      if(COMPARE(sum(RATES),SMALL_NUMBER)) CALC_REACTIONS = .FALSE.

      RETURN
      END FUNCTION CALC_REACTIONS

      END MODULE STIFF_CHEM
