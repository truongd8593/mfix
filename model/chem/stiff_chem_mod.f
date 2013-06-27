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
! Flag indicating if cell IJK is own by myPE.
      LOGICAL, dimension(:), allocatable :: notOwner

! ODEPACK Controlling parameters:
!---------------------------------------------------------------------//
! Dimension of ODEs solved in stiff solver.
      INTEGER :: ODE_DIMN_all
! Dimension of ODEs solved in stiff solver for gas phase only.
      INTEGER :: ODE_DIMN_g

! Dimension of ODEs solved in stiff solver for gas phase only.
      INTEGER :: NEQ_DIMN

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

      use mpi_utility   
      use stiff_chem_dbg

      implicit none

! Passed Variables:
!----------------------------------------------------------------------!
! Time integral length.
      DOUBLE PRECISION, intent(IN) :: ODE_DT
! Error Flag
      INTEGER, intent(OUT) :: iErr

! Local Variables:
!----------------------------------------------------------------------!
! Error flag -> Integration failed in one or more fluid cells.
      LOGICAL :: lErr_l   ! local
      LOGICAL :: gErr_l   ! global

! Fluid Cell Index
      INTEGER :: IJK
! The maximum number of ODEs to solve.
      DOUBLE PRECISION, dimension(ODE_DIMN_all) :: ODE_VARS

! (1) :: Number of ODEs
! (2) :: Fluid cell index (IJK) passed into ODEPACK
! (:) :: Flag for solving solid phase M density and species equations.
      INTEGER, dimension(NEQ_DIMN) :: lNEQ
! Start time for integration
      DOUBLE PRECISION :: lT
! Stop time for integration
      DOUBLE PRECISION :: lTOUT
! Indicates type of Error control.
      INTEGER :: lITOL
! Relative error tolerance paramter.
      DOUBLE PRECISION :: lRTOL(1)
! Absolue error tolerance parameter. (Dimension (ODE_DIMN))
      DOUBLE PRECISION :: lATOL(ODE_DIMN_all)
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

! The number of attempts of a specific fluid cell.
      INTEGER :: lAtps

      LOGICAL :: lReset

      INCLUDE 'function.inc'

      lErr_l = .FALSE.

      IF(FULL_LOG) CALL INIT_ODE_STATS()

      IJK_LP: DO IJK = IJKSTART3, IJKEND3
         IF(notOwner(IJK)) cycle IJK_LP
         IF(FLUID_AT(IJK)) THEN

            lAtps = 0
            lReset = .FALSE.

! Forced restset of tolerance values.
            lRTOL = ODE_RTOL
            lATOL = ODE_ATOL

! Increment the attempts counter.
  50        lAtps = lAtps + 1

! Forced restset to initial values.
            lT    = 0.0d0
            lTOUT = ODE_DT
            lITOL = ODE_ITOL
            lLRW  = ODE_LRW
            lLIW  = ODE_LIW
            lJT   = ODE_JT

! Fixed parameters
            lITASK  = 1
            lISTATE = 1
            lIOPT   = 1

! Calculate the number of ODEs to solve.
            CALL CALC_ODE_COEFF(lNEQ, IJK)

! Clear the work arrays.
            IWORK = 0
            RWORK = ZERO

! The maximum number of internal steps ODEPACK may use to integrate over
! the time interval. The default value is 500.
            IWORK(6) = 500000

            IF(CALC_REACTIONS(IJK)) THEN


!               write(*,"(3x,'Reaction in ',I4)") IJK

! Map MFIX variables to ODE variables.
               CALL mapMFIXtoODE(lNEQ, ODE_VARS)

! Store a copy of the original field variables. This allows for these
! values to be 'reset' in the event that the stiff solver fails.
               CALL ODE_UPDATE_OLD(ODE_DIMN_all, ODE_VARS)

! Clear the error flag.
 100           iErr = 0

! Integrate flow field variables to incorporate reactions.
               CALL DLSODA(STIFF_CHEM_RRATES, lNEQ, ODE_VARS, lT,      &
                  lTOUT, lITOL, lRTOL, lATOL, lITASK, lISTATE, lIOPT,  &
                  RWORK, lLRW, IWORK, lLIW, lJAC, lJT)

! Verify that the results are well defined.
               CALL CHECK_ODE_DATA(NEQ_DIMN, lNEQ, ODE_DIMN_all,       &
                  ODE_VARS, lISTATE, iErr)

! Successfully Integrated ODEs.
               IF(iErr == 0) THEN
                  lReset = .FALSE.
! Additional integration steps are needed (lT < lTOUT).
               ELSEIF(iErr == -1) THEN
! Reste the state flag and keep integrating.
                     lISTATE = 2
                     goto 100

! Too much accuracy was requested.
               ELSEIF(iErr == -2) THEN
                  IF(lAtps < 3) THEN
! Write to the error log file.
                     IF(ODE_DEBUG_LEVEL >= 2) CALL WRITE_ODE_LOG(iErr, &
                        NEQ_DIMN, lNEQ, ODE_DIMN_all, ODE_VARS)
! Reset the ODE variable array.
                     CALL ODE_RESET(ODE_DIMN_all, ODE_VARS, lAtps)
! Reset the field variables.
                     CALL mapODEtoMFIX(lNEQ, ODE_VARS)
! Loosen the convergence criteria and try again.
                     lRTOL = ODE_RTOL*10.0d0
                     lATOL = ODE_ATOL*10.0d0
                     goto 50
                  ELSE
! Write to the error log file.
                     IF(ODE_DEBUG_LEVEL >= 1) CALL WRITE_ODE_LOG(iErr, &
                        NEQ_DIMN, lNEQ, ODE_DIMN_all, ODE_VARS)
! Set the flag to reset the field variables to the initial values.
                     lReset = .TRUE.
                  ENDIF
! All other errors.
               ELSE
! Tighten the convergence criteria and try again.
                  IF(lAtps < 3) THEN
! Write to the error log file.
                     IF(ODE_DEBUG_LEVEL >= 2) CALL WRITE_ODE_LOG(iErr, &
                        NEQ_DIMN, lNEQ, ODE_DIMN_all, ODE_VARS)
! Reset the ODE variable array.
                     CALL ODE_RESET(ODE_DIMN_all, ODE_VARS, lAtps)
! Rest the filed variables to their original values.
               CALL mapODEtoMFIX(lNEQ, ODE_VARS)
! Reduce the tolerances and try again.
                     lRTOL = ODE_RTOL/(10.0d0**lAtps)
                     lATOL = ODE_ATOL/(10.0d0**lAtps)
                     goto 50
                  ELSE
! Write to the error log file.
                     IF(ODE_DEBUG_LEVEL >= 1) CALL WRITE_ODE_LOG(iErr, &
                        NEQ_DIMN, lNEQ, ODE_DIMN_all, ODE_VARS)
! Set the flag to reset the field variables to the initial values.
                     lReset = .TRUE.
                  ENDIF

               ENDIF ! IF(iErr == 0)

! Reset the field variables.
               if(lReset) CALL ODE_RESET(ODE_DIMN_all, ODE_VARS, lAtps)
! Store the results in the field variables.
               CALL mapODEtoMFIX(lNEQ, ODE_VARS)
! Collect solver stats.
               IF(FULL_LOG) CALL UPDATE_ODE_STATS(lNEQ, NEQ_DIMN,      &
                  IWORK(11), ODE_DIMN_all, lAtps)


            ENDIF  ! EndIF CALC_REACTIONS
         ENDIF  ! IF(CALC_REACTIONS(IJK))
      END DO IJK_LP ! End Loop over fluiod Cells, IJK

!      gErr_l = .FALSE.
!      CALL GLOBAL_ALL_OR(lErr_l, gErr_l)
!      IF(gErr_l) CALL WRITE_VTU_FILE


      CALL FINALIZE_STIFF_SOLVER()

      IF(FULL_LOG) CALL WRITE_ODE_STATS()

      iErr = 0

      RETURN
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
      SUBROUTINE mapMFIXtoODE(lNEQ, VARS)


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
      INTEGER, intent(in) :: lNEQ(NEQ_DIMN)
! ODE Variables.
      DOUBLE PRECISION, intent(out) :: VARS(ODE_DIMN_all)

! Local Variables:
!----------------------------------------------------------------------!
! Loop indicies:
      INTEGER :: M  ! phase
      INTEGER :: N  ! species
! ODE Equation Counter
      INTEGER :: Node

! Fluid cell index
      INTEGER :: IJK

! Wrapper functions for solids phase volume fraction.
      INCLUDE 'ep_s1.inc'
      INCLUDE 'ep_s2.inc'

! Initialize.
      Node = 1
      IJK = lNEQ(2)
      VARS = 0.0d0

! Gas phase density.
      VARS(Node) = RO_G(IJK);              Node = Node + 1
! Gas phase temperature.
      VARS(Node) = T_G(IJK);               Node = Node + 1
! Gas phase species mass fractions.
      DO N=1,NMAX(0)
         VARS(Node) = X_G(IJK,N);          Node = Node + 1
      ENDDO

! Solids temperature.
      DO M = 1, MMAX
         VARS(Node) = T_S(IJK,M);          Node = Node + 1
      ENDDO

      DO M = 1, MMAX
         IF(lNEQ(2+M) == 1) THEN
! Solids volume fraction.
            VARS(Node) = ROP_S(IJK,M);     Node = Node + 1
! Solids phase species mass fractions.
            DO N=1,NMAX(M)
               VARS(Node) = X_S(IJK,M,N);  Node = Node + 1
            ENDDO
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
      SUBROUTINE mapODEtoMFIX(lNEQ, VARS)

      use constant, only : GAS_CONST
      use fldvar,   only : EP_g, RO_g, T_g, X_g, P_g
      use fldvar,   only : ROP_S, T_s, X_s, D_p
      use param1,   only : ONE
      use physprop, only : NMAX, C_pg, MW_g, MW_MIX_g
      use physprop, only : MMAX, C_ps, MW_s, RO_s, RO_sv

      use compar,   only : myPE, PE_IO

      implicit none

      INTEGER, intent(in) :: lNEQ(NEQ_DIMN)
      DOUBLE PRECISION, intent(in) :: VARS(ODE_DIMN_all)

! Fluid Cell index.
      INTEGER :: IJK
      INTEGER :: L, M, N, Node

      INTEGER :: countNaN
      LOGICAL :: writeMsg

      IJK = lNEQ(2)

      countNaN = 0
      writeMsg = .FALSE.
      NaN_lp: do l=1, ode_dimn_all
         if(isNaN(VARS(l))) then
            countNaN = countNan + 1
            writeMsg = .TRUE.
         endif
      enddo NaN_lp


      if(writeMsg) then
         write(*,"(3x,'From MapODEtoMFIX: NaNs Found! :: ',3(3x,I4))") myPE, IJK, countNaN
         if(countNaN < ode_dimn_all) then
            do l=1, ode_dimn_all
               if(isNaN(VARS(l))) write(*,"(5x,' NaN in Var ',I2)") l
            enddo

            write(*,"(3x,'ODE Flags: ',I3)")lNEQ(1)
            DO M=1, MMAX
               write(*,"(3x,'Phase ',I3,': ',I3,3x,'ROPs: ',g18.6)")M, lNEQ(2+M), ROP_s(IJK,M)
            enddo
         endif
         stop
      endif


      Node = 1

! Gas phase density.
      RO_G(IJK) = VARS(Node);                     Node = Node + 1
! Gas phase temperature.
      T_G(IJK) = VARS(Node);                      Node = Node + 1
! Gas phase species mass fractions.
      DO N=1,NMAX(0)
         X_G(IJK,N) = VARS(Node);                 Node = Node + 1
      ENDDO

! Solids temperature.
      DO M = 1, MMAX
         IF(ROP_s(IJK,M) > 1.0d-8) &
            T_S(IJK,M) = VARS(Node);              Node = Node + 1
      ENDDO   

! Only map back what was calculated.
      DO M = 1, MMAX
         IF(lNEQ(2+M) == 1) THEN
! Solids volume fraction. (Constant Solids Density)
            ROP_S(IJK,M) = VARS(Node);            Node = Node + 1
! Solids phase species mass fractions.
            DO N=1,NMAX(M)
               X_S(IJK,M,N) = VARS(Node);         Node = Node + 1
            ENDDO
         ENDIF
      ENDDO   

! Calculate the gas volume fraction from solids volume fractions. Only
! update it's value if the solids equations are being solved.
      IF(sum(lNEQ(3:)) > 0) EP_G(IJK) = &
         ONE - sum(ROP_S(IJK,1:MMAX)/RO_S(1:MMAX))

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

      use rxns,   only : NO_OF_RXNS

      implicit none

      INTEGER, intent(in) :: IJK

      DOUBLE PRECISION :: RATES(NO_OF_RXNS)

      DOUBLE PRECISION, parameter :: rLimit = 1.0d-8

! Initialize
      RATES = 0.0d0

! Calculate user defined reaction rates.
      CALL USR_RATES(IJK, RATES)

! If there is little to no reaction in the cell, then set the ODE
! Time to zero to avoid calling the stiff solver.
      CALC_REACTIONS = .TRUE.
      if(maxval(RATES) < rLimit) CALC_REACTIONS = .FALSE.


      RETURN
      END FUNCTION CALC_REACTIONS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: CALC_DIMN_ODE                                          !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!  Author: J.Musser                                   Date: 07-Feb-13  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_ODE_COEFF(lNEQ, IJK)

      use fldvar,     only : ROP_s
      use physprop,   only : MMAX, RO_sv
      use run,        only : SPECIES_EQ

      implicit none

      INTEGER, intent(in)  :: IJK
      INTEGER, intent(out) :: lNEQ(NEQ_DIMN)

      INTEGER :: M

      LOGICAL :: USE_SOLIDS_ODEs

      INCLUDE 'ep_s1.inc'
      INCLUDE 'ep_s2.inc'

! Initialize.
      USE_SOLIDS_ODEs = .FALSE.
      lNEQ(2) = IJK
      lNEQ(3:) = 0

!       write(*,"(/3x,'IJK: ',I4)") IJK
!       DO M=1,MMAX
!          write(*,"(3x,'Phase: ',I2)") M
!          write(*,"(3x,'EP_s: ',g18.6)") EP_s(IJK,M) 
!          write(*,"(3x,'ROP_s:',g18.6)") ROP_S(IJK,M)
!          write(*,"(3x,'RO_SV:',g18.6)") RO_SV(IJK,M)
!       ENDDO

! If there is little to no solids in the cell, then set the ODE
! dimension to gas phase only.
      DO M=1, MMAX
         IF(SPECIES_EQ(M) .AND. (EP_s(IJK,M) > 1.0d-6)) THEN
            lNEQ(2+M) = 1
            USE_SOLIDS_ODEs = .TRUE.
         ENDIF
      ENDDO

      IF(USE_SOLIDS_ODEs)THEN
         lNEQ(1) = ODE_DIMN_all
      ELSE
         lNEQ(1) = ODE_DIMN_g
         lNEQ(3:) = 0
      ENDIF


      RETURN
      END SUBROUTINE CALC_ODE_COEFF


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: FINALIZE_STIFF_SOLVER                                  !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!  Author: J.Musser                                   Date: 07-Feb-13  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE FINALIZE_STIFF_SOLVER

      use constant, only : GAS_CONST
      use fldvar,   only : EP_g, RO_g, T_g, X_g, P_g
      use fldvar,   only : ROP_S, T_s, X_s
      use physprop, only : NMAX
      use physprop, only : MMAX

      use param1,   only : ONE
      use physprop, only : MW_g, MW_MIX_g



      use compar
      use mpi_utility
      use sendrecv

      implicit none

      INTEGER :: IJK  ! Fluid Cell index.
      INTEGER :: M    ! Solids phase index
      INTEGER :: N    ! Species index


      CALL send_recv(RO_G,2)
      CALL send_recv(T_G,2)

      DO N=1,NMAX(0)
         CALL send_recv(X_G(:,N),2)
      ENDDO

      DO M = 1, MMAX
! Solids temperature.
         CALL send_recv(T_S(:,M),2)
! Solids volume fraction. (Constant Solids Density)
         CALL send_recv(ROP_S(:,M),2)
! Solids phase species mass fractions.
         DO N=1,NMAX(M)
            CALL send_recv(X_S(:,M,N),2)
         ENDDO
      ENDDO   

      DO IJK = ijkStart3, ijkEnd3
! Calculate the mixture molecular weight.
         MW_MIX_G(IJK) = sum(X_G(IJK,1:NMAX(0))/MW_g(1:NMAX(0)))
         MW_MIX_G(IJK) = ONE/MW_MIX_G(IJK)
! Calculate the gas phase pressure.
         P_G(IJK) = (RO_G(IJK)*GAS_CONST*T_G(IJK))/MW_MIX_G(IJK)
      ENDDO

      RETURN
      END SUBROUTINE FINALIZE_STIFF_SOLVER

      END MODULE STIFF_CHEM
