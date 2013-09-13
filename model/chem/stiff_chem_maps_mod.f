      MODULE STIFF_CHEM_MAPS

      PRIVATE

! Variable Access:
!---------------------------------------------------------------------//

! Subroutine Access:
!---------------------------------------------------------------------//
      PUBLIC :: mapMFIXtoODE,   &
                mapODEtoMFIX

      contains



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: mapODEtoMFIX                                           !
!                                                                      !
!  Purpose: This is a driver routine for mapping variables stored in   !
!  the ODE array back to MFIX field variables.                         !
!                                                                      !
!  Author: J.Musser                                   Date: 07-Feb-13  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE mapMFIXtoODE(lnD, lNEQ, loD, lVars)

! Global Variables:
!---------------------------------------------------------------------//
      use fldvar, only: ROP_g
      use fldvar, only: ROP_s
      use fldvar, only: T_g
      use fldvar, only: T_s
      use fldvar, only: X_g
      use fldvar, only: X_s

      use physprop, only: NMAX
      use physprop, only: MMAX


      implicit none


! Passed Variables:
!---------------------------------------------------------------------//
! Passed array dimensions
      INTEGER, intent(in) :: lnD  ! lNEQ
      INTEGER, intent(in) :: loD  ! lVars

! (1) Number of ODEs to be solve
! (2) Fluid cell index
      INTEGER, intent(in) :: lNEQ(lnD)
! Array of dependent variable initial values.
      DOUBLE PRECISION, intent(out)  :: lVars(loD)


! Local Variables:
!---------------------------------------------------------------------//
! Fluid cell index
      INTEGER :: IJK
! Loop indicies:
      INTEGER :: M  ! phase
      INTEGER :: N  ! species
! ODE Equation Counter
      INTEGER :: Node


! Initialize.
      IJK = lNEQ(2)

      lVars = 0.0d0

      Node = 1

! Gas phase density.
      lVars(Node) = ROP_G(IJK);             Node = Node + 1
! Gas phase temperature.
      lVars(Node) = T_G(IJK);               Node = Node + 1
! Gas phase species mass fractions.
      DO N=1,NMAX(0)
         lVars(Node) = X_G(IJK,N);          Node = Node + 1
      ENDDO

! Solids temperature.
      DO M = 1, MMAX
         lVars(Node) = T_S(IJK,M);          Node = Node + 1
      ENDDO

      DO M = 1, MMAX
         IF(lNEQ(2+M) == 1) THEN
! Solids volume fraction.
            lVars(Node) = ROP_S(IJK,M);     Node = Node + 1
! Solids phase species mass fractions.
            DO N=1,NMAX(M)
               lVars(Node) = X_S(IJK,M,N);  Node = Node + 1
            ENDDO
         ENDIF
      ENDDO   

      RETURN
      END SUBROUTINE mapMFIXtoODE

!#######################################################################################################################



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: mapODEtoMFIX                                           !
!                                                                      !
!  Purpose: This is a driver routine for mapping variables stored in   !
!  the ODE array back to MFIX field variables.                         !
!                                                                      !
!  Author: J.Musser                                   Date: 07-Feb-13  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE mapODEtoMFIX(lnD, lNEQ, loD, lVars)

! Global Variables:
!---------------------------------------------------------------------//
      use fldvar,   only : ROP_g
      use fldvar,   only : ROP_S
      use fldvar,   only : T_g
      use fldvar,   only : T_s
      use fldvar,   only : X_g
      use fldvar,   only : X_s

      use physprop, only : MMAX
      use physprop, only : NMAX

      use compar,   only : myPE
      use compar,   only : PE_IO


      implicit none


! Passed Variables:
!---------------------------------------------------------------------//
! (1) Number of ODEs to be solve
! (2) Fluid cell index
      INTEGER, intent(in) :: lnD
      INTEGER, intent(in) :: lNEQ(lnD)

! Array of dependent variable initial values.
      INTEGER, intent(in) :: loD
      DOUBLE PRECISION, intent(in)  :: lVars(loD)


! Local Variables:
!---------------------------------------------------------------------//
! Fluid Cell index.
      INTEGER :: IJK
      INTEGER :: L

! Loop indicies:
      INTEGER :: M    ! phase
      INTEGER :: N    ! species
      INTEGER :: Node ! ODE Equation Counter

      INTEGER :: countNaN
      LOGICAL :: writeMsg


      IJK = lNEQ(2)

      countNaN = 0
      writeMsg = .FALSE.
      NaN_lp: do l=1, loD
         if(lVars(l).NE.lVars(l)) then
            countNaN = countNan + 1
            writeMsg = .TRUE.
         endif
      enddo NaN_lp

      if(writeMsg) then
         write(*,"(3x,'From MapODEtoMFIX: NaNs Found! :: ',3(3x,I4))") &
            myPE, IJK, countNaN

         if(countNaN < loD) then
            do l=1, loD
               if(lVars(l).NE.lVars(l))                                &
                  write(*,"(5x,' NaN in Var ',I2)") l
            enddo
         endif
      endif

! Directly map the ODE values into the field variable names.
!-----------------------------------------------------------------------
! Initialize the loop counter for ODEs.
      Node = 1

! Gas phase density.
      ROP_G(IJK) = lVars(Node);                    Node = Node + 1
! Gas phase temperature.
      T_G(IJK) = lVars(Node);                      Node = Node + 1
! Gas phase species mass fractions.
      DO N=1,NMAX(0)
         X_G(IJK,N) = lVars(Node);                 Node = Node + 1
      ENDDO

! Solids temperature.
      DO M = 1, MMAX
         IF(ROP_s(IJK,M) > 1.0d-8) &
            T_S(IJK,M) = lVars(Node);              Node = Node + 1
      ENDDO   

! Only map back what was calculated.
      DO M = 1, MMAX
         IF(lNEQ(2+M) == 1) THEN
! Solids volume fraction. (Constant Solids Density)
            ROP_S(IJK,M) = lVars(Node);            Node = Node + 1
! Solids phase species mass fractions.
            DO N=1,NMAX(M)
               X_S(IJK,M,N) = lVars(Node);         Node = Node + 1
            ENDDO
         ENDIF
      ENDDO   

!      IF(VARIABLE_DENSITY) THEN

!      ELSE
         CALL mapODEtoMFIX_ROs0
!      ENDIF

      RETURN

      contains
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: mapODEtoMFIX                                           !
!                                                                      !
!  Purpose: Finish mapping ODE result into MFIX field variables. This  !
!           routine is for constant solids density.                    !
!                                                                      !
!  Author: J.Musser                                   Date: 07-Feb-13  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE mapODEtoMFIX_ROs0  ! Constant solids density

      use constant, only : GAS_CONST

      use fldvar,   only : P_g
      use fldvar,   only : RO_g

      use fldvar,   only : EP_g
      use physprop, only : MW_g
      use physprop, only : MW_MIX_g
      use physprop, only : MW_s

      use param1,   only : ONE
      use param1,   only : LARGE_NUMBER
      use param1,   only : SMALL_NUMBER

      use physprop, only : RO_s
      use physprop, only : RO_sv

!      use physprop, only : C_pg
!      use physprop, only : C_ps



      implicit none

! Calculate the gas volume fraction from solids volume fractions. Only
! update it's value if the solids equations are being solved.
      IF(sum(lNEQ(3:)) > 0) EP_G(IJK) = &
         ONE - sum(ROP_S(IJK,1:MMAX)/RO_S(1:MMAX))

! Gas phase bulk density is updated within the stiff solver (lVar(1)).
! Now that the gas phase volume fraction is updated, the gas phase 
! density can be backed out. RO_g * EP_g = ROP_g
      IF(EP_g(IJK) > small_number) THEN
         RO_g(IJK) = ROP_g(IJK) / EP_g(IJK)
      ELSE
! This case shouldn't happen, however 'LARGE_NUMBER' is used to aid
! in tracking errors should this somehow become and issue.
         RO_g(IJK) = LARGE_NUMBER
      ENDIF

! Calculate the mixture molecular weight.
      MW_MIX_G(IJK) = sum(X_G(IJK,1:NMAX(0))/MW_g(1:NMAX(0)))
      MW_MIX_G(IJK) = ONE/MW_MIX_G(IJK)

! Calculate the gas phase pressure.
      P_G(IJK) = (RO_G(IJK)*GAS_CONST*T_G(IJK))/MW_MIX_G(IJK)

      RETURN
      END SUBROUTINE mapODEtoMFIX_ROs0

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: mapODEtoMFIX_ROs                                       !
!                                                                      !
!  Purpose: Finish mapping ODE result into MFIX field variables. This  !
!           routine is for constant solids density.                    !
!                                                                      !
!  Author: J.Musser                                   Date: 07-Feb-13  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE mapODEtoMFIX_ROs  ! Variable solids density
      END SUBROUTINE mapODEtoMFIX_ROs  ! Variable solids density


      END SUBROUTINE mapODEtoMFIX
      END MODULE STIFF_CHEM_MAPS
