!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: CHECK_BC_DEM                                             !
! Author: J.Musser                                    Date: 01-Mar-14  !
!                                                                      !
! Purpose: Determine if BCs are "DEFINED" and that they contain the    !
! minimum amount of geometry data.                                     !
!                                                                      !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
      SUBROUTINE CHECK_BC_DEM(M_TOT)

! Global Variables:
!---------------------------------------------------------------------//
! User specified BC 
      use bc, only: BC_TYPE
! User specifed: BC geometry
      use bc, only: BC_EP_s
! Solids phase identifier
      use run, only: SOLIDS_MODEL
! Number of DEM inlet/outlet BCs detected.
      use des_bc, only: DEM_BCMI, DEM_BCMO
!
      use des_bc, only: DEM_BC_MI_MAP
      use des_bc, only: DEM_BC_MO_MAP
! Global Parameters:
!---------------------------------------------------------------------//
! The max number of BCs.
      use param, only: DIMENSION_BC
! Parameter constants
      use param1, only: ZERO

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager


      IMPLICIT NONE


! Passed Arguments:
!---------------------------------------------------------------------//
! Total number of solids phases.
      INTEGER, INTENT(in) :: M_TOT

! Local Variables:
!---------------------------------------------------------------------//
! loop/variable indices
      INTEGER :: BCV, M
!......................................................................!


! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_BC_DEM")

! Initialize
      DEM_BCMI = 0
      DEM_BCMO = 0

! Loop over all BCs looking for DEM solids inlets/outlets
      DO BCV = 1, DIMENSION_BC 

         SELECT CASE (TRIM(BC_TYPE(BCV)))

! Determine the number of mass inlets that contain DEM solids.
         CASE ('MASS_INFLOW')
            M_LP: DO M=1,M_TOT
               IF(SOLIDS_MODEL(M)=='DEM' .AND.                         &
                  BC_EP_s(BCV,M) > ZERO) THEN
                  DEM_BCMI = DEM_BCMI + 1
                  DEM_BC_MI_MAP(DEM_BCMI) = BCV
                  EXIT M_LP
               ENDIF
            ENDDO M_LP

! Count the number of pressure outflows.
         CASE ('P_OUTFLOW', 'CG_PO')
            DEM_BCMO = DEM_BCMO + 1
            DEM_BC_MO_MAP(DEM_BCMO) = BCV

         END SELECT

      ENDDO

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE CHECK_BC_DEM
