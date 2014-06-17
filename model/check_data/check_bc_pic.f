!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
! minimum amount of geometry data.                                     !
!                                                                      !
! Subroutine: CHECK_BC_PIC                                             !
! Author: R. Garg                                     Date: 11-Jun-14  !
!                                                                      !
! Purpose: Determine if BCs are "DEFINED" and that they contain the    !
! minimum amount of geometry data.                                     !
!                                                                      !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
      SUBROUTINE CHECK_BC_PIC(M_TOT)

! Global Variables:
!---------------------------------------------------------------------//
! User specified BC 
      use bc, only: BC_TYPE
! User specifed: BC geometry
      use bc, only: BC_EP_s
! Use specified flag for ignoring PO BC for discrete solids 
      USE bc, only: BC_PO_APPLY_TO_DES
! PIC model specific BC region specification. 
      USE bc, only: BC_PIC_MI_CONST_NPC, BC_PIC_MI_CONST_STATWT
! Solids phase identifier
      use run, only: SOLIDS_MODEL
! Number of PIC inlet/outlet BCs detected.
      use pic_bc, only: PIC_BCMI, PIC_BCMO
!
      use pic_bc, only: PIC_BCMI_MAP
      use pic_bc, only: PIC_BCMO_MAP
! Global Parameters:
!---------------------------------------------------------------------//
! The max number of BCs.
      use param, only: DIMENSION_BC
! Parameter constants
      use param1, only: ZERO, UNDEFINED


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
      INTEGER :: BCV, M, BCV_I

! Temp logical variables for checking constant npc and statwt specification
      LOGICAL :: CONST_NPC, CONST_STATWT
!......................................................................!


! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_BC_PIC")

! Initialize
      PIC_BCMI = 0
      PIC_BCMO = 0

! Loop over all BCs looking for PIC solids inlets/outlets
      DO BCV = 1, DIMENSION_BC

         SELECT CASE (TRIM(BC_TYPE(BCV)))

! Determine the number of mass inlets that contain PIC solids.
         CASE ('MASS_INFLOW')
            M_LP: DO M=1,M_TOT
               IF(SOLIDS_MODEL(M)=='PIC' .AND.                         &
                  BC_EP_s(BCV,M) > ZERO) THEN
                  PIC_BCMI = PIC_BCMI + 1
                  PIC_BCMI_MAP(PIC_BCMI) = BCV
                  EXIT M_LP
               ENDIF
            ENDDO M_LP

! Count the number of pressure outflows.
         CASE ('P_OUTFLOW')
            IF(BC_PO_APPLY_TO_DES(BCV)) then
               PIC_BCMO = PIC_BCMO + 1
               PIC_BCMO_MAP(PIC_BCMO) = BCV
            ENDIF

! Flag CG_MI as an error if PIC solids are present.
         CASE ('CG_MI')
            DO M=1,M_TOT
               IF(SOLIDS_MODEL(M)=='PIC') THEN
                  IF(BC_EP_s(BCV,M) /= UNDEFINED .AND.                 &
                     BC_EP_s(BCV,M) > ZERO) THEN
                     WRITE(ERR_MSG,1000) trim(iVar('BC_TYPE',BCV)),    &
                        'GC_MI'
                     CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
                  ENDIF
               ENDIF
            ENDDO

         CASE ('CG_PO')
            WRITE(ERR_MSG,1000) trim(iVar('BC_TYPE',BCV)), 'GC_PO'
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

         CASE ('MASS_OUTFLOW', 'OUTFLOW', 'P_INFLOW')
            WRITE(ERR_MSG,1000) trim(iVar('BC_TYPE',BCV)),             &
               trim(BC_TYPE(BCV))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

         END SELECT

      ENDDO

      CALL FINL_ERR_MSG

   
 1000 FORMAT('Error 1000: Unsupported boundary condition specified ',  &
         'with',/'PIC simulation: ',A,' = ',A,/'Please correct the ',&
         'mfix.dat file.')


! Loop over all MI BC's for data consistency checks 
      DO BCV_I = 1, PIC_BCMI

! Get the user defined BC ID.
         BCV = PIC_BCMI_MAP(BCV_I)
      
         DO M=1,M_TOT
            IF(SOLIDS_MODEL(M)=='PIC' .AND.                         &
                 BC_EP_s(BCV,M) > ZERO) THEN
               CONST_NPC    = (BC_PIC_MI_CONST_NPC   (BCV, M) .ne. 0)
               CONST_STATWT = (BC_PIC_MI_CONST_STATWT(BCV, M) .ne. ZERO  )
               IF(CONST_NPC.and.CONST_STATWT) then
                  WRITE(ERR_MSG, 1100) BCV, M 
                  CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
               ENDIF
               
               IF(.not.CONST_NPC.and.(.not.CONST_STATWT)) then 
                  WRITE(ERR_MSG, 1101) BCV, M 
                  CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
               ENDIF
               

            ENDIF
         ENDDO
         
      ENDDO

 1100 FORMAT('Error 1100: In MPPIC model for BC # ',i5, & 
      ' and solid phase # ', i5, /, & 
      'Non zero Values specified for both ', &
      'BC_PIC_MI_CONST_NPC and BC_PIC_MI_CONST_STATWT.', /, &
      'Choose between constant number of parcels per cell or ', &
      'constant statistical weight', /, &
      'See MFIX readme',/'Please correct the data file.')


 1101 FORMAT('Error 1101: In MPPIC model for BC # ',i5, & 
      ' and solid phase # ', i5, /, & 
      'A non-zero value not specified for ', &
      'BC_PIC_MI_CONST_NPC or BC_PIC_MI_CONST_STATWT. ', /, &
      'Choose between constant number of parcels per cell or ', &
      'constant statistical weight', /, &
      'See MFIX readme',/'Please correct the data file.')

      END SUBROUTINE CHECK_BC_PIC
