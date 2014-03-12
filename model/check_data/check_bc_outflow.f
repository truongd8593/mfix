!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: CHECK_BC_OUTFLOW                                         !
! Author: J.Musser                                    Date: 01-Mar-14  !
!                                                                      !
! Purpose: Provided a detailed error message when the sum of volume    !
!                                                                      !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
      SUBROUTINE CHECK_BC_OUTFLOW(M_TOT, SKIP, BCV)

      USE param, only: DIM_M
      USE param1, only: UNDEFINED
      USE param1, only: ZERO
      use physprop, only: RO_g0

      use discretelement
      use physprop

      use bc

      use error_manager

      IMPLICIT NONE

! loop/variable indices
      INTEGER, INTENT(in) :: BCV
      INTEGER, INTENT(in) :: M_TOT
      LOGICAL, INTENT(in) :: SKIP(DIM_M)

      INTEGER :: M
      DOUBLE PRECISION :: SUM_EP

      CALL INIT_ERR_MSG("CHECK_BC_OUTFLOW")


! if bc_ep_g is defined at a PO, MO or O boundary, then the sum of ep_g
! and ep_s at the boundary may not equal one given the following code
! in the subroutine set_outflow (see code for details). therefore if
! bc_ep_g is is defined and any solids are present, provide the user
! with a warning.
      IF (BC_EP_G(BCV) /= UNDEFINED) THEN
         SUM_EP = BC_EP_G(BCV)

! Unclear how the discrete solids volume fraction can be dictated at 
! the boundary, so it is currently prevented!
         IF (DISCRETE_ELEMENT) THEN
            WRITE(ERR_MSG, 1101) trim(iVar('BC_EP_g',BCV))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

 1101 FORMAT('Error 1101: ',A,' should not be specified for outflow ', &
         'BCs',/'with DEM/PIC runs. Please correct the mfix.dat file.')

         ELSE

! by this point the code has checked that no discrete solids are present
! otherwise the routine will have exited
            DO M = 1, M_TOT 
               IF(BC_ROP_S(BCV,M) == UNDEFINED) THEN
                  IF(BC_EP_G(BCV) == ONE) THEN
! what does it mean to force the bulk density to zero at the
! boundary? (again does this value matter anyway)
                     BC_ROP_S(BCV,M) = ZERO

                  ELSEIF(M_TOT == 1 ) THEN
! no discrete solids are present so a bulk density can be defined from 
! 1-bc_ep_g even for hybrid model
                     BC_ROP_S(BCV,M) = (ONE - BC_EP_G(BCV))*RO_S0(M)
                  ELSE
! bc_ep_g is defined but some bc_rop_s(m) are undefined.
! in this ep_p in the outflow boundary will be based on the user defined
! value of bc_ep_g, while rop_s would become based on the adjacent fluid
! cell. consequently, no check ensures the result is consistent with
! a requirement for ep_g+ep_s=1.
                     WRITE(ERR_MSG, 1102) trim(BC_TYPE(BCV)), &
                        trim(iVar('BC_EP_g',BCV))

                     CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

 1102 FORMAT('Warning 1102: Volume fraction may not sum to one for ',  &
         A,/'when ',A,' is defined.')

                  ENDIF
               ENDIF  ! end if(bc_rop_s(bcv,m) == undefined)
! by this point bc_rop_s should either be defined or mfix exited
! therefore we can check that sum of void fraction and solids volume 
! fractions
               SUM_EP = SUM_EP + BC_ROP_S(BCV,M)/RO_S0(M) 
            ENDDO

         ENDIF   ! end if/else (.not.discrete_element)
      ENDIF

      CALL FINL_ERR_MSG

      RETURN

 1000 FORMAT('Error 1000: Required input not specified: ',A,/'Please ',&
         'correct the mfix.dat file.')

      END SUBROUTINE CHECK_BC_OUTFLOW



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: CHECK_BC_P_OUTFLOW                                       !
! Author: J.Musser                                    Date: 01-Mar-14  !
!                                                                      !
! Purpose: Provided a detailed error message when the sum of volume    !
!                                                                      !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
      SUBROUTINE CHECK_BC_P_OUTFLOW(M_TOT, SKIP, BCV)

      USE param, only: DIM_M
      USE param1, only: UNDEFINED
      USE param1, only: ZERO
      use physprop, only: RO_g0
      use bc, only: BC_P_g

      use error_manager

      IMPLICIT NONE

! loop/variable indices
      INTEGER, INTENT(in) :: BCV
      INTEGER, INTENT(in) :: M_TOT
      LOGICAL, INTENT(in) :: SKIP(DIM_M)


      CALL INIT_ERR_MSG("CHECK_BC_P_OUTFLOW")


      IF (BC_P_G(BCV) == UNDEFINED) THEN 
         WRITE(ERR_MSG,1000) 'BC_P_g', BCV 
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

      ELSEIF (BC_P_G(BCV)<=ZERO .AND. RO_G0==UNDEFINED) THEN 
         WRITE(ERR_MSG, 1100) BCV, BC_P_G(BCV)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF 

 1100 FORMAT('Error 1100: Pressure must be greater than zero for ',    &
         'compressible flow',/' >>>  BC_P_g(',I3,') = ',A,/'Please ',  &
         'correct the mfix.dat file.')

! Clean up and return.
      CALL FINL_ERR_MSG

      RETURN

 1000 FORMAT('Error 1000: Required input not specified: ',A,/'Please ',&
         'correct the mfix.dat file.')

      END SUBROUTINE CHECK_BC_P_OUTFLOW


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: CHECK_BC_MASS_OUTFLOW                                    !
! Author: J.Musser                                    Date: 01-Mar-14  !
!                                                                      !
! Purpose: Provided a detailed error message when the sum of volume    !
!                                                                      !
! Comments:                                                            !
!     The velocities at the outflow face are fixed and the momentum    !
!     equations are not solved in the outflow cells. Since the flow    !
!     is out of the domain none of the other scalars should need to    !
!     be specified (e.g., mass fractions, void fraction, etc.,).       !
!     Such values will become defined according to their adjacent      !
!     fluid cell                                                       !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_BC_MASS_OUTFLOW(M_TOT, SKIP, BCV)

      USE param 
      USE param1 
      USE geometry
      USE fldvar
      USE physprop
      USE run
      USE bc
      USE indices
      USE funits 
      USE scalars
      USE compar
      USE sendrecv
      USE discretelement
      USE mfix_pic
      USE cutcell

      use error_manager

      IMPLICIT NONE

! error flag
      LOGICAL :: ERROR
! loop/variable indices
      INTEGER, intent(in) :: BCV
      INTEGER, intent(in) :: M_TOT
      LOGICAL, intent(in) :: SKIP(DIM_M)

! Loop variable
      INTEGER :: M


      CALL INIT_ERR_MSG("CHECK_BC_MASS_OUTFLOW")


      IF(BC_DT_0(BCV) == UNDEFINED) THEN 
         WRITE(ERR_MSG, 1000) trim(iVar('BC_DT_0',BCV))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

! This check probably needs changed.
      IF(RO_G0 == UNDEFINED .AND. (BC_P_G(BCV) == UNDEFINED .OR.       &
         BC_T_G(BCV) == UNDEFINED) .AND.BC_MASSFLOW_G(BCV) /= ZERO) THEN

         IF(BC_PLANE(BCV)=='W' .OR. BC_PLANE(BCV)=='E') THEN
            IF(BC_U_G(BCV) /= ZERO) THEN 
               WRITE(ERR_MSG, 1100) BCV, 'BC_U_g'
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)  
            ENDIF 
         ELSEIF(BC_PLANE(BCV)=='N' .OR. BC_PLANE(BCV)=='S') THEN 
            IF(BC_V_G(BCV) /= ZERO) THEN 
               WRITE(ERR_MSG, 1100) BCV, 'BC_V_g' 
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)  
            ENDIF 
         ELSEIF (BC_PLANE(BCV)=='T' .OR. BC_PLANE(BCV)=='B') THEN 
            IF(BC_W_G(BCV) /= ZERO) THEN 
               WRITE(ERR_MSG, 1100)  BCV, 'BC_W_g' 
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)  
            ENDIF
         ENDIF
      ENDIF   ! end if/else (ro_g0 /=undefined)

 1100 FORMAT('Error 1100: Invalid mass outflow boundary condition: ',  &
         I3,/'RO_g0, BC_P_g, and BC_T_g are UNDEFINED and ',A,' is ',  &
         'non-zero',/'Please correct the mfix.dat file.')


      CALL FINL_ERR_MSG

      RETURN

 1000 FORMAT(/1X,70('*')//' From: CHECK_DATA_07',/' Message: ',A,'(',I2,&
         ') not specified',/1X,70('*')/) 

      END SUBROUTINE CHECK_BC_MASS_OUTFLOW
