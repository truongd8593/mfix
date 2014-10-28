!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CALC_d_e                                                !
!  Author: M. Syamlal                                 Date: 21-JUN-96  !
!                                                                      !
!  Purpose: calculate coefficients linking velocity correction to      !
!  pressure correction -- East                                         !
!                                                                      !
!  Notes: MFIX convention: center coeff is negative, hence:            !
!                            (-A_M(IJK,0,M)) > or = 0                  !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE CALC_D_E(A_M, VXF_GS, VXF_SS, D_E, IER)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE fldvar
      USE geometry
      USE indices
      USE physprop
      USE run
      USE scales
      USE compar
      USE sendrecv
      USE cutcell
      USE qmom_kinetic_equation
      USE discretelement
      USE fun_avg
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Error index
      INTEGER, INTENT(INOUT) :: IER
! Septadiagonal matrix A_m.  The center coefficient is negative.
      DOUBLE PRECISION, INTENT(IN) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Volume x average at momentum cell centers
      DOUBLE PRECISION, INTENT(IN) :: VxF_gs(DIMENSION_3, DIMENSION_M)
! Volume x average at momentum cell centers Solid-Solid Drag
      DOUBLE PRECISION, INTENT(IN) :: VxF_ss(DIMENSION_3, DIMENSION_LM)
! Coefficients for pressure correction equation. These coefficients are
! initialized to zero in the subroutine time_march before the time loop.
! Note that the D_e coefficients for phases M>0 are generally not used
! unless the solids phase has close_packed=F, in which case the d_e
! coefficients for that phase are employed in a 'mixture' pressure
! correction equation and for correcting velocities
      DOUBLE PRECISION, INTENT(INOUT) :: d_e(DIMENSION_3, 0:DIMENSION_M)
!-----------------------------------------------
! Local variables
      INTEGER :: LMMAX
! Temporary variable to store A matrix for local manipulation
      DOUBLE PRECISION :: AM0(DIMENSION_3, 0:DIMENSION_M)

      LOGICAL :: ANY_SOLIDS_X_MOMENTUM
!-----------------------------------------------

! Initializing
      AM0(:,:) = A_M(:,0,:)

      IF(DES_CONTINUUM_COUPLED) THEN
         lMMAX = merge(1, DES_MMAX, DES_INTERP_ON)
         AM0(:,0) = AM0(:,0) - sum(VXF_GDS(:,1:lMMAX))
         IF (DES_CONTINUUM_HYBRID) &
            AM0(:,:) = AM0(:,:) - sum(VXF_SDS(:,:,1:lMMAX))
      ENDIF

      ANY_SOLIDS_X_MOMENTUM = any(MOMENTUM_X_EQ(1:MMAX))

      IF(MOMENTUM_X_EQ(0)) THEN
         IF(ANY_SOLIDS_X_MOMENTUM) THEN
            CALL CALC_D_E_GAS_AND_SOLIDS(AM0, VXF_GS, VXF_SS, D_E, IER)
         ELSE
            CALL CALC_D_E_GAS_ONLY(AM0, VXF_GS, VXF_SS, D_E, IER)
         ENDIF

      ELSEIF(ANY_SOLIDS_X_MOMENTUM) THEN
         CALL CALC_D_E_SOLIDS_ONLY(AM0, VXF_GS, VXF_SS, D_E, IER)
      ENDIF

      RETURN

      END SUBROUTINE CALC_D_E




!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CALC_D_E_GAS_AND_SOLIDS                                 !
!  Author: M. Syamlal                                 Date: 21-JUN-96  !
!                                                                      !
!  Purpose: calculate coefficients linking velocity correction to      !
!  pressure correction -- East                                         !
!                                                                      !
!  Notes: MFIX convention: center coeff is negative, hence:            !
!                            (-A_M(IJK,0,M)) > or = 0                  !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE CALC_D_E_GAS_AND_SOLIDS(AM0, VXF_GS, VXF_SS, D_E, IER)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE fldvar
      USE geometry
      USE indices
      USE physprop
      USE run
      USE scales
      USE compar
      USE sendrecv
      USE cutcell
      USE qmom_kinetic_equation
      USE discretelement
      USE fun_avg
      USE functions

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Error index
      INTEGER, INTENT(INOUT) :: IER
! Temporary variable to store A matrix for local manipulation
      DOUBLE PRECISION, INTENT(IN) :: AM0(DIMENSION_3, 0:DIMENSION_M)
! Volume x average at momentum cell centers
      DOUBLE PRECISION, INTENT(IN) :: VxF_gs(DIMENSION_3, DIMENSION_M)
! Volume x average at momentum cell centers Solid-Solid Drag
      DOUBLE PRECISION, INTENT(IN) :: VxF_ss(DIMENSION_3, DIMENSION_LM)
! Coefficients for pressure correction equation. These coefficients are
! initialized to zero in the subroutine time_march before the time loop.
! Note that the D_e coefficients for phases M>0 are generally not used
! unless the solids phase has close_packed=F, in which case the d_e
! coefficients for that phase are employed in a 'mixture' pressure
! correction equation and for correcting velocities
      DOUBLE PRECISION, INTENT(INOUT) :: d_e(DIMENSION_3, 0:DIMENSION_M)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Usual Indices
      INTEGER :: LM, M, L, LPL, LP
      INTEGER :: DM
      INTEGER :: INN
      INTEGER :: I, IJK, IJKE
! Average solid volume fraction at momentum cell centers
      DOUBLE PRECISION :: EPSA(DIMENSION_M)
! Average volume fraction at momentum cell centers
      DOUBLE PRECISION :: EPGA
! local alias for face area
      DOUBLE PRECISION :: AREA_FACE
! sum of All Solid-Gas VolxDrag
      DOUBLE PRECISION :: SUM_VXF_GS
! sum of Solid M - All other Solid drag
      DOUBLE PRECISION :: SUM_VXF_SS(DIMENSION_M)
! sum of Solid L - All other Solid VolxDrag but M
      DOUBLE PRECISION :: SUM_VXF_SS_wt_M
! numerator and demoninator needed for evaluation of pressure correction
! coefficient for GAS phase
      DOUBLE PRECISION :: DEN_MGas, NUM_MGas
! component of the total numerator needed for evaluation of pressure
! correction coefficient for SOLIDS phase M. specifically, for when the
! sum L is over solids phase only
      DOUBLE PRECISION :: NUM_MSol_LSol(DIMENSION_M)
! component of the total denominator needed for evaluation of pressure
! correction coefficient for SOLIDS phase M.  specifically, for when the
! sum L is over solids phases only
      DOUBLE PRECISION :: DEN_MSol_LSol(DIMENSION_M)
! component of the total numerator needed for evaluation of pressure
! correction coefficient for SOLIDS phase M.  specifically, for when the
! sum L is over gas phase only
     DOUBLE PRECISION :: NUM_MSol_LGas
! component of the total denominator needed for evaluation of pressure
! correction coefficient for SOLIDS phase M.  specifically, for when the
! sum L is over gas phase only
      DOUBLE PRECISION :: DEN_MSol_LGas
!-----------------------------------------------

! Gas (m=0) and at least one solids phase (m>0) x-momentum equation
! are solved
! ---------------------------------------------------------------->>>
!!$omp   parallel do private(I,IJK, IJKE, EPGA, EPSA, AREA_FACE, &
!!$omp&  M, L, Lp, LpL, LM, SUM_VXF_SS, SUM_VXF_SS_wt_M, SUM_VXF_GS, &
!!$omp&  DEN_MGas, NUM_MGas, DEN_MSol_LGas, DEN_MSol_LSol,&
!!$omp&  NUM_MSol_LGas, NUM_MSol_LSol),&
!!$omp&  schedule(static)
  DO IJK = ijkstart3, ijkend3

! impermeable
     IF (IP_AT_E(IJK) .OR. MFLOW_AT_E(IJK)) THEN
        DO M= 0, MMAX
          D_E(IJK,M) = ZERO
        ENDDO
     ELSE

        AREA_FACE = merge(ONE, AYZ(IJK), CARTESIAN_GRID)
        I = I_OF(IJK)
        IJKE = EAST_OF(IJK)
        EPGA = AVG_X(EP_G(IJK),EP_G(IJKE),I)

        SUM_VXF_GS = ZERO
        DO M= 1, MMAX
          EPSA(M) = AVG_X(EP_S(IJK,M),EP_S(IJKE,M),I)
! Gas - All Solids VolxDrag summation
          SUM_VXF_GS = SUM_VXF_GS + VXF_GS(IJK,M)
          SUM_VXF_SS(M) = ZERO
          DO L = 1, MMAX
            IF (L .NE. M) THEN
              LM = FUNLM(L,M)
! Solids M - All other Solids VolxDrag summation
              SUM_VXF_SS(M) = SUM_VXF_SS(M) + VXF_SS(IJK,LM)
            ENDIF
          ENDDO
        ENDDO

! calculating DEN_MGas and NUM_MGas
        DEN_MGas  = ZERO
        NUM_MGas = ZERO
        DO M= 1, MMAX
          IF (MOMENTUM_X_EQ(M)) THEN
            NUM_MGas = NUM_MGas + (EPSA(M)*VXF_GS(IJK,M)/&
              ((-AM0(IJK,M)) + VXF_GS(IJK,M) + SUM_VXF_SS(M) + &
              SMALL_NUMBER))
            DEN_MGas = DEN_MGas + (VXF_GS(IJK,M)* &
              ((-AM0(IJK,M))+SUM_VXF_SS(M))/&
              ((-AM0(IJK,M))+VXF_GS(IJK,M)+SUM_VXF_SS(M)+&
              SMALL_NUMBER))
          ELSE
            DEN_MGas = DEN_MGas + VXF_GS(IJK,M)
          ENDIF
        ENDDO

! Model B
! -------------------------------------------------------------------
        IF (MODEL_B) THEN
! Linking velocity correction coefficient to pressure - GAS Phase
          IF ( (-AM0(IJK,0))>SMALL_NUMBER .OR. &
                    DEN_MGas>SMALL_NUMBER ) THEN
             D_E(IJK,0) = P_SCALE*AREA_FACE/((-AM0(IJK,0))+DEN_MGas)
          ELSE
             D_E(IJK,0) = ZERO
          ENDIF
! Linking velocity correction coefficient to pressure - SOLIDs Phase
          DO M = 1, MMAX
            IF (MOMENTUM_X_EQ(M) .AND.((-AM0(IJK,M))>SMALL_NUMBER .OR. &
                                       VXF_GS(IJK,M)>SMALL_NUMBER)) THEN
              D_E(IJK,M) = D_E(IJK,0)*(VXF_GS(IJK,M)/&
                           ((-AM0(IJK,M))+VXF_GS(IJK,M)))
            ELSE
              D_E(IJK,M) = ZERO
            ENDIF
          ENDDO

! Model A
! -------------------------------------------------------------------
        ELSE

! Linking velocity correction coefficient to pressure - GAS Phase
          IF ( (-AM0(IJK,0))>SMALL_NUMBER .OR. &
                    DEN_MGas>SMALL_NUMBER ) THEN
             D_E(IJK,0) = P_SCALE*AREA_FACE*(EPGA+NUM_MGas)/&
                          ((-AM0(IJK,0))+DEN_MGas)
          ELSE
             D_E(IJK,0) = ZERO
          ENDIF


          DO M= 1, MMAX
! calculating NUM_MSol_LGas and DEN_MSol_LGas
            NUM_MSol_LGas = VXF_GS(IJK,M)*EPGA/&
               ((-AM0(IJK,0))+SUM_VXF_GS+SMALL_NUMBER)
            DEN_MSol_LGas = VXF_GS(IJK,M)*( &
               ((-AM0(IJK,0)) + (SUM_VXF_GS - VXF_GS(IJK,M)))/&
               ((-AM0(IJK,0)) + SUM_VXF_GS + SMALL_NUMBER) )

! calculating NUM_MSol_LSol and DEN_MSol_LSol
            NUM_MSol_LSol(M) = ZERO
            DEN_MSol_LSol(M)  = ZERO
            DO L = 1, MMAX
              IF (L .NE. M) THEN
                LM = FUNLM(L,M)
                SUM_VXF_SS_wt_M = ZERO
                DO Lp = 1, MMAX
                  IF ( (Lp .NE. L) .AND. (Lp .NE. M) ) THEN
                    LpL = FUNLM(Lp,L)
! Solids L - All other Solids VolxDrag but M summation
                    SUM_VXF_SS_wt_M = SUM_VXF_SS_wt_M + VXF_SS(IJK,LpL)
                  ENDIF
                ENDDO

                IF(MOMENTUM_X_EQ(L)) THEN
                  NUM_MSol_LSol(M) = NUM_MSol_LSol(M) + &
                     ( VXF_SS(IJK,LM)*EPSA(L)/&
                     ((-AM0(IJK,L))+VXF_GS(IJK,L)+SUM_VXF_SS(L)+&
                     SMALL_NUMBER) )
                  DEN_MSol_LSol(M) = DEN_MSol_LSol(M) + VXF_SS(IJK,LM)*(&
                     ((-AM0(IJK,L))+VXF_GS(IJK,L)+SUM_VXF_SS_wt_M)/&
                     ((-AM0(IJK,L))+VXF_GS(IJK,L)+SUM_VXF_SS(L)+ &
                     SMALL_NUMBER) )
                ELSE
                  DEN_MSol_LSol(M) = DEN_MSol_LSol(M) + VXF_SS(IJK,LM)
                ENDIF
              ENDIF  ! end if (l.ne.m)
            ENDDO  ! end do (l=1,mmax)

! Linking velocity correction coefficient to pressure - SOLIDs Phase
            IF (MOMENTUM_X_EQ(M) .AND. (-AM0(IJK,M)>SMALL_NUMBER .OR. &
                                      DEN_MSol_LGas>SMALL_NUMBER .OR. &
                                   DEN_MSol_LSol(M)>SMALL_NUMBER)) THEN
              D_E(IJK,M) = P_SCALE*AREA_FACE*&
                (EPSA(M) + NUM_MSol_LSol(M) + NUM_MSol_LGas)/&
                ((-AM0(IJK,M)) + DEN_MSol_LGas + DEN_MSol_LSol(M))
            ELSE
              D_E(IJK,M) = ZERO
            ENDIF
          ENDDO   ! end do (m=1,mmax)

        ENDIF    !end if/else branch Model_B/Model_A

     ENDIF   ! end if/else branch(ip_at_e(ijk) .or. mflow_at_e(ijk))
   ENDDO  ! end do loop (ijk=ijkstart3,ijkend3)
! end pass1 branch
! ----------------------------------------------------------------<<<


 RETURN
 END SUBROUTINE CALC_D_E_GAS_AND_SOLIDS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CALC_d_e                                                !
!  Author: M. Syamlal                                 Date: 21-JUN-96  !
!                                                                      !
!  Purpose: calculate coefficients linking velocity correction to      !
!  pressure correction -- East                                         !
!                                                                      !
!  Notes: The gas phase x-momentum equation is solved; NO solids phase !
!  x-momentum equations are solved.                                    !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE CALC_D_E_GAS_ONLY(AM0, VXF_GS, VXF_SS, D_E, IER)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE fldvar
      USE geometry
      USE indices
      USE physprop
      USE run
      USE scales
      USE compar
      USE sendrecv
      USE cutcell
      USE qmom_kinetic_equation
      USE discretelement
      USE fun_avg
      USE functions

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Error index
      INTEGER, INTENT(INOUT) :: IER
! Temporary variable to store A matrix for local manipulation
      DOUBLE PRECISION, INTENT(IN) :: AM0(DIMENSION_3, 0:DIMENSION_M)
! Volume x average at momentum cell centers
      DOUBLE PRECISION, INTENT(IN) :: VxF_gs(DIMENSION_3, DIMENSION_M)
! Volume x average at momentum cell centers Solid-Solid Drag
      DOUBLE PRECISION, INTENT(IN) :: VxF_ss(DIMENSION_3, DIMENSION_LM)
! Coefficients for pressure correction equation. These coefficients are
! initialized to zero in the subroutine time_march before the time loop.
! Note that the D_e coefficients for phases M>0 are generally not used
! unless the solids phase has close_packed=F, in which case the d_e
! coefficients for that phase are employed in a 'mixture' pressure
! correction equation and for correcting velocities
      DOUBLE PRECISION, INTENT(INOUT) :: d_e(DIMENSION_3, 0:DIMENSION_M)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Usual Indices
      INTEGER :: LM, M, L, LPL, LP
      INTEGER :: DM
      INTEGER :: INN
      INTEGER :: I, IJK, IJKE
! Average solid volume fraction at momentum cell centers
      DOUBLE PRECISION :: EPSA(DIMENSION_M)
! Average volume fraction at momentum cell centers
      DOUBLE PRECISION :: EPGA
! local alias for face area
      DOUBLE PRECISION :: AREA_FACE
! sum of All Solid-Gas VolxDrag
      DOUBLE PRECISION :: SUM_VXF_GS
! sum of Solid M - All other Solid drag
      DOUBLE PRECISION :: SUM_VXF_SS(DIMENSION_M)
! sum of Solid L - All other Solid VolxDrag but M
      DOUBLE PRECISION :: SUM_VXF_SS_wt_M
! numerator and demoninator needed for evaluation of pressure correction
! coefficient for GAS phase
      DOUBLE PRECISION :: DEN_MGas, NUM_MGas
! component of the total numerator needed for evaluation of pressure
! correction coefficient for SOLIDS phase M. specifically, for when the
! sum L is over solids phase only
      DOUBLE PRECISION :: NUM_MSol_LSol(DIMENSION_M)
! component of the total denominator needed for evaluation of pressure
! correction coefficient for SOLIDS phase M.  specifically, for when the
! sum L is over solids phases only
      DOUBLE PRECISION :: DEN_MSol_LSol(DIMENSION_M)
! component of the total numerator needed for evaluation of pressure
! correction coefficient for SOLIDS phase M.  specifically, for when the
! sum L is over gas phase only
     DOUBLE PRECISION :: NUM_MSol_LGas
! component of the total denominator needed for evaluation of pressure
! correction coefficient for SOLIDS phase M.  specifically, for when the
! sum L is over gas phase only
      DOUBLE PRECISION :: DEN_MSol_LGas
!-----------------------------------------------

! ---------------------------------------------------------------->>>
!!$omp   parallel do private(IJK, I, IJKE, EPGA, M, SUM_VXF_GS,&
!!$omp&  AREA_FACE), &
!!$omp&  schedule(static)
      DO IJK = ijkstart3, ijkend3

         IF (IP_AT_E(IJK) .OR. MFLOW_AT_E(IJK)) THEN   !impermeable
            D_E(IJK,0) = ZERO
            CYCLE
         ENDIF

         AREA_FACE = merge(ONE, AYZ(IJK), CARTESIAN_GRID)

         I = I_OF(IJK)
         IJKE = EAST_OF(IJK)
         EPGA = AVG_X(EP_G(IJK),EP_G(IJKE),I)

         SUM_VXF_GS = ZERO
         IF (.NOT. QMOMK) THEN
            DO M= 1, MMAX
! Gas - All Solids VolxDrag summation
               SUM_VXF_GS = SUM_VXF_GS + VXF_GS(IJK,M)
            ENDDO
         ELSE
            DO INN = 1, QMOMK_NN
               DO M = 1, MMAX
                  SUM_VXF_GS = SUM_VXF_GS + VOL_U(IJK)* &
                     AVG_X(QMOMK_F_GS(INN,IJK,M),QMOMK_F_GS(INN,IJKE,M),I)
               ENDDO
            ENDDO
         ENDIF

         IF ((-AM0(IJK,0))>SMALL_NUMBER .OR. SUM_VXF_GS>SMALL_NUMBER) THEN
            IF (MODEL_B) THEN
               D_E(IJK,0) = P_SCALE*AREA_FACE/((-AM0(IJK,0))+SUM_VXF_GS)
            ELSE
               D_E(IJK,0) = P_SCALE*AREA_FACE*EPGA/((-AM0(IJK,0))+SUM_VXF_GS)
            ENDIF   !end if/else branch Model_B/Model_A
         ELSE
            D_E(IJK,0) = ZERO
         ENDIF

      ENDDO   ! end do (ijk=ijkstart3,ijkend3)

      RETURN
      END SUBROUTINE CALC_D_E_GAS_ONLY


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CALC_D_E_SOLIDS_ONLY                                    !
!  Author: M. Syamlal                                 Date: 21-JUN-96  !
!                                                                      !
!  Purpose: Calculate coefficients linking velocity correction to      !
!  pressure correction -- East                                         !
!                                                                      !
!  Notes: The gas phase momentum equations are NOT solved but at least !
!  one solids phase momentum equation is solved.                       !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE CALC_D_E_SOLIDS_ONLY(AM0, VXF_GS, VXF_SS, D_E, IER)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE fldvar
      USE geometry
      USE indices
      USE physprop
      USE run
      USE scales
      USE compar
      USE sendrecv
      USE cutcell
      USE qmom_kinetic_equation
      USE discretelement
      USE fun_avg
      USE functions

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Error index
      INTEGER, INTENT(INOUT) :: IER
! Temporary variable to store A matrix for local manipulation
      DOUBLE PRECISION, INTENT(IN) :: AM0(DIMENSION_3, 0:DIMENSION_M)
! Volume x average at momentum cell centers
      DOUBLE PRECISION, INTENT(IN) :: VxF_gs(DIMENSION_3, DIMENSION_M)
! Volume x average at momentum cell centers Solid-Solid Drag
      DOUBLE PRECISION, INTENT(IN) :: VxF_ss(DIMENSION_3, DIMENSION_LM)
! Coefficients for pressure correction equation. These coefficients are
! initialized to zero in the subroutine time_march before the time loop.
! Note that the D_e coefficients for phases M>0 are generally not used
! unless the solids phase has close_packed=F, in which case the d_e
! coefficients for that phase are employed in a 'mixture' pressure
! correction equation and for correcting velocities
      DOUBLE PRECISION, INTENT(INOUT) :: d_e(DIMENSION_3, 0:DIMENSION_M)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! What phase momentum equation is activated?
!    pass1 is true when at least one solids phase momentum equation is solved
!    and the gas/fluid (m=0) phase momentum equation is solved
!    pass2 is true when at least one solids phase momentum equation is solved
!    but the gas/fluid (m=0) phase momentum equation is NOT solved
      LOGICAL :: Pass1, Pass2
! Usual Indices
      INTEGER :: LM, M, L, LPL, LP
      INTEGER :: DM
      INTEGER :: INN
      INTEGER :: I, IJK, IJKE
! Average solid volume fraction at momentum cell centers
      DOUBLE PRECISION :: EPSA(DIMENSION_M)
! Average volume fraction at momentum cell centers
      DOUBLE PRECISION :: EPGA
! local alias for face area
      DOUBLE PRECISION :: AREA_FACE
! sum of All Solid-Gas VolxDrag
      DOUBLE PRECISION :: SUM_VXF_GS
! sum of Solid M - All other Solid drag
      DOUBLE PRECISION :: SUM_VXF_SS(DIMENSION_M)
! sum of Solid L - All other Solid VolxDrag but M
      DOUBLE PRECISION :: SUM_VXF_SS_wt_M
! numerator and demoninator needed for evaluation of pressure correction
! coefficient for GAS phase
      DOUBLE PRECISION :: DEN_MGas, NUM_MGas
! component of the total numerator needed for evaluation of pressure
! correction coefficient for SOLIDS phase M. specifically, for when the
! sum L is over solids phase only
      DOUBLE PRECISION :: NUM_MSol_LSol(DIMENSION_M)
! component of the total denominator needed for evaluation of pressure
! correction coefficient for SOLIDS phase M.  specifically, for when the
! sum L is over solids phases only
      DOUBLE PRECISION :: DEN_MSol_LSol(DIMENSION_M)
! component of the total numerator needed for evaluation of pressure
! correction coefficient for SOLIDS phase M.  specifically, for when the
! sum L is over gas phase only
     DOUBLE PRECISION :: NUM_MSol_LGas
! component of the total denominator needed for evaluation of pressure
! correction coefficient for SOLIDS phase M.  specifically, for when the
! sum L is over gas phase only
      DOUBLE PRECISION :: DEN_MSol_LGas
!-----------------------------------------------

! at least one solids phase x-momentum equation is solved but the gas
! phase (m=0) x-momentum equation is not solved

! ---------------------------------------------------------------->>>
!!$omp    parallel do private(IJK, I, IJKE, EPSA, L, Lp, M, SUM_VXF_SS, &
!!$omp&   NUM_MSol_LSol, DEN_MSol_LSol, SUM_VXF_SS_wt_M, AREA_FACE), &
!!$omp&  schedule(static)
   DO IJK = ijkstart3, ijkend3
     IF (IP_AT_E(IJK) .OR. MFLOW_AT_E(IJK) .OR. MODEL_B) THEN
       DO M= 1, MMAX
         D_E(IJK,M) = ZERO
       ENDDO
     ELSE

       AREA_FACE = merge(ONE, AYZ(IJK), CARTESIAN_GRID)
       I = I_OF(IJK)
       IJKE = EAST_OF(IJK)

       DO M= 1, MMAX
         EPSA(M) = AVG_X(EP_S(IJK,M),EP_S(IJKE,M),I)
         SUM_VXF_SS(M) = ZERO
         DO L = 1, MMAX
           IF (L .NE. M) THEN
             LM = FUNLM(L,M)
! Solids M - All other Solids VolxDrag summation
             SUM_VXF_SS(M) = SUM_VXF_SS(M) + VXF_SS(IJK,LM)
           ENDIF
         ENDDO
       ENDDO

       DO M= 1, MMAX
! calculating NUM_MSol_LSol and DEN_MSol_LSol
         NUM_MSol_LSol(M) = ZERO
         DEN_MSol_LSol(M)  = ZERO
         DO L = 1, MMAX
           IF (L .NE. M) THEN
             LM = FUNLM(L,M)
             SUM_VXF_SS_wt_M = ZERO
             DO Lp = 1, MMAX
               IF ( (Lp .NE. L) .AND. (Lp .NE. M) ) THEN
                 LpL = FUNLM(Lp,L)
! Solids L - All other Solids VolxDrag but M summation
                 SUM_VXF_SS_wt_M = SUM_VXF_SS_wt_M + VXF_SS(IJK,LpL)
               ENDIF
             ENDDO
             IF (MOMENTUM_X_EQ(L)) THEN
               NUM_MSol_LSol(M) = NUM_MSol_LSol(M) + &
                  (VXF_SS(IJK,LM)*EPSA(L)/&
                  ((-AM0(IJK,L))+VXF_GS(IJK,L)+SUM_VXF_SS(L)+&
                  SMALL_NUMBER))
               DEN_MSol_LSol(M) = DEN_MSol_LSol(M) + VXF_SS(IJK,LM)*(&
                  ( (-AM0(IJK,L))+VXF_GS(IJK,L)+SUM_VXF_SS_wt_M )/&
                  ( (-AM0(IJK,L))+VXF_GS(IJK,L)+SUM_VXF_SS(L)+&
                  SMALL_NUMBER ))
             ELSE
               DEN_MSol_LSol(M) = DEN_MSol_LSol(M) + VXF_SS(IJK,LM)
             ENDIF
           ENDIF   ! end if (L.ne.M)
         ENDDO  ! end do (l=1,mmax)

! Linking velocity correction coefficient to pressure - SOLIDs Phase (Model_A only)
         IF (MOMENTUM_X_EQ(M) .AND. (-AM0(IJK,M)>SMALL_NUMBER .OR. &
                                   VXF_GS(IJK,M)>SMALL_NUMBER .OR. &
                                DEN_MSol_LSol(M)>SMALL_NUMBER)) THEN
           D_E(IJK,M) = P_SCALE*AREA_FACE*&
             (EPSA(M) + NUM_MSol_LSol(M))/&
             ((-AM0(IJK,M))+VXF_GS(IJK,M)+DEN_MSol_LSol(M))
         ELSE
           D_E(IJK,M) = ZERO
         ENDIF
       ENDDO  ! end do (m=1,mmax)


     ENDIF   ! end if/else branch(ip_at_e(ijk) .or. mflow_at_e(ijk) .or. model_b)
   ENDDO   ! end do (ijk=ijkstart3,ijkend3)
! end pass2 branch
! ----------------------------------------------------------------<<<


 RETURN
 END SUBROUTINE CALC_D_E_SOLIDS_ONLY

