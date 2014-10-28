!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_d_t                                                C
!  Author: M. Syamlal                                 Date: 21-JUN-96  C
!  Purpose: Calculate coefficients linking velocity correction to      C
!           pressure correction -- Top                                 C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CALC_D_T(A_M, VXF_GS, VXF_SS, D_T, IER)

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
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(IN) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Volume x average at momentum cell centers
      DOUBLE PRECISION, INTENT(IN) :: VxF_gs(DIMENSION_3, DIMENSION_M)
! Volume x average at momentum cell centers Solid-Solid Drag
      DOUBLE PRECISION, INTENT(IN) :: VxF_ss(DIMENSION_3, DIMENSION_LM)
! Coefficients for pressure correction equation. These coefficients are
! initialized to zero in the subroutine time_march before the time loop.
      DOUBLE PRECISION, INTENT(INOUT) :: d_t(DIMENSION_3, 0:DIMENSION_M)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Temporary variable to store A matrix for local manipulation
      DOUBLE PRECISION :: AM0(DIMENSION_3, 0:DIMENSION_M)
      LOGICAL :: ANY_SOLIDS_Z_MOMENTUM
      INTEGER :: lMMAX
!-----------------------------------------------

! Initializing
      AM0(:,:) = A_M(:,0,:)

      IF(DES_CONTINUUM_COUPLED) THEN 
         lMMAX = merge(1, DES_MMAX, DES_INTERP_ON)
         AM0(:,0) = AM0(:,0) - sum(VXF_GDS(:,1:lMMAX))
         IF (DES_CONTINUUM_HYBRID) &
            AM0(:,:) = AM0(:,:) - sum(VXF_SDS(:,:,1:lMMAX))
      ENDIF

      ANY_SOLIDS_Z_MOMENTUM = any(MOMENTUM_Z_EQ(1:MMAX))

! Determine which calculations are needed
      IF (MOMENTUM_Z_EQ(0)) THEN
         IF(ANY_SOLIDS_Z_MOMENTUM)THEN
            CALL CALC_D_T_GAS_AND_SOLIDS(AM0, VXF_GS, VXF_SS, D_T, IER)
         ELSE
            CALL CALC_D_T_GAS_ONLY(AM0, VXF_GS, VXF_SS, D_T, IER)
         ENDIF
      ELSEIF (ANY_SOLIDS_Z_MOMENTUM) THEN
         CALL CALC_D_T_SOLIDS_ONLY(AM0, VXF_GS, VXF_SS, D_T, IER)
      ENDIF

      RETURN
      END SUBROUTINE CALC_D_T

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_d_t                                                C
!  Author: M. Syamlal                                 Date: 21-JUN-96  C
!  Purpose: Calculate coefficients linking velocity correction to      C
!           pressure correction -- Top                                 C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CALC_D_T_GAS_AND_SOLIDS(AM0, VXF_GS, VXF_SS, D_T, IER)

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
      DOUBLE PRECISION, INTENT(INOUT) :: d_t(DIMENSION_3, 0:DIMENSION_M)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! What phase momentum equation is activated?
!    pass1 is true when at least one solids phase momentum equation is solved
!    and the gas/fluid (m=0) phase momentum equation is solved
!    pass2 is true when at least one solids phase momentum equation is solved
!    but the gas/fluid (m=0) phase momentum equation is NOT solved
! Usual Indices
      INTEGER :: LM, M, L, LPL, LP
      INTEGER :: DM
      INTEGER :: INN
      INTEGER :: IJK, IJKT, K
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

! initializing
! Gas (m=0) and at least one solids phase (m>0) z-momentum equation are
! solved
! ---------------------------------------------------------------->>>
!!$omp   parallel do private(K,IJK, IJKT, EPGA, EPSA, AREA_FACE, &
!!$omp&  M, L, Lp, LpL, LM, SUM_VXF_SS, SUM_VXF_SS_wt_M, SUM_VXF_GS,  &
!!$omp&  DEN_MGas, NUM_MGas, DEN_MSol_LGas, DEN_MSol_LSol, &
!!$omp&  NUM_MSol_LSol, NUM_MSol_LGas),&
!!$omp&  schedule(static)
   DO IJK = ijkstart3, ijkend3

     IF (IP_AT_T(IJK) .OR. MFLOW_AT_T(IJK)) THEN   !impermeable
        DO M= 0, MMAX
           D_T(IJK,M) = ZERO
        ENDDO
     ELSE

        AREA_FACE = merge(ONE, AXY(IJK), CARTESIAN_GRID)
        K = K_OF(IJK)
        IJKT = TOP_OF(IJK)
        EPGA = AVG_Z(EP_G(IJK),EP_G(IJKT),K)

        SUM_VXF_GS = ZERO
        DO M= 1, MMAX
          EPSA(M) = AVG_Z(EP_S(IJK,M),EP_S(IJKT,M),K)
! Gas - All Solids VolxDrag summation
          SUM_VXF_GS = SUM_VXF_GS + VXF_GS(IJK,M)
          SUM_VXF_SS(M) = ZERO
          DO L = 1, MMAX
            IF (L .NE. M) THEN
              LM = FUNLM(L,M)
! Solid M - All other Solids VolxDrag summation
              SUM_VXF_SS(M) = SUM_VXF_SS(M) + VXF_SS(IJK,LM)
            ENDIF
          ENDDO
        ENDDO

! calculating DEN_MGas and NUM_MGas
        DEN_MGas  = ZERO
        NUM_MGas = ZERO
        DO M= 1, MMAX
          IF(MOMENTUM_Z_EQ(M)) THEN
            NUM_MGas = NUM_MGas + (EPSA(M)*VXF_GS(IJK,M)/&
              ((-AM0(IJK,M))+VXF_GS(IJK,M)+SUM_VXF_SS(M)+&
              SMALL_NUMBER))
            DEN_MGas = DEN_MGas + (VXF_GS(IJK,M)*&
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
             D_T(IJK,0) = P_SCALE*AREA_FACE/((-AM0(IJK,0))+DEN_MGas)
          ELSE
             D_T(IJK,0) = ZERO
          ENDIF
! Linking velocity correction coefficient to pressure - SOLIDs Phase
          DO M = 1, MMAX
             IF ( MOMENTUM_Z_EQ(M) .AND. ((-AM0(IJK,M))>SMALL_NUMBER .OR. &
                                          VXF_GS(IJK,M)>SMALL_NUMBER)) THEN
               D_T(IJK,M) = D_T(IJK,0)*(VXF_GS(IJK,M)/&
                            ((-AM0(IJK,M))+VXF_GS(IJK,M)))
             ELSE
               D_T(IJK,M) = ZERO
             ENDIF
          ENDDO

! Model A
! -------------------------------------------------------------------
        ELSE

! Linking velocity correction coefficient to pressure - GAS Phase
          IF ((-AM0(IJK,0))>SMALL_NUMBER .OR. &
                   DEN_MGas>SMALL_NUMBER) THEN
            D_T(IJK,0) = P_SCALE*AREA_FACE*(EPGA+NUM_MGas)/&
                         ((-AM0(IJK,0))+DEN_MGas)
          ELSE
            D_T(IJK,0) = ZERO
          ENDIF

          DO M= 1, MMAX
! calculating NUM_MSol_LGas and DEN_MSol_LGas
            NUM_MSol_LGas = VXF_GS(IJK,M)*EPGA/&
               ((-AM0(IJK,0))+SUM_VXF_GS+SMALL_NUMBER)
            DEN_MSol_LGas = VXF_GS(IJK,M)*(&
               ((-AM0(IJK,0)) + (SUM_VXF_GS - VXF_GS(IJK,M)))/&
               ((-AM0(IJK,0)) + SUM_VXF_GS+SMALL_NUMBER))

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

                IF(MOMENTUM_Z_EQ(L)) THEN
                  NUM_MSol_LSol(M) = NUM_MSol_LSol(M) + &
                     (VXF_SS(IJK,LM)*EPSA(L)/&
                     ((-AM0(IJK,L))+VXF_GS(IJK,L)+SUM_VXF_SS(L)+&
                     SMALL_NUMBER) )
                  DEN_MSol_LSol(M) = DEN_MSol_LSol(M) + VXF_SS(IJK,LM)*(&
                     ((-AM0(IJK,L))+VXF_GS(IJK,L)+SUM_VXF_SS_wt_M )/&
                     ((-AM0(IJK,L))+VXF_GS(IJK,L)+SUM_VXF_SS(L)+&
                     SMALL_NUMBER) )
                ELSE
                  DEN_MSol_LSol(M) = DEN_MSol_LSol(M) + VXF_SS(IJK,LM)
                ENDIF
              ENDIF  ! end if (l.ne.m)
            ENDDO  ! end do (l=1,mmax)

! Linking velocity correction coefficient to pressure - SOLIDs Phase
            IF ( MOMENTUM_Z_EQ(M) .AND. (-AM0(IJK,M)>SMALL_NUMBER .OR. &
                                       DEN_MSol_LGas>SMALL_NUMBER .OR. &
                                    DEN_MSol_LSol(M)>SMALL_NUMBER)) THEN
              D_T(IJK,M) = P_SCALE*AREA_FACE*&
                 (EPSA(M) + NUM_MSol_LSol(M) + NUM_MSol_LGas)/&
                 ((-AM0(IJK,M))+DEN_MSol_LGas+DEN_MSol_LSol(M))
            ELSE
              D_T(IJK,M) = ZERO
            ENDIF

          ENDDO   ! end do (m=1,mmax)

        ENDIF    !end if/else branch Model_B/Model_A

     ENDIF   ! end if/else branch(ip_at_t(ijk) .or. mflow_at_t(ijk))
   ENDDO  ! end do loop (ijk=ijkstart3,ijkend3)
! end pass1 branch
! ----------------------------------------------------------------<<<


 RETURN
 END SUBROUTINE CALC_D_T_GAS_AND_SOLIDS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_d_t                                                C
!  Author: M. Syamlal                                 Date: 21-JUN-96  C
!  Purpose: Calculate coefficients linking velocity correction to      C
!           pressure correction -- Top                                 C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CALC_D_T_GAS_ONLY(AM0, VXF_GS, VXF_SS, D_T, IER)

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
      DOUBLE PRECISION, INTENT(INOUT) :: d_t(DIMENSION_3, 0:DIMENSION_M)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! What phase momentum equation is activated?
!    pass1 is true when at least one solids phase momentum equation is solved
!    and the gas/fluid (m=0) phase momentum equation is solved
!    pass2 is true when at least one solids phase momentum equation is solved
!    but the gas/fluid (m=0) phase momentum equation is NOT solved
! Usual Indices
      INTEGER :: LM, M, L, LPL, LP
      INTEGER :: DM
      INTEGER :: INN
      INTEGER :: IJK, IJKT, K
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

! the solids z-momentum equations are not solved; only gas phase (m=0)
! z-momentum equation is solved. this branch is where a coupled DEM
! simulation should be directed for proper evaluation of pressure
! correction terms.
! ---------------------------------------------------------------->>>
!!$omp   parallel do private(IJK, K, IJKT, EPGA, M, SUM_VXF_GS,&
!!$omp&  AREA_FACE), &
!!$omp&  schedule(static)
   DO IJK = ijkstart3, ijkend3
     IF (IP_AT_T(IJK) .OR. MFLOW_AT_T(IJK)) THEN   !impermeable
       D_T(IJK,0) = ZERO
     ELSE

       AREA_FACE = merge(ONE, AXY(IJK), CARTESIAN_GRID)
       K = K_OF(IJK)
       IJKT = TOP_OF(IJK)
       EPGA = AVG_Z(EP_G(IJK),EP_G(IJKT),K)

       SUM_VXF_GS = ZERO
       IF (.NOT. QMOMK) THEN
          DO M= 1, MMAX
! Gas - All Solids VolxDrag summation
            SUM_VXF_GS = SUM_VXF_GS + VXF_GS(IJK,M)
          ENDDO
       ELSE
          DO INN = 1, QMOMK_NN
             DO M = 1, MMAX
                SUM_VXF_GS = SUM_VXF_GS + &
                  AVG_Z(QMOMK_F_GS(INN,IJK,M),QMOMK_F_GS(INN,IJKT,M),K)*&
                  VOL_W(IJK)
             ENDDO
          ENDDO
       ENDIF

       IF ((-AM0(IJK,0))>SMALL_NUMBER .OR. SUM_VXF_GS>SMALL_NUMBER) THEN
         IF (MODEL_B) THEN
            D_T(IJK,0) = P_SCALE*AREA_FACE/((-AM0(IJK,0))+SUM_VXF_GS)
         ELSE
            D_T(IJK,0) = P_SCALE*AREA_FACE*EPGA/((-AM0(IJK,0))+SUM_VXF_GS)
         ENDIF   !end if/else branch Model_B/Model_A
       ELSE
         D_T(IJK,0) = ZERO
       ENDIF

     ENDIF   ! end if/else branch(ip_at_t(ijk) .or. mflow_at_t(ijk))
   ENDDO   ! end do (ijk=ijkstart3,ijkend3)
! end only gas z-momentum equation is solved branch
! ----------------------------------------------------------------<<<



 RETURN
 END SUBROUTINE CALC_D_T_GAS_ONLY


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_d_t                                                C
!  Author: M. Syamlal                                 Date: 21-JUN-96  C
!  Purpose: Calculate coefficients linking velocity correction to      C
!           pressure correction -- Top                                 C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CALC_D_T_SOLIDS_ONLY(AM0, VXF_GS, VXF_SS, D_T, IER)

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
      DOUBLE PRECISION, INTENT(INOUT) :: d_t(DIMENSION_3, 0:DIMENSION_M)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! What phase momentum equation is activated?
!    pass1 is true when at least one solids phase momentum equation is solved
!    and the gas/fluid (m=0) phase momentum equation is solved
!    pass2 is true when at least one solids phase momentum equation is solved
!    but the gas/fluid (m=0) phase momentum equation is NOT solved
! Usual Indices
      INTEGER :: LM, M, L, LPL, LP
      INTEGER :: DM
      INTEGER :: INN
      INTEGER :: IJK, IJKT, K
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

! at least one solids phase z-momentum equation is solved but the gas
! phase (m=0) z-momentum equation is not solved
! ---------------------------------------------------------------->>>
!!$omp    parallel do private(IJK, K, IJKT, EPSA, L, Lp, M, SUM_VXF_SS, &
!!$omp&   NUM_MSol_LSol, DEN_MSol_LSol, SUM_VXF_SS_wt_M,AREA_FACE), &
!!$omp&  schedule(static)
   DO IJK = ijkstart3, ijkend3
     IF (IP_AT_T(IJK) .OR. MFLOW_AT_T(IJK) .OR. MODEL_B) THEN
       DO M= 1, MMAX
         D_T(IJK,M) = ZERO
       ENDDO
     ELSE

       AREA_FACE = merge(ONE, AXY(IJK), CARTESIAN_GRID)
       K = K_OF(IJK)
       IJKT = TOP_OF(IJK)

       DO M= 1, MMAX
         EPSA(M) = AVG_Z(EP_S(IJK,M),EP_S(IJKT,M),K)
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
             IF (MOMENTUM_Z_EQ(L)) THEN
               NUM_MSol_LSol(M) = NUM_MSol_LSol(M) + &
                 (VXF_SS(IJK,LM)*EPSA(L)/&
                 ((-AM0(IJK,L))+VXF_GS(IJK,L)+SUM_VXF_SS(L)+&
                 SMALL_NUMBER))
               DEN_MSol_LSol(M) = DEN_MSol_LSol(M) + VXF_SS(IJK,LM)*(&
                 ((-AM0(IJK,L))+VXF_GS(IJK,L)+SUM_VXF_SS_wt_M )/&
                 ((-AM0(IJK,L))+VXF_GS(IJK,L)+SUM_VXF_SS(L)+&
                 SMALL_NUMBER ))
              ELSE
                DEN_MSol_LSol(M) = DEN_MSol_LSol(M) + VXF_SS(IJK,LM)
              ENDIF
           ENDIF  ! end if (l.ne.m)
         ENDDO  ! end do (l=1,mmax)

! Linking velocity correction coefficient to pressure - SOLIDs Phase (Model_A only)
         IF (MOMENTUM_Z_EQ(M) .AND. (-AM0(IJK,M)>SMALL_NUMBER .OR. &
                                   VXF_GS(IJK,M)>SMALL_NUMBER .OR. &
                                DEN_MSol_LSol(M)>SMALL_NUMBER)) THEN
           D_T(IJK,M) = P_SCALE*AREA_FACE*&
             (EPSA(M) + NUM_MSol_LSol(M))/&
             ((-AM0(IJK,M))+VXF_GS(IJK,M)+DEN_MSol_LSol(M))
         ELSE
           D_T(IJK,M) = ZERO
         ENDIF
       ENDDO  ! end do (m=1,mmax)

     ENDIF   ! end if/else branch(ip_at_t(ijk) .or. mflow_at_t(ijk) .or. model_b)
   ENDDO   ! end do (ijk=ijkstart3,ijkend3)
! end pass2 branch
! ----------------------------------------------------------------<<<

 RETURN
 END SUBROUTINE CALC_D_T_SOLIDS_ONLY
