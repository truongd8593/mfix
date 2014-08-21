!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_d_e                                                C
!  Purpose: calculate coefficients linking velocity correction to      C
!           pressure correction -- East                                C
!                                                                      C
!  Notes: MFIX convention: center coeff is negative, hence:            C
!                            (-A_M(IJK,0,M)) > or = 0                  C
!         The EAST face area is AYZ                                    C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-JUN-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number: 2                                                  C
!  Purpose: allow multiparticle in D_E calculation in accounting for   C
!     the averaged Solid-Solid drag and multi-solid-particle drag      C
!  Author: S. Dartevelle, LANL                        Date: 28-FEB-04  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number: 3                                                  C
!  Purpose: more accurate formulas in the D_s(M) for Model_A           C
!  Author: S. Dartevelle, LANL                        Date: 17-MAR-04  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number: 4                                                  C
!  Purpose: To incorporate Cartesian grid modifications                C
!  Author: Jeff Dietiker                              Date: 01-Jul-09  C
!                                                                      C
!  Revision Number: 5                                                  C
!  Purpose: Incorporation of QMOM for the solution of the particle     C
!           kinetic equation                                           C
!  Author: Alberto Passalacqua - Fox Research Group   Date: 02-Dec-09  C
!                                                                      C
!  Literature/Document References:                                     C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

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
! Number of solids to use as LOOP counter
      INTEGER :: LC_MMAX
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
! Temporary variable to store A matrix for local manipulation
      DOUBLE PRECISION :: AM0(DIMENSION_3, 0:DIMENSION_M)
!-----------------------------------------------
! Include statement functions
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'
!-----------------------------------------------

! Initializing
      Pass1 = .FALSE.
      Pass2 = .FALSE.
      AM0(:,:) = A_M(:,0,:)

      IF (DES_CONTINUUM_HYBRID) THEN
! modify center coefficient to account for contribution from discrete phases.
! recall the center coefficient is negative.              
         DO IJK = ijkstart3, ijkend3              
            DO DM=1,DES_MMAX
               AM0(IJK,0) = AM0(IJK,0) - VXF_GDS(IJK,DM)
            ENDDO
            DO M = 1,MMAX
               DO DM=1,DES_MMAX
                  AM0(IJK,M) = AM0(IJK,M) - VXF_SDS(IJK,M,DM)
               ENDDO
            ENDDO
         ENDDO
      ELSE
! Loop over the the number of DES phases if this the is coupled
! DEM or PIC simulation. Otherwise, use the number of TFM solids.
         LC_MMAX = merge(DES_MMAX, MMAX, DES_CONTINUUM_COUPLED)
      ENDIF

! Determine which calculations are needed. This is a first pass filter
! that works well for1 solids phase. 
      DO M = 1, MMAX
         IF (MOMENTUM_X_EQ(0) .AND. MOMENTUM_X_EQ(M)) THEN
! there is at least one solids phase x-momentum equation and the gas phase
! x-momentum equation                  
            Pass1 = .TRUE.
            GOTO 10
         ELSEIF (MOMENTUM_X_EQ(M)) THEN
! there is at least one solids phase x-momentum equation but the gas phase
! x-momentum equation is NOT solved
            Pass2 = .TRUE.
            GOTO 10
         ENDIF
      ENDDO

   10 CONTINUE

 IF (Pass1) THEN
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


 ELSEIF (MOMENTUM_X_EQ(0)) THEN
! the solids x-momentum equations are not solved; only gas phase (m=0)
! x-momentum equation is solved.  this branch is where a coupled DEM
! simulation should be directed for proper evaluation of pressure
! correction terms
! ---------------------------------------------------------------->>>
!!$omp   parallel do private(IJK, I, IJKE, EPGA, M, SUM_VXF_GS,&
!!$omp&  AREA_FACE), &
!!$omp&  schedule(static)
   DO IJK = ijkstart3, ijkend3
     IF (IP_AT_E(IJK) .OR. MFLOW_AT_E(IJK)) THEN   !impermeable
       D_E(IJK,0) = ZERO
     ELSE

       AREA_FACE = merge(ONE, AYZ(IJK), CARTESIAN_GRID)
       I = I_OF(IJK)
       IJKE = EAST_OF(IJK)
       EPGA = AVG_X(EP_G(IJK),EP_G(IJKE),I)

       SUM_VXF_GS = ZERO
       IF (.NOT. QMOMK) THEN
          DO M= 1, LC_MMAX
! Gas - All Solids VolxDrag summation
            SUM_VXF_GS = SUM_VXF_GS + VXF_GS(IJK,M)
          ENDDO
       ELSE
          DO INN = 1, QMOMK_NN
             DO M = 1, MMAX
               SUM_VXF_GS = SUM_VXF_GS + &
                  AVG_X(QMOMK_F_GS(INN,IJK,M),QMOMK_F_GS(INN,IJKE,M),I)*&
                  VOL_U(IJK)
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

     ENDIF   ! end if/else branch(ip_at_e(ijk) .or. mflow_at_e(ijk))
   ENDDO   ! end do (ijk=ijkstart3,ijkend3)
! end only momentum_x_eq(0) is solved branch
! ----------------------------------------------------------------<<<


 ELSEIF (Pass2) THEN
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

 ENDIF   ! end if/else branch (pass1,momentum_x_eq(0),pass2)

 RETURN
 END SUBROUTINE CALC_D_E


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_d_n                                                C
!  Purpose: calculate coefficients linking velocity correction to      C
!           pressure correction -- North                               C
!                                                                      C
!  Notes: MFIX convention: center coeff is negative, hence             C
!                            (-A_M(IJK,0,M)) > or = 0                  C
!         The NORTH face area is AXZ                                   C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-JUN-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number: 2                                                  C
!  Purpose: allow multiparticle in D_N calculation in accounting for   C
!     the averaged Solid-Solid drag                                    C
!  Author: S. Dartevelle, LANL                        Date: 28-FEB-04  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number: 3                                                  C
!  Purpose: more accurate formulas in the D_s(M) for Model_A           C
!  Author: S. Dartevelle, LANL                        Date: 17-MAR-04  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number: 4                                                  C
!  Purpose: To incorporate Cartesian grid modifications                C
!  Author: Jeff Dietiker                              Date: 01-Jul-09  C
!                                                                      C
!  Revision Number: 5                                                  C
!  Purpose: Incorporation of QMOM for the solution of the particle     C
!           kinetic equation                                           C
!  Author: Alberto Passalacqua - Fox Research Group   Date: 02-Dec-09  C
!                                                                      C
!  Literature/Document References:                                     C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CALC_D_N(A_M, VXF_GS, VXF_SS, D_N, IER) 

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
! initialized to zero in the subroutine time_march before the time loop 
      DOUBLE PRECISION, INTENT(INOUT) :: d_n(DIMENSION_3, 0:DIMENSION_M)
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
      INTEGER :: J, IJK, IJKN
! Number of solids phases to use as loop counter      
      INTEGER :: LC_MMAX
! Average volume fraction at momentum cell centers
      DOUBLE PRECISION :: EPGA 
! Average solid volume fraction at momentum cell centers
      DOUBLE PRECISION :: EPSA(DIMENSION_M) 
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
! Temporary variable to store A matrix for local manipulation
      DOUBLE PRECISION :: AM0(DIMENSION_3, 0:DIMENSION_M)
!-----------------------------------------------
! Include statement functions
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'
!-----------------------------------------------

! initializing
      Pass1 = .FALSE.
      Pass2 = .FALSE.
      AM0(:,:) = A_M(:,0,:)

      IF (DES_CONTINUUM_HYBRID) THEN
         DO IJK = ijkstart3, ijkend3              
            DO DM=1,DES_MMAX
               AM0(IJK,0) = AM0(IJK,0) - VXF_GDS(IJK,DM)
            ENDDO
            DO M = 1,MMAX
               DO DM=1,DES_MMAX
                  AM0(IJK,M) = AM0(IJK,M) - VXF_SDS(IJK,M,DM)
               ENDDO
            ENDDO
         ENDDO
      ELSE
! Loop over the the number of DES phases if this the is coupled
! DEM or PIC simulation. Otherwise, use the number of TFM solids.
         LC_MMAX = merge(DES_MMAX, MMAX, DES_CONTINUUM_COUPLED)
      ENDIF


! Determine which calculations are needed      
      DO M = 1, MMAX
         IF (MOMENTUM_Y_EQ(0) .AND. MOMENTUM_Y_EQ(M)) THEN
! there is at least one solids phase y-momentum equation and the gas phase
! y-momentum equation                  
            Pass1 = .TRUE.
            GOTO 10
         ELSEIF (MOMENTUM_Y_EQ(M)) THEN
! there is at least one solids phase y-momentum equation but the gas phase
! y-momentum equation is NOT solved
            Pass2 = .TRUE.
            GOTO 10
         ENDIF
      ENDDO

   10 CONTINUE

 IF (Pass1) THEN    
! Gas (m=0) and at least one solids phase (m>0) y-momentum equation are
! solved
! ---------------------------------------------------------------->>>
!!$omp   parallel do private(J,IJK, IJKN, EPGA, EPSA, AREA_FACE, &
!!$omp&  M, L, Lp, LpL, LM, SUM_VXF_SS, SUM_VXF_SS_wt_M, SUM_VXF_GS, &
!!$omp&  DEN_MGas, NUM_MGas, DEN_MSol_LGas, DEN_MSol_LSol ,&
!!$omp&  NUM_MSol_LSol, NUM_MSol_LGas),&
!!$omp&  schedule(static)
   DO IJK = ijkstart3, ijkend3
     IF (IP_AT_N(IJK) .OR. MFLOW_AT_N(IJK)) THEN   !impermeable
        DO M= 0, MMAX
           D_N(IJK,M) = ZERO
        ENDDO
     ELSE

        AREA_FACE = merge(ONE, AXZ(IJK), CARTESIAN_GRID)
        J = J_OF(IJK)
        IJKN = NORTH_OF(IJK)
        EPGA = AVG_Y(EP_G(IJK),EP_G(IJKN),J)

        SUM_VXF_GS = ZERO
        DO M= 1, MMAX
          EPSA(M) = AVG_Y(EP_S(IJK,M),EP_S(IJKN,M),J)
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
        NUM_MGas = ZERO        
        DEN_MGas  = ZERO
        DO M= 1, MMAX
          IF(MOMENTUM_Y_EQ(M)) THEN
             NUM_MGas = NUM_MGas + (EPSA(M)*VXF_GS(IJK,M)/&
                ((-AM0(IJK,M))+VXF_GS(IJK,M)+SUM_VXF_SS(M)+&
                SMALL_NUMBER) )
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
            D_N(IJK,0) = P_SCALE*AREA_FACE/((-AM0(IJK,0))+DEN_MGas)
          ELSE
            D_N(IJK,0) = ZERO
          ENDIF
! Linking velocity correction coefficient to pressure - SOLIDs Phase
          DO M = 1, MMAX
            IF(MOMENTUM_Y_EQ(M) .AND. ((-AM0(IJK,M))>SMALL_NUMBER .OR. &
                                       VXF_GS(IJK,M)>SMALL_NUMBER)) THEN
              D_N(IJK,M) = D_N(IJK,0)*(VXF_GS(IJK,M)/&
                           ((-AM0(IJK,M))+VXF_GS(IJK,M)))
            ELSE
              D_N(IJK,M) = ZERO
            ENDIF
          ENDDO

! Model A
! -------------------------------------------------------------------
        ELSE
 
! Linking velocity correction coefficient to pressure - GAS Phase
          IF ( (-AM0(IJK,0))>SMALL_NUMBER .OR. &
                    DEN_MGas>SMALL_NUMBER ) THEN
            D_N(IJK,0) = P_SCALE*AREA_FACE*(EPGA+NUM_MGas)/&
                         ((-AM0(IJK,0))+DEN_MGas ) 
          ELSE
            D_N(IJK,0) = ZERO
          ENDIF

          DO M= 1, MMAX
! calculating NUM_MSol_LGas and DEN_MSol_LGas
            NUM_MSol_LGas = VXF_GS(IJK,M)*EPGA/&
              ((-AM0(IJK,0))+SUM_VXF_GS+SMALL_NUMBER)
            DEN_MSol_LGas = VXF_GS(IJK,M)*( &
              ((-AM0(IJK,0))+(SUM_VXF_GS - VXF_GS(IJK,M)))/&
              ((-AM0(IJK,0))+SUM_VXF_GS+SMALL_NUMBER) )
      
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

                IF(MOMENTUM_Y_EQ(L)) THEN
                  NUM_MSol_LSol(M) = NUM_MSol_LSol(M) + &
                     (VXF_SS(IJK,LM)*EPSA(L)/&
                     ((-AM0(IJK,L))+VXF_GS(IJK,L)+SUM_VXF_SS(L)+&
                     SMALL_NUMBER )) 
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
            IF (MOMENTUM_Y_EQ(M) .AND. (-AM0(IJK,M)>SMALL_NUMBER .OR. &
                                      DEN_MSol_LGas>SMALL_NUMBER .OR. &
                                   DEN_MSol_LSol(M)>SMALL_NUMBER)) THEN
              D_N(IJK,M) = P_SCALE*AREA_FACE*&
                (EPSA(M) + NUM_MSol_LSol(M) + NUM_MSol_LGas)/&
                ((-AM0(IJK,M))+DEN_MSol_LGas+DEN_MSol_LSol(M))
            ELSE
              D_N(IJK,M) = ZERO
            ENDIF
          ENDDO   ! end do (m=1,mmax)

        ENDIF    !end if/else branch Model_B/Model_A 

     ENDIF   ! end if/else branch(ip_at_n(ijk) .or. mflow_at_n(ijk))
   ENDDO  ! end do loop (ijk=ijkstart3,ijkend3)
! end pass1 branch
! ----------------------------------------------------------------<<<


 ELSEIF (MOMENTUM_Y_EQ(0)) THEN 
! the solids y-momentum equations are not solved; only gas phase (m=0)
! y-momentum equation is solved. this branch is where a coupled DEM 
! simulation should be directed for proper evaluation of pressure 
! correction terms
! ---------------------------------------------------------------->>>
!!$omp   parallel do private(IJK, J, IJKN, EPGA, M, SUM_VXF_GS,&
!!$omp&  AREA_FACE), &
!!$omp&  schedule(static)
   DO IJK = ijkstart3, ijkend3
     IF (IP_AT_N(IJK) .OR. MFLOW_AT_N(IJK)) THEN   !impermeable
       D_N(IJK,0) = ZERO
     ELSE

       AREA_FACE = merge(ONE, AXZ(IJK), CARTESIAN_GRID)
       J = J_OF(IJK)
       IJKN = NORTH_OF(IJK)
       EPGA = AVG_Y(EP_G(IJK),EP_G(IJKN),J)

       SUM_VXF_GS = ZERO
       IF (.NOT. QMOMK) THEN
          DO M= 1, LC_MMAX
! Gas - All Solids VolxDrag summation
            SUM_VXF_GS = SUM_VXF_GS + VXF_GS(IJK,M)
          ENDDO
       ELSE
          DO INN = 1, QMOMK_NN
             DO M = 1, MMAX
               SUM_VXF_GS = SUM_VXF_GS + &
                  AVG_Y(QMOMK_F_GS(INN,IJK,M),QMOMK_F_GS(INN,IJKN,M),J)*&
                  VOL_V(IJK)
             ENDDO
           ENDDO
       ENDIF

       IF ((-AM0(IJK,0))>SMALL_NUMBER .OR. SUM_VXF_GS>SMALL_NUMBER) THEN
         IF (MODEL_B) THEN
           D_N(IJK,0) = P_SCALE*AREA_FACE/((-AM0(IJK,0))+SUM_VXF_GS)
         ELSE
           D_N(IJK,0) = P_SCALE*AREA_FACE*EPGA/((-AM0(IJK,0))+SUM_VXF_GS)
         ENDIF   !end if/else branch Model_B/Model_A 
       ELSE
         D_N(IJK,0) = ZERO
       ENDIF

     ENDIF   ! end if/else branch(ip_at_n(ijk) .or. mflow_at_n(ijk))
   ENDDO   ! end do (ijk=ijkstart3,ijkend3)
! end only gas y-momentum equation is solved branch
! ----------------------------------------------------------------<<<


 ELSEIF (Pass2) THEN
! at least one solids phase y-momentum equation is solved but the gas 
! phase (m=0) y-momentum equation is not solved
! ---------------------------------------------------------------->>>
!!$omp    parallel do private(IJK, J, IJKN, EPSA, L, Lp, M, SUM_VXF_SS, &
!!$omp&   NUM_MSol_LSol, DEN_MSol_LSol, SUM_VXF_SS_wt_M, AREA_FACE), &
!!$omp&  schedule(static)
   DO IJK = ijkstart3, ijkend3
     IF (IP_AT_N(IJK) .OR. MFLOW_AT_N(IJK) .OR. MODEL_B) THEN
       DO M= 1, MMAX
         D_N(IJK,M) = ZERO
       ENDDO
     ELSE

       AREA_FACE = merge(ONE, AXZ(IJK), CARTESIAN_GRID)
       J = J_OF(IJK)
       IJKN = NORTH_OF(IJK)

       DO M= 1, MMAX
         EPSA(M) = AVG_Y(EP_S(IJK,M),EP_S(IJKN,M),J)     
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
             IF(MOMENTUM_Y_EQ(L)) THEN             
               NUM_MSol_LSol(M) = NUM_MSol_LSol(M) + &
                  (VXF_SS(IJK,LM)*EPSA(L)/&
                  ( (-AM0(IJK,L))+VXF_GS(IJK,L)+SUM_VXF_SS(L)+&
                  SMALL_NUMBER ))
               DEN_MSol_LSol(M) = DEN_MSol_LSol(M) + VXF_SS(IJK,LM)*(&
                  ( (-AM0(IJK,L))+VXF_GS(IJK,L)+SUM_VXF_SS_wt_M )/&
                  ( (-AM0(IJK,L))+VXF_GS(IJK,L)+SUM_VXF_SS(L)+&
                  SMALL_NUMBER ))
             ELSE
               DEN_MSol_LSol(M) = DEN_MSol_LSol(M) + VXF_SS(IJK,LM)
             ENDIF
           ENDIF   ! end if (l.ne.m)
         ENDDO  ! end do (l=1,mmax)
         
! Linking velocity correction coefficient to pressure - SOLIDs Phase (Model_A only)         
         IF (MOMENTUM_Y_EQ(M) .AND. (-AM0(IJK,M)>SMALL_NUMBER .OR. &
                                   VXF_GS(IJK,M)>SMALL_NUMBER .OR. &
                                DEN_MSol_LSol(M)>SMALL_NUMBER)) THEN
           D_N(IJK,M) = P_SCALE*AREA_FACE*&
              (EPSA(M) + NUM_MSol_LSol(M))/&
              ((-AM0(IJK,M))+VXF_GS(IJK,M)+DEN_MSol_LSol(M))
         ELSE
           D_N(IJK,M) = ZERO
         ENDIF
       ENDDO  ! end do (m=1,mmax)

     ENDIF   ! end if/else branch(ip_at_n(ijk) .or. mflow_at_n(ijk) .or. model_b)
   ENDDO   ! end do (ijk=ijkstart3,ijkend3)
! end pass2 branch
! ----------------------------------------------------------------<<<

 ENDIF   ! end if/else branch (pass1,momentum_y_eq(0),pass2)

 RETURN
 END SUBROUTINE CALC_D_N


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_d_t                                                C
!  Purpose: Calculate coefficients linking velocity correction to      C
!           pressure correction -- Top                                 C
!                                                                      C
!  Notes: MFIX convention: center coeff is negative, hence             C
!                            (-A_M(IJK,0,M)) > or = 0                  C
!         The TOP face area is AXY                                     C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-JUN-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number: 2                                                  C
!  Purpose: allow multiparticle in D_T calculation in accounting for   C
!           the averaged Solid-Solid drag                              C
!  Author: S. Dartevelle, LANL                        Date: 28-FEB-04  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number: 3                                                  C
!  Purpose: more accurate formulas in the D_s(M) for Model_A           C
!  Author: S. Dartevelle, LANL                        Date: 17-MAR-04  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number: 4                                                  C
!  Purpose: To incorporate Cartesian grid modifications                C
!  Author: Jeff Dietiker                              Date: 01-Jul-09  C
!                                                                      C
!  Revision Number: 5                                                  C
!  Purpose: Incorporation of QMOM for the solution of the particle     C
!           kinetic equation                                           C
!  Author: Alberto Passalacqua - Fox Research Group   Date: 02-Dec-09  C
!                                                                      C
!  Literature/Document References:                                     C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
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
      INTEGER :: IJK, IJKT, K
! Number of solids to use as loop counter
      INTEGER :: LC_MMAX
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
! Temporary variable to store A matrix for local manipulation
      DOUBLE PRECISION :: AM0(DIMENSION_3, 0:DIMENSION_M)
!-----------------------------------------------
! Include statement functions
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'
!-----------------------------------------------

! initializing
      Pass1 = .FALSE.
      Pass2 = .FALSE.
      AM0(:,:) = A_M(:,0,:)

      IF (DES_CONTINUUM_HYBRID) THEN
         DO IJK = ijkstart3, ijkend3              
            DO DM=1,DES_MMAX
               AM0(IJK,0) = AM0(IJK,0) - VXF_GDS(IJK,DM)
            ENDDO
            DO M = 1,MMAX
               DO DM=1,DES_MMAX
                  AM0(IJK,M) = AM0(IJK,M) - VXF_SDS(IJK,M,DM)
               ENDDO
            ENDDO
         ENDDO
      ELSE
! Loop over the the number of DES phases if this the is coupled
! DEM or PIC simulation. Otherwise, use the number of TFM solids.
         LC_MMAX = merge(DES_MMAX, MMAX, DES_CONTINUUM_COUPLED)
      ENDIF


! Determine which calculations are needed      
      DO M = 1, MMAX
         IF (MOMENTUM_Z_EQ(0) .AND. MOMENTUM_Z_EQ(M)) THEN
! there is at least one solids phase z-momentum equation and the gas phase
! z-momentum equation
            Pass1 = .TRUE.
            GOTO 10
         ELSEIF (MOMENTUM_Z_EQ(M)) THEN
! there is at least one solids phase z-momentum equation but the gas phase
! z-momentum equation is NOT solved
            Pass2 = .TRUE.
            GOTO 10
         ENDIF
      ENDDO

   10 CONTINUE

 IF (Pass1) THEN
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


 ELSEIF (MOMENTUM_Z_EQ(0)) THEN
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
          DO M= 1, LC_MMAX
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


 ELSEIF (Pass2) THEN
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

 ENDIF  ! end if/else branch (pass1,momentum_z_eq(0),pass2)


 RETURN
 END SUBROUTINE CALC_D_T


