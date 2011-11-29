!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_d_e(A_m, VxF_gs, VxF_ss, d_e, IER)                C
!  Purpose: calculte coefficients linking velocity correction to       C
!           pressure correction -- East                                C
!                                                                      C
!  Note:  MFIX convention: center coeff is negative, hence:            C
!                            (-A_M(IJK,0,M)) > or = 0                  C
!         The EAST face area is AYZ                                    C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-JUN-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number: 2                                                  C
!  Purpose: allow multiparticle in D_E calculation in accounting for   C
!           the averaged Solid-Solid drag and multi-solid-particle dragC
!                                                                      C
!  Author: S. Dartevelle, LANL                        Date: 28-FEB-04  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number: 3                                                  C
!  Purpose: more accurate formulas in the D_s(M) for Model_A           C
!                                                                      C
!  Author: S. Dartevelle, LANL                        Date: 17-MAR-04  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number: 4                                                  C
!  Purpose: To incorporate Cartesian grid modifications                C
!  Author: Jeff Dietiker                              Date: 01-Jul-09  C
!                                                                      C
!  Revision Number: 5                                                  C
!  Purpose: Incorporation of QMOM for the solution of the particle     C
!  kinetic equation                                                    C
!  Author: Alberto Passalacqua - Fox Research Group   Date: 02-Dec-09  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CALC_D_E(A_M, VXF_GS, VXF_SS, D_E, IER)
!
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
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
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
      USE cutcell
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================

! QMOMK - Alberto Passalacqua
      USE qmom_kinetic_equation
! QMOMK - End

      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!
!                      Error index
      INTEGER          IER
!
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
!
!                      Volume x average at momentum cell centers
      DOUBLE PRECISION VxF_gs(DIMENSION_3, DIMENSION_M)
!                      
      DOUBLE PRECISION d_e(DIMENSION_3, 0:DIMENSION_M)
!
!                      Average volume fraction at momentum cell centers
      DOUBLE PRECISION EPGA                                        !S. Dartevelle, LANL, Feb.2004
!          Volume x average at momentum cell centers Solid-Solid Drag
      DOUBLE PRECISION VxF_ss(DIMENSION_3, DIMENSION_LM)           !S. Dartevelle, LANL, Feb.2004
!          Usual Indices
      INTEGER           LM, M, L, LpL, Lp, I, IJK, IJKE, IN        !S. Dartevelle, LANL, Feb.2004
!          Average solid volume fraction at momentum cell centers
      DOUBLE PRECISION  EPSA(DIMENSION_M)                          !S. Dartevelle, LANL, Feb.2004
!          ratio of drag and A0 or A_solid and other sum of Solid-Solid drag
      DOUBLE PRECISION  other_ratio_1, FOA1                        !S. Dartevelle, LANL, Feb.2004
!          sum of All Solid-Gas VolxDrag
      DOUBLE PRECISION  SUM_VXF_GS                                 !S. Dartevelle, LANL, Feb.2004
!          What phase momentum equation is activated?
      LOGICAL  Pass1, Pass2                                        !S. Dartevelle, LANL, Feb.2004
!          sum of Solid M - All other Solid drag
      DOUBLE PRECISION  SUM_VXF_SS(DIMENSION_M)                    !S. Dartevelle, LANL, Feb.2004
!          numerator needed for solid phase M time EPSA
      DOUBLE PRECISION  numeratorxEP(DIMENSION_M)                  !S. Dartevelle, LANL, Feb.2004
!          denominator needed for solid phase M
      DOUBLE PRECISION  denominator(DIMENSION_M)                   !S. Dartevelle, LANL, Feb.2004
!          denominator needed for solid phase M
      DOUBLE PRECISION  other_denominator(DIMENSION_M)             !S. Dartevelle, LANL, Feb.2004
!          sum of Solid L - All other Solid VolxDrag but M
      DOUBLE PRECISION  SUM_VXF_SS_wt_M                            !S. Dartevelle, LANL, Feb.2004
!          total solids volume fraction
      DOUBLE PRECISION  EPStmp
!-----------------------------------------------
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'


   Pass1 = .FALSE.     !initialization
   Pass2 = .FALSE.
    DO M = 1, MMAX
         if (MOMENTUM_X_EQ(0) .AND. MOMENTUM_X_EQ(M)) then
           Pass1 = .TRUE.    !we have at least one solid phase X-momentum equation
           GO TO 10          !with the gas phase X-momentum equation
         elseif (MOMENTUM_X_EQ(M)) then
           Pass2 = .TRUE.    !we have at least one solid phase X-momentum equation
           GO TO 10          !but the gas phase X-momentum is not solved
         endif
    END DO

   10 CONTINUE


 IF (Pass1) THEN    !Gas and at least one solid phase X-momentum equation
!$omp   parallel do private(I,IJK, IJKE, EPGA, EPSA,EPStmp, numeratorxEP, &
!$omp&  M, L, Lp, LpL, LM, other_ratio_1, FOA1, SUM_VXF_SS, &
!$omp&  SUM_VXF_SS_wt_M, SUM_VXF_GS, other_denominator, denominator ),&
!$omp&  schedule(static)
  DO IJK = ijkstart3, ijkend3
     IF (IP_AT_E(IJK) .OR. MFLOW_AT_E(IJK)) THEN   !impermeable
        DO M= 0, MMAX
         D_E(IJK,M) = ZERO
        END DO
     ELSE
        I = I_OF(IJK)
        IJKE = EAST_OF(IJK)
        EPGA = AVG_X(EP_G(IJK),EP_G(IJKE),I)

        SUM_VXF_GS = ZERO
        DO M= 1, MMAX
          EPSA(M) = AVG_X(EP_S(IJK,M),EP_S(IJKE,M),I)
          SUM_VXF_GS = SUM_VXF_GS + VXF_GS(IJK,M)              !Gas - All Solids VolxDrag summation
          SUM_VXF_SS(M) = ZERO
          do L = 1, MMAX
            IF (L .NE. M) THEN
              LM = FUNLM(L,M)
              SUM_VXF_SS(M) = SUM_VXF_SS(M) + VXF_SS(IJK,LM) !Solid M - All other Solids VolxDrag summation
            ENDIF
          end do
        END DO

        other_ratio_1  = ZERO
        DO M= 1, MMAX
           other_ratio_1 = other_ratio_1 +&
                      ( VXF_GS(IJK,M)*( (-A_M(IJK,0,M))+SUM_VXF_SS(M) )/&
                                      ( (-A_M(IJK,0,M))+VXF_GS(IJK,M)+SUM_VXF_SS(M)+SMALL_NUMBER )&
                      )
        END DO

        IF (MODEL_B) THEN   !Model B
          !Linking velocity correction coefficient to pressure - GAS Phase
          if ( ((-A_M(IJK,0,0))>SMALL_NUMBER) .OR.  (other_ratio_1>SMALL_NUMBER) ) then
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            IF(.NOT.CARTESIAN_GRID) THEN
               D_E(IJK,0) = P_SCALE*AYZ(IJK)/( (-A_M(IJK,0,0))+other_ratio_1 )
            ELSE
               D_E(IJK,0) = P_SCALE         /( (-A_M(IJK,0,0))+other_ratio_1 )
            ENDIF
! Original terms
!             D_E(IJK,0) = P_SCALE*AYZ(IJK)/( (-A_M(IJK,0,0))+other_ratio_1 )
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
          else
             D_E(IJK,0) = ZERO
          endif
          !Linking velocity correction coefficient to pressure - SOLID Phase
          DO M = 1, MMAX
            if ( MOMENTUM_X_EQ(M) .AND. (  ((-A_M(IJK,0,M))>SMALL_NUMBER) .OR. &
                                                     (VXF_GS(IJK,M)>SMALL_NUMBER) ) )   then
               D_E(IJK,M) = D_E(IJK,0)*(&
                                VXF_GS(IJK,M)/((-A_M(IJK,0,M))+VXF_GS(IJK,M))&
                            )
            else
               D_E(IJK,M) = ZERO
            endif
          END DO
        ELSE                !Model A
          FOA1 = ZERO
          DO M= 1, MMAX
            FOA1 = FOA1 + (EPSA(M)*VXF_GS(IJK,M)/&
                             ( (-A_M(IJK,0,M))+VXF_GS(IJK,M)+SUM_VXF_SS(M)+SMALL_NUMBER )&
                          )
            other_denominator(M) = VXF_GS(IJK,M)*( ((-A_M(IJK,0,0))+(SUM_VXF_GS - VXF_GS(IJK,M)))/&
                                                   ((-A_M(IJK,0,0))+SUM_VXF_GS+SMALL_NUMBER)&
                                                 )
            numeratorxEP(M) = ZERO
            denominator(M)  = ZERO
            do L = 1, MMAX
              IF (L .NE. M) THEN
                LM = FUNLM(L,M)
                numeratorxEP(M) = numeratorxEP(M) + (&
                                                      VXF_SS(IJK,LM)*EPSA(L)/&
                                                      ( (-A_M(IJK,0,L))+VXF_GS(IJK,L)+SUM_VXF_SS(L)+SMALL_NUMBER )&
                                                    ) 
                SUM_VXF_SS_wt_M = ZERO
                do Lp = 1, MMAX
                  if ( (Lp .NE. L) .AND. (Lp .NE. M) ) then
                    LpL = FUNLM(Lp,L)
                    SUM_VXF_SS_wt_M = SUM_VXF_SS_wt_M + VXF_SS(IJK,LpL)   !Solid L - All other Solids VolxDrag but M summation
                  endif
                end do
                denominator(M) = denominator(M) + (&
                                                     VXF_SS(IJK,LM)*(&
                                                       ( (-A_M(IJK,0,L))+VXF_GS(IJK,L)+SUM_VXF_SS_wt_M )/&
                                                       ( (-A_M(IJK,0,L))+VXF_GS(IJK,L)+SUM_VXF_SS(L)+SMALL_NUMBER )&
                                                     )&
                                                  )
              ENDIF
            end do
          END DO
          !Linking velocity correction coefficient to pressure - GAS Phase
          if ( ((-A_M(IJK,0,0))>SMALL_NUMBER) .OR.  (other_ratio_1>SMALL_NUMBER) ) then
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            IF(.NOT.CARTESIAN_GRID) THEN
               D_E(IJK,0) = P_SCALE*AYZ(IJK)*(EPGA+FOA1)/( (-A_M(IJK,0,0))+other_ratio_1 )
            ELSE
               D_E(IJK,0) = P_SCALE         *(EPGA+FOA1)/( (-A_M(IJK,0,0))+other_ratio_1 )
            ENDIF
! Original terms
!            D_E(IJK,0) = P_SCALE*AYZ(IJK)*(EPGA+FOA1)/( (-A_M(IJK,0,0))+other_ratio_1 )
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
          else
            D_E(IJK,0) = ZERO
          endif
          !Linking velocity correction coefficient to pressure - SOLID Phase
          DO M = 1, MMAX
            if ( MOMENTUM_X_EQ(M) .AND. ( (-A_M(IJK,0,M)>SMALL_NUMBER)        .OR. &
                                                (other_denominator(M)>SMALL_NUMBER) .OR. &
                                                (denominator(M)>SMALL_NUMBER) ) )     then
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            IF(.NOT.CARTESIAN_GRID) THEN
               D_E(IJK,M) = P_SCALE*AYZ(IJK)*(&
                         ( EPSA(M) + numeratorxEP(M) + (VXF_GS(IJK,M)*EPGA/((-A_M(IJK,0,0))+SUM_VXF_GS+SMALL_NUMBER)) )/&
                         ( (-A_M(IJK,0,M))+other_denominator(M)+denominator(M) )&
                                                    )
            ELSE
               D_E(IJK,M) = P_SCALE           *(&
                         ( EPSA(M) + numeratorxEP(M) + (VXF_GS(IJK,M)*EPGA/((-A_M(IJK,0,0))+SUM_VXF_GS+SMALL_NUMBER)) )/&
                         ( (-A_M(IJK,0,M))+other_denominator(M)+denominator(M) )&
                                                    )

            ENDIF
! Original terms
!              D_E(IJK,M) = P_SCALE*AYZ(IJK)*(&
!                        ( EPSA(M) + numeratorxEP(M) + (VXF_GS(IJK,M)*EPGA/((-A_M(IJK,0,0))+SUM_VXF_GS+SMALL_NUMBER)) )/&
!                        ( (-A_M(IJK,0,M))+other_denominator(M)+denominator(M) )&
!                                                   )
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            else
              D_E(IJK,M) = ZERO
            endif
          END DO
        ENDIF    !end of Model_B/Model_A if then condition
     ENDIF
  ENDDO

 ELSE IF (MOMENTUM_X_EQ(0)) THEN    !the solid X-momentum equations are not solved
                                    !only gas phase X-momentum
!$omp   parallel do private(IJK, I, IJKE, EPGA, M, SUM_VXF_GS), &
!$omp&  schedule(static)
   DO IJK = ijkstart3, ijkend3
     IF (IP_AT_E(IJK) .OR. MFLOW_AT_E(IJK)) THEN   !impermeable
       D_E(IJK,0) = ZERO
     ELSE
       I = I_OF(IJK)
       IJKE = EAST_OF(IJK)
       EPGA = AVG_X(EP_G(IJK),EP_G(IJKE),I)

       ! QMOMK - Alberto Passalacqua 
       ! Added check for QMOMK
       SUM_VXF_GS = ZERO
       IF (.NOT. QMOMK) THEN
        DO M = 1, MMAX
         SUM_VXF_GS = SUM_VXF_GS + VXF_GS(IJK,M)              !Gas - All Solids VolxDrag summation
        END DO
       ELSE
        DO IN = 1, QMOMK_NN
          DO M = 1, MMAX           
            SUM_VXF_GS = SUM_VXF_GS + AVG_X(QMOMK_F_GS(IN,IJK,M),QMOMK_F_GS(IN,IJKE,M),I)*VOL_U(IJK)
          END DO
        END DO
       END IF
       ! QMOMK - End

       IF ( ((-A_M(IJK,0,0))>SMALL_NUMBER) .OR. (SUM_VXF_GS>SMALL_NUMBER) ) THEN
         IF (MODEL_B) THEN
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            IF(.NOT.CARTESIAN_GRID) THEN
               D_E(IJK,0) = P_SCALE*AYZ(IJK)/((-A_M(IJK,0,0))+SUM_VXF_GS)
            ELSE
               D_E(IJK,0) = P_SCALE         /((-A_M(IJK,0,0))+SUM_VXF_GS)
            ENDIF
! Original terms
!           D_E(IJK,0) = P_SCALE*AYZ(IJK)/((-A_M(IJK,0,0))+SUM_VXF_GS)
         ELSE

            IF(.NOT.CARTESIAN_GRID) THEN
               D_E(IJK,0) = P_SCALE*AYZ(IJK)*EPGA/((-A_M(IJK,0,0))+SUM_VXF_GS)
            ELSE
               D_E(IJK,0) = P_SCALE         *EPGA/((-A_M(IJK,0,0))+SUM_VXF_GS)
            ENDIF
! Original terms
!          D_E(IJK,0) = P_SCALE*AYZ(IJK)*EPGA/((-A_M(IJK,0,0))+SUM_VXF_GS)
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
         ENDIF
       ELSE
         D_E(IJK,0) = ZERO
       ENDIF
     ENDIF
   END DO

 ELSE IF (Pass2) THEN    !at least one solid phase X-momentum equation is solved
                         !but the gas phase X-momentum is not solved
!$omp    parallel do private(IJK, I, IJKE, EPSA,EPStmp, L, Lp, M, SUM_VXF_SS, &
!$omp&   numeratorxEP, denominator, SUM_VXF_SS_wt_M), &
!$omp&  schedule(static)
   DO IJK = ijkstart3, ijkend3
     IF (IP_AT_E(IJK) .OR. MFLOW_AT_E(IJK) .OR. MODEL_B) THEN
       DO M= 1, MMAX
         D_E(IJK,M) = ZERO
       END DO
     ELSE
       I = I_OF(IJK)
       IJKE = EAST_OF(IJK)

       DO M= 1, MMAX
         EPSA(M) = AVG_X(EP_S(IJK,M),EP_S(IJKE,M),I)
         SUM_VXF_SS(M) = ZERO
         do L = 1, MMAX
           IF (L .NE. M) THEN
             LM = FUNLM(L,M)
             SUM_VXF_SS(M) = SUM_VXF_SS(M) + VXF_SS(IJK,LM)  !Solid M - All other Solids VolxDrag summation
           ENDIF
         end do
       END DO
       DO M= 1, MMAX
         numeratorxEP(M) = ZERO
         denominator(M)  = ZERO
         do L = 1, MMAX
           IF (L .NE. M) THEN
             LM = FUNLM(L,M)
             numeratorxEP(M) = numeratorxEP(M) + (&
                                                    VXF_SS(IJK,LM)*EPSA(L)/&
                                                   ( (-A_M(IJK,0,L))+VXF_GS(IJK,L)+SUM_VXF_SS(L)+SMALL_NUMBER )&
                                                  )
             SUM_VXF_SS_wt_M = ZERO
             do Lp = 1, MMAX
               if ( (Lp .NE. L) .AND. (Lp .NE. M) ) then
                 LpL = FUNLM(Lp,L)
                 SUM_VXF_SS_wt_M = SUM_VXF_SS_wt_M + VXF_SS(IJK,LpL)     !Solid L - All other Solids VolxDrag but M summation
               endif
             end do
             denominator(M) = denominator(M) + (&
                                                  VXF_SS(IJK,LM)*(&
                                                   ( (-A_M(IJK,0,L))+VXF_GS(IJK,L)+SUM_VXF_SS_wt_M )/&
                                                   ( (-A_M(IJK,0,L))+VXF_GS(IJK,L)+SUM_VXF_SS(L)+SMALL_NUMBER )&
                                                                 )&
                                                )
           ENDIF
         end do
       END DO
!Linking velocity correction coefficient to pressure - SOLID Phase (Model_A only)
       DO M = 1, MMAX
         if ( MOMENTUM_X_EQ(M) .AND. ( (-A_M(IJK,0,M)>SMALL_NUMBER) .OR. &
                                        (VXF_GS(IJK,M)>SMALL_NUMBER) .OR. &
                                        (denominator(M)>SMALL_NUMBER) ) ) then
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            IF(.NOT.CARTESIAN_GRID) THEN
               D_E(IJK,M) = P_SCALE*AYZ(IJK)*( EPSA(M) + numeratorxEP(M) )/&
                                             ( (-A_M(IJK,0,M))+VXF_GS(IJK,M)+denominator(M) )
            ELSE
               D_E(IJK,M) = P_SCALE         *( EPSA(M) + numeratorxEP(M) )/&
                                             ( (-A_M(IJK,0,M))+VXF_GS(IJK,M)+denominator(M) )
            ENDIF
! Original terms
!           D_E(IJK,M) = P_SCALE*AYZ(IJK)*( EPSA(M) + numeratorxEP(M) )/&
!                                          ( (-A_M(IJK,0,M))+VXF_GS(IJK,M)+denominator(M) )
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
         else
           D_E(IJK,M) = ZERO
         endif
       END DO
     ENDIF
   END DO

 ENDIF

 RETURN
 END SUBROUTINE CALC_D_E


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_d_n(A_m, VxF_gs, VxF_ss, d_n, IER)                C
!  Purpose: calculte coefficients linking velocity correction to       C
!           pressure correction -- North                               C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-JUN-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Note:  MFIX convention: center coeff is negative, hence             C
!                            (-A_M(IJK,0,M)) > or = 0                  C
!         The NORTH face area is AXZ                                   C
!                                                                      C
!  Revision Number: 2                                                  C
!  Purpose: allow multiparticle in D_N calculation in accounting for   C
!                       the averaged Solid-Solid drag                  C
!  Author: S. Dartevelle, LANL                        Date: 28-FEB-04  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number: 3                                                  C
!  Purpose: more accurate formulas in the D_s(M) for Model_A           C
!                                                                      C
!  Author: S. Dartevelle, LANL                        Date: 17-MAR-04  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CALC_D_N(A_M, VXF_GS, VXF_SS, D_N, IER) 
!
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
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
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
      USE cutcell
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================

! QMOMK - Alberto Passalacqua
      USE qmom_kinetic_equation
! QMOMK - End

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
!
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
!
!                      Volume x average at momentum cell centers
      DOUBLE PRECISION VxF_gs(DIMENSION_3, DIMENSION_M)
!                      
      DOUBLE PRECISION d_n(DIMENSION_3, 0:DIMENSION_M)
!
!                      Average volume fraction at momentum cell centers
      DOUBLE PRECISION EPGA                                        !S. Dartevelle, LANL, Feb.2004
!          Volume x average at momentum cell centers Solid-Solid Drag
      DOUBLE PRECISION VxF_ss(DIMENSION_3, DIMENSION_LM)           !S. Dartevelle, LANL, Feb.2004
!          Usual Indices
      INTEGER           LM, M, L, LpL, Lp, J, IJK, IJKN, IN        !S. Dartevelle, LANL, Feb.2004
!          Average solid volume fraction at momentum cell centers
      DOUBLE PRECISION  EPSA(DIMENSION_M)                          !S. Dartevelle, LANL, Feb.2004
!          ratio of drag and A0 or A_solid and other sum of Solid-Solid drag
      DOUBLE PRECISION  other_ratio_1, FOA1                        !S. Dartevelle, LANL, Feb.2004
!          sum of All Solid-Gas VolxDrag
      DOUBLE PRECISION  SUM_VXF_GS                                 !S. Dartevelle, LANL, Feb.2004
!          What phase momentum equation is activated?
      LOGICAL  Pass1, Pass2                                        !S. Dartevelle, LANL, Feb.2004
!          sum of Solid M - All other Solid drag
      DOUBLE PRECISION  SUM_VXF_SS(DIMENSION_M)                    !S. Dartevelle, LANL, Feb.2004
!          numerator needed for solid phase M time EPSA
      DOUBLE PRECISION  numeratorxEP(DIMENSION_M)                  !S. Dartevelle, LANL, Feb.2004
!          denominator needed for solid phase M
      DOUBLE PRECISION  denominator(DIMENSION_M)                   !S. Dartevelle, LANL, Feb.2004
!          denominator needed for solid phase M
      DOUBLE PRECISION  other_denominator(DIMENSION_M)             !S. Dartevelle, LANL, Feb.2004
!          sum of Solid L - All other Solid VolxDrag but M
      DOUBLE PRECISION  SUM_VXF_SS_wt_M                            !S. Dartevelle, LANL, Feb.2004
!          total solids volume fraction
      DOUBLE PRECISION  EPStmp
!-----------------------------------------------
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'


   Pass1 = .FALSE.     !initialization
   Pass2 = .FALSE.
   DO M = 1, MMAX
     if (MOMENTUM_Y_EQ(0) .AND. MOMENTUM_Y_EQ(M)) then
       Pass1 = .TRUE.    !we have at least one solid phase Y-momentum equation
       GO TO 10          !with the gas phase Y-momentum equation
     elseif (MOMENTUM_Y_EQ(M)) then
       Pass2 = .TRUE.    !we have at least one solid phase Y-momentum equation
       GO TO 10          !but the gas phase Y-momentum is not solved
     endif
   END DO

   10 CONTINUE


 IF (Pass1) THEN    !Gas and at least one solid phases Y-momentum equation
!$omp   parallel do private(J,IJK, IJKN, EPGA, EPSA, EPStmp, numeratorxEP, &
!$omp&  M, L, Lp, LpL, LM, other_ratio_1, FOA1, SUM_VXF_SS, &
!$omp&  SUM_VXF_SS_wt_M, SUM_VXF_GS, other_denominator, denominator ),&
!$omp&  schedule(static)
   DO IJK = ijkstart3, ijkend3
     IF (IP_AT_N(IJK) .OR. MFLOW_AT_N(IJK)) THEN   !impermeable
       DO M= 0, MMAX
         D_N(IJK,M) = ZERO
       END DO
     ELSE
       J = J_OF(IJK)
       IJKN = NORTH_OF(IJK)
       EPGA = AVG_Y(EP_G(IJK),EP_G(IJKN),J)

       SUM_VXF_GS = ZERO
       DO M= 1, MMAX
         EPSA(M) = AVG_Y(EP_S(IJK,M),EP_S(IJKN,M),J)
         SUM_VXF_GS = SUM_VXF_GS + VXF_GS(IJK,M)              !Gas - All Solids VolxDrag summation
         SUM_VXF_SS(M) = ZERO
         do L = 1, MMAX
           IF (L .NE. M) THEN
             LM = FUNLM(L,M)
             SUM_VXF_SS(M) = SUM_VXF_SS(M) + VXF_SS(IJK,LM) !Solid M - All other Solids VolxDrag summation
           ENDIF
         end do
       END DO

       other_ratio_1  = ZERO
       DO M= 1, MMAX
         other_ratio_1 = other_ratio_1 +&
                      ( VXF_GS(IJK,M)*( (-A_M(IJK,0,M))+SUM_VXF_SS(M) )/&
                                      ( (-A_M(IJK,0,M))+VXF_GS(IJK,M)+SUM_VXF_SS(M)+SMALL_NUMBER )&
                      )
       END DO

       IF (MODEL_B) THEN   !Model B
          !Linking velocity correction coefficient to pressure - GAS Phase
         if ( ((-A_M(IJK,0,0))>SMALL_NUMBER) .OR.  (other_ratio_1>SMALL_NUMBER) ) then
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            IF(.NOT.CARTESIAN_GRID) THEN
               D_N(IJK,0) = P_SCALE*AXZ(IJK)/( (-A_M(IJK,0,0))+other_ratio_1 )
            ELSE
               D_N(IJK,0) = P_SCALE         /( (-A_M(IJK,0,0))+other_ratio_1 )
            ENDIF
! Original terms
!           D_N(IJK,0) = P_SCALE*AXZ(IJK)/( (-A_M(IJK,0,0))+other_ratio_1 )
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
         else
           D_N(IJK,0) = ZERO
         endif
         !Linking velocity correction coefficient to pressure - SOLID Phase
         DO M = 1, MMAX
           if ( MOMENTUM_Y_EQ(M) .AND. (  ((-A_M(IJK,0,M))>SMALL_NUMBER) .OR. &
                                                     (VXF_GS(IJK,M)>SMALL_NUMBER) ) )   then
             D_N(IJK,M) = D_N(IJK,0)*(&
                                       VXF_GS(IJK,M)/((-A_M(IJK,0,M))+VXF_GS(IJK,M))&
                                               )
           else
             D_N(IJK,M) = ZERO
           endif
         END DO
       ELSE                !Model A
         FOA1 = ZERO
         DO M= 1, MMAX
           FOA1 = FOA1 + (EPSA(M)*VXF_GS(IJK,M)/&
                             ( (-A_M(IJK,0,M))+VXF_GS(IJK,M)+SUM_VXF_SS(M)+SMALL_NUMBER )&
                          )
           other_denominator(M) = VXF_GS(IJK,M)*( ((-A_M(IJK,0,0))+(SUM_VXF_GS - VXF_GS(IJK,M)))/&
                                                   ((-A_M(IJK,0,0))+SUM_VXF_GS+SMALL_NUMBER)&
                                                 )
           numeratorxEP(M) = ZERO
           denominator(M)  = ZERO
           do L = 1, MMAX
             IF (L .NE. M) THEN
               LM = FUNLM(L,M)
               numeratorxEP(M) = numeratorxEP(M) + (&
                                                      VXF_SS(IJK,LM)*EPSA(L)/&
                                                      ( (-A_M(IJK,0,L))+VXF_GS(IJK,L)+SUM_VXF_SS(L)+SMALL_NUMBER )&
                                                    ) 
               SUM_VXF_SS_wt_M = ZERO
               do Lp = 1, MMAX
                 if ( (Lp .NE. L) .AND. (Lp .NE. M) ) then
                       LpL = FUNLM(Lp,L)
                   SUM_VXF_SS_wt_M = SUM_VXF_SS_wt_M + VXF_SS(IJK,LpL)     !Solid L - All other Solids VolxDrag but M summation
                 endif
               end do
               denominator(M) = denominator(M) + (&
                                                     VXF_SS(IJK,LM)*(&
                                                     ( (-A_M(IJK,0,L))+VXF_GS(IJK,L)+SUM_VXF_SS_wt_M )/&
                                                     ( (-A_M(IJK,0,L))+VXF_GS(IJK,L)+SUM_VXF_SS(L)+SMALL_NUMBER )&
                                                                    )&
                                                  )
             ENDIF
           end do
         END DO
         !Linking velocity correction coefficient to pressure - GAS Phase
         if ( ((-A_M(IJK,0,0))>SMALL_NUMBER) .OR.  (other_ratio_1>SMALL_NUMBER) ) then
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            IF(.NOT.CARTESIAN_GRID) THEN
               D_N(IJK,0) = P_SCALE*AXZ(IJK)*(EPGA+FOA1)/( (-A_M(IJK,0,0))+other_ratio_1 ) 
            ELSE
               D_N(IJK,0) = P_SCALE          *(EPGA+FOA1)/( (-A_M(IJK,0,0))+other_ratio_1 )
            ENDIF
! Original terms
!           D_N(IJK,0) = P_SCALE*AXZ(IJK)*(EPGA+FOA1)/( (-A_M(IJK,0,0))+other_ratio_1 )
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
         else
           D_N(IJK,0) = ZERO
         endif
         !Linking velocity correction coefficient to pressure - SOLID Phase
         DO M = 1, MMAX
           if ( MOMENTUM_Y_EQ(M) .AND. ( (-A_M(IJK,0,M)>SMALL_NUMBER)        .OR. &
                                                (other_denominator(M)>SMALL_NUMBER) .OR. &
                                                (denominator(M)>SMALL_NUMBER) ) )     then
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            IF(.NOT.CARTESIAN_GRID) THEN
               D_N(IJK,M) = P_SCALE*AXZ(IJK)*(&
                           ( EPSA(M) + numeratorxEP(M) + (VXF_GS(IJK,M)*EPGA/((-A_M(IJK,0,0))+SUM_VXF_GS+SMALL_NUMBER)) )/&
                           ( (-A_M(IJK,0,M))+other_denominator(M)+denominator(M) )&
                                                      )
            ELSE
               D_N(IJK,M) = P_SCALE         *(&
                         ( EPSA(M) + numeratorxEP(M) + (VXF_GS(IJK,M)*EPGA/((-A_M(IJK,0,0))+SUM_VXF_GS+SMALL_NUMBER)) )/&
                         ( (-A_M(IJK,0,M))+other_denominator(M)+denominator(M) )&
                                                      )
            ENDIF
! Original terms
!             D_N(IJK,M) = P_SCALE*AXZ(IJK)*(&
!                        ( EPSA(M) + numeratorxEP(M) + (VXF_GS(IJK,M)*EPGA/((-A_M(IJK,0,0))+SUM_VXF_GS+SMALL_NUMBER)) )/&
!                        ( (-A_M(IJK,0,M))+other_denominator(M)+denominator(M) )&
!                                                   )
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
           else
             D_N(IJK,M) = ZERO
           endif
         END DO
       ENDIF    !end of Model_B/Model_A if then condition
     ENDIF
   END DO

 ELSE IF (MOMENTUM_Y_EQ(0)) THEN    !the solid Y-momentum equations are not solved
                                    !only the gas phase
!$omp   parallel do private(IJK, J, IJKN, EPGA, M, SUM_VXF_GS), &
!$omp&  schedule(static)
   DO IJK = ijkstart3, ijkend3
     IF (IP_AT_N(IJK) .OR. MFLOW_AT_N(IJK)) THEN   !impermeable
       D_N(IJK,0) = ZERO
     ELSE
       J = J_OF(IJK)
       IJKN = NORTH_OF(IJK)
       EPGA = AVG_Y(EP_G(IJK),EP_G(IJKN),J)

       ! QMOMK - Alberto Passalacqua
       ! Added check for QMOMK
       SUM_VXF_GS = ZERO
       IF (.NOT. QMOMK) THEN
        DO M = 1, MMAX
         SUM_VXF_GS = SUM_VXF_GS + VXF_GS(IJK,M)              !Gas - All Solids VolxDrag summation
        END DO
       ELSE
        DO IN = 1, QMOMK_NN
          DO M = 1, MMAX           
            SUM_VXF_GS = SUM_VXF_GS + AVG_Y(QMOMK_F_GS(IN,IJK,M),QMOMK_F_GS(IN,IJKN,M),J)*VOL_V(IJK)
          END DO
        END DO
       END IF
       ! QMOMK - End

       IF ( ((-A_M(IJK,0,0))>SMALL_NUMBER) .OR. (SUM_VXF_GS>SMALL_NUMBER) ) THEN
         IF (MODEL_B) THEN
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            IF(.NOT.CARTESIAN_GRID) THEN
               D_N(IJK,0) = P_SCALE*AXZ(IJK)/((-A_M(IJK,0,0))+SUM_VXF_GS)
            ELSE
               D_N(IJK,0) = P_SCALE         /((-A_M(IJK,0,0))+SUM_VXF_GS)
            ENDIF
! Original terms
!           D_N(IJK,0) = P_SCALE*AXZ(IJK)/((-A_M(IJK,0,0))+SUM_VXF_GS)
         ELSE
            IF(.NOT.CARTESIAN_GRID) THEN
               D_N(IJK,0) = P_SCALE*AXZ(IJK)*EPGA/((-A_M(IJK,0,0))+SUM_VXF_GS)
            ELSE
               D_N(IJK,0) = P_SCALE         *EPGA/((-A_M(IJK,0,0))+SUM_VXF_GS)
            ENDIF
! Original terms
!              D_N(IJK,0) = P_SCALE*AXZ(IJK)*EPGA/((-A_M(IJK,0,0))+SUM_VXF_GS)
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
         ENDIF
       ELSE
         D_N(IJK,0) = ZERO
       ENDIF
     ENDIF
   END DO

 ELSE IF (Pass2) THEN    !at least one solid phase momentum Y-equation is solved
                         !but the gas phase Y-momentum is not solved
!$omp    parallel do private(IJK, J, IJKN, EPSA, EPStmp, L, Lp, M, SUM_VXF_SS, &
!$omp&   numeratorxEP, denominator, SUM_VXF_SS_wt_M), &
!$omp&  schedule(static)
   DO IJK = ijkstart3, ijkend3
     IF (IP_AT_N(IJK) .OR. MFLOW_AT_N(IJK) .OR. MODEL_B) THEN
       DO M= 1, MMAX
         D_N(IJK,M) = ZERO
       END DO
     ELSE
       J = J_OF(IJK)
       IJKN = NORTH_OF(IJK)

       DO M= 1, MMAX
         EPSA(M) = AVG_Y(EP_S(IJK,M),EP_S(IJKN,M),J)     
         SUM_VXF_SS(M) = ZERO
         do L = 1, MMAX
           IF (L .NE. M) THEN
             LM = FUNLM(L,M)
             SUM_VXF_SS(M) = SUM_VXF_SS(M) + VXF_SS(IJK,LM)  !Solid M - All other Solids VolxDrag summation
           ENDIF
         end do
       END DO
       DO M= 1, MMAX
         numeratorxEP(M) = ZERO
         denominator(M)  = ZERO
         do L = 1, MMAX
           IF (L .NE. M) THEN
             LM = FUNLM(L,M)
             numeratorxEP(M) = numeratorxEP(M) + (&
                                                    VXF_SS(IJK,LM)*EPSA(L)/&
                                                   ( (-A_M(IJK,0,L))+VXF_GS(IJK,L)+SUM_VXF_SS(L)+SMALL_NUMBER )&
                                                  )
             SUM_VXF_SS_wt_M = ZERO
             do Lp = 1, MMAX
               if ( (Lp .NE. L) .AND. (Lp .NE. M) ) then
                 LpL = FUNLM(Lp,L)
                 SUM_VXF_SS_wt_M = SUM_VXF_SS_wt_M + VXF_SS(IJK,LpL)     !Solid L - All other Solids VolxDrag but M summation
               endif
             end do
             denominator(M) = denominator(M) + (&
                                                  VXF_SS(IJK,LM)*(&
                                                   ( (-A_M(IJK,0,L))+VXF_GS(IJK,L)+SUM_VXF_SS_wt_M )/&
                                                   ( (-A_M(IJK,0,L))+VXF_GS(IJK,L)+SUM_VXF_SS(L)+SMALL_NUMBER )&
                                                                 )&
                                                )
           ENDIF
         end do
       END DO
       !Linking velocity correction coefficient to pressure - SOLID Phase (Model_A only)
       DO M = 1, MMAX
         if ( MOMENTUM_Y_EQ(M) .AND. ( (-A_M(IJK,0,M)>SMALL_NUMBER) .OR. &
                                        (VXF_GS(IJK,M)>SMALL_NUMBER) .OR. &
                                        (denominator(M)>SMALL_NUMBER) ) ) then
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            IF(.NOT.CARTESIAN_GRID) THEN
               D_N(IJK,M) = P_SCALE*AXZ(IJK)*( EPSA(M) + numeratorxEP(M) )/&
                                             ( (-A_M(IJK,0,M))+VXF_GS(IJK,M)+denominator(M) )
            ELSE
               D_N(IJK,M) = P_SCALE         *( EPSA(M) + numeratorxEP(M) )/&
                                             ( (-A_M(IJK,0,M))+VXF_GS(IJK,M)+denominator(M) )
            ENDIF
! Original terms
!           D_N(IJK,M) = P_SCALE*AXZ(IJK)*( EPSA(M) + numeratorxEP(M) )/&
!                                          ( (-A_M(IJK,0,M))+VXF_GS(IJK,M)+denominator(M) )
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
         else
           D_N(IJK,M) = ZERO
         endif
       END DO
     ENDIF
   END DO

 ENDIF

 RETURN
 END SUBROUTINE CALC_D_N


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_d_t(A_m, VxF_gs, VxF_ss, d_t, IER)                C
!  Purpose: calculte coefficients linking velocity correction to       C
!           pressure correction -- Top                                 C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-JUN-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Note:  MFIX convention: center coeff is negative, hence             C
!                            (-A_M(IJK,0,M)) > or = 0                  C
!         The TOP face area is AXY                                     C
!                                                                      C
!  Revision Number: 2                                                  C
!  Purpose: allow multiparticle in D_T calculation in accounting for   C
!                       the averaged Solid-Solid drag                  C
!  Author: S. Dartevelle, LANL                        Date: 28-FEB-04  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number: 3                                                  C
!  Purpose: more accurate formulas in the D_s(M) for Model_A           C
!                                                                      C
!  Author: S. Dartevelle, LANL                        Date: 17-MAR-04  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CALC_D_T(A_M, VXF_GS, VXF_SS, D_T, IER)
!
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
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
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
      USE cutcell
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================

! QMOMK - Alberto Passalacqua
      USE qmom_kinetic_equation
! QMOMK - End

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
!
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
!
!                      Volume x average at momentum cell centers
      DOUBLE PRECISION VxF_gs(DIMENSION_3, DIMENSION_M)
!                      
      DOUBLE PRECISION d_t(DIMENSION_3, 0:DIMENSION_M)
!
!                      Average volume fraction at momentum cell centers
      DOUBLE PRECISION EPGA                                        !S. Dartevelle, LANL, Feb.2004
!          Volume x average at momentum cell centers Solid-Solid Drag
      DOUBLE PRECISION VxF_ss(DIMENSION_3, DIMENSION_LM)           !S. Dartevelle, LANL, Feb.2004
!          Usual Indices
      INTEGER           LM, M, L, LpL, Lp, K, IJK, IJKT, IN        !S. Dartevelle, LANL, Feb.2004
!          Average solid volume fraction at momentum cell centers
      DOUBLE PRECISION  EPSA(DIMENSION_M)                          !S. Dartevelle, LANL, Feb.2004
!          ratio of drag and A0 or A_solid and other sum of Solid-Solid drag
      DOUBLE PRECISION  other_ratio_1, FOA1                        !S. Dartevelle, LANL, Feb.2004
!          sum of All Solid-Gas VolxDrag
      DOUBLE PRECISION  SUM_VXF_GS                                 !S. Dartevelle, LANL, Feb.2004
!          What phase momentum equation is activated?
      LOGICAL  Pass1, Pass2                                        !S. Dartevelle, LANL, Feb.2004
!          sum of Solid M - All other Solid drag
      DOUBLE PRECISION  SUM_VXF_SS(DIMENSION_M)                    !S. Dartevelle, LANL, Feb.2004
!          numerator needed for solid phase M time EPSA
      DOUBLE PRECISION  numeratorxEP(DIMENSION_M)                  !S. Dartevelle, LANL, Feb.2004
!          denominator needed for solid phase M
      DOUBLE PRECISION  denominator(DIMENSION_M)                   !S. Dartevelle, LANL, Feb.2004
!          denominator needed for solid phase M
      DOUBLE PRECISION  other_denominator(DIMENSION_M)             !S. Dartevelle, LANL, Feb.2004
!          sum of Solid L - All other Solid VolxDrag but M
      DOUBLE PRECISION  SUM_VXF_SS_wt_M                            !S. Dartevelle, LANL, Feb.2004
!          tmp variable for total solids volume fraction
      DOUBLE PRECISION  EPStmp
!-----------------------------------------------
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'


   Pass1 = .FALSE.     !initialization
   Pass2 = .FALSE.
   DO M = 1, MMAX
     if (MOMENTUM_Z_EQ(0) .AND. MOMENTUM_Z_EQ(M)) then
       Pass1 = .TRUE.    !we have at least one solid phase Z-momentum equation
       GO TO 10          !with the gas phase Z-momentum equation
     elseif (MOMENTUM_Z_EQ(M)) then
       Pass2 = .TRUE.    !we have at least one solid phase Z-momentum equation
       GO TO 10          !but the gas phase Z-momentum is not solved
     endif
   END DO

   10 CONTINUE


 IF (Pass1) THEN    !Gas and at least one solid phases Z-momentum equation
!$omp   parallel do private(K,IJK, IJKT, EPGA, EPSA, EPStmp, numeratorxEP, &
!$omp&  M, L, Lp, LpL, LM, other_ratio_1, FOA1, SUM_VXF_SS, &
!$omp&  SUM_VXF_SS_wt_M, SUM_VXF_GS, other_denominator, denominator ),&
!$omp&  schedule(static)
   DO IJK = ijkstart3, ijkend3
     IF (IP_AT_T(IJK) .OR. MFLOW_AT_T(IJK)) THEN   !impermeable
       DO M= 0, MMAX
         D_T(IJK,M) = ZERO
       END DO
     ELSE
       K = K_OF(IJK)
       IJKT = TOP_OF(IJK)
       EPGA = AVG_Z(EP_G(IJK),EP_G(IJKT),K)

       SUM_VXF_GS = ZERO
       DO M= 1, MMAX
         EPSA(M) = AVG_Z(EP_S(IJK,M),EP_S(IJKT,M),K)
         SUM_VXF_GS = SUM_VXF_GS + VXF_GS(IJK,M)              !Gas - All Solids VolxDrag summation
         SUM_VXF_SS(M) = ZERO
         do L = 1, MMAX
           IF (L .NE. M) THEN
             LM = FUNLM(L,M)
             SUM_VXF_SS(M) = SUM_VXF_SS(M) + VXF_SS(IJK,LM) !Solid M - All other Solids VolxDrag summation
           ENDIF
         end do
       END DO

       other_ratio_1  = ZERO
       DO M= 1, MMAX
         other_ratio_1 = other_ratio_1 +&
                      ( VXF_GS(IJK,M)*( (-A_M(IJK,0,M))+SUM_VXF_SS(M) )/&
                                      ( (-A_M(IJK,0,M))+VXF_GS(IJK,M)+SUM_VXF_SS(M)+SMALL_NUMBER )&
                      )
       END DO

       IF (MODEL_B) THEN   !Model B
         !Linking velocity correction coefficient to pressure - GAS Phase
         if ( ((-A_M(IJK,0,0))>SMALL_NUMBER) .OR.  (other_ratio_1>SMALL_NUMBER) ) then
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            IF(.NOT.CARTESIAN_GRID) THEN
               D_T(IJK,0) = P_SCALE*AXY(IJK)/( (-A_M(IJK,0,0))+other_ratio_1 )
            ELSE
               D_T(IJK,0) = P_SCALE         /( (-A_M(IJK,0,0))+other_ratio_1 )
            ENDIF
! Original terms
!           D_T(IJK,0) = P_SCALE*AXY(IJK)/( (-A_M(IJK,0,0))+other_ratio_1 )
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
         else
           D_T(IJK,0) = ZERO
         endif
         !Linking velocity correction coefficient to pressure - SOLID Phase
         DO M = 1, MMAX
           if ( MOMENTUM_Z_EQ(M) .AND. (  ((-A_M(IJK,0,M))>SMALL_NUMBER) .OR. &
                                                     (VXF_GS(IJK,M)>SMALL_NUMBER) ) )   then
             D_T(IJK,M) = D_T(IJK,0)*(&
                                       VXF_GS(IJK,M)/((-A_M(IJK,0,M))+VXF_GS(IJK,M))&
                                               )
           else
             D_T(IJK,M) = ZERO
           endif
         END DO
       ELSE                !Model A
         FOA1 = ZERO
         DO M= 1, MMAX
           FOA1 = FOA1 + (EPSA(M)*VXF_GS(IJK,M)/&
                             ( (-A_M(IJK,0,M))+VXF_GS(IJK,M)+SUM_VXF_SS(M)+SMALL_NUMBER )&
                          )
           other_denominator(M) = VXF_GS(IJK,M)*( ((-A_M(IJK,0,0))+(SUM_VXF_GS - VXF_GS(IJK,M)))/&
                                                   ((-A_M(IJK,0,0))+SUM_VXF_GS+SMALL_NUMBER)&
                                                 )
           numeratorxEP(M) = ZERO
           denominator(M)  = ZERO
           do L = 1, MMAX
             IF (L .NE. M) THEN
               LM = FUNLM(L,M)
               numeratorxEP(M) = numeratorxEP(M) + (&
                                                      VXF_SS(IJK,LM)*EPSA(L)/&
                                                      ( (-A_M(IJK,0,L))+VXF_GS(IJK,L)+SUM_VXF_SS(L)+SMALL_NUMBER )&
                                                    ) 
               SUM_VXF_SS_wt_M = ZERO
               do Lp = 1, MMAX
                 if ( (Lp .NE. L) .AND. (Lp .NE. M) ) then
                   LpL = FUNLM(Lp,L)
                   SUM_VXF_SS_wt_M = SUM_VXF_SS_wt_M + VXF_SS(IJK,LpL)     !Solid L - All other Solids VolxDrag but M summation
                 endif
               end do
               denominator(M) = denominator(M) + (&
                                                     VXF_SS(IJK,LM)*(&
                                                     ( (-A_M(IJK,0,L))+VXF_GS(IJK,L)+SUM_VXF_SS_wt_M )/&
                                                     ( (-A_M(IJK,0,L))+VXF_GS(IJK,L)+SUM_VXF_SS(L)+SMALL_NUMBER )&
                                                                    )&
                                                  )
             ENDIF
           end do
         END DO
         !Linking velocity correction coefficient to pressure - GAS Phase
         if ( ((-A_M(IJK,0,0))>SMALL_NUMBER) .OR.  (other_ratio_1>SMALL_NUMBER) ) then
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            IF(.NOT.CARTESIAN_GRID) THEN
               D_T(IJK,0) = P_SCALE*AXY(IJK)*(EPGA+FOA1)/( (-A_M(IJK,0,0))+other_ratio_1 )
            ELSE
               D_T(IJK,0) = P_SCALE         *(EPGA+FOA1)/( (-A_M(IJK,0,0))+other_ratio_1 )
            ENDIF
! Original terms
!           D_T(IJK,0) = P_SCALE*AXY(IJK)*(EPGA+FOA1)/( (-A_M(IJK,0,0))+other_ratio_1 )
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
         else
           D_T(IJK,0) = ZERO
         endif
         !Linking velocity correction coefficient to pressure - SOLID Phase
         DO M = 1, MMAX
           if ( MOMENTUM_Z_EQ(M) .AND. ( (-A_M(IJK,0,M)>SMALL_NUMBER)        .OR. &
                                                (other_denominator(M)>SMALL_NUMBER) .OR. &
                                                (denominator(M)>SMALL_NUMBER) ) )     then
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            IF(.NOT.CARTESIAN_GRID) THEN
               D_T(IJK,M) = P_SCALE*AXY(IJK)*(&
                          ( EPSA(M) + numeratorxEP(M) + (VXF_GS(IJK,M)*EPGA/((-A_M(IJK,0,0))+SUM_VXF_GS+SMALL_NUMBER)) )/&
                          ( (-A_M(IJK,0,M))+other_denominator(M)+denominator(M) )&
                                                     )
            ELSE
               D_T(IJK,M) = P_SCALE        *(&
                         ( EPSA(M) + numeratorxEP(M) + (VXF_GS(IJK,M)*EPGA/((-A_M(IJK,0,0))+SUM_VXF_GS+SMALL_NUMBER)) )/&
                         ( (-A_M(IJK,0,M))+other_denominator(M)+denominator(M) )&
                                                    )
            ENDIF
! Original terms
!             D_T(IJK,M) = P_SCALE*AXY(IJK)*(&
!                        ( EPSA(M) + numeratorxEP(M) + (VXF_GS(IJK,M)*EPGA/((-A_M(IJK,0,0))+SUM_VXF_GS+SMALL_NUMBER)) )/&
!                        ( (-A_M(IJK,0,M))+other_denominator(M)+denominator(M) )&
!                                                   )
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
           else
             D_T(IJK,M) = ZERO
           endif
         END DO
       ENDIF    !end of Model_B/Model_A if then condition
     ENDIF
   END DO

 ELSE IF (MOMENTUM_Z_EQ(0)) THEN    !the solid Z-momentum equations are not solved
                                    !only the gas phase Z-momentum is solved
!$omp   parallel do private(IJK, K, IJKT, EPGA, M, SUM_VXF_GS), &
!$omp&  schedule(static)
   DO IJK = ijkstart3, ijkend3
     IF (IP_AT_T(IJK) .OR. MFLOW_AT_T(IJK)) THEN   !impermeable
       D_T(IJK,0) = ZERO
     ELSE
       K = K_OF(IJK)
       IJKT = TOP_OF(IJK)
       EPGA = AVG_Z(EP_G(IJK),EP_G(IJKT),K)

       ! QMOMK - Alberto Passalacqua
       ! Added check for QMOMK
       SUM_VXF_GS = ZERO
       IF (.NOT. QMOMK) THEN
        DO M = 1, MMAX
         SUM_VXF_GS = SUM_VXF_GS + VXF_GS(IJK,M)              !Gas - All Solids VolxDrag summation
        END DO
       ELSE      
        DO IN = 1, QMOMK_NN
          DO M = 1, MMAX           
            SUM_VXF_GS = SUM_VXF_GS + AVG_Z(QMOMK_F_GS(IN,IJK,M),QMOMK_F_GS(IN,IJKT,M),K)*VOL_W(IJK)
          END DO
        END DO
       END IF
       ! QMOMK - End

       IF ( ((-A_M(IJK,0,0))>SMALL_NUMBER) .OR. (SUM_VXF_GS>SMALL_NUMBER) ) THEN
         IF (MODEL_B) THEN
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            IF(.NOT.CARTESIAN_GRID) THEN
               D_T(IJK,0) = P_SCALE*AXY(IJK)/((-A_M(IJK,0,0))+SUM_VXF_GS)
            ELSE
               D_T(IJK,0) = P_SCALE         /((-A_M(IJK,0,0))+SUM_VXF_GS)
            ENDIF
! Original terms
!           D_T(IJK,0) = P_SCALE*AXY(IJK)/((-A_M(IJK,0,0))+SUM_VXF_GS)
         ELSE
            IF(.NOT.CARTESIAN_GRID) THEN
               D_T(IJK,0) = P_SCALE*AXY(IJK)*EPGA/((-A_M(IJK,0,0))+SUM_VXF_GS)
            ELSE
               D_T(IJK,0) = P_SCALE         *EPGA/((-A_M(IJK,0,0))+SUM_VXF_GS)
            ENDIF
! Original terms
!           D_T(IJK,0) = P_SCALE*AXY(IJK)*EPGA/((-A_M(IJK,0,0))+SUM_VXF_GS)
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
         ENDIF
       ELSE
         D_T(IJK,0) = ZERO
       ENDIF
     ENDIF
   END DO

 ELSE IF (Pass2) THEN    !at least one solid phase momentum Z-equation is solved
                         !but the gas phase Z-momentum is not solved
!$omp    parallel do private(IJK, K, IJKT, EPSA, EPStmp, L, Lp, M, SUM_VXF_SS, &
!$omp&   numeratorxEP, denominator, SUM_VXF_SS_wt_M), &
!$omp&  schedule(static)
   DO IJK = ijkstart3, ijkend3
     IF (IP_AT_T(IJK) .OR. MFLOW_AT_T(IJK) .OR. MODEL_B) THEN
       DO M= 1, MMAX
         D_T(IJK,M) = ZERO
       END DO
     ELSE
       K = K_OF(IJK)
       IJKT = TOP_OF(IJK)

       DO M= 1, MMAX
         EPSA(M) = AVG_Z(EP_S(IJK,M),EP_S(IJKT,M),K)   
         SUM_VXF_SS(M) = ZERO
         do L = 1, MMAX
           IF (L .NE. M) THEN
             LM = FUNLM(L,M)
             SUM_VXF_SS(M) = SUM_VXF_SS(M) + VXF_SS(IJK,LM)  !Solid M - All other Solids VolxDrag summation
           ENDIF
         end do
       END DO
       DO M= 1, MMAX
         numeratorxEP(M) = ZERO
         denominator(M)  = ZERO
         do L = 1, MMAX
           IF (L .NE. M) THEN
             LM = FUNLM(L,M)
             numeratorxEP(M) = numeratorxEP(M) + (&
                                                    VXF_SS(IJK,LM)*EPSA(L)/&
                                                   ( (-A_M(IJK,0,L))+VXF_GS(IJK,L)+SUM_VXF_SS(L)+SMALL_NUMBER )&
                                                  )
             SUM_VXF_SS_wt_M = ZERO
             do Lp = 1, MMAX
               if ( (Lp .NE. L) .AND. (Lp .NE. M) ) then
                 LpL = FUNLM(Lp,L)
                 SUM_VXF_SS_wt_M = SUM_VXF_SS_wt_M + VXF_SS(IJK,LpL)     !Solid L - All other Solids VolxDrag but M summation
               endif
             end do
             denominator(M) = denominator(M) + (&
                                                  VXF_SS(IJK,LM)*(&
                                                   ( (-A_M(IJK,0,L))+VXF_GS(IJK,L)+SUM_VXF_SS_wt_M )/&
                                                   ( (-A_M(IJK,0,L))+VXF_GS(IJK,L)+SUM_VXF_SS(L)+SMALL_NUMBER )&
                                                                 )&
                                                )
           ENDIF
         end do
       END DO
       !Linking velocity correction coefficient to pressure - SOLID Phase (Model_A only)
       DO M = 1, MMAX
         if ( MOMENTUM_Z_EQ(M) .AND. ( (-A_M(IJK,0,M)>SMALL_NUMBER) .OR. &
                                        (VXF_GS(IJK,M)>SMALL_NUMBER) .OR. &
                                        (denominator(M)>SMALL_NUMBER) ) ) then
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            IF(.NOT.CARTESIAN_GRID) THEN
               D_T(IJK,M) = P_SCALE*AXY(IJK)*( EPSA(M) + numeratorxEP(M) )/&
                                             ( (-A_M(IJK,0,M))+VXF_GS(IJK,M)+denominator(M) )
            ELSE
             D_T(IJK,M) = P_SCALE         *( EPSA(M) + numeratorxEP(M) )/&
                                          ( (-A_M(IJK,0,M))+VXF_GS(IJK,M)+denominator(M) )
            ENDIF
! Original terms
!            D_T(IJK,M) = P_SCALE*AXY(IJK)*( EPSA(M) + numeratorxEP(M) )/&
!                                          ( (-A_M(IJK,0,M))+VXF_GS(IJK,M)+denominator(M) )
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
         else
           D_T(IJK,M) = ZERO
         endif
       END DO
     ENDIF
   END DO

 ENDIF

 RETURN
 END SUBROUTINE CALC_D_T

!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization 
!// 350 Changed do loop limits: 1,ijkmax2-> ijkstart3, ijkend3
