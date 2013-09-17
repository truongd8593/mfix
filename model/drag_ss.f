!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: DRAG_ss                                                 C
!  Purpose: This module computes the coefficient of drag between       C
!     two solids phases (M and L).                                     C
!                                                                      C
!  Author: M. Syamlal                                 Date: 20-MAY-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!     Gera, D., Syamlal, M., O'Brien T.J. 2004. International Journal  C
!        of Multiphase Flow, 30, p419-428.                             C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE DRAG_SS(L, M, IER) 

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE constant
      USE fldvar
      USE geometry
      USE indices
      USE physprop
      USE compar 
      USE sendrecv 
      USE drag
      USE discretelement
      IMPLICIT NONE
!-----------------------------------------------
! Dummy Arguments
!-----------------------------------------------
! Index of solids phases 
      INTEGER, INTENT(IN) :: L, M 
! Error index 
      INTEGER, INTENT(INOUT) :: IER 
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! indices 
      INTEGER :: I, IJK, IMJK, IJMK, IJKM
      INTEGER :: CM, DM 
! index for storing solids-solids drag coefficients in the upper 
! triangle of the matrix
      INTEGER :: LM 
! cell center value of U_sm, U_sl, V_sm, V_sl, W_sm, W_sl
      DOUBLE PRECISION :: USCM, USCL, VSCM, VSCL, WSCM, WSCL
! relative velocity between solids phase m and l 
      DOUBLE PRECISION :: VREL 
! particle diameters of phase M and phase L 
      DOUBLE PRECISION :: D_pm, D_pl 
! particle densities of phase M and phase L 
      DOUBLE PRECISION :: RO_M, RO_L 
! radial distribution function between phases M and L
      DOUBLE PRECISION :: G0_ML
! void fraction and solids volume fraction
      DOUBLE PRECISION :: EPg, EPs
! sum over all phases of ratio volume fraction over particle diameter
      DOUBLE PRECISION :: EPSoDP
! solid-solid drag coefficient 
      DOUBLE PRECISION :: lDss
!-----------------------------------------------
! External Functions
!-----------------------------------------------
      DOUBLE PRECISION, EXTERNAL :: G_0 
!-----------------------------------------------
! Include statement functions
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'
!-----------------------------------------------

      LM = FUNLM(L,M)
   
      DO IJK = ijkstart3, ijkend3

         IF (.NOT.WALL_AT(IJK)) THEN 
! Evaluate at all flow boundaries and fluid cells
! This is unlike the fluid-solid drag coefficient, which is only
! evluated in fluid cells and pressure inflow cells

            I = I_OF(IJK) 
            IMJK = IM_OF(IJK) 
            IJMK = JM_OF(IJK) 
            IJKM = KM_OF(IJK) 

! calculating velocity components at i, j, k (cell center)
            USCL = AVG_X_E(U_S(IMJK,L),U_S(IJK,L),I) 
            VSCL = AVG_Y_N(V_S(IJMK,L),V_S(IJK,L)) 
            WSCL = AVG_Z_T(W_S(IJKM,L),W_S(IJK,L)) 

            USCM = AVG_X_E(U_S(IMJK,M),U_S(IJK,M),I) 
            VSCM = AVG_Y_N(V_S(IJMK,M),V_S(IJK,M)) 
            WSCM = AVG_Z_T(W_S(IJKM,M),W_S(IJK,M)) 

! magnitude of solids-solids relative velocity
            VREL = SQRT((USCL - USCM)**2 + (VSCL - VSCM)**2 + &
                        (WSCL - WSCM)**2) 

! setting aliases for easy reference
            D_PM = D_P(IJK,M) 
            D_PL = D_P(IJK,L)
            RO_M = RO_S(IJK,M)
            RO_L = RO_S(IJK,L)

            IF (DES_CONTINUUM_HYBRID) THEN
! evaluating g0 - taken from G_0.f subroutine (lebowitz form)
! this section is needed to account for all solids phases until g0 for
! multiple solids types (i.e. discrete & continuum) can be addressed
! more effectively.
               EPSoDP = ZERO
               DO CM = 1, MMAX
                  EPS = EP_s(IJK, CM)
                  EPSoDP = EPSoDP + EPS / D_p(IJK,CM)
               ENDDO
               DO DM = 1, DES_MMAX
                  EPS = DES_ROP_S(IJK,DM)/DES_RO_S(DM)
                  EPSoDP = EPSoDP + EPS / DES_D_p0(DM)
               ENDDO
               EPg = EP_g(IJK)
               G0_ML = ONE/EPg + 3.0d0*EPSoDP*D_pM*D_PL / &
                  (EPg*EPg *(D_pM + D_pL))
            ELSE
               G0_ML = G_0(IJK,L,M)
            ENDIF

! determining the solids-solids 'drag coefficient'
            CALL DRAG_SS_SYAM(lDss,D_PM,D_PL,RO_M,RO_L,G0_ML,VREL)

            F_SS(IJK,LM) = lDss*ROP_S(IJK,M)*ROP_S(IJK,L)

! Gera: accounting for particle-particle drag due to enduring contact in a 
! close-packed system. Note the check for mmax >= 2 below is unnecessary
! since this routine will not be entered unless mmax >=2
            IF(CLOSE_PACKED(M) .AND. CLOSE_PACKED(L) .AND. &
              (MMAX >= 2))  F_SS(IJK,LM) = F_SS(IJK,LM) + &
               SEGREGATION_SLOPE_COEFFICIENT*P_star(IJK) 

          ELSE   ! elseif (.not.wall_at(ijk))

            F_SS(IJK,LM) = ZERO

         ENDIF   ! end if (.not.wall_at(ijk))

      ENDDO    ! end do (ijk=ijkstart3,ijkend3)
           
      RETURN  
      END SUBROUTINE DRAG_SS 

      

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: DRAG_SS_SYAM                                            C
!  Purpose: Calculate the solids-solids drag coefficient between a     C
!           continuous solids phase and discrete solids                C
!                                                                      C
!  Literature/Document References:                                     C
!     M. Syamlal. 1987. The particle-particle drag term in a           C
!        multiparticle model of fluidization. Technical Report.        C
!        DOE/MC/21353-2373. Office of Fossil Energy, Morgantown        C
!        Energy Technology Center, Morgantown, West Virginia.          C
!                                                                      C
!                                                                      C
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC

      SUBROUTINE DRAG_SS_SYAM(lDss,D_PM,D_PL,RO_M,RO_L, G0_ML, VREL)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param 
      USE param1
      USE constant
      use run
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! drag coefficient
      DOUBLE PRECISION, INTENT(OUT) :: ldss
! particle diameter of solids phase M
      DOUBLE PRECISION, INTENT(IN) :: D_PM
! particle diameter of solids phase L
      DOUBLE PRECISION, INTENT(IN) :: D_PL
! particle density of solids phase M
      DOUBLE PRECISION, INTENT(IN) :: RO_M
! particle density of solids phase L
      DOUBLE PRECISION, INTENT(IN) :: RO_L
! radial distribution function between phases M and L
      DOUBLE PRECISION, INTENT(IN) :: G0_ML
! relative velocity between solids phase m and l 
      DOUBLE PRECISION, INTENT(IN) :: VREL
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Sum of particle diameters 
      DOUBLE PRECISION :: DPSUM
! Intermediate calculation
      DOUBLE PRECISION :: const      
!-----------------------------------------------
! External functions
!-----------------------------------------------
!-----------------------------------------------
     
      DPSUM = D_PL + D_PM

      const = 3.d0*(ONE + C_E)*(PI/2.d0 + C_F*PI*PI/8.d0)*&
         DPSUM**2/(2.d0*PI*(RO_L*D_PL**3+RO_M*D_PM**3)) 

      ldss = const * G0_ML * VREL
         
      RETURN
      END SUBROUTINE DRAG_SS_SYAM

