!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: DRAG_ss                                                 C
!  Purpose: This module computes the coefficient of drag between       C
!           solids phase M and solids phase L                          C
!                                                                      C
!  Author: M. Syamlal                                 Date: 20-MAY-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!     M. Syamlal. 1987. The particle-particle drag term in a           C
!     multiparticle model of fluidization. Technical Report.           C
!     DOE/MC/21353-2373. Office of Fossil Energy, Morgantown Energy    C
!     Technology Center, Morgantown, West Virginia.                    C
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
! Relative velocity between solids phase m and l 
      DOUBLE PRECISION :: VREL 
! Indices 
      INTEGER :: I, IJK, IMJK, IJMK, IJKM 
! Index for storing solids-solids drag coefficients in the upper 
! triangle of the matrix
      INTEGER :: LM 
! Cell center value of U_sm, U_sl, V_sm, V_sl, W_sm, W_sl
      DOUBLE PRECISION :: USCM, USCL, VSCM, VSCL, WSCM, WSCL
! Particle diameters of phase M and phase L 
      DOUBLE PRECISION :: D_pm, D_pl 
! Sum of particle diameters 
      DOUBLE PRECISION :: DPSUM 
! Drag coefficient 
      DOUBLE PRECISION :: CONST 
!-----------------------------------------------
! External Functions
!-----------------------------------------------
      DOUBLE PRECISION, EXTERNAL :: G_0 
!-----------------------------------------------
! Include statement functions
!-----------------------------------------------
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
!-----------------------------------------------

      LM = FUNLM(L,M)
   
      DO IJK = ijkstart3, ijkend3  

         IF (.NOT.WALL_AT(IJK)) THEN 
            I = I_OF(IJK) 
            IMJK = IM_OF(IJK) 
            IJMK = JM_OF(IJK) 
            IJKM = KM_OF(IJK) 

! Calculate velocity components at i, j, k
            USCL = AVG_X_E(U_S(IMJK,L),U_S(IJK,L),I) 
            VSCL = AVG_Y_N(V_S(IJMK,L),V_S(IJK,L)) 
            WSCL = AVG_Z_T(W_S(IJKM,L),W_S(IJK,L)) 
            USCM = AVG_X_E(U_S(IMJK,M),U_S(IJK,M),I) 
            VSCM = AVG_Y_N(V_S(IJMK,M),V_S(IJK,M)) 
            WSCM = AVG_Z_T(W_S(IJKM,M),W_S(IJK,M)) 

! magnitude of solids-solids relative velocity
            VREL = SQRT((USCL - USCM)**2 + (VSCL - VSCM)**2 + (WSCL - WSCM)**2) 

            D_PM = D_P(IJK,M) 
            D_PL = D_P(IJK,L) 
            DPSUM = D_PL + D_PM        
       
            CONST = 3.d0*(ONE + C_E)*(PI/2.d0 + C_F*PI*PI/8.d0)*&
               DPSUM**2/(2.d0*PI*(RO_S(L)*D_PL**3+RO_S(M)*D_PM**3)) 
            F_SS(IJK,LM) = CONST*ROP_S(IJK,L)*ROP_S(IJK,M)*&
               G_0(IJK,L,M)*VREL 

! Gera: account for particle-particle drag due to enduring contact in a 
! close-packed system. Note the check for mmax >= 2 below is unnecessary
! since this routine will not be entered unless mmax >=2
            IF(CLOSE_PACKED(M) .AND. CLOSE_PACKED(L) .AND. &
              (MMAX >= 2))  F_SS(IJK,LM) = F_SS(IJK,LM) + &
               SEGREGATION_SLOPE_COEFFICIENT*P_star(IJK) 

         ENDIF  ! end if (.not.wall_at(ijk))

      ENDDO    ! end do (ijk=ijkstart3,ijkend3)
           
      RETURN  
      END SUBROUTINE DRAG_SS 




!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Purpose: This module computes a drag coefficient between two solids C
!           phases where one phase (M) is being modeled as a continuum C
!           and the other (L) as discrete particles.                   C
!           The solids-solids drag coefficient model developed by      C
!           Syamlal, 1987 (see drag_ss routine) is being used for      C
!           these computations.                                        C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

      SUBROUTINE DES_DRAG_SS(M, L, IER) 

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
      USE drag
      USE discretelement
      IMPLICIT NONE

!-----------------------------------------------
! Dummy Arguments
!-----------------------------------------------
! Index of continuum solids phase 
      INTEGER :: M 
! Index of discrete solids 'phase'
      INTEGER :: L
! Error index 
      INTEGER :: IER 
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! Relative velocity between solids phase m and l 
      DOUBLE PRECISION :: VREL 
! Indices 
      INTEGER :: I, IJK, IMJK, IJMK, IJKM 
! Continuous and Discrete solids phases indice, respectively
      INTEGER :: CM, DM
! Cell center value of U_sm, U_sl, V_sm, V_sl, W_sm, W_sl
      DOUBLE PRECISION :: USCM, USCL, VSCM, VSCL, WSCM, WSCL
! Particle diameters of phase M and phase L 
      DOUBLE PRECISION :: D_pm, D_pl 
! Sum of particle diameters 
      DOUBLE PRECISION :: DPSUM 
! Drag coefficient 
      DOUBLE PRECISION :: CONST 
! Solids volume fraction, and void fraction
      DOUBLE PRECISION :: EPS, EPG
! Sum over all phases of ratio volume fraction over particle diameter
      DOUBLE PRECISION :: EPSoDP
! Local calculation of radial distribution function between solids
! phases M and L
      DOUBLE PRECISION :: G_0loc_ML 
!-----------------------------------------------
! External Functions
!-----------------------------------------------
      DOUBLE PRECISION , EXTERNAL :: G_0 
!-----------------------------------------------
! Include statement functions
!-----------------------------------------------      
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'
!-----------------------------------------------

      DO IJK = ijkstart3, ijkend3  

         IF (.NOT.WALL_AT(IJK)) THEN 
            I = I_OF(IJK) 
            IMJK = IM_OF(IJK) 
            IJMK = JM_OF(IJK) 
            IJKM = KM_OF(IJK) 

! Calculate velocity components at i, j, k
! Note no further manipulation is needed for the discrete phase quantities since
! these are already determined at cell centers
            USCL = DES_U_S(IJK,L)
            VSCL = DES_V_S(IJK,L)
            WSCL = DES_W_S(IJK,L)

            USCM = AVG_X_E(U_S(IMJK,M),U_S(IJK,M),I) 
            VSCM = AVG_Y_N(V_S(IJMK,M),V_S(IJK,M)) 
            WSCM = AVG_Z_T(W_S(IJKM,M),W_S(IJK,M)) 

! magnitude of solids-solids relative velocity
            VREL = SQRT((USCL - USCM)**2 + (VSCL - VSCM)**2 + (WSCL - WSCM)**2) 

            D_PM = D_P(IJK,M) 
            D_PL = DES_D_P0(L) 
            DPSUM = D_PL + D_PM        

! taken from G_0.f subroutine (lebowitz form)
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
            G_0loc_ML = ONE/EPg + 3.0d0*EPSoDP*D_pM*D_PL / &
               (EPg*EPg *(D_pM + D_pL))


            CONST = 3.d0*(ONE + C_E)*(PI/2.d0 + C_F*PI*PI/8.d0)*&
               DPSUM**2/(2.d0*PI*(DES_RO_S(L)*D_PL**3+RO_S(M)*D_PM**3)) 

            F_SDS(IJK,M,L) = CONST*ROP_S(IJK,M)*DES_ROP_S(IJK,L)*&
               G_0loc_ML*VREL 

! Gera: account for particle-particle drag due to enduring contact in a 
! close-packed system.
            IF(CLOSE_PACKED(M)) F_SDS(IJK,M,L) = F_SDS(IJK,M,L) + &
               SEGREGATION_SLOPE_COEFFICIENT*P_star(IJK) 

         ENDIF  ! end if (.not.wall_at(ijk))

      ENDDO    ! end do (ijk=ijkstart3,ijkend3)
           
      RETURN  
      END SUBROUTINE DES_DRAG_SS 



