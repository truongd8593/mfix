!  Note that routines G_0AVG, G_0, and DG_0DNU need to be modified to effect 
!  a change in the radial distribution function g_0.
!  The old routine G_0EP has been replaced with G_0CS, which used
!  only for Carnahan-Starling g_0.
      DOUBLE PRECISION FUNCTION G_0AVG (IJK1, IJK2, DIR, L, M1, M2) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE physprop
      USE fldvar
      USE geometry
      USE indices
      USE compar
      USE visc_s 
      USE constant
      USE run
      USE toleranc
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
! 
! 
!                      cell index
      INTEGER          IJK 
! 
!                      cell indices for neighboring cells (passed variable) 
      INTEGER          IJK1, IJK2 
! 
!                      Direction (X, Y, or Z) (passed variable)
      CHARACTER        DIR 
! 
!                      Direction index (i, j, or k) 
      INTEGER          L 
! 
!                      Solids phase index-1 (passed variable) 
      INTEGER          M1 
! 
!                      Solids phase index-2 (passed variable) 
      INTEGER          M2 
! 
!                      Solids phase index 
      INTEGER          Mx 
! 
!                      average solids volume fraction
      DOUBLE PRECISION EPS
! 
!                      average void fraction
      DOUBLE PRECISION EPg, EPg_STAR_AVG
! 
!                      Sum over m (EP_sm/D_pm) 
      DOUBLE PRECISION EPSoDP   
!
!                      Sum over EP_s
      DOUBLE PRECISION SUM_EPS
!
!                      average number density of solids phase mm
      DOUBLE PRECISION NU_MM
! 
!                      volume of a single particle of solids phase mm
      DOUBLE PRECISION VOLP
!
!                      quantity employed in rdf calculation
      DOUBLE PRECISION XI
!
!                      Average D_P for phase M1 and M2
      DOUBLE PRECISION DP_AVG_M1, DP_AVG_M2, DP_AVG
! 
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
! 
!                      averaging function 
      DOUBLE PRECISION , EXTERNAL :: avg_xyz
! 
!                      radial distribution function 
      DOUBLE PRECISION , EXTERNAL :: G_0, G_0CS
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'ep_s2.inc'

      IF(IJK1 == IJK2)THEN
        G_0AVG = G_0(IJK1, M1, M2) 
      ELSE
!
! Part modified to automatically use Lebowitz if more than one solids
! phase is used. sof June 3rd 2005
!
        IF(MMAX > 1) THEN
!
               IF (TRIM(RDF_TYPE) .EQ. 'LEBOWITZ') THEN
!
! Start Lebowitz (1964)
                    DP_AVG_M1 = HALF*(D_p(IJK1,M1)+D_p(IJK2,M1))
                    DP_AVG_M2 = HALF*(D_p(IJK1,M2)+D_p(IJK2,M2))
                    EPSoDP = ZERO
                    DO Mx = 1, MMAX
                       EPS = AVG_XYZ(EP_s(IJK1, Mx), EP_s(IJK2, Mx), DIR, L)
                       EPSoDP = EPSoDP + 2d0*EPS / (D_p(IJK1,Mx)+D_p(IJK2,Mx))
                    END DO
                    EPg = AVG_XYZ(EP_g(IJK1), EP_g(IJK2), DIR, L)
                    G_0AVG = ONE / EPg                                  &
                            + 3.0d0 * EPSoDP * DP_AVG_M1 * DP_AVG_M2    &
                            / (EPg*EPg *(DP_AVG_M1 + DP_AVG_M2))
! End Lebowitz (1964)
!
               ELSEIF (TRIM(RDF_TYPE) .EQ. 'MODIFIED_LEBOWITZ') THEN
! I didn't check these types of RDF (sof) contact J. Galvin for questions	
!                   Start Modified Lebowitz (1964)
                    DP_AVG_M1 = HALF*(D_p(IJK1,M1)+D_p(IJK2,M1))
                    DP_AVG_M2 = HALF*(D_p(IJK1,M2)+D_p(IJK2,M2))
                    EPSoDP = ZERO
                    SUM_EPS = ZERO
!	
                    DO Mx = 1, MMAX
                         EPS = AVG_XYZ(EP_s(IJK1, Mx), EP_s(IJK2, Mx), DIR, L)
                         DP_AVG = HALF*( D_p(IJK1,Mx)+D_p(IJK2,Mx) )
                         EPSoDP = EPSoDP + EPS/DP_AVG
                         SUM_EPS = SUM_EPS + EPS
                    END DO
!
!                   The following if statement is to prevent G_0 from becoming negative 
!                   when the solids volume fraction approachs a maximum packing.  This 
!                   may occur during non-converged iterations when the local solids volume
!                   fraction exceeds the maximum 
!                   Need to revisit this and attempt van Wachem implementation ?
!
                    EPg_STAR_AVG = AVG_XYZ(EP_star_array(IJK1), EP_star_array(IJK2), DIR, L) 
                    IF (SUM_EPS .GE. (ONE - EPg_STAR_AVG)) THEN
                         SUM_EPS = SUM_EPS - DIL_EP_s
                    ENDIF
!
                    G_0AVG = (ONE/(ONE-SUM_EPS/(ONE-EPg_STAR_AVG))) + (3.0d0) * &
                         ( (DP_AVG_M1*DP_AVG_M2)/(DP_AVG_M1+DP_AVG_M2) )*&
                         EPSoDP
!
!                   End Modified Lebowitz (1964)
!
               ELSEIF (TRIM(RDF_TYPE) .EQ. 'MANSOORI') THEN
!
!                   Start extended Carnahan & Starling see Garzo & Dufty (1999) for
!                   details.  This expression is equivalent to the rdf of
!                   Mansoori et al. (1971)
                    DP_AVG_M1 = HALF*(D_p(IJK1,M1)+D_p(IJK2,M1))
                    DP_AVG_M2 = HALF*(D_p(IJK1,M2)+D_p(IJK2,M2))
                    SUM_EPS = ZERO
                    XI = ZERO

                    DO Mx = 1, MMAX
                         EPS = AVG_XYZ(EP_s(IJK1, Mx), EP_s(IJK2, Mx), DIR, L)
                         SUM_EPS = SUM_EPS + EPS
!		
                         DP_AVG = HALF*( D_p(IJK1,Mx)+D_p(IJK2,Mx) )
                         VOLP = (PI/6.0D0)*DP_AVG**3

                         IF (DP_AVG .LE. ZERO) THEN
                              XI = XI + ZERO
                         ELSE
                              NU_MM = EPS/VOLP
                              XI = XI + NU_MM*DP_AVG*DP_AVG
                         ENDIF
                    ENDDO
!
                    XI = (PI/6.0D0)*XI
!
                    G_0AVG = (ONE/(ONE-SUM_EPS)) + (3.0D0)*&
                         ( (DP_AVG_M1*DP_AVG_M2)/(DP_AVG_M1+DP_AVG_M2) )*&
                         ( XI/((ONE-SUM_EPS)*(ONE-SUM_EPS)) ) + (2.0D0) * &
                         ( (DP_AVG_M1*DP_AVG_M2)/(DP_AVG_M1+DP_AVG_M2) ) * &
                         ( (DP_AVG_M1*DP_AVG_M2)/(DP_AVG_M1+DP_AVG_M2) ) * &
                         ( (XI*XI)/((ONE-SUM_EPS)*(ONE-SUM_EPS)*(ONE-SUM_EPS)) )
!	
!                   End extended Carnahan & Starling 
!
               ELSEIF (TRIM(RDF_TYPE) .EQ. 'MODIFIED_MANSOORI') THEN
!
!                   Start modified Mansoori et al. (1971) see van Wachem et al. (2001)
!                   for details.
                    DP_AVG_M1 = HALF*(D_p(IJK1,M1)+D_p(IJK2,M1))
                    DP_AVG_M2 = HALF*(D_p(IJK1,M2)+D_p(IJK2,M2))
                    SUM_EPS = ZERO
                    XI = ZERO
!
                    DO Mx = 1, MMAX
                         EPS = AVG_XYZ(EP_s(IJK1, Mx), EP_s(IJK2, Mx), DIR, L)
                         SUM_EPS = SUM_EPS + EPS
	
                         DP_AVG = HALF*( D_p(IJK1,Mx)+D_p(IJK2,Mx) )
                         VOLP = (PI/6.0D0)*DP_AVG**3
                         IF (DP_AVG .LE. ZERO) THEN
                              XI = XI + ZERO
                         ELSE
                              NU_MM = EPS/VOLP
                              XI = XI + NU_MM*DP_AVG*DP_AVG
                         ENDIF
                    ENDDO
!
                    XI = (PI/6.0D0)*XI
!
                    EPg_STAR_AVG = AVG_XYZ(EP_star_array(IJK1), EP_star_array(IJK2), DIR, L) 
                    IF (SUM_EPS .GE. (ONE-EPg_STAR_AVG) ) THEN
                         SUM_EPS = SUM_EPS - DIL_EP_s
                    ENDIF
!

                    G_0AVG = (ONE/(ONE-SUM_EPS/(ONE-EPg_STAR_AVG))) + (3.0D0)*&
                         ( (DP_AVG_M1*DP_AVG_M2)/(DP_AVG_M1+DP_AVG_M2) )*&
                         ( XI/((ONE-SUM_EPS/(ONE-EPg_STAR_AVG))*&
                         (ONE-SUM_EPS/(ONE-EPg_STAR_AVG))) ) + (2.0D0) * &
                         ( (DP_AVG_M1*DP_AVG_M2)/(DP_AVG_M1+DP_AVG_M2) ) * &
                         ( (DP_AVG_M1*DP_AVG_M2)/(DP_AVG_M1+DP_AVG_M2) ) * &
                         ( (XI*XI)/((ONE-SUM_EPS/(ONE-EPg_STAR_AVG))*&
                         (ONE-SUM_EPS/(ONE-EPg_STAR_AVG))*&
                         (ONE-SUM_EPS/(ONE-EPg_STAR_AVG))) )
!	
!                   End modified Mansoori et al. (1971)
!	
               ENDIF
!
        ELSE !for mmax == 1 use Carnahan-Starling
!
!  Start Carnahan-Starling
!    (Do not use when there are more than one granular phase)
!     This is the form of the radial distribution function
!     proposed by Carnahan & Starling

          EPS = AVG_XYZ(EP_s(IJK1, M1), EP_s(IJK2, M1), DIR, L)

          G_0AVG = G_0CS(EPS) 
!  End Carnahan-Starling
        ENDIF ! for mmax > 1
      ENDIF

      RETURN  
      END FUNCTION G_0AVG 


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: G_0 (IJK, M1, M2)                                      C
!  Purpose: Calculate radial distribution function at contact for a    C
!           mixture of spheres of different diameters                  C
!                                                                      C
!  Author: M. Syamlal                                 Date: 16-MAR-92  C
!  Reviewer: W. Rogers                                Date: 11-DEC-92  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!    Lebowitz, J.L., "Exact solution of generalized Percus-Yevick      C
!      equation for a mixture of hard spheres," The Physical Review,   C
!      A133, 895-899 (1964).                                           C
!    Iddir, Y.H., "Modeling of the multiphase mixture of particles     C
!	  using the kinetic theory approach," PhD Thesis, Illinois     C
!	  Institute of Technology, Chicago, Illinois, 2004:            C 
!	  chapter 2, equations 2-49 through 2-52.                      C
!    Mansoori et al. (1971)                                            C
!       This RDF expression is equivalent to that cited by Jenkins and C
!	   Mancini (1987) & Garzo and Dufty (1999).                    C
!                                                                      C
!  Variables referenced: EP_g, DP_s, IJK                               C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables: EPSoDP, Mx                                         C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      DOUBLE PRECISION FUNCTION G_0 (IJK, M1, M2) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE physprop
      USE fldvar
      USE geometry
      USE indices
      USE compar 
      USE visc_s
      USE constant
      USE run
      USE toleranc
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
 
! 
!                      Solids phase index-1 (passed variable) 
      INTEGER          M1 
! 
!                      Solids phase index-2 (passed variable) 
      INTEGER          M2 
! 
!                      Local solids phase index 
      INTEGER          Mx, MM
! 
!                      cell index 
      INTEGER          IJK 
! 
!                      average solids volume fraction
      DOUBLE PRECISION EPS
! 
!                      average void fraction
      DOUBLE PRECISION EPg 
! 
!                      Sum over m (EP_sm/D_pm) 
      DOUBLE PRECISION EPSoDP 
!
!                      Sum over EP_s
      DOUBLE PRECISION SUM_EPS
!
!                      number density of solids phase mm
      DOUBLE PRECISION NU_MM
! 
!                      volume of a single particle of solids phase mm
      DOUBLE PRECISION VOLP
!
!                      quantity employed in rdf calculation
      DOUBLE PRECISION XI
! 
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
!                      other radial distribution functions 
      DOUBLE PRECISION , EXTERNAL :: G_0CS 
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'ep_s2.inc'
!
      SUM_EPS = ZERO
      EPg = EP_G(IJK)
      DO MM = 1, MMAX
          EPS = EP_s(IJK, MM)
          SUM_EPS = SUM_EPS + EPS
      END DO
        
      IF(MMAX > 1) THEN
          IF (TRIM(RDF_TYPE) .EQ. 'LEBOWITZ') THEN
!
! Start Lebowitz (1964)
               EPSoDP = ZERO
               DO Mx = 1, MMAX
                  EPS = EP_s(IJK, Mx)
                  EPSoDP = EPSoDP + EPS / D_p(IJK,Mx)
               END DO
               EPg = EP_g(IJK)
               G_0 = ONE / EPg                                      &
                   + 3.0d0 * EPSoDP * D_p(IJK,M1) * D_p(IJK,M2)     &
                   / (EPg*EPg *(D_p(IJK,M1) + D_p(IJK,M2)))
! End Lebowitz (1964)
! 
          ELSEIF (TRIM(RDF_TYPE) .EQ. 'MODIFIED_LEBOWITZ') THEN
!
!              Start modified Lebowitz (1964)
               EPSoDP = ZERO
               SUM_EPS = ZERO
!
               DO MM = 1, MMAX
                    EPS = EP_s(IJK, MM)
                    EPSoDP = EPSoDP + (EPS/D_p(IJK,MM))
                    SUM_EPS = SUM_EPS + EPS
               END DO
!
!              The following if statement is to prevent G_0 from becoming negative 
!              when the solids volume fraction approachs a maximum packing.  This 
!              may occur during non-converged iterations when the local solids volume
!              fraction exceeds the maximum 
!
               IF (SUM_EPS .GE. (ONE-EP_star_array(IJK)) ) THEN
                    SUM_EPS = SUM_EPS - DIL_EP_s
               ENDIF

               G_0 = (ONE/(ONE-SUM_EPS/(ONE-EP_star_array(IJK)) )) + (3.d0) * &
                    ( (D_p(IJK,M1)*D_p(IJK,M2))/(D_p(IJK,M1)+D_p(IJK,M2)) )*&
                    EPSoDP
!
!              End modified Lebowitz (1964)
!
          ELSEIF (TRIM(RDF_TYPE) .EQ. 'MANSOORI') THEN
!
!              Start extended Carnahan & Starling see Garzo & Dufty (1999) for
!              details.  This expression is equivalent to the rdf of
!              Mansoori et al. (1971)
               SUM_EPS = ZERO
               XI = ZERO

               DO MM = 1, MMAX
                    EPS = EP_s(IJK, MM)
                    SUM_EPS = SUM_EPS + EPS
                    VOLP = (PI/6.0D0)*D_P(IJK,MM)**3.0
                    IF (D_P(IJK,MM) .LE. ZERO) THEN
                         XI = XI + ZERO
                    ELSE
                         NU_MM = EPS/VOLP
                         XI = XI + NU_MM*D_P(IJK,MM)*D_P(IJK,MM)
                    ENDIF
               ENDDO
!
               XI = (PI/6.0D0)*XI
               G_0 = (ONE/(ONE-SUM_EPS)) + (3.0D0)*&
                    ( (D_P(IJK,M1)*D_P(IJK,M2))/(D_P(IJK,M1)+D_P(IJK,M2)) )*&
                    ( XI/((ONE-SUM_EPS)*(ONE-SUM_EPS)) ) + (2.0D0) * &
                    ( (D_P(IJK,M1)*D_P(IJK,M2))/(D_P(IJK,M1)+D_P(IJK,M2)) ) * &
                    ( (D_P(IJK,M1)*D_P(IJK,M2))/(D_P(IJK,M1)+D_P(IJK,M2)) ) * &
                    ( (XI*XI)/((ONE-SUM_EPS)*(ONE-SUM_EPS)*(ONE-SUM_EPS)) )
!	
!              End extended Carnahan & Starling 
!
          ELSEIF (TRIM(RDF_TYPE) .EQ. 'MODIFIED_MANSOORI') THEN
!
!              Start modified Mansoori et al. (1971) see van Wachem et al. (2001)
!              for details
               SUM_EPS = ZERO
               XI = ZERO

               DO MM = 1, MMAX
                    EPS = EP_s(IJK, MM)
                    SUM_EPS = SUM_EPS + EPS
                    VOLP = (PI/6.0D0)*D_P(IJK,MM)**3.0
                    IF (D_P(IJK,MM) .LE. ZERO) THEN
                         XI = XI + ZERO
                    ELSE
                         NU_MM = EPS/VOLP
                         XI = XI + NU_MM*D_P(IJK,MM)*D_P(IJK,MM)
                    ENDIF
               ENDDO
!
               XI = (PI/6.0D0)*XI
!
               IF (SUM_EPS .GE. (ONE-EP_star_array(IJK)) ) THEN
                    SUM_EPS = SUM_EPS - DIL_EP_s
               ENDIF
!
               G_0 = (ONE/(ONE-SUM_EPS/(ONE-EP_star_array(IJK)))) + (3.0D0)*&
                    ( (D_P(IJK,M1)*D_P(IJK,M2))/(D_P(IJK,M1)+D_P(IJK,M2)) )*&
                    ( XI/((ONE-SUM_EPS/(ONE-EP_star_array(IJK)))*&
                    (ONE-SUM_EPS/(ONE-EP_star_array(IJK)))) ) + (2.0D0) * &
                    ( (D_P(IJK,M1)*D_P(IJK,M2))/(D_P(IJK,M1)+D_P(IJK,M2)) ) * &
                    ( (D_P(IJK,M1)*D_P(IJK,M2))/(D_P(IJK,M1)+D_P(IJK,M2)) ) * &
                    ( (XI*XI)/((ONE-SUM_EPS/(ONE-EP_star_array(IJK)))*&
                    (ONE-SUM_EPS/(ONE-EP_star_array(IJK)))*&
                    (ONE-SUM_EPS/(ONE-EP_star_array(IJK)))) )
!	
!              End modified Mansoori et al. (1971)
          ENDIF       
      ELSE
!
!  Start Carnahan-Starling
!    (Do not use when there are more than one granular phase)
!     This is the form of the radial distribution function
!     proposed by Carnahan & Starling
!
        G_0 = G_0CS(EP_S(IJK,M1)) 
!  End Carnahan-Starling
!        
      ENDIF ! for mmax > 1
      RETURN  
      END FUNCTION G_0 
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: DG_0DNU (EPs)                                          C
!  Purpose: Calculate derivative of radial distribution function at    C
!           w.r.t granular volume fraction                             C
!                                                                      C
!  Author: K. Agrawal                                 Date: 16-FEB-98  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      DOUBLE PRECISION FUNCTION DG_0DNU (EPS) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE physprop 
      USE fldvar 
      USE geometry 
      USE indices 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      DOUBLE PRECISION EPS 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!-----------------------------------------------
!
!
! Start Lebowitz (1964).  g_0 derivative is not needed for multiparticle
! simulations; so this value is set to zero.
        
      IF(MMAX > 1) THEN
        DG_0DNU = ZERO
! End Lebowitz (1964)
!
      ELSE
!  Start Carnahan-Starling derivative
!
!     Derivative of (G0) wrt EP_s
!        DG_0DNU = ONE/(ONE - EPS)**2 + 1.5D0*(ONE + EPS)*(ONE/(ONE - EPS))**3 + &
!           0.5D0*(EPS**2 + 2.d0*EPS)*(ONE/(ONE - EPS))**4 
!     This simpler formula is less expensive to compute, sof Apr-03-2008
        DG_0DNU = (2.5D0-EPS)/(ONE - EPS)**4
!  End Carnahan-Starling derivative
!
      ENDIF
      RETURN  
      END FUNCTION DG_0DNU 



!-------------- No changes required below ------


      DOUBLE PRECISION FUNCTION G_0CS (EPS) 
!       Carnahan-Starling radial distribution function
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE physprop 
      USE fldvar 
      USE geometry 
      USE indices 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      DOUBLE PRECISION EPS 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!-----------------------------------------------
!
!
!      G_0CS = ONE/(ONE - EPS) + 1.5D0*EPS*(ONE/(ONE - EPS))**2 + 0.5D0*EPS**2*(ONE/&
!         (ONE - EPS))**3 
!     This simpler formula is less expensive to compute, sof Apr-03-2008
      G_0CS = (ONE-0.5D0*EPS)/(ONE - EPS)**3
      RETURN  
      END FUNCTION G_0CS 
      
      DOUBLE PRECISION FUNCTION AVG_XYZ (V1, V2, DIR, L) 
!       Arithmetic average in x, y or z directions
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE geometry 
      USE indices 
      USE compar
      IMPLICIT NONE
! 
!                      cell index
      INTEGER          IJK
! 
!                      Direction (X, Y, or Z) (passed variable)
      CHARACTER        DIR 
! 
!                      Direction index (i, j, or k) 
      INTEGER          L 
! 
!                      The variables to be averaged 
      DOUBLE PRECISION V1, V2
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      DOUBLE PRECISION EPS

      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
!
      IF(DIR == 'X')THEN
  	AVG_XYZ = AVG_X(V1, V2, L)
	  
      ELSEIF(DIR == 'Y')THEN
	AVG_XYZ = AVG_Y(V1, V2, L)
	  
      ELSEIF(DIR == 'Z')THEN
 	AVG_XYZ = AVG_Z(V1, V2, L)
	  
      ELSE
  	CALL WRITE_ERROR('AVG_XYZ', 'Unkown direction', 1)
      ENDIF

      RETURN  
      END FUNCTION AVG_XYZ 
