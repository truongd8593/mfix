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
      INTEGER          Mx 
! 
!                      cell index 
      INTEGER          IJK 
! 
!                      Sum over m (EP_sm/D_pm) 
      DOUBLE PRECISION EPSoDP 
! 
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
!                      other radial distribution functions 
      DOUBLE PRECISION , EXTERNAL :: G_0EP 
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'ep_s2.inc'
!
! Start Lebowitz (1964)
!      EPSoDP = ZERO
!      DO 10 Mx = 1, MMAX
!        EPSoDP = EPSoDP + EP_s(IJK, Mx) / D_p(Mx)
!10    CONTINUE
!      G_0 = ONE / EP_g(IJK)
!     &      + 3.0 * EPSoDP * D_p(M1) * D_p(M2)
!     &      / (EP_g(IJK)*EP_g(IJK) *(D_p(M1) + D_p(M2)))
!
! End Lebowitz (1964)
!
!  Start Carnahan-Starling
!    (Do not use when there are more than one granular phase)
!     This is the form of the radial distribution function
!     proposed by Carnahan & Starling
!
      G_0 = G_0EP(EP_S(IJK,M1)) 
!  End Carnahan-Starling
!
      RETURN  
      END FUNCTION G_0 
!
      DOUBLE PRECISION FUNCTION G_0EP (EPS) 
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
!  Start Carnahan-Starling
      G_0EP = 1D0/(1 - EPS) + 1.5D0*EPS*(1D0/(1 - EPS))**2 + 0.5D0*EPS**2*(1D0/&
         (1 - EPS))**3 
!  End Carnahan-Starling
!
      RETURN  
      END FUNCTION G_0EP 
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
!  Start Carnahan-Starling derivative
!
!     Derivative of (G0) wrt EP_s
      DG_0DNU = 1D0/(1. - EPS)**2 + 1.5D0*(1. + EPS)*(1D0/(1. - EPS))**3 + &
         0.5D0*(EPS**2 + 2.*EPS)*(1D0/(1. - EPS))**4 
!  End Carnahan-Starling derivative
      RETURN  
      END FUNCTION DG_0DNU 
