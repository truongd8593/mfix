CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: G_0 (IJK, M1, M2)                                      C
C  Purpose: Calculate radial distribution function at contact for a    C
C           mixture of spheres of different diameters                  C
C                                                                      C
C  Author: M. Syamlal                                 Date: 16-MAR-92  C
C  Reviewer: W. Rogers                                Date: 11-DEC-92  C
C                                                                      C
C  Revision Number:                                                    C
C  Purpose:                                                            C
C  Author:                                            Date: dd-mmm-yy  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Literature/Document References:                                     C
C    Lebowitz, J.L., "Exact solution of generalized Percus-Yevick      C
C      equation for a mixture of hard spheres," The Physical Review,   C
C      A133, 895-899 (1964).                                           C
C  Variables referenced: EP_g, DP_s, IJK                               C
C  Variables modified: None                                            C
C                                                                      C
C  Local variables: EPSoDP, Mx                                         C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
C
      DOUBLE PRECISION FUNCTION G_0(IJK, M1, M2)
C
      IMPLICIT NONE
C
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE 'physprop.inc'
      INCLUDE 'fldvar.inc'
      INCLUDE 'geometry.inc'
      INCLUDE 'indices.inc'
C
C  Local variables
C
C                      Solids phase index-1 (passed variable)
      INTEGER          M1
C
C                      Solids phase index-2 (passed variable)
      INTEGER          M2
C
C                      Local solids phase index
      INTEGER          Mx
C
C                      cell index
      INTEGER          IJK
C
C                      Sum over m (EP_sm/D_pm)
      DOUBLE PRECISION EPSoDP
C
C                      other radial distribution functions
      DOUBLE PRECISION G_0Ep
C
C  Statement functions
C
      INCLUDE 'ep_s1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'ep_s2.inc'

C Start Lebowitz (1964)
c      EPSoDP = ZERO
c      DO 10 Mx = 1, MMAX
c        EPSoDP = EPSoDP + EP_s(IJK, Mx) / D_p(Mx)
c10    CONTINUE
c      G_0 = ONE / EP_g(IJK) 
c     &      + 3.0 * EPSoDP * D_p(M1) * D_p(M2)
c     &      / (EP_g(IJK)*EP_g(IJK) *(D_p(M1) + D_p(M2)))
C
C End Lebowitz (1964)

C  Start Carnahan-Starling 
C    (Do not use when there are more than one granular phase)
C     This is the form of the radial distribution function
C     proposed by Carnahan & Starling

      G_0 = G_0Ep(EP_s(IJK,M1)) 
C  End Carnahan-Starling 

      RETURN
      END

      DOUBLE PRECISION FUNCTION G_0Ep(EPs)
C
      IMPLICIT NONE
C
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
C
C                      solids volume fraction
      DOUBLE PRECISION EPs


C  Start Carnahan-Starling 
      G_0Ep =   1d0/(1-EPs)
     &      + 1.5d0*EPs* ((1d0/(1-EPs))**2)
     &      + 0.5d0*(EPs**2)*((1d0/(1-EPs))**3)   
C  End Carnahan-Starling 

      RETURN
      END

CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: DG_0DNU (EPs)                                          C
C  Purpose: Calculate derivative of radial distribution function at    C
C           w.r.t granular volume fraction                             C
C                                                                      C
C  Author: K. Agrawal                                 Date: 16-FEB-98  C
C  Reviewer:                                          Date:            C
C                                                                      C
C  Revision Number:                                                    C
C  Purpose:                                                            C
C  Author:                                            Date: dd-mmm-yy  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Literature/Document References:                                     C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
C
      DOUBLE PRECISION FUNCTION DG_0DNU(EPs)
C
      IMPLICIT NONE
C
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
C
C                      solids volume fraction
      DOUBLE PRECISION EPs
C
C  Start Carnahan-Starling derivative
C
C     Derivative of (G0) wrt EP_s
      DG_0DNU = 1d0/( (1.-EPs)**2 )
     &         + 1.5d0 * (1.+EPs) * ( (1d0/(1.-EPs))**3 )
     &         + 0.5d0 * (EPs**2 + 2.*EPs) * ( 1d0/(1.-EPs) )**4
C  End Carnahan-Starling derivative
      RETURN
      END 
