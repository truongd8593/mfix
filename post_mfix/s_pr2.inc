CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: Solids_pr                                              C
C  Purpose: To compute solids pressure and its inverse                 C
C                                                                      C
C  Author: M. Syamlal                                 Date: 17-FEB-93  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Revision Number:                                                    C
C  Purpose:                                                            C
C  Author:                                            Date: dd-mmm-yy  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Literature/Document References:                                     C
C                                                                      C
C  Variables referenced:                                               C
C  Variables modified:                                                 C
C                                                                      C
C  Local variables:                                                    C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
C
C
C     Solids pressure in plastic-flow stress formulation
       Neg_H(XXX)    = 1.0D25 * (MAX(ZERO, (EP_star - XXX )))**10
       INV_H(XXX)    = EP_star - (XXX/1.0D25)**0.1
       dEPodP_s(XXX) = ((XXX/1.0D25)**9)/1.0D24
       dPodEP_s(XXX) = 1.0D25 * 10. 
     &                  * (MAX(ZERO, (EP_star - (1. - XXX) )))**9
C
C     Solids pressure in plastic-flow stress formulation. Linear model
c       Neg_H(XXX)    = 1.0D7 * (MAX(ZERO, (EP_star - XXX )))
c       INV_H(XXX)    = EP_star - (XXX/1.0D7)
c       dEPodP_s(XXX) = 1.0D-7
c       dPodEP_s(XXX) = 1.0D7
C
C      Incompressible granular medium
c       Neg_H(XXX) = P_star(IJK)
c       INV_H(XXX) = EP_star
c       dEPodP_s(XXX) = ZERO
c       dPodEP_s(XXX) = 1.0D15 
