!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: Solids_pr                                              C
!  Purpose: To compute solids pressure and its inverse                 C
!                                                                      C
!  Author: M. Syamlal                                 Date: 17-FEB-93  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
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
!
!     Solids pressure in plastic-flow stress formulation
       Neg_H(XXX)    = 1.0D25 * (MAX(ZERO, (EP_star - XXX )))**10
       INV_H(XXX)    = EP_star - (XXX/1.0D25)**0.1
       dEPodP_s(XXX) = ((XXX/1.0D25)**9)/1.0D24
       dPodEP_s(XXX) = 1.0D25 * 10. &
                       * (MAX(ZERO, (EP_star - (1. - XXX) )))**9
!
!     Solids pressure in plastic-flow stress formulation. Linear model
!       Neg_H(XXX)    = 1.0D7 * (MAX(ZERO, (EP_star - XXX )))
!       INV_H(XXX)    = EP_star - (XXX/1.0D7)
!       dEPodP_s(XXX) = 1.0D-7
!       dPodEP_s(XXX) = 1.0D7
!
!      Incompressible granular medium
!       Neg_H(XXX) = P_star(IJK)
!       INV_H(XXX) = EP_star
!       dEPodP_s(XXX) = ZERO
!       dPodEP_s(XXX) = 1.0D15 
