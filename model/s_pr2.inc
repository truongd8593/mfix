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
!     See S_pr1.inc for the parameters
!     Solids pressure in plastic-flow stress formulation. Mfix default model
       ! function P_s(EP_g)
        Neg_H(XXX)    = a_ps * (MAX(ZERO, (EP_star - XXX )))**b_ps
       
       ! inverse of Neg_H.  function EP_g(P_s)
        INV_H(XXX)    = EP_star - (XXX/a_ps)**(1./dble(b_ps))

       ! Differentiate P_s w.r.t. EP_s.  function dP_s/dEP_s (EP_s)       
        dPodEP_s(XXX) = a_ps * dble(b_ps) &
                        * (MAX(ZERO, (EP_star - (1. - XXX) )))**(b_ps-1)

!
!     Solids pressure in plastic-flow stress formulation. Jackson's model
       ! function P_s(EP_g)
!        Neg_H(XXX)    = a_ps_jackson * (MAX(ZERO, (EP_star - XXX )/(XXX - EP_g_cp)))
       
       ! inverse of Neg_H.  function EP_g(P_s)
!        INV_H(XXX)    = (a_ps_jackson * EP_star + XXX*EP_g_cp)/(a_ps_jackson + XXX)

       ! Differentiate P_s w.r.t. EP_s.  function dP_s/dEP_s (EP_s)       
!        dPodEP_s(XXX) = a_ps_jackson * (EP_star - EP_g_cp)/(XXX - EP_g_cp)**2
