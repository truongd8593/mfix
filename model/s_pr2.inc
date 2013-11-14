!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: Solids_pr                                              C
!  Purpose: To compute solids pressure and its inverse                 C
!                                                                      C
!  Author: M. Syamlal                                 Date: 17-FEB-93  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number:  1                                                 C
!  Purpose:  allow SI units                                            C
!  Author: S. Dartevelle                              Date: 01-Jul-02  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Comments:                                                           C
!     See S_pr1.inc for the parameters                                 C
!     to_SI is a constant to change from CGS to SI (Ba --> Pa)         C
!        see set_constants.f and constant_mod.f                        C
!        if CGS: to_SI=ONE and                                         C
!        if SI : to_SI=0.1                                             C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

! Solids pressure in plastic-flow stress formulation. MFIX default model

! function P_s(EP_g,Ep_star)
      Neg_H(XXX,YYY) = to_SI*a_ps * (MAX(ZERO, (YYY - XXX )))**b_ps

! inverse of Neg_H.  function EP_g(P_s)
      INV_H(XXX,YYY) = YYY - ( XXX/(to_SI*a_ps) )**(ONE/dble(b_ps))

! Differentiate P_s w.r.t. EP_s.  function dP_s/dEP_s (EP_s)
      dPodEP_s(XXX,YYY) = to_SI*a_ps * dble(b_ps)*(MAX(ZERO, (YYY - (ONE - XXX) )))**(b_ps-1)


! Solids pressure in plastic-flow stress formulation. Jackson's model
! function P_s(EP_g,Ep_star)
      !Neg_H(XXX,YYY) = to_SI*a_ps_jackson * (MAX(ZERO, (YYY - XXX )/(XXX - EP_g_cp)))

! inverse of Neg_H.  function EP_g(P_s)
      !INV_H(XXX,YYY) = (to_SI*a_ps_jackson * YYY + XXX*EP_g_cp)/(to_SI*a_ps_jackson + XXX)

! Differentiate P_s w.r.t. EP_s.  function dP_s/dEP_s (EP_s)       
      !dPodEP_s(XXX,YYY) = to_SI*a_ps_jackson * (YYY - EP_g_cp)/(XXX - EP_g_cp)**2

