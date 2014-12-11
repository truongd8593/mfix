!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: devol.inc                                              C
!  Purpose: Common block containing data relating to devolatilization  C
!                                                                      C
!  Author: S. Venkatesan                              Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References: MGAS code, Wen, et al. (1982)       C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      MODULE usr

      use param, only: DIMENSION_C

      DOUBLE PRECISION DUMMY_DP

! Sherwood number
      DOUBLE PRECISION, DIMENSION(:, :),    ALLOCATABLE :: N_sh

! Proximate Analysis:
      DOUBLE PRECISION :: PAC  ! Char
      DOUBLE PRECISION :: PAV  ! Volatiles
      DOUBLE PRECISION :: PAM  ! Moisture
      DOUBLE PRECISION :: PAA  ! Ash


      END MODULE usr
