!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!  Module name: ps_mod.f                                               C
!                                                                      C
!  Purpose: Common block containing point source data.                 C
!                                                                      C
!  Author: J. Musser                                  Date: 10-Jun-13  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References: None                                C
!                                                                      C
!  Variables referenced: None                                          C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables: None                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      MODULE ps

      LOGICAL :: POINT_SOURCE

      INTEGER, DIMENSION(:), ALLOCATABLE :: POINT_SOURCES

      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: PS_VOLUME

      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: PS_VEL_MAG_G
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: PS_VEL_MAG_S

      END MODULE ps
