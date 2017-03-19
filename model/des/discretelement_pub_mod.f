! -*- f90 -*-
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!   Module name: DISCRETELEMENT                                        C
!   Purpose: DES mod file                                              C
!            Common Block containing DEM conditions                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

MODULE DISCRETELEMENT_PUB

!-----------------------------------------------
! Modules
!-----------------------------------------------
   USE param, only: dim_m
   IMPLICIT NONE


! DES - Continuum
   LOGICAL DISCRETE_ELEMENT
   LOGICAL DES_CONTINUUM_COUPLED


END MODULE DISCRETELEMENT_PUB
