!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: mflux_mod                                              C
!  Purpose: Module for mass fluxes and densities at faces               C
!                                                                      C
!  Author: M. Syamlal                                 Date: dd-mmm-yy  C
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
 
 
      MODULE mflux
 
 
      Use param
      Use param1
 
 
!
!                      x-component of gas mass flux
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  Flux_gE 
!
!                      x-component of solids mass flux
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  Flux_sE 
!
!                      y-component of gas mass flux
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  Flux_gN 
!
!                      y-component of solids mass flux
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  Flux_sN 
!
!                      z-component of gas mass flux
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  Flux_gT 
!
!                      z-component of solids mass flux
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  Flux_sT 
!
!                      Added mass Flux Components to be used for scalar eq.
!                      Note: added mass apply only to one solids phase (M=M_AM)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  Flux_gSE, Flux_sSE, Flux_gSN
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  Flux_sSN, Flux_gST, Flux_sST
!
!
!                      macroscopic gas density at east face
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  ROP_gE 
!
!                      macroscopic solids density at east face
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  ROP_sE 
!
!                      macroscopic gas density at north face
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  ROP_gN 
!
!                      macroscopic solids density at north face
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  ROP_sN 
!
!                      macroscopic gas density at top face
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  ROP_gT 
!
!                      macroscopic solids density at top face
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  ROP_sT 
 
 
! for GHD Theory
!                      x-component of solids total number density flux
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  Flux_nE 
!                      y-component of solids total number density flux
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  Flux_nN 
!                      z-component of solids total number density flux
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  Flux_nT 
 
! end GHD Theory modification
 

      END MODULE mflux
