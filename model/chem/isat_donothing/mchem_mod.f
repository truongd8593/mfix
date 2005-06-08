!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                         C  
!     Module name: MCHEM                                                  C
!     Purpose: CHEM mod file                                              C
!                                                                         C
!     Author: Nan Xie                                   Date: 02-Aug-04   C
!     Reviewer:                                         Date:             C
!                                                                         C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      
      MODULE MCHEM

      Use param


!
!                      Time for ISAT calculation
      DOUBLE PRECISION TIME_isat  
!
!                      Dimension of ODEs solved in ISAT or DI
      INTEGER          NSpec
!
!                      Source terms for all species
      DOUBLE PRECISION RXN_source_g(DIM_N_g) 
      DOUBLE PRECISION RXN_source_s(DIM_M, DIM_N_s)
!
!                      Controlling parameters for ISATAB
      INTEGER          idtab
      INTEGER          mode
      INTEGER          info(100)
!                      Dimension (50+2+DIM_N_G+DIM_M*DIM_N_s+3*DIM_M)
      DOUBLE PRECISION rinfo(1182)
!
!
!                      Controlling parameters for ISATAB
      INTEGER          ITOL, JT, IOPT
!
!============================================================
!     The part that user must modify
!     user define the tolerances for ODEPACK according ITOL
!
      DOUBLE PRECISION RTOL, ATOL(12)
!============================================================

      END MODULE MCHEM
