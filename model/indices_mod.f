!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: indices.inc                                            C
!  Purpose: Common block containing arrays for index computations      C
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
 
 
      MODULE indices
 
 
      Use param
      Use param1
 
 
!
!                      Increments used for index computation
      INTEGER          INCREMENT_FOR_n (MAX_CLASS)
      INTEGER          INCREMENT_FOR_s (MAX_CLASS)
      INTEGER          INCREMENT_FOR_e (MAX_CLASS)
      INTEGER          INCREMENT_FOR_w (MAX_CLASS)
      INTEGER          INCREMENT_FOR_t (MAX_CLASS)
      INTEGER          INCREMENT_FOR_b (MAX_CLASS)
      INTEGER          INCREMENT_FOR_im(MAX_CLASS)
      INTEGER          INCREMENT_FOR_ip(MAX_CLASS)
      INTEGER          INCREMENT_FOR_jm(MAX_CLASS)
      INTEGER          INCREMENT_FOR_jp(MAX_CLASS)
      INTEGER          INCREMENT_FOR_km(MAX_CLASS)
      INTEGER          INCREMENT_FOR_kp(MAX_CLASS)
!
!                      Store LM index values
      INTEGER, DIMENSION(:, :), ALLOCATABLE ::           STORE_LM 
!
!                      Identification of the cell class
      INTEGER, DIMENSION(:), ALLOCATABLE ::           CELL_CLASS 
!
!                      I
      INTEGER, DIMENSION(:), ALLOCATABLE ::           I_OF 
!
!                      J
      INTEGER, DIMENSION(:), ALLOCATABLE ::           J_OF 
!
!                      K
      INTEGER, DIMENSION(:), ALLOCATABLE ::           K_OF 
!
!                      Store (I - 1)'s
      INTEGER, DIMENSION(:), ALLOCATABLE ::           Im1 
!
!                      Store (I + 1)'s
      INTEGER, DIMENSION(:), ALLOCATABLE ::           Ip1 
!
!                      Store (J - 1)'s
      INTEGER, DIMENSION(:), ALLOCATABLE ::           Jm1 
!
!                      Store (J + 1)'s
      INTEGER, DIMENSION(:), ALLOCATABLE ::           Jp1 
!
!                      Store (K - 1)'s
      INTEGER, DIMENSION(:), ALLOCATABLE ::           Km1 
!
!                      Store (K + 1)'s
      INTEGER, DIMENSION(:), ALLOCATABLE ::           Kp1 
!
 
 
!!HPF$ distribute CELL_CLASS(block)
!!HPF$ distribute I_OF(block)
!!HPF$ distribute J_OF(block)
!!HPF$ distribute K_OF(block)

      END MODULE indices
