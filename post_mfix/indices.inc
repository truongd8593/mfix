CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: indices.inc                                            C
C  Purpose: Common block containing arrays for index computations      C
C                                                                      C
C  Author: M. Syamlal                                 Date: dd-mmm-yy  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Revision Number:                                                    C
C  Purpose:                                                            C
C  Author:                                            Date: dd-mmm-yy  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Literature/Document References: None                                C
C                                                                      C
C  Variables referenced: None                                          C
C  Variables modified: None                                            C
C                                                                      C
C  Local variables: None                                               C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
C
C                      Increments used for index computation
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
C
C                      Store LM index values 
      INTEGER          STORE_LM (DIMENSION_M, DIMENSION_M)
C
C                      Identification of the cell class
      INTEGER          CELL_CLASS (DIMENSION_3)
C
C                      I
      INTEGER          I_OF (DIMENSION_3)
C
C                      J
      INTEGER          J_OF (DIMENSION_3)
C
C                      K
      INTEGER          K_OF (DIMENSION_3)
C
C                      Store (I - 1)'s
      INTEGER          Im1 (DIMENSION_I)
C
C                      Store (I + 1)'s
      INTEGER          Ip1 (DIMENSION_I)
C
C                      Store (J - 1)'s
      INTEGER          Jm1 (DIMENSION_J)
C
C                      Store (J + 1)'s
      INTEGER          Jp1 (DIMENSION_J)
C
C                      Store (K - 1)'s
      INTEGER          Km1 (DIMENSION_K)
C
C                      Store (K + 1)'s
      INTEGER          Kp1 (DIMENSION_K)
C
      COMMON /INDICES/
     & INCREMENT_FOR_n, INCREMENT_FOR_s, INCREMENT_FOR_e,
     & INCREMENT_FOR_w, INCREMENT_FOR_t, INCREMENT_FOR_b,
     & INCREMENT_FOR_im, INCREMENT_FOR_ip, INCREMENT_FOR_jm,
     & INCREMENT_FOR_jp, INCREMENT_FOR_km, INCREMENT_FOR_kp,
     & CELL_CLASS, STORE_LM,
     & Im1, Ip1, Jm1, Jp1, Km1, Kp1,
     & I_OF, J_OF, K_OF
