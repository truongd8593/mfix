CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: drag.inc                                               C
C  Purpose: Common block containing drag arrays                        C
C                                                                      C
C  Author: M. Syamlal                                 Date: 20-MAY-92  C
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
C                      Gas-solids drag
      DOUBLE PRECISION F_gs(DIMENSION_3, DIMENSION_M)
C
C                      Solids-solids drag
      DOUBLE PRECISION F_ss(DIMENSION_3, DIMENSION_LM)
C
      COMMON /DRAG_DP/
     &                  F_gs, F_ss
C
C      COMMON /DRAG_C/
C     &                    
