CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: IS.INC                                                 C
C  Purpose: Include file for internal surface specifications           C
C                                                                      C
C  Author: M. Syamlal                                 Date: 10-OCT-92  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Revision Number:                                                    C
C  Purpose:                                                            C
C  Author:                                            Date: dd-mmm-yy  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Literature/Document References:                                     C
C                                                                      C
C  Variables referenced:                                               C
C  Variables modified:                                                 C
C                                                                      C
C  Local variables:                                                    C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
C
C                      x coordinate of the west face of a region where
C                      internal surfaces are specified
      DOUBLE PRECISION IS_X_w (DIMENSION_IS) 
C
C                      x coordinate of the east face of a region where
C                      internal surfaces are specified
      DOUBLE PRECISION IS_X_e (DIMENSION_IS)
C
C                      y coordinate of the south face of a region where
C                      internal surfaces are specified
      DOUBLE PRECISION IS_Y_s (DIMENSION_IS)
C
C                      y coordinate of the north face of a region where
C                      internal surfaces are specified
      DOUBLE PRECISION IS_Y_n (DIMENSION_IS)
C
C                      z coordinate of the bottom face of a region where
C                      internal surfaces are specified
      DOUBLE PRECISION IS_Z_b (DIMENSION_IS)
C
C                      z coordinate of the top face of a region where
C                      internal surfaces are specified
      DOUBLE PRECISION IS_Z_t (DIMENSION_IS)
C
C                      i index of the west face of a region where
C                      internal surfaces are specified
      INTEGER          IS_I_w (DIMENSION_IS) 
C
C                      i index of the east face of a region where
C                      internal surfaces are specified
      INTEGER          IS_I_e (DIMENSION_IS)
C
C                      j index of the south face of a region where
C                      internal surfaces are specified
      INTEGER          IS_J_s (DIMENSION_IS)
C
C                      j index of the north face of a region where
C                      internal surfaces are specified
      INTEGER          IS_J_n (DIMENSION_IS)
C
C                      k index of the bottom face of a region where
C                      internal surfaces are specified
      INTEGER          IS_K_b (DIMENSION_IS)
C
C                      k index of the top face of a region where
C                      internal surfaces are specified
      INTEGER          IS_K_t (DIMENSION_IS)
C
C                      Type of internal surface: IMPERMEABLE - no gas or solids
C                      flow through the surface; SEMIPERMEABLE - only gas flows
C                      through the surface
      CHARACTER*16     IS_TYPE (DIMENSION_IS)
C
C                      Logical variable to determine whether an IS is defined
       LOGICAL         IS_DEFINED (DIMENSION_IS)
C
C                      Character variable with values E, N, and T
C                      to determine the IS plane of a flow cell
       CHARACTER       IS_PLANE (DIMENSION_IS)
C
C                      Permeability coefficients for semipermeable internal
C                      surface: 1- Darcy coefficient
C                               2- Inertial resistance factor
      DOUBLE PRECISION IS_PC (DIMENSION_IS, 2)
C
C                      Solids velocity at the semipermeable surface
      DOUBLE PRECISION IS_VEL_s (DIMENSION_IS, DIMENSION_M)
C
      COMMON / IS_DP /
     &       IS_X_w, IS_X_e, IS_Y_s, IS_Y_n, IS_Z_b, IS_Z_t, IS_PC,
     &       IS_VEL_s
C
      COMMON / IS_I /
     &       IS_I_w, IS_I_e, IS_J_s, IS_J_n, IS_K_b, IS_K_t
C
      COMMON / IS_L /
     &       IS_DEFINED
C
      COMMON / IS_C /
     &       IS_TYPE, IS_PLANE
