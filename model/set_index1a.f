!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SET_INDEX1A(I, J, K, IJK, IMJK, IPJK, IJMK, IJPK,      C
!                            IJKM, IJKP, IJKW, IJKE, IJKS, IJKN,       C
!                            IJKB, IJKT)                               C
!  Purpose: Set the indices of the neighbors of cell ijk (brute force) C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-JAN-92  C
!  Reviewer:M. Syamlal, S. Venkatesan, P. Nicoletti,  Date: 29-JAN-92  C
!           W. Rogers                                                  C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: Modify index computations for K for setting periodic       C
!           boundary conditions in a cylindrical geometry where z goes C
!           from 0 to 2 pi                                             C
!  Author: M. Syamlal                                 Date: 10-MAR-92  C
!  Revision Number: 2                                                  C
!  Purpose:  Calculate only the nearest neighbor indices.( for code    C
!            optimization)                                             C
!  Author: M. Syamlal                                 Date: 23-SEP-92  C
!  Reviewer: M. Syamlal                               Date: 11-DEC-92  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: I, J, K, IJK                                  C
!                                                                      C
!  Variables modified: IJKM, IJMK, IMJK, IPJK, IJPK, IJKP, IJKW, IJKE, C
!                      IJKS, IJKN, IJKB, IJKT                          C
!                                                                      C
!  Local variables: None                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE SET_INDEX1A(I, J, K, IJK, IMJK, IPJK, IJMK, IJPK, IJKM, IJKP, &
         IJKW, IJKE, IJKS, IJKN, IJKB, IJKT) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE physprop
      USE geometry
      USE constant
      USE fldvar
      USE indices
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER I, J, K, IJK, IMJK, IPJK, IJMK, IJPK, IJKM, IJKP, IJKW, IJKE, &
         IJKS, IJKN, IJKB, IJKT 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      LOGICAL :: COMPARE 
!-----------------------------------------------
      INCLUDE 'function.inc'
!
!
      IMJK = BOUND_FUNIJK(I - 1,J,K) 
      IPJK = BOUND_FUNIJK(I + 1,J,K) 
      IJMK = BOUND_FUNIJK(I,J - 1,K) 
      IJPK = BOUND_FUNIJK(I,J + 1,K) 
      IJKM = BOUND_FUNIJK(I,J,K - 1) 
      IJKP = BOUND_FUNIJK(I,J,K + 1) 
!
! Modify indices as needed to allow cyclic boundary conditions
!
      IF (CYCLIC_X) THEN 
         IF (I == IMAX1) IPJK = BOUND_FUNIJK(IMIN1,J,K) 
         IF (I == IMIN1) IMJK = BOUND_FUNIJK(IMAX1,J,K) 
      ENDIF 
      IF (CYCLIC_Y) THEN 
         IF (J == JMAX1) IJPK = BOUND_FUNIJK(I,JMIN1,K) 
         IF (J == JMIN1) IJMK = BOUND_FUNIJK(I,JMAX1,K) 
      ENDIF 
      IF (CYCLIC_Z) THEN 
         IF (K == KMAX1) IJKP = BOUND_FUNIJK(I,J,KMIN1) 
         IF (K == KMIN1) IJKM = BOUND_FUNIJK(I,J,KMAX1) 
      ENDIF 
!
!  IJKW
!
      IF (WALL_AT(IMJK)) THEN 
         IJKW = IJK 
      ELSE 
         IJKW = IMJK 
      ENDIF 
!
!  IJKE
!
      IF (WALL_AT(IPJK)) THEN 
         IJKE = IJK 
      ELSE 
         IJKE = IPJK 
      ENDIF 
!
!  IJKS
!
      IF (WALL_AT(IJMK)) THEN 
         IJKS = IJK 
      ELSE 
         IJKS = IJMK 
      ENDIF 
!
!  IJKN
!
      IF (WALL_AT(IJPK)) THEN 
         IJKN = IJK 
      ELSE 
         IJKN = IJPK 
      ENDIF 
!
!  IJKB
!
      IF (WALL_AT(IJKM)) THEN 
         IJKB = IJK 
      ELSE 
         IJKB = IJKM 
      ENDIF 
!
!  IJKT
!
      IF (WALL_AT(IJKP)) THEN 
         IJKT = IJK 
      ELSE 
         IJKT = IJKP 
      ENDIF 
!
      RETURN  
      END SUBROUTINE SET_INDEX1A 
