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
      USE compar
      USE fldvar
      USE indices
      USE boundfunijk
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
!//09/02/99 - Modified calls to BOUND_FUNIJK to have a self consistent formulation - Sreekanth
!
      IMJK = UNDEFINED_I
      IPJK = UNDEFINED_I
      IJMK = UNDEFINED_I
      IJPK = UNDEFINED_I
      IJKM = UNDEFINED_I
      IJKP = UNDEFINED_I
      IJKW = UNDEFINED_I
      IJKE = UNDEFINED_I
      IJKS = UNDEFINED_I
      IJKN = UNDEFINED_I
      IJKB = UNDEFINED_I
      IJKT = UNDEFINED_I


      IF(IM1(I).NE.UNDEFINED_I) THEN
        IMJK = BOUND_FUNIJK(IM1(I),J,K) 
!
!  IJKW
!
        IF (WALL_AT(IMJK).and..not.cyclic_at(IMJK)) THEN 
           IJKW = IJK 
        ELSE 
           IJKW = IMJK 
        ENDIF 
      ENDIF

      IF(IP1(I).NE.UNDEFINED_I) THEN
        IPJK = BOUND_FUNIJK(IP1(I),J,K) 
!
!  IJKE
!
        IF (WALL_AT(IPJK).and..not.cyclic_at(IPJK)) THEN 
           IJKE = IJK 
        ELSE 
           IJKE = IPJK 
        ENDIF 
      ENDIF

      IF(JM1(J).NE.UNDEFINED_I) THEN
        IJMK = BOUND_FUNIJK(I,JM1(J),K) 
!
!  IJKS
!
        IF (WALL_AT(IJMK).and..not.cyclic_at(IJMK)) THEN 
           IJKS = IJK 
        ELSE 
           IJKS = IJMK 
        ENDIF 
      ENDIF

      IF(JP1(J).NE.UNDEFINED_I) THEN
        IJPK = BOUND_FUNIJK(I,JP1(J),K) 
!
!  IJKN
!
        IF (WALL_AT(IJPK).and..not.cyclic_at(IJPK)) THEN 
           IJKN = IJK 
        ELSE 
           IJKN = IJPK 
        ENDIF 
      ENDIF

      IF(KM1(K).NE.UNDEFINED_I) THEN
        IJKM = BOUND_FUNIJK(I,J,KM1(K)) 
!
!  IJKB
!
        IF (WALL_AT(IJKM).and..not.cyclic_at(IJKM)) THEN 
           IJKB = IJK 
        ELSE 
           IJKB = IJKM 
        ENDIF 
      ENDIF

      IF(KP1(K).NE.UNDEFINED_I) THEN
        IJKP = BOUND_FUNIJK(I,J,KP1(K)) 
!
!  IJKT
!
        IF (WALL_AT(IJKP).and..not.cyclic_at(IJKP)) THEN 
           IJKT = IJK 
        ELSE 
           IJKT = IJKP 
        ENDIF 
      ENDIF
!
      RETURN  
      END SUBROUTINE SET_INDEX1A 

