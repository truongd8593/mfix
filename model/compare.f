!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: COMPARE(V1, V2)                                        C
!  Purpose: Return .TRUE. if values V1 and V2 are nearly equal         C
!                                                                      C
!  Author: M. Syamlal                                 Date: 29-JUL-92  C
!  Reviewer: W. Rogers                                Date: 11-DEC-92  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      LOGICAL FUNCTION COMPARE (V1, V2) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE toleranc 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!                      Values to be compared
      DOUBLE PRECISION V1, V2 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!-----------------------------------------------
!
!
      IF (ABS(V1) <= SMALL_NUMBER) THEN 
         IF (ABS(V2) <= SMALL_NUMBER) THEN 
            COMPARE = .TRUE. 
         ELSE 
            COMPARE = .FALSE. 
         ENDIF 
      ELSE 
         IF (ABS(V2/V1 - ONE) <= TOL_COM) THEN 
            COMPARE = .TRUE. 
         ELSE 
            COMPARE = .FALSE. 
         ENDIF 
      ENDIF 
      RETURN  
      END FUNCTION COMPARE 
!
!
      LOGICAL FUNCTION IS_SMALL (V, TOL) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE toleranc 
      USE geometry
      USE indices
      USE compar      !//d
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Tolerance value for small
      DOUBLE PRECISION TOL
!
!                      Field vriable array
      DOUBLE PRECISION, DIMENSION(DIMENSION_3) :: V 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IJK
!-----------------------------------------------
      INCLUDE 'function.inc'
!
      IS_SMALL = .FALSE. 
      DO IJK = IJKSTART3, IJKEND3
         IF (FLUID_AT(IJK)) THEN 
            IF (ABS(V(IJK)) > TOL) RETURN  
         ENDIF 
      END DO 
      IS_SMALL = .TRUE. 
!
      RETURN  
      END FUNCTION IS_SMALL 
