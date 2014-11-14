!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: LOCATION(L2, XMIN, DX)                                 C
!  Purpose: Find the cell center location in X, Y, or Z direction for  C
!           the given index L2.                                        C
!                                                                      C
!  Author: M. Syamlal                                 Date: 01-SEP-92  C
!  Reviewer: M. Syamlal                               Date: 11-DEC-92  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: None                                          C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables: L                                                  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      DOUBLE PRECISION FUNCTION LOCATION (L2, XMIN, DX)
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Index for which the location is required
      INTEGER          L2
!
!                      Starting location of the coordinate
      DOUBLE PRECISION XMIN
!
!                      Cell sizes (DX, DY, or DZ)
!//EFD Nov/11 avoid using dx(*)
!//      DOUBLE PRECISION DX(*)
      DOUBLE PRECISION DX(0:L2)
!
!  Local variables
!
!                      Index
      INTEGER          L
!-----------------------------------------------
!
      LOCATION = XMIN - HALF*DX(1)
      L = 2
      IF (L2 - 1 > 0) THEN

!//EFD      since indexing of dx starts from 0
!//         using DX(1:(L2-1)) instead of DX(:,L2)
!//         LOCATION = LOCATION + SUM(HALF*(DX(:L2-1)+DX(2:L2)))

         LOCATION = LOCATION + SUM(HALF*(DX(1:(L2-1))+DX(2:L2)))
         L = L2 + 1

      ENDIF
      RETURN
      END FUNCTION LOCATION
