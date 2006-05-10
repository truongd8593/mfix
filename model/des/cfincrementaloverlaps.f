!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFINCREMENTALOVERLAPS(Vno, Vtan, OVERLP_N, OVERLP_T)   C
!  Purpose: DES - calculate incremental overlaps between particles     C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CFINCREMENTALOVERLAPS(Vno, Vtan, OVERLP_N, OVERLP_T)

      USE discretelement        
      IMPLICIT NONE

      DOUBLE PRECISION OVERLP_N, OVERLP_T
      DOUBLE PRECISION Vno, Vtan

!-----------------------------------------------------------------------

      OVERLP_N = Vno * DTSOLID
      OVERLP_T = Vtan * DTSOLID

      RETURN
      END SUBROUTINE CFINCREMENTALOVERLAPS


