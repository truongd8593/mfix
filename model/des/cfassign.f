!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFASSIGN(PARTS)                                        C
!  Purpose: DES - To assign values for particles                       C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CFASSIGN(PARTS)

      Use discretelement
      IMPLICIT NONE

      INTEGER L, PARTS
      DOUBLE PRECISION Dp, pi, we, PVOL

!     
!---------------------------------------------------------------------
!     Assignments
!---------------------------------------------------------------------
!     

      DO L = 1, PARTS
      
         Dp = 2*DES_RADIUS(L) 
         pi = 22.0D0/7.0D0
         PVOL= (4.0D0/3.0D0)*pi*DES_RADIUS(L)**3
         PMASS(L) = PVOL*RO_Sol(L)
         MOI(L) = (PMASS(L)*Dp**2.0D0)/10.0D0 

      END DO
      

      RETURN
      END SUBROUTINE CFASSIGN



