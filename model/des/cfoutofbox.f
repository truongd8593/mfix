!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFOUTOFBOS(LL)                                         C
!  Purpose: DES - Check to make sure particles have not travered out   C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CFOUTOFBOX(LL)

      Use discretelement
      IMPLICIT NONE

      INTEGER LL, W, E, B, T, S, N
      DOUBLE PRECISION WT, OMEGA, A, ASINOMEGAT

!     
!---------------------------------------------------------------------
!     Checking if a particle is inside the box
!---------------------------------------------------------------------
!     
      IF(DES_F.NE.0.0) THEN
         OMEGA = 2*22*DES_F/7
         A = DES_GAMMA*9.81/(OMEGA*OMEGA)
      ELSE
         A = 0.0
         OMEGA = 0.0
      END IF

      ASINOMEGAT = A*SIN(OMEGA*CALLED*DTSOLID)
!     ASINOMEGAT = 0.0
      W = 0
      E = 0
      B = 0
      T = 0
      S = 0
      N = 0
      WT = DES_RADIUS(LL) 
!     WT = 0.0
      OUTOFBOX = 0 

      IF(DES_POS_NEW(1,LL).LT.(WX1-WT)) THEN
         W = 1
      ELSE
         W = 0
      END IF

      IF(DES_POS_NEW(1,LL).GT.(EX2+WT)) THEN
         E = 1
      ELSE
         E = 0
      END IF

      IF(DES_POS_NEW(2,LL).LT.(BY1+ASINOMEGAT-WT)) THEN
         B = 1
      ELSE
         B = 0
      END IF

      IF(DES_POS_NEW(2,LL).GT.(TY2+WT)) THEN
         T = 1
      ELSE
         T = 0
      END IF

      IF(DES_POS_NEW(3,LL).LT.(SZ1-WT)) THEN
         S = 1
      ELSE
         S = 0
      END IF

      IF(DES_POS_NEW(3,LL).GT.(NZ2+WT)) THEN
         N = 1
      ELSE
         N = 0
      END IF

      OUTOFBOX = W + 2*E + 4*B + 8*T + 16*S + 32*N

      RETURN
      END SUBROUTINE CFOUTOFBOX



