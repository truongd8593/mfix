!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: PRESSURE_DROP(PARTS)                                   C
!  Purpose: DES - Calculte the pressure drop across various points in  C
!           the domain                                                 C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE PRESSURE_DROP(PARTS)

      USE discretelement
      USE param
      USE param1
      USE parallel
      USE fldvar
      USE run
      USE geometry
      USE matrix
      USE indices
      USE physprop
      USE drag
      USE constant
      USE compar
      IMPLICIT NONE      

      INTEGER J0, J1, J2, J3
      INTEGER I, J, K, IJKL, PARTS, IJK
      DOUBLE PRECISION P0, P1, P2, P3,  DP1, DP2, DP3, DP4

      INCLUDE 'function.inc'

!     PRINT *,'PRESSURE DROP'

      J0 = 2
      J1 = 10
      J2 = 11
      J3 = 44

      P0 = 0.0
      P1 = 0.0
      P2 = 0.0
      P3 = 0.0

      DO IJKL = 1, IJKMAX2

         I = I_OF(IJKL)
         J = J_OF(IJKL)
         K = K_OF(IJKL)

         IF((I.GT.1).AND.(I.LT.IMAX2)) THEN
            IF(J.EQ.J0) THEN
               P0 = P0 + P_G(IJKL)/IMAX
            ELSE IF(J.EQ.J1) THEN
               P1 = P1 + P_G(IJKL)/IMAX 
            ELSE IF(J.EQ.J2) THEN
               P2 = P2 + P_G(IJKL)/IMAX 
            ELSE IF(J.EQ.J3) THEN
               P3 = P3 + P_G(IJKL)/IMAX
            END IF
         END IF

         DP1 = P1 - P0
         DP2 = P2 - P0 
         DP3 = P3 - P0
         DP4 = P3 - P1	
         
      END DO

      OPEN (UNIT=13, FILE='des_pressure.out', STATUS='REPLACE')
      WRITE(13,*) TIME, P0, P1, P2, P3	
      OPEN (UNIT=15, FILE='des_pressure_drop.out', STATUS='REPLACE')
      WRITE(15,*) TIME, DP1, DP2, DP3, DP4	
      OPEN (UNIT=17, FILE='des_granular_temp.out', STATUS='REPLACE')
      WRITE(17,*)' '
      WRITE(17,*)'ZONE T="',TIME,'"'
      DO IJKL = 1, IJKMAX2
         I = I_OF(IJKL)
         J = J_OF(IJKL)
         K = K_OF(IJKL)
         WRITE(17,*) I, J, DES_THETA(IJKL,1)
      END DO

      RETURN
      END SUBROUTINE PRESSURE_DROP


