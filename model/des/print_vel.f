!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: PRINT_VEL(TIM)                                         C
!  Purpose: DES - Print particle velocities as desired                 C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE PRINT_VEL(TIM)

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
      IMPLICIT NONE

      INTEGER L, LL, I, J, K, KK, IJKL
      DOUBLE PRECISION TIM

      PRINT *,'PRINT_VEL', TSTOP, IMAX2, JMAX2
      STOP

!     IF(TIM.GE.(TSTOP)) THEN

      K = IMAX2/2
      PRINT *,' U_G' 
      DO IJKL = 1, IJKMAX2
         I = I_OF(IJKL)
         J = J_OF(IJKL)
         IF(I.EQ.K) THEN
            PRINT *, (J-1)*DY(J), U_G(IJKL)
         END IF
      END DO

      K = JMAX2/2
      PRINT *, ' ' 
      PRINT *,' V_G'
      DO IJKL = 1, IJKMAX2
         I = I_OF(IJKL)
         J = J_OF(IJKL)
         IF(J.EQ.K) THEN
            PRINT *, (I-1)*DX(I), V_G(IJKL)
         END IF
      END DO

!     END IF
      
      RETURN
      END SUBROUTINE PRINT_VEL


