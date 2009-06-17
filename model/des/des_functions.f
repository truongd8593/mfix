!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: DES_DOTPRDCT                                           C
!  Purpose: Calculate Dot Product                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      DOUBLE PRECISION FUNCTION DES_DOTPRDCT(XX,YY) 
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1  
      USE discretelement
      IMPLICIT NONE
 
      INTEGER II
      DOUBLE PRECISION DOTP, XX(DIMN), YY(DIMN) 
 
      DOTP = ZERO

      DO II = 1, DIMN
         DOTP = DOTP + XX(II)*YY(II)
      ENDDO
      DES_DOTPRDCT = DOTP

      RETURN  
      END FUNCTION DES_DOTPRDCT 


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: DES_CROSSPRDCT                                         C
!  Purpose: Calculate Dot Product                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE DES_CROSSPRDCT (AA, XX,YY) 
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1  
      USE discretelement
      IMPLICIT NONE
 
      DOUBLE PRECISION AA(DIMN), XX(DIMN), YY(DIMN) 
 
      IF(DIMN.EQ.3) THEN
         AA(1) = XX(2)*YY(3) - XX(3)*YY(2) 
         AA(2) = XX(3)*YY(1) - XX(1)*YY(3) 
         AA(3) = XX(1)*YY(2) - XX(2)*YY(1)
      ELSE
         AA(1) = - XX(1)*YY(2) 
         AA(2) = XX(1)*YY(1)  
      ENDIF

      RETURN  
      END SUBROUTINE DES_CROSSPRDCT


