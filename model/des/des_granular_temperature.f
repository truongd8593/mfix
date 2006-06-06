!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: DES_GRANULAR_TEMPERATURE(GTC, FAC)                     C
!  Purpose: DES - Calculate the DES granular temperature               C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE DES_GRANULAR_TEMPERATURE(GTC, FAC)

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
      
      INTEGER I, J, K, IJK
      INTEGER NT, M
      INTEGER GTC, FAC
      DOUBLE PRECISION AVE_VEL_X(DIMENSION_3,MMAX)
      DOUBLE PRECISION AVE_VEL_Y(DIMENSION_3,MMAX)
      DOUBLE PRECISION AVE_VEL_Z(DIMENSION_3,MMAX), TEMP

        INCLUDE 'function.inc'

      IF(GTC.EQ.1) THEN
         AVE_VEL_X(:,:) = ZERO
         AVE_VEL_Y(:,:) = ZERO
         IF(DIMN.EQ.3) THEN
            AVE_VEL_Z(:,:) = ZERO
         END IF
      END IF   
  
      DO IJK = 1, DIMENSION_3
         DO M = 1, MMAX
            AVE_VEL_X(IJK,M) = AVE_VEL_X(IJK,M) + DES_U_s(IJK,M)/FAC
            AVE_VEL_Y(IJK,M) = AVE_VEL_Y(IJK,M) + DES_V_s(IJK,M)/FAC
            IF(DIMN.EQ.3) THEN  
               AVE_VEL_Z(IJK,M) = AVE_VEL_Z(IJK,M) + DES_W_s(IJK,M)/FAC  
            END IF
         END DO
      END DO

      IF(GTC.EQ.FAC) THEN
        DO IJK = 1, DIMENSION_3
            DO M = 1, MMAX
               TEMP = 0.0
               TEMP = (AVE_VEL_X(IJK,M)-DES_U_s(IJK,M))**2 
               TEMP = TEMP + (AVE_VEL_Y(IJK,M)-DES_V_s(IJK,M))**2
               IF(DIMN.EQ.3) THEN 
                  TEMP = TEMP + (AVE_VEL_Z(IJK,M)-DES_W_s(IJK,M))**2     
               END IF 
               DES_THETA(IJK,M) = TEMP/3
            END DO
         END DO

         OPEN (UNIT=17,FILE='des_granular_temp.out',STATUS='REPLACE')
         WRITE(17,*)' '
         WRITE(17,*)'T="',S_TIME,'"'
         DO IJK = IJKSTART3, IJKEND3
            I = I_OF(IJK)
            J = J_OF(IJK)
            K = k_OF(IJK)
            WRITE(17,*) IJK, I, J, K, DES_THETA(IJK,1)
         END DO
      END IF

      RETURN
      END SUBROUTINE DES_GRANULAR_TEMPERATURE

