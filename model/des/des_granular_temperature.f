!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: DES_GRANULAR_TEMPERATURE(NC)                           C
!  Purpose: DES - Calculate the DES granular temperature               C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE DES_GRANULAR_TEMPERATURE(NC)

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
      
      INTEGER I, J, K, IJKL, IJK
      INTEGER NT, ITER, M, AVE_STEPS
      INTEGER NC, NCOUNT
      DOUBLE PRECISION P0, P1, P2, P3,  DP1, DP2, DP3, DP4, TEMP
      DOUBLE PRECISION XV(1,1,1) !(IJKMAX2,100,MMAX)
      DOUBLE PRECISION YV(1,1,1) !(IJKMAX2,100,MMAX)
      DOUBLE PRECISION ZV(1,1,1) !(IJKMAX2,100,MMAX) 
      DOUBLE PRECISION AVE_VEL_X(1,1) !(IJKMAX2,MMAX)
      DOUBLE PRECISION AVE_VEL_Y(1,1) !(IJKMAX2,MMAX)
      DOUBLE PRECISION AVE_VEL_Z(1,1) !(IJKMAX2,MMAX)

        INCLUDE 'function.inc'

      AVE_STEPS = FACTOR 
      DO IJK = 1, IJKMAX2
         DO M = 1, MMAX
            XV(IJK,NC,M) = U_S(IJK,M)
            YV(IJK,NC,M) = V_S(IJK,M)
            IF(DIMN.EQ.3) THEN
               ZV(IJK,NC,M) = W_S(IJK,M)
            END IF
         END DO
      END DO
      
      IF(NC.EQ.AVE_STEPS) THEN

         DO IJK = 1, IJKMAX2
            DO M = 1, MMAX
               AVE_VEL_X(IJK,M) = 0.0
               AVE_VEL_Y(IJK,M) = 0.0
               IF(DIMN.EQ.3) THEN
                  AVE_VEL_Z(IJK,M) = 0.0
               END IF
            END DO
         END DO
         
         DO NT = 1, NC
            DO IJK = 1, IJKMAX2
               DO M = 1, MMAX
                  AVE_VEL_X(IJK,M) = AVE_VEL_X(IJK,M) + XV(IJK,NT,M)/FACTOR
                  AVE_VEL_Y(IJK,M) = AVE_VEL_Y(IJK,M) + YV(IJK,NT,M)/FACTOR
                  IF(DIMN.EQ.3) THEN  
                     AVE_VEL_Z(IJK,M) = AVE_VEL_Z(IJK,M) + ZV(IJK,NT,M)/FACTOR  
                  END IF
               END DO
            END DO
         END DO

         DO IJK = 1, IJKMAX2
            DO M = 1, MMAX
               TEMP = 0.0
               TEMP = (AVE_VEL_X(IJK,M)-U_S(IJK,M))**2 
               TEMP = TEMP + (AVE_VEL_Y(IJK,M)-V_S(IJK,M))**2
               IF(DIMN.EQ.3) THEN 
                  TEMP = TEMP + (AVE_VEL_Z(IJK,M)-W_S(IJK,M))**2     
               END IF 
               DES_THETA(IJK,M) = TEMP/3
            END DO
         END DO

         DO NT = 1, NC-1
            DO IJK = 1, IJKMAX2
               DO M = 1,MMAX
                  XV(IJK,NT,M) = XV(IJK,NT+1,M)
                  YV(IJK,NT,M) = YV(IJK,NT+1,M)  
                  IF(DIMN.EQ.3) THEN
                     ZV(IJK,NT,M) = ZV(IJK,NT+1,M)
                  END IF  
               END DO
            END DO
         END DO

      END IF

      RETURN
      END SUBROUTINE DES_GRANULAR_TEMPERATURE

