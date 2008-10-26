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

      SUBROUTINE DES_GRANULAR_TEMPERATURE

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
      USE sendrecv

      IMPLICIT NONE
      
      INTEGER I, J, K, IJK
      INTEGER NT, M, NP, IPART, NPG
      INTEGER GTC, FAC
      DOUBLE PRECISION  TEMP
      
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'fun_avg2.inc'
      
      DO IJK = 1, DIMENSION_3
         DO M = 1, MMAX
            AVE_VEL_X(IJK,M) = DES_U_s(IJK,M)
            AVE_VEL_Y(IJK,M) = DES_V_s(IJK,M)
            IF(DIMN.EQ.3) THEN  
               AVE_VEL_Z(IJK,M) = DES_W_s(IJK,M)
            END IF
         END DO
      END DO
      
      DO IJK = 1, DIMENSION_3
         
         IF(FLUID_AT(IJK)) THEN
            i = i_of(ijk)
            j = j_of(ijk)
            k = k_of(ijk)
            
            If (ASSOCIATED(PIC(I,J,K)%p)) then
               NPG = SIZE(PIC(I,J,K)%p)
            Else
               NPG = 0
            Endif
            
            DO IPART = 1, NPG 
               NP = PIC(I,J,K)%p(IPART)
               M = PIJK(NP,5)
               
               TEMP = ZERO
               TEMP = (DES_VEL_NEW(NP,1)-DES_U_s(IJK,M))**2 
               TEMP = TEMP + (DES_VEL_NEW(NP,2)-DES_V_s(IJK,M))**2
               IF(DIMN.EQ.3) THEN 
                  TEMP = TEMP + (DES_VEL_NEW(NP,3)-DES_W_s(IJK,M))**2     
               END IF
               DES_THETA(IJK,M) = TEMP/3.0d0
            end DO
         end IF
         
      end DO


!         OPEN (UNIT=17,FILE='des_granular_temp.out',STATUS='REPLACE')
!         WRITE(17,*)' '
!         WRITE(17,*)'T="',S_TIME,'"'
!         DO IJK = IJKSTART3, IJKEND3
!            I = I_OF(IJK)
!            J = J_OF(IJK)
!            K = k_OF(IJK)
!            WRITE(17,*) IJK, I, J, K, DES_THETA(IJK,1)
!         END DO

      RETURN
      END SUBROUTINE DES_GRANULAR_TEMPERATURE

