!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: DES_GRANULAR_TEMPERATURE
!>  Purpose: DES - Calculate the DES granular temperature               
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
!-----------------------------------------------
! Local Variables
!-----------------------------------------------      
      INTEGER I, J, K, IJK
      INTEGER M, NP, IPART, NPG
      DOUBLE PRECISION  TEMP
!-----------------------------------------------      
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'fun_avg2.inc'
      
      DO IJK = IJKSTART3, IJKEND3
         IF(FLUID_AT(IJK)) THEN
            DO M = 1, MMAX
               AVE_VEL_X(IJK,M) = DES_U_s(IJK,M)
               AVE_VEL_Y(IJK,M) = DES_V_s(IJK,M)
               IF(DIMN.EQ.3) THEN  
                  AVE_VEL_Z(IJK,M) = DES_W_s(IJK,M)
               ENDIF
            ENDDO
         ENDIF
      ENDDO

      DO IJK = ijkstart3, ijkend3 
         IF(FLUID_AT(IJK)) THEN
            I = I_OF(IJK)
            J = J_OF(IJK)
            K = K_OF(IJK)
            
            IF (ASSOCIATED(PIC(I,J,K)%p)) THEN
               NPG = SIZE(PIC(I,J,K)%p)
            ELSE
               NPG = 0
            ENDIF
            
            TEMP = ZERO
            DO IPART = 1, NPG 
               NP = PIC(I,J,K)%p(IPART)
               M = PIJK(NP,5)
               
               TEMP = TEMP + (DES_VEL_NEW(NP,1)-DES_U_s(IJK,M))**2 
               TEMP = TEMP + (DES_VEL_NEW(NP,2)-DES_V_s(IJK,M))**2
               IF(DIMN.EQ.3) THEN 
                  TEMP = TEMP + (DES_VEL_NEW(NP,3)-DES_W_s(IJK,M))**2 
               ENDIF
            ENDDO
            IF(NPG>0)  DES_THETA(IJK,M) = TEMP/(3.0d0 * DFLOAT(NPG))
         ENDIF
         
      ENDDO


!         OPEN (UNIT=17,FILE='des_granular_temp.out',STATUS='REPLACE')
!         WRITE(17,*)' '
!         WRITE(17,*)'T="',S_TIME,'"'
!         DO IJK = IJKSTART3, IJKEND3
!           IF(FLUID_AT(IJK)) THEN
!             I = I_OF(IJK)
!             J = J_OF(IJK)
!             K = k_OF(IJK)
!             WRITE(17,*) IJK, I, J, K, DES_THETA(IJK,1)
!           ENDIF
!         END DO

      RETURN
      END SUBROUTINE DES_GRANULAR_TEMPERATURE

