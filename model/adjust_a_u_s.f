!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: ADJUST_A_U_s                                            C
!  Purpose: Handle the special case of the center coefficient in       C
!           U_s momentum eq. becoming zero.                            C
!                                                                      C
!  Author: M. Syamlal                                 Date:  2-AUG-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE ADJUST_A_U_S(A_M, B_M, IER) 

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE matrix 
      USE fldvar
      USE physprop
      USE geometry
      USE run
      USE indices
      USE compar
      USE sendrecv  
      IMPLICIT NONE
!-----------------------------------------------
! Dummy Arguments
!-----------------------------------------------
! Error index
      INTEGER, INTENT(INOUT) :: IER
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3, 0:DIMENSION_M)
!-----------------------------------------------
! Local Variables 
!-----------------------------------------------
! Indices
      INTEGER :: I, IP, IJK, IJKE, IMJK
! Phase index
      INTEGER :: M
!-----------------------------------------------
! Include statement functions
!-----------------------------------------------
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
!-----------------------------------------------

      DO M = 1, MMAX  
         IF(TRIM(KT_TYPE) /= 'GHD' .OR. (TRIM(KT_TYPE) == 'GHD' .AND. M==MMAX)) THEN
            IF (MOMENTUM_X_EQ(M)) THEN 

!!$omp     parallel do private(I, IP, IJK, IJKE, IMJK )

               DO IJK = ijkstart3, ijkend3
                  IF (ABS(A_M(IJK,0,M)) < SMALL_NUMBER) THEN 
                     A_M(IJK,E,M) = ZERO 
                     A_M(IJK,W,M) = ZERO 
                     A_M(IJK,N,M) = ZERO 
                     A_M(IJK,S,M) = ZERO 
                     A_M(IJK,T,M) = ZERO 
                     A_M(IJK,B,M) = ZERO 
                     A_M(IJK,0,M) = -ONE 
                     IF (B_M(IJK,M) < ZERO) THEN 
                        IJKE = EAST_OF(IJK) 
                        IP = IP1(I_OF(IJK)) 
                        IF (ROP_S(IJKE,M) > SMALL_NUMBER) THEN 
                           B_M(IJK,M) = SQRT((-B_M(IJK,M)/&
                              (ROP_S(IJKE,M)*AVG_X_E(ONE,ZERO,IP)*&
                              AYZ_U(IJK)))) 
                        ELSE 
                           B_M(IJK,M) = ZERO 
                        ENDIF 
                     ELSEIF (B_M(IJK,M) > ZERO) THEN 
                        I = I_OF(IJK) 
                        IMJK = IM_OF(IJK) 
                        IF (ROP_S(IJK,M) > SMALL_NUMBER) THEN 
                           B_M(IJK,M) = SQRT(B_M(IJK,M)/(ROP_S(IJK,M)*&
                              AVG_X_E(ZERO,ONE,I)*AYZ_U(IMJK))) 
                        ELSE 
                           B_M(IJK,M) = ZERO 
                        ENDIF 
                     ENDIF 
                  ENDIF    ! end if (abs(a_m(ijk,0,m))<small_number)
               ENDDO    ! end do loop (ijk=ijkstart3,ijkend3)

            ENDIF   ! end if (momentum_x_eq(m))  
         ENDIF   ! end if check for GHD Theory
      ENDDO   ! end do loop (m=1,mmax)

      RETURN  
      END SUBROUTINE ADJUST_A_U_S
       

