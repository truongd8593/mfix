!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: GET_BC_AREA                                            C
!  Purpose: Compute area of boundary surfaces                          C
!                                                                      C
!  Author: M. Syamlal                                 Date: 29-JUL-92  C
!  Reviewer: W. Rogers                                Date: 11-DEC-92  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE GET_BC_AREA 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE geometry
      USE bc
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
! 
!                      BC number 
      INTEGER :: BCV 
!
!                      I, J, and K
      INTEGER ::  I, J, K 
!-----------------------------------------------
!



      DO BCV = 1, DIMENSION_BC 
         IF (BC_DEFINED(BCV)) THEN 
            BC_AREA(BCV) = ZERO 
            IF (BC_PLANE(BCV) == 'W') THEN 
               I = BC_I_W(BCV) 
               DO K = BC_K_B(BCV), BC_K_T(BCV) 
                  J = BC_J_S(BCV) 
                  IF (BC_J_N(BCV) - BC_J_S(BCV) + 1 > 0) THEN 
                     BC_AREA(BCV) = BC_AREA(BCV) + SUM(DY(BC_J_S(BCV):BC_J_N(BCV))*&
                        X_E(I-1)*DZ(K)) 
                     J = BC_J_N(BCV) + 1 
                  ENDIF 
               END DO 
            ELSE IF (BC_PLANE(BCV) == 'E') THEN 
               I = BC_I_W(BCV) 
               DO K = BC_K_B(BCV), BC_K_T(BCV) 
                  J = BC_J_S(BCV) 
                  IF (BC_J_N(BCV) - BC_J_S(BCV) + 1 > 0) THEN 
                     BC_AREA(BCV) = BC_AREA(BCV) + SUM(DY(BC_J_S(BCV):BC_J_N(BCV))*&
                        X_E(I)*DZ(K)) 
                     J = BC_J_N(BCV) + 1 
                  ENDIF 
               END DO 
            ELSE IF (BC_PLANE(BCV)=='S' .OR. BC_PLANE(BCV)=='N') THEN 
               J = BC_J_S(BCV) 
               DO K = BC_K_B(BCV), BC_K_T(BCV) 
                  I = BC_I_W(BCV) 
                  IF (BC_I_E(BCV) - BC_I_W(BCV) + 1 > 0) THEN 
                     BC_AREA(BCV) = BC_AREA(BCV) + SUM(DX(BC_I_W(BCV):BC_I_E(BCV))*&
                        X(BC_I_W(BCV):BC_I_E(BCV))*DZ(K)) 
                     I = BC_I_E(BCV) + 1 
                  ENDIF 
               END DO 
            ELSE IF (BC_PLANE(BCV)=='B' .OR. BC_PLANE(BCV)=='T') THEN 
               K = BC_K_B(BCV) 
               DO J = BC_J_S(BCV), BC_J_N(BCV) 
                  I = BC_I_W(BCV) 
                  IF (BC_I_E(BCV) - BC_I_W(BCV) + 1 > 0) THEN 
                     BC_AREA(BCV) = BC_AREA(BCV) + SUM(DX(BC_I_W(BCV):BC_I_E(BCV))*&
                        DY(J)) 
                     I = BC_I_E(BCV) + 1 
                  ENDIF 
               END DO 
            ENDIF 
         ENDIF 
      END DO 
      RETURN  
      END SUBROUTINE GET_BC_AREA 
