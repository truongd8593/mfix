!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: ZERO_NORM_VEL                                          C
!  Purpose: Set the velocity component normal to a wall to zero        C
!                                                                      C
!  Author: M. Syamlal                                 Date: 14-MAY-92  C
!  Reviewer: S. Venkatesan                            Date: 11-DEC-92  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: IJK                                           C
!  Variables modified: U_g, U_s, V_g, V_s, W_g, W_s                    C
!                                                                      C
!  Local variables: None                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE ZERO_NORM_VEL 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE geometry
      USE physprop
      USE fldvar
      USE indices
      USE is
      USE compar        !//d
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
      INTEGER          ISV
!      
!                      Indicies
      INTEGER          IJK,  IMJK, IJMK, IJKM, M
!
!-----------------------------------------------
      INCLUDE 'function.inc'
!
!!$omp  parallel do private( IMJK, IJMK, IJKM)
      DO IJK = 1, IJKMAX2 
         IF (.NOT.WALL_AT(IJK)) THEN 
            IF (IP_AT_E(IJK)) U_G(IJK) = ZERO 
            IF (IP_AT_N(IJK)) V_G(IJK) = ZERO 
            IF (IP_AT_T(IJK)) W_G(IJK) = ZERO 
         ELSE 
!-----------------------------------------------------------------------
            IMJK = IM_OF(IJK) 
            IJMK = JM_OF(IJK) 
            IJKM = KM_OF(IJK) 
!-----------------------------------------------------------------------
            U_G(IJK) = ZERO 
            V_G(IJK) = ZERO 
            W_G(IJK) = ZERO 
            IF (.NOT.(CYCLIC_AT(IJK) .AND. I_OF(IJK)==IMAX2)) U_G(IMJK) = ZERO 
            IF (.NOT.(CYCLIC_AT(IJK) .AND. J_OF(IJK)==JMAX2)) V_G(IJMK) = ZERO 
            IF (.NOT.(CYCLIC_AT(IJK) .AND. K_OF(IJK)==KMAX2)) W_G(IJKM) = ZERO 
         ENDIF 
      END DO 
      DO M = 1, MMAX 
!
!!$omp  parallel do private( ISV,  IMJK, IJMK, IJKM)
         DO IJK = 1, IJKMAX2 
            IF (.NOT.WALL_AT(IJK)) THEN 
               IF (IP_AT_E(IJK)) THEN 
                  U_S(IJK,M) = ZERO 
               ELSE IF (SIP_AT_E(IJK)) THEN 
                  ISV = FLAG_E(IJK) - 1000 
                  U_S(IJK,M) = IS_VEL_S(ISV,M) 
               ENDIF 
               IF (IP_AT_N(IJK)) THEN 
                  V_S(IJK,M) = ZERO 
               ELSE IF (SIP_AT_N(IJK)) THEN 
                  ISV = FLAG_N(IJK) - 1000 
                  V_S(IJK,M) = IS_VEL_S(ISV,M) 
               ENDIF 
               IF (IP_AT_T(IJK)) THEN 
                  W_S(IJK,M) = ZERO 
               ELSE IF (SIP_AT_T(IJK)) THEN 
                  ISV = FLAG_T(IJK) - 1000 
                  W_S(IJK,M) = IS_VEL_S(ISV,M) 
               ENDIF 
            ELSE 
!-----------------------------------------------------------------------
               IMJK = IM_OF(IJK) 
               IJMK = JM_OF(IJK) 
               IJKM = KM_OF(IJK) 
!-----------------------------------------------------------------------
               U_S(IJK,M) = ZERO 
               V_S(IJK,M) = ZERO 
               W_S(IJK,M) = ZERO 
               IF(.NOT.(CYCLIC_AT(IJK).AND.I_OF(IJK)==IMAX2))U_S(IMJK,M)=ZERO 
               IF(.NOT.(CYCLIC_AT(IJK).AND.J_OF(IJK)==JMAX2))V_S(IJMK,M)=ZERO 
               IF(.NOT.(CYCLIC_AT(IJK).AND.K_OF(IJK)==KMAX2))W_S(IJKM,M)=ZERO 
            ENDIF 
         END DO 
      END DO 
      RETURN  
      END SUBROUTINE ZERO_NORM_VEL 
