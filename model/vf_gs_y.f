!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: VF_gs_Y(F_gs, VxF_gs, IER)                             C
!  Purpose: Calculate the average drag coefficient at i, j+1/2, k and  C
!           multiply with v-momentum cell volume.                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 17-JUN-96  C
!  Reviewer:                                          Date:            C
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
      SUBROUTINE VF_GS_Y(F_GS, VXF_GS, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE geometry
      USE indices
      USE physprop
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
! 
! 
!                      Error index 
      INTEGER          IER 
! 
!                      Indices 
      INTEGER          J, IJK, IJKN 
! 
!                      Index of solids phases 
      INTEGER          M 
! 
!                      Drag array 
      DOUBLE PRECISION F_gs(DIMENSION_3, DIMENSION_M) 
! 
!                      Volume x Drag 
      DOUBLE PRECISION VxF_gs(DIMENSION_3, DIMENSION_M) 
! 
!-----------------------------------------------
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      DO M = 1, MMAX 
!


!$omp parallel do private(J,IJK,IJKN)
         DO IJK = 1, IJKMAX2 
            IF (.NOT.IP_AT_N(IJK)) THEN 
               J = J_OF(IJK) 
               IJKN = NORTH_OF(IJK) 
!
               VXF_GS(IJK,M) = AVG_Y(F_GS(IJK,M),F_GS(IJKN,M),J)*VOL_V(IJK) 
            ELSE 
               VXF_GS(IJK,M) = ZERO 
            ENDIF 
         END DO 
      END DO 
      RETURN  
      END SUBROUTINE VF_GS_Y 
