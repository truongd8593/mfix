!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: VF_gs_X(F_gs, VxF_gs, IER)                             C
!  Purpose: Calculate the average drag coefficient at i+1/2, j, k and  C
!           multiply with u-momentum cell volume.                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 20-MAY-96  C
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
      SUBROUTINE VF_GS_X(F_GS, VXF_GS, IER) 
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
      USE compar        !//d
      USE sendrecv      !// 400
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
      INTEGER          I, IJK, IJKE 
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
!// 350 1229 change do loop limits: 1,ijkmax2-> ijkstart3, ijkend3    

!$omp  parallel do private(I,IJK,IJKE)
         DO IJK = ijkstart3, ijkend3
            IF (.NOT.IP_AT_E(IJK)) THEN 
               I = I_OF(IJK) 
               IJKE = EAST_OF(IJK) 
!
               VXF_GS(IJK,M) = AVG_X(F_GS(IJK,M),F_GS(IJKE,M),I)*VOL_U(IJK) 
            ELSE 
               VXF_GS(IJK,M) = ZERO 
            ENDIF 
         END DO 
      END DO 

!//? Verify if the following COMM is necessary here (as inserted for fool proof)
!// 400 Communicate VXF_GS      
!!!!      call send_recv(VXF_GS,2)
      RETURN  
      END SUBROUTINE VF_GS_X 
