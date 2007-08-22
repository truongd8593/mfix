!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: VF_gs_X(VxF_gs, IER)                                   C
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
      SUBROUTINE VF_GS_X(VXF_GS, IER) 
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
      USE compar  
      USE drag  
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
!                      Volume x Drag 
      DOUBLE PRECISION VxF_gs(DIMENSION_3, DIMENSION_M) 
! 
!-----------------------------------------------
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'

      DO M = 1, MMAX
!$omp  parallel do private(I,IJK,IJKE)
         DO IJK = ijkstart3, ijkend3
            IF (.NOT.IP_AT_E(IJK)) THEN
               I = I_OF(IJK)
               IJKE = EAST_OF(IJK)
               ! it is not really needed to set VXF to zero, but reflects model in .pdf file.
	       VXF_GS(IJK,M) = ZERO !AVG_X(F_GS(IJK,M),F_GS(IJKE,M),I)*VOL_U(IJK)
            ELSE
               VXF_GS(IJK,M) = ZERO
            ENDIF
         END DO
      END DO

      RETURN  
      END SUBROUTINE VF_GS_X
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: VF_SS_X(VxF_ss, IER)                                   C
!  Purpose: Calculate the average Solid-Solid drag coefficient at      C
!           i+1/2, j, k and multiply with u-momentum cell volume.      C
!                                                                      C
!  Author: S. Dartevelle, LANL                        Date: 28-FEB-04  C
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
      SUBROUTINE VF_SS_X(VXF_SS, IER)
!
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
      USE compar
      USE drag
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
!                      Index of both solids phases 
      INTEGER          L, M, LM
! 
!                      Volume x Drag 
      DOUBLE PRECISION VxF_SS(DIMENSION_3, DIMENSION_LM)
! 
!-----------------------------------------------
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
!
  DO M = 1, MMAX
    DO L = 1, MMAX
      LM = FUNLM(L,M)
	  IF (L .NE. M) THEN
!$omp  parallel do private(I,IJK,IJKE)
            DO IJK = ijkstart3, ijkend3
              IF (.NOT.IP_AT_E(IJK)) THEN
                I = I_OF(IJK)
                IJKE = EAST_OF(IJK)
                VXF_SS(IJK,LM) = AVG_X(F_SS(IJK,LM),F_SS(IJKE,LM),I)*VOL_U(IJK)
              ELSE     !Impermeable wall
                VXF_SS(IJK,LM) = ZERO
              ENDIF
            END DO
      ENDIF
    END DO
  END DO
!
 RETURN
 END SUBROUTINE VF_SS_X

!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization
!// 350 Changed do loop limits: 1,ijkmax2-> ijkstart3, ijkend3
