!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: GAS_DRAG(A_M, B_M, VXF_GS, IER, UV, VV, WV)            C
!  Purpose: DES - Accounting for the equal and opposite drag force     C
!           on gas due to particles by introducing the drag            C      
!           as a source term                                           C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!     
      SUBROUTINE GAS_DRAG(A_M, B_M, VXF_GS, IER, UV, VV, WV)
!-----------------------------------------------
!     M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE matrix 
      USE scales 
      USE constant
      USE physprop
      USE fldvar
      USE visc_g
      USE rxns
      USE run
      USE toleranc 
      USE geometry
      USE indices
      USE is
      USE tau_g
      USE bc
      USE compar    
      USE sendrecv 
      USE discretelement
      USE drag 
      
      IMPLICIT NONE
!-----------------------------------------------
!     G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!     D u m m y   A r g u m e n t s
!-----------------------------------------------
!     
!     
!     Error index
      INTEGER          IER
!     
!     Indices
      INTEGER          IJK, NC
!     
!     Phase index
      INTEGER          M, UV, VV, WV
!     
      DOUBLE PRECISION VCOR
      DOUBLE PRECISION A_M(DIMENSION_3, -3:3, 0:DIMENSION_M)
      DOUBLE PRECISION B_M(DIMENSION_3, 0:DIMENSION_M)
      DOUBLE PRECISION VXF_GS(DIMENSION_3, DIMENSION_M)  
!     
!-----------------------------------------------
      INCLUDE 'function.inc'

      VCOR = 1
      
      IF(UV.EQ.1) THEN
         DO M = 1, MMAX
            DO IJK = 1, IJKMAX2
                 VCOR = VOL(IJK) !*DZ(K_OF(IJK))
                 A_M(IJK,0,0) = A_M(IJK,0,0) - (F_GS(IJK,M)*VCOR)
                 B_M(IJK,0) = B_M(IJK,0) - (F_GS(IJK,M)*U_S(IJK,M)*VCOR)
            END DO
         END DO


      ELSE IF(VV.EQ.1) THEN
         DO M = 1, MMAX
            DO IJK = 1, IJKMAX2
                 VCOR = VOL(IJK) !*DZ(K_OF(IJK))
                 A_M(IJK,0,0) = A_M(IJK,0,0) - (F_GS(IJK,M)*VCOR)
                 B_M(IJK,0) = B_M(IJK,0) - (F_GS(IJK,M)*V_S(IJK,M)*VCOR)
            END DO
         END DO


      ELSE IF(WV.EQ.1) THEN
         DO M = 1, MMAX
            DO IJK = 1, IJKMAX2
                 VCOR = VOL(IJK) !*DZ(K_OF(IJK))
                 A_M(IJK,0,0) = A_M(IJK,0,0) - (F_GS(IJK,M)*VCOR)
                 B_M(IJK,0) = B_M(IJK,0)  - (F_GS(IJK,M)*W_S(IJK,M)*VCOR)
            END DO
         END DO

      END IF

      RETURN
      END SUBROUTINE GAS_DRAG
