!TO DO:
!Check the formulation based on MCp.  The pressure correction
!should be based on all close-packed solids?
! Now PHASE_4_P_S(IJK) is undefined.
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_VOL_FR(P_star, RO_g, ROP_g, EP_g, ROP_s, IER)     C
!  Purpose: Calculate volume fractions of phases used for pressure     C
!           corrections.                                               C
!                                                                      C
!  Author: M. Syamlal                                 Date: 5-JUL-96   C
!  Reviewer:                                          Date:            C
!                                                                      C
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
      SUBROUTINE CALC_VOL_FR(P_STAR, RO_G, ROP_G, EP_G, ROP_S, IER) 
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
      USE indices
      USE physprop
      USE visc_s
      USE constant
      USE pgcor
      USE pscor
      USE compar 
      USE sendrecv 
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Solids pressure
      DOUBLE PRECISION P_star(DIMENSION_3)
!
!                      Gas density
      DOUBLE PRECISION RO_g(DIMENSION_3)
!
!                      Gas bulk density
      DOUBLE PRECISION ROP_g(DIMENSION_3)
!
!                      Gas volume fraction
      DOUBLE PRECISION EP_g(DIMENSION_3)
!
!                      solids bulk densities
      DOUBLE PRECISION ROP_s(DIMENSION_3, DIMENSION_M)
!
!                      error index
      INTEGER          IER
!
!                      volume fraction of close-packed region
      DOUBLE PRECISION EPcp
!
!                      sum of volume fractions
      DOUBLE PRECISION SUM
!
!                      Index of phase used for gas pressure correction
      INTEGER          Mf
!
!                      Indices
      INTEGER          IJK, M
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 's_pr1.inc'
      INCLUDE 'function.inc'
      INCLUDE 's_pr2.inc'
      INCLUDE 'ep_s2.inc'
!
      IER = 0 
!
!$omp  parallel do private( Mcp, EPcp, SUM, Mf, M) &
!$omp&  schedule(static)
      DO IJK = ijkstart3, ijkend3
         IF (FLUID_AT(IJK)) THEN 
!
!         bulk density of phase used for solids pr. correction
            IF (PHASE_4_P_S(IJK) /= UNDEFINED_I) THEN 
!            Mcp   = PHASE_4_P_s(IJK)
               EPCP = 1. - INV_H(P_STAR(IJK),EP_g_blend_end(ijk))
!
               SUM = ZERO 
               DO M = 1, MMAX 
                  IF (CLOSE_PACKED(M) .AND. M/=MCP) SUM = SUM + EP_S(IJK,M) 
               END DO 
               ROP_S(IJK,MCP) = (EPCP - SUM)*RO_S(MCP) 
            ENDIF 
!
!         bulk density of phase used for gas pr. correction
            MF = PHASE_4_P_G(IJK) 
!
            SUM = ZERO 
            IF (0 /= MF) THEN 
               EP_G(IJK) = ROP_G(IJK)/RO_G(IJK) 
               SUM = SUM + EP_G(IJK) 
            ENDIF 
            DO M = 1, MMAX 
               IF (M /= MF) SUM = SUM + EP_S(IJK,M) 
            END DO 
            IF (0 == MF) THEN 
               EP_G(IJK) = ONE - SUM 
!efd
               IF (EP_G(IJK) < ZERO) THEN
!$omp                  critical
                       IER = 1 
!$omp                  end critical
               ENDIF
               ROP_G(IJK) = EP_G(IJK)*RO_G(IJK) 
            ELSE 
               ROP_S(IJK,MF) = (ONE - SUM)*RO_S(MF) 
            ENDIF 
!
         ENDIF 
      END DO 

      CALL send_recv(EP_G, 2)
      CALL send_recv(ROP_G, 2)
      CALL send_recv(ROP_S, 2)
      
      RETURN  
      END SUBROUTINE CALC_VOL_FR 

!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization
!// 350 Changed do loop limits: 1,ijkmax2-> ijkstart3, ijkend3
!// 400 Added sendrecv module and send_recv calls for COMMunication
