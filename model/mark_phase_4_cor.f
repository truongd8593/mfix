!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: MARK_PHASE_4_COR(PHASE_4_P_g, PHASE_4_P_s, DO_CONT,    C
!     &                    DO_P_s, SWITCH_4_P_g, SWITCH_4_P_s, IER)    C
!  Purpose: For each cell mark the phase whose continuity will bed usedC
!           to form pressure correction equations.                     C
!                                                                      C
!  Author: M. Syamlal                                 Date: 19-JUN-96  C
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
      SUBROUTINE MARK_PHASE_4_COR(PHASE_4_P_G, PHASE_4_P_S, DO_CONT, MCP, &
         DO_P_S, SWITCH_4_P_G, SWITCH_4_P_S, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE geometry
      USE indices
      USE fldvar
      USE physprop
      USE constant
      USE compar       
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
! 
!                      Array marking the phase 
      INTEGER          PHASE_4_P_g(DIMENSION_3) 
! 
!                      Array marking the phase 
      INTEGER          PHASE_4_P_s(DIMENSION_3) 
! 
!                      whether continuity equation needs to be solved 
!                      in addition to pressure correction equation 
      LOGICAL          DO_CONT(0:DIMENSION_M) 
! 
!                      Index for close-packed solids phase 
      INTEGER          Mcp 
! 
!                      whether solids pressure correction is needed 
      LOGICAL          DO_P_s 
! 
!                      Whether different phases were used for gas 
!                      pressure correction 
      LOGICAL          SWITCH_4_P_g 
! 
!                      Whether different phases were used for solids 
!                      pressure correction 
      LOGICAL          SWITCH_4_P_s 
! 
!                      error index 
      INTEGER          IER 
! 
!                      Indices 
      INTEGER          IJK, M 
! 
!                      Index for second continuous fluid phase 
      INTEGER          Mf 
! 
!                      Count number true values 
      INTEGER          True_g, True_s 
! 
!                      check whether pressure switches are made 
      LOGICAL          SW_g(0:DIMENSION_M), SW_s(DIMENSION_M) 
! 
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'ep_s2.inc'
!
      MF = UNDEFINED_I 
      MCP = UNDEFINED_I 
      SW_G(0) = .FALSE. 
      DO_CONT(0) = .TRUE. 
      DO M = MMAX, 1, -1 
         IF (CLOSE_PACKED(M)) THEN 
            MCP = M 
         ELSE 
            MF = M 
         ENDIF 
         SW_G(M) = .FALSE. 
         SW_S(M) = .FALSE. 
         DO_CONT(M) = .TRUE. 
      END DO 

      DO IJK = ijkstart3, ijkend3
!
!
         IF (FLUID_AT(IJK)) THEN 
            IF (MF /= UNDEFINED_I) THEN 
               IF (EP_G(IJK)/RO_G(IJK) > EP_S(IJK,MF)/RO_S(MF)) THEN 
                  PHASE_4_P_G(IJK) = 0 
                  SW_G(0) = .TRUE. 
               ELSE 
                  PHASE_4_P_G(IJK) = MF 
                  SW_G(MF) = .TRUE. 
               ENDIF 
            ELSE 
               PHASE_4_P_G(IJK) = 0 
               SW_G(0) = .TRUE. 
            ENDIF 
!
!
!         Solids phase with the highest conc. that can be close-packed
!         is marked.
!          IF(EP_g(IJK) .LE. EP_star) THEN
!            PHASE_4_P_s(IJK) = Mcp
!            IF(Mcp .NE. UNDEFINED_I)SW_s(Mcp) = .TRUE.
!          ELSE
            PHASE_4_P_S(IJK) = UNDEFINED_I       !to indicate no need for pressure correction 
!          ENDIF
         ELSE 
            PHASE_4_P_G(IJK) = UNDEFINED_I       !to indicate a non-fluid cell 
            PHASE_4_P_S(IJK) = UNDEFINED_I       !to indicate a non-fluid cell 
         ENDIF 
      END DO 
      
      TRUE_G = 0 
      TRUE_S = 0 
      IF (SW_G(0)) TRUE_G = TRUE_G + 1 
      DO M = 1, MMAX 
         IF (SW_G(M)) TRUE_G = TRUE_G + 1 
         IF (SW_S(M)) TRUE_S = TRUE_S + 1 
      END DO 
      IF (TRUE_G > 1) THEN 
         SWITCH_4_P_G = .TRUE. 
      ELSE 
         SWITCH_4_P_G = .FALSE. 
      ENDIF 
!
      IF (TRUE_S > 1) THEN 
         SWITCH_4_P_S = .TRUE. 
      ELSE 
         SWITCH_4_P_S = .FALSE. 
      ENDIF 
!
      IF (MCP == UNDEFINED_I) THEN 
         DO_P_S = .FALSE. 
      ELSE 
         DO_P_S = .TRUE. 
      ENDIF 
!
!      DO_P_s = .FALSE.
!
!      IF(True_s .EQ. 0)THEN
!        DO_P_s = .FALSE.
!      ELSE
!        DO_P_s = .TRUE.
!      ENDIF
!
!  If a phase was used for pressure correction and no other phases were
!  used for pressure correction, there is no need to solve its continuity
!  equation.
!
      IF (SW_G(0) .AND. TRUE_G==1) DO_CONT(0) = .FALSE. 
!      DO 120 M = 1, MMAX
!        IF(SW_g(M) .AND. True_g .EQ. 1)DO_CONT(M) = .FALSE.
!        IF(SW_s(M) .AND. True_s .EQ. 1)DO_CONT(M) = .FALSE.
!120   CONTINUE


      RETURN  
      END SUBROUTINE MARK_PHASE_4_COR 

!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization
!// 350 Changed do loop limits: 1,ijkmax2-> ijkstart3, ijkend3
!
