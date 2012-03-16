! TO DO:
! 1. Check the formulation based on MCP.  
! 2. The pressure correction should be based on sum of close-packed 
!    solids?

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_VOL_FR                                             C
!  Purpose: Calculate volume fractions of phases used in pressure      C
!           corrections.                                               C
!                                                                      C
!  Notes: see mark_phase_4_cor for more details                        C
!                                                                      C      
!  Author: M. Syamlal                                 Date: 5-JUL-96   C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CALC_VOL_FR(P_STAR, RO_G, ROP_G, EP_G, ROP_S, IER) 

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE run
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
! Dummy Arguments
!-----------------------------------------------
! Solids pressure
      DOUBLE PRECISION, INTENT(IN) :: P_star(DIMENSION_3)
! Gas density
      DOUBLE PRECISION, INTENT(INOUT) :: RO_g(DIMENSION_3)
! Gas bulk density
      DOUBLE PRECISION, INTENT(INOUT) :: ROP_g(DIMENSION_3)
! Gas volume fraction
      DOUBLE PRECISION, INTENT(INOUT) :: EP_g(DIMENSION_3)
! solids bulk densities
      DOUBLE PRECISION, INTENT(INOUT) :: ROP_s(DIMENSION_3, DIMENSION_M)
! error index
      INTEGER, INTENT(INOUT) :: IER
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! volume of particle type M for GHD theory
      DOUBLE PRECISION :: VOL_M
! volume fraction of close-packed region
      DOUBLE PRECISION :: EPcp
! volume fraction of solids phase
      DOUBLE PRECISION :: EPS      
! sum of volume fractions
      DOUBLE PRECISION :: SUMVF
! Whichever phase is given by MF becomes the phase whose volume fraction
! is corrected based on all other phases present.  Generally MF gets
! defined as 0 so that the gas phase void fraction becomes corrected
! based on the summation of the volume fractions of all other phases       
      INTEGER :: MF
! Whichever phase is given by MCPl becomes the phase whose volume
! fraction is corrected based on value of maximum close packing and all
! other solids phase that can close pack.  This is basically a local
! variable for MCP but only in cells where close packed conditions 
! exist.   
      INTEGER :: MCPl
! Index of solids phase
      INTEGER :: M
! Indices
      INTEGER :: IJK      
!-----------------------------------------------
! Include statement functions
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 's_pr1.inc'
      INCLUDE 'function.inc'
      INCLUDE 's_pr2.inc'
      INCLUDE 'ep_s2.inc'
!-----------------------------------------------

! initialization
      IER = 0 

!!$omp  parallel do private(MCPl, EPCP, SUMVF, MF, M) &
!!$omp&  schedule(static)

      DO IJK = ijkstart3, ijkend3
         IF (FLUID_AT(IJK)) THEN 

! calculate the volume fraction of the solids phase based on the volume
! fractions of all other solids phases and on the value of maximum
! packing when that solids phase continuity was not solved for the
! indicated cell.
!----------------------------------------------------------------->>>
            IF (PHASE_4_P_S(IJK) /= UNDEFINED_I) THEN 
! for the current cell check if a solids phase continuity was skipped.
! this value will either be undefined (no cell was skipped in any
! solids continuity equations) or be a solids phase index (indicated
! solids continuity was skipped) for a solids phase that does close pack
! (i.e., close_packed=T) and the cell exhibits close packing conditions.
! note: this branch is never entered given the existing version of
! mark_phase_4_cor.                    
               MCPl   = PHASE_4_P_s(IJK)

! finding the value of maximum close packing based on the expression for
! plastic pressure                    
               EPCP = 1. - INV_H(P_STAR(IJK),EP_g_blend_end(ijk))
! summing the solids volume fraction of all continuum solids phases 
! that are marked as close_packed except for the solids phase which
! is also marked by phase_4_p_s (skipping its own continuity)
               SUMVF = ZERO 
               DO M = 1, MMAX 
                  IF (CLOSE_PACKED(M) .AND. M/=MCPl) SUMVF = SUMVF + EP_S(IJK,M) 
               ENDDO 
               ROP_S(IJK,MCPl) = (EPCP - SUMVF)*RO_S(MCPl) 
            ENDIF 
!-----------------------------------------------------------------<<<


! calculate the volume fraction of the 'solids' phase based on the 
! volume fractions of all other phases (including gas) if the gas 
! continuity was solved rather than that phases own continuity. 
! if the gas continuity was not solved then calculate the void 
! fraction of the gas phase based on all other solids phases.
!----------------------------------------------------------------->>>
! for the current cell check if the gas phase continuity was solved for
! the gas phase while the solids phase continuity for the indicated
! solids phase was skipped. this value will either be 0 (gas phase
! continuity was not solved) or be a solids phase index (gas continuity
! was solved while the indicated solids phase continuity was skipped)
! for a solids phase that does not close pack (i.e., close_packed=F) and
! is in greater concentration than the gas phase.
! Note: MF will always be set to 0 here given the existing version of
! mark_phase_4_cor            
            MF = PHASE_4_P_G(IJK)  

            SUMVF = ZERO 
            IF (0 /= MF) THEN 
! if gas continuity was solved rather than the solids phase, then
! include the gas phase void fraction in the summation here.
               EP_G(IJK) = ROP_G(IJK)/RO_G(IJK) 
               SUMVF = SUMVF + EP_G(IJK) 
            ENDIF 

! modified for GHD theory
            IF(TRIM(KT_TYPE) == 'GHD') THEN
              ROP_S(IJK,MMAX) = ZERO  ! mixture density
              DO M = 1, SMAX 
! volume of particle M based on fixed diamter Dp0              
                 VOL_M = PI*D_P0(M)**3/6d0 
                 IF (M /= MF) THEN
                   SUMVF = SUMVF + EP_S(IJK,M) 
                   ROP_S(IJK,MMAX) = ROP_S(IJK,MMAX) + RO_S(M)*EP_S(IJK,M) 
                 ENDIF
              ENDDO 
            ELSE
! summing the solids volume fraction of all continuum solids phases 
! except for the continuum solids phase which was marked by phase_4_p_g
! (skipping its continuity while solving gas continuity) 
              DO M = 1, MMAX 
                 IF (M /= MF) SUMVF = SUMVF + EP_S(IJK,M) 
              ENDDO 
            ENDIF   ! end if/else trim(kt_type)=='ghd'


            IF (0 == MF) THEN 
! if no gas phase continuity was solved in the current cell then correct
! the void fraction of the gas phase based on the total solids volume
! fraction of all solids phases
               EP_G(IJK) = ONE - SUMVF 

               IF (EP_G(IJK) < ZERO) THEN
!!$omp              critical
                   IER = 1 
!!$omp              end critical
               ENDIF
               ROP_G(IJK) = EP_G(IJK)*RO_G(IJK) 
            ELSE 
! else correct the volume fraction of the solids phase that was marked
               ROP_S(IJK,MF) = (ONE - SUMVF)*RO_S(MF) 
            ENDIF 
!-----------------------------------------------------------------<<<

         ENDIF    ! end if (fluid_at(ijk))
      ENDDO   ! end do (ijk=ijkstart3,ijkend3)

      CALL send_recv(EP_G, 2)
      CALL send_recv(ROP_G, 2)
      CALL send_recv(ROP_S, 2)
      
      RETURN  
      END SUBROUTINE CALC_VOL_FR 


