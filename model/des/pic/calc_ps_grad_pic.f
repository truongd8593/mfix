!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Subroutine: CALC_PS_GRAD_PIC                                        !
!  Author: R. Garg                                                     !
!                                                                      !
!  Purpose:                                                            !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_PS_GRAD_PIC

      implicit none

      CALL CALC_PS_GRAD_PIC_GARG

      END SUBROUTINE CALC_PS_GRAD_PIC



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Subroutine: CALC_PS_GRAD_PIC                                        !
!  Author: R. Garg                                                     !
!                                                                      !
!  Purpose:                                                            !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_PS_GRAD_PIC_GARG

      USE param
      USE param1
      USE parallel
      USE physprop
      USE fldvar, only: ep_g, u_g, v_g, w_g
      USE run
      USE geometry
      USE indices
      USE bc
      USE compar
      USE sendrecv
      USE discretelement
      USE constant
      USE cutcell
      USE interpolation
      USE mfix_pic
      USE fun_avg
      USE functions

      implicit none

      ! general i, j, k indices
      INTEGER I, J, K, IJK, IPJK, IJPK, IJKP, IDIM

      integer :: I1, J1, K1


! Since EP_G is already shared across the processors, the pressure
! gradient calculation can be made a function call so that the extra
! communication of P_S can be avoided.

      DO IJK = IJKSTART3, IJKEND3
         PS_FORCE_PIC(IJK,:) = ZERO
         IF(.NOT.FLUID_AT(IJK)) CYCLE

         I = I_OF(IJK)
         J = J_OF(IJK)
         K = K_OF(IJK)

         IPJK = IP_OF(IJK)
         IJPK = JP_OF(IJK)
         IJKP = KP_OF(IJK)

! Calculate the solids pressure gradient at east face.
         IF(FLUID_AT(IPJK)) THEN
            PS_FORCE_PIC(IJK,1) = 2.0d0 *                              &
               (PIC_P_S(IPJK,1) - PIC_P_S(IJK,1)) /                    &
               (DX(I) + DX(I_OF(IPJK)))
         ELSE
            IF(PIC_P_S(IJK,1) > ZERO) THEN
               PS_FORCE_PIC(IJK,1) = 2.0d0*                            &
                  (PIC_P_S(IPJK,1) - PIC_P_S(IJK,1)) /                 &
                  (DX(I) + DX(I_OF(IPJK)))
            ELSE
               PS_FORCE_PIC(IJK,1) = ZERO
            ENDIF
         ENDIF

! Calculate the solids pressure graident at the north face.
         IF(FLUID_AT(IJPK)) THEN
            PS_FORCE_PIC(IJK,2) = 2.0d0*                               &
               (PIC_P_S(IJPK,1) - PIC_P_S(IJK,1)) /                    &
               (DY(J)+DY(J_OF(IJPK)))
         ELSE
            IF(PIC_P_S(IJK,1) > ZERO) THEN
               PS_FORCE_PIC(IJK,2) = 2.0d0*                            &
                  (PIC_P_S(IJPK,1) - PIC_P_S(IJK,1))/                  &
                  (DY(j)+Dy(j_of(ijpk)))
            ELSE
               PS_FORCE_PIC(IJK,2) = ZERO
            ENDIF
         ENDIF

! Calculate the solids pressure graident at the top face.
         IF(DO_K) THEN
            IF(FLUID_AT(IJKP)) THEN
               PS_FORCE_PIC(IJK,3) = 2.0d0*                            &
                  (PIC_P_S(IJKP,1) - PIC_P_S(IJK,1))/                  &
                  (DZ(K)+DZ(K_OF(IJKP)))
            ELSE
               IF(PIC_P_S(IJK,1).GT.ZERO) then
                  PS_FORCE_PIC(IJK,3) = 2.0d0*&
                     (PIC_P_S(IJKP,1) - PIC_P_S(IJK,1))/               &
                     (DZ(K)+DZ(K_OF(IJKP)))
               ELSE
                  PS_FORCE_PIC(IJK,3) = ZERO
               ENDIF
            ENDIF
         ENDIF
      ENDDO

! Compute the pressure gradients along the west domain bondary which
! is previously skipped.
      I1 = IMIN2
      IF(I1 == ISTART2) THEN
         DO K1 = KSTART3, KEND3
         DO J1 = JSTART3, JEND3
            IF (.NOT.IS_ON_myPE_owns(I1,J1,K1)) CYCLE
            IJK = FUNIJK(I1,J1,K1)
            IPJK = IP_OF(IJK)
            IF(PIC_P_S(IPJK,1) > ZERO) THEN
               PS_FORCE_PIC(IJK,1) = 2.0d0*                            &
                  (PIC_P_S(IPJK,1) - PIC_P_S(IJK,1)) /                 &
                  (DX(I1) + DX(I_OF(IPJK)))
            ELSE
               PS_FORCE_PIC(IJK,1) = ZERO
            ENDIF
         ENDDO
         ENDDO
      ENDIF

! Compute the pressure gradients along the south domain bondary which
! is previously skipped.
      J1 = JMIN2
      IF(J1 == JSTART2) THEN
         DO K1 = KSTART3, KEND3
         DO I1 = ISTART3, IEND3
            IF (.NOT.IS_ON_myPE_owns(I1,J1,K1)) CYCLE
            IJK = FUNIJK(I1,J1,K1)
            IJPK = JP_OF(IJK)
            IF(PIC_P_S(IJPK,1).GT.ZERO) then
               PS_FORCE_PIC(IJK,2) = 2.0d0*                            &
                  (PIC_P_S(IJPK,1) - PIC_P_S(IJK,1))/                  &
                  (DY(J) + DY(J_OF(IJPK)))
            ELSE
               PS_FORCE_PIC(IJK,2) = ZERO
            ENDIF
         ENDDO
         ENDDO
      ENDIF

! Compute the pressure gradients along the south domain bondary which
! is previously skipped.
      IF(DO_K) then
         K1 = KMIN2
         IF(K1 == KSTART2) THEN
            DO J1 = JSTART3, JEND3
            DO I1 = ISTART3, IEND3
               IF (.NOT.IS_ON_myPE_owns(I1,J1,K1)) CYCLE
               IJK = FUNIJK(I1,J1,K1)
               IJKP = KP_OF(IJK)
               IF(PIC_P_S(IJKP,1).GT.ZERO) THEN
                  PS_FORCE_PIC(IJK,3) = 2.0d0*                         &
                     (PIC_P_S(IJKP,1) - PIC_P_S(IJK,1))/               &
                     (DZ(K)+DZ(K_OF(IJKP)))
               ELSE
                  PS_FORCE_PIC(IJK, 3) = ZERO
               ENDIF
            ENDDO
            ENDDO
         ENDIF
      ENDIF

      DO IDIM = 1, merge(2,3,NO_K)
         CALL SEND_RECV(PS_FORCE_PIC(:,IDIM),1)
      ENDDO

      END SUBROUTINE CALC_PS_GRAD_PIC_GARG
