!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_trD_s(trD_s, IER)                                 C
!  Purpose: Calculate the trace of solids phase rate of strain tensor  C
!                                                                      C
!  Author: M. Syamlal                                 Date: 19-DEC-96  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: To incorporate Cartesian grid modifications                C
!  Author: Jeff Dietiker                              Date: 01-Jul-09  C
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
      SUBROUTINE CALC_TRD_S(TRD_S, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!     Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE geometry
      USE fldvar
      USE indices
      USE physprop
      USE compar
      USE sendrecv
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
      USE bc
      USE cutcell
      USE quadric
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Error index
      INTEGER          IER
!                      Indices
      INTEGER          I, J, K, IJK, IMJK, IJMK, IJKM, IM, M
!
!                      Strain rate tensor components for mth solids phase
      DOUBLE PRECISION trD_s(DIMENSION_3, DIMENSION_M)
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
      DOUBLE PRECISION :: DEL_H,Nx,Ny,Nz
      DOUBLE PRECISION :: dudx,dvdy,dwdz
      DOUBLE PRECISION :: Xw,Xe,Yn,Ys,Xc,Yc
      DOUBLE PRECISION :: Xi,Yi,Zi,Ui,Vi,Wi,Sx,Sy,Sz

      LOGICAL :: U_NODE_AT_E, U_NODE_AT_W
      LOGICAL :: V_NODE_AT_N, V_NODE_AT_S
      LOGICAL :: W_NODE_AT_T, W_NODE_AT_B
      INTEGER :: BCV
      CHARACTER(LEN=9) :: BCT
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================

!-----------------------------------------------
      INCLUDE 'function.inc'
!
!
      DO M = 1, MMAX 
!!$omp    parallel do private(ijk,i,j,k,im,imjk,ijmk,ijkm)
         DO IJK = ijkstart3, ijkend3
            IF (.NOT.WALL_AT(IJK)) THEN 
               I = I_OF(IJK) 
               J = J_OF(IJK) 
               K = K_OF(IJK) 
               IM = IM1(I) 
               IMJK = IM_OF(IJK) 
               IJMK = JM_OF(IJK) 
               IJKM = KM_OF(IJK) 

!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
               IF(.NOT.CUT_CELL_AT(IJK)) THEN

                  TRD_S(IJK,M) = (X_E(I)*U_S(IJK,M)-X_E(IM)*U_S(IMJK,M))*OX(I)*ODX&
                     (I) + (V_S(IJK,M)-V_S(IJMK,M))*ODY(J) + (W_S(IJK,M)-W_S(IJKM,&
                     M))*(OX(I)*ODZ(K)) 

               ELSE  ! CUT CELL

                  BCV = BC_ID(IJK)

                  IF(BCV > 0 ) THEN
                     BCT = BC_TYPE(BCV)
                  ELSE
                     BCT = 'NONE'
                  ENDIF

                  SELECT CASE (BCT)
                     CASE ('CG_NSW')
                        NOC_TRDS = .TRUE.
                     CASE ('CG_FSW')
                        NOC_TRDS = .FALSE.
                     CASE('CG_PSW')
                        IF(BC_HW_S(BCV,M)==UNDEFINED) THEN   ! same as NSW
                           NOC_TRDS = .TRUE.
                        ELSEIF(BC_HW_S(BCV,M)==ZERO) THEN    ! same as FSW
                           NOC_TRDS = .FALSE.
                        ELSE                                 ! partial slip
                           NOC_TRDS = .FALSE.
                        ENDIF
                     CASE ('CG_MI')
                         TRD_S(IJK,M) = ZERO
                        RETURN 
                     CASE ('CG_PO')
                         TRD_S(IJK,M) = ZERO
                        RETURN 
                  END SELECT 

                  IF(FLOW_AT(IJK)) THEN
                      TRD_S(IJK,M) = ZERO
                     RETURN 
                  ENDIF

!              du/dx

                  U_NODE_AT_E = ((.NOT.BLOCKED_U_CELL_AT(IJK)) .AND.(.NOT.WALL_U_AT(IJK)))
                  U_NODE_AT_W = ((.NOT.BLOCKED_U_CELL_AT(IMJK)).AND.(.NOT.WALL_U_AT(IMJK)))

                  IF(U_NODE_AT_E.AND.U_NODE_AT_W) THEN

                     Ui = HALF * (U_S(IJK,M) + U_S(IMJK,M))
                     Xi = HALF * (X_U(IJK) + X_U(IMJK))
                     Yi = HALF * (Y_U(IJK) + Y_U(IMJK))
                     Zi = HALF * (Z_U(IJK) + Z_U(IMJK))
                     Sx = X_U(IJK) - X_U(IMJK)
                     Sy = Y_U(IJK) - Y_U(IMJK)
                     Sz = Z_U(IJK) - Z_U(IMJK)

                     CALL GET_DEL_H(IJK,'SCALAR',Xi,Yi,Zi,Del_H,Nx,Ny,Nz)

                     IF(Sx /= ZERO) THEN
                        dudx =  (U_S(IJK,M) - U_S(IMJK,M))/Sx
                        IF(NOC_TRDS) dudx = dudx - (Ui/(Sx*DEL_H)*(Sy*Ny+Sz*Nz))  
                     ELSE
                        dudx = ZERO
                     ENDIF

                  ELSE
                     dudx = ZERO
                  ENDIF

!              dv/dy

                  V_NODE_AT_N = ((.NOT.BLOCKED_V_CELL_AT(IJK)) .AND.(.NOT.WALL_V_AT(IJK)))
                  V_NODE_AT_S = ((.NOT.BLOCKED_V_CELL_AT(IJMK)).AND.(.NOT.WALL_V_AT(IJMK)))

                  IF(V_NODE_AT_N.AND.V_NODE_AT_S) THEN

                     Vi = HALF * (V_S(IJK,M) + V_S(IJMK,M))
                     Xi = HALF * (X_V(IJK) + X_V(IJMK))
                     Yi = HALF * (Y_V(IJK) + Y_V(IJMK))
                     Zi = HALF * (Z_V(IJK) + Z_V(IJMK))
                     Sx = X_V(IJK) - X_V(IJMK)
                     Sy = Y_V(IJK) - Y_V(IJMK)
                     Sz = Z_V(IJK) - Z_V(IJMK)

                     CALL GET_DEL_H(IJK,'SCALAR',Xi,Yi,Zi,Del_H,Nx,Ny,Nz)

                     IF(Sy /= ZERO) THEN
                        dvdy =  (V_S(IJK,M) - V_S(IJMK,M))/Sy
                        IF(NOC_TRDS) dvdy = dvdy - (Vi/(Sy*DEL_H)*(Sx*Nx+Sz*Nz))  
                     ELSE
                        dvdy =  ZERO
                     ENDIF

                  ELSE
                     dvdy = ZERO
                  ENDIF

!              dw/dz

                  IF(NO_K) THEN 

                     dwdz = ZERO

                  ELSE
   
                     W_NODE_AT_T = ((.NOT.BLOCKED_W_CELL_AT(IJK)) .AND.(.NOT.WALL_W_AT(IJK)))
                     W_NODE_AT_B = ((.NOT.BLOCKED_W_CELL_AT(IJKM)).AND.(.NOT.WALL_W_AT(IJKM)))

                     IF(W_NODE_AT_T.AND.W_NODE_AT_B) THEN

                        Wi = HALF * (W_S(IJK,M) + W_S(IJKM,M))
                        Xi = HALF * (X_W(IJK) + X_W(IJKM))
                        Yi = HALF * (Y_W(IJK) + Y_W(IJKM))
                        Zi = HALF * (Z_W(IJK) + Z_W(IJKM))
                        Sx = X_W(IJK) - X_W(IJKM)   
                        Sy = Y_W(IJK) - Y_W(IJKM)
                        Sz = Z_W(IJK) - Z_W(IJKM)

                        CALL GET_DEL_H(IJK,'SCALAR',Xi,Yi,Zi,Del_H,Nx,Ny,Nz)

                        IF(Sz /= ZERO) THEN
                           dwdz =  (W_S(IJK,M) - W_S(IJKM,M))/Sz
                           IF(NOC_TRDS) dwdz = dwdz - (Wi/(Sz*DEL_H)*(Sx*Nx+Sy*Ny))  
                        ELSE
                           dwdz = ZERO
                        ENDIF

                     ELSE
                        dwdz = ZERO
                     ENDIF

                  ENDIF  ! NO_K

                  TRD_S(IJK,M) = dudx + dvdy + dwdz

               ENDIF   ! CUT CELL

! Original term:
!               TRD_S(IJK,M) = (X_E(I)*U_S(IJK,M)-X_E(IM)*U_S(IMJK,M))*OX(I)*ODX&
!                  (I) + (V_S(IJK,M)-V_S(IJMK,M))*ODY(J) + (W_S(IJK,M)-W_S(IJKM,&
!                  M))*(OX(I)*ODZ(K)) 
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            ENDIF 
         END DO 
      END DO 

      END SUBROUTINE CALC_TRD_S 

!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization
!// 350 Changed do loop limits: 1,ijkmax2-> ijkstart3, ijkend3
