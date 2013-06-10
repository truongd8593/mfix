!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SOURCE_ROP_s                                            C
!  Purpose: Determine source terms for solids continuity equation.     C
!                                                                      C
!  Notes: The off-diagonal coefficients are positive. The center       C
!         coefficient and the source vector are negative.              C
!         See conv_rop_s0 for additional details.                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 3 -JUL-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE SOURCE_ROP_S(A_M, B_M, M, IER) 

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE matrix 
      USE fldvar
      USE rxns
      USE run
      USE geometry
      USE indices
      USE pgcor
      USE pscor
      USE compar 
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Septadiagonal matrix A_m 
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m 
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3, 0:DIMENSION_M)
! solids phase index 
      INTEGER, INTENT(IN) :: M 
! Error index 
      INTEGER, INTENT(INOUT) :: IER 
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! DEL dot V 
      DOUBLE PRECISION :: DEL_V 
! Mass source 
      DOUBLE PRECISION :: Src 
! Indices 
      INTEGER :: I, J, K, IJK, IMJK, IJMK, IJKM 
! error message 
      CHARACTER*80     LINE 
!-----------------------------------------------
! Include statement functions
!-----------------------------------------------
      INCLUDE 'function.inc'
!-----------------------------------------------

!!$omp  parallel do private( I, J, K, IJK, IMJK, IJMK, IJKM,  DEL_V, &
!!$omp&  Src, LINE) &
!!$omp&  schedule(static)
      DO IJK = ijkstart3, ijkend3
         IF (FLUID_AT(IJK) .AND. PHASE_4_P_G(IJK)/=M .AND. &
             PHASE_4_P_S(IJK)/=M) THEN 

            I = I_OF(IJK) 
            J = J_OF(IJK) 
            K = K_OF(IJK) 
            IMJK = IM_OF(IJK) 
            IJMK = JM_OF(IJK) 
            IJKM = KM_OF(IJK) 

            DEL_V = U_S(IJK,M)*AYZ(IJK) - U_S(IMJK,M)*AYZ(IMJK) +&
                    V_S(IJK,M)*AXZ(IJK) - V_S(IJMK,M)*AXZ(IJMK) + &
                    W_S(IJK,M)*AXY(IJK) - W_S(IJKM,M)*AXY(IJKM) 

            IF (ROP_S(IJK,M) > ZERO) THEN 
               SRC = VOL(IJK)*ZMAX((-SUM_R_S(IJK,M)))/ROP_S(IJK,M) 
            ELSE 
               SRC = ZERO 
            ENDIF 

            A_M(IJK,0,M) = -(A_M(IJK,E,M)+A_M(IJK,W,M)+A_M(IJK,N,M)+&
                             A_M(IJK,S,M)+A_M(IJK,T,M)+A_M(IJK,B,M)+&
                             VOL(IJK)*ODT+ZMAX(DEL_V)+SRC) 
            B_M(IJK,M) = -(ROP_SO(IJK,M)*VOL(IJK)*ODT+&
                           ZMAX((-DEL_V))*ROP_S(IJK,M)+&
                           ZMAX(SUM_R_S(IJK,M))*VOL(IJK)) 

            IF (ABS(A_M(IJK,0,M)) < SMALL_NUMBER) THEN 
               IF (ABS(B_M(IJK,M)) < SMALL_NUMBER) THEN 
                  A_M(IJK,0,M) = -ONE            ! Equation is undefined. 
                  B_M(IJK,M) = -ROP_S(IJK,M)     ! Use existing value 
               ELSE 
!!$omp             critical
                  WRITE (LINE, '(A,I6,A,I1,A,G12.5)') 'Error: At IJK = ', IJK, &
                     ' M = ', M, ' A = 0 and b = ', B_M(IJK,M) 
                  CALL WRITE_ERROR ('SOURCE_ROP_s', LINE, 1) 
!!$omp             end critical
               ENDIF 
            ENDIF 
         ELSE
! set the value of rop_s in all wall and flow boundary cells to what is
! known for that cell
            A_M(IJK,E,M) = ZERO 
            A_M(IJK,W,M) = ZERO 
            A_M(IJK,N,M) = ZERO 
            A_M(IJK,S,M) = ZERO 
            A_M(IJK,T,M) = ZERO 
            A_M(IJK,B,M) = ZERO 
            A_M(IJK,0,M) = -ONE 
            B_M(IJK,M) = -ROP_S(IJK,M) 
         ENDIF 
      ENDDO

      RETURN  
      END SUBROUTINE SOURCE_ROP_S 

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: POINT_SOURCE_ROP_S                                      C
!  Purpose: Adds point sources to the solids continuity equation.      C
!                                                                      C
!                                                                      C
!  Author: J. Musser                                  Date: 10-JUN-13  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE POINT_SOURCE_ROP_S(B_M, M, IER) 

      use compar    
      use constant
      use geometry
      use indices
      use physprop
      use ps
      use run

! To be removed upon complete integration of point source routines.
      use bc
      use usr

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Vector b_m 
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3, 0:DIMENSION_M) 
! Solids phase index.
      INTEGER, INTENT(IN) :: M

! Error index 
      INTEGER, INTENT(INOUT) :: IER 

!----------------------------------------------- 
! Local Variables
!----------------------------------------------- 

! Indices 
      INTEGER :: IJK, I, J, K
      INTEGER :: BCV

! terms of bm expression
      DOUBLE PRECISION :: pSource

!-----------------------------------------------
! Include statement functions
!-----------------------------------------------
      INCLUDE 'function.inc'
!-----------------------------------------------
      BC_LP: do BCV = 50, DIMENSION_BC
         if(POINT_SOURCES(BCV) == 0) cycle BC_LP

         do k = BC_K_B(BCV), BC_K_T(BCV)
         do j = BC_J_S(BCV), BC_J_N(BCV)
         do i = BC_I_W(BCV), BC_I_E(BCV)

            ijk = funijk(i,j,k)
            if(fluid_at(ijk)) then

               pSource =  BC_MASSFLOW_S(BCV,M) * &
                  (VOL(IJK)/PS_VOLUME(BCV))


               if(pSource /= 0.0d0) write(*,"(3x,'Ok then... now what?')")


               B_M(IJK,M) = B_M(IJK,M) - pSource 
            endif
         enddo
         enddo
         enddo

      enddo BC_LP

      RETURN
      END SUBROUTINE POINT_SOURCE_ROP_S
