!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SOURCE_phi                                              C
!  Purpose: Determine source terms for phi eq. The terms               C
!     appear in the center coefficient and RHS vector. The center      C
!     coefficient and source vector are negative. The off-diagonal     C
!     coefficients are positive. S_p must be positive.                 C
!                                                                      C
!  Note: This routine is now restricted to Non-Negative scalers when   C
!     using deferred correction.                                       C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 30-APR-97  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE SOURCE_PHI(S_P, S_C, EP, PHI, M, A_M, B_M, IER) 

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE matrix 
      USE scales 
      USE physprop
      USE fldvar
      USE visc_s
      USE rxns
      USE run
      USE toleranc 
      USE geometry
      USE indices
      USE is
      USE tau_s
      USE compar  
      IMPLICIT NONE
!-----------------------------------------------
! Dummy Arguments
!-----------------------------------------------
! Source term on LHS.  Must be positive. 
      DOUBLE PRECISION, INTENT(IN) :: S_p(DIMENSION_3) 
! Source term on RHS 
      DOUBLE PRECISION, INTENT(IN) :: S_C(DIMENSION_3) 
! Phase volume fraction 
      DOUBLE PRECISION, INTENT(IN) :: EP(DIMENSION_3) 
! Phi 
      DOUBLE PRECISION, INTENT(IN) :: Phi(DIMENSION_3) 
! phase index
      INTEGER, INTENT(IN) :: M      
! Septadiagonal matrix A_m 
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M) 
! Vector b_m 
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3, 0:DIMENSION_M) 
! Error index 
      INTEGER, INTENT(INOUT) :: IER 
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices 
      INTEGER :: IJK, IMJK, IJMK, IJKM 
! error message 
      CHARACTER*80     LINE(2) 
!QX
      LOGICAL dilute_e, dilute_w, dilute_s, dilute_n, dilute_b, dilute_t
!-----------------------------------------------
! Include statement functions
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'
!-----------------------------------------------

!!$omp    parallel do private(IJK)
!Q. Xue normal case
      IF(.NOT. SOLID_RO_V) THEN
!QX

         DO IJK = ijkstart3, ijkend3 

            IF (FLUID_AT(IJK)) THEN 

! dilute flow
               IF (EP(IJK) <= DIL_EP_S) THEN 
                  A_M(IJK,E,M) = ZERO 
                  A_M(IJK,W,M) = ZERO 
                  A_M(IJK,N,M) = ZERO 
                  A_M(IJK,S,M) = ZERO 
                  A_M(IJK,T,M) = ZERO 
                  A_M(IJK,B,M) = ZERO 
                  A_M(IJK,0,M) = -ONE 
                  B_M(IJK,M) = ZERO 

! use the average boundary cell values to compute phi (sof, Aug 23 2005)
                  IMJK = IM_OF(IJK)
                  IJMK = JM_OF(IJK)
                  IJKM = KM_OF(IJK)
                  IF (EP(WEST_OF(IJK)) > DIL_EP_S .AND. &
                      .NOT.IS_AT_E(IMJK)) A_M(IJK,W,M) = ONE
                  IF (EP(EAST_OF(IJK)) > DIL_EP_S .AND. &
                      .NOT.IS_AT_E(IJK)) A_M(IJK,E,M) = ONE 
                  IF (EP(SOUTH_OF(IJK)) > DIL_EP_S .AND. &
                      .NOT.IS_AT_N(IJMK)) A_M(IJK,S,M) = ONE
                  IF (EP(NORTH_OF(IJK)) > DIL_EP_S .AND. &
                      .NOT.IS_AT_N(IJK)) A_M(IJK,N,M) = ONE
                  IF(.NOT. NO_K) THEN
                    IF (EP(BOTTOM_OF(IJK)) > DIL_EP_S .AND. &
                        .NOT.IS_AT_T(IJKM)) A_M(IJK,B,M) = ONE
                    IF (EP(TOP_OF(IJK)) > DIL_EP_S .AND. &
                        .NOT.IS_AT_T(IJK)) A_M(IJK,T,M) = ONE 
                  ENDIF 
                  
                  IF((A_M(IJK,W,M)+A_M(IJK,E,M)+&
                      A_M(IJK,S,M)+A_M(IJK,N,M)+ &
                      A_M(IJK,B,M)+A_M(IJK,T,M)) == ZERO) THEN
                     B_M(IJK,M) = -PHI(IJK)
                  ELSE
                     A_M(IJK,0,M) = -(A_M(IJK,E,M) + A_M(IJK,W,M) +&
                                      A_M(IJK,N,M) + A_M(IJK,S,M) +&
                                      A_M(IJK,T,M)+A_M(IJK,B,M))
                  ENDIF

! Normal case
               ELSE 

! Collect the terms
                  A_M(IJK,0,M) = -(A_M(IJK,E,M)+A_M(IJK,W,M)+&
                                   A_M(IJK,N,M)+A_M(IJK,S,M)+&
                                   A_M(IJK,T,M)+A_M(IJK,B,M)+S_P(IJK))

! B_m and A_m are corrected in case deferred corrections computes B_m > S_c
! see CONV_DIF_PHI_DC.
                  IF(B_M(IJK,M) < S_C(IJK) .OR. PHI(IJK) == ZERO) THEN
                     B_M(IJK,M) = -S_C(IJK)+B_M(IJK,M)
                  ELSE ! disable ELSE statememt if PHI can be negative
                     A_M(IJK,0,M) = A_M(IJK,0,M) - B_M(IJK,M)/PHI(IJK)
                     B_M(IJK,M) = -S_C(IJK)
                  ENDIF

               ENDIF   ! end if/else (ep_g(ijk)<=dil_ep_s)
            ENDIF    ! end if (fluid_at(ijk))
         ENDDO   ! end do (ijk=ijkstart3,ijkend3)

      ELSE  ! variable particle density

         DO IJK = ijkstart3, ijkend3 
!
            IF (FLUID_AT(IJK)) THEN 
!
               dilute_w = .false.
               dilute_e = .false.
               dilute_s = .false.
               dilute_n = .false.
               dilute_t = .false.
               dilute_b = .false.

               IMJK = IM_OF(IJK)
               IJMK = JM_OF(IJK)
               IJKM = KM_OF(IJK)
               IF (EP(WEST_OF(IJK)) <= DIL_EP_S .or. IS_AT_E(IMJK)) dilute_w = .true. 
               IF (EP(EAST_OF(IJK)) <= DIL_EP_S .or. IS_AT_E(IJK)) dilute_e = .true.
               IF (EP(SOUTH_OF(IJK)) <= DIL_EP_S .or. IS_AT_N(IJMK)) dilute_s = .true.
               IF (EP(NORTH_OF(IJK)) <= DIL_EP_S .or. IS_AT_N(IJK)) dilute_n = .true.
               IF(.NOT. NO_K) THEN
                  IF (EP(BOTTOM_OF(IJK)) <= DIL_EP_S .or. IS_AT_T(IJKM)) dilute_b = .true.
                  IF (EP(TOP_OF(IJK)) <= DIL_EP_S .or. IS_AT_T(IJK)) dilute_t = .true.
               ELSEIF(NO_K) then
                  dilute_b = .true.
                  dilute_t = .true.
               ENDIF

!             dilute flow
               IF(EP(IJK) <= DIL_EP_S .and. &
                    dilute_e .and. dilute_w .and. &
                    dilute_s .and. dilute_n .and. &
                    dilute_t .and. dilute_b ) then

                  A_M(IJK,E,M) = ZERO 
                  A_M(IJK,W,M) = ZERO 
                  A_M(IJK,N,M) = ZERO 
                  A_M(IJK,S,M) = ZERO 
                  A_M(IJK,T,M) = ZERO 
                  A_M(IJK,B,M) = ZERO 
                  A_M(IJK,0,M) = -ONE 
                  B_M(IJK,M) = ZERO 
!
! using the average boundary cell values to compute phi (sof, Aug 23 2005)
!
                  goto 100

                  IMJK = IM_OF(IJK)
   	       IJMK = JM_OF(IJK)
   	       IJKM = KM_OF(IJK)
                  IF (EP(WEST_OF(IJK)) > DIL_EP_S .AND. .NOT.IS_AT_E(IMJK)) A_M(IJK,W,M) = ONE
                  IF (EP(EAST_OF(IJK)) > DIL_EP_S .AND. .NOT.IS_AT_E(IJK)) A_M(IJK,E,M) = ONE 
                  IF (EP(SOUTH_OF(IJK)) > DIL_EP_S .AND. .NOT.IS_AT_N(IJMK)) A_M(IJK,S,M) = ONE
                  IF (EP(NORTH_OF(IJK)) > DIL_EP_S .AND. .NOT.IS_AT_N(IJK)) A_M(IJK,N,M) = ONE
                  IF(.NOT. NO_K) THEN
   	         IF (EP(BOTTOM_OF(IJK)) > DIL_EP_S .AND. .NOT.IS_AT_T(IJKM)) A_M(IJK,B,M) = ONE
                    IF (EP(TOP_OF(IJK)) > DIL_EP_S .AND. .NOT.IS_AT_T(IJK)) A_M(IJK,T,M) = ONE 
   	       ENDIF 
100     continue

!               
   	       IF((A_M(IJK,W,M)+A_M(IJK,E,M)+A_M(IJK,S,M)+A_M(IJK,N,M)+ &
   	           A_M(IJK,B,M)+A_M(IJK,T,M)) == ZERO) THEN
   	          B_M(IJK,M) = -PHI(IJK)              
   	       ELSE
   	         A_M(IJK,0,M) = -(A_M(IJK,E,M)+A_M(IJK,W,M)+A_M(IJK,N,M)+ &
                                     A_M(IJK,S,M)+A_M(IJK,T,M)+A_M(IJK,B,M))
   	       ENDIF
!
!             Normal case
               ELSE 
!
!               Collect the terms
                  A_M(IJK,0,M) = -(A_M(IJK,E,M)+A_M(IJK,W,M)+A_M(IJK,N,M)+A_M(IJK,&
                     S,M)+A_M(IJK,T,M)+A_M(IJK,B,M)+S_P(IJK))
! B_m and A_m are corrected in case deferred corrections computes B_m > S_c
! see CONV_DIF_PHI_DC.
                  IF(B_M(IJK,M) < S_C(IJK) .OR. PHI(IJK) == ZERO) THEN
                    B_M(IJK,M) = -S_C(IJK)+B_M(IJK,M)
                  ELSE ! disable ELSE statememt if PHI can be negative
                    A_M(IJK,0,M) = A_M(IJK,0,M) - B_M(IJK,M)/PHI(IJK)
                    B_M(IJK,M) = -S_C(IJK)
                  ENDIF
!			
                  IF (ABS(A_M(IJK,0,M)) < SMALL_NUMBER) THEN 
                     IF (ABS(B_M(IJK,M)) < SMALL_NUMBER) THEN 
                        A_M(IJK,0,M) = -ONE            ! Equation is undefined. 
                        B_M(IJK,M) = -PHI(IJK)     ! Use existing value 
                     ELSE 
!!$omp             critical
                        WRITE (LINE, '(A,I6,A,I1,A,G12.5)') 'Error: At IJK = ', IJK, &
                             ' M = ', M, ' A = 0 and b = ', B_M(IJK,M) 
                        CALL WRITE_ERROR ('SOURCE_PHI', LINE, 1) 
!!$omp             end critical
                     ENDIF
                  ENDIF
!end

   	         ENDIF 
            ENDIF 
         END DO 
      ENDIF



      RETURN  
      END SUBROUTINE SOURCE_PHI 

