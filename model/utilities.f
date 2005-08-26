      logical function isNan(x)
      double precision x
!
!                      To check whether x is NAN
      CHARACTER *80 notnumber
      isNan = .False.
      WRITE(notnumber,*) x
! To check for NaN's in x, see if x (a real number) contains a letter "n" or symbol "?"
! in which case it's a NaN (Not a Number)
!
      IF(INDEX(notnumber,'?') > 0 .OR.     &
         INDEX(notnumber,'n') > 0 .OR.     &
         INDEX(notnumber,'N') > 0 ) THEN
        isNan = .TRUE.
	RETURN
      ENDIF
      
      return
      end function isNan
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: MAX_VEL_INLET                                          C
!  Purpose: Find maximum velocity at inlets.                           C
!                                                                      C
!  Author: S. Benyahia                                Date: 26-AUG-05  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
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
      DOUBLE PRECISION FUNCTION MAX_VEL_INLET()
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
      USE param 
      USE param1 
      USE parallel 
      USE bc
      USE fldvar
      USE geometry
      USE physprop
      USE indices
      USE constant
      USE run
      USE compar 
      IMPLICIT NONE
!-----------------------------------------------
!
!  Include files defining common blocks here
!
!
!  Define local variables here
!
       INTEGER L, I, J, K, IJK, IJK2, M
!
!
!  Include files defining statement functions here
!
      INCLUDE 'ep_s1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'ep_s2.inc'
!
!
      MAX_VEL_INLET = ZERO  !initializing
!
      DO L = 1, DIMENSION_BC 
         IF (BC_DEFINED(L)) THEN 
	   IF (BC_TYPE(L) == 'MASS_INFLOW' .OR. BC_TYPE(L) == 'P_INFLOW') THEN

                  DO K = BC_K_B(L), BC_K_T(L) 
                     DO J = BC_J_S(L), BC_J_N(L) 
                        DO I = BC_I_W(L), BC_I_E(L) 
                           IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                           IJK = FUNIJK(I,J,K)

                           SELECT CASE (BC_PLANE(L))  
                           CASE ('S')  
                              IJK2 = JM_OF(IJK) 
                              IF( ABS(V_G(IJK2)) > MAX_VEL_INLET ) MAX_VEL_INLET = ABS(V_G(IJK2))
			      DO M = 1, MMAX
			        IF( ABS(V_s(IJK2, M)) > MAX_VEL_INLET ) MAX_VEL_INLET = ABS(V_s(IJK2, M))
			      ENDDO
                           CASE ('N')  
                              IF( ABS(V_G(IJK)) > MAX_VEL_INLET ) MAX_VEL_INLET = ABS(V_G(IJK))
			      DO M = 1, MMAX
			        IF( ABS(V_s(IJK, M)) > MAX_VEL_INLET ) MAX_VEL_INLET = ABS(V_s(IJK, M))
			      ENDDO  
                           CASE ('W')  
                              IJK2 = IM_OF(IJK) 
                              IF( ABS(U_G(IJK2)) > MAX_VEL_INLET ) MAX_VEL_INLET = ABS(U_G(IJK2))
			      DO M = 1, MMAX
			        IF( ABS(U_s(IJK2, M)) > MAX_VEL_INLET ) MAX_VEL_INLET = ABS(U_s(IJK2, M))
			      ENDDO
                           CASE ('E')  
                              IF( ABS(U_G(IJK)) > MAX_VEL_INLET ) MAX_VEL_INLET = ABS(U_G(IJK))
			      DO M = 1, MMAX
			        IF( ABS(U_s(IJK, M)) > MAX_VEL_INLET ) MAX_VEL_INLET = ABS(U_s(IJK, M))
			      ENDDO  
                           CASE ('B')  
                              IJK2 = KM_OF(IJK) 
                              IF( ABS(W_G(IJK2)) > MAX_VEL_INLET ) MAX_VEL_INLET = ABS(W_G(IJK2))
			      DO M = 1, MMAX
			        IF( ABS(W_s(IJK2, M)) > MAX_VEL_INLET ) MAX_VEL_INLET = ABS(W_s(IJK2, M))
			      ENDDO
                           CASE ('T')  
                              IF( ABS(W_G(IJK)) > MAX_VEL_INLET ) MAX_VEL_INLET = ABS(W_G(IJK))
			      DO M = 1, MMAX
			        IF( ABS(W_s(IJK, M)) > MAX_VEL_INLET ) MAX_VEL_INLET = ABS(W_s(IJK, M))
			      ENDDO
                           END SELECT 

			ENDDO
		     ENDDO
		  ENDDO
	     
	   ENDIF 
         ENDIF
      ENDDO 
!
      RETURN  
      END FUNCTION MAX_VEL_INLET
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CHECK_VEL_BOUND()                                      C
!  Purpose: Check velocities upper bound to be less than speed of soundC
!                                                                      C
!  Author: S. Benyahia                                Date: 25-AUG-05  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
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
      LOGICAL FUNCTION CHECK_VEL_BOUND () 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE fldvar
      USE bc
      USE geometry
      USE physprop
      USE indices
      USE run
      USE toleranc
      USE compar        
      USE mpi_utility   

      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER          M
! 
!                      Indices 
      INTEGER          IJK 
! 
!                      Approximation of speed of sound in air 
      DOUBLE PRECISION SPEED_SOUND
! 
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'ep_s2.inc'
!
!!$omp   parallel do private(IJK)
      CHECK_VEL_BOUND = .FALSE. !initialisation
      DO IJK = IJKSTART3, IJKEND3
         IF (FLUID_AT(IJK)) THEN 
            SPEED_SOUND = 200d0*DSQRT(T_G(IJK)) !this is an approximation for air
	    IF (UNITS == 'SI') SPEED_SOUND = 20d0*DSQRT(T_G(IJK)) 
!
! if no inlet velocity is specified, use the speed of sound as an upper limit
            IF(MAX_INLET_VEL == ZERO) MAX_INLET_VEL = SPEED_SOUND 
!
	    IF(ABS(U_G(IJK)) > MAX_INLET_VEL .OR. ABS(V_G(IJK)) > MAX_INLET_VEL .OR. &
	       ABS(W_G(IJK)) > MAX_INLET_VEL) CHECK_VEL_BOUND = .TRUE.
	    DO M = 1, MMAX
	      IF(ABS(U_S(IJK,M)) > MAX_INLET_VEL .OR. ABS(V_S(IJK,M)) > MAX_INLET_VEL .OR. &
	         ABS(W_S(IJK,M)) > MAX_INLET_VEL) CHECK_VEL_BOUND = .TRUE.
	    ENDDO
         ENDIF 
      END DO 
!
      RETURN  
      END FUNCTION CHECK_VEL_BOUND 

!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization
!// 350 Changed do loop limits: 1,ijkmax2-> ijkstart3, ijkend3
!// 400 Added mpi_utility module and other global reduction (sum) calls
