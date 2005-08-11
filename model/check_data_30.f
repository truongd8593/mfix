!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CHECK_DATA_30                                          C
!  Purpose: Check whether the sum of reaction rates is zero and the sumC
!           of mass fractions is 1.0                                   C
!           and EP_g >= EP_star. Set miscellaneous constants           C
!                                                                      C
!  Author: M. Syamlal                                 Date: 27-OCT-92  C
!  Reviewer: W. Rogers                                Date: 11-DEC-92  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables: ABORT, SUM                                         C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CHECK_DATA_30 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE toleranc 
      USE fldvar
      USE rxns
      USE visc_s
      USE visc_g
      USE geometry
      USE run
      USE constant
      USE physprop
      USE indices
      USE funits 
      USE compar 
      USE mpi_utility  
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!
!                      Indices
      INTEGER          I, J, K, IJK
!
!                      Solids phase
      INTEGER          M
!
!                      Species index
      INTEGER          N
!
!                      Logical variable to set, if there is an error
      LOGICAL          ABORT
!
!                      Logical variable to set, if there is a message
      LOGICAL          MESSAGE
!
!                      sum of mass fractions or reaction rates
      DOUBLE PRECISION SUM
!
!                      Logical variable to set, if there is a message
      LOGICAL          MESSAGE_X_g
!
!                      minimum sum of gas species mass fractions
      DOUBLE PRECISION SUM_MIN_g
!
!                      maximum sum of gas species mass fractions
      DOUBLE PRECISION SUM_MAX_g
!
!                      Location of minimum
      INTEGER          I_MIN_g, J_MIN_g, K_MIN_g
!
!                      Location of maximum
      INTEGER          I_MAX_g, J_MAX_g, K_MAX_g
!
!                      Distribution of the sum of gas species mass fractions
      INTEGER          COUNT_g(9)
!
!                      Total count
      INTEGER          SUM_COUNT
!
!                      Fractional distribution of the sum of species mass fr.
      DOUBLE PRECISION FR_COUNT(9)
!
!                      Logical variable to set, if there is a message
      LOGICAL          MESSAGE_X_s (DIMENSION_M)
!
!                      minimum sum of solids species mass fractions
      DOUBLE PRECISION SUM_MIN_s (DIMENSION_M)
!
!                      maximum sum of gas species mass fractions
      DOUBLE PRECISION SUM_MAX_s (DIMENSION_M)
!
!                      Location of minimum
      INTEGER          I_MIN_s(DIMENSION_M), J_MIN_s(DIMENSION_M),&
                       K_MIN_s(DIMENSION_M) 
!
!                      Location of maximum
      INTEGER          I_MAX_s(DIMENSION_M), J_MAX_s(DIMENSION_M),&
                       K_MAX_s(DIMENSION_M)
!
!                      Distribution of the sum of gas species mass fractions
      INTEGER          COUNT_s (DIMENSION_M, 9)
!
!                      Do-loop counter
      INTEGER          L, LM
!
!                      There is a discrepancy in rxn sums
      LOGICAL          MESSAGE_rxnsum
!
!                      maximum discrepancy in rxn sums
      DOUBLE PRECISION RXNSUM_MAX
!
!                      Location of maximum discrepancy in rxn sums
      INTEGER          I_RXNSUM_MAX, J_RXNSUM_MAX, K_RXNSUM_MAX
!
!                      Distribution of the sum of gas species mass fractions
      INTEGER          COUNT_RXNSUM0, COUNT_RXNSUM1
!
!                      There is a discrepancy in interphase mass transfer
      LOGICAL          MESSAGE_masstr(0:DIMENSION_M)
!
!                      maximum discrepancy in interphase mass transfer
      DOUBLE PRECISION masstr_MAX(0:DIMENSION_M)
!
!                      Location of maximum discrepancy in interphase mass transfer
      INTEGER          I_masstr_MAX(0:DIMENSION_M), J_masstr_MAX(0:DIMENSION_M),&
                       K_masstr_MAX(0:DIMENSION_M)
!
!                      number of cells with minor and major discrepancy
      INTEGER          COUNT_masstr0(0:DIMENSION_M), COUNT_masstr1(0:DIMENSION_M)
!
!                      errror  flag
      INTEGER          IER
!
!-----------------------------------------------
      INCLUDE 'function.inc'

! For DM parallel runs we redo these checks again from here so that all processors can write
! log files.  There is a goto statement at the end to start from statement no. 1.
1     MESSAGE_rxnsum = .false.
      RXNSUM_MAX     = ZERO
      COUNT_RXNSUM0  = 0
      COUNT_RXNSUM1  = 0

      DO L = 0, MMAX
        MESSAGE_masstr(L) = .false.
        masstr_MAX(L)     = ZERO
        COUNT_masstr0(L)  = 0
        COUNT_masstr1(L)  = 0
      END DO 
!

      MESSAGE_X_G = .FALSE. 
      I_MIN_G = 0 
      J_MIN_G = 0 
      K_MIN_G = 0 
      I_MAX_G = 0 
      J_MAX_G = 0 
      K_MAX_G = 0 
      SUM_MIN_G = ONE 
      SUM_MAX_G = ONE 
      COUNT_G = 0 
      L = 10 
      M = 1 
      IF (MMAX > 0) THEN 
         MESSAGE_X_S(:MMAX) = .FALSE. 
         I_MIN_S(:MMAX) = 0 
         J_MIN_S(:MMAX) = 0 
         K_MIN_S(:MMAX) = 0 
         I_MAX_S(:MMAX) = 0 
         J_MAX_S(:MMAX) = 0 
         K_MAX_S(:MMAX) = 0 
         SUM_MIN_S(:MMAX) = ONE 
         SUM_MAX_S(:MMAX) = ONE 
         COUNT_S(:MMAX,:9) = 0 
         L = 10 
         M = MMAX + 1 
      ENDIF 
      CALL START_LOG 
      ABORT = .FALSE. 
      MESSAGE = .FALSE. 

      DO K = KSTART2, KEND2 
         DO J = JSTART2, JEND2 
            DO I = ISTART2, IEND2 
               IJK = FUNIJK(I,J,K) 
               IF (.NOT.WALL_AT(IJK)) THEN 
	       
	         IF(FLOW_AT(IJK)) THEN
		 !  The diffusivities must be zero in inflow and outflow
		 !  cells
                   IF(MU_gt(IJK) /= ZERO) THEN
                     IF (.NOT.MESSAGE) THEN 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1000) TIME 
                        MESSAGE = .TRUE. 
                     ENDIF 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1140) I, J, K, MU_gt(IJK), 'MU_gt' 
                     ABORT = .TRUE. 
		   ENDIF
		   
                   IF(LAMBDA_gt(IJK) /= ZERO) THEN
                     IF (.NOT.MESSAGE) THEN 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1000) TIME 
                        MESSAGE = .TRUE. 
                     ENDIF 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1140) I, J, K, LAMBDA_gt(IJK), 'LAMBDA_gt' 
                     ABORT = .TRUE. 
		   ENDIF
		   
                   IF(K_g(IJK) /= ZERO) THEN
                     IF (.NOT.MESSAGE) THEN 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1000) TIME 
                        MESSAGE = .TRUE. 
                     ENDIF 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1140) I, J, K, K_g(IJK), 'K_g' 
                     ABORT = .TRUE. 
		   ENDIF
		   
		   DO N = 1, NMAX(0)
                     IF( DIF_g(IJK, N) /= ZERO) THEN
                       IF (.NOT.MESSAGE) THEN 
                         IF(DMP_LOG)WRITE (UNIT_LOG, 1000) TIME 
                         MESSAGE = .TRUE. 
                       ENDIF 
                       IF(DMP_LOG)WRITE (UNIT_LOG, 1140) I, J, K, DIF_g(IJK, N), 'DIF_g' 
                       ABORT = .TRUE. 
		     ENDIF
		   ENDDO
		   
		   DO M = 1, MMAX
                     IF(MU_s(IJK, M) /= ZERO) THEN
                       IF (.NOT.MESSAGE) THEN 
                          IF(DMP_LOG)WRITE (UNIT_LOG, 1000) TIME 
                          MESSAGE = .TRUE. 
                       ENDIF 
                       IF(DMP_LOG)WRITE (UNIT_LOG, 1140) I, J, K, MU_s(IJK, M), 'MU_s' 
                       ABORT = .TRUE. 
    		     ENDIF
		   
                     IF(LAMBDA_s(IJK, M) /= ZERO) THEN
                       IF (.NOT.MESSAGE) THEN 
                         IF(DMP_LOG)WRITE (UNIT_LOG, 1000) TIME 
                         MESSAGE = .TRUE. 
                       ENDIF 
                       IF(DMP_LOG)WRITE (UNIT_LOG, 1140) I, J, K, LAMBDA_s(IJK, M), 'LAMBDA_s' 
                       ABORT = .TRUE. 
	  	     ENDIF
		   
                     IF(K_s(IJK, M) /= ZERO) THEN
                       IF (.NOT.MESSAGE) THEN 
                         IF(DMP_LOG)WRITE (UNIT_LOG, 1000) TIME 
                         MESSAGE = .TRUE. 
                       ENDIF 
                       IF(DMP_LOG)WRITE (UNIT_LOG, 1140) I, J, K, K_s(IJK, M), 'K_s' 
                       ABORT = .TRUE. 
		     ENDIF
		   
		     DO N = 1, NMAX(M)
                       IF( DIF_s(IJK, M, N) /= ZERO) THEN
                         IF (.NOT.MESSAGE) THEN 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1000) TIME 
                           MESSAGE = .TRUE. 
                         ENDIF 
                         IF(DMP_LOG)WRITE (UNIT_LOG, 1140) I, J, K, DIF_s(IJK, M, N), 'DIF_s' 
                         ABORT = .TRUE. 
		       ENDIF
		     ENDDO
		   	   
		   ENDDO
		   
	   
		 ENDIF
!
!   Gas viscosity, conductivity and specific heat must be positive
!
                  IF (MU_G(IJK) < ZERO) THEN 
                     IF (.NOT.MESSAGE) THEN 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1000) TIME 
                        MESSAGE = .TRUE. 
                     ENDIF 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1110) I, J, K, MU_G(IJK) 
                     ABORT = .TRUE. 
                  ENDIF 
!
                  IF (MW_MIX_G(IJK) <= ZERO) THEN 
                     IF (.NOT.MESSAGE) THEN 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1000) TIME 
                        MESSAGE = .TRUE. 
                     ENDIF 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1111) I, J, K, MW_MIX_G(IJK) 
                     ABORT = .TRUE. 
                  ENDIF 
!
                  IF (ENERGY_EQ) THEN 
                     IF (K_G(IJK) < ZERO) THEN 
                        IF (.NOT.MESSAGE) THEN 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1000) TIME 
                           MESSAGE = .TRUE. 
                        ENDIF 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1120) I, J, K, K_G(IJK) 
                        ABORT = .TRUE. 
                     ENDIF 
!
                     IF (C_PG(IJK) <= ZERO) THEN 
                        IF (.NOT.MESSAGE) THEN 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1000) TIME 
                           MESSAGE = .TRUE. 
                        ENDIF 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1130) I, J, K, C_PG(IJK) 
                        ABORT = .TRUE. 
                     ENDIF 
                  ENDIF 
		   
		  DO N = 1, NMAX(0)
                    IF( DIF_g(IJK, N) < ZERO) THEN
                      IF (.NOT.MESSAGE) THEN 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1000) TIME 
                        MESSAGE = .TRUE. 
                      ENDIF 
                      IF(DMP_LOG)WRITE (UNIT_LOG, 1131) I, J, K, DIF_g(IJK, N) 
                      ABORT = .TRUE. 
		    ENDIF
		  ENDDO
		   
!
!  Sum of over all reaction rates for phases should be zero
!
                  SUM = SUM_R_G(IJK) 
                  IF (MMAX > 0) THEN 
                     DO M = 1, MMAX 
                        SUM = SUM + SUM_R_S(IJK,M) 
                     END DO 
                  ENDIF 
                  IF (ABS(SUM) > SMALL_NUMBER) THEN
                     MESSAGE_rxnsum = .true.
		     if(abs(sum) > abs(RXNSUM_MAX))then
		       RXNSUM_MAX   = sum
		       I_RXNSUM_MAX = i
		       J_RXNSUM_MAX = j
		       K_RXNSUM_MAX = k
		     endif
                     IF (.NOT.MESSAGE) THEN 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1000) TIME 
                        MESSAGE = .TRUE. 
                     ENDIF 
!                     IF(DMP_LOG)WRITE (UNIT_LOG, 1100) I, J, K, SUM 
                     IF (ABS(SUM) > TOL_COM) then
		       COUNT_RXNSUM1 = COUNT_RXNSUM1 + 1
		       ABORT = .TRUE.
		     else
		       COUNT_RXNSUM0 = COUNT_RXNSUM0 + 1
		     endif
                  ENDIF 
!
!  The net production of each phase should match the total mass transferred
!  from the other phases
!
                  DO L = 0, MMAX 
                     IF (L == 0) THEN 
                        SUM = SUM_R_G(IJK) 
                     ELSE 
                        SUM = SUM_R_S(IJK,L) 
                     ENDIF 
                     DO M = 0, MMAX 
                        IF (M > L) THEN 
                           LM = L + 1 + (M - 1)*M/2 
                           SUM = SUM - R_PHASE(IJK,LM) 
                        ELSE IF (L > M) THEN 
                           LM = M + 1 + (L - 1)*L/2 
                           SUM = SUM + R_PHASE(IJK,LM) 
                        ENDIF 
                     END DO 
                     IF (ABS(SUM) > SMALL_NUMBER) THEN
		        MESSAGE_masstr(L) = .true. 
		        if(abs(sum) > abs(masstr_MAX(L)))then
		          masstr_MAX(L)   = sum
		          I_masstr_MAX(L) = i
		          J_masstr_MAX(L) = j
		          K_masstr_MAX(L) = k
		        endif
                        IF (.NOT.MESSAGE) THEN 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1000) TIME 
                           MESSAGE = .TRUE. 
                        ENDIF 
!                        IF(DMP_LOG)WRITE (UNIT_LOG, 1101) I, J, K, L, SUM 
                        IF (ABS(SUM) > TOL_COM) then
		          COUNT_masstr1(L) = COUNT_masstr1(L) + 1
			  ABORT = .TRUE. 
			ELSE
		          COUNT_masstr0(L) = COUNT_masstr0(L) + 1
			ENDIF
                     ENDIF 
                  END DO 

!
!  R_gp, RoX_gc, R_sp, and RoX_sc must be non-negative
!
                  IF(SPECIES_EQ(0))THEN
                    DO N = 1, NMAX(0) 
                       IF (R_GP(IJK,N) < ZERO) THEN 
                         IF (.NOT.MESSAGE) THEN 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1000) TIME 
                           MESSAGE = .TRUE. 
                         ENDIF 
                         IF(DMP_LOG)WRITE (UNIT_LOG, 1103) I, J, K, R_GP(IJK,N), N 
                         ABORT = .TRUE. 
                       ENDIF 
                       IF (ROX_GC(IJK,N) < ZERO) THEN 
                         IF (.NOT.MESSAGE) THEN 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1000) TIME 
                           MESSAGE = .TRUE. 
                         ENDIF 
                         IF(DMP_LOG)WRITE (UNIT_LOG, 1104) I, J, K, ROX_GC(IJK,N), N 
                         ABORT = .TRUE. 
                       ENDIF 
                    END DO 
                  ENDIF

                  DO M = 1, MMAX 
                    IF(SPECIES_EQ(M))THEN
                      DO N = 1, NMAX(M) 
                        IF (R_SP(IJK,M,N) < ZERO) THEN 
                           IF (.NOT.MESSAGE) THEN 
                              IF(DMP_LOG)WRITE (UNIT_LOG, 1000) TIME 
                              MESSAGE = .TRUE. 
                           ENDIF 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1105) I, J, K, M, R_SP(IJK,M,N), N 
                           ABORT = .TRUE. 
                        ENDIF 
                        IF (ROX_SC(IJK,M,N) < ZERO) THEN 
                           IF (.NOT.MESSAGE) THEN 
                              IF(DMP_LOG)WRITE (UNIT_LOG, 1000) TIME 
                              MESSAGE = .TRUE. 
                           ENDIF 
                           IF(DMP_LOG)WRITE(UNIT_LOG,1106)I,J,K,M,ROX_SC(IJK,M,N),N 
                           ABORT = .TRUE. 
                        ENDIF 
                      END DO 
                    ENDIF
                  END DO 

!
!  Sum of gas mass fractions should be one
!
                  IF (SPECIES_EQ(0)) THEN 
                     SUM = ZERO 
                     N = 1 
                     IF (NMAX(0) > 0) THEN 
                        DO N = 1, NMAX(0) 
                           SUM = SUM + X_G(IJK,N) 
                        END DO 
                        N = NMAX(0) + 1 
                     ENDIF 
                     IF (ABS(ONE - SUM) > TOL_COM) THEN 
                        IF (.NOT.MESSAGE) THEN 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1000) TIME 
                           MESSAGE = .TRUE. 
                        ENDIF 
                        MESSAGE_X_G = .TRUE. 
                     ENDIF 
                     IF (SUM < SUM_MIN_G) THEN 
                        SUM_MIN_G = SUM 
                        I_MIN_G = I 
                        J_MIN_G = J 
                        K_MIN_G = K 
                     ELSE IF (SUM > SUM_MAX_G) THEN 
                        SUM_MAX_G = SUM 
                        I_MAX_G = I 
                        J_MAX_G = J 
                        K_MAX_G = K 
                     ENDIF 
                     IF (SUM < 0.9) THEN         ! < 0.9 
                        COUNT_G(1) = COUNT_G(1) + 1 
                     ELSE IF (SUM < 0.99) THEN   ! 0.9    - 0.99 
                        COUNT_G(2) = COUNT_G(2) + 1 
                     ELSE IF (SUM < 0.999) THEN  ! 0.99   - 0.999 
                        COUNT_G(3) = COUNT_G(3) + 1 
                     ELSE IF (SUM < 0.9999) THEN ! 0.999  - 0.9999 
                        COUNT_G(4) = COUNT_G(4) + 1 
                     ELSE IF (SUM < 1.0001) THEN ! 0.9999 - 1.0001 
                        COUNT_G(5) = COUNT_G(5) + 1 
                     ELSE IF (SUM < 1.001) THEN  ! 1.0001 - 1.001 
                        COUNT_G(6) = COUNT_G(6) + 1 
                     ELSE IF (SUM < 1.01) THEN   ! 1.001  - 1.01 
                        COUNT_G(7) = COUNT_G(7) + 1 
                     ELSE IF (SUM < 1.1) THEN    ! 1.01   - 1.1 
                        COUNT_G(8) = COUNT_G(8) + 1 
                     ELSE                        ! > 1.1 
                        COUNT_G(9) = COUNT_G(9) + 1 
                     ENDIF 
                  ENDIF 
!
!
                  DO M = 1, MMAX 
!
!  Sum of solids mass fractions should be one
!
                     IF (SPECIES_EQ(M)) THEN 
                        SUM = ZERO 
                        N = 1 
                        IF (NMAX(M) > 0) THEN 
                           DO N = 1, NMAX(M) 
                              SUM = SUM + X_S(IJK,M,N) 
                           END DO 
                           N = NMAX(M) + 1 
                        ENDIF 
                        IF (ROP_S(IJK,M) /= ZERO) THEN 
                           IF (ABS(ONE - SUM) > TOL_COM) THEN 
                              IF (.NOT.MESSAGE) THEN 
                                 IF(DMP_LOG)WRITE (UNIT_LOG, 1000) TIME 
                                 MESSAGE = .TRUE. 
                              ENDIF 
                              MESSAGE_X_S(M) = .TRUE. 
                           ENDIF 
                           IF (SUM < SUM_MIN_S(M)) THEN 
                              SUM_MIN_S(M) = SUM 
                              I_MIN_S(M) = I 
                              J_MIN_S(M) = J 
                              K_MIN_S(M) = K 
                           ELSE IF (SUM > SUM_MAX_S(M)) THEN 
                              SUM_MAX_S(M) = SUM 
                              I_MAX_S(M) = I 
                              J_MAX_S(M) = J 
                              K_MAX_S(M) = K 
                           ENDIF 
                           IF (SUM < 0.9) THEN   ! < 0.9 
                              COUNT_S(M,1) = COUNT_S(M,1) + 1 
!                                                ! 0.9    - 0.99
                           ELSE IF (SUM < 0.99) THEN 
                              COUNT_S(M,2) = COUNT_S(M,2) + 1 
!                                                ! 0.99   - 0.999
                           ELSE IF (SUM < 0.999) THEN 
                              COUNT_S(M,3) = COUNT_S(M,3) + 1 
!                                                ! 0.999  - 0.9999
                           ELSE IF (SUM < 0.9999) THEN 
                              COUNT_S(M,4) = COUNT_S(M,4) + 1 
!                                                ! 0.9999 - 1.0001
                           ELSE IF (SUM < 1.0001) THEN 
                              COUNT_S(M,5) = COUNT_S(M,5) + 1 
!                                                ! 1.0001 - 1.001
                           ELSE IF (SUM < 1.001) THEN 
                              COUNT_S(M,6) = COUNT_S(M,6) + 1 
!                                                ! 1.001  - 1.01
                           ELSE IF (SUM < 1.01) THEN 
                              COUNT_S(M,7) = COUNT_S(M,7) + 1 
!                                                ! 1.01   - 1.1
                           ELSE IF (SUM < 1.1) THEN 
                              COUNT_S(M,8) = COUNT_S(M,8) + 1 
                           ELSE                  ! > 1.1 
                              COUNT_S(M,9) = COUNT_S(M,9) + 1 
                           ENDIF 
                        ENDIF 
                     ENDIF 
!
!   Solids viscosity, conductivity and specific heat must be positive
!
                     IF (MU_S(IJK,M) < ZERO) THEN 
                        IF (.NOT.MESSAGE) THEN 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1000) TIME 
                           MESSAGE = .TRUE. 
                        ENDIF 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1310) I, J, K, M, MU_S(IJK,M) 
                        ABORT = .TRUE. 
                     ENDIF 
!
                     IF (ENERGY_EQ) THEN 
                        IF (K_S(IJK,M) < ZERO) THEN 
                           IF (.NOT.MESSAGE) THEN 
                              IF(DMP_LOG)WRITE (UNIT_LOG, 1000) TIME 
                              MESSAGE = .TRUE. 
                           ENDIF 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1320) I, J, K, M, K_S(IJK,M) 
                           ABORT = .TRUE. 
                        ENDIF 
!
                        IF (C_PS(IJK,M) <= ZERO) THEN 
                           IF (.NOT.MESSAGE) THEN 
                              IF(DMP_LOG)WRITE (UNIT_LOG, 1000) TIME 
                              MESSAGE = .TRUE. 
                           ENDIF 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1330) I, J, K, M, C_PS(IJK,M) 
                           ABORT = .TRUE. 
                        ENDIF 
!
                        IF (T_S(IJK,M)<=TMIN .OR. T_S(IJK,M)>=TMAX) THEN 
                           IF (.NOT.MESSAGE) THEN 
                              IF(DMP_LOG)WRITE (UNIT_LOG, 1000) TIME 
                              MESSAGE = .TRUE. 
                           ENDIF 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1410) I, J, K, M, T_S(IJK,M) 
                           ABORT = .TRUE. 
                        ENDIF 
                     ENDIF 
		     
		     DO N = 1, NMAX(M)
                       IF( DIF_s(IJK, M, N) < ZERO) THEN
                         IF (.NOT.MESSAGE) THEN 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1000) TIME 
                           MESSAGE = .TRUE. 
                         ENDIF 
                         IF(DMP_LOG)WRITE (UNIT_LOG, 1331) I, J, K, M, DIF_s(IJK, M, N) 
                         ABORT = .TRUE. 
		       ENDIF
		     ENDDO
                  END DO 
                  IF (ENERGY_EQ) THEN 
                     IF (T_G(IJK)<=TMIN .OR. T_G(IJK)>=TMAX) THEN 
                        IF (.NOT.MESSAGE) THEN 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1000) TIME 
                           MESSAGE = .TRUE. 
                        ENDIF 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1400) I, J, K, T_G(IJK) 
                        ABORT = .TRUE. 
                     ENDIF 
                  ENDIF 
!
               ENDIF 
            END DO 
         END DO 
      END DO
       
      call global_all_or(MESSAGE_rxnsum)
      IF (MESSAGE_rxnsum) THEN 
	 call global_all_sum(COUNT_RXNSUM0)
	 call global_all_sum(COUNT_RXNSUM1)
         IF(DMP_LOG)WRITE (UNIT_LOG, 1415) COUNT_RXNSUM0, COUNT_RXNSUM1, RXNSUM_MAX, &
	                        I_RXNSUM_MAX, J_RXNSUM_MAX, K_RXNSUM_MAX 
      ENDIF 

      DO L = 0, MMAX
        call global_all_or(MESSAGE_masstr(L))
        IF (MESSAGE_masstr(L)) THEN 
	   call global_all_sum(COUNT_masstr0(L))
	   call global_all_sum(COUNT_masstr1(L))
           IF(DMP_LOG)WRITE (UNIT_LOG, 1420) L, COUNT_masstr0(L), COUNT_masstr1(L),&
	         masstr_MAX(L), &
	         I_masstr_MAX(L), J_masstr_MAX(L), K_masstr_MAX(L) 
        ENDIF 
      ENDDO
      
      call global_all_or(MESSAGE_X_G)
      IF (MESSAGE_X_G) THEN 
	 call global_all_sum(COUNT_G)
         SUM_COUNT = 0 
         DO L = 1, 9 
            SUM_COUNT = SUM_COUNT + COUNT_G(L) 
         END DO 
         IF(DMP_LOG)WRITE (UNIT_LOG, 1430) 
         FR_COUNT = DBLE(COUNT_G)/DBLE(SUM_COUNT) 
         L = 10 
         IF(DMP_LOG)WRITE (UNIT_LOG, 1432) (COUNT_G(L),FR_COUNT(L),L=1,9) 
         IF(DMP_LOG)WRITE (UNIT_LOG, 1434) SUM_MIN_G, I_MIN_G, J_MIN_G, K_MIN_G 
         IF(DMP_LOG)WRITE (UNIT_LOG, 1436) SUM_MAX_G, I_MAX_G, J_MAX_G, K_MAX_G 
      ENDIF 
!
      DO M = 1, MMAX 
         call global_all_or(MESSAGE_X_S(M))
         IF (MESSAGE_X_S(M)) THEN 
            call global_all_sum(COUNT_S(M,:))
            SUM_COUNT = 0 
            DO L = 1, 9 
               SUM_COUNT = SUM_COUNT + COUNT_S(M,L) 
            END DO 
            L = 10 
            IF(DMP_LOG)WRITE (UNIT_LOG, 1440) M 
            FR_COUNT = DBLE(COUNT_S(M,:))/DBLE(SUM_COUNT) 
            L = 10 
            IF(DMP_LOG)WRITE (UNIT_LOG, 1442) (COUNT_S(M,L),FR_COUNT(L),L=1,9) 
            IF(DMP_LOG)WRITE (UNIT_LOG, 1444) SUM_MIN_S(M), I_MIN_S(M), J_MIN_S(M), &
               K_MIN_S(M) 
            IF(DMP_LOG)WRITE (UNIT_LOG, 1446) SUM_MAX_S(M), I_MAX_S(M), J_MAX_S(M), &
               K_MAX_S(M) 
         ENDIF 
      END DO 
      IF (MESSAGE .AND. DMP_LOG)WRITE (UNIT_LOG, 1500) 
!
      CALL END_LOG 
      IF (ABORT) then
        if(.not.dmp_log)then
	  !enable dmp_log, open logfile, and redo the check (goto 1), so that this PE can
	  !  write the error message before aborting
          call open_pe_log(ier)
	  dmp_log = .true.
	  goto 1 

	else
          CALL MFIX_EXIT(myPE)
	endif
      endif
!
      RETURN  
 1000 FORMAT(/1X,70('*')//' From: CHECK_DATA_30',5X,'Time = ',G12.5,/&
         ' Message: One or more of the following errors detected:',/&
         '   1. Discrepancies in the reaction rates.',/&
         '   2. Viscosity, MW, conductivity, or specific heat < zero.',/&
         '   3. The sum of mass fractions is not equal to one.',/&
         '   4. Temperatures at the upper or lower bound.',/&
         '   5. The rate of production of phases (SUM_R_g or SUM_R_s)',/&
         '      and the interphase mass transfer rates (R_Phase) are',/&
         '      inconsistent (in subroutine RRATES).',/4X,'I',T14,'J',T24,'K',&
         T34,'M',T45,'Value') 
 1100 FORMAT(1X,I4,T11,I4,T21,I4,T41,G12.5,'  Sum of rates .NE. 0') 
 1101 FORMAT(1X,I4,T11,I4,T21,I4,T31,I4,T41,G12.5,'  See message 5') 
 1103 FORMAT(1X,I4,T11,I4,T21,I4,T41,G12.5,'  R_gp < 0 for N=',I2) 
 1104 FORMAT(1X,I4,T11,I4,T21,I4,T41,G12.5,'  RoX_gc<0 for N=',I2) 
 1105 FORMAT(1X,I4,T11,I4,T21,I4,T31,I4,T41,G12.5,'  R_sp < 0 for N=',I2) 
 1106 FORMAT(1X,I4,T11,I4,T21,I4,T31,I4,T41,G12.5,'  RoX_sc<0 for N=',I2) 
 1110 FORMAT(1X,I4,T11,I4,T21,I4,T41,G12.5,'  MU_g .LT. 0') 
 1111 FORMAT(1X,I4,T11,I4,T21,I4,T41,G12.5,'  MW_MIX_g .LE. 0') 
 1120 FORMAT(1X,I4,T11,I4,T21,I4,T41,G12.5,'  K_g .LT. 0') 
 1130 FORMAT(1X,I4,T11,I4,T21,I4,T41,G12.5,'  C_pg .LE. 0') 
 1131 FORMAT(1X,I4,T11,I4,T21,I4,T41,G12.5,'  DIF_g .LT. 0') 
 1140 FORMAT(1X,I4,T11,I4,T21,I4,T41,G12.5,A,' .NE. 0 in a flow boundary') 
 1200 FORMAT(1X,I4,T11,I4,T21,I4,T41,G12.5,'  Sum of X_g .NE. 1') 
 1300 FORMAT(1X,I4,T11,I4,T21,I4,T31,I4,T41,G12.5,'  Sum of X_s .NE. 1') 
 1310 FORMAT(1X,I4,T11,I4,T21,I4,T31,I4,T41,G12.5,'  MU_s .LT. 0') 
 1320 FORMAT(1X,I4,T11,I4,T21,I4,T31,I4,T41,G12.5,'  K_s .LT. 0') 
 1330 FORMAT(1X,I4,T11,I4,T21,I4,T31,I4,T41,G12.5,'  C_ps .LE. 0') 
 1331 FORMAT(1X,I4,T11,I4,T21,I4,T31,I4,T41,G12.5,'  DIF_s .LT. 0') 
 1400 FORMAT(1X,I4,T11,I4,T21,I4,T41,G12.5,'  T_g .EQ. TMIN or TMAX') 
 1410 FORMAT(1X,I4,T11,I4,T21,I4,T31,I4,T41,G12.5,'  T_s .EQ. TMIN or TMAX') 
 1415 FORMAT(//1X,'Sum of all the reaction rates is not zero!',/,1X,&
         'Number of cells with discrepancy < error tolerance = ',I5,/,1X,&
         'Number of cells with discrepancy > error tolerance = ',I5,/,1X,&
         'Maximum discrepancy = ',G12.5,/,1X,&
         'Location of maximum discrepancy: I = ',I4, '  J = ', I4, '  K = ', I4&
	 ) 
 1420 FORMAT(//1X,'Production of phase ', I2, ' not equal to total mass transfer from other phases!',/,1X,&
         'Number of cells with discrepancy < error tolerance = ',I5,/,1X,&
         'Number of cells with discrepancy > error tolerance = ',I5,/,1X,&
         'Maximum discrepancy = ',G12.5,/,1X,&
         'Location of maximum discrepancy: I = ',I4, '  J = ', I4, '  K = ', I4&
	 ) 
 1430 FORMAT(//1X,'Statistics of sum of gas species mass fraction',/,1X,&
         'Sum of X_g',9X,'No of cells',2X,'Distribution') 
 1432 FORMAT(1X,'<0.9',T20,I4,T33,G12.5,/,1X,'0.9    - 0.99',T20,I4,T33,G12.5,/&
         ,1X,'0.99   - 0.999',T20,I4,T33,G12.5,/,1X,'0.999  - 0.9999',T20,I4,&
         T33,G12.5,/,1X,'0.9999 - 1.0001',T20,I4,T33,G12.5,/,1X,&
         '1.0001 - 1.001',T20,I4,T33,G12.5,/,1X,'1.001  - 1.01',T20,I4,T33,&
         G12.5,/,1X,'1.01   - 1.1',T20,I4,T33,G12.5,/,1X,'>1.1',T20,I4,T33,&
         G12.5) 
 1434 FORMAT(/1X,'Minimum sum of X_g=',G12.5,'  at I=',I4,'  J=',I4,'  K=',I4) 
 1436 FORMAT(/1X,'Maximum sum of X_g=',G12.5,'  at I=',I4,'  J=',I4,'  K=',I4) 
 1440 FORMAT(//1X,'Statistics of sum of solids (',I2,') species mass fraction',&
         /,1X,'Sum of X_s',7X,'No of cells',2X,'Distribution') 
 1442 FORMAT(1X,'<0.9',T20,I4,T33,G12.5,/,1X,'0.9    - 0.99',T20,I4,T33,G12.5,/&
         ,1X,'0.99   - 0.999',T20,I4,T33,G12.5,/,1X,'0.999  - 0.9999',T20,I4,&
         T33,G12.5,/,1X,'0.9999 - 1.0001',T20,I4,T33,G12.5,/,1X,&
         '1.0001 - 1.001',T20,I4,T33,G12.5,/,1X,'1.001  - 1.01',T20,I4,T33,&
         G12.5,/,1X,'1.01   - 1.1',T20,I4,T33,G12.5,/,1X,'>1.1',T20,I4,T33,&
         G12.5) 
 1444 FORMAT(/1X,'Minimum sum of X_s=',G12.5,'  at I=',I4,'  J=',I4,'  K=',I4) 
 1446 FORMAT(/1X,'Maximum sum of X_s=',G12.5,'  at I=',I4,'  J=',I4,'  K=',I4) 
 1500 FORMAT(/1X,70('*')/) 
      END SUBROUTINE CHECK_DATA_30 

!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization
!// 350 change do loop limits: 1,kmax2->kmin3,kmax3      
