!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CHECK_DATA_09                                          C
!  Purpose: Check chemical reactions specifications                    C
!                                                                      C
!  Author: M. Syamlal                                 Date: 10-MAR-98  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!

      SUBROUTINE CHECK_DATA_09 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE geometry
      USE fldvar
      USE physprop
      USE run
      USE rxns
      USE indices
      USE funits 
      USE compar
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
!             error flag
      LOGICAL ERROR
!
!             loop/variable indices
      INTEGER L, M, N, ID, pos, neg, mp, mn

      DOUBLE PRECISION SUM, M_m(0:DIMENSION_M)
!
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      LOGICAL , EXTERNAL :: COMPARE 
!-----------------------------------------------
!
!
!
!     Create arrays for referencing species
      ID = 0 
      DO N = 1, NMAX(0) 
         ID = ID + 1 
         SPECIES_N2IDG(N) = ID 
         SPECIES_ID2N(ID,1) = 0 
         SPECIES_ID2N(ID,2) = N 
         MW_ALL(ID) = MW_G(N) 
      END DO 
      DO M = 1, MMAX 
         DO N = 1, NMAX(M) 
            ID = ID + 1 
            SPECIES_N2IDS(M,N) = ID 
            SPECIES_ID2N(ID,1) = M 
            SPECIES_ID2N(ID,2) = N 
            MW_ALL(ID) = MW_S(M,N) 
         END DO 
      END DO 
      N_ALL = ID 
!
!
!
      DO L = 1, NO_OF_RXNS 
!
!       Is a reaction defined in data file?
         IF (.NOT.GOT_RXN(L)) THEN 
            WRITE (UNIT_LOG, 1000) L, RXN_NAME(L) 
            call mfix_exit(myPE)  
         ENDIF 
!
!       Are molecular weights and stoichiometry consistent?
         SUM = ZERO 
         IF (N_ALL > 0) THEN 
            STOICHXMW(L,:N_ALL) = STOICH(L,:N_ALL)*MW_ALL(:N_ALL) 
            DO ID = 1, N_ALL 
	       IF(STOICH(L,ID) /= ZERO .AND. MW_ALL(ID) == UNDEFINED)THEN
                 WRITE (UNIT_LOG, 1001) ID
		 call mfix_exit(myPE) 
	       ENDIF
               SUM = SUM + STOICHXMW(L,ID) 
            END DO 
         ENDIF 
         IF (.NOT.COMPARE(SUM,ZERO)) THEN 
            WRITE (UNIT_LOG, 1010) L, RXN_NAME(L) 
            call mfix_exit(myPE)  
         ENDIF 
!
         IF (GOT_RATE(L)) THEN 
!
!
!         Phase index must be between 0 and MMAX
!
            IF (RATE_M4T(L)<0 .OR. RATE_M4T(L)>MMAX) THEN 
               WRITE (UNIT_LOG, 1012) L, RXN_NAME(L), RATE_M4T(L) 
               call mfix_exit(myPE)  
            ENDIF 
!
!
!         Preexponential factor should be positive
            IF (RATE_FAC(L,1) < ZERO) THEN 
               WRITE (UNIT_LOG, 1014) L, RXN_NAME(L), RATE_FAC(L,1) 
               call mfix_exit(myPE)  
            ENDIF 
!
!
!         Activation temperature should be positive
            IF (RATE_FAC(L,3) < ZERO) THEN 
               WRITE (UNIT_LOG, 1016) L, RXN_NAME(L), RATE_FAC(L,3) 
               call mfix_exit(myPE)  
            ENDIF 
!
         ENDIF 
!
!       Determine interphase exchanges
         IF (MMAX + 1 > 0) THEN 
            R_TEMP(L,:MMAX,:MMAX) = UNDEFINED 
         ENDIF 
         M_M(0) = ZERO 
         N = 1 
         IF (NMAX(0) > 0) THEN 
            DO N = 1, NMAX(0) 
               M_M(0) = M_M(0) + STOICHXMW(L,N) 
            END DO 
            N = NMAX(0) + 1 
         ENDIF 
         DO M = 1, MMAX 
            M_M(M) = ZERO 
            DO N = 1, NMAX(M) 
               ID = SPECIES_N2IDS(M,N) 
               M_M(M) = M_M(M) + STOICHXMW(L,ID) 
            END DO 
         END DO 
         POS = 0 
         NEG = 0 
         DO M = 0, MMAX 
            IF (M_M(M) > ZERO) THEN 
               POS = POS + 1 
               MP = M 
            ELSE IF (M_M(M) < ZERO) THEN 
               NEG = NEG + 1 
               MN = M 
            ENDIF 
         END DO 
         IF (POS == 1) THEN 
            DO M = 0, MMAX 
               IF (M /= MP) R_TEMP(L,MP,M) = -M_M(M) 
            END DO 
         ELSE IF (NEG == 1) THEN 
            DO M = 0, MMAX 
               IF (M /= MN) R_TEMP(L,MN,M) = -M_M(M) 
            END DO 
         ELSE 
            IF (POS/=0 .AND. NEG/=0) THEN 
               WRITE (UNIT_LOG, 1020) L, RXN_NAME(L) 
               call mfix_exit(myPE)
	    ELSE   !no interphase transfer
	      R_TEMP(L,:MMAX,:MMAX)  = ZERO
            ENDIF 
         ENDIF 
      END DO 
      RETURN  
 1000 FORMAT(/1X,70('*')//' From: CHECK_DATA_09',/' Error: ',&
         'Reaction scheme for reaction ',&
         I2,' (',A,') not specified',/1X,70('*')/) 
 1001 FORMAT(/1X,70('*')//' From: CHECK_DATA_09',/' Error: ',&
         'Undefined molecular weight for ID: ',I2,/,&
	 '  Total NMAX may be less than the items in SPECIES_NAME.',&
	 /1X,70('*')/) 
 1010 FORMAT(/1X,70('*')//' From: CHECK_DATA_09',/' Error: ',&
         'Stoichiometry for reaction ',I2,' (',A,')',/1X,&
         'is not consistent with molecular weights',/1X,70('*')/) 
 1012 FORMAT(/1X,70('*')//' From: CHECK_DATA_09',/' Error: ',&
         'Phase index (m) for reaction ',I2,' (',A,')',/1X,&
         'is not in 0-MMAX.  m = ',I2/1X,70('*')/) 
 1014 FORMAT(/1X,70('*')//' From: CHECK_DATA_09',/' Error: ',&
         'Pre-exponential factor for reaction ',I2,' (',A,')',/1X,&
         'is less than 0.  A = ',G12.5/1X,70('*')/) 
 1016 FORMAT(/1X,70('*')//' From: CHECK_DATA_09',/' Error: ',&
         'Activation temperature for reaction ',I2,' (',A,')',/1X,&
         'is less than 0.  E/R = ',G12.5/1X,70('*')/) 
 1020 FORMAT(/1X,70('*')//' From: CHECK_DATA_09',/' Error: ',&
         'Interphase exchange for reaction ',I2,' (',A,')',/1X,&
         'cannot be determined unambiguously.',/1X,70('*')/) 
      END SUBROUTINE CHECK_DATA_09 

!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization
!// 990 Replace STOP with exitMPI to terminate all processors
