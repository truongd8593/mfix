!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CHECK_DATA_08                                          C
!  Purpose: Check internal surface specifications, and convert         C
!           physical locations to i, j, k's.                           C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-OCT-92  C
!  Reviewer: W. Rogers                                Date: 11-DEC-92  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: IS_X_w, IS_X_e, IS_Y_s, IS_Y_n, IS_Z_b        C
!                        IS_Z_t, IS_I_w, IS_I_e, IS_J_s, IS_J_n        C
!                        IS_K_b, IS_K_t, IMAX2, JMAX2, KMAX2           C
!  Variables modified: IS_DEFINED                                      C
!                                                                      C
!  Local variables: ISV, I, J, K, IJK, VALID_IS_TYPE                  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CHECK_DATA_08 
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
      USE is
      USE indices
      USE funits 
      USE compar   !//d
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: DIM_ISTYPE = 4 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!
!             error flag
      LOGICAL ERROR
!
!             loop/variable indices
      INTEGER ISV , I , J , K , IJK
!
!             valid internal surface types
      CHARACTER*16, DIMENSION(1:DIM_ISTYPE) ::VALID_IS_TYPE = (/&
           'IMPERMEABLE     ', 'IP              ',&
	   'SEMIPERMEABLE   ', 'SP              '&
	    /)
      DOUBLE PRECISION SUM
!-----------------------------------------------
      INCLUDE 'function.inc'
      
!//? DID NOT PERFORMED A TEST WITH INTERNAL SURF. YET

!
! DETERMINE WHICH INTERNAL SURFACE INDICES HAVE VALUES
!
      L50: DO ISV = 1, DIMENSION_IS 
         IS_DEFINED(ISV) = .FALSE. 
         IF (IS_X_W(ISV) /= UNDEFINED) IS_DEFINED(ISV) = .TRUE. 
         IF (IS_X_E(ISV) /= UNDEFINED) IS_DEFINED(ISV) = .TRUE. 
         IF (IS_Y_S(ISV) /= UNDEFINED) IS_DEFINED(ISV) = .TRUE. 
         IF (IS_Y_N(ISV) /= UNDEFINED) IS_DEFINED(ISV) = .TRUE. 
         IF (IS_Z_B(ISV) /= UNDEFINED) IS_DEFINED(ISV) = .TRUE. 
         IF (IS_Z_T(ISV) /= UNDEFINED) IS_DEFINED(ISV) = .TRUE. 
         IF (IS_I_W(ISV) /= UNDEFINED_I) IS_DEFINED(ISV) = .TRUE. 
         IF (IS_I_E(ISV) /= UNDEFINED_I) IS_DEFINED(ISV) = .TRUE. 
         IF (IS_J_S(ISV) /= UNDEFINED_I) IS_DEFINED(ISV) = .TRUE. 
         IF (IS_J_N(ISV) /= UNDEFINED_I) IS_DEFINED(ISV) = .TRUE. 
         IF (IS_K_B(ISV) /= UNDEFINED_I) IS_DEFINED(ISV) = .TRUE. 
         IF (IS_K_T(ISV) /= UNDEFINED_I) IS_DEFINED(ISV) = .TRUE. 
         IF (IS_DEFINED(ISV)) THEN 
            IF (IS_X_W(ISV)==UNDEFINED .AND. IS_I_W(ISV)==UNDEFINED_I) THEN 
               IF (NO_I) THEN 
                  IS_X_W(ISV) = ZERO 
               ELSE 
                  WRITE (UNIT_LOG, 1000) 'IS_X_w and IS_I_w ', ISV 
                  call mfix_exit(myPE)  
               ENDIF 
            ENDIF 
            IF (IS_X_E(ISV)==UNDEFINED .AND. IS_I_E(ISV)==UNDEFINED_I) THEN 
               IF (NO_I) THEN 
                  IS_X_E(ISV) = XLENGTH 
               ELSE 
                  WRITE (UNIT_LOG, 1000) 'IS_X_e and IS_I_e ', ISV 
                  call mfix_exit(myPE)		  
               ENDIF 
            ENDIF 
            IF (IS_Y_S(ISV)==UNDEFINED .AND. IS_J_S(ISV)==UNDEFINED_I) THEN 
               IF (NO_J) THEN 
                  IS_Y_S(ISV) = ZERO 
               ELSE 
                  WRITE (UNIT_LOG, 1000) 'IS_Y_s and IS_J_s ', ISV 
                  call mfix_exit(myPE)		  
               ENDIF 
            ENDIF 
            IF (IS_Y_N(ISV)==UNDEFINED .AND. IS_J_N(ISV)==UNDEFINED_I) THEN 
               IF (NO_J) THEN 
                  IS_Y_N(ISV) = YLENGTH 
               ELSE 
                  WRITE (UNIT_LOG, 1000) 'IS_Y_n and IS_J_n ', ISV 
                  call mfix_exit(myPE)
               ENDIF 
            ENDIF 
            IF (IS_Z_B(ISV)==UNDEFINED .AND. IS_K_B(ISV)==UNDEFINED_I) THEN 
               IF (NO_K) THEN 
                  IS_Z_B(ISV) = ZERO 
               ELSE 
                  WRITE (UNIT_LOG, 1000) 'IS_Z_b and IS_K_b ', ISV 
                  call mfix_exit(myPE)
               ENDIF 
            ENDIF 
            IF (IS_Z_T(ISV)==UNDEFINED .AND. IS_K_T(ISV)==UNDEFINED_I) THEN 
               IF (NO_K) THEN 
                  IS_Z_T(ISV) = ZLENGTH 
               ELSE 
                  WRITE (UNIT_LOG, 1000) 'IS_Z_t and IS_K_t ', ISV 
                  call mfix_exit(myPE)
               ENDIF 
            ENDIF 
            DO I = 1, DIM_ISTYPE 
               IF (VALID_IS_TYPE(I) == IS_TYPE(ISV)) THEN 
                  IF (MOD(I,2) == 0) IS_TYPE(ISV) = VALID_IS_TYPE(I-1) 
                  CYCLE  L50 
               ENDIF 
               IF (VALID_IS_TYPE(I) == IS_TYPE(ISV)(3:16)) THEN 
                  IF (MOD(I,2) == 0) IS_TYPE(ISV)(3:16) = VALID_IS_TYPE(I-1) 
                  IF (IS_TYPE(ISV)(1:1)/='X' .AND. IS_TYPE(ISV)(1:1)/='Y' .AND. &
                     IS_TYPE(ISV)(1:1)/='Z') THEN 
                     WRITE (UNIT_LOG, 1000) ISV, IS_TYPE(ISV)(1:1) 
                     CALL MFIX_EXIT 
                  ENDIF 
                  CYCLE  L50 
               ENDIF 
            END DO 
            WRITE (UNIT_LOG, 1001) ISV, IS_TYPE(ISV) 
            WRITE (UNIT_LOG, 1002) VALID_IS_TYPE 
            call mfix_exit(myPE)  
         ENDIF 
      END DO L50 
      CALL GET_IS 
!
      DO ISV = 1, DIMENSION_IS 
         IF (IS_DEFINED(ISV)) THEN 
!
!  Check whether IS_PC is defined
!
            IF (IS_TYPE(ISV)=='SEMIPERMEABLE' .OR. IS_TYPE(ISV)(3:15)==&
               'SEMIPERMEABLE') THEN 
               IF (IS_PC(ISV,1) == UNDEFINED) THEN 
                  WRITE (UNIT_LOG, 1005) 'IS_PC', ISV 
                  call mfix_exit(myPE)  
               ENDIF 
               IF (IS_PC(ISV,1) == ZERO) THEN 
                  WRITE (UNIT_LOG, 1006) 'IS_PC', ISV 
                  call mfix_exit(myPE)  
               ENDIF 
               IF (IS_PC(ISV,2) == UNDEFINED) THEN 
                  WRITE (UNIT_LOG, 1010) 'IS_PC', ISV 
                  call mfix_exit(myPE)  
               ENDIF 
            ENDIF 
         ELSE 
!
! make sure IS_PC is not specified for undefined IS
!
            IF (IS_PC(ISV,1) /= LARGE_NUMBER) THEN 
               WRITE (UNIT_LOG, 1200) 'IS_PC', ISV 
               call mfix_exit(myPE)  
            ENDIF 
            IF (IS_PC(ISV,2) /= ZERO) THEN 
               WRITE (UNIT_LOG, 1210) 'IS_PC', ISV 
               call mfix_exit(myPE)  
            ENDIF 
         ENDIF 
      END DO 
      RETURN  
 1000 FORMAT(/1X,70('*')//' From: CHECK_DATA_08',/&
         ' Message: The IS_TYPE prefix should be X, Y, or Z for IS # = ',I2,/&
         '   Direction prefix of IS_TYPE = ',A) 
 1001 FORMAT(/1X,70('*')//' From: CHECK_DATA_08',/&
         ' Message: Illegal IS_TYPE for IS # = ',I2,/'   IS_TYPE = ',A,/&
         '  Valid IS_TYPE are: ') 
 1002 FORMAT(5X,A16) 
 1005 FORMAT(/1X,70('*')//' From: CHECK_DATA_08',/' Message: ',A,'(',I2,&
         ',1) not specified',/1X,70('*')/) 
 1006 FORMAT(/1X,70('*')//' From: CHECK_DATA_08',/' Message: ',A,'(',I2,&
         ',1) should not be zero',/1X,70('*')/) 
 1010 FORMAT(/1X,70('*')//' From: CHECK_DATA_08',/' Message: ',A,'(',I2,&
         ',2) not specified',/1X,70('*')/) 
 1200 FORMAT(/1X,70('*')//' From: CHECK_DATA_08',/' Message: ',A,'(',I2,&
         ',1) specified',' for an undefined IS location',/1X,70('*')/) 
 1210 FORMAT(/1X,70('*')//' From: CHECK_DATA_08',/' Message: ',A,'(',I2,&
         ',2) specified',' for an undefined IS location',/1X,70('*')/) 
      END SUBROUTINE CHECK_DATA_08 
