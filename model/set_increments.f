!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SET_INCREMENTS                                         C
!  Purpose: The purpose of this module is to create increments to be   C
!           stored in the array STORE_INCREMENT which will be added    C
!           to cell index ijk to find the effective indices of its     C
!           neighbors. These increments are found using the 'class'    C
!           of cell ijk. The class is determined based on the          C
!           neighboring cell type, i.e. wall or fluid.                 C
!                                                                      C
!  Author: M. Syamlal, W. Rogers                      Date: 10-DEC-91  C
!  Reviewer:M. Syamlal, S. Venkatesan, P. Nicoletti,  Date: 29-JAN-92  C
!           W. Rogers                                                  C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: Change the second index of STORE_INCREMENT to PARAMETER's  C
!  Author: M. Syamlal                                 Date: 18-FEB-92  C
!  Revision Number: 2                                                  C
!  Purpose: Do STORE_INCREMENT calculations for wall cells also.       C
!  Author: M. Syamlal                                 Date: 14-MAY-92  C
!  Revision Number: 3                                                  C
!  Purpose: Do STORE_LM calculations                                   C
!  Author: M. Syamlal                                 Date: 23-JUL-92  C
!  Revision Number: 4                                                  C
!  Purpose: Split the 2-D STORE_INCREMENT array into several 1-D       C
!           INCREMENT_FOR_xx arrays for speeding-up the program        C
!  Author: M.Syamlal                                  Date: 18-SEP-92  C
!  Reviewer: M. Syamlal                               Date: 11-DEC-92  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: IMAX2, JMAX2, KMAX2, IJKMAX2, IJKN,           C
!                        IJKS, IJKE, IJKW, IJKT, IJKB, IJKEE, IJKNE,   C
!                        IJKSE, IJKNW, IJKNN, IJKBE, IJKBN, IJKTT,     C
!                        IJKTE, IJKTW, IJKTN, IJKTS                    C
!                                                                      C
!  Variables modified: INCREMENT_FOR_xx, CELL_CLASS, DENOTE_CLASS, IJK,C
!                      I, J, K, STORE_LM                               C
!                                                                      C
!  Local variables: IC, ICLASS, DENOTE_CLASS, L                        C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE SET_INCREMENTS 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE indices
      USE geometry
      USE compar
      USE physprop
      USE fldvar
      USE funits 
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
      INTEGER          I, J, K, IJK, IMJK, IPJK, IJMK, IJPK, IJKM, IJKP,&
                       IJKW, IJKE, IJKS, IJKN, IJKB, IJKT
!                             DO-loop index, ranges from 1 to ICLASS
      INTEGER                 IC
!
!                      Index for the solids phase.
      INTEGER          M
!
!                             Local DO-loop index
      INTEGER                 L
!
!                             Array index denoting a cell class, it is a
!                             column of the array STORE_INCREMENTS
      INTEGER                 ICLASS
!
!                             Array of sum of increments to make the class
!                             determination faster. 
      INTEGER                 DENOTE_CLASS(MAX_CLASS)
!
!-----------------------------------------------
      INCLUDE 'function.inc'

!//09/02/99 Initialize the default values to Undefined_I

      IP1(:) = UNDEFINED_I
      IM1(:) = UNDEFINED_I
      JP1(:) = UNDEFINED_I
      JM1(:) = UNDEFINED_I
      KP1(:) = UNDEFINED_I
      KM1(:) = UNDEFINED_I
!
      DO I = ISTART2, IEND2
         IF (CYCLIC_X.AND.NODESI.EQ.1) THEN 
            IP1(I) = IMAP(IMAP(I)+1)
            IM1(I) = IMAP(IMAP(I)-1)
	 ELSE
            IM1(I) = MAX(ISTART3,I - 1) 
            IP1(I) = MIN(IEND3,I + 1) 
         ENDIF 
      END DO 
      DO J = JSTART2, JEND2
         IF (CYCLIC_Y.AND.NODESJ.EQ.1) THEN 
            JP1(J) = JMAP(JMAP(J)+1)
            JM1(J) = JMAP(JMAP(J)-1)
	 ELSE
            JM1(J) = MAX(JSTART3,J - 1)
            JP1(J) = MIN(JEND3,J + 1)
         ENDIF 
      END DO 
      DO K = KSTART2, KEND2
         IF (CYCLIC_Z.AND.NODESK.EQ.1.AND.DO_K) THEN 
            KP1(K) = KMAP(KMAP(K)+1)
            KM1(K) = KMAP(KMAP(K)-1)
	 ELSE
            KM1(K) = MAX(KSTART3,K - 1)
            KP1(K) = MIN(KEND3,K + 1)
         ENDIF 
      END DO 
      
      ICLASS = 0 
!
!     Loop over all cells
      DO K = KSTART3, KEND3
         DO J = JSTART3, JEND3
            DO I = ISTART3, IEND3
!
               IJK = FUNIJK(I,J,K)               !Find value of IJK 
!
!        Fill I, J, K arrays
               I_OF(IJK) = I 
               J_OF(IJK) = J 
               K_OF(IJK) = K 
!
            END DO 
         END DO 
      END DO 

!     Loop over all cells (minus the ghost layers)
      DO K = KSTART2, KEND2
         DO J = JSTART2, JEND2
            L100: DO I = ISTART2, IEND2
!
               IJK = FUNIJK(I,J,K)               !Find value of IJK
!

!          Find the the effective cell-center indices for all neighbor cells
               CALL SET_INDEX1A (I, J, K, IJK, IMJK, IPJK, IJMK, IJPK, IJKM, &
                  IJKP, IJKW, IJKE, IJKS, IJKN, IJKB, IJKT) 
!
               ICLASS = ICLASS + 1               !Increment the ICLASS counter 
               IF (ICLASS > MAX_CLASS) THEN 
                  WRITE (UNIT_LOG, 1000) MAX_CLASS 
                  CALL MFIX_EXIT 
               ENDIF 
               INCREMENT_FOR_N(ICLASS) = IJKN - IJK 
               INCREMENT_FOR_S(ICLASS) = IJKS - IJK 
               INCREMENT_FOR_E(ICLASS) = IJKE - IJK 
               INCREMENT_FOR_W(ICLASS) = IJKW - IJK 
               INCREMENT_FOR_T(ICLASS) = IJKT - IJK 
               INCREMENT_FOR_B(ICLASS) = IJKB - IJK 
               INCREMENT_FOR_IM(ICLASS) = IMJK - IJK 
               INCREMENT_FOR_IP(ICLASS) = IPJK - IJK 
               INCREMENT_FOR_JM(ICLASS) = IJMK - IJK 
               INCREMENT_FOR_JP(ICLASS) = IJPK - IJK 
               INCREMENT_FOR_KM(ICLASS) = IJKM - IJK 
               INCREMENT_FOR_KP(ICLASS) = IJKP - IJK 
!
               DENOTE_CLASS(ICLASS) = INCREMENT_FOR_N(ICLASS) + INCREMENT_FOR_S&
                  (ICLASS) + INCREMENT_FOR_E(ICLASS) + INCREMENT_FOR_W(ICLASS)&
                   + INCREMENT_FOR_T(ICLASS) + INCREMENT_FOR_B(ICLASS) + &
                  INCREMENT_FOR_IM(ICLASS) + INCREMENT_FOR_IP(ICLASS) + &
                  INCREMENT_FOR_JM(ICLASS) + INCREMENT_FOR_JP(ICLASS) + &
                  INCREMENT_FOR_KM(ICLASS) + INCREMENT_FOR_KP(ICLASS) 
!
               CELL_CLASS(IJK) = ICLASS 
!
!
!          Place the cell in a class based on its DENOTE_CLASS(ICLASS) value
               DO IC = 1, ICLASS - 1             !Loop over previous and present classes 
!                                                !IF a possible match in cell types
                  IF (DENOTE_CLASS(ICLASS) == DENOTE_CLASS(IC)) THEN 
!                                                !is found, compare all increments
                     IF (INCREMENT_FOR_N(ICLASS) /= INCREMENT_FOR_N(IC)) CYCLE  
                     IF (INCREMENT_FOR_S(ICLASS) /= INCREMENT_FOR_S(IC)) CYCLE  
                     IF (INCREMENT_FOR_E(ICLASS) /= INCREMENT_FOR_E(IC)) CYCLE  
                     IF (INCREMENT_FOR_W(ICLASS) /= INCREMENT_FOR_W(IC)) CYCLE  
                     IF (INCREMENT_FOR_T(ICLASS) /= INCREMENT_FOR_T(IC)) CYCLE  
                     IF (INCREMENT_FOR_B(ICLASS) /= INCREMENT_FOR_B(IC)) CYCLE  
                     IF (INCREMENT_FOR_IM(ICLASS) /= INCREMENT_FOR_IM(IC)) &
                        CYCLE  
                     IF (INCREMENT_FOR_IP(ICLASS) /= INCREMENT_FOR_IP(IC)) &
                        CYCLE  
                     IF (INCREMENT_FOR_JM(ICLASS) /= INCREMENT_FOR_JM(IC)) &
                        CYCLE  
                     IF (INCREMENT_FOR_JP(ICLASS) /= INCREMENT_FOR_JP(IC)) &
                        CYCLE  
                     IF (INCREMENT_FOR_KM(ICLASS) /= INCREMENT_FOR_KM(IC)) &
                        CYCLE  
                     IF (INCREMENT_FOR_KP(ICLASS) /= INCREMENT_FOR_KP(IC)) &
                        CYCLE  
                     CELL_CLASS(IJK) = IC        !Assign cell to a class 
                     ICLASS = ICLASS - 1 
                     CYCLE  L100                 !Go to next cell 
                  ENDIF 
               END DO 
            END DO L100 
         END DO 
      END DO 
      DO M = 1, MMAX 
         DO L = M, MMAX 
            IF (L == M) THEN 
               STORE_LM(L,M) = 0 
            ELSE 
               STORE_LM(L,M) = M + (L - 2)*(L - 1)/2 
               STORE_LM(M,L) = M + (L - 2)*(L - 1)/2 
            ENDIF 
         END DO 
      END DO 
      RETURN  
!
!     WRITE FOLLOWING IF THERE IS AN ERROR IN MODULE
 1000 FORMAT(/70('*')//'From: SET_INCREMENTS'/'Message: The number of',&
         'classes has exceeded the maximum allowed (',I3,').  Increase',&
         'MAX_CLASS in PARAM1.INC') 
!
      END SUBROUTINE SET_INCREMENTS 
