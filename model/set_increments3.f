!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SET_INCREMENTS3                                         C
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
      SUBROUTINE SET_INCREMENTS3 
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
!                      error Index
      INTEGER          ier
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
      LOGICAL          TRUEI, TRUEJ, TRUEK
!-----------------------------------------------
      INCLUDE 'function.inc'
      INCLUDE 'function3.inc'

!// Initialize the default values to Undefined_I

      IP1_3(:) = UNDEFINED_I
      IM1_3(:) = UNDEFINED_I
      JP1_3(:) = UNDEFINED_I
      JM1_3(:) = UNDEFINED_I
      KP1_3(:) = UNDEFINED_I
      KM1_3(:) = UNDEFINED_I
!
      DO I = ISTART4, IEND4
 	 TRUEI = .NOT.(I.EQ.IMIN4.OR.I.EQ.IMIN3.OR.I.EQ.IMIN2&
	 .OR.I.EQ.IMAX4.OR.I.EQ.IMAX3.OR.I.EQ.IMAX2)
         IF (CYCLIC_X.AND.NODESI.EQ.1.AND.DO_I.AND.TRUEI) THEN
            IP1_3(I) = IMAP_C(IMAP_C(I)+1)
            IM1_3(I) = IMAP_C(IMAP_C(I)-1)
	 ELSE
            IM1_3(I) = MAX(ISTART4,I - 1) 
            IP1_3(I) = MIN(IEND4,I + 1) 
         ENDIF 
	 
      END DO 
      DO J = JSTART4, JEND4
 	 TRUEJ = .NOT.(J.EQ.JMIN4.OR.J.EQ.JMIN3.OR.J.EQ.JMIN2&
	 .OR.J.EQ.JMAX4.OR.J.EQ.JMAX3.OR.J.EQ.JMAX2)
         IF (CYCLIC_Y.AND.NODESJ.EQ.1.AND.DO_J.AND.TRUEJ) THEN 
            JP1_3(J) = JMAP_C(JMAP_C(J)+1)
            JM1_3(J) = JMAP_C(JMAP_C(J)-1)
	 ELSE
            JM1_3(J) = MAX(JSTART4,J - 1)
            JP1_3(J) = MIN(JEND4,J + 1)
         ENDIF 
      END DO 
      DO K = KSTART4, KEND4
 	 TRUEK = .NOT.(K.EQ.KMIN4.OR.K.EQ.KMIN3.OR.K.EQ.KMIN2&
	 .OR.K.EQ.KMAX4.OR.K.EQ.KMAX3.OR.K.EQ.KMAX2)
         IF (CYCLIC_Z.AND.NODESK.EQ.1.AND.DO_K.AND.TRUEK) THEN
            KP1_3(K) = KMAP_C(KMAP_C(K)+1)
            KM1_3(K) = KMAP_C(KMAP_C(K)-1)
	 ELSE
            KM1_3(K) = MAX(KSTART4,K - 1)
            KP1_3(K) = MIN(KEND4,K + 1)
         ENDIF 
      END DO 
      
      ICLASS = 0 
!
!     Loop over all cells
      DO K = KSTART4, KEND4
         DO J = JSTART4, JEND4
            DO I = ISTART4, IEND4
               IJK = FUNIJK3(I,J,K)               !Find value of IJK 
!
!        Fill I, J, K arrays
               I3_OF(IJK) = I
               J3_OF(IJK) = J 
               K3_OF(IJK) = K 
            END DO 
         END DO 
      END DO 


!     Loop over all cells (minus the ghost layers)
      DO K = KSTART4, KEND4
         DO J = JSTART4, JEND4
            L100: DO I = ISTART4, IEND4
               IJK = FUNIJK3(I,J,K)               !Find value of IJK
!

!          Find the the effective cell-center indices for all neighbor cells
               CALL SET_INDEX1A3 (I, J, K, IJK, IMJK, IPJK, IJMK, IJPK, IJKM, &
                  IJKP, IJKW, IJKE, IJKS, IJKN, IJKB, IJKT) 
!
               ICLASS = ICLASS + 1               !Increment the ICLASS counter 
               IF (ICLASS > MAX_CLASS) THEN 
	          IF(.not.DMP_LOG)call open_pe_log(ier)
                  WRITE (UNIT_LOG, 1000) MAX_CLASS 
                  CALL MFIX_EXIT(myPE) 
               ENDIF 
               INCREMENT3_FOR_IM(ICLASS) = IMJK - IJK 
               INCREMENT3_FOR_IP(ICLASS) = IPJK - IJK 
               INCREMENT3_FOR_JM(ICLASS) = IJMK - IJK 
               INCREMENT3_FOR_JP(ICLASS) = IJPK - IJK 
               INCREMENT3_FOR_KM(ICLASS) = IJKM - IJK 
               INCREMENT3_FOR_KP(ICLASS) = IJKP - IJK 
!
               DENOTE_CLASS(ICLASS) =  &
                  INCREMENT3_FOR_IM(ICLASS) + INCREMENT3_FOR_IP(ICLASS) + &
                  INCREMENT3_FOR_JM(ICLASS) + INCREMENT3_FOR_JP(ICLASS) + &
                  INCREMENT3_FOR_KM(ICLASS) + INCREMENT3_FOR_KP(ICLASS) 
!
               CELL_CLASS3(IJK) = ICLASS 
!
!
!          Place the cell in a class based on its DENOTE_CLASS(ICLASS) value
               DO IC = 1, ICLASS - 1             !Loop over previous and present classes 
!                                                !IF a possible match in cell types
                  IF (DENOTE_CLASS(ICLASS) == DENOTE_CLASS(IC)) THEN 
!                                                !is found, compare all increments
                     IF (INCREMENT3_FOR_IM(ICLASS) /= INCREMENT3_FOR_IM(IC)) &
                        CYCLE  
                     IF (INCREMENT3_FOR_IP(ICLASS) /= INCREMENT3_FOR_IP(IC)) &
                        CYCLE  
                     IF (INCREMENT3_FOR_JM(ICLASS) /= INCREMENT3_FOR_JM(IC)) &
                        CYCLE  
                     IF (INCREMENT3_FOR_JP(ICLASS) /= INCREMENT3_FOR_JP(IC)) &
                        CYCLE  
                     IF (INCREMENT3_FOR_KM(ICLASS) /= INCREMENT3_FOR_KM(IC)) &
                        CYCLE  
                     IF (INCREMENT3_FOR_KP(ICLASS) /= INCREMENT3_FOR_KP(IC)) &
                        CYCLE  
                     CELL_CLASS3(IJK) = IC        !Assign cell to a class 
                     ICLASS = ICLASS - 1 
                     CYCLE  L100                 !Go to next cell 
                  ENDIF 
               END DO 
            END DO L100 
         END DO 
      END DO 
      RETURN  
!
!     WRITE FOLLOWING IF THERE IS AN ERROR IN MODULE
 1000 FORMAT(/70('*')//'From: SET_INCREMENTS3'/'Message: The number of',&
         'classes has exceeded the maximum allowed (',I3,').  Increase',&
         'MAX_CLASS in PARAM1.INC') 
!
      END SUBROUTINE SET_INCREMENTS3 

!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization
!// 020 New local variables for parallelization : TRUEI, TRUEJ, TRUEK
!// 350 Changed do loop limits: 1,kmax2->kmin3,kmax3      
