CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: SET_INCREMENTS                                         C
C  Purpose: The purpose of this module is to create increments to be   C
C           stored in the array STORE_INCREMENT which will be added    C
C           to cell index ijk to find the effective indices of its     C
C           neighbors. These increments are found using the 'class'    C
C           of cell ijk. The class is determined based on the          C
C           neighboring cell type, i.e. wall or fluid.                 C
C                                                                      C
C  Author: M. Syamlal, W. Rogers                      Date: 10-DEC-91  C
C  Reviewer:M. Syamlal, S. Venkatesan, P. Nicoletti,  Date: 29-JAN-92  C
C           W. Rogers                                                  C
C                                                                      C
C  Revision Number: 1                                                  C
C  Purpose: Change the second index of STORE_INCREMENT to PARAMETER's  C
C  Author: M. Syamlal                                 Date: 18-FEB-92  C
C  Revision Number: 2                                                  C
C  Purpose: Do STORE_INCREMENT calculations for wall cells also.       C
C  Author: M. Syamlal                                 Date: 14-MAY-92  C
C  Revision Number: 3                                                  C
C  Purpose: Do STORE_LM calculations                                   C
C  Author: M. Syamlal                                 Date: 23-JUL-92  C
C  Revision Number: 4                                                  C
C  Purpose: Split the 2-D STORE_INCREMENT array into several 1-D       C
C           INCREMENT_FOR_xx arrays for speeding-up the program        C
C  Author: M.Syamlal                                  Date: 18-SEP-92  C
C  Reviewer: M. Syamlal                               Date: 11-DEC-92  C
C                                                                      C
C  Literature/Document References:                                     C
C                                                                      C
C  Variables referenced: IMAX2, JMAX2, KMAX2, IJKMAX2, IJKN,           C
C                        IJKS, IJKE, IJKW, IJKT, IJKB, IJKEE, IJKNE,   C
C                        IJKSE, IJKNW, IJKNN, IJKBE, IJKBN, IJKTT,     C
C                        IJKTE, IJKTW, IJKTN, IJKTS                    C
C                                                                      C
C  Variables modified: INCREMENT_FOR_xx, CELL_CLASS, DENOTE_CLASS, IJK,C
C                      I, J, K, STORE_LM                               C
C                                                                      C
C  Local variables: IC, ICLASS, DENOTE_CLASS, L                        C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
C
      SUBROUTINE SET_INCREMENTS
C
      IMPLICIT NONE
C
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE 'indices.inc'
      INCLUDE 'geometry.inc'
      INCLUDE 'physprop.inc'
      INCLUDE 'fldvar.inc'
      INCLUDE 'funits.inc'
C
C Local Variables
C
C                      Indices
      INTEGER          I, J, K, IJK, IMJK, IPJK, IJMK, IJPK, IJKM, IJKP,
     &                 IJKW, IJKE, IJKS, IJKN, IJKB, IJKT
C                             DO-loop index, ranges from 1 to ICLASS
      INTEGER                 IC
C
C                      Index for the solids phase.
      INTEGER          M
C
C                             Local DO-loop index
      INTEGER                 L
C
C                             Array index denoting a cell class, it is a
C                             column of the array STORE_INCREMENTS
      INTEGER                 ICLASS
C
C                             Array of sum of increments to make the class
C                             determination faster. 
      INTEGER                 DENOTE_CLASS(MAX_CLASS)
C
C
      INCLUDE 'function.inc'
C
C
      DO 10 I = 1, IMAX2
        Im1(I) = MAX(1, (I-1))
        Ip1(I) = MIN(IMAX2, (I+1))
        IF(CYCLIC_X) THEN
            IF(I .EQ. IMAX1)Ip1(I) = IMIN1
            IF(I .EQ. 2)Im1(I) = IMAX1
        ENDIF
10    CONTINUE
C
      DO 20 J = 1, JMAX2
        Jm1(J) = MAX(1, (J-1))
        Jp1(J) = MIN(JMAX2, (J+1))
        IF(CYCLIC_Y) THEN
            IF(J .EQ. JMAX1)Jp1(J) = JMIN1
            IF(J .EQ. 2)Jm1(J) = JMAX1
        ENDIF
20    CONTINUE
C
      DO 30 K = 1, KMAX2
        Km1(K) = MAX(1, (K-1))
        Kp1(K) = MIN(KMAX2, (K+1))
        IF(CYCLIC_Z) THEN
            IF(K .EQ. KMAX1)Kp1(K) = KMIN1
            IF(K .EQ. 2)Km1(K) = KMAX1
        ENDIF
30    CONTINUE
C
C     Initialize the ICLASS flag
      ICLASS = 0
C
C     Loop over all cells
      DO 110 K = 1, KMAX2
      DO 105 J = 1, JMAX2
      DO 100 I = 1, IMAX2
C
         IJK = FUNIJK(I, J, K)  !Find value of IJK
C
C        Fill I, J, K arrays
         I_OF(IJK) = I
         J_OF(IJK) = J
         K_OF(IJK) = K
C
C          Find the the effective cell-center indices for all neighbor cells  
           CALL SET_INDEX1A(I, J, K, IJK, IMJK, IPJK, IJMK, IJPK,
     &                       IJKM, IJKP, IJKW, IJKE, IJKS, IJKN,
     &                       IJKB, IJKT) 
C
           ICLASS = ICLASS + 1           !Increment the ICLASS counter
           IF(ICLASS .GT. MAX_CLASS) THEN
             WRITE(UNIT_LOG, 1000) MAX_CLASS
             CALL EXIT
           ENDIF
           INCREMENT_FOR_n (ICLASS) = IJKN - IJK
           INCREMENT_FOR_s (ICLASS) = IJKS - IJK  
           INCREMENT_FOR_e (ICLASS) = IJKE - IJK
           INCREMENT_FOR_w (ICLASS) = IJKW - IJK
           INCREMENT_FOR_t (ICLASS) = IJKT - IJK
           INCREMENT_FOR_b (ICLASS) = IJKB - IJK
           INCREMENT_FOR_im(ICLASS) = IMJK - IJK
           INCREMENT_FOR_ip(ICLASS) = IPJK - IJK
           INCREMENT_FOR_jm(ICLASS) = IJMK - IJK
           INCREMENT_FOR_jp(ICLASS) = IJPK - IJK
           INCREMENT_FOR_km(ICLASS) = IJKM - IJK
           INCREMENT_FOR_kp(ICLASS) = IJKP - IJK
C
           DENOTE_CLASS(ICLASS) = 
     &       INCREMENT_FOR_n (ICLASS) +
     &       INCREMENT_FOR_s(ICLASS) + 
     &       INCREMENT_FOR_e(ICLASS) + 
     &       INCREMENT_FOR_w (ICLASS) +
     &       INCREMENT_FOR_t (ICLASS) +
     &       INCREMENT_FOR_b (ICLASS) +
     &       INCREMENT_FOR_im(ICLASS) +
     &       INCREMENT_FOR_ip(ICLASS) +
     &       INCREMENT_FOR_jm(ICLASS) +
     &       INCREMENT_FOR_jp(ICLASS) +
     &       INCREMENT_FOR_km(ICLASS) +
     &       INCREMENT_FOR_kp(ICLASS)
C
           CELL_CLASS(IJK) = ICLASS
C
C
C          Place the cell in a class based on its DENOTE_CLASS(ICLASS) value
           DO 50 IC = 1, ICLASS-1      !Loop over previous and present classes
             IF (DENOTE_CLASS(ICLASS)    
     &           .EQ. DENOTE_CLASS(IC)) THEN !IF a possible match in cell types 
               IF(INCREMENT_FOR_n (ICLASS) .NE.!is found, compare all increments 
     &            INCREMENT_FOR_n (IC))GOTO 50
               IF(INCREMENT_FOR_s(ICLASS) .NE.
     &            INCREMENT_FOR_s(IC))GOTO 50 
               IF(INCREMENT_FOR_e(ICLASS) .NE.
     &            INCREMENT_FOR_e(IC))GOTO 50 
               IF(INCREMENT_FOR_w (ICLASS) .NE.
     &           INCREMENT_FOR_w (IC))GOTO 50
               IF(INCREMENT_FOR_t (ICLASS) .NE.
     &           INCREMENT_FOR_t (IC))GOTO 50
               IF(INCREMENT_FOR_b (ICLASS) .NE.
     &           INCREMENT_FOR_b (IC))GOTO 50
               IF(INCREMENT_FOR_im(ICLASS) .NE.
     &           INCREMENT_FOR_im(IC))GOTO 50
               IF(INCREMENT_FOR_ip(ICLASS) .NE.
     &           INCREMENT_FOR_ip(IC))GOTO 50
               IF(INCREMENT_FOR_jm(ICLASS) .NE.
     &           INCREMENT_FOR_jm(IC))GOTO 50
               IF(INCREMENT_FOR_jp(ICLASS) .NE.
     &           INCREMENT_FOR_jp(IC))GOTO 50
               IF(INCREMENT_FOR_km(ICLASS) .NE.
     &           INCREMENT_FOR_km(IC))GOTO 50
               IF(INCREMENT_FOR_kp(ICLASS) .NE.
     &           INCREMENT_FOR_kp(IC))GOTO 50
               CELL_CLASS(IJK) = IC        !Assign cell to a class
               ICLASS = ICLASS - 1
               GOTO 100                    !Go to next cell
             ENDIF
50         CONTINUE
100   CONTINUE
105   CONTINUE
110   CONTINUE
C
C  Do STORE_LM calculations
C
      DO 210 M = 1, MMAX
      DO 200 L = M, MMAX
        IF(L .EQ. M) THEN
          STORE_LM(L, M) = 0
        ELSE
          STORE_LM(L, M) = M + (L-2) * (L-1) / 2
          STORE_LM(M, L) = M + (L-2) * (L-1) / 2
        ENDIF
200   CONTINUE
210   CONTINUE
      RETURN
C
C     WRITE FOLLOWING IF THERE IS AN ERROR IN MODULE
1000  FORMAT(/70('*')//'From: SET_INCREMENTS'/'Message: The number of',
     & 'classes has exceeded the maximum allowed (',I3,').  Increase',
     & 'MAX_CLASS in PARAM.INC')
C
      END
