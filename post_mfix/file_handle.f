CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: IN_BIN_512                                             C
C  Purpose: read an array in chunks of 512 bytes    (DP WORDS)         C
C                                                                      C
C  Author: P. Nicoletti                               Date: 02-JAN-92  C
C  Reviewer: P. Nicoletti, W. Rogers, M. Syamlal      Date: 24-JAN-92  C
C                                                                      C
C  Revision Number:                                                    C
C  Purpose:                                                            C
C  Author:                                            Date: dd-mmm-yy  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Literature/Document References:                                     C
C                                                                      C
C  Variables referenced:                                               C
C  Variables modified:                                                 C
C                                                                      C
C  Local variables: NWORDS, DS, L, NSEG, NREM, LC, N1, N2              C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
C
      SUBROUTINE IN_BIN_512(IUNIT,ARRAY,N,NEXT_REC)
C
      IMPLICIT NONE
C
      INCLUDE 'machine.inc'
C
C passed arguments
C
C                      array to write out
      DOUBLE PRECISION ARRAY(*)
C
C                      output unit number
      INTEGER          IUNIT
C
C                      number of elements in ARRAY
      INTEGER          N
C
C                      next record number in direct access output file
      INTEGER          NEXT_REC
C
C local variables
C
C                      number of words for 512 bytes
      INTEGER          NWORDS
C
C                      loop counter
      INTEGER          L
C
C                      number of full 512 byte segments need to write N
C                      double precision words
      INTEGER          NSEG
C
C                      number of double precision words in the partially
C                      filled last record
      INTEGER          NREM
C
C                      loop counter
      INTEGER          LC
C
C                      write out array elements N1 to N2
      INTEGER          N1 , N2
C
      NWORDS = NWORDS_DP
      IF (N.LE.NWORDS) THEN
         READ (IUNIT,REC=NEXT_REC) (ARRAY(L),L=1,N)
         NEXT_REC = NEXT_REC + 1
         RETURN
      END IF
C
      NSEG = N / NWORDS
      NREM = MOD(N,NWORDS)
      N1 = 1
      N2 = NWORDS
C
C read the full 512 byte segments
C
      DO 100 LC = 1,NSEG
         READ (IUNIT,REC=NEXT_REC) (ARRAY(L),L=N1,N2)
         N1 = N1 + NWORDS
         N2 = N2 + NWORDS
         NEXT_REC = NEXT_REC + 1
100   CONTINUE
C
C read the partially filled last record
C
      IF (NREM.NE.0) THEN
         READ (IUNIT,REC=NEXT_REC) (ARRAY(L),L=N1,N)
         NEXT_REC = NEXT_REC + 1
      END IF
C
      RETURN
      END
CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: IN_BIN_512I                                            C
C  Purpose: read an array in chunks of 512 bytes (INTEGER WORDS)       C
C                                                                      C
C  Author: P. Nicoletti                               Date: 02-JAN-92  C
C  Reviewer: P. Nicoletti, W. Rogers, M. Syamlal      Date: 24-JAN-92  C
C                                                                      C
C  Revision Number:                                                    C
C  Purpose:                                                            C
C  Author:                                            Date: dd-mmm-yy  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Literature/Document References:                                     C
C                                                                      C
C  Variables referenced:                                               C
C  Variables modified:                                                 C
C                                                                      C
C  Local variables: NWORDS, L, NSEG, NREM, LC, N1, N2                  C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
C
      SUBROUTINE IN_BIN_512I(IUNIT,ARRAY,N,NEXT_REC)
C
      IMPLICIT NONE
C
      INCLUDE 'machine.inc'
C
C passed arguments
C
C                      array to write out
      INTEGER          ARRAY(*)
C
C                      output unit number
      INTEGER          IUNIT
C
C                      number of elements in ARRAY
      INTEGER          N
C
C                      next record number in direct access output file
      INTEGER          NEXT_REC
C
C local variables
C
C                      number of words for 512 bytes (nwords * 4 = 512)
      INTEGER          NWORDS
C
C                      loop counter
      INTEGER          L
C
C                      number of full 512 byte segments need to write N
C                      double precision words
      INTEGER          NSEG
C
C                      number of double precision words in the partially
C                      filled last record
      INTEGER          NREM
C
C                      loop counter
      INTEGER          LC
C
C                      write out array elements N1 to N2
      INTEGER          N1 , N2
C
      NWORDS = NWORDS_I
      IF (N.LE.NWORDS) THEN
         READ (IUNIT,REC=NEXT_REC) (ARRAY(L),L=1,N)
         NEXT_REC = NEXT_REC + 1
         RETURN
      END IF
C
      NSEG = N / NWORDS
      NREM = MOD(N,NWORDS)
      N1 = 1
      N2 = NWORDS
C
C read the full 512 byte segments
C
      DO 100 LC = 1,NSEG
         READ (IUNIT,REC=NEXT_REC) (ARRAY(L),L=N1,N2)
         N1 = N1 + NWORDS
         N2 = N2 + NWORDS
         NEXT_REC = NEXT_REC + 1
100   CONTINUE
C
C read the partially filled last record
C
      IF (NREM.NE.0) THEN
         READ (IUNIT,REC=NEXT_REC) (ARRAY(L),L=N1,N)
         NEXT_REC = NEXT_REC + 1
      END IF
C
      RETURN
      END
CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: IN_BIN_512R                                            C
C  Purpose: read in an array in chunks of 512 bytes (REAL    WORDS)    C
C                                                                      C
C  Author: P. Nicoletti                               Date: 02-JAN-92  C
C  Reviewer: P. Nicoletti, W. Rogers, M. Syamlal      Date: 24-JAN-92  C
C                                                                      C
C  Revision Number:                                                    C
C  Purpose:                                                            C
C  Author:                                            Date: dd-mmm-yy  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Literature/Document References:                                     C
C                                                                      C
C  Variables referenced:                                               C
C  Variables modified:                                                 C
C                                                                      C
C  Local variables: NWORDS, L, NSEG, NREM, LC, N1, N2                  C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
C
      SUBROUTINE IN_BIN_512R(IUNIT,ARRAY,N,NEXT_REC)
C
      IMPLICIT NONE
C
      INCLUDE 'machine.inc'
C
C passed arguments
C
C                      array to write out
      REAL             ARRAY(*)
C
C                      output unit number
      INTEGER          IUNIT
C
C                      number of elements in ARRAY
      INTEGER          N
C
C                      next record number in direct access output file
      INTEGER          NEXT_REC
C
C local variables
C
C                      number of words for 512 bytes (nwords * 4 = 512)
      INTEGER          NWORDS
C
C                      loop counter
      INTEGER          L
C
C                      number of full 512 byte segments need to write N
C                      double precision words
      INTEGER          NSEG
C
C                      number of double precision words in the partially
C                      filled last record
      INTEGER          NREM
C
C                      loop counter
      INTEGER          LC
C
C                      write out array elements N1 to N2
      INTEGER          N1 , N2
C
      NWORDS = NWORDS_R
      IF (N.LE.NWORDS) THEN
         READ (IUNIT,REC=NEXT_REC) (ARRAY(L),L=1,N)
         NEXT_REC = NEXT_REC + 1
         RETURN
      END IF
C
      NSEG = N / NWORDS
      NREM = MOD(N,NWORDS)
      N1 = 1
      N2 = NWORDS
C
C write out the full 512 byte segments
C
      DO 100 LC = 1,NSEG
         READ (IUNIT,REC=NEXT_REC) (ARRAY(L),L=N1,N2)
         N1 = N1 + NWORDS
         N2 = N2 + NWORDS
         NEXT_REC = NEXT_REC + 1
100   CONTINUE
C
C write out the partially filled last record
C
      IF (NREM.NE.0) THEN
         READ (IUNIT,REC=NEXT_REC) (ARRAY(L),L=N1,N)
         NEXT_REC = NEXT_REC + 1
      END IF
C
      RETURN
      END
CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: IN_BIN_R(IUNIT,ARRAY,IJKMAX2,NEXT_REC)                 C
C  Purpose: read in a time-dependent restart variable (REAL)           C
C                                                                      C
C  Author: P. Nicoletti                               Date: 13-DEC-91  C
C  Reviewer: P. Nicoletti, W. Rogers, M. Syamlal      Date: 24-JAN-92  C
C                                                                      C
C  Revision Number:                                                    C
C  Purpose:                                                            C
C  Author:                                            Date: dd-mmm-yy  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Literature/Document References:                                     C
C                                                                      C
C  Variables referenced:                                               C
C  Variables modified:                                                 C
C                                                                      C
C  Local variables: LC, ARRAY_REAL                                     C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
C
      SUBROUTINE IN_BIN_R(IUNIT,ARRAY,IJKMAX2,NEXT_REC)
C
      IMPLICIT NONE
      INCLUDE  'param.inc'
      INCLUDE 'param1.inc'
C
C passed arguments
C
C                      double precision array to write out
      DOUBLE PRECISION ARRAY(*)
C
C                      unit number to write to
      INTEGER          IUNIT
C
C                      record pointer into file IUNIT
      INTEGER          NEXT_REC
C
C                      number of indices in ARRAY to write out
      INTEGER          IJKMAX2
C
C local variables
C
C                      single precision version of ARRAY
      REAL             ARRAY_REAL(DIMENSION_3)
C
C                      loop counter
      INTEGER          LC
C
      CALL IN_BIN_512R (IUNIT,ARRAY_REAL,IJKMAX2,NEXT_REC)
C
      DO 100 LC = 1,IJKMAX2
         ARRAY(LC) = ARRAY_REAL(LC)
100   CONTINUE
C
      RETURN
      END
CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: OPEN_FILEP(RUN_NAME, RUN_TYPE, NO_FILES)               C
C  Purpose: open all the files for this run (modified for POST_MFIX)   C
C                                                                      C
C  Author: P. Nicoletti                               Date: 12-DEC-91  C
C  Reviewer: P. Nicoletti, W. Rogers, M. Syamlal      Date: 24-JAN-92  C
C                                                                      C
C  Revision Number:                                                    C
C  Purpose:                                                            C
C  Author:                                            Date: dd-mmm-yy  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Literature/Document References:                                     C
C                                                                      C
C  Variables referenced:                                               C
C  Variables modified:                                                 C
C                                                                      C
C  Local variables: EXT, FILE_NAME, LC, NB                             C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
C
      LOGICAL FUNCTION OPEN_FILEP(RUN_NAME, RUN_TYPE, NO_FILES)
C
      IMPLICIT NONE
C
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE 'post3d.inc'
      INCLUDE 'machine.inc'
      INCLUDE 'funits.inc'
C
C passed arguments
C
C                   run_name (as specified in input file)
      CHARACTER*(*) RUN_NAME
C
C                   run type
      CHARACTER*(*) RUN_TYPE
C
C local variables
C
C                   extension to filename
      CHARACTER     EXT*4
C
C                   run_name + extension
      CHARACTER     FILE_NAME*64
C
C                   Loop counter
      INTEGER       LC
C
C                   index to first blank character in run_name
      INTEGER       NB
C
C                   Number of files to be opened
      INTEGER       NO_FILES
C
      OPEN_FILEP = .FALSE.
C
C DETERMINE THE FIRST BLANK CHARCATER IN RUN_NAME
C
      DO 100 LC = 1,LEN(RUN_NAME)
         IF (RUN_NAME(LC:LC).EQ.' ') THEN
            NB = LC
            GOTO 125
         END IF
100   CONTINUE
      WRITE (*,*) 'RUN_NAME TOOOOOOO LOOOONG'
      RETURN
C
125   IF (NB+4.GT.LEN(FILE_NAME)) THEN
         WRITE (*,*) 'RUN_NAME TOOOOOOO LOOOONG'
         RETURN
      END IF
C
C  Open RES file
C
      EXT = '.RES'
      FILE_NAME          = ' '
      FILE_NAME(1:NB-1)  = RUN_NAME(1:NB-1)
      FILE_NAME(NB:NB+3) = EXT(1:4)
      IF(RUN_TYPE .EQ. 'NEW')THEN
        OPEN (UNIT=UNIT_RES,FILE=FILE_NAME,STATUS='NEW',RECL=OPEN_N1,
     &      ACCESS='DIRECT',FORM='UNFORMATTED',ERR=300)
      ELSE
        OPEN (UNIT=UNIT_RES,FILE=FILE_NAME,STATUS='OLD',RECL=OPEN_N1,
     &      ACCESS='DIRECT',FORM='UNFORMATTED',ERR=300)
      ENDIF
      IF(NO_FILES .EQ. 0) THEN
        OPEN_FILEP = .TRUE.
        RETURN
      ENDIF
C
C Open SPx files
C
      EXT = '.SPx'
      DO 200 LC = 1,N_SPX
        WRITE (EXT(4:4),1000) LC
        FILE_NAME          = ' '
        FILE_NAME(1:NB-1)  = RUN_NAME(1:NB-1)
        FILE_NAME(NB:NB+3) = EXT(1:4)
        IF(RUN_TYPE .EQ. 'NEW') THEN
          OPEN (UNIT=UNIT_SPX+LC,FILE=FILE_NAME,STATUS='NEW',
     &          RECL=OPEN_N1,ACCESS='DIRECT',FORM='UNFORMATTED',ERR=150)
        ELSE
          OPEN (UNIT=UNIT_SPX+LC,FILE=FILE_NAME,STATUS='OLD',
     &          RECL=OPEN_N1,ACCESS='DIRECT',FORM='UNFORMATTED',ERR=150)
        ENDIF
        SPX_OPEN(LC) = .TRUE.
        GOTO 200
150     WRITE (*,*) 'ERROR OPENING FILE: ', FILE_NAME
        SPX_OPEN(LC) = .FALSE.
200   CONTINUE
      OPEN_FILEP = .TRUE.
      RETURN
300   WRITE (*,*) 'ERROR OPENING FILE: ', FILE_NAME
      RETURN
1000  FORMAT(I1)
      END
CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: OUT_ARRAY (ARRAY, MESSAGE)                             C
C  Purpose: print out a 3D array to standard output                    C
C                                                                      C
C  Author: P.Nicoletti                                Date: 02-DEC-91  C
C  Reviewer: W. Rogers, M. Syamlal, S. Venkatesan     Date: 31-JAN-92  C
C                                                                      C
C  Revision Number:                                                    C
C  Purpose:                                                            C
C  Author:                                            Date: dd-mmm-yy  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Literature/Document References:                                     C
C                                                                      C
C  Variables referenced: KMAX2                                         C
C  Variables modified: None                                            C
C                                                                      C
C  Local variables: POINTER, LK                                        C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
C
      SUBROUTINE OUT_ARRAY (ARRAY,MESSAGE)
C
      IMPLICIT NONE
C
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE 'geometry.inc'
      INCLUDE 'fldvar.inc'
      INCLUDE 'physprop.inc'
      INCLUDE 'indices.inc'
      INCLUDE 'funits.inc'
C
C passed arguments
C
C                       array to print out
      DOUBLE PRECISION  ARRAY(*)
C
C                       message to print out
      CHARACTER*(*)     MESSAGE
C
C local variables
C
C                       pointer into array (points to start of a k-plane)
      INTEGER           IJK
C
C                       loop counter
      INTEGER           L
C
      INCLUDE 'function.inc'
C
      DO 100 L = 1,KMAX2
         IJK = FUNIJK (1,1,L)
         WRITE (UNIT_OUT,1100) MESSAGE , L
         CALL OUT_ARRAY_K (ARRAY(IJK))
100   CONTINUE
C
1100  FORMAT(/,1X,A,' at K = ' ,I4,/)
C
      RETURN
      END
CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: OUT_ARRAY_C (ARRAY,MESSAGE)                            C
C  Purpose: print out a 3D array to standard output (character)        C
C                                                                      C
C  Author: P.Nicoletti                                Date: 10-JAN-92  C
C  Reviewer: W. Rogers, M. Syamlal, S. Venkatesan     Date: 31-JAN-92  C
C                                                                      C
C  Revision Number:                                                    C
C  Purpose:                                                            C
C  Author:                                            Date: dd-mmm-yy  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Literature/Document References:                                     C
C                                                                      C
C  Variables referenced: KMAX2                                         C
C  Variables modified: None                                            C
C                                                                      C
C  Local variables: POINTER, LK                                        C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
C
      SUBROUTINE OUT_ARRAY_C(ARRAY,MESSAGE)
C
      IMPLICIT NONE
C
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE 'geometry.inc'
      INCLUDE 'fldvar.inc'
      INCLUDE 'physprop.inc'
      INCLUDE 'indices.inc'
      INCLUDE 'funits.inc'
C
C passed arguments
C
C                       array to print out
      CHARACTER*3       ARRAY(*)
C
C                       message to print out
      CHARACTER*(*)     MESSAGE
C
C local variables
C
C                       pointer into array (points to start of a k-plane)
      INTEGER           IJK
C
C                       loop counter
      INTEGER           L
C
      INCLUDE 'function.inc'
C
      DO 100 L = 1,KMAX2
         IJK = FUNIJK (1,1,L)
         WRITE (UNIT_OUT,1100) MESSAGE , L
         CALL OUT_ARRAY_KC (ARRAY(IJK))
100   CONTINUE
C
1100  FORMAT(/,1X,A,' at K = ' ,I4,/)
C
      RETURN
      END
CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: OUT_ARRAY_K(ARRAY)                                     C
C  Purpose: print out a 2D (constant k-plane) array to standard output C
C                                                                      C
C  Author: P.Nicoletti                                Date: 02-DEC-91  C
C  Reviewer: W. Rogers, M. Syamlal, S. Venkatesan     Date: 31-JAN-92  C
C                                                                      C
C  Revision Number:                                                    C
C  Purpose:                                                            C
C  Author:                                            Date: dd-mmm-yy  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Literature/Document References:                                     C
C                                                                      C
C  Variables referenced: IMAX2, JMAX2                                  C
C  Variables modified:                                                 C
C                                                                      C
C  Local variables: NCOL, NTAB, L1, L2, L3, IFORM1, IFORM2, IJ1, IJ2   C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE OUT_ARRAY_K (ARRAY)
C
      IMPLICIT NONE
C
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE 'geometry.inc'
      INCLUDE 'fldvar.inc'
      INCLUDE 'physprop.inc'
      INCLUDE 'indices.inc'
      INCLUDE 'funits.inc'
C
C passed arguments
C
C                      2D array to print out
      DOUBLE PRECISION ARRAY(*)
C
C local variables
C
C                      number of columns to print out across the page
      INTEGER          NCOL
C
C                      number of tables the 2D array must be split into
C                      for printing
      INTEGER          NTAB
C
C                      loop indices
      INTEGER          L1x, L2x, L3, IJK
C
C                      start and end 'I' for current table
      INTEGER          IFORM1 , IFORM2
C
C                      start 'IJ' and end 'IJ' for a given 'J' to print out
      INTEGER          IJ1 , IJ2
C
      INCLUDE 'function.inc'
C
C NOTE:  IF NCOL IS CHANGED TO A NUMBER GREATER THAN 30, THEN THE "30"
C        IN FORMATS 5050 AND 5100 MUST BE CHANGED TO THAT NUMBER.
C
      NCOL = 10
      NTAB = IMAX2 / NCOL   +  1
      IF ( MOD(IMAX2,NCOL).EQ.0 ) NTAB = NTAB - 1
C
      DO 100 L1x = 1,NTAB
         IFORM1 = 1 + NCOL*(L1x-1)
         IFORM2 = NCOL * L1x
         IFORM2 = MIN(IFORM2,IMAX2)
         WRITE (UNIT_OUT,5050) (L3,L3=IFORM1,IFORM2)
         DO 50 L2x = JMAX2,1,-1
            IJ1 = FUNIJK(IFORM1,L2x,1)
            IJ2 = FUNIJK(IFORM2,L2x,1)
            WRITE (UNIT_OUT,5100) L2x , (ARRAY(L3),L3=IJ1,IJ2)
50       CONTINUE
100   CONTINUE
C
5050  FORMAT (3X,'J',3X,'I=',3X,10(I3,9X))
5100  FORMAT (1X,I3,3X,10(1PE12.4))
      RETURN
      END
CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: OUT_ARRAY_KC (ARRAY)                                   C
C  Purpose: print out a 2D (constant k-plane) array to standard output C
C           (character)                                                C
C                                                                      C
C  Author: P.Nicoletti                                Date: 02-DEC-91  C
C  Reviewer: W. Rogers, M. Syamlal, S. Venkatesan     Date: 31-JAN-92  C
C                                                                      C
C  Revision Number:                                                    C
C  Purpose:                                                            C
C  Author:                                            Date: dd-mmm-yy  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Literature/Document References:                                     C
C                                                                      C
C  Variables referenced: IMAX2, JMAX2                                  C
C  Variables modified: None                                            C
C                                                                      C
C  Local variables: NCOL, NTAB, L1, L2, L3, IFORM1, IFORM2, IJ1, IJ2   C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE OUT_ARRAY_KC (ARRAY)
C
      IMPLICIT NONE
C
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE 'geometry.inc'
      INCLUDE 'fldvar.inc'
      INCLUDE 'physprop.inc'
      INCLUDE 'indices.inc'
      INCLUDE 'funits.inc'
C
C passed arguments
C
C                      2D array to print out
      CHARACTER*3      ARRAY(*)
C
C local variables
C
C                      number of columns to print out across the page
      INTEGER          NCOL
C
C                      number of tables the 2D array must be split into
C                      for printing
      INTEGER          NTAB
C
C                      loop indices
      INTEGER          L1x, L2x, L3, IJK
C
C                      start and end 'I' for current table
      INTEGER          IFORM1 , IFORM2
C
C                      start 'IJ' and end 'IJ' for a given 'J' to print out
      INTEGER          IJ1 , IJ2
C
      INCLUDE 'function.inc'
C
C NOTE:  IF NCOL IS CHANGED TO A NUMBER GREATER THAN 30, THEN THE "30"
C        IN FORMATS 5050 AND 5100 MUST BE CHANGED TO THAT NUMBER.
C
      NCOL = 30
      NTAB = IMAX2 / NCOL   +  1
      IF ( MOD(IMAX2,NCOL).EQ.0 ) NTAB = NTAB - 1
C
      DO 100 L1x = 1,NTAB
         IFORM1 = 1 + NCOL*(L1x-1)
         IFORM2 = NCOL * L1x
         IFORM2 = MIN(IFORM2,IMAX2)
         WRITE (UNIT_OUT,5050) (L3,L3=IFORM1,IFORM2)
         DO 50 L2x = JMAX2,1,-1
            IJ1 = FUNIJK(IFORM1,L2x,1)
            IJ2 = FUNIJK(IFORM2,L2x,1)
            WRITE (UNIT_OUT,5100) L2x , (ARRAY(L3),L3=IJ1,IJ2)
50       CONTINUE
100   CONTINUE
C
5050  FORMAT (3X,'J',3X,'I=',3X,30(I3,1X))
5100  FORMAT (1X,I3,8X,30(A3,1X))
      RETURN
      END
CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: OUT_BIN_512                                            C
C  Purpose: write an array in chunks of 512 bytes    (DP WORDS)        C
C                                                                      C
C  Author: P. Nicoletti                               Date: 02-JAN-92  C
C  Reviewer: P. Nicoletti, W. Rogers, M. Syamlal      Date: 24-JAN-92  C
C                                                                      C
C  Revision Number:                                                    C
C  Purpose:                                                            C
C  Author:                                            Date: dd-mmm-yy  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Literature/Document References:                                     C
C                                                                      C
C  Variables referenced:                                               C
C  Variables modified:                                                 C
C                                                                      C
C  Local variables: NWORDS, DS, L, NSEG, NREM, LC, N1, N2              C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
C
      SUBROUTINE OUT_BIN_512(IUNIT,ARRAY,N,NEXT_REC)
C
      IMPLICIT NONE
C
      INCLUDE 'machine.inc'
C
C passed arguments
C
C                      array to write out
      DOUBLE PRECISION ARRAY(*)
C
C                      output unit number
      INTEGER          IUNIT
C
C                      number of elements in ARRAY
      INTEGER          N
C
C                      next record number in direct access output file
      INTEGER          NEXT_REC
C
C local variables
C
C                      number of words for 512 bytes
      INTEGER          NWORDS
C
C                      loop counter
      INTEGER          L
C
C                      number of full 512 byte segments need to write N
C                      double precision words
      INTEGER          NSEG
C
C                      number of double precision words in the partially
C                      filled last record
      INTEGER          NREM
C
C                      loop counter
      INTEGER          LC
C
C                      write out array elements N1 to N2
      INTEGER          N1 , N2
C
      NWORDS = NWORDS_DP
      IF (N.LE.NWORDS) THEN
         WRITE (IUNIT,REC=NEXT_REC) (ARRAY(L),L=1,N)
         NEXT_REC = NEXT_REC + 1
         RETURN
      END IF
C
      NSEG = N / NWORDS
      NREM = MOD(N,NWORDS)
      N1 = 1
      N2 = NWORDS
C
C read the full 512 byte segments
C
      DO 100 LC = 1,NSEG
         WRITE (IUNIT,REC=NEXT_REC) (ARRAY(L),L=N1,N2)
         N1 = N1 + NWORDS
         N2 = N2 + NWORDS
         NEXT_REC = NEXT_REC + 1
100   CONTINUE
C
C read the partially filled last record
C
      IF (NREM.NE.0) THEN
         WRITE (IUNIT,REC=NEXT_REC) (ARRAY(L),L=N1,N)
         NEXT_REC = NEXT_REC + 1
      END IF
C
      RETURN
      END
CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: OUT_BIN_512I                                           C
C  Purpose: write out an array in chunks of 512 bytes (INTEGER WORDS)  C
C                                                                      C
C  Author: P. Nicoletti                               Date: 02-JAN-92  C
C  Reviewer: P. Nicoletti, W. Rogers, M. Syamlal      Date: 24-JAN-92  C
C                                                                      C
C  Revision Number:                                                    C
C  Purpose:                                                            C
C  Author:                                            Date: dd-mmm-yy  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Literature/Document References:                                     C
C                                                                      C
C  Variables referenced:                                               C
C  Variables modified:                                                 C
C                                                                      C
C  Local variables: NWORDS, L, NSEG, NREM, LC, N1, N2                  C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
C
      SUBROUTINE OUT_BIN_512I(IUNIT,ARRAY,N,NEXT_REC)
C
      IMPLICIT NONE
C
      INCLUDE 'machine.inc'
C
C passed arguments
C
C                      array to write out
      INTEGER          ARRAY(*)
C
C                      output unit number
      INTEGER          IUNIT
C
C                      number of elements in ARRAY
      INTEGER          N
C
C                      next record number in direct access output file
      INTEGER          NEXT_REC
C
C local variables
C
C                      number of words for 512 bytes (nwords * 4 = 512)
      INTEGER          NWORDS
C
C                      loop counter
      INTEGER          L
C
C                      number of full 512 byte segments need to write N
C                      double precision words
      INTEGER          NSEG
C
C                      number of double precision words in the partially
C                      filled last record
      INTEGER          NREM
C
C                      loop counter
      INTEGER          LC
C
C                      write out array elements N1 to N2
      INTEGER          N1 , N2
C
      NWORDS = NWORDS_I
      IF (N.LE.NWORDS) THEN
         WRITE (IUNIT,REC=NEXT_REC) (ARRAY(L),L=1,N)
         NEXT_REC = NEXT_REC + 1
         RETURN
      END IF
C
      NSEG = N / NWORDS
      NREM = MOD(N,NWORDS)
      N1 = 1
      N2 = NWORDS
C
C write out the full 512 byte segments
C
      DO 100 LC = 1,NSEG
         WRITE (IUNIT,REC=NEXT_REC) (ARRAY(L),L=N1,N2)
         N1 = N1 + NWORDS
         N2 = N2 + NWORDS
         NEXT_REC = NEXT_REC + 1
100   CONTINUE
C
C write out the partially filled last record
C
      IF (NREM.NE.0) THEN
         WRITE (IUNIT,REC=NEXT_REC) (ARRAY(L),L=N1,N)
         NEXT_REC = NEXT_REC + 1
      END IF
C
      RETURN
      END
CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: OUT_BIN_512R                                           C
C  Purpose: write out an array in chunks of 512 bytes (REAL    WORDS)  C
C                                                                      C
C  Author: P. Nicoletti                               Date: 02-JAN-92  C
C  Reviewer: P. Nicoletti, W. Rogers, M. Syamlal      Date: 24-JAN-92  C
C                                                                      C
C  Revision Number:                                                    C
C  Purpose:                                                            C
C  Author:                                            Date: dd-mmm-yy  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Literature/Document References:                                     C
C                                                                      C
C  Variables referenced:                                               C
C  Variables modified:                                                 C
C                                                                      C
C  Local variables: NWORDS, L, NSEG, NREM, LC, N1, N2                  C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
C
      SUBROUTINE OUT_BIN_512R(IUNIT,ARRAY,N,NEXT_REC)
C
      IMPLICIT NONE
C
      INCLUDE 'machine.inc'
C
C passed arguments
C
C                      array to write out
      REAL             ARRAY(*)
C
C                      output unit number
      INTEGER          IUNIT
C
C                      number of elements in ARRAY
      INTEGER          N
C
C                      next record number in direct access output file
      INTEGER          NEXT_REC
C
C local variables
C
C                      number of words for 512 bytes (nwords * 4 = 512)
      INTEGER          NWORDS
C
C                      loop counter
      INTEGER          L
C
C                      number of full 512 byte segments need to write N
C                      double precision words
      INTEGER          NSEG
C
C                      number of double precision words in the partially
C                      filled last record
      INTEGER          NREM
C
C                      loop counter
      INTEGER          LC
C
C                      write out array elements N1 to N2
      INTEGER          N1 , N2
C
      NWORDS = NWORDS_R
      IF (N.LE.NWORDS) THEN
         WRITE (IUNIT,REC=NEXT_REC) (ARRAY(L),L=1,N)
         NEXT_REC = NEXT_REC + 1
         RETURN
      END IF
C
      NSEG = N / NWORDS
      NREM = MOD(N,NWORDS)
      N1 = 1
      N2 = NWORDS
C
C write out the full 512 byte segments
C
      DO 100 LC = 1,NSEG
         WRITE (IUNIT,REC=NEXT_REC) (ARRAY(L),L=N1,N2)
         N1 = N1 + NWORDS
         N2 = N2 + NWORDS
         NEXT_REC = NEXT_REC + 1
100   CONTINUE
C
C write out the partially filled last record
C
      IF (NREM.NE.0) THEN
         WRITE (IUNIT,REC=NEXT_REC) (ARRAY(L),L=N1,N)
         NEXT_REC = NEXT_REC + 1
      END IF
C
      RETURN
      END
CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: OUT_BIN_R                                              C
C  Purpose: write out a time-dependent restart variable (REAL)         C
C                                                                      C
C  Author: P. Nicoletti                               Date: 13-DEC-91  C
C  Reviewer: P. Nicoletti, W. Rogers, M. Syamlal      Date: 24-JAN-92  C
C                                                                      C
C  Revision Number:                                                    C
C  Purpose:                                                            C
C  Author:                                            Date: dd-mmm-yy  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Literature/Document References:                                     C
C                                                                      C
C  Variables referenced:                                               C
C  Variables modified:                                                 C
C                                                                      C
C  Local variables: LC, ARRAY_REAL                                     C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
C
      SUBROUTINE OUT_BIN_R(IUNIT,ARRAY,IJKMAX2,NEXT_REC)
C
      IMPLICIT NONE
      INCLUDE  'param.inc'
      INCLUDE 'param1.inc'
C
C passed arguments
C
C                      double precision array to write out
      DOUBLE PRECISION ARRAY(*)
C
C                      unit number to write to
      INTEGER          IUNIT
C
C                      record pointer into file IUNIT
      INTEGER          NEXT_REC
C
C                      number of indices in ARRAY to write out
      INTEGER          IJKMAX2
C
C local variables
C
C                      single precision version of ARRAY
      REAL             ARRAY_REAL(DIMENSION_3)
C
C                      loop counter
      INTEGER          LC
C
      DO 100 LC = 1,IJKMAX2
         ARRAY_REAL(LC) = SNGL(ARRAY(LC))
100   CONTINUE
C
      CALL OUT_BIN_512R (IUNIT,ARRAY_REAL,IJKMAX2,NEXT_REC)
C
      RETURN
      END
