CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: RES_FROM_SPX(AT_EOF,READ_SPX,REC_POINTER, TIME_REAL)   C
C  Purpose: Create a RES file from the corresponding SPX files         C
C           NOTE : USES N_SPX                                          C
C                                                                      C
C  Author: P. Nicoletti                               Date: 12-MAR-92  C
C  Reviewer:                                                           C
C                                                                      C
C  Revision Number:                                                    C
C  Purpose:                                                            C
C  Author:                                            Date: dd-mmm-yy  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Literature/Document References:                                     C
C                                                                      C
C  Variables referenced: IJKMAX2                                       C
C  Variables modified: IJK                                             C
C                                                                      C
C  Local variables: TIME_FOR_RES, L, PRINTED_MESS, ERROR               C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE RES_FROM_SPX(AT_EOF,READ_SPX,REC_POINTER, TIME_REAL)
C
      IMPLICIT NONE
C
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE 'geometry.inc'
      INCLUDE 'indices.inc'
      INCLUDE 'run.inc'
      INCLUDE 'funits.inc'
      INCLUDE 'post3d.inc'
      INCLUDE 'physprop.inc'
      INCLUDE 'fldvar.inc'
C
      REAL     TIME_FOR_RES, TIME_FOUND
      LOGICAL  AT_EOF(*), READ_SPX(*)
      INTEGER  REC_POINTER(*), REC_POINTER_t(N_SPX)
      INTEGER  NSTEP_1 , ERROR_CODE
      CHARACTER*1 IANS 
      REAL     TIME_REAL(*)
      LOGICAL  PRINTED_MESS(N_SPX)
      LOGICAL  ERROR
C
      INTEGER L, IJK
C
      ERROR = .FALSE.
C
      WRITE (*,'(A,$)') 'Enter time to retrieve from Spx files > '
      READ  (*,*) TIME_FOR_RES
C
      DO 10 L = 1,N_SPX
         READ_SPX(L) = .FALSE.
         REC_POINTER(L) = 4
         AT_EOF(L) = .FALSE.
         PRINTED_MESS(L) = .FALSE.
10    CONTINUE
      DO 20 L = 1,N_SPX
         READ_SPX(L) = .TRUE.
         IF(L .GT. 1)READ_SPX(L-1)=.FALSE.
         CALL SEEK_TIME(READ_SPX, TIME_FOR_RES, REC_POINTER, TIME_FOUND)
         REC_POINTER_t(L) =  REC_POINTER(L)
         IF(TIME_FOUND .NE. TIME_FOR_RES) THEN
           WRITE(*,'(A, I2, A, G12.5)')
     &       ' Did not find time in SPX file ' , L,
     &      '.  Time found = ',TIME_FOUND
           ERROR = .TRUE.
         ELSE
           WRITE (*,'(A, I2)') ' Found time in SPX file ', L
           PRINTED_MESS(L) = .TRUE.
         ENDIF
20    CONTINUE
      DO 30 L = 1,N_SPX
         READ_SPX(L) = .TRUE.
30    CONTINUE
      CALL READ_SPX1(READ_SPX,REC_POINTER_t,AT_EOF, TIME_REAL,NSTEP_1)
C
      ERROR_CODE = 0
      IF (ERROR) THEN
         IF (MMAX.GT.1 .AND. .NOT.PRINTED_MESS(5)) THEN
            WRITE (*,*) ' '
c            WRITE (*,*) ' MMAX > 1 and did not find EP_s'
c            WRITE (*,*) ' can not create restart file'
c            RETURN
         END IF
         IF (.NOT.PRINTED_MESS(5) .AND. .NOT.PRINTED_MESS(1)) THEN
            WRITE (*,*) ' '
c            WRITE (*,*) ' Did not find EP_s or EP_g'
c            WRITE (*,*) ' can not create restart file'
c            RETURN
         END IF
         IF (.NOT.PRINTED_MESS(1) .AND. PRINTED_MESS(5)) THEN
            WRITE (*,*) ' '
            WRITE (*,*) ' Did not find  EP_g   will use EP_s'
            WRITE (*,*) ' '
         END IF
         IF (.NOT.PRINTED_MESS(5) .AND. PRINTED_MESS(1)) THEN
            WRITE (*,*) ' '
            WRITE (*,*) ' Did not find  EP_s   will use EP_g'
            ERROR_CODE = 1
            WRITE (*,*) ' '
         END IF
         WRITE (*,*) ' '
         WRITE (*,*) ' Create RESTART file anyway ? (Y for yes)'
         READ  (*,'(1A1)') IANS
         IF (IANS .NE. 'Y' .AND. IANS .NE. 'y') RETURN
         WRITE (*,'(A,$)') ' Time and DT ? '
         READ(*,*)TIME, DT

      END IF
C
C
C
      WRITE(*,*)
      WRITE(*,'(A)')' The old restart file will be over-written.'
      WRITE(*,'(A,G12.5,A)')' The records in SPx files after time = ',
     & TIME_FOR_RES, ' will be irrecoverably lost!'
      WRITE(*,'(A,$)')' Press Y to over write RES and SPX files '
      READ(*,'(1a1)') IANS
      IF (IANS .NE. 'Y' .AND. IANS .NE. 'y') RETURN
C
C CHANGE THE POINTER RECORD IN THE SPx files to correspond to TIME =
C TIME_FOR_RES.  Times greater than TIME_FOR_RES will be overwritten
C when MFIX is continued.  In SPX files with TIME_FOUND > TIME_FOR_RES
C time is reset as TIME_FOR_RES
C
      DO L = 1,N_SPX
        IF(ERROR .AND. PRINTED_MESS(L) ) THEN
          READ (UNIT_SPX+L,REC=REC_POINTER(L)) TIME_FOUND , NSTEP_1
        ENDIF
        IF(.NOT. PRINTED_MESS(L))THEN
          READ (UNIT_SPX+L,REC=REC_POINTER(L)) TIME_FOUND , NSTEP
          IF(TIME_FOUND .GT. TIME_FOR_RES) THEN
            WRITE (UNIT_SPX+L,REC=REC_POINTER(L))TIME_FOR_RES, NSTEP_1
          ENDIF
        ENDIF
        REC_POINTER(L) =  REC_POINTER(L) + NUM_REC(L)
        WRITE (UNIT_SPX+L,REC=3) REC_POINTER(L), NUM_REC(L)
      END DO
C
C set the time and cycle to the appropriate values in the RESTART file
C
      TIME  = TIME_FOR_RES
      NSTEP = NSTEP_1
C
C MAKE SURE THAT THE VOLUME FRACTIONS ADD UP TO 1.0 ... MIGHT NOT DUE TO
C ROUND OFF ERROR IN GOING FROM DOUBLE TO SINGLE PRECISION
C
      IF (ERROR_CODE.EQ.0) THEN
c         DO IJK = 1,IJKMAX2
c            CALL CALC_EP_g(IJK)
c         END DO
c      ELSE
         if(mmax .gt. 1)Then
           WRITE(*,'(A)')' Cannot update SP5 file. Modify Post_mfix'
           return
         endif
         DO IJK = 1,IJKMAX2
            IF(EP_g(IJK) .NE. UNDEFINED) THEN
              ROP_s(IJK,1) = (1.0 - EP_g(IJK)) * RO_s(1)
            ELSE
              ROP_s(IJK,1) = ZERO
            ENDIF
         END DO
      END IF
C
C WRITE OUT THE RESTART FILE
C
      CALL WRITE_RES1
C
      WRITE(*,*)
      WRITE(*,'(A,A,$)')' RES and SPx files over written.',
     &  ' Press any key to continue.'
      READ(*,'(1a1)') IANS
      RETURN
      END
