      PROGRAM MAKEFILE_MFIX
!
!  Author: P. Nicoletti / M. Syamlal
!
! create a make file for MFIX
!
! assumptions :
!
!    1) module ABC is located in file  abc_mod.f
!    2) no non-module source file contains the string :
!               "_mod."    in it's name
!    3) list of all FORTRAN source files located in 
!       the file : files.lis    (one file per line)
!
      CHARACTER   FILENAME*60 , FN(400)*60 , IN(500)*60 , NOEXT*60
      CHARACTER   WoPATH*60   , f90_mods(400)*60 , fname*60
      CHARACTER   BACK*1      , inuse(500)*60
!
      BACK = CHAR(92)
!
      WRITE(*,*)
      WRITE(*,*) 'Determining dependencies...'
      WRITE(*,*) 
!
      OPEN (UNIT=10,FILE='files.lis',STATUS='OLD')            
!
C loop thru  files.lis ...
c
c     modules only     : in f90_mods(i) ... i = 1,num_mods
c     all files        : in fn(i)       ... i = 1,nfiles
C
      NFILES   = 0
      num_mods = 0
      DO WHILE(.TRUE.)
         READ (10,'(A)',END=100) fname
         NFILES = NFILES + 1
	 fn(nfiles) = fname
	 i = index(fname,'_mod.')
	 if (i.ne.0) then
	    num_mods = num_mods + 1
	    f90_mods(num_mods) = fname
	 end if
      END DO
100   CLOSE (UNIT=10)

C
C  Initialization
C
      OPEN (UNIT=30,FILE='mfix_u.make',STATUS='UNKNOWN')
      OPEN (UNIT=31,FILE='mfix_l.make',STATUS='UNKNOWN')
C
C  mfix.exe dependencies ( all object files + blas90.a )
C
      WRITE (30,'(A,A)')  '.$(FORTRAN_EXT).$(OBJ_EXT):'
      WRITE (30,'(A,A)')  '	$(FORTRAN_CMD) $(FORT_FLAGS) $<'
      WRITE (30,*) ' '
      WRITE (30,'(A,A)')  'mfix.exe : ' , BACK
      WRITE (31,'(A,A)')  '.$(FORTRAN_EXT).$(OBJ_EXT):'
      WRITE (31,'(A,A)')  '	$(FORTRAN_CMD) $(FORT_FLAGS) $<'
      WRITE (31,*) ' '
      WRITE (31,'(A,A)')  'mfix.exe : ' , BACK
c
      DO I = 1,NFILES
         CALL NO_EXT(FN(I),NOEXT,N)
         CALL NO_PATH(NOEXT, WoPATH, N, N1)
         WRITE (30,'(A,A,A,A)') '    ' , WoPATH(1:N1),'.$(OBJ_EXT) ',
     &                                     BACK
         WRITE (31,'(A,A,A,A)') '    ' , WoPATH(1:N1),'.$(OBJ_EXT) ',
     &                                     BACK
      END DO
      write (30,'(a)') '    blas90.a '
      write (31,'(a)') '    blas90.a '
C
C  mfix.exe link statement ( all object files + $(LIB_FLAGS) )
C
      WRITE (30,'(A,A)') '	$(LINK_CMD) $(LINK_FLAGS) ',BACK
      WRITE (31,'(A,A)') '	$(LINK_CMD) $(LINK_FLAGS) ',BACK
      DO I = 1,NFILES
         CALL NO_EXT(FN(I),NOEXT,N)
         CALL NO_PATH(NOEXT, WoPATH, N, N1)
         WRITE (30,'(A,A,A,A)') '    ' , WoPATH(1:N1),
     &                               '.$(OBJ_EXT) ',  BACK
         WRITE (31,'(A,A,A,A)') '    ' , WoPATH(1:N1),
     &                               '.$(OBJ_EXT) ',  BACK
      END DO
      WRITE (30,*) ' -o mfix.exe $(LIB_FLAGS)'
      WRITE (30,*) ' '            
      WRITE (31,*) ' -o mfix.exe $(LIB_FLAGS)'
      WRITE (31,*) ' ' 
c
c source code dependencies   ... blas90.a
c
      write (30,'(a)') 'blas90.a : BLAS.o'
      write (30,'(a)') '	ar cr blas90.a BLAS.o'
      write (30,'(a)') 'BLAS.o : BLAS.F'           
      WRITE (30,'(A)')'	$(FORTRAN_CMD) $(FORT_FLAGS) BLAS.F'
      write (31,'(a)') 'blas90.a : BLAS.o'
      write (31,'(a)') '	ar cr blas90.a BLAS.o'
      write (31,'(a)') 'BLAS.o : BLAS.F'           
      WRITE (31,'(A)')'	$(FORTRAN_CMD) $(FORT_FLAGS) BLAS.F'
C
C  source code Dependencies  ... MODULES
C
      DO I = 1,num_mods
         CALL NO_EXT(F90_MODS(I),NOEXT,N)
         CALL NO_PATH(NOEXT, WoPATH, N, N1)
         CALL GET_INC_FILES(NOEXT,N,IN,NINC)
	 call get_use_files(noext,n,inuse,nuse)
         j = index(WoPATH(1:N1),'_mod')
         IF (NINC+NUSE.EQ.0) THEN
            call make_upper_case(WoPATH(1:j-1),j-1)
            WRITE (30,'(A,A,A,A)') 
     &                WoPATH(1:J-1),'.mod : ' , NOEXT(1:N),'.f '
            call make_lower_case(WoPATH(1:j-1),j-1)
            WRITE (31,'(A,A,A,A)') 
     &                WoPATH(1:J-1),'.mod : ' , NOEXT(1:N),'.f '
         ELSE
            call make_upper_case(WoPATH(1:j-1),j-1)
            WRITE (30,'(A,A,A,A,A)') 
     &          WoPATH(1:J-1),'.mod : ' , NOEXT(1:N),'.f ',BACK
            call make_lower_case(WoPATH(1:j-1),j-1)
            WRITE (31,'(A,A,A,A,A)') 
     &          WoPATH(1:J-1),'.mod : ' , NOEXT(1:N),'.f ',BACK
            CALL UNIQUE_FILES(INUSE,NUSE)
            DO II = 1,NUSE
               call get_num_chars(inuse(ii),j)
               if (ii.eq.nuse .and. ninc.eq.0) then
                  call make_upper_case(inuse(ii)(1:j),j)
                  WRITE (30,*) '           ' , INUSE(II)(1:j),
     &                                                 '.mod '
                  call make_lower_case(inuse(ii)(1:j),j)
                  WRITE (31,*) '           ' , INUSE(II)(1:j),
     &                                                 '.mod '
               else
                  call make_upper_case(inuse(ii)(1:j),j)
                  WRITE (30,*) '           ' , INUSE(II)(1:j),
     &                                                 '.mod ',BACK
                  call make_lower_case(inuse(ii)(1:j),j)
                  WRITE (31,*) '           ' , INUSE(II)(1:j),
     &                                                 '.mod ',BACK
               end if
            END DO
            CALL UNIQUE_FILES(IN,NINC)
            DO II = 1,NINC
               if (ii.eq.ninc) then
                  WRITE (30,*) '           ' , IN(II) 
                  WRITE (31,*) '           ' , IN(II) 
               else
                  WRITE (30,*) '           ' , IN(II) , ' ' , BACK
                  WRITE (31,*) '           ' , IN(II) , ' ' , BACK
               end if
            END DO
         END IF

         WRITE (30,'(A,A,A)')
     &       '	$(FORTRAN_CMD) $(FORT_FLAGS) ', NOEXT(1:N),'.f '
         WRITE (31,'(A,A,A)')
     &       '	$(FORTRAN_CMD) $(FORT_FLAGS) ', NOEXT(1:N),'.f '
      END DO
C
C  source code Dependencies (object files)
C
      DO I = 1,NFILES
	 j = index(fn(i),'_mod.')
	 if (j.ne.0) goto 5189
         CALL NO_EXT(FN(I),NOEXT,N)
         CALL NO_PATH(NOEXT, WoPATH, N, N1)
         CALL GET_INC_FILES(NOEXT,N,IN,NINC)
	 call get_use_files(noext,n,inuse,nuse)
         IF (NINC+NUSE.EQ.0) THEN
            WRITE (30,'(A,A,A,A,A)') 
     &                WoPATH(1:N1),'.$(OBJ_EXT) : ' , NOEXT(1:N),'.f '
            WRITE (31,'(A,A,A,A,A)') 
     &                WoPATH(1:N1),'.$(OBJ_EXT) : ' , NOEXT(1:N),'.f '
         ELSE
            WRITE (30,'(A,A,A,A,A)') 
     &          WoPATH(1:N1),'.$(OBJ_EXT) : ' , NOEXT(1:N),'.f ',BACK
            WRITE (31,'(A,A,A,A,A)') 
     &          WoPATH(1:N1),'.$(OBJ_EXT) : ' , NOEXT(1:N),'.f ',BACK
            CALL UNIQUE_FILES(INUSE,NUSE)
            DO II = 1,NUSE
               call get_num_chars(inuse(ii),j)
               if (ii.eq.nuse .and. ninc.eq.0) then
                  call make_upper_case(inuse(ii)(1:j),j)
                  WRITE (30,*) '           ' , INUSE(II)(1:j),
     &                                                 '.mod '
                  call make_lower_case(inuse(ii)(1:j),j)
                  WRITE (31,*) '           ' , INUSE(II)(1:j),
     &                                                 '.mod '
               else
                  call make_upper_case(inuse(ii)(1:j),j)
                  WRITE (30,*) '           ' , INUSE(II)(1:j),
     &                                                 '.mod ',BACK
                  call make_lower_case(inuse(ii)(1:j),j)
                  WRITE (31,*) '           ' , INUSE(II)(1:j),
     &                                                 '.mod ',BACK
               end if
            END DO
            CALL UNIQUE_FILES(IN,NINC)
            DO II = 1,NINC
               if (ii.eq.ninc) then
                  WRITE (30,*) '           ' , IN(II) 
                  WRITE (31,*) '           ' , IN(II) 
               else
                  WRITE (30,*) '           ' , IN(II) , ' ' , BACK
                  WRITE (31,*) '           ' , IN(II) , ' ' , BACK
               end if
            END DO
         END IF

!
! for files in subdirectories compilation rule needs 
!        to be stated explicitly
!
! special processing of files : mark_phase_4.f  and
!                               calc_col_fr.f
!
!     for PGI compilers (needs -O1 optimization flag)
!
         IF(N1 .LT. N) THEN      
           WRITE (30,'(A,A,A)')
     &       '	$(FORTRAN_CMD) $(FORT_FLAGS) ', NOEXT(1:N),'.f '
           WRITE (31,'(A,A,A)')
     &       '	$(FORTRAN_CMD) $(FORT_FLAGS) ', NOEXT(1:N),'.f '
         ENDIF
         IF(fn(i).eq.'mark_phase_4_cor.f' .or.
     &      fn(i).eq.'calc_vol_fr.f') THEN      
           WRITE (30,'(A,A,A)')
     &       '	$(FORTRAN_CMD) $(FORT_FLAGS2) ', NOEXT(1:N),'.f '
           WRITE (31,'(A,A,A)')
     &       '	$(FORTRAN_CMD) $(FORT_FLAGS2) ', NOEXT(1:N),'.f '
         ENDIF
 5189    continue
      END DO
!
      close (unit=30)
      close (unit=31)
      STOP
      END
!
      SUBROUTINE NO_EXT(FN,NOEXT,N)
      CHARACTER  FN*60 , NOEXT*60
!
      NOEXT = ' '
      DO I = 1,60
         IF (FN(I:I).EQ.'.') THEN
             N = I - 1
             NOEXT = FN(1:N)
         END IF
      END DO
!
      RETURN
      END
!
      SUBROUTINE NO_PATH(NOEXT, WoPATH, N, N1)
      CHARACTER  WoPATH*60 , NOEXT*60
!
!  strip the path name from file name
!
      WoPath = ' '
      N1 = 0
      DO I = N, 1, -1
         IF (NOEXT(I:I).EQ.'/') THEN
             N1 = I
             GOTO 10
         END IF
      END DO
10    DO I = 1, N - N1
        WoPATH(I:I) = NOEXT(I+N1:I+N1)
      END DO
      N1 = N - N1
!
      RETURN
      END
!
      SUBROUTINE GET_INC_FILES(NOEXT,N,IN,NINC)
      CHARACTER*60   NOEXT , IN(*) , F1 , LINE*80
!
      F1 = NOEXT(1:N)
      F1(N+1:N+5)='.f   '
      OPEN (UNIT=20,FILE=F1,STATUS='OLD')
      NINC = 1
      DO WHILE(.TRUE.)
         READ (20,'(A)',END=999) LINE
         CALL FIND_INC(LINE,IN(NINC),NINC)
      END DO
999   NINC = NINC - 1
      close (unit=20)
      RETURN
      END
!
      SUBROUTINE FIND_INC(LINE,FNAME,NINC)
      CHARACTER*60 LINE,FNAME , QUOTE*1
!
      IF (LINE(1:1).EQ.'C') RETURN
      IF (LINE(1:1).EQ.'c') RETURN
      IF (LINE(1:1).EQ.'!') RETURN
      QUOTE = CHAR(39)
      DO I = 1,45
         IF (LINE(I:I+6).EQ.'INCLUDE' .OR.
     &       LINE(I:I+6).EQ.'include') then
             DO J = I+7,50
                IF (LINE(J:J).EQ.QUOTE) THEN
                   IQ1 = J
                   GOTO 100
                END IF
             END DO
             write (*,*) line
	     goto 101
100          DO J = IQ1+1,60
                IF (LINE(J:J).EQ.QUOTE) THEN
                   IQ2 = J
                   GOTO 200
                END IF
             END DO
             write (*,*) line
             STOP ' ERROR 2'
101          continue
         END IF
      END DO
      RETURN
200   CONTINUE
      FNAME = LINE(IQ1+1:IQ2-1)
      NINC = NINC + 1
      RETURN
      END
!
      subroutine unique_files(in,ninc)
      character*60 in(*) , out(500)
!
      if (ninc.le.1) return
      out(1) = in(1)
      next = 1
      do 100 i = 2,ninc
         do j = 1,i-1
            if (in(i).eq.in(j)) goto 100
         end do
         next = next + 1
         out(next) = in(i)
100   continue
!
      do i = 1,next
         in(i) = out(i)
      end do
      ninc = next
!
      return
      end
!
!
      SUBROUTINE GET_USE_FILES(NOEXT,N,INUSE,NUSE)
      CHARACTER*60   NOEXT , INUSE(*) , F1 , LINE*80
!
      F1 = NOEXT(1:N)
      F1(N+1:N+5)='.f   '
      OPEN (UNIT=20,FILE=F1,STATUS='OLD')
      NUSE = 1
      DO WHILE(.TRUE.)
         READ (20,'(A)',END=999) LINE
         CALL FIND_USE(LINE,INUSE(NUSE),NUSE)
      END DO
999   NUSE = NUSE - 1
      close (unit=20)
      RETURN
      END
!
      SUBROUTINE FIND_USE(LINE,USENAME,NUSE)
      CHARACTER*60 LINE,USENAME
!
      IF (LINE(1:1).EQ.'C') RETURN
      IF (LINE(1:1).EQ.'c') RETURN
      IF (LINE(1:1).EQ.'!') RETURN
      usename = ' '
      DO I = 1,20
         IF (LINE(I:I+3).EQ.'USE ' .OR.
     &       LINE(I:I+3).EQ.'Use ') then
             iq1 = i + 4
100          DO J = IQ1+1,60
                IF (LINE(J:J).EQ.' ' .or.LINE(J:J).EQ.',') THEN
                   IQ2 = J
                   GOTO 200
                END IF
             END DO
         END IF
      END DO
      RETURN
200   CONTINUE
      USENAME = LINE(IQ1:IQ2-1)
      NUSE = NUSE + 1
      RETURN
      END
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: MAKE_UPPER_CASE (LINE_STRING,MAXCOL)                   C
!  Purpose: change lowercase characters to uppercase in input line     C
!                                                                      C
!  Author: P.Nicoletti                                Date: 26-NOV-91  C
!  Reviewer: M.SYAMLAL, W.ROGERS, P.NICOLETTI         Date: 24-JAN-92  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: None                                          C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables: A_UP, A_LO, Z_LO, A_DIFF, INT_C, L                 C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE MAKE_UPPER_CASE(LINE_STRING, MAXCOL) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                   input line to change to uppercase
      CHARACTER*(*) LINE_STRING
!
!                   number of characters to look at in LINE_STRING
      INTEGER       MAXCOL
!
! local variables:
!
!                   ICHAR value for UPPERCASE A
      INTEGER       A_UP
!
!                   ICHAR value for lowercase a
      INTEGER       A_LO
!
!                   ICHAR value for lowercase z
      INTEGER       Z_LO
!
!                   ICHAR differnce between lower and uppercase letters
      INTEGER       A_DIFF
!
!                   holds ICHAR value of current character
      INTEGER       INT_C
!
!                   loop index
      INTEGER       L
!-----------------------------------------------
!
!
      A_UP = ICHAR('A') 
      A_LO = ICHAR('a') 
      Z_LO = ICHAR('z') 
      A_DIFF = A_LO - A_UP 
!
      DO L = 1, MAXCOL 
         INT_C = ICHAR(LINE_STRING(L:L)) 
         IF (A_LO.le.INT_C .AND. INT_C.le.Z_LO) THEN 
            INT_C = INT_C - A_DIFF 
            LINE_STRING(L:L) = CHAR(INT_C) 
         ENDIF 
      END DO 
      RETURN  
      END 
!
      subroutine get_num_chars(string,n)
      implicit none
      character*(*) string
      integer       i,n
!
      do i = 1,len(string)
         if (string(i:i).eq.' ') then
            n = i-1
            return
         end if
      end do
      n = len(string)
      return
      end
!
      SUBROUTINE MAKE_LOWER_CASE(LINE_STRING, MAXCOL) 
!
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                   input line to change to lowercase
      CHARACTER*(*) LINE_STRING
!
!                   number of characters to look at in LINE_STRING
      INTEGER       MAXCOL
!
! local variables:
!
!                   ICHAR value for UPPERCASE A
      INTEGER       A_UP
!
!                   ICHAR value for lowercase a
      INTEGER       A_LO
!
!                   ICHAR value for UPPERCASE Z
      INTEGER       Z_UP
!
!                   ICHAR differnce between lower and uppercase letters
      INTEGER       A_DIFF
!
!                   holds ICHAR value of current character
      INTEGER       INT_C
!
!                   loop index
      INTEGER       L
!-----------------------------------------------
!
!
      A_UP = ICHAR('A') 
      A_LO = ICHAR('a') 
      Z_UP = ICHAR('Z') 
      A_DIFF = A_LO - A_UP 
!
      DO L = 1, MAXCOL 
         INT_C = ICHAR(LINE_STRING(L:L)) 
         IF (A_UP.le.INT_C .AND. INT_C.le.Z_UP) THEN 
            INT_C = INT_C + A_DIFF 
            LINE_STRING(L:L) = CHAR(INT_C) 
         ENDIF 
      END DO 
      RETURN  
      END 
