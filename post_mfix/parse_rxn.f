CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: PARSE_RXN(LINE, LMAX)                                  C
C  Purpose: Parse input line                                           C
C                                                                      C
C  Author: P. Nicoletti                               Date: 30-JUN-97  C
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
C  Local variables:                                                    C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
C
      SUBROUTINE PARSE_RXN(LINE, LMAX)
c
      implicit none
      include 'param.inc'
      include 'param1.inc'
c
      include 'parse.inc'
      include 'rxns.inc'
      integer LMAX
c
      character*(*) line
      character*25 field(max_fields)
      character*10 rxn_id
      logical isRate, isDH, isRxn, get_rxn_id
      integer get_rxn_no
C           
      INTEGER nf
C
C           
      INTEGER ier , num_left , species_index
      INTEGER i, is
c
      real value, do_arith

C
c 
      ier = 0
c
      do i = 1,max_fields
	FIELD(i) = ' '
      end do
C
C  Preprocess the string
C
      CALL REMOVE_CHAR('(', LINE, LMAX)         !remove (
      CALL REMOVE_CHAR(')', LINE, LMAX)         !remove )

c
C  Look for a RXN keyword
C
c     Does it start with the rate key word?
      If( IsRate(LINE, LMAX) )THEN

        IF(READING_RATE)CALL CLOSE_READING_RATE
        IF(READING_RXN)CALL CLOSE_READING_RXN

c       yes -- set state to rate -- new rate expresion
        READING_RATE = .TRUE.
        FOUND_m4T    = .FALSE.
        FOUND_PREEXP = .FALSE.
        FOUND_TEXP   = .FALSE.
        FOUND_ACTEMP = .FALSE.

        CALL REMOVE_RATE(LINE, LMAX)

        IF(GET_RXN_ID(LINE, LMAX, RXN_ID))THEN   !rxn identifier found
          RXN_NO = GET_RXN_NO(RXN_ID)
          IF(GOT_RATE(RXN_NO))THEN
            WRITE(*,1010)RXN_ID, LINE(1:LMAX)
            STOP
          ELSE
            GOT_RATE(RXN_NO) = .TRUE.
          ENDIF

        ELSE         !rxn identifier not found -- error               

          WRITE(*,1000)LINE(1:LMAX)
          STOP

        ENDIF


C
c     Does it start with the DH (delta H) keyword ?
      ELSEIF( IsDH(LINE, LMAX) )THEN
        IF(READING_RATE)CALL CLOSE_READING_RATE
        IF(READING_RXN)CALL CLOSE_READING_RXN

        CALL REMOVE_DH(LINE, LMAX)

        IF(GET_RXN_ID(LINE, LMAX, RXN_ID))THEN   !rxn identifier found
          RXN_NO = GET_RXN_NO(RXN_ID)
        ELSE         !rxn identifier not found -- error               

          WRITE(*,1000)LINE(1:LMAX)
          STOP
        ENDIF

        CALL GET_FIELDS(LINE, nf, FIELD , num_left , ier)
c

        is = 1

C       Read enthalpy change due to rxn
        Delta_H(RXN_NO) = DO_ARITH(field(is), len(field(is)), ier)
        if(IER .ne. 0)then
          WRITE(*,1060)LINE(1:LMAX)
          STOP
        endif

        Return

C
c     If none of the above is it the reaction scheme?
      ELSEIF( IsRxn(LINE, LMAX) )THEN

        IF(READING_RATE)CALL CLOSE_READING_RATE
        IF(READING_RXN)CALL CLOSE_READING_RXN

c       rxn identifier found -- set state to rxn -- new rxn
        READING_Rxn = .TRUE.
        FOUND_RHS   = .FALSE.
        backward_rxn = .false.
        IF(GET_RXN_ID(LINE, LMAX, RXN_ID))THEN   !rxn identifier found

          RXN_NO = GET_RXN_NO(RXN_ID)
          IF(GOT_RXN(RXN_NO))THEN
            WRITE(*,1020)RXN_ID, LINE(1:LMAX)
            STOP
          ELSE
            GOT_RXN(RXN_NO) = .TRUE.
          ENDIF
         
        ELSE
c         rxn identifier not found -- error
          WRITE(*,1000)LINE(1:LMAX)
          STOP
        ENDIF
      ELSE
        IF(.NOT.READING_RXN .AND. .NOT.READING_RATE)THEN
          WRITE(*,1050)LINE(1:LMAX)
          STOP
        ENDIF

      ENDIF
C
c  Parse the reaction string and find the stoichiometry
c
C
      IF(READING_RXN)THEN

        CALL REMOVE_CHAR(' ', LINE, LMAX)         !remove spaces
        CALL REMOVE_CHAR('	', LINE, LMAX)  !remove tabs
        CALL REMOVE_CHAR('-', LINE, LMAX)         !remove -
c
c       detect backward rxn specification
        if(index(LINE(1:LMAX), '<') .ne. 0)then
          backward_rxn = .true.
          CALL REMOVE_CHAR('<', LINE, LMAX)         !remove <
        endif

        CALL GET_FIELDS(LINE, nf, FIELD , num_left , ier)
c
 	do i = 1, num_left
	   call get_coef(field(i),species_index,value, ier)

           if(ier .eq. 1)then
             WRITE(*,1051)field(i), LINE(1:LMAX)
             STOP
           elseif(ier .eq. 2)then
             WRITE(*,1052)field(i), LINE(1:LMAX)
             STOP
           endif

           STOICH(RXN_NO, species_index) = 
     &            STOICH(RXN_NO, species_index) - value  !reactant

	end do
	do i = num_left+ 1,nf
	   call get_coef(field(i),species_index,value, ier)

           if(ier .eq. 1)then
             WRITE(*,1051)field(i), LINE(1:LMAX)
             STOP
           elseif(ier .eq. 2)then
             WRITE(*,1052)field(i), LINE(1:LMAX)
             STOP
           endif

           STOICH(RXN_NO, species_index) =  
     &        STOICH(RXN_NO, species_index) + value   !product

	end do

      ENDIF

c
c  Parse the rate string and determine the rate expression
c
      IF(READING_RATE)THEN
        CALL REMOVE_CHAR('[', LINE, LMAX)         !remove [
        CALL REMOVE_CHAR(']', LINE, LMAX)         !remove ]
        CALL REMOVE_CHAR('+', LINE, LMAX)         !remove +
        CALL GET_FIELDS(LINE, nf, FIELD , num_left , ier)
c

        is = 1

C
C       Which temperature to use in the rate expression?
        if(.not.found_m4T .and. is .le. nf)then
          RATE_m4T(RXN_NO) = DO_ARITH(field(is), len(field(is)), ier)
          if(IER .ne. 0)then
            WRITE(*,1060)LINE(1:LMAX)
            STOP
          endif
          found_m4T = .true.
          is = is + 1
        endif

C
C       Preexponential factor
        if(.not.found_preexp .and. is .le. nf)then
          RATE_FAC(RXN_NO, 1) = DO_ARITH(field(is), len(field(is)), ier)
          if(IER .ne. 0)then
            WRITE(*,1062)LINE(1:LMAX)
            STOP
          endif
          found_preexp = .true.
          is = is + 1
        endif

C
C       Tempearture exponent
        if(.not.found_texp .and. is .le. nf)then
          RATE_FAC(RXN_NO, 2) = DO_ARITH(field(is), len(field(is)), ier)
          if(IER .ne. 0)then
            WRITE(*,1064)LINE(1:LMAX)
            STOP
          endif
          found_texp = .true.
          is = is + 1
        endif

C
C       Activation temperature
        if(.not.found_actemp .and. is .le. nf)then
          RATE_FAC(RXN_NO, 3) = DO_ARITH(field(is), len(field(is)), ier)
          if(IER .ne. 0)then
            WRITE(*,1066)LINE(1:LMAX)
            STOP
          endif
          found_actemp = .true.
          is = is + 1
        endif

	do i = is, nf
	   call get_exponent(field(i),species_index,value, ier)

           if(ier .eq. 1)then
             WRITE(*,1051)field(i), LINE(1:LMAX)
             STOP
           elseif(ier .eq. 2)then
             WRITE(*,1052)field(i), LINE(1:LMAX)
             STOP
           endif

           RATE_EXP(RXN_NO, species_index) =  value
	end do
      ENDIF
      Return
1000  FORMAT(/1X,70('*')//' From: PARSE_RXN',
     &/' Message: No reaction id for the rate expression: ',/9X,A
     & ,/1X, 70('*')/)

1010  FORMAT(/1X,70('*')//' From: PARSE_RXN',
     &/' Message: Duplicate rate expression for: ', A,/9X,A
     & ,/1X, 70('*')/)

1020  FORMAT(/1X,70('*')//' From: PARSE_RXN',
     &/' Message: Duplicate reaction scheme for : ', A,/9X,A
     & ,/1X, 70('*')/)

1050  FORMAT(/1X,70('*')//' From: PARSE_RXN',
     &/' Message: No reaction or rate label found: ',/9X,A
     & ,/1X, 70('*')/)

1051  FORMAT(/1X,70('*')//' From: PARSE_RXN',
     &/' Error: Undefined Species: ', A ,
     & /9X, 'in: ',A
     & ,/1X, 70('*')/)

1052  FORMAT(/1X,70('*')//' From: PARSE_RXN',
     &/' Error: Error reading coefficient: ', A ,
     & /9X, 'in: ',A
     & ,/1X, 70('*')/)

1060  FORMAT(/1X,70('*')//' From: PARSE_RXN',
     &/' Message: The first term in a rate expression should be',
     & ' the index for phase temperature. ',/9X,A
     & ,/1X, 70('*')/)

1062  FORMAT(/1X,70('*')//' From: PARSE_RXN',
     &/' Message: The second term in a rate expression should be',
     & ' the preexponential factor. ',/9X,A
     & ,/1X, 70('*')/)

1064  FORMAT(/1X,70('*')//' From: PARSE_RXN',
     &/' Message: The third term in a rate expression should be',
     & ' the temperature exponent. ',/9X,A
     & ,/1X, 70('*')/)

1066  FORMAT(/1X,70('*')//' From: PARSE_RXN',
     &/' Message: The fourth term in a rate expression should be',
     & ' the activation temperature. ',/9X,A
     & ,/1X, 70('*')/)

      END


      LOGICAL FUNCTION IsRate(LINE, LMAX)

      IMPLICIT none

      character*(*) LINE
      integer LMAX
      integer k, kr

      IsRate  = .false.
      k = index(LINE, ':')
      if(k .gt. 0)then
        kr = INDEX(LINE(1:k), 'RATE')
        IF( kr .NE. 0 .AND. (k-kr) .GT. 4)IsRate  = .true.
      endif

      return

      end





      LOGICAL FUNCTION IsDH(LINE, LMAX)

      IMPLICIT none

      character*(*) LINE
      integer LMAX
      integer k, kr

      IsDH  = .false.
      k = index(LINE, ':')
      if(k .gt. 0)then
        kr = INDEX(LINE(1:k), 'DH')
        IF( kr .NE. 0 .AND. (k-kr) .GT. 2)IsDH  = .true.
      endif

      return

      end



      LOGICAL FUNCTION IsRxn(LINE, LMAX)

      IMPLICIT none

      character*(*) LINE
      integer LMAX
      integer k

      IsRxn  = .false.
      k = index(LINE, ':')
      if(k .gt. 1)IsRxn  = .true.

      return

      end

      
      
      subroutine close_reading_rxn
      implicit none
      include 'param.inc'
      include 'param1.inc'
      include 'parse.inc'
      include 'rxns.inc'

      character*10 rxn_id
      integer l, rxn_no_back
      integer GET_RXN_NO

      reading_rxn = .false.
c
c     create a backward rxn stoichiometry, if specified
c
      if(backward_rxn)then
c
c       create a backward reaction name by appending '<'
        rxn_id = rxn_name(rxn_no)
        l = index(rxn_id, ' ')        !append a < or
        if(l .eq. 0)l = len(rxn_id)   !substitute the last char 
        rxn_id(l:l) = '<'
        rxn_no_back = GET_RXN_NO(RXN_ID)
c
c       copy the inverse stoichiometry of forward rxn
        do l = 1, DIMENSION_N_all
          STOICH(rxn_no_back, l) = - STOICH(rxn_no, l)
        enddo
        GOT_RXN(RXN_NO_back) = .TRUE.
        GOT_RATE(RXN_NO_back) = .FALSE.
      endif

      return
      end


      
      subroutine close_reading_rate
      implicit none
      include 'parse.inc'
      reading_rate = .false.
      return
      end

      Subroutine Remove_Rate(LINE, LMAX)
c      Remove the keyword RATE from the line
      IMPLICIT none

      character*(*) LINE
      integer LMAX
      integer k, kr

      k = index(LINE, ':')
      kr = INDEX(LINE(1:k), 'RATE')
      LINE(kr:LMAX-4) = LINE(kr+4:LMAX)
      LINE(LMAX-4:LMAX) = ' '
      LMAX = LMAX - 4

      return

      end

      

      Subroutine Remove_DH(LINE, LMAX)
c      Remove the keyword DH from the line
      IMPLICIT none

      character*(*) LINE
      integer LMAX
      integer k, kr

      k = index(LINE, ':')
      kr = INDEX(LINE(1:k), 'DH')
      LINE(kr:LMAX-2) = LINE(kr+2:LMAX)
      LINE(LMAX-2:LMAX) = ' '
      LMAX = LMAX - 2

      return

      end

      
      LOGICAL FUNCTION get_rxn_id(LINE, LMAX, rxn_id)

      IMPLICIT none

      character*(*) LINE
      character*10 rxn_id
      integer LMAX
      integer k

      get_rxn_id  = .false.
  
      k = index(LINE, ':')
      CALL REMOVE_CHAR(' ', LINE(1:k), k)         !remove spaces

      k = index(LINE, ':')
      CALL REMOVE_CHAR('	', LINE(1:k), k)  !remove tabs

      k = index(LINE, ':')
      if(k .gt. 1)get_rxn_id  = .true.

      rxn_id = LINE(1:k-1)
      LINE(1:LMAX-k) = LINE(k+1:LMAX)
      LINE(LMAX-k+1:LMAX) = ' '

      return

      end

      
      INTEGER FUNCTION get_rxn_no(rxn_id)

      IMPLICIT none
      include 'param.inc'
      include 'param1.inc'
      include 'rxns.inc'

      character*10 rxn_id

      integer l

      do l = 1, no_of_rxns
        if(rxn_id .eq. rxn_name(l))then     !existing reaction
          get_rxn_no = l
          return
        end if
      end do

      no_of_rxns = no_of_rxns  + 1         ! new reaction
      if(no_of_rxns .gt. DIMENSION_RXN) then
        write(*,1010)DIMENSION_RXN
        stop        
      endif
      rxn_name(no_of_rxns) = rxn_id
      got_rxn(no_of_rxns) = .false.
      got_rate(no_of_rxns) = .false.
      get_rxn_no = no_of_rxns

      do l = 1, DIMENSION_N_all
        STOICH(get_rxn_no, l)    = ZERO
        Rate_exp(get_rxn_no, l)  = ZERO
      end do

      return
1010  FORMAT(/1X,70('*')//' From: get_rxn_no',
     &/' Message: Number of reactions must not exceed ',I5,
     &'.  See param.inc.',
     &  /1X, 70('*')/)

      end

      


C
C ****************************************************************
C ************************** GET_FIELDS **************************
C ****************************************************************
C
      SUBROUTINE GET_FIELDS(LINE,nf,FIELD, num_left,ier)
C
      implicit none
      include 'parse.inc'
c
c
      character*(*) line 
      character*(*) field(*)
C
C          
      INTEGER nf
C
C             Return Flags
      INTEGER ier , ier_0
C
C           
      INTEGER I , n , k , ind_gr , num_left


      ier = 0     !initialize error flag
      nf  = 0     !initialize field counter
      I   = 0     !set initial character location in line string
c
      
      ind_gr = index(line,'>')       ! index of ">"
      if(ind_gr .gt. 0)then
        found_rhs = .true.
      elseif(.not.found_rhs)then
        ind_gr = len(LINE)
      endif

C     Read all the fields in the line
c
      num_left = 0
      do k = 1,max_fields
         CALL GET_NEXT_FIELD(LINE,I,N,FIELD(k),ind_gr,num_left,IER_0)
	   if (ier_0.eq.-1) goto 999 !read last field
	   nf = nf + 1
	   i = n
	end do
c
      ier = -1                    ! too many fields
c
c
  999 IF (nf .eq. 0) ier = -1     ! no fields
      RETURN
	END
C
C ****************************************************************
C ************************** GET_NEXT_FIELD **********************
C ****************************************************************
C
      SUBROUTINE GET_NEXT_FIELD(LINE,ISTART,IEND,FIELD,IND_GR,
     *                    NUM_LEFT,IER)
c
      implicit none
c
      character*(*)   line , field
      integer         istart , iend , n , i , i1 , i2 , IND_GR
      integer         i2a , i2b , i2c , i2d, ier , NUM_LEFT
c
      n = len(line)            !number of characters in line string
      do i = istart+1,n
c        find the first character that's not a blank,+,>, and tab
         if (line(i:i).ne.' ' .and. line(i:i).ne.'+' .and.
     *       line(i:i).ne.'	' .and. line(i:i) .ne. '>') then
                i1 = i
                goto 50
         end if
      end do
      ier = -1  !couldn't find another character, error
	return
 50   continue
      
      i2a = index(line(i1:n),' ')    !location of next space
      i2b = index(line(i1:n),'+')    !location of next +
      i2c = index(line(i1:n),'	') !location of next tab
      i2d = index(line(i1:n),'>')    !location of next >
c
      i2 = i2a
      if (i2b.gt.0) i2 = min(i2,i2b)
      if (i2c.gt.0) i2 = min(i2,i2c)
	if (i2d.gt.0) i2 = min(i2,i2d)
c
c     found last field in the line
      if (i2.eq.0) then
         field = line(i1:n)
         iend  = n
	   ier = 0
	   if (i1.lt.ind_gr) num_left = num_left + 1    ! actually ...
         return                                       ! this would be an error
      end if
c
c     specify end of field and field characters
      iend  = i1+i2-2         !set end of field
      field = line(i1:iend)
	if (i1.lt.ind_gr) num_left = num_left + 1  
	ier = 0
c
      return
      end
c
c ********************************************************************
c **************************** get_coef *****************************
c ********************************************************************
c
	subroutine get_coef(field,species_index,value,ier)
c
      implicit none
      include 'param.inc'
      include 'param1.inc'
      include 'rxns.inc'
c
      character*(*) field
      character*25 tmp_field
	integer      species_index
	real         value, do_arith
      integer compare_string

      integer i , k, species_len, field_len, ier
      ier = 0
c
      do i = 1, DIMENSION_N_all
	   k = compare_string(field,species_name(i))
	   if (k.gt.0) then

              species_index = i 
	      tmp_field = field(1:k-1)     ! needed for PC - why ?? Don't know
	      if (tmp_field.eq.' ') then   ! special processing if blank
	         value = 1.0
	         return
	      end if
              value = do_arith(tmp_field, k-1, ier)
              if(ier .ne. 0)ier=2
	      return
            endif
      end do
c
      ier = 1
      species_index = -1    ! error
      return
      end

c
c ********************************************************************
c **************************** get_exponent *********************
c ********************************************************************
c
      subroutine get_exponent(field,species_index,value, ier)
c
      implicit none
      include 'param.inc'
      include 'param1.inc'
      include 'rxns.inc'
c
      character*(*) field
      character*25 tmp_field
	integer      species_index
	real         value, do_arith
      integer compare_string

      integer i , k, ke, species_len, field_len, ier

      ier = 0
c
      ke = index(field,'^')
      if (ke .eq. 0) then   ! special processing if blank
	  value = 1.0
      else
          field_len = len(field)
	  tmp_field = field(ke+1:field_len)
          field(ke:field_len) = ' '
          value = do_arith(tmp_field, field_len-ke, ier)
          if(ier .ne. 0)ier = 2
      end if

      do i = 1, DIMENSION_N_all
	   k = compare_string(field,species_NAME(i))
	   if (k .eq. 1) then
              species_index = i 
	      return
            endif
      end do
c
      ier = 1
      species_index = -1    ! error
      return
      end

c
C *************************************************************
c ************************* compare_string *********************
c *************************************************************
c
      integer function compare_string(line, pattern)
      implicit none
      character*(*) line, pattern
      integer line_len, fill, lpat, lnext, i
      fill = index(pattern, ' ') - 1

      lpat = index(line, pattern(1:fill))
c
c     only a space is allowed after the pattern
      if(lpat .ne. 0) then
        line_len = len(line)
        lnext = lpat+fill
        if(line_len .ge. lnext)then
          if(line(lnext:lnext) .ne. ' ')lpat = 0
        endif
      endif

      compare_string = lpat
      return
      end

c
C *************************************************************
c ************************* remove_char *********************
c *************************************************************
c
      subroutine remove_char(ch, line, nl)
c
      implicit none
c
      character ch
      character*(*) line
      integer   nl , i

	do i = nl-1, 1, -1
	   if (line(i:i) .eq. ch ) then
	      line(i:nl-1) = line(i+1:nl)
	      line(nl:nl) = ' '
	   end if
	end do
c
      return
      end

c
c
c  
      REAL FUNCTION DO_ARITH(LINE, LMAX, ier)

      IMPLICIT NONE

      CHARACTER*(*) LINE
      CHARACTER OPERATION
      CHARACTER*80 SUB_STR
      INTEGER LMAX, ier
      INTEGER L, LSUB
      REAL VALUE, SUB_VALUE

      ier = 0
      DO_ARITH = 986754321.e23
      VALUE = 1.0
      OPERATION = '*'
      LSUB = 1
      DO 100 L = 1, LMAX
        IF(LINE(L:L) .EQ. '*' .OR. LINE(L:L) .EQ. '/')THEN

          IF(LSUB .EQ. 1)THEN
            WRITE(*,1015)LINE(1:LMAX)
            ier = 1
            RETURN
          ENDIF

          READ(SUB_STR(1:LSUB-1),*, ERR=900)SUB_VALUE

          IF(OPERATION .EQ. '*')THEN
            VALUE = VALUE * SUB_VALUE
          ELSEIF(OPERATION .EQ. '/')THEN
            VALUE = VALUE / SUB_VALUE
          ENDIF

          LSUB = 1

          OPERATION = LINE(L:L)

        ELSEIF(LINE(L:L) .EQ. ' ')THEN
        ELSEIF((ICHAR(LINE(L:L)) .GE. 48 .AND.              !0
     &          ICHAR(LINE(L:L)) .LE. 57       ) .OR.       !9
     &          ICHAR(LINE(L:L)) .EQ. 43         .OR.       !+
     &          ICHAR(LINE(L:L)) .EQ. 45         .OR.       !-
     &          ICHAR(LINE(L:L)) .EQ. 46         .OR.       !.
     &          ICHAR(LINE(L:L)) .EQ. 69         .OR.       !E
     &          ICHAR(LINE(L:L)) .EQ. 101            ) THEN !e
          SUB_STR(LSUB:LSUB) = LINE(L:L)
          LSUB = LSUB + 1
        ELSE
          WRITE(*,1015)LINE(1:LMAX)
          ier = 1
          RETURN
        ENDIF
100   CONTINUE

      READ(SUB_STR(1:LSUB-1),*, ERR=900)SUB_VALUE

      IF(OPERATION .EQ. '*')THEN
        VALUE = VALUE * SUB_VALUE
      ELSEIF(OPERATION .EQ. '/')THEN
        VALUE = VALUE / SUB_VALUE
      ENDIF
      
      DO_ARITH = VALUE
      RETURN

900   WRITE(*, 1010)SUB_STR(1:LSUB-1)
      IER = 2
      RETURN

1010  FORMAT(/1X,70('*')//' From: DO_ARITH',
     &/' Message: Error reading the input string: ',/9X,A
     & ,/1X, 70('*')/)

1015  FORMAT(/1X,70('*')//' From: DO_ARITH',
     &/' Message: Invalid operator in the input string: ',/9X,A
     & ,/1X, 70('*')/)
      END

