!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: PARSE_RXN(LINE, LMAX)                                  C
!  Purpose: Parse input line                                           C
!                                                                      C
!  Author: P. Nicoletti                               Date: 30-JUN-97  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE PARSE_RXN(LINE, LMAX) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parse 
      USE rxns
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER LMAX 
      CHARACTER LINE*(*) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NF, IER, NUM_LEFT, SPECIES_INDEX, I, IS 
      REAL :: VALUE 
      CHARACTER, DIMENSION(MAX_FIELDS) :: FIELD*25 
      CHARACTER :: RXN_ID*10 
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      INTEGER , EXTERNAL :: GET_RXN_NO 
      REAL , EXTERNAL :: DO_ARITH 
      LOGICAL , EXTERNAL :: ISRATE, ISDH, ISRXN, GET_RXN_ID 
!-----------------------------------------------
!
!
!
!
!
      IER = 0 
!
      DO I = 1, MAX_FIELDS 
         FIELD(I) = ' ' 
      END DO 
!
!  Preprocess the string
!
      CALL REMOVE_CHAR ('(', LINE, LMAX)         !remove ( 
      CALL REMOVE_CHAR (')', LINE, LMAX)         !remove ) 
!
!
!  Look for a RXN keyword
!
!     Does it start with the rate key word?
      IF (ISRATE(LINE,LMAX)) THEN 
!
         IF (READING_RATE) CALL CLOSE_READING_RATE 
         IF (READING_RXN) CALL CLOSE_READING_RXN 
!
!       yes -- set state to rate -- new rate expresion
         READING_RATE = .TRUE. 
         FOUND_M4T = .FALSE. 
         FOUND_PREEXP = .FALSE. 
         FOUND_TEXP = .FALSE. 
         FOUND_ACTEMP = .FALSE. 
!
         CALL REMOVE_RATE (LINE, LMAX) 
!
         IF (GET_RXN_ID(LINE,LMAX,RXN_ID)) THEN  !rxn identifier found 
            RXN_NO = GET_RXN_NO(RXN_ID) 
            IF (GOT_RATE(RXN_NO)) THEN 
               WRITE (*, 1010) RXN_ID, LINE(1:LMAX) 
               STOP  
            ELSE 
               GOT_RATE(RXN_NO) = .TRUE. 
            ENDIF 
!
         ELSE                                    !rxn identifier not found -- error 
!
            WRITE (*, 1000) LINE(1:LMAX) 
            STOP  
!
         ENDIF 
!
!
!
!     Does it start with the DH (delta H) keyword ?
      ELSE IF (ISDH(LINE,LMAX)) THEN 
         IF (READING_RATE) CALL CLOSE_READING_RATE 
         IF (READING_RXN) CALL CLOSE_READING_RXN 
!
         CALL REMOVE_DH (LINE, LMAX) 
!
         IF (GET_RXN_ID(LINE,LMAX,RXN_ID)) THEN  !rxn identifier found 
            RXN_NO = GET_RXN_NO(RXN_ID) 
         ELSE                                    !rxn identifier not found -- error 
!
            WRITE (*, 1000) LINE(1:LMAX) 
            STOP  
         ENDIF 
!
         CALL GET_FIELDS (LINE, NF, FIELD, NUM_LEFT, IER) 
!
!
         IS = 1 
!
!       Read enthalpy change due to rxn
         DELTA_H(RXN_NO) = DO_ARITH(FIELD(IS),LEN(FIELD(IS)),IER) 
         IF (IER /= 0) THEN 
            WRITE (*, 1060) LINE(1:LMAX) 
            STOP  
         ENDIF 
!
         RETURN  
!
!
!     If none of the above is it the reaction scheme?
      ELSE IF (ISRXN(LINE,LMAX)) THEN 
!
         IF (READING_RATE) CALL CLOSE_READING_RATE 
         IF (READING_RXN) CALL CLOSE_READING_RXN 
!
!       rxn identifier found -- set state to rxn -- new rxn
         READING_RXN = .TRUE. 
         FOUND_RHS = .FALSE. 
         BACKWARD_RXN = .FALSE. 
         IF (GET_RXN_ID(LINE,LMAX,RXN_ID)) THEN  !rxn identifier found 
!
            RXN_NO = GET_RXN_NO(RXN_ID) 
            IF (GOT_RXN(RXN_NO)) THEN 
               WRITE (*, 1020) RXN_ID, LINE(1:LMAX) 
               STOP  
            ELSE 
               GOT_RXN(RXN_NO) = .TRUE. 
            ENDIF 
!
         ELSE 
!         rxn identifier not found -- error
            WRITE (*, 1000) LINE(1:LMAX) 
            STOP  
         ENDIF 
      ELSE 
         IF ( .NOT.READING_RXN .AND.  .NOT.READING_RATE) THEN 
            WRITE (*, 1050) LINE(1:LMAX) 
            STOP  
         ENDIF 
!
      ENDIF 
!
!  Parse the reaction string and find the stoichiometry
!
!
      IF (READING_RXN) THEN 
!
         CALL REMOVE_CHAR (' ', LINE, LMAX)      !remove spaces 
         CALL REMOVE_CHAR ('	', LINE, LMAX)      !remove tabs 
         CALL REMOVE_CHAR ('-', LINE, LMAX)      !remove - 
!
!       detect backward rxn specification
         IF (INDEX(LINE(1:LMAX),'<') /= 0) THEN 
            BACKWARD_RXN = .TRUE. 
            CALL REMOVE_CHAR ('<', LINE, LMAX)   !remove < 
         ENDIF 
!
         CALL GET_FIELDS (LINE, NF, FIELD, NUM_LEFT, IER) 
!
         DO I = 1, NUM_LEFT 
            CALL GET_COEF (FIELD(I), SPECIES_INDEX, VALUE, IER) 
!
            IF (IER == 1) THEN 
               WRITE (*, 1051) FIELD(I), LINE(1:LMAX) 
               STOP  
            ELSE IF (IER == 2) THEN 
               WRITE (*, 1052) FIELD(I), LINE(1:LMAX) 
               STOP  
            ENDIF 
!
!                                                !reactant
            STOICH(RXN_NO,SPECIES_INDEX) = STOICH(RXN_NO,SPECIES_INDEX) - VALUE 
!
         END DO 
         DO I = NUM_LEFT + 1, NF 
            CALL GET_COEF (FIELD(I), SPECIES_INDEX, VALUE, IER) 
!
            IF (IER == 1) THEN 
               WRITE (*, 1051) FIELD(I), LINE(1:LMAX) 
               STOP  
            ELSE IF (IER == 2) THEN 
               WRITE (*, 1052) FIELD(I), LINE(1:LMAX) 
               STOP  
            ENDIF 
!
!                                                !product
            STOICH(RXN_NO,SPECIES_INDEX) = STOICH(RXN_NO,SPECIES_INDEX) + VALUE 
!
         END DO 
!
      ENDIF 
!
!
!  Parse the rate string and determine the rate expression
!
      IF (READING_RATE) THEN 
         CALL REMOVE_CHAR ('[', LINE, LMAX)      !remove [ 
         CALL REMOVE_CHAR (']', LINE, LMAX)      !remove ] 
         CALL REMOVE_CHAR ('+', LINE, LMAX)      !remove + 
         CALL GET_FIELDS (LINE, NF, FIELD, NUM_LEFT, IER) 
!
!
         IS = 1 
!
!
!       Which temperature to use in the rate expression?
         IF ( .NOT.FOUND_M4T .AND. IS<=NF) THEN 
            RATE_M4T(RXN_NO) = DO_ARITH(FIELD(IS),LEN(FIELD(IS)),IER) 
            IF (IER /= 0) THEN 
               WRITE (*, 1060) LINE(1:LMAX) 
               STOP  
            ENDIF 
            FOUND_M4T = .TRUE. 
            IS = IS + 1 
         ENDIF 
!
!
!       Preexponential factor
         IF ( .NOT.FOUND_PREEXP .AND. IS<=NF) THEN 
            RATE_FAC(RXN_NO,1) = DO_ARITH(FIELD(IS),LEN(FIELD(IS)),IER) 
            IF (IER /= 0) THEN 
               WRITE (*, 1062) LINE(1:LMAX) 
               STOP  
            ENDIF 
            FOUND_PREEXP = .TRUE. 
            IS = IS + 1 
         ENDIF 
!
!
!       Tempearture exponent
         IF ( .NOT.FOUND_TEXP .AND. IS<=NF) THEN 
            RATE_FAC(RXN_NO,2) = DO_ARITH(FIELD(IS),LEN(FIELD(IS)),IER) 
            IF (IER /= 0) THEN 
               WRITE (*, 1064) LINE(1:LMAX) 
               STOP  
            ENDIF 
            FOUND_TEXP = .TRUE. 
            IS = IS + 1 
         ENDIF 
!
!
!       Activation temperature
         IF ( .NOT.FOUND_ACTEMP .AND. IS<=NF) THEN 
            RATE_FAC(RXN_NO,3) = DO_ARITH(FIELD(IS),LEN(FIELD(IS)),IER) 
            IF (IER /= 0) THEN 
               WRITE (*, 1066) LINE(1:LMAX) 
               STOP  
            ENDIF 
            FOUND_ACTEMP = .TRUE. 
            IS = IS + 1 
         ENDIF 
!
         DO I = IS, NF 
            CALL GET_EXPONENT (FIELD(I), SPECIES_INDEX, VALUE, IER) 
!
            IF (IER == 1) THEN 
               WRITE (*, 1051) FIELD(I), LINE(1:LMAX) 
               STOP  
            ELSE IF (IER == 2) THEN 
               WRITE (*, 1052) FIELD(I), LINE(1:LMAX) 
               STOP  
            ENDIF 
!
            RATE_EXP(RXN_NO,SPECIES_INDEX) = VALUE 
         END DO 
      ENDIF 
      RETURN  
 1000 FORMAT(/1X,70('*')//' From: PARSE_RXN',/&
         ' Message: No reaction id for the rate expression: ',/9X,A,/1X,70('*')&
         /) 
!
 1010 FORMAT(/1X,70('*')//' From: PARSE_RXN',/&
         ' Message: Duplicate rate expression for: ',A,/9X,A,/1X,70('*')/) 
!
 1020 FORMAT(/1X,70('*')//' From: PARSE_RXN',/&
         ' Message: Duplicate reaction scheme for : ',A,/9X,A,/1X,70('*')/) 
!
 1050 FORMAT(/1X,70('*')//' From: PARSE_RXN',/&
         ' Message: No reaction or rate label found: ',/9X,A,/1X,70('*')/) 
!
 1051 FORMAT(/1X,70('*')//' From: PARSE_RXN',/' Error: Undefined Species: ',A,/&
         9X,'in: ',A,/1X,70('*')/) 
!
 1052 FORMAT(/1X,70('*')//' From: PARSE_RXN',/&
         ' Error: Error reading coefficient: ',A,/9X,'in: ',A,/1X,70('*')/) 
!
 1060 FORMAT(/1X,70('*')//' From: PARSE_RXN',/&
         ' Message: The first term in a rate expression should be',&
         ' the index for phase temperature. ',/9X,A,/1X,70('*')/) 
!
 1062 FORMAT(/1X,70('*')//' From: PARSE_RXN',/&
         ' Message: The second term in a rate expression should be',&
         ' the preexponential factor. ',/9X,A,/1X,70('*')/) 
!
 1064 FORMAT(/1X,70('*')//' From: PARSE_RXN',/&
         ' Message: The third term in a rate expression should be',&
         ' the temperature exponent. ',/9X,A,/1X,70('*')/) 
!
 1066 FORMAT(/1X,70('*')//' From: PARSE_RXN',/&
         ' Message: The fourth term in a rate expression should be',&
         ' the activation temperature. ',/9X,A,/1X,70('*')/) 
!
      END SUBROUTINE PARSE_RXN 
!
!
      LOGICAL FUNCTION ISRATE (LINE, LMAX) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parse 
      USE rxns 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER LMAX 
      CHARACTER LINE*(*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K, KR 
!-----------------------------------------------
!
!
      ISRATE = .FALSE. 
      K = INDEX(LINE,':') 
      IF (K > 0) THEN 
         KR = INDEX(LINE(1:K),'RATE') 
         IF (KR/=0 .AND. K-KR>4) ISRATE = .TRUE. 
      ENDIF 
!
      RETURN  
!
      END FUNCTION ISRATE 
!
!
!
!
!
      LOGICAL FUNCTION ISDH (LINE, LMAX) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parse 
      USE rxns 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER LMAX 
      CHARACTER LINE*(*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K, KR 
!-----------------------------------------------
!
!
      ISDH = .FALSE. 
      K = INDEX(LINE,':') 
      IF (K > 0) THEN 
         KR = INDEX(LINE(1:K),'DH') 
         IF (KR/=0 .AND. K-KR>2) ISDH = .TRUE. 
      ENDIF 
!
      RETURN  
!
      END FUNCTION ISDH 
!
!
!
      LOGICAL FUNCTION ISRXN (LINE, LMAX) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parse 
      USE rxns 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER LMAX 
      CHARACTER LINE*(*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K 
!-----------------------------------------------
!
!
      ISRXN = .FALSE. 
      K = INDEX(LINE,':') 
      IF (K > 1) ISRXN = .TRUE. 
!
      RETURN  
!
      END FUNCTION ISRXN 
!
!
!
      SUBROUTINE CLOSE_READING_RXN 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parse 
      USE rxns
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
      INTEGER :: L, RXN_NO_BACK 
      CHARACTER :: RXN_ID*10 
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      INTEGER , EXTERNAL :: GET_RXN_NO 
!-----------------------------------------------
!
!
      READING_RXN = .FALSE. 
!
!     create a backward rxn stoichiometry, if specified
!
      IF (BACKWARD_RXN) THEN 
!
!       create a backward reaction name by appending '<'
         RXN_ID = RXN_NAME(RXN_NO) 
         L = INDEX(RXN_ID,' ')                   !append a < or 
         IF (L == 0) L = LEN(RXN_ID)             !substitute the last char 
         RXN_ID(L:L) = '<' 
         RXN_NO_BACK = GET_RXN_NO(RXN_ID) 
!
!       copy the inverse stoichiometry of forward rxn
         STOICH(RXN_NO_BACK,:DIM_N_ALL)=-STOICH(RXN_NO,:DIM_N_ALL) 
         GOT_RXN(RXN_NO_BACK) = .TRUE. 
         GOT_RATE(RXN_NO_BACK) = .FALSE. 
      ENDIF 
!
      RETURN  
      END SUBROUTINE CLOSE_READING_RXN 
!
!
!
      SUBROUTINE CLOSE_READING_RATE 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parse 
      USE rxns 
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!-----------------------------------------------
      READING_RATE = .FALSE. 
      RETURN  
      END SUBROUTINE CLOSE_READING_RATE 
!
      SUBROUTINE REMOVE_RATE(LINE, LMAX) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parse 
      USE rxns 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER LMAX 
      CHARACTER LINE*(*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K, KR 
!-----------------------------------------------
!
!
      K = INDEX(LINE,':') 
      KR = INDEX(LINE(1:K),'RATE') 
      LINE(KR:LMAX-4) = LINE(KR+4:LMAX) 
      LINE(LMAX-4:LMAX) = ' ' 
      LMAX = LMAX - 4 
!
      RETURN  
!
      END SUBROUTINE REMOVE_RATE 
!
!
!
      SUBROUTINE REMOVE_DH(LINE, LMAX) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parse 
      USE rxns 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER LMAX 
      CHARACTER LINE*(*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K, KR 
!-----------------------------------------------
!
!
      K = INDEX(LINE,':') 
      KR = INDEX(LINE(1:K),'DH') 
      LINE(KR:LMAX-2) = LINE(KR+2:LMAX) 
      LINE(LMAX-2:LMAX) = ' ' 
      LMAX = LMAX - 2 
!
      RETURN  
!
      END SUBROUTINE REMOVE_DH 
!
!
      LOGICAL FUNCTION GET_RXN_ID (LINE, LMAX, RXN_ID) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parse 
      USE rxns 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER LMAX 
      CHARACTER LINE*(*), RXN_ID*10 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K 
!-----------------------------------------------
!
!
      GET_RXN_ID = .FALSE. 
!
      K = INDEX(LINE,':') 
      CALL REMOVE_CHAR (' ', LINE(1:K), K)       !remove spaces 
!
      K = INDEX(LINE,':') 
      CALL REMOVE_CHAR ('	', LINE(1:K), K)       !remove tabs 
!
      K = INDEX(LINE,':') 
      IF (K > 1) GET_RXN_ID = .TRUE. 
!
      RXN_ID = LINE(1:K-1) 
      LINE(1:LMAX-K) = LINE(K+1:LMAX) 
      LINE(LMAX-K+1:LMAX) = ' ' 
!
      RETURN  
!
      END FUNCTION GET_RXN_ID 
!
!
      INTEGER FUNCTION GET_RXN_NO (RXN_ID) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parse 
      USE rxns
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER RXN_ID*10 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: L 
!-----------------------------------------------
!
!
!
      DO L = 1, NO_OF_RXNS 
         IF (RXN_ID == RXN_NAME(L)) THEN         !existing reaction 
            GET_RXN_NO = L 
            RETURN  
         ENDIF 
      END DO 
!
      NO_OF_RXNS = NO_OF_RXNS + 1                ! new reaction 
      IF (NO_OF_RXNS > DIMENSION_RXN) THEN 
         WRITE (*, 1010) DIMENSION_RXN 
         STOP  
      ENDIF 
      RXN_NAME(NO_OF_RXNS) = RXN_ID 
      GOT_RXN(NO_OF_RXNS) = .FALSE. 
      GOT_RATE(NO_OF_RXNS) = .FALSE. 
      GET_RXN_NO = NO_OF_RXNS 
!
      STOICH(GET_RXN_NO,:DIM_N_ALL) = ZERO 
      RATE_EXP(GET_RXN_NO,:DIM_N_ALL) = ZERO 
!
      RETURN  
 1010 FORMAT(/1X,70('*')//' From: get_rxn_no',/&
         ' Message: Number of reactions must not exceed ',I5,&
         '.  See param.inc.',/1X,70('*')/) 
!
      END FUNCTION GET_RXN_NO 
!
!
!
!
!
! ****************************************************************
! ************************** GET_FIELDS **************************
! ****************************************************************
!
      SUBROUTINE GET_FIELDS(LINE, NF, FIELD, NUM_LEFT, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parse 
      USE rxns 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER NF, NUM_LEFT, IER 
      CHARACTER LINE*(*) 
      CHARACTER, DIMENSION(*) :: FIELD*(*) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IER_0, I, N, K, IND_GR 
!-----------------------------------------------
!
!
!
!
!
!             Return Flags
!
!
!
!
      IER = 0                                    !initialize error flag 
      NF = 0                                     !initialize field counter 
      I = 0                                      !set initial character location in line string 
!
!
      IND_GR = INDEX(LINE,'>')                   ! index of ">" 
      IF (IND_GR > 0) THEN 
         FOUND_RHS = .TRUE. 
      ELSE IF (.NOT.FOUND_RHS) THEN 
         IND_GR = LEN(LINE) 
      ENDIF 
!
!     Read all the fields in the line
!
      NUM_LEFT = 0 
      DO K = 1, MAX_FIELDS 
         CALL GET_NEXT_FIELD (LINE, I, N, FIELD(K), IND_GR, NUM_LEFT, IER_0) 
         IF (IER_0 == (-1)) GO TO 999            !read last field 
         NF = NF + 1 
         I = N 
      END DO 
!
      IER = -1                                   ! too many fields 
!
!
  999 CONTINUE 
      IF (NF == 0) IER = -1 
      RETURN                                     ! no fields 
      END SUBROUTINE GET_FIELDS 
!
! ****************************************************************
! ************************** GET_NEXT_FIELD **********************
! ****************************************************************
!
      subroutine get_next_field(line,istart,iend,field,ind_gr,num_left,ier) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      use param 
      use param1 
      use parse 
      use rxns 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer istart, iend, ind_gr, num_left, ier 
      character line*(*), field*(*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: n, i, i1, i2, i2a, i2b, i2c, i2d 
!-----------------------------------------------
!
!
      n = len(line)                              !number of characters in line string 
      do i = istart + 1, n 
!        find the first character that's not a blank,+,>, and tab
         if (line(i:i)/=' ' .and. line(i:i)/='+' .and. line(i:i)/='	' .and. &
            line(i:i)/='>') then 
            i1 = i 
            go to 50 
         endif 
      end do 
      ier = -1                                   !couldn't find another character, error 
      return  
   50 continue 
      i2a = index(line(i1:n),' ')                !location of next space 
      i2b = index(line(i1:n),'+')                !location of next + 
      i2c = index(line(i1:n),'	')                !location of next tab 
      i2d = index(line(i1:n),'>')                !location of next > 
!
      i2 = i2a 
      if (i2b > 0) i2 = min(i2,i2b) 
      if (i2c > 0) i2 = min(i2,i2c) 
      if (i2d > 0) i2 = min(i2,i2d) 
!
!     found last field in the line
      if (i2 == 0) then 
         field = line(i1:n) 
         iend = n 
         ier = 0 
         if (i1 < ind_gr) num_left = num_left + 1! actually ... 
         return                                  ! this would be an error 
      endif 
!
!     specify end of field and field characters
      iend = i1 + i2 - 2                         !set end of field 
      field = line(i1:iend) 
      if (i1 < ind_gr) num_left = num_left + 1 
      ier = 0 
!
      return  
      end subroutine get_next_field 
!
! ********************************************************************
! **************************** get_coef *****************************
! ********************************************************************
!
      SUBROUTINE GET_COEF(FIELD, SPECIES_INDEX, VALUE, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parse 
      USE rxns
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER SPECIES_INDEX, IER 
      REAL VALUE 
      CHARACTER FIELD*(*) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, K, SPECIES_LEN, FIELD_LEN 
      CHARACTER :: TMP_FIELD*25 
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      INTEGER , EXTERNAL :: COMPARE_STRING 
      REAL , EXTERNAL :: DO_ARITH 
!-----------------------------------------------
!
!
      IER = 0 
!
      DO I = 1, DIM_N_ALL 
         K = COMPARE_STRING(FIELD,SPECIES_NAME(I)) 
         IF (K > 0) THEN 
!
	    IER = 0
            SPECIES_INDEX = I 
            TMP_FIELD = FIELD(1:K-1)             ! needed for PC - why ?? Don't know 
            IF (TMP_FIELD == ' ') THEN           ! special processing if blank 
               VALUE = 1.0 
               RETURN  
            ENDIF 
            VALUE = DO_ARITH(TMP_FIELD,K - 1,IER) 
!            IF (IER /= 0) IER = 2 
            IF (IER == 0)  RETURN  
         ENDIF 
      END DO 
!
      IER = 1 
      SPECIES_INDEX = -1                         ! error 
      RETURN  
      END SUBROUTINE GET_COEF 
!
!
! ********************************************************************
! **************************** get_exponent *********************
! ********************************************************************
!
      SUBROUTINE GET_EXPONENT(FIELD, SPECIES_INDEX, VALUE, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parse 
      USE rxns
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER SPECIES_INDEX, IER 
      REAL VALUE 
      CHARACTER FIELD*(*) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, K, KE, SPECIES_LEN, FIELD_LEN 
      CHARACTER :: TMP_FIELD*25 
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      INTEGER , EXTERNAL :: COMPARE_STRING 
      REAL , EXTERNAL :: DO_ARITH 
!-----------------------------------------------
!
!
!
      IER = 0 
!
      KE = INDEX(FIELD,'^') 
      IF (KE == 0) THEN                          ! special processing if blank 
         VALUE = 1.0 
      ELSE 
         FIELD_LEN = LEN(FIELD) 
         TMP_FIELD = FIELD(KE+1:FIELD_LEN) 
         FIELD(KE:FIELD_LEN) = ' ' 
         VALUE = DO_ARITH(TMP_FIELD,FIELD_LEN - KE,IER) 
         IF (IER /= 0) IER = 2 
      ENDIF 
!
      DO I = 1, DIM_N_ALL 
         K = COMPARE_STRING(FIELD,SPECIES_NAME(I)) 
         IF (K == 1) THEN 
            SPECIES_INDEX = I 
            RETURN  
         ENDIF 
      END DO 
!
      IER = 1 
      SPECIES_INDEX = -1                         ! error 
      RETURN  
      END SUBROUTINE GET_EXPONENT 
!
!
! *************************************************************
! ************************* compare_string *********************
! *************************************************************
!
      integer function compare_string (line, pattern) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      use param 
      use param1 
      use parse 
      use rxns 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character line*(*), pattern*(*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: line_len, fill, lpat, lnext, i 
!-----------------------------------------------
      fill = index(pattern,' ') - 1 
!
      lpat = index(line,pattern(1:fill)) 
!
!     only a space is allowed after the pattern
      if (lpat /= 0) then 
         line_len = len(line) 
         lnext = lpat + fill 
         if (line_len >= lnext) then 
            if (line(lnext:lnext) /= ' ') lpat = 0 
         endif 
      endif 
!
      compare_string = lpat 
      return  
      end function compare_string 
!
!
! *************************************************************
! ************************* remove_char *********************
! *************************************************************
!
      subroutine remove_char(ch, line, nl) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      use param 
      use param1 
      use parse 
      use rxns 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer nl 
      character ch, line*(*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i 
!-----------------------------------------------
!
!
      do i = nl - 1, 1, -1 
         if (line(i:i) == ch) then 
            line(i:nl-1) = line(i+1:nl) 
            line(nl:nl) = ' ' 
         endif 
      end do 
!
      return  
      end subroutine remove_char 
!
!
!
!
      REAL FUNCTION DO_ARITH (LINE, LMAX, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parse 
      USE rxns 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER LMAX, IER 
      CHARACTER LINE*(*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: L, LSUB 
      REAL :: VALUE, SUB_VALUE 
      CHARACTER :: OPERATION, SUB_STR*80 
!-----------------------------------------------
!
!
      IER = 0 
      DO_ARITH = 986754321.E23 
      VALUE = 1.0 
      OPERATION = '*' 
      LSUB = 1 
      DO L = 1, LMAX 
         IF (LINE(L:L)=='*' .OR. LINE(L:L)=='/') THEN 
!
            IF (LSUB == 1) THEN 
               WRITE (*, 1015) LINE(1:LMAX) 
               IER = 1 
               RETURN  
            ENDIF 
!
            READ (SUB_STR(1:LSUB-1), *, ERR=900) SUB_VALUE 
!
            IF (OPERATION == '*') THEN 
               VALUE = VALUE*SUB_VALUE 
            ELSE IF (OPERATION == '/') THEN 
               VALUE = VALUE/SUB_VALUE 
            ENDIF 
!
            LSUB = 1 
!
            OPERATION = LINE(L:L) 
!
         ELSE IF (LINE(L:L) == ' ') THEN 
!                                                !e
         ELSE IF (ICHAR(LINE(L:L))>=48 .AND. ICHAR(LINE(L:L))<=57 .OR. ICHAR(&
               LINE(L:L))==43 .OR. ICHAR(LINE(L:L))==45 .OR. ICHAR(LINE(L:L))==&
               46 .OR. ICHAR(LINE(L:L))==69 .OR. ICHAR(LINE(L:L))==101) THEN 
            SUB_STR(LSUB:LSUB) = LINE(L:L) 
            LSUB = LSUB + 1 
         ELSE 
!            WRITE (*, 1015) LINE(1:LMAX) 
            IER = 1 
            RETURN  
         ENDIF 
      END DO 
      READ (SUB_STR(1:LSUB-1), *, ERR=900) SUB_VALUE 
!
      IF (OPERATION == '*') THEN 
         VALUE = VALUE*SUB_VALUE 
      ELSE IF (OPERATION == '/') THEN 
         VALUE = VALUE/SUB_VALUE 
      ENDIF 
!
      DO_ARITH = VALUE 
      RETURN  
!
  900 CONTINUE 
      WRITE (*, 1010) SUB_STR(1:LSUB-1) 
      IER = 2 
      RETURN  
!
 1010 FORMAT(/1X,70('*')//' From: DO_ARITH',/&
         ' Message: Error reading the input string: ',/9X,A,/1X,70('*')/) 
!
 1015 FORMAT(/1X,70('*')//' From: DO_ARITH',/&
         ' Message: Invalid operator in the input string: ',/9X,A,/1X,70('*')/) 
      END FUNCTION DO_ARITH 
