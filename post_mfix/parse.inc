C                      Starting and ending strings for arithmetic operation
      CHARACTER        START_STR*2, END_STR, RXN_BLK*4, END_BLK*3
C
      PARAMETER       (START_STR='@(')
      PARAMETER       (END_STR = ')')
      PARAMETER       (RXN_BLK = 'RXNS')
      PARAMETER       (END_BLK = 'END')

      integer max_fields
      parameter (max_fields = 10)
      integer rxn_no
c
      LOGICAL         READING_RXN, READING_RATE, FOUND_RHS
      LOGICAL         FOUND_M4T, FOUND_PREEXP, FOUND_TEXP, FOUND_ACTEMP
      LOGICAL         backward_rxn
      COMMON/PARSE_L/ READING_RXN, READING_RATE, FOUND_RHS,
     &                FOUND_M4T, FOUND_PREEXP, FOUND_TEXP, FOUND_ACTEMP,
     &                backward_rxn
      COMMON/PARSE_I/ rxn_no
