      MODULE parse


      Use param
      Use param1


!                      Starting and ending strings for arithmetic operation
      CHARACTER        START_STR*2, END_STR, RXN_BLK*4, END_BLK*3
!
      PARAMETER       (START_STR='@(')
      PARAMETER       (END_STR = ')')
      PARAMETER       (RXN_BLK = 'RXNS')
      PARAMETER       (END_BLK = 'END')
 
      integer max_fields
      parameter (max_fields = 10)
      integer rxn_no
!
      LOGICAL         READING_RXN, READING_RATE, FOUND_RHS
      LOGICAL         FOUND_M4T, FOUND_PREEXP, FOUND_TEXP, FOUND_ACTEMP
      LOGICAL         backward_rxn


      END MODULE parse                                                                           
