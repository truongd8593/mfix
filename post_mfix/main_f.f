      PROGRAM POST_MFIX
C
      INCLUDE 'xforms.inc'
C
      DO_XFORMS = .FALSE.
C
      CALL F_INIT
C
      STOP
      END
C
      SUBROUTINE  CHECK_INTER
      ENTRY       ADD_TO_RESULTS_BROWSER
      ENTRY       ADD_TO_SPX_BR
      ENTRY       SPX_TIME_SELECTED
      ENTRY       SPX_DESELECT_TIME
      ENTRY       GET_PTX_G
      RETURN
      END
C
