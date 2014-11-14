      PROGRAM POST_MFIX
!
      INCLUDE 'xforms.inc'
!
      DO_XFORMS = .FALSE.
!
      CALL F_INIT
!
      STOP
      END
!
      SUBROUTINE  CHECK_INTER(inter)
      ENTRY       ADD_TO_RESULTS_BROWSER(line)
      ENTRY       ADD_TO_SPX_BR(sel,line)
      ENTRY       SPX_TIME_SELECTED(i1,i2)
      ENTRY       SPX_DESELECT_TIME(i1)
      ENTRY       GET_PTX_G
!
      integer          :: inter , i1 , i2
      character(len=*) :: line
      logical          :: sel
!
      RETURN
      END
!
