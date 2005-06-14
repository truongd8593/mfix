      logical function isNan(x)
      double precision x
!
!                      To check whether x is NAN
      CHARACTER *80 notnumber
      isNan = .False.
      WRITE(notnumber,*) x
! To check for NaN's in x, see if x (a real number) contains a letter "n" or symbol "?"
! in which case it's a NaN (Not a Number)
!
      IF(INDEX(notnumber,'?') > 0 .OR.     &
         INDEX(notnumber,'n') > 0 .OR.     &
         INDEX(notnumber,'N') > 0 ) THEN
        isNan = .TRUE.
	RETURN
      ENDIF
      
      return
      end function isNan
