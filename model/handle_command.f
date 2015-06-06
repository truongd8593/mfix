SUBROUTINE handle_command(outbuffer, buffer, buffer_length)

  USE, INTRINSIC :: ISO_C_BINDING
  USE RUN
  IMPLICIT NONE

  CHARACTER(KIND=C_CHAR, LEN=1024), INTENT(IN)  :: buffer
  INTEGER(C_INT), INTENT(IN) :: buffer_length
  CHARACTER(KIND=C_CHAR, LEN=1024), INTENT(OUT) :: outbuffer

!  INTERACTIVE_MODE=.TRUE.
!  INTERUPT=.TRUE.

  WRITE(outbuffer,"(3(A),F10.5,1(A))")"You said '",buffer(1:buffer_length),"' at TIME = ", TIME, CHAR(0)

END SUBROUTINE handle_command
