      SUBROUTINE handle_command(outbuffer, buffer, buffer_length) BIND ( C )

      USE, INTRINSIC :: ISO_C_BINDING
      USE RUN
      IMPLICIT NONE

      CHARACTER(KIND=C_CHAR, LEN=1), INTENT(IN)  :: buffer(1024)
      INTEGER(C_INT), INTENT(IN) :: buffer_length
      CHARACTER(KIND=C_CHAR, LEN=1) :: outbuffer(1024)

      CHARACTER(len=1024) :: lMSG

      INTEGER :: LC, MLEN

!  INTERACTIVE_MODE=.TRUE.
!  INTERUPT=.TRUE.

! Copy over the formatted GUI message
      MLEN = 0
      DO LC = 1, 1024
         lMSG(LC:LC) = BUFFER(LC)
         MLEN = MLEN + 1
         IF(BUFFER(LC) == CHAR(00)) EXIT
      ENDDO
! Null terminate the string.


      WRITE(*,"(3x,'at Time ',f9.6,' BUFFER: ',I4)") TIME, MLEN
      WRITE(*,"(3x,'You said: ')",ADVANCE="NO")
      DO LC=1,1024
        WRITE(*,"(A1)",ADVANCE="NO") lMSG(LC:LC)
      ENDDO
      WRITE(*,"('@')")




      print *,"buffer: ",buffer

      stop
!      print *,"You said '",buffer(1:buffer_length),"' at TIME = ", TIME

!      WRITE(outbuffer,"(1(A))") "handled command."

!       WRITE(outchar,"(3(A),F10.5,1(A))")"You said '",buffer(1:buffer_length),"' at TIME = ", TIME, CHAR(0)

      END SUBROUTINE handle_command
