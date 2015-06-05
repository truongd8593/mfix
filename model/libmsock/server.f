      SUBROUTINE server(sockfd)
  
      USE, INTRINSIC :: ISO_C_BINDING
      USE MSockets
      USE RUN
      IMPLICIT NONE
  
      INTEGER(C_INT), INTENT(IN) :: sockfd
      CHARACTER(KIND=C_CHAR, LEN=1024) :: buffer,outbuffer
      INTEGER(C_INT) :: length, error
  
      IF(sockfd<0) STOP "Failed to open server socket"
  
      WRITE(*,*) "Client connected!"
  
      INTERACTIVE_MODE=.TRUE.
      INTERUPT=.TRUE.
  
      DO
         length=sockGets(sockfd,buffer,INT(LEN_TRIM(buffer),C_SIZE_T))
         IF(length<0) EXIT
         IF(buffer(1:length) == "quit") EXIT
         IF(length>0) THEN
            WRITE(*,*) "Client said: ", buffer(1:length)
            INTERUPT = .NOT.INTERUPT
            WRITE(*,*) "INTERACTIVE_MODE: ",INTERACTIVE_MODE,"  INTERUPT = ",INTERUPT
         ENDIF
 
         WRITE(outbuffer,"(3(A),F10.5)")"You just said ",buffer(1:length)," at TIME = ", TIME
         length=sockPuts(sockfd,outbuffer)
      END DO
      WRITE(*,*) "Client disconnected"
  
      error=CloseSocket(sockfd)
      IF(error<0) STOP "Failed to close server socket"

      END SUBROUTINE server
