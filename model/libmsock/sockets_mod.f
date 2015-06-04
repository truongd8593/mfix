MODULE MSockets
   ! Interface to libmsock
   ! A library for TCP/IP client-server applications on Unix by Muhammad A Muquit
   ! http://www.muquit.com/muquit/software/libmsock/libmsock.html
   USE, INTRINSIC :: ISO_C_BINDING
   IMPLICIT NONE
   PUBLIC
   
!   INTEGER, PARAMETER :: C_SIZE_T=C_INT
   
   INTERFACE
      ! Interfaces for the main functions
      
      ! int close(int fildes)
      ! This one is not in libmsock technically but is needed to close sockets:
      FUNCTION CloseSocket(fildes) BIND(C,NAME="close")  RESULT(error)
         IMPORT
         INTEGER(C_INT), VALUE :: fildes
         INTEGER(C_INT) :: error
      END FUNCTION          
      ! int ServerSocket(u_short port,int max_servers);
      FUNCTION ServerSocket(port,max_servers) BIND(C,NAME="ServerSocket")  RESULT(sockfd)
         IMPORT
         INTEGER(C_SHORT), VALUE :: port
         INTEGER(C_INT), VALUE :: max_servers
         INTEGER(C_INT) :: sockfd
      END FUNCTION
      ! int ClientSocket(char *netaddress,u_short port);
      FUNCTION ClientSocket(netaddress,port) BIND(C,NAME="ClientSocket")  RESULT(sockfd)
         IMPORT
         CHARACTER(C_CHAR), DIMENSION(*), INTENT(IN) :: netaddress
         INTEGER(C_SHORT), VALUE :: port
         INTEGER(C_INT) :: sockfd
      END FUNCTION
      ! int sockGets(int sockfd,char *str,size_t count);
      FUNCTION sockGets(sockfd,str,count) BIND(C,NAME="sockGets")  RESULT(length)
         IMPORT
         INTEGER(C_INT), VALUE :: sockfd
         CHARACTER(C_CHAR), DIMENSION(*), INTENT(IN) :: str
         INTEGER(C_SIZE_T), VALUE :: count
         INTEGER(C_INT) :: length ! Bytes read
      END FUNCTION
      ! int sockRead(int sockfd,char *str,size_t count);
      FUNCTION sockRead(sockfd,str,count) BIND(C,NAME="sockRead")  RESULT(error)
         IMPORT
         INTEGER(C_INT), VALUE :: sockfd
         INTEGER(C_SIZE_T), VALUE :: count
         CHARACTER(C_CHAR), DIMENSION(count), INTENT(OUT) :: str
         INTEGER(C_INT) :: error
      END FUNCTION
      ! int sockPuts(int sockfd,char *str);
      FUNCTION sockPuts(sockfd,str) BIND(C,NAME="sockPuts")  RESULT(error)
         IMPORT
         INTEGER(C_INT), VALUE :: sockfd
         CHARACTER(C_CHAR), DIMENSION(*), INTENT(IN) :: str ! NULL terminated
         INTEGER(C_INT) :: error
      END FUNCTION
      ! int sockWrite(int sockfd,char *str,size_t count);
      FUNCTION sockWrite(sockfd,str,count) BIND(C,NAME="sockWrite")  RESULT(error)
         IMPORT
         INTEGER(C_INT), VALUE :: sockfd
         INTEGER(C_SIZE_T), VALUE :: count
         CHARACTER(C_CHAR), DIMENSION(count), INTENT(IN) :: str
         INTEGER(C_INT) :: error
      END FUNCTION
   END INTERFACE
      
END MODULE

