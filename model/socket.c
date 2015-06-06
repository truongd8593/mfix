/*
Code taken from the example on Wikipedia:  http://en.wikipedia.org/wiki/Select_(Unix)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
 
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
 
#include <sys/select.h>
#include <fcntl.h>
#include <unistd.h>
#include <err.h>
#include <errno.h>
 
#define PORT "8888"
 
extern void handle_command_(char*, char*, int*);

/* function prototypes */
void die(const char*);
void check_socket();
void init_socket();
 
int sockfd, new, maxfd, on = 1, nready, ii;
struct addrinfo *res0, *res, hints;
struct timeval timeout;
char buffer[BUFSIZ];
char outbuffer[BUFSIZ];
fd_set master, readfds;
int error;
ssize_t nbytes;

/*
// for testing

int main(int argc, char **argv) {
  init_socket();

  while(1) {
      printf("DO SOMETHING USEFUL\n");
      sleep(1);
      check_socket();
  }
  return 0;
}*/

void init_socket() {
 
  (void)memset(&hints, '\0', sizeof(struct addrinfo));
 
  hints.ai_family = AF_INET;
  hints.ai_socktype = SOCK_STREAM;
  hints.ai_protocol = IPPROTO_TCP;
  hints.ai_flags = AI_PASSIVE;
 
  if(0 != (error = getaddrinfo(NULL, PORT, &hints, &res0)))
    errx(EXIT_FAILURE, "%s", gai_strerror(error));
 
  for(res = res0; res; res = res->ai_next)
    {
      if(-1 == (sockfd = socket(res->ai_family, res->ai_socktype, res->ai_protocol)))
	{
	  perror("socket()");
	  continue;
	}
 
      if(-1 == (setsockopt(sockfd, SOL_SOCKET, SO_REUSEADDR, (char*)&on, sizeof(int))))
	{
	  perror("setsockopt()");
	  continue;
	}
 
      if(-1 == (bind(sockfd, res->ai_addr, res->ai_addrlen)))
	{
	  perror("bind");
	  continue;
	}
 
      break;
 
    }
 
  if(-1 == sockfd)
    exit(EXIT_FAILURE);
 
  freeaddrinfo(res0);
 
  if(-1 == (listen(sockfd, 32)))
    die("listen()");
 
  if(-1 == (fcntl(sockfd, F_SETFD, O_NONBLOCK)))
    die("fcntl()");
 
  FD_ZERO(&master);
  FD_ZERO(&readfds);
 
  FD_SET(sockfd, &master);
 
  maxfd = sockfd;
 
  timeout.tv_sec = 0;
  timeout.tv_usec = 0;

}

    void check_socket() {

      memcpy(&readfds, &master, sizeof(master));
 
      (void)printf("running select()\n");
 
      if(-1 == (nready = select(maxfd+1, &readfds, NULL, NULL, &timeout)))
	die("select()");
 
      (void)printf("Number of ready descriptor: %d\n", nready);
 
      for(ii=0; ii<=maxfd && nready>0; ii++)
	{
	  if(FD_ISSET(ii, &readfds))
	    {
	      nready--;
 
	      if(ii == sockfd)
		{
		  (void)printf("Trying to accept() new connection(s)\n");
 
		  if(-1 == (new = accept(sockfd, NULL, NULL)))
		    {
		      if(EWOULDBLOCK != errno)
			perror("accept()");
 
		      break;
		    }
 
		  else
		    {
 
		      if(-1 == (fcntl(new, F_SETFD, O_NONBLOCK)))
			die("fcntl()");
 
		      FD_SET(new, &master);
 
		      if(maxfd < new)
			maxfd = new;
		    }
		}
 
	      else
		{
		  (void)printf("recv() data from one of descriptors(s)\n");
 
		  nbytes = recv(ii, buffer, sizeof(buffer), 0);
		  if(nbytes <= 0)
		    {
		      if(EWOULDBLOCK != errno)
			perror("recv()");
 
		      break;
		    }
 
		  (void)printf("%zi bytes received.\n", nbytes);

		  buffer[nbytes] = '\0';
		  printf("Client said:  %s", buffer);
		  fwrite(buffer,sizeof(char), nbytes, stdout);

		  handle_command_(outbuffer, buffer, &nbytes);

		  if ((nbytes = send(ii, outbuffer, strlen(outbuffer), 0)) == -1) {		  
		    perror("send");
		  }
 
		  (void)printf("%zi bytes sent.\n", nbytes);
 
		  close(ii);
		  FD_CLR(ii, &master);
 
		}
	    }
 
	}
 
}
 
void die(const char *msg)
{
  perror(msg);
  exit(EXIT_FAILURE);
}
