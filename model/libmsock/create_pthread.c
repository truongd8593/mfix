/*
 * C function to call pthread_create
 */

#include <pthread.h>
#include <stdio.h>
#include <msock.h>

extern void server_();
extern double simulation_time;

void *thread_function(void *pp);

void create_server() {

  pthread_t thread;
  long return_status;

  return_status = pthread_create(&thread, NULL, thread_function, NULL);

  printf("pthread_create exit status %d",return_status);

}

/*
**  htget - echo servers. prints out whatever a client writes.
**
**  Limitations and Comments:
**  an example of msock library.
**
**  Development History:
**      who                  when           why
**      ma_muquit@fccc.edu   Jul-10-1999    first cut
*/

void *thread_function(void *pp) 
{
    int
        connected,
        rc,
        port,
        client_port,
        nn,
        sock_fd;
    
    char
        szclient_host[64],
        szclient_ip[64],
        szbuf[BUFSIZ],
        szhost[64],
        szpage[1024];

    int argc = 2;
    char *argv[2];
    argv[1] = "8088";

    if (argc < 2)
    {
        (void) fprintf(stderr,"usage: %s port\n",argv[0]);
        return NULL;
    }
    port=atoi(argv[1]);

    /* open socket and start listening */
    sock_fd=ServerSocket((u_short) port,5);
    if (sock_fd == -1)
    {
        (void) fprintf(stderr,"Failed to open socket\n");
        return NULL;
    }

    /* we will be here only when a client connects */

    /* find out the name/ip of the client */
    rc=getPeerInfo(sock_fd,szclient_host,szclient_ip,&client_port);
    if (rc == 0)
    {
        (void) fprintf(stderr,"================================\n");
        (void) fprintf(stderr,"Client hostname %s\n",szclient_host);
        (void) fprintf(stderr,"Client IP: %s\n",szclient_ip);
        (void) fprintf(stderr,"Client port: %d\n",client_port);
        (void) fprintf(stderr,"================================\n");
    }

    server_(&sock_fd);
    return NULL;

    connected=1;
    /* server the client */
    while (connected)
    {
        rc=sockGets(sock_fd,szbuf,sizeof(szbuf));
        if (rc < 0)
            break;
        if (rc)
        {
            if (strcmp(szbuf,"quit") == 0)
            {
                (void) fprintf(stderr,"closing connection with %s ...",
                               szclient_host);
                break;
            }
            (void) fprintf(stdout,"%s\n",szbuf);
            (void) fflush(stdout);

	    //nn = sprintf(szbuf,"Current time is: %g",simulation_time);

	    //rc=sockPuts(sock_fd,szbuf); // rc should equal nn

        }
    }
    

    (void) fprintf(stderr," <closed>\n\n");
    close(sock_fd);
    return(0);
}
