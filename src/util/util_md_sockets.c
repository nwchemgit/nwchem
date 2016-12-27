/*
 $Id$
 */

#include <sys/types.h>
#include <stdio.h>
#include <unistd.h>
#if !defined(__MINGW32__)
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#endif

#if defined(CRAY) && !defined(__crayx1)
#define create_server_socket_ CREATE_SERVER_SOCKET
#define create_client_socket_ CREATE_CLIENT_SOCKET
#define client_socket_write_ CLIENT_SOCKET_WRITE
#define write_socket_ WRITE_SOCKET
#define close_socket_ CLOSE_SOCKET
#endif


/*
 Routine to create a server socket
 */
int create_server_socket(char *host, int port)
{
#if defined(CATAMOUNT) || defined(__MINGW32__)
    GA_Error("system calls do not work on this machine", 0);
#else
  int server_sockfd, client_sockfd;
  int server_len;
#ifdef NOSOCKLEN
  int client_len;
#else
  socklen_t client_len;
#endif
  struct sockaddr_in server_address;
  struct sockaddr_in client_address;

  server_sockfd= socket(AF_INET,SOCK_STREAM,0);
  server_address.sin_family=AF_INET;
  server_address.sin_addr.s_addr=inet_addr(host);
  server_address.sin_port=htons((ushort) port);
  server_len=sizeof(server_address);
  bind(server_sockfd,(struct sockaddr *)&server_address,server_len);
  listen(server_sockfd,5);
  client_sockfd=accept(server_sockfd,(struct sockaddr *)&client_address,&client_len);
  return client_sockfd;
#endif
}

int create_server_socket_(char *s, int *port)
{
  char host[255]; 
  int i=0;
  while(*s!=' ' && i<255){ host[i]=*s; s++; i++; }; host[i]=0;
  return create_server_socket(host,*port);
}

/*
 Routine to create a client socket to the specified host and port
 Return value is the socket file descriptor when successful, -1 otherwise
 */
int create_client_socket(char *host, int port)
{
#if defined(CATAMOUNT) || defined(__MINGW32__)
    GA_Error("system calls do not work on this machine", 0);
#else
  int sockfd, len, result;
  struct sockaddr_in address;
  /*  int size = 10; */

  sockfd=socket(AF_INET,SOCK_STREAM,0);
  address.sin_family=AF_INET;
  address.sin_addr.s_addr=inet_addr(host);
  address.sin_port=htons((ushort) port);
  len=sizeof(address);
  /*  setsockopt(sockfd,SOL_SOCKET,SO_RCVBUF,(char *) &size,sizeof size);
  setsockopt(sockfd,SOL_SOCKET,SO_SNDBUF,(char *) &size,sizeof size); */
  result=connect(sockfd,(struct sockaddr *)&address,len);

  if(result==-1) return -1;
  return sockfd;
#endif
}

int create_client_socket_(char *s, int *port)
{
  char host[255]; 
  int i=0;

  while(*s!=' ' && i<255){ host[i]=*s; s++; i++; }; host[i]=0;
  return create_client_socket(host,*port);
}


long client_socket_write_(int *sockfd, char *value, int *vlen)
{
  long nb;
  nb=write(*sockfd,value,*vlen);
  return nb;
}

void close_socket_(int *sockfd)
{
  close(*sockfd);
}





