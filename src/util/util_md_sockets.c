/*
 $Id: util_md_sockets.c,v 1.1 1999-08-10 14:27:39 d3j191 Exp $
 */

#include <sys/types.h>
#include <sys/socket.h>
#include <stdio.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <unistd.h>

int client_socket_open_(int port)
{
  int sockfd;
  int len;
  struct sockaddr_in address;
  int result;

  sockfd=socket(AF_INET,SOCK_STREAM,0);
  address.sin_family=AF_INET;
  address.sin_addr.s_addr=inet_addr("127.0.0.1");
  address.sin_port=3333;
  len=sizeof(address);

  result=connect(sockfd,(struct sockaddr *)&address,len);

  if(result==-1)
    { perror("Unable to connect to NWChem socket"); return sockfd; }

  return sockfd;
}

void client_socket_write_(int sockfd, char *value)
{
  write(sockfd,&value,1);
}

void client_socket_close_(int sockfd)
{
  close(sockfd);
}
