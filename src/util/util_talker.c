/*
 $Id$
*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <sys/types.h>
#if !defined(__MINGW32__)
#include <sys/socket.h>
#include <sys/un.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <netdb.h>
#endif
#include "typesf2c.h"



#if defined(CRAY) || defined(CRAY_T3D)
#include <fortran.h>
#if !defined(__crayx1)
#define USE_FCD
#endif
#endif

#if (defined(CRAY)  || defined(CRAY_T3D) || defined(WIN32)) &&!defined(__crayx1) &&!defined(__MINGW32__)
#define util_talker_ util_talker
#define util_talker_close_ util_talker_close
#define util_talker_write_ util_talker_write
#define util_talker_read_ util_talker_read
#endif

void FATR util_talker_
#if defined(USE_FCD)
( const _fcd fcd_addr_name,
 Integer *inet,
 Integer *n1,
 Integer *portin,
 Integer *sockout)
{
    char *addr_name = _fcdtocp(fcd_addr_name);

#else
(addr_name,inet,n1,portin,sockout)
char	addr_name[];
Integer *inet;
Integer	*n1;
Integer *portin;
Integer *sockout;
{

#endif

#if defined(__MINGW32__)
        perror("util_talker: not coded for this architecture");
        exit(1);
#else
    int sock = 0, valread;
    int na   = ((int) *n1);

    addr_name[na]   = 0;
    addr_name[na+1] = 0;

    if (*inet>0)
    {
        struct sockaddr_in serv_addr; 
        int port = ((int) *portin);
        
        printf("util_talker: addr_name=%s\n",addr_name);
        printf("util_talker: port=%d\n",port);

        if ((sock = socket(AF_INET, SOCK_STREAM, 0)) < 0) 
        { 
            printf("\nutil_talker:Socket creation error \n"); 
            exit(1);
        } 
        serv_addr.sin_family = AF_INET; 
        serv_addr.sin_port = htons(port); 

        // Convert IPv4 and IPv6 addresses from text to binary form 
        if(inet_pton(AF_INET, addr_name, &serv_addr.sin_addr)<=0)  
        { 
            printf("\nutil_talker:Invalid address/ Address not supported \n"); 
            exit(1);
        } 
        if (connect(sock, (struct sockaddr *)&serv_addr, sizeof(serv_addr)) < 0) 
        { 
            printf("\nutil_talker:Connection Failed \n"); 
            exit(1);
        } 
    }
    else
    {
        struct sockaddr_un serv_addr;
        memset(&serv_addr, 0, sizeof(struct sockaddr_un));
        serv_addr.sun_family = AF_UNIX;
        strcpy(serv_addr.sun_path, "/tmp/ipi_");
        strcpy(serv_addr.sun_path+9, addr_name);
        if ((sock = socket(AF_UNIX, SOCK_STREAM, 0)) < 0)
        {
            printf("\nutil_talker:Socket creation error \n");
            exit(1);
        }
        if (connect(sock, (struct sockaddr *)&serv_addr, sizeof(serv_addr)) < 0)
        {
            printf("\nutil_talker:Failed to connect to UNIX socket \n");
            exit(1);
        }
    }

    printf("util_talker: sockid=%d\n",sock);
    *sockout = ((Integer) sock);
#endif
}




void FATR util_talker_close_
#if defined(USE_FCD)
(Integer *socket1)
{
#else
(socket1)
Integer *socket1;
{
#endif
#if defined(__MINGW32__)
        perror("util_talker: not coded for this architecture");
        exit(1);
#else
   int sock = ((int) *socket1);
   close(sock);
#endif
}


void FATR util_talker_write_
#if defined(USE_FCD)
(Integer *socket1,
 const _fcd fcd_buffer,
 Integer *n1)
{
   char *buffer = _fcdtocp(fcd_buffer);
#else
(socket1,buffer,n1)
Integer *socket1;
char    buffer[];
Integer *n1;
{
#endif
#if defined(__MINGW32__)
        perror("util_talker_write: not coded for this architecture");
        exit(1);
#else
   int sock = ((int) *socket1);
   int nbuf = ((int) *n1);
   int it=0;
   int n = 0;
   while (((n += write(sock, &(buffer[n]), nbuf-n)) < nbuf) && (it<100)) {
      buffer[n] = 0;
      ++it;
   }
   //send(sock,buffer,nbuf,0);
   //printf("util_talker_write: sock=%d and n=%d nbuf=%d\n",sock,n,nbuf);
   //printf("util_talker_write: sock=%d and n=%d buf=%s\n",sock,nbuf,buffer);

#endif
}





void FATR util_talker_read_
#if defined(USE_FCD)
(Integer *socket1,
 const _fcd fcd_buffer,
 Integer *n1)
{   
   char *buffer = _fcdtocp(fcd_buffer);
#else
(socket1,buffer,n1)
Integer *socket1;
char    buffer[];
Integer *n1;
{       
#endif
#if defined(__MINGW32__)
        perror("util_talker_read: not coded for this architecture");
        exit(1);
#else   
   int sock = ((int) *socket1);
   int nbuf = ((int) *n1);
   int valread=0;
   int readit=0;
   bzero(buffer,nbuf+2);
   while ((((valread += read(sock,&(buffer[valread]),nbuf-valread))) < nbuf) && (readit<100))
   {
      buffer[valread] = 0;
      ++readit;
   }
   
    //printf("util_talker_read: received %d bytes with nbuf=%d from sock:%d\n", valread,nbuf,sock);
    //printf("util_talker_read: message received: %s\n\n",buffer);

#endif
}







