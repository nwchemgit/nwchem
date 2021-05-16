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



#if defined(CRAY)
#include <fortran.h>
#if !defined(__crayx1)
#define USE_FCD
#endif
#endif

#if (defined(CRAY) || defined(WIN32)) &&!defined(__crayx1) &&!defined(__MINGW32__)
#define nwpw_talker_ nwpw_talker
#define nwpw_talker_close_ nwpw_talker_close
#define nwpw_talker_write_ nwpw_talker_write
#define nwpw_talker_read_ nwpw_talker_read
#endif

void FATR nwpw_talker_
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
        perror("nwpw_talker: not coded for this architecture");
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
        
        printf("nwpw_talker: addr_name=%s\n",addr_name);
        printf("nwpw_talker: port=%d\n",port);

        if ((sock = socket(AF_INET, SOCK_STREAM, 0)) < 0) 
        { 
            printf("\nnwpw_talker:Socket creation error \n"); 
            exit(1);
        } 
        serv_addr.sin_family = AF_INET; 
        serv_addr.sin_port = htons(port); 

        // Convert IPv4 and IPv6 addresses from text to binary form 
        if(inet_pton(AF_INET, addr_name, &serv_addr.sin_addr)<=0)  
        { 
            printf("\nnwpw_talker:Invalid address/ Address not supported \n"); 
            exit(1);
        } 
        if (connect(sock, (struct sockaddr *)&serv_addr, sizeof(serv_addr)) < 0) 
        { 
            printf("\nnwpw_talker:Connection Failed \n"); 
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
            printf("\nnwpw_talker:Socket creation error \n");
            exit(1);
        }
        if (connect(sock, (struct sockaddr *)&serv_addr, sizeof(serv_addr)) < 0)
        {
            printf("\nnwpw_talker:Failed to connect to UNIX socket \n");
            exit(1);
        }
    }

    printf("nwpw_talker: sockid=%d\n",sock);
    *sockout = ((Integer) sock);
#endif
}




void FATR nwpw_talker_close_
#if defined(USE_FCD)
(Integer *socket1)
{
#else
(socket1)
Integer *socket1;
{
#endif
#if defined(__MINGW32__)
        perror("nwpw_talker: not coded for this architecture");
        exit(1);
#else
   int sock = ((int) *socket1);
   close(sock);
#endif
}


void FATR nwpw_talker_write_
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
        perror("nwpw_talker_write: not coded for this architecture");
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
   //printf("nwpw_talker_write: sock=%d and n=%d nbuf=%d\n",sock,n,nbuf);
   //printf("nwpw_talker_write: sock=%d and n=%d buf=%s\n",sock,nbuf,buffer);

#endif
}





void FATR nwpw_talker_read_
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
        perror("nwpw_talker_read: not coded for this architecture");
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
   
    //printf("nwpw_talker_read: received %d bytes with nbuf=%d from sock:%d\n", valread,nbuf,sock);
    //printf("nwpw_talker_read: message received: %s\n\n",buffer);

#endif
}






#if (defined(CRAY) || defined(WIN32)) &&!defined(__crayx1) &&!defined(__MINGW32__)
#define nwpw_listener_ nwpw_listener
#endif

void FATR nwpw_listener_
#if defined(USE_FCD)
( const _fcd fcd_addr_name,
 Integer *n1,
 Integer *portin,
 const _fcd fcd_mesg,
 Integer *n3)
{
    char *addr_name = _fcdtocp(fcd_addr_name);
    char *port_name = _fcdtocp(fcd_port_name);
    char *buffer    = _fcdtocp(fcd_mesg);

#else
(addr_name,n1,portin,buffer,n3)
char    addr_name[];
Integer *n1;
Integer *portin;
char    buffer[];
Integer *n3;
{

#endif

#if defined(__MINGW32__)
        perror("nwpw_listener: not coded for this architecture");
        exit(1);
#else
    int sock = 0, valread;
    struct sockaddr_in serv_addr;
    int port = ((int) *portin);
    int nbuf = ((int) *n3);

    if ((sock = socket(AF_INET, SOCK_STREAM, 0)) < 0)
    {
        printf("\nnwpw_listener: Socket creation error \n");
        exit(1);
    }
    serv_addr.sin_family = AF_INET;
    serv_addr.sin_port = htons(port);

    // Convert IPv4 and IPv6 addresses from text to binary form 
    if(inet_pton(AF_INET, addr_name, &serv_addr.sin_addr)<=0)
    {
        printf("\nnwpw_listener:Invalid address/ Address not supported \n");
        exit(1);
    }
    if (connect(sock, (struct sockaddr *)&serv_addr, sizeof(serv_addr)) < 0)
    {
        printf("\nnwpw_listener:Connection Failed \n");
        exit(1);
    }
    valread = read(sock,buffer,nbuf); 

    printf("nwpw_listener: received %d bytes to %s:%d\n", nbuf, addr_name,port);
    printf("nwpw_listener: message received: %s\n\n",buffer);
    close(sock);
#endif

}



