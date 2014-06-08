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
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <netdb.h>
#include "typesf2c.h"



#if defined(CRAY) || defined(CRAY_T3D)
#include <fortran.h>
#if !defined(__crayx1)
#define USE_FCD
#endif
#endif

#if (defined(CRAY) &&!defined(__crayx1)) || defined(CRAY_T3D) || defined(WIN32)
#define nwpw_talker_ nwpw_talker
#endif

void FATR nwpw_talker_
#if defined(USE_FCD)
( const _fcd fcd_addr_name,
 Integer *n1,
 const _fcd fcd_port_name,
 Integer *n2,
 const _fcd fcd_mesg,
 Integer *n3)
{
    char *addr_name = _fcdtocp(fcd_addr_name);
    char *port_name = _fcdtocp(fcd_port_name);
    char *mest      = _fcdtocp(fcd_mesg);

#else
(addr_name,n1,port_name,n2,mesg,n3)
char	addr_name[];
Integer	*n1;
char	port_name[];
Integer	*n2;
char	mesg[];
Integer	*n3;
{

#endif

    int sockfd;
    struct addrinfo hints, *servinfo, *p;
    int rv;
    int numbytes;


    memset(&hints, 0, sizeof(hints));
    hints.ai_family = AF_UNSPEC;
    hints.ai_socktype = SOCK_DGRAM;

    if ((rv = getaddrinfo(addr_name, port_name, &hints, &servinfo)) != 0) {
        fprintf(stderr, "getaddrinfo: %s\n", gai_strerror(rv));
    }

    // loop through all the results and make a socket
    for(p = servinfo; p != NULL; p = p->ai_next) {
        if ((sockfd = socket(p->ai_family, p->ai_socktype,
                p->ai_protocol)) == -1) {
            perror("talker: socket");
            continue;
        }

        break;
    }

    if (p == NULL) {
        fprintf(stderr, "talker: failed to bind socket\n");
    }

    if ((numbytes = sendto(sockfd, mesg, strlen(mesg), 0,
             p->ai_addr, p->ai_addrlen)) == -1) {
        perror("talker: sendto");
        exit(1);
    }

    freeaddrinfo(servinfo);

    printf("nwpw_talker: sent %d bytes to %s:%s\n", numbytes, addr_name,port_name);
    close(sockfd);

}

