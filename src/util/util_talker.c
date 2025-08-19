#include <stdio.h>
#include <stdlib.h>
#include <math.h>
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

void FATR util_talker_(char addr_name[], Integer * inet, Integer * n1, Integer * portin, Integer * sockout, Integer * max_retries_in, Integer * retry_delay_seconds_in)
{
#if defined(__MINGW32__)
    perror("util_talker: not coded for this architecture");
    exit(1);
#else
    int sock = 0;
    int max_retries = ((int) *max_retries_in);
    int retry_delay_seconds = ((int) *retry_delay_seconds_in);
    int retries = 0;

    if (max_retries <= 0){
        perror("Got an invalid (<=0) number of retries, resetting to 30\n");
        max_retries = 30;
    }
    if (retry_delay_seconds <= 0){
        perror("Got an invalid (<=0) delay for retries, resetting to 2\n");
        retry_delay_seconds = 2;
    }


    if (*inet > 0)
    {
        // --- TCP/IP Socket Logic ---
        struct sockaddr_in serv_addr;
        int port = ((int) *portin);
        int na   = ((int) *n1);
        char addr_name_buf[256]; // Use a safe buffer
        strncpy(addr_name_buf, addr_name, na);
        addr_name_buf[na] = '\0';

        printf("util_talker: Attempting to connect to TCP server %s:%d\n", addr_name_buf, port);

        if ((sock = socket(AF_INET, SOCK_STREAM, 0)) < 0)
        {
            perror("util_talker: TCP socket creation error");
            exit(1);
        }
        serv_addr.sin_family = AF_INET;
        serv_addr.sin_port = htons(port);

        if(inet_pton(AF_INET, addr_name_buf, &serv_addr.sin_addr) <= 0)
        {
            printf("\nutil_talker: Invalid address/Address not supported\n");
            exit(1);
        }

        // --- Retry Loop for TCP ---
        while (connect(sock, (struct sockaddr *)&serv_addr, sizeof(serv_addr)) < 0)
        {
            if (++retries >= max_retries) {
                fprintf(stderr, "\nutil_talker: TCP connection failed after %d attempts: %s. Exiting.\n", max_retries, strerror(errno));
                exit(1);
            }
            printf("util_talker: Connection failed. Retrying in %d second(s)... (%d/%d)\n",
                   retry_delay_seconds, retries, max_retries);
            sleep(retry_delay_seconds);
        }
    }
    else
    {
        // --- UNIX Domain Socket Logic ---
        struct sockaddr_un serv_addr;
        int na = ((int)*n1);
        char addr_name_buf[256];
        strncpy(addr_name_buf, addr_name, na);
        addr_name_buf[na] = '\0';

        memset(&serv_addr, 0, sizeof(struct sockaddr_un));
        serv_addr.sun_family = AF_UNIX;
        // Construct the full socket path, same as the server
        snprintf(serv_addr.sun_path, sizeof(serv_addr.sun_path), "/tmp/ipi_%s",
                 addr_name_buf);

        printf("util_talker: Attempting to connect to UNIX socket %s\n",
               serv_addr.sun_path);

        if ((sock = socket(AF_UNIX, SOCK_STREAM, 0)) < 0) {
            perror("util_talker: UNIX socket creation error");
            exit(1);
        }

        // --- Retry Loop for UNIX ---
        while (connect(sock, (struct sockaddr *)&serv_addr, sizeof(serv_addr)) <
               0) {
            if (++retries >= max_retries) {
                fprintf(stderr,
                        "\nutil_talker: UNIX socket connection failed after %d "
                        "attempts: %s. Exiting.\n",
                        max_retries, strerror(errno));
                exit(1);
            }
            printf("util_talker: Connection failed. Retrying in %d second(s)... (%d/%d)\n",
                   retry_delay_seconds, retries, max_retries);
            sleep(retry_delay_seconds);
        }
    }

    printf("util_talker: Connection successful! Socket ID: %d\n", sock);
    *sockout = ((Integer)sock);
#endif
}




void FATR util_talker_close_(Integer * socket1)
{
#if defined(__MINGW32__)
   perror("util_talker: not coded for this architecture");
   exit(1);
#else
   int sock = ((int) *socket1);
   close(sock);
#endif
}


void FATR util_talker_write_(Integer * socket1, char buffer[], Integer * n1)
{
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





void FATR util_talker_read_(Integer * socket1, char buffer[], Integer * n1)
{
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







