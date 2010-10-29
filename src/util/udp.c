/*
 $Id$
 */

#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <stdio.h>

int udp_server(int *sock, int *port)
/*
  Create a datagram socket in the internet domain (i.e. udp)
  and bind it to a port. Return the socket file descriptor and
  the port number in the argument list.
  Returned is 0 if all is OK or -1 if an error is detected.
*/
{
  int len;
  struct sockaddr_in name;

  if((*sock = socket(AF_INET, SOCK_DGRAM, 0)) < 0)
    return -1;

  name.sin_family = AF_INET;
  name.sin_addr.s_addr = INADDR_ANY;
  name.sin_port = 0;
  if (bind(*sock, (struct sockaddr *) &name, sizeof name) < 0)
    return -1;

  len = sizeof(name);
  if (getsockname(*sock, (struct sockaddr *) &name, &len) < 0)
    return -1;

  *port = ntohs(name.sin_port);
  return 0;
}

int udp_recv(int sock, char *buf, int lenbuf)
/*
  Read data from the unconnected udp socket sock. lenbuf is
  the size of buf.
  Returned is the amount of data actually read or -1 if
  an error was detected. Note that udp datagrams are forced
  to be small (typicall <=2048 bytes) by the protocol.
*/
{
  int zero = 0;

  return recvfrom(sock, buf, lenbuf, 0, 
		  (struct sockaddr *) NULL, &zero);
}

int udp_send(const char *hostname, int port, const char *buf, int lenbuf)
/*
  Send a udp packet of lenbuf bytes from buffer buf to the
  specified port and host. 
  Returned is the number of bytes actually sent or -1 if
  any error was detected.
*/
{
  int sock, len;
  struct sockaddr_in name;
  struct hostent *hp, *gethostbyname();

  if((sock = socket(AF_INET, SOCK_DGRAM, 0)) < 0)
    return -1;

  if((hp = gethostbyname(hostname)) == (struct hostent *) 0)
    return -1;

  bcopy((char *)hp->h_addr, (char *)&name.sin_addr, hp->h_length);
  name.sin_family = AF_INET;
  name.sin_port = htons((u_short) port);

  len = sendto(sock, buf, lenbuf, 0,
               (struct sockaddr *) &name, sizeof name);

  (void) close(sock);

  return len;
}
