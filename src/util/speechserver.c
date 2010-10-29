/*
 $Id$
 */

#include <stdio.h>

/*
  Test udp server program ... 

    1) Create udp socket and bind to address ... print
       out the port number so that client can use it

    2) Sit in infinite loop doing recvfrom on socket ...
       data is expected to be null terminated strings.
       If the string is "exit" then exit.
*/

int main()
{
  int sock, port, len;
  char buf[2048];

  if (udp_server(&sock, &port) != 0) {
    perror("udp_server");
    return 1;
  }

  (void) printf("UDP server port=%d\n",port);

  while ((len = udp_recv(sock, buf, 2048)) >= 0) {
    (void) printf("UDP server: len=%d, data='%s'\n", len, buf);
    if (strcmp(buf, "exit") == 0) {
      (void) close(sock);
      return 0;
    }
  }

  perror("recvfrom failed");
  (void) printf("UDP server: error from udp_recv, len=%d\n",len);

  (void) close(sock);
  return 1;
}
