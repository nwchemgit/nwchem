/*
 $Id$
 */

#include <stdio.h>

/*
  Simple udp socket client. Sends third command line argument
  as a udp datagram to the server specified by the host/port
  combination of the first two arguments.

  1) Invoke udp_get to run in the background on the host
     that U want to use as the server.

  2) Udp_get prints out the port number that the client (udp_put)
     needs.

  3) Invoke udp_put multiple times with the format

     udp_put hostname port 'data'

     where hostname is the host where udp_get is running.

     e.g. udp_put titan3 1182 'Hello there'

  4) When U get bored either specify "exit" as the data or
     just kill udp_get.

*/

static void Error(string, integer)
  char *string;
  int integer;
{
  (void) fflush(stdout);
  (void) fprintf(stderr, "!! error: %s %d\n", string, integer);
  perror("system message");
  exit(1);
}

int main(argc, argv)
  int argc;
  char **argv;
{
  int port, lenbuf;
  char *hostname, *buf;

  if (argc != 4)
    Error("usage: udp_put hostname port 'Message to server'", argc);

  hostname = argv[1];
  port = atoi(argv[2]);
  buf = argv[3];
  lenbuf = strlen(buf) + 1;

  if (udp_send(hostname, port, buf, lenbuf) != lenbuf) 
    Error("udp_put: udp_send failed",-1);

  return 0;
}
