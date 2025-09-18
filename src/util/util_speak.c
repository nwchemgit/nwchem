#include <stdio.h>
#include <string.h>
#ifdef WIN32
#include "typesf2c.h"
#endif

typedef long integer;		/* FORTRAN integer type */

#ifdef SPEECH
extern int udp_send(const char *hostname, int port, const char *buf, int);

static char *hostname;		/* Speech server hostname */
static int port;		/* Speech server port number */
#endif

void util_speak_init(const char *host, const int p)
/*
  Get info to connect to remote speech server
  */
{
#ifdef SPEECH
    if (hostname) free(hostname); /* In case of multiple inits */
    if (!(hostname = strdup(host)))
	return;

    port = p;
#endif
}

void util_speak(const char *string)
/*
  Application interface to speech server

  String is a null terminated character string
  */
{
#ifdef SPEECH
    if (!hostname) return;	/* Not initialized */

    (void) udp_send(hostname, port, string, strlen(string));
#endif
}

void util_speak_(const char *string, int len)
{
#ifdef SPEECH
    char buf[256];

    if (!fortchar_to_string(string, len, buf, sizeof(buf)))
	return;

    util_speak(buf);
#endif
}

void util_speak_init_(const char *host, integer *fp, int len)
{
#ifdef SPEECH
    char buf[256];
    int p = (int) *fp;

    if (!fortchar_to_string(host, len, buf, sizeof(buf)))
	return;
    
    util_speak_init(buf, p);
#endif
}

