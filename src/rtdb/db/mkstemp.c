/*$Id$*/
#include <sys/types.h>
#include <stdio.h>
#include <string.h>
#include <fcntl.h>

#ifdef IPSC
extern long mynode();
#define PID ((int) mynode())
#else
extern int getpid(void);
#define PID getpid()
#endif

int mkstemp(char *template)
/*
  Crude version of mkstemp
*/
{
  int flags = O_RDWR | O_CREAT;
  int len = (int) strlen(template);
  char *alphabet="ABCDEFGHIJKLMNOPQRSTUVWXYZ";
  int pid = PID;
  int letter, fd;

  if (pid > 100000)
    return -1;

  for (letter=0; letter<26; letter++) {
    (void) sprintf(template+len-6, "%05d%c", pid, alphabet[letter]);

    if ((fd = open(template, flags, 0660)) >= 0)
      return fd;
  }

  return -1;
}

