/*$Id$*/
#include <stdio.h>
#include <sys/types.h>
#include <unistd.h>
#include <fcntl.h>

extern off_t lseek(int, off_t, int);
extern int open(const char *, int, int);
extern int write(int, const char *, int);
extern int read(int, char *, int);

#define MAX_FD 256
static off_t length[MAX_FD];
static off_t position[MAX_FD];

int ipsc_open(char *name, int flags, int mode)
{
  int fd;

  if (flags & O_APPEND) {
    (void) fprintf(stderr, "ipsc_open: append not supported\n");
    return -1;
  }

  if ((fd = open(name, flags, mode)) < 0)
    return fd;

  position[fd] = 0;

  if ((length[fd] = lseek(fd, (off_t) 0, SEEK_END)) < 0) {
    (void) fprintf(stderr, "ipsc_lseek: setting length failed\n");
    (void) fflush(stderr);
    return -1;
  }

  (void) lseek(fd, (off_t) 0, SEEK_SET);

#ifdef IPSC_IO_DEBUG
  (void) fprintf(stderr, "ipsc_open: %s fd=%d, length=%d\n", name, fd, length[fd]);
#endif

  return fd;
}

int ipsc_close(fd)
{
  length[fd] = position[fd] = 0;

  return close(fd);
}

off_t ipsc_lseek(int fd, off_t offset, int whence)
{
  static char buf[1024];	/* static variables set to zero */

  /*
    On the CFS Intel fails to return zeroes as the result of a read
    from a hole in a file.  ipsc_lseek(), ipsc_write() , ipsc_read(),
    ipsc_open, ensure that holes are explicitly filled in with zeroes.
    */
  
  if (whence != SEEK_SET) {
    (void) fprintf(stderr, "ipsc_lseek: only seek_set supported\n");
    (void) fflush(stderr);
    return (off_t) -1;
  }
  
#ifdef IPSC_IO_DEBUG
  /* Temporary sanity check */
  
  if (length[fd] != lseek(fd, (off_t) 0, SEEK_END)) {
    (void) fprintf(stderr, "ipsc_lseek: length is inconsistent\n");
    (void) fflush(stderr);
    return -1;
  }
#endif
  
  if (length[fd] < offset) {
    /* About to position beyond position EOF ... fill in the gap */

    int ndo = offset - length[fd];
    
    if (length[fd] != lseek(fd, (off_t) 0, SEEK_END)) {
      (void) fprintf(stderr, "ipsc_lseek: length is inconsistent\n");
      (void) fflush(stderr);
      return -1;
    }
    
#ifdef IPSC_IO_DEBUG
    fprintf(stderr, "ipsc_lseek: padding fd=%d from %d to %d\n",
	    fd, length[fd], offset);
    fflush(stdout);
#endif

    while (ndo) {
      int nwrite = (ndo < sizeof(buf)) ? ndo : sizeof(buf);
      if (write(fd, buf, nwrite) != nwrite) {
	(void) fprintf(stderr, "ipsc_lseek: write of zeroes failed\n");
	(void) fflush(stderr);
	return -1;
      }
      ndo -= nwrite;
    }
  
    length[fd] = offset;
  }

  position[fd] = offset;

  return lseek(fd, offset, whence);
}

int ipsc_write(int fd, char *buf, int nbytes)
{
  position[fd] += nbytes;

  if (position[fd] > length[fd]) length[fd] = position[fd];

  return write(fd, buf, nbytes);
}

int ipsc_read(int fd, char *buf, int nbytes)
{
  position[fd] += nbytes;

  if (position[fd] > length[fd]) position[fd] = length[fd];

  return read(fd, buf, nbytes);
}
