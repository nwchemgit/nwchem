/* $Id: util_mygabcast.c 25519 2014-04-24 22:23:18Z d3y133 $ */
/* routine drop i/o caches */
#define _XOPEN_SOURCE 600
#include <unistd.h>
#include <fcntl.h>
#include "ga.h"
#include "macdecls.h"
#include "typesf2c.h"

Integer fortchar_to_string(const char *, Integer, char *, const Integer);

void FATR util_fadvise_(const char *fort_fname, double *offset, Integer *length, int flen){
    char buf[1024];
    int code, tmp;
#if defined(LINUX) || defined(LINUX64)
    if (!fortchar_to_string(fort_fname, flen, buf, sizeof(buf)))
      GA_Error("util_fadvise: fortchar_to_string failed for fname",0);
    
    int fd = open(buf, O_RDWR);
    
    printf(" fadvise fd %d fname %s\n", fd, fort_fname);
    
    (void) fdatasync(fd);
    
    if(posix_fadvise( fd, (size_t) offset,(size_t) length,POSIX_FADV_DONTNEED)!=0) perror("fadvise");
    close(fd);
#endif
}
