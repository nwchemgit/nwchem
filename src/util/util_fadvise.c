/* $Id$ */
/* routine drop i/o caches */
#define _FILE_OFFSET_BITS 64
#include "typesf2c.h"
#ifdef __linux__
#define _XOPEN_SOURCE 600
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#endif
#if  defined(POSIX_FADV_DONTNEED)
#include "ga.h"
#include "macdecls.h"

Integer fortchar_to_string(const char *, Integer, char *, const Integer);
void FATR util_fadvise(const char *, int);

void FATR util_fadvise_dontneed_(const char *fort_fname,  int flen){
    char buf[1024];
    if (!fortchar_to_string(fort_fname, flen, buf, sizeof(buf)))
      GA_Error("util_fadvise: fortchar_to_string failed for fname",0);

  (void) FATR util_fadvise(buf, POSIX_FADV_DONTNEED);
}

void FATR util_fadvise_noreuse_(const char *fort_fname,  int flen){
    char buf[1024];
    if (!fortchar_to_string(fort_fname, flen, buf, sizeof(buf)))
      GA_Error("util_fadvise: fortchar_to_string failed for fname",0);

  (void) FATR util_fadvise(buf, POSIX_FADV_NOREUSE);
}



void FATR util_fadvise(const char *buf, int mode){
    
    int fd = open(buf, O_RDWR);
    
    //    printf(" fadvise fd %d fname %s\n", fd, fort_fname);
    struct stat fd_stat ;

    if ( fstat( fd, &fd_stat ) < 0 ) {
        perror( "Could not stat file: " );
        return;
    }
    off_t offset = 0;
    off_t length = fd_stat.st_size;

    
    /*    (void) fdatasync(fd);*/
    if(posix_fadvise( fd,  offset, length, mode)!=0) perror("fadvise");
    close(fd);

}
#else
/* stubs */
void FATR util_fadvise_dontneed_(const char *fort_fname,  int flen){
}
void FATR util_fadvise_noreuse_(const char *fort_fname,  int flen){
}
#endif
