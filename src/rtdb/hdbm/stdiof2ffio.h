


extern FILE	*stdiof_ffopen (const char *path, const char *mode);
extern int	stdiof_ffclose(FILE *stream);
extern size_t	stdiof_ffread (void *buf, size_t recsz, size_t numrec, FILE *stream);
extern size_t	stdiof_ffwrite(const void *buf, size_t recsz, size_t numrec, FILE *stream);
extern int	stdiof_ffseek(FILE *stream, long offset, int whence);
extern long	stdiof_fftell(FILE *stream);
extern int	stdiof_fflush(FILE *stream);
extern void	stdiof_ffrewind(FILE *stream);



#define	fopen	stdiof_ffopen
#define	fclose	stdiof_ffclose
#define	fread	stdiof_ffread
#define	fwrite	stdiof_ffwrite
#define	fseek	stdiof_ffseek
#define	ftell	stdiof_fftell
#define	rewind	stdiof_ffrewind
/* $Id$ */
