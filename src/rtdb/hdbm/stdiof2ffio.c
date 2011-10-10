#include <stdio.h>
#include <fcntl.h>
#include <ffio.h>
#include <values.h>		/* MAXINT */

FILE *
stdiof_ffopen (const char *path, const char *mode)
{
	int	fd;
	int	oflag;
	long	cbits	= 0;
	int	cblks	= 0;
	struct ffsw	ffsw;
	static int	Trace_Open = MAXINT;

/*
 *	Check the environment variable FFIO_STDIO_TRACE_OPEN for its
 *	existence.  If it has been declared, then issue a msg to
 *	stderr about the file being opened via FFIO.
 */
	if ( Trace_Open == MAXINT ) {
		Trace_Open = (getenv("FFIO_STDIO_TRACE_OPEN")) ? 1 : 0;
	}

/*
 *	Based on the mode of the fopen, translate that into
 *	the oflags for the ffopen.
 */
	oflag = 0;
	switch ( mode[0] ) {
	case 'r':	oflag |= O_RDONLY;			break;
	case 'w':	oflag |= O_WRONLY | O_CREAT | O_TRUNC;	break;
	case 'a':	oflag |= O_WRONLY | O_CREAT | O_APPEND;	break;
	default:
		return ( NULL );
	}
	if ( (mode[1] == '+') || (mode[2] == '+')){
		oflag &= ~(O_RDONLY | O_WRONLY);
		oflag |= O_RDWR;
	}

/*
 *	Open the file.
 */
#ifndef	JUNK
	fd = ffopen (path, oflag, 0666);
#else
	oflag |= O_DIRECT;
	fd = ffopens(path, oflag, 0666, cbits, cblks, &ffsw, "cache:512:4");
#endif
	if ( Trace_Open == 1 ) {
		fprintf(stderr
			,"STDIO_FFOPEN: %3d = ffopen([%s], 0x%.4x, %03o);\n"
			,fd,path,oflag,0666);
	}
	if ( fd < 0 ) {
		return ( NULL );
	}


	return ( (FILE *)fd );
}



int
stdiof_ffclose(FILE *stream)
{
	int	fd;
	int	rc;

	fd = (long)stream;

	rc = ffclose(fd);

	return ( rc );
}




size_t
stdiof_ffread (void *buf, size_t recsz, size_t numrec, FILE *stream)
{
	int	fd;
	size_t	reqdbytes;
	size_t	repval;
	size_t	recdone;

	fd = (long)stream;	/* redfine it to be what it actually is */
	reqdbytes = recsz * numrec;

	repval = ffread(fd, buf, reqdbytes);
	recdone = repval / recsz;

	return ( recdone );
}




size_t
stdiof_ffwrite(const void *buf, size_t recsz, size_t numrec, FILE *stream)
{
	int	fd;
	size_t	reqdbytes;
	size_t	repval;
	size_t	recdone;

	fd = (long)stream;	/* redfine it to be what it actually is */
	reqdbytes = recsz * numrec;

	repval = ffwrite(fd, (void *)buf, reqdbytes);
	recdone = repval / recsz;

	return ( recdone );
}



int
stdiof_ffseek(FILE *stream, long offset, int whence)
{
	int	fd;
	int	rc;

	fd = (long)stream;

	rc = ffseek(fd, offset, whence);

	return ( (rc<0) ? -1 : 0  );
}



int
stdiof_fflush(FILE *stream)
{
	int	fd;
	int	rc;
	struct ffsw	ffsw;

	fd = (long)stream;

	rc = ffflush(fd, &ffsw);

	return ( rc );
}



void
stdiof_ffrewind(FILE *stream)
{
	int	fd;
	int	rc;
	struct ffsw	ffsw;

	fd = (long)stream;

	rc = ffflush(fd, &ffsw);
	rc = ffseek (fd, 0L, SEEK_SET);

}



long
stdiof_fftell(FILE *stream)
{
	int	fd;
	off_t	cpos;

	fd = (long)stream;

	cpos = ffseek(fd, 0L, SEEK_CUR);
	
	return ( (long)cpos );

}
/* $Id$ */
