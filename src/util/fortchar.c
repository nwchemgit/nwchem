/*$Id: fortchar.c,v 1.6 1999-11-13 03:15:13 bjohnson Exp $*/
/* Name munging to handle the various conventions for Fortran-C interfacing */
#if (defined(CRAY_T3D) || defined(ARDENT) || defined(WIN32))
#   define FCSND_  FCSND
#   define FCRCV_  FCRCV
#else
#   if (defined(AIX) || defined(NEXT) || defined(HPUX)) && !defined(EXTNAME)
#      define FCSND_  fcsnd
#      define FCRCV_  fcrcv
#   else
#      define FCSND_  fcsnd_
#      define FCRCV_  fcrcv_
#   endif
#endif

#define FC_MAXLEN 256 /* Length of message buffer.  Longer Fortran string
			 are sent in several messages. */

#include "sndrcv.h"
#include <string.h>
#ifdef CRAY_T3D
#include <fortran.h> /* Required for Fortran-C string interface on Crays */
#endif /* CRAY_T3D */
#ifdef WIN32
#include "typesf2c.h"
#endif

#define MAX(a,b) (((a) >= (b)) ? (a) : (b))
#define MIN(a,b) (((a) <= (b)) ? (a) : (b))
#define ABS(a) (((a) >= 0) ? (a) : (-(a)))

/**
 ** public routines
 **/


#if defined(CRAY_T3D) || defined(USE_FCD)
void FATR FCSND_(type, fcd, node, sync)
     long *type;
     _fcd fcd;
     long *node;
     long *sync;
#else /* CRAY_T3D */
void FCSND_(type, fstring, node, sync, flength)
     long *type;
     long *node;
     long *sync;
     char *fstring;
     long flength; /* Length of fstring, implicitly passed by FORTRAN */
#endif /* CRAY_T3D */
{
    char cstring[FC_MAXLEN];
    long fpos, clength, lenbuf=sizeof cstring;
#if defined(CRAY_T3D) || defined(USE_FCD)
    char	*fstring;	/* FORTRAN string */
    long	flength;	/* length of fstring */

    fstring = _fcdtocp(fcd);
    flength = _fcdlen(fcd);
#endif /* CRAY_T3D */
	    
    /* If the Fortran string is too long it is sent in chunks */
    fpos = 0;
    while (fpos < flength ) {
	clength = MIN(flength-fpos, lenbuf);
	strncpy(cstring, &fstring[fpos], (size_t) clength);
	fpos += clength;
	
	/* Our last message must be NULL-terminated.  We can add the
	   termination if there is room in the current buffer or send
	   this one and prepare a termination to be sent alone */

	if ( fpos == flength && clength < lenbuf ) {
	    clength++;
	    cstring[clength] = '\0';
	} else if ( fpos == flength && clength == lenbuf ) {
	    SND_(type, cstring, &clength, node, sync);
	    clength = 1;
	    cstring[0] = '\0';
	}

	SND_(type, cstring, &clength, node, sync);
    }
}

#if defined(CRAY_T3D) || defined(USE_FCD)
void FATR FCRCV_(type, fcd, flength, nodeselect, nodefrom, sync)
     long *type;
     _fcd fcd;
     long *flength;
     long *nodeselect;
     long *nodefrom;
     long *sync;
#else /* CRAY_T3D */
void FCRCV_(type, fstring, flength, nodeselect, nodefrom, sync, fsize)
     long *type;
     char *fstring;
     long *flength;
     long *nodeselect;
     long *nodefrom;
     long *sync;
     long fsize; /* Length of fstring, implicitly passed by FORTRAN */
#endif
{
    char cstring[FC_MAXLEN];
    long i, clength, lenbuf=sizeof cstring;
#if defined(CRAY_T3D) || defined(USE_FCD)
    char	*fstring;	/* FORTRAN string */
    long	fsize;	        /* length of fstring */

    fstring = _fcdtocp(fcd);
    fsize = _fcdlen(fcd);
#endif /* CRAY_T3D */
	    
    /* Collect the text, which may be broken into multiple messages */
    *flength = 0;
    do {
	RCV_(type, cstring, &lenbuf, &clength, nodeselect, nodefrom, sync);

	strncpy(&fstring[*flength], cstring,
		(size_t) MIN(fsize-*flength, clength) ); 
	*flength += MIN(fsize-*flength, clength);

	} while ( clength == lenbuf && cstring[clength] != '\0');

    /* The length includes the final NULL which Fortran doesn't want
       UNLESS we filled fstring up before getting the NULL */
    if ( *flength != fsize) (*flength)--;
    
    /* fill the rest of fstring with blanks --including the terminating NULL */
    for ( i = *flength; i < fsize; i++) 	fstring[i] = ' ';
}
