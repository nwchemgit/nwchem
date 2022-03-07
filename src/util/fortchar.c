/*$Id$*/
/* Name munging to handle the various conventions for Fortran-C interfacing */
#if (defined(ARDENT) || defined(WIN32)) && !defined(__MINGW32__)
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

#include <tcgmsg.h>
#include <string.h>
#ifdef WIN32
#include "typesf2c.h"
#endif

#define MAX(a,b) (((a) >= (b)) ? (a) : (b))
#define MIN(a,b) (((a) <= (b)) ? (a) : (b))
#define ABS(a) (((a) >= 0) ? (a) : (-(a)))

/**
 ** public routines
 **/


void FCSND_(long *type, char *fstring, long *node, long *sync,
            long flength) /* Length of fstring, implicitly passed by FORTRAN */
{
    char cstring[FC_MAXLEN];
    long fpos, clength, lenbuf=sizeof(cstring);
	    
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
	    tcg_snd(*type, cstring, clength, *node, *sync);
	    clength = 1;
	    cstring[0] = '\0';
	}

	tcg_snd(*type, cstring, clength, *node, *sync);
    }
}

void FCRCV_( long *type, char *fstring, long *flength,
             long *nodeselect, long *nodefrom, long *sync,
             long fsize) /* Length of fstring, implicitly passed by FORTRAN */
{
    char cstring[FC_MAXLEN];
    long i, clength, lenbuf=sizeof(cstring);
	    
    /* Collect the text, which may be broken into multiple messages */
    *flength = 0;
    do {
	tcg_rcv(*type, cstring, lenbuf, &clength, *nodeselect, nodefrom, *sync);

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
