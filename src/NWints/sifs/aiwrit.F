      subroutine aiwrit( unit, buffer, buflen )
c
c----------------------------------------------------------
c this is a primitive i/o routine to be used by all higher
c level routines requiring i/o of this type. machine i/o
c peculiarites are meant to reside in these routines only.
c----------------------------------------------------------
c
c routine name: aiwrit
c version: 1.0                        date: 10/13/88
c author: eric stahlberg - ohio state university
c purpose: this routine will write a record to an asynchronous file
c          buffer length is in real*8 words.
c          if an error occurs, bummer is called with an
c          encoded error number: ioerr*100+unit
c parameters:
c     unit:    unit number for file to be opened
c     buffer:  buffer to read record into
c     buflen:  length of buffer to be read
c
      implicit integer(a-z)
c
      integer unit,buflen
      real*8 buffer(buflen)
c
c     # bummer error types
      integer wrnerr, nfterr, faterr
      parameter (wrnerr=0, nfterr=1, faterr=2)
c
      integer ioerr
c
#if defined(DEBUG)
      print *, 'aiwrit: unit=, buflen=',unit,buflen
#endif 
#if (defined(DEBUG) && defined(OSU))
      print *, (buffer(i), i=1, min(100,buflen) )
#endif 
#if defined(REFERENCE)
c      asynchronous i/o for this case is not implmented because
c      calling program calls in wrong order. should be fixed
c      when calling programs are rewritten.
      write (unit=unit,iostat=ioerr,err=900,wait=.false.) buffer
#else 
c
c write regular file
c
      write (unit=unit,iostat=ioerr,err=900) buffer
#endif 
      return
c
 900  continue
      call bummer('i/o error in airead, (ioerr*100+unit)=',
     & (ioerr*100+unit), faterr )
      end
