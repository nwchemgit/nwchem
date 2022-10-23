      subroutine aiopen( unit, filnam, length )
c
c----------------------------------------------------------
c this is a primitive i/o routine to be used by all higher
c level routines requiring i/o of this type. machine i/o
c peculiarites are meant to reside in these routines only.
c----------------------------------------------------------
c
c routine name: aiopen
c version: 1.0                        date: 10/13/88
c author: eric stahlberg - ohio state university
c purpose: this routine will open a file for asynchronous i/o.
c          record length is in real*8 words for machines
c          which require this value.
c          if an error occurs, bummer is called with an
c          encoded error number: ioerr*100+unit
c parameters:
c     filnam: external file name to associate with unit
c     unit:    unit number for file to be opened
c     length:  length of direct access records in real*8 words
c
      implicit integer(a-z)
c
      integer unit,length
      character*(*) filnam
c
c     # bummer error types
      integer wrnerr, nfterr, faterr
      parameter (wrnerr=0, nfterr=1, faterr=2)
c
      integer ioerr
c
#if defined(DEBUG)
      print *, 'aiopen: unit=,filnam=,length=',
     + unit,filnam,' ',length
#endif 
#if defined(REFERENCE)
c      asynchronous i/o for this case is not implmented because
c      calling program calls in wrong order. should be fixed
c      when calling programs are rewritten.
      open(unit=unit,form='unformatted',status='unknown',
     + sync='asynchronous',file=filnam,iostat=ioerr,err=900)
#else 
c
c open regular file
c
      open(unit=unit,access='sequential',form='unformatted',status=
     +  'unknown',file=filnam,iostat=ioerr,err=900)
#endif 
      return
c
 900  continue
      call bummer('i/o error in openda, (ioerr*100+unit)=',
     & (ioerr*100+unit), faterr )
      end
