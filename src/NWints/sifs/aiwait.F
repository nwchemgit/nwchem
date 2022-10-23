      subroutine aiwait( unit )
c
c----------------------------------------------------------
c this is a primitive i/o routine to be used by all higher
c level routines requiring i/o of this type. machine i/o
c peculiarites are meant to reside in these routines only.
c----------------------------------------------------------
c
c routine name: aiwait
c version: 1.0                        date: 10/13/88
c author: eric stahlberg - ohio state university
c purpose: this routine will wait for an asynchronous file
c
c parameters:
c     unit:    unit number for file to be opened
c
      implicit integer(a-z)
c
      integer unit,ioerr
c
c     # bummer error types
      integer wrnerr, nfterr, faterr
      parameter (wrnerr=0, nfterr=1, faterr=2)
#if defined(DEBUG)
      print *, 'aiwait: unit=',unit
#endif 
c
c
#if defined(REFERENCE)
      wait(unit,iostat=ioerr)
      if (ierr.ne.0) goto 900
#endif 
      return
c
 900  call bummer('aiwait: i/o error, (ioerr*100+unit)=',
     & (ioerr*100+unit), faterr )
      end
