      subroutine seqrbf( unit, buffer, buflen )
c
c----------------------------------------------------------
c this is a primitive i/o routine to be used by all higher
c level routines requiring i/o of this type. machine i/o
c peculiarites are meant to reside in these routines only.
c----------------------------------------------------------
c
c routine name: seqrbf
c version: 1.0                        date: 8/24/88
c author: eric stahlberg - ohio state university
c purpose: this routine will read a given vector from a
c          sequential file. the length of the vector is given
c          in minimum necessary real*8 words to access the entire
c          vector
c parameters:
c     unit:   unit number to read from
c     buffer: buffer to transfer information to
c     buflen: number of real*8 words to read from file
c
      implicit integer(a-z)
c
      integer unit, buflen
      real*8 buffer ( buflen )
c
#if defined(DEBUG)
      print *,'seqrbf: unit=, buflen=',unit,buflen
#endif 
c
      read (unit) buffer
c
#if (defined(DEBUG) && defined(OSU))
      print *, (buffer(i), i=1, min(100,buflen) )
#endif 
      return
      end
