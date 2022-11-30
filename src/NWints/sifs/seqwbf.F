      subroutine seqwbf( unit, buffer, buflen )
c
c----------------------------------------------------------
c this is a primitive i/o routine to be used by all higher
c level routines requiring i/o of this type. machine i/o
c peculiarites are meant to reside in these routines only.
c----------------------------------------------------------
c
c routine name: seqwbf
c version: 1.0                        date: 8/24/88
c author: eric stahlberg - ohio state university
c purpose: this routine will write a given vector to a
c          sequential file. the length of the vector is given
c          in minimum necessary real*8 words to access the entire
c          vector
c parameters:
c     unit:   unit number to be written
c     buffer: buffer to be written
c     buflen: number of real*8 words to write to file
c
      implicit integer(a-z)
c
      integer unit, buflen
      real*8 buffer ( buflen )
c
#if defined(DEBUG)
      print *,'seqwbf: unit=, buflen=',unit,buflen
#endif 
#if (defined(DEBUG) && defined(OSU))
      print *, (buffer(i),i=1, min(100,buflen) )
#endif 
c
      write (unit) buffer
      return
      end
