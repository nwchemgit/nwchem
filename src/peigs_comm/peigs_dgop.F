*
* $Id$
*
#define SIZE_OF_DOUBLE sizeof(double)
c     
      integer function peigs_dgop(ptype, xbuf, plen, op, scratch)
      integer ptype
      integer plen
      double precision xbuf(plen)
      double precision scratch(plen)
c     
      external op
      integer mxcombv1
      external mxcombv1
c     
      integer nprocs, iii
      integer proclist(512)
      integer nnodes
      external nnodes
c     
      nprocs = nnodes()
c     
      do i = 1, nprocs
         proclist(i) = i-1
      enddo
c     
      peigs_dgop = mxcombv1(xbuf, op, 8,
     $     plen, nprocs, proclist, ptype, scratch)
c     
c     
      return
      end

      
