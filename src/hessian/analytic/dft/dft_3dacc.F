      subroutine dft_3dacc(g_fock,buf,imat,imat2,
     ,     ifirst,ilast,jfirst,jlast,ldin)
C$Id$
      implicit none
      integer g_fock
      integer imat,imat2
      integer ldin
      integer ifirst,ilast,jfirst,jlast
      double precision buf(*)
c      
      integer lo(3)  !array of starting indices for ga   [input]
      integer hi(3)  !  array of ending indices for ga     [input]
      integer ld(2) !array specifying leading dimensions for buffer 
      ld(1)=imat2-imat+1
      ld(2)=ldin
      lo(1)=imat
      hi(1)=imat2
      lo(2)=ifirst
      hi(2)=ilast
      lo(3)=jfirst
      hi(3)=jlast
      call nga_acc(g_fock, lo, hi, buf, ld, 1d0)
      return
      end
      subroutine dft_3dget(g_fock,buf,imat,
     ,     ifirst,ilast,jfirst,jlast,ldin)
      implicit none
#include "errquit.fh"
      integer g_fock
      integer imat
      integer ldin
      integer ifirst,ilast,jfirst,jlast
      double precision buf(*)
c      
      integer lo(3)  !array of starting indices for ga   [input]
      integer hi(3)  !  array of ending indices for ga     [input]
      integer ld(2) !array specifying leading dimensions for buffer 
#ifdef DEBUG
      integer ityp, ndim, dims(3)
      call nga_inquire(g_fock, ityp, ndim, dims)
      if(ndim.ne.3) call errquit(' 3d3d3d ',ndim, GA_ERR)
      write(6,*) ' 3dg imat ',imat
      if(imat.lt.1) call errquit( ' imat neg ',imat, UNKNOWN_ERR)
      if(imat.gt.dims(1)) call errquit(' imat gt dim1 ',imat,
     &       UNKNOWN_ERR)
      if(ilast.gt.dims(2)) call errquit(' il gt dim2 ',dims(2),
     &       UNKNOWN_ERR)
      if(jlast.gt.dims(3)) call errquit(' jl gt dim3 ',dims(3),
     &       UNKNOWN_ERR)
#endif
      ld(1)=1  
      ld(2)=ldin
      lo(1)=imat
      hi(1)=imat
      lo(2)=ifirst
      hi(2)=ilast
      lo(3)=jfirst
      hi(3)=jlast
      call nga_get(g_fock, lo, hi, buf, ld)
      return
      end
