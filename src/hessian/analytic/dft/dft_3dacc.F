      subroutine dft_3dacc(g_fock,buf,imat,
     ,     ifirst,ilast,jfirst,jlast,ldin)
C$Id: dft_3dacc.F,v 1.2 2002-01-14 20:40:53 edo Exp $
      implicit none
      integer g_fock
      integer imat
      integer ldin
      integer ifirst,ilast,jfirst,jlast
      double precision buf(*)
c      
      integer lo(3)  !array of starting indices for ga   [input]
      integer hi(3)  !  array of ending indices for ga     [input]
      integer ld(2) !array specifying leading dimensions for buffer 
      ld(1)=1  
      ld(2)=ldin
      lo(1)=imat
      hi(1)=imat
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
      if(ndim.ne.3) call errquit(' 3d3d3d ',ndim)
      write(6,*) ' 3dg imat ',imat
      if(imat.lt.1) call errquit( ' imat neg ',imat)
      if(imat.gt.dims(1)) call errquit(' imat gt dim1 ',imat)
      if(ilast.gt.dims(2)) call errquit(' il gt dim2 ',dims(2))
      if(jlast.gt.dims(3)) call errquit(' jl gt dim3 ',dims(3))
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
