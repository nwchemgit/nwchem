      subroutine dft_3dacc(g_fock,buf,imat,
     ,     ifirst,ilast,jfirst,jlast,ldin)
C$Id: dft_3dacc.F,v 1.1 2001-11-16 02:10:55 edo Exp $
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
