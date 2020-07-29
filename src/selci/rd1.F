*
* $Id$
*
C> \ingroup selci
C> @{
C>
C> \brief Routines to read columbus integrals
C>
      subroutine selci_rd1(iflint,h,norbs)
      implicit real*8 (a-h,o-z)
      parameter(lenbuf=2047, n1eabf=1636, n2eabf=1364)
      common/blkin/buffer(lenbuf),ibuf(n2eabf*4)
      dimension h(*)
c
c     get the 1 electron hamiltonian (t+v)
c
c     assume that we are positioned just after the header record
c
c     skip over overlap matrix
c
      call selci_skip1(iflint,lenbuf,buffer)
c
c  t integrals:
c
      nnorbs = norbs*(norbs+1)/2
      call dfill(nnorbs,0.0d0,h,1)
      call selci_rdlab1(iflint,h,lenbuf,n1eabf,buffer,ibuf,norbs)
c
c  v integrals:
c
      call selci_rdlab1(iflint,h,lenbuf,n1eabf,buffer,ibuf,norbs)
c
      end
      subroutine selci_skip1(ntape,lbuf,buf)
c
c skip over a triangle of 1-e integrals on disk
c
      implicit real*8(a-h,o-z)
      real*8 buf(lbuf)
      integer info(4)
c
10    read (ntape) buf
      call selci_ulab16(buf(lbuf),info,4)
      last=info(3)
      if (last.eq.0) goto 10
      return
      end
      subroutine selci_rdlab1(ntape,s,lbuf,nipbuf,buf,ibuf,nbft)
c
c read a triangle of 1-e integrals from disk
c
      implicit real*8(a-h,o-z)
#include "errquit.fh"
      real*8 buf(lbuf),s(nbft*(nbft+1)/2)
      integer ibuf(nipbuf*2)
      integer info(4)
c
      nnbft = nbft*(nbft+1)/2
c
10    read (ntape) buf
      call selci_ulab16(buf(lbuf),info,4)
      last=info(3)
      num=info(4)
      call selci_ulab8(buf(nipbuf+1),ibuf,2*num)
      ic2 = 1
      do 20 ic = 1,num
        ii = ibuf(ic2)
        jj = ibuf(ic2+1)
        ic2 = ic2 + 2
        i = max(ii,jj)
        j = min(ii,jj)
        if (i.gt.nbft) goto 20
        ij = i*(i-1)/2 + j
        if (ij.le.0 .or. ij.gt.nnbft) then
          write(6,*) ' rdlab1 ',i,j,ij,nbft,nnbft
c         call abort
          call errquit(' invalid triangle address ',ij, INPUT_ERR)
        endif
        s(ij) = s(ij) + buf(ic)
20    continue
      if (last.eq.0) goto 10
      return
      end
      subroutine selci_colint(ntape,q,iq,jq,kq,lq,nw)
c
c     read an unpack a buffer of two electron integrals and labels
c
      implicit real*8(a-h,o-z)
      parameter(lenbuf=2047, n1eabf=1636, n2eabf=1364)
      common/blkin/buffer(lenbuf),ibuf(n2eabf*4)
      integer info(4)
      dimension q(n2eabf),iq(n2eabf),jq(n2eabf),kq(n2eabf),lq(n2eabf)
      save last
      data last/0/
c
      if (last.ne.0) then
         nw = 0
         return
      endif
      read(ntape) buffer
      call selci_ulab16(buffer(lenbuf),info,4)
      last = info(3)
      num = info(4)
      nnum = (num+1)/2
      call selci_ulab8(buffer(n2eabf+1),ibuf,4*num)
      i4 = 1
      do 20 i = 1,num
         q(i) = buffer(i)
         iq(i) = ibuf(i4  )
         jq(i) = ibuf(i4+1)
         kq(i) = ibuf(i4+2)
         lq(i) = ibuf(i4+3)
          i4 = i4 + 4
20    continue
c
      nw = num
c
      end
      subroutine selci_ulab8(p,u,nuw)
c
c  unpack integral labels from p(*) into u(*,*).
c
c  p(*) = packed array (working precision in the calling program).
c  u(*) = unpacked array.  u(1: ((nuw+7)/8)*8 ) are referenced.
c  nuw  = number of unpacked integral labels.
c
c  written by ron shepard.
c  version: 14-mar-89
c
      implicit integer(a-z)
*mdc*if cray
*      real*8 p(*)
*      integer u(nuw)
*c
*      nuw8=((nuw+7)/8)*8
*      call unpack(p,8,u,nuw8)
*mdc*elseif fps
*      real*8 p(*)
*      integer u(nuw)
*c
*      npw=(nuw+7)/8
*      call viup8(p,1,u,1,npw)
*mdc*elseif vax
*c  byte order is just for consistency with existing conventions. -rls
*      byte p(8,*)
*      integer u(8,*)
*      intrinsic zext
*c  use the zext() function for 8-bit unsigned integer moves.
*      nuw8=(nuw+7)/8
*      do 10 i=1,nuw8
*         u(1,i)=zext(p(2,i))
*         u(2,i)=zext(p(1,i))
*         u(3,i)=zext(p(4,i))
*         u(4,i)=zext(p(3,i))
*         u(5,i)=zext(p(6,i))
*         u(6,i)=zext(p(5,i))
*         u(7,i)=zext(p(8,i))
*         u(8,i)=zext(p(7,i))
*10    continue
*mdc*elseif sunos3
*c  logical*1 does not work on sun OS3, so use character*1 kludge.
*      character*1 p(8,*)
*      character*1 u(4,8,*)
*c
*      call izero(nuw,u,1)
*      do 20 i=1,((nuw+7)/8)
*         u(4,1,i)=p(1,i)
*         u(4,2,i)=p(2,i)
*         u(4,3,i)=p(3,i)
*         u(4,4,i)=p(4,i)
*         u(4,5,i)=p(5,i)
*         u(4,6,i)=p(6,i)
*         u(4,7,i)=p(7,i)
*         u(4,8,i)=p(8,i)
*20    continue
*mdc*else
c  general byte-addressable 32-bit integer machines.
      logical*1 p(8,*)
      logical*1 u(4,8,*)
c
c  to determine the little-endian or big-endian addressing convention.
      integer longw
      integer*2 shortw(2)
#ifdef F77      
      equivalence (longw,shortw)
#endif      
      save longw, shortw
c
      call ifill(nuw,0,u,1)
      longw=1
#ifndef F77      
      shortw=transfer(longw,1)
#endif            
      if(shortw(1).eq.1)then
c        ...little-endian.
         do 10 i=1,((nuw+7)/8)
            u(1,1,i)=p(8,i)
            u(1,2,i)=p(7,i)
            u(1,3,i)=p(6,i)
            u(1,4,i)=p(5,i)
            u(1,5,i)=p(4,i)
            u(1,6,i)=p(3,i)
            u(1,7,i)=p(2,i)
            u(1,8,i)=p(1,i)
10       continue
      else
c        ...big-endian.
         do 20 i=1,((nuw+7)/8)
            u(4,1,i)=p(1,i)
            u(4,2,i)=p(2,i)
            u(4,3,i)=p(3,i)
            u(4,4,i)=p(4,i)
            u(4,5,i)=p(5,i)
            u(4,6,i)=p(6,i)
            u(4,7,i)=p(7,i)
            u(4,8,i)=p(8,i)
20       continue
      endif
*mdc*endif
      return
      end
      subroutine selci_ulab16(p,u,nuw)
c
c  unpack one-electron integral labels from p(*) into u(*,*).
c
c  p(*) = packed array (working precision in the calling program).
c  u(*) = unpacked array.  u(1: ((nuw+3)/4)*4 ) are referenced.
c  nuw  = number of unpacked integral labels.
c
c  written by ron shepard.
c  version: 14-mar-89
c
      implicit integer(a-z)
*mdc*if cray
*      real*8 p(*)
*      integer u(nuw)
*c
*      nuw16=((nuw+3)/4)*4
*      call unpack(p,16,u,nuw16)
*mdc*elseif fps
*      real*8 p(*)
*      integer u(nuw)
*c
*      npw=(nuw+3)/4
*      call viup16(p,1,u,1,npw)
*mdc*elseif vax
*c  byte order is just for consistency with existing conventions. -rls
*      integer*2 p(4,*)
*      integer u(4,*)
*      intrinsic zext
*c  use the zext() function for 16-bit unsigned integer moves.
*      nuw16=(nuw+3)/4
*      do 10 i=1,nuw16
*         u(1,i)=zext(p(1,i))
*         u(2,i)=zext(p(2,i))
*         u(3,i)=zext(p(3,i))
*         u(4,i)=zext(p(4,i))
*10    continue
*mdc*else
c  general byte-addressable 32-bit integer machines.
      integer*2 p(4,*)
      integer*2 u(2,4,*)
c
c  to determine the little-endian or big-endian addressing convention.
      integer longw
      integer*2 shortw(2)
#ifdef F77      
      equivalence (longw,shortw)
#endif      
      save longw, shortw
c
      call ifill(nuw,0,u,1)
      longw=1
#ifndef F77      
      shortw=transfer(longw,1)
#endif            
      if(shortw(1).eq.1)then
c        ...little-endian.
         do 10 i=1,((nuw+3)/4)
            u(1,1,i)=p(4,i)
            u(1,2,i)=p(3,i)
            u(1,3,i)=p(2,i)
            u(1,4,i)=p(1,i)
10       continue
      else
c        ...big-endian.
         do 20 i=1,((nuw+3)/4)
            u(2,1,i)=p(1,i)
            u(2,2,i)=p(2,i)
            u(2,3,i)=p(3,i)
            u(2,4,i)=p(4,i)
20       continue
      endif
*mdc*endif
      return
      end
C> @}
