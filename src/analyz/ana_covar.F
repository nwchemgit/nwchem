      subroutine ana_edinit()
c
c $Id: ana_covar.F,v 1.5 2000-03-21 22:23:50 d3j191 Exp $
c
      implicit none
c
#include "ana_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "global.fh"
c
c     create nsel x nsel global array covariance matrix
c
      if(.not.ga_create(mt_dbl,3*nsel,3*nsel,'covar',1,3*nsel,ga_cov))
     + call errquit('Failed to create global array covar',me)
      call ga_zero(ga_cov)
      if(.not.ga_create(mt_dbl,3*nsel,3*nsel,'vector',1,3*nsel,ga_vec))
     + call errquit('Failed to create global array vector',me)
c
c     allocate memory
c
      if(.not.ma_push_get(mt_int,3*nsel,'ndx',l_ndx,i_ndx))
     + call errquit('Failed to allocate memory for ndx',me)
      if(.not.ma_push_get(mt_dbl,3*nsel,'cov',l_cov,i_cov))
     + call errquit('Failed to allocate memory for cov',me)
c
      if(me.eq.0) then
      write(*,1000) 3*nsel
 1000 format(' Memory allocated for covariance ',i5)
      endif
c
      return
      end
      subroutine ana_edfinal()
c
c $Id: ana_covar.F,v 1.5 2000-03-21 22:23:50 d3j191 Exp $
c
      implicit none
c
#include "ana_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "global.fh"
c
c     deallocate memory
c
      if(.not.ma_pop_stack(l_cov))
     + call errquit('Failed to deallocate memory for cov',me)
      if(.not.ma_pop_stack(l_ndx))
     + call errquit('Failed to deallocate memory for ndx',me)
c
c     destroy global array covariance matrix
c
      if(.not.ga_destroy(ga_vec))
     + call errquit('Failed to destroy global array vector',me)
      if(.not.ga_destroy(ga_cov))
     + call errquit('Failed to destroy global array covar',me)
c
      return
      end
      subroutine ana_covar(ndx,cov,isel,x,xa)
c
      implicit none
c
#include "ana_common.fh"
#include "global.fh"
c
      integer ndx(nsel)
      real*8 cov(nsel)
      real*8 x(nsa,3),xa(nsa,3)
      integer isel(nsa)
c
      integer i,j,k,ii,ia,ib,ix
      integer il,ih,jl,jh
c
      j=0
      do 1 i=1,nsa
      if(isel(i).ne.0) then
      j=j+1
      ndx(j)=i
      endif
    1 continue
c
      call ga_distribution(ga_cov,me,il,ih,jl,jh)
c
      do 2 i=il,ih
      ia=ndx((i+2)/3)
      ix=mod(i+2,3)+1
      ii=0
      do 3 j=1,nsel
      ib=ndx(j)
      do 4 k=1,3
      ii=ii+1
      cov(ii)=(x(ia,ix)-xa(ia,ix))*(x(ib,k)-xa(ib,k))
      if(ia.eq.1.and.ib.eq.1.and.ix.eq.1.and.k.eq.1) then
c      write(*,'(7f12.6)') x(ia,ix),xa(ia,ix),
c     + x(ib,k),xa(ib,k),cov(ii)
      endif
    4 continue
    3 continue
      call ga_acc(ga_cov,i,i,jl,jh,cov,1,1.0)
    2 continue
c
      return
      end
      subroutine ana_projec(ndx,vct,isel,x,xa,val,xp)
c
      implicit none
c
#include "ana_common.fh"
#include "global.fh"
c
      integer ndx(nsel)
      real*8 vct(nsel),val
      real*8 x(nsa,3),xa(nsa,3),xp(nsa,3)
      integer isel(nsa)
c
      integer i,j,k
c
      val=0.0d0
      k=0
      do 1 i=1,nsel
      do 2 j=1,3
      k=k+1
      val=(x(ndx(i),j)-xa(ndx(i),j))*vct(k)
      xp(ndx(i),j)=val*vct(k)+xa(ndx(i),j)
    2 continue
    1 continue
c
      return
      end
