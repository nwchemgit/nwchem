      subroutine ana_edinit()
c
c $Id: ana_covar.F,v 1.6 2000-03-24 05:41:35 d3j191 Exp $
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
      if(.not.ga_create(mt_dbl,3*nsel,3*nsel,'covar',0,3*nsel,ga_cov))
     + call errquit('Failed to create global array covar',me)
      call ga_zero(ga_cov)
      if(.not.ga_create(mt_dbl,3*nsel,3*nsel,'vector',0,3*nsel,ga_vec))
     + call errquit('Failed to create global array vector',me)
c
c     allocate memory
c
      if(.not.ma_push_get(mt_int,3*nsel,'ndx',l_ndx,i_ndx))
     + call errquit('Failed to allocate memory for ndx',me)
      if(.not.ma_push_get(mt_dbl,3*nsel,'cov',l_cov,i_cov))
     + call errquit('Failed to allocate memory for cov',me)
      if(.not.ma_push_get(mt_int,3*nsel,'ord',l_ord,i_ord))
     + call errquit('Failed to allocate memory for ord',me)
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
      implicit none
c
#include "ana_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "global.fh"
c
c     deallocate memory
c
      if(.not.ma_pop_stack(l_ord))
     + call errquit('Failed to deallocate memory for ord',me)
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
#include "msgids.fh"
c
      integer ndx(*)
      real*8 cov(*)
      real*8 x(msa,3),xa(msa,3)
      integer isel(msa)
c
      integer i,j,ia,ja,ix,jx
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
      do 3 j=jl,jh
      ja=ndx((j+2)/3)
      jx=mod(j+2,3)+1
      cov(j-jl+1)=(x(ia,ix)-xa(ia,ix))*(x(ja,jx)-xa(ja,jx))
    3 continue
      call ga_acc(ga_cov,i,i,jl,jh,cov,1,1.0)
    2 continue
c
      call ga_sync()
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
      real*8 vct(nsel),val,pnorm,vval
      real*8 x(nsa,3),xa(nsa,3),xp(nsa,3)
      integer isel(nsa)
c
      integer i,j,k
c
      val=0.0d0
      pnorm=0.0d0
      k=0
      do 1 i=1,nsel
      vval=0.0d0
      do 2 j=1,3
      k=k+1
      val=val+(x(ndx(i),j)-xa(ndx(i),j))*vct(k)
      vval=vval+(x(ndx(i),j)-xa(ndx(i),j))*vct(k)
      pnorm=pnorm+vct(k)*vct(k)
    2 continue
    1 continue
      val=val/sqrt(pnorm)
c
      k=0
      do 3 i=1,nsel
      do 4 j=1,3
      k=k+1
      xp(ndx(i),j)=val*vct(k)+xa(ndx(i),j)
    4 continue
    3 continue
c
      return
      end
