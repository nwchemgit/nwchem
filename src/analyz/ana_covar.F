      subroutine ana_edinit()
c
c $Id$
c
      implicit none
c
#include "ana_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "global.fh"
c
      if(me.eq.0) then
      write(*,1000) (ma_sizeof(mt_dbl,3*nsel,mt_byte)+
     + ma_sizeof(mt_int,6*nsel,mt_byte))/104857,
     + ma_inquire_avail(mt_byte)/1048576,
     + ma_sizeof(mt_dbl,18*nsel*nsel/np,mt_byte)/1048576,
     + ga_memory_avail()/1048576
 1000 format(/,' Memory   Required    Available',//,
     + '   ma   ',2i10,' MB',/,
     + '   ga   ',2i10,' MB',/)
      endif
      print*,'nsel=',nsel
c
c     create nsel x nsel global array covariance matrix
c
      if(.not.ga_create(mt_dbl,3*nsel,3*nsel,'covar',0,3*nsel,ga_cov))
     + call md_abort('Failed to create global array covar',me)
      call ga_zero(ga_cov)
      if(.not.ga_create(mt_dbl,3*nsel,3*nsel,'vector',0,3*nsel,ga_vec))
     + call md_abort('Failed to create global array vector',me)
c
c     allocate memory
c
      if(.not.ma_push_get(mt_int,3*nsel,'ndx',l_ndx,i_ndx))
     + call md_abort('Failed to allocate memory for ndx',me)
      if(.not.ma_push_get(mt_dbl,3*nsel,'cov',l_cov,i_cov))
     + call md_abort('Failed to allocate memory for cov',me)
      if(.not.ma_push_get(mt_int,3*nsel,'ord',l_ord,i_ord))
     + call md_abort('Failed to allocate memory for ord',me)
c
      if(me.eq.0) then
      write(*,1002) 3*nsel
 1002 format(' Memory allocated for covariance ',i5)
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
     + call md_abort('Failed to deallocate memory for ord',me)
      if(.not.ma_pop_stack(l_cov))
     + call md_abort('Failed to deallocate memory for cov',me)
      if(.not.ma_pop_stack(l_ndx))
     + call md_abort('Failed to deallocate memory for ndx',me)
c
c     destroy global array covariance matrix
c
      if(.not.ga_destroy(ga_vec))
     + call md_abort('Failed to destroy global array vector',me)
      if(.not.ga_destroy(ga_cov))
     + call md_abort('Failed to destroy global array covar',me)
c
      return
      end
      subroutine ana_covar(ndx,cov,isel,x,xa)
c
      implicit none
c
#include "ana_params.fh"
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
c      if(i.eq.1.and.j.eq.1) then
c      write(*,'(5i5,2f12.6,2i5,2f12.6)') me,i,j,
c     + ia,ix,xa(ia,ix),x(ia,ix),
c     + ja,jx,xa(ja,jx),x(ja,jx)
c      endif
    3 continue
      call ga_acc(ga_cov,i,i,jl,jh,cov,1,one)
    2 continue
c
      call ga_sync()
      return
      end
      subroutine ana_projec(ndx,vct,isel,x,xa,mx,val,xp)
c
      implicit none
c
#include "ana_common.fh"
#include "global.fh"
#include "mafdecls.fh"
c
      integer mx
      integer ndx(nsel)
      real*8 vct(nsel),val,pnorm,vval
      real*8 x(mx,3),xa(mx,3),xp(nsa,3)
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
      if(valmin.gt.valmax) then
      call ana_wrtpdb(lfnmax,filmax,byte_mb(i_snam),isel,xp,
     + isel,xp,.true.,.false.,int_mb(i_ndxw))
      valmax=val
      call ana_wrtpdb(lfnmin,filmin,byte_mb(i_snam),isel,xp,
     + isel,xp,.true.,.false.,int_mb(i_ndxw))
      valmin=val
      elseif(val.lt.valmin) then
      call ana_wrtpdb(lfnmin,filmin,byte_mb(i_snam),isel,xp,
     + isel,xp,.true.,.false.,int_mb(i_ndxw))
      valmin=val
      elseif(val.gt.valmax) then
      call ana_wrtpdb(lfnmax,filmax,byte_mb(i_snam),isel,xp,
     + isel,xp,.true.,.false.,int_mb(i_ndxw))
      valmax=val
      endif
c
      return
      end
