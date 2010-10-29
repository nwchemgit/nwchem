      subroutine argos_pme_fft(nd1,nd2,nd3,isgn,x1,x2,
     $     map1,map2,knode,itype)
c $Id$
      implicit none
c
#include "argos_pme_common.fh"
#include "mafdecls.fh"
c
      integer nd1,nd2,nd3,isign,map1(*),map2(*),itype
      real*8 x1(*),x2(*)
      integer knode(ngz)
c
      integer ignr
c
c     3D-FFT interface
c
c     nd1,nd2,nd3  : dimension of x1 and x2
c     isign        : 1=forward, -1=reverse, -2=reverse without normalization
c     x1           : input slabwise distributed complex array
c     x2           : output slabwise distributed complex array
c     map1         : id of processor owning i-th plane
c     map2         : id of processor owning i-th plane after transpose
c     liwork,iwork : integer scratch space
c     lrwork,rwork : real*8 scratch space
c     itype        : 0=standard 3d-fft, 1=specific 3d-fft if available
c
      integer i,ndim,jstat,isgn
      real*8 scale
      integer liwrk,lrwrk,i_rwrk,l_rwrk,i_iwrk,l_iwrk
c
c     IBM-SP2 specific
c
#if defined(ESSL)
      integer icntxt,ip(40)
#endif
c
      ndim=nd1*nd2*nd3
      isign=isgn
      scale=one
      if(isgn.eq.-1) scale=one/dble(ndim)
      if(isgn.lt.-1) isign=-1
c
c     brand name 3d-fft routines
c
#if defined(ESSL) && defined(LAPI)
      if(itype.gt.1) then
      ip(1)=1
      ip(2)=1
      ip(20)=nd1
      ip(21)=nd2
      ip(22)=nd1
      ip(23)=nd2
c
      call timer_start(206)
      if(lpnode) call pdcft3(x1,x2,nd1,nd2,nd3,isign,scale,icntxt,ip)
      call timer_stop(206)
c
      return
      endif
#endif
      if(itype.gt.1)
     + call md_abort('argos_pme_fft: 3D-pFFt no implemented',itype)
c
c     generic 3d-fft routine
c
      do 1 i=1,ngz
      map1(i)=knode(i)
      map2(i)=knode(i)
    1 continue
c
      if(ngzloc.gt.0) then
      do 2 i=1,nd1*nd2*mgz*2
      x2(i)=x1(i)
    2 continue
c
      call wrkspc2(1,ngx,ngy,ngz,map1,map2,liwrk,lrwrk,ignr)
cx
      liwrk=2*liwrk
      lrwrk=2*lrwrk
cx
c
      if(.not.ma_push_get(mt_dbl,lrwrk,'rwork',l_rwrk,i_rwrk))
     + call md_abort('Failed to allocate rwork',0)
c
      if(.not.ma_push_get(mt_int,liwrk,'iwork',l_iwrk,i_iwrk))
     + call md_abort('Failed to allocate iwork',0)
c
      call timer_start(206)
      call pfft3d(isign,nd1,nd2,nd3,x2,map1,map2,
     + lrwrk,dbl_mb(i_rwrk),liwrk,int_mb(i_iwrk),jstat)
      call timer_stop(206)
c
      if(.not.ma_verify_allocator_stuff())
     +     call md_abort('FFT buffer problems after',me)
c      call ffflush(6)
c
      if(.not.ma_pop_stack(l_iwrk))
     + call md_abort('Failed to deallocate iwork',0)
      if(.not.ma_pop_stack(l_rwrk))
     + call md_abort('Failed to deallocate rwork',0)
c
      if(isign.eq.-1) then
      do 3 i=1,nd1*nd2*mgz*2
      x2(i)=x2(i)*scale
    3 continue
      endif
c
      if(jstat.ne.0)
     + call md_abort('argos_pme_fft: fft failed on node',me)
      endif
c
      return
      end
 
