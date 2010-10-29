      subroutine pme_start(a,m,imffti,nodpmi,
     + ngxi,ngyi,ngzi,nwm,nwa,nsa,ictrl,nbgeti)
c $Id$
      implicit none
c
      integer nodpmi,ngxi,ngyi,ngzi,imffti,nwm,nwa,nsa,m
      integer ictrl,nbgeti
      real*8 a
c
      integer mem
c
#include "pme_common.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      me=ga_nodeid()
      np=ga_nnodes()
c
      pi=four*atan(one)
      twopi=two*pi
      alpha=a
      morder=m
      nodpme=nodpmi
      icntrl=ictrl
      nbget=nbgeti
c
      imfft=imffti
      lcorr(1)=.false.
      lcorr(2)=.false.
      lcorr(3)=.false.
c
      mwa=nwa
      mwm=nwm
      msa=nsa
c
      ngx=ngxi
      ngy=ngyi
      ngz=ngzi
      ngmax=max(ngx,ngy,ngz)
      mgx=ngx
      mgy=ngy
      mgz=(ngz/nodpme)+1
      ngrx=ngx+morder
      ngry=ngy+morder
      ngrz=ngz
      mgrx=ngrx
      mgry=ngry
      mgrz=mgz
      if(nodpme*mgz.lt.ngz)
     + call md_abort('Inconsistent number of pme nodes',0)
c
c     allocate memory for spline coefficients
c
      if(.not.ma_push_get(mt_dbl,3*ngmax,'bmod',l_bmod,i_bmod))
     + call md_abort('Failed to allocate bmod',0)
      if(.not.ma_push_get(mt_dbl,ngmax,'barr',l_barr,i_barr))
     + call md_abort('Failed to allocate barr',0)
c
c     allocate memory for the grids
c
      mem=2*mgx*mgy*mgz
      if(.not.ma_push_get(mt_dbl,mem,'grid1',l_grd1,i_grd1))
     + call md_abort('Failed to allocate grid1',0)
      if(.not.ma_push_get(mt_dbl,mem,'grid2',l_grd2,i_grd2))
     + call md_abort('Failed to allocate grid2',0)
      if(.not.ma_push_get(mt_dbl,mgrx*mgry*mgz,'grid',l_grd,i_grd))
     + call md_abort('Failed to allocate grid',0)
      if(nbget.gt.0) then
      if(.not.ma_push_get(mt_dbl,mgrx*mgry*mgz,'gridt',l_grdt,i_grdt))
     + call md_abort('Failed to allocate gridt',0)
      if(.not.ma_push_get(mt_int,mgz,'gridh',l_grdh,i_grdh))
     + call md_abort('Failed to allocate gridh',0)
      else
      if(.not.ma_push_get(mt_dbl,1,'gridt',l_grdt,i_grdt))
     + call md_abort('Failed to allocate gridt',0)
      if(.not.ma_push_get(mt_int,1,'gridh',l_grdh,i_grdh))
     + call md_abort('Failed to allocate gridh',0)
      endif
c
c     allocate memory for node arrays
c
      if(.not.ma_push_get(mt_int,np,'lnode',l_lnod,i_lnod))
     + call md_abort('Failed to allocate lnode',0)
      if(.not.ma_push_get(mt_int,np,'kfrom',l_kfr,i_kfr))
     + call md_abort('Failed to allocate kfrom',0)
      if(.not.ma_push_get(mt_int,np,'kto',l_kto,i_kto))
     + call md_abort('Failed to allocate kto',0)
      if(.not.ma_push_get(mt_int,np,'nodep',l_nodp,i_nodp))
     + call md_abort('Failed to allocate nodep',0)
c
      if(.not.ma_push_get(mt_int,ngz,'knode',l_knod,i_knod))
     + call md_abort('Failed to allocate knode',0)
c
      if(.not.ma_push_get(mt_log,mgz*np,'lsldo',l_lsld,i_lsld))
     + call md_abort('Failed to allocate ldldo',0)
      if(.not.ma_push_get(mt_int,4*mgz*np,'slmax',l_slmax,i_slmax))
     + call md_abort('Failed to allocate slmax',0)
c
c     allocate memory for fft node maps
c
      mem=3*max(ngx,ngy,ngz)
      if(.not.ma_push_get(mt_int,mem,'map1',l_map1,i_map1))
     + call md_abort('Failed to allocate map1',0)
      if(.not.ma_push_get(mt_int,mem,'map2',l_map2,i_map2))
     + call md_abort('Failed to allocate map2',0)
c
c     allocate memory for slab arrays
c
      mem=mgz*np
      if(.not.ma_push_get(mt_int,mem,'slab',l_slab,i_slab))
     + call md_abort('Failed to allocate slab',0)
c
c     create the global arrays for the slabs
c
      call ga_create_list(np*mgrx,mgry,'sl',mgrx,mgry,mgz,lslab)
c
c     create the flag array
c
      if(.not.ga_create(mt_int,np,5,'flg',1,5,ga_flg))
     + call md_abort('Failed to allocate global array flg',0)
      call pme_flag(0,0,me)
      call pme_flag(1,0,me)
c
c     select the fft
c
      call pme_select(zero)
c
c     test the fft
c
      if(lpnode)
     + call pme_fftest(dbl_mb(i_grd1),dbl_mb(i_grd2),int_mb(i_nodp))
c
      call pme_coeff(dbl_mb(i_bmod),dbl_mb(i_barr))
c
      return
      end
      subroutine pme_finish
c
      implicit none
c
#include "pme_common.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      if(.not.ga_destroy(ga_flg))
     + call md_abort('Failed to deallocate flg',0)
c
      call ga_destroy_list(mgz,lslab)
c
      if(.not.ma_pop_stack(l_slab))
     + call md_abort('Failed to deallocate slab',0)
      if(.not.ma_pop_stack(l_map2))
     + call md_abort('Failed to deallocate map2',0)
      if(.not.ma_pop_stack(l_map1))
     + call md_abort('Failed to deallocate map1',0)
      if(.not.ma_pop_stack(l_slmax))
     + call md_abort('Failed to deallocate slmax',0)
      if(.not.ma_pop_stack(l_lsld))
     + call md_abort('Failed to deallocate lsldo',0)
      if(.not.ma_pop_stack(l_knod))
     + call md_abort('Failed to deallocate knode',0)
      if(.not.ma_pop_stack(l_nodp))
     + call md_abort('Failed to deallocate nodep',0)
      if(.not.ma_pop_stack(l_kto))
     + call md_abort('Failed to deallocate kto',0)
      if(.not.ma_pop_stack(l_kfr))
     + call md_abort('Failed to deallocate kfrom',0)
      if(.not.ma_pop_stack(l_lnod))
     + call md_abort('Failed to deallocate lnode',0)
      if(.not.ma_pop_stack(l_grdh))
     + call md_abort('Failed to deallocate gridh',0)
      if(.not.ma_pop_stack(l_grdt))
     + call md_abort('Failed to deallocate gridt',0)
      if(.not.ma_pop_stack(l_grd))
     + call md_abort('Failed to deallocate grid',0)
      if(.not.ma_pop_stack(l_grd2))
     + call md_abort('Failed to deallocate grid2',0)
      if(.not.ma_pop_stack(l_grd1))
     + call md_abort('Failed to deallocate grid1',0)
c
      if(.not.ma_pop_stack(l_barr))
     + call md_abort('Failed to deallocate barr',0)
      if(.not.ma_pop_stack(l_bmod))
     + call md_abort('Failed to deallocate bmod',0)
c
      return
      end
      subroutine pme_select(value)
c
      implicit none
c
#include "pme_common.fh"
#include "mafdecls.fh"
c
      real*8 value
      integer i_trgt,l_trgt
c
      if(.not.ma_push_get(mt_dbl,np,'target',l_trgt,i_trgt))
     + call md_abort('Failed to allocate target',0)
c
      call pme_selnod(value,dbl_mb(i_trgt),int_mb(i_nodp),
     + int_mb(i_knod),int_mb(i_kfr),int_mb(i_kto),int_mb(i_map1))
c
      if(.not.ma_pop_stack(l_trgt))
     + call md_abort('Failed to deallocate target',0)
c
      return
      end
      subroutine pme_selnod(value,target,nodep,knode,kfrom,kto,map)
c
      implicit none
c
#include "pme_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      real*8 value,target(np)
      integer nodep(nodpme),knode(ngz),kfrom(np),kto(np),map(np)
c
      integer i,j,k,npls,nlnd
#if defined(ESSL) && defined(LAPI)
      integer number,npr,npc,myr,myc
      integer icntxt,lfftok
#endif
c
c     select nodes with minimum value to be used for pme/fft
c
      do 1 i=1,np
      target(i)=zero
    1 continue
      target(me+1)=value
      call ga_dgop(mrp_d01,target,np,'+')
c
      do 2 i=1,np
      nodep(i)=i
    2 continue
      do 3 i=1,np-1
      do 4 j=i+1,np
      if(target(i).lt.target(j)) then
      k=nodep(i)
      nodep(i)=nodep(j)
      nodep(j)=k
      endif
    4 continue
    3 continue
      do 5 i=1,np
      nodep(i)=nodep(i)-1
    5 continue
      do 6 i=1,nodpme-1
      do 7 j=i+1,nodpme
      if(nodep(i).gt.nodep(j)) then
      k=nodep(i)
      nodep(i)=nodep(j)
      nodep(j)=k
      endif
    7 continue
    6 continue
c
      call ga_brdcst(mrp_i01,nodep,nodpme*ma_sizeof(mt_int,1,mt_byte),0)
c
c     for node i: kfrom(i) is the first grid point in the z direction
c     for node i: kto(i)   is the last  grid point in the z direction
c
      do 8 i=1,np
      kfrom(i)=0
      kto(i)=0
    8 continue
c
      npls=(ngz/nodpme)+1
      nlnd=mod(ngz,nodpme)
      kfrom(nodep(1)+1)=1
      do 9 i=1,nodpme
      if(i.eq.nlnd+1) npls=npls-1
      kto(nodep(i)+1)=kfrom(nodep(i)+1)+npls-1
      if(i.lt.nodpme) kfrom(nodep(i+1)+1)=kto(nodep(i)+1)+1
    9 continue
c
      if(kto(nodep(nodpme)+1).ne.ngz)
     + call md_abort('Error in pme_select',me)
      kto(nodep(nodpme)+1)=ngz
c
c     for point i in the z direction: knode(i) is the owning node
c     ngzloc is the number of points in the z direction on this node
c
      j=1
      ngzloc=0
      do 10 i=1,ngz
      if(i.gt.kto(nodep(j)+1)) j=j+1
      knode(i)=nodep(j)
      if(nodep(j).eq.me) ngzloc=ngzloc+1
   10 continue
c
      lpnode=.false.
      do 11 i=1,nodpme
      if(me.eq.nodep(i)) lpnode=.true.
   11 continue
c
#if defined(ESSL) && defined(LAPI)
c
c     setup process grid for pessl grid
c
      if(imfft.eq.2) then
c
      if(lfftok.gt.0) then
      call blacs_gridexit(icntxt)
      else
      call blacs_pinfo(myr,myc)
      if(me.ne.myr) call md_abort('Identity crisis',me)
      if(np.ne.myc) call md_abort('Node count inconsistent',me)
      endif
c
      do 12 i=1,np
      map(i)=-1
   12 continue
      number=1
      do 13 i=1,ngz
      if(map(number).lt.0) then
      map(number)=knode(i)
      elseif(map(number).ne.knode(i)) then
      number=number+1
      map(number)=knode(i)
      endif
   13 continue
      if(number.ne.nodpme) call md_abort('Node assignment problem',me)
c
      call blacs_get(0,0,icntxt)
c
      call blacs_gridmap(icntxt,map,1,1,nodpme)
      call blacs_gridinfo(icntxt,npr,npc,myr,myc)
c
      if(lpnode.and.(myr.ge.npr.or.myc.ge.npc)) then
      call md_abort('arg_fft3d: fft initialization failed',me)
      endif
c
      endif
c
#endif
c
      return
      end
      subroutine pme_coeff(bmod,barr)
c
      implicit none
c
#include "pme_common.fh"
c
      real*8 bmod(ngmax,3),barr(ngmax)
c
      real*8 w,arr(25),darr(25)
      integer i
c
      w=zero
c
      call pme_splfil(w,arr,darr)
c
      do 1 i=1,ngmax
      barr(i)=zero
    1 continue
c
      do 2 i=2,morder+1
      barr(i)=arr(i-1)
    2 continue
c
      call pme_dftmod(bmod(1,1),barr,ngx)
      call pme_dftmod(bmod(1,2),barr,ngy)
      call pme_dftmod(bmod(1,3),barr,ngz)
c
      return
      end
      subroutine pme_splfil(w,arr,darr)
c      
      implicit none
c
#include "pme_common.fh"
c      
      real*8 arr(morder),darr(morder),w
c      
      integer i,j
      real*8 rinv
c
      arr(morder)=zero
      arr(2)=w
      arr(1)=one-w
      do 1 i=3,morder-1
      rinv=one/dble(i-1)
      arr(i)=rinv*w*arr(i-1)
      do 2 j=1,i-2
      arr(i-j)=rinv*((w+dble(j))*arr(i-j-1)+(dble(i-j)-w)*arr(i-j))
    2 continue
      arr(1)=rinv*(one-w)*arr(1)
    1 continue
      darr(1)=-arr(1)
      do 3 j=2,morder
      darr(j)=arr(j-1)-arr(j)
    3 continue
      rinv=one/dble(morder-1)
      arr(morder)=rinv*w*arr(morder-1)
      do 4 j=1,morder-2
      arr(morder-j)=rinv*((w+dble(j))*arr(morder-j-1)+
     + (dble(morder-j)-w)*arr(morder-j))
    4 continue
      arr(1)=rinv*(one-w)*arr(1)
c
      return
      end
      subroutine pme_dftmod(bmod,barr,ng)
c
      implicit none
c
#include "pme_common.fh"
c
      integer ng
      real*8 bmod(ng),barr(ng)
      integer i,j
      real*8 sum1,sum2,arg
c
      do 1 i=1,ng
      sum1=zero
      sum2=zero
      do 2 j=1,ng
      arg=twopi*dble((i-1)*(j-1))/dble(ng)
      sum1=sum1+barr(j)*cos(arg)
      sum2=sum2+barr(j)*sin(arg)
    2 continue
      bmod(i)=sum1*sum1+sum2*sum2
    1 continue
      do 3 i=1,ng
      if(bmod(i).lt.small) bmod(i)=half*(bmod(i-1)+bmod(i+1))
    3 continue
c
      return
      end
      subroutine pme_init()
c
      implicit none
c
#include "pme_common.fh"
#include "global.fh"
c
      call ga_zero_list(mgz,lslab)
      call pme_flag(0,0,me)
      call pme_flag(1,0,me)
      call ga_sync()
c
      return
      end
