      subroutine pme_chgrid(ipset,xw,nwmloc,xs,isq,isgm,nsaloc,epme)
c $Id$
      implicit none
c
#include "pme_common.fh"
#include "mafdecls.fh"
c
      integer ipset,nwmloc,nsaloc
      real*8 xw(mwm,3,mwa),xs(msa,3),epme
      integer isq(msa),isgm(msa)
c
      call timer_start(24)
c
      call pme_init()
c
      if(.not.lcorr(ipset)) then
      call cf_pmecor(ipset,xw,nwmloc,xs,isq,isgm,nsaloc,ecorr(ipset))
      lcorr(ipset)=.true.
      endif
      call timer_stop(24)
c
      call pme_grid(ipset,xw,nwmloc,xs,isq,isgm,nsaloc,
     + int_mb(i_nodp),int_mb(i_knod),int_mb(i_kfr),int_mb(i_kto),
     + int_mb(i_lnod),log_mb(i_lsld),int_mb(i_slmax),
     + dbl_mb(i_grd1),dbl_mb(i_grd2),dbl_mb(i_grd),
     + int_mb(i_grdh),dbl_mb(i_grdt))
c
      epme=ecorr(ipset)
c
      if(lpnode) then
c
      call timer_start(26)
c
      call pme_fft(ngx,ngy,ngz,-2,dbl_mb(i_grd2),dbl_mb(i_grd1),
     + int_mb(i_map1),int_mb(i_map2),int_mb(i_knod),imfft)
c
      call timer_stop(26)
      call timer_start(27)
c
      call cf_epme(int_mb(i_knod),int_mb(i_kfr),int_mb(i_kto),
     + dbl_mb(i_bmod),dbl_mb(i_grd1),mgz,epme)
c
      epme=epme+ecorr(ipset)
c
      call timer_stop(27)
      call timer_start(28)
c
      call pme_sync(nodpme,int_mb(i_nodp))
c
      call timer_stop(28)
      call timer_start(31)
c
      call pme_fft(ngx,ngy,ngz,1,dbl_mb(i_grd1),dbl_mb(i_grd2),
     + int_mb(i_map1),int_mb(i_map2),int_mb(i_knod),imfft)
c
      call timer_stop(31)
      call timer_start(32)
c
      call pme_putgrd(dbl_mb(i_grd2),dbl_mb(i_grd))
c
      call timer_stop(32)
c
c     synchronize pme nodes
c
      call timer_start(33)
c
      call pme_sync(nodpme,int_mb(i_nodp))
c
      call timer_stop(33)
c
c     set flag on all nodes to indicate that potential
c     grid calculation has completed
c
      call timer_start(34)
c
      if(me.eq.int_mb(i_nodp)) call pme_flag(0,1,-1)
c
      call timer_stop(34)
      endif
c
      call timer_start(34)
      if(icntrl.eq.1) call ga_sync()
c
      call timer_stop(34)
c
      return
      end
      subroutine pme_energy(ipset,xw,nwmloc,xs,isq,isgm,nsaloc,epme)
c
      implicit none
c
#include "pme_common.fh"
#include "mafdecls.fh"
c
      integer ipset,nwmloc,nsaloc
      real*8 xw(mwm,3,mwa),xs(msa,3),epme
      integer isq(msa),isgm(msa)
c
      call timer_start(24)
      call pme_init()
c
      if(.not.lcorr(ipset)) then
      call cf_pmecor(ipset,xw,nwmloc,xs,isq,isgm,nsaloc,ecorr(ipset))
      lcorr(ipset)=.true.
      endif
      call timer_stop(24)
c
      call pme_grid(ipset,xw,nwmloc,xs,isq,isgm,nsaloc,
     + int_mb(i_nodp),int_mb(i_knod),int_mb(i_kfr),int_mb(i_kto),
     + int_mb(i_lnod),log_mb(i_lsld),int_mb(i_slmax),
     + dbl_mb(i_grd1),dbl_mb(i_grd2),dbl_mb(i_grd),
     + int_mb(i_grdh),dbl_mb(i_grdt))
c
      epme=ecorr(ipset)
c
      if(.not.lpnode) return
c
      call timer_start(28)
c
      call pme_fft(ngx,ngy,ngz,-2,dbl_mb(i_grd2),dbl_mb(i_grd1),
     + int_mb(i_map1),int_mb(i_map2),int_mb(i_knod),imfft)
c
      call timer_stop(28)
      call timer_start(29)
c
      call cf_epme(int_mb(i_knod),int_mb(i_kfr),int_mb(i_kto),
     + dbl_mb(i_bmod),dbl_mb(i_grd1),mgz,epme)
c
      epme=epme+ecorr(ipset)
c
      call timer_stop(29)
c
      return
      end
      subroutine pme_grid(ipset,xw,nwmloc,xs,isq,isgm,nsaloc,
     + nodep,knode,kfrom,kto,lnode,lsldo,maxsl,gridr,gridc,grid,
     + lgridh,gridt)
c
      implicit none
c
#include "pme_common.fh"
#include "global.fh"
#include "mafdecls.fh"
c
      integer ipset,nwmloc,nsaloc
      real*8 xw(mwm,3,mwa),xs(msa,3)
      integer isq(msa),isgm(msa)
      logical lsldo(np,mgz)
      integer maxsl(4,np,mgz)
      integer nodep(nodpme),knode(ngz),kfrom(np),kto(np),lnode(np)
      real*8 gridr(mgx,mgy,mgz),gridc(2,mgx,mgy,mgz)
      real*8 grid(mgrx,mgry,mgrz)
      integer lgridh(mgz)
      real*8 gridt(mgrx,mgry,mgrz)
c
      integer node,inode,ngznon
      integer i,j,k,ix,iy,ig,jg,minx,maxx,miny,maxy,minz,maxz
      integer ilsl,ihsl,jlsl,jhsl
      logical lset,lgrid
c
      call timer_start(25)
c
      lset=.true.
c
      if(nbget.gt.0) then
      do 321 i=1,mgz
      lgridh(i)=0
  321 continue
      endif
c
c     loop over all nodes
c
      do 1 inode=1,np
      node=me+inode
      if(node.gt.np) node=node-np
      node=node-1
      lnode(node+1)=0
c
c     determine the number of planes on node
c
      ngznon=0
      if(kto(node+1).gt.0) ngznon=kto(node+1)-kfrom(node+1)+1
      if(ngznon.eq.0) goto 1
c
      call cf_chgrid(lset,ipset,xw,nwmloc,xs,isq,isgm,nsaloc,
     + grid,ngznon,lgrid,knode,kfrom,lsldo,node,
     + minx,maxx,miny,maxy,minz,maxz)
c
      if(lgrid) then
      lnode(node+1)=1
c
      call timer_stop(25)
      call timer_start(26)
c
      if(nbget.eq.0) then
      do 21 k=1,ngznon
      if(lsldo(node+1,k)) then
      call ga_distribution(lslab(k),node,ilsl,ihsl,jlsl,jhsl)
      call ga_acc(lslab(k),ilsl+minx-1,ilsl+maxx-1,jlsl+miny-1,
     + jlsl+maxy-1,grid(minx,miny,k),mgrx,one)
      maxsl(1,node+1,k)=minx
      maxsl(2,node+1,k)=maxx
      maxsl(3,node+1,k)=miny
      maxsl(4,node+1,k)=maxy
      else
      maxsl(1,node+1,k)=0
      maxsl(2,node+1,k)=0
      maxsl(3,node+1,k)=0
      maxsl(4,node+1,k)=0
      endif
   21 continue
      else
      do 224 i=1,mgz
      if(lgridh(i).gt.0) call ga_nbwait(lgridh(i))
      lgridh(i)=0
  224 continue
      do 221 k=1,ngznon
      if(lsldo(node+1,k)) then
      do 222 iy=miny,maxy 
      do 223 ix=minx,maxx
      gridt(ix,iy,k)=grid(ix,iy,k)
  223 continue
  222 continue
      call ga_distribution(lslab(k),node,ilsl,ihsl,jlsl,jhsl)
      call ga_nbacc(lslab(k),ilsl+minx-1,ilsl+maxx-1,jlsl+miny-1,
     + jlsl+maxy-1,gridt(minx,miny,k),mgrx,one,lgridh(k))
      maxsl(1,node+1,k)=minx
      maxsl(2,node+1,k)=maxx
      maxsl(3,node+1,k)=miny
      maxsl(4,node+1,k)=maxy
      else
      maxsl(1,node+1,k)=0
      maxsl(2,node+1,k)=0
      maxsl(3,node+1,k)=0
      maxsl(4,node+1,k)=0
      lgridh(k)=0
      endif
  221 continue
      endif
      call timer_stop(26)
      call timer_start(25)
c
      endif
c
    1 continue
c
      if(nbget.gt.0) then
      do 333 i=1,mgz
      if(lgridh(i).gt.0) call ga_nbwait(lgridh(i))
  333 continue
      endif
c
c     pmenodes should wait until all nodes have finished contributions to the grid
c
      call timer_stop(25)
      call timer_start(27)
c
c     the following pme_flag/pme_wait/pme_flag combination may be
c     replaced by a ga_sync() at loss of performance
c
      call pme_flag(1,-1,-1)
      if(lpnode) call pme_wait(1,np)
c
c     retrieve local slab from global array space
c
      call timer_stop(27)
      call timer_start(26)
      if(ngzloc.gt.0) then
      do 23 k=ngzloc,1,-1
      call ga_distribution(lslab(k),me,ilsl,ihsl,jlsl,jhsl)
      call ga_get(lslab(k),ilsl,ilsl+ngrx-1,jlsl,jlsl+ngry-1,
     + grid(1,1,k),mgrx)
      do 125 j=1,ngy
      do 126 i=1,ngx
      gridr(i,j,k)=zero
  126 continue
  125 continue
      do 123 jg=1,ngry
      do 124 ig=1,ngrx
      i=ig
      j=jg
      if(i.gt.ngx) i=i-ngx
      if(j.gt.ngy) j=j-ngy
      gridr(i,j,k)=gridr(i,j,k)+grid(ig,jg,k)
  124 continue
  123 continue
      do 24 j=ngy,1,-1
      do 25 i=ngx,1,-1
      gridc(1,i,j,k)=gridr(i,j,k)
      gridc(2,i,j,k)=zero
   25 continue
   24 continue
   23 continue
      endif
      call timer_stop(26)
c
      return
      end
      subroutine pme_putgrd(gridc,grid)
c
      implicit none
c
#include "pme_common.fh"
#include "global.fh"
c
      real*8 gridc(2,mgx,mgy,*),grid(mgrx,mgry)
c
      integer i,j,k,m,n,ilsl,ihsl,jlsl,jhsl
c
      do 1 k=1,ngzloc
      do 2 i=1,ngx
      do 3 j=1,ngy
      grid(i,j)=gridc(1,i,j,k)
    3 continue
    2 continue
      do 4 m=1,morder
      do 5 i=1,ngx
      grid(i,ngy+m)=grid(i,m)
    5 continue
      do 6 j=1,ngy
      grid(ngx+m,j)=grid(m,j)
    6 continue
      do 7 n=1,morder
      grid(ngx+m,ngy+n)=grid(m,n)
    7 continue
    4 continue
      call ga_distribution(lslab(k),me,ilsl,ihsl,jlsl,jhsl)
      call ga_put(lslab(k),ilsl,ilsl+ngrx-1,jlsl,jlsl+ngry-1,
     + grid(1,1),mgrx)
    1 continue
c
      return
      end
      subroutine pme_getgrd(node,nsl,lsldo,maxsl,grid)
c
      implicit none
c
#include "pme_common.fh"
#include "global.fh"
c
      integer node,nsl
      logical lsldo(np,mgz)
      integer maxsl(4,np,mgz)
      real*8 grid(mgrx,mgry,*)
c
      integer k,ilsl,ihsl,jlsl,jhsl
      integer minx,maxx,miny,maxy
c
      do 1 k=1,nsl
      if(lsldo(node+1,k)) then
      minx=maxsl(1,node+1,k)
      maxx=maxsl(2,node+1,k)
      miny=maxsl(3,node+1,k)
      maxy=maxsl(4,node+1,k)
c      minx=1
c      maxx=ngx
c      miny=1
c      maxy=ngy
      call ga_distribution(lslab(k),node,ilsl,ihsl,jlsl,jhsl)
      call ga_get(lslab(k),ilsl+minx-1,ilsl+maxx-1,jlsl+miny-1,
     + jlsl+maxy-1,grid(minx,miny,k),mgrx)
      endif
    1 continue
c
      return
      end
      subroutine pme_forces(fw,nwmloc,fs,isq,isgm,nsaloc)
c
      implicit none
c
#include "pme_common.fh"
#include "mafdecls.fh"
c
      integer nwmloc,nsaloc
      real*8 fw(mwm,3,mwa),fs(msa,3)
      integer isq(msa),isgm(msa)
c
      call pme_for(fw,nwmloc,fs,isq,isgm,nsaloc,int_mb(i_lnod),
     + int_mb(i_kfr),int_mb(i_kto))
c
      return
      end
      subroutine pme_for(fw,nwmloc,fs,isq,isgm,nsaloc,lnode,kfrom,kto)
c
      implicit none
c
#include "pme_common.fh"
#include "mafdecls.fh"
c
      integer nwmloc,nsaloc
      real*8 fw(mwm,3,mwa),fs(msa,3)
      integer isq(msa),isgm(msa)
      integer lnode(np),kfrom(np),kto(np)
c
      integer inode,node,ngznon
c
c     wait until potential grid is available
c
      call timer_start(39)
      call pme_wait(0,1)
      call timer_stop(39)
      call timer_start(40)
c
      do 1 inode=1,np
      node=me+inode
      if(node.gt.np) node=node-np
      node=node-1
c
c     node can be skipped if it did not contribute to the charge grid
c
      if(lnode(node+1).eq.0) goto 1
c
c     determine the number of planes on non-local node
c
      ngznon=0
      if(kto(node+1).gt.0) ngznon=kto(node+1)-kfrom(node+1)+1
      if(ngznon.eq.0) goto 1
c
      call timer_stop(40)
      call timer_start(41)
      call pme_getgrd(node,ngznon,log_mb(i_lsld),int_mb(i_slmax),
     + dbl_mb(i_grd))
      call timer_stop(41)
      call timer_start(40)
c
      call cf_fpme(fw,nwmloc,fs,isq,isgm,nsaloc,dbl_mb(i_grd),mgz,
     + int_mb(i_knod),int_mb(i_kfr),log_mb(i_lsld),node)
c
    1 continue
c
      call timer_stop(40)
c
      return
      end
      subroutine pme_sync(number,nodes)
c
      implicit none
c
      integer number
      integer nodes(number)
c
      external mxcombv1,pme_synop
      integer mxcombv1,pme_synop
c
      integer ires
      real*8 datum,buffer
c
      datum=0.0d0
      ires=mxcombv1(datum,pme_synop,8,1,number,nodes,123,buffer)
c
      return
      end
      integer function pme_synop(idat,jdat,iteml,nitems)
c
      integer idat,jdat,iteml,nitems
c
      idat=0
      jdat=0
      pme_synop=0
c
      return
      end
      subroutine pme_flag(iflag,icolor,inode)
c
      implicit none
c
#include "pme_common.fh"
c
      integer iflag,icolor,inode
      integer il,ih,jl,jh,node,icolr
c
c     set flag iflag on node inode to color icolor:
c
c     icolor>=0 set flag to icolor
c     icolor<0  accumulate iabs(icolor) to flag
c
c     inode>=0  set/accumulate flag on node inode
c     inode<0   set/accumulate flag on all nodes
c
      icolr=iabs(icolor)
c
      if(inode.lt.0) then
      do 1 node=0,np-1
      call ga_distribution(ga_flg,node,il,ih,jl,jh)
      if(icolor.ge.0) then
      call ga_put(ga_flg,il,ih,jl+iflag,jl+iflag,icolr,1)
      else
      call ga_acc(ga_flg,il,ih,jl+iflag,jl+iflag,icolr,1,1)
      endif
    1 continue
      else
c
c     set flag on single node
c
      call ga_distribution(ga_flg,inode,il,ih,jl,jh)
      if(icolor.ge.0) then
      call ga_put(ga_flg,il,ih,jl+iflag,jl+iflag,icolr,1)
      else
      call ga_acc(ga_flg,il,ih,jl+iflag,jl+iflag,icolr,1,1)
      endif
      endif
c
      return
      end
      subroutine pme_wait(iflag,icolor)
c
      implicit none
c
#include "pme_common.fh"
#include "msgids.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      real*8 timer_wall
      external timer_wall
c
      integer iflag,icolor
      integer il,ih,jl,jh,icolr
c
      call ga_distribution(ga_flg,me,il,ih,jl,jh)
    1 continue
      call ga_get(ga_flg,il,ih,jl+iflag,jl+iflag,icolr,1)
      if(icolr.ne.icolor) goto 1
c
      return
      end
