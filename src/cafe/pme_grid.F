      subroutine pme_chgrid(ipset,xw,nwmloc,xs,isq,nsaloc,epme)
c $Id: pme_grid.F,v 1.9 1999-09-08 16:49:56 d3j191 Exp $
      implicit none
c
#include "pme_common.fh"
#include "mafdecls.fh"
c
      integer ipset,nwmloc,nsaloc
      real*8 xw(mwm,3,mwa),xs(msa,3),epme
      integer isq(msa)
c
      call timer_start(11)
      if(.not.lcorr(ipset)) then
      call cf_pmecor(ipset,xw,nwmloc,xs,isq,nsaloc,ecorr(ipset))
      lcorr(ipset)=.true.
      endif
c
      call pme_grid(ipset,xw,nwmloc,xs,isq,nsaloc,
     + int_mb(i_nodp),int_mb(i_knod),int_mb(i_kfr),int_mb(i_kto),
     + int_mb(i_lnod),log_mb(i_lsld),dbl_mb(i_grd1),dbl_mb(i_grd2))
c
      epme=ecorr(ipset)
      call timer_stop(11)
c
      if(.not.lpnode) return
c
      call timer_start(12)
c
      call pme_fft(ngx,ngy,ngz,-2,dbl_mb(i_grd2),dbl_mb(i_grd1),
     + int_mb(i_map1),int_mb(i_map2),int_mb(i_knod),imfft)
c
      call cf_epme(int_mb(i_knod),int_mb(i_kfr),int_mb(i_kto),
     + dbl_mb(i_bmod),dbl_mb(i_grd1),mgz,epme)
c
      epme=epme+ecorr(ipset)
c
      call pme_sync(nodpme,int_mb(i_nodp))
c
      call pme_fft(ngx,ngy,ngz,1,dbl_mb(i_grd1),dbl_mb(i_grd2),
     + int_mb(i_map1),int_mb(i_map2),int_mb(i_knod),imfft)
c
      call pme_putgrd(dbl_mb(i_grd2))
c
c     set flag on all nodes to indicate that potential
c     grid calculation has completed
c
      call pme_flag(1,-1)
c
      call timer_stop(12)
c
      return
      end
      subroutine pme_grid(ipset,xw,nwmloc,xs,isq,nsaloc,
     + nodep,knode,kfrom,kto,lnode,lsldo,gridr,gridc)
c
      implicit none
c
#include "pme_common.fh"
#include "global.fh"
#include "mafdecls.fh"
c
      integer ipset,nwmloc,nsaloc
      real*8 xw(mwm,3,mwa),xs(msa,3)
      integer isq(msa)
      logical lsldo(np,mgz)
      integer nodep(nodpme),knode(ngz),kfrom(np),kto(np),lnode(np)
      real*8 gridr(mgx,mgy,mgz),gridc(2,mgx,mgy,mgz)
c
      integer node,inode,ngznon
      integer i,j,k,minx,maxx,miny,maxy,minz,maxz
      integer ilsl,ihsl,jlsl,jhsl
      logical lset,lgrid
c
      lset=.true.
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
      call cf_chgrid(lset,ipset,xw,nwmloc,xs,isq,nsaloc,
     + gridr,ngznon,lgrid,knode,kfrom,lsldo,node,
     + minx,maxx,miny,maxy,minz,maxz)
c
      if(lgrid) then
      lnode(node+1)=1
c
      do 21 k=1,ngznon
      if(lsldo(node+1,k)) then
      call ga_distribution(lslab(k),node,ilsl,ihsl,jlsl,jhsl)
      call ga_acc(lslab(k),ilsl+minx-1,ilsl+maxx-1,jlsl+miny-1,
     + jlsl+maxy-1,gridr(minx,miny,k),mgx,one)
      endif
   21 continue
c
      endif
c
    1 continue
c
c     pmenodes should wait until all nodes have finished contributions to the grid
c
      call pme_sync(nodpme,nodep)
c
c     retrieve local slab from global array space
c
      if(ngzloc.gt.0) then
      do 23 k=ngzloc,1,-1
      call ga_distribution(lslab(k),me,ilsl,ihsl,jlsl,jhsl)
      call ga_get(lslab(k),ilsl,ilsl+ngx-1,jlsl,jlsl+ngy-1,
     + gridr(1,1,k),mgx)
      do 24 j=ngy,1,-1
      do 25 i=ngx,1,-1
      gridc(1,i,j,k)=gridr(i,j,k)
      gridc(2,i,j,k)=zero
   25 continue
   24 continue
   23 continue
      endif
c
      return
      end
      subroutine pme_putgrd(grid)
c
      implicit none
c
#include "pme_common.fh"
#include "global.fh"
c
      real*8 grid(2,mgx,mgy,*)
c
      integer k,ilsl,ihsl,jlsl,jhsl
c
      do 1 k=1,ngzloc
      call ga_distribution(lslab(k),me,ilsl,ihsl,jlsl,jhsl)
      call ga_put(lslab(k),ilsl,ilsl+2*ngx-1,jlsl,jlsl+ngy-1,
     + grid(1,1,1,k),2*mgx)
    1 continue
c
      return
      end
      subroutine pme_getgrd(node,nsl,grid)
c
      implicit none
c
#include "pme_common.fh"
#include "global.fh"
c
      integer node,nsl
      real*8 grid(2,mgx,mgy,*)
c
      integer k,ilsl,ihsl,jlsl,jhsl
c
      do 1 k=1,nsl
      call ga_distribution(lslab(k),node,ilsl,ihsl,jlsl,jhsl)
      call ga_get(lslab(k),ilsl,ilsl+2*ngx-1,jlsl,jlsl+ngy-1,
     + grid(1,1,1,k),2*mgx)
    1 continue
c
      return
      end
      subroutine pme_forces(fw,nwmloc,fs,isq,nsaloc)
c
      implicit none
c
#include "pme_common.fh"
#include "mafdecls.fh"
c
      integer nwmloc,nsaloc
      real*8 fw(mwm,3,mwa),fs(msa,3)
      integer isq(msa)
c
      call pme_for(fw,nwmloc,fs,isq,nsaloc,int_mb(i_lnod),
     + int_mb(i_kfr),int_mb(i_kto))
c
      return
      end
      subroutine pme_for(fw,nwmloc,fs,isq,nsaloc,lnode,kfrom,kto)
c
      implicit none
c
#include "pme_common.fh"
#include "mafdecls.fh"
c
      integer nwmloc,nsaloc
      real*8 fw(mwm,3,mwa),fs(msa,3)
      integer isq(msa)
      integer lnode(np),kfrom(np),kto(np)
c
      integer inode,node,ngznon
c
c     wait until potential grid is available
c
      call timer_start(15)
      call pme_wait(1)
      call timer_stop(15)
      call timer_start(16)
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
      call pme_getgrd(node,ngznon,dbl_mb(i_grd1))
c
      call cf_fpme(fw,nwmloc,fs,isq,nsaloc,dbl_mb(i_grd1),mgz,
     + int_mb(i_knod),int_mb(i_kfr),log_mb(i_lsld),node)
c
    1 continue
c
      call timer_stop(16)
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
      subroutine pme_flag(icolor,inode)
c
      implicit none
c
#include "pme_common.fh"
c
      integer icolor,inode
      integer il,ih,jl,jh,node,icolr
c
c     set flag on node:
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
      call ga_put(ga_flg,il,ih,jl,jh,icolr,1)
      else
      call ga_acc(ga_flg,il,ih,jl,jh,icolr,1,1)
      endif
    1 continue
      else
c
c     set flag on single node
c
      call ga_distribution(ga_flg,inode,il,ih,jl,jh)
      if(icolor.ge.0) then
      call ga_put(ga_flg,il,ih,jl,jh,icolr,1)
      else
      call ga_acc(ga_flg,il,ih,jl,jh,icolr,1,1)
      endif
      endif
c
      return
      end
      subroutine pme_wait(icolor)
c
      implicit none
c
#include "pme_common.fh"
c
      integer icolor
      integer il,ih,jl,jh,icolr
c
c     wait until flag color on this node is icolor
c
      call ga_distribution(ga_flg,me,il,ih,jl,jh)
    1 continue
      call ga_get(ga_flg,il,ih,jl,jh,icolr,1)
      if(icolr.ne.icolor) goto 1
c
      return
      end

