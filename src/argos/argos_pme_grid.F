      subroutine argos_pme_grid(ipset,xw,nwmloc,xs,isq,isgm,nsaloc,
     + nodep,knode,kfrom,kto,lnode,lsldo,maxsl,gridr,gridc,grid,
     + lgridh,gridt)
c
      implicit none
c
#include "argos_pme_common.fh"
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
      call argos_cafe_chgrid(lset,ipset,xw,nwmloc,xs,isq,isgm,nsaloc,
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
c     the following argos_pme_flag/argos_pme_wait/argos_pme_flag combination may be
c     replaced by a ga_sync() at loss of performance
c
      call argos_pme_flag(1,-1,-1)
      if(lpnode) call argos_pme_wait(1,np)
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
c $Id$
