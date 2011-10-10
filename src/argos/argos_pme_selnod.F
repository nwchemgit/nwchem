      subroutine argos_pme_selnod(value,target,nodep,knode,
     + kfrom,kto,map)
c
      implicit none
c
#include "argos_pme_common.fh"
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
     + call md_abort('Error in argos_pme_select',me)
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
c $Id$
