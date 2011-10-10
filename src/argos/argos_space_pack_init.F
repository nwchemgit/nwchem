      subroutine argos_space_pack_init(isl,numsa,iwl,numwm)
c
      implicit none
c
#include "argos_space_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "bitops.fh"
c
      integer isl(msa,mis2),numsa,iwl(mwm,miw2),numwm
c
      integer i,j,k,intsiz
c
      intsiz=8*ma_sizeof(mt_int,1,mt_byte)-1
c
      do 1 i=1,mis2
      misl(i)=0
    1 continue
c
      do 2 j=1,mis2
      do 3 i=1,numsa
      misl(j)=max(misl(j),isl(i,j))
    3 continue
    2 continue
      misl(lsdyn)=max(15,misl(lsdyn))
      misl(lsbox)=nbtot
      misl(lsnod)=np
c
      call ga_igop(msp_24,misl,mis2,'max')
c
      do 4 i=1,mis2
      k=misl(i)
      do 5 j=1,intsiz
      k=rshift(k,1)
      if(k.eq.0) then
      nbits(i)=j+1
      goto 4
      endif
    5 continue
      nbits(i)=intsiz
    4 continue
c
      npack=1
      ipack(1)=1
      k=nbits(1)
      do 6 i=2,mis2
      if(k+nbits(i).lt.intsiz) then
      ipack(npack)=ipack(npack)+1
      k=k+nbits(i)
      else
      npack=npack+1
      ipack(npack)=1
      k=nbits(i)
      endif
    6 continue
c
      do 7 i=1,mis2
      misl(i)=2**nbits(i)-1
    7 continue
c
      do 8 i=1,miw2
      miwl(i)=0
    8 continue
c
      do 9 j=1,miw2
      do 10 i=1,numwm
      miwl(j)=max(miwl(j),iwl(i,j))
   10 continue
    9 continue
      miwl(lwdyn)=15
      miwl(lwbox)=nbtot
      miwl(lwnod)=np
c
      call ga_igop(msp_25,miwl,miw2,'max')
c
      do 11 i=1,miw2
      k=miwl(i)
      do 12 j=1,intsiz
      k=rshift(k,1)
      if(k.eq.0) then
      nbitw(i)=j+1
      goto 11
      endif
   12 continue
      nbitw(i)=intsiz
   11 continue
c
      npackw=1
      ipackw(1)=1
      k=nbitw(1)
      do 13 i=2,miw2
      if(k+nbitw(i).lt.intsiz) then
      ipackw(npackw)=ipackw(npackw)+1
      k=k+nbitw(i)
      else
      npackw=npackw+1
      ipackw(npackw)=1
      k=nbitw(i)
      endif
   13 continue
c
      do 14 i=1,miw2
      miwl(i)=2**nbitw(i)-1
   14 continue
c
      if(me.eq.0) then
      write(lfnout,1001)
 1001 format(/,' DATA PACKING',/)
      write(lfnout,1002) ' Packing solute in  ',npack,
     + ' integers : ',(ipack(i),i=1,npack)
 1002 format(a,i2,a,15i4)
      write(lfnout,1003) ' significant bits : ',
     + (nbits(i),i=1,mis2)
 1003 format(14x,a,15i4)
      write(lfnout,1002) ' Packing solvent in ',npackw,
     + ' integers : ',(ipackw(i),i=1,npackw)
      write(lfnout,1003) ' significant bits : ',
     + (nbitw(i),i=1,miw2)
      endif
c
      if(npack.eq.mis2) npack=0
      if(npackw.eq.miw2) npackw=0
c      npackw=0
c
      if(npack.gt.0) then
      if(.not.ma_push_get(mt_int,msa*npack,'lpack',l_pack,i_pack))
     + call md_abort('Failed to allocate memory for pack',0)
      endif
      if(npackw.gt.0) then
      if(.not.ma_push_get(mt_int,mwm*npackw,'lpackw',l_packw,i_packw))
     + call md_abort('Failed to allocate memory for packw',0)
      endif
c
      call argos_space_pack_ga(numsa,isl,int_mb(i_pack),
     + numwm,iwl,int_mb(i_packw))
c
      return
      end
c $Id$
