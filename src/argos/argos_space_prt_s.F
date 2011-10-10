      subroutine argos_space_prt_s(filtop,lfntop,
     + ipl,isl,islp,xs,fs,nb,ibs,bxs,bfs,
     + npener,esa)
c
      implicit none
c
#include "argos_space_common.fh"
#include "msgids.fh"
#include "global.fh"
c
      integer lfntop,nb,npener
      character*70 filtop
      integer ipl(mbox,mip2),isl(msa,mis2),islp(msa,npack),ibs(nb)
      real*8 xs(msa,3),fs(msa,3)
      real*8 bxs(nb,3),bfs(nb,3)
      real*8 esa(nsa,2)
c
      integer i,j,k,number,icyc,ncyc,node,nums,nsan,idyn,ism,iss
      integer ilp,ihp,jlp,jhp,ili,ihi,jli,jhi,ils,ihs,jls,jhs
      integer naw,nbw,nhw,ndw,now,ntw,nnw,nat,nqt,idum,npars
      character*1 cdum
      character*16 cat
c
      call ga_distribution(ga_ip,me,ilp,ihp,jlp,jhp)
      call ga_get(ga_ip,ilp,ihp,jlp,jhp,ipl,mbox)
      nsan=ipl(2,2)
c
      if(nsan.gt.0) then
      call ga_distribution(ga_s,me,ils,ihs,jls,jhs)
      call ga_put(ga_s,ils,ils+nsan-1,jls+6,jls+8,fs,msa)
      endif
c
      call ga_sync()
c
      if(me.ne.0) return
c
      if(npener.le.0) then
      write(lfnout,1007)
 1007 format(//,' Solute coordinates and forces',//,
     + '   mol  segment     atom  dt     x       y       z ',
     + '             fx          fy          fz',/)
      else
      write(lfnout,1008)
 1008 format(//,' Solute coordinates, forces and energies',//,
     + '   mol  segment     atom  dt     x       y       z ',
     + '             fx          fy          fz        enb',
     + '          eb          ep',/)
      endif
c
      open(unit=lfntop,file=filtop(1:index(filtop,' ')-1),
     + status='old',form='formatted')
      rewind(lfntop)
      do 101 i=1,4
      read(lfntop,2001) cdum
 2001 format(a1)
  101 continue
      read(lfntop,2002) npars
      read(lfntop,2002) nat
      read(lfntop,2002) nqt
      read(lfntop,2002) nseq
 2002 format(i5)
      read(lfntop,2001) cdum
      do 1102 k=1,npars
      do 102 i=1,nat
      read(lfntop,2001) cdum
      do 103 j=i,nat
      read(lfntop,2001) cdum
  103 continue
  102 continue
 1102 continue
      do 104 i=1,nqt*npars
      read(lfntop,2001) cdum
  104 continue
      do 4104 i=1,nseq
      read(lfntop,2001) cdum
 4104 continue
      read(lfntop,2003) naw,nbw,nhw,ndw,now,ntw,nnw
 2003 format(5i7,2i10)
      read(lfntop,2001) cdum
      do 105 i=1,naw
      read(lfntop,2001) cdum
  105 continue
      do 106 i=1,nbw*(npars+1)
      read(lfntop,2001) cdum
  106 continue
      do 107 i=1,nhw*(npars+1)
      read(lfntop,2001) cdum
  107 continue
      do 108 i=1,ndw*(npars+1)
      read(lfntop,2001) cdum
  108 continue
      do 109 i=1,now*(npars+1)
      read(lfntop,2001) cdum
  109 continue
      if(ntw.gt.0) then
      read(lfntop,2004) (idum,i=1,ntw)
      read(lfntop,2004) (idum,i=1,ntw)
 2004 format(11i7)
      endif
      if(nnw.gt.0) then
      read(lfntop,2005) (idum,i=1,nnw)
      read(lfntop,2005) (idum,i=1,nnw)
 2005 format(11i7)
      endif
      read(lfntop,2001) cdum
      do 204 i=1,npars
      read(lfntop,2001) cdum
  204 continue
c
      if(nsa.gt.0) then
      number=0
      ncyc=nsa/nb+1
      nums=nb
      do 6 icyc=1,ncyc
      if(nsa-number.lt.nums) nums=nsa-number
      do 7 node=np-1,0,-1
      call ga_distribution(ga_ip,node,ilp,ihp,jlp,jhp)
      call ga_get(ga_ip,ilp,ihp,jlp,jhp,ipl,mbox)
      nsan=ipl(2,2)
      if(nsan.gt.0) then
      call ga_distribution(ga_is,node,ili,ihi,jli,jhi)
      if(npack.eq.0) then
      call ga_get(ga_is,ili,ili+nsan-1,jli,jli+lsdyn-1,isl,msa)
      else
      call ga_get(ga_is,ili,ili+nsan-1,jli,jli+npack-1,islp,msa)
      call argos_space_unpack(nsan,isl,islp)
      endif
      call ga_distribution(ga_s,node,ils,ihs,jls,jhs)
      call ga_get(ga_s,ils,ils+nsan-1,jls,jls+2,xs,msa)
      call ga_get(ga_s,ils,ils+nsan-1,jls+6,jls+8,fs,msa)
      do 8 i=1,nsan
      j=isl(i,lsgan)-number
      if(j.gt.0.and.j.le.nums) then
      bxs(j,1)=xs(i,1)
      bxs(j,2)=xs(i,2)
      bxs(j,3)=xs(i,3)
      bfs(j,1)=fs(i,1)
      bfs(j,2)=fs(i,2)
      bfs(j,3)=fs(i,3)
      ibs(j)=isl(i,lsdyn)
      endif
    8 continue
      endif
    7 continue
      do 9 i=1,nums
      read(lfntop,2009) cat,ism,iss
 2009 format(a16,3x,2i7)
c 2009 format(a16,3x,2i7)
      idyn=iand(ibs(i),3)
      if(npener.le.0) then
      write(lfnout,1009) ism,iss,cat,idyn,
     + (bxs(i,j),j=1,3),(bfs(i,j),j=1,3)
 1009 format(2i5,':',a16,i1,1x,3f8.4,3x,3f12.3)
      else
      write(lfnout,1010) ism,iss,cat,idyn,
     + (bxs(i,j),j=1,3),(bfs(i,j),j=1,3),esa(i,1),esa(i,2),
     + esa(i,1)+esa(i,2)
 1010 format(2i5,':',a16,i1,1x,3f8.4,3x,3f12.3,3(2x,1pe10.3))
      endif
    9 continue
    6 continue
      endif
c
      close(unit=lfntop,status='keep')
c
      return
      end
c $Id$
