      subroutine argos_space_initip(ibownr,ipl)
c
      implicit none
c
#include "argos_space_common.fh"
#include "global.fh"
c
      integer ibownr(maxbox,3),ipl(mbox,mip2)
      integer i,j,ibx,iby,ibz,nrbox,nrnod,nbiown
      integer il,ih,jl,jh
c
      do 1 j=1,mip2
      do 2 i=1,30
      ipl(i,j)=0
    2 continue
    1 continue
c
      if(iand(idebug,1).eq.1) then
      write(lfndbg,8003) nbx,nby,nbz,maxbox
 8003 format('ibownr in argos_space_initip',4i5)
      write(lfndbg,8004) (ibownr(ibx,1),ibx=1,nbx)
 8004 format('ibownr x',/,(20i5))
      write(lfndbg,8005) (ibownr(iby,2),iby=1,nby)
 8005 format('ibownr y',/,(20i5))
      write(lfndbg,8006) (ibownr(ibz,3),ibz=1,nbz)
 8006 format('ibownr z',/,(20i5))
      call util_flush(lfndbg)
      if(iand(idebug,2).eq.2) write(lfndbg,8002)
 8002 format('boxlist')
      endif
      nbiown=0
      do 3 ibx=1,nbx
      do 4 iby=1,nby
      do 5 ibz=1,nbz
      nrbox=((ibz-1)*nby+iby-1)*nbx+ibx-1
      nrnod=(ibownr(ibz,3)*npy+ibownr(iby,2))*npx+ibownr(ibx,1)
      if(iand(idebug,2).eq.2) then
      write(lfndbg,8001) ibx,iby,ibz,nrbox,nrnod
 8001 format(5i5)
      call util_flush(lfndbg)
      endif
      if(me.eq.nrnod) then
      nbiown=nbiown+1
      ipl(30+nbiown,1)=nrbox
      ipl(30+nbiown,2)=0
      ipl(30+nbiown,3)=0
      ipl(30+nbiown,4)=0
      ipl(30+nbiown,5)=0
      endif
    5 continue
    4 continue
    3 continue
      ipl(1,1)=nbiown
c
      call ga_distribution(ga_ip,me,il,ih,jl,jh)
      call ga_put(ga_ip,il,ih,jl,jh,ipl,mbox)
c
      if(iand(idebug,2).eq.2) then
      write(lfndbg,8000) (i,ipl(30+i,1),i=1,ipl(1,1))
 8000 format('ipl',/,(2i5))
      call util_flush(lfndbg)
      endif
c
      return
      end
c $Id$
