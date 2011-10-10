      subroutine argos_space_lbbl_ind(lbbl,ipl,jpl)
c
      implicit none
c
#include "argos_space_common.fh"
#include "global.fh"
c
      integer lbbl(mbbl,mbb2),ipl(mbox,mip2),jpl(mbox,mip2)
c
      integer ibbl,jproc,jcell,icell,jcproc,jccell,iccell
      integer i,j,indx,jndx,il,ih,jl,jh,indexw,indexs,nndexw,nndexs
      character*80 string
c
      jcproc=-1
      jccell=-1
      iccell=-1
      indexw=nwmloc+1
      indexs=nsaloc+1
      nndexw=indexw
      nndexs=indexs
c
c      if(iand(idebug,16).eq.16) then
c      write(lfndbg,3000) ((lbbl(ibbl,i),i=1,3),ibbl=1,nbbl)
c 3000 format(3i5)
c      endif
      do 1 ibbl=1,nbbl
c
      jproc=lbbl(ibbl,1)
      jcell=lbbl(ibbl,2)
      icell=lbbl(ibbl,3)
      lbbl(ibbl, 5)=0
      lbbl(ibbl, 6)=0
      if(iand(idebug,16).eq.16) then
      write(lfndbg,222) ibbl,nbbl,jproc,jcell,icell,jccell,iccell
  222 format(7i5)
      endif
      do 234 j=1,3
  234 continue
c
c     get indices for icell on the local node
c
      if(iccell.ne.icell) then
      indx=0
      do 2 i=1,ipl(1,1)
      if(ipl(30+i,1).eq.icell) indx=30+i
    2 continue
      if(indx.eq.0) then
      do 232 j=1,3
  232 continue
      write(*,1000) me,(ipl(30+i,1),i=1,ipl(1,1))
 1000 format(/,'Cells on proc',i5,' :',t20,20i5,/,(t20,20i5))
      write(string,1001) icell,me
 1001 format('SP0001: Could not find icell',i5,' on proc',i5)
      call md_abort(string,0)
      endif
      iccell=icell
      endif
      lbbl(ibbl, 7)=ipl(indx,2)
      lbbl(ibbl, 8)=ipl(indx,3)
      lbbl(ibbl, 9)=ipl(indx,4)
      lbbl(ibbl,10)=ipl(indx,5)
c
c     get indices for jcell if on the local node
c
      if(jproc.eq.me) then
c
      if(jccell.ne.jcell) then
      jndx=0
      do 3 i=1,ipl(1,1)
      if(ipl(30+i,1).eq.jcell) jndx=30+i
    3 continue
      if(jndx.eq.0) then
      write(*,1000) me,(ipl(30+i,1),i=1,ipl(1,1))
      write(string,1002) jcell,me
 1002 format('SP0001: Could not find jcell',i5,' on proc',i5)
      call md_abort(string,0)
      endif
      jccell=jcell
      endif
      lbbl(ibbl,11)=ipl(jndx,2)
      lbbl(ibbl,12)=ipl(jndx,3)
      lbbl(ibbl,13)=ipl(jndx,4)
      lbbl(ibbl,14)=ipl(jndx,5)
c
c     get indices for jcell if on the remote node
c
      else
c
      if(jproc.ne.jcproc) then
      call ga_distribution(ga_ip,jproc,il,ih,jl,jh)
      call ga_get(ga_ip,il,ih,jl,jh,jpl,mbox)
      jcproc=jproc
      endif
c
      if(jcell.ne.jccell) then
      jndx=0
      do 4 i=1,jpl(1,1)
      if(jpl(30+i,1).eq.jcell) jndx=30+i
    4 continue
      if(jndx.eq.0) then
      if(iand(idebug,16).eq.16) then
      write(lfndbg,1003) jcell,jproc
 1003 format('Could not find cell',i5,' on proc',i5,
     + ' in argos_space_lbbl')
      write(lfndbg,1004) jproc,(jpl(30+i,1),i=1,jpl(1,1))
 1004 format('Cell list obtained from proc',i5,' is',/,(20i5))
      call util_flush(lfndbg)
      endif
c
      write(*,1005) jproc,(jpl(30+i,1),i=1,jpl(1,1))
 1005 format(/,'Cells on proc',i5,':',t20,20i5,/,(t20,20i5))
      write(string,1006) jcell,jproc
 1006 format('SP0002: Could not find remote jcell',i5,' on proc',i5)
      call md_abort(string,me)
      endif
      jccell=jcell
c
      if(nbget.ne.0) then
      indexw=nndexw
      indexs=nndexs
      if(jpl(jndx,2).gt.0) nndexw=indexw+jpl(jndx,3)-jpl(jndx,2)+1
      if(jpl(jndx,4).gt.0) nndexs=indexs+jpl(jndx,5)-jpl(jndx,4)+1
      endif
c
      endif
c
      lbbl(ibbl, 5)=indexw
      lbbl(ibbl, 6)=indexs
      lbbl(ibbl,11)=jpl(jndx,2)
      lbbl(ibbl,12)=jpl(jndx,3)
      lbbl(ibbl,13)=jpl(jndx,4)
      lbbl(ibbl,14)=jpl(jndx,5)
c
      endif
c
    1 continue
c
      if(nndexw.gt.mwm) then
      write(string,1007) mwm,nndexw
 1007 format('Increase dimension mwm from ',i7,' to ',i7)
      call md_abort(string,me)
      endif
      if(nndexs.gt.msa) then
      write(string,1008) msa,nndexs
 1008 format('Increase dimension msa from ',i7,' to ',i7)
      call md_abort(string,me)
      endif
c
      return
      end
c $Id$
