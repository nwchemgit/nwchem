      logical function cf_hopping(lpbc,lpbcs,stimei,
     + isl,issgm,isgan,isq3,ishop,xs,nsaloc)
c
c $Id$
c
      implicit none
c
#include "cf_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      logical cf_hop
      external cf_hop
c
      logical lpbc,lpbcs
      integer isl(msa,mis2)
      integer issgm(msa),nsaloc,isgan(msa),ishop(msa),isq3(msa)
      real*8 xs(msa,3)
      real*8 stimei
c
      integer i_itmp,l_itmp,i_dtmp,l_dtmp
c
      stime=stimei
      cf_hopping=.false.
c
      if(lpair.and.nhop.gt.0) then
      lhop=.true.
c
      if(.not.ma_push_get(mt_int,np,'itmp',l_itmp,i_itmp))
     + call md_abort('Failed to allocate itmp',me)
c
      call cf_hoplist(issgm,nsaloc,int_mb(i_itmp),int_mb(i_lda),
     + dbl_mb(i_rda))
c
      if(.not.ma_pop_stack(l_itmp))
     + call md_abort('Failed to deallocate itmp',me)
c
      cf_hopping=.true.
      elseif(lhop) then
c
c     evaluate the E12
c
      call cf_uhop(lpbc,lpbcs,int_mb(i_lseq),isl,issgm,isgan,isq3,
     + xs,nsaloc,
     + dbl_mb(i_chg),int_mb(i_lda),dbl_mb(i_rda),dbl_mb(i_uda))
c
      if(.not.ma_push_get(mt_dbl,nldat,'dtmp',l_dtmp,i_dtmp))
     + call md_abort('Failed to allocate dtmp',me)
      if(.not.ma_push_get(mt_int,nldat,'itmp',l_itmp,i_itmp))
     + call md_abort('Failed to allocate itmp',me)
c
      cf_hopping=cf_hop(int_mb(i_lseq),issgm,
     + int_mb(i_lda),dbl_mb(i_rda),dbl_mb(i_uda),
     + int_mb(i_itmp),dbl_mb(i_dtmp),isgan,ishop,xs,nsaloc,
     + int_mb(i_lsthop),dbl_mb(i_timhop))
c
      if(.not.ma_pop_stack(l_itmp))
     + call md_abort('Failed to deallocate itmp',me)
      if(.not.ma_pop_stack(l_dtmp))
     + call md_abort('Failed to deallocate dtmp',me)
c
      endif
c
      return
      end
      subroutine cf_hoplist(issgm,nsaloc,nhopl,lda,rda)
c
      implicit none
c
#include "cf_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "util.fh"
c
      integer issgm(msa),nsaloc,nhopl(np)
      integer lda(16,*)
      real*8 rda(11,*)
c
      integer i,j,ioff
c
      do 1 i=1,np
      nhopl(i)=0
    1 continue
      nhopl(me+1)=nlda
c      write(*,'(i5,a,i5)') me,' nlda= ',nlda
c
      call ga_igop(mcf_76,nhopl,np,'+')
c
      do 2 i=2,np
      nhopl(i)=nhopl(i)+nhopl(i-1)
    2 continue
c      write(*,'(i5,a,10i5)') me,' nhopl ',(nhopl(i),i=1,np)
c
c
      ioff=0
      if(me.gt.0) then
      ioff=nhopl(me)
      do 3 i=1,nlda
      do 4 j=1,16
      lda(j,i+ioff)=lda(j,i)
    4 continue
      do 5 j=1,8
      rda(j,i+ioff)=rda(j,i)
    5 continue
    3 continue
      do 6 i=1,nhopl(me)
      do 7 j=1,16
      lda(j,i)=0
    7 continue
      do 8 j=1,8
      rda(j,i)=0.0d0
    8 continue
    6 continue
      endif
      if(me.lt.np-1) then
      do 9 i=nhopl(me+1)+1,nhopl(np)
      do 10 j=1,16
      lda(j,i)=0
   10 continue
      do 11 j=1,8
      rda(j,i)=0.0d0
   11 continue
    9 continue
      endif
c
      nldat=nhopl(np)
      call ga_igop(mcf_77,lda,16*nldat,'+')
c
c      do 322 i=1,nldat
c      write(*,'(18i5)') me,i,(lda(j,i),j=1,16)
c  322 continue
c
      do 201 i=1,nldat
      do 202 j=1,nsaloc
      if(lda(3,i).eq.issgm(j)) then
      if(lda(12,i).eq.0) lda(12,i)=j
      lda(13,i)=j
      endif
      if(lda(8,i).eq.issgm(j)) then
      if(lda(15,i).eq.0) lda(15,i)=j
      lda(16,i)=j
      endif
  202 continue
  201 continue
c
      if(me.gt.0) then
      do 221 j=1,nldat
      do 222 i=1,11
      lda(i,j)=0
  222 continue
      lda(14,j)=0
  221 continue
      endif
c
      call ga_igop(mcf_77,lda,16*nldat,'+')
      call ga_dgop(mcf_78,rda,11*nldat,'+')
c
c      do 2222 i=1,nlda
c      write(*,'(18i5)') me,i,(lda(j,i),j=1,16)
c 2222 continue
c
      nlda=0
      do 12 i=1,nldat
      do 13 j=1,nlda
      if(lda(1,i).eq.lda(1,j).and.lda(6,i).eq.lda(6,j)) goto 12
   13 continue
      nlda=nlda+1
      do 14 j=1,16
      lda(j,nlda)=lda(j,i)
   14 continue
      do 15 j=1,8
      rda(j,nlda)=rda(j,i)
   15 continue
   12 continue
      nldat=nlda
c
c
      if(me.eq.0) then
      if(util_print('qhop',print_high)) then
      write(lfnhop,110) nldat
  110 format(/,'Number of donor-acceptor pairs is ',i5,/,
     + '---------- Donor --------- --------- Acceptor ----',
     + ' --- Donor ---  ---Acceptor---',/,
     + ' glob  loc  sgm prot  off glob  loc  sgm prot  off',
     + '  own from   to  own from   to',/)
      do 212 i=1,nldat
      write(lfnhop,111) (lda(j,i),j=1,16),(rda(j,i),j=1,3)
  111 format(16i5,3f12.6)
  212 continue
      endif
      endif
c
      return
      end
      logical function cf_hop(lseq,issgm,lda,rda,uda,ndx,pda,
     + isgan,ishop,xs,nsaloc,lsthop,timhop)
c
      implicit none
c
#include "cf_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "util.fh"
#include "bitops_decls.fh"
#include "bitops_funcs.fh"
c
      integer lseq(mseq),issgm(msa)
      integer lda(16,*),nsaloc,isgan(msa),ishop(msa)
      real*8 rda(11,*),uda(4,*),pda(*),xs(msa,3),x(3)
      integer ndx(*)
      integer lsthop(2,*)
      real*8 timhop(*)
c
      integer i,j,k,ndxrel
      real*8 dprob,drand,e12,e120
c
      call ga_dgop(mcf_79,uda,4*nldat,'+')
c
      if(me.eq.0) then
      if(util_print('qhop',print_high)) then
      write(lfnhop,110) stime 
  110 format(/,' Time ',f12.6,/,
     + ' QHOP  ----donor---   --acceptor-- ',
     + '     pre-hop energies        post-hop energies  ',/,
     + ' numb atom  res prot atom  res prot',
     + '    donor      acceptor     donor      acceptor ',
     + '     e12        e120        distance     angle   probability',/)
      endif
      do 1 i=1,nldat
c      write(*,'(i5,8f12.6)') i,(rda(j,i),j=1,8)
      call qhop_prob(lda(1,i),lda(3,i),lseq(lda(3,i)),
     + lda(6,i),lda(8,i),lseq(lda(8,i)),uda(1,i),rda(4,i),pda(i),e120)
c      write(*,'(i5,8f12.6)') i,(rda(j,i),j=1,8)
c      pda(i)=0.0d0
c
      if(util_print('qhop',print_high)) then
      e12=uda(3,i)+uda(4,i)-uda(1,i)-uda(2,i)
      write(lfnhop,111) i,lda(1,i),lda(3,i),lseq(lda(3,i)),
     + lda(6,i),lda(8,i),lseq(lda(8,i)),
     + (uda(j,i),j=1,4),e12,e120*4.184,rda(4,i),rda(5,i),pda(i)
  111 format(7i5,8f12.6,f12.9)
      endif
c
    1 continue
c
      do 2 i=1,nldat
      ndx(i)=i
    2 continue
c
      do 3 i=1,nldat-1
      do 4 j=i+1,nldat
      if(pda(ndx(i)).lt.pda(ndx(j))) then
      k=ndx(i)
      ndx(i)=ndx(j)
      ndx(j)=k
      endif
    4 continue
    3 continue
c
      do 5 i=1,nldat
      drand=util_random(0)
      dprob=pda(ndx(i))
c      write(*,'(a,2f12.6)') 'RAND/PROB ',drand,dprob
      if(drand.lt.dprob) then
      write(lfnhop,'(/,a,i5,a,i5,a,f12.6/)')
     + ' PROTON HOP FROM ',lda(3,ndx(i)),' TO ',lda(8,ndx(i)),
     + ' AT TIME ',stime
      pda(ndx(i))=1.0d0
      do 6 j=i+1,nldat
      if(lda(3,ndx(i)).eq.lda(3,ndx(j))) pda(ndx(j))=0.0d0
      if(lda(8,ndx(i)).eq.lda(8,ndx(j))) pda(ndx(j))=0.0d0
      if(lda(3,ndx(i)).eq.lda(8,ndx(j))) pda(ndx(j))=0.0d0
      if(lda(8,ndx(i)).eq.lda(3,ndx(j))) pda(ndx(j))=0.0d0
    6 continue
      else
      pda(ndx(i))=0.0d0
      endif
    5 continue
c
      endif
c
      call ga_brdcst(mcf_80,pda,nldat*ma_sizeof(mt_dbl,1,mt_byte),0)
c
      cf_hop=.false.
      do 7 i=1,nldat
      if(pda(i).gt.0.5d0) then
      nhops=nhops+1
      lsthop(1,nhops)=lda(3,i)
      lsthop(2,nhops)=lda(8,i)
      timhop(nhops)=stime
      lseq(lda(3,i))=lda(4,i)
      lseq(lda(8,i))=lda(9,i)
      do 8 j=1,nsaloc
      if(lda(1,i).eq.isgan(j)) then
      call qhop_dsite(isgan(j),issgm(j),ndxrel)
c      xs(j,1)=0.9d0*xs(j+ndxrel,1)+0.1d0*xs(j,1)
c      xs(j,2)=0.9d0*xs(j+ndxrel,2)+0.1d0*xs(j,2)
c      xs(j,3)=0.9d0*xs(j+ndxrel,3)+0.1d0*xs(j,3)
      ishop(j)=ior(ishop(j),1)
      if(lda(5,i).ne.0) then
      x(1)=xs(j,1)
      x(2)=xs(j,2)
      x(3)=xs(j,3)
      xs(j,1)=xs(j+lda(5,i),1)
      xs(j,2)=xs(j+lda(5,i),2)
      xs(j,3)=xs(j+lda(5,i),3)
      xs(j+lda(5,i),1)=x(1)
      xs(j+lda(5,i),2)=x(2)
      xs(j+lda(5,i),3)=x(3)
      k=ishop(j)
      ishop(j)=ishop(j+lda(5,i))
      ishop(j+lda(5,i))=k
      endif
      endif
      if(lda(6,i).eq.isgan(j)) then
      xs(j,1)=rda(1,i)
      xs(j,2)=rda(2,i)
      xs(j,3)=rda(3,i)
      ishop(j)=ieor(ishop(j),1)
      if(lda(10,i).ne.0) then
      x(1)=xs(j,1)
      x(2)=xs(j,2)
      x(3)=xs(j,3)
      xs(j,1)=xs(j+lda(10,i),1)
      xs(j,2)=xs(j+lda(10,i),2)
      xs(j,3)=xs(j+lda(10,i),3)
      xs(j+lda(10,i),1)=x(1)
      xs(j+lda(10,i),2)=x(2)
      xs(j+lda(10,i),3)=x(3)
      k=ishop(j)
      ishop(j)=ishop(j+lda(10,i))
      ishop(j+lda(10,i))=k
      endif
      endif
    8 continue
      cf_hop=.true.
      endif
    7 continue
c
      return
      end
