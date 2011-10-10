      logical function argos_cafe_hop(lseq,issgm,lda,rda,uda,ndx,pda,
     + isgan,ishop,xs,nsaloc,lsthop,timhop)
c
      implicit none
c
#include "argos_cafe_common.fh"
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
      argos_cafe_hop=.false.
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
      argos_cafe_hop=.true.
      endif
    7 continue
c
      return
      end
c $Id$
