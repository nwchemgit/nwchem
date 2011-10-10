      subroutine argos_cafe_rdpmf(lfnout,lfntop,ipmf,jpmf,rpmf)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      integer lfnout,lfntop
      integer ipmf(numpmf,8),jpmf(numpmf,4,npmfa)
      real*8 rpmf(3,6,numpmf)
c
      integer i,j,k,m,lenint,lendbl,lenlog,idpmf
c
      lenint=ma_sizeof(mt_int,1,mt_byte)
      lendbl=ma_sizeof(mt_dbl,1,mt_byte)
      lenlog=ma_sizeof(mt_log,1,mt_byte)
c
      nbias=0
      npmfi=0
      idpmf=1
      if(me.eq.0) then
      write(lfnout,1100)
 1100 format(//,' Potentials of mean force',//,
     + '  pmf  terms',/)
      do 1 i=1,numpmf
      read(lfntop,1000) ipmf(i,1),ipmf(i,6),ipmf(i,7),
     + (rpmf(j,1,i),rpmf(j,3,i),j=1,2),
     + (ipmf(i,j),j=2,5),ipmf(i,8)
 1000 format(3i5,2f12.6,2e12.5,5i5)
      if(idpmf.ne.ipmf(i,8)) then
      write(lfnout,1200) idpmf,npmfi
 1200 format(2i5)
      npmfi=1
      idpmf=ipmf(i,8)
      else
      npmfi=npmfi+1
      endif
      rpmf(2,1,i)=facpmf*rpmf(2,1,i)
      rpmf(2,3,i)=facpmf*rpmf(2,3,i)
      if(ipmf(i,8).eq.0) ipmf(i,8)=1
      if(ipmf(i,6).ne.0) nbias=nbias+1
      do 2 j=1,2
      rpmf(j,2,i)=rpmf(j,1,i)
      rpmf(j,4,i)=rpmf(j,3,i)-rpmf(j,2,i)
      rpmf(j,5,i)=rpmf(j,2,i)
      rpmf(j,6,i)=rpmf(j,3,i)
    2 continue
      m=0
      if(ipmf(i,1).eq.1) m=2
      if(ipmf(i,1).eq.2) m=3
      if(ipmf(i,1).eq.3) m=4
      if(ipmf(i,1).eq.4) m=4
      if(ipmf(i,1).eq.5) m=2
      if(ipmf(i,1).eq.6) m=1
      if(ipmf(i,1).eq.7) m=0
      if(ipmf(i,1).eq.8) m=1
      if(ipmf(i,1).eq.9) m=2
      if(ipmf(i,1).eq.10) m=1
      if(ithint) then
      if(m.eq.1) ith(18)=.true.
      if(m.eq.2) ith(18)=.true.
      if(m.eq.3) ith(20)=.true.
      if(ipmf(i,1).eq.3) ith(21)=.true.
      if(ipmf(i,1).eq.4) ith(22)=.true.
      endif
      if(m.eq.0) call md_abort('Error in pmf from top',m)
      do 3 j=1,m
      read(lfntop,1001) (jpmf(i,j,k),k=1,ipmf(i,j+1))
 1001 format(10i7)
    3 continue
      if(ipmf(i,1).eq.3.or.ipmf(i,1).eq.4) then
      do 4 k=1,6
      do 5 j=2,1,-1
      rpmf(j+1,k,i)=rpmf(j,k,i)
      rpmf(j,k,i)=one
    5 continue
    4 continue
      else
      do 6 k=1,6
      rpmf(3,k,i)=zero
    6 continue
      endif
    1 continue
      write(lfnout,1200) idpmf,npmfi
      endif
c
      call ga_sync()
      call ga_brdcst(mcf_60,ipmf,8*numpmf*lenint,0)
      call ga_brdcst(mcf_61,jpmf,4*numpmf*npmfa*lenint,0)
      call ga_brdcst(mcf_62,rpmf,18*numpmf*lendbl,0)
      if(ithint) call ga_brdcst(mcf_30,ith,24*lenlog,0)
c
      npmfi=0
      do 7 i=1,numpmf
      npmfi=max(npmfi,ipmf(i,8))
    7 continue
c
      return
      end
c $Id$
