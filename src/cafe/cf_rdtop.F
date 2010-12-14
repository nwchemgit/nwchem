      subroutine cf_rdtop(lfntop,filtop,snam)
c $Id$
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "msgids.fh"
c
      integer lfntop
      character*(*) filtop
c
      character*80 title(3)
      character*10 topdat,toptim
      real*8 releps
      character*16 stemp,snam(nsatot)
      integer natyps,nqtyps
      real*8 rdata(24),vdata(24,4),qdata(24,3),pdata(4,24)
      integer idata(24)
      character*6 cdata(24)
      integer i,ia,iq,j,k,l,m,m2,m3,m4,itemp(24),i_i,i_j,l_i,l_j
      integer naw,nbw,nhw,ndw,now,ntw,nnw
      integer nhp
      integer nas,nbt,nhs,nds,nos,nts,nxs
      integer numbl,numbi,numbd,iqfr,iqto,nsad,maxpro
      character*10 string
c
      real*8 cf_charge
      external cf_charge
c
      if(.not.ma_verify_allocator_stuff()) then
      call md_abort('ERROR IN MEMORY USE AT START RDTOP',0)
      endif
c
      natyps = 0
      nqtyps = 0
      nbw = 0
      nhw = 0
      ndw = 0
      now = 0
      ntw = 0
      nnw = 0
      naw = 0
      nas = 0
      nbt = 0
      nhs = 0
      nds = 0
      nos = 0
      nts = 0
      nxs = 0
      mset = 0
      nparms = 0
      nhop = 0
      nhp = 0
c
      if(me.eq.0) then
c
      open(unit=lfntop,file=filtop(1:index(filtop,' ')-1),
     + form='formatted',status='old',err=9999)
c
      nwc=0
      nsc=0
      nmult(1)=0
      nmult(2)=0
      nmult(3)=0
      nmult(4)=0
      totchg=0.0d0
      nhop=0
c
 1000 format(a)
      read(lfntop,1001,end=9997,err=9998) title
 1001 format(a80)
      read(lfntop,1002,end=9997,err=9998) topdat,toptim,ffield
 1002 format(12x,3a10)
      nhp=2
      if(ffield(1:6).eq.'charmm') nhp=4
      read(lfntop,1003,end=9997,err=9998) nparms
      mset=nparms
      if(lfree) mset=6
      read(lfntop,1003,end=9997,err=9998) natyps
      read(lfntop,1003,end=9997,err=9998) nqtyps
      read(lfntop,1003,end=9997,err=9998) nseq
 1003 format(i5)
      read(lfntop,1004,end=9997,err=9998) q14fac,releps
 1004 format(2f12.6)
      if(q14fac.lt.small) q14fac=one
      if(releps.lt.small) releps=one
      qfac=sqrt(1.389354428d+02/releps)
      call cf_inita(natyps,4,nqtyps,4)
      do 1 i=1,natyps
      read(lfntop,1005,end=9997,err=9998)
     + (idata(k),cdata(k),rdata(k),k=1,nparms)
 1005 format(5x,i5,1x,a6,f12.6)
      call cf_para(i,cdata,rdata,idata)
  101 continue
    1 continue
      if(.not.ma_verify_allocator_stuff()) then
      call md_abort('ERROR IN MEMORY IN RDTOP ATOM TYPES',0)
      endif
      do 2 i=1,natyps
      do 3 j=i,natyps
      read(lfntop,1006,end=9997,err=9998)
     + ((vdata(k,l),l=1,4),k=1,nparms)
 1006 format(10x,4e12.5)
      call cf_parv(i,j,vdata)
    3 continue
    2 continue
      if(.not.ma_verify_allocator_stuff()) then
      call md_abort('ERROR IN MEMORY IN RDTOP VDWAALS',0)
      endif
      do 4 i=1,nqtyps
      read(lfntop,1007,end=9997,err=9998)
     + ((qdata(k,l),l=1,3),k=1,nparms)
 1007 format(5x,f12.6,e12.5,f12.6)
      call cf_parq(i,qdata)
    4 continue
      if(.not.ma_verify_allocator_stuff()) then
      call md_abort('ERROR IN MEMORY IN RDTOP CHARGE TYPES',0)
      endif
      do 4004 i=1,nseq
      read(lfntop,4008,end=9997,err=9998) string,maxpro
 4008 format(a10,33x,i5)
      call cf_parseq(i,maxpro)
 4004 continue
      read(lfntop,1008,end=9997,err=9998) naw,nbw,nhw,ndw,now,ntw,nnw
 1008 format(5i7,2i10)
      call cf_initb(1,1,nbw,2,nhw,nhp,ndw,3,now,3,ntw,2,nnw,2,naw)
      read(lfntop,1009,end=9997,err=9998) nas,nbt,nhs,nds,nos,nts,nxs
 1009 format(5i7,3i10,2i5)
      call cf_initb(2,nas,nbt,2,nhs,nhp,nds,3,nos,3,nts,2,nxs,2,0)
      do 5 i=1,naw
      read(lfntop,1010,end=9997,err=9998) wnam(i),ia,iq
 1010 format(a16,25x,2i5)
      call cf_parwiq(i,ia,iq)
    5 continue
      idata(4)=0
      do 6 i=1,nbw
      read(lfntop,1011,end=9997,err=9998) (idata(j),j=1,3)
 1011 format(3i7)
      read(lfntop,2011,end=9997,err=9998)
     + ((pdata(k,j),k=1,2),j=1,nparms)
 2011 format(f12.6,e12.5)
      if(noshak.eq.1.or.noshak.eq.3) idata(3)=0
      call cf_parbnd(1,i,idata,pdata)
      if(idata(3).gt.0) nwc=nwc+1
    6 continue
      idata(5)=0
      do 7 i=1,nhw
      read(lfntop,1012,end=9997,err=9998) (idata(j),j=1,3)
 1012 format(3i7)
      do 70 j=1,nparms
      read(lfntop,2012,end=9997,err=9998) (pdata(k,j),k=1,nhp)
 2012 format(2(f10.6,e12.5))
   70 continue
      call cf_parang(1,i,idata,pdata)
    7 continue
      idata(6)=0
      do 8 i=1,ndw
      read(lfntop,1013,end=9997,err=9998) (idata(j),j=1,4)
 1013 format(4i7)
      read(lfntop,2013,end=9997,err=9998)
     + (itemp(j),(pdata(k,j),k=2,3),j=1,nparms)
 2013 format(i3,f10.6,e12.5)
      do 2008 j=1,nparms
      pdata(1,j)=dble(itemp(j))
 2008 continue
      call cf_pardih(1,i,idata,pdata)
    8 continue
      do 9 i=1,now
      read(lfntop,1014,end=9997,err=9998) (idata(j),j=1,4)
 1014 format(4i7)
      read(lfntop,2014,end=9997,err=9998)
     + ((pdata(k,j),k=2,3),j=1,nparms)
 2014 format(3x,f10.6,e12.5)
      do 2009 j=1,nparms
      pdata(1,j)=0.0d0
 2009 continue
      call cf_parimp(1,i,idata,pdata)
    9 continue
      if(.not.ma_push_get(mt_int,max(ntw,nnw,nts,nxs),'i',l_i,i_i))
     + call md_abort('Failed to allocate temp array',0)
      if(.not.ma_push_get(mt_int,max(ntw,nnw,nts,nxs),'j',l_j,i_j))
     + call md_abort('Failed to allocate temp array',0)
      if(ntw.gt.0) then
      read(lfntop,1015,end=9997,err=9998) (int_mb(i_i+i-1),i=1,ntw)
 1015 format(11i7)
      read(lfntop,1015,end=9997,err=9998) (int_mb(i_j+i-1),i=1,ntw)
      call cf_ndxtrd(1,int_mb(i_i),int_mb(i_j),ntw)
      endif
      if(nnw.gt.0) then
      read(lfntop,1016,end=9997,err=9998) (int_mb(i_i+i-1),i=1,nnw)
 1016 format(11i7)
      read(lfntop,1016,end=9997,err=9998) (int_mb(i_j+i-1),i=1,nnw)
      call cf_ndxxcl(1,int_mb(i_i),int_mb(i_j),nnw)
      endif
      read(lfntop,8888) string
 8888 format(a)
      read(lfntop,1024,end=9997,err=9998) (ewc(i),i=1,nparms)
 1024 format(f12.6)
      if(.not.ma_verify_allocator_stuff()) then
      call md_abort('ERROR IN MEMORY IN RDTOP SOLVENT',0)
      endif
      do 10 i=1,nas
      read(lfntop,1017,end=9997,err=9998) stemp,j,k,l,m,m2,m3,m4
 1017 format(a16,i3,2i7,19x,3i5,5x,i5)
c 1017 format(a16,3i5,15x,3i5,5x,i5)
      if(i.eq.1) then
      iqfr=min(m,m2,m3)
      iqto=max(m,m2,m3)
      else
      iqfr=min(iqfr,m,m2,m3)
      iqto=max(iqto,m,m2,m3)
      endif
      write(snam(i),'(a5,a5,i6)') stemp(1:5),stemp(11:15),l
c      snam(i)=stemp
      if(j.gt.nsf) nsf=j
c      if(j.gt.nsm) nsm=j
      totchg=totchg+cf_charge(m)
      if(k.gt.msm)
     + call md_abort('Topology and Restart are incompatible',0)
      if(m4.lt.0) nhop=nhop+1
   10 continue
      if(scaleq.ge.zero) call cf_scaleq(iqfr,iqto)
      idata(4)=0
      do 11 i=1,nbt
      read(lfntop,1018,end=9997,err=9998) (idata(j),j=1,3)
 1018 format(3i7)
      read(lfntop,2018,end=9997,err=9998)
     + ((pdata(k,j),k=1,2),j=1,nparms)
 2018 format(f12.6,e12.5)
      if(noshak.eq.2.or.noshak.eq.3) idata(3)=0
      call cf_parbnd(2,i,idata,pdata)
      if(idata(3).gt.0) nsc=nsc+1
   11 continue
      idata(5)=0
      do 12 i=1,nhs
      read(lfntop,1019,end=9997,err=9998) (idata(j),j=1,3)
 1019 format(3i7)
      do 120 j=1,nparms
      read(lfntop,2019,end=9997,err=9998) (pdata(k,j),k=1,nhp)
 2019 format(2(f10.6,e12.5))
  120 continue
      call cf_parang(2,i,idata,pdata)
   12 continue
      idata(6)=0
      do 13 i=1,nds
      read(lfntop,1020,end=9997,err=9998) (idata(j),j=1,4)
 1020 format(4i7)
      read(lfntop,2020,end=9997,err=9998)
     + (itemp(j),(pdata(k,j),k=2,3),j=1,nparms)
 2020 format(i3,f10.6,e12.5)
      do 3013 j=1,nparms
      pdata(1,j)=dble(abs(itemp(j)))
      pdata(3,j)=abs(pdata(3,j))
 3013 continue
      call cf_pardih(2,i,idata,pdata)
   13 continue
      do 14 i=1,nos
      read(lfntop,1021,end=9997,err=9998) (idata(j),j=1,4)
 1021 format(4i7)
      read(lfntop,2021,end=9997,err=9998)
     + ((pdata(k,j),k=2,3),j=1,nparms)
 2021 format(3x,f10.6,e12.5)
      do 3014 j=1,nparms
      pdata(1,j)=dble(abs(itemp(j)))
      pdata(3,j)=abs(pdata(3,j))
 3014 continue
      call cf_parimp(2,i,idata,pdata)
   14 continue
      if(nts.gt.0) then
      read(lfntop,1022,end=9997,err=9998) (int_mb(i_i+i-1),i=1,nts)
 1022 format(11i7)
      read(lfntop,1022,end=9997,err=9998) (int_mb(i_j+i-1),i=1,nts)
      endif
      call cf_ndxtrd(2,int_mb(i_i),int_mb(i_j),nts)
      if(nxs.gt.0) then
      read(lfntop,1023,end=9997,err=9998) (int_mb(i_i+i-1),i=1,nxs)
 1023 format(11i7)
      read(lfntop,1023,end=9997,err=9998) (int_mb(i_j+i-1),i=1,nxs)
      endif
      call cf_ndxxcl(2,int_mb(i_i),int_mb(i_j),nxs)
      if(.not.ma_pop_stack(l_j)) call md_abort('Dealloc failed',0)
      if(.not.ma_pop_stack(l_i)) call md_abort('Dealloc failed',0)
c
      close(unit=lfntop)
c
      if(.not.ma_verify_allocator_stuff()) then
      call md_abort('ERROR IN MEMORY IN RDTOP',1)
      endif
c
      endif
c
      if(np.gt.1) then
c
      numbl=ma_sizeof(mt_log,1,mt_byte)
      numbi=ma_sizeof(mt_int,1,mt_byte)
      numbd=ma_sizeof(mt_dbl,1,mt_byte)
c
c     broadcast dimensions
c
      idata(1)=natyps
      idata(2)=nqtyps
      idata(3)=nbw
      idata(4)=nhw
      idata(5)=ndw
      idata(6)=now
      idata(7)=ntw
      idata(8)=nnw
      idata(9)=naw
      idata(10)=nas
      idata(11)=nbt
      idata(12)=nhs
      idata(13)=nds
      idata(14)=nos
      idata(15)=nts
      idata(16)=nxs
      idata(17)=mset
      idata(18)=nparms
      idata(19)=nhop
      idata(20)=nhp
      call ga_brdcst(mcf_01,idata,20*numbi,0)
      natyps=idata(1)
      nqtyps=idata(2)
      nbw=idata(3)
      nhw=idata(4)
      ndw=idata(5)
      now=idata(6)
      ntw=idata(7)
      nnw=idata(8)
      naw=idata(9)
      nas=idata(10)
      nbt=idata(11)
      nhs=idata(12)
      nds=idata(13)
      nos=idata(14)
      nts=idata(15)
      nxs=idata(16)
      mset=idata(17)
      nparms=idata(18)
      nhop=idata(19)
      nhp=idata(20)
c
c     initialize on nodes other than 0
c
      if(me.ne.0) then
      call cf_inita(natyps,4,nqtyps,4)
      call cf_initb(1,1,nbw,2,nhw,nhp,ndw,3,now,3,ntw,2,nnw,2,naw)
      call cf_initb(2,nas,nbt,2,nhs,nhp,nds,3,nos,3,nts,2,nxs,2,0)
      endif
c
c     broadcast force field parameters
c
c      call ga_brdcst(mcf_02,byte_mb(i_nam),16*mat,0)
c
      call ga_brdcst(mcf_02,int_mb(i_typ),nparms*mat*numbi,0)
      call ga_brdcst(mcf_03,dbl_mb(i_mas),mset*mat*numbd,0)
      call ga_brdcst(mcf_04,int_mb(i_num),nparms*mat*numbi,0)
      call ga_brdcst(mcf_05,dbl_mb(i_vdw),mset*mat*mat*map*numbd,0)
      call ga_brdcst(mcf_06,dbl_mb(i_chg),mset*mqt*mqp*numbd,0)
      call ga_brdcst(mcf_07,int_mb(i_iwa),mwa*numbi,0)
      call ga_brdcst(mcf_08,int_mb(i_iwq),mwa*numbi,0)
      call ga_brdcst(mcf_09,int_mb(i_ibnd(1)),4*mbt(1)*numbi,0)
      call ga_brdcst(mcf_10,int_mb(i_ibnd(2)),4*mbt(2)*numbi,0)
      call ga_brdcst(mcf_11,dbl_mb(i_bnd(1)),mset*mbp(1)*mbt(1)*numbd,0)
      call ga_brdcst(mcf_12,dbl_mb(i_bnd(2)),mset*mbp(2)*mbt(2)*numbd,0)
      call ga_brdcst(mcf_13,int_mb(i_iang(1)),5*mht(1)*numbi,0)
      call ga_brdcst(mcf_14,int_mb(i_iang(2)),5*mht(2)*numbi,0)
      call ga_brdcst(mcf_15,dbl_mb(i_ang(1)),mset*mhp(1)*mht(1)*numbd,0)
      call ga_brdcst(mcf_16,dbl_mb(i_ang(2)),mset*mhp(2)*mht(2)*numbd,0)
      call ga_brdcst(mcf_17,int_mb(i_idih(1)),6*mdt(1)*numbi,0)
      call ga_brdcst(mcf_18,int_mb(i_idih(2)),6*mdt(2)*numbi,0)
      call ga_brdcst(mcf_19,dbl_mb(i_dih(1)),mset*mdp(1)*mdt(1)*numbd,0)
      call ga_brdcst(mcf_20,dbl_mb(i_dih(2)),mset*mdp(2)*mdt(2)*numbd,0)
      call ga_brdcst(mcf_21,int_mb(i_iimp(1)),6*mit(1)*numbi,0)
      call ga_brdcst(mcf_22,int_mb(i_iimp(2)),6*mit(2)*numbi,0)
      call ga_brdcst(mcf_23,dbl_mb(i_imp(1)),mset*mip(1)*mit(1)*numbd,0)
      call ga_brdcst(mcf_24,dbl_mb(i_imp(2)),mset*mip(2)*mit(2)*numbd,0)
      call ga_brdcst(mcf_25,int_mb(i_itrd(1)),2*(mtt(1)+1)*numbi,0)
      call ga_brdcst(mcf_26,int_mb(i_itrd(2)),2*(mtt(2)+1)*numbi,0)
      call ga_brdcst(mcf_27,int_mb(i_ixcl(1)),2*(mxt(1)+1)*numbi,0)
      call ga_brdcst(mcf_28,int_mb(i_ixcl(2)),2*(mxt(2)+1)*numbi,0)
      call ga_brdcst(mcf_29,nmult,4*numbi,0)
      call ga_brdcst(mcf_30,ith,24*numbl,0)
      call ga_brdcst(mcf_31,ip2,24*numbl,0)
      call ga_brdcst(mcf_32,ip3,24*numbl,0)
      call ga_brdcst(mcf_56,qfac,numbd,0)
c
      do 15 i=1,nsatot
      stemp=snam(i)
      call util_char_ga_brdcst(mcf_66,stemp,0)
      if(me.ne.0) snam(i)=stemp
   15 continue
c      call ga_brdcst(mcf_66,byte_mb(i_snam),16*nsatot,0)
c
c      if(lanal) then
c      call ana_select(byte_mb(i_snam))
c      call ana_initx()
c      endif
c
      itemp(1)=nwc
      itemp(2)=nsc
      itemp(3)=nsf
      call ga_brdcst(mcf_33,itemp,3*ma_sizeof(mt_int,1,mt_byte),0)
      nwc=itemp(1)
      nsc=itemp(2)
      nsf=itemp(3)
c
      endif
c
      mmult=2*nmult(1)+3*nmult(2)+4*(nmult(3)+nmult(4))
      mmuli=nmult(1)+nmult(2)+nmult(3)+nmult(4)
      if(mmult.gt.0) then
      if(.not.ma_push_get(mt_int,mmuli,'ixmul',l_ixmul,i_ixmul))
     + call md_abort('Failed to allocate memory for ixmul',0)
      if(.not.ma_push_get(mt_int,4*mmult,'imul',l_imul,i_imul))
     + call md_abort('Failed to allocate memory for imul',0)
      if(.not.ma_push_get(mt_dbl,3*mmult,'xmul',l_xmul,i_xmul))
     + call md_abort('Failed to allocate memory for xmul',0)
      if(.not.ma_push_get(mt_dbl,3*mmult,'fmul',l_fmul,i_fmul))
     + call md_abort('Failed to allocate memory for fmul',0)
      call cf_lstmul(int_mb(i_ixmul),int_mb(i_imul),
     + mbt(2),int_mb(i_ibnd(2)),mht(2),int_mb(i_iang(2)),
     + mdt(2),int_mb(i_idih(2)),mit(2),int_mb(i_iimp(2)))
      endif
c
      factmw=zero
      factms=zero
      factmp=zero
      if(noshak.eq.2.or.noshak.eq.3) then
      nsad=3*nsa-3*ndums
      else
      nsad=3*nsa-2*ndums
      endif
      if(nwm*(3*nwa-nwc)-3*islow.gt.0)
     + factmw=two/(rgas*dble(nwm*(3*nwa-nwc)-3*islow))
      if(nsad-nsc-3*islow.gt.0)
     + factms=two/(rgas*dble(nsad-nsc-3*islow))
      if(nwm*(3*nwa-nwc)+nsad-nsc-3*islow.gt.0)
     + factmp=two/(rgas*dble(nwm*(3*nwa-nwc)+nsad-nsc-3*islow))
c
      if(.not.ma_verify_allocator_stuff()) then
      call md_abort('ERROR IN MEMORY USE AFTER RDTOP',0)
      endif
c
      return
c
 9997 continue
      call md_abort('EOF reading topology file',0)
 9998 continue
      call md_abort('Error reading topology file',0)
 9999 continue
      call md_abort('Failed to open topology file',0)
      return
      end
      subroutine cf_topol_init(lfnout)
c
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      integer numbi,lfnout
c
      numbi=ma_sizeof(mt_int,1,mt_byte)
c
      if(me.eq.0) then
      iffld=1
      if(ffield(1:6).eq.'charmm') iffld=2
      if(ffield(1:6).eq.'gromos') iffld=3
      endif
      call ga_brdcst(mcf_69,iffld,numbi,0)
c
c      if(iffld.eq.1) then
c      if(me.eq.0) write(lfnout,1001)
c 1001 format(' Force field is amber')
c      endif
c      if(iffld.eq.2) then
c      if(me.eq.0) write(lfnout,1002)
c 1002 format(' Force field is charmm')
c      endif
c      if(iffld.eq.3) then
c      if(me.eq.0) write(lfnout,1003)
c 1003 format(' Force field is gromos')
c      endif
c
      return
      end
