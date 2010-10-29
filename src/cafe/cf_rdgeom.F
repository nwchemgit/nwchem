      subroutine cf_rdgeom(snam)
c $Id$
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "msgids.fh"
#include "rtdb.fh"
#include "geom.fh"
c
      real*8 releps,c(3),q,w
      character*16 tag,snam(nsatot)
      integer natyps,nqtyps,igeom
c
      real*8 rdata(9),vadata(6),vrdata(6)
      integer idata(16)
      character*6 cdata(3)
      integer i,j,itemp(4)
      integer naw,nbw,nhw,ndw,now,ntw,nnw
      integer nas,nbt,nhs,nds,nos,nts,nxs
      integer numbl,numbi,numbd
c
      logical lstate
c
      lstate=rtdb_parallel(.false.)
c
      if(me.eq.0) then
c
      nwc=0
      nsc=0
      nmult(1)=0
      nmult(2)=0
      nmult(3)=0
      nmult(4)=0
c
      q14fac=one
      releps=one
      qfac=sqrt(1.389354428d+02/releps)
      natyps=nsatot
      nqtyps=nsatot
c
      call cf_inita(natyps,4,nqtyps,4)
c
      if(.not.geom_create(igeom,'geometry'))
     + call md_abort('cf_rdgeom: geom_create failed',0)
      if(.not.geom_rtdb_load(irtdb,igeom,'geometry'))
     + call md_abort('cf_rdgeom: geom_rtdb_load failed',0)
      if(.not.geom_ncent(igeom,nsa))
     + call md_abort('cf_rdgeom: geom_ncent failed',0)
c
      do 1 i=1,natyps
      if(.not.geom_cent_get(igeom,i,tag,c,q))
     + call md_abort('cf_rdgeom: geom_cent_get failed',i)
      if(.not.geom_mass_get(igeom,i,w))
     + call md_abort('cf_rdgeom: geom_mass_get failed',i)
      snam(i)='      '//tag(1:10)
      do 2 j=1,3
      cdata(j)=tag(1:6)
      rdata(j)=w
      idata(j)=int(q)
      vadata(2*j-1)=q
      vadata(2*j)=0.0d0
    2 continue
      call cf_para(i,cdata(1),rdata(1),idata(1))
      call cf_parq(i,vadata)
    1 continue
      if(.not.geom_destroy(igeom))
     + call md_abort('arg_rdqdim: geom_destroy failed',0)
c
      do 3 i=1,6
      vadata(i)=0.0d0
      vrdata(i)=0.0d0
    3 continue
      do 4 i=1,natyps
      do 5 j=i,natyps
      call cf_parv(i,j,vadata,vrdata)
    5 continue
    4 continue
c
      naw=0
      nbw=0
      nhw=0
      ndw=0
      now=0
      ntw=0
      nnw=0
      call cf_initb(1,1,nbw,2,nhw,2,ndw,3,now,3,ntw,2,nnw,2,naw)
c
      nas=nsatot
      nbt=0
      nhs=0
      nds=0
      nos=0
      nts=0
      nxs=0
      call cf_initb(2,nas,nbt,2,nhs,2,nds,3,nos,3,nts,2,nxs,2,0)
c
      ewc(1)=0.0d0
      ewc(2)=0.0d0
      ewc(3)=0.0d0
c
      endif
c
      lstate=rtdb_parallel(.true.)
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
      call ga_brdcst(mcf_01,idata,16*numbi,0)
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
c
c     initialize on nodes other than 0
c
      if(me.ne.0) then
      call cf_inita(natyps,4,nqtyps,4)
      call cf_initb(1,1,nbw,2,nhw,2,ndw,3,now,3,ntw,2,nnw,2,naw)
      call cf_initb(2,nas,nbt,2,nhs,2,nds,3,nos,3,nts,2,nxs,2,0)
      endif
c
c     broadcast force field parameters
c
c      call ga_brdcst(mcf_02,byte_mb(i_nam),16*mat,0)
c
      call ga_brdcst(mcf_02,int_mb(i_typ),3*mat,0)
      call ga_brdcst(mcf_03,dbl_mb(i_mas),6*mat*numbd,0)
      call ga_brdcst(mcf_04,int_mb(i_num),3*mat*numbi,0)
      call ga_brdcst(mcf_05,dbl_mb(i_vdw),6*mat*mat*map*numbd,0)
      call ga_brdcst(mcf_06,dbl_mb(i_chg),6*mqt*mqp*numbd,0)
      call ga_brdcst(mcf_07,int_mb(i_iwa),mwa*numbi,0)
      call ga_brdcst(mcf_08,int_mb(i_iwq),mwa*numbi,0)
      call ga_brdcst(mcf_09,int_mb(i_ibnd(1)),4*mbt(1)*numbi,0)
      call ga_brdcst(mcf_10,int_mb(i_ibnd(2)),4*mbt(2)*numbi,0)
      call ga_brdcst(mcf_11,dbl_mb(i_bnd(1)),6*mbp(1)*mbt(1)*numbd,0)
      call ga_brdcst(mcf_12,dbl_mb(i_bnd(2)),6*mbp(2)*mbt(2)*numbd,0)
      call ga_brdcst(mcf_13,int_mb(i_iang(1)),5*mht(1)*numbi,0)
      call ga_brdcst(mcf_14,int_mb(i_iang(2)),5*mht(2)*numbi,0)
      call ga_brdcst(mcf_15,dbl_mb(i_ang(1)),6*mhp(1)*mht(1)*numbd,0)
      call ga_brdcst(mcf_16,dbl_mb(i_ang(2)),6*mhp(2)*mht(2)*numbd,0)
      call ga_brdcst(mcf_17,int_mb(i_idih(1)),6*mdt(1)*numbi,0)
      call ga_brdcst(mcf_18,int_mb(i_idih(2)),6*mdt(2)*numbi,0)
      call ga_brdcst(mcf_19,dbl_mb(i_dih(1)),6*mdp(1)*mdt(1)*numbd,0)
      call ga_brdcst(mcf_20,dbl_mb(i_dih(2)),6*mdp(2)*mdt(2)*numbd,0)
      call ga_brdcst(mcf_21,int_mb(i_iimp(1)),6*mit(1)*numbi,0)
      call ga_brdcst(mcf_22,int_mb(i_iimp(2)),6*mit(2)*numbi,0)
      call ga_brdcst(mcf_23,dbl_mb(i_imp(1)),6*mip(1)*mit(1)*numbd,0)
      call ga_brdcst(mcf_24,dbl_mb(i_imp(2)),6*mip(2)*mit(2)*numbd,0)
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
      do 6 i=1,nsatot
      tag=snam(i)
      call util_char_ga_brdcst(mcf_66,tag,0)
      if(me.ne.0) snam(i)=tag
    6 continue
c
c      call ga_brdcst(mcf_66,byte_mb(i_snam),16*nsatot,0)
c
c      if(lanal) then
c      call ana_select(byte_mb(i_snam))
c      call ana_initx()
c      endif
c
      itemp(1)=nwc
      itemp(2)=nsc
      call ga_brdcst(mcf_33,itemp,2*ma_sizeof(mt_int,1,mt_byte),0)
      nwc=itemp(1)
      nsc=itemp(2)
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
      if(nwm*(3*nwa-nwc)-3*islow.gt.0)
     + factmw=two/(rgas*dble(nwm*(3*nwa-nwc)-3*islow))
      if(3*nsa-nsc-3*islow.gt.0)
     + factms=two/(rgas*dble(3*nsa-nsc-3*islow))
      if(nwm*(3*nwa-nwc)+3*nsa-nsc-3*islow.gt.0)
     + factmp=two/(rgas*dble(nwm*(3*nwa-nwc)+3*nsa-nsc-3*islow))
c
      return
      end
