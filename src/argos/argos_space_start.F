      subroutine argos_space_start(lout,lfntop,filtop,lfnrst,filrst,
     + lsyn,filsyn,nsyn,
     + rsht,rlng,rinp,rsgmi,
     + npxi,npyi,npzi,nbxi,nbyi,nbzi,
     + npbt,nbxt,boxt,vlatt,lpbc,
     + nwmi,mwmi,nwai,mwai,nsfi,msfi,nsmi,msmi,nsai,msai,
     + ldb,lbp,fld,lpol,lfre,temp,tempw,temps,iguidi,
     + ldbg,idbg,prjct,mbbi,nseri,isld,irset,ictrl,nseqi,
     + i_lseqi,ndumsi,nbgeti,npreci,madbxi)
c
c
c $Id$
c
      implicit none
c
#include "argos_space_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "global.fh"
c
      logical argos_space_rdrest,argos_space_diffbb
      external argos_space_rdrest,argos_space_diffbb
c
      integer lfntop,lfnrst,lsyn,nsyn,ldbg,idbg,lout
      character*255 filtop,filrst,filsyn
      real*8 rsht,rlng,rinp,boxt(3),vlatt(3,3),rtemp(3)
      real*8 temp,tempw,temps,fld,rsgmi
      integer itemp(2)
      integer npxi,npyi,npzi,nbxi,nbyi,nbzi,isld,irset,nseqi,i_lseqi
      integer nwmi,mwmi,nwai,mwai,nsfi,msfi,nsmi,msmi,nsai,msai
      integer lenscr,ldb,lbp,npbt,nbxt,iguidi,mbbi,nseri,ictrl,ndumsi
      integer nbgeti,npreci,madbxi
      logical lpol,lfre,ignore,lpbc
      character*80 prjct
c
      integer i,j,l_f,i_f
c
      if(.not.ma_verify_allocator_stuff()) then
      call md_abort('ERROR IN MEMORY USE AT START ARGOS_SPACE_START',0)
      endif
c
      me=ga_nodeid()
      np=ga_nnodes()
c
      project=prjct
c
      idebug=idbg
      lfndbg=ldbg
      icntrl=ictrl
      lfnout=lout
c
      lfnsyn=lsyn
      nfsync=nsyn
      isload=isld
      ireset=irset
      nbget=nbgeti
      nprec=npreci
      ibget=nbget
      madbox=madbxi
c
      rshort=rsht
      rlong=max(rsht,rlng)
      rbox=rinp
      ntwin=1
      if(rlng.gt.rsht) ntwin=2
      ltwin=rlng.gt.rsht
c
      loadb=ldb
      lbpair=lbp
      factld=fld
c
      lpola=lpol
      lfree=lfre
c
      nbx=nbxi
      nby=nbyi
      nbz=nbzi
      nbxin=nbxi
      nbyin=nbyi
      nbzin=nbzi
c
      npx=npxi
      npy=npyi
      npz=npzi
c
      mbbl=0
      nbbdif=-1
c
      mwm=mwmi
      msa=max(msai,msmi)
      mbbreq=mbbi
      mbblp=0
      nserie=nseri
c
      iguide=iguidi
c
      npack=0
      npackw=0
c
      rsgm=rsgmi
c
      call argos_space_nrnode()
c
      call argos_space_dimens(lfnrst,filrst)
c
      call argos_space_alloc()
c
      call argos_space_decomp(int_mb(i_iown),dbl_mb(i_boxs),
     + int_mb(i_buren),int_mb(i_bindex))
c
      call argos_space_initip(int_mb(i_iown),int_mb(i_ipl))
c
      call argos_space_numbb(int_mb(i_iown),dbl_mb(i_boxs))
c
      if(argos_space_diffbb(dbl_mb(i_boxs),int_mb(i_rng)))
     + call argos_space_listbb(int_mb(i_iown),dbl_mb(i_boxs),
     + int_mb(i_bb))
c
      call argos_space_alloc2()
c
      if(ireset.eq.0) then
      ignore=argos_space_rdrest(lfnrst,filrst,dbl_mb(i_boxs))
      endif
      if(np.gt.1) then
      call ga_brdcst(msp_02,nsm,ma_sizeof(mt_int,1,mt_byte),0)
      endif
      msm=max(1,nsm)
      if(.not.ma_push_get(mt_dbl,msm*3,'xscr',l_xscr,i_xscr))
     + call md_abort('Failed to allocate xscr',0)
c
      lenscr=3*max(mwm*mwa,msa)
      if(.not.ma_push_get(mt_dbl,lenscr,'x',l_x,i_x))
     + call md_abort('Failed to allocate x',0)
      if(.not.ma_push_get(mt_dbl,lenscr,'v',l_v,i_v))
     + call md_abort('Failed to allocate v',0)
      if(.not.ma_push_get(mt_dbl,lenscr,'f',l_f,i_f))
     + call md_abort('Failed to allocate f',0)
      if(.not.ma_push_get(mt_dbl,3*mwm*mwa,'r',l_r,i_r))
     + call md_abort('Failed to allocate r',0)
      lenscr=max(miw2*mwm,mis2*msa)
      if(.not.ma_push_get(mt_int,lenscr,'i',l_i,i_i))
     + call md_abort('Failed to allocate i',0)
      lenscr=ma_inquire_avail(mt_byte)/
     + ((9*mwa+3)*ma_sizeof(mt_dbl,1,mt_byte)+
     + (mis2+4)*ma_sizeof(mt_int,1,mt_byte))-1
      if(.not.ma_push_get(mt_dbl,lenscr*3*mwa,'bx',l_bx,i_bx))
     + call md_abort('Failed to allocate bx',0)
      if(.not.ma_push_get(mt_dbl,lenscr*3*mwa,'bv',l_bv,i_bv))
     + call md_abort('Failed to allocate bv',0)
      if(.not.ma_push_get(mt_dbl,lenscr*3*mwa,'bf',l_bf,i_bf))
     + call md_abort('Failed to allocate bf',0)
      if(.not.ma_push_get(mt_dbl,lenscr*3,'br',l_br,i_br))
     + call md_abort('Failed to allocate br',0)
      if(.not.ma_push_get(mt_int,lenscr*max(mis2,2),'bi',l_bi,i_bi))
     + call md_abort('Failed to allocate bi',0)
      if(.not.ma_push_get(mt_int,lenscr,'n',l_n,i_n))
     + call md_abort('Failed to allocate n',0)
c
      call ga_sync()
c
      if(.not.ma_verify_allocator_stuff()) then
      call md_abort
     + ('ERROR IN MEMORY USE BEFORE RDRST IN ARGOS_SPACE_START',0)
      endif
c
      call argos_space_rdrst(lfnrst,filrst,lfntop,filtop,
     + temp,tempw,temps,int_mb(i_ipl),
     + dbl_mb(i_x),dbl_mb(i_v),dbl_mb(i_f),dbl_mb(i_r),int_mb(i_i),
     + dbl_mb(i_x),dbl_mb(i_v),dbl_mb(i_f),dbl_mb(i_xscr),int_mb(i_i),
     + dbl_mb(i_bx),dbl_mb(i_bv),dbl_mb(i_bf),dbl_mb(i_br),
     + int_mb(i_bi),lenscr,
     + dbl_mb(i_bx),dbl_mb(i_bv),dbl_mb(i_bf),int_mb(i_bi),
     + lenscr,int_mb(i_n),
     + int_mb(i_iown),dbl_mb(i_boxs),int_mb(i_lseq),int_mb(i_sndx))
c
      if(.not.ma_verify_allocator_stuff()) then
      call md_abort
     + ('ERROR IN MEMORY USE AFTER RDRST IN ARGOS_SPACE_START',0)
      endif
c
      if(.not.ma_pop_stack(l_n))
     + call md_abort('Failed to deallocate n',0)
      if(.not.ma_pop_stack(l_bi))
     + call md_abort('Failed to deallocate bi',0)
      if(.not.ma_pop_stack(l_br))
     + call md_abort('Failed to deallocate br',0)
      if(.not.ma_pop_stack(l_bf))
     + call md_abort('Failed to deallocate bf',0)
      if(.not.ma_pop_stack(l_bv))
     + call md_abort('Failed to deallocate bv',0)
      if(.not.ma_pop_stack(l_bx))
     + call md_abort('Failed to deallocate bx',0)
      if(.not.ma_pop_stack(l_i))
     + call md_abort('Failed to deallocate i',0)
      if(.not.ma_pop_stack(l_r))
     + call md_abort('Failed to deallocate r',0)
      if(.not.ma_pop_stack(l_f))
     + call md_abort('Failed to deallocate f',0)
      if(.not.ma_pop_stack(l_v))
     + call md_abort('Failed to deallocate v',0)
      if(.not.ma_pop_stack(l_x))
     + call md_abort('Failed to deallocate x',0)
c
      if(np.gt.1) then
      rtemp(1)=temp
      rtemp(2)=tempw
      rtemp(3)=temps
      call ga_brdcst(msp_01,rtemp,ma_sizeof(mt_dbl,3,mt_byte),0)
      temp=rtemp(1)
      tempw=rtemp(2)
      temps=rtemp(3)
      itemp(1)=nsm
      itemp(2)=nsf
      call ga_brdcst(msp_02,itemp,ma_sizeof(mt_int,2,mt_byte),0)
      nsm=itemp(1)
      nsf=itemp(2)
      endif
c
      nwmi=nwm
      nwai=nwa
      nsmi=nsm
      nsai=nsa
c
      mwmi=mwm
      mwai=mwa
      msai=msa
c
      npbt=npbtyp
      nbxt=nbxtyp
c
      do 1 i=1,3
      boxt(i)=box(i)
      do 2 j=1,3
      vlatt(i,j)=vlat(i,j)
    2 continue
    1 continue
c
      nsfi=nsf
      msfi=max(1,nsf)
      msmi=max(1,nsm)
c
      third=1.0d0/3.0d0
      nldup=-2
      ipairf=-1
      ipairt=-1
      lpipo=.false.
c
      lpbc=npbtyp.gt.0
c
      if(me.eq.0) then
      if(nfsync.gt.0) then
      open(unit=lfnsyn,file=filsyn(1:index(filsyn,' ')-1),
     + status='unknown')
      write(lfnsyn,3000) np
 3000 format(i5)
      endif
      endif
c
      nseqi=nseq
      i_lseqi=i_lseq
c
      ndumsi=ndums
c
      if(.not.ma_verify_allocator_stuff()) then
      call md_abort('ERROR IN MEMORY USE AT EXIT ARGOS_SPACE_START',0)
      endif
c
      return
      end
