      subroutine sp_start(lfntop,filtop,lfnrst,filrst,lsyn,filsyn,nsyn,
     + rsht,rlng,rinp,
     + npxi,npyi,npzi,nbxi,nbyi,nbzi,
     + npbt,nbxt,boxt,vlatt,lpbc,
     + nwmi,mwmi,nwai,mwai,nsfi,msfi,nsmi,msmi,nsai,msai,
     + ldb,lbp,fld,lpol,lfre,temp,tempw,temps,lqd,iguidi)
c
c
c $Id: sp_init.F,v 1.40 1999-10-26 19:32:29 d3j191 Exp $
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "global.fh"
c
      logical sp_rdrest,sp_diffbb
      external sp_rdrest,sp_diffbb
c
      integer lfntop,lfnrst,lsyn,nsyn
      character*255 filtop,filrst,filsyn
      real*8 rsht,rlng,rinp,boxt(3),vlatt(3,3),rtemp(3)
      real*8 temp,tempw,temps,fld
      integer npxi,npyi,npzi,nbxi,nbyi,nbzi
      integer nwmi,mwmi,nwai,mwai,nsfi,msfi,nsmi,msmi,nsai,msai
      integer lenscr,ldb,lbp,npbt,nbxt,iguidi
      logical lpol,lfre,ignore,lpbc,lqd
c
      integer i,j
c
      me=ga_nodeid()
      np=ga_nnodes()
c
      lfnsyn=lsyn
      nfsync=nsyn
c
      rshort=rsht
      rlong=max(rsht,rlng)
      rbox=rinp
      ntwin=1
      if(rlng.gt.rsht) ntwin=2
      ltwin=rlng.gt.rsht
      lnode0=lqmd
      lqmd=lqd
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
c
      npx=npxi
      npy=npyi
      npz=npzi
c
      mbbl=0
      nbbdif=-1
c
      mwm=mwmi
      msa=msai
c
      iguide=iguidi
c
      call sp_nrnode()
c
      call sp_dimens(lfnrst,filrst)
c
      call sp_alloc()
c
      call sp_decomp(int_mb(i_iown),dbl_mb(i_boxs))
c
      call sp_initip(int_mb(i_iown),int_mb(i_ipl))
c
      call sp_numbb(int_mb(i_iown),dbl_mb(i_boxs))
c
      if(sp_diffbb(dbl_mb(i_boxs),int_mb(i_rng)))
     + call sp_listbb(int_mb(i_iown),dbl_mb(i_boxs),int_mb(i_bb))
c
      ignore=sp_rdrest(lfnrst,filrst)
c
      lenscr=3*max(mwm*mwa,msa)
      if(.not.ma_push_get(mt_dbl,lenscr,'x',l_x,i_x))
     + call errquit('Failed to allocate x',0)
      if(.not.ma_push_get(mt_dbl,lenscr,'v',l_v,i_v))
     + call errquit('Failed to allocate v',0)
      if(.not.ma_push_get(mt_dbl,3*mwm*mwa,'r',l_r,i_r))
     + call errquit('Failed to allocate r',0)
      lenscr=max(miw2*mwm,mis2*msa)
      if(.not.ma_push_get(mt_int,lenscr,'i',l_i,i_i))
     + call errquit('Failed to allocate i',0)
      lenscr=ma_inquire_avail(mt_byte)/
     + ((6*mwa+3)*ma_sizeof(mt_dbl,1,mt_byte)+
     + (mis2+4)*ma_sizeof(mt_int,1,mt_byte))-1
      if(.not.ma_push_get(mt_dbl,lenscr*3*mwa,'bx',l_bx,i_bx))
     + call errquit('Failed to allocate bx',0)
      if(.not.ma_push_get(mt_dbl,lenscr*3*mwa,'bv',l_bv,i_bv))
     + call errquit('Failed to allocate bv',0)
      if(.not.ma_push_get(mt_dbl,lenscr*3,'br',l_br,i_br))
     + call errquit('Failed to allocate br',0)
      if(.not.ma_push_get(mt_int,lenscr*max(mis2,2),'bi',l_bi,i_bi))
     + call errquit('Failed to allocate bi',0)
      if(.not.ma_push_get(mt_int,lenscr,'n',l_n,i_n))
     + call errquit('Failed to allocate n',0)
c
      call ga_sync()
      call sp_rdrst(lfnrst,filrst,lfntop,filtop,
     + temp,tempw,temps,int_mb(i_ipl),
     + dbl_mb(i_x),dbl_mb(i_v),dbl_mb(i_r),int_mb(i_i),
     + dbl_mb(i_x),dbl_mb(i_v),int_mb(i_i),
     + dbl_mb(i_bx),dbl_mb(i_bv),dbl_mb(i_br),int_mb(i_bi),lenscr,
     + dbl_mb(i_bx),dbl_mb(i_bv),int_mb(i_bi),lenscr,int_mb(i_n),
     + int_mb(i_iown),dbl_mb(i_boxs))
c
      if(.not.ma_pop_stack(l_n))
     + call errquit('Failed to deallocate n',0)
      if(.not.ma_pop_stack(l_bi))
     + call errquit('Failed to deallocate bi',0)
      if(.not.ma_pop_stack(l_br))
     + call errquit('Failed to deallocate br',0)
      if(.not.ma_pop_stack(l_bv))
     + call errquit('Failed to deallocate bv',0)
      if(.not.ma_pop_stack(l_bx))
     + call errquit('Failed to deallocate bx',0)
      if(.not.ma_pop_stack(l_i))
     + call errquit('Failed to deallocate i',0)
      if(.not.ma_pop_stack(l_r))
     + call errquit('Failed to deallocate r',0)
      if(.not.ma_pop_stack(l_v))
     + call errquit('Failed to deallocate v',0)
      if(.not.ma_pop_stack(l_x))
     + call errquit('Failed to deallocate x',0)
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
      if(np.gt.1) then
      rtemp(1)=temp
      rtemp(2)=tempw
      rtemp(3)=temps
      call ga_brdcst(msp_01,rtemp,ma_sizeof(mt_dbl,3,mt_byte),0)
      temp=rtemp(1)
      tempw=rtemp(2)
      temps=rtemp(3)
      call ga_brdcst(msp_02,nsf,ma_sizeof(mt_int,1,mt_byte),0)
      endif
c
      nsfi=nsf
      msfi=max(1,nsf)
      msmi=max(1,nsm)
c
      third=1.0d0/3.0d0
      nldup=-1
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
      call ga_sync()
c
      return
      end
      subroutine sp_dimens(lfnrst,filrst)
c
      implicit none
c
#include "sp_common.fh"
#include "msgids.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      integer lfnrst
      character*255 filrst
c
      character*1 cdum
      integer i,j,jdum,ibx,iby,ibz,itemp(9)
      real*8 rtemp(4)
c
      mbbl=0
c
      if(me.eq.0) then
c
      open(unit=lfnrst,file=filrst(1:index(filrst,' ')-1),
     + status='old',form='formatted',err=9999)
      rewind(lfnrst)
c
      do 1 i=1,4
      read(lfnrst,1001) cdum
 1001 format(a1)
    1 continue
      read(lfnrst,1002) npbtyp,nbxtyp
 1002 format(2i5)
      read(lfnrst,1004) ((vlat(i,j),j=1,3),i=1,3)
 1004 format(3f12.6)
      box(1)=vlat(1,1)
      box(2)=vlat(2,2)
      box(3)=vlat(3,3)
      read(lfnrst,1003) jdum
 1003 format(40x,i5)
      read(lfnrst,1001) cdum
      if(jdum.ne.0) then
      read(lfnrst,1001) cdum
      endif
      read(lfnrst,1005) nwm,nwa,nsm,nsa,nwmc
 1005 format(5i10)
      close(unit=lfnrst,status='keep')
c
      bsize=max(rshort,half*rlong,rbox)
      bxmin=bsize
      bymin=bsize
      bzmin=bsize
c
      if(nbx*nby*nbz.lt.np) then
      nbx=int(box(1)/bsize)
      nby=int(box(2)/bsize)
      nbz=int(box(3)/bsize)
      endif
      if(nbx.eq.0) nbx=1
      if(nby.eq.0) nby=1
      if(nbz.eq.0) nbz=1
c
      if(nbx.lt.npx) then
      i=npx/nbx
      if(mod(npx,nbx).gt.0) i=i+1
      bxmin=min(bsize,box(1)/dble(i))
      nbx=npx
      endif
      if(nby.lt.npy) then
      i=npy/nby
      if(mod(npy,nby).gt.0) i=i+1
      bymin=min(bsize,box(2)/dble(i))
      nby=npy
      endif
      if(nbz.lt.npz) then
      i=npz/nbz
      if(mod(npz,nbz).gt.0) i=i+1
      bzmin=min(bsize,box(3)/dble(i))
      nbz=npz
      endif
c
      endif
c
      if(np.gt.1) then
      itemp(1)=nwm
      itemp(2)=nwa
      itemp(3)=nsm
      itemp(4)=nsa
      itemp(5)=nbx
      itemp(6)=nby
      itemp(7)=nbz
      itemp(8)=npbtyp
      itemp(9)=nbxtyp
      rtemp(1)=bsize
      rtemp(2)=1.001d0*bxmin
      rtemp(3)=1.001d0*bymin
      rtemp(4)=1.001d0*bzmin
      call ga_brdcst(msp_03,itemp,ma_sizeof(mt_int,9,mt_byte),0)
      call ga_brdcst(msp_04,rtemp,ma_sizeof(mt_dbl,4,mt_byte),0)
      call ga_brdcst(msp_05,box,ma_sizeof(mt_dbl,3,mt_byte),0)
      call ga_brdcst(msp_06,vlat,ma_sizeof(mt_dbl,9,mt_byte),0)
      nwm=itemp(1)
      nwa=itemp(2)
      nsm=itemp(3)
      nsa=itemp(4)
      nbx=itemp(5)
      nby=itemp(6)
      nbz=itemp(7)
      npbtyp=itemp(8)
      nbxtyp=itemp(9)
      bsize=rtemp(1)
      bxmin=rtemp(2)
      bymin=rtemp(3)
      bzmin=rtemp(4)
      endif
c
      boxh(1)=half*box(1)
      boxh(2)=half*box(2)
      boxh(3)=half*box(3)
c
      maxbox=max(nbx,nby,nbz)
      nbtot=nbx*nby*nbz
      lpbc0=nbx.eq.1.or.nby.eq.1.or.nbz.eq.1.or.
     + npx.eq.1.or.npy.eq.1.or.npz.eq.1.or.lpbc9
c
      mbox=30
      do 2 ibx=1,nbx
      do 3 iby=1,nby
      do 4 ibz=1,nbz
      if(me.eq.((((ibz-1)*npz)/nbz)*npy+(((iby-1)*npy)/nby))*npx
     + +((ibx-1)*npx)/nbx) mbox=mbox+1
    4 continue
    3 continue
    2 continue
      mbxloc=mbox-30
c
      if(np.gt.1) call ga_igop(msp_07,mbox,1,'max')
c
      mwa=max(1,nwa)
      msa=max(1,msa,(mbox-30+madbox)*((nwm*nwa+nsa)/nbtot+1)+1)
      mwm=max(1,mwm,msa/mwa+1)
c
      msa=min(msa,2*nsa+1)
      mwm=min(mwm,2*nwm+1)
c
      if(lnode0) then
      msa=nsa+1
      mwm=nwm+1
      endif
c
      return
 9999 call errquit('Failed to open restart file',0)
      return
      end
      subroutine sp_alloc
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      integer isize
c
      if(.not.ma_push_get(mt_int,3*maxbox,'owner',l_iown,i_iown))
     + call errquit('Failed to allocate owner',0)
      if(.not.ma_push_get(mt_dbl,3*maxbox,'bxsiz',l_boxs,i_boxs))
     + call errquit('Failed to allocate bxsiz',0)
      if(.not.ma_push_get(mt_int,6*maxbox,'ibxrg',l_boxr,i_boxr))
     + call errquit('Failed to allocate ibxrg',0)
      if(.not.ma_push_get(mt_int,6*maxbox,'rng',l_rng,i_rng))
     + call errquit('Failed to allocate rng',me)
c
      if(.not.ma_push_get(mt_int,mip2*mbox,'ipl',l_ipl,i_ipl))
     + call errquit('Failed to allocate ipl',0)
      if(.not.ma_push_get(mt_int,mip2*mbox,'jpl',l_jpl,i_jpl))
     + call errquit('Failed to allocate jpl',0)
c
      if(.not.ga_create(mt_int,np*mbox,mip2,'ip',mbox,mip2,ga_ip))
     + call errquit('Failed to create global array ip',0)
      if(.not.ga_create(mt_int,np*mwm,miw2,'iw',mwm,miw2,ga_iw))
     + call errquit('Failed to create global array iw',0)
      isize=6+12*mwa
      if(lpola) isize=6+18*mwa
      if(lpola.and.lfree) isize=6+30*mwa
      if(.not.ga_create(mt_dbl,np*mwm,isize,'w',mwm,isize,ga_w))
     + call errquit('Failed to create global array w',0)
      if(.not.ga_create(mt_int,np*msa,mis2,'is',msa,mis2,ga_is))
     + call errquit('Failed to create global array is',0)
      isize=39
      if(lpola) isize=45
      if(lpola.and.lfree) isize=57
      if(.not.ga_create(mt_dbl,np*msa,isize,'s',msa,isize,ga_s))
     + call errquit('Failed to create global array s',0)
c
      if(.not.ga_create(mt_int,np*mwm,1,'iwz',mwm,1,ga_iwz))
     + call errquit('Failed to create global array iwz',0)
      if(.not.ga_create(mt_int,np*msa,1,'isz',msa,1,ga_isz))
     + call errquit('Failed to create global array isz',0)
c
      return
      end
      subroutine sp_free()
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      if(.not.ga_destroy(ga_isz))
     + call errquit('Error ga_destroy isz',ga_isz)
      if(.not.ga_destroy(ga_iwz))
     + call errquit('Error ga_destroy iwz',ga_iwz)
c
      if(.not.ga_destroy(ga_s))
     + call errquit('Error ga_destroy s',ga_s)
      if(.not.ga_destroy(ga_is))
     + call errquit('Error ga_destroy is',ga_is)
      if(.not.ga_destroy(ga_w))
     + call errquit('Error ga_destroy w',ga_w)
      if(.not.ga_destroy(ga_iw))
     + call errquit('Error ga_destroy iw',ga_iw)
      if(.not.ga_destroy(ga_ip))
     + call errquit('Error ga_destroy iw',ga_ip)
c
      if(.not.ma_pop_stack(l_jpl))
     + call errquit('Failed to deallocate jpl',0)
      if(.not.ma_pop_stack(l_ipl))
     + call errquit('Failed to deallocate ipl',0)
      if(.not.ma_pop_stack(l_rng))
     + call errquit('Failed to deallocate rng',0)
      if(.not.ma_pop_stack(l_boxr))
     + call errquit('Failed to deallocate boxr',0)
      if(.not.ma_pop_stack(l_boxs))
     + call errquit('Failed to deallocate boxs',0)
      if(.not.ma_pop_stack(l_iown))
     + call errquit('Failed to deallocate iown',0)
c
      return
      end
      subroutine sp_rdrst(lfnrst,filrst,lfntop,filtop,
     + temp,tempw,temps,ipl,xw,vw,xwcr,iwl,xs,vs,isl,
     + bxw,bvw,brw,ibw,nw,bxs,bvs,ibs,ns,ndx,
     + ibownr,boxsiz)
c
      implicit none
c
#include "sp_common.fh"
c
      integer sp_btop
      external sp_btop
c
      integer lfnrst,lfntop
      character*255 filrst,filtop
      integer nw,ns
      integer ibownr(maxbox,3)
      integer ipl(mbox,mip2),iwl(mwm,miw2),isl(msa,mis2)
      real*8 xw(mwm,3,mwa),xs(msa,3),xwcr(mwm,3)
      real*8 vw(mwm,3,mwa),vs(msa,3)
      real*8 bxw(nw,3,mwa),bxs(ns,3),brw(nw,3)
      real*8 bvw(nw,3,mwa),bvs(ns,3)
      integer ibw(nw,2),ibs(ns,mis2)
      real*8 boxsiz(maxbox,3)
      integer ndx(nw)
      real*8 temp,tempw,temps
c
      character*1 cdum
      real*8 rdum,cgx,cgy,cgz,xt(3)
      integer i,j,k,idum,jdum,kdum,number,ncyc,numw,nums
      integer icyc,ibx,iby,ibz,ipx,ipy,ipz,node,new,nold
      integer ilw,ihw,jlw,jhw,ils,ihs,jls,jhs
      integer ili,ihi,jli,jhi,ilp,ihp,jlp,jhp
c
      integer nat,nqt,naw,nbw,nhw,ndw,now,ntw,nnw,ipp,ips,ix
      real*8 boxi(3)
#if defined(CRAY_T3D) || defined(CRAY_T3E)
      real*8 arg_nint
      external arg_nint
#endif
c
      lpbc9=.true.
      nsf=0
c
      if(me.eq.0) then
c
      open(unit=lfnrst,file=filrst(1:index(filrst,' ')-1),
     + status='old',form='formatted',err=9999)
      rewind(lfnrst)
c
      open(unit=lfntop,file=filtop(1:index(filtop,' ')-1),
     + status='old',form='formatted',err=9899)
      rewind(lfntop)
c
      do 101 i=1,5
      read(lfntop,2001,end=9897,err=9898) cdum
 2001 format(a1)
  101 continue
      read(lfntop,2002,end=9897,err=9898) nat,nqt
 2002 format(i5)
      do 102 i=1,nat
      read(lfntop,2001,end=9897,err=9898) cdum
      do 103 j=i,nat
      read(lfntop,2001,end=9897,err=9898) cdum
      read(lfntop,2001,end=9897,err=9898) cdum
  103 continue
  102 continue
      do 104 i=1,nqt
      read(lfntop,2001,end=9897,err=9898) cdum
  104 continue
      read(lfntop,2003,end=9897,err=9898) naw,nbw,nhw,ndw,now,ntw,nnw
 2003 format(5i7,2i10)
      read(lfntop,2001,end=9897,err=9898) cdum
      do 105 i=1,naw
      read(lfntop,2001,end=9897,err=9898) cdum
  105 continue
      do 106 i=1,nbw
      read(lfntop,2001,end=9897,err=9898) cdum
      read(lfntop,2001,end=9897,err=9898) cdum
  106 continue
      do 107 i=1,nhw
      read(lfntop,2001,end=9897,err=9898) cdum
      read(lfntop,2001,end=9897,err=9898) cdum
  107 continue
      do 108 i=1,ndw
      read(lfntop,2001,end=9897,err=9898) cdum
      read(lfntop,2001,end=9897,err=9898) cdum
  108 continue
      do 109 i=1,now
      read(lfntop,2001,end=9897,err=9898) cdum
      read(lfntop,2001,end=9897,err=9898) cdum
  109 continue
      if(ntw.gt.0) then
      read(lfntop,2004,end=9897,err=9898) (idum,i=1,ntw)
      read(lfntop,2004,end=9897,err=9898) (idum,i=1,ntw)
 2004 format(11i7)
      endif
      if(nnw.gt.0) then
      read(lfntop,2005,end=9997,err=9998) (idum,i=1,nnw)
      read(lfntop,2005,end=9997,err=9998) (idum,i=1,nnw)
 2005 format(11i7)
      endif
      read(lfntop,2001,end=9897,err=9898) cdum
c
      do 1 i=1,4
      read(lfnrst,1001,end=9997,err=9998) cdum
 1001 format(a1)
    1 continue
      read(lfnrst,1002,end=9997,err=9998) npbtyp,nbxtyp
 1002 format(i5,i5)
      read(lfnrst,1003,end=9997,err=9998) ((vlat(i,j),j=1,3),i=1,3)
 1003 format(3f12.6)
      read(lfnrst,1004,end=9997,err=9998) jdum
 1004 format(40x,i5)
      read(lfnrst,1005,end=9997,err=9998) temp,tempw,temps
 1005 format(3f12.6)
      do 2 i=1,3
      box(i)=vlat(i,i)
      boxh(i)=half*box(i)
      boxi(i)=one/box(i)
      do 3 j=1,3
      vlati(i,j)=vlat(i,j)
    3 continue
    2 continue
      call matinv(vlati,3,3)
      if(jdum.ne.0) then
      read(lfnrst,1001,end=9997,err=9998) cdum
      endif
      read(lfnrst,1006,end=9997,err=9998) idum
 1006 format(70x,i5)
      if(idum.gt.0) then
      read(lfnrst,1007,end=9997,err=9998) idum,jdum,kdum
 1007 format(3i5)
      read(lfnrst,1008,end=9997,err=9998) (rdum,i=1,idum)
      read(lfnrst,1008,end=9997,err=9998) (rdum,i=1,jdum)
      read(lfnrst,1008,end=9997,err=9998) (rdum,i=1,kdum)
 1008 format(4e20.12)
      endif
c
      if(nwm.gt.0) then
      number=0
      ncyc=nwm/nw+1
      numw=nw
      do 4 icyc=1,ncyc
      if(nwm-number.lt.numw) numw=nwm-number
      do 44 i=1,numw
      read(lfnrst,1009,end=9997,err=9998)
     + ((bxw(i,j,k),j=1,3),(bvw(i,j,k),j=1,3),k=1,nwa)
 1009 format(2x,6f13.8)
      read(lfnrst,1010,end=9997,err=9998) ibw(i,2),(brw(i,k),k=1,3)
 1010 format(i1,1x,3f13.8)
   44 continue
      do 5 i=1,numw
      cgx=zero
      cgy=zero
      cgz=zero
      do 6 k=1,nwa
      cgx=cgx+bxw(i,1,k)
      cgy=cgy+bxw(i,2,k)
      cgz=cgz+bxw(i,3,k)
    6 continue
      ibx=0
      iby=0
      ibz=0
      if(nbxtyp.ne.1) then
      xt(1)=cgx
      xt(2)=cgy
      xt(3)=cgz
      else
      xt(1)=box(1)*(vlati(1,1)*cgx+vlati(1,2)*cgy+vlati(1,3)*cgz)
      xt(2)=box(2)*(vlati(2,1)*cgx+vlati(2,2)*cgy+vlati(2,3)*cgz)
      xt(3)=box(3)*(vlati(3,1)*cgx+vlati(3,2)*cgy+vlati(3,3)*cgz)
      endif
      do 7 j=1,nbx-1
      if(xt(1)/nwa+boxh(1).gt.boxsiz(j,1)) ibx=j
    7 continue
      do 8 j=1,nby-1
      if(xt(2)/nwa+boxh(2).gt.boxsiz(j,2)) iby=j
    8 continue
      do 9 j=1,nbz-1
      if(xt(3)/nwa+boxh(3).gt.boxsiz(j,3)) ibz=j
    9 continue
c
      if(npbtyp.gt.0) then
      if(ibx.ge.nbx) ibx=ibx-nbx
      if(iby.ge.nby) iby=iby-nby
      if(ibx.lt.0) ibx=ibx+nbx
      if(iby.lt.0) iby=iby+nby
      if(npbtyp.eq.1) then
      if(ibz.ge.nbz) ibz=ibz-nbz
      if(ibz.lt.0) ibz=ibz+nbz
      else
      if(ibz.ge.nbz) ibz=nbz-1
      if(ibz.lt.0) ibz=0
      endif
      else
      if(ibx.ge.nbx) ibx=nbx-1
      if(iby.ge.nby) iby=nby-1
      if(ibz.ge.nbz) ibz=nbz-1
      if(ibx.lt.0) ibx=0
      if(iby.lt.0) iby=0
      if(ibz.lt.0) ibz=0
      endif
      ipx=ibownr(ibx+1,1)
      ipy=ibownr(iby+1,2)
      ipz=ibownr(ibz+1,3)
c
      ndx(i)=(ipz*npy+ipy)*npx+ipx
      ibw(i,1)=(ibz*nby+iby)*nbx+ibx
      if(lnode0) then
      ndx(i)=0
      ibw(i,1)=0
      endif
    5 continue
      do 10 node=0,np-1
      new=0
      do 11 i=1,numw
      if(ndx(i).eq.node) then
      new=new+1
      iwl(new,lwgmn)=number+i
      iwl(new,lwnod)=node
      do 124 j=1,3
      do 12 k=1,nwa
      xw(new,j,k)=bxw(i,j,k)
      vw(new,j,k)=bvw(i,j,k)
   12 continue
      xwcr(new,j)=brw(i,j)
  124 continue
      iwl(new,lwbox)=ibw(i,1)
      iwl(new,lwdyn)=5*ibw(i,2)
      endif
   11 continue
c
      if(new.gt.0) then
      call ga_distribution(ga_ip,node,ilp,ihp,jlp,jhp)
      call ga_get(ga_ip,ilp,ilp,jlp,jhp,ipl,mbox)
      nold=ipl(1,2)
      if(nold+new.gt.mwm) call errquit('Dimension mwm too small',0)
      ipl(1,2)=ipl(1,2)+new
      call ga_put(ga_ip,ilp,ilp,jlp,jhp,ipl,mbox)
      call ga_distribution(ga_iw,node,ili,ihi,jli,jhi)
      ili=ili+nold
      ihi=ili+new-1
      call ga_put(ga_iw,ili,ihi,jli,jhi,iwl,mwm)
      call ga_distribution(ga_w,node,ilw,ihw,jlw,jhw)
      ilw=ilw+nold
      ihw=ilw+new-1
      call ga_put(ga_w,ilw,ihw,jlw,jlw+3*mwa-1,xw(1,1,1),mwm)
      call ga_put(ga_w,ilw,ihw,jlw+3*mwa,jlw+6*mwa-1,vw(1,1,1),mwm)
      call ga_put(ga_w,ilw,ihw,jlw+6*mwa,jlw+6*mwa+2,xwcr(1,1),mwm)
      endif
c
   10 continue
      number=number+numw
    4 continue
      endif
c
      if(nsa.gt.0) then
      ipp=-1
      ips=-1
      number=0
      ncyc=nsa/ns+1
      nums=ns
      do 13 i=1,ncyc
      if(nsa-number.lt.nums) nums=nsa-number
      read(lfnrst,1011,end=9997,err=9998)
     + (ibs(k,11),(bxs(k,j),j=1,3),(bvs(k,j),j=1,3),k=1,nums)
 1011 format(i1,1x,6f13.8)
      read(lfntop,2009,end=9897,err=9898) ((ibs(k,j),j=1,10),k=1,nums)
 2009 format(16x,10i5)
      do 14 k=1,nums
      if(nsf.lt.ibs(k,1)) nsf=ibs(k,1)
      ibx=0
      iby=0
      ibz=0
      do 15 ix=1,3
      xt(ix)=bxs(k,ix)
      if(npbtyp.ne.0) then
      if(abs(xt(ix)).gt.boxh(ix)) then
#if defined(CRAY_T3D) || defined(CRAY_T3E)
      xt(ix)=xt(ix)-arg_nint(xt(ix)*boxi(ix))*box(ix)
#else
      xt(ix)=xt(ix)-nint(xt(ix)*boxi(ix))*box(ix)
#endif
      endif
      endif
   15 continue
      do 16 j=1,nbx-1
      if(xt(1)+boxh(1).gt.boxsiz(j,1)) ibx=j
   16 continue
      do 17 j=1,nby-1
      if(xt(2)+boxh(2).gt.boxsiz(j,2)) iby=j
   17 continue
      do 18 j=1,nbz-1
      if(xt(3)+boxh(3).gt.boxsiz(j,3)) ibz=j
   18 continue
c
      if(npbtyp.gt.0) then
      if(ibx.ge.nbx) ibx=ibx-nbx
      if(iby.ge.nby) iby=iby-nby
      if(ibx.lt.0) ibx=ibx+nbx
      if(iby.lt.0) iby=iby+nby
      if(npbtyp.eq.1) then
      if(ibz.ge.nbz) ibz=ibz-nbz
      if(ibz.lt.0) ibz=ibz+nbz
      else
      if(ibz.ge.nbz) ibz=nbz-1
      if(ibz.lt.0) ibz=0
      endif
      else
      if(ibx.ge.nbx) ibx=nbx-1
      if(iby.ge.nby) iby=nby-1
      if(ibz.ge.nbz) ibz=nbz-1
      if(ibx.lt.0) ibx=0
      if(iby.lt.0) iby=0
      if(ibz.lt.0) ibz=0
      endif
      ipx=ibownr(ibx+1,1)
      ipy=ibownr(iby+1,2)
      ipz=ibownr(ibz+1,3)
      if(ibs(k,3).ne.ips) then
      ips=ibs(k,3)
      ipp=(ibz*nby+iby)*nbx+ibx
      endif
      ndx(k)=ipp
      if(lnode0) ndx(k)=0
   14 continue
      do 19 node=0,np-1
      new=0
      do 20 k=1,nums
      if(sp_btop(ndx(k),ibownr).eq.node) then
      new=new+1
      isl(new,lsgan)=number+k
      isl(new,lsfrc)=ibs(k,1)
      isl(new,lsmol)=ibs(k,2)
      isl(new,lssgm)=ibs(k,3)
      isl(new,lsgrp)=ibs(k,4)
      isl(new,lspgr)=ibs(k,5)
      isl(new,lsatt)=ibs(k,6)
      isl(new,lsct1)=ibs(k,7)
      isl(new,lsct2)=ibs(k,8)
      isl(new,lsct3)=ibs(k,9)
      isl(new,lssss)=ibs(k,10)
      isl(new,lsdyn)=5*ibs(k,11)
      xs(new,1)=bxs(k,1)
      xs(new,2)=bxs(k,2)
      xs(new,3)=bxs(k,3)
      vs(new,1)=bvs(k,1)
      vs(new,2)=bvs(k,2)
      vs(new,3)=bvs(k,3)
      endif
   20 continue
c
      if(new.gt.0) then
      call ga_distribution(ga_ip,node,ilp,ihp,jlp,jhp)
      call ga_get(ga_ip,ilp,ilp+1,jlp,jhp,ipl,mbox)
      nold=ipl(2,2)
      if(nold+new.gt.msa)
     + call errquit('Dimension msa too small',nold+new)
      ipl(2,2)=ipl(2,2)+new
      call ga_put(ga_ip,ilp,ilp+1,jlp,jhp,ipl,mbox)
      call ga_distribution(ga_is,node,ili,ihi,jli,jhi)
      ili=ili+nold
      ihi=ili+new-1
      call ga_put(ga_is,ili,ihi,jli,jhi,isl,msa)
      call ga_distribution(ga_s,node,ils,ihs,jls,jhs)
      ils=ils+nold
      ihs=ils+new-1
      call ga_put(ga_s,ils,ihs,jls,jls+2,xs(1,1),msa)
      call ga_put(ga_s,ils,ihs,jls+3,jls+5,vs(1,1),msa)
      endif
c
   19 continue
      number=number+nums
   13 continue
      endif
c
      close(unit=lfnrst,status='keep')
      close(unit=lfntop,status='keep')
c
      endif
c
      return
 9897 call errquit('EOF encountered on topology file',0)
 9898 call errquit('Error reading topology file',0)
 9899 call errquit('Error opening topology file',0)
 9997 call errquit('EOF encountered on restart file',0)
 9998 call errquit('Error reading restart file',0)
 9999 call errquit('Error opening restart file',0)
      return
      end
      integer function sp_btop(ibox,ibownr)
c
      implicit none
c
#include "sp_common.fh"
c
      integer ibox,ibownr(maxbox,3)
      integer iboxx,iboxy,iboxz
c
      iboxx=mod(ibox,nbx)
      iboxy=mod((ibox-iboxx)/nbx,nby)
      iboxz=((ibox-iboxx)/nbx-iboxy)/nby
      sp_btop=(ibownr(iboxz+1,3)*npy+ibownr(iboxy+1,2))*npx+
     + ibownr(iboxx+1,1)
c
      return
      end
      subroutine sp_decomp(ibownr,boxsiz)
c
      implicit none
c
#include "sp_common.fh"
#include "msgids.fh"
#include "global.fh"
c
      integer ibownr(maxbox,3)
      real*8 boxsiz(maxbox,3)
c
      integer ibx,iby,ibz,i,ix,jx,iy,jy,iz,jz,jnode,nrnod,nbiown
c
c     check dimensions of ibownr
c
      if(nbx.gt.maxbox.or.nby.gt.maxbox.or.nbz.gt.maxbox)
     + call errquit('Dimension maxbox too small',0)
c
c     determine the node dimension for each sub box
c
      do 1 ibx=1,nbx
      ibownr(ibx,1)=((ibx-1)*npx)/nbx
      boxsiz(ibx,1)=(box(1)*dble(ibx))/dble(nbx)
    1 continue
      do 2 iby=1,nby
      ibownr(iby,2)=((iby-1)*npy)/nby
      boxsiz(iby,2)=(box(2)*dble(iby))/dble(nby)
    2 continue
      do 3 ibz=1,nbz
      ibownr(ibz,3)=((ibz-1)*npz)/nbz
      boxsiz(ibz,3)=(box(3)*dble(ibz))/dble(nbz)
    3 continue
c
c     determine neighboring nodes and store in neighb(27,2)
c     such that:
c
c     neighb(n,1) is the n-th neighbor of node me
c     neighb(n,2) is the node of which node me is the n-th neighbor
c
c     a value of -1 indicates that such node does not exist
c
      do 10 i=1,27
      neighb(i,1)=-1
      neighb(i,2)=-1
 10   continue
c
      do 11 ix=1,3
      jx=mex+ix-2
      if(npbtyp.gt.0) then
      if(npx.gt.2.and.jx.lt.0) jx=jx+npx
      if(npx.gt.2.and.jx.ge.npx) jx=jx-npx
      endif
      if(jx.ge.0.and.jx.lt.npx) then
      do 12 iy=1,3
      jy=mey+iy-2
      if(npbtyp.gt.0) then
      if(npy.gt.2.and.jy.lt.0) jy=jy+npy
      if(npy.gt.2.and.jy.ge.npy) jy=jy-npy
      endif
      if(jy.ge.0.and.jy.lt.npy) then
      do 13 iz=1,3
      jz=mez+iz-2
      if(npbtyp.eq.1) then
      if(npz.gt.2.and.jz.lt.0) jz=jz+npz
      if(npz.gt.2.and.jz.ge.npz) jz=jz-npz
      endif
      if(jz.ge.0.and.jz.lt.npz) then
      jnode=((jz*npy)+jy)*npx+jx
      neighb(3*(3*(ix-1)+(iy-1))+iz,1)=jnode
      neighb(3*(3*(3-ix)+(3-iy))+4-iz,2)=jnode
      endif
   13 continue
      endif
   12 continue
      endif
   11 continue
c
      nbiown=30
      do 20 ibx=1,nbx
      do 21 iby=1,nby
      do 22 ibz=1,nbz
      nrnod=(ibownr(ibz,3)*npy+ibownr(iby,2))*npx+ibownr(ibx,1)
      if(me.eq.nrnod) nbiown=nbiown+1
   22 continue
   21 continue
   20 continue
c
c     broadcast the maximum number of sub-boxes per node to all nodes
c
      if(nbiown.gt.mbox) call errquit('Error in mbox',0)
c
      return
      end
      subroutine sp_nrnode
c
      implicit none
c
#include "sp_common.fh"
c
      integer ix,iy,iz,npt,i,j,k,l
c
c     this routine distributes the available processes in the
c     Cartesian directions
c
c     npx : number of processes in x direction
c     npy : number of processes in y direction
c     npz : number of processes in z direction
c
c     determine node dimensions
c
      if(npx*npy*npz.ne.np) then
      if(npx+npy+npz.gt.0) 
     + call errquit('Specified npx*npy*npz ne np',0)
c
      npt=0
      do 1 i=1,np
      do 2 j=i,np
      do 3 k=j,np
      l=i*j*k
      if(l.eq.np) then
      if(l.eq.npt) then
      if(k.gt.npz) goto 3
      if(i+j+k.lt.npx+npy+npz) then
      npx=i
      npy=j
      npz=k
      endif
      goto 3
      else
      npt=np
      npx=i
      npy=j
      npz=k
      endif
      endif
    3 continue
    2 continue
    1 continue
c
      if(npx*npy*npz.ne.np) call errquit('nrnode: code error',0)
      endif
c
c     determine processor location of me
c
      do 4 ix=1,npx
      do 5 iy=1,npy
      do 6 iz=1,npz
      if(me.eq.((iz-1)*npy+(iy-1))*npx+(ix-1)) then
      mex=ix-1
      mey=iy-1
      mez=iz-1
      endif
    6 continue
    5 continue
    4 continue
c
      return
      end
      subroutine sp_initip(ibownr,ipl)
c
      implicit none
c
#include "sp_common.fh"
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
      nbiown=0
      do 3 ibx=1,nbx
      do 4 iby=1,nby
      do 5 ibz=1,nbz
      nrbox=((ibz-1)*nby+iby-1)*nbx+ibx-1
      nrnod=(ibownr(ibz,3)*npy+ibownr(iby,2))*npx+ibownr(ibx,1)
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
      return
      end
      subroutine sp_numbb(ibownr,boxsiz)
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "global.fh"
c
      integer ibownr(maxbox,3)
      real*8 boxsiz(maxbox,3)
      logical lside,leven
c
      integer ibx,iby,ibz,ipx,ipy,ipz,ibox,inode
      integer jbx,jby,jbz,kbx,kby,kbz,jpx,jpy,jpz
      integer ilx,ihx,ily,ihy,ilz,ihz
      integer jbox,jnode,mbblb,mbblp
      real*8 dx,dxtmp,dy,dytmp,dz,dztmp,dist2
      character*255 string
c
c     determine size of box-box pairlist
c     ----------------------------------
c
      mbblp=mbbl
      mbbl=0
      mbblb=0
      do 1 ibx=0,nbx-1
      ipx=ibownr(ibx+1,1)
      do 2 iby=0,nby-1
      ipy=ibownr(iby+1,2)
      do 3 ibz=0,nbz-1
      ipz=ibownr(ibz+1,3)
      ibox=(ibz*nby+iby)*nbx+ibx
      inode=(ipz*npy+ipy)*npx+ipx
      if(inode.eq.me) then
      do 4 jbx=0,nbx-1
      kbx=jbx-ibx
      jpx=ibownr(jbx+1,1)
c
      if(ibx.le.jbx) then
      ilx=ibx
      ihx=jbx
      else
      ilx=jbx
      ihx=ibx
      endif
c
      dx=zero
      if(ibx.ne.jbx) then
      dx=boxsiz(ihx,1)-boxsiz(ilx+1,1)
      if(npbtyp.gt.0) then
      dxtmp=zero
      if(ilx.gt.0) dxtmp=boxsiz(ilx,1)
      if(ihx.lt.nbx-1) dxtmp=dxtmp-boxsiz(ihx+1,1)+box(1)
      if(dxtmp.lt.dx) dx=dxtmp
      if(kbx.gt.0.and.kbx.gt.iabs(kbx-nbx)) kbx=kbx-nbx
      if(kbx.lt.0.and.-kbx.gt.iabs(kbx+nbx)) kbx=kbx+nbx
      endif
      endif
c
      do 5 jby=0,nby-1
      kby=jby-iby
      jpy=ibownr(jby+1,2)
c
      if(iby.le.jby) then
      ily=iby
      ihy=jby
      else
      ily=jby
      ihy=iby
      endif
c
      dy=zero
      if(iby.ne.jby) then
      dy=boxsiz(ihy,2)-boxsiz(ily+1,2)
      if(npbtyp.gt.0) then
      dytmp=zero
      if(ily.gt.0) dytmp=boxsiz(ily,2)
      if(ihy.lt.nby-1) dytmp=dytmp-boxsiz(ihy+1,2)+box(2)
      if(dytmp.lt.dy) dy=dytmp
      if(kby.gt.0.and.kby.gt.iabs(kby-nby)) kby=kby-nby
      if(kby.lt.0.and.-kby.gt.iabs(kby+nby)) kby=kby+nby
      endif
      endif
c
      do 6 jbz=0,nbz-1
      kbz=jbz-ibz
      jpz=ibownr(jbz+1,3)
c
      if(ibz.le.jbz) then
      ilz=ibz
      ihz=jbz
      else
      ilz=jbz
      ihz=ibz
      endif
c
      dz=zero
      if(ibz.ne.jbz) then
      dz=boxsiz(ihz,3)-boxsiz(ilz+1,3)
      if(npbtyp.eq.1) then
      dztmp=zero
      if(ilz.gt.0) dztmp=boxsiz(ilz,3)
      if(ihz.lt.nbz-1) dztmp=dztmp-boxsiz(ihz+1,3)+box(3)
      if(dztmp.lt.dz) dz=dztmp
      if(kbz.gt.0.and.kbz.gt.iabs(kbz-nbz)) kbz=kbz-nbz
      if(kbz.lt.0.and.-kbz.gt.iabs(kbz+nbz)) kbz=kbz+nbz
      endif
      endif
c
      jbox=(jbz*nby+jby)*nbx+jbx
      jnode=(jpz*npy+jpy)*npx+jpx
c
c     determine orientation jbox in relation to ibox
c
c     lside is true if
c
c     i: 0  j: 0  k:  +
c     i: 0  j:  + k:-0+
c     i:  + j:-0+ k:-0+
c
      lside=(kbx.eq.0.and.kby.eq.0.and.kbz.ge.0)
     + .or.(kbx.eq.0.and.kby.gt.0) .or. kbx.gt.0
c
c     determine if ibox is identical to jbox
c
c     lsame=kbx.eq.0.and.kby.eq.0.and.kbz.eq.0
c
c     determine if difference in box numbers is even or odd
c
      leven=2*(iabs(ibox-jbox)/2).eq.iabs(ibox-jbox)
c
c     calculate the distance between the two boxes
c
      if(nbxtyp.eq.1) then
      dist2=
     + (vlat(1,1)*dx/box(1)+vlat(1,2)*dy/box(2)+vlat(1,3)*dz/box(3))**2+
     + (vlat(2,1)*dx/box(1)+vlat(2,2)*dy/box(2)+vlat(2,3)*dz/box(3))**2+
     + (vlat(3,1)*dx/box(1)+vlat(3,2)*dy/box(2)+vlat(3,3)*dz/box(3))**2
      else
      dist2=dx*dx+dy*dy+dz*dz
      endif
c
c     keep half of the box pairs
c
c     this test also appears in sp_numbb
c     any changes need to be made in both routines
c
      if((inode.eq.jnode.and.ibox.ge.jbox).or.(inode.ne.jnode.and.
     + ((lside.and.leven).or.(.not.lside.and..not.leven)))) then
c
c     keep only those within maximum cutoff distance
c
      if(rlong*rlong.gt.dist2) mbbl=mbbl+1
      endif
      if(rlong*rlong.gt.dist2) mbblb=mbblb+1
 6    continue
 5    continue
 4    continue
      endif
 3    continue
 2    continue
 1    continue
c
      mbbl=mbbl+1
      if(mbblb+1.gt.mbbl) mbbl=mbblb+1
c
c     allocate memory for the box-box list
c
      if(mbblp.eq.0) then
      if(np.gt.1) call ga_igop(msp_20,mbbl,1,'max')
      if(.not.ma_push_get(mt_int,15*mbbl,'bb',l_bb,i_bb))
     + call errquit('Failed to allocate memory for bb',0)
      else
      if(mbbl.gt.mbblp) then
      write(string,1111) mbblp,mbbl
 1111 format('error: lbbl increase from ',i6,' to ',i6)
      call errquit(string,me)
c      call errquit('lbbl increased beyond allocated memory',me)
      endif
      endif
c
      return
      end
      subroutine sp_listbb(ibownr,boxsiz,lbbl)
c
      implicit none
c
#include "sp_common.fh"
c
      integer ibownr(maxbox,3),lbbl(mbbl,mbb2)
      real*8 boxsiz(maxbox,3)
      logical lside,leven
c
      integer ibx,iby,ibz,ipx,ipy,ipz,ibox,inode
      integer jbx,jby,jbz,kbx,kby,kbz,jpx,jpy,jpz
      integer ilx,ihx,ily,ihy,ilz,ihz
      integer i,j,jbox,jnode,ltemp
      real*8 dx,dxtmp,dy,dytmp,dz,dztmp,dist2
c
c     Construction of the box-box pairlist
c
      nbbl=0
      do 1 ibx=0,nbx-1
      ipx=ibownr(ibx+1,1)
      do 2 iby=0,nby-1
      ipy=ibownr(iby+1,2)
      do 3 ibz=0,nbz-1
      ipz=ibownr(ibz+1,3)
      ibox=(ibz*nby+iby)*nbx+ibx
      inode=(ipz*npy+ipy)*npx+ipx
      if(inode.eq.me) then
      do 4 jbx=0,nbx-1
      kbx=jbx-ibx
      jpx=ibownr(jbx+1,1)
c
      if(ibx.le.jbx) then
      ilx=ibx
      ihx=jbx
      else
      ilx=jbx
      ihx=ibx
      endif
c
      dx=zero
      if(ibx.ne.jbx) then
      dx=boxsiz(ihx,1)-boxsiz(ilx+1,1)
      if(npbtyp.gt.0) then
      dxtmp=zero
      if(ilx.gt.0) dxtmp=boxsiz(ilx,1)
      if(ihx.lt.nbx-1) dxtmp=dxtmp-boxsiz(ihx+1,1)+box(1)
      if(dxtmp.lt.dx) dx=dxtmp
      if(kbx.gt.0.and.kbx.gt.iabs(kbx-nbx)) kbx=kbx-nbx
      if(kbx.lt.0.and.-kbx.gt.iabs(kbx+nbx)) kbx=kbx+nbx
      endif
      endif
c
      do 5 jby=0,nby-1
      kby=jby-iby
      jpy=ibownr(jby+1,2)
c
      if(iby.le.jby) then
      ily=iby
      ihy=jby
      else
      ily=jby
      ihy=iby
      endif
c
      dy=zero
      if(iby.ne.jby) then
      dy=boxsiz(ihy,2)-boxsiz(ily+1,2)
      if(npbtyp.gt.0) then
      dytmp=zero
      if(ily.gt.0) dytmp=boxsiz(ily,2)
      if(ihy.lt.nby-1) dytmp=dytmp-boxsiz(ihy+1,2)+box(2)
      if(dytmp.lt.dy) dy=dytmp
      if(kby.gt.0.and.kby.gt.iabs(kby-nby)) kby=kby-nby
      if(kby.lt.0.and.-kby.gt.iabs(kby+nby)) kby=kby+nby
      endif
      endif
c
      do 6 jbz=0,nbz-1
      kbz=jbz-ibz
      jpz=ibownr(jbz+1,3)
c
      if(ibz.le.jbz) then
      ilz=ibz
      ihz=jbz
      else
      ilz=jbz
      ihz=ibz
      endif
c
      dz=zero
      if(ibz.ne.jbz) then
      dz=boxsiz(ihz,3)-boxsiz(ilz+1,3)
      if(npbtyp.eq.1) then
      dztmp=zero
      if(ilz.gt.0) dztmp=boxsiz(ilz,3)
      if(ihz.lt.nbz-1) dztmp=dztmp-boxsiz(ihz+1,3)+box(3)
      if(dztmp.lt.dz) dz=dztmp
      if(kbz.gt.0.and.kbz.gt.iabs(kbz-nbz)) kbz=kbz-nbz
      if(kbz.lt.0.and.-kbz.gt.iabs(kbz+nbz)) kbz=kbz+nbz
      endif
      endif
c
      jbox=(jbz*nby+jby)*nbx+jbx
      jnode=(jpz*npy+jpy)*npx+jpx
c
c     determine orientation jbox in relation to ibox
c
c     lside is true if
c
c     i: 0  j: 0  k:  +
c     i: 0  j:  + k:-0+
c     i:  + j:-0+ k:-0+
c
      lside=(kbx.eq.0.and.kby.eq.0.and.kbz.ge.0)
     + .or.(kbx.eq.0.and.kby.gt.0) .or. kbx.gt.0
c
c     determine if ibox is identical to jbox
c
c     lsame=kbx.eq.0.and.kby.eq.0.and.kbz.eq.0
c
c     determine if difference in box numbers is even or odd
c
      leven=2*(iabs(ibox-jbox)/2).eq.iabs(ibox-jbox)
c
c     calculate the distance between the two boxes
c
      if(nbxtyp.eq.1) then
      dist2=
     + (vlat(1,1)*dx/box(1)+vlat(1,2)*dy/box(2)+vlat(1,3)*dz/box(3))**2+
     + (vlat(2,1)*dx/box(1)+vlat(2,2)*dy/box(2)+vlat(2,3)*dz/box(3))**2+
     + (vlat(3,1)*dx/box(1)+vlat(3,2)*dy/box(2)+vlat(3,3)*dz/box(3))**2
      else
      dist2=dx*dx+dy*dy+dz*dz
      endif
c
c     keep half of the box pairs
c
c     this test also appears in sp_numbb
c     any changes need to be made in both routines
c
      if((inode.eq.jnode.and.ibox.ge.jbox).or. (inode.ne.jnode.and.
     + ((lside.and.leven).or.(.not.lside.and..not.leven)))) then
c
c     keep only those within maximum cutoff distance
c
      if(rlong*rlong.gt.dist2) then
      nbbl=nbbl+1
      if(nbbl.gt.mbbl) call errquit('Box-box list too small',mbbl)
      lbbl(nbbl,1)=jnode
      lbbl(nbbl,2)=jbox
      lbbl(nbbl,3)=ibox
      lbbl(nbbl,4)=0
      endif
      endif
 6    continue
 5    continue
 4    continue
      endif
 3    continue
 2    continue
 1    continue
      npprev=0
c
      nbbloc=0
      do 7 i=1,nbbl-1
      do 8 j=i+1,nbbl
      if((lbbl(i,1).ne.me.and.lbbl(j,1).eq.me).or.
     + (lbbl(i,1).gt.lbbl(j,1).and.lbbl(i,1).ne.me).or.
     + (lbbl(i,1).eq.lbbl(j,1).and.lbbl(i,2).gt.lbbl(j,2)).or.
     + (lbbl(i,1).eq.lbbl(j,1).and.lbbl(i,2).eq.lbbl(j,2).and.
     + lbbl(i,3).gt.lbbl(j,3))) then
      ltemp=lbbl(i,1)
      lbbl(i,1)=lbbl(j,1)
      lbbl(j,1)=ltemp
      ltemp=lbbl(i,2)
      lbbl(i,2)=lbbl(j,2)
      lbbl(j,2)=ltemp
      ltemp=lbbl(i,3)
      lbbl(i,3)=lbbl(j,3)
      lbbl(j,3)=ltemp
      endif
    8 continue
      if(lbbl(i,1).eq.me) nbbloc=i
    7 continue
      if(lbbl(nbbl,1).eq.me) nbbloc=nbbl
c
      return
      end

      subroutine sp_finish()
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
c
      if(.not.ma_pop_stack(l_bb))
     + call errquit('Failed to deallocate memory for bb',0)
c
      call sp_free()
c
      return
      end
      subroutine sp_initf(fw,fs,llng,iwz,isz,lpair)
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      real*8 fw(mwm,3,mwa,2),fs(msa,3,2)
      integer iwz(mwm),isz(msa)
      logical llng,lpair
c
      integer i,j,k,l,m,il,ih,jl,jh
c
      llong=llng
c
      m=1
      if(llong) m=2
c
      do 1 l=1,m
      if(nwm.gt.0) then
      do 2 k=1,mwa
      do 3 j=1,3
      do 4 i=1,mwm
      fw(i,j,k,l)=zero
    4 continue
    3 continue
    2 continue
      endif
      if(nsa.gt.0) then
      do 5 j=1,3
      do 6 i=1,msa
      fs(i,j,l)=zero
    6 continue
    5 continue
      endif
    1 continue
c
      if(nwm.gt.0) then
      call ga_distribution(ga_w,me,il,ih,jl,jh)
      call ga_put(ga_w,il,ih,jl+6*mwa+3,jl+9*mwa+2,fw,mwm)
      if(llong) call ga_put(ga_w,il,ih,jl+9*mwa+3,jl+12*mwa+2,
     + fw(1,1,1,2),mwm)
      endif
      if(nsa.gt.0) then
      call ga_distribution(ga_s,me,il,ih,jl,jh)
      call ga_put(ga_s,il,ih,jl+6,jl+8,fs,msa)
      if(llong) call ga_put(ga_s,il,ih,jl+9,jl+11,fs(1,1,2),msa) 
      endif
c
      if(lpair) then
      do 7 i=1,mwm
      iwz(i)=0
    7 continue
      do 8 i=1,msa
      isz(i)=0
    8 continue
      call ga_zero(ga_iwz)
      call ga_zero(ga_isz)
      endif
c
      return
      end
      subroutine sp_final(fw,fs,lpair,iwz,isz)
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      real*8 fw(mwm,3,mwa,2),fs(msa,3,2)
      logical lpair
      integer iwz(mwm),isz(msa)
c
      integer i,il,ih,jl,jh
c
      if(np.gt.0) then
      if(nwm.gt.0) then
      call ga_distribution(ga_w,me,il,ih,jl,jh)
      call ga_acc(ga_w,il,ih,jl+6*mwa+3,jl+9*mwa+2,fw,mwm,one)
      if(llong) call ga_acc(ga_w,il,ih,jl+9*mwa+3,jl+12*mwa+2,
     + fw(1,1,1,2),mwm,one)
      call ga_get(ga_w,il,ih,jl+6*mwa+3,jl+9*mwa+2,fw,mwm)
      if(ltwin) call ga_get(ga_w,il,ih,jl+9*mwa+3,jl+12*mwa+2,
     + fw(1,1,1,2),mwm)
      endif
      if(nsa.gt.0) then
      call ga_distribution(ga_s,me,il,ih,jl,jh)
      call ga_acc(ga_s,il,ih,jl+6,jl+8,fs,msa,one)
      if(llong) call ga_acc(ga_s,il,ih,jl+9,jl+11,fs(1,1,2),msa,one)
      call ga_get(ga_s,il,ih,jl+6,jl+8,fs,msa)
      if(ltwin) call ga_get(ga_s,il,ih,jl+9,jl+11,fs(1,1,2),msa)
      endif
      endif
c
      if(lpair) then
      if(nwm.gt.0) then
      call ga_distribution(ga_iwz,me,il,ih,jl,jh)
      call ga_acc(ga_iwz,il,ih,1,1,iwz,mwm,1)
      call ga_get(ga_iwz,il,ih,1,1,iwz,mwm)
      do 1 i=1,nwmloc
      iwz(i)=min(1,iwz(i))
    1 continue
      endif
      if(nsa.gt.0) then
      call ga_distribution(ga_isz,me,il,ih,jl,jh)
      call ga_acc(ga_isz,il,ih,1,1,isz,msa,1)
      call ga_get(ga_isz,il,ih,1,1,isz,msa)
      do 2 i=1,nsaloc
      isz(i)=min(1,isz(i))
    2 continue
      endif
      endif
c
      return
      end
      subroutine sp_wrtrst(lfnrst,filrst,lveloc,
     + pres,temp,tempw,temps,iwl,xw,vw,xwcr,isl,xs,vs)
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      integer lfnrst
      character*255 filrst
      logical lveloc
      integer iwl(mwm,miw2),isl(msa,mis2)
      real*8 pres,temp,tempw,temps
      real*8 xw(mwm,3,mwa),vw(mwm,3,mwa),xwcr(mwm,3)
      real*8 xs(msa,3),vs(msa,3)
c
      integer lenscr
c
      lenscr=ma_inquire_avail(mt_byte)/
     + ((6*mwa+3)*ma_sizeof(mt_dbl,1,mt_byte)+
     + (mis2+4)*ma_sizeof(mt_int,1,mt_byte))-1
      if(.not.ma_push_get(mt_dbl,lenscr*3*mwa,'bx',l_bx,i_bx))
     + call errquit('Failed to allocate bx',0)
      if(.not.ma_push_get(mt_dbl,lenscr*3*mwa,'bv',l_bv,i_bv))
     + call errquit('Failed to allocate bv',0)
      if(.not.ma_push_get(mt_dbl,lenscr*3,'br',l_br,i_br))
     + call errquit('Failed to allocate br',0)
      if(.not.ma_push_get(mt_int,lenscr*max(mis2,2),'bi',l_bi,i_bi))
     + call errquit('Failed to allocate bi',0)
      if(.not.ma_push_get(mt_int,lenscr,'n',l_n,i_n))
     + call errquit('Failed to allocate n',0)
c
      call sp_wtrst(lfnrst,filrst,lveloc,pres,temp,tempw,temps,
     + iwl,xw,vw,xwcr,isl,xs,vs,int_mb(i_ipl),
     + lenscr,int_mb(i_bi),dbl_mb(i_bx),dbl_mb(i_bv),dbl_mb(i_br),
     + int_mb(i_bi),dbl_mb(i_bx),dbl_mb(i_bv))
c
      if(.not.ma_pop_stack(l_n))
     + call errquit('Failed to deallocate n',0)
      if(.not.ma_pop_stack(l_bi))
     + call errquit('Failed to deallocate bi',0)
      if(.not.ma_pop_stack(l_br))
     + call errquit('Failed to deallocate br',0)
      if(.not.ma_pop_stack(l_bv))
     + call errquit('Failed to deallocate bv',0)
      if(.not.ma_pop_stack(l_bx))
     + call errquit('Failed to deallocate bx',0)
c
      return
      end
      subroutine sp_wtrst(lfnrst,filrst,lveloc,pres,temp,tempw,temps,
     + iwl,xw,vw,xwcr,isl,xs,vs,
     + ipl,nb,ibw,bxw,bvw,brw,ibs,bxs,bvs)
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      integer lfnrst,nb
      character*255 filrst
      logical lveloc
      real*8 pres,temp,tempw,temps
      integer iwl(mwm,miw2),isl(msa,mis2)
      real*8 xw(mwm,3,mwa),vw(mwm,3,mwa),xwcr(mwm,3)
      real*8 xs(msa,3),vs(msa,3)
      integer ipl(mbox,mip2),ibw(nb),ibs(nb)
      real*8 bxw(nb,3,mwa),bvw(nb,3,mwa),brw(nb,3)
      real*8 bxs(nb,3),bvs(nb,3)
c
      integer i,j,k,node,ncyc,icyc,numw,nums,number,nwmn,nsan
      integer ilp,ihp,jlp,jhp,ili,ihi,jli,jhi,ilw,ihw,jlw,jhw
      integer ils,ihs,jls,jhs
      character*10 rdate,rtime
      integer idyn
c
      if(me.eq.0) then
c
      call swatch(rdate,rtime)
c
      rewind(lfnrst)
      write(lfnrst,1000)
 1000 format('Restart file',/,' ',/,' ')
      write(lfnrst,1001) 3.3,rdate,rtime
 1001 format(f12.6,2a10)
      write(lfnrst,1002) npbtyp,nbxtyp,((vlat(i,j),j=1,3),i=1,3)
 1002 format(2i5,/,(3f12.6))
      write(lfnrst,1003) pres
 1003 format(1pe12.5)
      write(lfnrst,1004) temp,tempw,temps
 1004 format(3f12.6)
      write(lfnrst,1005) nwm,nwa,nsm,nsa,nwmc,0,0,0,0
 1005 format(7i10,2i5)
c
      if(nwm.gt.0) then
      number=0
      ncyc=nwm/nb+1
      numw=nb
      do 1 icyc=1,ncyc
      if(nwm-number.lt.numw) numw=nwm-number
      do 2 node=np-1,0,-1
      call ga_distribution(ga_ip,node,ilp,ihp,jlp,jhp)
      call ga_get(ga_ip,ilp,ihp,jlp,jhp,ipl,mbox)
      nwmn=ipl(1,2)
      if(nwmn.gt.0) then
      call ga_distribution(ga_iw,node,ili,ihi,jli,jhi)
      call ga_get(ga_iw,ili,ili+nwmn-1,jli,jli+lwdyn-1,iwl,mwm)
      call ga_distribution(ga_w,node,ilw,ihw,jlw,jhw)
      call ga_get(ga_w,ilw,ilw+nwmn-1,jlw,jlw+3*mwa-1,xw,mwm)
      if(lveloc)
     + call ga_get(ga_w,ilw,ilw+nwmn-1,jlw+3*mwa,jlw+6*mwa-1,vw,mwm)
      call ga_get(ga_w,ilw,ilw+nwmn-1,jlw+6*mwa,jlw+6*mwa+2,xwcr,mwm)
      do 3 i=1,nwmn
      j=iwl(i,lwgmn)-number
      if(j.gt.0.and.j.le.numw) then
      do 4 k=1,nwa
      bxw(j,1,k)=xw(i,1,k)
      bxw(j,2,k)=xw(i,2,k)
      bxw(j,3,k)=xw(i,3,k)
      bvw(j,1,k)=vw(i,1,k)
      bvw(j,2,k)=vw(i,2,k)
      bvw(j,3,k)=vw(i,3,k)
    4 continue
      brw(j,1)=xwcr(i,1)
      brw(j,2)=xwcr(i,2)
      brw(j,3)=xwcr(i,3)
      ibw(j)=iwl(i,lwdyn)
      endif
    3 continue
      endif
    2 continue
      do 5 i=1,numw
      if(lveloc) then
      write(lfnrst,1006) ((bxw(i,j,k),j=1,3),(bvw(i,j,k),j=1,3),k=1,nwa)
      else
      write(lfnrst,1006) ((bxw(i,j,k),j=1,3),(zero,j=1,3),k=1,nwa)
      endif
 1006 format(2x,6f13.8)
      idyn=iand(ibw(i),3)
      write(lfnrst,1007) idyn,(brw(i,k),k=1,3)
 1007 format(i1,1x,3f13.8)
    5 continue
    1 continue
      endif
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
      call ga_get(ga_is,ili,ili+nsan-1,jli,jli+lsdyn-1,isl,msa)
      call ga_distribution(ga_s,node,ils,ihs,jls,jhs)
      call ga_get(ga_s,ils,ils+nsan-1,jls,jls+2,xs,msa)
      if(lveloc) call ga_get(ga_s,ils,ils+nsan-1,jls+3,jls+5,vs,msa)
      do 8 i=1,nsan
      j=isl(i,lsgan)-number
      if(j.gt.0.and.j.le.nums) then
      bxs(j,1)=xs(i,1)
      bxs(j,2)=xs(i,2)
      bxs(j,3)=xs(i,3)
      bvs(j,1)=vs(i,1)
      bvs(j,2)=vs(i,2)
      bvs(j,3)=vs(i,3)
      ibs(j)=isl(i,lsdyn)
      endif
    8 continue
      endif
    7 continue
      do 9 i=1,nums
      idyn=iand(ibs(i),3)
      if(lveloc) then
      write(lfnrst,1008) idyn,(bxs(i,j),j=1,3),(bvs(i,j),j=1,3)
      else
      write(lfnrst,1008) idyn,(bxs(i,j),j=1,3),(zero,j=1,3)
      endif
 1008 format(i1,1x,6f13.8)
    9 continue
    6 continue
      endif
c
      endif
c
      return
 9999 continue
      call errquit('Failed to open restart file',me)
      return
      end
      subroutine sp_print(lfnout)
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
c
      integer lfnout
c
      integer i_lcnt,l_lcnt
c
      if(me.eq.0) then
      write(lfnout,1000)
 1000 format(/,' DOMAIN DECOMPOSITION',/)
c
      write(lfnout,1001) np,npx,npy,npz
 1001 format(' Processor count ',i5,' =',i5,' x',i5,' x',i5)
      write(lfnout,1002) nbx*nby*nbz,nbx,nby,nbz
 1002 format(' Box count       ',i5,' =',i5,' x',i5,' x',i5)
c
      if(mod(nbx,npx)+mod(nby,npy)+mod(nbz,npz).ne.0) then
      write(lfnout,1003)
 1003 format(/,' WARNING: Inefficient distribution of boxes over ',
     + 'processors')
      endif
      endif
c
      if(.not.ma_push_get(mt_int,3*np,'lcnt',l_lcnt,i_lcnt))
     + call errquit('Failed to allocate memory for lcnt',0)
c
      call sp_prtcnt(lfnout,int_mb(i_lcnt))
c
      if(.not.ma_pop_stack(l_lcnt))
     + call errquit('Failed to deallocate lcnt',0)
c
      return
      end
      subroutine sp_prtcnt(lfnout,lcnt)
c
      implicit none
c
#include "sp_common.fh"
#include "msgids.fh"
#include "global.fh"
c
      integer lfnout,lcnt(np,3)
c
      integer i,j
c
      do 1 i=1,np
      lcnt(i,1)=0
      lcnt(i,2)=0
      lcnt(i,3)=0
    1 continue
c
      lcnt(me+1,1)=mbxloc
      lcnt(me+1,2)=nwmloc*nwa
      lcnt(me+1,3)=nsaloc
c
      if(np.gt.1) call ga_igop(msp_08,lcnt,3*np,'+')
c
      if(me.eq.0) then
      write(lfnout,1000)
 1000 format(/,' Initial distribution p:b(w+s)',/)
      write(lfnout,1001) (i-1,(lcnt(i,j),j=1,3),i=1,np)
 1001 format(4(3x,i4,':',i5,'(',i7,'+',i7,')'))
      endif
c
      return
      end
      subroutine sp_wtrest(lfn)
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
c
      integer lfn
c
      integer i_ltmp,l_ltmp
c
      if(.not.ma_push_get(mt_int,(mbbl+1)*mbb2,'ltmp',l_ltmp,i_ltmp))
     + call errquit('Failed to allocate memory for ltmp',0)
      call sp_wrest(lfn,int_mb(i_bb),int_mb(i_ltmp),mbbl+1)
      if(.not.ma_pop_stack(l_ltmp))
     + call errquit('Failed to deallocate ltmp',0)
c
      return
      end
      subroutine sp_wrest(lfn,lbbl,ltemp,mdim)
c
      implicit none
c
#include "sp_common.fh"
#include "msgids.fh"
#include "global.fh"
c
      integer lfn,mdim
      integer lbbl(mbbl,mbb2),ltemp(mdim,mbb2)
c
      integer i,j,k
c
      if(me.eq.0) then
      write(lfn,1000)
 1000 format('restart space')
      write(lfn,1001) np,mbbl
 1001 format(2i7)
      endif
c
      do 1 i=1,np
      if(i.eq.me+1) then
      do 2 k=1,mbb2
      ltemp(1,k)=0
      do 3 j=1,nbbl
      ltemp(j+1,k)=lbbl(j,k)
    3 continue
    2 continue
      ltemp(1,1)=nbbl
      else
      do 4 k=1,mbb2
      do 5 j=1,mdim
      ltemp(j,k)=0
    5 continue
    4 continue
      endif
c
      call ga_igop(msp_09,ltemp,mdim*mbb2,'+')
c
      if(me.eq.0) then
      write(lfn,1002) i-1,ltemp(1,1)
 1002 format(2i7)
      do 6 j=1,ltemp(1,1)
      write(lfn,1003) (ltemp(j+1,k),k=1,mbb2)
 1003 format(10i7)
    6 continue
      endif
c
    1 continue
c
      return
      end
      logical function sp_rdrest(lfn,fil)
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
c
      integer lfn
      character*255 fil
c
      logical sp_rrest
      external sp_rrest
c
      integer i_ltmp,l_ltmp
c
      if(.not.ma_push_get(mt_int,(mbbl+1)*mbb2,'ltmp',l_ltmp,i_ltmp))
     + call errquit('Failed to allocate memory for ltmp',0)
      sp_rdrest=sp_rrest(lfn,fil,int_mb(i_bb),int_mb(i_ltmp),mbbl+1)
      if(.not.ma_pop_stack(l_ltmp))
     + call errquit('Failed to deallocate ltmp',0)
c
      return
      end
      logical function sp_rrest(lfn,fil,lbbl,ltemp,mdim)
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      integer lfn,mdim
      character*255 fil
      integer lbbl(mbbl,mbb2),ltemp(mdim,mbb2)
c
      integer i,j,k,node,npp,nbl,mbblp,nbytes
      character*13 string
c
      if(me.eq.0) then
      open(unit=lfn,file=fil(1:index(fil,' ')-1),
     + status='old',form='formatted',err=9999)
      rewind(lfn)

    1 continue
      npp=0
      read(lfn,1000,end=9997) string
 1000 format(a13)
      if(string.ne.'restart space') goto 1
      read(lfn,1001) npp,mbblp
 1001 format(2i7)
      if(mbbl.ne.mbblp) npp=0
 9997 continue
      do 2 i=1,mbb2
      ltemp(1,i)=0
    2 continue
      endif
c
      nbytes=ma_sizeof(mt_int,1,mt_byte)
      call ga_brdcst(msp_10,npp,nbytes,0)
c
      if(np.ne.npp) goto 9998
c
      nbytes=mdim*mbb2*ma_sizeof(mt_int,1,mt_byte)
c
      do 3 i=1,np
c
      if(me.eq.0) then
      read(lfn,1002) node,nbl
 1002 format(2i7)
      ltemp(1,1)=node
      ltemp(1,2)=nbl
      do 4 j=1,nbl
      read(lfn,1003) (ltemp(j+1,k),k=1,mbb2)
 1003 format(10i7)
    4 continue
      endif
c
      call ga_brdcst(msp_11,ltemp,nbytes,0)
c
      if(ltemp(1,1).eq.me) then
      nbbl=ltemp(1,2)
      do 5 k=1,mbb2
      do 6 j=1,nbbl
      lbbl(j,k)=ltemp(j+1,k)
    6 continue
    5 continue
      endif
c
    3 continue
c
      sp_rrest=.true.
      return
 9998 continue
      sp_rrest=.false.
      return
 9999 continue
      call errquit('Failed to open restart file',0)
      return
      end
      subroutine sp_wrttrj(lfntrj,lxw,lvw,lxs,lvs,
     + stime,pres,temp,tempw,temps,iwl,xw,vw,xwcr,isl,xs,vs)
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      integer lfntrj
      logical lxw,lvw,lxs,lvs
      integer iwl(mwm,miw2),isl(msa,mis2)
      real*8 stime,pres,temp,tempw,temps
      real*8 xw(mwm,3,mwa),vw(mwm,3,mwa),xwcr(mwm,3)
      real*8 xs(msa,3),vs(msa,3)
c
      integer lenscr
c
      lenscr=ma_inquire_avail(mt_byte)/
     + ((6*mwa+3)*ma_sizeof(mt_dbl,1,mt_byte)+
     + (mis2+4)*ma_sizeof(mt_int,1,mt_byte))-1
      if(.not.ma_push_get(mt_dbl,lenscr*3*mwa,'bx',l_bx,i_bx))
     + call errquit('Failed to allocate bx',0)
      if(.not.ma_push_get(mt_dbl,lenscr*3*mwa,'bv',l_bv,i_bv))
     + call errquit('Failed to allocate bv',0)
      if(.not.ma_push_get(mt_dbl,lenscr*3,'br',l_br,i_br))
     + call errquit('Failed to allocate br',0)
      if(.not.ma_push_get(mt_int,lenscr*max(mis2,2),'bi',l_bi,i_bi))
     + call errquit('Failed to allocate bi',0)
      if(.not.ma_push_get(mt_int,lenscr,'n',l_n,i_n))
     + call errquit('Failed to allocate n',0)
c
      call sp_wttrj(lfntrj,lxw,lvw,lxs,lvs,
     + stime,pres,temp,tempw,temps,
     + iwl,xw,vw,xwcr,isl,xs,vs,int_mb(i_ipl),
     + lenscr,int_mb(i_bi),dbl_mb(i_bx),dbl_mb(i_bv),dbl_mb(i_br),
     + int_mb(i_bi),dbl_mb(i_bx),dbl_mb(i_bv))
c
      if(.not.ma_pop_stack(l_n))
     + call errquit('Failed to deallocate n',0)
      if(.not.ma_pop_stack(l_bi))
     + call errquit('Failed to deallocate bi',0)
      if(.not.ma_pop_stack(l_br))
     + call errquit('Failed to deallocate br',0)
      if(.not.ma_pop_stack(l_bv))
     + call errquit('Failed to deallocate bv',0)
      if(.not.ma_pop_stack(l_bx))
     + call errquit('Failed to deallocate bx',0)
c
      return
      end
      subroutine sp_wttrj(lfntrj,lxw,lvw,lxs,lvs,
     + stime,pres,temp,tempw,temps,iwl,xw,vw,xwcr,isl,xs,vs,
     + ipl,nb,ibw,bxw,bvw,brw,ibs,bxs,bvs)
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      integer lfntrj,nb
      logical lxw,lvw,lxs,lvs
      real*8 stime,pres,temp,tempw,temps
      integer iwl(mwm,miw2),isl(msa,mis2)
      real*8 xw(mwm,3,mwa),vw(mwm,3,mwa),xwcr(mwm,3)
      real*8 xs(msa,3),vs(msa,3)
      integer ipl(mbox,mip2),ibw(nb),ibs(nb)
      real*8 bxw(nb,3,mwa),bvw(nb,3,mwa),brw(nb,3)
      real*8 bxs(nb,3),bvs(nb,3)
c
      integer i,j,k,node,ncyc,icyc,numw,nums,number,nwmn,nsan
      integer ilp,ihp,jlp,jhp,ili,ihi,jli,jhi,ilw,ihw,jlw,jhw
      integer ils,ihs,jls,jhs
      character*10 rdate,rtime
c
      if(me.eq.0) then
c
      call swatch(rdate,rtime)
c
      write(lfntrj,1000)
 1000 format('frame')
      write(lfntrj,1001) stime,temp,pres,rdate,rtime
 1001 format(2f12.6,1pe12.5,1x,2a10)
      write(lfntrj,1002) ((vlat(i,j),j=1,3),i=1,3)
 1002 format(9f12.6)
      write(lfntrj,1003) lxw,lvw,lxs,lvs,nwm,nwa,nsa
 1003 format(4l1,3i10)
c
      if((lxw.or.lvw).and.nwm.gt.0) then
      number=0
      ncyc=nwm/nb+1
      numw=nb
      do 1 icyc=1,ncyc
      if(nwm-number.lt.numw) numw=nwm-number
      do 2 node=np-1,0,-1
      call ga_distribution(ga_ip,node,ilp,ihp,jlp,jhp)
      call ga_get(ga_ip,ilp,ihp,jlp,jhp,ipl,mbox)
      nwmn=ipl(1,2)
      if(nwmn.gt.0) then
      call ga_distribution(ga_iw,node,ili,ihi,jli,jhi)
      call ga_get(ga_iw,ili,ili+nwmn-1,jli,jli+lwdyn-1,iwl,mwm)
      call ga_distribution(ga_w,node,ilw,ihw,jlw,jhw)
      call ga_get(ga_w,ilw,ilw+nwmn-1,jlw,jlw+3*mwa-1,xw,mwm)
      if(lvw) 
     + call ga_get(ga_w,ilw,ilw+nwmn-1,jlw+3*mwa,jlw+6*mwa-1,vw,mwm)
      do 3 i=1,nwmn
      j=iwl(i,lwgmn)-number
      if(j.gt.0.and.j.le.numw) then
      do 4 k=1,nwa
      bxw(j,1,k)=xw(i,1,k)
      bxw(j,2,k)=xw(i,2,k)
      bxw(j,3,k)=xw(i,3,k)
    4 continue
      if(lvw) then
      do 5 k=1,nwa
      bvw(j,1,k)=vw(i,1,k)
      bvw(j,2,k)=vw(i,2,k)
      bvw(j,3,k)=vw(i,3,k)
    5 continue
      endif
      endif
    3 continue
      endif
    2 continue
      if(lvw) then
      do 6 i=1,numw
      write(lfntrj,1004) ((bxw(i,j,k),j=1,3),(bvw(i,j,k),j=1,3),k=1,nwa)
 1004 format(6f8.3)
    6 continue
      else
      do 7 i=1,numw
      write(lfntrj,1005) ((bxw(i,j,k),j=1,3),k=1,nwa)
 1005 format(3f8.3)
    7 continue
      endif
    1 continue
      endif
c
      if((lxs.or.lvs).and.nsa.gt.0) then
      number=0
      ncyc=nsa/nb+1
      nums=nb
      do 8 icyc=1,ncyc
      if(nsa-number.lt.nums) nums=nsa-number
      do 9 node=np-1,0,-1
      call ga_distribution(ga_ip,node,ilp,ihp,jlp,jhp)
      call ga_get(ga_ip,ilp,ihp,jlp,jhp,ipl,mbox)
      nsan=ipl(2,2)
      if(nsan.gt.0) then
      call ga_distribution(ga_is,node,ili,ihi,jli,jhi)
      call ga_get(ga_is,ili,ili+nsan-1,jli,jli+lsdyn-1,isl,msa)
      call ga_distribution(ga_s,node,ils,ihs,jls,jhs)
      call ga_get(ga_s,ils,ils+nsan-1,jls,jls+2,xs,msa)
      if(lvs) call ga_get(ga_s,ils,ils+nsan-1,jls+3,jls+5,vs,msa)
      do 10 i=1,nsan
      j=isl(i,lsgan)-number
      if(j.gt.0.and.j.le.nums) then
      bxs(j,1)=xs(i,1)
      bxs(j,2)=xs(i,2)
      bxs(j,3)=xs(i,3)
      if(lvs) then
      bvs(j,1)=vs(i,1)
      bvs(j,2)=vs(i,2)
      bvs(j,3)=vs(i,3)
      endif
      endif
   10 continue
      endif
    9 continue
      if(lvs) then
      do 11 i=1,nums
      write(lfntrj,1006) (bxs(i,j),j=1,3),(bvs(i,j),j=1,3)
 1006 format(6f8.3)
   11 continue
      else
      do 12 i=1,nums
      write(lfntrj,1007) (bxs(i,j),j=1,3)
 1007 format(3f8.3)
   12 continue
      endif
    8 continue
      endif
c
      endif
c
      return
      end
      subroutine sp_wrtmro(lfnmro,stime,pres,temp,tempw,temps,
     + iwl,xw,vw,xwcr,isl,xs,vs)
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      integer lfnmro
      integer iwl(mwm,miw2),isl(msa,mis2)
      real*8 stime,pres,temp,tempw,temps
      real*8 xw(mwm,3,mwa),vw(mwm,3,mwa),xwcr(mwm,3)
      real*8 xs(msa,3),vs(msa,3)
c
      call sp_wtmro(lfnmro,stime,pres,temp,tempw,temps,
     + iwl,xw,vw,xwcr,isl,xs,vs,int_mb(i_ipl))
c
      return
      end
      subroutine sp_wtmro(lfnmro,stime,pres,temp,tempw,temps,
     + iwl,xw,vw,xwcr,isl,xs,vs,ipl)
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      integer lfnmro
      real*8 stime,pres,temp,tempw,temps
      integer iwl(mwm,miw2),isl(msa,mis2)
      real*8 xw(mwm,3,mwa),vw(mwm,3,mwa),xwcr(mwm,3)
      real*8 xs(msa,3),vs(msa,3)
      integer ipl(mbox,mip2)
c
      integer j,k,l,nwmn,nsan,node,ilp,ihp,jlp,jhp,ili,ihi,jli,jhi
      integer ilw,ihw,jlw,jhw,ils,ihs,jls,jhs
c
      write(lfnmro) nwm,nwa,nsa,stime,temp,pres,vlat
c
      do 1 node=np-1,0,-1
      call ga_distribution(ga_ip,node,ilp,ihp,jlp,jhp)
      call ga_get(ga_ip,ilp,ihp,jlp,jhp,ipl,mbox)
      write(lfnmro) ((ipl(j,k),j=1,mbox),k=1,mip2)
      nwmn=ipl(1,2)
      nsan=ipl(2,2)
      if(nwmn.gt.0) then
      call ga_distribution(ga_iw,node,ili,ihi,jli,jhi)
      call ga_get(ga_iw,ili,ili+nwmn-1,jli,jhi,iwl,mwm)
      call ga_distribution(ga_w,node,ilw,ihw,jlw,jhw)
      call ga_get(ga_w,ilw,ilw+nwmn-1,jlw,jlw+3*mwa-1,xw,mwm) 
      call ga_get(ga_w,ilw,ilw+nwmn-1,jlw+3*mwa,jlw+6*mwa-1,vw,mwm)
      write(lfnmro) ((iwl(j,k),j=1,nwmn),k=1,miw2)
      write(lfnmro) (((xw(j,k,l),j=1,nwmn),k=1,3),l=1,nwa)
      write(lfnmro) (((vw(j,k,l),j=1,nwmn),k=1,3),l=1,nwa)
      endif
      if(nsan.gt.0) then
      call ga_distribution(ga_is,node,ili,ihi,jli,jhi)
      call ga_get(ga_is,ili,ili+nsan-1,jli,jhi,isl,msa)
      call ga_distribution(ga_s,node,ils,ihs,jls,jhs)
      call ga_get(ga_s,ils,ils+nsan-1,jls,jls+2,xs,msa)
      call ga_get(ga_s,ils,ils+nsan-1,jls+3,jls+5,vs,msa)
      write(lfnmro) ((isl(j,k),j=1,nsan),k=1,mis2)
      write(lfnmro) ((xs(j,k),j=1,nsan),k=1,3)
      write(lfnmro) ((vs(j,k),j=1,nsan),k=1,3)
      endif
    1 continue
c
      return
      end
      subroutine sp_rdmri(lfnmri,stime,pres,temp,tempw,temps,
     + iwl,xw,vw,xwcr,isl,xs,vs)
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      integer lfnmri
      integer iwl(mwm,miw2),isl(msa,mis2)
      real*8 stime,pres,temp,tempw,temps
      real*8 xw(mwm,3,mwa),vw(mwm,3,mwa),xwcr(mwm,3)
      real*8 xs(msa,3),vs(msa,3)
c
      call sp_rmri(lfnmri,stime,pres,temp,tempw,temps,
     + iwl,xw,vw,xwcr,isl,xs,vs,int_mb(i_ipl))
c
      return
      end
      subroutine sp_rmri(lfnmri,stime,pres,temp,tempw,temps,
     + iwl,xw,vw,xwcr,isl,xs,vs,ipl)
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "global.fh"
c
      integer lfnmri
      real*8 stime,pres,temp,tempw,temps
      integer iwl(mwm,miw2),isl(msa,mis2)
      real*8 xw(mwm,3,mwa),vw(mwm,3,mwa),xwcr(mwm,3)
      real*8 xs(msa,3),vs(msa,3)
      integer ipl(mbox,mip2)
c
      integer i,j,k,l,nwmn,nsan,node,ilp,ihp,jlp,jhp,ili,ihi,jli,jhi
      integer ilw,ihw,jlw,jhw,ils,ihs,jls,jhs
      integer ltemp(3)
      real*8 rtemp(12)
c
      if(me.eq.0) then
      read(lfnmri) ltemp,rtemp
      do 1 node=np-1,0,-1
      read(lfnmri) ((ipl(j,k),j=1,mbox),k=1,mip2)
      nwmn=ipl(1,2)
      nsan=ipl(2,2)
      call ga_distribution(ga_ip,node,ilp,ihp,jlp,jhp)
      call ga_put(ga_ip,ilp,ihp,jlp,jhp,ipl,mbox)
      if(nwmn.gt.0) then
      read(lfnmri) ((iwl(j,k),j=1,nwmn),k=1,miw2)
      read(lfnmri) (((xw(j,k,l),j=1,nwmn),k=1,3),l=1,nwa)
      read(lfnmri) (((vw(j,k,l),j=1,nwmn),k=1,3),l=1,nwa)
      call ga_distribution(ga_iw,node,ili,ihi,jli,jhi)
      call ga_put(ga_iw,ili,ili+nwmn-1,jli,jhi,iwl,mwm)
      call ga_distribution(ga_w,node,ilw,ihw,jlw,jhw)
      call ga_put(ga_w,ilw,ilw+nwmn-1,jlw,jlw+3*mwa-1,xw,mwm) 
      call ga_put(ga_w,ilw,ilw+nwmn-1,jlw+3*mwa,jlw+6*mwa-1,vw,mwm)
      endif
      if(nsan.gt.0) then
      read(lfnmri) ((isl(j,k),j=1,nsan),k=1,mis2)
      read(lfnmri) ((xs(j,k),j=1,nsan),k=1,3)
      read(lfnmri) ((vs(j,k),j=1,nsan),k=1,3)
      call ga_distribution(ga_is,node,ili,ihi,jli,jhi)
      call ga_put(ga_is,ili,ili+nsan-1,jli,jhi,isl,msa)
      call ga_distribution(ga_s,node,ils,ihs,jls,jhs)
      call ga_put(ga_s,ils,ils+nsan-1,jls,jls+2,xs,msa)
      call ga_put(ga_s,ils,ils+nsan-1,jls+3,jls+5,vs,msa)
      endif
    1 continue
      endif
c
      call ga_brdcst(msp_12,ltemp,3*ma_sizeof(mt_int,1,mt_byte),0)
      call ga_brdcst(msp_13,rtemp,12*ma_sizeof(mt_dbl,1,mt_byte),0)
      nwm=ltemp(1)
      nwa=ltemp(2)
      nsa=ltemp(3)
      stime=rtemp(1)
      temp=rtemp(2)
      pres=rtemp(3)
      k=3
      do 2 i=1,3
      do 3 j=1,3
      k=k+1
      vlat(i,j)=rtemp(k)
    3 continue
    2 continue
      call ga_distribution(ga_ip,me,ilp,ihp,jlp,jhp)
      call ga_get(ga_ip,ilp,ihp,jlp,jhp,ipl,mbox)
      nwmloc=ipl(1,2)
      nsaloc=ipl(2,2)
      if(nwmloc.gt.0) then
      call ga_distribution(ga_iw,me,ili,ihi,jli,jhi)
      call ga_get(ga_iw,ili,ili+nwmloc-1,jli,jhi,iwl,mwm)
      call ga_distribution(ga_w,me,ilw,ihw,jlw,jhw)
      call ga_get(ga_w,ilw,ilw+nwmloc-1,jlw,jlw+3*mwa-1,xw,mwm) 
      call ga_get(ga_w,ilw,ilw+nwmloc-1,jlw+3*mwa,jlw+6*mwa-1,vw,mwm)
      endif
      if(nsaloc.gt.0) then
      call ga_distribution(ga_is,me,ili,ihi,jli,jhi)
      call ga_get(ga_is,ili,ili+nsaloc-1,jli,jhi,isl,msa)
      call ga_distribution(ga_s,me,ils,ihs,jls,jhs)
      call ga_get(ga_s,ils,ils+nsaloc-1,jls,jls+2,xs,msa)
      call ga_get(ga_s,ils,ils+nsaloc-1,jls+3,jls+5,vs,msa)
      endif
c
      return
      end
      subroutine sp_fix()
c
      implicit none
c
      return
      end
